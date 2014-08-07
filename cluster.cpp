/*
 * Compile command : 
 * g++ -I/home/mridul/seqan/seqan-trunk/core/include -I/home/mridul/seqan/seqan-trunk/extras/include/ -lrt -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -std=c++0x match_reads.cpp
 */

#include <iostream>
#include <string>
#include <cstdint>
#include <fstream>
#include <unordered_map>
#include <vector>

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
using namespace seqan;

#define BOUND 5
#define THRESHOLD 10
using namespace std;

typedef String<char> TSequence;                 // sequence type
typedef Align<TSequence,ArrayGaps> TAlign;      // align type

class llist {
    public:
    llist* next;
    char side;
    int pos;
    int entrynum;
    llist* cluster;
    llist() {
        side='c';
    }
};

#define RANK(x) (2*x->entrynum + (x->side == 'l' ? 0 : 1) + 1) // ranks are 1, 2, ...
typedef unordered_map<int64_t,llist* > HashTable;


// ------------ disjoint sets begin -----------------

#include <boost/pending/disjoint_sets.hpp>
#include <boost/pending/property.hpp>

class Element
{
public:
    explicit
    Element(int n) : mSomeInt(n) { }

    int someInt() const { return mSomeInt; }

    // The following class members are specific for the disjoint_sets
    // implementation
    size_t dsID;
    size_t dsRank;
    size_t dsParent;

private:
    int mSomeInt;
};

inline bool
operator==(Element const& lhs, Element const& rhs)
{
    return lhs.someInt() == rhs.someInt();
}

inline bool
operator!=(Element const& lhs, Element const& rhs)
{
    return ! operator==(lhs, rhs);
}

class Parent
{
public:
    Parent(std::vector<Element>& e) : mElements(e) { }
    std::vector<Element>& mElements;
};

class Rank
{
public:
    Rank(std::vector<Element>& e) : mElements(e) { }
    std::vector<Element>& mElements;
};

namespace boost {
template <>
struct property_traits<Rank*>
{
    typedef size_t value_type;
};
}

inline Element const&
get(Parent* pa, Element const& k)
{
    return pa->mElements.at(k.dsParent);
}

inline void
put(Parent* pa, Element k, Element const& val)
{
    pa->mElements.at(k.dsID).dsParent = val.dsID;
}

inline size_t const&
get(Rank*, Element const& k)
{
    return k.dsRank;
}

inline void
put(Rank* pa, Element k, size_t const& val)
{
    pa->mElements.at(k.dsID).dsRank = val;
}

void
printElements(std::vector<Element>& elements)
{
    std::cout << "Elements:            ";
    for (size_t i = 0; i < elements.size(); ++i)
    {
        std::cout << std::setw(4) << elements[i].someInt();
    }
    std::cout << std::endl;
    std::cout << "Set representatives: ";
    for (size_t i = 0; i < elements.size(); ++i)
    {
        std::cout << std::setw(4) << elements[i].dsParent;
    }
    std::cout << std::endl;
    std::cout << "ID                 : ";
    for (size_t i = 0; i < elements.size(); ++i)
    {
        std::cout << std::setw(4) << elements[i].dsID;
    }
    std::cout << std::endl;
}

inline bool
compareByParent(Element const& lhs, Element const& rhs)
{
    return lhs.dsParent < rhs.dsParent;
}

inline bool
compareBySomeInt(Element const& lhs, Element const& rhs)
{
    return lhs.someInt() < rhs.someInt();
}

typedef boost::disjoint_sets<Rank*, Parent*> DisjointSets;

// ------------ disjoint sets end -------------


vector<string> left_reads,right_reads;

bool match_reads(llist* x,llist* y) {
    
    TSequence seql1 = left_reads[x->entrynum];
    TSequence seql2 = left_reads[y->entrynum];
    TSequence seqr1 = right_reads[x->entrynum];
    TSequence seqr2 = right_reads[y->entrynum];

    TAlign alignl;
    resize(rows(alignl), 2);
    assignSource(row(alignl,0),seql1);
    assignSource(row(alignl,1),seql2);
    int scorel = globalAlignment(alignl, Score<int,Simple>(0,-1,-1),-BOUND,BOUND);
    cout << scorel << endl;
    cout << alignl << endl;

    TAlign alignr;
    resize(rows(alignr), 2);
    assignSource(row(alignr,0),seqr1);
    assignSource(row(alignr,1),seqr2);
    int scorer = globalAlignment(alignr, Score<int,Simple>(0,-1,-1),-BOUND,BOUND);
    cout << scorer << endl;
    cout << alignr << endl;
    if(scorel<=THRESHOLD && scorer <=THRESHOLD) {
        return true;
    }

    return false;

}


void bin_reads(const HashTable &hashtab, const std::vector<Element> &elements, DisjointSets &ds) {

    for ( auto it = hashtab.begin(); it != hashtab.end(); ++it ) {
        cout << it->first << ":";
        for (auto x = it->second; (x != NULL) && (x->next != NULL); x = x->next) {
          int xrank = RANK(x);
          for (auto y = x->next; y != NULL; y = y->next) {
            int yrank = RANK(y);
            if (match_reads(x, y)) {
              ds.union_set(elements[xrank], elements[yrank]);
            }
          }
        }
    }
}


int main(int argc, char*  argv[]) {
    int k,skipN=0,num=0;
    cin >> k; //k-mer size // not greater than 32
    int64_t slid,mask=(1<<(2*k))-1;
    HashTable hashtab;
    ifstream f1(argv[1]),f2(argv[2]);
    string p1,p2;
    int count=0, len;

    while(getline(f1,p1) && getline(f2,p2)) {
        //f1.push_back(p1);
        //f2.push_back(p2);
        len=p1.length();
        if(count==1) {
            left_reads.push_back(p1);
            right_reads.push_back(p2);
            slid=0;
            for(int i=0;i<len-k+1;i++) {
                if(skipN>0) {
                    skipN--;
                    continue;
                }
                if(p1[i]=='A') {
                    slid=slid<<2;
                }
                else if(p1[i]=='C') {
                    slid=slid<<2;
                    slid+=1;
                }
                else if(p1[i]=='T') {
                    slid=slid<<2;
                    slid+=2;
                }
                else if (p1[i]=='G') {
                    slid=slid<<2;
                    slid+=3;
                }
                else {
                    skipN=k;
                    continue;
                }

                if(i>=k-1) {
                    slid&=mask;
                    llist* info = new llist();
                    info->side='l';
                    info->pos=i-k+1;
                    info->entrynum=num;
                    info->cluster = info;
                    llist** point = &hashtab[slid];
                    if(*point==NULL) {
                        *point=info;
                        info->next=NULL;
                    }
                    else {
                        info->next=*point;
                        hashtab[slid]=info;
                        //point->next=info;
                    }
                }
            }

            slid=0;
            skipN=0;
            len=p2.length();
            for(int i=0;i<len-k+1;i++) {
                if(skipN>0) {
                    skipN--;
                    continue;
                }
                if(p2[i]=='A') {
                    slid=slid<<2;
                }
                else if(p2[i]=='C') {
                    slid=slid<<2;
                    slid+=1;
                }
                else if(p2[i]=='T') {
                    slid=slid<<2;
                    slid+=2;
                }
                else if(p2[i]=='G') {
                    slid=slid<<2;
                    slid+=3;
                }
                else {
                    skipN=k;
                    continue;
                }
                if(i>=k-1) {
                    slid&=mask;
                    llist* info = new llist();
                    info->side='r';
                    info->pos=i-k+1;
                    info->entrynum=num;
                    info->cluster = info;
                    llist** point = &hashtab[slid];
                    if(*point==NULL) {
                        *point=info;
                        info->next=NULL;
                    }
                    else {
                        info->next=*point;
                        hashtab[slid]=info;
                    }
                }
            }
            num++;
        }

        count+=1;
        count%=4;
    }
    // HASH TABLE CREATED
    llist* temp;
    for ( auto it = hashtab.begin(); it != hashtab.end(); ++it ) {
        cout << it->first << ":";
        temp=it->second;
        while(temp!=NULL) {
            cout <<"--> [ " << temp->entrynum << ","<< temp->side << "," << temp->pos <<" ] " ;
            temp=temp->next;
        }
        cout << "\n";
    }

    for ( auto it = hashtab.begin(); it != hashtab.end(); ++it ) {
        temp=it->second;
        while(temp!=NULL) {

            temp=temp->next;
        }
    }

    
    std::vector<Element> elements;
    elements.reserve(2*num);
    for (size_t i = 0; i < elements.capacity(); ++i)
    {
        elements.push_back(Element(i+1));
    }

    for (size_t i = 0; i < elements.size(); ++i)
    {
        elements[i].dsID = i+1;
    }

    Rank rank(elements);
    Parent parent(elements);

    DisjointSets ds(&rank, &parent);

    for (size_t i = 0; i < elements.size(); ++i)
    {
        ds.make_set(elements.at(i));
    }

    bin_reads(hashtab, elements, ds);

    std::cout << "Found " << ds.count_sets(elements.begin(), elements.end()) << " sets:" << std::endl;
    printElements(elements);

    ds.compress_sets(elements.begin(), elements.end());

    std::cout << std::endl << "After path compression:" << std::endl;
    printElements(elements);

    std::sort(elements.begin(), elements.end(), compareByParent);

    std::cout << std::endl << "After path compression and sorting by parent:" << std::endl;
    printElements(elements);

    std::cout << std::endl << "Now we can iterate through all elements of each set separately using the indices:" << std::endl;
    {
        size_t first = 0;
        while (first < elements.size())
        {
            size_t currentParent = elements.at(first).dsParent;
            size_t last = first;
            while (last < elements.size() && elements.at(last).dsParent == currentParent)
            {
                ++last;
            }
            std::cout << "\tRange: [" << first << "," << last << "). Sorted elements: ";
            for (size_t i = first; i < last; ++i)
            {
                std::cout << elements.at(i).someInt() << " ";
            }
            std::cout << std::endl;
            first = last;
        }
    }
    

    return 0;
}
