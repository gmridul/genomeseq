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

#include <boost/pending/disjoint_sets.hpp>
#include <boost/unordered/unordered_set.hpp>


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

typedef std::vector<int> VecInt;
typedef boost::disjoint_sets<int*,int*> DisjointSets;


bool match_reads(llist* x,llist* y) {//TSequence s1,TSequence s2,int k,int thold) {
    
     TSequence seql1 =left[x->entrynum];
    TSequence seql2 =left[y->entrynum];
    TSequence seqr1 =right[x->entrynum];
    TSequence seqr2 =right[y->entrynum];

    TAlign alignl;
    resize(rows(alignl), 2);
    assignSource(row(alignl,0),seql1);
    assignSource(row(alignl,1),seql2);
    int scorel = globalAlignment(alignl, Score<int,Simple>(0,-1,-1),-BOUND,BOUND);

    TAlign alignr;
    resize(rows(alignr), 2);
    assignSource(row(alignr,0),seqr1);
    assignSource(row(alignr,1),seqr2);
    int scorer = globalAlignment(alignr, Score<int,Simple>(0,-1,-1),-BOUND,BOUND);
    cout << score << endl;
    cout << align << endl;
    if(scorel<=THRESHOLD && scorer <=THRESHOLD) {
        return true;
    }

    return false;

}


void bin_reads(const HashTable &hashtab, DisjointSets &ds) {

    for ( auto it = hashtab.begin(); it != hashtab.end(); ++it ) {
        cout << it->first << ":";
        for (auto x = it->second; x != NULL; x = x->next) {
          int rank = RANK(x);
          if (!ds.find_set(rank)) { // find_set returns 0 if not found
            ds.make_set(rank); 
          }
        }
        for (auto x = it->second; (x != NULL) && (x->next != NULL); x = x->next) {
          int xrank = RANK(x);
          for (auto y = x->next; y != NULL; y = y->next) {
            int yrank = RANK(y);
            if (match_reads(x, y)) {
              ds.union_set(xrank, yrank);
            }
          }
        }
    }
}


int main(int argc, char*  argv[]) {
    int k,skipN=0,num=1;
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
        }

        count+=1;
        count%=4;
        num++;
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

    
    VecInt rank (10);
    VecInt parent (10);
    DisjointSets ds(&rank[0], &parent[0]);
    

    return 0;
}
