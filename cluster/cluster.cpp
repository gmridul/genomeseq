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
//using namespace seqan;

#include "dswrapper.hpp"

#define BOUND 2
#define THRESHOLD -10
using namespace std;

typedef seqan::String<char> TSequence;                     // sequence type
typedef seqan::Align<TSequence,seqan::ArrayGaps> TAlign;   // align type

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

vector<string> left_reads,right_reads;

int match_seqs(const string &s1, const string &s2) {
    TSequence seq1 = s1, seq2 = s2;
    TAlign align;
    seqan::resize(rows(align), 2);
    seqan::assignSource(row(align, 0), seq1);
    seqan::assignSource(row(align, 1), seq2);
    int score = seqan::globalAlignment(align, seqan::Score<int,seqan::Simple>(0,-1,-1), -BOUND, BOUND);
    cout << "score " << score << endl;
    cout << align << endl;
    return score;
}

bool match_reads(llist* x,llist* y) {
    if (match_seqs(left_reads[x->entrynum], left_reads[y->entrynum]) < THRESHOLD) return false;
    if (match_seqs(right_reads[x->entrynum], right_reads[y->entrynum]) < THRESHOLD) return false;
    return true;
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

extern "C" {
#include "clustal-omega.h"
}

/**************************************************************************************************************************************************************/

int myAlign(mseq_t *prMSeq, opts_t *prOpts) {

    /* HMM
     * structs with pseudocounts etc; one for each HMM infile, i.e.
     * index range: 0..iHMMInputFiles */
    hmm_light *prHMMs = NULL;

    /* MSA order in which nodes/profiles are to be merged/aligned
       (order of nodes in guide tree (left/right)*/
    int *piOrderLR = NULL;

    /* weights per sequence */
    double *pdSeqWeights = NULL;

    double dAlnScore;

    assert(NULL != prMSeq);
    /*
       if (TRUE == prOpts->bPileup){
       PileUp(prMSeq, prOpts->rHhalignPara, prOpts->iClustersizes);
       return 0;
       }
       */
    if (OK != AlignmentOrder(&piOrderLR, &pdSeqWeights, prMSeq,
                prOpts->iPairDistType,
                prOpts->pcDistmatInfile, prOpts->pcDistmatOutfile,
                prOpts->iClusteringType, prOpts->iClustersizes, 
                prOpts->pcGuidetreeInfile, prOpts->pcGuidetreeOutfile, prOpts->pcClustfile, 
                prOpts->bUseMbed, prOpts->bPercID)) {
        Log(&rLog, LOG_ERROR, "AlignmentOrder() failed. Cannot continue");
        return -1;
    }

    /* Progressive alignment of input sequences. Order defined by
     * branching of guide tree (piOrderLR). Use optional
     * background HMM information (prHMMs[0..prOpts->iHMMInputFiles-1])
     *
     */
    dAlnScore = HHalignWrapper(prMSeq, piOrderLR, pdSeqWeights,
            2*prMSeq->nseqs -1/* nodes */,
            prHMMs, prOpts->iHMMInputFiles, -1, prOpts->rHhalignPara);
    Log(&rLog, LOG_VERBOSE,
            "Alignment score for first alignment = %f", dAlnScore);        

    if (NULL != piOrderLR) {
        free(piOrderLR);
        piOrderLR=NULL;
    }
    if (NULL != pdSeqWeights) {
        free(pdSeqWeights);
        pdSeqWeights=NULL;
    }

    return 0;
}
/*************************************************************************************************************************/
void tree_to_align(char* t) {
    int no_of_seq=4;
    string fname,tname;
    fname="read_file";tname="guide_tree_example";
    //cin >> fname >> tname;
    cout << fname<<tname <<endl;
    string tree;
    char * writable = new char[fname.size() + 1];
    std::copy(fname.begin(), fname.end(), writable);
    writable[fname.size()] = '\0'; // don't forget the terminating 0
    mseq_t *prMSeq = (mseq_t*)CKCALLOC(1,sizeof(mseq_t));
    prMSeq->seqtype  = SEQTYPE_DNA;
    prMSeq->aligned  = false;
    prMSeq->filename = (char *)fname.c_str();
    prMSeq->seqtype =SQFILE_FASTA; 
    prMSeq->nseqs    = 0;
    prMSeq->seq =  (char **) CKREALLOC(prMSeq->seq, (prMSeq->nseqs+1) * sizeof(char *));
    ReadSequences(prMSeq, writable, SEQTYPE_DNA, SQFILE_FASTA, false,false, 10000, 300);
    for(int i=0;i<prMSeq->nseqs;i++) {
        printf("%s\n",prMSeq->seq[i]);
    }

    /* alignment options to use */
    opts_t rAlnOpts;
    LogDefaultSetup(&rLog);
    SetDefaultAlnOpts(&rAlnOpts);
    rAlnOpts.pcGuidetreeInfile = t;

    if(myAlign(prMSeq, &rAlnOpts)) {
        Log(&rLog, LOG_FATAL, "A fatal error happended during the alignment process");
    }

    if (WriteAlignment(prMSeq, NULL, MSAFILE_A2M, 1000, TRUE)) {
        Log(&rLog, LOG_FATAL, "Could not save alignment");
    } 
    FreeMSeq(&prMSeq);

}


int main(int argc, char*  argv[]) {
    int no_of_seq;
    //    cin >> no_of_seq;

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
