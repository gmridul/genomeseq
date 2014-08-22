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

#include <boost/pending/disjoint_sets.hpp>
#include <boost/unordered/unordered_set.hpp>


#define BOUND 2
#define THRESHOLD -10
using namespace std;


typedef seqan::String<char> TSequence;                 // sequence type
typedef seqan::Align<TSequence,seqan::ArrayGaps> TAlign;      // align type


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

vector<string> left_reads,right_reads;

bool match_reads(llist* x,llist* y) {
    
    TSequence seql1 = left_reads[x->entrynum];
    TSequence seql2 = left_reads[y->entrynum];
    TSequence seqr1 = right_reads[x->entrynum];
    TSequence seqr2 = right_reads[y->entrynum];

    TAlign alignl;
    seqan::resize(rows(alignl), 2);
    seqan::assignSource(row(alignl,0),seql1);
    seqan::assignSource(row(alignl,1),seql2);
    int scorel = seqan::globalAlignment(alignl, seqan::Score<int,seqan::Simple>(0,-1,-1),-BOUND,BOUND);
    cout << "left score " << scorel << endl;
    cout << alignl << endl;

    if(scorel < THRESHOLD) return false;

    TAlign alignr;
    seqan::resize(rows(alignr), 2);
    seqan::assignSource(row(alignr,0),seqr1);
    seqan::assignSource(row(alignr,1),seqr2);
    int scorer = seqan::globalAlignment(alignr, seqan::Score<int,seqan::Simple>(0,-1,-1),-BOUND,BOUND);
    cout <<" right score " <<  scorer << endl;
    cout << alignr << endl;

    if(scorer < THRESHOLD) return false;

    return true;

}


void bin_reads(const HashTable &hashtab, DisjointSets &ds) {

    for ( auto it = hashtab.begin(); it != hashtab.end(); ++it ) {
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
            if (x->entrynum!=y->entrynum && match_reads(x, y)) {
              ds.union_set(xrank, yrank);
            }
          }
        }
    }
}

extern "C" {
#include "clustal-omega.h"
}
/**************************************************************************************************************************************************************/
mseq_t *prMSeq;
int myAlign(mseq_t *prMSeq, mseq_t *prMSeqProfile, opts_t *prOpts) {
   
    /* HMM
     */
    /* structs with pseudocounts etc; one for each HMM infile, i.e.
     * index range: 0..iHMMInputFiles */
    hmm_light *prHMMs = NULL;

    /* MSA order in which nodes/profiles are to be merged/aligned
       (order of nodes in guide tree (left/right)*/
    int *piOrderLR = NULL;

    /* weights per sequence */
    double *pdSeqWeights = NULL;

    /* Iteration
     */
    int iIterationCounter = 0;
    double dAlnScore;
    /* last dAlnScore for iteration */
    double dLastAlnScore = -666.666;

    

    assert(NULL != prMSeq);
    if (NULL != prMSeqProfile) {
        assert(TRUE == prMSeqProfile->aligned);
    }



/*    if (TRUE == prOpts->bPileup){
        PileUp(prMSeq, prOpts->rHhalignPara, prOpts->iClustersizes);
        return 0;
    }
     */

#if 0
    Log(&rLog, LOG_WARN, "Using a sequential alignment order.");
    SequentialAlignmentOrder(&piOrderLR, prMSeq->nseqs);
#else

    if (OK != AlignmentOrder(&piOrderLR, &pdSeqWeights, prMSeq,
                prOpts->iPairDistType,
                prOpts->pcDistmatInfile, prOpts->pcDistmatOutfile,
                prOpts->iClusteringType, prOpts->iClustersizes, 
                prOpts->pcGuidetreeInfile, prOpts->pcGuidetreeOutfile, prOpts->pcClustfile, 
                prOpts->bUseMbed, prOpts->bPercID)) {
        Log(&rLog, LOG_ERROR, "AlignmentOrder() failed. Cannot continue");
        return -1;
    }
#endif

    /* if max-hmm-iter is set < 0 then do not perform alignment 
     * there is a problem/feature(?) that the un-aligned sequences are output 
     */
    if (prOpts->iMaxHMMIterations < 0){
        Log(&rLog, LOG_VERBOSE,
            "iMaxHMMIterations < 0 (%d), will not perform alignment", prOpts->iMaxHMMIterations);
        return 0;
    }


    /* Progressive alignment of input sequences. Order defined by
     * branching of guide tree (piOrderLR). Use optional
     * background HMM information (prHMMs[0..prOpts->iHMMInputFiles-1])
     *
     */
    dAlnScore = HHalignWrapper(prMSeq, piOrderLR, pdSeqWeights,
                               2*prMSeq->nseqs -1/* nodes */,
                               prHMMs, prOpts->iHMMInputFiles, -1, prOpts->rHhalignPara);
    dLastAlnScore = dAlnScore;
    Log(&rLog, LOG_VERBOSE,
        "Alignment score for first alignment = %f", dAlnScore);        

    /* ------------------------------------------------------------
     * prMSeq is aligned now. Now start iterations if requested and save the
     * alignment at the very end.
     * ------------------------------------------------------------ */

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
/**********************************************************************************************************************************************************/
void tree_to_align(char* t) {
/* the multiple sequence structure */

    /* for openmp: number of threads to use */
    int iThreads = 0;
    /* alignment options to use */
    opts_t rAlnOpts;
    int iAux;
    LogDefaultSetup(&rLog);
    SetDefaultAlnOpts(&rAlnOpts);
    rAlnOpts.pcGuidetreeInfile = t;
    if(myAlign(prMSeq, (mseq_t *)NULL, &rAlnOpts)) {
        Log(&rLog, LOG_FATAL, "A fatal error happended during the alignment process");
    }
    if (WriteAlignment(prMSeq, NULL, MSAFILE_A2M, 1000, false)) {
        Log(&rLog, LOG_FATAL, "Could not save alignment");
    } 
    FreeMSeq(&prMSeq);

    Log(&rLog, LOG_INFO, "Successfull program exit");
}


int main(int argc, char*  argv[]) {
    int no_of_seq;
    cin >> no_of_seq;
    string fname;
    cin >> fname;
    char * writable = new char[fname.size() + 1];
    std::copy(fname.begin(), fname.end(), writable);
    writable[fname.size()] = '\0'; // don't forget the terminating 0
    prMSeq->seqtype  = SEQTYPE_DNA;
    prMSeq->aligned  = false;
    prMSeq->filename = writable;
    prMSeq->seqtype =SQFILE_FASTA; 
    prMSeq->nseqs    = 0;
    prMSeq->seq =  (char **) CKREALLOC(prMSeq->seq, (prMSeq->nseqs+1) * sizeof(char *));
    ReadSequences(prMSeq, writable, SEQTYPE_DNA, SQFILE_FASTA, false,false, 10000, 300);
        
    /*
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
    
    VecInt rank (10);
    VecInt parent (10);
    DisjointSets ds(&rank[0], &parent[0]);

    bin_reads(hashtab, ds);
    */

    return 0;
}
