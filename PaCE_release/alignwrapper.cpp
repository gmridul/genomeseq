#include "alignwrapper.hpp"
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cstring>

#define MIN_SCORE -2147483648
void print_sequences(mseq_t *prMSeq) {
    for(int i=0;i<prMSeq->nseqs;i++) {
        printf("%s\n",prMSeq->seq[i]);
    }
}

void setup_sequences(mseq_t *prMSeq, const Reads& reads, int num_leaves, uint *leaf_ids) {
    prMSeq->seqtype  = SEQTYPE_DNA;
    prMSeq->aligned  = false;
    prMSeq->filename = NULL;
    prMSeq->nseqs    = num_leaves;
    prMSeq->seq =  (char **) CKMALLOC((prMSeq->nseqs) * sizeof(char *));
    prMSeq->orig_seq =  (char **) CKMALLOC((prMSeq->nseqs) * sizeof(char *));
    prMSeq->sqinfo = (SQINFO *)CKREALLOC(prMSeq->sqinfo, (prMSeq->nseqs+1) * sizeof(SQINFO));
    for (int i=0; i<num_leaves; ++i) {
        uint read_id = leaf_ids[i];
        prMSeq->seq[i] = CkStrdup(reads[read_id].seq.c_str());;
        prMSeq->orig_seq[i] = CkStrdup(prMSeq->seq[i]);;
        prMSeq->sqinfo[i].flags = SQINFO_NAME | SQINFO_ID;;
        sprintf(prMSeq->sqinfo[i].name, "%s", reads[read_id].name.c_str());
        sprintf(prMSeq->sqinfo[i].id, "%d", i);
    }
}

std::ofstream dataf;
int serial_num;
void process_alignment(int num_reads, char **aligned_reads) {

    for(int i=0;i<num_reads;i++) {
        printf("%s\n",aligned_reads[i]);
    }
    // do other stuff
    int count[5]={0,0,0,0,0};
    int align_read_len = std::strlen(aligned_reads[0]);
    dataf.open("cluster.dat", /*std::ofstream::out | std::ofstream::app*/std::ios_base::app);
    for(int i=0;i<align_read_len;i++) {
        for(int j=0;j<num_reads;j++) {
            if(aligned_reads[j][i]=='A') count[0]++;
            else if(aligned_reads[j][i]=='C') count[1]++;
            else if(aligned_reads[j][i]=='T') count[2]++;
            else if(aligned_reads[j][i]=='G') count[3]++;
            else if(aligned_reads[j][i]=='-') count[4]++;
        }
        std::sort(count,count+5);
        dataf << serial_num++ << " ";
        //std::cout << count[4] << " " << count[3] << " "<< count[2] << " " << count[1] << " "<<count[0] << "\n";
        for(int j=3;j>=0;j--) {
            dataf << (100*((float)(count[j+1]-count[j])))/num_reads << " ";
        }
        dataf << "\n";
        for(int j=0;j<5;j++) count[j]=0;
    }
    
    //std::stringstream ss;
    //ss << cluster_num;
    dataf <<serial_num++ <<  " 0 0 0 0\n";
    dataf.close();
}


void align_cluster(const Reads &reads, int num_leaves, uint *left, uint *right, uint *leaf_ids) {
    LogDefaultSetup(&rLog);
    mseq_t *prMSeq = (mseq_t*)CKCALLOC(1,sizeof(mseq_t));
    setup_sequences(prMSeq, reads, num_leaves, leaf_ids);
    print_sequences(prMSeq);

    float* leftLength  = (float *)malloc(num_leaves*sizeof(float));
    float* rightLength = (float *)malloc(num_leaves*sizeof(float));
    char** leafNames   = (char **)malloc(num_leaves*sizeof(char*));
    for (int i=0; i<num_leaves; i++) {
        leaf_ids[i] = i;
        leftLength[i] = rightLength[i] = 1.0;
        leafNames[i] = prMSeq->sqinfo[i].name;
    }

    tree_t *tree = (tree_t *)malloc(sizeof(tree_t));
    // root node is always numbered 0
    MuscleTreeCreate(tree, num_leaves, 0, left, right, leftLength, rightLength, leaf_ids, leafNames);
    //MuscleTreeToFile(stdout, tree);
    int *piOrderLR = NULL;
    TraverseTree(&piOrderLR, tree, prMSeq);
    hhalign_para rHhalignPara;
    SetDefaultHhalignPara(&rHhalignPara);
    double dAlnScore = HHalignWrapper(prMSeq, piOrderLR, NULL/*pdSeqWeights*/,
            2*num_leaves-1/* nodes */, NULL/*prHMMs*/, 0/*iHMMInputFiles*/, -1, rHhalignPara);
    //printf("Align score = %f\n", dAlnScore);
    //process_alignment(prMSeq->nseqs, prMSeq->seq);
    FreeMuscleTree(tree);
    FreeMSeq(&prMSeq);
    if (NULL != piOrderLR) CkFree(piOrderLR, __FUNCTION__, __LINE__);
    free(leftLength);
    free(rightLength);
    free(leafNames);
}
