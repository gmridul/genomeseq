#include "alignwrapper.hpp"
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cstring>

#define MIN_SCORE -2147483648

int match_seqs(const std::string &s1, const std::string &s2) {
    //std::cout << "aligning " << s1 << " and " << s2 << std::endl;
    std::cout.flush();
    TSequence seq1 = s1, seq2 = s2;
    TAlign align;
    seqan::resize(rows(align), 2);
    seqan::assignSource(row(align, 0), seq1);
    seqan::assignSource(row(align, 1), seq2);
    int score = seqan::globalAlignment(align, seqan::Score<int,seqan::Simple>(0,-1,-1), -BOUND, BOUND);
    //if (score > MIN_SCORE) {
    //    std::cout << "score " << score << std::endl;
    //    //std::cout << align << std::endl;
    //}
    return score;
}

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
        prMSeq->seq[i] = CkStrdup(reads[read_id].c_str());;
        prMSeq->orig_seq[i] = CkStrdup(prMSeq->seq[i]);;
        prMSeq->sqinfo[i].flags = SQINFO_NAME | SQINFO_ID;;
        sprintf(prMSeq->sqinfo[i].name, "read-%d", read_id);
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
    int mini[4]={num_reads,num_reads,num_reads,num_reads},maxi[4]={0,0,0,0};
    int count[5]={0,0,0,0,0};
    int align_read_len = std::strlen(aligned_reads[0]);
    for(int i=0;i<align_read_len;i++) {
        for(int j=0;j<num_reads;j++) {
            if(aligned_reads[j][i]=='A') count[0]++;
            else if(aligned_reads[j][i]=='C') count[1]++;
            else if(aligned_reads[j][i]=='T') count[2]++;
            else if(aligned_reads[j][i]=='G') count[3]++;
            else if(aligned_reads[j][i]=='-') count[4]++;
        }
        std::sort(count,count+5);
        
        for(int j=0;j<4;j++) {
            maxi[j]=std::max(maxi[j],count[j+1]-count[j]);
            mini[j]=std::min(mini[j],count[j+1]-count[j]);
        }
        for(int j=0;j<5;j++) count[j]=0;
    }
    
    //std::stringstream ss;
    //ss << cluster_num;
    dataf.open("cluster.dat", /*std::ofstream::out | std::ofstream::app*/std::ios_base::app);
    for(int i=3;i>=0;i--) {
        dataf << serial_num << " " << (mini[i]*100/num_reads) << " " << (maxi[i]*100/num_reads) << "\n";
        serial_num++;
        serial_num%=4;
    }
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
    printf("Align score = %f\n", dAlnScore);
    process_alignment(prMSeq->nseqs, prMSeq->seq);
    FreeMuscleTree(tree);
    FreeMSeq(&prMSeq);
    if (NULL != piOrderLR) CkFree(piOrderLR, __FUNCTION__, __LINE__);
    free(leftLength);
    free(rightLength);
    free(leafNames);
}

/*******************************************************************************************************************/
/*  Routines from this point onwards are only for test purposes */
/*******************************************************************************************************************/
void testMuscleTree() {
    uint left[] = {0, 2, 10, 9, 5, 7, 13, 16};
    uint right[] = {1, 3, 4, 11, 6, 8, 14, 16};
    float leftLength[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    float rightLength[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    uint leafIds[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    char* leafNames[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8"};
    const char* seqs[] = {"AAA", "AC", "AAT", "AAG", "AAAT", "ACAT", "AAT", "ATAA", "AATA"};
    mseq_t *prMSeq = (mseq_t*)CKCALLOC(1,sizeof(mseq_t));
    prMSeq->seqtype  = SEQTYPE_DNA;
    prMSeq->aligned  = false;
    prMSeq->filename = (char *)"read_file";
    prMSeq->nseqs    = 9;
    prMSeq->seq =  (char **) CKMALLOC((prMSeq->nseqs) * sizeof(char *));
    prMSeq->orig_seq =  (char **) CKMALLOC((prMSeq->nseqs) * sizeof(char *));
    for (int i=0; i<prMSeq->nseqs; ++i) {
        prMSeq->seq[i] = CkStrdup(seqs[i]);;
        prMSeq->orig_seq[i] = CkStrdup(seqs[i]);;
    }
    for(int i=0;i<prMSeq->nseqs;i++) {
        printf("%s\n",prMSeq->seq[i]);
    }
    LogDefaultSetup(&rLog);
    prMSeq->sqinfo = (SQINFO *)CKREALLOC(prMSeq->sqinfo, (prMSeq->nseqs+1) * sizeof(SQINFO));
    {
        tree_t *tree = (tree_t *)malloc(sizeof(tree_t));
        MuscleTreeCreate(tree, 9, 3, left, right, leftLength, rightLength, leafIds, leafNames);
        //TreeValidate(tree);
        MuscleTreeToFile(stdout, tree);
        int *piOrderLR = NULL;
        TraverseTree(&piOrderLR, tree, prMSeq);
        hhalign_para rHhalignPara;
        SetDefaultHhalignPara(&rHhalignPara);
        double dAlnScore = HHalignWrapper(prMSeq, piOrderLR, NULL/*pdSeqWeights*/, 2*prMSeq->nseqs -1/* nodes */, NULL/*prHMMs*/, 0/*iHMMInputFiles*/, -1, rHhalignPara);
        Log(&rLog, LOG_VERBOSE, "Alignment score for first alignment = %f", dAlnScore);        
        for(int i=0;i<prMSeq->nseqs;i++) {
            printf("%s\n",prMSeq->seq[i]);
        }
        if (WriteAlignment(prMSeq, NULL, MSAFILE_A2M, 1000, TRUE)) {
            Log(&rLog, LOG_FATAL, "Could not save alignment");
        } 
        FreeMuscleTree(tree);
        std::cout << "\n\n DONE WITH FIRST TREE \n\n";
    }
    {
        tree_t *tree = (tree_t *)malloc(sizeof(tree_t));
        MuscleTreeCreate(tree, 9, 6, left, right, leftLength, rightLength, leafIds, leafNames);
        std::cout << "\n before tree validate\n";
        //TreeValidate(tree);
        std::cout << "\n after tree validate\n";
        MuscleTreeToFile(stdout, tree);
        printf("\n Sequences before calling TraverseTree\n");
        for(int i=0;i<prMSeq->nseqs;i++) {
            printf("%s\n",prMSeq->seq[i]);
        }
        int *piOrderLR = NULL;
        TraverseTree(&piOrderLR, tree, prMSeq);
        printf("\n Sequences after calling TraverseTree\n");
        for(int i=0;i<prMSeq->nseqs;i++) {
            printf("%s\n",prMSeq->seq[i]);
        }
        hhalign_para rHhalignPara;
        SetDefaultHhalignPara(&rHhalignPara);
        printf("\n Sequences before calling HHalignWrapper\n");
        for(int i=0;i<prMSeq->nseqs;i++) {
            printf("%s\n",prMSeq->seq[i]);
        }

        double dAlnScore = HHalignWrapper(prMSeq, piOrderLR, NULL/*pdSeqWeights*/, 2*prMSeq->nseqs -1/* nodes */, NULL/*prHMMs*/, 0/*iHMMInputFiles*/, -1, rHhalignPara);
        Log(&rLog, LOG_VERBOSE, "Alignment score for first alignment = %f", dAlnScore);        
        if (WriteAlignment(prMSeq, NULL, MSAFILE_A2M, 1000, TRUE)) {
            Log(&rLog, LOG_FATAL, "Could not save alignment");
        } 
        printf("\n Sequences after calling HHalignWrapper\n");
        for(int i=0;i<prMSeq->nseqs;i++) {
            printf("%s\n",prMSeq->seq[i]);
        }

        FreeMuscleTree(tree);
        std::cout << "\n\n DONE WITH SECOND TREE \n\n";
    }
    tree_to_align();
}


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
void tree_to_align() {
    std::string fname,tname;
    fname="read_file";tname="guide_tree_example";
    //cin >> fname >> tname;
    std::cout << fname<<tname <<std::endl;
    mseq_t *prMSeq = (mseq_t*)CKCALLOC(1,sizeof(mseq_t));
    prMSeq->seqtype  = SEQTYPE_DNA;
    prMSeq->aligned  = false;
    prMSeq->filename = (char *)fname.c_str();
    prMSeq->nseqs    = 0;
    prMSeq->seq =  (char **) CKREALLOC(prMSeq->seq, (prMSeq->nseqs+1) * sizeof(char *));
    ReadSequences(prMSeq, (char *)fname.c_str(), SEQTYPE_DNA, SQFILE_FASTA, false,false, 10000, 300);
    for(int i=0;i<prMSeq->nseqs;i++) {
        printf("%s\n",prMSeq->seq[i]);
    }

    /* alignment options to use */
    opts_t rAlnOpts;
    LogDefaultSetup(&rLog);
    SetDefaultAlnOpts(&rAlnOpts);
    rAlnOpts.pcGuidetreeInfile = (char *)tname.c_str();

    if(myAlign(prMSeq, &rAlnOpts)) {
        Log(&rLog, LOG_FATAL, "A fatal error happended during the alignment process");
    }

    if (WriteAlignment(prMSeq, NULL, MSAFILE_A2M, 1000, TRUE)) {
        Log(&rLog, LOG_FATAL, "Could not save alignment");
    } 
    FreeMSeq(&prMSeq);

}

