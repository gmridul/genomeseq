#include "alignwrapper.hpp"

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

void setup_sequences(mseq_t *prMSeq, int nreads, Reads& left_reads, Reads& right_reads) {
    prMSeq->seqtype  = SEQTYPE_DNA;
    prMSeq->aligned  = false;
    prMSeq->filename = "";
    prMSeq->nseqs    = 2*nreads;
    prMSeq->seq =  (char **) CKMALLOC((prMSeq->nseqs) * sizeof(char *));
    prMSeq->orig_seq =  (char **) CKMALLOC((prMSeq->nseqs) * sizeof(char *));
    for (int i=0; i<nreads; ++i) {
        prMSeq->seq[2*i] = CkStrdup(left_reads[i].c_str());;
        prMSeq->orig_seq[2*i] = CkStrdup(left_reads[i].c_str());;
        prMSeq->seq[2*i+1] = CkStrdup(right_reads[i].c_str());;
        prMSeq->orig_seq[2*i+1] = CkStrdup(right_reads[i].c_str());;
    }
    for(int i=0;i<prMSeq->nseqs;i++) {
        printf("%s\n",prMSeq->seq[i]);
    }
    LogDefaultSetup(&rLog);
    prMSeq->sqinfo = (SQINFO *)CKREALLOC(prMSeq->sqinfo, (prMSeq->nseqs+1) * sizeof(SQINFO));
}



void align_cluster(int n, int ni, uint *left, uint *right, float *leftLength, float *rightLength, uint *leafIds, char **leafNames, uint root, mseq_t *prMSeq) {
    tree_t *tree = (tree_t *)malloc(sizeof(tree_t));
    std::cout << "n=" << n << ", root=" << root << "\n";
    MuscleTreeCreate(tree, n, root, left, right, leftLength, rightLength, leafIds, leafNames);
    MuscleTreeToFile(stdout, tree);
    int *piOrderLR = NULL;
    TraverseTree(&piOrderLR, tree, prMSeq);
    /*
    for (int i=0; i<(2*n-1); i++) {
        register int riAux = DIFF_NODE * i;
        std::cout << piOrderLR[riAux + LEFT_NODE] << " " << piOrderLR[riAux + RGHT_NODE] << " " << piOrderLR[riAux + PRNT_NODE] << "\n";
    }
    */
    hhalign_para rHhalignPara;
    SetDefaultHhalignPara(&rHhalignPara);
    double dAlnScore = HHalignWrapper(prMSeq, piOrderLR, NULL/*pdSeqWeights*/, 2*n-1/* nodes */, NULL/*prHMMs*/, 0/*iHMMInputFiles*/, -1, rHhalignPara);
    Log(&rLog, LOG_VERBOSE, "Alignment score for first alignment = %f", dAlnScore);        
    if (WriteAlignment(prMSeq, NULL, MSAFILE_A2M, 1000, TRUE)) {
        Log(&rLog, LOG_FATAL, "Could not save alignment");
    } 
    FreeMuscleTree(tree);
}


void testSetupSequences(mseq_t *prMSeq) {
    char* seqs[] = {"AAA", "AC", "AAT", "AAG", "AAAT", "ACAT", "AAT", "ATAA", "AATA"};
    prMSeq->seqtype  = SEQTYPE_DNA;
    prMSeq->aligned  = false;
    prMSeq->filename = "read_file";
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
}

void testMuscleTree() {
    uint left[] = {0, 2, 10, 9, 5, 7, 13, 16};
    uint right[] = {1, 3, 4, 11, 6, 8, 14, 16};
    float leftLength[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    float rightLength[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    uint leafIds[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    char* leafNames[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8"};
    char* seqs[] = {"AAA", "AC", "AAT", "AAG", "AAAT", "ACAT", "AAT", "ATAA", "AATA"};
    mseq_t *prMSeq = (mseq_t*)CKCALLOC(1,sizeof(mseq_t));
    prMSeq->seqtype  = SEQTYPE_DNA;
    prMSeq->aligned  = false;
    prMSeq->filename = "read_file";
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

