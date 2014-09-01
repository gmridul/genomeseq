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

/*************************************************************************************************************************/
void tree_to_align(char* t) {
    int no_of_seq=4;
    std::string fname,tname;
    fname="read_file";tname="guide_tree_example";
    //cin >> fname >> tname;
    std::cout << fname<<tname <<std::endl;
    std::string tree;
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

