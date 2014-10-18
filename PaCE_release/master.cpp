  /***************************************************************************
  *
  * If you are a U.S. Governmental agency or a non-profit educational
  * institution, you may use this Product, royalty free, if you agree to
  * the below terms.  If you agree to these terms, contact Dr. Srinivas Aluru
  * at aluru@iastate.edu with the following information:
  * (i)   that you agree with the terms of this license;
  * (ii)  identification of the U.S. Governmental agency or non-profit
  * educational institution agreeing to this license; and
  * (iii) the contact name, address, phone number, and email at your place of
  * business.
  * All other interested parties should contact licensing@iastate.edu.
  * ________________________________________________________________________
  *                 Research License for PaCE Software ("Product")
  * Copyright ¨ 2003 Iowa State University Research Foundation, Inc.
  * All rights reserved
  * This material is based on work supported by the National Science
  * Foundation under Grant Nos. EIA-0130861 and ACI-0203782.
  * READ THIS LICENSE AGREEMENT CAREFULLY BEFORE USING THIS PRODUCT. BY
  * USING THIS PRODUCT YOU INDICATE YOUR ACCEPTANCE OF THE TERMS OF THIS
  * LICENSE. THESE TERMS APPLY TO YOU AND ANY SUBSEQUENT USER-LICENSEE
  * OF THIS PRODUCT.
  *
  * Iowa State University Research Foundation, Inc. ("Licensor") retains the
  * ownership of this copy and any subsequent copies of the Product. Licensor
  * grants to Licensee, a U.S. Governmental agency or a non-profit educational
  * institution ("Licensee"), a non-exclusive, royalty free, non-transferable
  * license to use the copy of the Product in accordance with the terms and
  * conditions of this License Agreement.
  * 1. Permitted Uses.  Licensee may:
  * a) use the Product solely for Licensee's own internal research purposes.
  * b) alter, modify, or adapt the Product for Licensee's own internal
  * research purposes.
  * 2. Prohibited Uses.  Licensee may not:
  * a) transfer, distribute, lease or sub-license the Product.
  * b) use the Product for any commercial purpose.
  * c) charge for access or viewing of the Product.
  * d) alter, modify, or adapt the Product or documentation, or portions
  * thereof except as provided herein.
  *
  * Limited Warranty and Liability:
  * LICENSOR MAKES NO WARRANTY OR REPRESENTATION, EITHER EXPRESS OR IMPLIED,
  * WITH RESPECT TO THE PRODUCT, INCLUDING ITS QUALITY, PERFORMANCE,
  * MERCHANTABILITY, OR FITNESS FOR A PARTICULAR PURPOSE OR THAT THE USE OF
  * THE PRODUCT OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS,
  * TRADEMARKS, OR OTHER RIGHTS. THE PRODUCT PROVIDED HEREUNDER IS ON AN
  * "AS IS" BASIS.  IN NO EVENT WILL LICENSOR BE LIABLE FOR DIRECT, INDIRECT,
  * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR
  * INABILITY TO USE THE PRODUCT OR DOCUMENTATION, EVEN IF ADVISED OF THE
  * POSSIBILITY OF SUCH DAMAGES.
  *
  * General:
  * Licensor retains all rights not expressly granted herein. Nothing in this
  * License Agreement constitutes a waiver of Licensor's rights under United
  * States copyright law. This License and Licensee's right to use the Product
  * automatically terminate without notice from Licensor if Licensee fails to
  * comply with any provision of this License Agreement, or any terms and
  * conditions associated with the transfer of this Product. Upon termination,
  * you will remove and destroy all copies of the Product. This License
  * Agreement is governed by the laws of the State of Iowa.
  *
  *****************************************************************************/





/*   master.c  :  This file supports the code for the master and slave model */
/*                 The contents of this master slave is specific to this   */
/*                 application though the framework of the master and slave */
/*                 could be used generically                                */
 
 #define _XOPEN_SOURCE 500


/* #define MAXBATCHES 100 : obsolete Planned */
#define PACKBUF 15000000  
#define MBUFSIZE 2500000
/* #define BUCKBUF 120000000 : obsolete Planned */

#define M0 0
#define M1 1
#define M2 2
#define M3 3
#define M4 4
#define M5 5
#define M6 6
#define M7 7
#define M8 8
#define M9 9
#define M10 10
#define M11 11
#define M12 12
#define M13 13
#define M14 14
#define MFINAL 15

#define S0 0
#define S1 1
#define S2 2
#define S3 3
#define S4 4
#define S5 5
#define S6 6
#define S7 7
#define S8 8
#define S9 9
#define S10 10
#define S11 11
#define S12 12
#define S13 13
#define SFINAL 14

#define YES 1
#define NO 0


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>
#include <fcntl.h>
#include <math.h>
#include <errno.h>
#include "mastercode.h"
#include "keys.h"
#include "cfg.h"
#include "err.h"
#include "stree.h"
#include "debug.h"


/* global variables */ 

int N=0;
int rank,p;
MPI_Comm myComm;
float Kfactor=0;
int iClusters=0;
int iSingletons=0;
int LOADPERPROC=0;          /* Number of units per proc at a give time*/
int *containArr;
int iContained=0;
unsigned int iPairsAccepted=0,iPairsRejected=0;
unsigned long tOnlyFetch=0;
int numbuckets=0;
int batchsize=10;
char pbuf[PACKBUF],rpbuf[PACKBUF];
char estFile[FNAME];
char CFGFile[FNAME];

int g_startEST,g_endEST,g_inumLenLessK,g_isizLenLessK;
int g_iLongest=0;
unsigned long g_startOffset=0,g_endOffset=0;
unsigned int myFragsLen=0;
unsigned int g_iNoLocalFrags=0,g_iExpecKeys=0;
unsigned int g_iFragKeys=0;
struct est *AllESTs=NULL,*g_localESTs=NULL; 
struct estgi *allEstGis=NULL;
struct FileOffset *AllESTlocator;
int *AllESTlen;
struct fkp *fragkey;
struct fkp *final_fragkeys;
unsigned int g_iSortedKeys=0;
unsigned int cmask32;
const short intSize = sizeof(int)*8;
int lastNode=-1;
int mindep=30;
MPI_Datatype pairtype,respairtype;
int g_iMaxPairsPerNode=0;
int g_iMaxPairsNode=-1;
int *deviArray,*m5times;
unsigned int g_treeSpace=0;
int keysFree=0,fkeysFree=0;
unsigned int g_maxBuckStrSize=0;

long g_partSize=0;

int *g_dataPartition=NULL;
int *g_itStarters=NULL;
char *g_recv_str=NULL;

unsigned int *g_TreeStringPtr=NULL;
char *g_TreeStrings=NULL;
char *g_TreeStringMarker=NULL;
unsigned int g_TreeStringsSize=0;

//unsigned long t_fseek=0,t_sfseek=0;
//unsigned long t_fscanf=0,t_sfscanf=0;
unsigned long t_pread=0,t_spread=0;
unsigned int n_callsToSeqtoMem=0;
unsigned int n_bytesReadFromIO=0;

int g_IntNodesProcessed=0;

unsigned int g_TotalSeqLen=0,g_AvgSeqLen=0;

int g_matchscore=0;
int g_bTranscriptsTogether= 0;  /* flag for clustering transcript isoforms together - default is 0 */
int g_bReportSplicedCandidates = 0;  /* flag for reporting spliced candidates - default is 0 */
int g_bReportMaximalPairs= 0;  /* flag for reporting maximal pairs - default is 0 */
int g_bReportAcceptedPairs= 0;  /* flag for reporting accepted pairs - default is 0 */
int g_bReportMaximalSubstrings= 0;  /* flag for reporting maximal substrings (slave side)- default is 0 */
int g_bReportGeneratedPairs= 0;  /* flag for reporting all pairs generated (slave side)- default is 0 */
int g_bDumpClustersMidway = 0;  /* flag for reporting incremental cluster file update */

int g_ReportPairsCountUnit=1; /* basic unit for counting the number of pairs */

int g_iKeep_Mbuf_Full = 0; /* How much to keep the mbuf full ? */
int g_bMPI_Block_Sends = 0; /* Block Sends or not */ 

char g_sOutputFolder[MAXFILENAMESIZE];

int EndToEndScoreRatioThreshold,EndtoEndAlignLenThreshold,MaxScoreRatioThreshold,TranscriptCoverageThreshold;
int MaxStringsInABucket=100000;
unsigned int g_BucketSizeThreshold=100000000;
FILE *fpSpliced=NULL;
FILE *fpMaximalPairs=NULL;
FILE *fpAcceptedPairs=NULL;
FILE *fpMaximalSubstrings=NULL;
FILE *fpGeneratedPairs=NULL;
FILE *fpLargeMerges=NULL;

//FILE **g_fp;
int *g_fd;
struct FileSize *global_fs=NULL;
int nFilesToRead=0;

int g_bOutputLargeMerges = 0;
int LargeClusterThreshold=500;

extern int errnofpMaximalSubstrings;

void debugSync(MPI_Comm mycomm);
void printIntArr(char *name,int *arr,int size);
void checkInRange(int a,int lb,int ub,char *tag);

#include "stree2.h"


void reportClustersMidway(struct mastVar *);
struct FileSize *getInputFileStats(char *,int *);


int main(int argc,char *argv[])
{
   int k=0;
   int gRank,gP;
  
   double start,end,slot1,slot2;
   double tTotalTime=0,tLoadTime=0,tSlideTime=0,tBuckTime=0,tBuildTime=0,tSortNodesTime=0,tReadAgain=0,tSlaveTime=0;

   int splitkey,slaveReady=0,q=0,wrong=0;


   struct ufind *ESTcluster;

   void bucketing(struct fkp *fkey, int size,int k);
   void master(struct ufind *ESTcluster);
   void slave(int mastRank);
   int retrieveK();
   void doCloneMates(struct ufind *,int );



   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&gRank);
   MPI_Comm_size(MPI_COMM_WORLD,&gP); 

   if(gRank==0) {
	printf(" ******************************************************************************\n");
	printf(" If you are a U.S. Governmental agency or a non-profit educational   \n");
	printf(" institution, you may use this Product, royalty free, if you agree to   \n");
	printf(" the below terms.  If you agree to these terms, contact Dr. Srinivas Aluru   \n");
	printf(" at aluru@iastate.edu with the following information:  \n");
	printf(" (i)   that you agree with the terms of this license;  \n");
	printf(" (ii)  identification of the U.S. Governmental agency or non-profit   \n");
	printf(" educational institution agreeing to this license; and  \n");
	printf(" (iii) the contact name, address, phone number, and email at your place of   \n");
	printf(" business.  \n");
	printf(" All other interested parties should contact licensing@iastate.edu.      \n");
	printf(" ________________________________________________________________________  \n");
	printf("                 Research License for PaCE Software (\"Product\")  \n");
	printf(" Copyright ¨ 2003 Iowa State University Research Foundation, Inc.   \n");
	printf(" All rights reserved   \n");
	printf(" This material is based on work supported by the National Science   \n");
	printf(" Foundation under Grant Nos. EIA-0130861 and ACI-0203782.  \n");
	printf(" READ THIS LICENSE AGREEMENT CAREFULLY BEFORE USING THIS PRODUCT. BY   \n");
	printf(" USING THIS PRODUCT YOU INDICATE YOUR ACCEPTANCE OF THE TERMS OF THIS   \n");
	printf(" LICENSE. THESE TERMS APPLY TO YOU AND ANY SUBSEQUENT USER-LICENSEE   \n");
	printf(" OF THIS PRODUCT.   \n");
	printf(" Iowa State University Research Foundation, Inc. (\"Licensor\") retains the   \n");
	printf(" ownership of this copy and any subsequent copies of the Product. Licensor   \n");
	printf(" grants to Licensee, a U.S. Governmental agency or a non-profit educational   \n");
	printf(" institution (\"Licensee\"), a non-exclusive, royalty free, non-transferable   \n");
	printf(" license to use the copy of the Product in accordance with the terms and   \n");
	printf(" conditions of this License Agreement.   \n");
	printf(" 1. Permitted Uses.  Licensee may:  \n");
	printf(" a) use the Product solely for Licensee's own internal research purposes.  \n");
	printf(" b) alter, modify, or adapt the Product for Licensee's own internal   \n");
	printf(" research purposes.  \n");
	printf("    \n");
	printf(" 2. Prohibited Uses.  Licensee may not:  \n");
	printf(" a) transfer, distribute, lease or sub-license the Product.  \n");
	printf(" b) use the Product for any commercial purpose.  \n");
	printf(" c) charge for access or viewing of the Product.  \n");
	printf(" d) alter, modify, or adapt the Product or documentation, or portions   \n");
	printf(" thereof except as provided herein.  \n");
	printf(" Limited Warranty and Liability:  \n");
	printf(" LICENSOR MAKES NO WARRANTY OR REPRESENTATION, EITHER EXPRESS OR IMPLIED,   \n");
	printf(" WITH RESPECT TO THE PRODUCT, INCLUDING ITS QUALITY, PERFORMANCE,  \n");
	printf(" MERCHANTABILITY, OR FITNESS FOR A PARTICULAR PURPOSE OR THAT THE USE OF   \n");
	printf(" THE PRODUCT OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS,   \n");
	printf(" TRADEMARKS, OR OTHER RIGHTS. THE PRODUCT PROVIDED HEREUNDER IS ON AN   \n");
	printf(" \"AS IS\" BASIS.  IN NO EVENT WILL LICENSOR BE LIABLE FOR DIRECT, INDIRECT,   \n");
	printf(" SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR   \n");
	printf(" INABILITY TO USE THE PRODUCT OR DOCUMENTATION, EVEN IF ADVISED OF THE   \n");
	printf(" POSSIBILITY OF SUCH DAMAGES.   \n");
	printf(" General:  \n");
	printf(" Licensor retains all rights not expressly granted herein. Nothing in this   \n");
	printf(" License Agreement constitutes a waiver of Licensor's rights under United   \n");
	printf(" States copyright law. This License and Licensee's right to use the Product   \n");
	printf(" automatically terminate without notice from Licensor if Licensee fails to   \n");
	printf(" comply with any provision of this License Agreement, or any terms and   \n");
	printf(" conditions associated with the transfer of this Product. Upon termination,   \n");
	printf(" you will remove and destroy all copies of the Product. This License   \n");
	printf(" Agreement is governed by the laws of the State of Iowa.   \n");
	printf(" ******************************************************************************\n\n\n");
   }

   int pLimit;

   pLimit = gP;
   if(argc==5) {
   	pLimit = atoi(argv[4]);
	if(pLimit<2 || pLimit > gP) {
		if(rank==0) printf("Warning: pLimit reset to %d\n",gP);
		pLimit = gP;
	}
   }

   if(pLimit<2) {
     printf("Error : There has to be at least two processors!\n");
     MPI_Finalize();
     return 0;
   }

   
   q=p=pLimit-1;    /* p is both the number of slaves and the rank of the master */
   
   if(gRank<p) splitkey=0;
   else if(gRank==p) splitkey=1;
   else splitkey=2; // dummy processor class 
     
   MPI_Comm_split(MPI_COMM_WORLD,splitkey,gRank,&myComm);
   MPI_Comm_rank(myComm,&rank);
   MPI_Comm_size(myComm,&p);

   if(gRank>=pLimit) {
   	printf("Rank=%d: Processor retired.\n",gRank);
     MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
   }

   /* Added : Ananth 6/10/2004 */
   if(argc<2) {
        if(gRank==0)  {
           printf("Usage : PaCE {EST filename} {#ESTfrags} [Optional: PaCE.cfg file location] [Optional: Number of processors to run] \n");
	}
	MPI_Finalize();
        return 0;
   }
   N = atoi(argv[2]);
   strcpy(estFile,argv[1]);
   /* Added PaCE.cfg file as an optional 3rd parameter : Ananth 6/10/2004 */
   strcpy(CFGFile,"PaCE.cfg"); /* default it to PaCE.cfg in the cwd*/
   if(argc==4) {
   	strcpy(CFGFile,argv[3]);
   }
   /* End of Change : Ananth 6/10/2004*/
   
   if(splitkey==0)
   {
      start =  MPI_Wtime();
      k = retrieveK(); 
      /* STEP 1 : Load roughly N/p ESTs into memory */

      readAllESTs(estFile,N,k);


      slot2 =  MPI_Wtime();
      tLoadTime = slot2 - start;
#ifdef debug
      printf("Rank=%d : Time for phase: Load<%d>\n",rank,(int) tLoadTime); fflush(stdout);
      printf("Rank=%d : Space After Loading<%d> \n",rank,spaceSoFar); fflush(stdout);
#endif


     
      /* STEP 2 : Compute Key for all ESTs and rev complements */
      /*fprintf(ofp,"Rank=%d: windowSize=%d\n",rank,k);*/
      slot1 = slot2;

      slideWindowK(k);

      slot2 =  MPI_Wtime();
      tSlideTime= slot2 - slot1;
#ifdef debug
      printf("Rank=%d : Time for phase: Slide<%d>\n",rank,(int) tSlideTime); fflush(stdout);
      printf("Rank=%d : Space After Sliding<%d> \n",rank,spaceSoFar); fflush(stdout);
#endif

      /* STEP 3 : Sort the keys in parallel*/ 
      slot1 = slot2;

      bucketing(fragkey,g_iFragKeys,k); 

      slot2 =  MPI_Wtime();
      tBuckTime= slot2 - slot1;
      spaceSoFar-=keysFree;
#ifdef debug
      printf("Rank=%d : Time for phase: Bucket<%d>\n",rank,(int) tBuckTime); fflush(stdout);
      printf("Rank=%d : Space After Bucketing<%d> \n",rank,spaceSoFar); fflush(stdout);
#endif


      /* Step 4 : Build the Suffix tree in parts */
      slot1 = slot2;
      buildTree(k); 
      slot2 =  MPI_Wtime();
      tBuildTime= slot2 - slot1;
      spaceSoFar-=fkeysFree;
#ifdef debug
      printf("Rank=%d : Time for phase: Build<%d>\n",rank,(int) tBuildTime); fflush(stdout);
      printf("Rank=%d : Space After GST Construction<%d> \n",rank,spaceSoFar); fflush(stdout);
#endif


      /* Step 5 : Locally sort all the local nodes */
      slot1 = slot2;
      localSortNodes();  
      slot2 =  MPI_Wtime();
      tSortNodesTime= slot2 - slot1;
#ifdef debug
      printf("Rank=%d : Time for phase: SortNodes<%d>\n",rank,(int) tSortNodesTime); fflush(stdout);
      printf("Rank=%d : Space After SortNodes<%d> \n",rank,spaceSoFar); fflush(stdout);
#endif


      /* Step 6 : Read Local ESTs into memory again */
      slot1 = slot2;
      readTreeStrings(estFile,N); 
      slot2 =  MPI_Wtime();
      tReadAgain= slot2 - slot1;
#ifdef debug
      printf("Rank=%d : Time for phase: ReadAgain<%d>\n",rank,(int) tReadAgain); fflush(stdout);
      printf("Rank=%d : Space After ReadAgain<%d> \n",rank,spaceSoFar); fflush(stdout);
#endif

      /* minlen = retrieveMinLen(); */
      /*generatePairs(minlen); this was for testing purpose  */

      slaveReady=1;
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      //MPI_Comm_size(MPI_COMM_WORLD,&p);
	 p = pLimit;
      
      /* MPI_Send(&slaveReady,1,MPI_INT,q,1000,MPI_COMM_WORLD); */
      
      slot1 = MPI_Wtime();
      slave(q);
      end =  MPI_Wtime();
      tSlaveTime = end - slot1;
      tTotalTime = end - start;
#ifdef debug
      printf("Rank=%d : Time for phase: Slave<%d>, total<%d>\n",rank,(int) tSlaveTime, (int) tTotalTime); fflush(stdout);
      printf("Rank=%d : Space After MasterSlave<%d> \n",rank,spaceSoFar); fflush(stdout);
#endif

      /* printf("Slave=%d : Time Load<%d>+Slide<%d>+Bucketing<%d>+Build<%d>+SortNodes<%d>+Slave<%d>= %d\n",rank,(int) tLoadTime,(int) tSlideTime,(int) tBuckTime,(int) tBuildTime,(int) tSortNodesTime,(int) tSlaveTime,(int) tTotalTime);*/
      printf("Time taken by slave %d : Load<%d> + Preprocessing Phase<%d> + Clustering Phase<%d>= %d secs\n",rank,(int) tLoadTime,(int) tSlideTime+ (int) tBuckTime+ (int) tBuildTime + (int) tSortNodesTime + (int) tReadAgain, (int) tSlaveTime,(int) tTotalTime);
      fflush(stdout);
      /*fclose(ofp);*/
   } 
   else
   {


          printf("Master processor is %d\n",rank); fflush(stdout);
       /* Initialize cluster and process RAFL */
       ESTcluster = MakeSet(N);
          printf("Master: %d Clusters Initialized..\n",N); fflush(stdout);
       doCloneMates(ESTcluster,N);
          printf("Master: Clone Mates Initialized..\n"); fflush(stdout);

       MPI_Barrier(MPI_COMM_WORLD);
       MPI_Comm_rank(MPI_COMM_WORLD,&rank);
       //MPI_Comm_size(MPI_COMM_WORLD,&p);
	  p = pLimit;
       if(rank==q) {
          printf("Master is %d\n",rank);
	  	wrong=0;
         	if(!wrong)  master(ESTcluster);
       } 
   } 
   
#ifdef debug
   printf("Rank=%d: Finalizing\n",rank);
#endif
   MPI_Finalize();
   return 0;
} /* end main */


/*******************************************************************/
/***************START OF MASTER STATE ROUTINES *********************/


   void initMaster(struct mastVar *st)
   {
       int i;
       struct pair dumpair;
       struct resPair dumrespair;
       struct dynParams dumDynParams;
       char s[200];

	  void initBuf(struct cBuffer *,int );
	  void initWaitQ(struct waitq *,int );
       FILE *fileOpenForWrite(char *fname);

       #ifdef debug
          /*printf("Master: initMaster\n");*/
       #endif
       st->nSlaves=rank;
       /*st->UFcluster = MakeSet(N); */ /* init clusters */
       LOADPERPROC = getLoadPerProc();
       /* initBuf(&(st->MBUF),LOADPERPROC*(st->nSlaves)*MAXBATCHES);*/
       /*initBuf(&(st->MBUF),MAXNEWPAIRS*(st->nSlaves));*/
       initBuf(&(st->MBUF),MBUFSIZE);

       needs = sizeof(int)*(st->nSlaves);
       st->active = (int *)malloc(needs);
       checkAlloc(st->active,"active");
       for(i=0;i<(st->nSlaves);i++) st->active[i]=1;
       st->ACTIVE = st->nSlaves;
       st->TORECV = st->nSlaves;

       initWaitQ(&st->waitQ,st->nSlaves + 1);

       needs = sizeof(struct pair)*LOADPERPROC;
       st->nxtBatch = (struct pair*)malloc(needs);
       checkAlloc(st->nxtBatch,"nxtBatch");

       needs = sizeof(struct resPair)*LOADPERPROC;
       st->newResults = (struct resPair*)malloc(needs);
       checkAlloc(st->newResults,"newResults") ;

       needs = sizeof(int)*N;
       containArr = (int *)malloc(needs);
       checkAlloc(containArr,"containArr");
       for(i=0;i<N;i++)
       {
          containArr[i]=-1;
       }

       buildPairMPI(&dumpair,&pairtype);
       buildresPairMPI(&dumrespair,&respairtype);
       st->send_req = MPI_REQUEST_NULL;
       st->P = st->R = st->W = 0;
       st->tSendTime = st->tIdleTime = 0;
       st->tsSendTime = st->tsIdleTime = 0;

       
       needs = sizeof(int)*st->nSlaves;
       m5times = (int *)malloc(needs);
       checkAlloc(m5times,"m5times");
       deviArray = (int *)malloc(needs);
       checkAlloc(deviArray,"deviArray");
       for(i=0;i<st->nSlaves;i++)
       {
          m5times[i]=deviArray[i]=0;
       }

 
       dumDynParams = getDynParams();
       g_matchscore = dumDynParams.match;
       g_bTranscriptsTogether = getTranscriptsTogether();
       g_bReportSplicedCandidates= getReportSplicedCandidates();
       g_bReportMaximalPairs = getReportMaximalPairs();
       g_bReportAcceptedPairs = getReportAcceptedPairs();
       EndtoEndAlignLenThreshold = getThreshold("EndtoEndAlignLenThreshold");
       EndToEndScoreRatioThreshold = getThreshold("EndToEndScoreRatioThreshold");
       MaxScoreRatioThreshold = getThreshold("MaxScoreRatioThreshold");
       TranscriptCoverageThreshold = getThreshold("TranscriptCoverageThreshold");
       g_ReportPairsCountUnit = getThreshold("ReportPairsCountUnit");
	  if(g_ReportPairsCountUnit<=0) g_ReportPairsCountUnit=1;

	  LargeClusterThreshold = getThreshold("LargeClusterThreshold");
	  g_bOutputLargeMerges = getOutputLargeMerges();


       g_iKeep_Mbuf_Full = getKeep_Mbuf_Full();
       g_bMPI_Block_Sends = getMPI_Block_Sends();

       getOutputFolder(g_sOutputFolder);

       if(g_bReportSplicedCandidates) {
	         strcpy(s,"");
		 sprintf(s,"%s/spliced_candidates.%d.%d.PaCE",g_sOutputFolder,N,rank);
		 fpSpliced=fileOpenForWrite(s);
       }

       if(g_bReportMaximalPairs) {
	         strcpy(s,"");
		 sprintf(s,"%s/maximal_pairs.%d.%d.PaCE",g_sOutputFolder,N,rank);
		 fpMaximalPairs=fileOpenForWrite(s);
       }

       if(g_bReportAcceptedPairs) {
	         strcpy(s,"");
	         sprintf(s,"%s/accepted_pairs.%d.%d.PaCE",g_sOutputFolder,N,rank);
	         fpAcceptedPairs=fileOpenForWrite(s);
       }

	  if(g_bOutputLargeMerges) { 
	  	strcpy(s,"");
		sprintf(s,"%s/large_merges.%d.%d.PaCE",g_sOutputFolder,N,rank);
		fpLargeMerges =fileOpenForWrite(s);
	  }


       g_bReportGeneratedPairs= getReportGeneratedPairs();

       g_bDumpClustersMidway = getDumpClustersMidway();
       st->numClusterDumpSteps=0;

       needs = sizeof(int)*(N);
       st->UFcluster_size= (int *)malloc(needs);
       checkAlloc(st->UFcluster_size,"UFcluster_size");
       for(i=0;i<N;i++) st->UFcluster_size[i]=1;	

   } /* end initMaster */

   void handlerM0(struct mastVar *st)
   {
      int src,position;
      MPI_Status status;
      /*double slot2,slot1;*/
	 struct timeval t_1,t_2;
      static int s_recvCounter=0;
      static int s_dotCounter=0;

      if(((++s_recvCounter)%10)==0) {
		printf("."); 
#ifdef BGL_debug
		/*if((st->ACTIVE<=(p-1) && st->ACTIVE>=4000) || (st->ACTIVE<=2))*/
		if((st->ACTIVE<=(p-1) && st->ACTIVE>=(p-50)) || (st->ACTIVE<=2))
			fflush(stdout);
#endif
		if((++s_dotCounter%60)==0) { 
		/*	if((st->ACTIVE<=10)) {*/
				printf(" Mbuf_size=%d\n",getBufSize(st->MBUF));
			/*} else printf("\n"); */
			fflush(stdout); 
			s_dotCounter=0;
		}
		s_recvCounter=0;
      }
      /*slot1 = MPI_Wtime();*/
      gettimeofday(&t_1,0);
      MPI_Recv(rpbuf,PACKBUF,MPI_PACKED,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      /*slot2 = MPI_Wtime();
      st->tIdleTime+=slot2-slot1;*/

      gettimeofday(&t_2,0);
      st->tIdleTime+= (t_2.tv_sec-t_1.tv_sec)*1000000 + (t_2.tv_usec-t_1.tv_usec);
      if(st->tIdleTime>=MAXULONGBY2) {
               st->tsIdleTime+=st->tIdleTime/1000000;
               st->tIdleTime=0;
      }

      src = status.MPI_SOURCE;

      if(status.MPI_TAG==100)
      {
        printf("Error reported from slave=%d : \n",src);
	fflush(stdout);
        terminate();
	return; 
      }
      position=0;
      MPI_Unpack(rpbuf,PACKBUF,&position,&st->P,1,MPI_INT,MPI_COMM_WORLD);
      MPI_Unpack(rpbuf,PACKBUF,&position,&st->R,1,MPI_INT,MPI_COMM_WORLD);
      if(st->P>0)  MPI_Unpack(rpbuf,PACKBUF,&position,st->newPairs,st->P,pairtype,MPI_COMM_WORLD);
      if(st->R>0)  MPI_Unpack(rpbuf,PACKBUF,&position,st->newResults,st->R,respairtype,MPI_COMM_WORLD);

      st->TORECV--;
      st->dst = src; /*reply should go to src*/ 

	 if(g_ReportPairsCountUnit<=1) {
      	st->numTotalPairs+=st->P;
	 }
	 else {
		if( (st->tmp_numTotalPairs + st->P) <= g_ReportPairsCountUnit)  {
			st->tmp_numTotalPairs += st->P;
		}
		else {
			int quo = (st->tmp_numTotalPairs + st->P) / g_ReportPairsCountUnit;
			int rem = (st->tmp_numTotalPairs + st->P) % g_ReportPairsCountUnit;
			st->numTotalPairs += quo;
			st->tmp_numTotalPairs = rem;
		}
	 }

      /*printf("Master: Received  P=%d R=%d from %d\n", st->P,st->R,src); fflush(stdout);
      for(i=0;i<st->R;i++) {
	      printf("Master:  "); dispresPair(st->newResults[i]); fflush(stdout);
      }*/
   } /* handlerM0 */

   void handlerM1(struct mastVar *st)
   {
   }

   void handlerM2(struct mastVar *st)
   {
      interpretResults(st->newResults,st->R,st->UFcluster,st->UFcluster_size);
      //st->numResultsRecv+=st->R;
	 if(g_ReportPairsCountUnit<=1) {
      	st->numResultsRecv+=st->R;
	 }
	 else {
		if( (st->tmp_numResultsRecv + st->R) <= g_ReportPairsCountUnit)  {
			st->tmp_numResultsRecv += st->R;
		}
		else {
			int quo = (st->tmp_numResultsRecv + st->R) / g_ReportPairsCountUnit;
			int rem = (st->tmp_numResultsRecv + st->R) % g_ReportPairsCountUnit;
			st->numResultsRecv += quo;
			st->tmp_numResultsRecv = rem;
		}
	 }
   }  
   
   void handlerM3(struct mastVar *st)
   {
   }

   void handlerM4(struct mastVar *st)
   {
      interpretResults(st->newResults,st->R,st->UFcluster,st->UFcluster_size);
      //st->numResultsRecv+=st->R;
	 if(g_ReportPairsCountUnit<=1) {
      	st->numResultsRecv+=st->R;
	 }
	 else {
		if( (st->tmp_numResultsRecv + st->R) <= g_ReportPairsCountUnit)  {
			st->tmp_numResultsRecv += st->R;
		}
		else {
			int quo = (st->tmp_numResultsRecv + st->R) / g_ReportPairsCountUnit;
			int rem = (st->tmp_numResultsRecv + st->R) % g_ReportPairsCountUnit;
			st->numResultsRecv += quo;
			st->tmp_numResultsRecv = rem;
		}
	 }
   }

   void handlerM5(struct mastVar *st) {
	int canhold=0,delta=1;
	int Pprime;
	int i=0;

	if(g_bReportGeneratedPairs) {
		st->E=LOADPERPROC;
		st->W=0;
		return;
	}
#ifdef Report_Only_Contained
	for(i=0;i<st->P;i++) {
		containArr[st->newPairs[i].f1]=st->newPairs[i].f2;
	}
	st->E=LOADPERPROC;
	st->W=0;
	return;
#endif

	if(st->ACTIVE<=0) {
	    	// this case is not possible
           printf("Master Error: Case ACTIVE=%d and STATE=M5\n",st->ACTIVE);
		 fflush(stdout);
		 terminate();
	}

	Pprime = insertCBuf(st->newPairs,st->P,&st->MBUF,st->UFcluster);
	if(Pprime<0) { 
		printf("Master Error :State M5 : insertCBuf failed \n");
		fflush(stdout);
		terminate();
		return;
	}

      /* printf("Assigned %d out of %d from %d \n",Pprime,st->P,st->dst); fflush(stdout);*/
	canhold = (st->MBUF.totSize - st->MBUF.currSize)/st->ACTIVE;

	if(Pprime>0) {
		delta = st->P/Pprime;
      	/* calculate expected value from source next time */
      	st->E = (delta*LOADPERPROC*st->nSlaves)/st->ACTIVE; /* for load balancing */
	}
	else {
		st->E = canhold;
	}
	/*printf("Master : st->E = min{%d,%d,%d} is expected\n",st->E,canhold,MAXNEWPAIRS/st->ACTIVE);
	fflush(stdout);*/
	if(st->E > canhold) st->E = canhold;
	if(st->E > (MAXNEWPAIRS/st->ACTIVE)) st->E = MAXNEWPAIRS/st->ACTIVE;

	//if(st->ACTIVE==1) {
		// populate as much in mbuf as possible
		// this will override everything above
	//	st->E = canhold/g_iKeep_Mbuf_Full;
	//	printf("Master: DEBUG: demanding: active=1, st->E = %d\n",st->E);
	//	fflush(stdout);
	//}

      st->W = prepareNextBatch(&st->MBUF,st->nxtBatch,st->UFcluster);
      /* printf("Assigned E=%d and prepared W=%d\n",st->E,st->W); fflush(stdout);*/

      
      deviArray[st->dst] = ( (deviArray[st->dst]*m5times[st->dst]) + ((Pprime-LOADPERPROC)*100/LOADPERPROC) )/ (m5times[st->dst]+1);
      m5times[st->dst]++;
      
   }

   void handlerM6(struct mastVar *st)
   {
      int position=0;
      MPI_Status sendstat;
	 struct timeval t_1,t_2; 

      /*slot1 = MPI_Wtime();*/
      gettimeofday(&t_1,0);
	 if(!g_bMPI_Block_Sends) MPI_Wait(&st->send_req,&sendstat);
      position=0;
      MPI_Pack(&st->W,1,MPI_INT,pbuf,PACKBUF,&position,MPI_COMM_WORLD);
      MPI_Pack(&st->E,1,MPI_INT,pbuf,PACKBUF,&position,MPI_COMM_WORLD);

      /* Flush outstanding sends when running on BlueGene/L
         This will ensure message is sent immediately rather than
         made to wait until the Wait() call on BlueGene/L
      */
      if(!g_bMPI_Block_Sends)
      	MPI_Isend(pbuf,position,MPI_PACKED,st->dst,1,MPI_COMM_WORLD,&st->send_req);  
	 else
	#ifdef Use_Ssend
      	MPI_Ssend(pbuf,position,MPI_PACKED,st->dst,1,MPI_COMM_WORLD);  
	#else
      	MPI_Send(pbuf,position,MPI_PACKED,st->dst,1,MPI_COMM_WORLD);  
	#endif

      st->TORECV++;
      /*slot2 = MPI_Wtime();*/
      /*st->tSendTime+=slot2-slot1;*/
	 gettimeofday(&t_2,0);
	 st->tSendTime+= (t_2.tv_sec-t_1.tv_sec)*1000000 + (t_2.tv_usec-t_1.tv_usec);
	 if(st->tSendTime>=MAXULONGBY2) {
		st->tsSendTime+=st->tSendTime/1000000;
		st->tSendTime=0;
	 }
       /* printf("Master: Sent W=%d, E=%d to %d\n",st->W,st->E,st->dst); fflush(stdout);*/
    
   }

   void handlerM7( struct mastVar *st)
   {
      int position=0;
      MPI_Status sendstat;
	 struct timeval t_1,t_2;
     
      /*double slot1,slot2;*/
      
      /*slot1 = MPI_Wtime();*/
	 gettimeofday(&t_1,0);


	 if(!g_bMPI_Block_Sends) MPI_Wait(&st->send_req,&sendstat);

      position=0;
      MPI_Pack(&(st->W),1,MPI_INT,pbuf,PACKBUF,&position,MPI_COMM_WORLD);
      MPI_Pack(&(st->E),1,MPI_INT,pbuf,PACKBUF,&position,MPI_COMM_WORLD);
      MPI_Pack(st->nxtBatch,st->W,pairtype,pbuf,PACKBUF,&position,MPI_COMM_WORLD);

      /* Flush outstanding sends when running on BlueGene/L
         This will ensure message is sent immediately rather than
         made to wait until the Wait() call on BlueGene/L
      */
	 if(!g_bMPI_Block_Sends)
      	MPI_Isend(pbuf,position,MPI_PACKED,st->dst,1,MPI_COMM_WORLD,&(st->send_req) );  
	 else
	#ifdef Use_Ssend
		MPI_Ssend(pbuf,position,MPI_PACKED,st->dst,1,MPI_COMM_WORLD);
	#else
      	MPI_Send(pbuf,position,MPI_PACKED,st->dst,1,MPI_COMM_WORLD);  
	#endif

      /*slot2 = MPI_Wtime();
      st->tSendTime+=slot2-slot1;*/
	 gettimeofday(&t_2,0);
	 st->tSendTime+= (t_2.tv_sec-t_1.tv_sec)*1000000 + (t_2.tv_usec-t_1.tv_usec);
	 if(st->tSendTime>=MAXULONGBY2) {
		st->tsSendTime+=st->tSendTime/1000000;
		st->tSendTime=0;
	 }
      
	 /* printf("Master: Sent W=%d, E=%d to %d\n",st->W,st->E,st->dst); fflush(stdout);*/
      (st->TORECV)++;
      //st->numPairsAssigned=st->numPairsAssigned + st->W;
	 if(g_ReportPairsCountUnit<=1) {
      	st->numPairsAssigned += st->W;
	 }
	 else {
		if( (st->tmp_numPairsAssigned + st->W) <= g_ReportPairsCountUnit)  {
			st->tmp_numPairsAssigned += st->W;
		}
		else {
			int quo = (st->tmp_numPairsAssigned + st->W) / g_ReportPairsCountUnit;
			int rem = (st->tmp_numPairsAssigned + st->W) % g_ReportPairsCountUnit;
			st->numPairsAssigned += quo;
			st->tmp_numPairsAssigned = rem;
		}
	 }
   }

   void handlerM8(struct mastVar *st)
   {
      if(!deQWaitQ(&st->waitQ,&st->dst)) 
      { 
         printf("Error :State M8 Master : deQWaitQ failed \n");
	 fflush(stdout);
         terminate();
	 return;
      }
      st->E=0;
      st->W = prepareNextBatch(&st->MBUF,st->nxtBatch,st->UFcluster);
   }  

   void handlerM9(struct mastVar *st)
   {


       if(st->active[st->dst]!=-1) {
 	    st->active[st->dst]=-1;
            st->ACTIVE--;
            printf("Master : %d made passive with ACTIVE=%d Mbuf_size=%d\n",
				st->dst,st->ACTIVE,getBufSize(st->MBUF)); fflush(stdout);
	    if(g_bDumpClustersMidway && 
		((st->ACTIVE==(3*p/4)) || (st->ACTIVE==(p/2)) || (st->ACTIVE==(p/4)) || (st->ACTIVE<2)) ) {
		reportClustersMidway(st);
	    }
       }
       st->E=0;
       st->W = prepareNextBatch(&st->MBUF,st->nxtBatch,st->UFcluster);
   } 

   void handlerM10(struct mastVar *st)
   {
       enQWaitQ(&st->waitQ,st->dst);
   } 
    
   void handlerM11(struct mastVar *st)
   {
   }

   void handlerM12(struct mastVar *st)
   {
      int src,position;
      MPI_Status status;
      /*double slot1,slot2;*/
      struct timeval t_1,t_2;

      /*slot1 = MPI_Wtime();*/
	 gettimeofday(&t_1,0);
      MPI_Recv(rpbuf,PACKBUF,MPI_PACKED,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      gettimeofday(&t_2,0);
      st->tIdleTime+= (t_2.tv_sec-t_1.tv_sec)*1000000 + (t_2.tv_usec-t_1.tv_usec);
      if(st->tIdleTime>=MAXULONGBY2) {
               st->tsIdleTime+=st->tIdleTime/1000000;
               st->tIdleTime=0;
      }

      /*slot2 = MPI_Wtime();
      st->tIdleTime+=slot2-slot1;*/

      src = status.MPI_SOURCE;

      if(status.MPI_TAG==100) {
		printf("Error reported from slave=%d : \n",src); 
		fflush(stdout);
		terminate();
		return; 
      }
      position=0;
      MPI_Unpack(rpbuf,PACKBUF,&position,&st->P,1,MPI_INT,MPI_COMM_WORLD);
      MPI_Unpack(rpbuf,PACKBUF,&position,&st->R,1,MPI_INT,MPI_COMM_WORLD);
      MPI_Unpack(rpbuf,PACKBUF,&position,st->newPairs,st->P,pairtype,MPI_COMM_WORLD);
      MPI_Unpack(rpbuf,PACKBUF,&position,st->newResults,st->R,respairtype,MPI_COMM_WORLD);
      st->TORECV--;
      interpretResults(st->newResults,st->R,st->UFcluster,st->UFcluster_size);

	 if(g_ReportPairsCountUnit<=1) {
      	st->numResultsRecv+=st->R;
	 }
	 else {
		if( (st->tmp_numResultsRecv + st->R) <= g_ReportPairsCountUnit)  {
			st->tmp_numResultsRecv += st->R;
		}
		else {
			int quo = (st->tmp_numResultsRecv + st->R) / g_ReportPairsCountUnit;
			int rem = (st->tmp_numResultsRecv + st->R) % g_ReportPairsCountUnit;
			st->numResultsRecv += quo;
			st->tmp_numResultsRecv = rem;
		}
	 }

	 if(g_ReportPairsCountUnit<=1) {
      	st->numTotalPairs+=st->P;
	 }
	 else {
		if( (st->tmp_numTotalPairs + st->P) <= g_ReportPairsCountUnit)  {
			st->tmp_numTotalPairs += st->P;
		}
		else {
			int quo = (st->tmp_numTotalPairs + st->P) / g_ReportPairsCountUnit;
			int rem = (st->tmp_numTotalPairs + st->P) % g_ReportPairsCountUnit;
			st->numTotalPairs += quo;
			st->tmp_numTotalPairs = rem;
		}
	 }
   }
 
   void handlerM13(struct mastVar *st)
   {
      int i,tempSize=0,position=0;
      MPI_Status sendstat;
      /*double slot1,slot2;*/
	 struct timeval t_1,t_2;

      for(i=0;i<st->nSlaves;i++) /* send Tag 0 for termination of all slaves */
      { 
	 	gettimeofday(&t_1,0);
          /*slot1 = MPI_Wtime();*/

		if(!g_bMPI_Block_Sends)
          	MPI_Wait(&st->send_req,&sendstat);

          position=0;
          MPI_Pack(&tempSize,1,MPI_INT,pbuf,PACKBUF,&position,MPI_COMM_WORLD);

      /* Flush outstanding sends when running on BlueGene/L
         This will ensure message is sent immediately rather than
         made to wait until the Wait() call on BlueGene/L
      */
		if(!g_bMPI_Block_Sends)
          	MPI_Isend(pbuf,position,MPI_PACKED,i,0,MPI_COMM_WORLD,&st->send_req);  
		else
          #ifdef Use_Ssend
          	MPI_Ssend(pbuf,position,MPI_PACKED,i,0,MPI_COMM_WORLD);  
     	#else
          	MPI_Send(pbuf,position,MPI_PACKED,i,0,MPI_COMM_WORLD);  
		#endif

          /*slot2 = MPI_Wtime();
          st->tSendTime+=slot2-slot1;*/
	 
	 	gettimeofday(&t_2,0);
	     st->tSendTime+= (t_2.tv_sec-t_1.tv_sec)*1000000 + (t_2.tv_usec-t_1.tv_usec);
	 	if(st->tSendTime>=MAXULONGBY2) {
			st->tsSendTime+=st->tSendTime/1000000;
			st->tSendTime=0;
	 	}
          #ifdef debug
             printf("Master: Termination sent to %d\n",i); 
          #endif
      }
   } 

   void handlerM14(struct mastVar *st)
   {
      int src,position;
      MPI_Status status;
      /*double slot1,slot2;*/
 	 struct timeval t_1,t_2;
      int i=0;

      for(i=0;i<st->nSlaves;i++) {

      	/*slot1 = MPI_Wtime();*/
		gettimeofday(&t_1,0);
      	MPI_Recv(rpbuf,PACKBUF,MPI_PACKED,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      	/*slot2 = MPI_Wtime();
      	st->tIdleTime+=slot2-slot1;*/

      	gettimeofday(&t_2,0);
      	st->tIdleTime+= (t_2.tv_sec-t_1.tv_sec)*1000000 + (t_2.tv_usec-t_1.tv_usec);
      	if(st->tIdleTime>=MAXULONGBY2) {
               st->tsIdleTime+=st->tIdleTime/1000000;
               st->tIdleTime=0;
      	}

      	src = status.MPI_SOURCE;

     	if(status.MPI_TAG==100) {
        	printf("Error reported from slave=%d : \n",src);
		fflush(stdout);
        	terminate();
		return; 
      	}
      	position=0;
      	MPI_Unpack(rpbuf,PACKBUF,&position,&st->P,1,MPI_INT,MPI_COMM_WORLD);
      	MPI_Unpack(rpbuf,PACKBUF,&position,&st->R,1,MPI_INT,MPI_COMM_WORLD);
	if(st->P>0) MPI_Unpack(rpbuf,PACKBUF,&position,st->newPairs,st->P,pairtype,MPI_COMM_WORLD);
      	if(st->R>0) MPI_Unpack(rpbuf,PACKBUF,&position,st->newResults,st->R,respairtype,MPI_COMM_WORLD);
      	st->TORECV--;
        	interpretResults(st->newResults,st->R,st->UFcluster,st->UFcluster_size);
	   	if(g_ReportPairsCountUnit<=1) {
      		st->numResultsRecv+=st->R;
	   	}
	   	else {
			if( (st->tmp_numResultsRecv + st->R) <= g_ReportPairsCountUnit)  {
				st->tmp_numResultsRecv += st->R;
			}
			else {
				int quo = (st->tmp_numResultsRecv + st->R) / g_ReportPairsCountUnit;
				int rem = (st->tmp_numResultsRecv + st->R) % g_ReportPairsCountUnit;
				st->numResultsRecv += quo;
				st->tmp_numResultsRecv = rem;
			}
	   	}

      	if(st->R>0) {
		printf("Master: Info: handlerM14: Last minute results st->R=%d from slave rank %d\n",
		st->R,src);
		fflush(stdout);
      	}
      } /* end for*/
   }

   void reportClusters(struct mastVar *st)
   {
       int i;
       void displayNonContainArr(int *containArr,int N,struct estgi *);
	  struct FileSize *getInputFileStats(char *,int *);

       /*sprintf(dumpfile,"%s.%d.%d.temp",estFile,N,rank);
       fp=fopen(dumpfile,"w");*/

       /*
	* Commented out temporarily for performance studies
	* Ananth 6/22/2004
	*/ 

	// for the master processor
	if(!global_fs || nFilesToRead<=0) global_fs = getInputFileStats(estFile,&nFilesToRead);

#ifndef PERF_STUDY
	if(allEstGis==NULL) {
       		needs = sizeof(struct estgi)*N;
       		allEstGis= (struct estgi*)malloc(needs);
       		checkAlloc(allEstGis,"allEstGis");

       		loadESTGIs(allEstGis,N,global_fs,nFilesToRead);
	}


       iSingletons=0;
       printf("Final Individual Clusters:\n");
       iClusters = dispAllClusters(st->UFcluster,N,allEstGis,&iSingletons);  
       iContained=0;
       for(i=0;i<N;i++)
       {
          if(containArr[i]>=0) iContained++;
       }
       printf("\n");

       displayNonContainArr(containArr,N,allEstGis);
#endif
   
#ifdef debug
       printf("Master: Memory utilization := %d bytes\n",spaceSoFar);
#endif
       printf("Master: Number of ESTs := %d \n",N);

	 if(g_ReportPairsCountUnit<=1) {
       	printf("Master: #Pairs Input := %u, #Pairs Selected :=%u, #Pairs Accepted=%u, #Results =%u\n",
	  		st->numTotalPairs+LOADPERPROC*2*rank,
			st->numPairsAssigned+LOADPERPROC*2*rank,
			iPairsAccepted,st->numResultsRecv);
	 }
	 else {
       	printf("Master: #Pairs Input := %u units, #Pairs Selected :=%u units, #Pairs Accepted=%u, #Results =%u units (UNIT: %d pairs)\n",
	  		st->numTotalPairs + (st->tmp_numTotalPairs+LOADPERPROC*2*rank)/g_ReportPairsCountUnit,
	  		st->numPairsAssigned + (st->tmp_numPairsAssigned+LOADPERPROC*2*rank)/g_ReportPairsCountUnit,
			iPairsAccepted, //iPairsRejected, st->numResultsRecv);
	  		st->numResultsRecv + st->tmp_numResultsRecv/g_ReportPairsCountUnit,
			g_ReportPairsCountUnit);
	 }
	
       printf("Master: #Clusters Output:= %d #Singletons=%d \n",iClusters,iSingletons);
       printf("Master: #Contained ESTs:= %d \n",iContained);
 

    } /* end reportClusters*/

   void reportClustersMidway(struct mastVar *st)
   {
       int i;
	// for the master processor
	if(!global_fs || nFilesToRead<=0) global_fs = getInputFileStats(estFile,&nFilesToRead);

#ifndef PERF_STUDY
	if(allEstGis==NULL) {
       		needs = sizeof(struct estgi)*N;
       		allEstGis= (struct estgi*)malloc(needs);
       		checkAlloc(allEstGis,"allEstGis");

       		loadESTGIs(allEstGis,N,global_fs,nFilesToRead);
	}

       ++(st->numClusterDumpSteps);
       printf("Individual Clusters Midway Step = %d :\n",st->numClusterDumpSteps);
       iSingletons=0;
       iClusters = dispAllClustersMidway(st->UFcluster,N,allEstGis,
					&iSingletons,st->numClusterDumpSteps);  

	/* No need to display contain arrays */
       /*iContained=0;
       for(i=0;i<N;i++)
       {
          if(containArr[i]>=0) iContained++;
       }
       printf("\n");
       displayNonContainArr(containArr,N,allEstGis);*/
#endif
   
#ifdef debug
       /* printf("Master: Memory utilization := %d bytes\n",spaceSoFar);*/
#endif
       printf("Master: (Midway %d Active %d ): Number of ESTs := %d \n",st->numClusterDumpSteps,st->ACTIVE,N);

	  if(g_ReportPairsCountUnit<=1) {
       	printf("Master: (Midway %d Active %d ): #Pairs Input := %u , #Pairs Selected :=%u, #Pairs Accepted=%u, #Results =%u\n",
			st->numClusterDumpSteps,st->ACTIVE,
	  		st->numTotalPairs+LOADPERPROC*2*rank,
			st->numPairsAssigned+LOADPERPROC*2*rank,
			iPairsAccepted,st->numResultsRecv);
	  }
	  else {
       	printf("Master: (Midway %d Active %d ): #Pairs Input := %u units, #Pairs Selected :=%u units, #Pairs Accepted=%u, #Results =%u units (UNIT: %d pairs)\n",
			st->numClusterDumpSteps,st->ACTIVE,
	  		st->numTotalPairs + (st->tmp_numTotalPairs+LOADPERPROC*2*rank)/g_ReportPairsCountUnit,
	  		st->numPairsAssigned + (st->tmp_numPairsAssigned+LOADPERPROC*2*rank)/g_ReportPairsCountUnit,
			iPairsAccepted, //iPairsRejected, st->numResultsRecv);
	  		st->numResultsRecv + st->tmp_numResultsRecv/g_ReportPairsCountUnit,
			g_ReportPairsCountUnit);
	  }

       printf("Master: (Midway %d Active %d ): #Clusters Output:= %d #Singletons=%d \n",st->numClusterDumpSteps,st->ACTIVE,iClusters,iSingletons);

    } /* end reportClustersMidway */


    void initMastVar(struct mastVar *st)
    {
       st->P = st->R = st->E = st->W = 0;
       st->ACTIVE = st->TORECV = 0;
       st->nSlaves = 0;
       st->numTotalPairs = st->numPairsAssigned = st->numResultsRecv = 0;
       st->tmp_numTotalPairs = st->tmp_numPairsAssigned = st->tmp_numResultsRecv = 0;
       st->waitQ.q = 0;
    }
   /*******************************************************************/
   /***********   START OF MASTER CODE  *******************************/

void master(struct ufind *ESTcluster)
{
   struct mastVar starvars;
   int state=M0; /* current state of master */
   double start,end,tDiff=0;
   int i=0,avgdevi=0;


#ifdef debug
   printf("Master: Start of master\n"); fflush(stdout);
#endif
   start = MPI_Wtime();

   initMastVar(&starvars);
   initMaster(&starvars);
   starvars.UFcluster = ESTcluster;

   while(state!=MFINAL) 
   {
     /*printf("Master: STATE=M%d Mbuf_size = %d TORECV = %d\n",
	*	state,getBufSize(starvars.MBUF),starvars.TORECV); 
     * fflush(stdout); 
	*/
     switch (state)
     {
       case M0 : handlerM0(&starvars);
                 if(starvars.P>0 && starvars.R>0) state=M4;
                 if(starvars.P>0 && starvars.R==0) state=M3;
                 if(starvars.P==0 && starvars.R>0) state=M2;
                 if(starvars.P==0 && starvars.R==0) state=M1;
	 	 break;

       case M1 : handlerM1(&starvars);
		 state=M9;
		 break;

       case M2 : handlerM2(&starvars);
                 state=M9;
		 break;		 

       case M3 : handlerM3(&starvars);
	         state=M5;
		 break;
      
       case M4 : handlerM4(&starvars);
		 state=M5;
		 break; 

       case M5 : handlerM5(&starvars);
	         if(starvars.W >0)  state=M7;
		 else state=M6;
	         break;

       case M6 : handlerM6(&starvars);
		 state=M0;
		 break;

       case M7 : handlerM7(&starvars);
		 if(!emptyBuf(starvars.MBUF) && !emptyWaitQ(starvars.waitQ)) state=M8;
		 else state=M0;
		 break;
		
       case M8 : handlerM8(&starvars);
		 state=M7;
		 break;

       case M9 : handlerM9(&starvars);
	         if(starvars.ACTIVE==0 && starvars.W==0) state=M11;
	         else 
		 if(starvars.W==0) state=M10;
		 else  state=M7;
		 break;
    
       case M10: handlerM10(&starvars);
                 state=M0;
		 break;
		  
       case M11: handlerM11(&starvars);
		 if(starvars.TORECV>0) state=M12;
                 else state=M13;
                 break;
       
       case M12: handlerM12(&starvars);
	 	 state=M11;
	  	 break;

       case M13: handlerM13(&starvars);
		 state=M14;
		 break;

       case M14: handlerM14(&starvars);
		 state=MFINAL;
		 break;

       default : printf("Error :Undefined State %d reached by Master\n",state);
	         fflush(stdout);
                 terminate();
	         return; 
		 
     } /* end switch */
   } /* end while */

   end = MPI_Wtime();
   #ifdef debug
      printf("Master: End of master\n"); fflush(stdout);
   #endif
   tDiff = end - start;
   printf("Master: Send Time := %u sec(s), Idle Time := %u sec(s)\n ",
   		starvars.tsSendTime + starvars.tSendTime/1000000, 
		starvars.tsIdleTime + starvars.tIdleTime/1000000);
   printf("Master: Total Time Taken for Master-Slave := %d sec(s)\n ",
   		(int) tDiff);

   for(i=0;i<starvars.nSlaves;i++) 
   {
     avgdevi+=deviArray[i];
     #ifdef debug
         /*printf("devi[%d]=%d\n",i,deviArray[i]);*/
     #endif
   }
   #ifdef debug
      printf("Avg Deviation = %d\n",avgdevi/starvars.nSlaves);
   #endif

    reportClusters(&starvars);
    

   if(fpSpliced) fclose(fpSpliced);
   if(fpMaximalPairs) fclose(fpMaximalPairs);
   if(fpMaximalSubstrings) fclose(fpMaximalSubstrings);

   free(starvars.MBUF.pairs);
   free(starvars.nxtBatch);
   free(starvars.newResults);
   free(starvars.UFcluster);
   free(starvars.UFcluster_size);
   free(starvars.active);
   free(starvars.waitQ.q);

   /***********   END OF MASTER CODE  *******************************/
   /*******************************************************************/
} /* end master */



/***************************************************************/
/* SLAVE ROUTINES  */

   void reportErr(char *msg,struct slaveVar *sl)
   {
       int position,tempSize=0;
       printf("Rank=%d: Error :%s\n",rank,msg); fflush(stdout);
       position=0;
       MPI_Pack(&tempSize,1,MPI_INT,pbuf,PACKBUF,&position,MPI_COMM_WORLD);
       MPI_Send(pbuf,position,MPI_PACKED,sl->mastRank,100,MPI_COMM_WORLD);  
       return;
   }
     

   int inRange(int fragid)
   {
     if(fragid>=g_startEST && fragid<=g_endEST) return 1;
     return 0;
   }

   void collectStringsForWork(struct slaveVar *sl,int *allfds,int *strIds)
   {
     int i=0,j=0;
     int jfid=0;
     long byteOffEST=0;
     char frag[MAXFRAGSIZE];
     struct pair *currWork;
     int W;
     char **workPairs;
     int wit=0;
	int fId; // input file identifier
	int flen;
     void getSeqtoMem(int fd, long byteOff,int bytes,char *frag);

     currWork=sl->currWork;
     W=sl->W;
     workPairs=sl->workPairs;

     if(W<=0) return ;


     j=0;
     for(i=0;i<W;i++) {
        	jfid=currWork[i].f1;
		if(g_TreeStringPtr[jfid]<=g_TreeStringsSize) {
			/* Link locally */
			workPairs[wit++]=&(g_TreeStrings[g_TreeStringPtr[jfid]]);
		}
		else {
			/* Go for IO */
			fId = AllESTlocator[jfid].fId;
           	byteOffEST=AllESTlocator[jfid].Off;
			flen=AllESTlen[jfid];

           	//getESTtoMem(allfps[fId],byteOffEST,frag);
           	getSeqtoMem(allfds[fId],byteOffEST,flen,frag);
           	workPairs[wit] = (char *)malloc(sizeof(char)*(strlen(frag)+1));
           	strcpy(workPairs[wit],frag);
			wit++;
		}
		strIds[j++]=jfid;

        	jfid=currWork[i].f2;
        	if(g_TreeStringPtr[jfid]<=g_TreeStringsSize) {
                /* Link locally */
                workPairs[wit++]=&(g_TreeStrings[g_TreeStringPtr[jfid]]);
        	}
        	else {
                /* Go for IO */
			fId = AllESTlocator[jfid].fId;
           	byteOffEST=AllESTlocator[jfid].Off;
           	//getESTtoMem(allfps[fId],byteOffEST,frag);
			flen=AllESTlen[jfid];
           	getSeqtoMem(allfds[fId],byteOffEST,flen,frag);
               workPairs[wit] = (char *)malloc(sizeof(char)*(strlen(frag)+1));
               strcpy(workPairs[wit],frag);
               wit++;
        	}
        	strIds[j++]=jfid;
     }
   } /* end collectStringsForWork */

   void cleanupStringsPostWork(struct slaveVar *sl,int *strIds)
   {
     int i;
     int count;

     count=2*sl->W;

     if(!strIds || count<=0) return ;

     for(i=0;i<count;i++)
     {
        if(g_TreeStringPtr[strIds[i]]<=g_TreeStringsSize) continue;

        if(sl->workPairs[i]) free(sl->workPairs[i]);
        sl->workPairs[i]=NULL;
     } 
     
   } /* end cleanupStringsPostWork */

   void initSlaveVar(struct slaveVar *sl)
   {
     sl->P = sl->R = sl->oE = sl->nE = sl->W =0;
     sl->terminationSent=NO;
     sl->gotit=0;
     sl->numPairsGenerated=0;
     sl->tIdleTime = sl->tSendTime = 0;
     sl->tsIdleTime = sl->tsSendTime = 0;
     sl->SBUF.pairs = 0;
   }

   void initSlave(struct slaveVar *sl)
   {
      struct pair dumpair;
      struct resPair dumrespair;
      char s[200];
      FILE *fileOpenForWrite(char *fname);
      int i=0;
	 void initBuf(struct cBuffer *,int );
	 int *OpenInputFileDescriptors(struct FileSize *gfs,int nFiles);
      

      #ifdef debug 
         /*printf("Rank=%d: initSlave() \n",rank);*/
      #endif
  
      LOADPERPROC = getLoadPerProc();
      #ifdef debug 
        if(rank==0) printf("Rank %d : Loadperproc=%d \n",rank,LOADPERPROC);
      #endif
      sl->param = getDynParams();
      Kfactor = sl->param.match*(1-sl->param.threshold)/(sl->param.match-sl->param.gap) ;
      mindep = retrieveMinLen();

      needs = sizeof(struct pair)*LOADPERPROC;
      sl->currWork = (struct pair *) malloc(needs);
      if(!sl->currWork) AlertErr("initSlave_currWork");
      
      needs = sizeof(struct resPair)*LOADPERPROC;
      sl->resToSend= (struct resPair *) malloc(needs);
      if(!sl->resToSend) AlertErr("initSlave_resToSend");

      needs = sizeof(char *)*2*LOADPERPROC;
      sl->workPairs = (char **) malloc(needs);
      if(!sl->workPairs) AlertErr("initSlave_workPairs");
      for(i=0;i<2*LOADPERPROC;i++) sl->workPairs[i]=NULL;

      needs = sizeof(int)*2*LOADPERPROC;
      sl->stringsForWork = (int *) malloc(needs);
      if(!sl->stringsForWork) AlertErr("initSlave_stringsForWork");

      sl->P=sl->R=sl->W=0;
      sl->oE=sl->nE=LOADPERPROC;
      initBuf(&sl->SBUF,MAXBUFSIZE);
      sl->terminationSent=NO;
      buildPairMPI(&dumpair,&pairtype);
      buildresPairMPI(&dumrespair,&respairtype);
      sl->send_req = MPI_REQUEST_NULL;
      sl->recv_req = MPI_REQUEST_NULL;
      
       getOutputFolder(g_sOutputFolder);
       if(rank==0) printf("Info: Output Folder is %s\n",g_sOutputFolder);

       g_bReportMaximalSubstrings= getReportMaximalSubstrings();
      if(g_bReportMaximalSubstrings) {
	         strcpy(s,"");
	   	 sprintf(s,"%s/maximal_substrings.%d.%d.PaCE",g_sOutputFolder,N,rank);
		 fpMaximalSubstrings=fileOpenForWrite(s);
      }

       g_bReportGeneratedPairs= getReportGeneratedPairs();
      if(g_bReportGeneratedPairs) {
	         strcpy(s,"");
	   	 sprintf(s,"%s/Gen_pairs.%d.%d.PaCE",g_sOutputFolder,N,rank);
		 fpGeneratedPairs=fileOpenForWrite(s);
      }

	 // open file pointer handles for all input FASTA files
	 //g_fp = OpenInputFilePointers(global_fs,nFilesToRead);
	 g_fd = OpenInputFileDescriptors(global_fs,nFilesToRead);
   } /* end init */



   void handlerS0(struct slaveVar *sl)
   {
      int i;
      struct dynResult result;
      int extractNextPairs(struct pair *batch,int k,struct cBuffer *SBUF);
      int generateNextPairsInBuf(struct cBuffer *SBUF,int k);

      #ifdef debug
         /*printf("Rank=%d : GOING TO EXTRACTED\n",rank);  fflush(stdout);*/
      #endif
      sl->W = extractNextPairs(sl->currWork,LOADPERPROC,&sl->SBUF); /* K1 */
      sl->numPairsGenerated=sl->numPairsGenerated + sl->W;

      collectStringsForWork(sl,g_fd,sl->stringsForWork);

      for(i=0;i<sl->W;i++)
      {
         result.typeAlign='X';

         processPair(&result,sl->currWork[i],sl->param,
			sl->workPairs[2*i],sl->workPairs[2*i+1]);
         if(result.typeAlign=='X')
         {
             reportErr("Process Pair failed !",sl);
             return;
         }
         sl->resToSend[i].maxBorderScore = result.maxBorderScore;
         sl->resToSend[i].idealBorderScore = result.idealBorderScore;
         sl->resToSend[i].maxScore = result.maxScore;
         sl->resToSend[i].idealScore = result.idealScore;
         sl->resToSend[i].maxScoreCoverage= result.coverage;
         sl->resToSend[i].f1 = sl->currWork[i].f1;
         sl->resToSend[i].f2 = sl->currWork[i].f2;
         sl->resToSend[i].type = result.typeAlign;
         sl->resToSend[i].atag = result.atag;
	 /* printf("Rank=%d: processed result   ",rank); fflush(stdout);
	 dispresPair(sl->resToSend[i]);*/
      } /* end for */

      cleanupStringsPostWork(sl,sl->stringsForWork);
      
      sl->R=sl->W;
      sl->W = extractNextPairs(sl->currWork,LOADPERPROC,&sl->SBUF);  /* K2 */
      sl->numPairsGenerated=sl->numPairsGenerated + sl->W;
      generateNextPairsInBuf(&sl->SBUF,LOADPERPROC);  /* K3 */
   } /* handlerS0 */

   void handlerS1(struct slaveVar *sl)
   {
      int extractNextPairs(struct pair *batch,int k,struct cBuffer *SBUF);

      sl->P = extractNextPairs(sl->newPairs,sl->nE,&sl->SBUF); 
      sl->numPairsGenerated=sl->numPairsGenerated + sl->P;
   }

   void handlerS2(struct slaveVar *sl)
   {
      int position=0;
      MPI_Status sendstat;
      /*double slot1,slot2;*/
	 struct timeval t_1,t_2;
      int i;
      
      /*slot1 = MPI_Wtime();*/
	 gettimeofday(&t_1,0);

	 if(!g_bMPI_Block_Sends)
      	MPI_Wait(&sl->send_req,&sendstat);

      position=0;
      MPI_Pack(&sl->P,1,MPI_INT,pbuf,PACKBUF,&position,MPI_COMM_WORLD);
      MPI_Pack(&sl->R,1,MPI_INT,pbuf,PACKBUF,&position,MPI_COMM_WORLD);
      if(sl->P >0) MPI_Pack(sl->newPairs,sl->P,pairtype,pbuf,PACKBUF,&position,MPI_COMM_WORLD);
      if(sl->R >0) MPI_Pack(sl->resToSend,sl->R,respairtype,pbuf,PACKBUF,&position,MPI_COMM_WORLD);

      /* Flush outstanding sends when running on BlueGene/L
         This will ensure message is sent immediately rather than
         made to wait until the Wait() call on BlueGene/L
      */
	 if(!g_bMPI_Block_Sends)
      	MPI_Isend(pbuf,position,MPI_PACKED,sl->mastRank,0,MPI_COMM_WORLD,&sl->send_req);  
	 else
     #ifdef Use_Ssend
      	MPI_Ssend(pbuf,position,MPI_PACKED,sl->mastRank,0,MPI_COMM_WORLD);  
     #else
      	MPI_Send(pbuf,position,MPI_PACKED,sl->mastRank,0,MPI_COMM_WORLD);  
     #endif

	 gettimeofday(&t_2,0);
	 sl->tSendTime+= (t_2.tv_sec-t_1.tv_sec)*1000000 + (t_2.tv_usec-t_1.tv_usec);
	 if(sl->tSendTime>=MAXULONGBY2) {
			sl->tsSendTime+=sl->tSendTime/1000000;
			sl->tSendTime=0;
	 }

      /*slot2 = MPI_Wtime();
      sl->tSendTime+=slot2-slot1;*/
      /*printf("Slave : %d Sent P=%d, R=%d to master\n",rank,sl->P,sl->R); fflush(stdout);
      for(i=0;i<sl->R;i++) {
	      printf("Slave=%d:  ",rank); dispresPair(sl->resToSend[i]); fflush(stdout);
      } */
   }  
   
   void handlerS3(struct slaveVar *sl)
   {
      handlerS2(sl);
   }

   void handlerS4(struct slaveVar *sl)
   {
      handlerS2(sl);
   }

   void handlerS5(struct slaveVar *sl)
   {
      handlerS2(sl);
   }

   void handlerS6(struct slaveVar *sl)
   {
      sl->R=0;
   }

   void handlerS7(struct slaveVar *sl )
   {
      int i;
      struct dynResult result;

      collectStringsForWork(sl,g_fd,sl->stringsForWork);

      for(i=0;i<sl->W;i++)
      {
         result.typeAlign='X';
         /*processPair(&result,sl->currWork[i],AllESTs,sl->param);*/
         processPair(&result,sl->currWork[i],sl->param,
			sl->workPairs[2*i],sl->workPairs[2*i+1]);
         if(result.typeAlign=='X')
         {
             reportErr("Process Pair failed !",sl);
             return;
         }
         sl->resToSend[i].maxBorderScore = result.maxBorderScore;
         sl->resToSend[i].idealBorderScore = result.idealBorderScore;
         sl->resToSend[i].maxScore = result.maxScore;
         sl->resToSend[i].idealScore = result.idealScore;
         sl->resToSend[i].maxScoreCoverage= result.coverage;
         sl->resToSend[i].f1 = sl->currWork[i].f1;
         sl->resToSend[i].f2 = sl->currWork[i].f2;
         sl->resToSend[i].type = result.typeAlign;
         sl->resToSend[i].atag = result.atag;
	 /*printf("Rank=%d: processed result  ",rank); fflush(stdout);
	 dispresPair(sl->resToSend[i]); */

      } /* end for */
      cleanupStringsPostWork(sl,sl->stringsForWork);
      sl->R = sl->W;
      /*printf("Slave : %d S7 processed %d\n",rank,sl->R); fflush(stdout);*/
   }

   void handlerS8(struct slaveVar *sl)
   {
      int position=0;
      /*double slot1,slot2; */
	 struct timeval t_1,t_2;
     
      /*slot1 = MPI_Wtime(); */
	 gettimeofday(&t_1,0);
      MPI_Recv(rpbuf,PACKBUF,MPI_PACKED,sl->mastRank,MPI_ANY_TAG,MPI_COMM_WORLD,&sl->recv_status);
      /*slot2 = MPI_Wtime();
      sl->tIdleTime+=slot2-slot1;*/

      gettimeofday(&t_2,0);
      sl->tIdleTime+= (t_2.tv_sec-t_1.tv_sec)*1000000 + (t_2.tv_usec-t_1.tv_usec);
      if(sl->tIdleTime>=MAXULONGBY2) {
               sl->tsIdleTime+=sl->tIdleTime/1000000;
               sl->tIdleTime=0;
      }

      if(sl->recv_status.MPI_TAG==0) sl->terminationSent=YES;
      else
      { 
         position=0;
         MPI_Unpack(rpbuf,PACKBUF,&position,&sl->W,1,MPI_INT,MPI_COMM_WORLD);
         MPI_Unpack(rpbuf,PACKBUF,&position,&sl->nE,1,MPI_INT,MPI_COMM_WORLD);
         MPI_Unpack(rpbuf,PACKBUF,&position,sl->currWork,sl->W,pairtype,MPI_COMM_WORLD);
         /* printf("Slave: %d Received W=%d, E=%d\n",rank,sl->W,sl->nE); fflush(stdout);*/
      }
   }  

   void handlerS9(struct slaveVar *sl)
   {
      MPI_Irecv(rpbuf,PACKBUF,MPI_PACKED,sl->mastRank,MPI_ANY_TAG,MPI_COMM_WORLD,&sl->recv_req);
   } 

   void handlerS10(struct slaveVar *sl)
   {
      MPI_Test(&sl->recv_req,&sl->gotit,&sl->recv_status);
   } 
    
   void handlerS11(struct slaveVar *sl)
   {
      int generateNextPairsInBuf(struct cBuffer *SBUF,int k);
      generateNextPairsInBuf(&sl->SBUF,sl->oE);  /* K3 */
   }

   void handlerS12(struct slaveVar *sl)
   {
      int position;
      
      if(sl->recv_status.MPI_TAG==0) sl->terminationSent=YES;
      else
      {
         position=0;
         MPI_Unpack(rpbuf,PACKBUF,&position,&sl->W,1,MPI_INT,MPI_COMM_WORLD);
         MPI_Unpack(rpbuf,PACKBUF,&position,&sl->nE,1,MPI_INT,MPI_COMM_WORLD);
         MPI_Unpack(rpbuf,PACKBUF,&position,sl->currWork,sl->W,pairtype,MPI_COMM_WORLD);
         sl->oE=sl->nE; 
          /*printf("Slave: %d Received W=%d, E=%d\n",rank,sl->W,sl->nE); fflush(stdout);*/
      }
   }

   void handlerS13(struct slaveVar *sl)
   {
      int tempSize=0,pSendSize=0,position=0;
      MPI_Status sendstat;
      /*double slot1,slot2;*/
	 struct timeval t_1,t_2;
	 void CloseInputFileDescriptors(int *fds,int nFiles);
      	
	 gettimeofday(&t_1,0);
      /*slot1 = MPI_Wtime();*/
      
	 if(!g_bMPI_Block_Sends) MPI_Wait(&sl->send_req,&sendstat);

      /* Can the following Send be removed ?? */
      position=0;
      MPI_Pack(&sl->P,1,MPI_INT,pbuf,PACKBUF,&position,MPI_COMM_WORLD);
      MPI_Pack(&sl->R,1,MPI_INT,pbuf,PACKBUF,&position,MPI_COMM_WORLD);
      if(sl->P >0) 
	 	MPI_Pack(sl->newPairs,sl->P,pairtype,pbuf,PACKBUF,&position,MPI_COMM_WORLD);
      if(sl->R >0) 
	 	MPI_Pack(sl->resToSend,sl->R,respairtype,pbuf,PACKBUF,&position,MPI_COMM_WORLD);

      if(sl->P>0 || sl->R>0) {
	printf("Rank=%d: Info: handlerS13: sl->P=%d and sl->R=%d\n",rank,sl->P,sl->R);
	fflush(stdout);
      }
      /*MPI_Isend(pbuf,position,MPI_PACKED,sl->mastRank,0,MPI_COMM_WORLD,&sl->send_req);  */

      /* Flush outstanding sends when running on BlueGene/L
         This will ensure message is sent immediately rather than
         made to wait until the Wait() call on BlueGene/L
      */
     #ifdef Use_Ssend
          MPI_Ssend(pbuf,position,MPI_PACKED,sl->mastRank,0,MPI_COMM_WORLD);
     #else
          MPI_Send(pbuf,position,MPI_PACKED,sl->mastRank,0,MPI_COMM_WORLD);
     #endif

      /*slot2 = MPI_Wtime();
      sl->tSendTime+=slot2-slot1;*/

      gettimeofday(&t_2,0);
      sl->tSendTime+= (t_2.tv_sec-t_1.tv_sec)*1000000 + (t_2.tv_usec-t_1.tv_usec);
      if(sl->tSendTime>=MAXULONGBY2) {
               sl->tsSendTime+=sl->tSendTime/1000000;
               sl->tSendTime=0;
      }
	 
	 //CloseInputFilePointers(g_fs,nFilesToRead);
	 CloseInputFileDescriptors(g_fd,nFilesToRead);
   }

   /*******************************************************************/
   /***********   START OF SLAVE CODE  *******************************/
void slave(int mastRank)
{
   int state=S0; /* current state of slave */
   double start,end,tDiff=0;
   struct slaveVar slv;

   /* printf("Slave %d: Start of Slave\n",rank); fflush(stdout); */
   start = MPI_Wtime();

   initSlaveVar(&slv);
   slv.mastRank=mastRank;
   initSlave(&slv);
   while(state!=SFINAL) 
   {
     /*if(rank==0 && state!=S10) {printf("Rank %d: STATE=S%d\n",rank,state); fflush(stdout);}*/
     switch (state)
     {
       case S0 : handlerS0(&slv);
    		 state=S1;
	 	 break;

       case S1 : handlerS1(&slv);
		 if(slv.P==0 && slv.R==0) state=S2;
		 if(slv.P>0 && slv.R==0) state=S3;
		 if(slv.P==0 && slv.R>0) state=S4;
		 if(slv.P>0 && slv.R>0) state=S5;
		 break;

       case S2 : handlerS2(&slv);
		 if(slv.W>0) state=S7;
                 else state=S6;
		 break;		 

       case S3 : handlerS3(&slv);
		 if(slv.W>0) state=S7;
                 else state=S6;
		 break;
      
       case S4 : handlerS4(&slv);
		 if(slv.W>0) state=S7;
                 else state=S6;
		 break; 

       case S5 : handlerS5(&slv);
		 if(slv.W>0) state=S7;
                 else state=S6;
	         break;

       case S6 : handlerS6(&slv);
		 if(slv.P > 0) state=S9;
      		 else state=S8;
		 break;

       case S7 : handlerS7(&slv);
		 if(slv.P > 0) state=S9;
      		 else state=S8;
		 break;
		
       case S8 : handlerS8(&slv);
                 if(slv.terminationSent==YES) state=S13;
                 else state=S1;
		 break;

       case S9 : handlerS9(&slv);
                 state=S10;
		 break;
    
       case S10: handlerS10(&slv);
                 if(slv.gotit) state=S12;
                 else if(!fullBuf(slv.SBUF) && (getBufSize(slv.SBUF)<MAXPAIRS) ) state=S11;
                 else state=S10;
		 break;
		  
       case S11: handlerS11(&slv);
                 state=S10;
                 break;
       
       case S12: handlerS12(&slv);
                 if(slv.terminationSent==YES) state=S13;
                 else state=S1;
                 break;

       case S13: handlerS13(&slv);
                 state=SFINAL;
                 break;

       default : printf("Error : Slave %d: Undefined State %d reached \n",rank,state);
	         fflush(stdout);
                 terminate();
	         return; 
     } /* end switch */
   } /* end while */

   free(slv.SBUF.pairs);
   free(slv.currWork);
   free(slv.resToSend);
   free(slv.workPairs);
   free(slv.stringsForWork);

   end = MPI_Wtime();
   tDiff = end - start;

   if(fpMaximalSubstrings) fclose(fpMaximalSubstrings);
   if(fpGeneratedPairs) {
        fprintf(fpGeneratedPairs,"%s",g_printpairstr);
	fclose(fpGeneratedPairs);
   }

     /*printf("rank: End of slave\n",rank); fflush(stdout);*/
     printf("Rank=%d: Clustering Phase IO stat: %u sequences (%u bytes) read: pread %u s \n",
		rank,n_callsToSeqtoMem,n_bytesReadFromIO,
		t_spread+(t_pread/1000000));
     printf("Rank=%d: Send Time := %u sec(s), Idle Time := %u sec(s)\n ",
			rank,slv.tsSendTime + slv.tSendTime/1000000, 
			slv.tsIdleTime + slv.tIdleTime/1000000);
     printf("Rank=%d: Total Time Taken for Master-Slave := %d sec(s)\n ",
			rank,(int) tDiff);
     printf("Rank=%d: Total Pairs generated = %d, Max Pairs per node=%d in node ID %d\n",
			rank,slv.numPairsGenerated,g_iMaxPairsPerNode,g_iMaxPairsNode);


   /***********   END OF SLAVE CODE  *******************************/
   /*******************************************************************/

} /* end slave */




/* Exception handler from the master proc */

void killAll()
{
  int i;
  MPI_Request send_req;
 

  for(i=1;i<p;i++)
  {
      MPI_Send(NULL,0,pairtype,i,0,MPI_COMM_WORLD);
  }
  MPI_Finalize();
  exit(0);
} /* end killAll */


/* readAllESTs : Reads the local portion of ests into a global array AllESTs */
void readAllESTs(char *estdata,int N,int k)
{
   int i=0;
   int f1;
   unsigned long ifsize=0;
   unsigned long iOffset=0;
   struct loadRes ldRes;
   int maxlocal=0;
   int *recv=NULL;
   int *recv_counts=NULL,*recv_disp=NULL,prev=0;
   //long *localESTlocator=NULL;
   struct FileOffset *localESTlocator=NULL;
   int *localESTlen=NULL;
   int mysted[2];

   // make these three global variables
   struct FileOffset *fo_boundaries=NULL;

   struct FileOffset *determineByteBoundaries(char *,struct FileSize **,int *);
   void dispLocator(long *locator,int n,char *datafile);
   void printIntArr(char *name,int *arr,int size);

   if(rank==0) {
	   printf("Total number of input sequences = %d\n",N);
	   fflush(stdout);
   }

   thisPhaseSpace=spaceSoFar;

   // assign this inside determineByteBoundaries
   // g_partSize = fabs(ifsize/p);

   /* Ananth: 6/25/2004 
    * Added to deterministically set the start and end fragment
    * read numbers without ambiguity 
    * Each proc sets its starting byte offset to the byte offset in
    * the file pointer value at which its corresponding new fasta
    * sequence starts.
    * It then communicates this value to the proc at its left
    * and sets the value it gets from its right to its end offset
    * */

   //determineByteBoundaries(estdata,&g_startOffset,&g_endOffset);
   fo_boundaries = determineByteBoundaries(estdata,&global_fs,&nFilesToRead);

   /*maxlocal=(3*N)/p;*/
   maxlocal=N;  
   /* Change 11/11/2003:  Cannot guarantee how many sequences will be uploaded locally */

   needs = sizeof(struct FileOffset)*(maxlocal);
   localESTlocator = (struct FileOffset *)malloc(maxlocal*sizeof(struct FileOffset));
   if(localESTlocator) { thisPhaseSpace+=needs; }
   else AlertErr("localESTlocator");

   needs = sizeof(int)*(maxlocal);
   localESTlen = (int *)calloc(maxlocal,sizeof(int));
   if(localESTlen) { thisPhaseSpace+=needs; }
   else AlertErr("localESTlen");

   needs = sizeof(struct est)*(maxlocal);
   g_localESTs = (struct est *)malloc(needs);
   /*if(g_localESTs) { thisPhaseSpace+=needs; spaceSoFar+=needs;}*/
   if(g_localESTs) { thisPhaseSpace+=needs; }
   else AlertErr("g_localESTs");
   
   //g_iNoLocalFrags = loadESTs(g_localESTs,N,estdata,g_startOffset,g_endOffset,
   //								k,localESTlocator,localESTlen,&ldRes);
   g_iNoLocalFrags = loadESTs(g_localESTs,N,estdata,fo_boundaries,global_fs,nFilesToRead,
   								k,localESTlocator,localESTlen,&ldRes);


   if(g_iNoLocalFrags>maxlocal) {
     printf("Rank=%d: Error: #Sequences loaded locally <%d> exceeds its limit <%d>\n",rank,g_iNoLocalFrags,maxlocal);
     fflush(stdout);
     terminate();
   }

   /* Change 11/11/2003: Reallocate all local storage to g_iNoLocalFrags*/
   localESTlocator = (struct FileOffset *)realloc(localESTlocator,
   								g_iNoLocalFrags*sizeof(struct FileOffset));
   localESTlen= (int *)realloc(localESTlen,g_iNoLocalFrags*sizeof(int));
   g_localESTs= (struct est *)realloc(g_localESTs,g_iNoLocalFrags*sizeof(struct est));
   spaceSoFar+=g_iNoLocalFrags*sizeof(struct est);
   
   if(!localESTlocator || !localESTlen || !g_localESTs) {
	  printf("Rank=%d: readAllESTs: Error: Reallocation error for localESTlocator or localESTlen or g_localESTs.\n",rank); 
	fflush(stdout);
	terminate();
   }


   /* Each proc has iNoLocalFrags frags in their local array
   Do an Allgather of this to calculate and store each frag's frag#
   accordingly
   A parallel prefix can also be done here : Planned
   */

   recv = (int *)malloc(sizeof(int)*p);
   MPI_Allgather(&g_iNoLocalFrags,1,MPI_INT,recv,1,MPI_INT,myComm);
   ldRes.startEST=computeStartFragId(recv);
   ldRes.endEST = ldRes.startEST + g_iNoLocalFrags - 1;

   if(sumUpAllFrags(recv)!=N) {
     printf("Warning: Rank =%d Calculated number of total fragments %d does NOT match supplied %d  <%d,%d>\n",
			rank,sumUpAllFrags(recv),N,ldRes.startEST,ldRes.endEST);
     fflush(stdout);
     terminate();
   }

   /* Gather all byte offset array from all procs */

   needs = sizeof(struct FileOffset)*N;
   AllESTlocator = (struct FileOffset *)malloc(N*sizeof(struct FileOffset));
   if(AllESTlocator) { thisPhaseSpace+=needs; spaceSoFar+=needs; }
   else AlertErr("AllESTlocator");

   recv_counts=(int *)malloc(sizeof(int)*p);
   for(i=0;i<p;i++)
   {
      recv_counts[i]=recv[i];
   }

   recv_disp=(int *)malloc(sizeof(int)*p);
   prev=0;
   for(i=0;i<p;i++)
   {
      recv_disp[i]=prev;
      prev+=recv[i];
   }

   MPI_Datatype mpidt_fo;
   struct FileOffset temp_fo;
	
   buildFileOffsetMPI(&temp_fo,&mpidt_fo);

   MPI_Allgatherv(localESTlocator,g_iNoLocalFrags,mpidt_fo,
   				AllESTlocator,recv_counts,recv_disp,mpidt_fo,
				myComm);
   /*dispLocator(AllESTlocator,N,estdata);*/

   /* Gather all len array from all procs */

   needs = sizeof(int)*(N);
   AllESTlen= (int *)calloc(N,sizeof(int));
   if(AllESTlen) { thisPhaseSpace+=needs; spaceSoFar+=needs; }
   else AlertErr("AllESTlen");

   MPI_Allgatherv(localESTlen,g_iNoLocalFrags,MPI_INT,
   				AllESTlen,recv_counts,recv_disp,MPI_INT,
				myComm);
   /*dispLocator(AllESTlocator,N,estdata);*/

   g_startEST = ldRes.startEST;
   g_endEST = ldRes.endEST;
   myFragsLen = ldRes.myFragsLen;
   g_inumLenLessK = ldRes.numLenLessK;
   g_isizLenLessK = ldRes.sizLenLessK;
   g_iLongest = ldRes.Longest;

   // this may overflow for more than 4 million sequences (assuming average length 1Kbp)
   g_TotalSeqLen = 0;
   g_AvgSeqLen = 0;
   unsigned int avgSeqLen_temp=0;
   for(i=0;i<N;i++) {
   	g_TotalSeqLen += AllESTlen[i];
	avgSeqLen_temp += AllESTlen[i];
	if(i%10000==0) { // for every 10,000 sequences update average (to prevent overflow)
		g_AvgSeqLen += avgSeqLen_temp/N;
		avgSeqLen_temp=0;
	}
   }
   g_AvgSeqLen += avgSeqLen_temp/N;
   //g_AvgSeqLen = g_TotalSeqLen/N;


   /* A2A: Need to gather the start and end sequences ranges for all processors */
   g_dataPartition=(int *)malloc(sizeof(int)*2*p);
   mysted[0]=g_startEST;
   mysted[1]=g_endEST;
   MPI_Allgather(mysted,2,MPI_INT,g_dataPartition,2,MPI_INT,myComm);
   /*printIntArr("g_dataPartition",g_dataPartition,2*p);*/


   if(recv!=NULL) free(recv);
   if(recv_counts!=NULL) free(recv_counts);
   if(recv_disp!=NULL) free(recv_disp);
   if(localESTlocator!=NULL) free(localESTlocator);
   if(localESTlen!=NULL) free(localESTlen);  

   if (fo_boundaries) free(fo_boundaries);

   #ifdef debug
     printf("Rank=%d : Longest Frag %d myFragsLen=%u\n",rank,g_iLongest,myFragsLen); 
     printf("Rank=%d : NumLess %d out of %d\n",rank,g_inumLenLessK,g_iNoLocalFrags); 
     printf("Rank=%d: START EST = %d, END EST = %d\n",rank,g_startEST,g_endEST);
	if(rank==0)
		printf("Rank=%d : Total Seq Len %u Average Seq Len %u\n",rank,g_TotalSeqLen,g_AvgSeqLen);

     fflush(stdout);
   #endif


} /* end readAllESTs */


/* getESTtoMem : gets the EST frag from disk to memory (frag)  */
/*******
void getESTtoMem(FILE *fp,unsigned long byteOff,char *frag)
{
  struct timeval t_1,t_2;


   if(byteOff<0) 
   {
      strcpy(frag,"");
      return;
   }

   n_callsToESTtoMem++;
   gettimeofday(&t_1,0);
   if(fseek(fp,byteOff,SEEK_SET)!=0)
   {
     printf("Rank=%d: Error : EST data file does not have the fragment with offset %d \n",rank,byteOff);
     fflush(stdout);
     terminate();
   }
   gettimeofday(&t_2,0);
   t_fseek+= (t_2.tv_sec-t_1.tv_sec)*1000000 + (t_2.tv_usec-t_1.tv_usec);
   if(t_fseek>=MAXULONGBY2) {
		t_sfseek+=t_fseek/1000000;
		t_fseek=0;
   }	

   fscanf(fp,"%[^\n]%\n",frag);
   gettimeofday(&t_1,0);
   t_fscanf+= (t_1.tv_sec-t_2.tv_sec)*1000000 + (t_1.tv_usec-t_2.tv_usec);
   if(t_fscanf>=MAXULONGBY2) {
		t_sfscanf+=t_fscanf/1000000;
		t_fscanf=0;
   }	
   n_bytesReadFromIO+=strlen(frag)+1;
}
**************/


void dispLocator(long *locator,int n,char *datafile)
{
  FILE *fp,*wfp;
  int i;
  char s[100],frag[MAXFRAGSIZE];

  sprintf(s,"%s.%d.frags",datafile,rank);
  wfp = fopen(s,"w");
  fp = fopen(datafile,"r");

  for(i=0;i<n;i++)
  {
   if(locator[i]<0)
   {
     fprintf(wfp,"%d:\n\n",i);
     continue;
   }

   if(fseek(fp,locator[i],SEEK_SET)!=0)
   {
     printf("Rank=%d: Error : EST data file does not have the fragment with offset %d \n",rank,locator[i]);
     fflush(stdout);
     terminate();
   }
   fscanf(fp,"%[^\n]%\n",frag);
   fprintf(wfp,"%d:\n%s\n",i,frag);
  }

  fclose(fp);
  fclose(wfp);
} /* end dispLocator */

int retrieveK()
{
   FILE *fp;
   char word[50];
   int k=10;
  
   fp = fopen(CFGFile,"r");
   while(!feof(fp))
   {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"window"))
    {
       fscanf(fp,"%s",word);
       fclose(fp);
       return atoi(word);
    }
   }
  fclose(fp);
  return k;
} /* end retrieveK */


int retrieveBatchSize()
{
   FILE *fp;
   char word[50];
   int k=10;
  
   fp = fopen(CFGFile,"r");
   while(!feof(fp))
   {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"BatchSize"))
    {
       fscanf(fp,"%s",word);
       fclose(fp);
       return atoi(word);
    }
   }
  fclose(fp);
  return k;
} /* end retrieveBatchSize */

int retrieveMinLen()
{
   FILE *fp;
   char word[50];
   int k=30;
  
   fp = fopen(CFGFile,"r");
   while(!feof(fp))
   {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"MinLen"))
    {
       fscanf(fp,"%s",word);
       fclose(fp);
       return atoi(word);
    }
   }
  fclose(fp);
  return k;
} /* end retrieveMinLen */

/* This function slides the window of length K over the local frags */
/* Local Frags Include : AllESTs[g_startEST.....g_endEST]  in the conceptual global array*/

void slideWindowK(int k)
{

  int i=0,iFragLen=0;
  struct est rFragi;
  
  void computeKeyAndStore(struct est fragi,int k,int id,char dir);

  thisPhaseSpace=spaceSoFar;

  g_iExpecKeys=(myFragsLen-g_isizLenLessK)*2 - (g_iNoLocalFrags-g_inumLenLessK)*2*(k-1);
  
  needs = sizeof(struct fkp)*g_iExpecKeys;
  fragkey = (struct fkp *) malloc(needs);
  if(fragkey) { thisPhaseSpace+=needs; spaceSoFar+=needs;}
  else AlertErr("fragkey");
  keysFree=needs; /* this will be freed later after bucketing */

  cmask32 = (1<<(intSize-1)) + (1<<(intSize-2));

  needs = sizeof(char)*(MAXFRAGSIZE+1);
  rFragi.s = (char *) malloc(needs);
  if(rFragi.s) { thisPhaseSpace+=needs; }
  else AlertErr("rFragi.s");

  for(i=0;i<g_iNoLocalFrags;i++)
  {

    iFragLen= strlen(g_localESTs[i].s);
    if(iFragLen<k)
    {   /* this means k is too big or this frag is too small! */
        continue;
    }
    /* Compute key, insert into a global array (key#,frag#,pos#) structures */
    computeKeyAndStore(g_localESTs[i],k,i+g_startEST,'D');


    /* Hash the Rev Comp too on the fly */
    strcpy(rFragi.s,"");
    getRComp(g_localESTs[i].s,rFragi.s);

    computeKeyAndStore(rFragi,k,i+g_startEST,'P');

  }  /* end for   */
 
  if(rFragi.s!=NULL) /* Ananth Change: for AiX/Sri 1/29/2004 */
    free(rFragi.s);
  /* A2A change:
   * Use this array for GST construction 
   *for(i=0;i<g_iNoLocalFrags;i++) {
	*  if(g_localESTs[i].s!=NULL) 
	*   free(g_localESTs[i].s);
    *}
    * if(g_localESTs!=NULL) free(g_localESTs); */

  /* Planned: subtract even the space allocated for the local sequences from spaceSoFar */
  /* Rochester:
   * spaceSoFar-=sizeof(struct est)*g_iNoLocalFrags;
   * spaceSoFar-=myFragsLen;
   * */

#ifdef debug
  /*printf(" Rank=%d: #FragKeys <Expected=%u,Generated=%u>   \n",rank,g_iExpecKeys,g_iFragKeys); */
#endif

    /* dispLocalFragKey();  */
} /* end slideWindowK() */


/* Computes and Stores the hash key for the fragment */
void computeKeyAndStore(struct est fragi,int k,int id,char dir)
{
    struct keyS thisKey;
    int iKey2;
    short j=0,iFirst,iLast,iFragLen;
    char frag[MAXFRAGSIZE];
    
    int Nposition=-1,currCharVal;

    
    iFragLen=strlen(fragi.s);

    sscanf(fragi.s,"%[A-Z]",frag);
    /* Compute key for the first window and then enter the loop*/
    thisKey.a=0;
    thisKey.b=0;
    j=0;
    while(j<k)
    {
	    /* Ananth: Nchange 24 Jun 2003 */
      currCharVal = getValue(frag[j++]);
      if(currCharVal==-1) { /* this char is a N or X  */
	Nposition = j-1;
	currCharVal=0;   /* N or X char value changed from -1 to 0, but i wont generate and store the key containing this char-pos*/
      }

      addToBase4(&thisKey,currCharVal); /* key*4 + currValue */
    }

    /* Now insert the first triplet*/
    if(Nposition<0) {
       fragkey[g_iFragKeys].fragid = id;
       fragkey[g_iFragKeys].dir = dir; 
       copyKey(&(fragkey[g_iFragKeys].key2),thisKey);
       fragkey[g_iFragKeys++].pos = 0;
    }
	    /* Ananth change end : Nchange 24 Jun 2003 */

    /* initialize first and last ptrs of the k-size window  */
    iFirst = 0;
    iLast = iFirst+k-1;
    iFragLen=strlen(frag);
    while(iLast<(iFragLen-1))
    {
      iFirst++;

	    /* Ananth: Nchange 24 Jun 2003 */
      currCharVal = getValue(frag[++iLast]);
      if(currCharVal==-1) { /* this char is a N  or X*/
		Nposition = iLast;
		currCharVal=0;   /* N or X char value changed from -1 to 0, but i wont generate and store the key containing this char-pos*/
      }
      iKey2 = currCharVal;
	    /* Ananth change end : Nchange 24 Jun 2003 */

      if(iKey2==-2) {
         printf("Warning: &%s,%c,%c&\n",frag,frag[iFirst-1],frag[iLast]);
	 iKey2=0;
      }

      /* iKey= ( (iKey-iKey1*iPowerK) << 2) + iKey2;   */
      addNewKey(&thisKey,iKey2,k);
    
      if(iFirst>Nposition) {
        fragkey[g_iFragKeys].fragid = id;
        fragkey[g_iFragKeys].dir = dir; 
        copyKey(&(fragkey[g_iFragKeys].key2),thisKey);
        fragkey[g_iFragKeys++].pos = iFirst;
		/*if(fragkey[g_iFragKeys-1].key2==0 && rank==0) {
			 printf("Rank=%d: Key is Zero at %d, id=%d, dir=%c,Nposition=%d\n",rank,iFirst,id,dir,Nposition);
		}*/
      }
      /*g_iFragKeys++; Ananth:  9 Sep 2003: Nchange bugfix*/
    }
    /* After the above execution l-k+1 more key triplets would have
     been inserted into the fragkey array due to this fragment
     where l is the length of frag#i
         */
} /* end computeKeyAndStore*/

 
int getValue(char c)
{
   double ran;
   int x,t;
   int lower=0,upper=3;


      switch(c)
      {
        case 'A' :  return 0;
        case 'C' :  return 1;
        case 'G' :  return 2;
        case 'T' :  return 3;
        case 'X' :      /*  X treated as N - Ananth: noN change */
        case 'N' :  return -1;
		    ran = drand48();
 		    x = (int) floor(ran*(upper-lower+1));
		    t = lower+x;
                    return t;
        case 'R' :  return 0;
        case 'Y' :  return 3;
        case 'K' :  return 2;
        case 'M' :  return 1;
        case 'S' :  return 2;
        case 'W' :  return 0;
        case 'B' :  return 2;
        case 'V' :  return 1;
        case 'D' :  return 2;
        case 'H' :  return 0;
        case 'U' :  return 3;
        default  :
                    printf("\nrank =%d : %c : %dError in decoding of Frag.\n",rank,c,c);
                    fflush(stdout);
                    return 0;
      }
      return -2;
} /* end getValue */

/* Copy KeyS Structure -the 2 32bit words of the key*/
void copyKey(unsigned int *b,struct keyS q)
{
  /**a = q.a; */
  *b = q.b;
} /* End copyKey() */

/* addToBase4 : Adds the base4 number new_val after
                shifting the existing p by one base 4 position
*/

void addToBase4(struct keyS *p,int new_val)
{
 unsigned int d;
 void what(unsigned int,unsigned int);

 d = (p->b) & cmask32;
 p->b=(p->b<<2) | new_val;
 d = d>>30;
 p->a=(p->a<<2) | d;
 /* what(p->a,p->b);*/

} /* End addToBase4() */

void what(unsigned int a,unsigned int b)
{
  printf("(%u*2P32)+%u\n",a,b);
} /* End what() */

/* addNewKey : Moves the sliding window k by one position to the right.
               This eliminates the left most base4 number
               and includes the next new_val right most base4 number
*/

void addNewKey(struct keyS *p,int new_val,int k)
{
  int l=0;
  unsigned int d;

  l=(k-1)<<1;

  if(k<=16)
  {
    p->b = ((p->b & ((1<<l)-1))<<2) | new_val;
  }
  else
  {
    l=l-32;
    p->a = ( p->a & ((1<<l)-1))<<2;
    d = (p->b & cmask32)>>30;
    p->b = (p->b<<2) | new_val;
    p->a = p->a | d;
  }

} /* End addNewKey() */


void dispLocalFragKey()
{
  int i;
  printf("\nrank =%d : ****START***%d*********\n",rank,g_iFragKeys);
  for(i=0;i<g_iFragKeys;i++)
  {
                   /*  key fragid   pos */
  printf("rank =%d : %u  %d  %c %d \n",rank,fragkey[i].key2,fragkey[i].fragid,fragkey[i].dir,fragkey[i].pos);
  }
  printf("\nrank =%d : *************END*****************\n",rank);
}

void dispfragkeys(struct fkp *s,int len)
{
  int i=0;
  printf("Rank : %d Number of Frag Keys =%d ** \n",rank,len);
  for(i=0;i<len;i++)
  {
     if(i%100 == 0) printf("%d %u <%d,%c> %d\n",i,s[i].key2,s[i].fragid,s[i].dir,s[i].pos);

  }
}  /* end dispfragkeys*/



/* buildTree : This parses the local sorted buckets of keys and
               builds the corresponding portions of GST in parts  
*/

void buildTree(int k)
{

  int i=0,j=0,bucsize=0;
  unsigned int k2;
  int *bucket;
 
  char frag[MAXFRAGSIZE];
  
  int start=0,end=0;
  struct stnode *bucroot;
  int distinct_strings=0,maxdist=0;
  int *buckStrings=NULL,buckit=0,jfid=0;
  long byteOffEST=0;
  FILE *fp;
  int er=0;
  struct fkp *buckeyes,*newfragkeys=NULL;
  unsigned int currBuckSize=0,maxBuckSize=0,g_maxBuckSize=0;
  int total_distinct_strings=0;
  unsigned int  gmaxBuckStrKeys=0;
  unsigned int maxBuckStrSize=0,currBuckStrSize=0,maxBuckStrKeys=0;


  int *send_buf=NULL,sendBufCap=0;
  int s=0,t=0;
  int n_A2A_rounds=0,n_A2A_iter=0;
  int mysted[2],upto=0;
  struct timeval t_beginA2A,t_endA2A;
  unsigned long t_A2A=0,t_sA2A=0;

  struct timeval t_start,t_end;
  struct timeval t_beginLoop,t_endLoop,t_endTemp;
  unsigned long t_ESTtoMemTimer=0,t_sESTtoMemTimer=0;
  unsigned long t_buildBucketsTimer=0,t_sbuildBucketsTimer=0;
  unsigned long t_reallocKeysBuckets=0;
  unsigned long t_whileLoopTimer=0;
  unsigned int u_currBuckStrSize;
  unsigned int start_index,end_index;
  unsigned int u_i,u_j,u_k;
  int bVoidedBucket;
  

  void doA2A(int *,int ,int ,int);

  thisPhaseSpace=spaceSoFar;
  g_treeSpace=0;



  /* first create N-size AllESTs array for ESTs */

  needs = sizeof(struct est)*N; 
  AllESTs = (struct est *) malloc(needs);
  if(AllESTs) { thisPhaseSpace+=needs;  }
  else AlertErr("AllESTs");
  
  for(i=0;i<N;i++)  AllESTs[i].s=NULL;
  mysted[0]=g_dataPartition[2*rank];
  mysted[1]=g_dataPartition[2*rank+1];
  for(i=mysted[0];i<=mysted[1];i++) AllESTs[i].s=g_localESTs[i-mysted[0]].s;
  
  /* Void out buckets with more than BucketSizeThreshold strings */

  //MaxStringsInABucket = getThreshold("MaxStringsInABucket");
  //g_BucketSizeThreshold = (unsigned int) (MaxStringsInABucket * (g_AvgSeqLen+1) );
  printf("Rank=%d: Subtree leaves threshold = %u bytes\n",rank,g_BucketSizeThreshold);
  //fflush(stdout);


	// new v10 - do it inplace without a temp_frag_keys

  u_currBuckStrSize = 0;
  start_index = 0;
  end_index = 0;
  u_k=0;
  for(u_i=0;u_i<g_iSortedKeys;u_i++) {
	if(u_i==0) k2=final_fragkeys[u_i].key2;
	if(k2!=final_fragkeys[u_i].key2) {
		if(u_currBuckStrSize < g_BucketSizeThreshold) {
			for(u_j=start_index;u_j<=end_index;u_j++) {
				if(u_k<u_j) copyFKP(&(final_fragkeys[u_k++]),final_fragkeys[u_j]);
				else {
					if(u_k>u_j) {
						printf("Serious Warning: Something wrong while inplace copying\n");
						fflush(stdout);
					}
				}
			}
		}
		else {
			printf("Rank=%d: Skipping bucket <%u, %u> with %u \n",rank,start_index,end_index,u_currBuckStrSize);
		}
		start_index = u_i;
		end_index = u_i;
		u_currBuckStrSize = 0;
	}
	k2=final_fragkeys[u_i].key2;
	u_currBuckStrSize+=(unsigned int) AllESTlen[final_fragkeys[u_i].fragid]+1;
	end_index=u_i;
  }
  if(u_currBuckStrSize < g_BucketSizeThreshold) {
	for(u_j=start_index;u_j<=end_index;u_j++) {
		if(u_k<u_j) copyFKP(&(final_fragkeys[u_k++]),final_fragkeys[u_j]);
		else {
			if(u_k>u_j) {
				printf("Serious Warning: Something wrong while inplace copying\n");
				fflush(stdout);
			}
		}
	}
  }

  printf("Rank=%d: Adjusted SortedKeys to %u (from %u) \n",rank,u_k,g_iSortedKeys);
  g_iSortedKeys = u_k;
  if(g_iSortedKeys>0) {
  	final_fragkeys = (struct fkp *) realloc(final_fragkeys,sizeof(struct fkp)*g_iSortedKeys);
  	if(!final_fragkeys) AlertErr("final_fragkeys_adjusted");
  }

 
  /*  Ananth: 19 Oct 2003 - New Space Saving Scheme:  
   *  (Re)Allocate final_fragkeys to 2*g_iSortedKeys entries and then
   *  Move the present set of keys to the second half
   *  Initialize Nodes pointer to the index 0, and
   *  Initialize final_fragkeys to the index g_iSortedKeys
   *  The main idea is to keep creating nodes in the local tree (DFS) from index 0,
   *  as we replenish the keys from left to right as buckets
   *  This method takes advantage of the fact that the struct stnode and struct keys
   *  take 12 bytes each.  Thus the method is subject to change if any of the above 
   *  structures change.
   * */
    
  g_currMaxNodes=2*g_iSortedKeys;
  #ifdef debug
    printf("Rank=%d: Initial expected Number of nodes and leaves %u \n",rank,g_currMaxNodes);
    fflush(stdout);
  #endif
  
  gettimeofday(&t_start,0);

    if(g_currMaxNodes>0) {
  	needs = sizeof(struct fkp )*g_currMaxNodes; 
  	final_fragkeys= (struct fkp *) realloc(final_fragkeys,needs);
  	if(final_fragkeys) { thisPhaseSpace+=needs/2; }
  	else AlertErr("final_fragkeys-Nodes");
    }
    else {
	    printf("Rank=%d: Warning! Local GST empty !!\n",rank);
	    fflush(stdout);
	    final_fragkeys=NULL;
    }
  newfragkeys=&(final_fragkeys[g_currMaxNodes/2]);
 
  needs=sizeof(char)*N;
  g_TreeStringMarker = (char *)malloc(needs);
  if(g_TreeStringMarker) {
	thisPhaseSpace+=needs;
	spaceSoFar+=needs;
  }
  else AlertErr("buildTree_g_TreeStringMarker");

  for(i=0;i<N;i++) g_TreeStringMarker[i]='0';
  

  maxBuckSize=currBuckSize=0;
  maxBuckStrSize=currBuckStrSize=0;
  maxBuckStrKeys=0;
  j=0;
  for(i=0;i<g_iSortedKeys;i++) {

	g_TreeStringMarker[final_fragkeys[i].fragid]='1';

	if(i==0) k2=final_fragkeys[i].key2;
	if(k2!=final_fragkeys[i].key2) {
		  if(currBuckSize>maxBuckSize) maxBuckSize=currBuckSize;
		  if(currBuckStrSize>maxBuckStrSize) {
			  	maxBuckStrSize=currBuckStrSize;
				maxBuckStrKeys=currBuckSize;
		  }
		  currBuckSize=0;
		  currBuckStrSize=0;
	}
	currBuckSize++;
	k2=final_fragkeys[i].key2;
	if(final_fragkeys[i].fragid<0) {
		printf("Rank=%d: Error: Fragid is <0 in %d\n",rank,i);
		fflush(stdout);
		terminate();
	}
	currBuckStrSize+=(unsigned int) AllESTlen[final_fragkeys[i].fragid]+1;
  	copyFKP(&(newfragkeys[j++]),final_fragkeys[i]);
  }
  if(currBuckSize>maxBuckSize) maxBuckSize=currBuckSize;
  if(currBuckStrSize>maxBuckStrSize) {
			  	maxBuckStrSize=currBuckStrSize;
				maxBuckStrKeys=currBuckSize;
  }

  assert(sizeof(struct stnode) < sizeof(struct fkp));

  Nodes = (struct stnode *) final_fragkeys; 
  final_fragkeys=&(final_fragkeys[g_currMaxNodes/2]);

  gettimeofday(&t_end,0);
  t_reallocKeysBuckets+= (t_end.tv_sec-t_start.tv_sec)*1000000 + (t_end.tv_usec-t_start.tv_usec);
  if(checkEndian()==LSBinvalid)  terminate();

  MPI_Allreduce(&maxBuckSize,&g_maxBuckSize,1,MPI_UNSIGNED,MPI_MAX,myComm);
  MPI_Allreduce(&maxBuckStrKeys,&gmaxBuckStrKeys,1,MPI_UNSIGNED,MPI_MAX,myComm);
  MPI_Allreduce(&maxBuckStrSize,&g_maxBuckStrSize,1,MPI_UNSIGNED,MPI_MAX,myComm);

  #ifdef debug
    printf("Rank=%d: Local vs. Global maxBuckStrKeys %u : %u\n",rank,maxBuckStrKeys,gmaxBuckStrKeys);
    printf("Rank=%d: Local vs. Global maxBuckStrSize %u : %u\n",rank,maxBuckStrSize,g_maxBuckStrSize);
    fflush(stdout);
  #endif

  /*n_A2A_iter=determineA2AIterations(g_maxBuckSize);*/
  n_A2A_iter=determineA2AIterations(gmaxBuckStrKeys);
  #ifdef debug
    if(rank==0) {
	printf("Rank=%d: Planned number of A2A rounds = %d\n",rank,n_A2A_iter);
    	fflush(stdout);
    }
  #endif
  sendBufCap=N;
  if(g_maxBuckSize<sendBufCap) sendBufCap=g_maxBuckSize;
  if(sendBufCap) {
  	  needs = sizeof(int)*sendBufCap;
	  send_buf = (int *)malloc(sizeof(int)*sendBufCap);
	  if(!send_buf) AlertErr("send_buf");
	  else thisPhaseSpace+=needs;

  }
  else {
	  send_buf = (int *)malloc(sizeof(int)*p);
  }

  g_iNodes=0;
  n_A2A_rounds=0;

  
  gettimeofday(&t_beginLoop,0);
  i=0; 
  while(i<g_iSortedKeys) {
     distinct_strings=0;
	g_TransientMemory=0;

     start=i;
     k2 = final_fragkeys[i].key2;
     j=i+1;
     while(k2==final_fragkeys[j].key2) {
        j++;
        if(j<g_iSortedKeys) continue;
        break;
     }
     end=j-1;
     bucsize=end-start+1;

     needs = sizeof(struct fkp)*bucsize;
     buckeyes = (struct fkp *) malloc(needs);
     if(!buckeyes) AlertErr("buckeyes");
	g_TransientMemory+=needs;

     t=0;
     for(j=start;j<=end;j++) {
	     copyFKP(&(buckeyes[t++]),final_fragkeys[j]);
     }


     /* A2A change:
 	* Get strings into local memory through communication instead of IO*/
	 if(start==g_itStarters[n_A2A_rounds]) {
		n_A2A_rounds++;
		upto=g_iSortedKeys-1;
		if(n_A2A_rounds < n_A2A_iter) upto=g_itStarters[n_A2A_rounds];
		if(upto<0) { /* last round */ 
			upto=g_iSortedKeys-1;
		}
		gettimeofday(&t_beginA2A,0);

		doA2A(send_buf,start,upto,n_A2A_rounds-1);
		gettimeofday(&t_endA2A,0);
     		t_A2A+= (t_endA2A.tv_sec-t_beginA2A.tv_sec)*1000000 + (t_endA2A.tv_usec-t_beginA2A.tv_usec);
     		if(t_A2A>=MAXULONGBY2) {
	     		t_sA2A+=t_A2A/1000000;
	     		t_A2A=0;
     		}
	 } /* end if */

     needs=bucsize*sizeof(int);
     bucket=(int *)calloc(bucsize,sizeof(int));
	if(!bucket) AlertErr("BuildTree_bucket");
	g_TransientMemory+=needs;

     for(t=0;t<bucsize-1;t++) { bucket[t]=t+1; }
     bucket[t]=-1;
      
     /* SpaceSaver: create node on demand (realloc feature) */
     /* bucroot = createNode(&Nodes,&g_currMaxNodes,g_iNodes);*/
     bucroot = &(Nodes[g_iNodes]);

     buildNode(bucroot,k);
     g_iNodes++;
/*printf("Rank=%d: buildBucket %d to %d\n",rank,0,bucsize-1);
fflush(stdout);*/

     gettimeofday(&t_start,0);
     buildBucket(0,bucsize-1,0,bucket,bucroot,buckeyes);
     gettimeofday(&t_end,0);
/*printf("Rank=%d: buildBucket %d to %d over\n",rank,0,bucsize-1);
fflush(stdout);*/
     t_buildBucketsTimer+= (t_end.tv_sec-t_start.tv_sec)*1000000 + (t_end.tv_usec-t_start.tv_usec);
     if(t_buildBucketsTimer>=MAXULONGBY2) {
	     t_sbuildBucketsTimer+=t_buildBucketsTimer/1000000;
	     t_buildBucketsTimer=0;
     }
   
     numbuckets++;
     i=end+1;

     if(bucket!=NULL) free(bucket);
     if(buckeyes!=NULL) free(buckeyes);
	g_TransientMemory=0;
  
  } /* end while */


  for(i=n_A2A_rounds;i<n_A2A_iter;i++) {
		gettimeofday(&t_beginA2A,0);
		doA2A(send_buf,0,-1,i);
		gettimeofday(&t_endA2A,0);
     		t_A2A+= (t_endA2A.tv_sec-t_beginA2A.tv_sec)*1000000 + (t_endA2A.tv_usec-t_beginA2A.tv_usec);
     		if(t_A2A>=MAXULONGBY2) {
	     		t_sA2A+=t_A2A/1000000;
	     		t_A2A=0;
     		}
  }
  t_sA2A+=t_A2A/1000000;
  gettimeofday(&t_endLoop,0);
  t_whileLoopTimer= (t_endLoop.tv_sec-t_beginLoop.tv_sec)+ (t_endLoop.tv_usec-t_beginLoop.tv_usec)/1000000;

  g_iGSTscratchSpace=0;
		
  if(send_buf) free(send_buf);
  if(g_recv_str) free(g_recv_str);
  if(g_itStarters) free(g_itStarters);

  /* update with residual timers */
  t_sESTtoMemTimer+=t_ESTtoMemTimer/1000000;
  t_sbuildBucketsTimer+=t_buildBucketsTimer/1000000;

  /* Nodes = (struct stnode *)realloc(Nodes,sizeof(struct stnode)*g_iNodes);*/
  spaceSoFar+=sizeof(struct stnode)*g_iNodes;
  g_treeSpace+=sizeof(struct stnode)*g_iNodes;

   for(i=0;i<N;i++)  {
	if(i<mysted[0] || i>mysted[1]) AllESTs[i].s=NULL;
   }

   for(i=0;i<g_iNoLocalFrags;i++)  {
   	if(g_localESTs[i].s) free(g_localESTs[i].s); 
        g_localESTs[i].s=NULL;
   }
   if(g_localESTs) free(g_localESTs);

   if(AllESTs) free(AllESTs);

   needs=sizeof(int)*2*N;
   g_NodeMarker = (int *)malloc(needs);
   if(g_NodeMarker) { spaceSoFar+=needs;}
   else AlertErr("g_NodeMarker");

  for(i=0;i<2*N;i++) g_NodeMarker[i]=-1;  /* initialize */

  /*if(rank==3) dispNodesDetails();   */
  #ifdef debug 
     printf("Rank=%d: Work done =%u bytes\n",rank,g_iNumCharsVisited);
     printf("Rank=%d: Max distinct strings per bucket=%d\n",rank,maxdist);
     printf("Rank=%d: Time for A2A %u s \n",rank,t_sA2A);
     if(t_sESTtoMemTimer>0) {

     }
     else {

     }
     printf("Rank=%d: Time for building all buckets =%u s\n",rank,t_sbuildBucketsTimer);
     printf("Rank=%d: Total number of distinct strings = %d\n",rank,total_distinct_strings);
     printf("Rank=%d: Time for the build tree while loop=%u s\n",rank,t_whileLoopTimer);
  #endif

#ifdef debug
  printf("Rank=%d: Number of Buckets=%d,Nodes=%u,Leaves=%u\n",rank,numbuckets,g_iNodes,g_iLeaves); 
  printf("Rank=%d: Tree Space taken (in bytes) = %u, (nodes=%u, lset=%u)\n",rank,g_treeSpace,sizeof(struct stnode)*g_iNodes,g_treeSpace-sizeof(struct stnode)*g_iNodes);
  printf("Rank=%d: Tree Space Peak Usage (in bytes) = %d\n",rank,thisPhaseSpace);
  
  fflush(stdout); 
#endif


} /* end buildTree*/





/* localSortNodes : Locally sort all tbe nodes depthwise in the decreasing
                    order 
*/

void localSortNodes()
{
    int i;
    /*struct depnode *deps;*/
    

    thisPhaseSpace=spaceSoFar;

    /*
    * needs = sizeof(struct depnode)*g_iNodes;
    * deps = (struct depnode *) malloc(needs);
    * if(deps) { thisPhaseSpace+=needs;}
    * else AlertErr("deps"); 
    */

    /*
    * needs = sizeof(short)*g_iNodes;
    * depthArray=(short *)malloc(needs);
    * if(depthArray) { thisPhaseSpace+=needs; }
    * else AlertErr("depthArray"); 
    */

    needs = sizeof(int)*g_iNodes;
    if(g_iNodes>0) {
	    depnodes = (int *) malloc(needs);
    	    if(depnodes) { thisPhaseSpace+=needs; spaceSoFar+=needs;}
    	    else AlertErr("depnodes"); 
    } else {
	    depnodes=NULL;
    }

    /* prepare array of depth for sorting next phase */
    for(i=0;i<g_iNodes;i++)
    {
      depnodes[i] = i;
      /* depthArray[i] = Nodes[i].depth;*/
    }  /* end for */   

  
    /* quicksortNodes(0,g_iNodes-1,deps);*/
    nodesort(depnodes,Nodes,g_iNodes);

    /*for(i=0;i<g_iNodes;i++)
    * {
    *  depnodes[i]=deps[i].id;
    * }*/

    /*if(depthArray!=NULL) free(depthArray);*/

} /* end localSortNodes */


int extractNextPairs(struct pair *batch,int k,struct cBuffer *SBUF)
{
   int i=0,j=0;
   static int outOfPairs=0;
   int generateNextPairsInBuf(struct cBuffer *SBUF,int k);

   if(outOfPairs) return 0;
   while(i<k)
   {
      if(emptyBuf(*SBUF))
      {
	    j = generateNextPairsInBuf(SBUF,k-i);
         if(j<=0) {
             outOfPairs=1;
             break;
         }
      }
      deQBuf(SBUF,&(batch[i]));
      i++;
   }
   return i;
} /* end extractNextPairs */



/* generateNextPairsInBuf : generates next set of k or the pairs of a single processed node
                        whichever is minimum.  
		       Returns the number of pairs fed in the SBUF
                       returns -1 on error 
*/

int generateNextPairsInBuf(struct cBuffer *SBUF,int k)
{
  int j=0,generated=0,r=0; 
  int done=0;
  int recalK=0;

  /* static variables */
  static int i=0,bFirst=1,outOfPairs=0;
  static unsigned int s_pairsCount=0;
  

  if(outOfPairs) return 0;

  if(bFirst)
  {
     /* use depnodes and Nodes arrays to process in descending order of depth */
     bFirst=0;
     i=g_iNodes-1;
  } 

  if( (getBufSize(*SBUF)+k)>=MAXPAIRS )
  {  /* prevent overrun of the buffer */
     recalK=MAXPAIRS-getBufSize(*SBUF);
  }
  else
  {
     recalK=k;
  }
  
  if(recalK<1) return 0;  /* Newly added : Ananth 05/03/2003 */

  while(i>=0)
  {
     done=0;
     j=depnodes[i];
     lastNode=j;
  
     if((Nodes[j].depth)<mindep) 
     {
       outOfPairs=1; 
       break;
     }

     /* printf("Rank=%d: Node being processed =%d %d\n",rank,j,j==Nodes[j].rmost); fflush(stdout);  */

     if(j==Nodes[j].rmost) 
     {
       r=processLeaf(j,SBUF,recalK-generated,&done);
       if(r<0)
       {
         printf("Slave: %d Error generateNext returning -1 \n",rank); fflush(stdout);
         return -1;
       }
     }
     else
     {
       r=processNode(j,SBUF,recalK-generated,&done);
       if(r<0)
       {
         printf("Slave: %d Error generateNext returning -1 \n",rank); fflush(stdout);
         return -1;
       }
     }

     s_pairsCount+=r;
     if(done) {
     	if(g_iMaxPairsPerNode<s_pairsCount) {
		g_iMaxPairsPerNode=s_pairsCount;
		g_iMaxPairsNode=j;
		s_pairsCount=0;
     	}
        i--;
     }

     generated+=r;

     /* if(r>0) { printf("Rank=%d: Node #Pairs %d\n",rank,r); fflush(stdout); }*/
     if(generated>=recalK) break;

  }  /* end while */
  return generated;
} /* end generateNextPairsInBuf */


void bucketing(struct fkp *fkey, int size,int k)
{
   int i=0;
      unsigned int it=0,it2=0;
   unsigned int j=0;
   int totBucks=0;
   unsigned int *buck=NULL;
   unsigned int *rbuck=NULL,totKeys=0,myTotKeys=0;
   int zone=0;
   unsigned int parti=0;
   int start=0,end=0;
   unsigned int temp=0;
   int mysted[2],*sted;
   int *send_count=NULL,*recv_count=NULL;
   struct fkp *send_buf=NULL,*recv_buf=NULL,n1;
   int *send_disp=NULL,*recv_disp=NULL,*recvG=NULL;
   int disp,r=0,s=0;
   unsigned int iGoingToRecv=0;
   MPI_Datatype sFKPtype;
   struct bucketkeys *buckeys=NULL,*bky=NULL,*bkeys2=NULL;
   int curred,len=0,t=0;
   int sumbkeysspace=0,sumbkeysnext=0;
    unsigned   int max_buck_keys=0; 
         int max_buck=-1;
   int MaxStringsInABucket;
  
   char bsize_file[200];

   void dispLocalBucks(int *buck,int totBuc);

   thisPhaseSpace=spaceSoFar;

   if(k>16)
   {
        printf("Error : Change bucketing code to work for key1 and key2 fields \n");
       return;
   }

   totBucks = 1<<(2*k);  /* this is 4 the power k */

   needs = sizeof(unsigned int)*totBucks;
   buck = (unsigned int *)malloc(needs);
   if(buck) { thisPhaseSpace+=needs; /*spaceSoFar+=needs;*/  
              #ifdef debug 
                    /*printf("Rank:=%d Inside bucketing: buck=%d\n",rank,needs); */
              #endif
            }
   else AlertErr("buck");

   needs = sizeof(unsigned int)*totBucks;
   rbuck = (unsigned int *)malloc(needs);
   if(rbuck) {  thisPhaseSpace+=needs; /* spaceSoFar+=needs; */
              #ifdef debug 
                /*printf("Rank:=%d Inside bucketing: rbuck=%d\n",rank,needs); */
              #endif
             }
   else AlertErr("rbuck");

   for(i=0;i<totBucks;i++)
   {
      buck[i]=0;
      rbuck[i]=0;
   }

   for(i=0;i<size;i++)
   {
      j=fkey[i].key2;
      buck[j]++ ;
   }

   MPI_Allreduce(buck,rbuck,totBucks,MPI_UNSIGNED,MPI_SUM,myComm);

   /* if(rank==0) dispLocalBucks(rbuck,totBucks); */
   /*sprintf(bsize_file,"bsize.%d",rank);
   bsize_fp=fopen(bsize_file,"w"); */

  /* Void out buckets with more than BucketSizeThreshold strings */

  MaxStringsInABucket = getThreshold("MaxStringsInABucket");
  g_BucketSizeThreshold = (unsigned int) (MaxStringsInABucket * (g_AvgSeqLen+1) );
  printf("Rank=%d: Bucket Size Threshold = %u bytes\n",rank,g_BucketSizeThreshold);
  fflush(stdout);

  char *OverFullBucket=NULL;
  
  OverFullBucket = (char *)malloc(sizeof(char)*totBucks);
   for(i=0;i<totBucks;i++) OverFullBucket[i] = '0';


   totKeys=0;
   parti=0;
   for(i=0;i<totBucks;i++)
   {
   	if(rbuck[i]>=g_BucketSizeThreshold) {
		OverFullBucket[i] = '1';
		printf("Rank=%d: OverFullBucket[%u] bucket\n",rank,i);
		continue;
	}
	if(totKeys>=1000000) {
		parti+=totKeys/p; // running average
		totKeys=0;
	}
     totKeys+=rbuck[i];
      if(max_buck_keys<rbuck[i]) {
            max_buck_keys=rbuck[i];
               max_buck=i;
         }
   }
   parti+=totKeys/p;
      //if(rank==0) {
         printf("Rank=%d: max_buck_keys=%u , max_buck=%d\n",rank,max_buck_keys,max_buck);
      //}

   #ifdef debug
      //printf("Rank=%d totKeys=%u, parti=%u\n",rank,totKeys,parti);
      printf("Rank=%d parti=%u\n",rank,parti);
   #endif
   fflush(stdout);
   zone = 0;
   temp =0;
   i=0;
   while(zone<rank && i<totBucks)
   {
   	if(OverFullBucket[i]=='0')  {
     	temp+=rbuck[i];
     	if(temp>=parti) {
        		zone++;
        		temp=0;
		}
     }
     i++;
   }
   start = i;
   myTotKeys=0;
   i=start;
   while(i<totBucks)
   {
   	if(OverFullBucket[i]=='0') myTotKeys+=rbuck[i];
     /*fprintf(bsize_fp,"%d	%d\n",i,rbuck[i]);*/
	if(rank==1022) {
     	//printf("Rank=%d: rbuck	%d	%d\n",rank,i,rbuck[i]);
	}
     i++;
     if( rank<(p-1) && myTotKeys>=parti)  break;
   }

   end = i-1;
   #ifdef debug
   /*printf("Rank=%d: start=%d, end=%d, myTotKeys=%d \n",rank,start,end,myTotKeys);
   fflush(stdout);*/
   #endif

   mysted[0]=start;
   mysted[1]=end;

   needs = sizeof(int)*p*2;
   sted = (int *)malloc(needs);
   if(sted) { thisPhaseSpace+=needs; /* spaceSoFar+=needs;*/}
   else AlertErr("sted");

   MPI_Allgather(mysted,2,MPI_INT,sted,2,MPI_INT,myComm);

if(rank==0)
   for(i=0;i<(2*p);i=i+2)
   {
      printf("Rank=%d : bucket start=%d, end=%d\n",i/2,sted[i],sted[i+1]);
         fflush(stdout);
   }

   needs = sizeof(struct bucketkeys)*totBucks;
   buckeys= (struct bucketkeys *)malloc(needs);
   if(buckeys) { thisPhaseSpace+=needs; /* spaceSoFar+=needs;*/  
              #ifdef debug
                 /*printf("Rank:=%d Inside bucketing: buckeys=%d\n",rank,needs); */
              #endif
             }
   else AlertErr("buckeys");
  
   sumbkeysspace=0;
   for(i=0;i<totBucks;i++)
   {
   	 if(OverFullBucket[i]=='1' || buck[i] <= 0)
      {
         buckeys[i].bk=NULL;
         buckeys[i].it=0;
         continue;
      }

      needs = sizeof(struct fkp)*buck[i];
      if(buck[i]<=0) buckeys[i].bk=NULL;
      else {
      	buckeys[i].bk = (struct fkp *)malloc(needs);
      	if(buckeys[i].bk) { thisPhaseSpace+=needs; /* spaceSoFar+=needs;*/  sumbkeysspace+=needs; }
      	else AlertErr("buckeys[i].bk");
      }
      buckeys[i].it=0;
   }

   #ifdef debug
      /* printf("Rank:=%d Inside bucketing: buckeys[].bk(sum)=%d\n",rank,sumbkeysspace); */
   #endif

   for(i=0;i<size;i++)
   {
      j=fkey[i].key2;
	 if(OverFullBucket[j]=='1') continue;
      bky = &(buckeys[j]);
	 if(!(bky->bk)) continue;
      copyFKP(&(bky->bk[bky->it]),fkey[i]);
      (bky->it)++;
   }

   /* get ready for an All-To-AllV communication */
   needs = sizeof(int)*p;
   send_count = (int *) malloc(needs);
   if(send_count) { thisPhaseSpace+=needs; /* spaceSoFar+=needs; */}
   else AlertErr("send_count");

   needs = sizeof(int)*p;
   recv_count = (int *) malloc(needs);
   if(recv_count) { thisPhaseSpace+=needs; /* spaceSoFar+=needs; */}
   else AlertErr("recv_count");


   /* SpaceSaver: just re-use or overwrite fkey for send_buf */
   send_buf = fkey;

 /*  needs = sizeof(struct fkp)*size;
   send_buf = (struct fkp *) malloc(needs);
   if(send_buf) { thisPhaseSpace+=needs; 
                  #ifdef debug
                     printf("Rank:=%d Inside bucketing: send_buf=%d\n",rank,needs); 
                  #endif 
                }
   else AlertErr("send_buf"); */

   needs = sizeof(int)*p;
   send_disp = (int *) malloc(needs);
   if(send_disp) { thisPhaseSpace+=needs; /* spaceSoFar+=needs; */}
   else AlertErr("send_disp");

   needs = sizeof(int)*p;
   recv_disp = (int *) malloc(needs);
   if(recv_disp) { thisPhaseSpace+=needs; /* spaceSoFar+=needs; */}
   else AlertErr("recv_disp");

   for(i=0;i<p;i++) {
      send_count[i] = 0;
      recv_count[i] = 0;
   }

   /* Fill up the send_buf     */
   curred=sted[1];
   r=1;   /* for sted */
   s=0; /*  for send_count */
   k=0;  /* for the send_buf*/
   for(i=0;i<totBucks;i++)
   {
      if(i>curred)
      {
         s++;
         r=r+2;
         curred=sted[r];
      }
	 if(OverFullBucket[i]=='1') continue;

      len = buck[i];
      if(len<=0) continue;
      bky = &(buckeys[i]);
      for(t=0;t<len;t++)
      {
        send_count[s]++;
        copyFKP(&send_buf[k++],bky->bk[t]);
      }

   }  /* end for */

   /*  SpaceSaver: Free buckeys here */ 
   for(i=0;i<totBucks;i++)
   {
     if(buckeys[i].bk!=NULL) free(buckeys[i].bk);
   }
   if(buckeys!=NULL) /* Ananth Change: for AiX/Sri 1/29/2004 */
     free(buckeys);

   if(buck!=NULL) /* Ananth Change: for AiX/Sri 1/29/2004 */
     free(buck);

   if(OverFullBucket) free(OverFullBucket);

   /* Fill up send_disp*/
   send_disp[0]=0;
   disp=send_count[0];
   for(i=1;i<p;i++)
   {
     send_disp[i]=disp;
     disp+= send_count[i];


   }

   /* Fill up recv_disp and recv_count using A2A instead of Allgather*/
   /*MPI_Allgather(send_count,p,MPI_INT,recvG,p,MPI_INT,myComm);*/

   MPI_Alltoall(send_count,1,MPI_INT,
   			 recv_count,1,MPI_INT,
			 myComm);

   iGoingToRecv=0;
   for(i=0;i<p;i++) {
     iGoingToRecv+=recv_count[i];
   }
   /* Fill up recv_disp */
   disp=0;
   for(i=0;i<p;i++) {
     recv_disp[i]=disp;
     disp+= recv_count[i];
   }
   printf("Rank=%d: bucketing: B4 A2Av iGoingToRecv=%u\n", 
   		rank,iGoingToRecv);
   fflush(stdout);

   if(iGoingToRecv) {
        needs = sizeof(struct fkp)*iGoingToRecv;
        recv_buf= (struct fkp *) malloc(needs);
   }
   else  {
        needs = sizeof(struct fkp)*p; 
	   /* To make sure Alltoallv does not have a null recv_buf */;
        recv_buf= (struct fkp *) malloc(needs);
   }
   if(recv_buf) { thisPhaseSpace+=needs; spaceSoFar+=needs; }
   else AlertErr("recv_buf");
   fkeysFree=needs;


   buildFKPMPI(&n1,&sFKPtype);

   /* Do all-to-all communication*/
    MPI_Alltoallv(send_buf,send_count,send_disp,sFKPtype,recv_buf,
			recv_count,recv_disp,sFKPtype,myComm); 


   if(send_buf!=NULL) /* Ananth Change: for AiX/Sri 1/29/2004 */
     free(send_buf); 
  
   needs = sizeof(struct bucketkeys)*totBucks;
   bkeys2= (struct bucketkeys *)malloc(needs);
   if(bkeys2) {
	   /* thisPhaseSpace+=needs; spaceSoFar+=needs; */ 
	   /* Do not count even thisPhaseSpace because buckeys got released*/
           #ifdef debug
            /*printf("Rank:=%d Inside bucketing: bkeys2=%d\n",rank,needs); */
           #endif 
   }
   else AlertErr("bkeys2");

   sumbkeysnext=0;
   for(i=0;i<totBucks;i++)
   {
        if(i>=start && i<=end)
        {
            needs = sizeof(struct fkp)*rbuck[i];
	    if(rbuck[i]<=0) {
            	bkeys2[i].bk=NULL;
            	bkeys2[i].it=0;
	    }
            else {
            	bkeys2[i].bk= (struct fkp *)malloc(needs);
             	if(!bkeys2[i].bk)  AlertErr("bkeys2[i].bk");
            	bkeys2[i].it=0;
	    }
        }
        else
        {
            bkeys2[i].bk=NULL;
            bkeys2[i].it=0;
        }
   }
 
   for(i=0;i<iGoingToRecv;i++)
   {
      j = recv_buf[i].key2;
      bky = &(bkeys2[j]);
      copyFKP(&(bky->bk[bky->it]),recv_buf[i]);
      (bky->it)++;
   }

  
   final_fragkeys=recv_buf;
   r=0;
   for(i=start;i<=end;i++)
   {
      len = rbuck[i];
      if(len<=0) continue;
      bky = &(bkeys2[i]);
      for(t=0;t<len;t++)
      {
         copyFKP(&(final_fragkeys[r++]),bky->bk[t]);
      }
   }
   g_iSortedKeys = r;
   /*printf("Rank=%d: bucketing(): sorted keys=%d\n",rank,g_iSortedKeys); 
   fflush(stdout);*/

   if(rbuck!=NULL) /* Ananth Change: for AiX/Sri 1/29/2004 */
     free(rbuck);
   if(sted!=NULL) /* Ananth Change: for AiX/Sri 1/29/2004 */
     free(sted);
   if(send_count!=NULL) /* Ananth Change: for AiX/Sri 1/29/2004 */
     free(send_count);
   if(recv_count!=NULL) /* Ananth Change: for AiX/Sri 1/29/2004 */
     free(recv_count);
   if(recv_disp!=NULL) /* Ananth Change: for AiX/Sri 1/29/2004 */
     free(recv_disp);
   if(recvG!=NULL) /* Ananth Change: for AiX/Sri 1/29/2004 */
     free(recvG);
   if(send_disp!=NULL) /* Ananth Change: for AiX/Sri 1/29/2004 */
     free(send_disp);

   for(i=0;i<totBucks;i++) {
     if(bkeys2[i].bk!=NULL) free(bkeys2[i].bk);
   }
   if(bkeys2!=NULL) /* Ananth Change: for AiX/Sri 1/29/2004 */
     free(bkeys2);
   /*free(fkey);*/
  
} /* end bucketing */

void dispLocalBucks(int *buck,int totBuc)
{
   int i=0;
   printf("Slave=%d: Local Bucket\n",rank);

   for(i=0;i<totBuc;i++)
   {
        if(buck[i]>0)
        printf("<%d,%d>\n",i,buck[i]);
   }
} /* end dispLocalBucks */



void displayNonContainArr(int *containArr,int N,struct estgi *allEstGis)
{
    int i;
    FILE *confp,*nconfp;
    char s1[200],s2[200];
    

      
    strcpy(s1,"");
    sprintf(s1,"%s/ContainedESTs.%d.PaCE",g_sOutputFolder,N);
    confp = fopen(s1,"w");

    strcpy(s2,"");
    sprintf(s2,"%s/NonContainedESTs.%d.PaCE",g_sOutputFolder,N);
    nconfp = fopen(s2,"w");
     

     for(i=0;i<N;i++)
     {
         if(containArr[i]>=0) 
	 { 
		  /*strcpy(h1,""); 
		  strcpy(h2,""); 
                  sscanf(allEstGis[i].gi,"%*[^|]%*[|]%[^|]%*s",h1);
                  sscanf(allEstGis[containArr[i]].gi,"%*[^|]%*[|]%[^|]%*s",h2); 
		  fprintf(confp,"%s    %s\n",h1,h2);*/
		  fprintf(confp,"{%s} IN  {%s}\n",allEstGis[i].gi,allEstGis[containArr[i]].gi);
	 }
	 else
		  fprintf(nconfp,"%s\n",allEstGis[i].gi);
     }
     fclose(confp);
     fclose(nconfp);
	  
} /* end displayNonContainArr */



FILE *fileOpenForWrite(char *fname) {
	FILE *fp=NULL;
	int er;

	fp= fopen(fname,"w");
	er=errno;
	#ifdef debug
	if(!fp) {
		printf("Error: File %s could not be open for writing \n",fname);
		printf("Error: errno:%d, error: %s\n",er,strerror(er));
		printErr(er);
		fflush(stdout);
		terminate();
	}
	#endif
	return fp;
} /* end fileOpenForWrite */

/* Ananth: 6/25/2004 
 * Added to deterministically set the start and end fragment
 * read numbers without ambiguity 
 * Each proc sets its starting byte offset to the byte offset in
 * the file pointer value at which its corresponding new fasta
 * sequence starts.
 * It then communicates this value to the proc at its left
 * and sets the value it gets from its right to its end offset
 * Return Value:
 *	Begin and End File Offset entries for my rank
 * */

struct FileOffset *determineByteBoundaries(
					char *estfile,
					struct FileSize **global_fs,
					int *nFilesToRead) {

	FILE *fp;
	char frag[MAXFRAGSIZE];
	long aftOffset=0,b4Offset=0;
	long lastOffset=0;
	int from,to;
	MPI_Status sendstatus,recvstatus;
	MPI_Request sendreq;
	struct FileSize *getInputFileStats(char *,int *);
	void InitializeFileOffset(struct FileOffset *,int );

	// first take all input files and determine iOffset and partSize

	*global_fs = getInputFileStats(estfile,nFilesToRead);

	struct FileOffset *local_fo=NULL;
	struct FileOffset *fo;
	
	local_fo = (struct FileOffset *) malloc(sizeof(struct FileOffset)*2);
	InitializeFileOffset(local_fo,2);

	// local_fo[0] is begin and local_fo[1] is end

	long rem,curr_pos;
	int which_proc;
	int fo_it,nfile;
	struct FileSize *fs=NULL;

	which_proc=0;
	nfile=0;
	rem = g_partSize;
	curr_pos = 0;
	fo_it=0;

	// initial assignment of beginning offset
	while(1) {
		if(rank==0) {
			local_fo[0].fId = nfile; //==0
			local_fo[0].Off = 0L;
			break;
		}
		if(rank<which_proc) break;
		fs = &((*global_fs)[nfile]);
		if((curr_pos+rem) <= fs->size ) {
			curr_pos += rem;
			which_proc++;
			if(rank==which_proc) {
				//printf("Rank=%d: rank==which_proc: assigning %d, %ld \n",
				//		rank,nfile,curr_pos);
				local_fo[0].fId = nfile;
				local_fo[0].Off = curr_pos;
				break;
			}
			if(curr_pos==fs->size) {
				nfile++;
				curr_pos=0;
			}
			rem = g_partSize;
			continue;
		}
		else {
			//printf("Rank=%d: Moving to next file %d: %ld, %ld, %ld \n",
			//		rank,nfile+1,rem,fs->size,curr_pos);
			rem = rem - (fs->size - curr_pos);
			nfile++;
			curr_pos=0;
			continue;
		}
	} /* end while*/

	//printf("Rank=%d: determineByteBoundary: Initial Start %ld of file (%s) : (file size %ld)\n",
		//		rank,local_fo[0].Off,(*global_fs)[local_fo[0].fId].name,
		//		(*global_fs)[local_fo[0].fId].size);
	//fflush(stdout);

	// adjust the beginOff 

	char *file;
	int bBeginOffSet;
	int fId;
	long beginOff;

	bBeginOffSet=0;
	fId=local_fo[0].fId;
	beginOff = local_fo[0].Off;
	

	while(1) {
		// loop variants: fid, beginOff
		if( fId >= *nFilesToRead ) {
			printf("Warning: Rank=%d: to be assigned no part of input !\n",rank);
			fflush(stdout);
			fId = -1;
			beginOff=0;
			break;
		}

		fs = &((*global_fs)[fId]);
		file = fs->name;
		fp = fopen(file,"r");
		if(fseek(fp,beginOff,SEEK_SET)!=0) {
			printf("Rank=%d: Error: determineByteBoundary fseek failed (%ld, %ld)!\n",
					rank, beginOff,fs->size);
			fflush(stdout);
			terminate();
		}

  		while(!feof(fp)) {
			strcpy(frag,"");
			b4Offset= ftell(fp);
			lastOffset=b4Offset;
			fscanf(fp,"%[^\n]%*\n",frag);
			aftOffset= ftell(fp);
			if(b4Offset==aftOffset) { 
				b4Offset= ftell(fp);
				fseek(fp,b4Offset+1L,SEEK_SET);
				continue;
			}
			if(strcmp(frag,"")==0) continue;
			if(strstr(frag,tag)==0)  /* continue until tag > is seen */
				continue;

			/* This fix was crucial in letting read byte boundaries from
		 	* the input fasta file properly */
			if(frag[0]!='>') continue;

			bBeginOffSet=1;
			// adjust begin file offset to start of the fasta header
			beginOff = b4Offset; 

			break;
		} // inner while
		fclose(fp);
		if(bBeginOffSet) break;

		fId++;
		beginOff=0;
		
	} // outer while

	local_fo[0].fId = fId;
	local_fo[0].Off = beginOff;

	int endFileId;
	long endOff=0;

	// communicate beginOff of rank and assign it to endOff of rank-1
	from=rank+1;
	to=rank-1;

	if(rank>0) MPI_Isend(&fId,1,MPI_INT,to,0,myComm,&sendreq);
	if(rank<(p-1)) MPI_Recv(&endFileId,1,MPI_INT,from,0,myComm,&recvstatus);
	if(rank>0) MPI_Wait(&sendreq,&sendstatus);

	if(rank>0) MPI_Isend(&beginOff,1,MPI_LONG,to,0,myComm,&sendreq);
	if(rank<(p-1)) MPI_Recv(&endOff,1,MPI_LONG,from,0,myComm,&recvstatus);
	if(rank>0) MPI_Wait(&sendreq,&sendstatus);

	if(rank==(p-1)) {
		if(local_fo[0].fId>=0) {
			endFileId = *nFilesToRead-1;
			endOff = (*global_fs)[*nFilesToRead-1].size + 1L;
		}
		else {
			endFileId = -1;
			endOff = 0;
		}
	}

	local_fo[1].fId = endFileId;
	local_fo[1].Off= endOff;

	printf("Rank=%d: determineByteBoundaries: File: [%d:%ld] , [%d:%ld]\n",
				rank,local_fo[0].fId,local_fo[0].Off,
				local_fo[1].fId,local_fo[1].Off);
	fflush(stdout);

	return local_fo;
} /* end determineByteBoundaries*/

struct FileSize *getInputFileStats(char *ftemplate,int *nFiles) {
	int i;
	char f[MAXFILENAMESIZE];
	FILE *fp;
	struct FileSize *fs=NULL;
	int f1;

	*nFiles=0;

	i=0;
	while(1) {
		sprintf(f,ftemplate,i++);
		if(strcmp(f,ftemplate)==0) {
			// there is only one input file
			i++;
			break;
		}
		fp = fopen(f,"r");
		if(!fp) break;
		fclose(fp);
	}
	*nFiles=i-1;

	fs = (struct FileSize *)malloc(sizeof(struct FileSize)*(*nFiles));
	for(i=0;i<*nFiles;i++) {
		sprintf(f,ftemplate,i);
		f1 = open(f,O_RDONLY);
		strcpy(fs[i].name,f);
   		fs[i].size = lseek(f1,0L,SEEK_END);
		close(f1);
	}

	g_partSize=0;
	for(i=0;i<(*nFiles);i++) g_partSize += fabs(fs[i].size/p);


	if(rank==0) {
		printf("Rank=%d: getInputFileStats: Number of input files = %d (partSize = %ld)\n",rank,*nFiles,g_partSize);
		fflush(stdout);
	}

	return fs;

} /* end getInputFileStats */

void InitializeFileOffset(struct FileOffset *a,int n) {
	int i;
	for(i=0;i<n;i++) {
		a[i].fId = -1;
		a[i].Off = 0L;
	}
	
} /* InitializeFileSizeOffset */



int min(int a,int b) {
	if(a<=b) return a;
	return b;
} /* end min */


void debugSync(MPI_Comm mycomm) {
	printf("Rank=%d: DebugSync...\n",rank); fflush(stdout);
	MPI_Barrier(mycomm);
	MPI_Abort(MPI_COMM_WORLD,1);
}


void printIntArr(char *name,int *arr,int size) {
	int i=0;
	printf("Rank=%d: %s<%d>  ",rank,name,size);
	for(i=0;i<size;i++) printf("[%d], ",arr[i]);
	printf("\n");
	fflush(stdout);

} /* end printIntArr*/


void checkInRange(int a,int lb,int ub,char *tag) {
	if(a>ub || a<lb) {
			printf("Rank=%d: Error: Id %d out of range <%d,%d> for %s\n",rank,a,lb,ub,tag);
			fflush(stdout);
			terminate();
	}
} /* end checkInRange */


int determineA2AIterations(unsigned int g_maxBuckSize) {

	unsigned int upto=0,nextI=0;
	int i=0,j=0;
	int n_iter=0;
	int n_A2A_iter=0;
	const unsigned int MaxScratchSpace=256000000; /* 256MB */
        int *markString=NULL;
	unsigned int uit=0,valid_upto=0;
	int singler=0;
	char *markStarters=NULL;
	unsigned int n_spaceForScratch=0;
	int fid=0,len=0;
	unsigned int recentStart=0,currKey=0;
	
	if(g_maxBuckStrSize>=MaxScratchSpace/2) {
		printf("Rank=%d: Warning: determineA2A: g_maxBuckStrSize %u greater than half of ScratchSpace %u\n",rank,g_maxBuckStrSize,MaxScratchSpace/2);
		fflush(stdout);
	}


	markString=(int *)malloc(sizeof(int)*N);
	if(!markString) AlertErr("determineA2A_markString");
	
	for(i=0;i<N;i++) markString[i]=-1;

	markStarters=(char *)malloc(sizeof(char)*g_iSortedKeys);
	if(!markStarters) AlertErr("determineA2A_markStarters");

	for(uit=0;uit<g_iSortedKeys;uit++) markStarters[uit]='0';

	nextI=0;
	while(nextI<g_iSortedKeys) {
		markStarters[nextI]='1';
		singler=0;
		n_iter++;
		upto=nextI+g_maxBuckSize;
	     if(upto>=g_iSortedKeys) {
			upto=g_iSortedKeys-1; 
			break;
		}
		while(final_fragkeys[upto].key2 == final_fragkeys[upto-1].key2 
				&& upto>nextI) {
			upto--;
		}
		if(upto!=nextI) upto--;
		else {
			upto=nextI+1;
			while((upto+1)<g_iSortedKeys && 
				final_fragkeys[upto].key2==final_fragkeys[upto+1].key2) {
				upto++;
			}
			singler=1;
		}

		/* Now calculate based on the string lengths 
		   in the current span
		*/
		if(!singler) {
			n_spaceForScratch=0;
			recentStart=nextI; 
			uit=valid_upto=nextI;
			currKey=final_fragkeys[uit].key2;
			for(uit=nextI;uit<=upto;uit++) {
				if(final_fragkeys[uit].key2!=currKey) {
					currKey=final_fragkeys[uit].key2;
					recentStart=uit;
				}
				fid=final_fragkeys[uit].fragid;
				len=AllESTlen[fid]+1;
				if(markString[fid]==n_iter) {
				     valid_upto=uit;
					continue;
 				}
				markString[fid]=n_iter;
				if( (n_spaceForScratch+len)>=MaxScratchSpace) {
					valid_upto=recentStart-1;
					printf("Rank=%d: determineA2A: readjusting space.. <%u, %u, %u>\n",
							rank,nextI,valid_upto,n_spaceForScratch);
					break;
				}
				valid_upto=uit;
				n_spaceForScratch+=len;
			
			} /* end for uit */	
			upto=valid_upto;
			if(upto<nextI) {
				printf("Rank=%d: Warning: determineA2A: upto %u less than nextI %u \n",
							rank,upto,nextI);
				fflush(stdout);
			}
		}

	 	/* Now goto next hop*/
		nextI=upto+1;
	} /* end while nextI */

        /*printf("Rank=%d: Projected number of iterations=%d\n",rank,n_iter);*/
   	MPI_Allreduce(&n_iter,&n_A2A_iter,1,MPI_INT,MPI_MAX,myComm);
        if(rank==0) printf("Rank=%d: MPI_Reduced number of iterations=%d\n",rank,n_A2A_iter);

        needs=sizeof(int)*n_A2A_iter;
	g_itStarters=(int *)malloc(sizeof(int)*n_A2A_iter);
	if(!g_itStarters) AlertErr("g_itStarters");
	else thisPhaseSpace+=needs;
	for(i=0;i<n_A2A_iter;i++) g_itStarters[i]=-1;

	i=0;
	for(uit=0;uit<g_iSortedKeys;uit++) {
		if(markStarters[uit]!='1')  continue;
		if(i>=n_A2A_iter) {
			AlertErr("determineA2A_startoutofbound");
		}
		g_itStarters[i++]= (int) uit;
		if(rank==-1) {
			printf("Rank=%d: itStarter[%d]=%d\n",rank,i-1,g_itStarters[i-1]);
			fflush(stdout);
		}
		
		if(i>1 && final_fragkeys[g_itStarters[i-1]].key2==final_fragkeys[g_itStarters[i-2]].key2)  {
			
			printf("Rank=%d: Warning: itStart %u ",rank,g_itStarters[i-1]);
			fflush(stdout);
		}
	}

	free(markString);
	free(markStarters);
	/*printIntArr("g_itStarters",g_itStarters,n_A2A_iter);*/
	return n_A2A_iter;
} /* end determineA2AIteration */


  /* A2A change:
   * Add this function to do 2 A2A 
   */

void doA2A(int *send_buf,int start,int end,int n_A2A_rounds) {

  	char *stringMarker=NULL;
  	int *recv_buf=NULL;
  	int *send_disp=NULL,*recv_disp=NULL;
  	int *send_count=NULL,*recv_count=NULL;
  	int sendBufCap=0,sendBufSize=0,recvBufSize=0,sendStrSize=0;
     unsigned int iGoingToRecv=0;
  	int upto=0;
  	char *tempStr=NULL;
  	int s=0,t=0;
  	int currEnd=0,countForThisProc=0,offset=0,fid=0;
  	int disp=0;
  	int mysted[2],b_inRangeFirst=1,b_inRange=0;
	int i=0,j=0;
	unsigned int it=0;


	int left,right;
	int *strlen_ToExpect=NULL;
	int nxt=0;
	void doSendRecv(int ,int ,
                int *,int *,int *,
                int *,int *,int *,
                char *,int *,int *);
  
	if(rank==0) {
		printf("Rank=%d: A2A Round %d : A2A getting ready <%d,%d>...\n",rank,n_A2A_rounds,start,end);
		fflush(stdout);
	}
	MPI_Barrier(myComm);

  
 	g_iGSTscratchSpace=N*sizeof(int); /* upperbound on send_buf size*/
 	g_iGSTscratchSpace+=4*p*sizeof(int); /* for the count and disp arrays */

  	send_count = (int *)malloc(sizeof(int)*p);
  	recv_count = (int *)malloc(sizeof(int)*p);
  	send_disp = (int *)malloc(sizeof(int)*p);
  	recv_disp = (int *)malloc(sizeof(int)*p);

  	if(!send_count || !send_disp || !recv_count || !recv_disp) 
	  AlertErr("A2A_send_recv_count_disp");

  	strlen_ToExpect = (int *)malloc(sizeof(int)*p);
  	if(!strlen_ToExpect) AlertErr("strlen_ToExpect");
  	for(j=0;j<p;j++) strlen_ToExpect[j]=0;

  	/* Second A2A */
	if(g_recv_str) free(g_recv_str);
	g_recv_str=NULL;

  	stringMarker = (char *)malloc(sizeof(char)*N);
	g_iGSTscratchSpace+=N;

  	mysted[0]=g_dataPartition[2*rank];
  	mysted[1]=g_dataPartition[2*rank+1];

	for(j=0;j<N;j++) {
		if(j<mysted[0] || j>mysted[1]) AllESTs[j].s=NULL;
		stringMarker[j]='0';
	}
	for(j=start;j<=end;j++) {
		stringMarker[final_fragkeys[j].fragid]='1';
	}
		  
  	for(j=0;j<p;j++) {
	  	send_count[j]=recv_count[j]=0;
	  	send_disp[j]=recv_disp[j]=0;
  	}
		
	/* Prepare of what string to ask whom */
	/* First A2A */
	s=0;
	t=0;
	for(j=0;j<N;j++) {
		if(stringMarker[j]!='1') continue;
		if((j>=mysted[0] && j<=mysted[1])) b_inRange=1;
		else b_inRange=0;

		if(!b_inRange || (b_inRange && b_inRangeFirst)) {
			checkInRange(t,0,p-1,"send_buf_init_t");
			currEnd=g_dataPartition[2*t+1];
		}

		if(j>currEnd) {
			t++;
			while(t<p) {
				currEnd=g_dataPartition[2*t+1];
				if(j<=currEnd) break;
				send_count[t++]=0;
			}
			if(t>=p) {
				printf("Rank=%d: Error: j=%d, currEnd=%d, range exceeded!\n",
						rank,j,currEnd);
				fflush(stdout);
				terminate();
			}
		}
			
		if((j>=mysted[0] && j<=mysted[1])) {
			if(b_inRangeFirst) {
				send_count[t]=0;
				t++;
			}
			b_inRangeFirst=0;
			continue;
		}

		/* j<=currEnd */

		send_buf[s++]=j;
		send_count[t]++;
		strlen_ToExpect[t]+=AllESTlen[j]+1;
	} 
	sendBufSize=s;
	b_inRangeFirst=1;

	send_disp[0]=0;
	disp=send_count[0];
	for(t=1;t<p;t++) {
		send_disp[t]=disp;
		disp+= send_count[t];
	}

  	/* First A2A */
  	recv_buf=NULL;

        /*MPI_Allgather(send_count,p,MPI_INT,recvG,p,MPI_INT,myComm);*/
	   MPI_Alltoall(send_count,1,MPI_INT,
	   			 recv_count,1,MPI_INT,
				 myComm);

        iGoingToRecv=0;
        for(t=0;t<p;t++) {
                iGoingToRecv+=recv_count[t];
        }

	   disp=0;
        for(i=0;i<p;i++) {
                recv_disp[i]=disp;
                disp+= recv_count[i];
        }

        needs=sizeof(int)*iGoingToRecv;
        if(iGoingToRecv) recv_buf=(int *)malloc(sizeof(int)*iGoingToRecv);
        else {
                needs=sizeof(int)*p;
                recv_buf=(int *)malloc(sizeof(int)*p); 
			 /* Because A2Av needs a non-null recv_buf */
        }
        if(!recv_buf) AlertErr("A2Av_1_recv_buf");
        recvBufSize=iGoingToRecv;
        if(recv_buf) { g_iGSTscratchSpace+=needs; }

        if(recv_buf==NULL) {
                printf("Rank=%d: Warning: A2Av_1: Recv buffer is empty!!\n",rank);
                fflush(stdout);
        }

     MPI_Alltoallv(send_buf,send_count,send_disp,MPI_INT,
				recv_buf,recv_count,recv_disp,MPI_INT,
				myComm);


        /* Ask */
	/* Second A2A */

        /* Prepare the Send Buffer in the following way:
	*	- do p-Send/Recvs as a p-permutation
	*	- for i=1 to p-1
	*		Send/Recv data between (rank) and ((rank+i)%p) processors
	*/ 

	/* First allocate and initialize g_recv_str */
	iGoingToRecv=0;
  	for(j=0;j<p;j++) iGoingToRecv+=strlen_ToExpect[j];

   	printf("Rank=%d: doA2Av_2: iGoingToRecv=%u bytes\n", rank,iGoingToRecv);

        if(iGoingToRecv>0) {
                needs = sizeof(char)*iGoingToRecv;
                g_recv_str = (char *) malloc(needs);
        }
        else  {
                printf("Rank=%d: Warning: A2Av_2: Recv buffer is empty!!\n",rank);
                fflush(stdout);
                /* To make sure Send/Recv does not have a null recv_buf */;
                needs = sizeof(char)*p;
                g_recv_str = (char *) malloc(needs);
        }
        if(g_recv_str) { g_iGSTscratchSpace+=needs; }
        else AlertErr("A2Av_2_g_recv_str");

        /* Initialize the g_recv_str with null entries (character=\0)*/
        for(it=0;it<iGoingToRecv;it++) g_recv_str[it]='\0';

	/* pointers into the A2Av_1 arrays 
	*   Beware: Misuse of recv_count and recv_disp for nonA2Av purposes!!
	*/
	t=0;
	for(i=0;i<p;i++) {
		if(recv_count[i]>0) recv_disp[i]=t; 
		else recv_disp[i]=-1;
		t+=recv_count[i]; 
		
	}
	/* Now deal with the send_str as p-permutation send/recv */
        nxt=0;
	for(i=1;i<=p;i++) {
		right=(rank+i)%p;
		left=(rank-i+p)%p;
		doSendRecv(left,right,recv_count,recv_disp,recv_buf,
				send_count,send_disp,send_buf,
				g_recv_str,&nxt,strlen_ToExpect);
		
        }
	

		 
	if(recv_buf) free(recv_buf); 
	if(send_count) free(send_count);
	if(recv_count) free(recv_count);
	if(send_disp) free(send_disp);
	if(recv_disp) free(recv_disp);
	if(stringMarker) free(stringMarker);
	if(strlen_ToExpect) free(strlen_ToExpect);

	printf("Rank=%d: A2A Round %d : A2A ready for GST construction (scratch space %u bytes)\n",rank,n_A2A_rounds,g_iGSTscratchSpace+g_iSendRecvScratchSpace);
	fflush(stdout);

} /* end doA2A */



/* Change for PaCEA2A_v5 */

/* readTreeStrings: Reads a subset of strings from the Fasta file 
*	that are represented in the local GST, and 
*	with the constraint that the sum of the local strings
*	cannot exceed g_scratchMemory.
*/

void readTreeStrings(char *fastafile,int N) {
   int i=0;
   struct loadRes ldRes;
   int iNewLocalFrags=0;

   unsigned int scratchMemory=0;
   char *testStr=NULL;
   char *thisLset=NULL;
   int upto=-1,maxupto;
   int temp=0;
   int b_exceeded=0;
   unsigned int max_allowed_memory=0;


   thisPhaseSpace=spaceSoFar;
   scratchMemory = getScratchMemory();

   max_allowed_memory = scratchMemory + myFragsLen ;
   printf("Rank=%d: readTreeStrings: Tree memory Capacity = %u bytes (%u MB)\n",
   				rank,max_allowed_memory,max_allowed_memory/1000000);



   /* Already g_TreeStringMarker has all the tree strings 
   * 	So sum up their length and re-mark from first id only those 
   * 	which add up to a max of scratchMemory.
   */
 
   upto=-1;
   b_exceeded=0;
   g_TreeStringsSize=0;
   for(i=0;i<N;i++) {
	if(g_TreeStringMarker[i]=='1' && !b_exceeded) {
		temp=AllESTlen[i]+1;
		if((g_TreeStringsSize+temp)>max_allowed_memory) {
			g_TreeStringMarker[i]='0';
			b_exceeded=1;
		}
		g_TreeStringsSize+=temp;
		upto=i;
	} else {
		g_TreeStringMarker[i]='0';
	}
   }
   maxupto=upto;
   if(!b_exceeded) {
	/* There is still more space to fill !!
	*	So load as many strings as possible, 
	*	even those which are not in GST
	*	until max_allowed_memory is reached
	*/
        for(i=0;i<N;i++) {
		if(g_TreeStringMarker[i]=='1') continue;
		temp=AllESTlen[i]+1;
		if((g_TreeStringsSize+temp)>max_allowed_memory) break;
		g_TreeStringsSize+=temp;
		g_TreeStringMarker[i]='1'; /* toggle it*/
		upto=i;
	}
   }
   if(maxupto<upto) maxupto=upto;
   
   /* Scan g_TreeStringMarker until "upto" and load those strings
   * 	into memory from IO
   *    Caution: This scheme will have to be modified for >= 4 GB 
   *    to be loaded in the memory.
   */

   needs=N*sizeof(unsigned int);
   g_TreeStringPtr =(unsigned int *)malloc(needs);
   if(!g_TreeStringPtr) AlertErr("readTreeStrings_g_TreeStringPtr");
   else { 
	thisPhaseSpace+=needs;
	spaceSoFar+=needs;
   }
   /* Initialize the string index pointers to out of bounds on
    *	the g_TreeStrings array.
   */
   for(i=0;i<N;i++) g_TreeStringPtr[i]=g_TreeStringsSize+1;

   needs=g_TreeStringsSize*sizeof(char);
   g_TreeStrings =(char *)malloc(needs);
   if(!g_TreeStrings) AlertErr("readTreeStrings_g_TreeStrings");
   else { 
	thisPhaseSpace+=needs;
	spaceSoFar+=needs;
   }

   iNewLocalFrags = loadTreeStrings(fastafile,global_fs,nFilesToRead,
   							g_TreeStringPtr,g_TreeStrings,
							g_TreeStringMarker,maxupto);

   if(iNewLocalFrags<=0) {
     printf("Rank=%d: Some Serious Warning in readTreeStrings: No Strings loaded\n",rank);
     fflush(stdout);
   }
  
   if(g_TreeStringMarker) free(g_TreeStringMarker);

} /* end readTreeStrings */
/* End of Change for PaCEA2A_v5 */


/* Change for PaCEA2A_v6 */
/* doSendRecv:
*	Sends string data rank->to and 
* 	receives string data from->rank,
*	all as per request in the A2Av_1 exchange.
*	Also set AllESTs[].s in the process.
*	PS: send_count,_disp,_buf and recv_count,_disp,_buf buffers 
*		are all as per A2Av_1 exchange.
*/

void doSendRecv(int from,int to,	
		int *recv_count,int *recv_disp,int *recv_buf,
		int *send_count,int *send_disp,int *send_buf,
		char *recv_str,int *nxt,int *strlen_ToExpect){


  char *send_str=NULL;
  int i,j,t;
  int sendStrSize=0;
  int startoff,endoff;
  int n_lenRecvd=0;
  int index=0,fid;
  char *tempStr=NULL;
  int s,offset;
  MPI_Status recv_status;

  if(recv_disp[to]<0) { 
	/* right proc has not requested any string from myrank */
	send_str=NULL;
	sendStrSize=0;
  }
  else {
	/* populate send_buf as per request from right proc */
	startoff=g_dataPartition[2*rank];
	endoff=g_dataPartition[2*rank+1];
	index=recv_disp[to];
        sendStrSize=0;
        for(t=0;t<recv_count[to];t++) {
                checkInRange(recv_buf[index],startoff,
                        endoff,"recv_buf[index]");
                sendStrSize+=AllESTlen[recv_buf[index]]+1;
		index++;
        }
        needs=sendStrSize;
        if(sendStrSize>0) {
                send_str=(char *)malloc(sendStrSize);
                if(!send_str) AlertErr("A2Av_2_send_str");
        }
	else {
		printf("Rank=%d: Warning: doSendRecv: sendStrSize is %d inspite of checking before!!\n",
				rank,sendStrSize);
		fflush(stdout);
		terminate();
	}
	s=0;
	index=recv_disp[to];
        for(t=0;t<recv_count[to];t++) {
		/*if(rank==0 && to==1) {
			printf("Rank=%d: to %d: id=%d, of length %d (%d) \n",
				rank,to,recv_buf[index],strlen(AllESTs[recv_buf[index]].s),
							AllESTlen[recv_buf[index]]);
			fflush(stdout);
		}*/
		tempStr=AllESTs[recv_buf[index]].s;
		index++;
                strcpy(&(send_str[s]),tempStr);
                s+=strlen(&(send_str[s]))+1;
        }
	if(s!=sendStrSize) {
		printf("Rank=%d: Error: doSendRecv: sendStrSize %d unmatching %d send_str copied!!\n",
				rank,sendStrSize,s);
		fflush(stdout);
	}
  }

  /*printf("Rank=%d: to %d: sendStrSize=%d and strlen_ToExpect[%d]=%d\n",
		rank,to,sendStrSize,from,strlen_ToExpect[from]); 
  fflush(stdout);*/

  /* &recv_str[*nxt] is the recv_buf for MPI_Sendrecv
  *	however, recvStrSize=strlen_ToExpect[left]
  */

  MPI_Sendrecv(send_str,sendStrSize,MPI_CHAR,to,0,
		       &(recv_str[*nxt]),strlen_ToExpect[from],MPI_CHAR,from,0,
			myComm,&recv_status);


  /* fix the AllESTlen pointers */
  n_lenRecvd=0;
  j=*nxt;
  for(t=send_disp[from];t<(send_disp[from]+send_count[from]);t++) {
                fid=send_buf[t];
                checkInRange(fid,0,N-1,"send_buf[t]");
                tempStr=&(recv_str[j]);
                if(AllESTlen[fid]!=strlen(tempStr)) {
                        printf("Rank=%d: Error: A2Av_2_sendrecv string %d should have length %d (not %d %s) \n",
				  rank,fid,AllESTlen[fid],strlen(tempStr),tempStr);
                        fflush(stdout);
                        terminate();
                }
                AllESTs[fid].s=&(recv_str[j]);
                n_lenRecvd+=strlen(tempStr)+1;
                j+=strlen(tempStr)+1;
  }

  if(n_lenRecvd!=strlen_ToExpect[from]) {
                        printf("Rank=%d: Error: A2Av_2_sendrecv received str %d mismatch w.r.t strlen_ToExpect %d\n",
				  rank,n_lenRecvd,strlen_ToExpect[from]);
                        fflush(stdout);
			terminate();
  }
	
  *nxt+=strlen_ToExpect[from];

  if(send_str) free(send_str);
  if(g_iSendRecvScratchSpace>sendStrSize) g_iSendRecvScratchSpace=sendStrSize;
  
} /* end doSendRecv*/
/* End of Change for PaCEA2A_v6 */

FILE **OpenInputFilePointers(struct FileSize *gfs,int nFiles) {
	int i;
	FILE **fps;
	fps = (FILE **)malloc(sizeof(FILE *)*nFilesToRead);
	for(i=0;i<nFiles;i++) {
		struct FileSize *fs_ptr;
		fs_ptr = &(gfs[i]);
		fps[i] = fopen(fs_ptr->name,"r");
		if(!fps[i]) {
			printf("Rank=%d: Error OpenInputFilePointers: File %s is possibly empty or corrupted!\n",
				rank,fs_ptr->name);
			fflush(stdout);
			terminate();
		}
	}
	return fps;
}

void CloseInputFilePointers(FILE **fps,int nFiles) {
	int i;
	for(i=0;i<nFiles;i++) if(fps[i]) fclose(fps[i]);
}


int *OpenInputFileDescriptors(struct FileSize *gfs,int nFiles) {
	int i;
	int *fds;
	fds = (int *)malloc(sizeof(int)*nFilesToRead);
	for(i=0;i<nFiles;i++) {
		struct FileSize *fs_ptr;
		fs_ptr = &(gfs[i]);
		fds[i] = open(fs_ptr->name,O_RDONLY);
		if(fds[i]<0) {
			printf("Rank=%d: Error OpenInputFileDescriptors: File %s is possibly empty or corrupted!\n",
				rank,fs_ptr->name);
			fflush(stdout);
			terminate();
		}
	}
	return fds;
}

void CloseInputFileDescriptors(int *fds,int nFiles) {
	int i;
	for(i=0;i<nFiles;i++) if(fds[i]>=0) close(fds[i]);
}

/* getSeqtoMem : reads a sequence from disk  */
/* Ananth 2007 - Jaric' suggestion: replace fseek+fscanf by pread */

void getSeqtoMem(int fd,long byteOff,int bytes,char *frag) {
  struct timeval t_1,t_2;

   if(byteOff<0) 
   {
      strcpy(frag,"");
      return;
   }

   n_callsToSeqtoMem++;

   gettimeofday(&t_1,0);
   if( pread(fd,frag,bytes,(off_t) byteOff) < bytes) {
     printf("Rank=%d: Error : getSeqtoMem failed to IO read sequence from %d offset! \n",rank,byteOff);
     fflush(stdout);
     terminate();
   }

   frag[bytes]='\0';

   gettimeofday(&t_2,0);
   t_pread+= (t_2.tv_sec-t_1.tv_sec)*1000000 + (t_2.tv_usec-t_1.tv_usec);
   if(t_pread>=MAXULONGBY2) {
		t_spread+=t_pread/1000000;
		t_pread=0;
   }	

   n_bytesReadFromIO+=strlen(frag)+1;
} /* end getSeqtoMem */
