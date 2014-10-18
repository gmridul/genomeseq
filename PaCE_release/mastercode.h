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



#ifndef __MASTERCODE_H__
#define __MASTERCODE_H__



#include <mpi.h>
#include "uFind.h"
#include "dynamic.h"
#include "est.h"

#define MATCHLINESIZE 100
#define FNAME 200
#define MAXNEWPAIRS 200000
#define MAXLOADPERPROC 50
/* max pairs that could be stored on a slave after generation before dispatch to master */
#define MAXBUFSIZE 200100 
#define MAXPAIRS 200000
/* Will try to keep the master work buffer as full as possible */
/*#define KEEP_MBUF_FULL*/ 
/* Parameterize using Keep_Mbuf_Full in PaCE.cfg */
/*  eg.,
	Keep_Mbuf_Full 0	=> dont keep it full
	Keep_Mbuf_Full 10	=> keep it 1/10th full
	Keep_Mbuf_Full 5	=> dont keep it 1/5th full
*/

#define Use_Ssend

/*#define BLUEGENEL*/
#define BGL_debug
#define PrintPairBatchSize 100  /* in number of lines */
#define PrintPairOneLine 1000   /* each line is of max size this*/

/*  This is a framework for a very simple application.  */
/*  Input : An array of integers of size say N                */
/*  Output : An array of integer pairs 'pairs' such that      */
/*           pairs[i] has the ith ordered pair from the list of */
/*           N choose 2 and has (x,y) such that x<y if i is even*/
/*           and (y,x) if i is odd.  Hence output size is N Choose 2 */      
/*  Unit Operation for the slaves :   Take m pairs at a time,   */
/*           If index of a pair is odd then swap the elements' positions */


/* Application Specific   */
/* This is the data type of the output array */
struct pair{
 char dummy ;
 int f1;
 int f2;
 int l1;
 int l2;
 int p1;
 int p2;
 int d;
 float e;
 char type;
};


struct resPair{
 /*char dummy;*/
 int maxBorderScore;
 int idealBorderScore;
 int maxScore;
 int idealScore;
 int maxScoreCoverage;
 int f1;
 int f2;
 char type;
 char atag;    /* if type='C' : atag is contained in the other string  */
               /* if type='S' : atag is left of the other string       */
};

/* master receives this structure from slave everytime */
/* this is not used now : Obsolete: Planned*/
/*struct work {
  struct pair newPairs[MAXNEWPAIRS];
  struct resPair res[MAXLOADPERPROC];
}; */ 


struct cBuffer {
  struct pair *pairs;
  int front;
  int rear;
  int currSize;
  int totSize;
};

/* This is not used now */
/*  struct slaveBuf {
  struct pair pairs[MAXNEWPAIRS];
  int start;
  int end;
  int currSize;
  int totSize;
}; */

struct waitq {
  int *q;
  int front;
  int rear;
  int currSize;
  int totSize;
};

struct mastVar {
  /* STATE VARIABLES */
  int P;
  int R;
  int E;
  int W;
  struct cBuffer MBUF;
  int ACTIVE;
  int TORECV;
  struct ufind *UFcluster;
  int *UFcluster_size;
  int nSlaves;
  int *active;
  int dst;
  struct waitq waitQ;
  struct pair *nxtBatch;
  struct pair newPairs[MAXNEWPAIRS];
  struct resPair *newResults;
  /* OTHER VARIABLES */
  unsigned int numTotalPairs;
  unsigned int numPairsAssigned;
  unsigned int numResultsRecv;
  /* Temp variables: accumulative pairs count */
  int tmp_numTotalPairs;
  int tmp_numPairsAssigned;
  int tmp_numResultsRecv;

  unsigned long tIdleTime;
  unsigned long tSendTime;
  unsigned long tsIdleTime;
  unsigned long tsSendTime;
  int numClusterDumpSteps;
  MPI_Request send_req;
};

struct slaveVar {
   /*  STATE VARIABLES */
   int P;
   int R;
   int oE;
   int nE;
   int W;
   struct cBuffer SBUF;
   struct dynParams param;
   struct pair newPairs[MAXNEWPAIRS];
   struct pair *currWork;
   struct resPair *resToSend;
   int terminationSent;
   int gotit;
   char **workPairs;
   int *stringsForWork;
   /*  OTHER VARIABLES */
   int mastRank;
   int numPairsGenerated;
   unsigned long tIdleTime;
   unsigned long tSendTime;
   unsigned long tsIdleTime;
   unsigned long tsSendTime;
   MPI_Request send_req;
   MPI_Request recv_req;
   MPI_Status recv_status;
   FILE *estFP;
};


extern int rank,p;
extern int N;
extern int numPairsProvided;
extern int numPairsServiced;
extern int numTotalPairs;
extern int LOADPERPROC;
extern int *containArr;
extern unsigned int iPairsAccepted;
extern unsigned long tOnlyFetch;
extern int g_matchscore;
extern int g_bTranscriptsTogether;
extern int g_bReportSplicedCandidates,g_bReportMaximalPairs,g_bReportAcceptedPairs;
extern int g_bReportGeneratedPairs;
extern int EndToEndScoreRatioThreshold,EndtoEndAlignLenThreshold,MaxScoreRatioThreshold,TranscriptCoverageThreshold;
extern FILE *fpSpliced;
extern FILE *fpMaximalPairs;
extern FILE *fpAcceptedPairs;
extern FILE *fpMaximalSubstrings;
extern FILE *fpGeneratedPairs;
extern char g_sOutputFolder[];
 
extern long g_partSize;
extern char g_printpairstr[];

extern int LargeClusterThreshold;
extern int g_bOutputLargeMerges;
extern int *UFcluster_size;
extern FILE *fpLargeMerges;
extern MPI_Comm myComm;

void processPair(struct dynResult *,struct pair rBuf,struct dynParams,char *,char *);
/*void master();
void slave(int);*/
void killAll();

void readAllESTs(char *,int,int);
void slideWindowK(int);
void buildTree(int);
void localSortNodes();
void generatePairs(int);
void readMyESTs(char *,int,int);
void readTreeStrings(char *fastafile,int N);


/*  This method incorporates the results just after they are obtained */
/*  from a slave and just before the slave is alloted the next batch of */
/*  work load.  This is Application data specific */
void interpretResults1(struct resPair *batch,int count,struct ufind *UF,int *);


/*  This method incorporates the results just after  */
/*  a slave is alloted the next batch of */
/*  work load.  This is Application data specific */

/* void interpretResults2(int data[],int size,struct pair *pairs); */



/*  This method fetches the next batch of data to be processed   */
/*  The logic of checking whether there is data available in the */
/*  initialized output array can be used in a generic way  although */
/*  the data types and way to populate the data in the array may vary */
/*  with applications                                               */
int prepareNextBatch(struct cBuffer *,struct pair *,struct ufind *);

/*  This is an Application Specific method to display the result  */
void dispPairs(struct pair *pairs, int size);
void dispPairsDetail(struct pair *p, int size);

void buildPairMPI(struct pair *t, MPI_Datatype *dt);
void buildresPairMPI(struct resPair *t, MPI_Datatype *dt);

int getLoadPerProc();
unsigned int getScratchMemory();

int dispAllClusters(struct ufind *uf,int size,struct estgi*,int *singletons);

int insertCBuf(struct pair *newPrs,int nPairs,struct cBuffer *cBuf,struct ufind *UFcluster);

void copyPair(struct pair *a,struct pair b);
void getESTtoMem(FILE *fp,unsigned long byteOff,char *frag);

int getBufSize(struct cBuffer buf);
void interpretResults(struct resPair *batch,int count,struct ufind *UF,int *UF_size);
int deQWaitQ(struct waitq *buf,int *ret);
int enQWaitQ(struct waitq *buf,int slave);
int dispAllClustersMidway(struct ufind *uf, int size, struct estgi *allEstGis, int *singletons,	int step);
int emptyBuf(struct cBuffer buf);
int emptyWaitQ(struct waitq buf);
int retrieveMinLen();
int fullBuf(struct cBuffer buf);
int determineA2AIterations(unsigned int g_maxBuckSize);
int deQBuf(struct cBuffer *buf,struct pair *ret);
void printAPair(struct pair p,FILE *fp);
int enQBuf(struct cBuffer *buf,struct pair pr);

#endif
