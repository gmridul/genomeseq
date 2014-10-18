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





/*          mastercode.c:  Application specific methods require 
 *          complete implementation by the user of this code */



#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include "mastercode.h"
#include "err.h"
#include "debug.h"

const char *comment="#";
char g_printpairstr[PrintPairBatchSize*PrintPairOneLine];

int *UFcluster_size;


void processPair(struct dynResult *result,struct pair rBuf,struct dynParams param,char *s1,char *s2) {
   alignPair(result,rBuf.f1,rBuf.f2,s1,s2,rBuf.p1,rBuf.p2,rBuf.l1,rBuf.l2,rBuf.type,rBuf.d,param);
} /* end processPair */



/*  Applicaton specific: Displays an array of pair structures 
*/ 

void dispPairs(struct pair *pairs, int size)
{
 int i;
 for(i=0;i<size;i++)
  {
     printf("(%d,%d) ",pairs[i].f1,pairs[i].f2);
  }
 printf("\n");
} /* end dispPairs */


void dispPairsDetail(struct pair *p, int size)
{
 int i;
 for(i=0;i<size;i++)
  {
     printf("(%d,%d- %d,%d,%d,%d,%d,%c) ",p[i].f1,p[i].f2,p[i].p1,p[i].p2,p[i].l1,p[i].l2,p[i].d,p[i].type);
  }
 printf("\n");
} /* end dispPairsDetail */



/*
    interpretResults : Update the cluster based on the ratio of 
                        bandScore Vs idealScore and cluster
			based on gene or transcript homology
*/     
void interpretResults(struct resPair *batch,int count,struct ufind *UF,int *UF_size) {
  int i ;
  int mbScore,ibScore,mScore,iScore;
  float borderratio,maxscoreratio;
  float borderalignlen=0;
  int b_geneHomology,b_transcriptHomology,acceptCondition =0;
  int size1,size2;
  void printAccepted(int f1,int f2, char type,int size1,int size2,FILE *fp);
  void dispresPair(struct resPair st);

  for(i=0;i<count;i++)
  {
     /*printf("Master: interpretResults(%d,%d)\n",batch[i].f1,batch[i].f2); fflush(stdout);*/
     if(Find(UF,batch[i].f1)==Find(UF,batch[i].f2)) continue;

    /* printf("Master: interpretResults(%d,%d) checking ... \n",batch[i].f1,batch[i].f2); fflush(stdout);
     dispresPair(batch[i]);*/
     mbScore = batch[i].maxBorderScore;
     ibScore = batch[i].idealBorderScore;
     mScore = batch[i].maxScore;
     iScore = batch[i].idealScore;

   /*    printf("Master: #%d,%d,%d,%d#\n",mbScore,ibScore,mScore,iScore); fflush(stdout);*/

     borderratio = abs(ibScore-mbScore)*100/ibScore;
     borderalignlen = ibScore/g_matchscore; /* change for borderalign fix to reflect min-align-len*/
     /*borderalignlen = iScore/g_matchscore;*/

     maxscoreratio = abs(iScore-mScore)*100/iScore;

     if(borderratio<=EndToEndScoreRatioThreshold && borderalignlen>=EndtoEndAlignLenThreshold) b_transcriptHomology=1;
     else b_transcriptHomology=0;

     if(maxscoreratio<=MaxScoreRatioThreshold && batch[i].maxScoreCoverage>=TranscriptCoverageThreshold) b_geneHomology=1;
     else b_geneHomology=0;


     if(g_bTranscriptsTogether)
     {
        /* clusters alternative transcript isoforms together */
     	if( b_transcriptHomology || b_geneHomology ) acceptCondition=1;
	else acceptCondition=0;

	if(g_bReportSplicedCandidates && !b_transcriptHomology && b_geneHomology) {
		fprintf(fpSpliced,"Spliced_Candidate  %d         %d\n",batch[i].f1,batch[i].f2);
	}

     }
     else
     {
        /* intends separates alternative transcript isoforms */
     	if(b_transcriptHomology) acceptCondition=1;
	else acceptCondition=0;
     }

     if(acceptCondition) {
     	iPairsAccepted++;
		size1=UF_size[Find(UF,batch[i].f1)];	
		size2=UF_size[Find(UF,batch[i].f2)];	

		if(g_bOutputLargeMerges) {
			if((size1 >= LargeClusterThreshold) && (size2 >= LargeClusterThreshold)) {
				printAccepted(batch[i].f1,batch[i].f2,batch[i].type,
					size1,size2,fpLargeMerges);
			}
		}



	
       	Union(UF,batch[i].f1,batch[i].f2);
		UF_size[Find(UF,batch[i].f1)]=size1+size2;
       	/* Detect and mark containments */
       	/* Ananth change - Jul 1 2003, make this containment condition more stringent 
		* This is currently hardcoded to 0% but parameterize this */

		if(g_bReportAcceptedPairs)  {
			printAccepted(batch[i].f1,batch[i].f2,batch[i].type,
				size1,size2,fpAcceptedPairs);
		}

       	/*if(batch[i].type=='C') */ 
       	if(batch[i].type=='C' && borderratio<=0) {
         		if(batch[i].atag=='m')  containArr[batch[i].f1]=batch[i].f2;
         		else containArr[batch[i].f2]=batch[i].f1;
       	}		
       	/* printf("Rank=%d: Pair(%d,%d) Accept with %d%% and type=%c,%c\n",rank,batch[i].f1,batch[i].f2,ratio,batch[i].type,batch[i].atag); 
       	fflush(stdout); */
       
       	continue;
     }

     //iPairsRejected++;
     /* printf("Rank=%d: Pair(%d,%d) Reject with %d%% type=%c,%c \n",rank,batch[i].f1,batch[i].f2,ratio,batch[i].type,batch[i].atag); 
     fflush(stdout); */
     
  }

} /* end interpretResults1 */

/* insertCBuf : master loop calls this to insert into cBuffer */
/*              returns the number of inserted pairs out of nPairs */
/* 	  	This modifies only cBuf->end and cBuf->currSize */

int insertCBuf(struct pair *newPrs,int nPairs,struct cBuffer *cBuf,struct ufind *UFcluster)
{
  int i=0,ass=0;
  int f1,f2;
  void printAPair(struct pair p,FILE *fp);
  int enQBuf(struct cBuffer *buf,struct pair pr);

  ass=0;
  for(i=0;i<nPairs;i++)
  {
    f1 = newPrs[i].f1;
    f2 = newPrs[i].f2; 
    if( (f1<0 || f1>=N) || (f2<0 || f2>=N)) {
    	printf("Rank=%d: Error: insertCBuf f1 or f2 out of range %d %d\n",rank,f1,f2);
	fflush(stdout);
    //   MPI_Abort(MPI_COMM_WORLD,1);
    	continue;
    }
    if(f1==f2) continue;
    if(g_bReportMaximalPairs) printAPair(newPrs[i],fpMaximalPairs);
    
    if(Find(UFcluster,f1)==Find(UFcluster,f2))
     continue;
    if(containArr[f1]>=0 || containArr[f2]>=0)
     continue;

    if(enQBuf(cBuf,newPrs[i])) ass++;
    else return -1;
  }
  return ass;
} /* end insertCBuf */


/* copyPair : copies pair b to *a */
void copyPair(struct pair *a,struct pair b)
{
  a->f1=b.f1;
  a->f2=b.f2;
  a->l1=b.l1;
  a->l2=b.l2;
  a->p1=b.p1;
  a->p2=b.p2;
  a->d=b.d;
  a->e=b.e;
  a->type=b.type;
} /* end copyPair */


/* prepareNextBatch : prepares next batch of pairs + the expected pairs for */
/*                    a src.  It modifies cBuf->start and cBuf->currSize  */
/*	              returns #assigned pairs from buffer to nxtWork */
/*		               -1 otherwise 		*/

int prepareNextBatch(struct cBuffer *cBuf,struct pair *batch,struct ufind *UFcluster)
{
   int j=0;
   int f1,f2;
   int emptyBuf(struct cBuffer buf);
   int deQBuf(struct cBuffer *buf,struct pair *ret);

   if(emptyBuf(*cBuf)) return 0;
   j=0;
   while(!emptyBuf(*cBuf)&& j<LOADPERPROC)
   {
       deQBuf(cBuf,&(batch[j]));

      /* Check if work is necessary even while distributing the pairs
	*/
       f1=batch[j].f1;
       f2=batch[j].f2;

       if(Find(UFcluster,f1)==Find(UFcluster,f2)) {
		/* No need to distribute this pair.. and do not incr. j
		  	because batch[j] will be overwritten or not used if last
		*/
		continue;
	}
	/* similarly do not allocate if either of the fragments is contained*/
        if(containArr[f1]>=0 || containArr[f2]>=0) continue;
	
       j++;
   } /* end while */
   return j;
} /* end prepareNextBatch */


/* Build MPI_Datatype for struct pair*/
void buildPairMPI(struct pair *t, MPI_Datatype *dt)
{
 int block_lengths[10];
 MPI_Datatype type_list[10];
 MPI_Aint saddr,addr;
 MPI_Aint disp[10];
 int i=0;

 for(i=0;i<10;i++)
   block_lengths[i] = 1;

 type_list[0]=MPI_CHAR;
 for(i=1;i<8;i++)
   type_list[i] = MPI_INT;

 type_list[8] = MPI_FLOAT;
 type_list[9] = MPI_CHAR;

 disp[0]=0;
 MPI_Address(&(t->dummy),&saddr);
 MPI_Address(&(t->f1),&addr);
 disp[1]=addr-saddr;
 MPI_Address(&(t->f2),&addr);
 disp[2]=addr-saddr;
 MPI_Address(&(t->l1),&addr);
 disp[3]=addr-saddr;
 MPI_Address(&(t->l2),&addr);
 disp[4]=addr-saddr;
 MPI_Address(&(t->p1),&addr);
 disp[5]=addr-saddr;
 MPI_Address(&(t->p2),&addr);
 disp[6]=addr-saddr;
 MPI_Address(&(t->d),&addr);
 disp[7]=addr-saddr;
 MPI_Address(&(t->e),&addr);
 disp[8]=addr-saddr;
 MPI_Address(&(t->type),&addr);
 disp[9]=addr-saddr;


 MPI_Type_struct(10,block_lengths,disp,type_list,dt);
 MPI_Type_commit(dt);
} /* end buildPairMPI */


/* Build MPI_Datatype for struct resPair*/
void buildresPairMPI(struct resPair *t, MPI_Datatype *dt)
{
 int block_lengths[9];
 MPI_Datatype type_list[9];
 MPI_Aint saddr,addr;
 MPI_Aint disp[9];
 int i=0;

 for(i=0;i<9;i++)
   block_lengths[i] = 1;

 /*type_list[0]=MPI_CHAR;*/
 for(i=0;i<=6;i++)
   type_list[i] = MPI_INT;

 type_list[7] = MPI_CHAR;
 type_list[8] = MPI_CHAR;

 disp[0]=0;
 /*MPI_Address(&t->dummy,&saddr);*/
 MPI_Address(&t->maxBorderScore,&saddr);
 /*disp[1]=addr-saddr;*/
 MPI_Address(&t->idealBorderScore,&addr);
 disp[1]=addr-saddr;
 MPI_Address(&t->maxScore,&addr);
 disp[2]=addr-saddr;
 MPI_Address(&t->idealScore,&addr);
 disp[3]=addr-saddr;
 MPI_Address(&t->maxScoreCoverage,&addr);
 disp[4]=addr-saddr;
 MPI_Address(&t->f1,&addr);
 disp[5]=addr-saddr;
 MPI_Address(&t->f2,&addr);
 disp[6]=addr-saddr;
 MPI_Address(&t->type,&addr);
 disp[7]=addr-saddr;
 MPI_Address(&t->atag,&addr);
 disp[8]=addr-saddr;

 MPI_Type_struct(9,block_lengths,disp,type_list,dt);
 MPI_Type_commit(dt);
} /* end buildresPairMPI */

/* Build MPI_Datatype for struct work - Planned: Obsolete*/
/*int buildWorkMPI(struct work *t, MPI_Datatype *dt)
{
 int block_lengths[4];
 MPI_Datatype type_list[4];
 MPI_Aint saddr,addr;
 MPI_Aint disp[4];
 int i=0,n;
 struct pair n1;
 struct resPair n2;
 MPI_Datatype pType,rType;

 n=4;
 buildPairMPI(&n1,&pType);
 buildresPairMPI(&n2,&rType);

 block_lengths[0] = MAXNEWPAIRS;
 block_lengths[1] = MAXLOADPERPROC;
 block_lengths[2] = 1;
 block_lengths[3] = 1;

 type_list[0] = pType;
 type_list[1] = rType;
 type_list[2] = MPI_INT;
 type_list[3] = MPI_INT;

 disp[0]=0;
 MPI_Address(&(t->newPairs),&saddr);
 MPI_Address(&(t->res),&addr);
 disp[1]=addr-saddr;
 disp[2]=addr-saddr;
 disp[3]=addr-saddr;

 MPI_Type_struct(n,block_lengths,disp,type_list,dt);
 MPI_Type_commit(dt);
} */

int getLoadPerProc()
{
   FILE *fp;
   char word[50];

   fp = fopen(CFGFile,"r");
   while(!feof(fp))
   {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"LOADPERPROC"))
    {
       fscanf(fp,"%s",word);
       return atoi(word);
    }
    else
    {
       continue;
    }

   } /* end while */
 
  fclose(fp);
  printf("Error: %s does not have LOADPERPROC.\n",CFGFile);
  terminate();
  return -1;
} /* end getLoadPerProc */


/* Displays all the individual clusters
   Returns the number of individual clusters (including singletons)
   and the number of singletons is returned separately in the parameter passed
*/

int dispAllClusters(struct ufind *uf,int size,struct estgi *allEstGis,int *singletons)
{
    int i=0,iClusters=0;
    struct clusters *clust;
    struct link *temp, *temp_fix;
    int root;
    char s[200];
    FILE *fp,*fp1;

    /* add code for the alternative splicing testing */
    /*FILE *zfp;
    char sz[20];
    int zeroRoot=0;*/
    /* end of change */

    needs=sizeof(struct clusters)*size;
    clust = (struct clusters *)malloc(needs);
    if(clust) spaceSoFar+=needs;
    else AlertErr("clusters");

    for(i=0;i<size;i++)
    {
      clust[i].count=1;
      clust[i].singleton='1';
      clust[i].head=NULL;
    }

    
    for(i=0;i<size;i++)
    {
       root = Find(uf,i);
       if(i!=root)
       {
          /* add code for the alternative splicing testing */
          /*if(i==0) zeroRoot= root; */
          /* end of change */

          needs=sizeof(struct link);
          temp = (struct link *)malloc(needs);
          if(temp) spaceSoFar+=needs;
          else AlertErr("temp");

          temp->id=i;
          temp->next=clust[root].head;
          clust[root].head=temp;
          clust[root].count++;
          clust[root].singleton='0';
          clust[i].singleton = '0';
          clust[i].head = NULL;
       }
    } /* end for */

#ifndef PERF_STUDY
    sprintf(s,"%s/estClustSize.%d.%d.PaCE",g_sOutputFolder,size,rank);
    fp1 = fopen(s,"w");
	 
    sprintf(s,"%s/estClust.%d.%d.PaCE",g_sOutputFolder,size,rank); 
    fp = fopen(s,"w");
#endif

    /* add code for the alternative splicing testing */
      /*printf("printing in zeroRoot 0a\n"); fflush(stdout);
    sprintf(sz,"zeroCluster.%d.%d",size,rank);
    zfp = fopen(sz,"w");
      printf("printing in zeroRoot 0b\n"); fflush(stdout);*/
    /* end of change */

    iClusters=0;
    *singletons=0;
    for(i=0;i<size;i++)
    {
       if(clust[i].singleton=='0' && clust[i].head!=NULL)
       {

#ifndef PERF_STUDY
         fprintf(fp,"\n{Cluster#}   %d\n",iClusters);
         fprintf(fp,"\n  {Member#}  %s\n",allEstGis[i].gi);
#endif
         temp = clust[i].head;
         while(temp!=NULL)
         {
#ifndef PERF_STUDY
            fprintf(fp,"\n  {Member#}  %s\n",allEstGis[temp->id].gi);
#endif
            temp = temp->next;
         }
#ifndef PERF_STUDY
         fprintf(fp1,"%d %d\n",iClusters,clust[i].count);
#endif
         iClusters++;
       }
       else
       {
        if(clust[i].singleton=='1')
        {
#ifndef PERF_STUDY
          fprintf(fp,"\n{Cluster#}   %d\n",iClusters);
          fprintf(fp,"\n  {Member#}  %s\n",allEstGis[i].gi);
          fprintf(fp1,"%d %d\n",iClusters,clust[i].count);
#endif

          iClusters++;
          (*singletons)++;
        }

       }
    } /* end for 2 */

    fclose(fp);
    fclose(fp1);

    for(i=0;i<size;i++) {
	if(clust[i].head) { 
	    temp = clust[i].head;
	    while (temp->next != 0) {
		temp_fix = temp->next;
		free(temp);
		temp = temp_fix;
	    }
	    free(temp);
//	    free(clust[i].head);
	}
    }
    if(clust) free(clust);

    return iClusters;

} /* end dispAllClusters */

/* Displays all the individual clusters at some arbitrary point of the program
   Returns the number of individual clusters (including singletons)
   and the number of singletons is returned separately in the parameter passed
*/

int dispAllClustersMidway(struct ufind *uf,
				int size,
				struct estgi *allEstGis,
				int *singletons,	
				int step) {
    int i=0,iClusters=0;
    struct clusters *clust;
    struct link *temp;
    struct link *temp2;
    int root;
    char s[200];
    FILE *fp;

    /* add code for the alternative splicing testing */
    /*FILE *zfp;
    char sz[20];
    int zeroRoot=0;*/
    /* end of change */

    needs=sizeof(struct clusters)*size;
    clust = (struct clusters *)malloc(needs);
    if(clust) spaceSoFar+=needs;
    else AlertErr("clusters");

    for(i=0;i<size;i++)
    {
      clust[i].count=1;
      clust[i].singleton='1';
      clust[i].head=NULL;
    }

    
    for(i=0;i<size;i++)
    {
       root = Find(uf,i);
       if(i!=root)
       {
          /* add code for the alternative splicing testing */
          /*if(i==0) zeroRoot= root; */
          /* end of change */

          needs=sizeof(struct link);
          temp = (struct link *)malloc(needs);
          if(temp) spaceSoFar+=needs;
          else AlertErr("temp");

          temp->id=i;
          temp->next=clust[root].head;
          clust[root].head=temp;
          clust[root].count++;
          clust[root].singleton='0';
          clust[i].singleton = '0';
          clust[i].head = NULL;
       }
    } /* end for */

#ifndef PERF_STUDY
 /*   sprintf(s,"%s/estClustSize.%d.%d.PaCE.%d",g_sOutputFolder,size,rank,step);
    fp1 = fopen(s,"w");*/
	 
    sprintf(s,"%s/estClust.%d.%d.PaCE.%d",g_sOutputFolder,size,rank,step); 
    fp = fopen(s,"w");
#endif

    /* add code for the alternative splicing testing */
      /*printf("printing in zeroRoot 0a\n"); fflush(stdout);
    sprintf(sz,"zeroCluster.%d.%d",size,rank);
    zfp = fopen(sz,"w");
      printf("printing in zeroRoot 0b\n"); fflush(stdout);*/
    /* end of change */

    iClusters=0;
    *singletons=0;
    for(i=0;i<size;i++)
    {
       if(clust[i].singleton=='0' && clust[i].head!=NULL)
       {

#ifndef PERF_STUDY
         fprintf(fp,"\n{Cluster#}   %d\n",iClusters);
         fprintf(fp,"\n  {Member#}  %s\n",allEstGis[i].gi);
#endif
         temp = clust[i].head;
         while(temp!=NULL)
         {
#ifndef PERF_STUDY
            fprintf(fp,"\n  {Member#}  %s\n",allEstGis[temp->id].gi);
#endif
     /*       temp = temp->next;*/
            temp2 = temp->next;
	    if(temp) free(temp);
	    temp = temp2;
         }
#ifndef PERF_STUDY
         /*fprintf(fp1,"%d %d\n",iClusters,clust[i].count);*/
#endif
         iClusters++;
       }
       else
       {
        if(clust[i].singleton=='1')
        {
#ifndef PERF_STUDY
          fprintf(fp,"\n{Cluster#}   %d\n",iClusters);
          fprintf(fp,"\n  {Member#}  %s\n",allEstGis[i].gi);
          /*fprintf(fp1,"%d %d\n",iClusters,clust[i].count);*/
#endif

          iClusters++;
          (*singletons)++;
        }

       }
    } /* end for 2 */

    fclose(fp);
    /*fclose(fp1);*/

    /*for(i=0;i<size;i++) {
	if(clust[i].head) free(clust[i].head);
    }*/
    if(clust) free(clust);

    return iClusters;
    
} /* end dispAllClusters */


/*******************************************************/
/* Routines to maintain a circular buffer of struct pairs */
/* Will be used for both master and slave buffers */
/*******************************************************/

void initBuf(struct cBuffer *buf,int size)
{
/*
  if(buf->pairs) { buf->pairs=NULL;} 
*/
    
/*     Comment:
       free(buf->pairs); Ananth - fix for AiX 10/22/2003
       good joke Ananth - Jarek ;-)
*/

  needs = sizeof(struct pair)*size;
  buf->pairs = (struct pair *)malloc(needs);
  checkAlloc(buf->pairs,"buf->pairs");
  buf->totSize = size;
  buf->currSize = 0;
  buf->front=0;
  buf->rear=size-1;
} /* end initBuffer */


int emptyBuf(struct cBuffer buf)
{
   if( ((buf.rear+1)%buf.totSize)==buf.front)  return 1;
   return 0;
} /* end emptyBuf */

int fullBuf(struct cBuffer buf)
{
   if( ( ( ((buf.rear+1)%buf.totSize)+1) % buf.totSize)==buf.front)  return 1;
   return 0;
} /* end fullBuf*/

int enQBuf(struct cBuffer *buf,struct pair pr)
{
   if(fullBuf(*buf))
   {
	printf("Rank: %d Error: Buffer Full with size %d\n",rank,buf->currSize); fflush(stdout);
        return 0; 
   }
   buf->rear = (buf->rear +1) % buf->totSize;
   copyPair(&(buf->pairs[buf->rear]),pr);
   (buf->currSize)++;
   return 1;
} /* end enQBuf */


int deQBuf(struct cBuffer *buf,struct pair *ret)
{
   if(emptyBuf(*buf)) 
   {
	printf("Buffer Empty with size %d\n",buf->currSize);
        return 0; 
   }

   copyPair(ret,buf->pairs[buf->front]);
   buf->front = (buf->front + 1)%buf->totSize;
   (buf->currSize)--;
   return 1; 
} /* end deQBuf */

int getBufSize(struct cBuffer buf)
{
   return buf.currSize;
}

/*******************************************************/
/* Routines to maintain a circular buffer for the WAITQ */
/*******************************************************/

void initWaitQ(struct waitq *buf,int size)
{
  int i=0;
  if(buf->q) free(buf->q);

  needs = sizeof(int)*size;
  buf->q= (int *)malloc(needs);
  checkAlloc(buf->q,"buf->q");

  buf->totSize = size;
  buf->currSize = 0;
  buf->front=0;
  buf->rear=size-1;
  for(i=0;i<size;i++) buf->q[i]=-1;
} /* end initWaitQ*/

int emptyWaitQ(struct waitq buf)
{
   if( ((buf.rear+1)%buf.totSize)==buf.front)  return 1;
   return 0;
} /* end emptyBuf */

int fullWaitQ(struct waitq buf)
{
   if( ( ( ((buf.rear+1)%buf.totSize)+1) % buf.totSize)==buf.front)  return 1;
   return 0;
} /* end fullWaitQ */

int enQWaitQ(struct waitq *buf,int slave)
{
   if(fullWaitQ(*buf))  
   {
	printf("Error: WAITQ Full with size %d\n",buf->currSize);
        return 0; 
   }
   buf->rear = (buf->rear +1) % buf->totSize;
   buf->q[buf->rear] = slave;
   (buf->currSize)++;
   return 1;
} /* end enQWaitQ */


int deQWaitQ(struct waitq *buf,int *ret)
{
   if(emptyWaitQ(*buf)) 
   {
	printf("WAITQ Empty with size %d\n",buf->currSize);
        return 0; 
   }
   *ret = buf->q[buf->front];
   buf->front = (buf->front + 1)%buf->totSize;
   (buf->currSize)--;
   return 1; 
} /* end deQWaitQ*/

/*******************************************************/


void dispresPair(struct resPair st)
{

	printf("#ResPair:<%d,%d>: %d,%d,%d,%d,%d\n",st.f1,st.f2,st.maxBorderScore,st.idealBorderScore,st.maxScore,st.idealScore,st.maxScoreCoverage);
	fflush(stdout); 
} /* end dispresPair */


void printAPairInBatches(struct pair p,FILE *fp) {
	static short bFirst=1;
	static int nlines=0;
	char temp[PrintPairOneLine];

	if(bFirst) {
		fprintf(fp,"#Maximal_pair:	f1	f2	p1	p2	len	orientation \n\n\n");
		bFirst=0;
		strcpy(g_printpairstr,"");
		nlines=0;
	}
	/*if(fp) fprintf(fp,"maximal_pair:	%d	%d	%d	%d	%d	%c \n",p.f1,p.f2,p.p1,p.p2,p.l1,p.type);*/
	strcpy(temp,"");
	sprintf(temp,"%d	%d	%d	%d	%d	%c\n",p.f1,p.f2,p.p1,p.p2,p.l1,p.type);
	strcat(g_printpairstr,temp);
	if(fp && (++nlines)>=PrintPairBatchSize) {
		/*fwrite(g_printpairstr,sizeof(char),g_lenPrintPairStr,fp);*/
		fprintf(fp,"%s",g_printpairstr);
		strcpy(g_printpairstr,"");
		nlines=0;
	}
} /* end printAPairInBatches */


void printAPair(struct pair p,FILE *fp) {
	static short bFirst=1;
	if(bFirst) {
		fprintf(fp,"#Maximal_pair:	f1	f2	p1	p2	len	orientation \n\n\n");
		bFirst=0;
	}
	/*if(fp) fprintf(fp,"maximal_pair:	%d	%d	%d	%d	%d	%c \n",p.f1,p.f2,p.p1,p.p2,p.l1,p.type);*/
	if(fp) fprintf(fp,"%d %d %d %d %d %c\n",p.f1,p.f2,p.p1,p.p2,p.l1,p.type);
} /* end printAPair */


void printAccepted(int f1,int f2, char type,int size1,int size2,FILE *fp) {
	static short bFirst=1;
	if(bFirst) {
		fprintf(fp,"#Accepted_pair:   f1   f2   type size1     size2     merged_size\n\n\n");
		bFirst=0;
	}
	if(fp) {
		fprintf(fp,"accepted_pair:	%d	%d	%c	%d	%d	%d\n",
			f1,f2,type,size1,size2,size1+size2);
	}
} /* end printAccepted */



void printMaximalSubstring(struct pair p,int id,FILE *fp) {
	static int s_id=-1;

	if(id!=s_id) {
		fprintf(fp,"\n<%d,%d,%d>       ",p.f1,p.p1,p.l1);
		s_id=id;
	}
	fprintf(fp,"%d       %d       ",p.f1,p.f2);
	fflush(fp);
} /* end printMaximalSubstring */


unsigned int getScratchMemory() {
   FILE *fp;
   char word[50];
   unsigned int scratchMem=0;

   fp = fopen(CFGFile,"r");
   while(!feof(fp)) {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"ScratchMemory")) {
       fscanf(fp,"%s",word);
       /* Convert from MB (Mega Bytes) to Bytes*/
       scratchMem = (unsigned int) (1000000* atoi(word)); 

       if(scratchMem<=0 || scratchMem<g_partSize) {
		scratchMem = g_partSize;
	}
	if(rank==0){
		printf("Rank=%d: Info: Scratch Memory from config file = %u bytes\n",
				rank,scratchMem);
	}
	return scratchMem;
    }
    else continue;

   } /* end while */
 
  fclose(fp);
  printf("Error: %s does not have ScratchMemory.\n",CFGFile);
  terminate();
  return 0;
} /* end getScratchMemory*/
