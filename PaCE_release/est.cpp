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


#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include "est.h"
#include "err.h"
#include "uFind.h"
#include "debug.h"

const char *tag=">";

/* REturn the Reverse Comp of the frag in rFrag */
void getRComp(char *frag, char *rFrag)
{
    int i=0,iLen=0;
    char c;

    iLen=strlen(frag);

    for(i=iLen-1;i>=0;i--)
    {
      switch(frag[i])
      {
        case 'A' :  c='T'; break;
        case 'T' :  c='A'; break;
        case 'G' :  c='C'; break;
        case 'C' :  c='G'; break;
        case 'N' :  c='N'; break;
        case 'U' :  c='A'; break;
        case 'R' :  c='Y'; break;
        case 'Y' :  c='R'; break;
        case 'K' :  c='M'; break;
        case 'M' :  c='K'; break;
        case 'S' :  c='W'; break;
        case 'W' :  c='S'; break;
        case 'B' :  c='V'; break;
        case 'V' :  c='B'; break;
        case 'D' :  c='H'; break;
        case 'H' :  c='D'; break;
        default  :
                   /* printf("\nError: getRComp : Frag=%s=%d Error Detected in format here.\n",frag,iLen);
                    return; */
                    c='N'; break;
      } /* end switch */

      rFrag[iLen-1-i]=c;
    } /* end for*/
    rFrag[iLen]='\0';
} /* end getRComp*/



int loadESTs(struct est *localEsts,int N,char *estfile,
				struct FileOffset *loc_fo, struct FileSize *glob_fs, int nFilesToRead,
				int k,struct FileOffset *ESTloc,int *ESTlen,struct loadRes *res) {
	FILE *fp;
	int bFirst=1;
	int iReadingFrag=0,iNoLocalFrags=0;
	char frag[MAXFRAGSIZE],fragstr[MAXFRAGSIZE],tempfrag[MAXFRAGSIZE];
	int fraglen=0;
	int thislen=0;
	char lastHeader[MAXFRAGSIZE];

	long b4Offset=0,aftOffset=0;
	long lastOffset=0;

	struct timeval t_1,t_2;
	unsigned long t_IO=0;
	unsigned long t_sIO=0;
	long n_bytesIORead=0;

	int computeStartFragId(int *recv);
	int sumUpAllFrags(int *rec);

	if(!glob_fs || !loc_fo) {
		printf("Rank=%d: loadESTs Error: global_fs or local_fo empty !\n",rank);
		fflush(stdout);
		terminate();
	}

  	res->startEST = 0;
  	res->endEST = 0;
	res->myFragsLen=0;
	res->Longest = 0;
	res->numLenLessK=0;
	res->sizLenLessK=0;
	iNoLocalFrags=0;


	int fId;
	int bfId,efId;
	long bOff,eOff;
	long iOff;

	bfId=loc_fo[0].fId;
	bOff=loc_fo[0].Off;

	efId=loc_fo[1].fId;
	eOff=loc_fo[1].Off;

	fId = bfId;
	iOff = bOff;

	if(bfId<0 || efId<0 ) {
		printf("Rank=%d: Warning loadESTs: Empty list of local sequences !\n",rank);
		fflush(stdout);
		return 0;
	}
		

  	while(1){
		// loop variants: fId and iOff

		if(fId>=nFilesToRead) {
			//printf("Rank=%d: Warning: Last file read without prior notice! \n",rank);
			break;
		}
		
		struct FileSize *g_fs;

		g_fs = &(glob_fs[fId]);
		fp = fopen(g_fs->name,"r");
  		fseek(fp,iOff,SEEK_SET);
  		if(feof(fp) || fp==NULL) {
     		printf("Rank=%d: Error : EST data file (%d, %ld) possibly empty or corrupted \n",
								rank,fId,iOff);
			fflush(stdout);
		     terminate();
		}

		bFirst=1;
		iReadingFrag=0;
		int bDone=0;
		strcpy(fragstr,"");
		strcpy(lastHeader,"");

		//printf("Rank=%d: Start reading file %d with %ld\n",rank,fId,ftell(fp));

		while(!feof(fp)) {

    			/* this will read even the frag boundary */
    			strcpy(frag,"");
    			b4Offset= ftell(fp);
			if(fId==efId && b4Offset>=eOff) {
				bDone=1;
				break;
			}
			
    			lastOffset=b4Offset;
    			//if(lastOffset>=endOffset) break;

    			gettimeofday(&t_1,0);
    			fscanf(fp,"%[^\n]%*\n",frag);
    			gettimeofday(&t_2,0);

    			aftOffset= ftell(fp);
    			if(b4Offset==aftOffset || (strcmp(frag,"")==0) ) {
       			if(b4Offset==aftOffset) {
    					b4Offset= ftell(fp);
	       			fseek(fp,b4Offset+1L,SEEK_SET);
       			}
       			continue;
    			}

    			fraglen=strlen(frag);

    			t_IO+= (t_2.tv_sec-t_1.tv_sec)*1000000 + (t_2.tv_usec-t_1.tv_usec);
			if(t_IO>=MAXULONGBY2) {
                	t_sIO+=t_IO/1000000;
                	t_IO=0;
    			}
    			n_bytesIORead+=strlen(frag)+1;

    			if(bFirst) {
       			if(!strstr(frag,tag)) {
	  				printf("Rank=%d: loadESTs(): First line not a header frag!##%s##\n",
							rank,frag);
	  				fflush(stdout);
          			continue;
       			}
    			} /* end if */

    			if(strstr(frag,tag)) {
       			strcpy(lastHeader,frag);
       			/* this means its the next frag - so make sure you store the prev frag*/
       			if(iReadingFrag) {
         				thislen = strlen(fragstr);
         				needs = sizeof(char)*(thislen+1);
         				localEsts[iNoLocalFrags].s = (char *) malloc(needs);
         				if(localEsts[iNoLocalFrags].s) { 
						thisPhaseSpace+=needs; 
						spaceSoFar+=needs;
					}
         				else AlertErr("localEsts[].s");

         				strcpy(localEsts[iNoLocalFrags].s,fragstr);
         				if(ESTlen) ESTlen[iNoLocalFrags]=strlen(fragstr);
         				iNoLocalFrags++;
         				if(thislen>res->Longest) res->Longest=thislen;
         				if(thislen<k) {
            				res->numLenLessK+=1;
            				res->sizLenLessK+=thislen;
         				}

         				if(ESTloc) {
						ESTloc[iNoLocalFrags].fId = fId;
						ESTloc[iNoLocalFrags].Off =ftell(fp);   
					}
					/* set to next fragid start*/   
       			}
       			else {
          			if(!bFirst) {
             				needs = sizeof(char);
             				localEsts[iNoLocalFrags].s = (char *) malloc(needs);
             				if(localEsts[iNoLocalFrags].s) { 
							thisPhaseSpace+=needs; spaceSoFar+=needs;
						}
             				else AlertErr("localEsts[].s");

             				strcpy(localEsts[iNoLocalFrags].s,"");
             				if(ESTlen) ESTlen[iNoLocalFrags]=0;
             				iNoLocalFrags++;
			          	res->numLenLessK+=1;

						/* set to -1 to avoid consulting the file at all for empty frag*/   
             				if(ESTloc) {
							ESTloc[iNoLocalFrags].fId = -1;
							ESTloc[iNoLocalFrags].Off = 0;
						}
          			}
          			else {
               			/*printf("Rank=%d: First Header =%s\n",rank,frag); fflush(stdout);   */
               			bFirst=0;
               			if(ESTloc) {
							ESTloc[iNoLocalFrags].fId = fId;
							ESTloc[iNoLocalFrags].Off = ftell(fp);
						}
          			}
       			}
       			strcpy(fragstr,"");
       			iReadingFrag=0;
    			}
    			else {
        			/* fragment string - so store it temporarily */
        			strcpy(tempfrag,"");
        			sscanf(frag,"%[A-Za-z]",tempfrag);
        			strcat(fragstr,tempfrag);
        			iReadingFrag=1;
        			res->myFragsLen+=fraglen;  
    			}

		} // end inner while
		//printf("Rank=%d: Finish reading file %d with %ld\n",rank,fId,ftell(fp));
		fclose(fp);

  		/* Write out the last fragment */
       	if(iReadingFrag) {
         		thislen = strlen(fragstr);
         		needs = sizeof(char)*(thislen+1);
         		localEsts[iNoLocalFrags].s = (char *) malloc(needs);
         		if(localEsts[iNoLocalFrags].s) { 
				thisPhaseSpace+=needs; 
				spaceSoFar+=needs;
			}
         		else AlertErr("localEsts[].s");

         		strcpy(localEsts[iNoLocalFrags].s,fragstr);
         		if(ESTlen) ESTlen[iNoLocalFrags]=strlen(fragstr);
         		if(thislen>res->Longest) res->Longest=thislen;
         		if(thislen<k) {
            		res->numLenLessK+=1;
            		res->sizLenLessK+=thislen;
         		}
         		iNoLocalFrags++;
       	}


		if(bDone==1) break;
		// read from the next file
		fId++;
		iOff=0;

	} // end outer while

	printf("Rank=%d: Local IO stat: %d sequences (%ld bytes) read: fscanf %u s \n",
		rank,iNoLocalFrags,n_bytesIORead,
	t_sIO+(t_IO/1000000));
	printf("Rank=%d: Number of fragments =%d \n",rank,iNoLocalFrags);

  	return iNoLocalFrags;
}/* end loadESTs */


/* this calculates the global frag id of the starting fragment in the local array*/
int computeStartFragId(int *recv)
{
  int i=0;
  int iFragsB4Me=0;

  for(i=0;i<rank;i++)
  {
    iFragsB4Me+=recv[i];
  }

  return iFragsB4Me;

} /* end computeStartFragId*/

/* keeps track of the total number of fragments in the global array*/
int sumUpAllFrags(int *rec)
{
  int i=0;
  int totalFrags=0;

  for(i=0;i<p;i++)
  {
    /*if(rank==0)  printf("Recv[%d]=%d\n",i,rec[i]);*/
    totalFrags+= rec[i];
  }

  return totalFrags;
} /*  end sumUpAllFrags*/


void dispESTs(struct est *ESTs, int size,int N)
{
 int i;
  FILE *fp;
 char s[20];

  sprintf(s,"est#.%d.%d",N,rank);
  fp = fopen(s,"w");
  for(i=0;i<size;i++)
  {
     fprintf(fp,"%d:\n%s\n",i,ESTs[i].s);
  }
  fprintf(fp,"\n");
  fclose(fp);
} /* end dispESTs */


void dispESTGIs(struct estgi *ESTgis, int size)
{
 int i;
  FILE *fp;
 char s[20];

  sprintf(s,"estgi#.%d",size);
  fp = fopen(s,"w");
  for(i=0;i<size;i++)
  {
     fprintf(fp,"%d %s\n",i,ESTgis[i].gi);
  }
  fprintf(fp,"\n");
  fclose(fp);
} /* end dispESTGIs */


/* loadESTGIs :  This function maps each EST to its header.  
 * PS:  This need not be a GI header. This will work for any FASTA format header.
 * This function is used in reporting the final set of clusters and outputing in CF format. */

void loadESTGIs(struct estgi *allEstGis,int N,
				struct FileSize *glob_fs,int nFilesToRead) {
  	FILE *fp;
  	int iNoLocalFrags=0;
  	char frag[MAXFRAGSIZE];
  	char firstword[MAXFRAGSIZE];


	if(!glob_fs || nFilesToRead<=0) {
		printf("Rank=%d: loadESTGIs: global_fs empty !\n",rank);
		fflush(stdout);
		terminate();
	}

	iNoLocalFrags=0;

	int fId;
	long iOff;

	fId = 0;
	iOff=0;

  	while(1){
		// loop variants: fId and iOff

		if(fId>=nFilesToRead) {
			break;
		}
		
		struct FileSize *g_fs;

		g_fs = &(glob_fs[fId]);
		fp = fopen(g_fs->name,"r");
  		fseek(fp,iOff,SEEK_SET);
  		if(feof(fp) || fp==NULL) {
     		printf("Rank=%d: Error : EST data file (%d, %ld) possibly empty or corrupted \n",
								rank,fId,iOff);
			fflush(stdout);
		     terminate();
		}

		int bDone=0;
		int bReading=0;
		int bFirst=1;

		while(!feof(fp)) {

    			/* this will read even the frag boundary */
    			strcpy(frag,"");
    			fscanf(fp,"%[^\n]%\n",frag);
    			sscanf(frag,"%s",firstword);
    			if(bFirst) {
       			/* first iteration - make sure read from the next frag header */
       			if(!strstr(frag,tag)) {
          			printf("Warning in input file from GI: First line not a header \n");
          			continue;
       			}
       			bFirst =0;
    			}

    			if(strstr(frag,tag)) {
       			strcpy(allEstGis[iNoLocalFrags].gi,firstword);
       			iNoLocalFrags++;
    			}
    			else continue;

  		} /* inner while end */
  		fclose(fp);
		fId++;
		iOff=0;
	} // outer while end


  	if(iNoLocalFrags!=N) {
    		printf("Master Error: Supplied number <%d> of ESTs does not match with the number in the input file <%d>\n",N,iNoLocalFrags);
    	//	terminate();
  	}

  	#ifdef debug
  		printf("Rank=%d: Number of fragments =%d \n",rank,iNoLocalFrags);
  		fflush(stdout);
  	#endif

}/* end loadESTGIs */



/* Loads just the GI numbers */
/* loadGIs :  This function maps each EST to its GI number in its header.  
 * PS:  This works for only restricted cases where the FASTA header is of the format >gi|878787|...
 * This function is used while updating the set of clusters with clone pair information */


void loadGIs(int *allEstGis,int N,char *estfile,int *smallest,int *largest)
{
  FILE *fp;
  int bFirst=1;
  int iNoLocalFrags=0;
  char frag[MAXFRAGSIZE],header[MAXFRAGSIZE];
  int thisgi=0;


  fp = fopen(estfile,"r");
  if(feof(fp) || fp==NULL)
  {
     printf("Error : EST data file possibly empty or corrupted.\n");
     exit(1);
  }

  *smallest=SOMEBIGNUMBER;
  *largest=0;

  while(!feof(fp))
  {
    /* this will read even the frag boundary */
    fscanf(fp,"%[^\n]%\n",frag);
     /* printf("Frag=%s\n",frag); 
     fflush(stdout);*/
    if(bFirst)
    {
       /* first iteration - make sure read from the next frag header */
       if(!strstr(frag,tag))
       {
          printf("Error in input file from GI: First line not a header \n");
          exit(1);
       }
       bFirst =0;
    }

    if(strstr(frag,tag))
    {
       sscanf(frag,"%*[^|]%*[|]%[^|]%*s",header);
       thisgi=atoi(header);
       allEstGis[iNoLocalFrags]=thisgi;
       iNoLocalFrags++;

       if(thisgi<*smallest) *smallest=thisgi;
       if(thisgi>*largest) *largest=thisgi;

    }
    else
    {
       continue;
    }
  } /* while end */

  if(iNoLocalFrags!=N)
  {
    printf("Master Error: Supplied number <%d> of ESTs does not match with the number in the input file <%d>\n",N,iNoLocalFrags);
    exit(1);
  }

  fclose(fp);
  #ifdef debug
  printf("Rank=%d: Number of fragments =%d \n",rank,iNoLocalFrags);
  fflush(stdout);
  #endif

}/* end loadGIs */

void doClonePairs(struct ufind *estcluster,int N,char *estFile)
{
       int *allEstGis,*rmapGi;
       int i;
       int smallest=SOMEBIGNUMBER,largest=0;
       int rmapsize=0;
       FILE *fp=NULL;
       char word[50];
       int lastGi=-1,thisgi=0,f1=0,f2=0;
       char cpfile[200];
       void getCloneMatesFile(char *);


       getCloneMatesFile(cpfile);
       if(!strcmp(cpfile,"")) return;

       fp=fopen(cpfile,"r");
       if(!fp) 
       {
         printf("Master Warning: Clone Pairs file %s is possibly missing!\n",cpfile);
         return;
       }

       needs = sizeof(int)*N;
       allEstGis= (int *)malloc(needs);
       checkAlloc(allEstGis,"allEstGis in doClonePairs");

       loadGIs(allEstGis,N,estFile,&smallest,&largest);

       if(largest<smallest) 
       {
          printf("Master: Error in clonePairs\n");
          exit(1);
       }
   
       rmapsize=largest-smallest+1;
       needs = sizeof(int)*rmapsize;
       rmapGi = (int *)malloc(needs);
       checkAlloc(allEstGis,"rmapGi in doClonePairs");

       for(i=0;i<rmapsize;i++)  rmapGi[i]=-1; /* invalid GI */

       /* reverse Map Gi numbers to their EST index */
       for(i=0;i<N;i++)  rmapGi[allEstGis[i]-smallest]=i;

       lastGi=-1;
       while(!feof(fp))
       {
          fscanf(fp,"%s",word);
          thisgi=atoi(word);
          if(thisgi<=0) 
          {
             lastGi=-1;
             continue;
          }
          
          if(lastGi<0)
          {
             lastGi=thisgi; continue;
          }
   
          if(lastGi<smallest || lastGi>largest) continue;
          if(thisgi<smallest || thisgi>largest) continue;

          f1 = rmapGi[lastGi-smallest];
          f2 = rmapGi[thisgi-smallest];
          if(f1<0 || f2<0) continue;

          Union(estcluster,f1,f2);
       } /* end while */

       if(allEstGis!=NULL) /* Ananth Change: for AiX/Sri 1/29/2004 */
         free(allEstGis);
       if(rmapGi!=NULL) /* Ananth Change: for AiX/Sri 1/29/2004 */
         free(rmapGi);
 
} /* end doClonePairs */



/* Build MPI_Datatype for struct FileOffset */
int buildFileOffsetMPI(struct FileOffset *t, MPI_Datatype *dt) {
 int block_lengths[2];
 MPI_Datatype type_list[2];
 MPI_Aint saddr,addr;
 MPI_Aint disp[2];

 block_lengths[0]=block_lengths[1]=1;

 type_list[0]=MPI_INT;
 type_list[1]=MPI_LONG;

 disp[0]=0;
 MPI_Address(&t->fId,&saddr);
 MPI_Address(&t->Off,&addr);
 disp[1]=addr-saddr;

 MPI_Type_struct(2,block_lengths,disp,type_list,dt);
 MPI_Type_commit(dt);
 return 1;
}

int loadTreeStrings(char *fastafile,
			 struct FileSize *glob_fs, int nFilesToRead,
			unsigned int *treeStringPtr,
			char *treeStrings,
			char *treeStringMarker,
			int maxupto) {

	FILE *fp;
	int iNoLocalFrags=0;
	char frag[MAXFRAGSIZE],fragstr[MAXFRAGSIZE],tempfrag[MAXFRAGSIZE];
	int fraglen=0;

	long b4Offset=0,aftOffset=0;
	long lastOffset=0;

	struct timeval t_1,t_2;
	unsigned long t_TreeScan=0;
	unsigned long t_sTreeScan=0;
	unsigned int n_bytesTreeReadFromIO=0;
	unsigned int str_it=0;
	int iLoadedStrings=0;


	if(!glob_fs || nFilesToRead<=0) {
		printf("Rank=%d: loadTreeStrings Error: global_fs empty !\n",rank);
		fflush(stdout);
		terminate();
	}

	iNoLocalFrags=-1;

	int fId;
	long iOff;

	fId = 0;
	iOff=0;
	str_it=0;
	iLoadedStrings=0;

  	while(1){
		// loop variants: fId and iOff

		if(fId>=nFilesToRead) {
			//printf("Rank=%d: Warning: Last file read without prior notice! \n",rank);
			break;
		}
		
		struct FileSize *g_fs;

		g_fs = &(glob_fs[fId]);
		fp = fopen(g_fs->name,"r");
  		fseek(fp,iOff,SEEK_SET);
  		if(feof(fp) || fp==NULL) {
     		printf("Rank=%d: Error : EST data file (%d, %ld) possibly empty or corrupted \n",
								rank,fId,iOff);
			fflush(stdout);
		     terminate();
		}

		int bDone=0;
		strcpy(fragstr,"");
		int bReading=0;
		int bFirst=1;

		//printf("Rank=%d: Start reading file %d with %ld\n",rank,fId,ftell(fp));

		while(!feof(fp)) {

    			/* this will read even the frag boundary */
    			strcpy(frag,"");
    			b4Offset= ftell(fp);
			
    			lastOffset=b4Offset;

    			gettimeofday(&t_1,0);
    			fscanf(fp,"%[^\n]%*\n",frag);
    			gettimeofday(&t_2,0);

    			aftOffset= ftell(fp);
    			if(b4Offset==aftOffset || (strcmp(frag,"")==0) ) {
       			if(b4Offset==aftOffset) {
    					b4Offset= ftell(fp);
	       			fseek(fp,b4Offset+1L,SEEK_SET);
       			}
       			continue;
    			}

    			fraglen=strlen(frag);

    			t_TreeScan+= (t_2.tv_sec-t_1.tv_sec)*1000000 + (t_2.tv_usec-t_1.tv_usec);
			if(t_TreeScan>=MAXULONGBY2) {
                	t_sTreeScan+=t_sTreeScan/1000000;
                	t_TreeScan=0;
    			}
    			n_bytesTreeReadFromIO+=strlen(frag)+1;

    			if(bFirst) {
       			if(!strstr(frag,tag)) {
	  				printf("Rank=%d: loadTreeStrings(): First line of %s not a header frag!##%s##\n",
							rank,g_fs->name,frag);
	  				fflush(stdout);
          			continue;
       			}
    			} /* end if */

			bFirst=0;

    			if(strstr(frag,tag)) {
       			/* this means its the next frag - so make sure you store the prev frag*/
       			if(bReading==1) {
					if(treeStringMarker[iNoLocalFrags]=='1') {
                			iLoadedStrings++;
	          			strcpy(&(treeStrings[str_it]),fragstr);
     	     			treeStringPtr[iNoLocalFrags]=str_it;
          				str_it+=strlen(fragstr)+1;
          				if(strlen(fragstr)<=0)  { /* empty frag*/
               				printf("Rank=%d: Warning_1:  loadTreeStrings(): Empty Sequence %d (maxupto=%d)!\n",
                    				rank,iNoLocalFrags,maxupto);
          				}
       				}
				}
       			iNoLocalFrags++;
				bReading=0;
       			strcpy(fragstr,"");
       			if(iNoLocalFrags>maxupto) {
					bDone=1;
					break; // out of the inner loop
				}
			}
    			else {
       			if(treeStringMarker[iNoLocalFrags]!='1') continue;
	         		/* fragment string - so store it temporarily */
				strcpy(tempfrag,"");
				sscanf(frag,"%[A-Za-z]",tempfrag);
				strcat(fragstr,tempfrag);
				bReading=1;
    			}

		} // end inner while
		//printf("Rank=%d: Finish reading file %d with %ld\n",rank,fId,ftell(fp));
		fclose(fp);

		if(bDone==1) break; // out of the outer loop

  		/* Write out the last fragment if necessary */
     	if(bReading==1) {
			if(treeStringMarker[iNoLocalFrags]=='1') {
		     	iLoadedStrings++;
				strcpy(&(treeStrings[str_it]),fragstr);
				treeStringPtr[iNoLocalFrags]=str_it;
				str_it+=strlen(fragstr)+1;
				if(strlen(fragstr)<=0)  { /* empty frag*/ 
					printf("Rank=%d: Warning_2:  loadTreeStrings(): Empty Sequence %d (maxupto=%d)!\n", 		
								rank,iNoLocalFrags,maxupto);
				}
			}
		}

		// read from the next file
		fId++;
		iOff=0;

	} // end outer while

	iNoLocalFrags++;

	printf("Rank=%d: loadTreeStrings: IO stat: %d sequences (%u bytes) read: fscanf %u s \n",
				rank,iNoLocalFrags,n_bytesTreeReadFromIO,
				t_sTreeScan+(t_TreeScan/1000000));
	printf("Rank=%d: loadTreeStrings: Number of fragments =%d occupying %u bytes \n",
				rank,iLoadedStrings,str_it);
	fflush(stdout);

	return iNoLocalFrags;

}/* end loadTreeStrings*/

