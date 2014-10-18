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
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include <sys/time.h>

#define MAXFRAGSIZE 16000
#define OUTFILEEXTN 100
#define CUTOFFRATIO 20  /* 1/<CUTOFFRATIO> fraction of a frag has to come 
                           covered by leading/trailing A/Ts portion 
                           inorder to be chopped off */


int g_iNoLocalFrags=0;         /* keeps track of Number of local fragments */
int nStripped=0;

const char *tag = ">";


int main(int argc,char **argv)
{
void processFile(char *,int);

if(argc==3)
   processFile(argv[1],atoi(argv[2]));
else
   printf("Usage: preprocessPaCE {FASTA data file} {1: prune PolyA/T, 0 otherwise} \n");
} /*  end main*/

void processFile(char *fName,int stripFlag)
{
   int iOffset;
   struct timeval start,end;   
   int tDiff=0;
   char *outfile;
   FILE *rfp,*wfp;

   void parseFile(FILE *,FILE *,int);

   /* Start the timer */
   gettimeofday(&start,0);
   
   rfp=fopen(fName,"r");
   iOffset = 0;

   outfile = (char *) malloc(sizeof(char)*(strlen(fName)+OUTFILEEXTN));
   sprintf(outfile,"preprocessed.PaCE",fName);
   printf("Output file= %s\n",outfile);
   wfp = fopen(outfile,"w");
 
   if( fseek(rfp,iOffset,SEEK_SET)==0)
   {
     /* Step.1 : */
     /* Read the file now */
     parseFile(rfp,wfp,stripFlag);

   }
   else
   {
     printf("Error in fseek\n");
   }

   fclose(rfp);
   fclose(wfp);
   
   /* Start the timer */
   gettimeofday(&end,0);

   tDiff = 0; /* in seconds  */

   tDiff = (end.tv_sec-start.tv_sec + (end.tv_usec-start.tv_usec)/1000000 );
   printf("Total Time Taken : %d secs\n\n ",tDiff);
 

} /* processFile end */


void parseFile(FILE *rfp,FILE *wfp,int stripFlag)
{
  int i=0,bFirst=1;
  int iReadingFrag=0;
  char frag[MAXFRAGSIZE],fragstr[MAXFRAGSIZE],tempfrag[MAXFRAGSIZE];
  int numChars=0,maxlen=0;
  char b4fragstr[MAXFRAGSIZE];
  void stripTs(char *);
  void toUpperFrag(char *s);

  strcpy(fragstr,"");
   
  while(!feof(rfp))
  {
    /* this will read even the frag boundary */
    fscanf(rfp,"%[^\n]%\n",frag);
     /* printf("Frag=%s\n",frag); */
     fflush(stdout);
    if(strcmp(frag,"")==0) continue;
    /*if(strstr(frag,"Set")) continue;*/
    if(bFirst)
    {
       /* first iteration - make sure read from the next frag header */
       if(!strstr(frag,tag))
       {  
          printf("Error in input file : First line not a header \n");
          exit(0);
       }
       bFirst =0;
    }
    
    if(strstr(frag,tag))
    {
        /* this means its the next frag - so make sure you store the 
         prev fragment
        */

        if(iReadingFrag)
        {
		if(stripFlag) {
	   		strcpy(b4fragstr,fragstr);
           		stripTs(fragstr);
	   		if(strcmp(b4fragstr,fragstr)) {
		   		nStripped++;
	   		}
		}
           	fprintf(wfp,"%s\n",fragstr);
	   	numChars+=strlen(fragstr);
           	if(maxlen<strlen(fragstr)) maxlen=strlen(fragstr);
           	g_iNoLocalFrags++;
        } 
        strcpy(fragstr,"");
        iReadingFrag=0;
       
        fprintf(wfp,"%s\n",frag);
    }
    else
    {
        
        /* fragment string - so store it temporarily */
        strcpy(tempfrag,"");
        sscanf(frag,"%[A-Za-z]",tempfrag);
	toUpperFrag(tempfrag);
        strcat(fragstr,tempfrag);
        iReadingFrag=1; 
    }

  } /* while end */
     

  if(!iReadingFrag)
  {
     printf("Warning : Input File has last sequence incomplete\n");
  }
  else
  {
    /* Write out the last fragment */
    if(stripFlag) {
   	strcpy(b4fragstr,fragstr);
    	stripTs(fragstr);
    	if(strcmp(b4fragstr,fragstr)) {
	   nStripped++;
    	}
    }
    fprintf(wfp,"%s\n",fragstr);
    numChars+=strlen(fragstr);
    if(maxlen<strlen(fragstr)) maxlen=strlen(fragstr);
    g_iNoLocalFrags++;
  }
 
  printf("Total Number of sequences = %d \n",g_iNoLocalFrags);
  printf("Total Number of bases = %d  bp \n",numChars);
  printf("Max length of a sequence = %d  bp \n",maxlen);
  printf("Average Number of bases per sequence = %d  bp \n",numChars/g_iNoLocalFrags);
  printf("Number of sequences modified = %d  \n",nStripped);
  fflush(stdout);
 
}/* parseFile */


void stripTs(char *frag)
{
  char temp[MAXFRAGSIZE],temp1[MAXFRAGSIZE],temp2[MAXFRAGSIZE],*str;
  int len,cutoff,i,j;
  int sawT =-1;  /* 0:A, 1:T, -1:none */

  len=strlen(frag);
  strcpy(temp,frag);
  cutoff = len/CUTOFFRATIO;
  while(1)
  {

   str = strstr(temp,"TTTTTTTTTT");
   if(!str) 
   {
      str= strstr(temp,"AAAAAAAAAA");
      if(!str) break;
      else sawT=0;
   }
   else  
   { 
      sawT=1; 
   }

   if((len-strlen(str))>=cutoff) break;
   strcpy(temp2,"");

   if(sawT==1)
      sscanf(str,"%*[T]%s",temp2); 
   else
      sscanf(str,"%*[A]%s",temp2); 
   
   strcpy(temp,temp2);
   /* printf("temp=%s\n",temp); fflush(stdout); */
  } /* end while 1 */

  sawT=-1;
  i=j=-1;
  strcpy(temp1,temp);
  while(1)
  {
   str = strstr(temp1,"AAAAAAAAAA");
   if(!str) 
   {
      str= strstr(temp1,"TTTTTTTTTT");
      if(!str) break;
      else sawT=1;
   }
   else  
   { 
      sawT=0; 
   }

   if(strlen(str)<=cutoff) 
   {
     i=strlen(temp);
     j=strlen(str);
     temp[i-j]='\0';     /* Truncate trailing AAAs */
     break; 
   }
   strcpy(temp2,"");
  
   if(sawT==0)
      sscanf(str,"%*[A]%s",temp2); 
   else
      sscanf(str,"%*[T]%s",temp2); 
      
   if(strlen(temp2)<=cutoff)
   {
     i=strlen(temp);
     j=strlen(str);
     temp[i-j]='\0';     /* Truncate trailing AAAs */
     break; 
   } 
   strcpy(temp1,temp2);
  } /* end while 2 */

  strcpy(frag,temp);
}


void toUpperFrag(char *s) {
	int i=0,l=0;
	l=strlen(s);
	for(i=0;i<l;i++)  if(s[i]>='a' && s[i]<='z')  s[i]=toupper(s[i]);
} /* end toUpperFrag */	
