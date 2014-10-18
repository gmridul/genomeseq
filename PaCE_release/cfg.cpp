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



#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "cfg.h"
#include "err.h"
#include "debug.h"
#include "dynamic.h"


extern int errno;

struct dynParams getDynParams()
{
   struct dynParams param;
   FILE *fp;
   char word[50];

   param.match=MNINT;  /* match */
   param.gap=MNINT;  /* continuous gap penalty */
   param.hgap=MNINT;  /* start gap penalty */
   param.mism=MNINT;  /* mismatch penalty */
   param.AlignmentWithN=MNINT;  
   param.threshold=0.5; /* Optimal score variation for calculating band width */

   fp = fopen(CFGFile,"r");
   while(!feof(fp))
   {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"match")) 
    {
       fscanf(fp,"%s",word);
       param.match = atoi(word);
    }

    if(!strcmp(word,"mismatch")) 
    {
       fscanf(fp,"%s",word);
       param.mism = atoi(word);
    }
    
    if(!strcmp(word,"gap")) 
    {
       fscanf(fp,"%s",word);
       param.gap = atoi(word);
    }

    if(!strcmp(word,"hgap")) 
    {
       fscanf(fp,"%s",word);
       param.hgap = atoi(word);
    }

    if(!strcmp(word,"EndToEndScoreRatioThreshold")) 
    {
       fscanf(fp,"%s",word);
       param.threshold = (100-atof(word))/100;
    }

    if(!strcmp(word,"AlignmentWithN")) 
    {
       fscanf(fp,"%s",word);
       param.AlignmentWithN= atoi(word);
    }

  } /* end while */
  
  fclose(fp);

  if(param.match==MNINT || param.mism==MNINT || param.gap==MNINT || param.hgap==MNINT || param.AlignmentWithN==MNINT ||
           param.threshold==(2.0))
  {
    printf(" Error: %s file missing or incomplete.\n",CFGFile);
    terminate();
    exit(1);
  }
 
  /* printf("Rank=%d: Load-Parameters : %d %d %d %d %f\n\n",rank,param.match,param.mism,param.gap,param.hgap,param.threshold);  */

  if(rank==0) dispParam(param);
  return param;

} /* end getDynParams */


int getThreshold(char *whatThreshold) {
   FILE *fp;
   char word[50];
   
   fp = fopen(CFGFile,"r");
   while(!feof(fp))
   {
    fscanf(fp,"%s",word);
    if(!strcmp(word,whatThreshold)) 
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
  printf("Error: %s does not have %s.\n",CFGFile,whatThreshold);
  terminate();
  exit(1);
} /* end getThreshold */

int getTranscriptsTogether()
{
   FILE *fp;
   char word[50];
   
   fp = fopen(CFGFile,"r");
   while(!feof(fp))
   {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"TranscriptsTogether")) 
    {
       fscanf(fp,"%s",word);
       if(atoi(word)<0 || atoi(word)>1) 
       {
          printf("Error: %s TranscriptTogether should be 0 or 1.\n",CFGFile);
          terminate();
          exit(1);
       }
       return atoi(word);
    }
    else
    {
       continue;
    }

   } /* end while */
  
  fclose(fp);
  printf("Serious Warning: %s does not have TranscriptsTogether.\n",CFGFile);
  return 0;
} /* end getTranscriptsTogether */

int getDumpClustersMidway()
{
   FILE *fp;
   char word[50];
   
   fp = fopen(CFGFile,"r");
   while(!feof(fp))
   {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"DumpClustersMidway")) 
    {
       fscanf(fp,"%s",word);
       if(atoi(word)<0 || atoi(word)>1) 
       {
          printf("Error: %s DumpClustersMidway should be 0 or 1.\n",CFGFile);
          terminate();
          exit(1);
       }
       return atoi(word);
    }
    else
    {
       continue;
    }

   } /* end while */
  
  fclose(fp);
  printf("Warning: %s does not have DumpClustersMidway.\n",CFGFile);
  return 0;
} /* end getDumpClustersMidway */

int getReportSplicedCandidates()
{
   FILE *fp;
   char word[50];
   
   fp = fopen(CFGFile,"r");
   while(!feof(fp))
   {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"ReportSplicedCandidates")) 
    {
       fscanf(fp,"%s",word);
       if(atoi(word)<0 || atoi(word)>1) 
       {
          printf("Error: %s ReportSplicedCandidates should be 0 or 1.\n",CFGFile);
          terminate();
          exit(1);
       }
       return atoi(word);
    }
    else
    {
       continue;
    }
   } /* end while */
  fclose(fp);
  printf("Warning: %s does not have ReportSplicedCandidates.\n",CFGFile);
  return 0;
} /* end getReportSplicedCandidates */

int getReportMaximalPairs()
{
   FILE *fp;
   char word[50];
   
   fp = fopen(CFGFile,"r");
   while(!feof(fp))
   {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"ReportMaximalPairs")) 
    {
       fscanf(fp,"%s",word);
       if(atoi(word)<0 || atoi(word)>1) 
       {
          printf("Error: %s ReportMaximalPairs should be 0 or 1.\n",CFGFile);
          terminate();
          exit(1);
       }
       return atoi(word);
    }
    else
    {
       continue;
    }
   } /* end while */
  fclose(fp);
  printf("Warning: %s does not have ReportMaximalPairs.\n",CFGFile);
  return 0;
} /* end getReportMaximalPairs*/

int getReportMaximalSubstrings()
{
   FILE *fp;
   char word[50];
   
   fp = fopen(CFGFile,"r");
   while(!feof(fp))
   {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"ReportMaximalSubstrings")) 
    {
       fscanf(fp,"%s",word);
       if(atoi(word)<0 || atoi(word)>1) 
       {
          printf("Error: %s ReportMaximalSubstrings should be 0 or 1.\n",CFGFile);
          terminate();
          exit(1);
       }
       return atoi(word);
    }
    else
    {
       continue;
    }
   } /* end while */
  fclose(fp);
  printf("Warning: %s does not have ReportMaximalSubstrings.\n",CFGFile);
  return 0;
} /* end getReportMaximalSubstrings*/

int getReportGeneratedPairs()
{
   FILE *fp;
   char word[50];
   
   fp = fopen(CFGFile,"r");
   while(!feof(fp))
   {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"ReportGeneratedPairs")) 
    {
       fscanf(fp,"%s",word);
       if(atoi(word)<0 || atoi(word)>1) 
       {
          printf("Error: %s ReportGeneratedPairs should be 0 or 1.\n",CFGFile);
          terminate();
          exit(1);
       }
      if(rank==0  && (atoi(word)==1)) 
		printf("Rank=%d: INFO: ReportGeneratedPairs SET TO 1\n",rank);

       return atoi(word);
    }
    else
    {
       continue;
    }
   } /* end while */
  fclose(fp);
  printf("Warning: %s does not have ReportGeneratedPairs.\n",CFGFile);
  return 0;
} /* end getReportGeneratedPairs*/


int getReportAcceptedPairs()
{
   FILE *fp;
   char word[50];
   
   fp = fopen(CFGFile,"r");
   while(!feof(fp))
   {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"ReportAcceptedPairs")) 
    {
       fscanf(fp,"%s",word);
       if(atoi(word)<0 || atoi(word)>1) 
       {
          printf("Error: %s ReportAcceptedPairs should be 0 or 1.\n",CFGFile);
          terminate();
          exit(1);
       }
       return atoi(word);
    }
    else
    {
       continue;
    }
   } /* end while */
  fclose(fp);
  printf("Error: %s does not have ReportAcceptedPairs.\n",CFGFile);
  terminate();
  exit(1);
} /* end getReportAcceptedPairs*/



void getCloneMatesFile(char *cf)
{
   FILE *fp;
   char word[1000];
   
   fp = fopen(CFGFile,"r");
   while(!feof(fp))
   {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"ClonePairsFile")) 
    {
       fscanf(fp,"%s",word);
       if(!strcmp(word,"None")) 
       {
	       strcpy(cf,"");
	       printf("Master: Info: No clone pairs information supplied\n");
       }
       else 
       {
	       strcpy(cf,word);
	       printf("Master: Info: Clone pairs information supplied is %s\n",word);
       }
       return;
    }
    else
    {
       continue;
    }

   } /* end while */
  
  fclose(fp);
  printf("Warning: %s does not have clonePairsFile field.\n",CFGFile);
  fflush(stdout);
  strcpy(cf,"");
  return;
} /* end getCloneMatesFile*/



void getOutputFolder(char *opfolder)
{
   FILE *fp;
   char word[1000];
   
   fp = fopen(CFGFile,"r");
   while(!feof(fp))
   {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"OutputFolder")) 
    {
       fscanf(fp,"%s",word);
       if(!strcmp(word,"None")) 
       {
	       strcpy(opfolder,".");
	       printf("Info: Output folder is current working directory\n");
       }
       else 
       {
	       strcpy(opfolder,word);
	       /*printf("Info: Output Folder is %s\n",word);*/
       }
       return;
    }
    else continue;

   } /* end while */
  
  fclose(fp);
  printf("Warning: %s does not have OutputFolder field.\n",CFGFile);
  fflush(stdout);
  strcpy(opfolder,"");
  return;
} /* end getOutputFolder */



int getKeep_Mbuf_Full() {
   FILE *fp;
   char word[50];
   
   fp = fopen(CFGFile,"r");
   while(!feof(fp)) {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"Keep_Mbuf_Full")) {
       fscanf(fp,"%s",word);
       if(atoi(word)<0) {
          printf("Error: %s Keep_Mbuf_Full should be >= 0.\n",CFGFile);
          terminate();
          exit(1);
       }
       return atoi(word);
    }
    else continue;
   } /* end while */
  fclose(fp);
  printf("Warning: %s does not have Keep_Mbuf_Full.\n",CFGFile);
  return 0; /* default  : means do not override st->E */
} /* end getKeepMbufFull */

/* Parameterize using MPI_Block_Sends in PaCE.cfg */
/*  eg.,
	MPI_Block_Sends 0	=> Can Use MPI_Isend wherever possible
	MPI_Block_Sends 1	=> Use MPI_Send in all places
*/

int getMPI_Block_Sends() {
   FILE *fp;
   char word[50];
   
   fp = fopen(CFGFile,"r");
   while(!feof(fp)) {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"MPI_Block_Sends")) {
       fscanf(fp,"%s",word);
       if(atoi(word)<0 || atoi(word)>1) {
          printf("Warning: %s  MPI_Block_Sends defaulting to 0.\n",CFGFile);
		return 0;
       }
       return atoi(word);
    }
    else continue;
   } /* end while */
  fclose(fp);
  printf("Warning: %s does not have MPI_Block_Sends.\n",CFGFile);
  return 0; /* default  */
} /* end get MPI_Block_Sends*/


int getOutputLargeMerges() {
   FILE *fp;
   char word[50];

   fp = fopen(CFGFile,"r");
   while(!feof(fp))
   {
    fscanf(fp,"%s",word);
    if(!strcmp(word,"OutputLargeMerges"))
    {
       fscanf(fp,"%s",word);
       if(atoi(word)<0 || atoi(word)>1)
       {
          printf("Error: %s OutputLargeMerges should be 0 or 1.\n",CFGFile);
          terminate();
          exit(1);
       }
       return atoi(word);
    }
    else
    {
       continue;
    }

   } /* end while */

  fclose(fp);
  printf("Serious Warning: %s does not have OutputLargeMerges.\n",CFGFile);
  return 0;
} /* end getOutputLargeMerges*/

