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
#include <errno.h>
#include "dynamic.h"
#include "err.h"
#include "debug.h"
#include "est.h"

const int MNINT=(-1*(1<<(sizeof(int)*8-1)-1));

int  match;  /* match */
int gap;  /* continuous gap penalty */
int hgap;  /* start gap penalty */
int mism;  /* mismatch penalty */
float threshold; /* Optimal score variation for calculating band width */
int AlignmentWithN;  /* aligning a base with N */

int g_dpf1=-1,g_dpf2=-1;
int g_dpp1=-1,g_dpp2=-1;
int g_dplen1=0,g_dplen2=0;
char g_dpdir='x';
char *g_dps1=NULL,*g_dps2=NULL;

char msg[200];

extern int errno;

/* Main function to align a pair of sequences */

void alignPair(struct dynResult *result,int f1,int f2,
				char *s1,char *s2,int p1,int p2,
				int len1, int len2,char type,int d,
				struct dynParams param) {
	int i,j;
  	int m,n,wp2;
  	char m1[MAXFRAGSIZE],n1[MAXFRAGSIZE],m2[MAXFRAGSIZE],n2[MAXFRAGSIZE];
  	struct bandResult area1,area2;
 
  	int midscore;
  	char s2work[MAXFRAGSIZE];
  	int calculateMiddleScore(char *a,char *b,int p1,int p2,int len);

  	/* store which pair is being aligned - that way pair info is 
	*available in case there is an error while aligning.
  	*/
  
  	g_dpf1=f1;
  	g_dpf2=f2;
  	g_dps1=s1;
  	g_dps2=s2;
  	g_dpp1=p1;
  	g_dpp2=p2;
  	g_dplen1=len1;
  	g_dplen2=len2;
  	g_dpdir=type;

  	sprintf(msg,"Pair: %d %d : pos %d %d : len %d %d : dir %c ",
		g_dpf1,g_dpf2,g_dpp1,g_dpp2,g_dplen1,g_dplen2,g_dpdir);	

  
  	match = param.match;
  	gap = param.gap;
  	hgap = param.hgap;
  	mism = param.mism;
  	AlignmentWithN = param.AlignmentWithN;
  	threshold = param.threshold;
 
  	if(match==MNINT || mism==MNINT || gap==MNINT || 
			hgap==MNINT || AlignmentWithN==MNINT || threshold>1.0) {
    		printf("Error: %s file missing or incomplete.\n",CFGFile);
    		return; 
  	}

	/* convention: always have s1 in direct form, and s2 reverse complemented if necessary */
 
  	if(type=='D') {
     	strcpy(s2work,s2);
     	wp2=p2; 
  	}
  	else if(type=='P') {
     	getRComp(s2,s2work);
     	wp2 = strlen(s2work)-p2-len2;
  	}
  	else {
     	printf("Rank=%d : Error : type<%c> has to be either P or D!\n",rank,type);
     	printf("Rank=%d : Error : p1=%d p2=%d l1=%d l2=%d!\n",rank,p1,p2,len1,len2);
     	return;
  	}

  	m = strlen(s1);
  	n = strlen(s2work);

  	reverse(s1,0,p1-1,m1);
  	reverse(s2work,0,wp2-1,n1);

  	/* printf("s1=\n%s\n,s2work=\n%s\n",s1,s2work); */

  	for(i=p1+len1,j=0;i<m;i++,j++)  m2[j]=s1[i];
  	m2[j]='\0';

  	for(i=wp2+len2,j=0;i<n;i++,j++)  n2[j]=s2work[i];
  	n2[j]='\0';

	area1 = bandAlign(m1,n1);
   	area2 = bandAlign(m2,n2); 
 

    	/* Ananth - Change for N,X in maximal match - 07/05/2003 */
    	/* Scan middle portion of maximal match for N or X and treat them as mismatch*/

    	if(len1!=len2) {  /* right now len1 will always equal len2*/
		printf("Rank=%d: Warning: len1!=len2: %d!=%d for pair %s\n",rank,len1,len2,msg);
		fflush(stdout);
    	}

    	/* midscore = calculateScore(len1,len2,d);*/
    	midscore = calculateMiddleScore(s1,s2work,p1,wp2,len1);
    	/* Ananth - End Change for N,X in maximal match - 07/05/2003 */


 	/* printf("midscore=%d\n",midscore);  */
  	buildResult(result,area1,area2,midscore,
				len1,len2,strlen(m1),strlen(n1),
                    strlen(m2),strlen(n2));
}   /* end alignPair */


/* bandAlign : Uses Banded Dynamic programming , with affine gap
               penalty and applies semiglobal alignment for ignoring
               end gaps.
*/

struct bandResult bandAlign(char *m,char *n) {
	int k; /* band width */ 
    	int l1,l2;
    	char xaxis = '\0',yaxis='\0';
    	char s_x[MAXFRAGSIZE],s_y[MAXFRAGSIZE],c_x,c_y;
    	struct rowptr a,b,c; // current row
    	struct rowptr prow_a,prow_b,prow_c; // previous row
    	int i,j,p;
    	int left[2],up[2],cross[2],curr[2];
    	struct bandResult bandres;
    	int maxscore,maxborder=0,tempmax=0;
    	int maxpos[2],maxpos1[2],maxborderpos[2],maxborderpos1[2];
    	char dum,maxin;
    	struct pth *path=NULL;
    	char stemp[2000];
    	int max_row_size;
    	int maxborderscore;

    	bandres.borderScore = MNINT;
    	bandres.mORn = 'e';
    	bandres.touches = 0;
    	bandres.maxScore = MNINT;
    	bandres.maxScoreLen = 0;

    	if(strlen(m)==0) {
     	bandres.borderScore=0;
       	bandres.maxScore=0;
       	bandres.mORn='n';
       	return bandres;
    	}

    	if(strlen(n)==0) {
       	bandres.borderScore=0;
       	bandres.maxScore=0;
       	bandres.mORn='m';
       	return bandres;
    	}
 
    	/* convention: l1<=l2 enforced */
    	if(strlen(m)<=strlen(n)) {
      	l1 = strlen(m);
      	l2 = strlen(n);
      	xaxis = 'm';
      	strcpy(s_x,m);
      	yaxis = 'n';
      	strcpy(s_y,n);
    	}
    	else {
      	l1 = strlen(n);
      	l2 = strlen(m);
      	xaxis = 'n';
      	strcpy(s_x,n);
      	yaxis = 'm';
      	strcpy(s_y,m);
    	}

 
    	/* calculate k  - slightly overestimate*/
    	/* k = l1 * (match*(1-threshold)/(match-gap)) + 1; */
    	k = l1 * Kfactor + 1;

	// maximum size of each row is 2k+1 
	// not all cells may be used though

	max_row_size = 2*k+1;
	a.ptr = (struct cell *)malloc(sizeof(struct cell)*max_row_size);
	b.ptr = (struct cell *)malloc(sizeof(struct cell)*max_row_size);
	c.ptr = (struct cell *)malloc(sizeof(struct cell)*max_row_size);
	prow_a.ptr = (struct cell *)malloc(sizeof(struct cell)*max_row_size);
	prow_b.ptr = (struct cell *)malloc(sizeof(struct cell)*max_row_size);
	prow_c.ptr = (struct cell *)malloc(sizeof(struct cell)*max_row_size);

	if(!a.ptr || !b.ptr || !c.ptr || !prow_a.ptr || !prow_b.ptr || !prow_c.ptr ) {
		sprintf(stemp,"bandAlign: %s (%s) Banded DP rows = %d bytes",msg,strerror(errno),6*sizeof(struct cell)*max_row_size);
		AlertErr(stemp);
	}

	initBandRow(&a,0,l1,l2,k,'a');
	initBandRow(&b,0,l1,l2,k,'b');
	initBandRow(&c,0,l1,l2,k,'c');


    	maxborderscore = MNINT;
    	maxscore = MNINT;

    	/* Do banded dynamic programming using affine gaps*/
    	for(i=1;i<=l1;i++) {
    		// copy row to prev_row

		copyRows(&prow_a,&a);  // prow_a = a & prow_size = row_size
		copyRows(&prow_b,&b);  // prow_a = a & prow_size = row_size
		copyRows(&prow_c,&c);  // prow_a = a & prow_size = row_size

		initBandRow(&a,i,l1,l2,k,'a');
		initBandRow(&b,i,l1,l2,k,'b');
		initBandRow(&c,i,l1,l2,k,'c');

       	for(j=0;j<a.n;j++) { // for all cells in the current (i^th) row
         		if(i<=k && j==0) continue;   
         		if(i==0 && j<=k) continue;   
           	/* because values has already been initialized */

         		bandToArr(i,j,k,curr);
         		if(curr[0]==l1 || curr[1]==l2) maxborder=1;
         		else maxborder=0;
         
         		/* now (curr[0],curr[1]) are the original co-ords */

         		/* left,up,cross are band co-ords */     
         		arrToBand(curr[0],curr[1]-1,k,left);
         		arrToBand(curr[0]-1,curr[1],k,up);
			arrToBand(curr[0]-1,curr[1]-1,k,cross); 
 
         		if( (cross[0]==-1) || (cross[1]==-1)) {
            		printf("Rank %d: Error : bandAlign : cross[i-1,j-1] out of bound for pair %s!\n",rank,msg);
            		return 	bandres;
         		}

         		/* compute arrays a,b,c for i,j */ 
         		c_x= s_x[curr[0]-1];
         		c_y= s_y[curr[1]-1];

         		p=findMatch(c_x,c_y,match,mism);

         		a.ptr[j].score = p + maxi(prow_a.ptr[cross[1]].score, 
                                      	prow_b.ptr[cross[1]].score,
                                      	prow_c.ptr[cross[1]].score
							   );


         		if( (left[0]!=-1) && (left[1]!=-1)) { 
             		b.ptr[j].score = maxi(addscore(hgap+gap,a.ptr[left[1]].score),
                                      addscore(gap,b.ptr[left[1]].score),
                                      addscore(hgap+gap,c.ptr[left[1]].score)
							   );
         		}
         
         		if( (up[0]!=-1) && (up[1]!=-1)) { 
             		c.ptr[j].score = maxi(addscore(hgap+gap,prow_a.ptr[up[1]].score),
                                      addscore(hgap+gap,prow_b.ptr[up[1]].score), 
                                      addscore(gap,prow_c.ptr[up[1]].score)
							   );
         		}

         		tempmax = maxi(a.ptr[j].score,b.ptr[j].score,c.ptr[j].score);

         		if(tempmax>maxscore) {
             		maxscore = tempmax;
             		maxpos[0]=i; 
             		maxpos[1]=j;
         		}

	 		if(maxborder) {
             		if(tempmax>=maxborderscore) {
             			maxborderscore = tempmax;
             			maxborderpos[0]=i; 
              			maxborderpos[1]=j;
             			maxin = dum;
	     		}
	 		} 

       	} /* end for j*/
	} /* end for i*/

	/* Assign the maximum score both inside and border of the band*/

    	bandres.borderScore=maxborderscore;
    	bandToArr(maxborderpos[0],maxborderpos[1],k,maxborderpos1);

    	bandres.maxScore=maxscore;
    	bandToArr(maxpos[0],maxpos[1],k,maxpos1);

    	if(maxborderpos1[0]==l1) {
      	bandres.mORn = yaxis;
      	bandres.touches = maxborderpos1[1];
    	}
    	else if(maxborderpos1[1]==l2) {
      	bandres.mORn = xaxis;
      	bandres.touches = maxborderpos1[0];
    	}


    	if(maxpos1[0]<maxpos1[1])  bandres.maxScoreLen=maxpos1[0];
    	else bandres.maxScoreLen=maxpos1[1];

 

     if(a.ptr!=NULL) free(a.ptr);
	else printf("Warning: cannot free a \n");

     if(b.ptr!=NULL) free(b.ptr);
	else printf("Warning: cannot free b \n");

     if(c.ptr!=NULL) free(c.ptr);
	else printf("Warning: cannot free c \n");

     if(prow_a.ptr!=NULL) free(prow_a.ptr);
	else printf("Warning: cannot free prow_a \n");

     if(prow_b.ptr!=NULL) free(prow_b.ptr);
	else printf("Warning: cannot free prow_b \n");

     if(prow_c.ptr!=NULL) free(prow_c.ptr);
	else printf("Warning: cannot free prow_c \n");

    	if(path!=NULL) free(path);
    	return bandres;
    
} /* end bandAlign */

/* Transformation functions : 2-D array <---> banded array */
void bandToArr(int i,int j,int k,int a[])
{

   if(i<=k)  a[1]=j;
   else      a[1]=j+i-k;
   a[0]=i; 
} /* end bandToArr */

void arrToBand(int i,int j,int k,int a[])
{
   
   if(!insideStrip(i,j,k)) 
   {
     a[0]=a[1]=-1;
     return;
   }

   if(i<=k)  a[1]=j;
   else      a[1]=j-(i-k);
   a[0]=i; 
} /* end arrToBand */

/* input i,j in original array co-ords - not that of the band */
short insideStrip(int i,int j,int k)
{
   if(abs(i-j)<=k) return 1;
   else return 0;
}  /* end insideStrip */

/* Band manipulation functions */

/* initBandRow : Initializes a row with the affine gap penalties
                 Assumes that l1<=l2
                 l1: #rows
                 l2: #cols
	eg., initBandRow(&a,0,l1,l2,k,'a');

	Returns: number of columns/cells initialized for this row row_i
*/
int initBandRow(struct rowptr *row,int row_i,int l1,int l2,int k,char type) {
    	int i,j,t; 
    	char stemp[1000];
    
    	if(l1>l2 || k>l1)  {
        printf("Rank=%d: Error: initBandRow : l1>l2 or k>l1 ! for pair %s\n",rank,msg);
        return 0;
    	}  

    	i = row_i; 
    	if( (i+k)<=l2)  j = k+1;
    	else j = l2-i+1;

    	if(i>=k) j+=k;
    	else	j+=i;

	if(j<=0) {
	  	printf("Rank=%d: Warning : Table row %c[%d] is empty for %s\n",rank,type,row_i,msg);
		row->n = 0;
		return 0;
	}

     row->n = j;

     for(t=0;t<j;t++) {
         switch(type) {
           case 'a' :
                      if(i==0) {
                		if(t==0)  row->ptr[t].score=0;
                        	else  row->ptr[t].score = MNINT;
                      }
                      else 
                      if(i<=k && t==0) row->ptr[t].score = MNINT;
                      else  row->ptr[t].score = 0;
                      break;
                         
           case 'b' :
                      if(i==0 && t!=0) {
                        row->ptr[t].score = hgap + gap*t;
                      }
                      else 
                      if(t==0) row->ptr[t].score = MNINT; 
                      else  row->ptr[t].score = 0;
                      break;

           case 'c' :
                      if((i==0) || t==(2*k) || (t==(j-1) && i<=k)) {
                        row->ptr[t].score = MNINT;
                      }
                      else 
                      if(i<=k && t==0) row->ptr[t].score = hgap + gap*i; 
                      else  row->ptr[t].score = 0;
                      break;

           default  :  
                      printf("Rank=%d: Error : initBandRow : 'type' should be one of these : a,b,c for pair %s\n",rank,msg);
				  row->n=0;
  		      	  return 0;
 

         } /* end switch */
	} /* end for */

    return row->n;
} /* end initBandRow */



int calculateScore(int len1,int len2,int d) {

  int big,diff;

  if(len1==len2)
  {
    if(d==0) return len1*match;
    if(d<0) return (-(d*mism)+(len1+d)*match); 
    return (d*gap + (len1-d)*match);
  }
 
  big=len1;
  diff = len1-len2;   

  if(len2>big)
  {
    big=len2;
    diff = len2-len1;
  }

  if(d>=0)  
     return (diff*gap + (d-diff)*mism + (big-d)*match);
 
   printf("Warning : Input has d (-)ve for -e option!\n"); 
   return (-(d*mism)+(big+d)*match); 

}/* end calculateScore */



int calculateMiddleScore(char *a,char *b,int p1,int p2,int len) {
  int i,j,l1,l2;
  int mcount=0,mismcount=0;
  void printSubString(char *a,int pos,int len);
  int bNorX(char);
  
  l1=strlen(a);  l2=strlen(b);
  if((p1+len-1)>=l1 || (p2+len-1)>=l2) {
     printf("Rank=%d: Warning: calculateMiddleScore: out of bounds for %s!\n",rank,msg); 
     fflush(stdout);
     return 0;
  }
  
  i=p1; j=p2;
  while( i<(p1+len) ) {
     if(a[i]==b[j]) {
          if(a[i]!='A' && a[i]!='C' && a[i]!='G' && a[i]!='T' ) mismcount++;
          else mcount++;
     }
     else {
          if(rank==-1 && !bNorX(a[i]) &&  !bNorX(b[j]) ) {
                   printf("Rank=%d: Warning: calculateMiddleScore: a[%d] b[%d] %c, %c!\n",rank,i,j,a[i],b[j]); fflush(stdout);
                   printf("Rank=%d: a: %s\n %d to %d\n",rank,a,p1,p1+len-1);
                   printf("Rank=%d: b: %s\n %d to %d\n",rank,b,p2,p2+len-1);
                   printf("Rank=%d:\n",rank);
                   printSubString(a,p1,len);
                   printf("Rank=%d:\n",rank);
                   printSubString(b,p2,len);
          }
          mismcount++;
     }
     i++;
     j++;
   } /* end while */
   return (mcount*match + mismcount*mism);
} /* end calculateMiddleScore */



/* This function pulls together the two alignment scores computed
*	on either side of the maximal match
*	Total score = area1.score + score for maximal match + area2.score
*/

void buildResult(struct dynResult *final,
				struct bandResult area1,struct bandResult area2,
                    int midscore,int l1,int l2,
				int m1,int n1,int m2, int n2 ) {
	int mlen,nlen,min;


   	final->maxBorderScore = area1.borderScore + midscore + area2.borderScore;
   	final->maxScore = area1.maxScore + midscore + area2.maxScore;

   	/* detect type of overlap */
   	if(area1.mORn==area2.mORn) {
     	final->typeAlign = 'C';
     	if(area1.mORn=='m') final->atag='n';  /* n is possibly contained in m */
     	else final->atag='m';                 /* m is possibly contained in n */
     	/* printf("We have a possible containment in %c! \n",area1.mORn); */
    
     	if(area1.mORn=='m') {
       		mlen = area1.touches+l1+area2.touches;
       		nlen = n1 + l2 + n2;

     	}
     	else {
       		nlen = area1.touches+l2+area2.touches;
       		mlen = m1 + l1 + m2;
     	}  

    	 	min=mlen;
     	if(nlen<min) min=nlen;

     	final->idealBorderScore = min*match;
   	}
   	else {
     	final->typeAlign = 'S';
     	/* printf("We have a prefix-suffix overlap! \n"); */
     	if(area1.mORn=='m') {
       		final->atag='m';               /* m is left of n */
       		mlen=area1.touches + l1 + m2;
       		nlen=n1 + l2 + area2.touches;
     	}    
     	else {
       		final->atag='n';               /* n is left of m */
       		nlen=area1.touches + l2 + n2;
       		mlen=m1 + l1 + area2.touches;
     	}    
   
     	min=mlen;
     	if(nlen<min) min=nlen;
     	final->idealBorderScore = min*match;
	}

   	final->idealScore = (area1.maxScoreLen+l1+area2.maxScoreLen)*match;
   	/* compute the coverage of the maxScore Alignment with the smaller string */
   	min=m1+l1+m2;
   	if(min>(n1+l2+n2)) min=n1+l2+n2;
   	final->coverage = ( (area1.maxScoreLen+l1+area2.maxScoreLen)*100)/min;

} /* end buildResult */



/* Returns the score after trying to match x and y input chars */
/* m - is match score  , mism - is mismatch score */

int findMatch(char x,char y,int m,int mism) {
   short forx=0,fory=0;
   char a,b;


    /* Ananth Change: Swapped 1st and 3rd check and included X check middle:  16/6/2003 */
    if(x=='N' || y=='N') return AlignmentWithN;
    if(x=='X' || y=='X') return AlignmentWithN;
    if(x==y) return m;
    /* Change over*/

   if(x=='A' || x=='G' || x=='T' || x=='C' || x=='U') forx=1;
   if(y=='A' || y=='G' || y=='T' || y=='C' || y=='U') fory=1;

   if(forx && fory) return mism;

   if(forx) 
   {
     a=y; b=x; 
   }
   else
   {
     a=x; b=y; 
   }

   /*if(rank==0) { printf("##<%c,%c>##",a,b); fflush(stdout); }*/

   switch(a)
   {
      case 'R' :  if(b=='G' || b=='A') return m;
                  return mism;
      case 'Y' :  if(b=='T' || b=='C') return m;
                  return mism;
      case 'K' :  if(b=='G' || b=='T') return m;
                  return mism;
      case 'M' :  if(b=='A' || b=='C') return m;
                  return mism;
      case 'S' :  if(b=='G' || b=='C') return m;
                  return mism;
      case 'W' :  if(b=='A' || b=='T') return m;
                  return mism;
      case 'B' :  if(b=='G' || b=='T'|| b=='C') return m;
                  return mism;
      case 'D' :  if(b=='G' || b=='T'|| b=='A') return m;
                  return mism;
      case 'H' :  if(b=='C' || b=='T'|| b=='A') return m;
                  return mism;
      case 'V' :  if(b=='C' || b=='G'|| b=='A') return m;
                  return mism;
      default  :  
		  printf("Rank=%d: Error : Unrecognized character<%c> with <%c> found in sequence for pair: %s\n",rank,a,b,msg);
                  return mism; 
   } /* end switch  */

} /* end findMatch */



/************** SUPPORT FUNCTIONS **********************/

void copyRows(struct rowptr *r1, struct rowptr *r2) {
	int i=0;

	for(i=0;i<r2->n;i++) r1->ptr[i].score = r2->ptr[i].score;
	r1->n = r2->n;
} // end copyRows

int addscore(int i,int score) {
  if(score<=MNINT) return MNINT;
  return i+score;
}

/* Reverse a string  */
void reverse(char *src,int start,int end,char *dst) {
  int i=0,j=0;

  for(i=end,j=0;i>=start;i--,j++)
  {
     dst[j]=src[i]; 
  }
  dst[j]='\0';

} /* end reverse */

int maxi(int r,int s,int t) {
  int max;
  max=r;
  if(s>r)   max=s;
  if(t>max) max=t;
  return max;
} /* end maxi */


void dispParam(struct dynParams p)
{
   printf("Dynamic programming parameters: \n");
   printf("Match=%d, MisMatch=%d, Gap opening=%d, Gap continuation =%d, Aligning with N =%d, Band width threshold=%f\n",p.match,p.mism,p.hgap,p.gap,p.AlignmentWithN,p.threshold);

   fflush(stdout);
}/* end dispParam */ 


void printSubString(char *a,int pos,int len) {
  int i=0,l;
  l=strlen(a);
  if(pos+len>l) return;
  for(i=pos;i<(pos+len);i++)  printf("%c",a[i]);
  printf("\n");
} /* end printSubString */
