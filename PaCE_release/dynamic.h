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
  * Copyright � 2003 Iowa State University Research Foundation, Inc.
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


#ifndef __DYMANIC_H__
#define __DYMANIC_H__

extern int  match;  /* match */
extern int gap;  /* continuous gap penalty */
extern int hgap;  /* start gap penalty */
extern int mism;  /* mismatch penalty */
extern float threshold; /* Optimal score variation for calculating band width */
extern int AlignmentWithN;  /* score/penalty for aligning a base with N */
extern int bandsize;  /* # cells in the band */

extern const int MNINT;

extern int p,rank;
extern float Kfactor;
extern char CFGFile[];


struct cell{
  int score;
  #ifdef AlignmentRetrace
  char from;  
  #endif
};


struct rowptr {
   struct cell *ptr;
   int n;
};

struct dynResult{
  int maxBorderScore;
  int idealBorderScore;
  char typeAlign;  /* C - containment, S - suffix-prefix */ 
  char atag;
  int  maxScore;
  int idealScore;
  int coverage; /* coverage of the max score alignlen over shorter string in percentage */
};

struct bandResult{
  int borderScore; /* best score on the border */
  char mORn;  /* m - if band align ended in m's side
                  else if on n's side then 'n', or if error then 'e' */ 
  int touches; /* length from the origin of the maxBorderScore */
  int maxScore;  /* could be anywhere in the band*/
  int maxScoreLen; /* minimum of the number of characters in the alignment till the maxScore */
};

struct pth {
  int in[2];
};


struct dynParams {
  int match;  /* match */
  int gap;  /* continuous gap penalty */
  int hgap;  /* start gap penalty */
  int mism;  /* mismatch penalty */
  int AlignmentWithN;  /* score for aligning a base with N*/
  float threshold; /* Optimal score variation for calculating band width */
};



void alignPair(struct dynResult *,int f1,int f2,char *s1,char *s2,int p1,int p2,int len1, int len2,char type,int d,struct dynParams); 
void reverse(char *src,int start,int end,char *dst);
struct bandResult bandAlign(char *,char *);
void bandToArr(int i,int j,int k,int []);
void arrToBand(int i,int j,int k,int []);
short insideStrip(int i,int j,int k);
#ifdef AlignmentRetrace
int maxi(int r,int s,int t,char *which);
#else
int maxi(int r,int s,int t);
#endif

int initBandRow(struct rowptr *row,int row_i,int l1,int l2,int k,char type);
void copyRows(struct rowptr *r1, struct rowptr *r2);


int addscore(int i,int score);
#ifdef AlignmentRetrace
int retrace(int *maxpos,char maxin,struct rowptr *a,struct rowptr *b,
                           struct rowptr *c,int k,struct pth path[],int p);
#endif
void dispPath(struct pth *path,int len);
void initPath(struct pth *path,int len);
void dispAlign(char *s1,char *s2,struct pth *path,int len);
void dispString(char *s,struct pth *path,int len,short ix);
int calculateScore(int len1,int len2,int Eval);
void buildResult(struct dynResult *,struct bandResult area1,struct bandResult area2,
                      int midscore,int l1,int l2,int m1,int n1,int m2, int n2);
void dispParam(struct dynParams);
int findMatch(char x,char y,int m,int mism);
#endif
