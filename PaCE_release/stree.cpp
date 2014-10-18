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
#include <math.h>
#include <assert.h>
#include "keys.h"
#include "stree.h"
#include "est.h"
#include "err.h"
#include "debug.h"

/*extern struct fkp *final_fragkeys;*/
extern struct est *AllESTs;
struct stnode *Nodes;
int *depnodes;
int *g_NodeMarker;
unsigned int g_iNodes=0;
unsigned int g_currMaxNodes=0;
unsigned int g_iLeaves=0;
unsigned int g_iGenPairs=0;
unsigned int g_iNumCharsVisited=0;
unsigned int g_iAvgLeafDepth=0;
int minDepth=10000,maxDepth=0;


char charToSig(char d)
{
   char c;

   switch(d)
   {
        case 'A' :  c='A'; break;
        case 'C' :  c='C'; break;
        case 'G' :  c='G'; break;
        case 'T' :  c='T'; break;
        case 'X' :      /*  X treated as N and will be treated as $ for suffix tree- Ananth:  Nchange  */
        case 'N' :  c='A'; break;
        case 'R' :  c='A'; break;
        case 'Y' :  c='T'; break;
        case 'K' :  c='G'; break;
        case 'M' :  c='C'; break;
        case 'S' :  c='G'; break;
        case 'W' :  c='A'; break;
        case 'B' :  c='G'; break;
        case 'V' :  c='C'; break;
        case 'D' :  c='G'; break;
        case 'H' :  c='A'; break;
        case 'U' :  c='T'; break;
        default  :  c='A';
   }
   return c;
} /* end charToSig */

char charToPSig(char d)
{
   char c;

   switch(d)
   {
        case 'A' :  c='T'; break;
        case 'C' :  c='G'; break;
        case 'G' :  c='C'; break;
        case 'T' :  c='A'; break;
        case 'X' :      /*  X treated as N and will be treated as $ for suffix tree- Ananth:  Nchange  */
        case 'N' :  c='T'; break;
        case 'R' :  c='T'; break;
        case 'Y' :  c='A'; break;
        case 'K' :  c='C'; break;
        case 'M' :  c='G'; break;
        case 'S' :  c='C'; break;
        case 'W' :  c='T'; break;
        case 'B' :  c='C'; break;
        case 'V' :  c='G'; break;
        case 'D' :  c='C'; break;
        case 'H' :  c='T'; break;
        case 'U' :  c='A'; break;
        default  :  c='T';
   }
   return c;
} /* end charToSig */



void buildBucket(int start,int end,int head,int *pbuck,struct stnode *broot,struct fkp *buckeyes)
{
   int i=0,j=0,t=0;
   struct fkp *key;
   int heads[Sigma],next[Sigma];
   int lftlen[Sigma],it[Sigma];
   int offset[Sigma];
   char c;
   short len=0,ind=0;
   short d=0;  /* tracks the depth for each edge going out of this broot */
   int id,comeon=0;
   int dep=0;
   struct stnode *leaf,*newnode;
   int onlyDollar;
   int position;
   int nlocalLists=0,lsetlen=0,nxtoffset=0;
   struct suff delimit[Sigma],dumsuf;
   struct suff *lsptr=NULL;
   int iamat=0;
   char getLeftChar(struct fkp);
   int bNorX(char c) ;
   int getlcharFromDelimit(int id);
   int goTillDifferent(int start,int end,int head,int *pbuck,int depSoFar,struct fkp *buckeyes);

   if(start>end) return;
   if(head== -1) return;

   /* if(rank==2) printf("g_iNodes=%d\n",g_iNodes);*/
  /* a.  Partition the initial bucket into A,C,G,T,$ lists */
   for(t=0;t<Sigma;t++) { heads[t]=-1; next[t]=-1; }
   i=head;
   while(i>=0)
   {
        j=start+i;
        key = &(buckeyes[j]);
        len=strlen(AllESTs[key->fragid].s);
        if(key->dir=='D')
        {
           /* if((key->pos+(broot->depth))<len)*/
	   iamat = key->pos+(broot->depth);
	   if(iamat<len && !bNorX(AllESTs[key->fragid].s[iamat]) )
              c=charToSig(AllESTs[key->fragid].s[iamat]);
           else
              c='$';
        }
        else
        {
           comeon=len-1-(key->pos)-(broot->depth);
           if( comeon>=0  && !bNorX(AllESTs[key->fragid].s[comeon]) )
              c=charToPSig(AllESTs[key->fragid].s[comeon]);
           else
              c='$';
        }

        ind = sigToInd(c);
        if(heads[ind]==-1)
        {
          next[ind]=i; 
          heads[ind]=next[ind];
        }
        else
        {
          pbuck[next[ind]]=i;
          next[ind]=i;
        }
        i=pbuck[i];
        g_iNumCharsVisited++;
   } /* end while */ 
   for(t=0;t<Sigma;t++) { if(heads[t]>=0) pbuck[next[t]]=-1; }

   /* if(rank==2 && numbuckets==1) dispBucket(start,end,heads,pbuck); */
   /* b. For each edge list a in Sigma :
              - go as far as you can differentiate for A,C,G,T lists .
              - create a new node or leaf where it differentiates or ends 
              - call buildbucket recursively for each of the new nodes 
              - for $ lists, merge the leafnode to the parent node broot.
   */

   for(i=0;i<Sigma;i++)
   {
      if(heads[i]<=-1) continue;   /* empty */
      for(t=0;t<Sigma;t++) lftlen[t]=0; 
      if(indToSig(i)=='$')
      {
	 for(t=heads[i];t>=0;t=pbuck[t])
         {
 	    /* for each string in the list determine the char left of pos */
            j=start+t;
            c = getLeftChar(buckeyes[j]);
            ind = sigToInd(c);
            lftlen[ind]++;    
            g_iNumCharsVisited++;
         } /* end for t */


         /* lset change : include code for array creation here */
         /* sum leftlen[] and create array for that much size + for delimiters*/
        
         nlocalLists=0;
         lsetlen=0;

         for(t=0;t<Sigma;t++)
         {
	    if(lftlen[t]>0)
            {
               lsetlen=lsetlen + lftlen[t] + 1;
               delimit[nlocalLists].fid=getDelimitId(indToLftSig(t));
               delimit[nlocalLists].pos=lftlen[t];
		 assert(delimit[nlocalLists].pos>=0);
               nlocalLists++;
	    }
 	 }
         
         needs = sizeof(struct suff)*lsetlen;
         broot->lset=createSA(lsetlen);
   	 if(broot->lset) spaceSoFar+=needs;
   	 else AlertErr("broot->lset in $");
         g_treeSpace+=needs;
         thisPhaseSpace+=needs;

         /* compute offset for start of lset of each character */

         lsptr=broot->lset;
         nxtoffset=nlocalLists; /* Starting index of the first lset list after the delimiters */
         for(t=0;t<Sigma;t++) offset[t]=-1;
         
         for(t=0;t<nlocalLists;t++)
         {
            offset[sigToInd(getlcharFromDelimit(delimit[t].fid))] = nxtoffset;
            nxtoffset+=delimit[t].pos;

            /* Also enter the delimiters into the lset array of the node */
            putSuff(lsptr,t,delimit[t]);
            /*copySuff(&(lsptr[t]),delimit[t]); */
         }


         for(t=0;t<Sigma;t++)
         {
            if(lftlen[t]<=0) continue;
            it[t]=0;  /* offset[]+it[] indicates where next string in lset[] for that char should go */
         }

         /* end change */

         for(t=heads[i];t>=0;t=pbuck[t])
         {
 	    /* for each string in the list determine the char left of pos */
            j=start+t;
            c = getLeftChar(buckeyes[j]);
            ind = sigToInd(c);
            id=2*(buckeyes[j].fragid);
            if(buckeyes[j].dir=='P') id++;
            position = buckeyes[j].pos;

            if(lsptr==NULL)
            {
		printf("Rank=%d: Error: curr is null when there is a left char\n",rank);
		return; 
            }

            /* left list change : offset[lchar]+it[lchar] should give the next index where 
                                  this suffix should go */

            dumsuf.fid=id;
            dumsuf.pos=position;

	     assert(dumsuf.fid>=0);
	     assert(dumsuf.pos>=0);
            /*  lsptr[offset[ind]+it[ind]].fid=id;
            lsptr[offset[ind]+it[ind]].pos=position; */
            putSuff(lsptr,offset[ind]+it[ind],dumsuf);
        
            it[ind]++;
         } /* end for t */

          /* if(rank==3) dispNode(broot);  */
         onlyDollar=1;
         for(t=0;t<Sigma;t++)
         {
            if( (indToSig(t)!='$') && (heads[t]>=0) ) onlyDollar=0;
         }
         if(onlyDollar==1)
         { /* broot is a leaf */
            g_iLeaves++;
            broot->rmost=g_iNodes-1;  /* point to itself */
         }
   
      } /* end if $ */
      else /* A,C,G,T lists */
      {
         /* printf("Rank=%d: B4 BuildBucket: heads[i]=%d,pbuck[]=%d\n",rank,heads[i],pbuck[heads[i]]); fflush(stdout); */
 	 if((pbuck[heads[i]])<0) /* only one frag - so build a leaf */
         {
            j=start+heads[i];
            key = &(buckeyes[j]);
            len=strlen(AllESTs[key->fragid].s);
    	    dep=broot->depth + (len-(key->pos)-(broot->depth)); 
                                 /* works for P and D */

            /* leaf = createNode(&Nodes,&g_currMaxNodes,g_iNodes);*/
            leaf = &(Nodes[g_iNodes]);
            buildNode(leaf,dep);
            g_iNodes++;
            g_iLeaves++;
            
            c = getLeftChar(buckeyes[j]);
            ind = sigToInd(c);
            id=2*(key->fragid);
            if(key->dir=='P') id++;
	    position = key->pos;

            needs = sizeof(struct suff)*2;
            leaf->lset=createSA(2);
   	    if(leaf->lset) spaceSoFar+=needs;
   	    else AlertErr("leaf->lset in other chars");
            g_treeSpace+=needs;
            thisPhaseSpace+=needs;
     
            lsptr = leaf->lset;

            dumsuf.fid = getDelimitId(c);
            dumsuf.pos =1;
	     assert(dumsuf.fid>=0);
	     assert(dumsuf.pos>=0);
            putSuff(lsptr,0,dumsuf);
          
            /*lsptr[0].fid = getDelimitId(c); 
	    lsptr[0].pos = 1; */

            dumsuf.fid = id;
            dumsuf.pos =position;
	     assert(dumsuf.fid>=0);
	     assert(dumsuf.pos>=0);
            putSuff(lsptr,1,dumsuf);

            /*
            if(rank==0 && (g_iNodes-1)==3608)
            {
              printf("**This 3608\n");
              dispSA(lsptr,0);
              dispSA(lsptr,1);

            }*/

            /* lsptr[1].fid = id; 
	    lsptr[1].pos = position; */

            leaf->rmost=g_iNodes-1;
            broot->rmost=leaf->rmost;
            g_iNumCharsVisited++;

         }
         else
         {
            d = goTillDifferent(start,end,heads[i],pbuck,(broot->depth)+1,buckeyes);   
 	    /* create a node with depth broot->depth+d+1 */  
            /* this could either be a node or a leaf depending on */
            /* where it split */             
           
            /*newnode = createNode(&Nodes,&g_currMaxNodes,g_iNodes);*/
            newnode = &(Nodes[g_iNodes]);
            buildNode(newnode,(broot->depth)+d+1);

            g_iNodes++;
            
            /* Recursively call buildBucket for this new node */
            buildBucket(start,end,heads[i],pbuck,newnode,buckeyes);
            broot->rmost=newnode->rmost;    
         }
      } /* end else A,C,G,T lists */
   } /* end for i */
} /* end buildBucket */


/* buildNode : Builds a node corresponding within a bucket */
void buildNode(struct stnode *node,int depth)
{
    node->depth=depth; 
    node->rmost=-10;
    node->lset=NULL; 
} /* end buildNode*/



/* createNode : creates a new node in the NodeArr array */
/*              If the NodeArr is full upto currMax, this function 
                re-allocates NodeArr to increase the max size of the array
                by #NodesIncreaseBy new node entries */

struct stnode *createNode(struct stnode **NodesArr,int *currMax,int newIndex)
{

  struct stnode *temp;

  if(newIndex>=*currMax)
  {
    needs = sizeof(struct stnode)*NodesIncreaseBy;
    temp = (struct stnode *)realloc((*NodesArr),sizeof(struct stnode)*(*currMax+NodesIncreaseBy));
    *NodesArr = temp;
    if(*NodesArr) thisPhaseSpace+=needs;
    else AlertErr("Nodes-realloc");

    *currMax = *currMax + NodesIncreaseBy;
    printf("Slave %d: currMax b4=%d, after=%d\n",rank,*currMax-NodesIncreaseBy,*currMax);
  }

  return &( (*NodesArr)[newIndex]);
} /* end createNode */



/* computeOffsetFromLset : Given an lset array, it computes the
                           offset (starting index) of each of the
                           character thats has lset elements in this lset.
 		           The function returns the number of distinct
                           lsets in this collection (one for each char). 
*/

int computeOffsetFromLset(struct suff *lset,int *offs)
{
  int i=0,nlsets=0;
  int nxtoffset=0;		
  int getlcharFromDelimit(int id);

  if(!offs) return 0;
  for(i=0;i<Sigma;i++) offs[i]=-1;

  if(!lset) return 0;
  for(i=0;getSuffFid(lset,i)>=2*N;i++) {} /* just count the number of lsets */
  nlsets=i;

  nxtoffset=nlsets;
  for(i=0;i<nlsets;i++)
  {
     offs[sigToInd(getlcharFromDelimit(getSuffFid(lset,i)))] = nxtoffset;
     nxtoffset+=getSuffPos(lset,i);
  }
 
  return nlsets;
} /* end computeOffsetFromLset */

/* computeBoundFromLset : Given an lset array, it computes the
                           offsets (starting and ending) of each of the
                           character thats has lset elements in this lset.
 		           The function returns the number of distinct
                           lsets in this collection (one for each char). 
*/

int computeBoundFromLset(struct suff *lset,int *start,int *end)
{
  int i=0,nlsets=0,indx=0;
  int nxtoffset=0;
  int getlcharFromDelimit(int id);

  if(!start || !end) return 0;
  for(i=0;i<Sigma;i++) { start[i]=end[i]=-1; }

  if(!lset) return 0;
  for(i=0;getSuffFid(lset,i)>=2*N;i++) {} /* just count the number of lsets */
  nlsets=i;

  nxtoffset=nlsets;
  for(i=0;i<nlsets;i++)
  {
     indx = sigToInd(getlcharFromDelimit(getSuffFid(lset,i)));
     start[indx] = nxtoffset;
     end[indx] = nxtoffset+getSuffPos(lset,i)-1;
     nxtoffset+=getSuffPos(lset,i);
  }
 
  return nlsets;
} /* end computeBoundFromLset */

/* getLsetLen : returns the sum of lengths of all char lsets under a node */
int getLsetLen(struct suff *lset)
{
   int i=0,lsetlen=0;

   if(lset==NULL) return 0;
   for(i=0;getSuffFid(lset,i)>=2*N;i++) {lsetlen+=getSuffPos(lset,i);} /* just sum the length of lsets */
   return lsetlen; 
} /*  end getLsetLen */


/* getnLsets : returns the number of distinct lsets in a node */
int getnLsets(struct suff *lset)
{
   int i=0;

   if(!lset) return 0;
   for(i=0;getSuffFid(lset,i)>=2*N;i++) {}
   return i; 
} /* end getnLsets */


/* compressLset : Removes delimiters that have len 0 and reallocates the lset array*/
struct suff *compressLset(struct suff *lset,int *compressedby)
{
  struct suff *newlset;
  int i=0,j=0;
  int nlsets=0,lsetlen=0;
  int temp;

  
  if(!lset) return NULL;
  for(i=0;getSuffFid(lset,i)>=2*N;i++) {}
  nlsets=i;

  j=nlsets-1;
  for(i=nlsets-1;i>=0;i--)
  {
    temp = getSuffPos(lset,i);
    if(temp<0) {
    	printf("Debug: getSuffPos returning -ve : %d, fid=%d, at lset[%d], j=%d, nlsets=%d\n", temp,lset[i].fid, i,j,nlsets);
	fflush(stdout);
	assert(temp>=0);
    }

    if(getSuffPos(lset,i)<=0) continue;
    /* if(j!=i)  copySuff(&(lset[j]),lset[i]); */
    if(j!=i)  copySuffix(lset,j,lset,i);
    j--;


    lsetlen+= getSuffPos(lset,i) + 1; /* +1 for the delimiter storage */
  }

  /* just enuf to move the main lset ptr to now point to j+1 index */
  newlset=lset;
  lset+=(j+1);
  *compressedby=(j+1);
  memmove(newlset,lset,sizeof(struct suff)*lsetlen);
  newlset=(struct suff *) realloc(newlset,sizeof(struct suff)*lsetlen);
  return newlset;
  
} /* end compressLset */

/* dispNode : displays a node and its lset contents */
short dispNode(struct stnode *node,int id)
{
   char str[10];
   int notempty=0;
   struct suff *lsets;   
   int offsets[Sigma];
   int nlsets=0;
   int i=0,j=0,len=0;
   int start=0,end=0;
   char lchar;
   int getlcharFromDelimit(int id);
 
   if(node->rmost==id) 
   { 
      strcpy(str,"Leaf"); 
      g_iAvgLeafDepth+=node->depth; 
      if(node->depth<minDepth) minDepth=node->depth;
      if(node->depth>maxDepth) maxDepth=node->depth;
   }
   else strcpy(str,"Node");

   printf("%s : %d depth=%d,rmost=%d\n",str,id,node->depth,node->rmost);

   /* make the change here - see if you still need the len parameter */
   lsets=node->lset;
   /* first compute offsets - start index of each char list */
   nlsets = computeOffsetFromLset(lsets,offsets);  
    
   for(i=0;i<nlsets;i++)
   {
        lchar = getlcharFromDelimit(getSuffFid(lsets,i));
        printf("<%d> Left[%c]->{",getSuffPos(lsets,i),lchar);
        len=getSuffPos(lsets,i);
        start = offsets[sigToInd(lchar)];
        end = start + getSuffPos(lsets,i) - 1;
	for(j=start;j<=end;j++)
        {
           notempty=1;
  	   printf("#%d,%d#,",getSuffFid(lsets,j),getSuffPos(lsets,j));
	}
        printf("} \n");
   }
   /* end change */   

    return notempty;
} /* end dispNode*/


void dispBucket(unsigned int start,unsigned int end,int sig[],int *bucket)
{
  int i;
  int j,t;


  for(i=0;i<Sigma;i++)
  {
      j=sig[i];
      t=0;
      /*  printf("%d,",j);  */
      if(j==-1) continue;
      t++; 
      while(bucket[j]!=-1)
      {
         /* printf("%d,",bucket[j]);  */
         j=bucket[j];
         t++;
      }
      /* printf("%d.",bucket[j]);  */
      printf("\nList<%d> %c is: <%d>\n ",end-start+1,indToSig(i),t);
  }
  fflush(stdout);
} /* end dispBucket */


void copySuff(struct suff *dst,struct suff src)
{
  dst->fid=src.fid;
  dst->pos=src.pos;
} /* end copySuff */

int sigToInd(char c)
{
        switch(c)
        {
          case 'A' :  return 0; break;
          case 'C' :  return 1; break;
          case 'G' :  return 2; break;
          case 'T' :  return 3; break;
          case '$' :  return 4; break;
          case 'B' :  return 4; break;
          default  :  return 0;
        } /* end switch */

} /* sigToInd */


char indToSig(int d)
{
       switch(d)
       {
         case 0 :  return 'A'; break;
         case 1 :  return 'C'; break;
         case 2 :  return 'G'; break;
         case 3 :  return 'T'; break;
         case 4 :  return '$'; break;
         default  :  return '$';
       } /* end switch */
} /* end indToSig */

char indToLftSig(int d)
{
       switch(d)
       {
         case 0 :  return 'A'; break;
         case 1 :  return 'C'; break;
         case 2 :  return 'G'; break;
         case 3 :  return 'T'; break;
         case 4 :  return 'B'; break;
         default  :  return 'B';
       } /* end switch */
} /* end indToLftSig */


int getDelimitId(char c)
{
       /*switch(c)
       {
         case 'A' :  return -1; break;
         case 'C' :  return -2; break;
         case 'G' :  return -3; break;
         case 'T' :  return -4; break;
         case 'B' :  return -5; break;
         default  :  return -5;
       } */
       switch(c)
       {
         case 'A' :  return 2*N; break;
         case 'C' :  return 2*N+1; break;
         case 'G' :  return 2*N+2; break;
         case 'T' :  return 2*N+3; break;
         case 'B' :  return 2*N+4; break;
         default  :  return 2*N+4;
       } 
} /* end getDelimitId */

int getlcharFromDelimit(int id)
{
      /* switch(id)
      {
	 case -1 : return 'A'; break;
	case -2 : return 'C'; break;
	case -3 : return 'G'; break;
	case -4 : return 'T'; break;
	case -5 : return 'B'; break;
	default : return 'B';  
      } */
      
      if(id==(2*N)) return 'A';
      if(id==(2*N+1)) return 'C';
      if(id==(2*N+2)) return 'G';
      if(id==(2*N+3)) return 'T';
      if(id==(2*N+4)) return 'B';
      return 'B';

} /* end getlcharFromDelimit */

char getLeftChar(struct fkp ky)
{
   char *s,c;
   int l;
   s=AllESTs[ky.fragid].s;
   l=strlen(s);

   if(ky.dir=='D')
   {
     if(ky.pos<=0) return 'B';
     c=s[ky.pos-1];
     if(c!='A' && c!='C' && c!='G' && c!='T') return 'B';
     return charToSig(s[ky.pos-1]);             
   }
   else
   {
     /* if(ky.pos >= (l-1)) return 'B'; */
     if(ky.pos <=0 ) return 'B'; 
     c=s[l-ky.pos];
     if(c!='A' && c!='C' && c!='G' && c!='T') return 'B';
     return charToPSig(s[l-ky.pos]);
   }
  
} /* end getLeftChar */


int goTillDifferent(int start,int end,int head,int *pbuck,int depSoFar,struct fkp *buckeyes)
{
  int i,j;
  int d=0;
  char c1,c2;
  struct fkp *key;
  int len=0,comeon=0;
  int iamat=0;
  int bNorX(char c) ;

  if(start>end) return 0;
  if(head== -1) return 0;

  while(1)   
  {
     i=head;
     j=start+i;
     key = &(buckeyes[j]);
     len=strlen(AllESTs[key->fragid].s);
     if(key->dir=='D')
     {
	iamat=key->pos+depSoFar+d; 
        if(iamat<len &&  !bNorX(AllESTs[key->fragid].s[iamat]))
        { 
            c1=charToSig(AllESTs[key->fragid].s[iamat]);
	}
        else
        {
           /* printf("String [%d,%c] has ended at %d+%d=%d\n",key->fragid,key->dir,key->pos,depSoFar+d,len); */
           return d; 
        }
     }
     else
     {
        comeon=len-1-(key->pos)-depSoFar-d;
        if( comeon>=0 && !bNorX(AllESTs[key->fragid].s[comeon]))
	{
           c1=charToPSig(AllESTs[key->fragid].s[comeon]);
	}
        else
        {
       /*          printf("String [%d,%c] has ended at %d\n",key->fragid,key->dir,comeon); */
           return d; 
        }
     }
     g_iNumCharsVisited++;

     i=pbuck[i]; 
     while(i>=0)
     {
        j=start+i;
        key = &(buckeyes[j]);
        len=strlen(AllESTs[key->fragid].s);
        if(key->dir=='D')
        {
	   iamat=key->pos+depSoFar+d;
	   if(iamat<len &&  !bNorX(AllESTs[key->fragid].s[iamat]))
	   {
              c2=charToSig(AllESTs[key->fragid].s[iamat]);
	   }
           else
           {
              /* printf("String [%d,%c] has ended at %d+%d=%d\n",key->fragid,key->dir,key->pos,depSoFar+d,len); */
              return d; 
           }
        }
        else
        {
           comeon=len-1-(key->pos)-depSoFar-d;
           if( comeon>=0 && !bNorX(AllESTs[key->fragid].s[comeon]) )
	   {
              c2=charToPSig(AllESTs[key->fragid].s[comeon]);
	   }
           else
           {
              /* printf("String [%d,%c] has ended at %d\n",key->fragid,key->dir,comeon); */
              return d; 
           }
        }
        if(c1!=c2)
        {
              /* printf("String [%d,%c] differed <%c,%c>\n",key->fragid,key->dir,c2,c1); */
              return d; 
        }

        i=pbuck[i];
        g_iNumCharsVisited++;
     } /* end while */ 
     d++;
  } /* end while 1 */
} /* end goTillDifferent */



void dispTree()
{
   int i;
   printf("Rank=%d: All Nodes Display \n",rank);
   for(i=0;i<g_iNodes;i++)
   {
        printf("<%d,%d,%d> ",i,Nodes[i].depth,Nodes[i].rmost); 
        if(i%5==4) printf("\n"); 
   }
} /* end dispTree */

void dispLeaves()
{
   int i;
   printf("Rank=%d: All Leaves Display \n",rank);
   for(i=0;i<g_iNodes;i++)
   {
        if(i!=Nodes[i].rmost) continue;
        if(!dispNode(&Nodes[i],i)) printf("Error: Leaf left lists probably empty for %d\n",i);     
   }
} /* end dispLeaves */


void dispNodesDetails()
{
  int i;
   printf("Rank=%d: All Nodes Display \n",rank);
   for(i=0;i<g_iNodes;i++)
   {
        dispNode(&Nodes[i],i);     
   }
  
} /* end dispNodesDetails */


/******************************************************************

 * Function to sort locally on a processor
 * Input :  Takes the lower and upper bound of list
 *          Input array to be sorted
 * Note  :  This takes in an array of nodedep structures and
 *          sorts the structures with respect to their depth
*******************************************************************/


void quicksortNodes(int lower, int upper, struct depnode *nodes)
{
        int par1,par2;
        int  down, up,down2,up2;
        int t,x;
        double ran;
        struct depnode pivotElem;
     	int i=0;
        void swapNodes(struct depnode *a,struct depnode *b);
        int cmpNodes(struct depnode a,struct depnode b);
        void copyNodes(struct depnode *,struct depnode);
   	int isAlessB(struct depnode a,struct depnode b);
	int isAequalB(struct depnode a,struct depnode b);

        if(rank==-1) 
        {
	   if(lower==0 && (upper==258))
           {
              for(i=lower;i<=upper;i++) printf("<%d,%d> ",i,nodes[i].depth);
 	      printf("\n");
           }
 	}
       /* generate a random pivot element */
        ran = drand48(); 
        x=(int) floor(ran*(upper-lower+1));
        t= lower + x; 
        swapNodes(&nodes[lower],&nodes[t]);
      
        up = upper;
        down = lower;

        if(lower>=upper)
                return; /*the array is sorted*/

        copyNodes(&pivotElem,nodes[lower]);

        while(down < up)
        {
                while((cmpNodes(nodes[down], pivotElem)) && down < up)
                {
                        down++;
                }
                while(!(cmpNodes(nodes[up], pivotElem)))
                        up--;
                if(down < up)
                {
                        swapNodes(&nodes[down],&nodes[up]);
                }
        } /* end of while*/

        swapNodes(&nodes[lower],&nodes[up]);
        par1=par2 = up;

        /* Adding from here */
        /* Partition lower..par-1 to [< | = ] portions */
        down2 = lower;
        up2 = par2;
    
        while(down2 < up2)
        {
                while((isAlessB(nodes[down2], pivotElem)) && down2 < up2)
                        down2++;
                while(isAequalB(nodes[up2], pivotElem) && down2 < up2)
                        up2--;
                if(down2 < up2)
                {
                        swapNodes(&nodes[down2],&nodes[up2]);
                }
        } 

        par1 = down2;
        
       /* if(rank==0) printf("<%d to %d> <%d to %d> t=%d\n",lower,par1-1,par2+1,upper,t); */
      
        if(lower <( par1-1))
                quicksortNodes(lower, par1-1, nodes);
        if((par2+1)<upper)
                quicksortNodes(par2+1, upper, nodes);

        return;
} /* end of quicksortNodes */


void swapNodes(struct depnode *a,struct depnode *b)
{
   struct depnode temp;
   void copyNodes(struct depnode *,struct depnode);

   copyNodes(&temp,*a);
   copyNodes(a,*b);
   copyNodes(b,temp);
    
} /* end of swapNodes*/

/* Returns  1 if a<=b
          0 if a>b
*/

int cmpNodes(struct depnode a,struct depnode b)
{
      if( (a.depth) <= (b.depth))
                return 1;
      return 0;
}  /* end of cmpNodes */

void copyNodes(struct depnode *a, struct depnode b)
{
   a->id = b.id;
   a->depth = b.depth;
}

/* isAlessB : returns 1 if a.depth<b.depth , 0 otherwise */

int isAlessB(struct depnode a,struct depnode b)
{
  if((a.depth)<(b.depth)) return 1;
  else return 0;
} /* end isAlessB */


/* isAequalB : returns 1 if a.depth == b.depth
                       0 otherwise 
*/
int isAequalB(struct depnode a,struct depnode b)
{
  if(a.depth==b.depth) return 1;
  else return 0;
} /* end isAequalB */


void dispParentChild(int id)
{
   int cid;

   printf("Parent:%d\n",id);
   dispNode(&(Nodes[id]),id);
   cid=id+1;
   while(1)
   {
      printf("Parent:%d, Child:%d\n",id,cid);
      dispNode(&(Nodes[cid]),cid);
      if(Nodes[id].rmost==Nodes[cid].rmost) break;
      cid=Nodes[cid].rmost+1;
   }

} /* end dispParentChild */


int bNorX(char c) {
	   if(c=='N' || c=='X') return 1;
	    return 0;
} /*  bNorX */




/* Counting sort implementation for sorting the nodes of the GST based on string depth*/

void nodesort(int *id,struct stnode *Nodes,int size) {
	int i;
        short curr_depth,next_depth;
	int curr_id,next_id;
	int *c;
	char *b;
	int min,max,range;
	int j,dst;

	if(id==NULL || Nodes==NULL || size<=0) {
	       printf("Warning: nodesort: Node ID array or Nodes array (GST) is empty !!\n");
       	       return;
	}

	min=max=Nodes[0].depth;

	b=(char *)malloc(size);
        for(i=0;i<size;i++) {
                curr_depth=Nodes[i].depth;
		if(max<curr_depth) max=curr_depth;
		if(min>curr_depth) min=curr_depth;
		b[i]='0';
	}
	if(max<0 || min<0) { 
		printf("nodesort: Error: Input string depth array has negative values.\n");
		return;
	}
	range=max-min+1;
	c=(int *)malloc(sizeof(int)*range);
	
	for(i=0;i<range;i++) c[i]=0;
	for(i=0;i<size;i++) c[Nodes[i].depth-min]++;
	/* c[0]=c[0];*/
	for(i=1;i<range;i++) c[i]=c[i-1]+c[i];; /* prefix sum as index */

	curr_id=id[0];
	curr_depth=Nodes[curr_id].depth;
	j=0;
	b[0]='1'; /* index 0 is taken - so mark it*/
	for(i=1;i<=size;i++) {
		dst = c[curr_depth-min]-1;

		next_id=id[dst];
		next_depth=Nodes[next_id].depth;

		id[dst]=curr_id;
		c[curr_depth-min]--;
		if(b[dst]=='0') {
			b[dst]='1'; /* Mark as taken */
			curr_depth=next_depth;
			curr_id=next_id;
			continue;
		}
		else { /* ignore next_depth - already taken care of*/
			while(j<size){
				if(b[j]=='0') break;
				/*printf("	moving to %d\n",j); fflush(stdout);*/
				j++;
			}
			if(j>=size) {
		                /*printf("Counting Sort Finished\n"); fflush(stdout);*/
			    free(b);
			    free(c);
				return;
			}

			curr_id=id[j];
			curr_depth=Nodes[curr_id].depth;
			b[j]='1';
		}
	} /* end for */

	free(b);
	free(c);
} /* end nodesort */
