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


/* This is a .h file for compilation purposes. */

/* processLeaf : generate pairs for the leaf identified by Nodes[id] */
/*#define Report_Only_Contained*/

void printMaximalSubstring(struct pair p,int id,FILE *fp);
void printAPairInBatches(struct pair p,FILE *fp);
int getSuffPos(struct suff *,int );
int enQBuf(struct cBuffer *,struct pair );
int computeBoundFromLset(struct suff *lset,int *start,int *end);
int computeOffsetFromLset(struct suff *,int *);
int getlcharFromDelimit(int );
int getLsetLen(struct suff *);
void incrSuffPos(struct suff *,int ,int);

int processLeaf(int id,struct cBuffer *sBuf,int howmany,int *done)
{
   int f1,f2,dept=0;
   short p1,p2;
   char dir;
   int minfrag;
   int pcount=0;
   struct pair newpair;
   int suspend_lLset=0,suspend_rLset=0,suspend_lsuff=0,suspend_rsuff=0;
   int lstart=0,lend=0;
   int rstart=0,rend=0;
   
   struct suff *lsets=NULL;
   int offsets[Sigma];
   int nlsets=0;
   char lLchar,rLchar;

   /* static function variables : */
   static int s_suspended=0;
   static int s_lLset=0,s_rLset=0;
   static int s_lsuff=0,s_rsuff=0;



   dept = Nodes[id].depth;

   if(s_suspended) 
   {
      s_suspended=0;
      suspend_lLset=suspend_rLset=suspend_lsuff=suspend_rsuff=1;
   }


   /* lset change from here */
   lsets=Nodes[id].lset;

   /* first compute offsets - start index of each char list */
   nlsets = computeOffsetFromLset(lsets,offsets);

   if(!suspend_lLset) s_lLset=0;
   else suspend_lLset=0;

   for(;s_lLset<nlsets;s_lLset++)
   {
        lLchar = getlcharFromDelimit(getSuffFid(lsets,s_lLset));
        lstart = offsets[sigToInd(lLchar)];
        lend = lstart + getSuffPos(lsets,s_lLset) - 1;
          
        if(!suspend_rLset)
        {
           if(lLchar=='B')
    	   {
              /* so this is a BLANK (B) left list - cross it with itself in addition */
              s_rLset=s_lLset;
	   }
           else
           {
              /* so this is a A,C,G,T list - so cross only with whatever is ahead */
              s_rLset=s_lLset + 1;
	   }
        }
	else suspend_rLset=0;  /* just remember the previous value */
	

        for(;s_rLset<nlsets;s_rLset++)
        {
             rLchar = getlcharFromDelimit(getSuffFid(lsets,s_rLset));
             rstart = offsets[sigToInd(rLchar)];
             rend = rstart + getSuffPos(lsets,s_rLset) - 1;
                
             if(!suspend_lsuff)  s_lsuff=lstart;
             else suspend_lsuff=0;

             for(;s_lsuff<=lend;s_lsuff++) 
	     {
                 f1 = getSuffFid(lsets,s_lsuff);
                 p1 = getSuffPos(lsets,s_lsuff);

                 if(!suspend_rsuff)  
		 {
	   	    if(lLchar==rLchar) s_rsuff=s_lsuff+1; /* for BLANK B */
                    else s_rsuff=rstart;
                 }
                 else suspend_rsuff=0;

		 while(s_rsuff<=rend)
 		 {
                    f2 = getSuffFid(lsets,s_rsuff);
                    p2 = getSuffPos(lsets,s_rsuff);
                    if(f1/2 == f2/2) { s_rsuff++; continue; }
                    if(f2%2==1 || f1%2==1) dir='P';
                    else dir='D';
                    minfrag=f1;
                    if(f2<minfrag) minfrag=f2;
                    if(minfrag%2==1) { s_rsuff++; continue; }
#ifdef Report_Only_Contained
		    if(p1==0 && AllESTlen[f1/2]==dept) { }
		    else 
		    if(p2==0 && AllESTlen[f2/2]==dept) { }
		    else {
			s_rsuff++;
			continue;
		    }

                    if(f1%2==1) {
                        p1 = AllESTlen[f1/2] - p1 - dept;
                    }
                    if(f2%2==1) {
                        p2 = AllESTlen[f2/2] - p2 - dept;
                    }

                  /* printf("\nrank=%d: <%d,%d,%d,%d,%d,%c> ",rank,dept,f1/2,f2/2,p1,p2,dir);
                    fflush(stdout); */
                    /* populate sBuf */
		    if(p1==0 && AllESTlen[f1/2]==dept) {
			/* f1 is contained in f2 */
                    	newpair.f1 = f1/2;
                    	newpair.f2 = f2/2;
                    	newpair.p1 = p1;
                    	newpair.p2 = p2;
		    }
		    else {
			/* f2 is contained in f1 */
                    	newpair.f1 = f2/2;
                    	newpair.f2 = f1/2;
                    	newpair.p1 = p2;
                    	newpair.p2 = p1;
		    }
                    newpair.l1 = dept;
                    newpair.l2 = dept;
                    newpair.d = 0;
                    newpair.e = 0;
                    newpair.type= dir;

#else
                    if(f1%2==1) {
                        p1 = AllESTlen[f1/2] - p1 - dept;
                    }
                    if(f2%2==1) {
                        p2 = AllESTlen[f2/2] - p2 - dept;
                    }

                  /* printf("\nrank=%d: <%d,%d,%d,%d,%d,%c> ",rank,dept,f1/2,f2/2,p1,p2,dir);
                    fflush(stdout); */
                    /* populate sBuf */
                    newpair.f1 = f1/2;
                    newpair.f2 = f2/2;
                    newpair.p1 = p1;
                    newpair.p2 = p2;
                    newpair.l1 = dept;
                    newpair.l2 = dept;
                    newpair.d = 0;
                    newpair.e = 0;
                    newpair.type= dir;
#endif
   				if( (newpair.f1<0 || newpair.f1>=N) || (newpair.f2<0 || newpair.f2>=N)) {
				          printf("Rank=%d: Error: processLeaf f1 or f2 out of range %d %d for node %d: other depth %d, p1,p2 %d,%d\n",rank,newpair.f1,newpair.f2,id,dept,p1,p2);
				          fflush(stdout);
				           //MPI_Abort(MPI_COMM_WORLD,1);
						s_rsuff++;
					    continue;
			        }

                    enQBuf(sBuf,newpair);
		    if(g_bReportMaximalSubstrings) printMaximalSubstring(newpair,id,fpMaximalSubstrings);
		    if(g_bReportGeneratedPairs) printAPairInBatches(newpair,fpGeneratedPairs);
                    g_iGenPairs++;
                    s_rsuff++;

                    pcount++;

                    if(pcount>=howmany)
                    {
                        /* save and suspend here to return later on */
                        s_suspended=1;
                        return pcount;
                    }
		     
		 } /* end for s_rsuff */
     	     } /* end for s_lsuff  */
              
	} /* end for s_rLset */
   } /* end for s_lLset */

   /* end of lset change */


   /* re-initialize all the static variables to be on the safer side */
   s_suspended=0;
   s_lLset=s_rLset=0;
   s_lsuff=s_rsuff=0;

   *done=1;
   return pcount;

} /* end processLeaf */


/* processNode : Processes an internal node to generate promising pairs on-demand */

int processNode(int id,struct cBuffer *sBuf,int howmany,int *done)
{
  
  int f1,f2;
  short p1,p2;
  int cid=0, dept=0,minfrag=-1;
  char dir;
  struct pair newpair;
  void dispParentChild(int id);
  struct suff *compressLset(struct suff *,int *);

  int x=0,y=0,z=0;
  int end=0,strt=0,pcount=0,resLeaf=0;
  int c_start[Sigma][Sigma];
  int c_end[Sigma][Sigma];
  int newlsetlen =0,compressedby=0;
  struct suff *newlset,*lsptr,*lsets;
  int nxtsuff=0;
  int lstart=0,lend=0,rstart=0,rend=0;
  struct suff dumsuf;

  int suspended_i=0,suspended_j=0,suspended_k=0;
  int suspended_l=0,suspended_m=0,suspended_n=0;
  int already_dup=0;

   
  /* static variable for processNode */

  static int s_leafdone=0,s_suspended=0;
  static int s_i=0,s_j=0,s_k=0,s_l=0,s_m=0,s_n=0;
  static int s_child[Sigma];
  static int s_nc=0;

  static int s_nstart[Sigma][Sigma];
  static int s_nend[Sigma][Sigma];

  dept = Nodes[id].depth;
 /* printf("Node depth =%d\n",dept);
  if(dept==83) dispNode(&(Nodes[id]),id); */

  /* 1.0 : first generate pairs out of the $ branch of the node */
  if(!s_leafdone) 
  {
      resLeaf=processLeaf(id,sBuf,howmany,&s_leafdone);
      if(!s_leafdone || resLeaf>=howmany )
      {
        return resLeaf;
      }
  }

  pcount=0;

  /* 1.0 : Just count here the number of children*/


  if(!s_suspended)
  {
    x=0;
    s_child[x]=id+1;
    while(1)   /* executes for each child of Nodes[id] */
    {
      /* new code here */
      cid = s_child[x];
       x++;
       if(Nodes[cid].rmost==Nodes[id].rmost)
       {
         s_child[x]=id; /* point to self bcos $ list is contained within node */
	 x++;
	 break;
	 
       }
       else
       {
         s_child[x] = Nodes[s_child[x-1]].rmost + 1;
	 continue;
       }
    
       /* new code ends here */
    } /* end while cid*/
 
    s_nc=x; /* number of children for Nodes[id] */

  }
  else 
  { /* suspended before - so release the suspended flag */
    s_suspended=0;
    suspended_i=suspended_j=suspended_m=1;
    suspended_k=suspended_l=suspended_n=1;
    already_dup=1;  /* init to 0 */
  }


   /* lset change starts here */
  /* 1.1 : mark the first occurences of strings of every child  and in the process
            build the new lset array for this parent internal node 
  */
  /* Do this only once per node - static */

   if(!already_dup)
   {
     g_IntNodesProcessed++;

     newlsetlen =0;
     for(x=0;x<s_nc;x++)
     {  
        newlsetlen+= getLsetLen(Nodes[s_child[x]].lset);
     }
     newlsetlen+=Sigma;  /* Add 5 for sigma delimiters - this will get realloc later */

    /*if(rank==0) { printf("Rank=%d: Stage 1 <%d,%d> \n",rank,sufsize,newlsetlen); fflush(stdout); } */
     needs=sizeof(struct suff)*newlsetlen;
     newlset = createSA(newlsetlen);
     if(!newlset) AlertErr("newlset");

     for(x=0;x<s_nc;x++) 
     {
       computeBoundFromLset(Nodes[s_child[x]].lset,c_start[x],c_end[x]);
     }
   
     for(nxtsuff=0;nxtsuff<Sigma;nxtsuff++)
     { /* initialize the delimiters of the newlset */
        dumsuf.fid = getDelimitId(indToLftSig(nxtsuff));
        dumsuf.pos=0;

        putSuff(newlset,nxtsuff,dumsuf);

     }

     /* initialize start and end offset matrix for child-lchar mapping */
     for(x=0;x<s_nc;x++)
     {
        for(y=0;y<Sigma;y++) 
        {
	  s_nstart[x][y]=s_nend[x][y]=-1;
        }
     }

     /* traverse each lset of each child node and remove duplicates and
      simulataneously copy the suffixes to the newlset array */
  
     for(x=0;x<Sigma;x++)
     {
        for(y=0;y<s_nc;y++)
        {
           strt=c_start[y][x];
           end=c_end[y][x];
           if(strt==-1 || end==-1) continue;  /* list empty */

           lsptr=Nodes[s_child[y]].lset;
           for(z=strt;z<=end;z++)
           {
              if(g_NodeMarker[getSuffFid(lsptr,z)]==id) continue; 
        
              g_NodeMarker[getSuffFid(lsptr,z)]=id;
              copySuffix(newlset,nxtsuff,lsptr,z);
              if(s_nstart[y][x]==-1) s_nstart[y][x]=nxtsuff;
              s_nend[y][x]=nxtsuff;
              nxtsuff++;

              incrSuffPos(newlset,x,1);

           } /* end for z*/
        } /* end for y children */
     } /* end for x Sigma */

     
     compressedby=0;
     newlset = compressLset(newlset,&compressedby);

     for(x=0;x<s_nc;x++)
     {
        for(y=0;y<Sigma;y++)  
        {
           if(s_nstart[x][y]==-1) continue;
           s_nstart[x][y]-=compressedby;
           s_nend[x][y]-=compressedby;
        }
     }

     /* free all children' lsets except this node's own lsets */
     for(x=0;x<s_nc;x++)
     {
	     if(Nodes[s_child[x]].lset!=NULL) free(Nodes[s_child[x]].lset); /* Ananth Change: for AiX/Sri 1/29/2004 */
     } 

     Nodes[id].lset = newlset;

    
   } /* end if already_dup */

   /* 1.2 Form pairs by cross-product */
   
   lsets=Nodes[id].lset;
   
   if(!suspended_i) s_i=0;
   else suspended_i=0;

   for(;s_i<s_nc;s_i++)
   {
     if(!suspended_j) s_j=0;
     else suspended_j=0;

     for(;s_j<Sigma;s_j++)
     {
        if(s_nstart[s_i][s_j]==-1) continue;

        lstart=s_nstart[s_i][s_j];
        lend=s_nend[s_i][s_j];
    
        if(!suspended_m) s_m=lstart;
        else suspended_m=0;

        for(;s_m<=lend;s_m++)
        {
          f1 = getSuffFid(lsets,s_m);
          p1 = getSuffPos(lsets,s_m);

          if(!suspended_k) s_k=s_i+1;
          else suspended_k=0;

          for(;s_k<s_nc;s_k++)
          {
             if(!suspended_l) s_l=0;
             else suspended_l=0;

             for(;s_l<Sigma;s_l++)
             {
                if(s_l==s_j && indToLftSig(s_l)!='B') continue;
                if(s_nstart[s_k][s_l]==-1) continue;
             
                rstart=s_nstart[s_k][s_l];
                rend=s_nend[s_k][s_l];

                if(!suspended_n) s_n=rstart;
                else suspended_n=0;

                while(s_n<=rend)
                {
                           f2 = getSuffFid(lsets,s_n);
                           p2 = getSuffPos(lsets,s_n);
		           if(f1/2 == f2/2) { s_n++; continue; }
                           if(f2%2==1 || f1%2==1) dir='P';
                           else dir='D';
                           minfrag=f1;
                           if(f2<minfrag) minfrag=f2;
                           if(minfrag%2==1) { s_n++; continue; }
#ifdef Report_Only_Contained
		    	   if(p1==0 && AllESTlen[f1/2]==dept) { }
		    	   else 
		           if(p2==0 && AllESTlen[f2/2]==dept) { }
		           else {
				s_n++;
				continue;
			   }

                           if(f1%2==1) {
                        	p1 = AllESTlen[f1/2] - p1 - dept;
                    	   }
                    	   if(f2%2==1) {
                        	p2 = AllESTlen[f2/2] - p2 - dept;
                    	   }

                           /* populate sBuf */
		           if(p1==0 && AllESTlen[f1/2]==dept) {
			   /* f1 is contained in f2 */
                    	   	newpair.f1 = f1/2;
                    	   	newpair.f2 = f2/2;
                    	   	newpair.p1 = p1;
                    	   	newpair.p2 = p2;
		    	   }
		           else {
				/* f2 is contained in f1 */
                    		newpair.f1 = f2/2;
                    		newpair.f2 = f1/2;
                    		newpair.p1 = p2;
                    		newpair.p2 = p1;
		    	   }
                    	   newpair.l1 = dept;
                           newpair.l2 = dept;
                           newpair.d = 0;
                           newpair.e = 0;
                           newpair.type= dir;
#else

	                   if(f1%2==1)
                           {
		              p1 = AllESTlen[f1/2] - p1 - dept;
		           }
	                   if(f2%2==1)
                           {
		              p2 = AllESTlen[f2/2] - p2 - dept;
		           }
                           /*printf("\nrank:%d <%d,%d,%d,%d,%d,%c> ",rank,dept,f1/2,f2/2,p1,p2,dir);  */
                           /* populate sBuf */
                           newpair.f1 = f1/2;
                           newpair.f2 = f2/2;
                           newpair.p1 = p1;
                           newpair.p2 = p2;
                           newpair.l1 = dept;
                           newpair.l2 = dept;
                           newpair.d = 0;
                           newpair.e = 0;
                           newpair.type= dir;
#endif
   				if( (newpair.f1<0 || newpair.f1>=N) || (newpair.f2<0 || newpair.f2>=N)) {
				          printf("Rank=%d: Error: processNode f1 or f2 out of range %d %d for node %d: other depth %d, p1,p2 %d,%d\n",rank,newpair.f1,newpair.f2,id,dept,p1,p2);
				          fflush(stdout);
				           //MPI_Abort(MPI_COMM_WORLD,1);
						s_n++;
					    continue;
			        }

                           enQBuf(sBuf,newpair);
		           if(g_bReportMaximalSubstrings) printMaximalSubstring(newpair,id,fpMaximalSubstrings);

		           if(g_bReportGeneratedPairs) printAPairInBatches(newpair,fpGeneratedPairs);
                           g_iGenPairs++;
                           pcount++;
                           s_n++;

			   if(dept==-1) {
			     printf("GeneratePair: <%d,%d> %d,%d,%d,%c\n",newpair.f1,newpair.f2,newpair.p1,newpair.p2,newpair.l1,newpair.type); fflush(stdout);
			   }

                           if( (pcount+resLeaf)>=howmany)
                           {
			       /* save and suspend here to return later on */
                               s_suspended=1; 
                               return (pcount+resLeaf);
 	                   }
                } /* end s_n while */
             } /* end s_l for */
          } /* end s_k for */
        } /* end s_m for */
     }  /* end s_j for */
   } /* end s_i for */
   

   /* end change */

   /* re-initialize static variables here */
   s_leafdone=0;
   s_suspended=0;
   s_i=s_j=s_k=s_l=s_m=s_n=0;
   s_nc=0;
  
   *done=1;
   return (pcount+resLeaf);

} /* end processNode */


