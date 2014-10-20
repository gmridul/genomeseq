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



#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "suff.h"

int g_endian=LSBinvalid;
   
                 /* 0 means LSB last endian */
                 /* 1 means LSB first endian */
                 /* 2 means Unidentified endian */

int checkEndian()
{
  int x=1,y=0;
  char *p,*q;

  y=0;
  q=(char *)&y;
  p=(char *)&x;
  p++;
  q++;
  memcpy(q,p,3);
  if(x==y) 
  {
    g_endian=LSBlast;
    return g_endian;
  }

  y=0;
  q=(char *)&y;
  p=(char *)&x;
  memcpy(q,p,3);
  if(x==y) 
  {
    g_endian=LSBfirst;
    return g_endian;
  }

  printf("Neither Big-Endian nor Little-Endian\n"); fflush(stdout);
  g_endian=LSBinvalid;
  return g_endian;

}


struct suff *createSA(int n) 
{
  struct suff *s;
  int len;
  
  len=n;

  s=(struct suff *)calloc(n,sizeof(struct suff));
  /*s=(char *)malloc(sizeof(char)*n*sufsize); */
  /*for(i=0;i<len;i++) s[i]=s[i]&0; */
  return s;
}

void putSuff(struct suff *SA,int index,struct suff sf)
{
	assert(sf.fid>=0 && sf.fid<(2*N+1000));
	assert(sf.pos>=0);
	SA[index].fid = sf.fid;
	SA[index].pos= sf.pos;
}

void getSuff(struct suff *sf,struct suff *SA,int index)
{
	sf->fid = SA[index].fid;
	sf->pos= SA[index].pos;
	assert(sf->fid>=0 && sf->fid<(2*N+1000));
	assert(sf->pos>=0);
}

int getSuffFid(struct suff *SA,int index)
{
  int fid;
  fid = SA[index].fid;
  if(fid<0 || fid>=(2*N+1000)) {
  	printf("Debug: getSuffFid fid is not in range: %d\n",fid);
	fflush(stdout);
  }
  assert(fid>=0 && fid<(2*N+1000));
  return fid; 
}

int getSuffPos(struct suff *SA,int index)
{
  int pos;
  pos = SA[index].pos;
  //assert(pos>=0);
  return pos; 
}

void dispSA(struct suff *SA,int index)
{
  printf("[%d %d]\n",SA[index].fid,SA[index].pos);
}

void copySuffix(struct suff *dSA,int dst,struct suff *sSA,int src)
{

   dSA[dst].fid = sSA[src].fid;
   dSA[dst].pos = sSA[src].pos;

}

void incrSuffPos(struct suff *SA,int index,int addby)
{

  SA[index].pos = SA[index].pos + addby;
  
}
