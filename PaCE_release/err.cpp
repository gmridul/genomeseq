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


#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
extern int errno;
#include "err.h"

int spaceSoFar=0,needs=0,effSpace=0,thisPhaseSpace=0,g_TransientMemory=0;
unsigned int g_iGSTscratchSpace=0;
unsigned int g_iSendRecvScratchSpace=0;
 
void checkAlloc(void *newSpace,char *msg)
{
   if(newSpace) spaceSoFar+=needs;
   else AlertErr(msg); 
} /* end check */

/* AlertErr : displays a malloc error occured */
void AlertErr(char *msg)
{
   /*printf("Rank=%d: Malloc Error occured while trying to allocate for <%s>\n",rank,msg);
   printf("Rank=%d: Memory Consumed So Far until start of this phase =%d\n",rank,spaceSoFar);
   printf("Rank=%d: Memory Consumed So Far This Phase =%d; \nMemory needed more = %d\n",rank,thisPhaseSpace+g_TransientMemory,needs);
   printf("Rank=%d: Crash Point = %d bytes\n",rank,thisPhaseSpace+needs+g_TransientMemory);
   */
   
   printf("Rank=%d: Malloc Error occured while allocating <%s>\n",rank,msg);
   printf("Rank=%d: Memory Consumed: Allocated<%d> + Scratch<%u> + NeededMore<%d>=%u bytes\n",
		rank,thisPhaseSpace+g_TransientMemory,g_iGSTscratchSpace+g_iSendRecvScratchSpace,needs,
		thisPhaseSpace+g_TransientMemory+g_iGSTscratchSpace+g_iSendRecvScratchSpace+needs);
   fflush(stdout);
   terminate();
  
} /* end AlertErr */

void terminate()
{
  /*MPI_Finalize();*/
  MPI_Abort(MPI_COMM_WORLD,0); /* Change: Planned: Version 2.1*/
  exit(1);
}

void printErr(int eno)
{
  int e;

  e=eno;
  if(eno==-1) eno=errno;
   
  switch(e)
  {
     case E2BIG:  printf("Error: Arg list too long\n"); fflush(stdout); break;
     case EACCES: printf("Error: Permission denied\n"); fflush(stdout); break;
     case EAGAIN: printf("Error: Resource temporarily unavailable\n"); fflush(stdout); break;
     case EBADF: printf("Error: Bad file descriptor\n"); fflush(stdout); break;
     case EBADMSG: printf("Error: Bad message\n"); fflush(stdout); break;
     case EBUSY: printf("Error: Resource busy\n"); fflush(stdout); break;
     /* case ECANCELED: printf("Error: Operation canceled\n"); fflush(stdout); break;*/
     case ECHILD: printf("Error: No child processes\n"); fflush(stdout); break;
     case EDEADLK: printf("Error: Resource deadlock avoided\n"); fflush(stdout); break;
     case EDOM: printf("Error: Domain error\n"); fflush(stdout); break;
     case EEXIST: printf("Error: File exists\n"); fflush(stdout); break;
     case EFAULT: printf("Error: Bad address\n"); fflush(stdout); break;
     case EFBIG: printf("Error: File too large\n"); fflush(stdout); break;
     case EINPROGRESS: printf("Error: Operation in progress\n"); fflush(stdout); break;
     case EINTR: printf("Error: Interrupted function call\n"); fflush(stdout); break;
     case EINVAL: printf("Error: Invalid argument \n"); fflush(stdout); break;
     case EIO: printf("Error: Input/output error \n"); fflush(stdout); break;
     case EISDIR: printf("Error: Is a directory \n"); fflush(stdout); break;
     case EMFILE: printf("Error: Too many open files \n"); fflush(stdout); break;
     case EMLINK: printf("Error: Too many links \n"); fflush(stdout); break;
     case EMSGSIZE: printf("Error: Inappropriate message buffer length \n"); fflush(stdout); break;
     case ENAMETOOLONG: printf("Error: Filename too long \n"); fflush(stdout); break;
     case ENFILE: printf("Error: Too many open files in system \n"); fflush(stdout); break;
     case ENODEV: printf("Error: No such device \n"); fflush(stdout); break;
     case ENOENT: printf("Error: No such file or directory\n"); fflush(stdout); break;
     case ENOEXEC: printf("Error: Exec format error \n"); fflush(stdout); break;
     case ENOLCK: printf("Error: No locks available \n"); fflush(stdout); break;
     case ENOMEM: printf("Error: Not enough space \n"); fflush(stdout); break;
     case ENOSPC: printf("Error: No space left on device \n"); fflush(stdout); break;
     case ENOSYS: printf("Error: Function not implemented \n"); fflush(stdout); break;
     case ENOTDIR: printf("Error: Not a directory \n"); fflush(stdout); break;
     case ENOTEMPTY: printf("Error: Directory not empty\n"); fflush(stdout); break;
     /*case ENOTSUP: printf("Error: Not supported \n"); fflush(stdout); break;*/
     case ENOTTY: printf("Error: Inappropriate I/O control operation\n"); fflush(stdout); break;
     case ENXIO: printf("Error: No such device or address \n"); fflush(stdout); break;
     case EPERM: printf("Error: Operation not permitted \n"); fflush(stdout); break;
     case EPIPE: printf("Error: Broken pipe \n"); fflush(stdout); break;
     case ERANGE: printf("Error: Result too large\n"); fflush(stdout); break;
     case EROFS: printf("Error:  Read-only file system \n"); fflush(stdout); break;
     case ESPIPE: printf("Error: Invalid seek\n"); fflush(stdout); break;
     case ESRCH: printf("Error: No such process\n"); fflush(stdout); break;
     case ETIMEDOUT: printf("Error: Operation timed out\n"); fflush(stdout); break;
     case EXDEV: printf("Error: Improper link\n"); fflush(stdout); break;


     default :    printf("Error: Some error %d\n",e); fflush(stdout); break;

  }
} /* end printErr */
