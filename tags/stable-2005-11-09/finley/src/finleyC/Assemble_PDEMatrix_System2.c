/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2003,2004,2005 -  All Rights Reserved              *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

/**************************************************************/

/*    assembles the system of numEq PDEs into the stiffness matrix S: */

/*     -div(A*grad u)-div(B*u)+C*grad u + D*u */

/*      -(A_{k,i,m,j} u_m,j)_i-(B_{k,i,m} u_m)_i+C_{k,m,j} u_m,j-D_{k,m} u_m  */

/*    u has numComp components. */

/*    Shape of the coefficients: */

/*      A = numEqu x numDim x numComp x numDim */
/*      B = numDim x numEqu x numComp  */
/*      C = numEqu x numDim x numComp  */
/*      D = numEqu x numComp  */


/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/


#include "Assemble.h"

/**************************************************************/

void  Finley_Assemble_PDEMatrix_System2(dim_t NS,dim_t numDim,dim_t numQuad,dim_t numEqu,dim_t numComp,
                                              double* S,double* DSDX, double* Vol,
                                              int NN, double* EM_S,
                                              double* A, bool_t extendedA,  
                                              double* B, bool_t extendedB, 
                                              double* C, bool_t extendedC, 
                                              double* D, bool_t extendedD ) {
   dim_t s,r,k,m,i,j,q;
   double rtmp;


   /**************************************************************/
   /*   process A: */
   /**************************************************************/
   if (NULL!=A) {
      if (extendedA) {
         for (s=0;s<NS;s++) {
           for (r=0;r<NS;r++) {
             for (k=0;k<numEqu;k++) {
               for (m=0;m<numComp;m++) {
                 for (i=0;i<numDim;i++) {
                   for (j=0;j<numDim;j++) {
                     for (q=0;q<numQuad;q++) {
                       EM_S[INDEX4(k,m,s,r,numEqu,numComp,NN)]+=Vol[q]*DSDX[INDEX3(s,i,q,NS,numDim)]*
                                              A[INDEX5(k,i,m,j,q,numEqu,numDim,numComp,numDim)]*DSDX[INDEX3(r,j,q,NS,numDim)];
                     }
                   }
                 }
               }
             }
           }
         }
      } else {
         for (s=0;s<NS;s++) {
           for (r=0;r<NS;r++) {
              for (i=0;i<numDim;i++) {
                 for (j=0;j<numDim;j++) {
                    rtmp=0;
                    for (q=0;q<numQuad;q++) rtmp+=Vol[q]*DSDX[INDEX3(s,i,q,NS,numDim)]*DSDX[INDEX3(r,j,q,NS,numDim)];
                    for (k=0;k<numEqu;k++) {
                       for (m=0;m<numComp;m++) {
                           EM_S[INDEX4(k,m,s,r,numEqu,numComp,NN)]+=rtmp*A[INDEX4(k,i,m,j,numEqu,numDim,numComp)];
                       }
                    }
                 }
              }
           }
         }
       }
   }
   /**************************************************************/
   /*   process B: */
   /**************************************************************/
   if (NULL!=B) {
      if (extendedB) {
         for (s=0;s<NS;s++) {
           for (r=0;r<NS;r++) {
             for (k=0;k<numEqu;k++) {
               for (m=0;m<numComp;m++) {
                 for (i=0;i<numDim;i++) {
                   for (q=0;q<numQuad;q++) {
                     EM_S[INDEX4(k,m,s,r,numEqu,numComp,NN)]+=Vol[q]*DSDX[INDEX3(s,i,q,NS,numDim)]*
                                                             B[INDEX4(k,i,m,q,numEqu,numDim,numComp)]*S[INDEX2(r,q,NS)];
                   }
                 }
               }
             }
           }
         }
      } else {
         for (s=0;s<NS;s++) {
           for (r=0;r<NS;r++) {
              for (i=0;i<numDim;i++) {
                  rtmp=0;
                  for (q=0;q<numQuad;q++) rtmp+=Vol[q]*DSDX[INDEX3(s,i,q,NS,numDim)]*S[INDEX2(r,q,NS)];
                  for (k=0;k<numEqu;k++) {
                     for (m=0;m<numComp;m++) {
                        EM_S[INDEX4(k,m,s,r,numEqu,numComp,NN)]+=rtmp*B[INDEX3(k,i,m,numEqu,numDim)];
                     }
                  }
              }
           }
         }
      }
   }
   /**************************************************************/
   /*   process C: */
   /**************************************************************/
   if (NULL!=C) {
     if (extendedC) {
         for (s=0;s<NS;s++) {
           for (r=0;r<NS;r++) {
             for (k=0;k<numEqu;k++) {
               for (m=0;m<numComp;m++) {
                 for (j=0;j<numDim;j++) {
                   for (q=0;q<numQuad;q++) {
                      EM_S[INDEX4(k,m,s,r,numEqu,numComp,NN)]+=Vol[q]*S[INDEX2(s,q,NS)]*
                                                               C[INDEX4(k,m,j,q,numEqu,numComp,numDim)]*DSDX[INDEX3(r,j,q,NS,numDim)];
                   }
                 }
               }
             }
           }
         }
     } else {
         for (s=0;s<NS;s++) {
           for (r=0;r<NS;r++) {
               for (j=0;j<numDim;j++) {
                  rtmp=0;
                  for (q=0;q<numQuad;q++) rtmp+=Vol[q]*S[INDEX2(s,q,NS)]*DSDX[INDEX3(r,j,q,NS,numDim)];
                  for (k=0;k<numEqu;k++) {
                     for (m=0;m<numComp;m++) {
                           EM_S[INDEX4(k,m,s,r,numEqu,numComp,NN)]+=rtmp*C[INDEX3(k,m,j,numEqu,numComp)];
                       }
                  }
               }
           }
         }
     }
   }
   /************************************************************* */
   /* process D */
   /**************************************************************/
   if (NULL!=D) {
     if (extendedD) {
         for (s=0;s<NS;s++) {
           for (r=0;r<NS;r++) {
             for (k=0;k<numEqu;k++) {
               for (m=0;m<numComp;m++) {
                 for (q=0;q<numQuad;q++) {
                   EM_S[INDEX4(k,m,s,r,numEqu,numComp,NN)]+=Vol[q]*S[INDEX2(s,q,NS)]*
                                                                 D[INDEX3(k,m,q,numEqu,numComp)]*S[INDEX2(r,q,NS)];
                 }
               }
             }
           }
         }
     } else {
         for (s=0;s<NS;s++) {
           for (r=0;r<NS;r++) {
               rtmp=0;
               for (q=0;q<numQuad;q++) rtmp+=Vol[q]*S[INDEX2(s,q,NS)]*S[INDEX2(r,q,NS)];
               for (k=0;k<numEqu;k++) {
                   for (m=0;m<numComp;m++) {
                     EM_S[INDEX4(k,m,s,r,numEqu,numComp,NN)]+=rtmp*D[INDEX2(k,m,numEqu)];
                  }
               }
           }
         }
     }
   }
}

/*
 * $Log$
 * Revision 1.3  2005/09/15 03:44:21  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.2.2.1  2005/09/07 06:26:17  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.2  2005/07/08 04:07:46  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.1  2005/06/29 02:34:47  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.2  2004/07/30 04:37:06  gross
 * escript and finley are linking now and RecMeshTest.py has been passed
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */