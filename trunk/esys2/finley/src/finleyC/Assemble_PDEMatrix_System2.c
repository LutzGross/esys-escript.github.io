/* $Id$ */

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

/*   Copyrights by ACcESS Australia, 2003 */
/*   author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/


#include "Common.h"
#include "Assemble.h"

/**************************************************************/

void  Finley_Assemble_PDEMatrix_System2(int NS,int numDim,int numQuad,int numEqu,int numComp,
                                              double* S,double* DSDX, double* Vol,
                                              int NN, double* EM_S,
                                              double* A, int extendedA,  
                                              double* B, int extendedB, 
                                              double* C, int extendedC, 
                                              double* D, int extendedD ) {
   int s,r,k,m,i,j,q;
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
 * Revision 1.1  2004/10/26 06:53:57  jgs
 * Initial revision
 *
 * Revision 1.2  2004/07/30 04:37:06  gross
 * escript and finley are linking now and RecMeshTest.py has been passed
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
