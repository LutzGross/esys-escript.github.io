/* $Id$ */

/**************************************************************/

/*    Updates the element matrices: */

/*    assembles the system of numEq components right hand side F */

/*     -div X + Y */

/*      -(X_{k,i})_i + Y_k */

/*    Shape of the coefficients: */

/*      X = numEqu x numDim   */
/*      Y = numEqu */


/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003 */
/*   author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"

/**************************************************************/

void Finley_Assemble_RHSMatrix_System(int NS,int numDim,int numQuad,int numEqu,
                                           double* S,double* DSDX,double* Vol,
                                           int NN, double* EM_F,
                                           double* X, int extendedX,  
                                           double* Y, int extendedY) {
   int s,k,q,i;
   double rtmp;

         /**************************************************************/
         /*   process X: */
         /**************************************************************/
         if (NULL!=X) {
            if (extendedX) {
               for (s=0;s<NS;s++) {
                  for (k=0;k<numEqu;k++) {
                     for (i=0;i<numDim;i++) {
                        for (q=0;q<numQuad;q++) {
                           EM_F[INDEX2(k,s,numEqu)]+=Vol[q]*DSDX[INDEX3(s,i,q,NS,numDim)]*X[INDEX3(k,i,q,numEqu,numDim)];
                        }
                     }
                  }
               }
            } else {
               for (s=0;s<NS;s++) {
                  for (i=0;i<numDim;i++) {
                     rtmp=0;
                     for (q=0;q<numQuad;q++) rtmp+=Vol[q]*DSDX[INDEX3(s,i,q,NS,numDim)];
                     for (k=0;k<numEqu;k++) EM_F[INDEX2(k,s,numEqu)]+=rtmp*X[INDEX2(k,i,numEqu)];
                  }
               }
            }
         }
         /**************************************************************/
         /*   process Y: */
         /**************************************************************/
         if (NULL!=Y) {
            if (extendedY) {
               for (s=0;s<NS;s++) {
                  for (k=0;k<numEqu;k++) {
                     for (q=0;q<numQuad;q++) EM_F[INDEX2(k,s,numEqu)]+=Vol[q]*S[INDEX2(s,q,NS)]*Y[INDEX2(k,q,numEqu)];
                  }
               }
             } else {
               for (s=0;s<NS;s++) {
                   rtmp=0;
                   for (q=0;q<numQuad;q++) rtmp+=Vol[q]*S[INDEX2(s,q,NS)];
                   for (k=0;k<numEqu;k++) EM_F[INDEX2(k,s,numEqu)]+=rtmp*Y[k];
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
