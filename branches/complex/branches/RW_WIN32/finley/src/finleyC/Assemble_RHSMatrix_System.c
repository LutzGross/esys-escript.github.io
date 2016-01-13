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

/*    Updates the element matrices: */

/*    assembles the system of numEq components right hand side F */

/*     -div X + Y */

/*      -(X_{k,i})_i + Y_k */

/*    Shape of the coefficients: */

/*      X = numEqu x numDim   */
/*      Y = numEqu */


/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "Assemble.h"

/**************************************************************/

void Finley_Assemble_RHSMatrix_System(dim_t NS,dim_t numDim,dim_t numQuad,dim_t numEqu,
                                           double* S,double* DSDX,double* Vol,
                                           dim_t NN, double* EM_F,
                                           double* X, bool_t extendedX,  
                                           double* Y, bool_t extendedY) {
   dim_t s,k,q,i;
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
