/* $Id$ */

/**************************************************************/

/*    Updates the element matrices: */

/*    assembles the system of numEq PDEs into the stiffness matrix S: */

/*     -div(A*grad u)-div(B*u)+C*grad u + D*u */

/*      -(A_{i,j} u_j)_i-(B_{i} u)_i+C_{j} u_j-D u */

/*    u has numComp components. */

/*    Shape of the coefficients: */

/*      A = numDim x numDim */
/*      B = numDim  */
/*      C = numDim */
/*      D = 1 */


/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003 */
/*   author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"
#include "Assemble.h"

/**************************************************************/
void Finley_Assemble_PDEMatrix_Single2(dim_t NS,dim_t numDim,dim_t numQuad,
                                              double* S,double* DSDX, double* Vol,
                                              dim_t NN, double* EM_S,
                                              double* A, bool_t extendedA,  
                                              double* B, bool_t extendedB, 
                                              double* C, bool_t extendedC, 
                                              double* D, bool_t extendedD ) {
   dim_t s,r,i,j,q;
   double rtmp;

   for (s=0;s<NS;s++) {
     for (r=0;r<NS;r++) {
        /**************************************************************/
         /*   process A: */
         /**************************************************************/
         if (NULL!=A) {
            if (extendedA) {
               for (i=0;i<numDim;i++) {
                 for (j=0;j<numDim;j++) {
                   for (q=0;q<numQuad;q++) {
                     EM_S[INDEX4(0,0,s,r,1,1,NN)]+=Vol[q]*DSDX[INDEX3(s,i,q,NS,numDim)]*
                                                        A[INDEX3(i,j,q,numDim,numDim)]*DSDX[INDEX3(r,j,q,NS,numDim)];
                   
                 }
               }
             }
            } else {
               for (i=0;i<numDim;i++) {
                 for (j=0;j<numDim;j++) {
                   if (A[INDEX2(i,j,numDim)]!=0) {
                      rtmp=0;
                      for (q=0;q<numQuad;q++) rtmp+=Vol[q]*DSDX[INDEX3(s,i,q,NS,numDim)]*DSDX[INDEX3(r,j,q,NS,numDim)];
                      EM_S[INDEX4(0,0,s,r,1,1,NN)]+=rtmp*A[INDEX2(i,j,numDim)];
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
              for (i=0;i<numDim;i++) {
                for (q=0;q<numQuad;q++) {
                  EM_S[INDEX4(0,0,s,r,1,1,NN)]+=Vol[q]*DSDX[INDEX3(s,i,q,NS,numDim)]*B[INDEX2(i,q,numDim)]*S[INDEX2(r,q,NS)];
                }
              }
           } else {
              for (i=0;i<numDim;i++) {
                if (B[i]!=0) {
                   rtmp=0;
                   for (q=0;q<numQuad;q++) rtmp+=Vol[q]*DSDX[INDEX3(s,i,q,NS,numDim)]*S[INDEX2(r,q,NS)];
                   EM_S[INDEX4(0,0,s,r,1,1,NN)]+=rtmp*B[i];
                }
              }
            }
         }
         /**************************************************************/
         /*   process C: */
         /**************************************************************/
         if (NULL!=C) {
           if (extendedC) {
             for (j=0;j<numDim;j++) {
               for (q=0;q<numQuad;q++) {
                 EM_S[INDEX4(0,0,s,r,1,1,NN)]+=Vol[q]*S[INDEX2(s,q,NS)]*C[INDEX2(j,q,numDim)]*DSDX[INDEX3(r,j,q,NS,numDim)];
               }
             }
           } else {
             for (j=0;j<numDim;j++) {
               if (C[j]!=0) {
                  rtmp=0;
                  for (q=0;q<numQuad;q++) rtmp+=Vol[q]*S[INDEX2(s,q,NS)]*DSDX[INDEX3(r,j,q,NS,numDim)];
                  EM_S[INDEX4(0,0,s,r,1,1,NN)]+=rtmp*C[j];
               }
             }
           }
         }
         /************************************************************* */
         /* process D */
         /**************************************************************/
         if (NULL!=D) {
           if (extendedD) {
             for (q=0;q<numQuad;q++) {
                EM_S[INDEX4(0,0,s,r,1,1,NN)]+=Vol[q]*S[INDEX2(s,q,NS)]*D[q]*S[INDEX2(r,q,NS)];
             }
           } else {
             if (D[0]!=0) {
                rtmp=0;
                for (q=0;q<numQuad;q++) rtmp+=Vol[q]*S[INDEX2(s,q,NS)]*S[INDEX2(r,q,NS)];
                EM_S[INDEX4(0,0,s,r,1,1,NN)]+=rtmp*D[0];
             }
           }
         }
     }
  }
}
/*
 * $Log$
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
