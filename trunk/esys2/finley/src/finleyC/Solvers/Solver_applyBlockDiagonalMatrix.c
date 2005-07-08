/* $Id$ */

/**************************************************************/

/* Finley: apply block diagonal matrix D: x=D*b               */

/* should be called within a parallel region                                              */
/* barrier synconization should be performed to make sure that the input vector available */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004, 2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Common.h"

/**************************************************************/


void Finley_Solver_applyBlockDiagonalMatrix(dim_t n_block,dim_t n,double* D,index_t* pivot,double* x,double* b) {
     dim_t i;
     if (n_block==1) {
         #pragma omp for private(i) schedule(static)
         for (i=0;i<n;i++) {
            x[i]=D[i]*b[i];
         }
     } else if (n_block==2) {
         #pragma omp for private(i) schedule(static)
         for (i=0;i<n;i++) {
            x[2*i  ]=D[4*i  ]*b[2*i]+D[4*i+2]*b[2*i+1];
            x[2*i+1]=D[4*i+1]*b[2*i]+D[4*i+3]*b[2*i+1];
         }
     } else if (n_block==3) {
         #pragma omp for private(i) schedule(static)
         for (i=0;i<n;i++) {
            x[3*i  ]=D[9*i  ]*b[3*i]+D[9*i+3]*b[3*i+1]+D[9*i+6]*b[3*i+2];
            x[3*i+1]=D[9*i+1]*b[3*i]+D[9*i+4]*b[3*i+1]+D[9*i+7]*b[3*i+2];
            x[3*i+2]=D[9*i+2]*b[3*i]+D[9*i+5]*b[3*i+1]+D[9*i+8]*b[3*i+2];
         }
     }
     return;
}

/*
 * $Log$
 * Revision 1.3  2005/07/08 04:08:00  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.2  2005/03/04 07:12:47  jgs
 * *** empty log message ***
 *
 * Revision 1.1.2.2  2005/06/29 02:35:00  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.2.1  2005/03/04 05:09:46  gross
 * a missing file from the ILU reimplementation
 *
 *
 */
