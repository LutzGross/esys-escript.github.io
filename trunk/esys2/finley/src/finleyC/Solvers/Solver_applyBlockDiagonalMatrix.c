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


void Finley_Solver_applyBlockDiagonalMatrix(int n_block,int n,double* D,maybelong* pivot,double* x,double* b) {
     int i;
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
 * Revision 1.2  2005/03/04 07:12:47  jgs
 * *** empty log message ***
 *
 * Revision 1.1.2.1  2005/03/04 05:09:46  gross
 * a missing file from the ILU reimplementation
 *
 *
 */
