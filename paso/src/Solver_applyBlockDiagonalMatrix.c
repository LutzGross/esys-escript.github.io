
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: apply block diagonal matrix D: x=D*b               */

/* should be called within a parallel region                                              */
/* barrier synconization should be performed to make sure that the input vector available */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004, 2005 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"

/**************************************************************/


void Paso_Solver_applyBlockDiagonalMatrix(dim_t n_block,dim_t n,double* D,index_t* pivot,double* x,double* b) {
     dim_t i;
     register dim_t i3,i9;
     register double b0,b1,b2,D00,D10,D20,D01,D11,D21,D02,D12,D22;

     if (n_block==1) {
         #pragma omp parallel for private(i) schedule(static)
         for (i=0;i<n;++i) {
            x[i]=D[i]*b[i];
         }
     } else if (n_block==2) {
         #pragma omp parallel for private(i,b0,b1,D00,D10,D01,D11,i3,i9) schedule(static)
         for (i=0;i<n;++i) {
            i3=2*i;
            i9=4*i;
            b0=b[i3];
            b1=b[i3+1];
            D00=D[i9  ];
            D10=D[i9+1];
            D01=D[i9+2];
            D11=D[i9+3];
            x[i3  ]=D00*b0+D01*b1;
            x[i3+1]=D10*b0+D11*b1;
         }
     } else if (n_block==3) {
         #pragma omp parallel for private(i,b0,b1,b2,D00,D10,D20,D01,D11,D21,D02,D12,D22,i3,i9) schedule(static)
         for (i=0;i<n;++i) {
            i3=3*i;
            i9=9*i;
            b0=b[i3];
            b1=b[i3+1];
            b2=b[i3+2];
            D00=D[i9  ];
            D10=D[i9+1];
            D20=D[i9+2];
            D01=D[i9+3];
            D11=D[i9+4];
            D21=D[i9+5];
            D02=D[i9+6];
            D12=D[i9+7];
            D22=D[i9+8];
            x[i3  ]=D00*b0+D01*b1+D02*b2;
            x[i3+1]=D10*b0+D11*b1+D12*b2;
            x[i3+2]=D20*b0+D21*b1+D22*b2;
         }
     }
     return;
}
