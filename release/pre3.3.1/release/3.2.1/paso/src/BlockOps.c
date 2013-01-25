
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "BlockOps.h"
#include "Paso.h"

/**************************************************************

    PASO block operations

**************************************************************/

void Paso_BlockOps_solveAll(dim_t n_block,dim_t n,double* D,index_t* pivot,double* x) {
     dim_t i;
     int failed=0;
     const dim_t block_size=n_block*n_block+17;
     (void)block_size;		/* silence the warning which happens when because var is sometimes ignored by a macro*/
     
     if (n_block==1) {
         #pragma omp parallel for private(i) schedule(static)
         for (i=0;i<n;++i) x[i]*=D[i];
     } else if (n_block==2) {
         #pragma omp parallel for private(i) schedule(static)
         for (i=0;i<n;++i) Paso_BlockOps_MViP_2(&D[4*i], &x[2*i]);

     } else if (n_block==3) {
         #pragma omp parallel for private(i) schedule(static)
         for (i=0;i<n;++i) Paso_BlockOps_MViP_3(&D[9*i], &x[3*i]);
     } else {

	#pragma omp parallel for private(i) schedule(static)
	for (i=0;i<n;++i) {
	   Paso_BlockOps_solve_N(n_block, &x[n_block*i], &D[block_size*i], &pivot[n_block*i], &failed);
	}
     }
     if (failed > 0) {
	Esys_setError(ZERO_DIVISION_ERROR, "Paso_BlockOps_solveAll: solution failed.");
     }
     return;
}
