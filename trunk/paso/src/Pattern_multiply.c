
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/



/************************************************************************************/

/* Paso: Pattern */

/************************************************************************************/
 
/* Author: Lutz Gross, l.gross@uq.edu.au */

/************************************************************************************/

#include "Paso.h"
#include "Pattern.h"

/************************************************************************************/

/* computes the pattern coming from matrix-matrix multiplication
*
**/

Paso_Pattern* Paso_Pattern_multiply(int type, Paso_Pattern* A, Paso_Pattern* B) {
  Paso_Pattern*out=NULL;
  index_t iptrA,iptrB;
  dim_t i,j,k;
  Paso_IndexListArray* index_list = Paso_IndexListArray_alloc(A->numOutput);
  
  #pragma omp parallel for private(i,iptrA,j,iptrB,k) schedule(static) 
  for(i = 0; i < A->numOutput; i++) {
     for(iptrA = A->ptr[i]; iptrA < A->ptr[i+1]; ++iptrA) {
      j = A->index[iptrA];
      for(iptrB = B->ptr[j]; iptrB < B->ptr[j+1]; ++iptrB) {
    	k = B->index[iptrB];
	Paso_IndexListArray_insertIndex(index_list,i,k);
     }
    }
  }
    
  out=Paso_Pattern_fromIndexListArray(0, index_list,0,B->numInput,0);

 /* clean up */
 Paso_IndexListArray_free(index_list);
  
return out;
}



/*
 * Computes the pattern  of C = A binary operation B for CSR matrices A,B
 *
 * Note: we do not check whether A_ij(op)B_ij=0
 *
 */
Paso_Pattern* Paso_Pattern_binop(int type, Paso_Pattern* A, Paso_Pattern* B) {
  Paso_Pattern *out=NULL;
  index_t iptrA,iptrB;
  dim_t i,j,k;

  Paso_IndexListArray* index_list = Paso_IndexListArray_alloc(A->numOutput);
  
  #pragma omp parallel for private(i,iptrA,j,iptrB,k) schedule(static) 
  for(i = 0; i < B->numOutput; i++){
    iptrA = A->ptr[i],
    iptrB = B->ptr[i];
    
    while (iptrA < A->ptr[i+1] && iptrB < B->ptr[i+1]) {
        j = A->index[iptrA];
        k = B->index[iptrB];
        if (j<k) {
	   Paso_IndexListArray_insertIndex(index_list,i,j);
           iptrA++;
        } else if (j>k) {
	   Paso_IndexListArray_insertIndex(index_list,i,k);
            iptrB++;
        } else if (j==k) {
	   Paso_IndexListArray_insertIndex(index_list,i,j);
            iptrB++;
            iptrA++;
        }
    }
    while(iptrA < A->ptr[i+1]) {
        j = A->index[iptrA];
	Paso_IndexListArray_insertIndex(index_list,i,j);
        iptrA++;
    }
    while(iptrB < B->ptr[i+1]) {
        k = B->index[iptrB];
	Paso_IndexListArray_insertIndex(index_list,i,k);
        iptrB++;
    }
  }
 
  out=Paso_Pattern_fromIndexListArray(0, index_list,0,A->numInput,0);


 /* clean up */
 Paso_IndexListArray_free(index_list);

  return out;
}
