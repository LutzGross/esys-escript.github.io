// @HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2011) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Trilinos_Util.h"

int Trilinos_Util_csrcsc(int n, int n2, int job, int ipos, double * a, 
           int *ja, int *ia, double *ao, int *jao, int *iao)
{
     
    int next, i, j, k;

/*
 -----------------------------------------------------------------------
 Compressed Sparse Row     to      Compressed Sparse Column 
 (transposition operation)   Not in place. 
 -----------------------------------------------------------------------
 Rectangular version.  n is number of rows of CSR matrix, 
                       n2 (input) is number of columns of CSC matrix. 
 -----------------------------------------------------------------------
 -- not in place -- 
 this subroutine transposes a matrix stored in a, ja, ia format. 
 --------------- 
 on entry: 
 ---------- 
 n	= number of rows of CSR matrix. 
 n2    = number of columns of CSC matrix. 
 job	= integer to indicate whether to fill the values (job.eq.0) of the 

         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1) 

 ipos  = starting position in ao, jao of the transposed matrix. 
        the iao array takes this into account (thus iao(1) is set to ipo
s.)
        Note: this may be useful if one needs to append the data structu
re
         of the transpose to that of A. In this case use for example 
                call csrcsc2 (n,n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
 	  for any other normal usage, enter ipos=1. 
 a	= real array of length nnz (nnz=number of nonzero elements in input 

         matrix) containing the nonzero elements. 
 ja	= integer array of length nnz containing the column positions 
 	  of the corresponding elements in a. 
 ia	= integer of size n+1. ia(k) contains the position in a, ja of 
 	  the beginning of the k-th row. 

 on return: 
 ---------- 
 output arguments: 
 ao	= real array of size nzz containing the "a" part of the transpose 
 jao	= integer array of size nnz containing the column indices. 
 iao	= integer array of size n+1 containing the "ia" index array of 
 	  the transpose. 
 -----------------------------------------------------------------------
*/
 // ----------------- compute lengths of rows of transp(A) ----------------
    for (i = 0; i <= n2; i++) iao[i] = 0;

    for (i = 0; i < n; i++) 
	for (k = ia[i]; k <ia[i+1]; k++) 
	    ++iao[ja[k]+1];

// ---------- compute pointers from lengths ------------------------------
 
    iao[0] = ipos;
    for (i = 0; i < n2; i++)
	iao[i+1] = iao[i] + iao[i+1];

// --------------- now do the actual copying -----------------------------
 
    for (i = 0; i < n; i++) {
	for (k = ia[i]; k <ia[i+1]; k++) {
	    j = ja[k];
	    next = iao[j];
	    if (job == 0) {
		ao[next] = a[k];
	    }
	    jao[next] = i;
	    iao[j] = next + 1;
	}
    }
 
  // -------------------------- reshift iao and leave ----------------------
 
    for (i=n2-1; i >= 0; i--) iao[i+1] = iao[i];
    iao[0] = ipos;
 return(0);
}

