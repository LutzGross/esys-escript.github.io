/* $Id$ */

/**************************************************************/

/*   Finley: system matrix pattern                            */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2004 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#ifndef INC_FINLEY_SYSTEMPATTERN
#define INC_FINLEY_SYSTEMPATTERN

#include "Finley.h"

/**************************************************************/

/* Typically INDEX_OFFSET and PTR_OFFSET are 1 if the matrix is
   handed over to a library written in FORTRAN. */

#define INDEX_OFFSET 0
#define PTR_OFFSET 0

typedef struct Finley_SystemMatrixPattern {
  dim_t n_ptr;
  dim_t n_index;
  index_t* ptr;
  index_t* index;
  dim_t len;
  dim_t reference_counter;

} Finley_SystemMatrixPattern;


/*  interfaces: */

Finley_SystemMatrixPattern* Finley_SystemMatrixPattern_alloc(dim_t, index_t*,index_t*);
Finley_SystemMatrixPattern* Finley_SystemMatrixPattern_reference(Finley_SystemMatrixPattern*);
void Finley_SystemMatrixPattern_dealloc(Finley_SystemMatrixPattern*);
int Finley_comparIndex(const void *,const void *);
Finley_SystemMatrixPattern* Finley_SystemMatrixPattern_getSubpattern(Finley_SystemMatrixPattern*,dim_t,index_t*,index_t*);
void Finley_SystemMatrixPattern_mis(Finley_SystemMatrixPattern* pattern_p, index_t* mis_marker);
Finley_SystemMatrixPattern* Finley_SystemMatrixPattern_unrollBlocks(Finley_SystemMatrixPattern*, dim_t,dim_t);

#endif /* #ifndef INC_FINLEY_SYSTEMPATTERN */

/*
 * $Log$
 * Revision 1.6  2005/07/08 04:07:57  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.2.6  2005/06/29 02:34:56  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.2.5  2005/03/15 07:23:55  gross
 * Finley's interface to the SCSL library can deal with systems of PDEs now. tests shows that the SCSL library cannot deal with problems with more then 200000 unknowns. problem has been reported to SGI.
 *
 * Revision 1.1.2.4  2005/03/02 23:35:06  gross
 * reimplementation of the ILU in Finley. block size>1 still needs some testing
 *
 * Revision 1.1.2.3  2005/02/18 02:27:31  gross
 * two function that will be used for a reimplementation of the ILU preconditioner
 *
 * Revision 1.1.2.2  2004/11/24 01:37:16  gross
 * some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
 *
 * Revision 1.1.2.1  2004/11/14 23:49:10  gross
 * the forgotten files
 *
 *
 */
