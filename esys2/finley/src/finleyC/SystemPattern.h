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
  maybelong n_ptr;
  maybelong n_index;
  maybelong* ptr;
  maybelong* index;
  size_t len;
  maybelong reference_counter;

} Finley_SystemMatrixPattern;


/*  interfaces: */

Finley_SystemMatrixPattern* Finley_SystemMatrixPattern_alloc(int, maybelong*,maybelong*);
Finley_SystemMatrixPattern* Finley_SystemMatrixPattern_reference(Finley_SystemMatrixPattern*);
void Finley_SystemMatrixPattern_dealloc(Finley_SystemMatrixPattern*);
int Finley_comparIndex(const void *,const void *);
Finley_SystemMatrixPattern* Finley_SystemMatrixPattern_getSubpattern(Finley_SystemMatrixPattern*,int,maybelong*,maybelong*);
void Finley_SystemMatrixPattern_mis(Finley_SystemMatrixPattern* pattern_p, maybelong* mis_marker);
Finley_SystemMatrixPattern* Finley_SystemMatrixPattern_unrollBlocks(Finley_SystemMatrixPattern*, maybelong,maybelong);

#endif /* #ifndef INC_FINLEY_SYSTEMPATTERN */
