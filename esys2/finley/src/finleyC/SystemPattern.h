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

#endif /* #ifndef INC_FINLEY_SYSTEMPATTERN */

/*
 * $Log$
 * Revision 1.2  2004/12/15 07:08:33  jgs
 * *** empty log message ***
 *
 * Revision 1.1.2.2  2004/11/24 01:37:16  gross
 * some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
 *
 * Revision 1.1.2.1  2004/11/14 23:49:10  gross
 * the forgotten files
 *
 *
 */
