/* $Id$ */
                                                                                                                                                                                                   
/**************************************************************/
                                                                                                                                                                                                   
/* Finley: Converting an element list into a matrix shape     */
                                                                                                                                                                                                   
/**************************************************************/
                                                                                                                                                                                                   
/* Copyrights by ACcESS Australia 2003,2004 */
/* Author: gross@access.edu.au */
                                                                                                                                                                                                   
/**************************************************************/

#ifndef INC_FINLEY_INDEXLIST
#define INC_FINLEY_INDEXLIST

#include "ElementFile.h" 
#include "System.h" 

/* structure to build system matrix */

#define INDEXLIST_LENGTH 85

typedef struct Finley_IndexList {
  index_t index[INDEXLIST_LENGTH];
  dim_t n;
  struct Finley_IndexList *extension;
} Finley_IndexList;

void Finley_IndexList_insertElements(Finley_IndexList*, Finley_ElementFile*,dim_t, index_t*,dim_t, index_t*);
void Finley_IndexList_insertIndex(Finley_IndexList*, index_t);
void Finley_IndexList_toArray(Finley_IndexList*, index_t*);
dim_t Finley_IndexList_count(Finley_IndexList*);
void Finley_IndexList_free(Finley_IndexList*);

#endif /* #ifndef INC_FINLEY_INDEXLIST */

/*
 * $Log$
 * Revision 1.5  2005/07/08 04:07:51  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.4  2004/12/15 07:08:32  jgs
 * *** empty log message ***
 * Revision 1.1.1.1.2.3  2005/06/29 02:34:50  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1.2.2  2004/11/24 01:37:13  gross
 * some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
 *
 *
 *
 */
