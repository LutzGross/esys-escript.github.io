
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
  maybelong index[INDEXLIST_LENGTH];
  maybelong n;
  struct Finley_IndexList *extension;
} Finley_IndexList;

void Finley_IndexList_insertElements(Finley_IndexList*, Finley_ElementFile*,int, maybelong*,int, maybelong*);
void Finley_IndexList_insertIndex(Finley_IndexList*, maybelong);
void Finley_IndexList_toArray(Finley_IndexList*, maybelong*);
int Finley_IndexList_count(Finley_IndexList*);
void Finley_IndexList_free(Finley_IndexList*);

#endif /* #ifndef INC_FINLEY_INDEXLIST */

/*
 * $Log$
 * Revision 1.4  2004/12/15 07:08:32  jgs
 * *** empty log message ***
 *
 *
 *
 */
