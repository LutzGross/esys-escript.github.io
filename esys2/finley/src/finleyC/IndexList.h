
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

#define INDEXLIST_LENGTH 200

typedef struct Finley_IndexList {
  int index[INDEXLIST_LENGTH];
  int n;
  struct Finley_IndexList *extension;
} Finley_IndexList;

void Finley_IndexList_insertElements(Finley_IndexList*, Finley_ElementFile*,int,int,maybelong*,int,int,maybelong*,int,Finley_SystemMatrixType);
void Finley_IndexList_insertIndex(Finley_IndexList*, maybelong);
void Finley_IndexList_toArray(Finley_IndexList*, maybelong*);
int Finley_IndexList_count(Finley_IndexList*);
void Finley_IndexList_free(Finley_IndexList*);

#endif /* #ifndef INC_FINLEY_INDEXLIST */

/*
 * $Log$
 * Revision 1.3  2004/12/15 03:48:45  jgs
 * *** empty log message ***
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
