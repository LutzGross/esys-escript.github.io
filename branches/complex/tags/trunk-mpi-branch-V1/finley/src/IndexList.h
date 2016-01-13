/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/
/**************************************************************/

/* Finley: Converting an element list into a matrix shape     */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#ifndef INC_FINLEY_INDEXLIST
#define INC_FINLEY_INDEXLIST

#include "Finley.h"
#include "ElementFile.h"
#include "Mesh.h"

/* structure to build system matrix */

#define INDEXLIST_LENGTH 85

typedef struct Finley_IndexList {
  index_t index[INDEXLIST_LENGTH];
  dim_t n;
  struct Finley_IndexList *extension;
} Finley_IndexList;

int Finley_IndexList_localToGlobal(Finley_NodeDistribution* dofDistribution, int localIndex);
void Finley_IndexList_insertElements(Finley_IndexList*, Finley_Mesh*, Finley_ElementFile*,dim_t, index_t*,dim_t, index_t*);
void Finley_IndexList_insertIndex(Finley_IndexList*, index_t);
void Finley_IndexList_toArray(Finley_IndexList*, index_t*);
dim_t Finley_IndexList_count(Finley_IndexList*);
void Finley_IndexList_free(Finley_IndexList*);

/* this will go */
#ifdef PASO_MPI
int Finley_IndexList_localToGlobal(Finley_NodeDistribution *dofDistribution, int localIndex);
#endif

#endif /* #ifndef INC_FINLEY_INDEXLIST */
