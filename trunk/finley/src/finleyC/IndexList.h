/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2003,2004,2005 -  All Rights Reserved              *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
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
 * Revision 1.6  2005/09/15 03:44:22  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.5.2.1  2005/09/07 06:26:19  gross
 * the solver from finley are put into the standalone package paso now
 *
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
