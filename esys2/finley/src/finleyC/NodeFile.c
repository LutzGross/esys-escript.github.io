/* $Id$ */
/**************************************************************/
/*                                                             */
/*   Finley: Mesh : NodeFile */
/*                                                             */
/*   allocates and deallocates node files                      */
/*                                                             */
/**************************************************************/

/*   Copyrights by ACcESS Australia 2003/04 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Finley.h"
#include "NodeFile.h"

/**************************************************************/

/*   allocates a node file to hold nodes */
/*   use Finley_NodeFile_allocTable to allocate the node table (Id,Coordinatess). */

Finley_NodeFile* Finley_NodeFile_alloc(int numDim){
  Finley_NodeFile *out;
  
  /*  allocate the return value */
  
  out=MEMALLOC(1,Finley_NodeFile);
  if (Finley_checkPtr(out)) return NULL;
  out->numNodes=0;
  out->numDegreesOfFreedom=0;
  out->reducedNumDegreesOfFreedom=0;
  out->reducedNumNodes=0;
  out->numDim=numDim;
  out->Id=NULL;
  out->Tag=NULL;
  out->Coordinates=NULL;
  out->degreeOfFreedom=NULL;
  out->reducedDegreeOfFreedom=NULL;
  out->toReduced=NULL;
  return out;
}

/*  deallocates a node file: */

void Finley_NodeFile_dealloc(Finley_NodeFile* in) {
  if (in!=NULL) {
     #ifdef Finley_TRACE
     printf("node file is deallocated.\n");
     #endif
     Finley_NodeFile_deallocTable(in);   
     MEMFREE(in);      
  }
}
/* 
* $Log$
* Revision 1.2  2004/12/14 05:39:30  jgs
* *** empty log message ***
*
* Revision 1.1.1.1.2.1  2004/11/24 01:37:14  gross
* some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/
