/* $Id$ */
/**************************************************************/

/*   Finley: ElementFile */

/*   allocates and deallocates element table                  */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$                                            */

/**************************************************************/

#include "Finley.h"
#include "ElementFile.h"

/****************************************************************************/

/*  allocates the element table within an element file to hold numElements: */

void Finley_ElementFile_allocTable(Finley_ElementFile* in,int numElements) {
  maybelong *Id2=NULL,*Nodes2=NULL,*Tag2=NULL,*Color2=NULL;
  maybelong numNodes,e,i;
  
  /*  allocate memory: */
  
  numNodes=(maybelong) in->ReferenceElement->Type->numNodes;
  Id2=(maybelong*) MEMALLOC(numElements*sizeof(maybelong));
  Nodes2=(maybelong*) MEMALLOC(numElements*numNodes*sizeof(maybelong));
  Tag2=(maybelong*) MEMALLOC(numElements*sizeof(maybelong));
  Color2=(maybelong*) MEMALLOC(numElements*sizeof(maybelong));
  
  /*  if fine, deallocate the old table and replace by new: */
  
  if (Finley_checkPtr(Id2) || Finley_checkPtr(Nodes2) || Finley_checkPtr(Tag2) || Finley_checkPtr(Color2)) {
    MEMFREE(Id2);
    MEMFREE(Nodes2);
    MEMFREE(Tag2);
    MEMFREE(Color2);
  } else { 
    Finley_ElementFile_deallocTable(in);
    in->numElements=numElements;
    in->Id=Id2;
    in->Nodes=Nodes2;
    in->Tag=Tag2;
    in->Color=Color2;

    /* this initialization makes sure that data are located on the right processor */

    #pragma omp parallel for private(e,i) schedule(static)
    for (e=0;e<numElements;e++) {
       for (i=0;i<numNodes;i++) in->Nodes[INDEX2(i,e,numNodes)]=-1;
       in->Id[e]=-1;
       in->Tag[e]=-1;
       in->Color[e]=-1;
    }
    in->numColors=0;
  }
  return;
}

/*  deallocates the element table within an element file: */

void Finley_ElementFile_deallocTable(Finley_ElementFile* in) {
  MEMFREE(in->Id);
  MEMFREE(in->Nodes);
  MEMFREE(in->Tag);
  MEMFREE(in->Color);
  in->numElements=0;
  in->numColors=0;
}

/* 
* $Log$
* Revision 1.1  2004/10/26 06:53:57  jgs
* Initial revision
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/
