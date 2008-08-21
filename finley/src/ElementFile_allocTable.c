
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/*   Finley: ElementFile */

/*   allocates and deallocates element table                  */

/**************************************************************/

#include "ElementFile.h"

/****************************************************************************/

/*  allocates the element table within an element file to hold numElements: */

void Finley_ElementFile_allocTable(Finley_ElementFile* in,dim_t numElements) 
{
  index_t *Id2=NULL,*Nodes2=NULL,*Tag2=NULL,*Color2=NULL;
  Paso_MPI_rank *Owner2=NULL;
  dim_t numNodes,e,i;

  Finley_resetError();
  /*  allocate memory: */ 
  numNodes=in->numNodes;
  Owner2=MEMALLOC(numElements,Paso_MPI_rank);
  Id2=MEMALLOC(numElements,index_t);
  Nodes2=MEMALLOC(numElements*in->numNodes,index_t);
  Tag2=MEMALLOC(numElements,index_t);
  Color2=MEMALLOC(numElements,index_t);
  
  /*  if fine, deallocate the old table and replace by new: */
  
  if ( Finley_checkPtr(Owner2) || Finley_checkPtr(Id2) || Finley_checkPtr(Nodes2) || 
                                                       Finley_checkPtr(Tag2) || Finley_checkPtr(Color2) ) {
    MEMFREE(Owner2);	
    MEMFREE(Nodes2);
    MEMFREE(Id2);
    MEMFREE(Tag2);
    MEMFREE(Color2);
  } else { 
    Finley_ElementFile_freeTable(in);
    in->Owner=Owner2;
    in->numElements=numElements;
    in->Id=Id2;
    in->Nodes=Nodes2;
    in->Tag=Tag2;
    in->Color=Color2;

    /* this initialization makes sure that data are located on the right processor */

    #pragma omp parallel for private(e,i) schedule(static)
    for (e=0;e<numElements;e++) {
       for (i=0;i<numNodes;i++) in->Nodes[INDEX2(i,e,numNodes)]=-1;
       in->Owner[e]=-1;
       in->Id[e]=-1;
       in->Tag[e]=-1;
       in->Color[e]=-1;
    }
    in->maxColor=-1;
    in->minColor=0;
  }
  return;
}

void Finley_ElementFile_setTagsInUse(Finley_ElementFile* in)
{
    index_t *tagsInUse=NULL;
    dim_t numTagsInUse;
    if (in !=NULL) {
       Finley_Util_setValuesInUse(in->Tag, in->numElements, &numTagsInUse, &tagsInUse, in->MPIInfo);
       if (Finley_noError()) {
          MEMFREE(in->tagsInUse);
          in->tagsInUse=tagsInUse;
          in->numTagsInUse=numTagsInUse;
       }
    }
}

/*  deallocates the element table within an element file: */

void Finley_ElementFile_freeTable(Finley_ElementFile* in) {
  MEMFREE(in->Owner);
  MEMFREE(in->Id);
  MEMFREE(in->Nodes);
  MEMFREE(in->Tag);
  MEMFREE(in->Color);
  MEMFREE(in->tagsInUse);
  in->numTagsInUse=0;
  in->numElements=0;
  in->maxColor=-1;
  in->minColor=0;
}
