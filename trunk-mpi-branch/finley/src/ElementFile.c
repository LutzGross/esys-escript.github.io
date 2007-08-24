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

/*   Finley: ElementFile */

/*   allocates an element file to hold elements of type id and with integration order order. */
/*   use Finley_Mesh_allocElementTable to allocate the element table (Id,Nodes,Tag,Owner). */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "ElementFile.h"

/**************************************************************/

Finley_ElementFile* Finley_ElementFile_alloc(ElementTypeId id, index_t order, index_t reduced_order, Paso_MPIInfo *MPIInfo)
{
  extern Finley_RefElementInfo Finley_RefElement_InfoList[];
  dim_t NQ, reduced_NQ;
  Finley_ElementFile *out;
  
  /*   get the number of quadrature nodes needed to achieve integration order order: */
  
  if (order<0) order=MAX(2*Finley_RefElement_InfoList[id].numOrder,0);
  if (reduced_order<0) reduced_order=MAX(2*(Finley_RefElement_InfoList[id].numOrder-1),0);
  NQ= Finley_RefElement_InfoList[id].getNumQuadNodes(order);
  reduced_NQ= Finley_RefElement_InfoList[id].getNumQuadNodes(reduced_order);
  if (! Finley_noError()) return NULL;
  
  /*  allocate the return value */
  
  out=MEMALLOC(1,Finley_ElementFile);
  if (Finley_checkPtr(out)) return NULL;
  out->order = order;
  out->reduced_order = reduced_order;
  out->ReferenceElement=NULL;
  out->LinearReferenceElement=NULL;
  out->ReferenceElementReducedOrder=NULL;
  out->LinearReferenceElementReducedOrder=NULL;
  out->numElements=0;
  out->Id=NULL;
  out->Nodes=NULL;
  out->Tag=NULL;
  out->Color=NULL;
  out->minColor=0;
  out->maxColor=-1;
  out->jacobeans=NULL;
  out->jacobeans_reducedQ=NULL;
  out->jacobeans_reducedS=NULL;
  out->jacobeans_reducedS_reducedQ=NULL;

  out->Owner=NULL;                
  out->MPIInfo = Paso_MPIInfo_getReference( MPIInfo );

  /*  allocate the reference element: */
  
  out->ReferenceElement=Finley_RefElement_alloc(id,NQ);
  out->jacobeans=Finley_ElementFile_Jacobeans_alloc(out->ReferenceElement);
  out->ReferenceElementReducedOrder=Finley_RefElement_alloc(id,reduced_NQ);
  out->jacobeans_reducedQ=Finley_ElementFile_Jacobeans_alloc(out->ReferenceElementReducedOrder);

  out->LinearReferenceElement=Finley_RefElement_alloc(Finley_RefElement_InfoList[id].LinearTypeId,NQ);
  out->jacobeans_reducedS=Finley_ElementFile_Jacobeans_alloc(out->LinearReferenceElement);
  out->LinearReferenceElementReducedOrder=Finley_RefElement_alloc(Finley_RefElement_InfoList[id].LinearTypeId,reduced_NQ);
  out->jacobeans_reducedS_reducedQ=Finley_ElementFile_Jacobeans_alloc(out->LinearReferenceElementReducedOrder);

  out->numNodes=out->ReferenceElement->Type->numNodes;

  if (! Finley_noError()) {
     Finley_ElementFile_free(out);
     return NULL;
  }
  return out;
}

/*  deallocates an element file: */

void Finley_ElementFile_free(Finley_ElementFile* in) {
  if (in!=NULL) {
     #ifdef Finley_TRACE
     if (in->ReferenceElement!=NULL) printf("element file for %s is deallocated.\n",in->ReferenceElement->Type->Name);
     #endif
     Finley_ElementFile_freeTable(in);   
     Finley_RefElement_dealloc(in->ReferenceElement);
     Finley_RefElement_dealloc(in->ReferenceElementReducedOrder);
     Finley_RefElement_dealloc(in->LinearReferenceElement);
     Finley_RefElement_dealloc(in->LinearReferenceElementReducedOrder);
     Finley_ElementFile_Jacobeans_dealloc(in->jacobeans);
     Finley_ElementFile_Jacobeans_dealloc(in->jacobeans_reducedS);
     Finley_ElementFile_Jacobeans_dealloc(in->jacobeans_reducedQ);
     Finley_ElementFile_Jacobeans_dealloc(in->jacobeans_reducedS_reducedQ);
     Paso_MPIInfo_free( in->MPIInfo );
     MEMFREE(in);      
  }
}
void Finley_ElementFile_setElementDistribution(Finley_ElementFile* in, dim_t* distribution) {
  dim_t local_num_elements,e,out,num_elements=0, size;
  Paso_MPI_rank myRank;
  if (in == NULL) {
      distribution[0]=num_elements;
  } else {
      if (in->MPIInfo->size>1) {
         num_elements=0;
         myRank=in->MPIInfo->rank;
         size=in->MPIInfo->size;
         #pragma omp parallel private(local_num_elements)
         {
            local_num_elements=0;
            #pragma omp for private(e)
            for (e=0;e<in->numElements;e++) {
               if (in->Owner[e] == myRank) local_num_elements++;
            }
            #pragma omp critical
            num_elements+=local_num_elements;
         }
         #ifdef PASO_MPI
           MPI_Allgather(&num_elements,1,MPI_INT,distribution,1,MPI_INT,in->MPIInfo->comm);
         #else
           distribution[0]=num_elements;
         #endif
      } else {
        distribution[0]=in->numElements;
      }
  }
}

dim_t Finley_ElementFile_getGlobalNumElements(Finley_ElementFile* in) {
  dim_t size, *distribution=NULL, out, p;
  if (in == NULL) {
      return 0;
  } else {
    size=in->MPIInfo->size;
    distribution=TMPMEMALLOC(size,dim_t);
    Finley_ElementFile_setElementDistribution(in,distribution);
    out=0;
    for (p=0;p<size;++p) out+=distribution[p];
    TMPMEMFREE(distribution);
    return out;
  }
}
dim_t Finley_ElementFile_getMyNumElements(Finley_ElementFile* in) {
  dim_t size, *distribution=NULL, out, p;
  if (in == NULL) {
      return 0;
  } else {
    size=in->MPIInfo->size;
    distribution=TMPMEMALLOC(size,dim_t);
    Finley_ElementFile_setElementDistribution(in,distribution);
    out=distribution[in->MPIInfo->rank];
    TMPMEMFREE(distribution);
    return out;
  }

}
index_t Finley_ElementFile_getFirstElement(Finley_ElementFile* in){
  dim_t size, *distribution=NULL, out, p;
  if (in == NULL) {
      return 0;
  } else {
    size=in->MPIInfo->size;
    distribution=TMPMEMALLOC(size,dim_t);
    Finley_ElementFile_setElementDistribution(in,distribution);
    out=0;
    for (p=0;p<in->MPIInfo->rank;++p) out+=distribution[p];
    TMPMEMFREE(distribution);
    return out;
  }
}
