
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Finley: Mesh */

/**************************************************************/

#include "Mesh.h"
#include "IndexList.h"

/**************************************************************/

/* returns a reference to the matrix pattern                  */

Paso_SystemMatrixPattern* Finley_getPattern(Finley_Mesh *mesh,bool_t reduce_row_order, bool_t reduce_col_order) {
   Paso_SystemMatrixPattern *out=NULL;
   Finley_resetError();
   /* make sure that the requested pattern is available */
   if (reduce_row_order) {
      if (reduce_col_order) {
         if (mesh->ReducedReducedPattern==NULL) mesh->ReducedReducedPattern=Finley_makePattern(mesh,reduce_row_order,reduce_col_order);
      } else {
         if (mesh->ReducedFullPattern==NULL) mesh->ReducedFullPattern=Finley_makePattern(mesh,reduce_row_order,reduce_col_order);
      }
   } else {
      if (reduce_col_order) {
         if (mesh->FullReducedPattern==NULL) mesh->FullReducedPattern=Finley_makePattern(mesh,reduce_row_order,reduce_col_order);
      } else {
         if (mesh->FullFullPattern==NULL) mesh->FullFullPattern=Finley_makePattern(mesh,reduce_row_order,reduce_col_order);
      }
   }
   if (Finley_noError()) {
      if (reduce_row_order) {
         if (reduce_col_order) {
            out=Paso_SystemMatrixPattern_getReference(mesh->ReducedReducedPattern);
         } else {
            out=Paso_SystemMatrixPattern_getReference(mesh->ReducedFullPattern);
         }
      } else {
         if (reduce_col_order) {
            out=Paso_SystemMatrixPattern_getReference(mesh->FullReducedPattern);
         } else {
            out=Paso_SystemMatrixPattern_getReference(mesh->FullFullPattern);
         }
      }
   }  
   return out;
}
Paso_SystemMatrixPattern* Finley_makePattern(Finley_Mesh *mesh,bool_t reduce_row_order, bool_t reduce_col_order) {
  double time0;
  Paso_SystemMatrixPattern* out=NULL;
  Paso_Pattern *main_pattern = NULL, *col_couple_pattern=NULL, *row_couple_pattern=NULL;
  Paso_Connector *col_connector, *row_connector;
  Finley_IndexList* index_list=NULL;
  Finley_NodeMapping *colMap=NULL, *rowMap=NULL;
  Paso_Distribution *colDistribution=NULL, *rowDistribution=NULL;
  
  index_t i;
  Finley_resetError();
  time0=Finley_timer();

  if (reduce_col_order) {
       colMap=mesh->Nodes->reducedDegreesOfFreedomMapping;
       colDistribution=mesh->Nodes->reducedDegreesOfFreedomDistribution;
       col_connector=mesh->Nodes->reducedDegreesOfFreedomConnector;

  } else {
       colMap=mesh->Nodes->degreesOfFreedomMapping;
       colDistribution=mesh->Nodes->degreesOfFreedomDistribution;
       col_connector=mesh->Nodes->degreesOfFreedomConnector;
  }
     
  if (reduce_row_order) {
      rowMap=mesh->Nodes->reducedDegreesOfFreedomMapping;
      rowDistribution=mesh->Nodes->reducedDegreesOfFreedomDistribution;
      row_connector=mesh->Nodes->reducedDegreesOfFreedomConnector;
  } else {
      rowMap=mesh->Nodes->degreesOfFreedomMapping;
      rowDistribution=mesh->Nodes->degreesOfFreedomDistribution;
      row_connector=mesh->Nodes->degreesOfFreedomConnector;
  }

  index_list=TMPMEMALLOC(rowMap->numTargets,Finley_IndexList);
  if (! Finley_checkPtr(index_list)) {
  
      #pragma omp parallel private(i)
      {
        #pragma omp for schedule(static)
        for(i=0;i<rowMap->numTargets;++i) {
             index_list[i].extension=NULL;
             index_list[i].n=0;
        }
        /*  insert contributions from element matrices into colums index index_list: */
        Finley_IndexList_insertElements(index_list,mesh->Elements,
                                        reduce_row_order,rowMap->target,reduce_col_order,colMap->target);
        Finley_IndexList_insertElements(index_list,mesh->FaceElements,
                                        reduce_row_order,rowMap->target,reduce_col_order,colMap->target);
        Finley_IndexList_insertElements(index_list,mesh->ContactElements,
                                        reduce_row_order,rowMap->target,reduce_col_order,colMap->target);
        Finley_IndexList_insertElements(index_list,mesh->Points,
                                        reduce_row_order,rowMap->target,reduce_col_order,colMap->target);

     }
 
     /* create pattern */
     main_pattern=Finley_IndexList_createPattern(0,Paso_Distribution_getMyNumComponents(rowDistribution),index_list,
                                                 0,Paso_Distribution_getMyNumComponents(colDistribution),
                                                 0);
     col_couple_pattern=Finley_IndexList_createPattern(0,Paso_Distribution_getMyNumComponents(rowDistribution),index_list,
                                                 Paso_Distribution_getMyNumComponents(colDistribution),colMap->numTargets,
                                                 -Paso_Distribution_getMyNumComponents(colDistribution));
     row_couple_pattern=Finley_IndexList_createPattern(Paso_Distribution_getMyNumComponents(rowDistribution),rowMap->numTargets,index_list,
                                                       0,Paso_Distribution_getMyNumComponents(colDistribution),
                                                       0);

     /* if everthing is in order we can create the return value */
     if (Finley_noError()) {
          out=Paso_SystemMatrixPattern_alloc(PATTERN_FORMAT_DEFAULT,
                                             rowDistribution,
                                             colDistribution,
                                             main_pattern,
                                             col_couple_pattern,
                                             row_couple_pattern,
                                             col_connector,
                                             row_connector);
     }
     /* clean up */
     if (index_list!=NULL) {
        #pragma omp parallel for private(i) 
        for(i=0;i<rowMap->numTargets;++i) Finley_IndexList_free(index_list[i].extension);
     }
     TMPMEMFREE(index_list);
     Paso_Pattern_free(main_pattern);
     Paso_Pattern_free(col_couple_pattern);
     Paso_Pattern_free(row_couple_pattern);
  }
  #ifdef Finley_TRACE
  printf("timing: mesh to matrix pattern: %.4e sec\n",Finley_timer()-time0);
  #endif
  Paso_MPIInfo_noError(mesh->MPIInfo);
  return out;
}
