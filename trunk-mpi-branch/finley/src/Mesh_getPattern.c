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

/* Finley: Mesh */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

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
            out=Paso_SystemMatrixPattern_reference(mesh->ReducedReducedPattern);
         } else {
            out=Paso_SystemMatrixPattern_reference(mesh->ReducedFullPattern);
         }
      } else {
         if (reduce_col_order) {
            out=Paso_SystemMatrixPattern_reference(mesh->FullReducedPattern);
         } else {
            out=Paso_SystemMatrixPattern_reference(mesh->FullFullPattern);
         }
      }
   }  
   return out;
}
Paso_SystemMatrixPattern* Finley_makePattern(Finley_Mesh *mesh,bool_t reduce_row_order, bool_t reduce_col_order) {
  double time0;
  dim_t i,n;
  index_t s,itmp,*rowLabel=NULL,*colLabel=NULL, *ptr=NULL, *index=NULL;
  Finley_IndexList* index_list=NULL;
  Finley_resetError();
  Finley_NodeDistribution *row_degreeOfFreedomDistribution;
  Finley_NodeDistribution *col_degreeOfFreedomDistribution;

  time0=Finley_timer();

#ifdef PASO_MPI
  if (reduce_col_order) {
       n=mesh->Nodes->reducedDegreeOfFreedomDistribution->numLocal;
       colLabel=mesh->Nodes->reducedDegreeOfFreedom;
  } else {
       n=mesh->Nodes->degreeOfFreedomDistribution->numLocal;
       colLabel=mesh->Nodes->degreeOfFreedom;
  }
     
  if (reduce_row_order) {
      n=mesh->Nodes->reducedDegreeOfFreedomDistribution->numLocal;
      rowLabel=mesh->Nodes->reducedDegreeOfFreedom;
  } else {
      n=mesh->Nodes->degreeOfFreedomDistribution->numLocal;
      rowLabel=mesh->Nodes->degreeOfFreedom;
  }
#else
  if (reduce_col_order) {
       n=mesh->Nodes->reducedNumDegreesOfFreedom;
       colLabel=mesh->Nodes->reducedDegreeOfFreedom;
  } else {
       n=mesh->Nodes->numDegreesOfFreedom;
       colLabel=mesh->Nodes->degreeOfFreedom;
  }
     
  if (reduce_row_order) {
      n=mesh->Nodes->reducedNumDegreesOfFreedom;
      rowLabel=mesh->Nodes->reducedDegreeOfFreedom;
  } else {
      n=mesh->Nodes->numDegreesOfFreedom;
      rowLabel=mesh->Nodes->degreeOfFreedom;
  }
#endif

  index_list=TMPMEMALLOC(n,Finley_IndexList);
  ptr=MEMALLOC(n+1,index_t);
  if (! (Finley_checkPtr(index_list) || Finley_checkPtr(ptr)) ) {
      #pragma omp parallel private(i,s,itmp)
      {
        #pragma omp for schedule(static)
        for(i=0;i<n;++i) {
             index_list[i].extension=NULL;
             index_list[i].n=0;
        }
        /*  insert contributions from element matrices into colums index index_list: */
	printf("ksteube mesh->Elements\n");
        Finley_IndexList_insertElements(index_list,mesh,mesh->Elements,
                                        reduce_row_order,rowLabel,reduce_col_order,colLabel);
	printf("ksteube mesh->FaceElements\n");
        Finley_IndexList_insertElements(index_list,mesh,mesh->FaceElements,
                                        reduce_row_order,rowLabel,reduce_col_order,colLabel);
	printf("ksteube mesh->ContactElements\n");
        Finley_IndexList_insertElements(index_list,mesh,mesh->ContactElements,
                                        reduce_row_order,rowLabel,reduce_col_order,colLabel);
	printf("ksteube mesh->Points\n");
        Finley_IndexList_insertElements(index_list,mesh,mesh->Points,
                                        reduce_row_order,rowLabel,reduce_col_order,colLabel);
	printf("ksteube done with 4 calls to Finley_IndexList_insertElements\n");
        /* get the number of connections per row */
        #pragma omp for schedule(static)
        for(i=0;i<n;++i) {
            ptr[i]=Finley_IndexList_count(&index_list[i]);
        }
        /* accumalate ptr */
        /* OMP */
        #pragma omp master 
        {
            s=0;
            for(i=0;i<n;++i) {
                 itmp=ptr[i];
                 ptr[i]=s;
                 s+=itmp;
            }
            ptr[n]=s;
            /* fill index */
            index=MEMALLOC(ptr[n],index_t);
        }
        #pragma omp barrier
        if (! Finley_checkPtr(index)) {
            #pragma omp for schedule(static)
            for(i=0;i<n;++i) Finley_IndexList_toArray(&index_list[i],&index[ptr[i]]);
        }
     }
  }
  /* all done */
  /* clean up */
  if (index_list!=NULL) {
    #pragma omp parallel for private(i) 
    for(i=0;i<n;++i) Finley_IndexList_free(index_list[i].extension);
  }
  TMPMEMFREE(index_list);
  #ifdef Finley_TRACE
  printf("timing: mesh to matrix pattern: %.4e sec\n",Finley_timer()-time0);
  #endif
  if (! Finley_noError()) {
    MEMFREE(ptr);
    MEMFREE(index);
    return NULL;
  } else {
    // Use a getReference method to get the DOF distributions to avoid memory leaks
    if (reduce_col_order) {
      col_degreeOfFreedomDistribution = mesh->Nodes->reducedDegreeOfFreedomDistribution;
    }
    else {
      col_degreeOfFreedomDistribution = mesh->Nodes->degreeOfFreedomDistribution;
    }
    if (reduce_row_order) {
      row_degreeOfFreedomDistribution = mesh->Nodes->reducedDegreeOfFreedomDistribution;
    }
    else {
      row_degreeOfFreedomDistribution = mesh->Nodes->degreeOfFreedomDistribution;
    }
    Paso_Distribution *row_distribution=Paso_Distribution_alloc(mesh->MPIInfo, row_degreeOfFreedomDistribution->vtxdist,1,0);
    Paso_Distribution *col_distribution=Paso_Distribution_alloc(mesh->MPIInfo, col_degreeOfFreedomDistribution->vtxdist,1,0);
    Paso_SystemMatrixPattern *pattern = Paso_SystemMatrixPattern_alloc(PATTERN_FORMAT_DEFAULT,
                                  row_distribution,
                                  col_distribution,
                                  ptr,index); /* add distr info */
    return pattern;
  }
}
