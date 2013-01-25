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
  Paso_SystemMatrixPattern *pattern = NULL;
  dim_t i,n, numHops;
  index_t s,itmp,*rowLabel=NULL,*colLabel=NULL, *ptr=NULL, *index=NULL, *hop=NULL;
  Finley_IndexList* index_list=NULL;
  Finley_resetError();
  Finley_NodeDistribution *row_degreeOfFreedomDistribution=NULL;
  Finley_NodeDistribution *col_degreeOfFreedomDistribution=NULL;
  Paso_Distribution *row_distribution=NULL;
  Paso_Distribution *col_distribution=NULL;
#ifndef PASO_MPI
  index_t tmp[2];
  Paso_MPIInfo *temp_mpi_info = Paso_MPIInfo_alloc( MPI_COMM_WORLD );
#endif

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
        Finley_IndexList_insertElements(index_list,mesh,mesh->Elements,
                                        reduce_row_order,rowLabel,reduce_col_order,colLabel);
        Finley_IndexList_insertElements(index_list,mesh,mesh->FaceElements,
                                        reduce_row_order,rowLabel,reduce_col_order,colLabel);
        Finley_IndexList_insertElements(index_list,mesh,mesh->ContactElements,
                                        reduce_row_order,rowLabel,reduce_col_order,colLabel);
        Finley_IndexList_insertElements(index_list,mesh,mesh->Points,
                                        reduce_row_order,rowLabel,reduce_col_order,colLabel);
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
  if (Finley_noError()) {
    // Use a getReference method to get the DOF distributions to avoid memory leaks
#ifndef PASO_MPI
    /* Create row and column distributions for serial matrix, i.e. not distributed */
    tmp[0] = 0;
    if (reduce_col_order) {
      tmp[1] = mesh->Nodes->reducedNumDegreesOfFreedom;
    } else {
      tmp[1] = mesh->Nodes->numDegreesOfFreedom;
    }
    col_distribution=Paso_Distribution_alloc(temp_mpi_info,tmp,1,0);
    if (reduce_row_order) {
      tmp[1] = mesh->Nodes->reducedNumDegreesOfFreedom;
    } else {
      tmp[1] = mesh->Nodes->numDegreesOfFreedom;
    }
    row_distribution=Paso_Distribution_alloc(temp_mpi_info,tmp,1,0);
    /* Allocate MPIInfo for serial job */
    Paso_MPIInfo_dealloc(temp_mpi_info);
#else
    if (reduce_col_order) {
      col_degreeOfFreedomDistribution = mesh->Nodes->reducedDegreeOfFreedomDistribution;
    } else {
      col_degreeOfFreedomDistribution = mesh->Nodes->degreeOfFreedomDistribution;
    }
    if (reduce_row_order) {
      row_degreeOfFreedomDistribution = mesh->Nodes->reducedDegreeOfFreedomDistribution;
    } else {
      row_degreeOfFreedomDistribution = mesh->Nodes->degreeOfFreedomDistribution;
    }
    row_distribution=Paso_Distribution_alloc(mesh->MPIInfo, row_degreeOfFreedomDistribution->vtxdist,1,0);
    col_distribution=Paso_Distribution_alloc(mesh->MPIInfo, col_degreeOfFreedomDistribution->vtxdist,1,0);
#endif
    
    if (Finley_noError()) {
       Paso_SystemMatrixPattern_makeHops(PATTERN_FORMAT_DEFAULT,col_distribution,ptr,index,&numHops,&hop);
       if (Finley_noError()) {
          pattern = Paso_SystemMatrixPattern_alloc(PATTERN_FORMAT_DEFAULT,
                                                   row_distribution,
                                                   col_distribution,
                                                   ptr,index,numHops, hop); 
          if (Finley_noError()) {
            #ifdef PASO_MPI
            pattern->output_node_distribution=Finley_NodeDistribution_getReference(row_degreeOfFreedomDistribution);
            pattern->input_node_distribution=Finley_NodeDistribution_getReference(col_degreeOfFreedomDistribution);
            #endif
          }
      }
    }
    Paso_Distribution_free(row_distribution);
    Paso_Distribution_free(col_distribution);
  }
  MEMFREE(hop);
  if (! Finley_noError()) {
    MEMFREE(ptr);
    MEMFREE(index);
  }
  return pattern;
}
