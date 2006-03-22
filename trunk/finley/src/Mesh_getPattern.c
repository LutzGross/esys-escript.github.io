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
  
  time0=Finley_timer();
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
        Finley_IndexList_insertElements(index_list,mesh->Elements,
                                        reduce_row_order,rowLabel,reduce_col_order,colLabel);
        Finley_IndexList_insertElements(index_list,mesh->FaceElements,
                                        reduce_row_order,rowLabel,reduce_col_order,colLabel);
        Finley_IndexList_insertElements(index_list,mesh->ContactElements,
                                        reduce_row_order,rowLabel,reduce_col_order,colLabel);
        Finley_IndexList_insertElements(index_list,mesh->Points,
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
  if (! Finley_noError()) {
    MEMFREE(ptr);
    MEMFREE(index);
    return NULL;
  } else {
    return Paso_SystemMatrixPattern_alloc(MATRIX_FORMAT_DEFAULT,n,ptr,index);
  }
}
/*
 * $Log$
 * Revision 1.5  2005/09/15 03:44:22  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.4  2005/09/01 03:31:35  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-01
 *
 * Revision 1.3.2.2  2005/09/07 06:26:19  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.3.2.1  2005/08/24 02:02:18  gross
 * timing output switched off. solver output can be swiched through getSolution(verbose=True) now.
 *
 * Revision 1.3  2005/07/08 04:07:52  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.2  2004/12/15 07:08:33  jgs
 * *** empty log message ***
 *
 * Revision 1.1.2.5  2005/06/30 01:53:55  gross
 * a bug in coloring fixed
 *
 * Revision 1.1.2.4  2005/06/29 02:34:51  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.2.3  2004/11/24 01:37:14  gross
 * some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
 *
 * Revision 1.1.2.1  2004/11/15 01:01:09  gross
 * and anotther missing file
 *
 */
