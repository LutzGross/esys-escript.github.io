/* $Id:$ */

/*******************************************************
 *
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/* Paso: FCTransportProblem                                          */

/**************************************************************/

/* Author: l.gross@uq.edu.au                                */

/**************************************************************/


#include "Paso.h"
#include "SolverFCT.h"
#include "PasoUtil.h"


/**************************************************************/

/* free all memory used by                                */
void Paso_FCTransportProblem_free(Paso_FCTransportProblem* in) {
     if (in!=NULL) {
        in->reference_counter--;
        if (in->reference_counter<=0) {

           Paso_SystemMatrix_free(in->transport_matrix);
           Paso_SystemMatrix_free(in->flux_matrix);
           Paso_MPIInfo_free(in->mpi_info);

           MEMFREE(in->u);
           MEMFREE(in->lumped_mass_matrix);
           MEMFREE(in->row_sum_flux_matrix);
           MEMFREE(in->transport_matrix_diagonal);
           MEMFREE(in->r_p);
           MEMFREE(in->r_n);
           MEMFREE(in->main_iptr);
           MEMFREE(in);
        }
    }
}

Paso_FCTransportProblem* Paso_FCTransportProblem_getReference(Paso_FCTransportProblem* in) {
     if (in!=NULL) {
        ++(in->reference_counter);
     }
     return in;
}    

Paso_SystemMatrix* Paso_FCTransportProblem_borrowTransportMatrix(Paso_FCTransportProblem* in) {
   return in->transport_matrix;
}

Paso_SystemMatrix* Paso_FCTransportProblem_borrowFluxMatrix(Paso_FCTransportProblem* in) {
    return in->flux_matrix;
}

double* Paso_FCTransportProblem_borrowLumpedMassMatrix(Paso_FCTransportProblem* in) {
    return in->lumped_mass_matrix;
}

dim_t Paso_FCTransportProblem_getTotalNumRows(Paso_FCTransportProblem* in) {
    return Paso_SystemMatrix_getTotalNumRows(in->transport_matrix);
}

Paso_FCTransportProblem* Paso_FCTransportProblem_alloc(double theta, double dt_max, Paso_SystemMatrixPattern *pattern, int block_size


) {
     Paso_SystemMatrixType matrix_type=MATRIX_FORMAT_DEFAULT+MATRIX_FORMAT_BLK1;  /* at the moment only block size 1 is supported */
     Paso_FCTransportProblem* out=NULL;
     dim_t n,i;
     index_t iptr,iptr_main,k;

     if ((theta<0.) || (theta >1.)) {
        Paso_setError(TYPE_ERROR,"Paso_FCTransportProblem_alloc: theta needs to be between 0. and. 1.");
        return NULL;
     }

     out=MEMALLOC(1,Paso_FCTransportProblem);
     if (Paso_checkPtr(out)) return NULL;

     out->theta=theta;
     out->dt_max=dt_max;
     out->valid_matrices=FALSE;
     out->transport_matrix=Paso_SystemMatrix_alloc(matrix_type,pattern,block_size,block_size);
     Paso_SystemMatrix_allocBuffer(out->transport_matrix);
     out->flux_matrix=Paso_SystemMatrix_alloc(matrix_type,pattern,block_size,block_size);
     out->mpi_info=Paso_MPIInfo_getReference(pattern->mpi_info);

     out->main_iptr=NULL;
     out->lumped_mass_matrix=NULL;
     out->row_sum_flux_matrix=NULL;
     out->transport_matrix_diagonal=NULL;
     out->r_p=NULL;
     out->r_n=NULL;

     if (Paso_noError()) {
         n=Paso_SystemMatrix_getTotalNumRows(out->transport_matrix);

         out->r_p=MEMALLOC(n,double);
         out->r_n=MEMALLOC(n,double);
         out->main_iptr=MEMALLOC(n,index_t);
         out->lumped_mass_matrix=MEMALLOC(n,double);
         out->row_sum_flux_matrix=MEMALLOC(n,double);
         out->transport_matrix_diagonal=MEMALLOC(n,double);
         out->u=MEMALLOC(n,double);

         if ( ! (Paso_checkPtr(out->r_p) || Paso_checkPtr(out->r_n) || Paso_checkPtr(out->main_iptr) || 
                 Paso_checkPtr(out->lumped_mass_matrix) || Paso_checkPtr(out->transport_matrix_diagonal) || Paso_checkPtr(out->row_sum_flux_matrix) || Paso_checkPtr(out->u)) ) {
             
             #pragma omp parallel for schedule(static) private(i)
             for (i = 0; i < n; ++i) {
                out->lumped_mass_matrix[i]=0.;
                out->row_sum_flux_matrix[i]=0.;
                out->u[i]=0.;
                out->r_p[i]=0.;
                out->r_n[i]=0.;
             }
             /* identify the main diagonals */
             #pragma omp parallel for schedule(static) private(i,iptr,iptr_main,k)
             for (i = 0; i < n; ++i) {
                    iptr_main=pattern->mainPattern->ptr[0]-1;
                    for (iptr=pattern->mainPattern->ptr[i];iptr<pattern->mainPattern->ptr[i+1]; iptr++) {
                          if (pattern->mainPattern->index[iptr]==i) {
                               iptr_main=iptr;
                               break;
                          }
                    }
                    out->main_iptr[i]=iptr_main;
                    if (iptr_main==pattern->mainPattern->ptr[0]-1)
                         Paso_setError(VALUE_ERROR, "Paso_FCTransportProblem_alloc: no main diagonal");
             }

      }

  }
  if (Paso_noError()) {
     return out;
  } else {
     Paso_FCTransportProblem_free(out);
     return NULL;
  }
} 

void Paso_FCTransportProblem_checkinSolution(Paso_FCTransportProblem* in, double* u) {
    dim_t i, n;
   
    n=Paso_FCTransportProblem_getTotalNumRows(in);
    #pragma omp parallel for schedule(static) private(i)
    for (i = 0; i < n; ++i) {
         in->u[i]=u[i];
    }
}
