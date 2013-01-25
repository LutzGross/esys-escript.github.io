
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
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
           Paso_SystemMatrix_free(in->mass_matrix);
           Paso_SystemMatrix_free(in->iteration_matrix);
           Paso_MPIInfo_free(in->mpi_info);
           Paso_Coupler_free(in->u_coupler);
           MEMFREE(in->u);
           MEMFREE(in->constraint_weights);
           MEMFREE(in->main_iptr);
           MEMFREE(in->lumped_mass_matrix);
           MEMFREE(in->main_diagonal_low_order_transport_matrix);
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
Paso_SystemMatrix* Paso_FCTransportProblem_borrowMassMatrix(Paso_FCTransportProblem* in) {
   return in->mass_matrix;
}

double* Paso_FCTransportProblem_borrowLumpedMassMatrix(Paso_FCTransportProblem* in) {
    return in->lumped_mass_matrix;
}

dim_t Paso_FCTransportProblem_getTotalNumRows(Paso_FCTransportProblem* in) {
    return Paso_SystemMatrix_getTotalNumRows(in->transport_matrix);
}

Paso_FCTransportProblem* Paso_FCTransportProblem_alloc(double theta, Paso_SystemMatrixPattern *pattern, int block_size) 
{
     Paso_SystemMatrixType matrix_type=MATRIX_FORMAT_DEFAULT+MATRIX_FORMAT_BLK1;  /* at the moment only block size 1 is supported */
     Paso_FCTransportProblem* out=NULL;
     Paso_SystemMatrixPattern *transport_pattern;
     dim_t n,i;
     index_t iptr,iptr_main;

     if ((theta<0.) || (theta >1.)) {
        Paso_setError(TYPE_ERROR,"Paso_FCTransportProblem_alloc: theta needs to be between 0. and. 1.");
        return NULL;
     }

     out=MEMALLOC(1,Paso_FCTransportProblem);
     if (Paso_checkPtr(out)) return NULL;
     out->reference_counter=0;
     out->theta=theta;
     out->dt_max=LARGE_POSITIVE_FLOAT;
     out->constraint_factor=sqrt(LARGE_POSITIVE_FLOAT);
     if (out->theta < 1.) {
            out->dt_factor=MIN(1./(1.-out->theta),DT_FACTOR_MAX);
     } else {
            out->dt_factor=DT_FACTOR_MAX;
     }
     out->valid_matrices=FALSE;
     out->transport_matrix=Paso_SystemMatrix_alloc(matrix_type,pattern,block_size,block_size);
     out->mass_matrix=Paso_SystemMatrix_alloc(matrix_type,pattern,block_size,block_size);
     out->iteration_matrix=NULL;
     out->u=NULL;
     out->constraint_weights=NULL;
     out->mpi_info=Paso_MPIInfo_getReference(pattern->mpi_info);
     out->main_iptr=NULL;
     out->lumped_mass_matrix=NULL;
     out->main_diagonal_low_order_transport_matrix=NULL;

     if (Paso_noError()) {
         n=Paso_SystemMatrix_getTotalNumRows(out->transport_matrix);
         transport_pattern=out->transport_matrix->pattern;
         out->u=MEMALLOC(n,double);
         out->constraint_weights=MEMALLOC(n,double);
         out->main_iptr=MEMALLOC(n,index_t);
         out->lumped_mass_matrix=MEMALLOC(n,double);
         out->main_diagonal_low_order_transport_matrix=MEMALLOC(n,double);
         out->u_coupler=Paso_Coupler_alloc(Paso_FCTransportProblem_borrowConnector(out),block_size);

         if ( ! (Paso_checkPtr(out->u) || Paso_checkPtr(out->main_iptr) || Paso_checkPtr(out->constraint_weights) || 
                 Paso_checkPtr(out->lumped_mass_matrix) || Paso_checkPtr(out->main_diagonal_low_order_transport_matrix)) && Paso_noError()  ) {
             
             #pragma omp parallel 
             {
                 #pragma omp for schedule(static) private(i)
                 for (i = 0; i < n; ++i) {
                    out->lumped_mass_matrix[i]=0.;
                    out->main_diagonal_low_order_transport_matrix[i]=0.;
                    out->u[i]=0.;
                 }
                 /* identify the main diagonals */
                 #pragma omp for schedule(static) private(i,iptr,iptr_main)
                 for (i = 0; i < n; ++i) {
                        iptr_main=transport_pattern->mainPattern->ptr[0]-1;
                        for (iptr=transport_pattern->mainPattern->ptr[i];iptr<transport_pattern->mainPattern->ptr[i+1]; iptr++) {
                              if (transport_pattern->mainPattern->index[iptr]==i) {
                                   iptr_main=iptr;
                                   break;
                              }
                        }
                        out->main_iptr[i]=iptr_main;
                        if (iptr_main==transport_pattern->mainPattern->ptr[0]-1)
                             Paso_setError(VALUE_ERROR, "Paso_FCTransportProblem_alloc: no main diagonal");
                 }
     
             }
      }

  }
  if (Paso_noError()) {
     out->reference_counter=1;
     return out;
  } else {
     Paso_FCTransportProblem_free(out);
     return NULL;
  }
} 

void Paso_FCTransportProblem_checkinSolution(Paso_FCTransportProblem* in, double* u) {
    dim_t i, n;
    double local_u_min,u_min;
    n=Paso_SystemMatrix_getTotalNumRows(in->transport_matrix);
    u_min=LARGE_POSITIVE_FLOAT;
    #pragma omp parallel private(local_u_min)
    {
         local_u_min=0.;
         #pragma omp for schedule(static) private(i)
         for (i=0;i<n;++i) {
              in->u[i]=u[i];
              local_u_min=MIN(local_u_min,u[i]);
         }
         #pragma omp critical 
         {
            u_min=MIN(u_min,local_u_min);
         }
    }
    #ifdef PASO_MPI
         local_u_min=u_min;
         MPI_Allreduce(&local_u_min,&u_min, 1, MPI_DOUBLE, MPI_MIN, in->mpi_info->comm);
    #endif
    if (u_min<0) {
       Paso_setError(VALUE_ERROR, "Paso_FCTransportProblem_checkinSolution: initial guess must be non-negative.");
    }
}
dim_t Paso_FCTransportProblem_getBlockSize(const Paso_FCTransportProblem* in)
{
   return in->transport_matrix->row_block_size;
}

Paso_Connector* Paso_FCTransportProblem_borrowConnector(const Paso_FCTransportProblem* in)
{
   return in->transport_matrix->pattern->col_connector;
}

index_t Paso_FCTransportProblem_getTypeId(const index_t solver,const index_t preconditioner, const index_t package,const  bool_t symmetry) 
{
   return MATRIX_FORMAT_DEFAULT + MATRIX_FORMAT_BLK1;
}

void Paso_FCTransportProblem_setUpConstraint(Paso_FCTransportProblem* fctp,  const double* q, const double factor)
{
   dim_t i, n;
   register double m, rtmp;
   double factor2= fctp->dt_factor * factor;
   n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);

   if ( fctp->valid_matrices ) {
      Paso_setError(VALUE_ERROR, "Paso_FCTransportProblem_insertConstraint: you must not insert a constraint is a valid system.");
      return;
   }
   if (factor<=0) {
      Paso_setError(VALUE_ERROR, "Paso_FCTransportProblem_insertConstraint: constraint_factor needs to be positive.");
      return;
   }

   
   #pragma omp for schedule(static) private(i,m,rtmp)
   for (i=0;i<n;++i) {
        if (q[i]>0) {
           m=fctp->mass_matrix->mainBlock->val[fctp->main_iptr[i]];
           rtmp=factor2 * (m == 0 ? 1 : m);
           fctp->constraint_weights[i]=rtmp;
           fctp->mass_matrix->mainBlock->val[fctp->main_iptr[i]]=m+rtmp;
        } else {
           fctp->constraint_weights[i]=0;
        }
   }
   fctp->constraint_factor=factor;
}

void Paso_FCTransportProblem_insertConstraint(Paso_FCTransportProblem* fctp,  const double* r,  double* source)
{
   dim_t i, n;
   n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);

   #pragma omp for schedule(static) private(i)
   for (i=0;i<n;++i) source[i]+=fctp->constraint_weights[i] * r[i];
}
