
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: TransportProblem (see Paso_TransportSolver_solve)   */                                      

/**************************************************************/

/* Author: l.gross@uq.edu.au                                */

/**************************************************************/


#include "Transport.h"
#include "PasoUtil.h"


/**************************************************************/

/* free all memory used by                                */

void Paso_TransportProblem_free(Paso_TransportProblem* in) {
     if (in!=NULL) {
        in->reference_counter--;
        if (in->reference_counter<=0) {
           Paso_SystemMatrix_free(in->transport_matrix);
           Paso_SystemMatrix_free(in->mass_matrix);
           Paso_SystemMatrix_free(in->iteration_matrix);
           Paso_MPIInfo_free(in->mpi_info);
           Paso_Coupler_free(in->u_coupler);
           MEMFREE(in->constraint_weights);
           MEMFREE(in->reactive_matrix);
           MEMFREE(in->main_diagonal_mass_matrix);
           MEMFREE(in->lumped_mass_matrix);
           MEMFREE(in->main_diagonal_low_order_transport_matrix);
           MEMFREE(in);
        }
    }
}

Paso_TransportProblem* Paso_TransportProblem_getReference(Paso_TransportProblem* in) {
     if (in!=NULL) {
        ++(in->reference_counter);
     }
     return in;
}    

Paso_SystemMatrix* Paso_TransportProblem_borrowTransportMatrix(Paso_TransportProblem* in) {
   return in->transport_matrix;
}
Paso_SystemMatrix* Paso_TransportProblem_borrowMassMatrix(Paso_TransportProblem* in) {
   return in->mass_matrix;
}

double* Paso_TransportProblem_borrowLumpedMassMatrix(Paso_TransportProblem* in) {
    return in->lumped_mass_matrix;
}

dim_t Paso_TransportProblem_getTotalNumRows(Paso_TransportProblem* in) {
    return Paso_SystemMatrix_getTotalNumRows(in->transport_matrix);
}

Paso_TransportProblem* Paso_TransportProblem_alloc(bool_t useBackwardEuler, Paso_SystemMatrixPattern *pattern, int block_size) 
{
     Paso_SystemMatrixType matrix_type=MATRIX_FORMAT_DEFAULT+MATRIX_FORMAT_BLK1;  /* at the moment only block size 1 is supported */
     Paso_TransportProblem* out=NULL;
     Paso_SystemMatrixPattern *transport_pattern;
     dim_t n,i;

     out=MEMALLOC(1,Paso_TransportProblem);
     if (Paso_checkPtr(out)) return NULL;
     out->reference_counter=0;
     out->useBackwardEuler=useBackwardEuler;
     out->dt_max=LARGE_POSITIVE_FLOAT;
     out->dt_failed=LARGE_POSITIVE_FLOAT;
/****************** REVISE ****************************/
     out->constraint_factor=sqrt(LARGE_POSITIVE_FLOAT);
     if (out->useBackwardEuler) {
            out->dt_factor=DT_FACTOR_MAX;
     } else {
            out->dt_factor=2.;
     }
/*****************************************************/
     out->valid_matrices=FALSE;
     out->transport_matrix=Paso_SystemMatrix_alloc(matrix_type,pattern,block_size,block_size,FALSE);
     out->mass_matrix=Paso_SystemMatrix_alloc(matrix_type,pattern,block_size,block_size,FALSE);
     out->iteration_matrix=NULL;
     out->constraint_weights=NULL;
     out->mpi_info=Paso_MPIInfo_getReference(pattern->mpi_info);

     out->lumped_mass_matrix=NULL;
     out->main_diagonal_low_order_transport_matrix=NULL;
     out->reactive_matrix=NULL;
     out->main_diagonal_mass_matrix=NULL;

     if (Paso_noError()) {
         n=Paso_SystemMatrix_getTotalNumRows(out->transport_matrix);
         transport_pattern=out->transport_matrix->pattern;
         out->constraint_weights=MEMALLOC(n,double);
         out->lumped_mass_matrix=MEMALLOC(n,double);
         out->reactive_matrix=MEMALLOC(n,double);;
         out->main_diagonal_mass_matrix=MEMALLOC(n,double);	 
         out->main_diagonal_low_order_transport_matrix=MEMALLOC(n,double);
         out->u_coupler=Paso_Coupler_alloc(Paso_TransportProblem_borrowConnector(out),block_size);

         if ( ! (Paso_checkPtr(out->constraint_weights) || 
	         Paso_checkPtr(out->reactive_matrix) || Paso_checkPtr(out->main_diagonal_mass_matrix) || 
                 Paso_checkPtr(out->lumped_mass_matrix) || Paso_checkPtr(out->main_diagonal_low_order_transport_matrix)) && Paso_noError()  ) 
	 {
                 #pragma omp parallel for schedule(static) private(i)
                 for (i = 0; i < n; ++i) {
                    out->lumped_mass_matrix[i]=0.;
                    out->main_diagonal_low_order_transport_matrix[i]=0.;
                 }
	 }
  }
  if (Paso_noError()) {
     out->reference_counter=1;
     return out;
  } else {
     Paso_TransportProblem_free(out);
     return NULL;
  }
} 

void Paso_TransportProblem_reset(Paso_TransportProblem* in) 
{
    Paso_SystemMatrix_setValues(in->transport_matrix, 0.);
    Paso_SystemMatrix_setValues(in->mass_matrix, 0.);
    Paso_solve_free(in->iteration_matrix);
    in->valid_matrices=FALSE;
}

dim_t Paso_TransportProblem_getBlockSize(const Paso_TransportProblem* in)
{
   return in->transport_matrix->row_block_size;
}

Paso_Connector* Paso_TransportProblem_borrowConnector(const Paso_TransportProblem* in)
{
   return in->transport_matrix->pattern->col_connector;
}

index_t Paso_TransportProblem_getTypeId(const index_t solver,const index_t preconditioner, const index_t package,const  bool_t symmetry, Paso_MPIInfo *mpi_info) 
{
   return MATRIX_FORMAT_DEFAULT + MATRIX_FORMAT_BLK1;
}

/****************** REVISE ****************************/
void Paso_TransportProblem_setUpConstraint(Paso_TransportProblem* fctp,  const double* q, const double factor)
{
   dim_t i;
   register double m, rtmp;
   double factor2= fctp->dt_factor * factor;
   const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
   const index_t* main_iptr=Paso_SparseMatrix_borrowMainDiagonalPointer(fctp->mass_matrix->mainBlock);
   
   if ( fctp->valid_matrices ) {
      Paso_setError(VALUE_ERROR, "Paso_TransportProblem_insertConstraint: you must not insert a constraint is a valid system.");
      return;
   }
   if (factor<=0) {
      Paso_setError(VALUE_ERROR, "Paso_TransportProblem_insertConstraint: constraint_factor needs to be positive.");
      return;
   }

   
   #pragma omp for schedule(static) private(i,m,rtmp)
   for (i=0;i<n;++i) {
        if (q[i]>0) {
           m=fctp->mass_matrix->mainBlock->val[main_iptr[i]];
           rtmp=factor2 * (m == 0 ? 1 : m);
           fctp->constraint_weights[i]=rtmp;
           fctp->mass_matrix->mainBlock->val[main_iptr[i]]=m+rtmp;
        } else {
           fctp->constraint_weights[i]=0;
        }
   }
   fctp->constraint_factor=factor;
}

void Paso_TransportProblem_insertConstraint(Paso_TransportProblem* fctp,  const double* r,  double* source)
{
   dim_t i, n;
   n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);

   #pragma omp for schedule(static) private(i)
   for (i=0;i<n;++i) source[i]+=fctp->constraint_weights[i] * r[i];
}
/*****************************************************/