
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

/* Paso: TransportProblem (see Paso_TransportSolver_solve)    */                                      

/**************************************************************/

/* Author: l.gross@uq.edu.au                                  */

/**************************************************************/


#include "Transport.h"
#include "PasoUtil.h"


/**************************************************************/

/* free all used memory                                       */

void Paso_TransportProblem_free(Paso_TransportProblem* in) {
     if (in!=NULL) {
        in->reference_counter--;
        if (in->reference_counter<=0) {
           Paso_SystemMatrix_free(in->transport_matrix);
           Paso_SystemMatrix_free(in->mass_matrix);
           Paso_SystemMatrix_free(in->iteration_matrix);
           Esys_MPIInfo_free(in->mpi_info);
           MEMFREE(in->constraint_mask);
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

Paso_TransportProblem* Paso_TransportProblem_alloc(Paso_SystemMatrixPattern *pattern, const int block_size) 
{
     Paso_SystemMatrixType matrix_type=MATRIX_FORMAT_DEFAULT+MATRIX_FORMAT_BLK1;  /* at the moment only block size 1 is supported */
     Paso_TransportProblem* out=NULL;
     dim_t n,i;
     
     out=MEMALLOC(1,Paso_TransportProblem);
     if (Esys_checkPtr(out)) return NULL;
     out->reference_counter=0; 
     out->dt_max_R=LARGE_POSITIVE_FLOAT;
     out->dt_max_T=LARGE_POSITIVE_FLOAT;
     out->valid_matrices=FALSE;
     out->transport_matrix=Paso_SystemMatrix_alloc(matrix_type,pattern,block_size,block_size,FALSE);
     out->mass_matrix=Paso_SystemMatrix_alloc(matrix_type,pattern,block_size,block_size,FALSE);
     out->iteration_matrix=NULL;
     out->constraint_mask=NULL;
     out->mpi_info=Esys_MPIInfo_getReference(pattern->mpi_info);

     out->lumped_mass_matrix=NULL;
     out->main_diagonal_low_order_transport_matrix=NULL;
     out->reactive_matrix=NULL;
     out->main_diagonal_mass_matrix=NULL;

     if (Esys_noError()) {
         n=Paso_SystemMatrix_getTotalNumRows(out->transport_matrix);
         out->constraint_mask=MEMALLOC(n,double); /* ? */
         out->lumped_mass_matrix=MEMALLOC(n,double);  /* ? */
         out->reactive_matrix=MEMALLOC(n,double); /* ? */
         out->main_diagonal_mass_matrix=MEMALLOC(n,double); /* ? */	 
         out->main_diagonal_low_order_transport_matrix=MEMALLOC(n,double); /* ? */


         if ( ! (Esys_checkPtr(out->constraint_mask) || 
	         Esys_checkPtr(out->reactive_matrix) || Esys_checkPtr(out->main_diagonal_mass_matrix) || 
                 Esys_checkPtr(out->lumped_mass_matrix) || Esys_checkPtr(out->main_diagonal_low_order_transport_matrix)) && Esys_noError()  ) 
	 {
                 #pragma omp parallel for schedule(static) private(i)
                 for (i = 0; i < n; ++i) {
                    out->lumped_mass_matrix[i]=0.;
                    out->main_diagonal_low_order_transport_matrix[i]=0.;
		    out->constraint_mask[i]=0.;
                 }
	 }
  }
  if (Esys_noError()) {
     out->reference_counter=1;
     return out;
  } else {
     Paso_TransportProblem_free(out);
     return NULL;
  }
} 

void Paso_TransportProblem_reset(Paso_TransportProblem* in) 
{
    const dim_t n = Paso_SystemMatrix_getTotalNumRows(in->transport_matrix);
    Paso_SystemMatrix_setValues(in->transport_matrix, 0.);
    Paso_SystemMatrix_setValues(in->mass_matrix, 0.);
    Paso_solve_free(in->iteration_matrix);
    Paso_zeroes( n, in->constraint_mask);
    in->valid_matrices=FALSE;
}


index_t Paso_TransportProblem_getTypeId(const index_t solver,const index_t preconditioner, const index_t package,const  bool_t symmetry, Esys_MPIInfo *mpi_info) 
{
   return MATRIX_FORMAT_DEFAULT + MATRIX_FORMAT_BLK1;
}

void Paso_TransportProblem_setUpConstraint(Paso_TransportProblem* fctp,  const double* q)
{
   dim_t i;
   const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);   
   if ( fctp->valid_matrices ) {
      Esys_setError(VALUE_ERROR, "Paso_TransportProblem_setUpConstraint: Cannot insert a constraint into a valid system.");
      return;
   }
   
   #pragma omp parallel for schedule(static) private(i)
   for (i=0;i<n;++i) {
        if (q[i]>0) {
	  fctp->constraint_mask[i]=1;
        } else {
           fctp->constraint_mask[i]=0;
        }
   }
   fctp->valid_matrices=FALSE;
}

void Paso_TransportProblem_insertConstraint(Paso_TransportProblem* fctp,  const double* r,  double* source)
{
   dim_t i;
   const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);

   #pragma omp parallel for schedule(static) private(i)
   for (i=0;i<n;++i) {
       if (fctp->constraint_mask[i] > 0) source[i] =  r[i];
   }
}

