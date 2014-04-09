
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/************************************************************************************/

/* Paso: TransportProblem (see Paso_TransportSolver_solve)    */

/************************************************************************************/

/* Author: l.gross@uq.edu.au                                  */

/************************************************************************************/


#include "Transport.h"
#include "PasoUtil.h"
#include "Preconditioner.h"


void Paso_TransportProblem_free(Paso_TransportProblem* in)
{
     if (in!=NULL) {
        in->reference_counter--;
        if (in->reference_counter<=0) {
           Esys_MPIInfo_free(in->mpi_info);
           delete[] in->constraint_mask;
           delete[] in->reactive_matrix;
           delete[] in->main_diagonal_mass_matrix;
           delete[] in->lumped_mass_matrix;
           delete[] in->main_diagonal_low_order_transport_matrix;
           delete in;
        }
    }
}

Paso_TransportProblem* Paso_TransportProblem_getReference(Paso_TransportProblem* in) {
     if (in!=NULL) {
        ++(in->reference_counter);
     }
     return in;
}

Paso_TransportProblem* Paso_TransportProblem_alloc(paso::SystemMatrixPattern_ptr pattern, const int block_size)
{
    /* at the moment only block size 1 is supported */
    paso::SystemMatrixType matrix_type = MATRIX_FORMAT_DEFAULT+MATRIX_FORMAT_BLK1;
    Paso_TransportProblem* out=NULL;
    dim_t n,i;

    out=new Paso_TransportProblem;
    out->reference_counter=0;
    out->dt_max_R=LARGE_POSITIVE_FLOAT;
    out->dt_max_T=LARGE_POSITIVE_FLOAT;
    out->valid_matrices=FALSE;
    out->transport_matrix.reset(new paso::SystemMatrix(matrix_type, pattern, block_size, block_size, false));
    out->mass_matrix.reset(new paso::SystemMatrix(matrix_type, pattern, block_size, block_size, false));
    out->constraint_mask=NULL;
    out->mpi_info=Esys_MPIInfo_getReference(pattern->mpi_info);

    out->lumped_mass_matrix=NULL;
    out->main_diagonal_low_order_transport_matrix=NULL;
    out->reactive_matrix=NULL;
    out->main_diagonal_mass_matrix=NULL;

    if (Esys_noError()) {
        n = out->transport_matrix->getTotalNumRows();
        out->constraint_mask=new double[n]; /* ? */
        out->lumped_mass_matrix=new double[n];  /* ? */
        out->reactive_matrix=new double[n]; /* ? */
        out->main_diagonal_mass_matrix=new double[n]; /* ? */
        out->main_diagonal_low_order_transport_matrix=new double[n]; /* ? */

        #pragma omp parallel for schedule(static) private(i)
        for (i = 0; i < n; ++i) {
            out->lumped_mass_matrix[i]=0.;
            out->main_diagonal_low_order_transport_matrix[i]=0.;
		    out->constraint_mask[i]=0.;
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
    const dim_t n = in->transport_matrix->getTotalNumRows();
    in->transport_matrix->setValues(0.);
    in->mass_matrix->setValues(0.);
    Paso_solve_free(in->iteration_matrix.get());
    paso::util::zeroes(n, in->constraint_mask);
    in->valid_matrices=FALSE;
}


index_t Paso_TransportProblem_getTypeId(const index_t solver,const index_t preconditioner, const index_t package,const  bool symmetry, Esys_MPIInfo *mpi_info)
{
   return MATRIX_FORMAT_DEFAULT + MATRIX_FORMAT_BLK1;
}

void Paso_TransportProblem_setUpConstraint(Paso_TransportProblem* fctp,  const double* q)
{
   dim_t i;
   const dim_t n = fctp->transport_matrix->getTotalNumRows();
   if (fctp->valid_matrices) {
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
   const dim_t n = fctp->transport_matrix->getTotalNumRows();

   #pragma omp parallel for schedule(static) private(i)
   for (i=0;i<n;++i) {
       if (fctp->constraint_mask[i] > 0) source[i] =  r[i];
   }
}

