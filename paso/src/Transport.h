
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


#ifndef INC_PASOTRANSPORT
#define INC_PASOTRANSPORT

#define DT_FACTOR_MAX 1000.

#include "SystemMatrix.h"
#include "Options.h"
#include "performance.h"
#include "Paso.h"

typedef struct Paso_TransportProblem {

    bool_t useBackwardEuler;

    bool_t valid_matrices;
/****************** REVISE ****************************/
    double dt_factor;  
    double dt_max; 
    double constraint_factor;  
    double* constraint_weights;
/*****************************************************/

    Paso_SystemMatrix * transport_matrix;
    Paso_SystemMatrix * mass_matrix;
    
    Paso_Coupler* u_coupler;
    Paso_SystemMatrix * iteration_matrix;
    double* main_diagonal_low_order_transport_matrix;
    double* lumped_mass_matrix;
    double* reactive_matrix;
    double* main_diagonal_mass_matrix;

    Paso_MPIInfo *mpi_info;
    dim_t reference_counter;

} Paso_TransportProblem;



PASO_DLL_API
Paso_TransportProblem* Paso_TransportProblem_getReference(Paso_TransportProblem* in);

PASO_DLL_API
Paso_TransportProblem* Paso_TransportProblem_alloc(bool_t useBackwardEuler, Paso_SystemMatrixPattern *pattern, int block_size);

PASO_DLL_API
dim_t Paso_TransportProblem_getBlockSize(const Paso_TransportProblem* in);

PASO_DLL_API
double Paso_TransportProblem_getSafeTimeStepSize(Paso_TransportProblem* in);

PASO_DLL_API
Paso_SystemMatrix* Paso_TransportProblem_borrowTransportMatrix(Paso_TransportProblem* in);

PASO_DLL_API
Paso_SystemMatrix* Paso_TransportProblem_borrowMassMatrix(Paso_TransportProblem* in);

PASO_DLL_API
void Paso_TransportProblem_solve(Paso_TransportProblem* fctp, double* u, double dt, double* u0, double* q, Paso_Options* options);

PASO_DLL_API
double* Paso_TransportProblem_borrowLumpedMassMatrix(Paso_TransportProblem* in);

PASO_DLL_API
dim_t Paso_TransportProblem_getTotalNumRows(Paso_TransportProblem* in);

PASO_DLL_API
void Paso_TransportProblem_free(Paso_TransportProblem* in);

PASO_DLL_API
void Paso_TransportProblem_reset(Paso_TransportProblem* in);

PASO_DLL_API
Paso_Connector* Paso_TransportProblem_borrowConnector(const Paso_TransportProblem* in);

PASO_DLL_API
index_t Paso_TransportProblem_getTypeId(const index_t solver,const index_t preconditioner, const index_t package,const  bool_t symmetry, Paso_MPIInfo *mpi_info);

PASO_DLL_API
void Paso_TransportProblem_insertConstraint(Paso_TransportProblem* fctp,  const double* r,  double* source);


PASO_DLL_API
void Paso_TransportProblem_setUpConstraint(Paso_TransportProblem* fctp,  const double* q, const double factor);

#define Paso_TransportProblem_borrowMainDiagonalPointer(_fct_) Paso_SparseMatrix_borrowMainDiagonalPointer((_fct_)->mass_matrix->mainBlock)
#define Paso_Transport_getTheta(_fct_) ( ( (_fct_)->useBackwardEuler ) ? 1. : 0.5 )

#endif /* #ifndef INC_PASOTRANSPORT */
