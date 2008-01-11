/* $Id: $ */

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

#ifndef INC_SOLVERFCT
#define INC_SOLVERFCT

#include "SystemMatrix.h"

typedef struct Paso_FCTransportProblem {

    double theta;
    double dt;
    double dt_max;
    bool_t valid_matrices;

    double* u;

    Paso_SystemMatrix * transport_matrix;
    Paso_SystemMatrix * flux_matrix;
    double* lumped_mass_matrix;
    double* row_sum_flux_matrix;
    double* transport_matrix_diagonal;

    dim_t num_colors;
    index_t *colorOf;
    index_t *main_iptr;
    
    Paso_MPIInfo *mpi_info;
    dim_t reference_counter;

} Paso_FCTransportProblem;

void Paso_FCTransportProblem_free(Paso_FCTransportProblem* in);
Paso_FCTransportProblem* Paso_FCTransportProblem_getReference(Paso_FCTransportProblem* in);
Paso_SystemMatrix* Paso_FCTransportProblem_borrowTransportMatrix(Paso_FCTransportProblem* in);
Paso_SystemMatrix* Paso_FCTransportProblem_borrowFluxMatrix(Paso_FCTransportProblem* in);
double* Paso_FCTransportProblem_borrowLumpedMassMatrix(Paso_FCTransportProblem* in);
dim_t Paso_FCTransportProblem_getTotalNumRows(Paso_FCTransportProblem* in);
Paso_FCTransportProblem* Paso_FCTransportProblem_alloc(double theta, double dt_max, Paso_SystemMatrixPattern *pattern, int block_size);
void Paso_FCTransportProblem_setAntiDiffusiveFlux(Paso_FCTransportProblem * fc, double * u, double *u_remote, double* fa);
void Paso_FCTransportProblem_addAdvectivePart(Paso_FCTransportProblem * fc, double alpha);
void Paso_SolverFCT_solve(Paso_FCTransportProblem* fctp, double* u, double dt, double* source, Paso_Options* options);
void Paso_FCTransportProblem_checkinSolution(Paso_FCTransportProblem* in, double* u) ;
void Paso_FCTransportProblem_setFlux(Paso_FCTransportProblem * fc, double * u, double* fa);

#endif /* #ifndef INC_SOLVERFCT */
