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
    double dt_max;
    bool_t valid_matrices;

    Paso_SystemMatrix * transport_matrix;
    Paso_SystemMatrix * mass_matrix;

    double* u;
    double u_min;

    /* x */
    index_t *main_iptr;
    Paso_SystemMatrix * iteration_matrix;
    double* main_diagonal_low_order_transport_matrix;
    double* lumped_mass_matrix;

    Paso_MPIInfo *mpi_info;
    dim_t reference_counter;

} Paso_FCTransportProblem;

Paso_FCTransportProblem* Paso_FCTransportProblem_getReference(Paso_FCTransportProblem* in);
Paso_FCTransportProblem* Paso_FCTransportProblem_alloc(double theta, Paso_SystemMatrixPattern *pattern, int block_size);
double Paso_FCTransportProblem_getSafeTimeStepSize(Paso_FCTransportProblem* in);
void Paso_FCTransportProblem_setLowOrderOperator(Paso_FCTransportProblem * fc);
Paso_SystemMatrix* Paso_FCTransportProblem_borrowTransportMatrix(Paso_FCTransportProblem* in);
Paso_SystemMatrix* Paso_FCTransportProblem_borrowMassMatrix(Paso_FCTransportProblem* in);
double* Paso_FCTransportProblem_borrowLumpedMassMatrix(Paso_FCTransportProblem* in);
dim_t Paso_FCTransportProblem_getTotalNumRows(Paso_FCTransportProblem* in);
void Paso_FCTransportProblem_free(Paso_FCTransportProblem* in);
void Paso_SolverFCT_solve(Paso_FCTransportProblem* fctp, double* u, double dt, double* source, Paso_Options* options);
void Paso_FCTransportProblem_checkinSolution(Paso_FCTransportProblem* in, double* u);
void Paso_FCTransportProblem_applyPreAntiDiffusionCorrection(Paso_SystemMatrix *f,const double* u);
void Paso_SolverFCT_setMuPaLuPbQ(double* out,const double* M, const  double* u,const  double a, Paso_SystemMatrix *L, const  double b,const double* Q);
void Paso_SolverFCT_setQs(const double* u,double* QN, double* QP, Paso_SystemMatrix *L);
void Paso_FCTransportProblem_updateAntiDiffusionFlux(const Paso_FCTransportProblem * fc, Paso_SystemMatrix *flux_matrix,const double a, const double b, const double* u);
void Paso_FCTransportProblem_setRs(const Paso_SystemMatrix *f,const double* lumped_mass_matrix,const double* QN,const double* QP,double* RN,double* RP);
void Paso_FCTransportProblem_addCorrectedFluxes(double* f,Paso_SystemMatrix *flux_matrix,const double* RN,const double* RP);

#endif /* #ifndef INC_SOLVERFCT */
