
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#ifndef INC_PASOFCT_SOLVER
#define INC_PASOFCT_SOLVER

#include "Transport.h"
#include "FluxLimiter.h"
#include "Solver.h"





typedef struct Paso_FCT_Solver {
      Paso_TransportProblem* transportproblem;
      Esys_MPIInfo *mpi_info;
      Paso_FCT_FluxLimiter* flux_limiter;
      index_t method;
      double omega;
      double dt;
      double *b;
      double *z;
      double *du;
      Paso_Coupler *u_coupler;
      Paso_Coupler *u_old_coupler; /* last time step */
      
} Paso_FCT_Solver;

PASO_DLL_API
Paso_FCT_Solver* Paso_FCT_Solver_alloc(Paso_TransportProblem *fctp, Paso_Options* options);

PASO_DLL_API
void Paso_FCT_Solver_free(Paso_FCT_Solver *in);

PASO_DLL_API
err_t Paso_FCT_Solver_update(Paso_FCT_Solver *fct_solver, double *u, double *u_old, Paso_Options* options, Paso_Performance *pp) ;


PASO_DLL_API
void Paso_FCT_setLowOrderOperator(Paso_TransportProblem * fc);

PASO_DLL_API
err_t Paso_FCT_Solver_updateNL(Paso_FCT_Solver *fct_solver, double* u, double *u_old, Paso_Options* options, Paso_Performance *pp) ;

PASO_DLL_API
err_t Paso_FCT_Solver_update_LCN(Paso_FCT_Solver *fct_solver, double * u, double *u_old, Paso_Options* options, Paso_Performance *pp) ;

void Paso_FCT_setAntiDiffusionFlux_linearCN(Paso_SystemMatrix *flux_matrix, const Paso_TransportProblem* fct, 
				            const double dt, const Paso_Coupler* u_tilde_coupler,  
				            const Paso_Coupler* u_old_coupler);

void Paso_FCT_setAntiDiffusionFlux_BE(Paso_SystemMatrix *flux_matrix,
                                      const Paso_TransportProblem* fct, 
				      const double dt,
			              const Paso_Coupler* u_coupler,  
				      const Paso_Coupler* u_old_coupler);

void Paso_FCT_setAntiDiffusionFlux_CN(Paso_SystemMatrix *flux_matrix,
                                      const Paso_TransportProblem* fct, 
				      const double dt,
			              const Paso_Coupler* u_coupler,  
				      const Paso_Coupler* u_old_coupler);

void Paso_FCT_Solver_initialize(const double dt, Paso_FCT_Solver *fct_solver, Paso_Options* options, Paso_Performance* pp) ;
double Paso_FCT_Solver_getSafeTimeStepSize(Paso_TransportProblem* fctp);
void Paso_FCT_Solver_setMuPaLu(double* out, const double* M, const Paso_Coupler* u_coupler, const double a, const Paso_SystemMatrix *L);

#define Paso_FCT_Solver_getTheta(_fct_) ( ( (_fct_)->method == PASO_BACKWARD_EULER ) ? 1. : 0.5 )

#endif /* #ifndef INC_PASOFCT_SOLVER */

