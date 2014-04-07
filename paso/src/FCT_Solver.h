
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
      paso::Coupler_ptr u_coupler;
      paso::Coupler_ptr u_old_coupler; /* last time step */
      
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

void Paso_FCT_setAntiDiffusionFlux_linearCN(paso::SystemMatrix_ptr flux_matrix, const Paso_TransportProblem* fct, 
				            const double dt, paso::const_Coupler_ptr u_tilde_coupler,  
				            paso::const_Coupler_ptr u_old_coupler);

void Paso_FCT_setAntiDiffusionFlux_BE(paso::SystemMatrix_ptr flux_matrix,
                                      const Paso_TransportProblem* fct, 
				      const double dt,
			              paso::const_Coupler_ptr u_coupler,  
				      paso::const_Coupler_ptr u_old_coupler);

void Paso_FCT_setAntiDiffusionFlux_CN(paso::SystemMatrix_ptr flux_matrix,
                                      const Paso_TransportProblem* fct, 
				      const double dt,
			          paso::const_Coupler_ptr u_coupler,  
				      paso::const_Coupler_ptr u_old_coupler);

void Paso_FCT_Solver_initialize(const double dt, Paso_FCT_Solver *fct_solver, Paso_Options* options, Paso_Performance* pp) ;
double Paso_FCT_Solver_getSafeTimeStepSize(Paso_TransportProblem* fctp);
void Paso_FCT_Solver_setMuPaLu(double* out, const double* M, paso::const_Coupler_ptr u_coupler, const double a, paso::const_SystemMatrix_ptr L);

#define Paso_FCT_Solver_getTheta(_fct_) ( ( (_fct_)->method == PASO_BACKWARD_EULER ) ? 1. : 0.5 )

#endif /* #ifndef INC_PASOFCT_SOLVER */

