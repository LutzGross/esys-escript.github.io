
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


#ifndef INC_PASOFCT
#define INC_PASOFCT

#include "Functions.h"
#include "Transport.h"


typedef struct Paso_FCTSolver {
      Paso_TransportProblem* transportproblem;
      Paso_SystemMatrix *flux_matrix_m;
      double dt;
      double* uTilde_n;
      double* QN_n;
      double* QP_n;
      double* RN_m;
      double* RP_m;
      Paso_Coupler *QN_n_coupler;
      Paso_Coupler *QP_n_coupler;
      Paso_Coupler *RN_m_coupler;
      Paso_Coupler *RP_m_coupler;
      Paso_Coupler *uTilde_n_coupler;
      Paso_Coupler *u_m_coupler;
} Paso_FCTSolver;


PASO_DLL_API
err_t Paso_FCTSolver_Function_call(Paso_Function * F,double* value, const double* arg, Paso_Performance *pp);

PASO_DLL_API
void Paso_FCTSolver_Function_free(Paso_Function * in);

PASO_DLL_API
Paso_Function* Paso_FCTSolver_Function_alloc(Paso_TransportProblem *tp, Paso_Options* options);
	       
PASO_DLL_API
err_t Paso_FCTSolver_solve(Paso_Function* F, double* u, double dt, Paso_Options* options, Paso_Performance *pp);

PASO_DLL_API
double Paso_FCTSolver_getSafeTimeStepSize(Paso_TransportProblem* fctp);

PASO_DLL_API
void Paso_FCTSolver_applyPreAntiDiffusionCorrection(Paso_SystemMatrix *f,const Paso_Coupler* u_coupler);

PASO_DLL_API
void Paso_FCTSolver_setQs(const Paso_Coupler* u_coupler,double* QN, double* QP, const Paso_SystemMatrix *L);

PASO_DLL_API
void Paso_FCTSolver_setAntiDiffusionFlux(const double dt, const Paso_TransportProblem * fc, Paso_SystemMatrix *flux_matrix, const Paso_Coupler* u_coupler);

PASO_DLL_API
void Paso_FCTSolver_setRs(const Paso_SystemMatrix *f,const double* lumped_mass_matrix,const Paso_Coupler* QN,const Paso_Coupler* QP,double* RN,double* RP);

PASO_DLL_API
void Paso_FCTSolver_addCorrectedFluxes(double* f,const Paso_SystemMatrix *flux_matrix,const Paso_Coupler* RN,const Paso_Coupler* RP);


PASO_DLL_API
void Paso_FCTSolver_setMuPaLu(double* out, const double* M, const Paso_Coupler* u_coupler, const double a, const Paso_SystemMatrix *L);

PASO_DLL_API
void Paso_FCTSolver_setLowOrderOperator(Paso_TransportProblem * fc);

PASO_DLL_API
void Paso_FCTSolver_setUp(Paso_TransportProblem* fctp, const double dt, const double *u, double* b, double* uTilde,
                          Paso_Coupler* uTilde_coupler, double *QN, Paso_Coupler* QN_coupler, double *QP, Paso_Coupler* QP_coupler,
                          Paso_Options* options, Paso_Performance* pp);

PASO_DLL_API
err_t Paso_FCTSolver_setUpRightHandSide(Paso_TransportProblem* fctp, const double dt, 
				        const double *u_m, Paso_Coupler* u_m_coupler,  double * z_m,
                                        Paso_SystemMatrix* flux_matrix, Paso_Coupler* uTilde_coupler, const double *b, 
                                        Paso_Coupler* QN_coupler, Paso_Coupler* QP_coupler,
                                        double *RN_m, Paso_Coupler* RN_m_coupler, double* RP_m, Paso_Coupler* RP_m_coupler,
                                        Paso_Performance* pp);


#endif /* #ifndef INC_PASOFCT */
