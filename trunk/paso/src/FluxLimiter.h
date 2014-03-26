
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


#ifndef INC_PASOFCTLIMITER
#define INC_PASOFCTLIMITER

#include "Transport.h"


typedef struct Paso_FCT_FluxLimiter {
      Paso_SystemMatrix *antidiffusive_fluxes;
      Esys_MPIInfo *mpi_info;
      double dt;
      double* u_tilde;
      double* MQ;   /* (M_C* Q_min, M_C* Q_max) */ 
      double* R;   /* (R-, R+) */
      /* paso::Coupler *MQ_coupler; */
      paso::Coupler *R_coupler;
      paso::Coupler *u_tilde_coupler;
      double*  borrowed_lumped_mass_matrix; /* borrowd reference */
} Paso_FCT_FluxLimiter;

#define Paso_FCT_FluxLimiter_getTotalNumRows(_f_) Paso_SystemMatrix_getTotalNumRows((_f_)->antidiffusive_fluxes)
#define Paso_FCT_FluxLimiter_getFluxPattern(_f_) ((_f_)->antidiffusive_fluxes->pattern)

PASO_DLL_API Paso_FCT_FluxLimiter* Paso_FCT_FluxLimiter_alloc(Paso_TransportProblem *fctp);
PASO_DLL_API void Paso_FCT_FluxLimiter_free(Paso_FCT_FluxLimiter * in);
PASO_DLL_API void Paso_FCT_FluxLimiter_setU_tilda(Paso_FCT_FluxLimiter* flux_limiter, const double *Mu_tilda);
PASO_DLL_API void Paso_FCT_FluxLimiter_addLimitedFluxes_Start(Paso_FCT_FluxLimiter* flux_limiter);
PASO_DLL_API void Paso_FCT_FluxLimiter_addLimitedFluxes_Complete(Paso_FCT_FluxLimiter* flux_limiter, double* b);

#endif /* #ifndef INC_PASOFCTLIMITER */
