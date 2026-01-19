
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __PASO_FLUXLIMITER_H__
#define __PASO_FLUXLIMITER_H__

#include "Transport.h"

namespace paso {


struct PASO_DLL_API FCT_FluxLimiter
{
    FCT_FluxLimiter(const_TransportProblem_ptr tp);
    ~FCT_FluxLimiter();

    inline dim_t getTotalNumRows() const
    {
        return antidiffusive_fluxes->getTotalNumRows();
    }

    inline SystemMatrixPattern_ptr getFluxPattern() const
    {
        return antidiffusive_fluxes->pattern;
    }

    void setU_tilde(const double* Mu_tilde);
    void addLimitedFluxes_Start();
    void addLimitedFluxes_Complete(double* b);

    SystemMatrix_ptr<double> antidiffusive_fluxes;
    escript::JMPI mpi_info;
    double dt;
    double* u_tilde;
    double* MQ;  // (M_C* Q_min, M_C* Q_max)
    double* R;   // (R-, R+)
    //Coupler_ptr MQ_coupler;
    Coupler_ptr<real_t> R_coupler;
    Coupler_ptr<real_t> u_tilde_coupler;
    double* borrowed_lumped_mass_matrix; // borrowed reference
};

} // namespace paso

#endif // __PASO_FLUXLIMITER_H__

