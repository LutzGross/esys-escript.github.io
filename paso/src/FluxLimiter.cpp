
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/


/****************************************************************************/

/* Paso: CFT: fluxlimiter for a transport problem
 *
 ****************************************************************************/

/* Author: l.gross@uq.edu.au */

/****************************************************************************/

#include "FluxLimiter.h"

namespace paso {

FCT_FluxLimiter::FCT_FluxLimiter(const_TransportProblem_ptr tp)
{
    const dim_t n = tp->transport_matrix->getTotalNumRows();
    const dim_t blockSize = tp->getBlockSize();

    mpi_info = tp->mpi_info;
    u_tilde = new double[n];
    MQ = new double[2*n];
    R = new double[2*n];

    R_coupler.reset(new Coupler<real_t>(tp->borrowConnector(), 2*blockSize, mpi_info));
    u_tilde_coupler.reset(new Coupler<real_t>(tp->borrowConnector(), blockSize, mpi_info));
    antidiffusive_fluxes.reset(new SystemMatrix<double>(
                tp->transport_matrix->type, tp->transport_matrix->pattern,
                tp->transport_matrix->row_block_size,
                tp->transport_matrix->col_block_size, true,
                tp->transport_matrix->getRowFunctionSpace(),
                tp->transport_matrix->getColumnFunctionSpace()));
    borrowed_lumped_mass_matrix = tp->lumped_mass_matrix;
}

FCT_FluxLimiter::~FCT_FluxLimiter()
{
    delete[] u_tilde;
    delete[] MQ;
    delete[] R;
}

// sets the predictor u_tilde from Mu_tilde by solving M_C * u_tilde = Mu_tilde
// and calculates the limiters QP and QN
void FCT_FluxLimiter::setU_tilde(const double* Mu_tilde)
{
    const real_t LARGE_POSITIVE_FLOAT = escript::DataTypes::real_t_max();
    const dim_t n = getTotalNumRows();
    const_SystemMatrixPattern_ptr pattern(getFluxPattern());

#pragma omp parallel for
    for (index_t i = 0; i < n; ++i) {
        const double m = borrowed_lumped_mass_matrix[i];
        if (m > 0) {
            u_tilde[i] = Mu_tilde[i]/m;
        } else {
            u_tilde[i] = Mu_tilde[i];
        }
    }

    // distribute u_tilde
    u_tilde_coupler->startCollect(u_tilde);

    // calculate
    //   MQ_P[i] = lumped_mass_matrix[i] * max_{j} (\tilde{u}[j]- \tilde{u}[i])
    //   MQ_N[i] = lumped_mass_matrix[i] * min_{j} (\tilde{u}[j]- \tilde{u}[i])

    // first we calculate the min and max of u_tilde in the main block
    // QP, QN are used to hold the result
#pragma omp parallel for
    for (index_t i = 0; i < n; ++i) {
        if (borrowed_lumped_mass_matrix[i] > 0) { // no constraint
            double u_min_i = LARGE_POSITIVE_FLOAT;
            double u_max_i = -LARGE_POSITIVE_FLOAT;
            #pragma ivdep
            for (index_t iptr_ij = pattern->mainPattern->ptr[i];
                         iptr_ij < pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
                const index_t j=pattern->mainPattern->index[iptr_ij];
                const double u_j = u_tilde[j];
                u_min_i = std::min(u_min_i, u_j);
                u_max_i = std::max(u_max_i, u_j);
            }
            MQ[2*i] = u_min_i;
            MQ[2*i+1] = u_max_i;

        } else {
            MQ[2*i  ] = LARGE_POSITIVE_FLOAT;
            MQ[2*i+1] = LARGE_POSITIVE_FLOAT;
        }
    }

    // complete distribute u_tilde
    u_tilde_coupler->finishCollect();
    const double* remote_u_tilde = u_tilde_coupler->borrowRemoteData();

    // now we look at the couple matrix and set the final value for QP, QN
#pragma omp parallel for
    for (index_t i = 0; i < n; ++i) {
        if (borrowed_lumped_mass_matrix[i] > 0) { // no constraint
            const double u_i = u_tilde[i];
            double u_min_i = MQ[2*i];
            double u_max_i = MQ[2*i+1];
            #pragma ivdep
            for (index_t iptr_ij = pattern->col_couplePattern->ptr[i];
                         iptr_ij < pattern->col_couplePattern->ptr[i+1];
                         iptr_ij++) {
                const index_t j = pattern->col_couplePattern->index[iptr_ij];
                const double u_j = remote_u_tilde[j];
                u_min_i = std::min(u_min_i, u_j);
                u_max_i = std::max(u_max_i, u_j);
            }
            MQ[2*i  ] = (u_min_i-u_i)*borrowed_lumped_mass_matrix[i];//M_C*Q_min
            MQ[2*i+1] = (u_max_i-u_i)*borrowed_lumped_mass_matrix[i];//M_C*Q_max
        }
    }
}

// starts to update a vector (not given yet) from the antidiffusion fluxes
// in flux_limiter->antidiffusive_fluxes (needs u_tilde and Q)
void FCT_FluxLimiter::addLimitedFluxes_Start()
{
    const dim_t n = getTotalNumRows();
    const_SystemMatrixPattern_ptr pattern(getFluxPattern());
    const double* remote_u_tilde = u_tilde_coupler->borrowRemoteData();
    SystemMatrix_ptr<double> adf(antidiffusive_fluxes);

#pragma omp parallel for
    for (dim_t i = 0; i < n; ++i) {
        double R_N_i = 1;
        double R_P_i = 1;
        if (borrowed_lumped_mass_matrix[i] > 0) { // no constraint
            const double u_tilde_i = u_tilde[i];
            double P_P_i = 0.;
            double P_N_i = 0.;
            const double MQ_min = MQ[2*i];
            const double MQ_max = MQ[2*i+1];
            #pragma ivdep
            for (index_t iptr_ij = pattern->mainPattern->ptr[i];
                         iptr_ij < pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
                const index_t j = pattern->mainPattern->index[iptr_ij];
                if (i != j ) {
                    const double f_ij=adf->mainBlock->val[iptr_ij];
                    const double u_tilde_j=u_tilde[j];
                    /* pre-limiter */
                    if (f_ij * (u_tilde_j-u_tilde_i) >= 0) {
                        adf->mainBlock->val[iptr_ij]=0;
                    } else {
                        if (f_ij <=0) {
                            P_N_i+=f_ij;
                        } else {
                            P_P_i+=f_ij;
                        }
                    }
                }
            }

            // now the couple matrix
            #pragma ivdep
            for (index_t iptr_ij = pattern->col_couplePattern->ptr[i];
                    iptr_ij < pattern->col_couplePattern->ptr[i+1]; ++iptr_ij) {
                const index_t j=pattern->col_couplePattern->index[iptr_ij];
                const double f_ij=adf->col_coupleBlock->val[iptr_ij];
                const double u_tilde_j=remote_u_tilde[j];
                // pre-limiter
                if (f_ij * (u_tilde_j-u_tilde_i) >= 0) {
                    adf->col_coupleBlock->val[iptr_ij]=0;
                } else {
                    if (f_ij <= 0) {
                        P_N_i+=f_ij;
                    } else {
                        P_P_i+=f_ij;
                    }
                }
            }
            /* finally the R+ and R- are calculated */
            if (P_N_i<0) R_N_i=std::min(1., MQ_min/P_N_i);
            if (P_P_i>0) R_P_i=std::min(1., MQ_max/P_P_i);
        }
        R[2*i]   = R_N_i;
        R[2*i+1] = R_P_i;
    }

    // now we kick off the distribution of the R's
    R_coupler->startCollect(R);
}


// completes the exchange of the R factors and adds the weighted
// antidiffusion fluxes to the residual b
void FCT_FluxLimiter::addLimitedFluxes_Complete(double* b)
{
    const dim_t n = getTotalNumRows();
    const_SystemMatrixPattern_ptr pattern(getFluxPattern());
    const_SystemMatrix_ptr<double> adf(antidiffusive_fluxes);
    const double* remote_R = R_coupler->finishCollect();

#pragma omp parallel for
    for (dim_t i = 0; i < n; ++i) {
        const double R_N_i = R[2*i];
        const double R_P_i = R[2*i+1];
        double f_i = b[i];

        #pragma ivdep
        for (index_t iptr_ij = pattern->mainPattern->ptr[i];
                     iptr_ij < pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
            const index_t j = pattern->mainPattern->index[iptr_ij];
            const double f_ij = adf->mainBlock->val[iptr_ij];
            const double R_P_j = R[2*j+1];
            const double R_N_j = R[2*j];
            const double rtmp=(f_ij>=0 ? std::min(R_P_i, R_N_j) : std::min(R_N_i, R_P_j));
            f_i += f_ij*rtmp;
        }
        #pragma ivdep
        for (index_t iptr_ij=pattern->col_couplePattern->ptr[i];
                     iptr_ij<pattern->col_couplePattern->ptr[i+1]; ++iptr_ij) {
            const index_t j = pattern->col_couplePattern->index[iptr_ij];
            const double f_ij = adf->col_coupleBlock->val[iptr_ij];
            const double R_P_j = remote_R[2*j+1];
            const double R_N_j = remote_R[2*j];
            const double rtmp=(f_ij>=0) ? std::min(R_P_i, R_N_j) : std::min(R_N_i, R_P_j);
            f_i += f_ij*rtmp;
        }
        b[i]=f_i;
    } // end of i loop
}

} // namespace paso

