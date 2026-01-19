
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
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

/* Paso: Transport solver with flux correction (L is row sum zero)
 *
 *   - Mv_t=Lv   v(0)=u
 *
 *  to return v(dt)
 *
*/
/****************************************************************************/

/* Author: l.gross@uq.edu.au */

/****************************************************************************/

#include "FCT_Solver.h"
#include "PasoUtil.h"
#include "Preconditioner.h"

#include <iostream>

#define MIN3(_arg1_,_arg2_,_arg3_) std::min(_arg1_, std::min(_arg2_,_arg3_))

namespace paso {

static const real_t LARGE_POSITIVE_FLOAT = escript::DataTypes::real_t_max();

FCT_Solver::FCT_Solver(const_TransportProblem_ptr tp, Options* options) :
    transportproblem(tp),
    omega(0),
    z(NULL),
    du(NULL)
{
    const dim_t blockSize = tp->getBlockSize();
    const dim_t n = tp->getTotalNumRows();
    mpi_info = tp->mpi_info;
    flux_limiter = new FCT_FluxLimiter(tp);
    b = new double[n];
    if (options->ode_solver == PASO_CRANK_NICOLSON || options->ode_solver == PASO_BACKWARD_EULER) {
        du = new double[n];
        z = new double[n];
    }
    u_coupler.reset(new Coupler<real_t>(tp->borrowConnector(), blockSize, mpi_info));
    u_old_coupler.reset(new Coupler<real_t>(tp->borrowConnector(), blockSize, mpi_info));

    if (options->ode_solver == PASO_LINEAR_CRANK_NICOLSON) {
        method = PASO_LINEAR_CRANK_NICOLSON;
    } else if (options->ode_solver == PASO_CRANK_NICOLSON) {
        method = PASO_CRANK_NICOLSON;
    } else if (options->ode_solver == PASO_BACKWARD_EULER) {
        method = PASO_BACKWARD_EULER;
    } else {
        throw PasoException("FCT_Solver: unknown integration scheme.");
    }
}

FCT_Solver::~FCT_Solver()
{
    delete flux_limiter;
    delete[] b;
    delete[] z;
    delete[] du;
}

// modifies the main diagonal of the iteration matrix to introduce new dt
void FCT_Solver::initialize(double _dt, Options* options, Performance* pp)
{
    const real_t EPSILON = escript::DataTypes::real_t_eps();
    const_TransportProblem_ptr fctp(transportproblem);
    const index_t* main_iptr = fctp->borrowMainDiagonalPointer();
    const dim_t n = fctp->transport_matrix->getTotalNumRows();
    const double theta = getTheta();
    omega = 1. / (_dt * theta);
    dim_t i;
    Options options2;

    solve_free(fctp->iteration_matrix.get());
    //   fctp->iteration_matrix[i,i]=m[i]/(dt theta) -l[i,i]
    dt = _dt;
    #pragma omp parallel for private(i)
    for (i = 0; i < n; ++i) {
        const double m_i = fctp->lumped_mass_matrix[i];
        const double l_ii = fctp->main_diagonal_low_order_transport_matrix[i];
        if ( m_i > 0 ) {
            fctp->iteration_matrix->mainBlock->val[main_iptr[i]] = m_i * omega - l_ii;
        } else {
            fctp->iteration_matrix->mainBlock->val[main_iptr[i]] = std::abs(m_i * omega - l_ii)/(EPSILON*EPSILON);
        }
    }

    // allocate preconditioner/solver
    options2.verbose = options->verbose;
    if (method == PASO_LINEAR_CRANK_NICOLSON) {
        options2.preconditioner = PASO_GS;
    } else {
        options2.preconditioner = PASO_JACOBI;
        //options2.preconditioner = PASO_GS;
    }
    options2.use_local_preconditioner = false;
    options2.sweeps = -1;

    Performance_startMonitor(pp, PERFORMANCE_PRECONDITIONER_INIT);
    fctp->iteration_matrix->setPreconditioner(&options2);
    Performance_stopMonitor(pp, PERFORMANCE_PRECONDITIONER_INIT);
}

// entry point for update procedures
SolverResult FCT_Solver::update(double* u, double* u_old, Options* options,
                                Performance* pp)
{
    SolverResult err_out = NoError;

    if (method == PASO_LINEAR_CRANK_NICOLSON) {
        err_out = updateLCN(u, u_old, options, pp);
    } else if (method == PASO_CRANK_NICOLSON) {
        err_out = updateNL(u, u_old, options, pp);
    } else if (method == PASO_BACKWARD_EULER) {
        err_out = updateNL(u, u_old, options, pp);
    } else {
        err_out = InputError;
    }
    return err_out;
}

/// linear crank-nicolson update
SolverResult FCT_Solver::updateLCN(double* u, double* u_old, Options* options,
                                   Performance* pp)
{
    dim_t sweep_max, i;
    double const RTOL = options->tolerance;
    const dim_t n = transportproblem->getTotalNumRows();
    SystemMatrix_ptr<double> iteration_matrix(transportproblem->iteration_matrix);
    const index_t* main_iptr = transportproblem->borrowMainDiagonalPointer();
    SolverResult errorCode = NoError;
    double norm_u_tilde;

    u_old_coupler->startCollect(u_old);
    u_old_coupler->finishCollect();

    // b[i]=m*u_tilde[i] = m u_old[i] + dt/2 sum_{j <> i} l_{ij}*(u_old[j]-u_old[i])
    //     = u_tilde[i]  = u_old[i] where constraint m<0.
    //  note that iteration_matrix stores the negative values of the
    //  low order transport matrix l. Therefore a=-dt*0.5 is used.

    setMuPaLu(b, u_old_coupler, -dt*0.5);
    /* solve for u_tilde : u_tilda = m^{-1} * b   */
    flux_limiter->setU_tilde(b);
    // u_tilde_connector is completed

    // calculate anti-diffusive fluxes for u_tilde
    setAntiDiffusionFlux_linearCN(flux_limiter->antidiffusive_fluxes);

    /* b_i += sum_{j} limitation factor_{ij} * antidiffusive_flux_{ij} */
    flux_limiter->addLimitedFluxes_Start();
    flux_limiter->addLimitedFluxes_Complete(b);

    util::scale(n, b, omega);
    // solve (m-dt/2*L) u = b in the form (omega*m-L) u = b * omega with omega*dt/2=1
   #pragma omp for private(i) schedule(static)
    for (i = 0; i < n; ++i) {
       if (!(transportproblem->lumped_mass_matrix[i] > 0)) {
         b[i] = flux_limiter->u_tilde[i]
             * transportproblem->iteration_matrix->mainBlock->val[main_iptr[i]];
       }
    }
    // initial guess is u<- -u + 2*u_tilde
    util::update(n, -1., u, 2., flux_limiter->u_tilde);

    sweep_max = std::max((int) (- 2 * log(RTOL)/log(2.)-0.5),1);
    norm_u_tilde = util::lsup(n, flux_limiter->u_tilde, flux_limiter->mpi_info);
    if (options->verbose) {
        std::cout << "FCT_Solver::updateLCN: u_tilde lsup = " << norm_u_tilde
            << " (rtol = " << RTOL*norm_u_tilde << ", max. sweeps = "
            << sweep_max << ")" << std::endl;
    }
    errorCode = Preconditioner_Smoother_solve_byTolerance(iteration_matrix,
            ((Preconditioner*)(iteration_matrix->solver_p))->gs, u, b, RTOL,
            &sweep_max, true);
    if (errorCode == NoError) {
        if (options->verbose)
            std::cout << "FCT_Solver::updateLCN: convergence after "
                << sweep_max << " Gauss-Seidel steps." << std::endl;
        errorCode = NoError;
    } else {
        if (options->verbose)
            std::cout << "FCT_Solver::updateLCN: Gauss-Seidel failed within "
                << sweep_max << " steps (rel. tolerance " << RTOL << ")."
                << std::endl;
        errorCode = MaxIterReached;
    }
    return errorCode;
}

SolverResult FCT_Solver::updateNL(double* u, double* u_old, Options* options,
                                  Performance* pp)
{
    // number of rates >=critical_rate accepted before divergence is triggered
    const dim_t num_critical_rates_max = 3;
    // expected value of convergence rate
    const double critical_rate = 0.95;

    const_TransportProblem_ptr fctp(transportproblem);
    dim_t i;
    double norm_u_tilde, norm_du=LARGE_POSITIVE_FLOAT, norm_du_old, rate=1.;
    const dim_t n = fctp->transport_matrix->getTotalNumRows();
    const double atol = options->absolute_tolerance;
    const double rtol = options->tolerance;
    const dim_t max_m = options->iter_max;
    dim_t m = 0, num_critical_rates = 0;
    SolverResult errorCode = NoError;
    bool converged=false, max_m_reached=false, diverged=false;
    /* //////////////////////////////////////////////////////////////////// */

    options->num_iter=0;
    u_old_coupler->startCollect(u_old);
    u_old_coupler->finishCollect();
    // prepare u_tilde and flux limiter
    if (method == PASO_BACKWARD_EULER) {
        // b[i]=m_i* u_old[i]
        #pragma omp for private(i) schedule(static)
        for (i = 0; i < n; ++i) {
            if (fctp->lumped_mass_matrix[i] > 0 ) {
                b[i]=u_old[i]* fctp->lumped_mass_matrix[i];
            } else {
                b[i]=u_old[i];
            }
        }
    } else {
       /* b[i]=m_i* u_old[i] + dt/2 sum_{j <> i} l_{ij}*(u_old[j]-u_old[i]) = m_i * u_tilde_i where m_i>0
        *     = u_old[i]  otherwise
        * note that iteration_matrix stores the negative values of the
        * low order transport matrix l. Therefore a=-dt*0.5 is used. */
        setMuPaLu(b, u_old_coupler, -dt*0.5);
    }
    flux_limiter->setU_tilde(b); // u_tilde = m^{-1} b */
    // u_tilde_connector is completed
    /************************************************************************/
    // calculate stopping criterion
    norm_u_tilde=util::lsup(n, flux_limiter->u_tilde, flux_limiter->mpi_info);
    const double ATOL = rtol * norm_u_tilde + atol;
    if (options->verbose)
        std::cout << "FCT_Solver::updateNL: iteration starts u_tilde lsup = "
            << norm_u_tilde << " (abs. tol = " << ATOL << ")" << std::endl;

    // u_old is an initial guess for u
    util::copy(n, u, u_old);

    while (!converged && !diverged && !max_m_reached) {
        u_coupler->startCollect(u);
        u_coupler->finishCollect();

        // set antidiffusive_flux_{ij} for u
        if (method == PASO_BACKWARD_EULER) {
             setAntiDiffusionFlux_BE(flux_limiter->antidiffusive_fluxes);
         } else {
             setAntiDiffusionFlux_CN(flux_limiter->antidiffusive_fluxes);
         }
         // start the calculation of the limitation factors_{fct_solver->ij}
         flux_limiter->addLimitedFluxes_Start(); // uses u_tilde

         /*
          * z_m[i]=b[i] - (m_i*u[i] - omega*sum_{j<>i} l_{ij} (u[j]-u[i]) ) where m_i>0
          *       ==b[i] - u[i] = u_tilda[i]-u[i] =0 otherwise
          *
          * omega = dt/2 or dt .
          *
          * note that iteration_matrix stores the negative values of the
          * low order transport matrix l. Therefore a=dt*theta is used.
          */
          if (method == PASO_BACKWARD_EULER) {
              setMuPaLu(z, u_coupler, dt);
          } else {
              setMuPaLu(z, u_coupler, dt/2);
          }

        util::update(n, -1., z, 1., b);  // z=b-z

        // z_i += sum_{j} limitation factor_{ij} * antidiffusive_flux_{ij}
        flux_limiter->addLimitedFluxes_Complete(z);

        // we solve (m/omega - L ) * du = z
        if (method == PASO_BACKWARD_EULER) {
            dim_t cntIter = options->iter_max;
            double tol= util::l2(n, z, fctp->mpi_info) ;

            if (m==0) {
                tol *= 0.5;
            } else {
                tol *= std::min(std::max(rate*rate, 1e-2), 0.5);
            }
            // use BiCGStab with Jacobi preconditioner ( m - omega * L )
            util::zeroes(n,du);
            errorCode = Solver_BiCGStab(fctp->iteration_matrix, z, du, &cntIter, &tol, pp);

            // errorCode = Solver_GMRES(fctp->iteration_matrix, z, du, &cntIter, &tol, 10, 2000, pp);
            if (options->verbose)
                std::cout << "FCT_Solver::updateNL: BiCGStab completed after "
                    << cntIter << " steps (residual = " << tol << ")." << std::endl;
            options->num_iter += cntIter;
            if (errorCode != NoError) break;
        } else {
            // just use the main diagonal of (m/omega - L )

            Preconditioner_Smoother_solve(fctp->iteration_matrix,
                ((Preconditioner*) (fctp->iteration_matrix->solver_p))->jacobi,
                du, z, 1, false);

            options->num_iter++;
        }

        util::update(n, 1., u, omega, du);
        norm_du_old = norm_du;
        norm_du = util::lsup(n, du, fctp->mpi_info);
        if (m == 0) {
            if (options->verbose)
                std::cout << "FCT_Solver::updateNL: step " << m+1
                    << ": increment = " << norm_du * omega << std::endl;
        } else {
            if (norm_du_old > 0.) {
                rate = norm_du/norm_du_old;
            } else if (norm_du <= 0.) {
                rate = 0.;
            } else {
                rate = LARGE_POSITIVE_FLOAT;
            }
            if (options->verbose)
                std::cout << "FCT_Solver::updateNL: step " << m+1
                    << ": increment= " << norm_du * omega << " (rate = "
                    << rate << ")" << std::endl;
            num_critical_rates += (rate<critical_rate ? 0 : 1);
            max_m_reached = (m>max_m);
            diverged = (num_critical_rates >= num_critical_rates_max);
            converged = (norm_du * omega <= ATOL);
        }
        m++;
    } // end of while loop
    if (errorCode == NoError) {
        if (converged) {
            if (options->verbose)
                std::cout << "FCT_Solver::updateNL: iteration is completed." << std::endl;
            errorCode = NoError;
        } else if (diverged) {
            if (options->verbose)
                std::cout << "FCT_Solver::updateNL: divergence." << std::endl;
            errorCode = Divergence;
        } else if (max_m_reached) {
            if (options->verbose)
                std::cout << "FCT_Solver::updateNL: maximum number of iteration steps reached." << std::endl;
            errorCode = MaxIterReached;
        }
    }
    return errorCode;
}


/*
 *  AntiDiffusionFlux:
 *
 *        f_{ij} = (m_{ij} - dt (1-theta) d_{ij}) (u_old[j]-u_old[i]) - (m_{ij} + dt theta d_{ij}) (u[j]-u[i])
 *
 *     m=fc->mass matrix
 *     d=artificial diffusion matrix = L - K = - fc->iteration matrix - fc->transport matrix (away from main diagonal)
 *
 *   for CN : theta = 0.5
 *   for BE : theta = 1.
 */

void FCT_Solver::setAntiDiffusionFlux_CN(SystemMatrix_ptr<double> flux_matrix)
{
    const double* u = u_coupler->borrowLocalData();
    const double* u_old = u_old_coupler->borrowLocalData();
    const double* remote_u = u_coupler->borrowRemoteData();
    const double* remote_u_old = u_old_coupler->borrowRemoteData();
    const double dt_half = dt/2;
    const_TransportProblem_ptr fct(transportproblem);
    const_SystemMatrixPattern_ptr pattern(fct->iteration_matrix->pattern);
    const dim_t n = fct->iteration_matrix->getTotalNumRows();

#pragma omp parallel for
    for (dim_t i = 0; i < n; ++i) {
        const double u_i = u[i];
        const double u_old_i = u_old[i];

        #pragma ivdep
        for (index_t iptr_ij = pattern->mainPattern->ptr[i];
                     iptr_ij < pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
            const index_t j   = pattern->mainPattern->index[iptr_ij];
            const double m_ij = fct->mass_matrix->mainBlock->val[iptr_ij];
            // this is in fact -d_ij
            const double d_ij = fct->transport_matrix->mainBlock->val[iptr_ij]+
                fct->iteration_matrix->mainBlock->val[iptr_ij];
            const double u_old_j = u_old[j];
            const double u_j = u[j];

            // (m_{ij} - dt (1-theta) d_{ij}) (u_old[j]-u_old[i]) - (m_{ij} + dt theta d_{ij}) (u[j]-u[i])
            flux_matrix->mainBlock->val[iptr_ij] =
                (m_ij+dt_half*d_ij)*(u_old_j-u_old_i) -
                        (m_ij-dt_half*d_ij)*(u_j-u_i);

        }
        #pragma ivdep
        for (index_t iptr_ij = pattern->col_couplePattern->ptr[i];
                   iptr_ij < pattern->col_couplePattern->ptr[i+1]; iptr_ij++) {
            const index_t j = pattern->col_couplePattern->index[iptr_ij];
            const double m_ij = fct->mass_matrix->col_coupleBlock->val[iptr_ij];
            // this is in fact -d_ij
            const double d_ij =
                fct->transport_matrix->col_coupleBlock->val[iptr_ij] +
                fct->iteration_matrix->col_coupleBlock->val[iptr_ij];
            const double u_old_j = remote_u_old[j];
            const double u_j = remote_u[j];
            flux_matrix->col_coupleBlock->val[iptr_ij] =
                (m_ij+dt_half*d_ij)*(u_old_j-u_old_i) -
                        (m_ij-dt_half*d_ij)*(u_j-u_i);
        }
    }
}

void FCT_Solver::setAntiDiffusionFlux_BE(SystemMatrix_ptr<double> flux_matrix)
{
    const double* u = u_coupler->borrowLocalData();
    const double* u_old = u_old_coupler->borrowLocalData();
    const double* remote_u = u_coupler->borrowRemoteData();
    const double* remote_u_old = u_old_coupler->borrowRemoteData();
    const_TransportProblem_ptr fct(transportproblem);
    const_SystemMatrixPattern_ptr pattern(fct->iteration_matrix->pattern);
    const dim_t n = fct->iteration_matrix->getTotalNumRows();

#pragma omp parallel for
    for (dim_t i = 0; i < n; ++i) {
        const double u_i = u[i];
        const double u_old_i = u_old[i];
        #pragma ivdep
        for (index_t iptr_ij = pattern->mainPattern->ptr[i];
                     iptr_ij < pattern->mainPattern->ptr[i+1]; iptr_ij++) {
            const index_t j = pattern->mainPattern->index[iptr_ij];
            const double m_ij = fct->mass_matrix->mainBlock->val[iptr_ij];
            // this is in fact -d_ij
            const double d_ij = fct->transport_matrix->mainBlock->val[iptr_ij]+
                fct->iteration_matrix->mainBlock->val[iptr_ij];
            const double u_old_j = u_old[j];
            const double u_j = u[j];

            flux_matrix->mainBlock->val[iptr_ij] =
                m_ij*(u_old_j-u_old_i) - (m_ij-dt*d_ij)*(u_j-u_i);
        }
        #pragma ivdep
        for (index_t iptr_ij = pattern->col_couplePattern->ptr[i];
                   iptr_ij < pattern->col_couplePattern->ptr[i+1]; iptr_ij++) {
            const index_t j = pattern->col_couplePattern->index[iptr_ij];
            const double m_ij = fct->mass_matrix->col_coupleBlock->val[iptr_ij];
            // this is in fact -d_ij
            const double d_ij =
                fct->transport_matrix->col_coupleBlock->val[iptr_ij] +
                fct->iteration_matrix->col_coupleBlock->val[iptr_ij];
            const double u_old_j = remote_u_old[j];
            const double u_j = remote_u[j];

            flux_matrix->col_coupleBlock->val[iptr_ij] =
                m_ij*(u_old_j-u_old_i) - (m_ij-dt*d_ij)*(u_j-u_i);
        }
    }
}

/* special version of the ant-diffusive fluxes for the linear Crank-Nicolson
 * scheme. In fact this is evaluated for u = 2*u_tilde - u_old which is the
 * predictor of the solution of the the stabilised problem at time dt using
 * the forward Euler scheme
 *
 * f_{ij} = (m_{ij} - dt/2 d_{ij}) (u_old[j]-u_old[i]) - (m_{ij} + dt/2 d_{ij}) (u[j]-u[i])
 *    =  (m_{ij} - dt/2 d_{ij}) * (u_old[j]-u_old[i]) - (m_{ij} + dt/2 d_{ij}) * ( 2*(u_tilde[j]-u_tilde[i]) - (u_old[j] -u_old [i]) )
 *    =  2*  m_{ij} * ( (u_old[j]-u_tilde[j] - (u_old[i]) - u_tilde[i]) ) - dt d_{ij} * (u_tilde[j]-u_tilde[i])
 *
 */
void FCT_Solver::setAntiDiffusionFlux_linearCN(SystemMatrix_ptr<double> flux_matrix)
{
    const_Coupler_ptr<real_t> u_tilde_coupler(flux_limiter->u_tilde_coupler);
    const double* u_tilde = u_tilde_coupler->borrowLocalData();
    const double* u_old = u_old_coupler->borrowLocalData();
    const double* remote_u_tilde = u_tilde_coupler->borrowRemoteData();
    const double* remote_u_old = u_old_coupler->borrowRemoteData();
    const_TransportProblem_ptr fct(transportproblem);
    const_SystemMatrixPattern_ptr pattern(fct->iteration_matrix->pattern);
    const dim_t n = fct->iteration_matrix->getTotalNumRows();

#pragma omp parallel for
    for (dim_t i = 0; i < n; ++i) {
        const double u_tilde_i = u_tilde[i];
        const double u_old_i = u_old[i];
        const double du_i = u_tilde_i - u_old_i;
        #pragma ivdep
        for (index_t iptr_ij = pattern->mainPattern->ptr[i];
                     iptr_ij < pattern->mainPattern->ptr[i+1]; iptr_ij++) {
            const index_t j = pattern->mainPattern->index[iptr_ij];
            const double m_ij = fct->mass_matrix->mainBlock->val[iptr_ij];
            // this is in fact -d_ij
            const double d_ij = fct->transport_matrix->mainBlock->val[iptr_ij]+
                fct->iteration_matrix->mainBlock->val[iptr_ij];
            const double u_tilde_j = u_tilde[j];
            const double u_old_j = u_old[j];
            const double du_j = u_tilde_j - u_old_j;

            flux_matrix->mainBlock->val[iptr_ij] = 2 * m_ij * ( du_i - du_j ) -
                                          dt * d_ij * ( u_tilde_i - u_tilde_j);
        }
        #pragma ivdep
        for (index_t iptr_ij=pattern->col_couplePattern->ptr[i];
                     iptr_ij<pattern->col_couplePattern->ptr[i+1]; iptr_ij++) {

            const index_t j = pattern->col_couplePattern->index[iptr_ij];
            const double m_ij = fct->mass_matrix->col_coupleBlock->val[iptr_ij];
            // this is in fact -d_ij
            const double d_ij =
                fct->transport_matrix->col_coupleBlock->val[iptr_ij] +
                fct->iteration_matrix->col_coupleBlock->val[iptr_ij];
            const double u_tilde_j = remote_u_tilde[j];
            const double u_old_j = remote_u_old[j];
            const double du_j = u_tilde_j - u_old_j;

            flux_matrix->col_coupleBlock->val[iptr_ij] =
                2*m_ij * ( du_i - du_j ) - dt * d_ij * (u_tilde_i - u_tilde_j);
        }
    }
}

/****************************************************************************/

double FCT_Solver::getSafeTimeStepSize(const_TransportProblem_ptr fctp)
{
    double dt_max = LARGE_POSITIVE_FLOAT;
    const dim_t n = fctp->transport_matrix->getTotalNumRows();

    // set low order transport operator
    setLowOrderOperator(boost::const_pointer_cast<TransportProblem>(fctp));

    // calculate time step size
    dt_max = LARGE_POSITIVE_FLOAT;
#pragma omp parallel
    {
        double dt_max_loc = LARGE_POSITIVE_FLOAT;
#pragma omp for schedule(static)
        for (dim_t i=0; i<n; ++i) {
            const double l_ii = fctp->main_diagonal_low_order_transport_matrix[i];
            const double m_i = fctp->lumped_mass_matrix[i];
            if (m_i > 0) {
                if (l_ii<0)
                    dt_max_loc = std::min(dt_max_loc,m_i/(-l_ii));
            }
        }
        #pragma omp critical
        {
            dt_max = std::min(dt_max,dt_max_loc);
        }
    }
#ifdef ESYS_MPI
    double dt_max_loc = dt_max;
    MPI_Allreduce(&dt_max_loc, &dt_max, 1, MPI_DOUBLE, MPI_MIN, fctp->mpi_info->comm);
#endif
    if (dt_max < LARGE_POSITIVE_FLOAT)
        dt_max *= 2.;

    return dt_max;
}

/* Creates the low order transport matrix and stores its negative values
 * into the iteration_matrix except for the main diagonal which is stored
 * separately.
 * If fc->iteration_matrix==NULL, fc->iteration_matrix is allocated
 *
 * a=transport_matrix
 * b= low_order_transport_matrix = - iteration_matrix
 * c=main diagonal low_order_transport_matrix
 * initialise c[i] mit a[i,i]
 *
 *    d_ij=max(0,-a[i,j],-a[j,i])
 *    b[i,j]=-(a[i,j]+d_ij)
 *    c[i]-=d_ij
 */

void FCT_Solver::setLowOrderOperator(TransportProblem_ptr fc)
{
    const index_t* main_iptr = fc->borrowMainDiagonalPointer();

    if (!fc->iteration_matrix.get()) {
        fc->iteration_matrix.reset(new SystemMatrix<double>(
                  fc->transport_matrix->type, fc->transport_matrix->pattern,
                  fc->transport_matrix->row_block_size,
                  fc->transport_matrix->col_block_size, true,
                  fc->transport_matrix->getRowFunctionSpace(),
                  fc->transport_matrix->getColumnFunctionSpace()));
    }

    const_SystemMatrixPattern_ptr pattern(fc->iteration_matrix->pattern);
    const dim_t n = fc->iteration_matrix->getTotalNumRows();
#pragma omp parallel for
    for (dim_t i = 0; i < n; ++i) {
        double sum = fc->transport_matrix->mainBlock->val[main_iptr[i]];
        //printf("sum[%d] = %e -> ", i, sum);

        // look at a[i,j]
        for (index_t iptr_ij=pattern->mainPattern->ptr[i];iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
            const index_t j = pattern->mainPattern->index[iptr_ij];
            const double rtmp1 = fc->transport_matrix->mainBlock->val[iptr_ij];
            if (j != i) {
                // find entry a[j,i]
                #pragma ivdep
                for (index_t iptr_ji=pattern->mainPattern->ptr[j]; iptr_ji<pattern->mainPattern->ptr[j+1]; ++iptr_ji) {

                    if (pattern->mainPattern->index[iptr_ji] == i) {
                        const double rtmp2=fc->transport_matrix->mainBlock->val[iptr_ji];
                        //printf("a[%d,%d]=%e\n",i,j,rtmp1);
                        //printf("a[%d,%d]=%e\n",j,i,rtmp2);
                        const double d_ij=-MIN3(0.,rtmp1,rtmp2);
                        fc->iteration_matrix->mainBlock->val[iptr_ij]=-(rtmp1+d_ij);
//printf("l[%d,%d]=%e\n",i,j,fc->iteration_matrix->mainBlock->val[iptr_ij]);
                        sum-=d_ij;
                        break;
                    }
                }
            }
        }
        for (index_t iptr_ij=pattern->col_couplePattern->ptr[i];iptr_ij<pattern->col_couplePattern->ptr[i+1]; ++iptr_ij) {
            const index_t j = pattern->col_couplePattern->index[iptr_ij];
            const double rtmp1 = fc->transport_matrix->col_coupleBlock->val[iptr_ij];
            // find entry a[j,i]
            #pragma ivdep
            for (index_t iptr_ji=pattern->row_couplePattern->ptr[j]; iptr_ji<pattern->row_couplePattern->ptr[j+1]; ++iptr_ji) {
                if (pattern->row_couplePattern->index[iptr_ji]==i) {
                    const double rtmp2=fc->transport_matrix->row_coupleBlock->val[iptr_ji];
                    const double d_ij=-MIN3(0.,rtmp1,rtmp2);
                    fc->iteration_matrix->col_coupleBlock->val[iptr_ij]=-(rtmp1+d_ij);
                    fc->iteration_matrix->row_coupleBlock->val[iptr_ji]=-(rtmp2+d_ij);
                    sum-=d_ij;
                    break;
                }
            }
        }
        // set main diagonal entry
        fc->main_diagonal_low_order_transport_matrix[i] = sum;
        //printf("%e\n", sum);
    }
}

/*
 * out_i=m_i u_i + a * \sum_{j <> i} l_{ij} (u_j-u_i) where m_i>0
 *       = u_i                                        where m_i<=0
 *
 */
void FCT_Solver::setMuPaLu(double* out, const_Coupler_ptr<real_t> coupler, double a)
{
    const_SystemMatrix_ptr<double> L(transportproblem->iteration_matrix);
    const double* M = transportproblem->lumped_mass_matrix;
    const_SystemMatrixPattern_ptr pattern(L->pattern);
    const double* u = coupler->borrowLocalData();
    const double* remote_u = coupler->borrowRemoteData();
    const dim_t n = L->getTotalNumRows();

#pragma omp parallel for
    for (dim_t i = 0; i < n; ++i) {
        if (M[i] > 0.) {
            out[i] = M[i]*u[i];
        } else {
            out[i] = u[i];
        }
    }
    if (std::abs(a) > 0) {
#pragma omp parallel for
        for (dim_t i = 0; i < n; ++i) {
            if (M[i] > 0.) {
                double sum = 0;
                const double u_i = u[i];
                #pragma ivdep
                for (index_t iptr_ij = pattern->mainPattern->ptr[i];
                             iptr_ij < pattern->mainPattern->ptr[i+1];
                             iptr_ij++) {
                    const index_t j = pattern->mainPattern->index[iptr_ij];
                    const double l_ij = L->mainBlock->val[iptr_ij];
                    sum += l_ij*(u[j]-u_i);
                }
                #pragma ivdep
                for (index_t iptr_ij = pattern->col_couplePattern->ptr[i];
                             iptr_ij < pattern->col_couplePattern->ptr[i+1];
                             iptr_ij++) {
                    const index_t j=pattern->col_couplePattern->index[iptr_ij];
                    const double l_ij = L->col_coupleBlock->val[iptr_ij];
                    sum += l_ij*(remote_u[j]-u_i);
                }
                out[i] += a*sum;
            }
        }
    }
}

} // namespace paso

