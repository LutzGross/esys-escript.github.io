
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/****************************************************************************/

/* Paso: AMG preconditioner  (local version)                  */

/****************************************************************************/

/* Author: artak@uq.edu.au, l.gross@uq.edu.au                 */

/****************************************************************************/

#define SHOW_TIMING 0

#include "Paso.h"
#include "Preconditioner.h"
#include "MergedSolver.h"
#include "Options.h"
#include "PasoUtil.h"
#include "MKL.h"

#include <iostream>

namespace {

double random_seed = .4142135623730951;

}

namespace paso {

inline double* createRandomVector(escript::const_Distribution_ptr dist)
{
    const index_t n_0 = dist->getFirstComponent();
    const index_t n_1 = dist->getLastComponent();
    const index_t n = dist->getGlobalNumComponents();
    const dim_t my_n = n_1 - n_0;
    double* out = new double[my_n];

#pragma omp parallel for schedule(static)
    for (index_t i = 0; i < my_n; ++i) {
        out[i] = fmod(random_seed * (n_0+i+1), 1.);
    }

    random_seed = fmod(random_seed * (n+1.7), 1.);
    return out;
}


void Preconditioner_AMG_free(Preconditioner_AMG* in)
{
    if (in!=NULL) {
        Preconditioner_Smoother_free(in->Smoother);
        Preconditioner_AMG_free(in->AMG_C);
        delete[] in->r;
        delete[] in->x_C;
        delete[] in->b_C;
        delete in->merged_solver;
        delete in;
    }
}

int Preconditioner_AMG_getMaxLevel(const Preconditioner_AMG* in)
{
    if (in->AMG_C == NULL) {
        return in->level;
    } else {
        return Preconditioner_AMG_getMaxLevel(in->AMG_C);
    }
}

double Preconditioner_AMG_getCoarseLevelSparsity(const Preconditioner_AMG* in)
{
    if (in->AMG_C == NULL) {
        if (in->A_C == NULL) {
            return 1.;
        } else {
            return in->A_C->getSparsity();
        }
    }
    return Preconditioner_AMG_getCoarseLevelSparsity(in->AMG_C);
}

dim_t Preconditioner_AMG_getNumCoarseUnknowns(const Preconditioner_AMG* in)
{
    if (in->AMG_C == NULL) {
        if (in->A_C == NULL) {
            return 0;
        } else {
            return in->A_C->getTotalNumRows();
        }
    }
    return Preconditioner_AMG_getNumCoarseUnknowns(in->AMG_C);
}

/****************************************************************************

   constructs AMG

*****************************************************************************/
Preconditioner_AMG* Preconditioner_AMG_alloc(SystemMatrix_ptr A, int level,
                                             Options* options)
{
    Preconditioner_AMG* out = NULL;
    const bool verbose = options->verbose;
    const dim_t my_n = A->mainBlock->numRows;
    const dim_t overlap_n = A->row_coupleBlock->numRows;
    const dim_t n = my_n + overlap_n;
    const dim_t n_block = A->row_block_size;
    const double sparsity = A->getSparsity();
    const dim_t global_n = A->getGlobalNumRows();

    // is the input matrix A suitable for coarsening?
    if (sparsity >= options->min_coarse_sparsity ||
         global_n <= options->min_coarse_matrix_size ||
         level > options->level_max) {
        if (verbose) {
            /*
            print stopping condition:
                - 'SPAR' = min_coarse_matrix_sparsity exceeded
                - 'SIZE' = min_coarse_matrix_size exceeded
                - 'LEVEL' = level_max exceeded
            */
            std::cout << "Preconditioner: AMG: termination of coarsening by ";

            if (sparsity >= options->min_coarse_sparsity)
                std::cout << "SPAR";

            if (global_n <= options->min_coarse_matrix_size)
                std::cout << "SIZE";

            if (level > options->level_max)
                std::cout << "LEVEL";

            std::cout << std::endl << "Preconditioner: AMG level " << level
                << " (limit = " << options->level_max << ") stopped."
                << " Sparsity = " << sparsity << " (limit = "
                << options->min_coarse_sparsity << "), unknowns = " << global_n
                << " (limit = " << options->min_coarse_matrix_size << ")"
                << std::endl;
        }
        return out;
    }

    ///// Start Coarsening /////

    double time0 = 0;
    // this is the table for strong connections combining mainBlock,
    // col_coupleBlock and row_coupleBlock
    const dim_t len_S = A->mainBlock->pattern->len +
                        A->col_coupleBlock->pattern->len +
                        A->row_coupleBlock->pattern->len +
                        A->row_coupleBlock->numRows*A->col_coupleBlock->numCols;

    dim_t* degree_S = new dim_t[n];
    index_t* offset_S = new index_t[n];
    index_t* S = new index_t[len_S];
    dim_t* degree_ST = new dim_t[n];
    index_t* offset_ST = new index_t[n];
    index_t* ST = new index_t[len_S];
    AMGBlockSelect* F_marker = new AMGBlockSelect[n];
    index_t* counter = new index_t[n];
    const double theta = options->coarsening_threshold;
    const double tau = options->diagonal_dominance_threshold;

    // make sure that corresponding values in the row_coupleBlock and
    // col_coupleBlock are identical
    A->copyColCoupleBlock();
    A->copyRemoteCoupleBlock(false);

    // set splitting of unknowns
    time0 = escript::gettime();
    if (n_block > 1) {
        Preconditioner_AMG_setStrongConnections_Block(A, degree_S, offset_S,
                                                      S, theta, tau);
    } else {
        Preconditioner_AMG_setStrongConnections(A, degree_S, offset_S, S,
                                                theta, tau);
    }
    Preconditioner_AMG_transposeStrongConnections(n, degree_S, offset_S, S,
                                                  n, degree_ST, offset_ST, ST);
    //A->extendedRowsForST(degree_ST, offset_ST, ST);

    Preconditioner_AMG_CIJPCoarsening(n, my_n, F_marker, degree_S, offset_S, S,
                                      degree_ST, offset_ST, ST,
                                      A->col_coupler->connector,
                                      A->col_distribution);

    // in BoomerAMG if interpolation is used FF connectivity is required
    //if (options->interpolation_method == PASO_CLASSIC_INTERPOLATION_WITH_FF_COUPLING)
    //    Preconditioner_AMG_enforceFFConnectivity(n, A->pattern->ptr, degree_S, S, F_marker);

    options->coarsening_selection_time = escript::gettime()-time0 +
                               std::max(0., options->coarsening_selection_time);
#pragma omp parallel for
    for (dim_t i = 0; i < n; ++i)
        F_marker[i] = (F_marker[i]==PASO_AMG_IN_F ? PASO_AMG_IN_C : PASO_AMG_IN_F);

    // count number of unknowns to be eliminated:
    dim_t my_n_F = util::cumsum_maskedTrue(my_n, counter, (int*)F_marker);
    const dim_t n_F = util::cumsum_maskedTrue(n, counter, (int*)F_marker);
    // collect my_n_F values on all processes, a direct solver should
    // be used if any my_n_F value is 0
    dim_t* F_set = new dim_t[A->mpi_info->size];
#ifdef ESYS_MPI
    MPI_Allgather(&my_n_F, 1, MPI_INT, F_set, 1, MPI_INT, A->mpi_info->comm);
#endif
    dim_t global_n_F = 0;
    bool F_flag = true;
    for (dim_t i=0; i<A->mpi_info->size; i++) {
        global_n_F += F_set[i];
        if (F_set[i] == 0)
            F_flag = false;
    }
    delete[] F_set;

    const dim_t global_n_C = global_n-global_n_F;
    if (verbose)
        std::cout << "Preconditioner: AMG (non-local) level " << level
            << ": " << global_n_F << " unknowns are flagged for"
            << " elimination. " << global_n_C << " left." << std::endl;

    //if (n_F == 0) { nasty case. a direct solver should be used!
    if (F_flag) {
        out = new Preconditioner_AMG;
        out->level = level;
        out->post_sweeps = options->post_sweeps;
        out->pre_sweeps  = options->pre_sweeps;
        out->r = NULL;
        out->x_C = NULL;
        out->b_C = NULL;
        out->AMG_C = NULL;
        out->Smoother = NULL;
        out->merged_solver = NULL;
        out->Smoother = Preconditioner_Smoother_alloc(A,
                (options->smoother == PASO_JACOBI), 0, verbose);

        if (global_n_C != 0) {
            index_t* mask_C = new index_t[n];
            index_t* rows_in_F = new index_t[n_F];
            // create mask of C nodes with value >-1, gives new id
            const dim_t n_C = util::cumsum_maskedFalse(n, mask_C,
                                                       (int*)F_marker);
            const dim_t my_n_C = my_n-my_n_F;
            // if nothing has been removed we have a diagonal dominant
            // matrix and we just run a few steps of the smoother

            out->x_C = new double[n_block*my_n_C];
            out->b_C = new double[n_block*my_n_C];
            out->r = new double[n_block*my_n];

            // creates index for F
#pragma omp parallel for
            for (dim_t i = 0; i < n; ++i) {
                if (F_marker[i])
                    rows_in_F[counter[i]] = i;
            }
            // get Prolongation
            time0 = escript::gettime();
            out->P = Preconditioner_AMG_getProlongation(A,
                    offset_S, degree_S, S, n_C, mask_C,
                    options->interpolation_method);

            // construct Restriction operator as transposed of
            // Prolongation operator:
            time0 = escript::gettime();
            out->R = Preconditioner_AMG_getRestriction(out->P);
            if (SHOW_TIMING)
                std::cout << "timing: level " << level
                    << ": getTranspose: " << escript::gettime()-time0
                    << std::endl;
            // construct coarse level matrix
            SystemMatrix_ptr A_C;
            time0 = escript::gettime();
            A_C = Preconditioner_AMG_buildInterpolationOperator(A, out->P, out->R);
            if (SHOW_TIMING)
                std::cout << "timing: level " << level
                    << ": construct coarse matrix: "
                    << escript::gettime()-time0 << std::endl;

            out->AMG_C = Preconditioner_AMG_alloc(A_C, level+1, options);
            out->A_C = A_C;
            if (out->AMG_C == NULL) {
                // merge the system matrix into 1 rank when
                // it's not suitable coarsening due to the
                // total number of unknowns being too small
                out->merged_solver = new MergedSolver(A_C, options);
            }
            delete[] mask_C;
            delete[] rows_in_F;
        }
    }
    delete[] counter;
    delete[] F_marker;
    delete[] degree_S;
    delete[] offset_S;
    delete[] S;
    delete[] degree_ST;
    delete[] offset_ST;
    delete[] ST;

    return out;
}


void Preconditioner_AMG_solve(SystemMatrix_ptr A,
                              Preconditioner_AMG* amg, double* x, double* b)
{
    const dim_t n = A->mainBlock->numRows * A->mainBlock->row_block_size;
    const dim_t post_sweeps=amg->post_sweeps;
    const dim_t pre_sweeps=amg->pre_sweeps;

    // presmoothing
    double time0 = escript::gettime();
    Preconditioner_Smoother_solve(A, amg->Smoother, x, b, pre_sweeps, false);
    time0 = escript::gettime()-time0;
    if (SHOW_TIMING)
        std::cout << "timing: level " << amg->level << ": Presmoothing: "
            << time0 << std::endl;
    // end of presmoothing

    time0=escript::gettime();
    // r <- b
    util::copy(n, amg->r, b);
    // r = r-Ax
    A->MatrixVector_CSR_OFFSET0(-1., x, 1., amg->r);
    // b_C = R*r
    amg->R->MatrixVector_CSR_OFFSET0(1., amg->r, 0., amg->b_C);
    time0 = escript::gettime()-time0;

    // coarse level solve
    if (amg->AMG_C == NULL) {
        time0 = escript::gettime();
        // A_C is the coarsest level
        amg->merged_solver->solve(amg->x_C, amg->b_C);
        if (SHOW_TIMING)
            std::cout << "timing: level " << amg->level << ": DIRECT SOLVER: "
                << escript::gettime()-time0 << std::endl;
    } else {
        // x_C = AMG(b_C)
        Preconditioner_AMG_solve(amg->A_C, amg->AMG_C, amg->x_C, amg->b_C);
    }

    time0 = time0+escript::gettime();
    // x = x + P*x_c
    amg->P->MatrixVector_CSR_OFFSET0(1., amg->x_C, 1., x);

    // postsmoothing:
    // solve Ax=b with initial guess x
    time0 = escript::gettime();
    Preconditioner_Smoother_solve(A, amg->Smoother, x, b, post_sweeps, true);
    time0 = escript::gettime()-time0;
    if (SHOW_TIMING)
        std::cout << "timing: level " << amg->level << ": Postsmoothing: "
            << time0 << std::endl;
}

/// theta = threshold for strong connections
/// tau = threshold for diagonal dominance
///
/// S_i={j \in N_i; i strongly coupled to j}
/// in the sense that |A_{ij}| >= theta * max_k |A_{ik}|
void Preconditioner_AMG_setStrongConnections(SystemMatrix_ptr A,
                                             dim_t* degree_S,
                                             index_t* offset_S, index_t* S,
                                             double theta, double tau)
{
    const dim_t my_n=A->mainBlock->numRows;
    const dim_t overlap_n=A->row_coupleBlock->numRows;
    index_t iptr, i;
    double *threshold_p=NULL;
    threshold_p = new double[2*my_n];

#pragma omp parallel for private(i,iptr) schedule(static)
    for (i=0;i<my_n;++i) {
        double max_offdiagonal = 0.;
        double sum_row=0;
        double main_row=0;
        dim_t kdeg=0;
        const index_t koffset=A->mainBlock->pattern->ptr[i]+A->col_coupleBlock->pattern->ptr[i];

        // collect information for row i:
        #pragma ivdep
        for (iptr=A->mainBlock->pattern->ptr[i];iptr<A->mainBlock->pattern->ptr[i+1]; ++iptr) {
            index_t j=A->mainBlock->pattern->index[iptr];
            const double fnorm=std::abs(A->mainBlock->val[iptr]);
            if(j != i) {
                max_offdiagonal = std::max(max_offdiagonal,fnorm);
                sum_row+=fnorm;
            } else {
                main_row=fnorm;
            }
        }

        #pragma ivdep
        for (iptr=A->col_coupleBlock->pattern->ptr[i];iptr<A->col_coupleBlock->pattern->ptr[i+1]; ++iptr) {
            const double fnorm = std::abs(A->col_coupleBlock->val[iptr]);

            max_offdiagonal = std::max(max_offdiagonal,fnorm);
            sum_row+=fnorm;
        }

        // inspect row i:
        {
            const double threshold = theta*max_offdiagonal;
            threshold_p[2*i+1]=threshold;
            if (tau*main_row < sum_row) { // no diagonal dominance
                threshold_p[2*i]=1;
                #pragma ivdep
                for (iptr=A->mainBlock->pattern->ptr[i];iptr<A->mainBlock->pattern->ptr[i+1]; ++iptr) {
                    const index_t j=A->mainBlock->pattern->index[iptr];
                    if(std::abs(A->mainBlock->val[iptr])>threshold && i!=j) {
                        S[koffset+kdeg] = j;
                        kdeg++;
                    }
                }
                #pragma ivdep
                for (iptr=A->col_coupleBlock->pattern->ptr[i];iptr<A->col_coupleBlock->pattern->ptr[i+1]; ++iptr) {
                    const index_t j=A->col_coupleBlock->pattern->index[iptr];
                    if(std::abs(A->col_coupleBlock->val[iptr])>threshold) {
                        S[koffset+kdeg] = j + my_n;
                        kdeg++;
                    }
                }
            } else {
                threshold_p[2*i]=-1;
            }
        }
        offset_S[i]=koffset;
        degree_S[i]=kdeg;
    }

    // now we need to distribute the threshold and the diagonal dominance
    // indicator
    if (A->mpi_info->size > 1) {
        const index_t koffset_0 =
            A->mainBlock->pattern->ptr[my_n]+A->col_coupleBlock->pattern->ptr[my_n]
            -A->mainBlock->pattern->ptr[0]-A->col_coupleBlock->pattern->ptr[0];

        Coupler_ptr<real_t> threshold_coupler(new Coupler<real_t>(A->row_coupler->connector, 2, A->mpi_info));
        threshold_coupler->startCollect(threshold_p);
        double* remote_threshold = threshold_coupler->finishCollect();

#pragma omp parallel for private(i,iptr) schedule(static)
        for (i=0; i<overlap_n; i++) {
            const double threshold = remote_threshold[2*i+1];
            dim_t kdeg=0;
            const index_t koffset=koffset_0+A->row_coupleBlock->pattern->ptr[i]+A->remote_coupleBlock->pattern->ptr[i];
            if (remote_threshold[2*i]>0) {
                #pragma ivdep
                for (iptr=A->row_coupleBlock->pattern->ptr[i];iptr<A->row_coupleBlock->pattern->ptr[i+1]; ++iptr) {
                  index_t j=A->row_coupleBlock->pattern->index[iptr];
                  if(std::abs(A->row_coupleBlock->val[iptr])>threshold) {
                     S[koffset+kdeg] = j ;
                     kdeg++;
                  }
                }

                #pragma ivdep
                for (iptr=A->remote_coupleBlock->pattern->ptr[i];iptr<A->remote_coupleBlock->pattern->ptr[i+1]; iptr++) {
                  index_t j=A->remote_coupleBlock->pattern->index[iptr];
                  if(std::abs(A->remote_coupleBlock->val[iptr])>threshold && i!=j) {
                      S[koffset+kdeg] = j + my_n;
                      kdeg++;
                  }
                }
            }
            offset_S[i+my_n]=koffset;
            degree_S[i+my_n]=kdeg;
        }
    }
    delete[] threshold_p;
}

// theta = threshold for strong connections
// tau = threshold for diagonal dominance
// S_i={j \in N_i; i strongly coupled to j}
// in the sense that |A_{ij}|_F >= theta * max_k |A_{ik}|_F
void Preconditioner_AMG_setStrongConnections_Block(SystemMatrix_ptr A,
                                                   dim_t* degree_S,
                                                   index_t* offset_S,
                                                   index_t* S,
                                                   double theta, double tau)
{
    const dim_t block_size=A->block_size;
    const dim_t my_n=A->mainBlock->numRows;
    const dim_t overlap_n=A->row_coupleBlock->numRows;
    double* threshold_p = new double[2*my_n];
    index_t iptr, i, bi;

#pragma omp parallel private(i,iptr,bi)
    {
        dim_t max_deg=0;

        #pragma omp for schedule(static)
        for (i=0; i<my_n; ++i)
            max_deg=std::max(max_deg, A->mainBlock->pattern->ptr[i+1]-A->mainBlock->pattern->ptr[i]
                                +A->col_coupleBlock->pattern->ptr[i+1]-A->col_coupleBlock->pattern->ptr[i]);

        double* rtmp = new double[max_deg];

        #pragma omp for schedule(static)
        for (i=0;i<my_n;++i) {
            double max_offdiagonal = 0.;
            double sum_row=0;
            double main_row=0;
            index_t rtmp_offset=-A->mainBlock->pattern->ptr[i];
            dim_t kdeg=0;
            const index_t koffset=A->mainBlock->pattern->ptr[i]+A->col_coupleBlock->pattern->ptr[i];

            /* collect information for row i: */
            for (iptr=A->mainBlock->pattern->ptr[i];iptr<A->mainBlock->pattern->ptr[i+1]; ++iptr) {
                const index_t j=A->mainBlock->pattern->index[iptr];
                double fnorm=0;
                #pragma ivdep
                for(bi=0;bi<block_size;++bi) {
                    const double rtmp2 = A->mainBlock->val[iptr*block_size+bi];
                    fnorm+=rtmp2*rtmp2;
                }
                fnorm=sqrt(fnorm);
                rtmp[iptr+rtmp_offset]=fnorm;

                if( j != i) {
                    max_offdiagonal = std::max(max_offdiagonal,fnorm);
                    sum_row+=fnorm;
                } else {
                    main_row=fnorm;
                }
            }
            rtmp_offset+=A->mainBlock->pattern->ptr[i+1]-A->col_coupleBlock->pattern->ptr[i];
            for (iptr=A->col_coupleBlock->pattern->ptr[i];iptr<A->col_coupleBlock->pattern->ptr[i+1]; ++iptr) {
                double fnorm=0;
                #pragma ivdep
                for(bi=0;bi<block_size;++bi) {
                    const double rtmp2 = A->col_coupleBlock->val[iptr*block_size+bi];
                    fnorm+=rtmp2*rtmp2;
                }
                fnorm=sqrt(fnorm);

                rtmp[iptr+rtmp_offset]=fnorm;
                max_offdiagonal = std::max(max_offdiagonal,fnorm);
                sum_row+=fnorm;
            }

            // inspect row i:
            {
                const double threshold = theta*max_offdiagonal;
                rtmp_offset=-A->mainBlock->pattern->ptr[i];

                threshold_p[2*i+1]=threshold;
                if (tau*main_row < sum_row) { /* no diagonal dominance */
                    threshold_p[2*i]=1;
                    #pragma ivdep
                    for (iptr=A->mainBlock->pattern->ptr[i];iptr<A->mainBlock->pattern->ptr[i+1]; ++iptr) {
                        const index_t j=A->mainBlock->pattern->index[iptr];
                        if(rtmp[iptr+rtmp_offset] > threshold && i!=j) {
                            S[koffset+kdeg] = j;
                            kdeg++;
                        }
                    }
                    rtmp_offset+=A->mainBlock->pattern->ptr[i+1]-A->col_coupleBlock->pattern->ptr[i];
                    #pragma ivdep
                    for (iptr=A->col_coupleBlock->pattern->ptr[i];iptr<A->col_coupleBlock->pattern->ptr[i+1]; ++iptr) {
                        const index_t j=A->col_coupleBlock->pattern->index[iptr];
                        if( rtmp[iptr+rtmp_offset] >threshold) {
                            S[koffset+kdeg] = j + my_n;
                            kdeg++;
                        }
                    }
                } else {
                    threshold_p[2*i]=-1;
                }
            }
            offset_S[i]=koffset;
            degree_S[i]=kdeg;
        }
        delete[] rtmp;
    } // parallel section

    // now we need to distribute the threshold and the diagonal dominance
    // indicator
    if (A->mpi_info->size > 1) {
        const index_t koffset_0 =
            A->mainBlock->pattern->ptr[my_n]+A->col_coupleBlock->pattern->ptr[my_n]
            -A->mainBlock->pattern->ptr[0]-A->col_coupleBlock->pattern->ptr[0];

        Coupler_ptr<real_t> threshold_coupler(new Coupler<real_t>(A->row_coupler->connector, 2, A->mpi_info));
        threshold_coupler->startCollect(threshold_p);
        double* remote_threshold = threshold_coupler->finishCollect();

        #pragma omp parallel for private(i,iptr) schedule(static)
        for (i=0; i<overlap_n; i++) {
            const double threshold2 = remote_threshold[2*i+1]*remote_threshold[2*i+1];
            dim_t kdeg=0;
            const index_t koffset = koffset_0+A->row_coupleBlock->pattern->ptr[i]+A->remote_coupleBlock->pattern->ptr[i];
            if (remote_threshold[2*i]>0) {
                #pragma ivdep
                for (iptr=A->row_coupleBlock->pattern->ptr[i];iptr<A->row_coupleBlock->pattern->ptr[i+1]; ++iptr) {
                    const index_t j=A->row_coupleBlock->pattern->index[iptr];
                    double fnorm2=0;
                    #pragma ivdep
                    for(bi=0;bi<block_size;++bi) {
                        const double rtmp2 = A->row_coupleBlock->val[iptr*block_size+bi];
                        fnorm2+=rtmp2*rtmp2;
                    }

                    if(fnorm2 > threshold2 ) {
                        S[koffset+kdeg] = j ;
                        kdeg++;
                    }
                }

                #pragma ivdep
                for (iptr=A->remote_coupleBlock->pattern->ptr[i];iptr<A->remote_coupleBlock->pattern->ptr[i+1]; ++iptr) {
                    const index_t j=A->remote_coupleBlock->pattern->index[iptr];
                    double fnorm2 = 0.;
                    #pragma ivdep
                    for (bi=0; bi<block_size; ++bi) {
                        const double v = A->remote_coupleBlock->val[iptr*block_size+bi];
                        fnorm2 += v*v;
                    }
                    if (fnorm2 > threshold2 && i != j) {
                        S[koffset+kdeg] = j + my_n;
                        kdeg++;
                    }
                }
            }
            offset_S[i+my_n]=koffset;
            degree_S[i+my_n]=kdeg;
        }
    }
    delete[] threshold_p;
}

void Preconditioner_AMG_transposeStrongConnections(dim_t n,
        const dim_t* degree_S, const index_t* offset_S, const index_t* S,
        const dim_t nT, dim_t* degree_ST, index_t* offset_ST, index_t* ST)
{
#pragma omp parallel for
    for (index_t i=0; i<nT; ++i) {
        degree_ST[i]=0;
    }
    for (index_t i=0; i<n ;++i) {
        for (dim_t p=0; p<degree_S[i]; ++p)
            degree_ST[S[offset_S[i]+p]]++;
    }
    dim_t len=0;
    for (index_t i=0; i<nT; ++i) {
        offset_ST[i]=len;
        len+=degree_ST[i];
        degree_ST[i]=0;
    }
    for (index_t i=0; i<n; ++i) {
        for (dim_t p=0; p<degree_S[i]; ++p) {
            const index_t j = S[offset_S[i]+p];
            ST[offset_ST[j]+degree_ST[j]]=i;
            degree_ST[j]++;
        }
    }
}

void Preconditioner_AMG_CIJPCoarsening(dim_t n, dim_t my_n,
                                       AMGBlockSelect* split_marker,
                                       const dim_t* degree_S,
                                       const index_t* offset_S,
                                       const index_t* S,
                                       const dim_t* degree_ST,
                                       const index_t* offset_ST,
                                       const index_t* ST,
                                       const_Connector_ptr col_connector,
                                       escript::const_Distribution_ptr col_dist)
{
    Coupler_ptr<real_t> w_coupler(new Coupler<real_t>(col_connector, 1, col_dist->mpi_info));
    double* w = new double[n];
    double* Status = new double[n];
    double* random = createRandomVector(col_dist);
    index_t* ST_flag = new index_t[offset_ST[n-1] + degree_ST[n-1]];
    dim_t i, numUndefined, iter=0;
    index_t iptr, jptr, kptr;

    #pragma omp parallel for private(i)
    for (i=0; i<my_n; ++i) {
        w[i]=degree_ST[i]+random[i];
        if (degree_ST[i] < 1) {
            Status[i]=-100; // F point
        } else {
            Status[i]=1; // status undefined
        }
    }

    #pragma omp parallel for private(i, iptr)
    for (i=0; i<n; ++i) {
        for (iptr =0 ; iptr < degree_ST[i]; ++iptr)  {
            ST_flag[offset_ST[i]+iptr] = 1;
        }
    }

    numUndefined = util::numPositives(col_dist->getMyNumComponents(), Status,
                                      col_dist->mpi_info);
    //printf("coarsening loop start: num of undefined rows = %d \n",numUndefined);
    iter=0;
    while (numUndefined > 0) {
        w_coupler->fillOverlap(n, w);

        // calculate the maximum value of neighbours following active strong
        // connections:
        //    w2[i]=MAX(w[k]) with k in ST[i] or S[i] and (i,k) connection
        //    is still active
        #pragma omp parallel for private(i, iptr)
        for (i=0; i<my_n; ++i) {
            if (Status[i] > 0) { // status is still undefined
                bool inD = true;
                const double wi=w[i];

                for (iptr=0; iptr < degree_S[i]; ++iptr) {
                    const index_t k=S[offset_S[i]+iptr];
                    const index_t* start_p = &ST[offset_ST[k]];
                    const index_t* where_p=(index_t*)bsearch(&i, start_p, degree_ST[k], sizeof(index_t), util::comparIndex);

                    if (ST_flag[offset_ST[k] + (index_t)(where_p-start_p)]>0) {
                        if (wi <= w[k]) {
                            inD = false;
                            break;
                        }
                    }
                }

                if (inD) {
                    for (iptr=0; iptr < degree_ST[i]; ++iptr) {
                        const index_t k=ST[offset_ST[i]+iptr];
                        if (ST_flag[offset_ST[i]+iptr] > 0) {
                            if (wi <= w[k]) {
                                inD = false;
                                break;
                            }
                        }
                    }
                }
                if (inD) {
                    Status[i]=0.; // is in D
                }
            }
        }

        w_coupler->fillOverlap(n, Status);

        /*   remove connection to D points :
               for each i in D:
                  for each j in S_i:
                     w[j]--
                     ST_tag[j,i]=-1
                  for each j in ST[i]:
                     ST_tag[i,j]=-1
                     for each k in ST[j]:
                        if k in ST[i]:
                           w[j]--;
                        ST_tag[j,k]=-1
        */
        // w is updated for local rows only
        {
            #pragma omp parallel for private(i, jptr)
            for (i=0; i<my_n; ++i) {
                for (jptr=0; jptr<degree_ST[i]; ++jptr) {
                    const index_t j=ST[offset_ST[i]+jptr];
                    if (Status[j] == 0. && ST_flag[offset_ST[i]+jptr] > 0) {
                        w[i]--;
                        ST_flag[offset_ST[i]+jptr] = -1;
                    }
                }
            }
            #pragma omp parallel for private(i, jptr)
            for (i=my_n; i<n; ++i) {
                for (jptr=0; jptr<degree_ST[i]; ++jptr) {
                    const index_t j = ST[offset_ST[i]+jptr];
                    if (Status[j] == 0.)
                        ST_flag[offset_ST[i]+jptr] = -1;
                }
            }

            for (i=0; i< n; ++i) {
                if (Status[i] == 0.) {
                    const index_t* start_p = &ST[offset_ST[i]];

                    for (jptr=0; jptr<degree_ST[i]; ++jptr) {
                        const index_t j=ST[offset_ST[i]+jptr];
                        ST_flag[offset_ST[i]+jptr]=-1;
                        for (kptr=0; kptr<degree_ST[j]; ++kptr) {
                            const index_t k=ST[offset_ST[j]+kptr];
                            // k in ST[i]?
                            if (bsearch(&k, start_p, degree_ST[i],
                                        sizeof(index_t), util::comparIndex)) {
                                if (ST_flag[offset_ST[j]+kptr] > 0) {
                                    if (j< my_n) {
                                        w[j]--;
                                    }
                                    ST_flag[offset_ST[j]+kptr]=-1;
                                }
                            }
                        }
                    }
                }
            }
        }

        // adjust status
        #pragma omp parallel for private(i)
        for (i=0; i< my_n; ++i) {
            if ( Status[i] == 0. ) {
                Status[i] = -10;   // this is now a C point
            } else if (Status[i] == 1. && w[i]<1.) {
                Status[i] = -100;  // this is now a F point
            }
        }

        i = numUndefined;
        numUndefined = util::numPositives(col_dist->getMyNumComponents(),
                                          Status, col_dist->mpi_info);
        if (numUndefined == i) {
            throw PasoException("AMG: Can NOT reduce numUndefined.");
        }

        iter++;
        //printf("coarsening loop %d: num of undefined rows = %d \n",iter, numUndefined);

    } // end of while loop

    // map to output
    w_coupler->fillOverlap(n, Status);
    #pragma omp parallel for private(i)
    for (i=0; i< n; ++i) {
        if (Status[i] > -50.) {
            split_marker[i]=PASO_AMG_IN_C;
        } else {
            split_marker[i]=PASO_AMG_IN_F;
        }
    }
    delete[] random;
    delete[] w;
    delete[] Status;
    delete[] ST_flag;
}

} // namespace paso

