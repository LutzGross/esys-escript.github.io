
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

/* Author: artak@uq.edu.au, l.gross@uq.edu.au l.gao@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "Options.h"
#include "PasoUtil.h"
#include "Preconditioner.h"
#include "MKL.h"
#include "UMFPACK.h"

#include <iostream>

#define SHOW_TIMING false
#define USE_TRANSPOSE true
#define SMALL_PANEL true

namespace paso {

void Preconditioner_LocalAMG_free(Preconditioner_LocalAMG* in)
{
    if (in!=NULL) {
        Preconditioner_LocalSmoother_free(in->Smoother);
        Preconditioner_LocalAMG_free(in->AMG_C);
        delete[] in->r;
        delete[] in->x_C;
        delete[] in->b_C;
        delete in;
    }
}

int Preconditioner_LocalAMG_getMaxLevel(const Preconditioner_LocalAMG* in)
{
    if (in->AMG_C == NULL) {
        return in->level;
    }
    return Preconditioner_LocalAMG_getMaxLevel(in->AMG_C);
}

double Preconditioner_LocalAMG_getCoarseLevelSparsity(const Preconditioner_LocalAMG* in)
{
    if (in->AMG_C == NULL) {
        if (in->A_C == NULL) {
            return 1.;
        } else {
            return in->A_C->getSparsity();
        }
    }
    return Preconditioner_LocalAMG_getCoarseLevelSparsity(in->AMG_C);
}

dim_t Preconditioner_LocalAMG_getNumCoarseUnknowns(const Preconditioner_LocalAMG* in)
{
    if (in->AMG_C == NULL) {
        if (in->A_C == NULL) {
            return 0;
        } else {
            return in->A_C->getTotalNumRows();
        }
   }
   return Preconditioner_LocalAMG_getNumCoarseUnknowns(in->AMG_C);
}

/*****************************************************************************

   Constructs AMG

******************************************************************************/
Preconditioner_LocalAMG* Preconditioner_LocalAMG_alloc(SparseMatrix_ptr A_p,
                                                 int level, Options* options)
{
    Preconditioner_LocalAMG* out=NULL;
    bool verbose=options->verbose;

    SparseMatrix_ptr A_C;
    const dim_t n = A_p->numRows;
    const dim_t n_block = A_p->row_block_size;
    AMGBlockSelect* F_marker = NULL;
    index_t *counter=NULL, *mask_C=NULL, *rows_in_F=NULL, *S=NULL, *degree_S=NULL;
    dim_t n_F=0, n_C=0, i;
    double time0=0;
    const double theta = options->coarsening_threshold;
    const double tau = options->diagonal_dominance_threshold;
    const double sparsity = A_p->getSparsity();
    const dim_t total_n = A_p->getTotalNumRows();

    // is the input matrix A suitable for coarsening
    if ( (sparsity >= options->min_coarse_sparsity) ||
            (total_n <= options->min_coarse_matrix_size) ||
            (level > options->level_max) ) {

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
            else if (total_n <= options->min_coarse_matrix_size)
                std::cout << "SIZE";
            else if (level > options->level_max)
                std::cout << "LEVEL";

            std::cout << std::endl << "Preconditioner: AMG level " << level
                << " (limit = " << options->level_max << ") stopped. "
                   "Sparsity = " << sparsity << " (limit = "
                << options->min_coarse_sparsity << "), unknowns = " << total_n
                << " (limit = " << options->min_coarse_matrix_size << ")"
                << std::endl;
        }
        return NULL;
    }

    // Start Coarsening

    F_marker = new AMGBlockSelect[n];
    counter = new index_t[n];
    degree_S = new dim_t[n];
    S = new index_t[A_p->pattern->len];
    /*
         set splitting of unknowns:
    */
    time0=escript::gettime();
    if (n_block>1) {
        Preconditioner_LocalAMG_setStrongConnections_Block(A_p, degree_S, S, theta,tau);
    } else {
        Preconditioner_LocalAMG_setStrongConnections(A_p, degree_S, S, theta,tau);
    }

    Preconditioner_LocalAMG_RungeStuebenSearch(n, A_p->pattern->ptr, degree_S, S, F_marker, options->usePanel);

    /* in BoomerAMG if interpolation is used FF connectivity is required: */
    if (options->interpolation_method == PASO_CLASSIC_INTERPOLATION_WITH_FF_COUPLING)
        Preconditioner_LocalAMG_enforceFFConnectivity(n, A_p->pattern->ptr, degree_S, S, F_marker);
    options->coarsening_selection_time=escript::gettime()-time0 + std::max(0., options->coarsening_selection_time);

    #pragma omp parallel for private(i) schedule(static)
    for (i = 0; i < n; ++i) F_marker[i]=((F_marker[i] ==  PASO_AMG_IN_F) ? PASO_AMG_IN_C : PASO_AMG_IN_F);

    /*
       count number of unknowns to be eliminated:
    */
    n_F=util::cumsum_maskedTrue(n,counter, (int*)F_marker);
    n_C=n-n_F;
    if (verbose)
        std::cout << "Preconditioner: AMG level " << level << ": "
            << n_F << " unknowns are flagged for elimination. "
            << n-n_F << " left." << std::endl;

    if ( n_F == 0 ) {  /* This is a nasty case. A direct solver should be used, return NULL */
        out = NULL;
    } else {
        out = new Preconditioner_LocalAMG;
        out->level = level;
        out->post_sweeps = options->post_sweeps;
        out->pre_sweeps  = options->pre_sweeps;
        out->r = NULL;
        out->x_C = NULL;
        out->b_C = NULL;
        out->AMG_C = NULL;
        out->Smoother=NULL;
        mask_C=new index_t[n];
        rows_in_F=new index_t[n_F];
        out->Smoother = Preconditioner_LocalSmoother_alloc(A_p, (options->smoother == PASO_JACOBI), verbose);

        if (n_C != 0) {
            /* if nothing has been removed we have a diagonal
             * dominant matrix and we just run a few steps of
             * the smoother */
            /* allocate helpers :*/
            out->x_C = new double[n_block*n_C];
            out->b_C = new double[n_block*n_C];
            out->r   = new double[n_block*n];

            /* creates index for F */
#pragma omp parallel private(i)
            {
#pragma omp for schedule(static)
                for (i = 0; i < n; ++i) {
                    if (F_marker[i])
                        rows_in_F[counter[i]]=i;
                }
            }
            // create mask of C nodes with value >-1 gives new id
            i=util::cumsum_maskedFalse(n, counter, (int*)F_marker);

#pragma omp parallel for private(i) schedule(static)
            for (i = 0; i < n; ++i) {
                if (F_marker[i]) {
                    mask_C[i]=-1;
                } else {
                    mask_C[i]=counter[i];;
                }
            }
            /*
              get Prolongation :
            */
            time0=escript::gettime();
            out->P=Preconditioner_LocalAMG_getProlongation(A_p,A_p->pattern->ptr, degree_S,S,n_C,mask_C, options->interpolation_method);
            if (SHOW_TIMING)
                std::cout << "timing: level " << level <<
                   ": getProlongation: " << escript::gettime()-time0
                   << std::endl;
            /*
               construct Restriction operator as transposed of Prolongation operator:
            */
            time0=escript::gettime();
            out->R = out->P->getTranspose();
            if (SHOW_TIMING)
                std::cout << "timing: level " << level
                    << ": SparseMatrix::getTranspose: "
                    << escript::gettime()-time0 << std::endl;
            /*
            construct coarse level matrix:
            */
            SparseMatrix_ptr Atemp;
            time0=escript::gettime();
            if (USE_TRANSPOSE)
                Atemp = SparseMatrix_MatrixMatrixTranspose(A_p,out->P,out->R);
            else
                Atemp = SparseMatrix_MatrixMatrix(A_p,out->P);
            A_C=SparseMatrix_MatrixMatrix(out->R, Atemp);
            if (SHOW_TIMING)
                std::cout << "timing: level " << level
                       << ": construct coarse matrix: "
                       << escript::gettime()-time0 << std::endl;

            /*
               construct coarser level:
            */
            out->AMG_C=Preconditioner_LocalAMG_alloc(A_C,level+1,options);
            if ( out->AMG_C == NULL ) {
                out->reordering = options->reordering;
                out->refinements = options->coarse_matrix_refinements;
                // no coarse level matrix has been constructed.
                // Use direct solver
#ifdef ESYS_HAVE_MKL
                out->A_C = A_C->unroll(MATRIX_FORMAT_BLK1 + MATRIX_FORMAT_OFFSET1);
                A_C.reset();
                out->A_C->solver_package = PASO_MKL;
                if (verbose)
                    std::cout << "Preconditioner: AMG: use MKL "
                      << "direct solver on the coarsest level "
                      << "(number of unknowns = "
                      << n_C*n_block << ")." << std::endl;
#elif defined ESYS_HAVE_UMFPACK
                out->A_C = A_C->unroll(MATRIX_FORMAT_BLK1 + MATRIX_FORMAT_CSC);
                A_C.reset();
                out->A_C->solver_package = PASO_UMFPACK;
                if (verbose)
                    std::cout << "Preconditioner: AMG: use "
                      << "UMFPACK direct solver on the "
                      << "coarsest level (number of unknowns = "
                      << n_C*n_block << ")." << std::endl;
#else
                out->A_C = A_C;
                out->A_C->solver_p = Preconditioner_LocalSmoother_alloc(out->A_C, (options->smoother == PASO_JACOBI), verbose);
                out->A_C->solver_package = PASO_SMOOTHER;
                if (verbose)
                    std::cout << "Preconditioner: AMG: use "
                      << "smoother on the coarsest level "
                      << "(number of unknowns = "
                      << n_C*n_block << ")." << std::endl;
#endif
            } else {
                // finally we set some helpers for the solver step
                out->A_C = A_C;
            }
        }
        delete[] mask_C;
        delete[] rows_in_F;
    }
    delete[] counter;
    delete[] F_marker;
    delete[] degree_S;
    delete[] S;
    return out;
}


void Preconditioner_LocalAMG_solve(SparseMatrix_ptr A,
                                   Preconditioner_LocalAMG* amg,
                                   double* x, double* b)
{
    const dim_t n = A->numRows * A->row_block_size;
    double time0=0;
    const dim_t post_sweeps=amg->post_sweeps;
    const dim_t pre_sweeps=amg->pre_sweeps;

    // presmoothing
    time0=escript::gettime();
    Preconditioner_LocalSmoother_solve(A, amg->Smoother, x, b, pre_sweeps, false);
    time0=escript::gettime()-time0;
    if (SHOW_TIMING)
        std::cout << "timing: level " << amg->level << ": Presmoothing: "
             << time0 << std::endl;;
    // end of presmoothing

    time0=escript::gettime();
    util::copy(n, amg->r, b);                            /*  r <- b */
    SparseMatrix_MatrixVector_CSR_OFFSET0(-1.,A,x,1.,amg->r); /*r=r-Ax*/
    SparseMatrix_MatrixVector_CSR_OFFSET0_DIAG(1.,amg->R,amg->r,0.,amg->b_C);  /* b_c = R*r  */
    time0=escript::gettime()-time0;
    /* coarse level solve */
    if (amg->AMG_C == NULL) {
            time0=escript::gettime();
            /*  A_C is the coarsest level */
            switch (amg->A_C->solver_package) {
               case (PASO_MKL):
                  MKL_solve(amg->A_C, amg->x_C,amg->b_C, amg->reordering, amg->refinements, SHOW_TIMING);
                  break;
               case (PASO_UMFPACK):
                  UMFPACK_solve(amg->A_C, amg->x_C,amg->b_C, amg->refinements, SHOW_TIMING);
                  break;
               case (PASO_SMOOTHER):
                  Preconditioner_LocalSmoother_solve(amg->A_C, reinterpret_cast<Preconditioner_LocalSmoother*>(amg->A_C->solver_p),amg->x_C,amg->b_C,pre_sweeps+post_sweeps, false);
                  break;
            }
            if (SHOW_TIMING)
                std::cout << "timing: level " << amg->level
                    << ": DIRECT SOLVER: " << escript::gettime()-time0 << std::endl;
    } else {
            Preconditioner_LocalAMG_solve(amg->A_C, amg->AMG_C,amg->x_C,amg->b_C); /* x_C=AMG(b_C)     */
    }
    time0=time0+escript::gettime();
    SparseMatrix_MatrixVector_CSR_OFFSET0_DIAG(1.,amg->P,amg->x_C,1.,x); /* x = x + P*x_c */

    /*postsmoothing*/

    /*solve Ax=b with initial guess x */
    time0=escript::gettime();
    Preconditioner_LocalSmoother_solve(A, amg->Smoother, x, b, post_sweeps, true);
    time0=escript::gettime()-time0;
    if (SHOW_TIMING)
         std::cout << "timing: level " << amg->level << ": Postsmoothing: "
             << time0 << std::endl;
    // end of postsmoothing
}

/* theta = threshold for strong connections */
/* tau = threshold for diagonal dominance */

/*S_i={j \in N_i; i strongly coupled to j}

in the sense that |A_{ij}| >= theta * max_k |A_{ik}|
*/

void Preconditioner_LocalAMG_setStrongConnections(SparseMatrix_ptr A,
                                          dim_t *degree_S, index_t *S,
                                          const double theta, const double tau)
{
    const dim_t n=A->numRows;
    index_t iptr, i;

#pragma omp parallel for private(i,iptr) schedule(static)
    for (i=0;i<n;++i) {
        double max_offdiagonal = 0.;
        double sum_row=0;
        double main_row=0;
        #pragma ivdep
        for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
            index_t j=A->pattern->index[iptr];
            const double fnorm=std::abs(A->val[iptr]);

            if(j != i) {
                max_offdiagonal = std::max(max_offdiagonal,fnorm);
                sum_row+=fnorm;
            } else {
                main_row=fnorm;
            }
        }
        {
            const double threshold = theta*max_offdiagonal;
            dim_t kdeg=0;
            if (tau*main_row < sum_row) { /* no diagonal dominance */
                #pragma ivdep
                for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
                    index_t j=A->pattern->index[iptr];
                    if (std::abs(A->val[iptr])>threshold && i!=j) {
                        S[A->pattern->ptr[i]+kdeg] = j;
                        kdeg++;
                    }
                }
            }
            degree_S[i]=kdeg;
        }
    }
}

/* theta = threshold for strong connections */
/* tau = threshold for diagonal dominance */
/*S_i={j \in N_i; i strongly coupled to j}

in the sense that |A_{ij}|_F >= theta * max_k |A_{ik}|_F
*/
void Preconditioner_LocalAMG_setStrongConnections_Block(SparseMatrix_ptr A,
                                                        dim_t* degree_S,
                                                        index_t* S,
                                                        double theta,
                                                        double tau)
{
    const dim_t n_block=A->row_block_size;
    const dim_t n=A->numRows;
    index_t iptr, i, bi;

#pragma omp parallel private(i, iptr,  bi)
    {
         dim_t max_deg=0; /* this is local on each thread */
         double* rtmp=NULL;

         #pragma omp for schedule(static)
         for (i=0;i<n;++i) max_deg=std::max(max_deg, A->pattern->ptr[i+1]-A->pattern->ptr[i]);

         rtmp=new double[max_deg];

         #pragma omp for schedule(static)
         for (i=0;i<n;++i) {

            double max_offdiagonal = 0.;
            double sum_row=0;
            double main_row=0;

            for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
               index_t j=A->pattern->index[iptr];
               double fnorm=0;
               #pragma ivdep
               for(bi=0;bi<n_block*n_block;++bi) {
                   const double rtmp2 = A->val[iptr*n_block*n_block+bi];
                   fnorm+=rtmp2*rtmp2;
               }
               fnorm=sqrt(fnorm);

               rtmp[iptr-A->pattern->ptr[i]]=fnorm;
               if( j != i) {
                  max_offdiagonal = std::max(max_offdiagonal,fnorm);
                  sum_row+=fnorm;
               } else {
                  main_row=fnorm;
               }
            }
            {
               const double threshold = theta*max_offdiagonal;
               dim_t kdeg=0;
               if (tau*main_row < sum_row) { /* no diagonal dominance */
                  #pragma ivdep
                  for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
                     index_t j=A->pattern->index[iptr];
                     if(rtmp[iptr-A->pattern->ptr[i]] > threshold && i!=j) {
                        S[A->pattern->ptr[i]+kdeg] = j;
                        kdeg++;
                     }
                  }
               }
               degree_S[i]=kdeg;
            }
        }
        delete[] rtmp;
    } // end of parallel region
}

/* the Runge Stueben coarsening algorithm: */
void Preconditioner_LocalAMG_RungeStuebenSearch(dim_t n,
        const index_t* offset_S, const dim_t* degree_S, const index_t* S,
        AMGBlockSelect* split_marker, bool usePanel)
{
    bool* notInPanel=NULL;
    index_t *lambda=NULL, *ST=NULL, *panel=NULL, lambda_max, lambda_k;
    dim_t i,k, p, q, *degree_ST=NULL, len_panel, len_panel_new;
    index_t j, itmp;

    // make sure that the return of util::arg_max is not pointing to nirvana
    if (n<=0)
        return;

    lambda=new index_t[n];
    degree_ST=new dim_t[n];
    ST=new index_t[offset_S[n]];
    if (usePanel) {
        if (!SMALL_PANEL) {
            notInPanel=new bool[n];
        }
        panel=new index_t[n];
    }

    /* initialize split_marker: */
    /* Those unknowns which are not influenced go into F, the rest is available for F or C */
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<n;++i) {
        degree_ST[i]=0;
        if (degree_S[i]>0) {
            lambda[i]=0;
            split_marker[i]=PASO_AMG_UNDECIDED;
        } else {
            split_marker[i]=PASO_AMG_IN_F;
            lambda[i]=-1;
        }
    }
    /* create transpose :*/
    for (i=0;i<n;++i) {
        for (p=0; p<degree_S[i]; ++p) {
           j=S[offset_S[i]+p];
           ST[offset_S[j]+degree_ST[j]]=i;
           degree_ST[j]++;
        }
    }

    /* lambda[i] = |undecided k in ST[i]| + 2 * |F-unknown in ST[i]| */
    #pragma omp parallel for private(i, j, p, itmp) schedule(static)
    for (i=0;i<n;++i) {
        if (split_marker[i]==PASO_AMG_UNDECIDED) {
            itmp=lambda[i];
            for (p=0; p<degree_ST[i]; ++p) {
                j=ST[offset_S[i]+p];
                if (split_marker[j]==PASO_AMG_UNDECIDED) {
                    itmp++;
                } else {  /* at this point there are no C points */
                    itmp+=2;
                }
            }
            lambda[i]=itmp;
        }
    }
    if (usePanel && !SMALL_PANEL) {
        #pragma omp parallel for private(i) schedule(static)
        for (i=0;i<n;++i)
            notInPanel[i]=true;
    }

    // start search
    i=util::arg_max(n,lambda);
    while (lambda[i] > -1) { // is there any undecided unknown?
        if (SMALL_PANEL) {
            do {
                len_panel=0;
                // the unknown i is moved to C
                split_marker[i]=PASO_AMG_IN_C;
                lambda[i]=-1; // lambda from unavailable unknowns is set to -1

                // all undecided unknowns strongly coupled to i are moved to F
                for (p=0; p<degree_ST[i]; ++p) {
                    j=ST[offset_S[i]+p];

                    if (split_marker[j]==PASO_AMG_UNDECIDED) {
                        split_marker[j]=PASO_AMG_IN_F;
                        lambda[j]=-1;
                        for (q=0; q<degree_S[j]; ++q) {
                            k=S[offset_S[j]+q];
                            if (split_marker[k]==PASO_AMG_UNDECIDED) {
                                lambda[k]++;
                                panel[len_panel]=k;
                                len_panel++;
                            }
                        }
                    }
                }
                for (p=0; p<degree_S[i]; ++p) {
                    j=S[offset_S[i]+p];
                    if (split_marker[j]==PASO_AMG_UNDECIDED) {
                        lambda[j]--;
                        panel[len_panel]=j;
                        len_panel++;
                    }
                }

                lambda_max=-1;
                for (q=0; q<len_panel; q++) {
                    k = panel[q];
                    j = lambda[k];
                    if (lambda_max < j) {
                        lambda_max = j;
                        i = k;
                    }
                }
            } while (len_panel>0);
        } else if (usePanel) {
            len_panel=0;
            do {
                // the unknown i is moved to C
                split_marker[i]=PASO_AMG_IN_C;
                lambda[i]=-1; // lambda from unavailable unknowns is set to -1

                // all undecided unknowns strongly coupled to i are moved to F
                for (p=0; p<degree_ST[i]; ++p) {
                    j=ST[offset_S[i]+p];
                    if (split_marker[j]==PASO_AMG_UNDECIDED) {
                        split_marker[j]=PASO_AMG_IN_F;
                        lambda[j]=-1;
                        for (q=0; q<degree_S[j]; ++q) {
                            k=S[offset_S[j]+q];
                            if (split_marker[k]==PASO_AMG_UNDECIDED) {
                                lambda[k]++;
                                if (notInPanel[k]) {
                                    notInPanel[k]=false;
                                    panel[len_panel]=k;
                                    len_panel++;
                                }
                            } // the unknown i is moved to C
                            split_marker[i]=PASO_AMG_IN_C;
                            lambda[i]=-1; // lambda from unavailable unknowns is set to -1
                        }
                    }
                }
                for (p=0; p<degree_S[i]; ++p) {
                    j=S[offset_S[i]+p];
                    if (split_marker[j]==PASO_AMG_UNDECIDED) {
                        lambda[j]--;
                        if (notInPanel[j]) {
                            notInPanel[j]=false;
                            panel[len_panel]=j;
                            len_panel++;
                        }
                    }
                }

                // consolidate panel
                // remove lambda[q]=-1
                lambda_max=-1;
                i=-1;
                len_panel_new=0;
                for (q=0; q<len_panel; q++) {
                    k=panel[q];
                    lambda_k=lambda[k];
                    if (split_marker[k]==PASO_AMG_UNDECIDED) {
                        panel[len_panel_new]=k;
                        len_panel_new++;

                        if (lambda_max == lambda_k) {
                            if (k<i) i=k;
                        } else if (lambda_max < lambda_k) {
                            lambda_max =lambda_k;
                            i=k;
                        }
                    }
                }
                len_panel=len_panel_new;
            } while (len_panel>0);
        } else {
            // the unknown i is moved to C
            split_marker[i]=PASO_AMG_IN_C;
            lambda[i]=-1; // lambda from unavailable unknowns is set to -1

            // all undecided unknowns strongly coupled to i are moved to F
            for (p=0; p<degree_ST[i]; ++p) {
                j=ST[offset_S[i]+p];
                if (split_marker[j]==PASO_AMG_UNDECIDED) {
                    split_marker[j]=PASO_AMG_IN_F;
                    lambda[j]=-1;

                    for (q=0; q<degree_S[j]; ++q) {
                        k=S[offset_S[j]+q];
                        if (split_marker[k]==PASO_AMG_UNDECIDED)
                            lambda[k]++;
                    }
                }
            }
            for (p=0; p<degree_S[i]; ++p) {
                j=S[offset_S[i]+p];
                if (split_marker[j]==PASO_AMG_UNDECIDED)
                    lambda[j]--;
            }
        }
        i=util::arg_max(n,lambda);
    }
    delete[] lambda;
    delete[] ST;
    delete[] degree_ST;
    if (usePanel) {
        delete[] panel;
        if (!SMALL_PANEL)
            delete[] notInPanel;
    }
}

/// ensures that two F nodes are connected via a C node
void Preconditioner_LocalAMG_enforceFFConnectivity(dim_t n,
            const index_t* offset_S, const dim_t* degree_S, const index_t* S,
            AMGBlockSelect* split_marker)
{
    dim_t i, p, q;
    // now we make sure that two (strongly) connected F nodes are (strongly)
    // connected via a C node.
    for (i=0; i<n; ++i) {
        if (split_marker[i]==PASO_AMG_IN_F && degree_S[i]>0) {
            for (p=0; p<degree_S[i]; ++p) {
                index_t j=S[offset_S[i]+p];
                if ( (split_marker[j]==PASO_AMG_IN_F)  && (degree_S[j]>0) )  {
                    // i and j are now two F nodes which are strongly connected
                    // is there a C node they share ?
                    index_t sharing=-1;
                    for (q=0; q<degree_S[i]; ++q) {
                        index_t k=S[offset_S[i]+q];
                        if (split_marker[k]==PASO_AMG_IN_C) {
                            index_t* where_k = (index_t*)bsearch(
                                    &k, &(S[offset_S[j]]), degree_S[j],
                                    sizeof(index_t), util::comparIndex);
                            if (where_k != NULL) {
                                sharing=k;
                                break;
                            }
                        }
                    }
                    if (sharing < 0) {
                        if (i < j) {
                            split_marker[j] = PASO_AMG_IN_C;
                        } else {
                            split_marker[i] = PASO_AMG_IN_C;
                            // no point to look any further as i is now a C node
                            break;
                        }
                    }
                }
            }
        }
    }
}

} // namespace paso

