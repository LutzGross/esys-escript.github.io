
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

#include "Solver.h"
#include "PasoUtil.h"

#include <iostream>

namespace paso {

SolverResult Solver_GMRES2(Function* F, const double* f0, const double* x0,
                           double* dx, dim_t* iter, double* tolerance,
                           Performance* pp)
{
    static double RENORMALIZATION_CONST=0.001;
    const dim_t l=(*iter)+1, iter_max=*iter;
    dim_t k=0, i, j;
    const dim_t n = F->getLen();
    const double rel_tol = *tolerance;
    double abs_tol, normf0, normv, normv2, hh, hr, nu, norm_of_residual = 0.;
    bool breakFlag = false, maxIterFlag = false, convergeFlag = false;

    if (n < 0 || iter_max<=0 || l<1 || rel_tol<0) {
        return InputError;
    }

    SolverResult status=NoError;

    double* h = new double[l*l];
    double** v = new double*[l];
    double* c = new double[l];
    double* s = new double[l];
    double* g = new double[l];
    double* work = new double[n];

    for (i=0; i<iter_max; i++)
        v[i]=NULL;

    util::zeroes(n,dx);

    /*
     *  the show begins:
     */
    normf0 = util::l2(n, f0, F->mpi_info);
    k = 0;
    convergeFlag = (std::abs(normf0)<=0);
    if (!convergeFlag) {
        abs_tol = rel_tol*normf0;
        std::cout << "GMRES2 initial residual norm " << normf0
            << " (rel. tol = " << rel_tol << ")" << std::endl;
        v[0] = new double[n];
        util::zeroes(n, v[0]);
        util::update(n, 1., v[0], -1./normf0, f0); // v = -1./normf0*f0
        g[0] = normf0;
        while (!breakFlag && !maxIterFlag && !convergeFlag && status==NoError) {
            k++;
            v[k]=new double[n];
            /*
             * call directional derivative function
             */
            F->derivative(v[k], v[k-1], f0, x0, work, pp);
            normv=util::l2(n,v[k],F->mpi_info);
            /*
             * Modified Gram-Schmidt
             */
            for (j=0; j<k; j++){
                hh = util::innerProduct(n,v[j],v[k],F->mpi_info);
                util::update(n,1.,v[k],(-hh),v[j]); // v[k]-hh*v[j]
                h[INDEX2(j,k-1,l)]=hh;
                //printf("%d :  %d = %e\n",k,j,hh);
            }
            normv2=util::l2(n,v[k],F->mpi_info);
            h[INDEX2(k,k-1,l)]=normv2;
            /*
             * reorthogonalize
             */
            if (!(normv + RENORMALIZATION_CONST*normv2 > normv)) {
                // printf("GMRES2: renormalization!");
                for (j=0; j<k; j++) {
                    hr = util::innerProduct(n,v[j],v[k],F->mpi_info);

                    h[INDEX2(j,k-1,l)]+=hr;
                    util::update(n,1.,v[k],(-hr),v[j]);
                }
                normv2=util::l2(n,v[k],F->mpi_info);
                h[INDEX2(k,k-1,l)]=normv2;
            }
            /*
             * watch out for happy breakdown
             */
            if (normv2 > 0.) {
                util::update(n,1./normv2,v[k],0.,v[k]); /* normalize v[k] */
            }
            /*
             * Form and store the information for the new Givens rotation
             */
            util::applyGivensRotations(k,&h[INDEX2(0,k-1,l)],c,s);

            /*
             * Don't divide by zero if solution has been found
             */
            g[k] = 0;
            nu = sqrt(h[INDEX2(k-1,k-1,l)]*h[INDEX2(k-1,k-1,l)]+h[INDEX2(k,k-1,l)]*h[INDEX2(k,k-1,l)]);
            if (nu > 0) {
                c[k-1]= h[INDEX2(k-1,k-1,l)]/nu;
                s[k-1]=-h[INDEX2(k,k-1,l)]/nu;
                h[INDEX2(k-1,k-1,l)]=c[k-1]*h[INDEX2(k-1,k-1,l)]-s[k-1]*h[INDEX2(k,k-1,l)];
                h[INDEX2(k,k-1,l)]=0;
                util::applyGivensRotations(2,&(g[k-1]),&(c[k-1]),&(s[k-1]));
            }
            norm_of_residual = fabs(g[k]);
            maxIterFlag = (k >= iter_max);
            convergeFlag = (norm_of_residual <= abs_tol);
            std::cout << "GMRES2 step " << k << ": residual " << fabs(g[k])
                << " (abs. tol = " << abs_tol << ")" << std::endl;
        }
    }

    // all done and ready for the forward substitution:
    for (i=k-1;i>=0;--i) {
        for (j=i+1;j<k;j++) {
            g[i]-=h[INDEX2(i,j,l)]*g[j];
        }
        g[i] /= h[INDEX2(i,i,l)];
        util::update(n, 1., dx, g[i], v[i]); // dx = dx+g[i]*v[i]
    }
    if ( v != NULL) {
        for (i=0; i<iter_max; i++)
            delete[] v[i];
    }
    delete[] h;
    delete[] v;
    delete[] c;
    delete[] s;
    delete[] g;
    delete[] work;
    *iter=k;
    *tolerance=norm_of_residual;
    return status;
}

} // namespace paso

