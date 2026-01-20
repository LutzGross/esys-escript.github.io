
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/

#include "Functions.h"
#include "PasoUtil.h"
#include "Solver.h"

namespace paso {

Function::Function(const escript::JMPI& mpiInfo) :
    mpi_info(mpiInfo)
{
}

Function::~Function()
{
}

SolverResult Function::derivative(double* J0w, const double* w, const double* f0,
                           const double* x0, double* setoff, Performance* pp)
{
    const real_t EPSILON = escript::DataTypes::real_t_eps();
    SolverResult err = NoError;
    dim_t i;
    double aw;
    const double epsnew = sqrt(EPSILON);
    double ttt, s=epsnew, local_s, norm_w=0.;
    const dim_t n = getLen();

    //double norm_x0=util::lsup(n,x0,mpi_info);
    norm_w=util::lsup(n, w, mpi_info);
    ttt=sqrt(EPSILON)*norm_w;
#pragma omp parallel private(local_s)
    {
        local_s=s;
#pragma omp for private(i, aw)
        for (i=0;i<n;++i) {
            aw=fabs(w[i]);
            if ( aw>ttt ) {
                local_s=std::max(local_s,fabs(x0[i])/aw);
            }
        }
#pragma omp critical
        {
            s=std::max(s,local_s);
        }
    }
#ifdef ESYS_MPI
    {
        double local_v[2], v[2];
        local_v[0]=s;
        local_v[1]=norm_w;
        MPI_Allreduce(local_v,v, 2, MPI_DOUBLE, MPI_MAX, mpi_info->comm);
        s=v[0];
        norm_w=v[1];
    }
#endif
    //printf("s ::  = %e, %e \n",s, norm_x0/norm_w);
    if (norm_w>0) {
        s=s*epsnew;
        //printf("s = %e\n",s);
        util::linearCombination(n,setoff,1.,x0,s,w);
        err = call(J0w, setoff, pp);
        if (err==NoError) {
            util::update(n,1./s,J0w,-1./s,f0); // J0w = (J0w - f0)/epsnew;
            //for (int i=0;i<n; i++) printf("df[%d]=%e %e\n",i,J0w[i],w[i]);
        }
    } else {
        util::zeroes(n,J0w);
    }
    return err;
}

} // namespace paso

