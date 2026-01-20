
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

// Some utility routines

#include "PasoUtil.h"

namespace paso {

namespace util {

int comparIndex(const void* index1, const void* index2)
{
    return *(index_t*)index1 - *(index_t*)index2;
}

bool isAny(dim_t N, const index_t* array, index_t value)
{
    bool out = false;
#pragma omp parallel for reduction(||:out)
    for (index_t i = 0; i < N; i++)
        out = out || (array[i]==value);
    return out;
}

dim_t numPositives(dim_t N, const double* x, escript::JMPI mpiInfo)
{
    dim_t myOut = 0;
    dim_t out = 0;

#pragma omp parallel
    {
        dim_t localOut = 0;
#pragma omp for schedule(static)
        for (index_t i = 0; i < N; i++) {
            if (x[i] > 0)
                localOut++;
        }
#pragma omp critical
        {
            myOut = myOut + localOut;
        }
    }
#ifdef ESYS_MPI
    MPI_Allreduce(&myOut, &out, 1, MPI_DIM_T, MPI_SUM, mpiInfo->comm);
#else
    out = myOut;
#endif
    return out;
}

index_t iMax(dim_t N, const index_t* array)
{
    const index_t INDEX_T_MIN = escript::DataTypes::index_t_min();
    index_t out = INDEX_T_MIN;

    if (N > 0) {
#pragma omp parallel
        {
            index_t out2 = INDEX_T_MIN;
#pragma omp for schedule(static)
            for (index_t i = 0; i < N; i++)
                out2 = std::max(out2, array[i]);
#pragma omp critical
            {
                out = std::max(out, out2);
            }
        }
    }
    return out;
}

index_t cumsum(dim_t N, index_t* array)
{
    index_t out = 0, tmp;
    dim_t i;
#ifdef _OPENMP
    const int num_threads = omp_get_max_threads();
#else
    const int num_threads = 1;
#endif

    if (num_threads > 1) {
#ifdef _OPENMP
        index_t* partial_sums = new index_t[num_threads];
#pragma omp parallel private(i,tmp)
        {
            index_t sum = 0;
            const int thread_num = omp_get_thread_num();
#pragma omp for schedule(static)
            for (i = 0; i < N; ++i)
                sum += array[i];

            partial_sums[thread_num]=sum;
#pragma omp barrier
#pragma omp master
            {
                out=0;
                for (i = 0; i < num_threads; ++i) {
                    tmp = out;
                    out += partial_sums[i];
                    partial_sums[i] = tmp;
                }
            }
#pragma omp barrier
            sum = partial_sums[thread_num];
#pragma omp for schedule(static)
            for (i = 0; i < N; ++i) {
                tmp = sum;
                sum += array[i];
                array[i] = tmp;
            }
        }
        delete[] partial_sums;
#endif // _OPENMP
    } else {
        for (i = 0; i < N; ++i) {
            tmp = out;
            out += array[i];
            array[i] = tmp;
        }
    }
    return out;
}

index_t cumsum_maskedTrue(dim_t N, index_t* array, int* mask)
{
    index_t out=0;
    dim_t i;
#ifdef _OPENMP
    const int num_threads=omp_get_max_threads();
#else
    const int num_threads=1;
#endif

    if (num_threads>1) {
#ifdef _OPENMP
        index_t* partial_sums=new index_t[num_threads];
        #pragma omp parallel private(i)
        {
            index_t sum=0, tmp;
            const int thread_num = omp_get_thread_num();
            #pragma omp for schedule(static)
            for (i=0;i<N;++i) {
                if (mask[i]) {
                    array[i] =1;
                    sum++;
                } else {
                    array[i] =0;
                }
            }
            partial_sums[thread_num]=sum;
            #pragma omp barrier
            #pragma omp master
            {
                out=0;
                for (i=0;i<num_threads;++i) {
                    tmp=out;
                    out+=partial_sums[i];
                    partial_sums[i]=tmp;
                }
            }
            #pragma omp barrier
            sum=partial_sums[thread_num];
            #pragma omp for schedule(static)
            for (i=0;i<N;++i) {
                if (mask[i]) {
                    tmp=sum;
                    sum+=array[i];
                    array[i]=tmp;
                } else {
                    array[i]=-1;
                }
            }
        }
        delete[] partial_sums;
#endif
    } else { /* num_threads=1 */
        for (i=0;i<N;++i) {
            if (mask[i]) {
                array[i]=out;
                out++;
            } else {
                array[i]=-1;
            }
        }
    }
    return out;
}

index_t cumsum_maskedFalse(dim_t N, index_t* array, int* mask)
{
    index_t out=0;
    dim_t i;
#ifdef _OPENMP
    const int num_threads=omp_get_max_threads();
#else
    const int num_threads=1;
#endif

    if (num_threads>1) {
#ifdef _OPENMP
        index_t* partial_sums=new index_t[num_threads];
        #pragma omp parallel private(i)
        {
            index_t sum = 0, tmp = 0;
            const int thread_num=omp_get_thread_num();
            #pragma omp for schedule(static)
            for (i=0;i<N;++i) {
                if (! mask[i]) {
                    array[i] =1;
                    sum++;
                } else {
                    array[i] =0;
                }
            }
            partial_sums[thread_num]=sum;
            #pragma omp barrier
            #pragma omp master
            {
                out=0;
                for (i=0;i<num_threads;++i) {
                    tmp=out;
                    out+=partial_sums[i];
                    partial_sums[i]=tmp;
                }
            }
            #pragma omp barrier
            sum=partial_sums[thread_num];
            #pragma omp for schedule(static)
            for (i=0;i<N;++i) {
                if (! mask[i]) {
                    tmp=sum;
                    sum+=array[i];
                    array[i]=tmp;
                } else {
                    array[i]=-1;
                }
            }
        }
        delete[] partial_sums;
#endif
    } else {
        for (i=0;i<N;++i) {
            if (! mask[i]) {
                array[i]=out;
                out++;
            } else {
                array[i]=-1;
            }
        }
    }
    return out;
}

index_t arg_max(dim_t n, dim_t* lambda)
{
    dim_t i;
    index_t max=-1;
    index_t argmax=-1;
    index_t lmax=-1;
    index_t li=-1;
#ifdef _OPENMP
    const int num_threads=omp_get_max_threads();
#else
    const int num_threads=1;
#endif

    if (n>0) {
        max=lambda[0];
        argmax=0;
        if (num_threads > 1) {
#pragma omp parallel private(i,lmax,li)
            {
                lmax=max;
                li=argmax;
                #pragma omp for schedule(static)
                for (i=0;i<n;++i) {
                    if(lmax<lambda[i]){
                        lmax=lambda[i];
                        li=i;
                    }
                }
                #pragma omp critical
                {
                    if (max < lmax)  {
                        max=lmax;
                        argmax=li;
                    } else if (max==lmax && argmax>li) {
                        argmax=li;
                    }
                }
            }
        } else {
            for (i=0;i<n;++i) {
                if(max<lambda[i]){
                    max=lambda[i];
                    argmax=i;
                }
            }
        }
    }
    return argmax;
}

void zeroes(dim_t n, double* x)
{
    dim_t i,local_n,rest,n_start,n_end,q;
#ifdef _OPENMP
    const int num_threads=omp_get_max_threads();
#else
    const int num_threads=1;
#endif

#pragma omp parallel for private(i,local_n,rest,n_start,n_end,q)
    for (i=0;i<num_threads;++i) {
        local_n=n/num_threads;
        rest=n-local_n*num_threads;
        n_start=local_n*i+std::min(i,rest);
        n_end=local_n*(i+1)+std::min(i+1,rest);
        #pragma ivdep
        for (q=n_start;q<n_end;++q) x[q]=0;
    }
}

void update(dim_t n, double a, double* x, double b, const double* y)
{
    dim_t i,local_n,rest,n_start,n_end,q;
#ifdef _OPENMP
    const int num_threads=omp_get_max_threads();
#else
    const int num_threads=1;
#endif

    #pragma omp parallel for private(i,local_n,rest,n_start,n_end,q)
    for (i=0;i<num_threads;++i) {
        local_n=n/num_threads;
        rest=n-local_n*num_threads;
        n_start=local_n*i+std::min(i,rest);
        n_end=local_n*(i+1)+std::min(i+1,rest);
        if (a==1.) {
            #pragma ivdep
            for (q=n_start;q<n_end;++q) {
                x[q]+=b*y[q];
            }
        } else if ((a==1.) && (b==1)) {
            #pragma ivdep
            for (q=n_start;q<n_end;++q) {
                x[q]+=y[q];
            }
        } else if ((std::abs(a)==0) && (std::abs(b)==0)) {
            #pragma ivdep
            for (q=n_start;q<n_end;++q) {
                x[q]=0.;
            }
        } else if ((std::abs(a)==0) && (b==1)) {
            #pragma ivdep
            for (q=n_start;q<n_end;++q) {
                x[q]=y[q];
            }
        } else if (std::abs(a)==0) {
            #pragma ivdep
            for (q=n_start;q<n_end;++q) {
                x[q]=b*y[q];
            }
        } else {
            #pragma ivdep
            for (q=n_start;q<n_end;++q) {
                x[q]=a*x[q]+b*y[q];
            }
        }
    }
}

void linearCombination(dim_t n, double* z, double a, const double* x,
                       double b, const double* y)
{
    dim_t i,local_n,rest,n_start,n_end,q;
#ifdef _OPENMP
    const int num_threads=omp_get_max_threads();
#else
    const int num_threads=1;
#endif

#pragma omp parallel for private(i,local_n,rest,n_start,n_end,q)
    for (i=0;i<num_threads;++i) {
        local_n=n/num_threads;
        rest=n-local_n*num_threads;
        n_start=local_n*i+std::min(i,rest);
        n_end=local_n*(i+1)+std::min(i+1,rest);
        if (((std::abs(a)==0) && (std::abs(b)==0)) || (y==NULL) || (x==NULL)) {
            #pragma ivdep
            for (q=n_start;q<n_end;++q) {
                z[q]=0.;
            }
        } else if ( ((std::abs(a)==0) && (std::abs(b)>0.)) || (x==NULL) )  {
            #pragma ivdep
            for (q=n_start;q<n_end;++q) {
                z[q]=b*y[q];
            }
        } else if (((std::abs(a)>0) && (std::abs(b)==0.)) || (y==NULL) ) {
            #pragma ivdep
            for (q=n_start;q<n_end;++q) {
                z[q]=a*x[q];
            }
        } else {
            #pragma ivdep
            for (q=n_start;q<n_end;++q) {
                z[q]=a*x[q]+b*y[q];
            }
        }
    }
}

double innerProduct(dim_t n, const double* x, const double* y,
                    escript::JMPI mpiinfo)
{
    dim_t i,local_n,rest,n_start,n_end,q;
    double my_out=0, local_out=0., out=0.;
#ifdef _OPENMP
    const int num_threads=omp_get_max_threads();
#else
    const int num_threads=1;
#endif

#pragma omp parallel for private(i,local_out,local_n,rest,n_start,n_end,q)
    for (i=0; i<num_threads; ++i) {
        local_out = 0;
        local_n = n/num_threads;
        rest = n-local_n*num_threads;
        n_start = local_n*i+std::min(i,rest);
        n_end = local_n*(i+1)+std::min(i+1,rest);
        #pragma ivdep
        for (q=n_start; q<n_end; ++q)
            local_out += x[q]*y[q];
#pragma omp critical
        {
            my_out+=local_out;
        }
    }
#ifdef ESYS_MPI
#pragma omp single
    {
        MPI_Allreduce(&my_out, &out, 1, MPI_DOUBLE, MPI_SUM, mpiinfo->comm);
    }
#else
       out=my_out;
#endif

    return out;
}

double lsup(dim_t n, const double* x, escript::JMPI mpiinfo)
{
    dim_t i,local_n,rest,n_start,n_end,q;
    double my_out=0., local_out=0., out=0.;
#ifdef _OPENMP
    const int num_threads=omp_get_max_threads();
#else
    const int num_threads=1;
#endif

#pragma omp parallel for private(i,local_n,rest,n_start,n_end,q, local_out)
    for (i=0; i<num_threads; ++i) {
        local_n=n/num_threads;
        rest=n-local_n*num_threads;
        n_start=local_n*i+std::min(i,rest);
        n_end=local_n*(i+1)+std::min(i+1,rest);
        local_out=0;
        for (q=n_start; q<n_end; ++q)
            local_out = std::max(std::abs(x[q]), local_out);
#pragma omp critical
        {
            my_out = std::max(my_out,local_out);
        }
    }
#ifdef ESYS_MPI
#pragma omp single
    {
        MPI_Allreduce(&my_out, &out, 1, MPI_DOUBLE, MPI_MAX, mpiinfo->comm);
    }
#else
    out = my_out;
#endif
    return out;
}

double l2(dim_t n, const double* x, escript::JMPI mpiinfo)
{
    double my_out = 0, out = 0.;
#ifdef _OPENMP
    const int num_threads=omp_get_max_threads();
#else
    const int num_threads=1;
#endif

#pragma omp parallel for
    for (dim_t i=0; i<num_threads; ++i) {
        const dim_t local_n = n/num_threads;
        const dim_t rest = n-local_n*num_threads;
        const dim_t n_start = local_n*i+std::min(i,rest);
        const dim_t n_end = local_n*(i+1)+std::min(i+1,rest);
        double local_out = 0.;
        for (dim_t q=n_start; q<n_end; ++q)
            local_out += x[q]*x[q];
#pragma omp critical
        {
            my_out += local_out;
        }
    }
#ifdef ESYS_MPI
    #pragma omp single
    {
        MPI_Allreduce(&my_out, &out, 1, MPI_DOUBLE, MPI_SUM, mpiinfo->comm);
    }
#else
    out = my_out;
#endif

    return sqrt(out);
}

void applyGivensRotations(dim_t n, double* v, const double* c, const double* s)
{
    #pragma ivdep
    for (dim_t i=0; i<n-1;++i) {
        const double w1 = c[i]*v[i]-s[i]*v[i+1];
        const double w2 = s[i]*v[i]+c[i]*v[i+1];
        v[i] = w1;
        v[i+1] = w2;
   }
}

} // namespace util
} // namespace paso

