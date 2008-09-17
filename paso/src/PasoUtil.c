/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/*   Some utility routines: */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004,2005 */

/**************************************************************/

#include "Common.h"
#include "PasoUtil.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/* returns true if array contains value */
bool_t Paso_Util_isAny(dim_t N,index_t* array,index_t value) {
   bool_t out=FALSE;
   dim_t i;
   #pragma omp parallel for private(i) schedule(static) reduction(||:out)
   for (i=0;i<N;i++) {
            out = out || (array[i]==value);
   }
   return out;
}

/**************************************************************/

/* calculates the cummultative sum in array and returns the total sum */

/**************************************************************/
index_t Paso_Util_cumsum(dim_t N,index_t* array) {
   index_t out=0,tmp;
   dim_t i;
   #ifdef _OPENMP
      index_t *partial_sums=NULL,sum;
      partial_sums=TMPMEMALLOC(omp_get_max_threads(),index_t);
      #pragma omp parallel private(sum,i,tmp)
      {
        sum=0;
        #pragma omp for schedule(static)
        for (i=0;i<N;++i) sum+=array[i];

        partial_sums[omp_get_thread_num()]=sum;
        #pragma omp barrier
        #pragma omp single
        {
          out=0;
          for (i=0;i<omp_get_max_threads();++i) {
             tmp=out;
             out+=partial_sums[i];
             partial_sums[i]=tmp;
           } 
        }
        sum=partial_sums[omp_get_thread_num()];
        #pragma omp for schedule(static)
        for (i=0;i<N;++i) {
          tmp=sum;
          sum+=array[i];
          array[i]=tmp;
        } 
      }
      TMPMEMFREE(partial_sums);
   #else 
      for (i=0;i<N;++i) {
         tmp=out;
         out+=array[i];
         array[i]=tmp;
      }
   #endif
   return out;
}
/*
 * set x to zero:
 */
void Paso_zeroes(const dim_t n, double* x) 
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
        n_start=local_n*i+MIN(i,rest);
        n_end=local_n*(i+1)+MIN(i+1,rest);
        #pragma ivdep
        for (q=n_start;q<n_end;++q) x[q]=0;
   }

}
/*
 * performs an update of the form x=a*x+b*y  where y and x are long vectors.
 * if b=0, y is not used.
 */
void Paso_Update(const dim_t n, const double a, double* x, const double b, const double* y) 
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
        n_start=local_n*i+MIN(i,rest);
        n_end=local_n*(i+1)+MIN(i+1,rest);
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
        } else if ((ABS(a)==0) && (ABS(b)==0)) {
            #pragma ivdep
            for (q=n_start;q<n_end;++q) {
              x[q]=0.;
            }
        } else if ((ABS(a)==0) && (ABS(b)==1)) {
            #pragma ivdep
            for (q=n_start;q<n_end;++q) {
              x[q]=y[q];
            }
        } else if (ABS(a)==0) {
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
void Paso_Copy(const dim_t n, double* out, const double* in) {
    Paso_LinearCombination(n,out, 1.,in,0., in);
}
/*
 * performs an update of the form z=a*x+b*y  where y and x are long vectors.
 * if a=0, x is not used.
 * if b=0, y is not used.
 */
void Paso_LinearCombination(const dim_t n, double*z, const double a,const double* x, const double b, const double* y)
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
        n_start=local_n*i+MIN(i,rest);
        n_end=local_n*(i+1)+MIN(i+1,rest);
        if (((ABS(a)==0) && (ABS(b)==0)) || (y==NULL) || (x==NULL)) {
            #pragma ivdep
            for (q=n_start;q<n_end;++q) {
              z[q]=0.;
            }
        } else if ( ((ABS(a)==0) && (ABS(b)>0.)) || (x==NULL) )  {
            #pragma ivdep
            for (q=n_start;q<n_end;++q) {
              z[q]=b*y[q];
            }
        } else if (((ABS(a)>0) && (ABS(b)==0.)) || (y==NULL) ) {
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



double Paso_InnerProduct(const dim_t n,const double* x, const double* y, Paso_MPIInfo* mpiinfo)
{
   dim_t i,local_n,rest,n_start,n_end,q;
   double my_out=0, local_out=0., out=0.;
   #ifdef _OPENMP
       const int num_threads=omp_get_max_threads();
   #else
       const int num_threads=1;
   #endif
   #pragma omp parallel for private(i,local_out,local_n,rest,n_start,n_end,q) reduction(+:my_out)
   for (i=0;i<num_threads;++i) {
        local_out=0;
        local_n=n/num_threads;
        rest=n-local_n*num_threads;
        n_start=local_n*i+MIN(i,rest);
        n_end=local_n*(i+1)+MIN(i+1,rest);
        #pragma ivdep
        for (q=n_start;q<n_end;++q) local_out+=x[q]*y[q];
        my_out+=local_out;
   }
   #ifdef PASO_MPI
      #pragma omp single
      {
          MPI_Allreduce(&my_out,&out, 1, MPI_DOUBLE, MPI_SUM, mpiinfo->comm);
       }
   #else
       out=my_out;
   #endif

   return out;

}
double Paso_l2(const dim_t n, const double* x, Paso_MPIInfo* mpiinfo) 
{
   dim_t i,local_n,rest,n_start,n_end,q;
   double my_out=0, local_out=0., out=0.;
   #ifdef _OPENMP
       const int num_threads=omp_get_max_threads();
   #else
       const int num_threads=1;
   #endif
   #pragma omp parallel for private(i,local_n,rest,n_start,n_end,q, local_out) reduction(+:my_out)
   for (i=0;i<num_threads;++i) {
        local_n=n/num_threads;
        rest=n-local_n*num_threads;
        n_start=local_n*i+MIN(i,rest);
        n_end=local_n*(i+1)+MIN(i+1,rest);
        local_out=0;
        for (q=n_start;q<n_end;++q) local_out+=x[q]*x[q];
        my_out+=local_out;
   }
   #ifdef PASO_MPI
      #pragma omp single
      {
          MPI_Allreduce(&my_out,&out, 1, MPI_DOUBLE, MPI_SUM, mpiinfo->comm);
       }
   #else
       out=my_out;
   #endif

   return sqrt(out);

}
/*  
 *  Apply a sequence of n-1 Givens rotations (c,s) to v of length n which assumed to be small.
 *
 */
void ApplyGivensRotations(const dim_t n,double* v,const double* c,const double* s)
{
    dim_t i;
    register double w1,w2;
    #pragma ivdep
    for (i=0; i<n-1;++i) {
        w1=c[i]*v[i]-s[i]*v[i+1];
        w2=s[i]*v[i]+c[i]*v[i+1];
        v[i]=w1;
        v[i+1]=w2;
   }
}

bool_t Paso_fileExists( const char* filename )
{
    FILE* fp = NULL;
    fp = fopen(filename,"r");
    if( fp != NULL ) {
        fclose(fp);
        return TRUE;
    } else {
        return FALSE;
    }
}
