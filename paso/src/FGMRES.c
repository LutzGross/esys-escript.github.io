/* $Id:$ */

/*******************************************************
 *
 *       Copyright 2008 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/
/*
*  Purpose
*  =======
*
*  FGMRES solves the non-linear system f0+J_0*d=0 
*  where f0=F(x0), J_0 is the jacobian of F at x0.
*
*  Convergence test: norm(f0+J_0*d)<=tolerance*norm(f0)
*
*  Arguments
*  =========
*
*   Paso_Function * F   evaluation of F (including any preconditioner)
*
*   x0                  (input) current point
*
*   f0                  (input) function F at current point x0
*
*   d                   (output) solution of f0+J0*d=0 with accuracy tolerance
*
*   iter                (input/output)
*                       On input, the maximum num_iterations to be performed.
*                       On output, actual number of num_iterations performed.
*   tolerance           (input/output) 
*                        On input, the allowable convergence measure for norm(f0+J0*d)/norm(f0)
*                        On output, the final value for norm(f0+J0*d)
*  return value
*
*         =SOLVER_NO_ERROR: Successful exit. approximate solution returned.
*         =SOLVER_MAXNUM_ITER_REACHED
*         =SOLVER_INPUT_ERROR Illegal parameter:
*         =SOLVER_BREAKDOWN: bad luck!
*         =SOLVER_MEMORY_ERROR : no memory available
*
*  ==============================================================
*/

#include "Common.h"
#include "Solver.h"
#ifdef _OPENMP
#include <omp.h>
#endif

err_t Paso_Solver_NLGMRES(
    Paso_Function * F,
    const double* f0,
    const double* x0,
    double * x,
    dim_t *iter,
    double* tolerance,
    Paso_Performance* pp) 
{
  double static RENORMALIZATION_CONST=0.001;
  dim_t l=(*iter)+1, iter_max=*iter,k=0,n=F->local_n,i,j;
  double rel_tol=*tolerance, abs_tol, normf0, normv, normv2, hh, hr, nu, norm_of_residual=0.;
  bool_t breakFlag=FALSE, maxIterFlag=FALSE, convergeFlag=FALSE;
  double *h=NULL, **v=NULL, *c=NULL,*s=NULL,*g=NULL, *work=NULL;
  err_t Status=SOLVER_NO_ERROR;

  /*     Test the input parameters. */
  if (n < 0 || iter_max<=0 || l<1 || rel_tol<0) {
    return SOLVER_INPUT_ERROR;
  }
  Paso_zeroes(n,x);
  /*     
   *  allocate memory: 
   */
  h=TMPMEMALLOC(l*l,double);
  v=TMPMEMALLOC(iter_max,double*);
  c=TMPMEMALLOC(l,double);
  s=TMPMEMALLOC(l,double);
  g=TMPMEMALLOC(l,double);
  work=TMPMEMALLOC(n,double);

  if (h==NULL || v ==NULL || c== NULL || s == NULL || g== NULL || work==NULL) {
     Status=SOLVER_MEMORY_ERROR;
  } else {
     for (i=0;i<iter_max;i++) v[i]=NULL;
     for (i=0;i<iter_max;i++) {
       v[i]=TMPMEMALLOC(n,double);
       if (v[i]==NULL) {
          Status=SOLVER_MEMORY_ERROR;
          break;
       }
     }
  }
  if ( Status ==SOLVER_NO_ERROR ) {
     /*     
      *  the show begins:
      */
      normf0=Paso_l2(n,f0,F->mpi_info);
      k=0;
      convergeFlag=(ABS(normf0)<=0);
      if (! convergeFlag) {
          abs_tol=rel_tol*normf0;
          Paso_zeroes(n,v[0]);
          Paso_Update(n,1.,v[0],-1./normf0,f0);
          g[0]=normf0;
          while(! (breakFlag || maxIterFlag || convergeFlag)) {
               k++;
               /*
               *      call directional derivative function
               */
               Paso_FunctionDerivative(v[k],v[k-1],F,f0,x0,work);
               normv=Paso_l2(n,v[k],F->mpi_info);
               /*
                * Modified Gram-Schmidt
                */
                for (j=0;j<k;j++){
                   hh=Paso_InnerProduct(n,v[j],v[k],F->mpi_info);
                   Paso_Update(n,1.,v[k],(-hh),v[j]);
                   h[INDEX2(j,k-1,l)]=hh;
                }
                normv2=Paso_l2(n,v[k],F->mpi_info);
                h[INDEX2(k,k-1,l)]=normv2;
                /*
                 * reorthogonalize
                 */
                if (! (normv + RENORMALIZATION_CONST*normv2 > normv)) {
                     for (j=0;j<k;j++){
                         hr=Paso_InnerProduct(n,v[j],v[k],F->mpi_info);
	                 h[INDEX2(j,k-1,l)]+=hr;
                         Paso_Update(n,1.,v[k],(-hr),v[j]);
                     }
                     normv2=Paso_l2(n,v[k],F->mpi_info);
                     h[INDEX2(k,k-1,l)]=normv2;
                 }
                /* 
                 *   watch out for happy breakdown
                 */
                if(normv2 > 0.) {
                   Paso_update(n,1./normv2,v[k],0.,v[k]);
                } 
                /*
                 *   Form and store the information for the new Givens rotation
                 */
                ApplyGivensRotations(iter,&h[INDEX2(0,k-1,l)],c,s);
                /*
                 *  Don't divide by zero if solution has  been found
                 */
                g[k]=0;
                nu=sqrt(h[INDEX2(k-1,k-1,l)]*h[INDEX2(k-1,k-1,l)]+h[INDEX2(k,k-1,l)]*h[INDEX2(k,k-1,l)]);
                if (nu>0) {
                    c[k-1]=h[INDEX2(k-1,k-1,l)]/nu;
                    s[k-1]=-h[k,k-1]/nu;
                    h[INDEX2(k-1,k-1,l)]=c[k-1]*h[INDEX2(k-1,k-1,l)]-s[k-1]*h[k,k-1];
                    h[INDEX2(k,k-1,l)]=0;
                    ApplyGivensRotations(2,&(g[k-1]),&(c[k-1]),&(s[k-1]));
                }
                norm_of_residual=abs(g[k]);
                maxIterFlag = (k>=iter_max);
                convergeFlag = (abs(g[k]) <= abs_tol);
printf("FGMRES step %d: error %e (tol=%d)\n",k,abs(g[k]),abs_tol);
          }
      }
      /*
       * all done and ready for the forward substitution:
      */
      
      for (i=0;i<k;i++) {
           for (j=0;j<i;j++) {
              g[i]-=h[INDEX2(j,i,l)]*g[j];
           }
           g[i]/=h[INDEX2(i,i,l)];
           Paso_Update(n,x[k],g[i],v[i]);
      }
 }
 /*     
  *  clean up:
  */
  if ( v !=NULL) {
    for (i=0;i<iter_max;i++) TMPMEMFREE(v);
  }
  TMPMEMFREE(h);
  TMPMEMFREE(v);
  TMPMEMFREE(c);
  TMPMEMFREE(s);
  TMPMEMFREE(g);
  TMPMEMFREE(work);
  *iter=k;
  *tolerance=norm_of_residual;
  return Status;
}

