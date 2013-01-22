
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


#include "Common.h"
#include "Solver.h"
#include "PasoUtil.h"
#ifdef _OPENMP
#include <omp.h>
#endif

err_t Paso_Solver_GMRES2(
    Paso_Function * F,
    const double* f0,
    const double* x0,
    double * dx,
    dim_t *iter,
    double* tolerance,
    Paso_Performance* pp) 
{
  static double RENORMALIZATION_CONST=0.001;
  const dim_t l=(*iter)+1, iter_max=*iter;
  dim_t k=0,i,j;
  const dim_t n=F->n;
  const double rel_tol=*tolerance;
  double abs_tol, normf0, normv, normv2, hh, hr, nu, norm_of_residual=0.;
  bool_t breakFlag=FALSE, maxIterFlag=FALSE, convergeFlag=FALSE;
  double *h=NULL, **v=NULL, *c=NULL,*s=NULL,*g=NULL, *work=NULL;
  err_t Status=SOLVER_NO_ERROR;

  /*     Test the input parameters. */
  if (n < 0 || iter_max<=0 || l<1 || rel_tol<0) {
    return SOLVER_INPUT_ERROR;
  }
  Paso_zeroes(n,dx);
  /*     
   *  allocate memory: 
   */
  h=TMPMEMALLOC(l*l,double);
  v=TMPMEMALLOC(l,double*);
  c=TMPMEMALLOC(l,double);
  s=TMPMEMALLOC(l,double);
  g=TMPMEMALLOC(l,double);
  work=TMPMEMALLOC(n,double);

  if (h==NULL || v ==NULL || c== NULL || s == NULL || g== NULL || work==NULL) {
     Status=SOLVER_MEMORY_ERROR;
  }
  if ( Status ==SOLVER_NO_ERROR ) {
     for (i=0;i<iter_max;i++) v[i]=NULL;
     /*     
      *  the show begins:
      */
      normf0=Paso_l2(n,f0,F->mpi_info);
      k=0;
      convergeFlag=(ABS(normf0)<=0);
      if (! convergeFlag) {
          abs_tol=rel_tol*normf0;
          printf("GMRES2 initial residual norm %e (rel. tol=%e)\n",normf0,rel_tol);
          v[0]=TMPMEMALLOC(n,double);
          if (v[0]==NULL) {
             Status=SOLVER_MEMORY_ERROR;
          } else {
             Paso_zeroes(n,v[0]);
             Paso_Update(n,1.,v[0],-1./normf0,f0); /* v = -1./normf0*f0 */
             g[0]=normf0;
          }
          while( (! (breakFlag || maxIterFlag || convergeFlag)) && (Status ==SOLVER_NO_ERROR) ) {
               k++;
               v[k]=TMPMEMALLOC(n,double);
               if (v[k]==NULL) {
                  Status=SOLVER_MEMORY_ERROR;
               } else {
                  /*
                  *      call directional derivative function
                  */
                  Paso_FunctionDerivative(v[k],v[k-1],F,f0,x0,work,pp);
                  normv=Paso_l2(n,v[k],F->mpi_info);
                  /*
                   * Modified Gram-Schmidt
                   */
                   for (j=0;j<k;j++){
                      hh=Paso_InnerProduct(n,v[j],v[k],F->mpi_info);
                      Paso_Update(n,1.,v[k],(-hh),v[j]); /* v[k]-hh*v[j] */
                      h[INDEX2(j,k-1,l)]=hh; 
 /* printf("%d :  %d = %e\n",k,j,hh); */ 
		   }
                   normv2=Paso_l2(n,v[k],F->mpi_info);
                   h[INDEX2(k,k-1,l)]=normv2;
                   /*
                    * reorthogonalize
                    */
                   if (! (normv + RENORMALIZATION_CONST*normv2 > normv)) {
/* printf("GMRES2: renormalization!"); */
                        for (j=0;j<k;j++){
                            hr=Paso_InnerProduct(n,v[j],v[k],F->mpi_info);

	                    h[INDEX2(j,k-1,l)]+=hr;
                            Paso_Update(n,1.,v[k],(-hr),v[j]);
                        }
                        normv2=Paso_l2(n,v[k],F->mpi_info);
                        h[INDEX2(k,k-1,l)]=normv2;
                    }
/*
{int p,q;
   for (p=0;p<k+1;p++){
   for (q=0;q<k+1;q++)printf("%e ",h[INDEX2(p,q,l)]);
    printf("\n");
   
   }
}
*/
                   /* 
                    *   watch out for happy breakdown
                    */
                   if(normv2 > 0.) {
                      Paso_Update(n,1./normv2,v[k],0.,v[k]); /* normalize v[k] */
                   } 
                   /*
                    *   Form and store the information for the new Givens rotation
                    */
                   ApplyGivensRotations(k,&h[INDEX2(0,k-1,l)],c,s);

                   /*
                    *  Don't divide by zero if solution has  been found
                    */
                   g[k]=0;
		   nu=sqrt(h[INDEX2(k-1,k-1,l)]*h[INDEX2(k-1,k-1,l)]+h[INDEX2(k,k-1,l)]*h[INDEX2(k,k-1,l)]);
                   if (nu>0) {
                       c[k-1]= h[INDEX2(k-1,k-1,l)]/nu;
                       s[k-1]=-h[INDEX2(k,k-1,l)]/nu;
                       h[INDEX2(k-1,k-1,l)]=c[k-1]*h[INDEX2(k-1,k-1,l)]-s[k-1]*h[INDEX2(k,k-1,l)];
                       h[INDEX2(k,k-1,l)]=0;
                       ApplyGivensRotations(2,&(g[k-1]),&(c[k-1]),&(s[k-1]));
                   }
                   norm_of_residual=fabs(g[k]);
                   maxIterFlag = (k>=iter_max);
                   convergeFlag = (norm_of_residual <= abs_tol);
                   printf("GMRES2 step %d: residual %e (abs. tol=%e)\n",k,fabs(g[k]),abs_tol);
              }
          }
      }
      /*
       * all done and ready for the forward substitution:
      */

      for (i=k-1;i>=0;--i) {
           for (j=i+1;j<k;j++) {
              g[i]-=h[INDEX2(i,j,l)]*g[j];
           }
           g[i]/=h[INDEX2(i,i,l)];
           Paso_Update(n,1.,dx,g[i],v[i]); /* dx = dx+g[i]*v[i] */ 
      }
 }
 /*     
  *  clean up:
  */
  if ( v !=NULL) {
    for (i=0;i<iter_max;i++) TMPMEMFREE(v[i]);
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

