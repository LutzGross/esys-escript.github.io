
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/*
*  Purpose
*  =======
*
*  GMRES solves the linear system A*x=b using the
*  truncated and restered GMRES method with preconditioning. 
*
*  Convergence test: norm( b - A*x )< TOL.
*
*  
*
*  Arguments
*  =========
*
*  r       (input/output) double array, dimension n.
*          On entry, residual of inital guess X
*
*  x       (input/output) double array, dimension n.
*          On input, the initial guess.
*
*  iter    (input/output) int
*          On input, the maximum num_iterations to be performed.
*          On output, actual number of num_iterations performed.
*
*  tolerance (input/output) DOUBLE PRECISION
*          On input, the allowable convergence measure for
*          norm( b - A*x )
*          On output, the final value of this measure.
*
*  Length_of_recursion (input) gives the number of residual to be kept in orthogonalization process
*
*  restart (input) If restart>0, iteration is resterted a after restart steps.
*
*  INFO    (output) int
*
*         =SOLVER_NO_ERROR: Successful exit. num_iterated approximate solution returned.
*         =SOLVER_MAXNUM_ITER_REACHED
*         =SOLVER_INPUT_ERROR Illegal parameter:
*         =SOLVER_BREAKDOWN: bad luck!
*         =SOLVER_MEMORY_ERROR : no memory available
*
*  ==============================================================
*/

#include "Common.h"
#include "SystemMatrix.h"
#include "Solver.h"
#ifdef _OPENMP
#include <omp.h>
#endif

err_t Paso_Solver_GMRES(
    Paso_SystemMatrix * A,
    double * r,
    double * x,
    dim_t *iter,
    double * tolerance,dim_t Length_of_recursion,dim_t restart,
    Paso_Performance* pp) {

   /* Local variables */

   #ifdef _OPENMP
       const int num_threads=omp_get_max_threads();
   #else
       const int num_threads=1;
   #endif
  double *AP,**X_PRES,**R_PRES,**P_PRES, *dots, *loc_dots;
  double *P_PRES_dot_AP,*R_PRES_dot_P_PRES,*BREAKF,*ALPHA;
  double R_PRES_dot_AP0,P_PRES_dot_AP0,P_PRES_dot_AP1,P_PRES_dot_AP2,P_PRES_dot_AP3,P_PRES_dot_AP4,P_PRES_dot_AP5,P_PRES_dot_AP6,R_PRES_dot_P,breakf0;
  double tol,Factor,sum_BREAKF,gamma,SC1,SC2,norm_of_residual=0,diff,L2_R,Norm_of_residual_global=0;
  double *save_XPRES, *save_P_PRES, *save_R_PRES,save_R_PRES_dot_P_PRES;
  dim_t maxit,Num_iter_global=0,num_iter_restart=0,num_iter;
  dim_t i,z,order,n, Length_of_mem, th, local_n , rest, n_start ,n_end;
  bool_t breakFlag=FALSE, maxIterFlag=FALSE, convergeFlag=FALSE,restartFlag=FALSE;
  err_t Status=SOLVER_NO_ERROR;

  
  /* adapt original routine parameters */
  n = Paso_SystemMatrix_getTotalNumRows(A);
  Length_of_mem=MAX(Length_of_recursion,0)+1;

  /*     Test the input parameters. */
  if (restart>0) restart=MAX(Length_of_recursion,restart);
  if (n < 0 || Length_of_recursion<=0) {
    return SOLVER_INPUT_ERROR;
  }

  /*     allocate memory: */
  X_PRES=TMPMEMALLOC(Length_of_mem,double*);
  R_PRES=TMPMEMALLOC(Length_of_mem,double*);
  P_PRES=TMPMEMALLOC(Length_of_mem,double*);
  loc_dots=TMPMEMALLOC(MAX(Length_of_mem+1,3),double);
  dots=TMPMEMALLOC(MAX(Length_of_mem+1,3),double);
  P_PRES_dot_AP=TMPMEMALLOC(Length_of_mem,double);
  R_PRES_dot_P_PRES=TMPMEMALLOC(Length_of_mem,double);
  BREAKF=TMPMEMALLOC(Length_of_mem,double);
  ALPHA=TMPMEMALLOC(Length_of_mem,double);
  AP=TMPMEMALLOC(n,double);
  if (AP==NULL || X_PRES ==NULL || R_PRES == NULL || P_PRES == NULL || 
         P_PRES_dot_AP == NULL || R_PRES_dot_P_PRES ==NULL || BREAKF == NULL || ALPHA == NULL || dots==NULL || loc_dots==NULL) {
     Status=SOLVER_MEMORY_ERROR;
  } else {
     for (i=0;i<Length_of_mem;i++) {
       X_PRES[i]=TMPMEMALLOC(n,double);
       R_PRES[i]=TMPMEMALLOC(n,double);
       P_PRES[i]=TMPMEMALLOC(n,double);
       if (X_PRES[i]==NULL || R_PRES[i]==NULL ||  P_PRES[i]==NULL) Status=SOLVER_MEMORY_ERROR;
     }
  }
  if ( Status ==SOLVER_NO_ERROR ) {

    /* now PRES starts : */
    maxit=*iter;
    tol=*tolerance;

     /* initialization */

      restartFlag=TRUE;
      num_iter=0;

      #pragma omp parallel for private(th,z,i,local_n, rest, n_start, n_end)
      for (th=0;th<num_threads;++th) {
        local_n=n/num_threads;
        rest=n-local_n*num_threads;
        n_start=local_n*th+MIN(th,rest);
        n_end=local_n*(th+1)+MIN(th+1,rest);
        memset(&AP[n_start],0,sizeof(double)*(n_end-n_start));
        for(i=0;i<Length_of_mem;++i) {
           memset(&P_PRES[i][n_start],0,sizeof(double)*(n_end-n_start));
           memset(&R_PRES[i][n_start],0,sizeof(double)*(n_end-n_start));
           memset(&X_PRES[i][n_start],0,sizeof(double)*(n_end-n_start));
        }
      }

      while (! (convergeFlag || maxIterFlag || breakFlag))  {
         if (restartFlag) {
             BREAKF[0]=PASO_ONE;
             #pragma omp parallel for private(th,z, local_n, rest, n_start, n_end)
             for (th=0;th<num_threads;++th) {
               local_n=n/num_threads;
               rest=n-local_n*num_threads;
               n_start=local_n*th+MIN(th,rest);
               n_end=local_n*(th+1)+MIN(th+1,rest);
               memcpy(&R_PRES[0][n_start],&r[n_start],sizeof(double)*(n_end-n_start));
               memcpy(&X_PRES[0][n_start],&x[n_start],sizeof(double)*(n_end-n_start));
             }
             num_iter_restart=0;
             restartFlag=FALSE;
         }
         ++num_iter;
         ++num_iter_restart;
         /* order is the dimension of the space on which the residual is minimized: */
         order=MIN(num_iter_restart,Length_of_recursion);
         /***                                                                 
         *** calculate new search direction P from R_PRES
         ***/
         Paso_Solver_solvePreconditioner(A,&P_PRES[0][0], &R_PRES[0][0]);
         /***                                                                 
         *** apply A to P to get AP 
         ***/
	 Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(PASO_ONE, A, &P_PRES[0][0],PASO_ZERO, &AP[0]);
         /***                                                                 
         ***** calculation of the norm of R and the scalar products of       
         ***   the residuals and A*P:                                        
         ***/
         memset(loc_dots,0,sizeof(double)*(order+1));
         #pragma omp parallel for private(th,z,i,local_n, rest, n_start, n_end, R_PRES_dot_P, P_PRES_dot_AP0, P_PRES_dot_AP1, P_PRES_dot_AP2, P_PRES_dot_AP3, P_PRES_dot_AP4, P_PRES_dot_AP5, P_PRES_dot_AP6)
         for (th=0;th<num_threads;++th) {
               local_n=n/num_threads;
               rest=n-local_n*num_threads;
               n_start=local_n*th+MIN(th,rest);
               n_end=local_n*(th+1)+MIN(th+1,rest);
               if (order==0) {
                   R_PRES_dot_P=PASO_ZERO;
                   #pragma ivdep
                   for (z=n_start; z < n_end; ++z) {
                       R_PRES_dot_P+=R_PRES[0][z]*P_PRES[0][z];
                   }
                   #pragma omp critical
                   {
                        loc_dots[0]+=R_PRES_dot_P;
                   }
               } else if (order==1) {
                   R_PRES_dot_P=PASO_ZERO;
                   P_PRES_dot_AP0=PASO_ZERO;
                   #pragma ivdep
                   for (z=n_start; z < n_end; ++z) {
                       R_PRES_dot_P+=R_PRES[0][z]*P_PRES[0][z];
                       P_PRES_dot_AP0+=P_PRES[0][z]*AP[z];
                   }
                   #pragma omp critical
                   {
                        loc_dots[0]+=R_PRES_dot_P;
                        loc_dots[1]+=P_PRES_dot_AP0;
                   }
               } else if (order==2) {
                   R_PRES_dot_P=PASO_ZERO;
                   P_PRES_dot_AP0=PASO_ZERO;
                   P_PRES_dot_AP1=PASO_ZERO;
                   #pragma ivdep
                   for (z=n_start; z < n_end; ++z) {
                       R_PRES_dot_P+=R_PRES[0][z]*P_PRES[0][z];
                       P_PRES_dot_AP0+=P_PRES[0][z]*AP[z];
                       P_PRES_dot_AP1+=P_PRES[1][z]*AP[z];
                   }
                   #pragma omp critical
                   {
                        loc_dots[0]+=R_PRES_dot_P;
                        loc_dots[1]+=P_PRES_dot_AP0;
                        loc_dots[2]+=P_PRES_dot_AP1;
                   }
               } else if (order==3) {
                   R_PRES_dot_P=PASO_ZERO;
                   P_PRES_dot_AP0=PASO_ZERO;
                   P_PRES_dot_AP1=PASO_ZERO;
                   P_PRES_dot_AP2=PASO_ZERO;
                   #pragma ivdep
                   for (z=n_start; z < n_end; ++z) {
                       R_PRES_dot_P+=R_PRES[0][z]*P_PRES[0][z];
                       P_PRES_dot_AP0+=P_PRES[0][z]*AP[z];
                       P_PRES_dot_AP1+=P_PRES[1][z]*AP[z];
                       P_PRES_dot_AP2+=P_PRES[2][z]*AP[z];
                   }
                   #pragma omp critical
                   {
                        loc_dots[0]+=R_PRES_dot_P;
                        loc_dots[1]+=P_PRES_dot_AP0;
                        loc_dots[2]+=P_PRES_dot_AP1;
                        loc_dots[3]+=P_PRES_dot_AP2;
                   }
               } else if (order==4) {
                   R_PRES_dot_P=PASO_ZERO;
                   P_PRES_dot_AP0=PASO_ZERO;
                   P_PRES_dot_AP1=PASO_ZERO;
                   P_PRES_dot_AP2=PASO_ZERO;
                   P_PRES_dot_AP3=PASO_ZERO;
                   #pragma ivdep
                   for (z=n_start; z < n_end; ++z) {
                       R_PRES_dot_P+=R_PRES[0][z]*P_PRES[0][z];
                       P_PRES_dot_AP0+=P_PRES[0][z]*AP[z];
                       P_PRES_dot_AP1+=P_PRES[1][z]*AP[z];
                       P_PRES_dot_AP2+=P_PRES[2][z]*AP[z];
                       P_PRES_dot_AP3+=P_PRES[3][z]*AP[z];
                   }
                   #pragma omp critical
                   {
                        loc_dots[0]+=R_PRES_dot_P;
                        loc_dots[1]+=P_PRES_dot_AP0;
                        loc_dots[2]+=P_PRES_dot_AP1;
                        loc_dots[3]+=P_PRES_dot_AP2;
                        loc_dots[4]+=P_PRES_dot_AP3;
                   }
               } else if (order==5) {
                   R_PRES_dot_P=PASO_ZERO;
                   P_PRES_dot_AP0=PASO_ZERO;
                   P_PRES_dot_AP1=PASO_ZERO;
                   P_PRES_dot_AP2=PASO_ZERO;
                   P_PRES_dot_AP3=PASO_ZERO;
                   P_PRES_dot_AP4=PASO_ZERO;
                   #pragma ivdep
                   for (z=n_start; z < n_end; ++z) {
                       R_PRES_dot_P+=R_PRES[0][z]*P_PRES[0][z];
                       P_PRES_dot_AP0+=P_PRES[0][z]*AP[z];
                       P_PRES_dot_AP1+=P_PRES[1][z]*AP[z];
                       P_PRES_dot_AP2+=P_PRES[2][z]*AP[z];
                       P_PRES_dot_AP3+=P_PRES[3][z]*AP[z];
                       P_PRES_dot_AP4+=P_PRES[4][z]*AP[z];
                   }
                   #pragma omp critical
                   {
                        loc_dots[0]+=R_PRES_dot_P;
                        loc_dots[1]+=P_PRES_dot_AP0;
                        loc_dots[2]+=P_PRES_dot_AP1;
                        loc_dots[3]+=P_PRES_dot_AP2;
                        loc_dots[4]+=P_PRES_dot_AP3;
                        loc_dots[5]+=P_PRES_dot_AP4;
                   }
               } else if (order==6) {
                   R_PRES_dot_P=PASO_ZERO;
                   P_PRES_dot_AP0=PASO_ZERO;
                   P_PRES_dot_AP1=PASO_ZERO;
                   P_PRES_dot_AP2=PASO_ZERO;
                   P_PRES_dot_AP3=PASO_ZERO;
                   P_PRES_dot_AP4=PASO_ZERO;
                   P_PRES_dot_AP5=PASO_ZERO;
                   #pragma ivdep
                   for (z=n_start; z < n_end; ++z) {
                       R_PRES_dot_P+=R_PRES[0][z]*P_PRES[0][z];
                       P_PRES_dot_AP0+=P_PRES[0][z]*AP[z];
                       P_PRES_dot_AP1+=P_PRES[1][z]*AP[z];
                       P_PRES_dot_AP2+=P_PRES[2][z]*AP[z];
                       P_PRES_dot_AP3+=P_PRES[3][z]*AP[z];
                       P_PRES_dot_AP4+=P_PRES[4][z]*AP[z];
                       P_PRES_dot_AP5+=P_PRES[5][z]*AP[z];
                   }
                   #pragma omp critical
                   {
                        loc_dots[0]+=R_PRES_dot_P;
                        loc_dots[1]+=P_PRES_dot_AP0;
                        loc_dots[2]+=P_PRES_dot_AP1;
                        loc_dots[3]+=P_PRES_dot_AP2;
                        loc_dots[4]+=P_PRES_dot_AP3;
                        loc_dots[5]+=P_PRES_dot_AP4;
                        loc_dots[6]+=P_PRES_dot_AP5;
                   }
               } else {
                  R_PRES_dot_P=PASO_ZERO;
                  P_PRES_dot_AP0=PASO_ZERO;
                  P_PRES_dot_AP1=PASO_ZERO;
                  P_PRES_dot_AP2=PASO_ZERO;
                  P_PRES_dot_AP3=PASO_ZERO;
                  P_PRES_dot_AP4=PASO_ZERO;
                  P_PRES_dot_AP5=PASO_ZERO;
                  P_PRES_dot_AP6=PASO_ZERO;
                  #pragma ivdep
                  for (z=n_start; z < n_end; ++z) {
                      R_PRES_dot_P+=R_PRES[0][z]*P_PRES[0][z];
                      P_PRES_dot_AP0+=P_PRES[0][z]*AP[z];
                      P_PRES_dot_AP1+=P_PRES[1][z]*AP[z];
                      P_PRES_dot_AP2+=P_PRES[2][z]*AP[z];
                      P_PRES_dot_AP3+=P_PRES[3][z]*AP[z];
                      P_PRES_dot_AP4+=P_PRES[4][z]*AP[z];
                      P_PRES_dot_AP5+=P_PRES[5][z]*AP[z];
                      P_PRES_dot_AP6+=P_PRES[6][z]*AP[z];
                  }
                  #pragma omp critical
                  {
                       loc_dots[0]+=R_PRES_dot_P;
                       loc_dots[1]+=P_PRES_dot_AP0;
                       loc_dots[2]+=P_PRES_dot_AP1;
                       loc_dots[3]+=P_PRES_dot_AP2;
                       loc_dots[4]+=P_PRES_dot_AP3;
                       loc_dots[5]+=P_PRES_dot_AP4;
                       loc_dots[6]+=P_PRES_dot_AP5;
                       loc_dots[7]+=P_PRES_dot_AP6;
                  }
                  for (i=7;i<order;++i) {
                       P_PRES_dot_AP0=PASO_ZERO;
                       #pragma ivdep
                       for (z=n_start; z < n_end; ++z) {
                             P_PRES_dot_AP0+=P_PRES[i][z]*AP[z];
                       }
                       #pragma omp critical
                       {
                           loc_dots[i+1]+=P_PRES_dot_AP0;
                       }
                  }
               }
         }
         #ifdef PASO_MPI
                MPI_Allreduce(loc_dots, dots, order+1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
                R_PRES_dot_P_PRES[0]=dots[0];
                memcpy(P_PRES_dot_AP,&dots[1],sizeof(double)*order);
         #else
                R_PRES_dot_P_PRES[0]=loc_dots[0];
                memcpy(P_PRES_dot_AP,&loc_dots[1],sizeof(double)*order);
         #endif
         R_PRES_dot_AP0=R_PRES_dot_P_PRES[0];

         /***   if sum_BREAKF is equal to zero a breakdown occurs.
          ***   iteration procedure can be continued but R_PRES is not the
          ***   residual of X_PRES approximation.
          ***/
         sum_BREAKF=0.;
         #pragma ivdep
         for (i=0;i<order;++i) sum_BREAKF +=BREAKF[i];
         breakFlag=!((ABS(R_PRES_dot_AP0) > PASO_ZERO) &&  (sum_BREAKF >PASO_ZERO));
         if (!breakFlag) {
            breakFlag=FALSE;
            /***
            ***** X_PRES and R_PRES are moved to memory:
            ***/
            Factor=0.;
            #pragma ivdep
            for (i=0;i<order;++i) {
                   ALPHA[i]=-P_PRES_dot_AP[i]/R_PRES_dot_P_PRES[i];
                   Factor+=BREAKF[i]*ALPHA[i];
            }

             save_R_PRES_dot_P_PRES=R_PRES_dot_P_PRES[Length_of_mem-1];
             save_R_PRES=R_PRES[Length_of_mem-1];
             save_XPRES=X_PRES[Length_of_mem-1];
             save_P_PRES=P_PRES[Length_of_mem-1];
             #pragma ivdep
             for (i=Length_of_mem-1;i>0;--i) {
                   BREAKF[i]=BREAKF[i-1];
                   R_PRES_dot_P_PRES[i]=R_PRES_dot_P_PRES[i-1];
                   R_PRES[i]=R_PRES[i-1];
                   X_PRES[i]=X_PRES[i-1];
                   P_PRES[i]=P_PRES[i-1];
             }
             R_PRES_dot_P_PRES[0]=save_R_PRES_dot_P_PRES;
             R_PRES[0]=save_R_PRES;
             X_PRES[0]=save_XPRES;
             P_PRES[0]=save_P_PRES;

             if (ABS(Factor)<=PASO_ZERO) {
                  Factor=1.;
                  BREAKF[0]=PASO_ZERO;
             } else {
                  Factor=1./Factor;
                  BREAKF[0]=PASO_ONE;
            }
            for (i=0;i<order;++i) ALPHA[i]*=Factor;
            /*                                                                 
            ***** update of solution X_PRES and its residual R_PRES:            
            ***                                                               
            ***   P is used to accumulate X and AP to accumulate R. X and R     
            ***   are still needed until they are put into the X and R memory   
            ***   R_PRES and X_PRES                                     
            ***                                                                 
            **/
            breakf0=BREAKF[0];
            #pragma omp parallel for private(th,z,i,local_n, rest, n_start, n_end)
            for (th=0;th<num_threads;++th) {
               local_n=n/num_threads;
               rest=n-local_n*num_threads;
               n_start=local_n*th+MIN(th,rest);
               n_end=local_n*(th+1)+MIN(th+1,rest);
               if (order==0) {
                   #pragma ivdep
                   for (z=n_start; z < n_end; ++z) {
                      R_PRES[0][z]= Factor*       AP[z];
                      X_PRES[0][z]=-Factor*P_PRES[1][z];
                   }
               } else if (order==1) {
                   #pragma ivdep
                   for (z=n_start; z < n_end; ++z) {
                       R_PRES[0][z]= Factor*       AP[z]+ALPHA[0]*R_PRES[1][z];
                       X_PRES[0][z]=-Factor*P_PRES[1][z]+ALPHA[0]*X_PRES[1][z];
                   }
               } else if (order==2) {
                   #pragma ivdep
                   for (z=n_start; z < n_end; ++z) {
                       R_PRES[0][z]= Factor*       AP[z]+ALPHA[0]*R_PRES[1][z]
                                                        +ALPHA[1]*R_PRES[2][z];
                       X_PRES[0][z]=-Factor*P_PRES[1][z]+ALPHA[0]*X_PRES[1][z]
                                                        +ALPHA[1]*X_PRES[2][z];
                   }
               } else if (order==3) {
                   #pragma ivdep
                   for (z=n_start; z < n_end; ++z) {
                       R_PRES[0][z]= Factor*AP[z]+ALPHA[0]*R_PRES[1][z]
                                                 +ALPHA[1]*R_PRES[2][z]
                                                 +ALPHA[2]*R_PRES[3][z];
                       X_PRES[0][z]=-Factor*P_PRES[1][z]+ALPHA[0]*X_PRES[1][z]
                                                        +ALPHA[1]*X_PRES[2][z]
                                                        +ALPHA[2]*X_PRES[3][z];
                   }
               } else if (order==4) {
                   #pragma ivdep
                   for (z=n_start; z < n_end; ++z) {
                       R_PRES[0][z]= Factor*AP[z]+ALPHA[0]*R_PRES[1][z]
                                                 +ALPHA[1]*R_PRES[2][z]
                                                 +ALPHA[2]*R_PRES[3][z]
                                                 +ALPHA[3]*R_PRES[4][z];
                       X_PRES[0][z]=-Factor*P_PRES[1][z]+ALPHA[0]*X_PRES[1][z]
                                                        +ALPHA[1]*X_PRES[2][z]
                                                        +ALPHA[2]*X_PRES[3][z]
                                                        +ALPHA[3]*X_PRES[4][z];
                   }
               } else if (order==5) {
                   #pragma ivdep
                   for (z=n_start; z < n_end; ++z) {
                       R_PRES[0][z]=Factor*AP[z]+ALPHA[0]*R_PRES[1][z]
                                                +ALPHA[1]*R_PRES[2][z]
                                                +ALPHA[2]*R_PRES[3][z]
                                                +ALPHA[3]*R_PRES[4][z]
                                                +ALPHA[4]*R_PRES[5][z];
                       X_PRES[0][z]=-Factor*P_PRES[1][z]+ALPHA[0]*X_PRES[1][z]
                                                        +ALPHA[1]*X_PRES[2][z]
                                                        +ALPHA[2]*X_PRES[3][z]
                                                        +ALPHA[3]*X_PRES[4][z]
                                                        +ALPHA[4]*X_PRES[5][z];
                    }
               } else if (order==6) {
                   #pragma ivdep
                   for (z=n_start; z < n_end; ++z) {
                      R_PRES[0][z]=Factor*AP[z]+ALPHA[0]*R_PRES[1][z]
                                               +ALPHA[1]*R_PRES[2][z]
                                               +ALPHA[2]*R_PRES[3][z]
                                               +ALPHA[3]*R_PRES[4][z]
                                               +ALPHA[4]*R_PRES[5][z]
                                               +ALPHA[5]*R_PRES[6][z];
                      X_PRES[0][z]=-Factor*P_PRES[1][z]+ALPHA[0]*X_PRES[1][z]
                                                       +ALPHA[1]*X_PRES[2][z]
                                                       +ALPHA[2]*X_PRES[3][z]
                                                       +ALPHA[3]*X_PRES[4][z]
                                                       +ALPHA[4]*X_PRES[5][z]
                                                       +ALPHA[5]*X_PRES[6][z];
                   }
               } else {
                   #pragma ivdep
                   for (z=n_start; z < n_end; ++z) {
                        R_PRES[0][z]=Factor*AP[z]+ALPHA[0]*R_PRES[1][z]
                                                 +ALPHA[1]*R_PRES[2][z]
                                                 +ALPHA[2]*R_PRES[3][z]
                                                 +ALPHA[3]*R_PRES[4][z]
                                                 +ALPHA[4]*R_PRES[5][z]
                                                 +ALPHA[5]*R_PRES[6][z]
                                                 +ALPHA[6]*R_PRES[7][z];
                        X_PRES[0][z]=-Factor*P_PRES[1][z]+ALPHA[0]*X_PRES[1][z]
                                                         +ALPHA[1]*X_PRES[2][z]
                                                         +ALPHA[2]*X_PRES[3][z]
                                                         +ALPHA[3]*X_PRES[4][z]
                                                         +ALPHA[4]*X_PRES[5][z]
                                                         +ALPHA[5]*X_PRES[6][z]
                                                         +ALPHA[6]*X_PRES[7][z];
                   }
                   for (i=7;i<order;++i) {
                      #pragma ivdep
                      for (z=n_start; z < n_end; ++z) {
                         R_PRES[0][z]+=ALPHA[i]*R_PRES[i+1][z];
                         X_PRES[0][z]+=ALPHA[i]*X_PRES[i+1][z];
                     }
                   }
               }
            }
            if (breakf0>0.) {
              /***
              ***** calculate gamma from min_(gamma){|R+gamma*(R_PRES-R)|_2}:
              ***/
              memset(loc_dots,0,sizeof(double)*3);
              #pragma omp parallel for private(th,z,i,local_n, rest, n_start, n_end,diff,SC1,SC2)
              for (th=0;th<num_threads;++th) {
                  local_n=n/num_threads;
                  rest=n-local_n*num_threads;
                  n_start=local_n*th+MIN(th,rest);
                  n_end=local_n*(th+1)+MIN(th+1,rest);
                  SC1=PASO_ZERO;
                  SC2=PASO_ZERO;
                  #pragma ivdep
                  for (z=n_start; z < n_end; ++z) {
                       diff=R_PRES[0][z]-r[z];
                       SC1+=diff*diff;
                       SC2+=diff*r[z];
                  }
                  #pragma omp critical
                  {
                       loc_dots[0]+=SC1;
                       loc_dots[1]+=SC2;
                  }
              }
              #ifdef PASO_MPI
                  MPI_Allreduce(loc_dots, dots, 2, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
                  SC1=dots[0];
                  SC2=dots[1];
              #else
                  SC1=loc_dots[0];
                  SC2=loc_dots[1];
              #endif
              gamma=(SC1<=PASO_ZERO) ? PASO_ZERO : -SC2/SC1;
              #pragma omp parallel for private(th,z,local_n, rest, n_start, n_end,diff,L2_R)
              for (th=0;th<num_threads;++th) {
                  local_n=n/num_threads;
                  rest=n-local_n*num_threads;
                  n_start=local_n*th+MIN(th,rest);
                  n_end=local_n*(th+1)+MIN(th+1,rest);
                  L2_R=PASO_ZERO;
                  #pragma ivdep
                  for (z=n_start; z < n_end; ++z) {
                      x[z]+=gamma*(X_PRES[0][z]-x[z]);
                      r[z]+=gamma*(R_PRES[0][z]-r[z]);
                      L2_R+=r[z]*r[z];
                  }
                  #pragma omp critical
                  {
                       loc_dots[2]+=L2_R;
                  }
              }
              #ifdef PASO_MPI
                  MPI_Allreduce(&loc_dots[2], &dots[2], 1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
                  L2_R=dots[2];
              #else
                  L2_R=loc_dots[2];
              #endif
              norm_of_residual=sqrt(L2_R);
              convergeFlag = (norm_of_residual <= tol);
              if (restart>0) restartFlag=(num_iter_restart >= restart);
            } else { 
              convergeFlag=FALSE;
              restartFlag=FALSE;
            }
            maxIterFlag = (num_iter >= maxit);
         }
        }
        /* end of iteration */
        Norm_of_residual_global=norm_of_residual;
	Num_iter_global=num_iter;
	if (maxIterFlag) { 
	       Status = SOLVER_MAXITER_REACHED;
	   } else if (breakFlag) {
	       Status = SOLVER_BREAKDOWN;
	}
    }
    for (i=0;i<Length_of_recursion;i++) {
       TMPMEMFREE(X_PRES[i]);
       TMPMEMFREE(R_PRES[i]);
       TMPMEMFREE(P_PRES[i]);
    }
    TMPMEMFREE(AP);
    TMPMEMFREE(X_PRES);
    TMPMEMFREE(R_PRES);
    TMPMEMFREE(P_PRES);
    TMPMEMFREE(P_PRES_dot_AP);
    TMPMEMFREE(R_PRES_dot_P_PRES);
    TMPMEMFREE(BREAKF);
    TMPMEMFREE(ALPHA);
    TMPMEMFREE(dots);
    TMPMEMFREE(loc_dots);
    *iter=Num_iter_global;
    *tolerance=Norm_of_residual_global;
  return Status;
}

