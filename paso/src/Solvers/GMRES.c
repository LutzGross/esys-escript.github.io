/* $Id$ */

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
    double * tolerance,dim_t Length_of_recursion,dim_t restart) {

  /* Local variables */

  double *AP,*X_PRES[MAX(Length_of_recursion,0)+1],*R_PRES[MAX(Length_of_recursion,0)+1],*P_PRES[MAX(Length_of_recursion,0)+1];
  double P_PRES_dot_AP[MAX(Length_of_recursion,0)],R_PRES_dot_P_PRES[MAX(Length_of_recursion,0)+1],BREAKF[MAX(Length_of_recursion,0)+1],ALPHA[MAX(Length_of_recursion,0)];
  double P_PRES_dot_AP0,P_PRES_dot_AP1,P_PRES_dot_AP2,P_PRES_dot_AP3,P_PRES_dot_AP4,P_PRES_dot_AP5,P_PRES_dot_AP6,R_PRES_dot_P,abs_RP,breakf0;
  double tol,Factor,sum_BREAKF,gamma,SC1,SC2,norm_of_residual,diff,L2_R,Norm_of_residual_global;
  double *save_XPRES, *save_P_PRES, *save_R_PRES,save_R_PRES_dot_P_PRES;
  dim_t maxit,Num_iter_global,num_iter_restart,num_iter;
  dim_t i,z,order;
  bool_t breakFlag=FALSE, maxIterFlag=FALSE, convergeFlag=FALSE,restartFlag=FALSE;
  err_t Status=SOLVER_NO_ERROR;

  /* adapt original routine parameters */

  dim_t n=A->num_cols * A-> col_block_size;
  dim_t Length_of_mem=MAX(Length_of_recursion,0)+1;

  /*     Test the input parameters. */
  if (restart>0) restart=MAX(Length_of_recursion,restart);
  if (n < 0 || Length_of_recursion<=0) {
    return SOLVER_INPUT_ERROR;
  }

  /*     allocate memory: */

  AP=TMPMEMALLOC(n,double);
  if (AP==NULL) Status=SOLVER_MEMORY_ERROR;
  for (i=0;i<Length_of_mem;i++) {
    X_PRES[i]=TMPMEMALLOC(n,double);
    R_PRES[i]=TMPMEMALLOC(n,double);
    P_PRES[i]=TMPMEMALLOC(n,double);
    if (X_PRES[i]==NULL || R_PRES[i]==NULL ||  P_PRES[i]==NULL) Status=SOLVER_MEMORY_ERROR;
  }
  if ( Status ==SOLVER_NO_ERROR ) {

    /* now PRES starts : */
    maxit=*iter;
    tol=*tolerance;

    #pragma omp parallel firstprivate(maxit,tol,convergeFlag,maxIterFlag,breakFlag) \
       private(num_iter,i,num_iter_restart,order,sum_BREAKF,gamma,restartFlag,norm_of_residual,abs_RP,breakf0,\
               save_XPRES,save_P_PRES,save_R_PRES,save_R_PRES_dot_P_PRES)
    {
      /* initialization */

      restartFlag=TRUE;
      num_iter=0;
      #pragma omp for private(z) schedule(static) nowait
      for (z=0; z < n; ++z) AP[z]=0;
      for(i=0;i<Length_of_mem;++i) {
        #pragma omp for private(z) schedule(static) nowait
        for (z=0; z < n; ++z) {
                   P_PRES[i][z]=0;
                   R_PRES[i][z]=0;
                   X_PRES[i][z]=0;
        }
      }

      while (! (convergeFlag || maxIterFlag || breakFlag))  {
         #pragma omp barrier
         if (restartFlag) {
             #pragma omp master
             BREAKF[0]=ONE;
             #pragma omp for private(z) schedule(static) nowait
             for (z=0; z < n; ++z) {
                R_PRES[0][z]=r[z];
                X_PRES[0][z]=x[z];
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
	 Paso_SystemMatrix_MatrixVector(ONE, A, &P_PRES[0][0],ZERO, &AP[0]);
         /***                                                                 
         ***** calculation of the norm of R and the scalar products of       
         ***   the residuals and A*P:                                        
         ***/
         if (order==0) {
            #pragma omp master
            R_PRES_dot_P=ZERO;
            #pragma omp barrier
            #pragma omp for private(z) reduction(+:R_PRES_dot_P) schedule(static)
            for (z=0;z<n;++z) 
                R_PRES_dot_P+=R_PRES[0][z]*P_PRES[0][z];
            #pragma omp master
            R_PRES_dot_P_PRES[0]=R_PRES_dot_P;
         } else if (order==1) {
            #pragma omp master
            {
               R_PRES_dot_P=ZERO;
               P_PRES_dot_AP0=ZERO;
            }
            #pragma omp barrier
            #pragma omp for private(z) reduction(+:R_PRES_dot_P,P_PRES_dot_AP0) schedule(static)
            for (z=0;z<n;++z) {
                R_PRES_dot_P+=R_PRES[0][z]*P_PRES[0][z];
                P_PRES_dot_AP0+=P_PRES[0][z]*AP[z];
            }
            #pragma omp master
            {
              P_PRES_dot_AP[0]=P_PRES_dot_AP0;
              R_PRES_dot_P_PRES[0]=R_PRES_dot_P;
            }
         } else if (order==2) {
            #pragma omp master
            {
              R_PRES_dot_P=ZERO;
              P_PRES_dot_AP0=ZERO;
              P_PRES_dot_AP1=ZERO;
            }
            #pragma omp barrier
            #pragma omp for private(z) reduction(+:R_PRES_dot_P,P_PRES_dot_AP0,P_PRES_dot_AP1) schedule(static)
            for (z=0;z<n;++z) {
                R_PRES_dot_P+=R_PRES[0][z]*P_PRES[0][z];
                P_PRES_dot_AP0+=P_PRES[0][z]*AP[z];
                P_PRES_dot_AP1+=P_PRES[1][z]*AP[z];
            }
            #pragma omp master
            {
              P_PRES_dot_AP[0]=P_PRES_dot_AP0;
              P_PRES_dot_AP[1]=P_PRES_dot_AP1;
              R_PRES_dot_P_PRES[0]=R_PRES_dot_P;
            }
         } else if (order==3) {
            #pragma omp master
            {
               R_PRES_dot_P=ZERO;
               P_PRES_dot_AP0=ZERO;
               P_PRES_dot_AP1=ZERO;
               P_PRES_dot_AP2=ZERO;
            }
            #pragma omp barrier
            #pragma omp for private(z) reduction(+:R_PRES_dot_P,P_PRES_dot_AP0,P_PRES_dot_AP1,P_PRES_dot_AP2) schedule(static)
            for (z=0;z<n;++z) {
                R_PRES_dot_P+=R_PRES[0][z]*P_PRES[0][z];
                P_PRES_dot_AP0+=P_PRES[0][z]*AP[z];
                P_PRES_dot_AP1+=P_PRES[1][z]*AP[z];
                P_PRES_dot_AP2+=P_PRES[2][z]*AP[z];
            }
            #pragma omp master
            {
              P_PRES_dot_AP[0]=P_PRES_dot_AP0;
              P_PRES_dot_AP[1]=P_PRES_dot_AP1;
              P_PRES_dot_AP[2]=P_PRES_dot_AP2;
              R_PRES_dot_P_PRES[0]=R_PRES_dot_P;
            }
         } else if (order==4) {
            #pragma omp master
            {
              R_PRES_dot_P=ZERO;
              P_PRES_dot_AP0=ZERO;
              P_PRES_dot_AP1=ZERO;
              P_PRES_dot_AP2=ZERO;
              P_PRES_dot_AP3=ZERO;
            }
            #pragma omp barrier
            #pragma omp for private(z) reduction(+:R_PRES_dot_P,P_PRES_dot_AP0,P_PRES_dot_AP1,P_PRES_dot_AP2,P_PRES_dot_AP3) schedule(static)
            for (z=0;z<n;++z) {
                R_PRES_dot_P+=R_PRES[0][z]*P_PRES[0][z];
                P_PRES_dot_AP0+=P_PRES[0][z]*AP[z];
                P_PRES_dot_AP1+=P_PRES[1][z]*AP[z];
                P_PRES_dot_AP2+=P_PRES[2][z]*AP[z];
                P_PRES_dot_AP3+=P_PRES[3][z]*AP[z];
            }
            #pragma omp master
            {
              P_PRES_dot_AP[0]=P_PRES_dot_AP0;
              P_PRES_dot_AP[1]=P_PRES_dot_AP1;
              P_PRES_dot_AP[2]=P_PRES_dot_AP2;
              P_PRES_dot_AP[3]=P_PRES_dot_AP3;
              R_PRES_dot_P_PRES[0]=R_PRES_dot_P;
            }
         } else if (order==5) {
            #pragma omp master
            {
              R_PRES_dot_P=ZERO;
              P_PRES_dot_AP0=ZERO;
              P_PRES_dot_AP1=ZERO;
              P_PRES_dot_AP2=ZERO;
              P_PRES_dot_AP3=ZERO;
              P_PRES_dot_AP4=ZERO;
            }
            #pragma omp barrier
            #pragma omp for private(z) reduction(+:R_PRES_dot_P,P_PRES_dot_AP0,P_PRES_dot_AP1,P_PRES_dot_AP2,P_PRES_dot_AP3,P_PRES_dot_AP4) schedule(static)
            for (z=0;z<n;++z) {
                R_PRES_dot_P+=R_PRES[0][z]*P_PRES[0][z];
                P_PRES_dot_AP0+=P_PRES[0][z]*AP[z];
                P_PRES_dot_AP1+=P_PRES[1][z]*AP[z];
                P_PRES_dot_AP2+=P_PRES[2][z]*AP[z];
                P_PRES_dot_AP3+=P_PRES[3][z]*AP[z];
                P_PRES_dot_AP4+=P_PRES[4][z]*AP[z];
            }
            #pragma omp master
            {
              P_PRES_dot_AP[0]=P_PRES_dot_AP0;
              P_PRES_dot_AP[1]=P_PRES_dot_AP1;
              P_PRES_dot_AP[2]=P_PRES_dot_AP2;
              P_PRES_dot_AP[3]=P_PRES_dot_AP3;
              P_PRES_dot_AP[4]=P_PRES_dot_AP4;
              R_PRES_dot_P_PRES[0]=R_PRES_dot_P;
            }
         } else if (order==6) {
            #pragma omp master
            {
              R_PRES_dot_P=ZERO;
              P_PRES_dot_AP0=ZERO;
              P_PRES_dot_AP1=ZERO;
              P_PRES_dot_AP2=ZERO;
              P_PRES_dot_AP3=ZERO;
              P_PRES_dot_AP4=ZERO;
              P_PRES_dot_AP5=ZERO;
            }
            #pragma omp barrier
            #pragma omp for private(z) reduction(+:R_PRES_dot_P,P_PRES_dot_AP0,P_PRES_dot_AP1,P_PRES_dot_AP2,P_PRES_dot_AP3,P_PRES_dot_AP4,P_PRES_dot_AP5) schedule(static)
            for (z=0;z<n;++z) {
                R_PRES_dot_P+=R_PRES[0][z]*P_PRES[0][z];
                P_PRES_dot_AP0+=P_PRES[0][z]*AP[z];
                P_PRES_dot_AP1+=P_PRES[1][z]*AP[z];
                P_PRES_dot_AP2+=P_PRES[2][z]*AP[z];
                P_PRES_dot_AP3+=P_PRES[3][z]*AP[z];
                P_PRES_dot_AP4+=P_PRES[4][z]*AP[z];
                P_PRES_dot_AP5+=P_PRES[5][z]*AP[z];
             }
            #pragma omp master
            {
              P_PRES_dot_AP[0]=P_PRES_dot_AP0;
              P_PRES_dot_AP[1]=P_PRES_dot_AP1;
              P_PRES_dot_AP[2]=P_PRES_dot_AP2;
              P_PRES_dot_AP[3]=P_PRES_dot_AP3;
              P_PRES_dot_AP[4]=P_PRES_dot_AP4;
              P_PRES_dot_AP[5]=P_PRES_dot_AP5;
              R_PRES_dot_P_PRES[0]=R_PRES_dot_P;
            }
         } else if (order>6) {
            #pragma omp master
            {
              R_PRES_dot_P=ZERO;
              P_PRES_dot_AP0=ZERO;
              P_PRES_dot_AP1=ZERO;
              P_PRES_dot_AP2=ZERO;
              P_PRES_dot_AP3=ZERO;
              P_PRES_dot_AP4=ZERO;
              P_PRES_dot_AP5=ZERO;
              P_PRES_dot_AP6=ZERO;
            }
            #pragma omp barrier
            #pragma omp for private(z) reduction(+:R_PRES_dot_P,P_PRES_dot_AP0,P_PRES_dot_AP1,P_PRES_dot_AP2,P_PRES_dot_AP3,P_PRES_dot_AP4,P_PRES_dot_AP5,P_PRES_dot_AP6) schedule(static)
            for (z=0;z<n;++z) {
                R_PRES_dot_P+=R_PRES[0][z]*P_PRES[0][z];
                P_PRES_dot_AP0+=P_PRES[0][z]*AP[z];
                P_PRES_dot_AP1+=P_PRES[1][z]*AP[z];
                P_PRES_dot_AP2+=P_PRES[2][z]*AP[z];
                P_PRES_dot_AP3+=P_PRES[3][z]*AP[z];
                P_PRES_dot_AP4+=P_PRES[4][z]*AP[z];
                P_PRES_dot_AP5+=P_PRES[5][z]*AP[z];
                P_PRES_dot_AP6+=P_PRES[6][z]*AP[z];
             }
            #pragma omp master
            {
              P_PRES_dot_AP[0]=P_PRES_dot_AP0;
              P_PRES_dot_AP[1]=P_PRES_dot_AP1;
              P_PRES_dot_AP[2]=P_PRES_dot_AP2;
              P_PRES_dot_AP[3]=P_PRES_dot_AP3;
              P_PRES_dot_AP[4]=P_PRES_dot_AP4;
              P_PRES_dot_AP[5]=P_PRES_dot_AP5;
              P_PRES_dot_AP[6]=P_PRES_dot_AP6;
              R_PRES_dot_P_PRES[0]=R_PRES_dot_P;

              P_PRES_dot_AP0=ZERO;
            }
            for (i=7;i<order;++i) {
                #pragma omp barrier
                #pragma omp for private(z) reduction(+:P_PRES_dot_AP0) schedule(static)
                for (z=0;z<n;++z) P_PRES_dot_AP0+=P_PRES[i][z]*AP[z];
                #pragma omp master
                {
                  P_PRES_dot_AP[i]=P_PRES_dot_AP0;
                  P_PRES_dot_AP0=ZERO;
                }
            }
         }
         /* this fixes a problem with the intel compiler */
         #pragma omp master
         P_PRES_dot_AP0=R_PRES_dot_P_PRES[0];
         #pragma omp barrier
         /***   if sum_BREAKF is equal to zero a breakdown occurs.
          ***   iteration procedure can be continued but R_PRES is not the
          ***   residual of X_PRES approximation.
          ***/
         sum_BREAKF=0.;
         for (i=0;i<order;++i) sum_BREAKF +=BREAKF[i];
         breakFlag=!((ABS(P_PRES_dot_AP0) > ZERO) &&  (sum_BREAKF >ZERO));
         if (!breakFlag) {
            breakFlag=FALSE;
            /***
            ***** X_PRES and R_PRES are moved to memory:
            ***/
            #pragma omp master
            {
               Factor=0.;
               for (i=0;i<order;++i) {
                   ALPHA[i]=-P_PRES_dot_AP[i]/R_PRES_dot_P_PRES[i];
                   Factor+=BREAKF[i]*ALPHA[i];
               }

               save_R_PRES_dot_P_PRES=R_PRES_dot_P_PRES[Length_of_mem-1];
               save_R_PRES=R_PRES[Length_of_mem-1];
               save_XPRES=X_PRES[Length_of_mem-1];
               save_P_PRES=P_PRES[Length_of_mem-1];
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

               if (ABS(Factor)<=ZERO) {
                  Factor=1.;
                  BREAKF[0]=ZERO;
               } else {
                  Factor=1./Factor;
                  BREAKF[0]=ONE;
               }
               for (i=0;i<order;++i) ALPHA[i]*=Factor;
            }
            /*                                                                 
            ***** update of solution X_PRES and its residual R_PRES:            
            ***                                                               
            ***   P is used to accumulate X and AP to accumulate R. X and R     
            ***   are still needed until they are put into the X and R memory   
            ***   R_PRES and X_PRES                                     
            ***                                                                 
            **/
            #pragma omp barrier
            breakf0=BREAKF[0];
            if (order==0) {
              #pragma omp for private(z) schedule(static)
              for (z=0;z<n;++z) {
                  R_PRES[0][z]= Factor*       AP[z];
                  X_PRES[0][z]=-Factor*P_PRES[1][z];
              }
            } else if (order==1) {
              #pragma omp for private(z) schedule(static)
              for (z=0;z<n;++z) {
                  R_PRES[0][z]= Factor*       AP[z]+ALPHA[0]*R_PRES[1][z];
                  X_PRES[0][z]=-Factor*P_PRES[1][z]+ALPHA[0]*X_PRES[1][z];
              }
            } else if (order==2) {
              #pragma omp for private(z) schedule(static)
              for (z=0;z<n;++z) {
                 R_PRES[0][z]= Factor*       AP[z]+ALPHA[0]*R_PRES[1][z]
                                                  +ALPHA[1]*R_PRES[2][z];
                 X_PRES[0][z]=-Factor*P_PRES[1][z]+ALPHA[0]*X_PRES[1][z]
                                                  +ALPHA[1]*X_PRES[2][z];
              }
            } else if (order==3) {
              #pragma omp for private(z) schedule(static)
              for (z=0;z<n;++z) {
                 R_PRES[0][z]= Factor*AP[z]+ALPHA[0]*R_PRES[1][z]
                                           +ALPHA[1]*R_PRES[2][z]
                                           +ALPHA[2]*R_PRES[3][z];
                 X_PRES[0][z]=-Factor*P_PRES[1][z]+ALPHA[0]*X_PRES[1][z]
                                                  +ALPHA[1]*X_PRES[2][z]
                                                  +ALPHA[2]*X_PRES[3][z];
              }
            } else if (order==4) {
              #pragma omp for private(z) schedule(static)
              for (z=0;z<n;++z) {
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
              #pragma omp for private(z) schedule(static)
              for (z=0;z<n;++z) {
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
              #pragma omp for private(z) schedule(static)
              for (z=0;z<n;++z) {
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
            } else if (order>6) {
              #pragma omp for private(z) schedule(static)
              for (z=0;z<n;++z) {
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
                #pragma omp for private(z) schedule(static)
                for (z=0;z<n;++z) {
                    R_PRES[0][z]+=ALPHA[i]*R_PRES[i+1][z];
                    X_PRES[0][z]+=ALPHA[i]*X_PRES[i+1][z];
                }
              }
            }
            if (breakf0>0.) {
              /***
              ***** calculate gamma from min_(gamma){|R+gamma*(R_PRES-R)|_2}:
              ***/
              #pragma omp master
              {
                 SC1=ZERO;
                 SC2=ZERO;
                 L2_R=ZERO;
              }
              #pragma omp barrier
              #pragma omp for private(z,diff) reduction(+:SC1,SC2) schedule(static)
              for (z=0;z<n;++z) {
                diff=R_PRES[0][z]-r[z];
                SC1+=diff*diff;
                SC2+=diff*r[z];
              }
              gamma=(SC1<=ZERO) ? ZERO : -SC2/SC1;
              #pragma omp for private(z) reduction(+:L2_R) schedule(static)
              for (z=0;z<n;++z) {
                 x[z]+=gamma*(X_PRES[0][z]-x[z]);
                 r[z]+=gamma*(R_PRES[0][z]-r[z]);
                 L2_R+=r[z]*r[z];
              }
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
        #pragma omp master 
        {
           Norm_of_residual_global=norm_of_residual;
	   Num_iter_global=num_iter;
	   if (maxIterFlag) { 
	       Status = SOLVER_MAXITER_REACHED;
	   } else if (breakFlag) {
	       Status = SOLVER_BREAKDOWN;
	  }
        }
    }  /* end of parallel region */
    TMPMEMFREE(AP);
    for (i=0;i<Length_of_recursion;i++) {
       TMPMEMFREE(X_PRES[i]);
       TMPMEMFREE(R_PRES[i]);
       TMPMEMFREE(P_PRES[i]);
    }
    *iter=Num_iter_global;
    *tolerance=Norm_of_residual_global;
  }
  return Status;
}

