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
*  r       (input/output) double array, dimension n_row.
*          On entry, residual of inital guess X
*
*  x       (input/output) double array, dimension n_col.
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
*  length_of_recursion (input) gives the number of residual to be kept in orthogonalization process
*
*  restart (input) If restart>0, iteration is resterted a after restart steps.
*
*  INFO    (output) int
*
*         =SOLVER_NO_ERROR: Successful exit. num_iterated approximate solution returned.
*         =SOLVER_MAXNUM_ITEr_REACHED
*         =SOLVER_INPUT_ERROR Illegal parameter:
*         =SOLVER_BREAKDOWN: If parameters RHO or OMEGA become smaller
*         =SOLVER_MEMORY_ERROR : If parameters RHO or OMEGA become smaller
*
*  ==============================================================
*/

static double ONE = 1.000000000000000;
static double ZERO = 0.000000000000000;

#include "Common.h"
#include "System.h"
#include "Solver.h"
#ifdef _OPENMP
#include <omp.h>
#endif

int Finley_Solver_GMRES(
    Finley_SystemMatrix * A,
    double * r,
    double * x,
    int *iter,
    double * tolerance,int length_of_recursion,int restart) {

  /* Local variables */

  double *x_PRES,*r_PRES,*P,*AP,*x_PRES_MEM[MAX(length_of_recursion,0)],*r_PRES_MEM[MAX(length_of_recursion,0)];
  double r_PRES_MEMdotAP[MAX(length_of_recursion,0)],L2_r_PRES_MEM[MAX(length_of_recursion,0)],BREAKF_MEM[MAX(length_of_recursion,0)];
  double r_PRESdotAP,r_PRES_MEMdotAP0,r_PRES_MEMdotAP1,r_PRES_MEMdotAP2,r_PRES_MEMdotAP3,r_PRES_MEMdotAP4,L2_r_PRES;
  double *tmp,tol,BREAKF,factor,GAMMA,SC1,SC2,norm_of_residual,diff,L2_R,norm_of_residual_global;
  int maxit,num_iter_global,num_iter_restart,num_iter;
  int i,z,OLDEST,ORDER;
  int breakFlag=FALSE, maxIterFlag=FALSE, convergeFlag=FALSE,restartFlag=FALSE;
  int status=SOLVER_NO_ERROR;

  /* adapt original routine parameters */

  int n_col=A->num_cols * A-> col_block_size;
  int n_row=A->num_rows * A-> row_block_size;

  /*     Test the input parameters. */
  if (restart>0) restart=MAX(length_of_recursion,restart);
  if (n_col < 0 || n_row<0 || length_of_recursion<=0) {
    return SOLVER_INPUT_ERROR;
  }

  /*     allocate memory: */

  x_PRES=TMPMEMALLOC(n_col,double);
  r_PRES=TMPMEMALLOC(n_row,double);
  P=TMPMEMALLOC(n_col,double);
  AP=TMPMEMALLOC(n_row,double);
  if (x_PRES==NULL || r_PRES==NULL || P==NULL || AP==NULL) status=SOLVER_MEMORY_ERROR;
  for (i=0;i<length_of_recursion;i++) {
    x_PRES_MEM[i]=TMPMEMALLOC(n_col,double);
    r_PRES_MEM[i]=TMPMEMALLOC(n_row,double);
    if (x_PRES_MEM[i]==NULL || r_PRES_MEM[i]==NULL) status=SOLVER_MEMORY_ERROR;
  }
  if ( status ==SOLVER_NO_ERROR ) {

    /* now PRES starts : */
    maxit=*iter;
    tol=*tolerance;
    norm_of_residual=tol;

    #pragma omp parallel firstprivate(maxit,tol) \
       private(num_iter,i,num_iter_restart,ORDER,OLDEST,BREAKF,factor,GAMMA,restartFlag,convergeFlag,maxIterFlag,breakFlag)
    {
      /* initialization */

      restartFlag=TRUE;
      num_iter=0;
      #pragma omp for private(z) schedule(static) nowait
      for (z=0; z < n_row; ++z) AP[z]=0;
      #pragma omp for private(z) schedule(static) nowait
      for (z=0; z < n_col; ++z) P[z]=0;
      for(i=0;i<length_of_recursion;++i) {
        #pragma omp for private(z) schedule(static) nowait
        for (z=0; z < n_row; ++z) r_PRES_MEM[i][z]=0;
        #pragma omp for private(z) schedule(static) nowait
        for (z=0; z < n_col; ++z) x_PRES_MEM[i][z]=0;
      }

      next:
         if (restartFlag) {
             #pragma omp for private(z) schedule(static) nowait
             for (z=0; z < n_row; ++z) r_PRES[z]=r[z];
             #pragma omp for private(z) schedule(static) nowait
             for (z=0; z < n_col; ++z) x_PRES[z]=x[z];
             num_iter_restart=0;
             OLDEST=length_of_recursion-1;
             BREAKF=ONE;
         }
         ++num_iter;
         ++num_iter_restart;
         /* ORDER is the dimension of the space on which the residual is minimized: */
         ORDER=MIN(num_iter_restart-1,length_of_recursion);
         /* OLDEST points to the oldest r_PRES and x_PRES in memory :*/
         OLDEST=(OLDEST==length_of_recursion-1) ? 0 : OLDEST+1;

         /***                                                                 
         *** calculate new search direction P from r_PRES
         ***/
         Finley_Solver_solvePreconditioner(A,&P[0], &r_PRES[0]);
         /***                                                                 
         *** apply A to P to get AP 
         ***/
	 Finley_RawScaledSystemMatrixVector(ONE, A, &P[0],ZERO, &AP[0]);
         /***                                                                 
         ***** calculation of the norm of R and the scalar products of       
         ***   the residuals and A*P:                                        
         ***/
         if (ORDER==0) {
            #pragma omp master
            {
               L2_r_PRES=ZERO;
               r_PRESdotAP=ZERO;
            }
            #pragma omp for private(z) reduction(+:L2_r_PRES,r_PRESdotAP)
            for (z=0;z<n_row;++z) {
                L2_r_PRES+=r_PRES[z]*r_PRES[z];
                r_PRESdotAP+=r_PRES[z]*AP[z];
            }
         } else if (ORDER==1) {
            #pragma omp master
            {
              L2_r_PRES=ZERO;
              r_PRESdotAP=ZERO;
              r_PRES_MEMdotAP0=ZERO;
            }
            #pragma omp for private(z) reduction(+:L2_r_PRES,r_PRESdotAP,r_PRES_MEMdotAP0)
            for (z=0;z<n_row;++z) {
                L2_r_PRES+=r_PRES[z]*r_PRES[z];
                r_PRESdotAP+=r_PRES[z]*AP[z];
                r_PRES_MEMdotAP0+=r_PRES_MEM[0][z]*AP[z];
            }
            #pragma omp master
            {
              r_PRES_MEMdotAP[0]=r_PRES_MEMdotAP0;
            }
   
         } else if (ORDER==2) {
            #pragma omp master
            {
               L2_r_PRES=ZERO;
               r_PRESdotAP=ZERO;
               r_PRES_MEMdotAP0=ZERO;
               r_PRES_MEMdotAP1=ZERO;
            }
            #pragma omp for private(z) reduction(+:L2_r_PRES,r_PRESdotAP,r_PRES_MEMdotAP0,r_PRES_MEMdotAP1)
            for (z=0;z<n_row;++z) {
                L2_r_PRES+=r_PRES[z]*r_PRES[z];
                r_PRESdotAP+=r_PRES[z]*AP[z];
                r_PRES_MEMdotAP0+=r_PRES_MEM[0][z]*AP[z];
                r_PRES_MEMdotAP1+=r_PRES_MEM[1][z]*AP[z];
            }
            #pragma omp master
            {
              r_PRES_MEMdotAP[0]=r_PRES_MEMdotAP0;
              r_PRES_MEMdotAP[1]=r_PRES_MEMdotAP1;
            }

         } else if (ORDER==3) {
            #pragma omp master
            {
              L2_r_PRES=ZERO;
              r_PRESdotAP=ZERO;
              r_PRES_MEMdotAP0=ZERO;
              r_PRES_MEMdotAP1=ZERO;
              r_PRES_MEMdotAP2=ZERO;
            }
            #pragma omp for private(z) reduction(+:L2_r_PRES,r_PRESdotAP,r_PRES_MEMdotAP0,r_PRES_MEMdotAP1,r_PRES_MEMdotAP2)
            for (z=0;z<n_row;++z) {
                L2_r_PRES+=r_PRES[z]*r_PRES[z];
                r_PRESdotAP+=r_PRES[z]*AP[z];
                r_PRES_MEMdotAP0+=r_PRES_MEM[0][z]*AP[z];
                r_PRES_MEMdotAP1+=r_PRES_MEM[1][z]*AP[z];
                r_PRES_MEMdotAP2+=r_PRES_MEM[2][z]*AP[z];
            }
            #pragma omp master
            {
              r_PRES_MEMdotAP[0]=r_PRES_MEMdotAP0;
              r_PRES_MEMdotAP[1]=r_PRES_MEMdotAP1;
              r_PRES_MEMdotAP[2]=r_PRES_MEMdotAP2;
            }
         } else if (ORDER==4) {
            #pragma omp master
            {
              L2_r_PRES=ZERO;
              r_PRESdotAP=ZERO;
              r_PRES_MEMdotAP0=ZERO;
              r_PRES_MEMdotAP1=ZERO;
              r_PRES_MEMdotAP2=ZERO;
              r_PRES_MEMdotAP3=ZERO;
            }
            #pragma omp for private(z) reduction(+:L2_r_PRES,r_PRESdotAP,r_PRES_MEMdotAP0,r_PRES_MEMdotAP1,r_PRES_MEMdotAP2,r_PRES_MEMdotAP3)
            for (z=0;z<n_row;++z) {
                L2_r_PRES+=r_PRES[z]*r_PRES[z];
                r_PRESdotAP+=r_PRES[z]*AP[z];
                r_PRES_MEMdotAP0+=r_PRES_MEM[0][z]*AP[z];
                r_PRES_MEMdotAP1+=r_PRES_MEM[1][z]*AP[z];
                r_PRES_MEMdotAP2+=r_PRES_MEM[2][z]*AP[z];
                r_PRES_MEMdotAP3+=r_PRES_MEM[3][z]*AP[z];
            }
            #pragma omp master
            {
              r_PRES_MEMdotAP[0]=r_PRES_MEMdotAP0;
              r_PRES_MEMdotAP[1]=r_PRES_MEMdotAP1;
              r_PRES_MEMdotAP[2]=r_PRES_MEMdotAP2;
              r_PRES_MEMdotAP[3]=r_PRES_MEMdotAP3;
            }
         } else if (ORDER==5) {
            #pragma omp master
            {
              L2_r_PRES=ZERO;
              r_PRESdotAP=ZERO;
              r_PRES_MEMdotAP0=ZERO;
              r_PRES_MEMdotAP1=ZERO;
              r_PRES_MEMdotAP2=ZERO;
              r_PRES_MEMdotAP3=ZERO;
              r_PRES_MEMdotAP4=ZERO;
            }
            #pragma omp for private(z) reduction(+:L2_r_PRES,r_PRESdotAP,r_PRES_MEMdotAP0,r_PRES_MEMdotAP1,r_PRES_MEMdotAP2,r_PRES_MEMdotAP3,r_PRES_MEMdotAP4)
            for (z=0;z<n_row;++z) {
                L2_r_PRES+=r_PRES[z]*r_PRES[z];
                r_PRESdotAP+=r_PRES[z]*AP[z];
                r_PRES_MEMdotAP0+=r_PRES_MEM[0][z]*AP[z];
                r_PRES_MEMdotAP1+=r_PRES_MEM[1][z]*AP[z];
                r_PRES_MEMdotAP2+=r_PRES_MEM[2][z]*AP[z];
                r_PRES_MEMdotAP3+=r_PRES_MEM[3][z]*AP[z];
                r_PRES_MEMdotAP4+=r_PRES_MEM[4][z]*AP[z];
             }
            #pragma omp master
            {
              r_PRES_MEMdotAP[0]=r_PRES_MEMdotAP0;
              r_PRES_MEMdotAP[1]=r_PRES_MEMdotAP1;
              r_PRES_MEMdotAP[2]=r_PRES_MEMdotAP2;
              r_PRES_MEMdotAP[3]=r_PRES_MEMdotAP3;
              r_PRES_MEMdotAP[4]=r_PRES_MEMdotAP4;
            }
         } else {
            #pragma omp master
            {
              L2_r_PRES=ZERO;
              r_PRESdotAP=ZERO;
              r_PRES_MEMdotAP0=ZERO;
            }
            #pragma omp for private(z) reduction(+:L2_r_PRES,r_PRESdotAP)
            for (z=0;z<n_row;++z) {
                L2_r_PRES+=r_PRES[z]*r_PRES[z];
                r_PRESdotAP+=r_PRES[z]*AP[z];
            }
            for (i=0;i<ORDER;++i) {
                #pragma omp for private(z) reduction(+:r_PRES_MEMdotAP0)
                for (z=0;z<n_row;++z) r_PRES_MEMdotAP0+=r_PRES_MEM[i][z]*AP[z];
                #pragma omp master
                {
                  r_PRES_MEMdotAP[i]=r_PRES_MEMdotAP0;
                  r_PRES_MEMdotAP0=ZERO;
                }
            }
         }
         /* if L2_r_PRES=0 there is a fatal breakdown */
         if (L2_r_PRES<=ZERO) {
            breakFlag=TRUE;
         } else {
            /***                                                                 
            ***** calculation of the weights for the update of X:               
            */                                                               
            #pragma omp master 
            {
               r_PRESdotAP*= -ONE/L2_r_PRES;
               for (i=0;i<ORDER;++i) r_PRES_MEMdotAP[i]*=-ONE/L2_r_PRES_MEM[i];
            }
            /*                                                                 
            ***** update of solution x_PRES and its residual r_PRES:            
            ***                                                               
            ***   P is used to accumulate X and AP to accumulate R. X and R     
            ***   are still needed until they are put into the X and R memory   
            ***   r_PRES_MEM and x_PRES_MEM                                     
            ***                                                                 
            **/
            if (ORDER==0) {
              #pragma omp for private(z)
              for (z=0;z<n_row;++z)
                  AP[z]+=r_PRESdotAP*r_PRES[z];
              #pragma omp for private(z)
              for (z=0;z<n_col;++z)
                  P[z]=-P[z]+r_PRESdotAP*x_PRES[z];
            } else if (ORDER==1) {
              #pragma omp for private(z)
              for (z=0;z<n_row;++z)
                 AP[z]+=r_PRESdotAP*r_PRES[z]+r_PRES_MEMdotAP[0]*r_PRES_MEM[0][z];
              #pragma omp for private(z)
              for (z=0;z<n_col;++z)
                P[z]=-P[z]+r_PRESdotAP*x_PRES[z]+r_PRES_MEMdotAP[0]*x_PRES_MEM[0][z];
            } else if (ORDER==2) {
              #pragma omp for private(z)
              for (z=0;z<n_row;++z)
                 AP[z]+=r_PRESdotAP*r_PRES[z]+r_PRES_MEMdotAP[0]*r_PRES_MEM[0][z]
                                             +r_PRES_MEMdotAP[1]*r_PRES_MEM[1][z];
              #pragma omp for private(z)
              for (z=0;z<n_col;++z)
                 P[z]=-P[z]+r_PRESdotAP*x_PRES[z]+r_PRES_MEMdotAP[0]*x_PRES_MEM[0][z]
                                            +r_PRES_MEMdotAP[1]*x_PRES_MEM[1][z];
            } else if (ORDER==3) {
              #pragma omp for private(z)
              for (z=0;z<n_row;++z)
                 AP[z]+=r_PRESdotAP*r_PRES[z]+r_PRES_MEMdotAP[0]*r_PRES_MEM[0][z]
                                             +r_PRES_MEMdotAP[1]*r_PRES_MEM[1][z]
                                             +r_PRES_MEMdotAP[2]*r_PRES_MEM[2][z];
              #pragma omp for private(z)
              for (z=0;z<n_col;++z)
                 P[z]=-P[z]+r_PRESdotAP*x_PRES[z]+r_PRES_MEMdotAP[0]*x_PRES_MEM[0][z]
                                            +r_PRES_MEMdotAP[1]*x_PRES_MEM[1][z]
                                            +r_PRES_MEMdotAP[2]*x_PRES_MEM[2][z];
            } else if (ORDER==4) {
              #pragma omp for private(z)
              for (z=0;z<n_row;++z)
                 AP[z]+=r_PRESdotAP*r_PRES[z]+r_PRES_MEMdotAP[0]*r_PRES_MEM[0][z]
                                             +r_PRES_MEMdotAP[1]*r_PRES_MEM[1][z]
                                             +r_PRES_MEMdotAP[2]*r_PRES_MEM[2][z]
                                             +r_PRES_MEMdotAP[3]*r_PRES_MEM[3][z];
              #pragma omp for private(z)
              for (z=0;z<n_col;++z)
                 P[z]=-P[z]+r_PRESdotAP*x_PRES[z]+r_PRES_MEMdotAP[0]*x_PRES_MEM[0][z]
                                            +r_PRES_MEMdotAP[1]*x_PRES_MEM[1][z]
                                            +r_PRES_MEMdotAP[2]*x_PRES_MEM[2][z]
                                            +r_PRES_MEMdotAP[3]*x_PRES_MEM[3][z];
            } else if (ORDER==5) {
              #pragma omp for private(z)
              for (z=0;z<n_row;++z)
                 AP[z]+=r_PRESdotAP*r_PRES[z]+r_PRES_MEMdotAP[0]*r_PRES_MEM[0][z]
                                             +r_PRES_MEMdotAP[1]*r_PRES_MEM[1][z]
                                             +r_PRES_MEMdotAP[2]*r_PRES_MEM[2][z]
                                             +r_PRES_MEMdotAP[3]*r_PRES_MEM[3][z]
                                             +r_PRES_MEMdotAP[4]*r_PRES_MEM[4][z];
              #pragma omp for private(z)
              for (z=0;z<n_col;++z)
                 P[z]=-P[z]+r_PRESdotAP*x_PRES[z]+r_PRES_MEMdotAP[0]*x_PRES_MEM[0][z]
                                            +r_PRES_MEMdotAP[1]*x_PRES_MEM[1][z]
                                            +r_PRES_MEMdotAP[2]*x_PRES_MEM[2][z]
                                            +r_PRES_MEMdotAP[3]*x_PRES_MEM[3][z]
                                            +r_PRES_MEMdotAP[4]*x_PRES_MEM[4][z];
            } else {

              #pragma omp for private(z)
              for (z=0;z<n_row;++z)
                AP[z]+=r_PRESdotAP*r_PRES[z];
                /* AP[z]+=r_PRESdotAP*r_PRES[z]; */
              #pragma omp for private(z)
              for (z=0;z<n_col;++z)
                P[z]=-P[z]+r_PRESdotAP*x_PRES[z];

              for (i=0;i<ORDER;++i) {
                #pragma omp for private(z)
                for (z=0;z<n_row;++z)
                    AP[z]+=r_PRES_MEMdotAP[i]*r_PRES_MEM[i][z];
                #pragma omp for private(z)
                for (z=0;z<n_col;++z)
                    P[z]+=r_PRES_MEMdotAP[i]*x_PRES_MEM[i][z];
              }
            }
            /***  factor scales AP and P to make it a residual and
            ***   as solution approximation
            ***
            ***   if factor is equal to zero a breakdown occurs. the
            ***   iteration procedure can be continued but r_PRES is not the
            ***   residual of x_PRES approximation.
            ***/
            factor=BREAKF*r_PRESdotAP;
            for (i=0;i<ORDER;++i) factor +=BREAKF_MEM[i]*r_PRES_MEMdotAP[i];
            /***
            ***** x_PRES and r_PRES are moved to memory:
            ***/
            #pragma omp master
            {
               L2_r_PRES_MEM[OLDEST]=L2_r_PRES;
               BREAKF_MEM[OLDEST]=BREAKF;
               tmp=x_PRES; x_PRES=x_PRES_MEM[OLDEST];x_PRES_MEM[OLDEST]=tmp;
               tmp=r_PRES; r_PRES=r_PRES_MEM[OLDEST];r_PRES_MEM[OLDEST]=tmp;
            }
            if (ABS(factor)<=ZERO) {
                /* in case of a break down: */
                BREAKF=ZERO;
                #pragma omp for private(z)
                for (z=0;z<n_row;++z) r_PRES[z]=AP[z];
                #pragma omp for private(z)
                for (z=0;z<n_col;++z) x_PRES[z]=P[z];
                /* is there any progress */
                breakFlag=TRUE;
                for (i=0;i<ORDER;++i) if (BREAKF_MEM[i]>ZERO) breakFlag=FALSE;
                convergeFlag = FALSE;
                maxIterFlag = FALSE;
            } else {
              BREAKF=ONE;
              breakFlag=FALSE;
              /***
              *** rescale P an AP and apply smooting:
              ***
              ***** calculate GAMMA from min_(GAMMA){|R+GAMMA*(r_PRES-R)|_2}:
              ***/
              #pragma omp master
              {
                 SC1=ZERO;
                 SC2=ZERO;
                 L2_R=ZERO;
              }
              factor=ONE/factor;
              #pragma omp for private(z,diff) reduction(+:SC1,SC2)
              for (z=0;z<n_row;++z) {
                r_PRES[z]=factor*AP[z];
                diff=r_PRES[z]-r[z];
                SC1+=diff*diff;
                SC2+=diff*r[z];
              }
              GAMMA=(SC1<=ZERO) ? ZERO : -SC2/SC1;
              #pragma omp for private(z)
              for (z=0;z<n_col;++z) {
                 x_PRES[z]=factor* P[z];
                 x[z]+=GAMMA*(x_PRES[z]-x[z]);
              }
              #pragma omp for private(z) reduction(+:L2_R)
              for (z=0;z<n_row;++z) {
                 r[z]+=GAMMA*(r_PRES[z]-r[z]);
                 L2_R+=r[z]*r[z];
              }
              norm_of_residual=sqrt(L2_R);
              convergeFlag = (norm_of_residual <= tol);
              maxIterFlag = (num_iter >= maxit);
            }
         } 
         if (restart>0) restartFlag=(num_iter_restart >= restart);
         if (! (convergeFlag || maxIterFlag || breakFlag))  goto next;
         /* end of iteration */
         #pragma omp master 
         {
           norm_of_residual_global=norm_of_residual;
	   num_iter_global=num_iter;
	   if (maxIterFlag) { 
	       status = SOLVER_MAXITER_REACHED;
	   } else if (breakFlag) {
	       status = SOLVER_BREAKDOWN;
	   }
        }
    }  /* end of parallel region */
    TMPMEMFREE(x_PRES);
    TMPMEMFREE(r_PRES);
    TMPMEMFREE(P);
    TMPMEMFREE(AP);
    for (i=0;i<length_of_recursion;i++) {
       TMPMEMFREE(x_PRES_MEM[i]);
       TMPMEMFREE(r_PRES_MEM[i]);
    }
    *iter=num_iter_global;
    *tolerance=norm_of_residual_global;
  }
  return status;
}
