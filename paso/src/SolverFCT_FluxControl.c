/* $Id:$ */

/*******************************************************
 *
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/* Paso: FluxControl                                          */

/**************************************************************/

/* Author: l.gross@uq.edu.au                                */

/**************************************************************/


#include "Paso.h"
#include "SolverFCT.h"
#include "PasoUtil.h"

#define FLUX_S(a,b) ((SIGN(a)+SIGN(b))/2.)
#define MINMOD(a,b) (FLUX_S(a,b)*MIN(ABS(a),ABS(b)))
#define SUPERBEE(a,b) (FLUX_S(a,b)*MAX(MIN(2*ABS(a),ABS(b)),MIN(ABS(a),2*ABS(b))))
#define MC(a,b) (FLUX_S(a,b)*MIN3(ABS((a)+(b))/2,2*ABS(a),2*ABS(b)))

#define FLUX_L(a,b) MC(a,b)  /* alter for other flux limiter */

#define FLUX_LIMITER(a) FLUX_L(1,a)

/**************************************************************/

/* adds A plus stabelising diffusion into the matrix B        */

/* d_ij=alpha*max(0,-a[i,j],-a[j,i])  */
/* b[i,j]+=alpha*(a[i,j]+d_ij)  */
/* b[j,i]+=alpha*(a[j,i]+d_ij)  */
/* b[i,i]-=alpha*d_ij  */
/* b[j,j]-=alpha*d_ij  */

void Paso_FCTransportProblem_addAdvectivePart(Paso_FCTransportProblem * fc, double alpha) {
  dim_t n,i;
  index_t color, iptr_ij,j,iptr_ji;
  register double d_ij, sum;

  if (fc==NULL) return;
  n=Paso_SystemMatrix_getTotalNumRows(fc->flux_matrix);

  #pragma omp parallel for private(i,iptr_ij,j,iptr_ji,d_ij,sum)  schedule(static)
  for (i = 0; i < n; ++i) {
     sum=alpha*fc->flux_matrix->mainBlock->val[fc->main_iptr[i]];
     for (iptr_ij=fc->flux_matrix->mainBlock->pattern->ptr[i];iptr_ij<fc->flux_matrix->mainBlock->pattern->ptr[i+1]; ++iptr_ij) {
         j=fc->flux_matrix->mainBlock->pattern->index[iptr_ij];
         if (j!=i) {
             /* find entry a[j,i] */
             for (iptr_ji=fc->flux_matrix->mainBlock->pattern->ptr[j]; iptr_ji<fc->flux_matrix->mainBlock->pattern->ptr[j+1]; ++iptr_ji) {
                if (fc->flux_matrix->mainBlock->pattern->index[iptr_ji]==i) {
                    d_ij=(-alpha)*MIN3(0.,fc->flux_matrix->mainBlock->val[iptr_ij],
                                          fc->flux_matrix->mainBlock->val[iptr_ji]);
                    fc->transport_matrix->mainBlock->val[iptr_ij]+=
                                         alpha*fc->flux_matrix->mainBlock->val[iptr_ij]+d_ij;
                    sum-=d_ij;
                    break;
                }
             }
         }
     }
     /* TODO process couple block */

     /* update main diagonal */
     fc->transport_matrix->mainBlock->val[fc->main_iptr[i]]+=sum;
  }

}

void Paso_FCTransportProblem_setFlux(Paso_FCTransportProblem * fc, double * u, double* fa) {
  /*
   *   sets fa=transport_matrix*u+anti-diffuison_flux(u)
   */

  double *remote_u=NULL;
  
  if (fc==NULL) return;
  Paso_SystemMatrix_startCollect(fc->transport_matrix,u);
  /* process main block */
  Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(1.,fc->transport_matrix->mainBlock,u,0.,fa);
  /* finish exchange */
  remote_u=Paso_SystemMatrix_finishCollect(fc->transport_matrix);
  /* process couple block */
  Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(1.,fc->transport_matrix->coupleBlock,remote_u,1.,fa);

  Paso_FCTransportProblem_setAntiDiffusiveFlux(fc,u,remote_u,fa); 
}
/**************************************************************/

/* adds antidiffusion to fa  
 
   P_p[i] + = sum_j min(0,a[i,j]) min(0,u[j]-u[i])   
   P_n[i] + = sum_j min(0,a[i,j]) max(0,u[j]-u[i])  
   Q_p[i] + = sum_j max(0,a[i,j]) max(0,u[j]-u[i])  
   Q_n[i] + = sum_j max(0,a[i,j]) min(0,u[j]-u[i])   
   d_ij=max(0,-a[i,j],-a[j,i])
   l_ji=max(0,a[j,i],a[j,i]-a[i,j])
   if a[j,i] >= a[i,j] and 0>a[i,j] : (i.e d_ij>0 and l_ji>=l_ij)
      r_ij = u[i]>u[j] ? Q_p[i]/P_p[i] : Q_n[i]/Q_n[i]
      f_ij =min(FLUX_LIMITER(r_ij)*d_ij,l_ji) (u[i]-u[j])=min(FLUX_LIMITER(r_ij)*a[i,j],a[j,i]-a[i,j]) (u[i]-u[j])
      fa[i]+=f_ij
      fa[j]-=f_ij

*/

void Paso_FCTransportProblem_setAntiDiffusiveFlux(Paso_FCTransportProblem * fc, double * u, double *u_remote, double* fa) {

  register double u_i,P_p,P_n,Q_p,Q_n,r_p,r_n, a_ij, d, u_j, r_ij, f_ij, a_ji, d_ij, sum;
  index_t color, iptr_ij,j,iptr_ji, i;
  dim_t n;


  if (fc==NULL) return;
  n=Paso_SystemMatrix_getTotalNumRows(fc->flux_matrix);


  #pragma omp parallel
  {
      /*
       * calculate the smootness sensors 
      */
      #pragma omp for schedule(static) private(i, u_i,P_p,P_n,Q_p,Q_n,iptr_ij,a_ij,j,d)
      for (i = 0; i < n; ++i) {
          u_i=u[i];
          P_p=0.;
          P_n=0.;
          Q_p=0.;
          Q_n=0.;
          #pragma ivdep
  	  for (iptr_ij=(fc->flux_matrix->mainBlock->pattern->ptr[i]);iptr_ij<(fc->flux_matrix->mainBlock->pattern->ptr[i+1]); ++iptr_ij) {
               a_ij=fc->flux_matrix->mainBlock->val[iptr_ij];
               j=fc->flux_matrix->mainBlock->pattern->index[iptr_ij];
               d=u[j]-u_i;
               if (a_ij<0.) {
                  if (d<0.) {
                      P_p+=a_ij*d;
                  } else {
                      P_n+=a_ij*d;
                  }
               } else {
                  if (d>0.) {
                    Q_p+=a_ij*d;
                  } else {
                    Q_n+=a_ij*d;
                  }
               }
  	  }
          #pragma ivdep
  	  for (iptr_ij=(fc->flux_matrix->coupleBlock->pattern->ptr[i]);iptr_ij<(fc->flux_matrix->coupleBlock->pattern->ptr[i+1]); ++iptr_ij) {
               a_ij=fc->flux_matrix->coupleBlock->val[iptr_ij];
               j=fc->flux_matrix->coupleBlock->pattern->index[iptr_ij];
               d=u_remote[j]-u_i;
               if (a_ij<0.) {
                   if (d<0.) {
                     P_p+=a_ij*d;
                   } else {
                     P_n+=a_ij*d;
                   }
               } else {
                  if (d>0.) {
                     Q_p+=a_ij*d;
                  } else {
                     Q_n+=a_ij*d;
                  }
               }
  	  }
          /* set the smoothness indicators */
          fc->r_p[i] = (P_p > 0.) ? FLUX_LIMITER(Q_p/P_p) : 0.;
          fc->r_n[i] = (P_n < 0.) ? FLUX_LIMITER(Q_n/P_n) : 0.;

      } /* end of row loop for smootheness indicator */

      /*
       * calculate antidiffusion
      */
      #pragma omp for schedule(static) private(i, u_i, sum, iptr_ij, a_ij, j, iptr_ji, a_ji,d_ij, u_j, r_ij, f_ij)
      for (i = 0; i < n; ++i) {
          u_i=u[i];
          sum=0;
          /* anti diffusive flux from main block */
          for (iptr_ij=fc->flux_matrix->mainBlock->pattern->ptr[i];iptr_ij<fc->flux_matrix->mainBlock->pattern->ptr[i+1]; ++iptr_ij) {
              a_ij=fc->flux_matrix->mainBlock->val[iptr_ij];
              j=fc->flux_matrix->mainBlock->pattern->index[iptr_ij];
              if ( i!=j ) {
                   /* find entry a[j,i] */
                   for (iptr_ji=fc->flux_matrix->mainBlock->pattern->ptr[j];iptr_ji<fc->flux_matrix->mainBlock->pattern->ptr[j+1]; ++iptr_ji) {
                      if (fc->flux_matrix->mainBlock->pattern->index[iptr_ji]==i) {
                           a_ji=fc->flux_matrix->mainBlock->val[iptr_ji];
                           d_ij=-MIN3(0,a_ji,a_ij);
                           if (d_ij > 0.) {
                                u_j=u[j];
                                if (a_ji >= a_ij) {
                                    r_ij = u_i>u_j ? fc->r_p[i] : fc->r_n[i];
                                    f_ij =MIN(r_ij*d_ij,a_ji+d_ij);
                                } else {
                                    r_ij = u_j>u_i ? fc->r_p[j] : fc->r_n[j];
                                    f_ij =MIN(r_ij*d_ij,a_ij+d_ij);
                                }
                                sum+=f_ij*(u_i-u_j);
                           }
                           break;
                      }
                   }
              }
          }
          /* anti diffusive flux from couple block */
          /* TODO */



          fa[i]+=sum;
      }
  } /* end of parallel block */
}
