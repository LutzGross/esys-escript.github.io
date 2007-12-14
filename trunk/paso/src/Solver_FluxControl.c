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

#define FLUX_L(a,b) SUPERBEE(a,b)  /* alter for other flux limiter */

#define FLUX_LIMITER(a) FLUX_L(a,1)

/**************************************************************/

/* free all memory used by FluxControl                                */

void Paso_Solver_FluxControl_free(Paso_Solver_FluxControl* in) {
     if (in!=NULL) {
        Paso_SystemMatrix_freeBuffer(in->matrix);
        Paso_SystemMatrix_free(in->matrix);
        MEMFREE(in->colorOf);
        MEMFREE(in->main_iptr);
        MEMFREE(in);
     }
}

/**************************************************************/

/*   constructs a flux control mechanism                      */

Paso_Solver_FluxControl* Paso_SolverFCT_getFluxControl(Paso_SystemMatrix * A) {

  Paso_Solver_FluxControl* out=NULL;
  dim_t n,i;
  index_t iptr,iptr_main,k;

  if (A==NULL) return out;
  n=Paso_SystemMatrix_getTotalNumRows(A);
  if (A->block_size!=1) {
        Paso_setError(TYPE_ERROR,"Paso_SolverFCT_getFluxControl: block size > 1 is not supported.");
        return NULL;
  }
  out=MEMALLOC(1,Paso_Solver_FluxControl);
  if (Paso_checkPtr(out)) return NULL;

  out->matrix=Paso_SystemMatrix_reference(A);
  out->colorOf=NULL;
  out->main_iptr=NULL;
  

  /* allocations: */  
  out->colorOf=MEMALLOC(n,index_t);
  out->main_iptr=MEMALLOC(n,index_t);
  if ( ! (Paso_checkPtr(out->colorOf) || Paso_checkPtr(out->main_iptr) ) ) {
      printf("Paso_SolverFCT_getFluxControl: Revise coloring!!\n");
      Paso_Pattern_color(A->mainBlock->pattern,&(out->num_colors),out->colorOf);
      Paso_SystemMatrix_allocBuffer(A);

      #pragma omp parallel for schedule(static) private(i,iptr,iptr_main,k)
      for (i = 0; i < n; ++i) {
        for (iptr=A->mainBlock->pattern->ptr[i];iptr<A->mainBlock->pattern->ptr[i+1]; ++iptr) {
             iptr_main=A->mainBlock->pattern->ptr[0]-1;
              for (iptr=A->mainBlock->pattern->ptr[i];iptr<A->mainBlock->pattern->ptr[i+1]; iptr++) {
                   if (A->mainBlock->pattern->index[iptr]==i) {
                        iptr_main=iptr;
                        break;
                   }
               }
               out->main_iptr[i]=iptr_main;
               if (iptr_main==A->mainBlock->pattern->ptr[0]-1)
                  Paso_setError(VALUE_ERROR, "Paso_SolverFCT_getFluxControl: no main diagonal");
           }
       }

  }
  if (Paso_noError()) {
     return out;
  } else {
     Paso_Solver_FluxControl_free(out);
     return NULL;
  }
} 

/**************************************************************/

/* adds A plus stabelising diffusion into the matrix B        */

/* d_ij=alpha*max(0,-a[i,j],-a[j,i])  */
/* b[i,j]+=alpha*(a[i,j]+d_ij)  */
/* b[j,i]+=alpha*(a[j,i]+d_ij)  */
/* b[i,i]-=alpha*d_ij  */
/* b[j,j]-=alpha*d_ij  */

void Paso_Solver_FluxControl_addDiffusion(Paso_Solver_FluxControl * fc, double alpha, Paso_SystemMatrix * B) {
  dim_t n,i;
  index_t color, iptr_ij,j,iptr_ji;
  register double d_ij;

  if (fc==NULL) return;
  n=Paso_SystemMatrix_getTotalNumRows(fc->matrix);
  /* TODO test - same pattern + block size */

  #pragma omp parallel private(color) 
  {
       /* process main block */
       for (color=0;color<fc->num_colors;++color) {
           #pragma omp for private(i,iptr_ij,j,iptr_ji,d_ij)  schedule(static)
           for (i = 0; i < n; ++i) {
               if (fc->colorOf[i]==color) {
                  for (iptr_ij=fc->matrix->mainBlock->pattern->ptr[i];iptr_ij<fc->matrix->mainBlock->pattern->ptr[i+1]; ++iptr_ij) {
                     j=fc->matrix->mainBlock->pattern->index[iptr_ij];
                     if (i<j) {
                        /* find entry a[j,i] */
                        for (iptr_ji=fc->matrix->mainBlock->pattern->ptr[i];iptr_ji<fc->matrix->mainBlock->pattern->ptr[j+1]-1; ++iptr_ji) {
                            if (fc->matrix->mainBlock->pattern->index[iptr_ji]==i) {
                                d_ij=(-alpha)*MIN3(0.,fc->matrix->mainBlock->val[iptr_ij],fc->matrix->mainBlock->val[iptr_ji]);
                                B->mainBlock->val[iptr_ij]+=alpha*fc->matrix->mainBlock->val[iptr_ij]+d_ij;
                                B->mainBlock->val[iptr_ji]+=alpha*fc->matrix->mainBlock->val[iptr_ji]+d_ij;
                                B->mainBlock->val[fc->main_iptr[i]]-=d_ij;
                                B->mainBlock->val[fc->main_iptr[j]]-=d_ij;
                                break;
                            }
                        }
                     }
                        
                  }
                  /* TODO process couple block */
               }
           }
           #pragma omp barrier
       }
  }

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

void Paso_Solver_FluxControl_setAntiDiffusiveFlux(Paso_Solver_FluxControl * fc, double * u, double* fa) {

  register double u_i,P_p,P_n,Q_p,Q_n,r_p,r_n, a_ij, d, u_j, r_ij, f_ij, a_ji;
  double *u_remote=NULL;
  index_t color, iptr_ij,j,iptr_ji, i;
  dim_t n;


  if (fc==NULL) return;
  n=Paso_SystemMatrix_getTotalNumRows(fc->matrix);
  /* exchange */
  Paso_SystemMatrix_startCollect(fc->matrix,u);
  u_remote=Paso_SystemMatrix_finishCollect(fc->matrix);

  #pragma omp parallel private(color) 
  {
       for (color=0;color<fc->num_colors;++color) {
           #pragma omp for schedule(static) private(i, u_i,P_p,P_n,Q_p,Q_n,r_p,r_n,iptr_ij,a_ij,d,j,iptr_ji, u_j, r_ij, f_ij, a_ji)
           for (i = 0; i < n; ++i) {
              if (fc->colorOf[i]==color) {
                  u_i=u[i];
                  /* gather the smoothness sensor */
                  P_p=0.;
                  P_n=0.;
                  Q_p=0.;
                  Q_n=0.;
                  #pragma ivdep
  	          for (iptr_ij=(fc->matrix->mainBlock->pattern->ptr[i]);iptr_ij<(fc->matrix->mainBlock->pattern->ptr[i+1]); ++iptr_ij) {
                      a_ij=fc->matrix->mainBlock->val[iptr_ij];
                      j=fc->matrix->mainBlock->pattern->index[iptr_ij];
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
  	          for (iptr_ij=(fc->matrix->coupleBlock->pattern->ptr[i]);iptr_ij<(fc->matrix->coupleBlock->pattern->ptr[i+1]); ++iptr_ij) {
                      a_ij=fc->matrix->coupleBlock->val[iptr_ij];
                      j=fc->matrix->coupleBlock->pattern->index[iptr_ij];
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
                  r_p = (P_p > 0.) ? FLUX_LIMITER(Q_p/P_p) : 0.;
                  r_n = (P_n < 0 ) ? FLUX_LIMITER(Q_n/P_n) : 0.;
                  /* anti diffusive flux from main block */
                  for (iptr_ij=fc->matrix->mainBlock->pattern->ptr[i];iptr_ij<fc->matrix->mainBlock->pattern->ptr[i+1]; ++iptr_ij) {
                     a_ij=fc->matrix->mainBlock->val[iptr_ij];
                     j=fc->matrix->mainBlock->pattern->index[iptr_ij];
                     if (a_ij < 0 && i!=j) {
                        /* find entry a[j,i] */
                        for (iptr_ji=fc->matrix->mainBlock->pattern->ptr[i];iptr_ji<fc->matrix->mainBlock->pattern->ptr[j+1]-1; ++iptr_ji) {
                            if (fc->matrix->mainBlock->pattern->index[iptr_ji]==i) {
                                a_ji=fc->matrix->mainBlock->val[iptr_ji];
                                if  (a_ji > a_ij || (a_ji == a_ij && j<i) ) {
                                    u_j=u[j];
                                    r_ij = u_i>u_j ? r_p : r_n;
                                    f_ij =MIN(r_ij*a_ij,a_ji-a_ij)*(u_i-u_j);
                                    fa[i]+=f_ij;
                                    fa[j]-=f_ij;
                                    break;
                                }
                            }
                        }
                     }
                  }
                  /* anti diffusive flux from couple block */

                  /* TODO */
              }
           }
           #pragma omp barrier
       }
  }
}
