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


/**************************************************************/

/* create the low order transport matrix and stores it negative values 
 * into the iteration_matrix accept the main diagonal which is stored seperately
 * if fc->iteration_matrix==NULL, fc->iteration_matrix is allocated 
 *
 * a=transport_matrix 
 * b= low_order_transport_matrix = - iteration_matrix
 * c=main diagonal low_order_transport_matrix 
 * initialize c[i] mit a[i,i] 
 *    d_ij=max(0,-a[i,j],-a[j,i])  
 *    b[i,j]=-(a[i,j]+d_ij)         
 *    c[i]-=d_ij                  
 */

void Paso_FCTransportProblem_setLowOrderOperator(Paso_FCTransportProblem * fc) {
  dim_t n=Paso_SystemMatrix_getTotalNumRows(fc->transport_matrix),i;
  index_t color, iptr_ij,j,iptr_ji;
  Paso_SystemMatrixPattern *pattern;
  register double d_ij, sum, rtmp1, rtmp2;

  if (fc==NULL) return;
  if (fc->iteration_matrix==NULL) {
      fc->iteration_matrix=Paso_SystemMatrix_alloc(fc->transport_matrix->type,
                                                   fc->transport_matrix->pattern,
                                                   fc->transport_matrix->row_block_size,
                                                   fc->transport_matrix->col_block_size);
  }

  if (Paso_noError()) {
      pattern=fc->iteration_matrix->pattern;
      n=Paso_SystemMatrix_getTotalNumRows(fc->iteration_matrix);
      #pragma omp parallel for private(i,sum,iptr_ij,j,iptr_ji,rtmp1, rtmp2,d_ij)  schedule(static)
      for (i = 0; i < n; ++i) {
          sum=fc->transport_matrix->mainBlock->val[fc->main_iptr[i]];
          /* look at a[i,j] */
          for (iptr_ij=pattern->mainPattern->ptr[i];iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
              j=pattern->mainPattern->index[iptr_ij];
              if (j!=i) {
                 /* find entry a[j,i] */
                 for (iptr_ji=pattern->mainPattern->ptr[j]; iptr_ji<pattern->mainPattern->ptr[j+1]; ++iptr_ji) {
                    if (pattern->mainPattern->index[iptr_ji]==i) {
                        rtmp1=fc->transport_matrix->mainBlock->val[iptr_ij];
                        rtmp2=fc->transport_matrix->mainBlock->val[iptr_ji];
                        d_ij=-MIN3(0.,rtmp1,rtmp2);
                        fc->iteration_matrix->mainBlock->val[iptr_ij]=-(rtmp1+d_ij);
                        sum-=d_ij;
                        break;
                    }
                 }
             }
         }
         /* TODO process couple block */
        
         /* set main diagonal entry */
         fc->main_diagonal_low_order_transport_matrix[i]=sum;
      }
  }
}
/*
 * out_i=m_i u_i + alpha \sum_{j <> i} l_{ij} (u_j-u_i) + beta q_i
 *
 */
void Paso_SolverFCT_setMuPaLuPbQ(double* out,
                                 const double* M, 
                                 const double* u,
                                 const double a,
                                 Paso_SystemMatrix *L,
                                 const double b,
                                 const double* Q) 
{
  dim_t n, i, j;
  Paso_SystemMatrixPattern *pattern;
  double *remote_u;
  register double sum, u_i, l_ij;
  register index_t iptr_ij;

  n=Paso_SystemMatrix_getTotalNumRows(L);

  if (ABS(a)>0) {
      Paso_SystemMatrix_startCollect(L,u);
  }
  #pragma omp parallel for private(i) schedule(static)
  for (i = 0; i < n; ++i) {
      out[i]=M[i]*u[i]+b*Q[i];
  } 
  if (ABS(a)>0) {
      pattern=L->pattern;
      remote_u=Paso_SystemMatrix_finishCollect(L);
      #pragma omp parallel for schedule(static) private(i, sum, u_i, iptr_ij, j, l_ij)
      for (i = 0; i < n; ++i) {
          sum=0;
          u_i=u[i];
          #pragma ivdep
  	  for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
               j=pattern->mainPattern->index[iptr_ij];
               l_ij=L->mainBlock->val[iptr_ij];
               sum+=l_ij*(u[j]-u_i);
          }
          #pragma ivdep
  	  for (iptr_ij=(pattern->couplePattern->ptr[i]);iptr_ij<pattern->couplePattern->ptr[i+1]; ++iptr_ij) {
               j=pattern->couplePattern->index[iptr_ij];
               l_ij=L->coupleBlock->val[iptr_ij];
               sum+=l_ij*(remote_u[j]-u_i);
          }
          out[i]+=a*sum;
      }
  }
}
/* 
 * calculate QP[i] max_{i in L->pattern[i]} (u[j]-u[i] ) 
 *           QN[i] min_{i in L->pattern[i]} (u[j]-u[i] ) 
 *
*/
void Paso_SolverFCT_setQs(const double* u,double* QN, double* QP, Paso_SystemMatrix *L) 
{
  dim_t n, i, j;
  Paso_SystemMatrixPattern *pattern;
  double *remote_u;
  register double u_i, u_min_i, u_max_i, u_j;
  register index_t iptr_ij;

  Paso_SystemMatrix_startCollect(L,u);
  n=Paso_SystemMatrix_getTotalNumRows(L);
  pattern=L->pattern;
  remote_u=Paso_SystemMatrix_finishCollect(L);
  #pragma omp parallel for schedule(static) private(i, u_i, u_min_i, u_max_i, j, u_j, iptr_ij)
  for (i = 0; i < n; ++i) {
     u_i=u[i];
     u_min_i=u_i;
     u_max_i=u_i;
     #pragma ivdep
     for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
         j=pattern->mainPattern->index[iptr_ij];
         u_j=u[j];
         u_min_i=MIN(u_min_i,u_j);
         u_max_i=MAX(u_max_i,u_j);
     }
     #pragma ivdep
     for (iptr_ij=(pattern->couplePattern->ptr[i]);iptr_ij<pattern->couplePattern->ptr[i+1]; ++iptr_ij) {
          j=pattern->couplePattern->index[iptr_ij];
          u_j=remote_u[j];
          u_min_i=MIN(u_min_i,u_j);
          u_max_i=MAX(u_max_i,u_j);
      }
      QN[i]=u_min_i-u_i;
      QP[i]=u_max_i-u_i;
  }
}

/*
 *
 *  f_{ij} = (m_{ij} - dt (1-theta) d_{ij}) (u_last[j]-u_last[i]) - (m_{ij} + dt theta d_{ij}) (u[j]-u[i])
 *         
 * m=fc->mass matrix
 * d=artifical diffusion matrix = L - K = - fc->iteration matrix - fc->transport matrix (away from main diagonal)
 */
void Paso_FCTransportProblem_setAntiDiffusionFlux(const double dt, const Paso_FCTransportProblem * fc, Paso_SystemMatrix *flux_matrix, const double* u, const double* u_last)
{
  dim_t n, j, i;
  double *remote_u=NULL, *remote_u_last=NULL;
  index_t iptr_ij;
  const double f1=- dt * (1.-fc->theta);
  const double f2=  dt * fc->theta;
  register double m_ij, d_ij, u_i, u_last_i, d_u_last, d_u;
  Paso_SystemMatrixPattern *pattern;
  Paso_SystemMatrix_startCollect(fc->iteration_matrix,u);
  Paso_SystemMatrix_startCollect(fc->iteration_matrix,u_last); 
  n=Paso_SystemMatrix_getTotalNumRows(fc->iteration_matrix);
  pattern=fc->iteration_matrix->pattern;
  remote_u=Paso_SystemMatrix_finishCollect(fc->iteration_matrix);
  remote_u_last=Paso_SystemMatrix_finishCollect(fc->iteration_matrix); 
  if ( (ABS(f1) >0 ) ) {
     if ( (ABS(f2) >0 ) ) {
        #pragma omp parallel for schedule(static) private(i, u_i, u_last_i, iptr_ij, j,m_ij,d_ij, d_u_last, d_u)
        for (i = 0; i < n; ++i) {
           u_i=u[i];
           u_last_i=u_last[i];
           #pragma ivdep
           for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
              j=pattern->mainPattern->index[iptr_ij];
              m_ij=fc->mass_matrix->mainBlock->val[iptr_ij];
              d_ij=-(fc->transport_matrix->mainBlock->val[iptr_ij]+
                                                fc->iteration_matrix->mainBlock->val[iptr_ij]);
              d_u=u[j]-u_i;
              d_u_last=u_last[j]-u_last_i;
              flux_matrix->mainBlock->val[iptr_ij]=(m_ij+f1*d_ij)*d_u_last- (m_ij+f2*d_ij)*d_u;
           }
           #pragma ivdep
           for (iptr_ij=(pattern->couplePattern->ptr[i]);iptr_ij<pattern->couplePattern->ptr[i+1]; ++iptr_ij) {
              j=pattern->couplePattern->index[iptr_ij];
              m_ij=fc->mass_matrix->coupleBlock->val[iptr_ij];
              d_ij=-(fc->transport_matrix->coupleBlock->val[iptr_ij]+
                                             fc->iteration_matrix->coupleBlock->val[iptr_ij]);
              d_u=remote_u[j]-u_i;
              d_u_last=remote_u_last[j]-u_last_i;
              flux_matrix->coupleBlock->val[iptr_ij]=(m_ij+f1*d_ij)*d_u_last- (m_ij+f2*d_ij)*d_u;
           }
        }
     } else {
        #pragma omp parallel for schedule(static) private(i, u_i, u_last_i, iptr_ij, j,m_ij,d_ij, d_u_last, d_u)
        for (i = 0; i < n; ++i) {
           u_i=u[i];
           u_last_i=u_last[i];
           #pragma ivdep
           for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
              j=pattern->mainPattern->index[iptr_ij];
              m_ij=fc->mass_matrix->mainBlock->val[iptr_ij];
              d_ij=-(fc->transport_matrix->mainBlock->val[iptr_ij]+
                                                fc->iteration_matrix->mainBlock->val[iptr_ij]);
              d_u=u[j]-u_i;
              d_u_last=u_last[j]-u_last_i;
              flux_matrix->mainBlock->val[iptr_ij]=(m_ij+f1*d_ij)*d_u_last-m_ij*d_u;
           }
           #pragma ivdep
           for (iptr_ij=(pattern->couplePattern->ptr[i]);iptr_ij<pattern->couplePattern->ptr[i+1]; ++iptr_ij) {
              j=pattern->couplePattern->index[iptr_ij];
              m_ij=fc->mass_matrix->coupleBlock->val[iptr_ij];
              d_ij=-(fc->transport_matrix->coupleBlock->val[iptr_ij]+
                                             fc->iteration_matrix->coupleBlock->val[iptr_ij]);
              d_u=remote_u[j]-u_i;
              d_u_last=remote_u_last[j]-u_last_i;
              flux_matrix->coupleBlock->val[iptr_ij]=(m_ij+f1*d_ij)*d_u_last-m_ij*d_u;
           }
        }
     }
  } else {
     if ( (ABS(f2) >0 ) ) {
        #pragma omp parallel for schedule(static) private(i, u_i, u_last_i, iptr_ij, j,m_ij,d_ij, d_u_last, d_u)
        for (i = 0; i < n; ++i) {
           u_i=u[i];
           u_last_i=u_last[i];
           #pragma ivdep
           for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
              j=pattern->mainPattern->index[iptr_ij];
              m_ij=fc->mass_matrix->mainBlock->val[iptr_ij];
              d_ij=-(fc->transport_matrix->mainBlock->val[iptr_ij]+
                                                fc->iteration_matrix->mainBlock->val[iptr_ij]);
              d_u=u[j]-u_i;
              d_u_last=u_last[j]-u_last_i;
              flux_matrix->mainBlock->val[iptr_ij]=m_ij*d_u_last- (m_ij+f2*d_ij)*d_u;
           }
           #pragma ivdep
           for (iptr_ij=(pattern->couplePattern->ptr[i]);iptr_ij<pattern->couplePattern->ptr[i+1]; ++iptr_ij) {
              j=pattern->couplePattern->index[iptr_ij];
              m_ij=fc->mass_matrix->coupleBlock->val[iptr_ij];
              d_ij=-(fc->transport_matrix->coupleBlock->val[iptr_ij]+
                                             fc->iteration_matrix->coupleBlock->val[iptr_ij]);
              d_u=remote_u[j]-u_i;
              d_u_last=remote_u_last[j]-u_last_i;
              flux_matrix->coupleBlock->val[iptr_ij]=m_ij*d_u_last- (m_ij+f2*d_ij)*d_u;
           }
        }
     } else {
        #pragma omp parallel for schedule(static) private(i, u_i, u_last_i, iptr_ij, j,m_ij,d_ij, d_u_last, d_u)
        for (i = 0; i < n; ++i) {
           u_i=u[i];
           u_last_i=u_last[i];
           #pragma ivdep
           for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
              j=pattern->mainPattern->index[iptr_ij];
              m_ij=fc->mass_matrix->mainBlock->val[iptr_ij];
              d_ij=-(fc->transport_matrix->mainBlock->val[iptr_ij]+
                                                fc->iteration_matrix->mainBlock->val[iptr_ij]);
              d_u=u[j]-u_i;
              d_u_last=u_last[j]-u_last_i;
              flux_matrix->mainBlock->val[iptr_ij]=m_ij*(d_u_last-d_u);
           }
           #pragma ivdep
           for (iptr_ij=(pattern->couplePattern->ptr[i]);iptr_ij<pattern->couplePattern->ptr[i+1]; ++iptr_ij) {
              j=pattern->couplePattern->index[iptr_ij];
              m_ij=fc->mass_matrix->coupleBlock->val[iptr_ij];
              d_ij=-(fc->transport_matrix->coupleBlock->val[iptr_ij]+
                                             fc->iteration_matrix->coupleBlock->val[iptr_ij]);
              d_u=remote_u[j]-u_i;
              d_u_last=remote_u_last[j]-u_last_i;
              flux_matrix->coupleBlock->val[iptr_ij]=m_ij*(d_u_last-d_u);
           }
        }
      }
  }
}
/*
 *
 * f_{ij} + = (a*m_{ij} + b* d_{ij}) (u[j]-u[i])
 *
 * m=fc->mass matrix
 * d=artifical diffusion matrix = L - K = - fc->iteration matrix - fc->transport matrix (away from main diagonal)
 */
void Paso_FCTransportProblem_updateAntiDiffusionFlux(const Paso_FCTransportProblem * fc, Paso_SystemMatrix *flux_matrix,const double a, const double b, const double* u)
{
  dim_t n, j, i;
  double *remote_u;
  index_t iptr_ij;
  const double b2=-b;
  register double m_ij, ml_ij, k_ij, u_i;
  Paso_SystemMatrixPattern *pattern;
  Paso_SystemMatrix_startCollect(fc->iteration_matrix,u);
  n=Paso_SystemMatrix_getTotalNumRows(fc->iteration_matrix);
  pattern=fc->iteration_matrix->pattern;
  remote_u=Paso_SystemMatrix_finishCollect(fc->iteration_matrix);
  if ( (ABS(a) >0 ) ) {
      if ( (ABS(b) >0 ) ) {
         #pragma omp parallel for schedule(static) private(i, u_i, iptr_ij, j,m_ij,k_ij,ml_ij)
         for (i = 0; i < n; ++i) {
            u_i=u[i];
            #pragma ivdep
            for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
                j=pattern->mainPattern->index[iptr_ij];
                m_ij=fc->mass_matrix->mainBlock->val[iptr_ij];
                k_ij=fc->transport_matrix->mainBlock->val[iptr_ij];
                ml_ij=fc->iteration_matrix->mainBlock->val[iptr_ij];
                flux_matrix->mainBlock->val[iptr_ij]+=(u[j]-u_i)*(a*m_ij+b2*(ml_ij+k_ij));
            }
            #pragma ivdep
            for (iptr_ij=(pattern->couplePattern->ptr[i]);iptr_ij<pattern->couplePattern->ptr[i+1]; ++iptr_ij) {
                j=pattern->couplePattern->index[iptr_ij];
                m_ij=fc->mass_matrix->coupleBlock->val[iptr_ij];
                k_ij=fc->transport_matrix->coupleBlock->val[iptr_ij];
                ml_ij=fc->iteration_matrix->coupleBlock->val[iptr_ij];
                flux_matrix->coupleBlock->val[iptr_ij]+=(remote_u[j]-u_i)*(a*m_ij+b2*(ml_ij+k_ij));
            }
         }
      } else {
         #pragma omp parallel for schedule(static) private(i, u_i, iptr_ij, j,m_ij)
         for (i = 0; i < n; ++i) {
            u_i=u[i];
            #pragma ivdep
            for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
                j=pattern->mainPattern->index[iptr_ij];
                m_ij=fc->mass_matrix->mainBlock->val[iptr_ij];
                flux_matrix->mainBlock->val[iptr_ij]+=(u[j]-u_i)*a*m_ij;
            }
            #pragma ivdep
            for (iptr_ij=(pattern->couplePattern->ptr[i]);iptr_ij<pattern->couplePattern->ptr[i+1]; ++iptr_ij) {
                j=pattern->couplePattern->index[iptr_ij];
                m_ij=fc->mass_matrix->coupleBlock->val[iptr_ij];
                flux_matrix->coupleBlock->val[iptr_ij]+=(remote_u[j]-u_i)*a*m_ij;
            }
         }


      }
  } else {
      if ( (ABS(b) >0 ) ) {
         #pragma omp parallel for schedule(static) private(i, u_i, iptr_ij, j,k_ij, ml_ij)
         for (i = 0; i < n; ++i) {
            u_i=u[i];
            #pragma ivdep
            for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
                j=pattern->mainPattern->index[iptr_ij];
                k_ij=fc->transport_matrix->mainBlock->val[iptr_ij];
                ml_ij=fc->iteration_matrix->mainBlock->val[iptr_ij];
                flux_matrix->mainBlock->val[iptr_ij]+=(u[j]-u_i)*b2*(ml_ij+k_ij);
            }
            #pragma ivdep
            for (iptr_ij=(pattern->couplePattern->ptr[i]);iptr_ij<pattern->couplePattern->ptr[i+1]; ++iptr_ij) {
                j=pattern->couplePattern->index[iptr_ij];
                k_ij=fc->transport_matrix->coupleBlock->val[iptr_ij];
                ml_ij=fc->iteration_matrix->coupleBlock->val[iptr_ij];
                flux_matrix->coupleBlock->val[iptr_ij]+=(remote_u[j]-u_i)*b2*(ml_ij+k_ij);
            }
         }

      } 
  }
}
/*
 *  apply pre flux-correction: f_{ij}:=0 if f_{ij}*(u_[i]-u[j])<=0
 *  
 */
void Paso_FCTransportProblem_applyPreAntiDiffusionCorrection(Paso_SystemMatrix *f,const double* u)
{
  dim_t n, i, j;
  Paso_SystemMatrixPattern *pattern;
  double *remote_u;
  register double u_i, f_ij;
  register index_t iptr_ij;

  n=Paso_SystemMatrix_getTotalNumRows(f);
  pattern=f->pattern;
  remote_u=Paso_SystemMatrix_finishCollect(f);
  #pragma omp parallel for schedule(static) private(i, u_i, iptr_ij, j, f_ij)
  for (i = 0; i < n; ++i) {
      u_i=u[i];
      #pragma ivdep
      for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
          j=pattern->mainPattern->index[iptr_ij];
          f_ij=f->mainBlock->val[iptr_ij];
          if (f_ij * (u_i-u[j]) <= 0) f->mainBlock->val[iptr_ij]=0;
      }
      #pragma ivdep
      for (iptr_ij=(pattern->couplePattern->ptr[i]);iptr_ij<pattern->couplePattern->ptr[i+1]; ++iptr_ij) {
          j=pattern->couplePattern->index[iptr_ij];
          f_ij=f->coupleBlock->val[iptr_ij];
          if (f_ij * (u_i-remote_u[j]) <= 0) f->coupleBlock->val[iptr_ij]=0;
      }
  }
}


void Paso_FCTransportProblem_setRs(const Paso_SystemMatrix *f,const double* lumped_mass_matrix,
                                   const double* QN,const double* QP,double* RN,double* RP)
{
  dim_t n, i, j;
  Paso_SystemMatrixPattern *pattern;
  register double f_ij, PP_i, PN_i;
  register index_t iptr_ij;

  n=Paso_SystemMatrix_getTotalNumRows(f);
  pattern=f->pattern;
  #pragma omp parallel for schedule(static) private(i, iptr_ij, j, f_ij, PP_i, PN_i)
  for (i = 0; i < n; ++i) {
      PP_i=0.;
      PN_i=0.;
      #pragma ivdep
      for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
          j=pattern->mainPattern->index[iptr_ij];
          if (i != j ) {
             f_ij=f->mainBlock->val[iptr_ij];
             if (f_ij <=0) {
                PN_i+=f_ij;
             } else {
                PP_i+=f_ij;
             }
          }
      }
      #pragma ivdep
      for (iptr_ij=(pattern->couplePattern->ptr[i]);iptr_ij<pattern->couplePattern->ptr[i+1]; ++iptr_ij) {
          j=pattern->couplePattern->index[iptr_ij];
          f_ij=f->coupleBlock->val[iptr_ij];
          if (f_ij <=0  ) {
                PN_i+=f_ij;
          } else {
                PP_i+=f_ij;
          }
      }
      if (PN_i<0) {
         RN[i]=MIN(1,QN[i]*lumped_mass_matrix[i]/PN_i);
      } else {
         RN[i]=1.;
      }
      if (PP_i>0) {
         RP[i]=MIN(1,QP[i]*lumped_mass_matrix[i]/PP_i);
      } else {
         RP[i]=1.;
      }
  }

}

void Paso_FCTransportProblem_addCorrectedFluxes(double* f,Paso_SystemMatrix *flux_matrix,const double* RN,const double* RP) 
{
  dim_t n, i, j;
  Paso_SystemMatrixPattern *pattern;
  double *remote_RN, *remote_RP;
  register double RN_i, RP_i, f_i, f_ij;
  register index_t iptr_ij;
  n=Paso_SystemMatrix_getTotalNumRows(flux_matrix);

  Paso_SystemMatrix_startCollect(flux_matrix,RN);
  pattern=flux_matrix->pattern;
  remote_RN=Paso_SystemMatrix_finishCollect(flux_matrix);
  #pragma omp parallel for schedule(static) private(i, RN_i, RP_i, iptr_ij, j, f_ij, f_i)
  for (i = 0; i < n; ++i) {
     RN_i=RN[i];
     RP_i=RP[i];
     f_i=0;
     #pragma ivdep
     for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
         j=pattern->mainPattern->index[iptr_ij];
         f_ij=flux_matrix->mainBlock->val[iptr_ij];
         if (f_ij >=0) {
              f_i+=f_ij*MIN(RP_i,RN[j]);
         } else {
              f_i+=f_ij*MIN(RN_i,RP[j]);
         }
     }
     #pragma ivdep
     for (iptr_ij=(pattern->couplePattern->ptr[i]);iptr_ij<pattern->couplePattern->ptr[i+1]; ++iptr_ij) {
          j=pattern->couplePattern->index[iptr_ij];
          f_ij=flux_matrix->coupleBlock->val[iptr_ij];
          if (f_ij >=0) {
              f_i+=f_ij*MIN(RP_i,remote_RN[j]);
          }
      }
      f[i]+=f_i;
  }
  Paso_SystemMatrix_startCollect(flux_matrix,RP);
  remote_RP=Paso_SystemMatrix_finishCollect(flux_matrix);
  #pragma omp parallel for schedule(static) private(i, RN_i, iptr_ij, j, f_ij, f_i)
  for (i = 0; i < n; ++i) {
     RN_i=RN[i];
     f_i=0;
     #pragma ivdep
     for (iptr_ij=(pattern->couplePattern->ptr[i]);iptr_ij<pattern->couplePattern->ptr[i+1]; ++iptr_ij) {
          j=pattern->couplePattern->index[iptr_ij];
          f_ij=flux_matrix->coupleBlock->val[iptr_ij];
          if (f_ij < 0) {
              f_i+=f_ij*MIN(RN_i,remote_RP[j]);
          }
     }
     f[i]+=f_i;
  }
}
