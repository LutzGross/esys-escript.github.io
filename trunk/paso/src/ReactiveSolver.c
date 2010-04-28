
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: reactive solver (D is a diagonal matrix) 
 *        
 *   - Mv_t=Dv+q   v(0)=u       
 *
 *  to return v(dt)
 *
*/

/**************************************************************/

/* Author: l.gross@uq.edu.au                                */

/**************************************************************/


#include "ReactiveSolver.h"
#include "PasoUtil.h"
#include "Solver.h"

err_t  Paso_ReactiveSolver_solve(Paso_ReactiveSolver* support, Paso_TransportProblem* fctp, double* u, double dt, double* source, Paso_Options* options, Paso_Performance *pp)
{
    return SOLVER_NO_ERROR;
}

Paso_ReactiveSolver* Paso_ReactiveSolver_alloc(Paso_TransportProblem* fctp)
{
    Paso_ReactiveSolver* out=NULL;
    out=MEMALLOC(1,Paso_ReactiveSolver);
    if (Paso_checkPtr(out)) return NULL;
    return out;
}

void Paso_ReactiveSolver_free(Paso_ReactiveSolver* in)
{
   if (in!=NULL) {
      MEMFREE(in);
   }
}
double Paso_ReactiveSolver_getSafeTimeStepSize(Paso_TransportProblem* fctp)
{
  double dt_max=LARGE_POSITIVE_FLOAT;



  double beta_0=0., beta_0_loc, betap, beta=EPSILON;
  dim_t i, n;
  int fail=0;
  double dt_max_loc;
  register double m_ii,m_i, d_ii;
  #ifdef PASO_MPI
     double rtmp[2], rtmp_loc[2];
  #endif
  n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix); 

  #pragma omp parallel private(beta_0_loc)
  {
      beta_0_loc=0.;
      #pragma omp for schedule(static) private(i,m_ii,m_i, d_ii) 
      for (i=0;i<n;++i) {
         m_i=fctp->lumped_mass_matrix[i];
         m_ii=fctp->main_diagonal_mass_matrix[i];
/* printf(" %d : %e %e -> %e\n", i, m_i, m_ii, m_i/m_ii); */
         if ( m_ii>0 ) {
            beta_0_loc=MAX(m_i/m_ii, beta_0_loc);
         } else {
            fail=1;
         }
      }
      #pragma omp critical 
      {
          beta_0=MAX(beta_0,beta_0_loc);
      }
  }
  #ifdef PASO_MPI
     rtmp_loc[0]=beta_0;
     rtmp_loc[1]=1.*((double) fail );
     MPI_Allreduce(rtmp_loc, rtmp, 2, MPI_DOUBLE, MPI_MAX, fctp->mpi_info->comm);
     beta_0=rtmp_loc[0];
     fail= (rtmp_loc[0] > 0) ? 1 : 0;
  #endif
  if (fail>0) {
     Paso_setError(TYPE_ERROR,"Paso_ReactiveSolver_getSafeTimeStepSize: non-positive main diagonal entry in mass matrix.");
  }
  
  if ( Paso_noError()) {
      if ( ( beta_0 < 2. ) && ( beta_0 >= 1. ) ) {
          beta=MAX(2*(beta_0-1.)/beta_0, EPSILON); 
       } else {
          Paso_setError(TYPE_ERROR,"Paso_ReactiveSolver_getSafeTimeStepSize: non-positive main diagonal entry in mass matrix.");
       }
  }
printf("beta = %e from beta_0 = %e\n",beta, beta_0-1);
  if ( Paso_noError()) {
     
     #pragma omp parallel private(dt_max_loc, betap)
     {
         betap=beta+1.;
         dt_max_loc=LARGE_POSITIVE_FLOAT;
         #pragma omp for schedule(static) private(i, m_ii, m_i, d_ii) 
         for (i=0;i<n;++i) {
             d_ii=fctp->reactive_matrix[i];
             m_i=fctp->lumped_mass_matrix[i];
             m_ii=fctp->main_diagonal_mass_matrix[i];

             if (d_ii > 0.)  dt_max_loc=MIN(dt_max_loc, (betap* m_ii - m_i)/d_ii );
         } 
         #pragma omp critical 
         {
             dt_max=MIN(dt_max,dt_max_loc);
        }
     }
     dt_max=1./(beta* Paso_Transport_getTheta(fctp));
      #ifdef PASO_MPI
         dt_max_loc=dt_max;
         MPI_Allreduce(&dt_max_loc, &dt_max, 1, MPI_DOUBLE, MPI_MIN, fctp->mpi_info->comm);
      #endif
  }
  return dt_max;
}
