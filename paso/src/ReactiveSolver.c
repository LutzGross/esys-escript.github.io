
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


    
err_t  Paso_ReactiveSolver_solve(Paso_ReactiveSolver* support, Paso_TransportProblem* fctp, double* u, const double dt, const double* source, Paso_Options* options, Paso_Performance *pp)
{
     const double EXP_LIM_MIN =PASO_RT_EXP_LIM_MIN;
     const double EXP_LIM_MAX =PASO_RT_EXP_LIM_MAX;
     index_t fail=0;
     register double d_ii, m_i, x_i, e_i, u_i, F_i;
     dim_t i;
     const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
    
     #pragma omp parallel for schedule(static) private(i, d_ii, m_i, x_i, e_i, u_i, F_i) 
     for (i=0;i<n;++i) {
        d_ii=fctp->reactive_matrix[i];
        m_i=fctp->lumped_mass_matrix[i];
        x_i=dt*d_ii/m_i;
	if (x_i >= EXP_LIM_MAX) {
	  fail=1;
	} else  {
	    F_i=source[i];
	    e_i=exp(x_i);
	    u_i=e_i*u[i];
	    if ( abs(x_i) > EXP_LIM_MIN) {
		u_i+=F_i/d_ii*(e_i-1.);
	    } else {
		u_i+=F_i*dt/m_i * (1. + x_i/2); /* second order approximation of ( exp(x_i)-1)/x_i */
	    }
	    u[i]=u_i;
	}
    }
    #ifdef PASO_MPI
    {
        index_t fail_loc = fail;
        MPI_Allreduce(&fail_loc, &fail, 1, MPI_INT, MPI_MAX, fctp->mpi_info->comm);
    }
    #endif
    if (fail < 0 ) {
        return SOLVER_DIVERGENCE;
    } else {
        return SOLVER_NO_ERROR;
    }
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

     const double EXP_LIM_MAX =PASO_RT_EXP_LIM_MAX;
     double dt_max=LARGE_POSITIVE_FLOAT, dt_max_loc;  
     index_t fail_loc, fail;
     dim_t i;
     const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
     register double d_ii,m_i;
     if (Paso_noError()) {
        /*
         *  calculate time step size:                                           
        */
        dt_max=LARGE_POSITIVE_FLOAT;
	fail=0;
        #pragma omp parallel private(dt_max_loc, fail_loc)
        {
               dt_max_loc=LARGE_POSITIVE_FLOAT;
	       fail_loc=0;
               #pragma omp for schedule(static) private(i,d_ii,m_i) 
               for (i=0;i<n;++i) {
                  d_ii=fctp->reactive_matrix[i];
                  m_i=fctp->lumped_mass_matrix[i];
		  if (m_i > 0) {
		      if (d_ii>0) dt_max_loc=MIN(dt_max_loc, m_i/d_ii);
		  } else {
		      fail_loc=-1;
		  }
               }
               #pragma omp critical 
               {
                  dt_max=MIN(dt_max,dt_max_loc);
		  fail=MIN(fail, fail_loc);
               }
        }
        #ifdef PASO_MPI
        {
	       double rtmp_loc[2], rtmp[2];
               rtmp_loc[0]=dt_max;
	       rtmp_loc[1]= (double) fail;
               MPI_Allreduce(rtmp_loc, rtmp, 2, MPI_DOUBLE, MPI_MIN, fctp->mpi_info->comm);
	       dt_max=rtmp[0];
	       fail = rtmp[1] < 0 ? -1 : 0;
	}
        #endif
        if (fail < 0 ) {
	   Paso_setError(VALUE_ERROR, "Paso_ReactiveSolver_getSafeTimeStepSize: negative mass term detected.");
	   return -1;
	} else {
	    if (dt_max<LARGE_POSITIVE_FLOAT) dt_max*=0.5*EXP_LIM_MAX;
	}
   }
   return dt_max;
}
void Paso_ReactiveSolver_initialize(const double dt, Paso_TransportProblem* fctp, Paso_Options* options)
{
}
