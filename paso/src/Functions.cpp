
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
#include "Functions.h"
#include "PasoUtil.h"
#include "Solver.h"
/*
 * numerical calculation of the directional derivative J0w if F at x0 in the direction w. f0 is the value of F at x0.
 * setoff is workspace
 */

err_t Paso_FunctionDerivative(double* J0w, const double* w, Paso_Function* F, const double *f0, const double *x0, double* setoff, Paso_Performance* pp) 
{
   err_t err=SOLVER_NO_ERROR;
   const dim_t n=F->n;
   dim_t i;
   register double aw;
   const double epsnew=sqrt(EPSILON);
   double /*norm_x0,*/ ttt, s=epsnew, local_s, norm_w=0.;
   
   /*norm_x0=Paso_lsup(n,x0,F->mpi_info);*/
   norm_w=Paso_lsup(n,w,F->mpi_info);
   ttt=sqrt(EPSILON)*norm_w;
   #pragma omp parallel private(local_s) 
   {
      local_s=s;
      #pragma omp for private(i, aw) 
      for (i=0;i<n;++i) {
	   aw=fabs(w[i]);
	   if ( aw>ttt ) {
	      local_s=MAX(local_s,fabs(x0[i])/aw);
	   }
      }
      #pragma omp critical
      {
            s=MAX(s,local_s);
      }
      
   }
   #ifdef ESYS_MPI
   {
       double local_v[2], v[2];
       local_v[0]=s;
       local_v[1]=norm_w;
       MPI_Allreduce(local_v,v, 2, MPI_DOUBLE, MPI_MAX, F->mpi_info->comm);
       s=v[0];
       norm_w=v[1];
   }
   #endif
/* printf("s ::  = %e, %e \n",s, norm_x0/norm_w); */
   if (norm_w>0) {
        s=s*epsnew;
/* printf("s = %e\n",s); */
        Paso_LinearCombination(n,setoff,1.,x0,s,w); 
        err=Paso_FunctionCall(F,J0w,setoff,pp);
        if (err==SOLVER_NO_ERROR) {
              Paso_Update(n,1./s,J0w,-1./s,f0); /* J0w = (J0w - f0)/epsnew; */
/*	      {
		int i;
		for (i=0;i<n; i++) printf("df[%d]=%e %e\n",i,J0w[i],w[i]);
	      }
*/
           }
   } else {
       Paso_zeroes(n,J0w);
   }
   return err;
}

/*
 * sets value=F(arg)
 *
 */
err_t Paso_FunctionCall(Paso_Function * F,double* value, const double* arg, Paso_Performance *pp) 
{ 
   if (F!=NULL) {
      switch(F->kind) {
          case LINEAR_SYSTEM:
               return Paso_Function_LinearSystem_call(F, value, arg,pp);
	       break;
          default:
               return SYSTEM_ERROR;
      }
   }
   /* Added by PGH, assume a null pointer is an error */
   return SYSTEM_ERROR;
}

/*
 * clear Paso_Function
 */
void Paso_Function_free(Paso_Function * F) {
   if (F!=NULL) {

      switch(F->kind) {
          case LINEAR_SYSTEM:
               Paso_Function_LinearSystem_free(F);
               break;
          default:
               MEMFREE(F);
      }
   }
}

