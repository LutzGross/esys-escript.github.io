
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "Common.h"
#include "Functions.h"
#include "PasoUtil.h"
#include "Solver.h"
/*
 * numerical calculation of the directional derivative J0w if F at x0 in the direction w. f0 is the value of F at x0.
 * setoff is workspace
 */

err_t Paso_FunctionDerivative(double* J0w, const double* w, Paso_Function* F, const double *f0, const double *x0, double* setoff, const bool_t w_is_normalized, Paso_Performance* pp) 
{
   err_t err=0;
   dim_t n=F->n;
   double norm_w,epsnew,norm_x0;
   epsnew=10.*sqrt(EPSILON);
   if (! w_is_normalized) {
        norm_w=Paso_l2(n,w,F->mpi_info);
   } else {
        norm_w=1.;
   }

   if (norm_w>0) {
       epsnew = epsnew/norm_w;
       norm_x0=Paso_l2(n,x0,F->mpi_info);
       if (norm_x0>0) epsnew*=norm_x0;
       Paso_LinearCombination(n,setoff,1.,x0,epsnew,w);
       err=Paso_FunctionCall(F,J0w,setoff,pp);
       if (err==NO_ERROR) {
           Paso_Update(n,1./epsnew,J0w,-1./epsnew,f0); /* J0w = (J0w - f0)/epsnew; */
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
          default:
               return SYSTEM_ERROR;
      }
   }
   /* Added by PGH, assume a null pointe is an error */
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
