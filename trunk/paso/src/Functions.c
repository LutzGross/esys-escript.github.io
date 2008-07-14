/* $Id:$ */

/*******************************************************
 *
 *       Copyright 2008 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#include "Common.h"
#include "Functions.h"
#include "PasoUtil.h"
/*
 * numerical calculation of the directional derivative J0w if F at x0 in the direction w. f0 is the value of F at x0.
 * setoff is workspace
 */

err_t Paso_FunctionDerivative(double* J0w, const double* w, Paso_Function* F, const double *f0, const double *x0, double* setoff)
{
   err_t err=0;
   dim_t n=F->n;
   double norm_w,epsnew,norm_x0;
   epsnew=10.*sqrt(EPSILON);
   norm_w=Paso_l2(n,w,F->mpi_info);

   if (norm_w>0) {
        Paso_zeroes(n,J0w);
   } else {
       epsnew = epsnew/norm_w;
       norm_x0=Paso_l2(n,x0,F->mpi_info);
       if (norm_x0>0) epsnew*=norm_x0;
       Paso_LinearCombination(n,setoff,1.,x0,epsnew,w);
       err=Paso_FunctionCall(F,J0w,setoff);
       if (err==NO_ERROR) {
           Paso_Update(n,1./epsnew,J0w,-1./epsnew,f0); /* J0w = (J0w - f0)/epsnew; */
       }
   }
   return err;
}

/*
 * sets value=F(arg)
 *
 */
err_t Paso_FunctionCall(Paso_Function * F,double* value, const double* arg) 
{ 
   if (F!=NULL) {
      switch(F->kind) {
          case LINEAR_SYSTEM:
               return Paso_Function_LinearSystem_call(F, value, arg);
          case FCT:
               return Paso_Function_FCT_call(F, value, arg);
          default:
               return SYSTEM_ERROR;
      }
   }
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
          case FCT:
               Paso_Function_FCT_free(F);
               break;
          default:
               MEMFREE(F);
      }
   }
}
