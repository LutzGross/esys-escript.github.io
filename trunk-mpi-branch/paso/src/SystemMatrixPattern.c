/* $Id$ */

/*
********************************************************************************
*               Copyright   2006 by ACcESS MNRF                                *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

/**************************************************************/

/* Paso: SystemMatrixPatternPattern */

/**************************************************************/
 
/* Copyrights by ACcESS Australia 2003, 2004,2005, 2006, 2007 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrixPattern.h"

/**************************************************************/

/* allocates a SystemMatrixPattern  */

Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_alloc(int type,
                                                         Paso_Distribution *output_distribution,
                                                         Paso_Distribution *input_distribution,
                                                         Paso_Pattern* mainPattern,
                                                         Paso_Pattern* couplePattern,
                                                         Paso_Coupler* coupler) 
{
  Paso_SystemMatrixPattern*out=NULL;


  Paso_resetError();

  if (mainPattern->type != type)  {
      Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: type of mainPattern does not match expected type.");
  }
  if (couplePattern->type != type)  {
      Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: type of couplePattern does not match expected type.");
  }
  if ( couplePattern->numOutput != mainPattern->numOutput) {
            Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: number of output for couple and main pattern don't match.");
  }
  if (mainPattern->numOutput !=  Paso_Distribution_getMyNumComponents(output_distribution)) {
      Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: number of output and given distribution don't match.");
  }
  if (mainPattern->numInput != Paso_Distribution_getMyNumComponents(input_distribution)) {
     Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: number of input for main pattern and number of send components in coupler don't match.");
  }
printf("couplePattern->numInput != coupler->recv->numSharedComponents %d %d\n", couplePattern->numInput, coupler->recv->numSharedComponents);
  if (couplePattern->numInput != coupler->recv->numSharedComponents) {
     Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: number of inputs for couple pattern and number of received components in coupler don't match.");
  }
  out=MEMALLOC(1,Paso_SystemMatrixPattern);
  if (Paso_checkPtr(out)) return NULL;
  out->type=type;
  out->reference_counter=1;
  out->mainPattern=Paso_Pattern_getReference(mainPattern);
  out->couplePattern=Paso_Pattern_getReference(couplePattern);
  out->coupler=Paso_Coupler_getReference(coupler);
  out->output_distribution=Paso_Distribution_getReference(output_distribution);
  out->input_distribution=Paso_Distribution_getReference(input_distribution);
  out->mpi_info= Paso_MPIInfo_getReference(coupler->mpi_info);
  #ifdef Paso_TRACE
  printf("Paso_SystemMatrixPattern_dealloc: system matrix pattern as been allocated.\n");
  #endif
  return out;
}

/* returns a reference to in */

Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_reference(Paso_SystemMatrixPattern* in) {
     if (in!=NULL) {
        ++(in->reference_counter);
     }
     return in;
}
  
/* deallocates a SystemMatrixPattern: */

void Paso_SystemMatrixPattern_free(Paso_SystemMatrixPattern* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<=0) {
        Paso_Pattern_free(in->mainPattern);
        Paso_Pattern_free(in->couplePattern);
        Paso_Coupler_free(in->coupler);
        Paso_Distribution_free(in->output_distribution);
        Paso_Distribution_free(in->input_distribution);
        Paso_MPIInfo_free(in->mpi_info);
        MEMFREE(in);
        #ifdef Paso_TRACE
        printf("Paso_SystemMatrixPattern_free: system matrix pattern as been deallocated.\n");
        #endif
     }
   }
}
