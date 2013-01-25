
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/*   Paso: system matrix pattern                            */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2004,2005,2007 */
/*   Author: gross@access.edu.au */

/**************************************************************/

#include "Distribution.h"

Paso_Distribution* Paso_Distribution_alloc( Paso_MPIInfo *mpi_info, 
                                            index_t *first_component,
                                            index_t m, index_t b) 
{
  int i;
  Paso_Distribution *out=NULL;
  out = MEMALLOC( 1, Paso_Distribution );
  if (Paso_checkPtr(out)) return NULL;
  out->mpi_info = Paso_MPIInfo_getReference(mpi_info);
  out->reference_counter = 0;
  out->first_component=NULL;

  out->first_component = MEMALLOC( (mpi_info->size)+1, index_t );
  if (Paso_checkPtr(out->first_component)) {
       Paso_Distribution_free(out);
       return NULL;
  }
  for (i=0; i<(mpi_info->size)+1; ++i) out->first_component[i]=m*first_component[i]+b;
  out->reference_counter++;
  return out;
}

void Paso_Distribution_free( Paso_Distribution *in )
{
  if (in != NULL) {
    --(in->reference_counter);
    if (in->reference_counter<=0) {
      Paso_MPIInfo_free( in->mpi_info );
      MEMFREE( in->first_component );
      MEMFREE( in );
    } 
  }
}

Paso_Distribution* Paso_Distribution_getReference( Paso_Distribution *in )
{
  if ( in != NULL) 
    in->reference_counter++;
  return in;
}

index_t Paso_Distribution_getFirstComponent(Paso_Distribution *in ) {
 if (in !=NULL) {
   return in->first_component[in->mpi_info->rank];
 } else {
   return 0;
 }
}
index_t Paso_Distribution_getLastComponent(Paso_Distribution *in ){
 if (in !=NULL) {
   return in->first_component[(in->mpi_info->rank)+1];
 } else {
   return 0;
 }
}

dim_t Paso_Distribution_getGlobalNumComponents(Paso_Distribution *in ){
   return Paso_Distribution_getMaxGlobalComponents(in)-Paso_Distribution_getMinGlobalComponents(in);
}
dim_t Paso_Distribution_getMyNumComponents(Paso_Distribution *in ) {
   return Paso_Distribution_getLastComponent(in)-Paso_Distribution_getFirstComponent(in);
}

dim_t Paso_Distribution_getMinGlobalComponents(Paso_Distribution *in ){
 if (in !=NULL) {
   return in->first_component[0];
 } else {
   return 0;
 }
}
dim_t Paso_Distribution_getMaxGlobalComponents(Paso_Distribution *in ){
 if (in !=NULL) {
   return in->first_component[in->mpi_info->size];
 } else {
   return 0;
 }
}

