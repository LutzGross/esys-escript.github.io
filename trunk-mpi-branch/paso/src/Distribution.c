/*
********************************************************************************
*               Copyright  2006,2007 by ACcESS MNRF                            *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0                     *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

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
  out->numComponents=first_component[mpi_info->size]-first_component[0];
  out->firstComponent=first_component[0];
  out->myNumComponents=first_component[(mpi_info->rank)+1]-first_component[mpi_info->rank];
  out->maxNumComponents=out->myNumComponents;
  for (i=0; i< mpi_info->size; ++i) out->maxNumComponents=MAX(out->maxNumComponents,first_component[i+1]-first_component[i]);
  out->myFirstComponent=first_component[mpi_info->rank];
  return out;
}

void Paso_Distribution_free( Paso_Distribution *in )
{
  index_t i;

  if ( in && !(--in->reference_counter) )
  {
    Paso_MPIInfo_dealloc( in->mpi_info );
    MEMFREE( in->first_component );
    MEMFREE( in );
  } 
}

Paso_Distribution* Paso_Distribution_getReference( Paso_Distribution *in )
{
  if ( in ) 
    in->reference_counter++;
  
  return in;
}
