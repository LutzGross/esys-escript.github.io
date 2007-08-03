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

#ifndef INC_PASO_DISTRIBUTION
#define INC_PASO_DISTRIBUTION

#include "Common.h"
#include "Paso_MPI.h"

/**************************************************** 
  describes the distribution of a vector stored
  on the local process
****************************************************/
struct Paso_Distribution
{
  index_t *first_component;  /* process i has nodes with global indices first_component[i+1] to first_component[i]. */
  dim_t reference_counter;
  Paso_MPIInfo *mpi_info;
};

typedef struct Paso_Distribution Paso_Distribution;

/***************************************
 Function prototypes 
**************************************/

Paso_Distribution*  Paso_Distribution_alloc( Paso_MPIInfo *mpi_info, index_t* first_component, index_t m, index_t b);
void                Paso_Distribution_free( Paso_Distribution *in );
Paso_Distribution*  Paso_Distribution_getReference( Paso_Distribution *in );
index_t Paso_Distribution_getFirstComponent(Paso_Distribution *in );
index_t Paso_Distribution_getLastComponent(Paso_Distribution *in );
dim_t Paso_Distribution_getGlobalNumComponents(Paso_Distribution *in );
dim_t Paso_Distribution_getMyNumComponents(Paso_Distribution *in );
dim_t Paso_Distribution_getMinGlobalComponents(Paso_Distribution *in );
dim_t Paso_Distribution_getMaxGlobalComponents(Paso_Distribution *in );
#endif
