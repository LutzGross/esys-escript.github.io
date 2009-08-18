
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


/**************************************************************/

/*   Paso: system matrix pattern                            */

/**************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

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


PASO_DLL_API
Paso_Distribution*  Paso_Distribution_alloc( Paso_MPIInfo *mpi_info, index_t* first_component, index_t m, index_t b);

PASO_DLL_API
void                Paso_Distribution_free( Paso_Distribution *in );

PASO_DLL_API
Paso_Distribution*  Paso_Distribution_getReference( Paso_Distribution *in );

PASO_DLL_API
index_t Paso_Distribution_getFirstComponent(Paso_Distribution *in );

PASO_DLL_API
index_t Paso_Distribution_getLastComponent(Paso_Distribution *in );

PASO_DLL_API
dim_t Paso_Distribution_getGlobalNumComponents(Paso_Distribution *in );

PASO_DLL_API
dim_t Paso_Distribution_getMyNumComponents(Paso_Distribution *in );

PASO_DLL_API
dim_t Paso_Distribution_getMinGlobalComponents(Paso_Distribution *in );

PASO_DLL_API
dim_t Paso_Distribution_getMaxGlobalComponents(Paso_Distribution *in );
#endif
