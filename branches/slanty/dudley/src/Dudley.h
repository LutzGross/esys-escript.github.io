
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

#ifndef INC_DUDLEY
#define INC_DUDLEY

/************************************************************************************/

/*    Dudley finite element solver */

/************************************************************************************/

#include "esysUtils/types.h"
#include "esysUtils/Esys_MPI.h"
#include "esysUtils/error.h"

/************************************************************************************/
/*#define Dudley_TRACE */
#define DUDLEY_UNKNOWN -1
#define DUDLEY_DEGREES_OF_FREEDOM 1
#define DUDLEY_NODES 3
#define DUDLEY_ELEMENTS 4
#define DUDLEY_FACE_ELEMENTS 5
#define DUDLEY_POINTS 6
#define DUDLEY_REDUCED_DEGREES_OF_FREEDOM 2
#define DUDLEY_REDUCED_NODES 14
#define DUDLEY_REDUCED_ELEMENTS 10
#define DUDLEY_REDUCED_FACE_ELEMENTS 11

/* status stuff */
typedef int Dudley_Status_t;
#define Dudley_increaseStatus(self) ((self)->status)++
#define DUDLEY_INITIAL_STATUS 0

/* error codes */

typedef Esys_ErrorCodeType Dudley_ErrorCodeType;

/* interfaces */

double Dudley_timer(void);
bool_t Dudley_checkPtr(void *);
void Dudley_resetError(void);
void Dudley_setError(Dudley_ErrorCodeType err, __const char *msg);
bool_t Dudley_noError(void);
Dudley_ErrorCodeType Dudley_getErrorType(void);
char *Dudley_getErrorMessage(void);
void Dudley_convertPasoError(void);
bool_t Dudley_MPI_noError(Esys_MPIInfo * mpi_info);
void Dudley_setTagsInUse(const index_t Tag, const dim_t numTags, dim_t * numTagsInUse, index_t ** tagsInUse,
			 Esys_MPIInfo * mpiinfo);

#endif				/* #ifndef INC_DUDLEY */
