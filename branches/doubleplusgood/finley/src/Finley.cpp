
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


/************************************************************************************/

/*    Finley finite element solver */

/************************************************************************************/

#include "Finley.h"
#include "esysUtils/error.h"

/* This function returns a time mark */
double Finley_timer(void) {
   return Esys_timer();
}

/* This function checks if the pointer ptr has a target. If not an
   error is raised and TRUE is returned. */
bool_t Finley_checkPtr(void* arg) {
   return Esys_checkPtr(arg);
}

/* reset the error to NO_ERROR */
void Finley_resetError(void) {
  Esys_resetError();
}

/* sets an error */
void Finley_setError(Finley_ErrorCodeType err,__const char* msg) {
  Esys_setError(err,msg);
}

/* checks if there is no error */
bool_t Finley_noError(void) {
   return Esys_noError();
}

/* return the error code */
Finley_ErrorCodeType Finley_getErrorType(void) {
    return Esys_getErrorType();
}

/* return the error message */
char* Finley_getErrorMessage(void) {
  return Esys_getErrorMessage();
}
/* return the error message */
void Finley_convertPasoError(void) {
  /* nothing has to be done here */
}

/* checks that there is no error across all processes in a communicator */
/* NOTE : does not guarantee consistency of error string on each process */
bool_t Finley_MPI_noError( Esys_MPIInfo *mpi_info )
{
    return Esys_MPIInfo_noError( mpi_info );
}


