
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*    Dudley finite element solver */

/**************************************************************/

#include "Dudley.h"

/* This function returns a time mark */
double Dudley_timer(void) {
   return Paso_timer();
}

/* This function checks if the pointer ptr has a target. If not an
   error is raised and TRUE is returned. */
bool_t Dudley_checkPtr(void* arg) {
   return Paso_checkPtr(arg);
}

/* reset the error to NO_ERROR */
void Dudley_resetError(void) {
  Paso_resetError();
}

/* sets an error */
void Dudley_setError(Dudley_ErrorCodeType err,__const char* msg) {
  Paso_setError(err,msg);
}

/* checks if there is no error */
bool_t Dudley_noError(void) {
   return Paso_noError();
}

/* return the error code */
Dudley_ErrorCodeType Dudley_getErrorType(void) {
    return Paso_getErrorType();
}

/* return the error message */
char* Dudley_getErrorMessage(void) {
  return Paso_getErrorMessage();
}
/* return the error message */
void Dudley_convertPasoError(void) {
  /* nothing has to be done here */
}

/* checks that there is no error accross all processes in a communicator */
/* NOTE : does not make guarentee consistency of error string on each process */
bool_t Dudley_MPI_noError( Paso_MPIInfo *mpi_info )
{
    return Paso_MPIInfo_noError( mpi_info );
}


