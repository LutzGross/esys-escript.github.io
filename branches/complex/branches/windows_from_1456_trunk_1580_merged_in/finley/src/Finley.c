
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

/*    Finley finite element solver */

/**************************************************************/

#include "Finley.h"

/* This function returns a time mark */
double Finley_timer(void) {
   return Paso_timer();
}

/* This function checks if the pointer ptr has a target. If not an
   error is raised and TRUE is returned. */
bool_t Finley_checkPtr(void* arg) {
   return Paso_checkPtr(arg);
}

/* reset the error to NO_ERROR */
void Finley_resetError(void) {
  Paso_resetError();
}

/* sets an error */
void Finley_setError(Finley_ErrorCodeType err,char* msg) {
  Paso_setError(err,msg);
}

/* checks if there is no error */
bool_t Finley_noError(void) {
   return Paso_noError();
}

/* return the error code */
Finley_ErrorCodeType Finley_getErrorType(void) {
    return Paso_getErrorType();
}

/* return the error message */
char* Finley_getErrorMessage(void) {
  return Paso_getErrorMessage();
}
/* return the error message */
void Finley_convertPasoError(void) {
  /* nothing has to be done here */
}

/* checks that there is no error accross all processes in a communicator */
/* NOTE : does not make guarentee consistency of error string on each process */
bool_t Finley_MPI_noError( Paso_MPIInfo *mpi_info )
{
    return Paso_MPIInfo_noError( mpi_info );
}


