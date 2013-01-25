
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*    Paso finite element solver */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003 */
/*   Version: $Id$ */

/**************************************************************/

#include "Paso.h"

#ifdef _OPENMP 
#include <omp.h>
#endif

#ifdef PASO_MPI
#include "mpi_C.h"
#else
#include <time.h>
#endif

Paso_ErrorCodeType Paso_ErrorCode_=NO_ERROR;
char Paso_ErrorMsg_[LenErrorMsg_MAX]={'\0'};

/* reset the error to NO_ERROR */
void Paso_resetError(void) {
  Paso_ErrorCode_=NO_ERROR;
}
                                                                                                                                                                                                     
/* sets an error */
void Paso_setError(Paso_ErrorCodeType err,__const char* msg) {
  size_t lenMsg=strlen(msg);
  if (Paso_noError()) {
printf("error set = %d %s\n",err,msg);
     Paso_ErrorCode_=err;
     strncpy(Paso_ErrorMsg_,msg,MIN(LenErrorMsg_MAX,lenMsg));
     Paso_ErrorMsg_[MIN(LenErrorMsg_MAX,lenMsg)]='\0';
  }
}
                                                                                                                                                                                                     
/* checks if there is no error */
bool_t Paso_noError(void) {
   Paso_ErrorCodeType err=Paso_getErrorType();
   /* return (err==NO_ERROR ||  err==WARNING);*/
   return (err==NO_ERROR);
}
/* This function checks if the pointer ptr has a target. If not an
   error is raised and TRUE is returned. */

bool_t Paso_checkPtr(void* ptr) {
   if (ptr==NULL) {
      Paso_setError(MEMORY_ERROR,"Out of memory.");
      return TRUE;
   } else {
      return FALSE;
   }
} 

/* This function returns a timer */
double Paso_timer(void) {
  double out;

#ifdef PASO_MPI
  out = MPI_Wtime();
#else
#ifdef _OPENMP 
  out=omp_get_wtime();
#else
  out=((double) clock())/CLOCKS_PER_SEC;
#endif
#endif
  return out;
}
#ifndef _OPENMP 
int omp_get_max_threads(void) {
  return 1;
}
#endif


/* return the error code */
Paso_ErrorCodeType Paso_getErrorType(void) {
   return Paso_ErrorCode_;
}

/* return the error message */
char* Paso_getErrorMessage(void) {
   return Paso_ErrorMsg_;
}
/**************************************************************/
