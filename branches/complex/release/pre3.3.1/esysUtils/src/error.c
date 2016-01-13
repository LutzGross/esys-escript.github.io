
/*****************************************************************************
*
* Copyright (c) 2010-2013 by University of Queensland
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


#include <string.h>
#include "error.h"

#ifdef _OPENMP 
#include <omp.h>
#endif

#ifdef ESYS_MPI
#include "mpi_C.h"
#else
#include <time.h>
#endif



#define MIN(X,Y) ((X)<(Y)?(X):(Y))



Esys_ErrorCodeType Esys_ErrorCode_=NO_ERROR;
char Esys_ErrorMsg_[LenErrorMsg_MAX]={'\0'};

/* reset the error to NO_ERROR */
void Esys_resetError(void) {
  Esys_ErrorCode_=NO_ERROR;
}
                                                                                                                                                                                                     
/* sets an error */
void Esys_setError(Esys_ErrorCodeType err,__const char* msg) {
  size_t lenMsg=strlen(msg);
  if (Esys_noError()) {
/* printf("error set = %d %s\n",err,msg); */
     Esys_ErrorCode_=err;
     strncpy(Esys_ErrorMsg_,msg,MIN(LenErrorMsg_MAX,lenMsg));
     Esys_ErrorMsg_[MIN(LenErrorMsg_MAX,lenMsg)]='\0';
  }
}
                                                                                                                                                                                                     
/* checks if there is no error */
bool_t Esys_noError(void) {
   Esys_ErrorCodeType err=Esys_getErrorType();
   /* return (err==NO_ERROR ||  err==WARNING);*/
   return (err==NO_ERROR);
}
/* This function checks if the pointer ptr has a target. If not an
   error is raised and TRUE is returned. */

bool_t Esys_checkPtr(void* ptr) {
   if (ptr==NULL) {
      Esys_setError(MEMORY_ERROR,"Out of memory.");
      return TRUE;
   } else {
      return FALSE;
   }
} 

/* This function returns a timer */
double Esys_timer(void) {
  double out;

#ifdef ESYS_MPI
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



/* return the error code */
Esys_ErrorCodeType Esys_getErrorType(void) {
   return Esys_ErrorCode_;
}

/* return the error message */
char* Esys_getErrorMessage(void) {
   return Esys_ErrorMsg_;
}
/************************************************************************************/
