
/*****************************************************************************
*
* Copyright (c) 2010-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


#include "error.h"
#include "Esys_MPI.h"

#include <stdio.h>	/* For FILENAME_MAX */
#include <string.h>
#include <time.h>
#ifdef _OPENMP 
#include <omp.h>
#endif

#define LenString_MAX FILENAME_MAX*2
#define LenErrorMsg_MAX LenString_MAX

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
bool Esys_noError(void)
{
   Esys_ErrorCodeType err=Esys_getErrorType();
   /* return (err==NO_ERROR ||  err==WARNING);*/
   return (err==NO_ERROR);
}

/* This function returns a timer */
double Esys_timer(void)
{
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
Esys_ErrorCodeType Esys_getErrorType(void)
{
   return Esys_ErrorCode_;
}

/* return the error message */
char* Esys_getErrorMessage(void)
{
   return Esys_ErrorMsg_;
}

