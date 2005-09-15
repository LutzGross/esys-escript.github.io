/* $Id$ */

/**************************************************************/

/*    Paso finite element solver */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003 */
/*   Version: $Id$ */

/**************************************************************/

#include "Paso.h"
#ifdef _OPENMP 
#include <omp.h>
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
void Paso_setError(Paso_ErrorCodeType err,char* msg) {
  size_t lenMsg=strlen(msg);
  if (Paso_noError()) {
     Paso_ErrorCode_=err;
     strncpy(Paso_ErrorMsg_,msg,MIN(LenErrorMsg_MAX,lenMsg));
     Paso_ErrorMsg_[MIN(LenErrorMsg_MAX,lenMsg)]='\0';
  }
}
                                                                                                                                                                                                     
/* checks if there is no error */
bool_t Paso_noError(void) {
   Paso_ErrorCodeType err=Paso_getErrorType();
   return (err==NO_ERROR ||  err==WARNING);
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
  #ifdef _OPENMP 
  out=omp_get_wtime();
  #else
  out=((double) clock())/CLOCKS_PER_SEC;
  #endif
  return out;
}

/* return the error code */
Paso_ErrorCodeType Paso_getErrorType(void) {
   return Paso_ErrorCode_;
}

/* return the error message */
char* Paso_getErrorMessage(void) {
   return Paso_ErrorMsg_;
}
                                                                                                                                                                                                     


/**************************************************************/


/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:38  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.3  2005/09/08 08:28:39  gross
 * some cleanup in savevtk
 *
 * Revision 1.1.2.2  2005/09/07 00:59:08  gross
 * some inconsistent renaming fixed to make the linking work.
 *
 * Revision 1.1.2.1  2005/09/05 06:29:47  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 * Revision 1.2  2005/07/08 04:07:50  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.1  2005/06/29 02:34:50  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.2  2004/07/02 04:21:13  gross
 * Paso C code has been included
 *
 * Revision 1.1.1.1  2004/06/24 04:00:40  johng
 * Initial version of eys using boost-python.
 *
 *
 */
