/* $Id$ */

/**************************************************************/

/*    Finley finite element solver */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003 */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"
#include "Finley.h"
#ifdef _OPENMP 
#include <omp.h>
#else 
#include <time.h>
#endif



/* This function checks if the pointer ptr has a target. If not an
   error is raised and TRUE is returned. */

int Finley_checkPtr(void* ptr) {
   if (ptr==NULL) {
      Finley_ErrorCode=MEMORY_ERROR;
      sprintf(Finley_ErrorMsg,"Out of memory.");
      return TRUE;
   } else {
      return FALSE;
   }
} 

/* This function returns a timer */
double Finley_timer(void) {
  double out;
  #ifdef _OPENMP 
  out=omp_get_wtime();
  #else
  out=((double) clock())/CLOCKS_PER_SEC;
  #endif
  return out;
}

/**************************************************************/


/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:57  jgs
 * Initial revision
 *
 * Revision 1.2  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 * Revision 1.1.1.1  2004/06/24 04:00:40  johng
 * Initial version of eys using boost-python.
 *
 *
 */
