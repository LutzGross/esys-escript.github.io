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

bool_t Finley_checkPtr(void* ptr) {
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
 * Finley C code has been included
 *
 * Revision 1.1.1.1  2004/06/24 04:00:40  johng
 * Initial version of eys using boost-python.
 *
 *
 */
