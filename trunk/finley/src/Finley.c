/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

/**************************************************************/

/*    Finley finite element solver */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

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


/**************************************************************/


/*
 * $Log$
 * Revision 1.3  2005/09/15 03:44:22  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.2.2.1  2005/09/07 06:26:18  gross
 * the solver from finley are put into the standalone package paso now
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
 * Finley C code has been included
 *
 * Revision 1.1.1.1  2004/06/24 04:00:40  johng
 * Initial version of eys using boost-python.
 *
 *
 */
