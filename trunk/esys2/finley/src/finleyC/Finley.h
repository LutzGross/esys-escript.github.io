/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2003,2004,2005 -  All Rights Reserved              *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/


#ifndef INC_FINLEY
#define INC_FINLEY

/**************************************************************/

/*    Finley finite element solver */

/**************************************************************/

/*   Version: $Id$ */

/**************************************************************/

#include "paso/Paso.h"

/**************************************************************/

#define FINLEY_DEGREES_OF_FREEDOM 1
#define FINLEY_REDUCED_DEGREES_OF_FREEDOM 2
#define FINLEY_NODES 3
#define FINLEY_ELEMENTS 4
#define FINLEY_FACE_ELEMENTS 5
#define FINLEY_POINTS 6
#define FINLEY_CONTACT_ELEMENTS_1 7
#define FINLEY_CONTACT_ELEMENTS_2 8

/* error codes */


typedef Paso_ErrorCodeType Finley_ErrorCodeType;

/* interfaces */

double Finley_timer(void);
bool_t Finley_checkPtr(void*);
void Finley_resetError(void);
void Finley_setError(Finley_ErrorCodeType err,char* msg);
bool_t Finley_noError(void);
Finley_ErrorCodeType Finley_getErrorType(void);
char* Finley_getErrorMessage(void);
void Finley_convertPasoError(void);

#endif /* #ifndef INC_FINLEY */

/*
 * $Log$
 * Revision 1.3  2005/09/15 03:44:22  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.2.2.1  2005/09/07 06:26:18  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.2  2005/07/08 04:07:51  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.1  2005/06/29 02:34:50  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.3  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 * Revision 1.2  2004/06/29 01:59:31  johng
 * ??
 *
 * Revision 1.1.1.1  2004/06/24 04:00:40  johng
 * Initial version of eys using boost-python.
 *
 *
 */
