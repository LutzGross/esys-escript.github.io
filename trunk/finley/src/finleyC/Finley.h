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

#define FINLEY_UNKNOWN -1
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

