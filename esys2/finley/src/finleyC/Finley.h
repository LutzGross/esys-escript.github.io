/* $Id$ */

#ifndef INC_FINLEY
#define INC_FINLEY

/**************************************************************/

/*    Finley finite element solver */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003 */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"

/**************************************************************/

/* Atom types:                                                                 */
#define FINLEY_DEGREES_OF_FREEDOM 1
#define FINLEY_REDUCED_DEGREES_OF_FREEDOM 2
#define FINLEY_NODES 3
#define FINLEY_ELEMENTS 4
#define FINLEY_FACE_ELEMENTS 5
#define FINLEY_POINTS 6
#define FINLEY_CONTACT_ELEMENTS_1 7
#define FINLEY_CONTACT_ELEMENTS_2 8



/* error codes */

enum Finley_ErrorCodeType {
  NO_ERROR,
  WARNING,
  VALUE_ERROR,
  TYPE_ERROR,
  MEMORY_ERROR,
  IO_ERROR,
  ZERO_DIVISION_ERROR,
  EOF_ERROR,
  FLOATING_POINT_ERROR,
  INDEX_ERROR,
  OS_ERROR,
  OVERFLOW_ERROR,
  SYSTEM_ERROR
};

/* interfaces */

extern enum Finley_ErrorCodeType Finley_ErrorCode;
extern char Finley_ErrorMsg[LenErrorMsg_MAX];


double Finley_timer(void);
int Finley_checkPtr(void*);

#endif /* #ifndef INC_FINLEY */

/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:57  jgs
 * Initial revision
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
