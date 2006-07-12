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


#ifndef INC_FINLEY
#define INC_FINLEY

/**************************************************************/

/*    Finley finite element solver */

/**************************************************************/

/*   Version: $Id$ */

/**************************************************************/

#include "paso/Paso.h"

/**************************************************************/
/*#define Finley_TRACE */
#define FINLEY_UNKNOWN -1
#define FINLEY_DEGREES_OF_FREEDOM 1
#define FINLEY_REDUCED_DEGREES_OF_FREEDOM 2
#define FINLEY_NODES 3
#define FINLEY_ELEMENTS 4
#define FINLEY_FACE_ELEMENTS 5
#define FINLEY_POINTS 6
#define FINLEY_CONTACT_ELEMENTS_1 7
#define FINLEY_CONTACT_ELEMENTS_2 8
/* these are not used yet */
#define FINLEY_REDUCED_NODES 9
#define FINLEY_REDUCED_ELEMENTS 10
#define FINLEY_REDUCED_FACE_ELEMENTS 11
#define FINLEY_REDUCED_CONTACT_ELEMENTS_1 12
#define FINLEY_REDUCED_CONTACT_ELEMENTS_2 13

/* status stuff */
typedef int Finley_Status_t;
#define Finley_increaseStatus(self) ((self)->status)++
#define FINLEY_INITIAL_STATUS 0

#ifdef PASO_MPI
#define FINLEY_INIT_ITEMSIZE (sizeof(double)*8)
#define FINLEY_NODE_TAG 0;
#define FINLEY_ELEMENT_TAG 10000;
extern int __g_nodeTag;
extern int __g_elementTag; 
#endif

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

#ifdef PASO_MPI
bool_t Finley_MPI_noError( Paso_MPIInfo *mpi_info );
#endif

#endif /* #ifndef INC_FINLEY */

