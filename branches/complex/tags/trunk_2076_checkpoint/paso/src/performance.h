
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: perfomance monitor interface using papi              */

/**************************************************************/

/* Copyrights by ACcESS Australia 2006 */
/* Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_PERFORMANCE
#define INC_PASO_PERFORMANCE

#define PERFORMANCE_UNMONITORED_EVENT -1
#define PERFORMANCE_NUM_EVENTS 10                         /* maximum number of events handeled by PAPI */

#define PERFORMANCE_ALL 0
#define PERFORMANCE_SOLVER 1
#define PERFORMANCE_PRECONDITIONER_INIT 2
#define PERFORMANCE_PRECONDITIONER 3
#define PERFORMANCE_MVM 4
#define PERFORMANCE_ASSEMBLAGE 5
#define PERFORMANCE_UNKNOWN 6                             /* more can be added here */
#define PERFORMANCE_NUM_MONITORS PERFORMANCE_UNKNOWN+1

#define PERFORMANCE_UNUSED -1
#define PERFORMANCE_CLOSED 0
#define PERFORMANCE_OPENED 1

#ifdef PAPI
#include <papi.h>
struct Paso_Performance {
    int event_set;                                                       /* papi event sets for the monitors */
    int num_events;                                                     /* number of events tracked by the monitors */
    int events[PERFORMANCE_NUM_EVENTS];                                 /* the event tracked by the monitors */
    long_long values[PERFORMANCE_NUM_MONITORS][PERFORMANCE_NUM_EVENTS]; /* counter accummulator */
    long_long cycles[PERFORMANCE_NUM_MONITORS];                         /* cycle accumulator */
    int set[PERFORMANCE_NUM_MONITORS];
};
#else
struct Paso_Performance {
    int none;
};
#endif
typedef struct Paso_Performance Paso_Performance;


void Performance_open(Paso_Performance* pp,int verbose);
int  Performance_getEventIndex(Paso_Performance* pp, int event_id);
void Performance_close(Paso_Performance* pp,int verbose);
void Performance_startMonitor(Paso_Performance* pp,int monitor);
void Performance_stopMonitor(Paso_Performance* pp,int monitor);

#endif

/*
 * $Log$
*/
