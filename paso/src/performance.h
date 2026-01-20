
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


/****************************************************************************/

/* Paso: perfomance monitor interface using PAPI                            */

/****************************************************************************/

/* Copyrights by ACcESS Australia 2006 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_PERFORMANCE_H__
#define __PASO_PERFORMANCE_H__

#ifdef ESYS_HAVE_PAPI
#include <papi.h>
#endif

namespace paso {

#define PERFORMANCE_UNMONITORED_EVENT -1
#define PERFORMANCE_NUM_EVENTS 10 // maximum number of events handled by PAPI

#define PERFORMANCE_ALL 0
#define PERFORMANCE_SOLVER 1
#define PERFORMANCE_PRECONDITIONER_INIT 2
#define PERFORMANCE_PRECONDITIONER 3
#define PERFORMANCE_MVM 4
#define PERFORMANCE_ASSEMBLAGE 5
#define PERFORMANCE_UNKNOWN 6  // more can be added here
#define PERFORMANCE_NUM_MONITORS PERFORMANCE_UNKNOWN+1

#define PERFORMANCE_UNUSED -1
#define PERFORMANCE_CLOSED 0
#define PERFORMANCE_OPENED 1

struct Performance
{
#ifdef ESYS_HAVE_PAPI
    /// PAPI event sets for the monitors
    int event_set;
    /// number of events tracked by the monitors
    int num_events;
    /// the events tracked by the monitors
    int events[PERFORMANCE_NUM_EVENTS];
    /// counter accumulator
    long_long values[PERFORMANCE_NUM_MONITORS][PERFORMANCE_NUM_EVENTS];
    /// cycle accumulator
    long_long cycles[PERFORMANCE_NUM_MONITORS];
    int set[PERFORMANCE_NUM_MONITORS];
#else
    int dummy;
#endif
};

void Performance_open(Performance* pp, int verbose);
int  Performance_getEventIndex(Performance* pp, int event_id);
void Performance_close(Performance* pp, int verbose);
void Performance_startMonitor(Performance* pp, int monitor);
void Performance_stopMonitor(Performance* pp, int monitor);

} // namespace paso

#endif // __PASO_PERFORMANCE_H__

