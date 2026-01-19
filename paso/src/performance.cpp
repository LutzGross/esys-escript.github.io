
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/


/****************************************************************************/

/* Paso: performance monitor interface using PAPI                           */

/****************************************************************************/

/* Copyrights by ACcESS Australia 2006 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "PasoException.h"
#include "performance.h"

namespace paso {

/// sets up the monitoring process
void Performance_open(Performance* pp, int verbose)
{
#ifdef ESYS_HAVE_PAPI
    #pragma omp single
    {
        pp->event_set = PAPI_NULL;
        // Initialize the PAPI library
        int retval = PAPI_library_init(PAPI_VER_CURRENT);
        if (retval != PAPI_VER_CURRENT && retval > 0) {
            throw PasoException("performance: PAPI library version mismatch.");
        } else if (retval < 0) {
            throw PasoException("performance: PAPI initialization error.");
        } else {
            if (PAPI_create_eventset(&(pp->event_set)) != PAPI_OK)
                throw PasoException("performance: PAPI event set up failed.");
        }
        // try to add various monitors
        pp->num_events=0;
        if (PAPI_add_event(pp->event_set, PAPI_FP_OPS) == PAPI_OK) {
            pp->events[pp->num_events]=PAPI_FP_OPS;
            pp->num_events++;
        }
        if (PAPI_add_event(pp->event_set, PAPI_L1_DCM) == PAPI_OK) {
            pp->events[pp->num_events]=PAPI_L1_DCM;
            pp->num_events++;
        }
        if (PAPI_add_event(pp->event_set, PAPI_L2_DCM) == PAPI_OK) {
            pp->events[pp->num_events]=PAPI_L2_DCM;
            pp->num_events++;
        }
        if (PAPI_add_event(pp->event_set, PAPI_L3_DCM) == PAPI_OK) {
            pp->events[pp->num_events]=PAPI_L3_DCM;
            pp->num_events++;
        }
        for (int i=0; i<PERFORMANCE_NUM_MONITORS; ++i) {
            pp->cycles[i] = 0;
            pp->set[i] = PERFORMANCE_UNUSED;
            for (int j=0; j<PERFORMANCE_NUM_EVENTS; ++j)
                pp->values[i][j] = 0.;
        }
        PAPI_start(pp->event_set);
    } // omp single
#endif // ESYS_HAVE_PAPI
}

/// find the index of an event in the list of monitored events
int Performance_getEventIndex(Performance* pp, int event_id)
{
#ifdef ESYS_HAVE_PAPI
    for (int i=0; i<pp->num_events; ++i)
        if (pp->events[i]==event_id)
            return i;
#endif
    return PERFORMANCE_UNMONITORED_EVENT;
}

/// shuts down the monitoring process
void Performance_close(Performance* pp, int verbose)
{
#ifdef ESYS_HAVE_PAPI
#pragma omp single
    {
        if (verbose) {
            int i_ops = Performance_getEventIndex(pp, PAPI_FP_OPS);
            int i_l1_miss = Performance_getEventIndex(pp, PAPI_L1_DCM);
            int i_l2_miss = Performance_getEventIndex(pp, PAPI_L2_DCM);
            int i_l3_miss = Performance_getEventIndex(pp, PAPI_L3_DCM);
            printf(" monitor               |");
            if (i_ops != PERFORMANCE_UNMONITORED_EVENT)
                printf("  flops/cycle |");
            if (i_ops != PERFORMANCE_UNMONITORED_EVENT && i_l1_miss!=PERFORMANCE_UNMONITORED_EVENT)
                printf(" L1 miss/flops |");
            if (i_ops!=PERFORMANCE_UNMONITORED_EVENT && i_l2_miss!=PERFORMANCE_UNMONITORED_EVENT)
                printf(" L2 miss/flops |");
            if (i_ops!=PERFORMANCE_UNMONITORED_EVENT &&  i_l3_miss!=PERFORMANCE_UNMONITORED_EVENT)
                printf(" L3 miss/flops |");
            printf("\n");
            for (int i=0; i<PERFORMANCE_NUM_MONITORS; ++i) {
                if (pp->set[i] != PERFORMANCE_UNUSED) {
                    switch(i) {
                        case PERFORMANCE_ALL:
                            printf(" over all              |");
                        break;
                        case PERFORMANCE_SOLVER:
                            printf(" solver                |");
                        break;
                        case PERFORMANCE_PRECONDITIONER_INIT:
                            printf(" init. preconditioner  |");
                        break;
                        case PERFORMANCE_PRECONDITIONER:
                            printf(" preconditioner        |");
                        break;
                        case PERFORMANCE_MVM:
                            printf(" matrix-vector product |");
                        break;
                        case PERFORMANCE_ASSEMBLAGE:
                            printf(" assemblage            |");
                        break;
                        default:
                            printf(" no name               |");
                        break;
                    }
                    if (pp->set[i]==PERFORMANCE_CLOSED) {
                        if (i_ops != PERFORMANCE_UNMONITORED_EVENT)
                            printf(" %12.5e |",(((double)(pp->values[i][i_ops]))/((double)(pp->cycles[i]))));
                        if (i_ops!=PERFORMANCE_UNMONITORED_EVENT && i_l1_miss!=PERFORMANCE_UNMONITORED_EVENT)
                            printf(" %13.5e |",(((double)(pp->values[i][i_l1_miss]))/((double)(pp->values[i][i_ops]))));
                        if (i_ops!=PERFORMANCE_UNMONITORED_EVENT && i_l2_miss!=PERFORMANCE_UNMONITORED_EVENT)
                            printf(" %13.5e |",(((double)(pp->values[i][i_l2_miss]))/((double)(pp->values[i][i_ops]))));
                        if (i_ops!=PERFORMANCE_UNMONITORED_EVENT && i_l3_miss!=PERFORMANCE_UNMONITORED_EVENT)
                            printf(" %13.5e |",(((double)(pp->values[i][i_l3_miss]))/((double)(pp->values[i][i_ops]))));
                    } else {
                        printf("not closed!!!");
                    }
                    printf("\n");
                }
            }
        }
        long_long values[PERFORMANCE_NUM_EVENTS];
        PAPI_stop(pp->event_set, values);
        PAPI_cleanup_eventset(pp->event_set);
        PAPI_destroy_eventset(&pp->event_set);
    }
#endif
}

/// switches on a monitor
void Performance_startMonitor(Performance* pp, int monitor)
{
#ifdef ESYS_HAVE_PAPI
#pragma omp barrier
#pragma omp single
    {
        long_long values[PERFORMANCE_NUM_EVENTS];
        // Start counting events in the Event Set
        PAPI_read(pp->event_set, values);
        for (int i=0; i<pp->num_events; ++i)
            pp->values[monitor][i] -= values[i];
        // set cycles
        pp->cycles[monitor] -= PAPI_get_real_cyc();
        pp->set[monitor] = PERFORMANCE_OPENED;
    }
#endif
}

/// switches off a monitor
void Performance_stopMonitor(Performance* pp,int monitor)
{
#ifdef ESYS_HAVE_PAPI
#pragma omp barrier
#pragma omp single
    {
        long_long values[PERFORMANCE_NUM_EVENTS];
        // Add the counters in the Event Set
        PAPI_read(pp->event_set, values);
        for (int i=0; i<pp->num_events; ++i)
            pp->values[monitor][i] += values[i];
        // set cycles
        pp->cycles[monitor] += PAPI_get_real_cyc();
        pp->set[monitor] = PERFORMANCE_CLOSED;
    }
#endif
}

} // namespace paso

