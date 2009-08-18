
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
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
/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "performance.h"

/*                                                            */
/*  sets up the the monitoring process                        */
/*                                                            */
void Performance_open(Paso_Performance* pp,int verbose) {
   #ifdef PAPI
      int i,j;
      #pragma omp single
      {
         pp->event_set=PAPI_NULL;
         /* Initialize the PAPI library */
         int retval = PAPI_library_init(PAPI_VER_CURRENT);
         if (retval != PAPI_VER_CURRENT && retval > 0) {
           Paso_setError(SYSTEM_ERROR,"Paso_performance: PAPI library version mismatch.");
         } else if (retval < 0) {
           Paso_setError(SYSTEM_ERROR,"Paso_performance: PAPI initialization error.");
         } else {
            if (PAPI_create_eventset(&(pp->event_set)) != PAPI_OK) 
               Paso_setError(SYSTEM_ERROR,"Paso_performance: PAPI event set set up failed.");
         }
         if (Paso_noError()) {
            /* try to add various monitors */
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
            for (i=0;i<PERFORMANCE_NUM_MONITORS;++i) {
               pp->cycles[i]=0;
               pp->set[i]=PERFORMANCE_UNUSED;
               for (j=0;j<PERFORMANCE_NUM_EVENTS;++j) pp->values[i][j]=0.;
            }
            PAPI_start(pp->event_set);
         }
      }
   #endif
}
/* find the index of an event in the list of monitored events */
int  Performance_getEventIndex(Paso_Performance* pp, int event_id) {
   #ifdef PAPI
     int i;
     for (i=0;i<pp->num_events;++i) 
        if (pp->events[i]==event_id) return i;
   #endif
   return PERFORMANCE_UNMONITORED_EVENT;
}
/*                                                            */
/*  shuts down the monitoring process                         */
/*                                                            */
void Performance_close(Paso_Performance* pp,int verbose) {
    #ifdef PAPI
      long_long values[PERFORMANCE_NUM_EVENTS];
      #pragma omp single
      {
        if (Paso_noError() && verbose) {
           int i;
           int i_ops=Performance_getEventIndex(pp,PAPI_FP_OPS);
           int i_l1_miss=Performance_getEventIndex(pp,PAPI_L1_DCM);
           int i_l2_miss=Performance_getEventIndex(pp,PAPI_L2_DCM);
           int i_l3_miss=Performance_getEventIndex(pp,PAPI_L3_DCM);
           printf(" monitor               |");
           if (i_ops!=PERFORMANCE_UNMONITORED_EVENT) printf("  flops/cycle |");
           if (i_ops!=PERFORMANCE_UNMONITORED_EVENT &&  i_l1_miss!=PERFORMANCE_UNMONITORED_EVENT) printf(" L1 miss/flops |");
           if (i_ops!=PERFORMANCE_UNMONITORED_EVENT &&  i_l2_miss!=PERFORMANCE_UNMONITORED_EVENT) printf(" L2 miss/flops |");
           if (i_ops!=PERFORMANCE_UNMONITORED_EVENT &&  i_l3_miss!=PERFORMANCE_UNMONITORED_EVENT) printf(" L3 miss/flops |");
           printf("\n");
           for (i=0;i<PERFORMANCE_NUM_MONITORS;++i) {
               if (pp->set[i]!=PERFORMANCE_UNUSED ) {
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
                  if (pp->set[i]==PERFORMANCE_CLOSED ) {
                     if (i_ops!=PERFORMANCE_UNMONITORED_EVENT) printf(" %12.5e |",(((double)(pp->values[i][i_ops]))/((double)(pp->cycles[i]))));
                     if (i_ops!=PERFORMANCE_UNMONITORED_EVENT &&  i_l1_miss!=PERFORMANCE_UNMONITORED_EVENT) 
                                                   printf(" %13.5e |",(((double)(pp->values[i][i_l1_miss]))/((double)(pp->values[i][i_ops]))));
                     if (i_ops!=PERFORMANCE_UNMONITORED_EVENT &&  i_l2_miss!=PERFORMANCE_UNMONITORED_EVENT) 
                                                   printf(" %13.5e |",(((double)(pp->values[i][i_l2_miss]))/((double)(pp->values[i][i_ops]))));
                     if (i_ops!=PERFORMANCE_UNMONITORED_EVENT &&  i_l3_miss!=PERFORMANCE_UNMONITORED_EVENT)
                                                   printf(" %13.5e |",(((double)(pp->values[i][i_l3_miss]))/((double)(pp->values[i][i_ops]))));
                  } else {
                     printf("not closed!!!");
                  }
                  printf("\n");
             }
          }
        }
        PAPI_stop(pp->event_set,values);
        PAPI_cleanup_eventset(pp->event_set);
        PAPI_destroy_eventset(&(pp->event_set));
      }
    #endif
}
/*                                                                     */
/*    switches on a monitor                                            */
/*                                                                     */
void Performance_startMonitor(Paso_Performance* pp,int monitor) {
    #ifdef PAPI
       int i;
       long_long values[PERFORMANCE_NUM_EVENTS];
       #pragma omp barrier
       #pragma omp single
       {
          /* Start counting events in the Event Set */
          PAPI_read(pp->event_set,values);
          for (i=0;i<pp->num_events;++i) pp->values[monitor][i]-=values[i];
          /* set cycles */
          pp->cycles[monitor]-=PAPI_get_real_cyc();
          pp->set[monitor]=PERFORMANCE_OPENED;
       }
    #endif
}
/*                                                                     */
/*    switches off a monitor                                           */
/*                                                                     */
void Performance_stopMonitor(Paso_Performance* pp,int monitor) {
    #ifdef PAPI
       int i;
       long_long values[PERFORMANCE_NUM_EVENTS];
       #pragma omp barrier
       #pragma omp single
       {
          /* Add the counters in the Event Set */
          PAPI_read(pp->event_set,values);
          for (i=0;i<pp->num_events;++i) pp->values[monitor][i]+=values[i];
          /* set cycles */
          pp->cycles[monitor]+=PAPI_get_real_cyc();
          pp->set[monitor]=PERFORMANCE_CLOSED;
       }
    #endif
}
