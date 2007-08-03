
#if 0

  A simple library to make it easy to use PAPI to time a single block of code.

  You can start and stop the timers as often as you like and the results will
  be accumulated.

  If you use the event PAPI_FP_OPS or NATV_FP_OPS_RETIRED your MFLOPS will be
  computed based automatically.

  Usage example:

      #include "blockpapi.h"

      /* Add a PAPI preset event */
      blockpapi_addEvent(PAPI_FP_OPS,			"total floating-point operations");

      /* Add a couple native events */
      blockpapi_addEvent(NATV_CPU_CYCLES,		"total cpu cycles");
      blockpapi_addEvent(NATV_BACK_END_BUBBLE_ALL,	"cycles stalled for any reason");

      blockpapi_start();
      /* Compute something here */
      blockpapi_stop();

      blockpapi_writeReport();
      blockpapi_writeSystemInfo();

  Compile with cc -DBLOCKPAPI file.c blockpapi.c -lpapi

#endif




#include "blockpapi.h"

#ifdef BLOCKPAPI
#include "papi.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_COUNTERS 100

int		 g_numCounters=0;
int		 g_events[MAX_COUNTERS];
long_long_t	 g_values[MAX_COUNTERS];
char		*g_descriptions[MAX_COUNTERS];
int		 g_numPasses = 0;
double		 g_totalMicroSecs = 0.0;
long_long_t	 g_numFpOps = -1;

void blockpapi_addEvent(int event, char *description) {
#ifdef BLOCKPAPI
  if (g_numCounters >= MAX_COUNTERS) {
    fprintf(stderr, "blockpapi_addEvent: too many counters, cannot add '%s'\n", description);
    exit(1);
  }
  if (g_numCounters > PAPI_num_counters()) {
    fprintf(stderr, "blockpapi_addEvent: CPU does not support this many counters, cannot add '%s'\n", description);
    exit(1);
  }
  g_values[g_numCounters] = (long_long_t) 0;
  g_events[g_numCounters] = event;
  g_descriptions[g_numCounters] = strdup(description);
  g_numCounters++;
#endif
}

void blockpapi_start() {
#ifdef BLOCKPAPI
  int retval;
  if (g_numCounters == 0) { return; }
  g_numPasses++;
  g_totalMicroSecs -= PAPI_get_real_usec();
  if ((retval = PAPI_start_counters(g_events, g_numCounters)) != PAPI_OK) {
    fprintf(stderr, "blockpapi_start: PAPI_start_counters failed\n");
    fprintf(stderr, "PAPI error %d: %s\n",retval,PAPI_strerror(retval));
    exit(1);
  }
#endif
}

void blockpapi_stop() {
#ifdef BLOCKPAPI
  long_long_t g_values_tmp[MAX_COUNTERS];
  if (g_numCounters == 0) { return; }
  g_totalMicroSecs += PAPI_get_real_usec();
  /* Add the the latest counts to the running total */
  if ( (PAPI_accum_counters(g_values, g_numCounters)) != PAPI_OK) {
    fprintf(stderr, "blockpapi_stop: PAPI_accum_counters failed");
    exit(1);
  }
  /* Stop the counters, toss the values since we accumulated them above */
  if ((PAPI_stop_counters(g_values_tmp, g_numCounters)) != PAPI_OK) {
    fprintf(stderr, "blockpapi_stop: PAPI_stop_counters failed");	
    exit(1);
  }
#endif
}

void blockpapi_writeReport() {
#ifdef BLOCKPAPI
  int i;
  if (g_numCounters == 0) { return; }
  printf("Blockpapi report: this block was executed %5d time(s) in %10.2lf seconds\n", g_numPasses, (double) g_totalMicroSecs/1000000.0);
  for (i=0; i<g_numCounters; i++) {
    /* blockpapi_writeReportLine(g_values[i], g_events[i], g_descriptions[i]); */
    char name[200];
    PAPI_event_code_to_name(g_events[i], name);
    printf("%15ld  %-25s   %s\n", g_values[i], name, g_descriptions[i]);
    /* If num FP ops available we'll save it to compute FLOPS */
    if (!strcmp(name, "PAPI_FP_OPS") || !strcmp(name, "FP_OPS_RETIRED")) {
      g_numFpOps = g_values[i];
    }
  }
  if (g_numFpOps > 0 && g_totalMicroSecs > 0.0) {
    printf("%15.2f  %-25s   %s\n", g_numFpOps/g_totalMicroSecs, "MFLOPS", "millions of floating point operations per second");
  }
#endif
}

long_long_t *blockpapi_getValues() {
  return(g_values);
}

void blockpapi_writeSystemInfo() {
#ifdef BLOCKPAPI
  int num_hwcntrs = 0;
  const PAPI_hw_info_t *hwinfo = NULL;
  if ((hwinfo = PAPI_get_hardware_info ()) != NULL) {
    printf ("This system has %d %s CPUs running at %.0f MHz.\n", hwinfo->totalcpus,
      hwinfo->model_string, hwinfo->mhz);
    printf ("They are revision %f ", hwinfo->revision);
  }
  if ((num_hwcntrs = PAPI_num_counters ()) < PAPI_OK) {
    printf ("and have no performance counters.\n");
    exit (1);
  }
  else {
    printf ("and have %d performance counters.\n", num_hwcntrs);
  }
#endif
}

