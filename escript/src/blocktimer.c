
/* TODO: note if has been initialized, give error only once if not */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <search.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "escript/blocktimer.h"

#ifdef PASO_MPI
#include "mpi.h"
#endif

#ifdef BLOCKTIMER
static char  *g_names[NUM_TIMERS];	/* Names of the timers */
static int    g_count[NUM_TIMERS];	/* How many times was the timer incremented? */
static double g_times[NUM_TIMERS];	/* The total time spent in the block */
static double g_start_time;		/* Start time for the entire program */
static int    g_initialized = 0;	/* Has the blocktimer been initialized? */
static int    g_end_computed = 0;	/* Has the end time been set? */
#endif /* BLOCKTIMER */

void
blocktimer_initialize()
{
#ifdef BLOCKTIMER
  int i;

  for (i=0; i<NUM_TIMERS; i++) {
    g_names[i] = (char *)NULL;
    g_times[i] = 0.0;
    g_count[i] = 0;
  }

  if (hcreate(NUM_TIMERS) == 0) {
    perror("hcreate");
    fprintf(stderr, "blocktimer_initialize: Could not initialize hash table\n");
    exit(1);
  }

  g_initialized = 1;

  g_start_time = blocktimer_time();

  /* Initialize timer for "entire program" to zero so it appears first in the report */
  blocktimer_increment("entire program", g_start_time);
  g_count[0] = 0; /* Reset counter for "entire program" to zero */
#endif /* BLOCKTIMER */
}

void
blocktimer_increment(char *name, double start_time)
{
#ifdef BLOCKTIMER
  int id;

  if (!g_initialized) { return; }

  id = blocktimer_getOrCreateTimerId(name);

  g_times[id] += blocktimer_time() - start_time;
  g_count[id] += 1;
#endif /* BLOCKTIMER */
}

int
blocktimer_getOrCreateTimerId(char *name)
{
  int id=0;
#ifdef BLOCKTIMER
  static int nextId = 0;		/* Next timer ID to assign */
  ENTRY item, *found_item;

  if (!g_initialized) { return(0); }

  /* Has a timer with 'name' already been defined? */
  item.key = name;
  item.data = (void *) NULL;
  found_item = hsearch(item, FIND);

  if (found_item != NULL) {	/* Already defined so retrieve it from the hash */
    /* Return the ID of the entry we found */
    int *idTmp = found_item->data;
    id = *idTmp;
  }
  else {			/* Not already defined so create one */
    /* malloc new int, can't use stack var or all items share same data */
    int *idTmp = (int *)malloc(sizeof(int));
    /* Enter the new name in the hash */
    if (nextId >= NUM_TIMERS) {
      fprintf(stderr, "blocktimer: exceeded limit of %d timers, increase NUM_TIMERS\n", NUM_TIMERS);
      exit(1);
    }
    *idTmp = nextId++;
    item.key = name;
    item.data = (void *) idTmp;
    hsearch(item, ENTER);
    id = *idTmp;
    strncpy(g_names[id], name, strlen(name));
  }

#endif /* BLOCKTIMER */
  return(id);
}

void
blocktimer_reportSortByName()
{
#ifdef BLOCKTIMER
  int i;

  if (!g_initialized) { return; }

  if (!g_end_computed) {
    blocktimer_increment("entire program", g_start_time);
    g_end_computed = 1;
  }
  printf("BlockTimer sorted by name (sorting TBD):\n");
  for(i=0; i<NUM_TIMERS; i++) {
    if (g_names[i] != (char *) NULL) {
      printf("	%7d %15.2lf   %s\n", g_count[i], g_times[i], g_names[i]);
    }
  }
#endif /* BLOCKTIMER */
}

void
blocktimer_reportSortByTime()
{
#ifdef BLOCKTIMER
  int i;

  if (!g_initialized) { return; }

  if (!g_end_computed) {
    blocktimer_increment("entire program", g_start_time);
    g_end_computed = 1;
  }
  printf("BlockTimer sorted by time (sorting TBD):\n");
  for(i=0; i<NUM_TIMERS; i++) {
    if (g_names[i] != (char *) NULL) {
      printf("	%7d %15.2lf seconds for %s\n", g_count[i], g_times[i], g_names[i]);
    }
  }
#endif /* BLOCKTIMER */
}

/* Copied from Paso_timer() */
double
blocktimer_time()
{
  double out=0.0;
#ifdef PASO_MPI
  out = MPI_Wtime();
#else
#ifdef _OPENMP
  out=omp_get_wtime();
#else
  out=((double) clock())/CLOCKS_PER_SEC;
#endif
#endif
  return(out);
}

