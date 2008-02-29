#include <stdio.h>
#include <search.h>
#include "system_dep.h"

/* Enable the block timer (or remove this and use -DBLOCKTIMER) */
/* # define BLOCKTIMER */

# define NUM_TIMERS 1024

ESCRIPT_DLL_API
void blocktimer_initialize();
ESCRIPT_DLL_API
void blocktimer_increment(char *name, double start_time);
ESCRIPT_DLL_API
int blocktimer_getOrCreateTimerId(char *name);
ESCRIPT_DLL_API
void blocktimer_reportSortByName();
ESCRIPT_DLL_API
void blocktimer_reportSortByTime();
ESCRIPT_DLL_API
double blocktimer_time();

