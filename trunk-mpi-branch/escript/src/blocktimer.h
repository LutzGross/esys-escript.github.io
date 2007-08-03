#include <stdio.h>
#include <search.h>

/* Enable the block timer (or remove this and use -DBLOCKTIMER) */
/* # define BLOCKTIMER */

# define NUM_TIMERS 1024

void blocktimer_initialize();
void blocktimer_increment(char *name, double start_time);
int blocktimer_getOrCreateTimerId(char *name);
void blocktimer_reportSortByName();
void blocktimer_reportSortByTime();
void blocktimer_reportSystemInfo();
double blocktimer_time();

