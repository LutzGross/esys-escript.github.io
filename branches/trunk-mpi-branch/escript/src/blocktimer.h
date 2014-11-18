
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

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
