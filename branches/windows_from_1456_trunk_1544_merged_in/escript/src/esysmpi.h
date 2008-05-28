
/* $Id: esysmpi.h 1306 2007-09-18 05:51:09Z ksteube $ */

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

#ifndef ESYSMPI_H
#define ESYSMPI_H

#ifdef PASO_MPI
#include "mpi.h"
#else 
#define MPI_Comm long
#endif

#endif 
