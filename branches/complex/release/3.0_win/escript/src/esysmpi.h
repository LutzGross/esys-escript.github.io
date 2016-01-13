
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


#ifndef ESYSMPI_H
#define ESYSMPI_H

#ifdef PASO_MPI
#include "mpi.h"
#else 
#define MPI_Comm long
#endif

#endif 
