// $Id$
/*
 ************************************************************
 *          Copyright 2006,2007 by ACcESS MNRF              *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

#ifndef ESYSMPI_H
#define ESYSMPI_H

#ifdef PASO_MPI
#include "mpi.h"
#else 
#define MPI_Comm long
#endif

#endif 
