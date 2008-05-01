
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

/*
    mpi_C.h

    Ensures that mpi++.h no C++ stuff leaks into Paso/Finley
*/


/*
#ifndef MPI_NO_CPPBIND
  #define MPI_NO_CPPBIND
  #include <mpi.h>
  #undef MPI_NO_CPPBIND
#else
  #include <mpi.h>
#endif
*/
#ifndef PASO_MPI_C
#define PASO_MPI_C
  #ifdef PASO_MPI
     #define MPICH_SKIP_MPICXX
     #include <mpi.h>
     #undef MPICH_SKIP_MPICXX
  #endif
#endif
