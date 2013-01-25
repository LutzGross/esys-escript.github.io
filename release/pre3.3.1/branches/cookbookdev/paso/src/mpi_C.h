
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


/*
    mpi_C.h

    Ensures that no C++ stuff leaks into Paso/Finley from mpi.h
*/

#ifdef PASO_MPI

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

#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#undef MPICH_SKIP_MPICXX
#undef OMPI_SKIP_MPICXX

#endif /* PASO_MPI_C */

#endif /* PASO_MPI */

