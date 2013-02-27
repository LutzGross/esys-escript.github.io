
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/*
    mpi_C.h

    Ensures that no C++ stuff leaks into Paso/Finley from mpi.h
*/

#ifdef ESYS_MPI

/*
#ifndef MPI_NO_CPPBIND
  #define MPI_NO_CPPBIND
  #include <mpi.h>
  #undef MPI_NO_CPPBIND
#else
  #include <mpi.h>
#endif
*/
#ifndef ESYS_MPI_C
#define ESYS_MPI_C

#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#undef MPICH_SKIP_MPICXX
#undef OMPI_SKIP_MPICXX

#endif /* ESYS_MPI_C */

#endif /* ESYS_MPI */

