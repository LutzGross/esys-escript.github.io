
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/****************************************************************************/

/* Paso: interface to the MUMPS library */

/****************************************************************************/

#ifndef __PASO_MUMPS_H__
#define __PASO_MUMPS_H__

#include "SparseMatrix.h"

#ifdef ESYS_HAVE_MUMPS
// TODO: is this needed? #pragma push_macro("MPI_COMM_WORLD")
#if defined(MPI_COMM_WORLD)
#undef MPI_COMM_WORLD    // breaks mumps_mpi.h, defined in escriptcore/src/EsysMPI.h
#endif
#include <mumps_mpi.h>
// TODO: is this needed? #pragma pop_macro("MPI_COMM_WORLD")
// #include <zmumps_c.h>
#include <dmumps_c.h>
#include <zmumps_c.h>
#define MUMPS_JOB_INIT -1
#define MUMPS_JOB_END -2
#define MUMPS_USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1] // macro s.t. indices match documentation

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
typedef HRESULT (CALLBACK* DMUMPS_C_FUNC_PTR)(DMUMPS_STRUC_C*);
typedef HRESULT (CALLBACK* ZMUMPS_C_FUNC_PTR)(ZMUMPS_STRUC_C*);
#undef NOMINMAX
#endif
#endif // ESYS_HAVE_MUMPS

namespace paso {

struct MUMPS_Handler {
    bool isComplex;
    bool verbose;
    double* rhs;
    cplx_t* crhs;
    std::stringstream ssExceptMsg;
#ifdef ESYS_HAVE_MUMPS
    MUMPS_INT myid;
    DMUMPS_STRUC_C id;
    ZMUMPS_STRUC_C zid;
#ifdef _WIN32 // workaround for dmumps/zdmumps dll clash
    HINSTANCE h_dmumps_c_dll;
    DMUMPS_C_FUNC_PTR dmumps_c;
    HINSTANCE h_zmumps_c_dll;
    ZMUMPS_C_FUNC_PTR zmumps_c;
#endif
#endif // ESYS_HAVE_MUMPS
};

void MUMPS_free(SparseMatrix* A);
void MUMPS_solve(SparseMatrix_ptr A, double* out, double* in,
                   dim_t numRefinements, bool verbose);
void MUMPS_solve(SparseMatrix_ptr A, cplx_t* out, cplx_t* in,
                   dim_t numRefinements, bool verbose);

} // namespace paso

#endif // __PASO_MUMPS_H__

