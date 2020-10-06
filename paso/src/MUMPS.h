
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
#include "Options.h"
#include "PasoException.h"

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
#undef NOMINMAX
#endif
#endif // ESYS_HAVE_MUMPS

namespace paso {

struct MUMPS_Handler_t {
    bool verbose;
    std::stringstream ssExceptMsg;
#ifdef ESYS_HAVE_MUMPS
    MUMPS_INT myid;
#ifdef _WIN32 // workaround for d/zmumps dll clash
    HINSTANCE h_mumps_c_dll;
#endif
#endif // ESYS_HAVE_MUMPS
};

template <typename T>
struct MUMPS_Handler : MUMPS_Handler_t {
    T* rhs;
};

template <typename T>
void MUMPS_free(SparseMatrix<T>* A);

template <typename T>
void MUMPS_solve(SparseMatrix_ptr<T> A, T* out, T* in, dim_t numRefinements, bool verbose);

template <typename T>
void MUMPS_print_list(const char* name, const T* vals, const int n, const int max_n=100);

std::ostream& operator<<(std::ostream& os, const cplx_t& c);

template <>
struct MUMPS_Handler<double> : MUMPS_Handler_t {
    double* rhs;
#ifdef ESYS_HAVE_MUMPS
    DMUMPS_STRUC_C id;
    typedef double mumps_t;
#ifdef _WIN32 // workaround for d/zmumps dll clash
    typedef HRESULT (CALLBACK* MUMPS_C_FUNC_PTR)(DMUMPS_STRUC_C*);
    MUMPS_C_FUNC_PTR mumps_c;
    const char* mumps_lib = "libdmumps";
    const char* mumps_proc = "dmumps_c";
#else
    void (*mumps_c)(DMUMPS_STRUC_C*) = &dmumps_c;
#endif
#endif // ESYS_HAVE_MUMPS
};

template <>
struct MUMPS_Handler<cplx_t> : MUMPS_Handler_t {
    cplx_t* rhs;
#ifdef ESYS_HAVE_MUMPS
    ZMUMPS_STRUC_C id;
    typedef ZMUMPS_COMPLEX mumps_t;
#ifdef _WIN32 // workaround for dmumps/zdmumps dll clash
    typedef HRESULT (CALLBACK* MUMPS_C_FUNC_PTR)(ZMUMPS_STRUC_C*);
    MUMPS_C_FUNC_PTR mumps_c;
    const char* mumps_lib = "libzmumps";
    const char* mumps_proc = "zmumps_c";
#else
    void (*mumps_c)(ZMUMPS_STRUC_C*) = &zmumps_c;
#endif
#endif // ESYS_HAVE_MUMPS
};

/// frees any MUMPS related data from the matrix
template <typename T>
void MUMPS_free(SparseMatrix<T>* A)
{
    if (A && A->solver_p) {
#ifdef ESYS_HAVE_MUMPS
        // Terminate instance.
        auto pt = static_cast<MUMPS_Handler<T>*>(A->solver_p);
        delete[] pt->rhs;
        pt->id.job = MUMPS_JOB_END;
        pt->mumps_c(&pt->id);
#ifdef _WIN32
        FreeLibrary(pt->h_mumps_c_dll);
#endif
        if (pt->myid == 0) {
            std::string message = pt->ssExceptMsg.str();
            if (!message.empty()) {
                // terminating with solve error message
                throw PasoException(message);
            }
        }
        MUMPS_INT ierr = MPI_Finalize();
        if (pt->verbose) {
            std::cout << "MUMPS: instance terminated." << std::endl;
        }
        delete pt;
#endif
        A->solver_p = NULL;
    }
}

/// calls the solver
template <typename T>
void MUMPS_solve(SparseMatrix_ptr<T> A, T* out, T* in, dim_t numRefinements, bool verbose)
{
#ifdef ESYS_HAVE_MUMPS
    if (! (A->type & (MATRIX_FORMAT_OFFSET1 + MATRIX_FORMAT_BLK1)) ) {
        throw PasoException("Paso: MUMPS requires CSR format with index offset 1 and block size 1.");
    }

    auto pt = reinterpret_cast<MUMPS_Handler<T>*>(A->solver_p);
    if (pt == NULL) {
        pt = new MUMPS_Handler<T>;
#ifdef _WIN32
        pt->h_mumps_c_dll = LoadLibrary(pt->mumps_lib);
        if (pt->h_mumps_c_dll == NULL) {
            std::stringstream ss;
            ss << "Paso: MUMPS LoadLibrary failed - \"" << pt->mumps_lib << "\".";
            throw PasoException(ss.str());
        }
        pt->mumps_c = (MUMPS_Handler<T>::MUMPS_C_FUNC_PTR)GetProcAddress(pt->h_mumps_c_dll, pt->mumps_proc);
        if (pt->mumps_c == NULL) {
            std::stringstream ss;
            ss << "Paso: MUMPS GetProcAddress failed - \"" << pt->mumps_proc << "\".";
            throw PasoException(ss.str());
        }
#endif
        A->solver_p = (void*) pt;
        A->solver_package = PASO_MUMPS;
        double time0 = escript::gettime();

        A->pattern->csrToHB(); // generate Harwell-Boeing format needed for MUMPS from CSR
        MUMPS_INT n = A->numRows;  // matrix order
        MUMPS_INT8 nnz = A->pattern->len;  // number non-zeros
        MUMPS_INT* irn = reinterpret_cast<MUMPS_INT*>(A->pattern->hb_row);  // row indices array
        MUMPS_INT* jcn = reinterpret_cast<MUMPS_INT*>(A->pattern->hb_col);  // col indices array
        pt->verbose = verbose;
        if (pt->verbose) {
            std::cout << "MUMPS in  ===>" << std::endl;
            std::cout << "n = " << n << std::endl;
            std::cout << "nnz = " << nnz << std::endl;
            MUMPS_print_list("val", A->val, nnz);
            MUMPS_print_list("in", in, n);
            MUMPS_print_list("ptr", A->pattern->ptr, n+1);
            MUMPS_print_list("index", A->pattern->index, nnz);
            MUMPS_print_list("hb_row", A->pattern->hb_row, nnz);
            MUMPS_print_list("hb_col", A->pattern->hb_col, nnz);
        }
        pt->rhs = new T[n];
        std::memcpy(pt->rhs, in, n*sizeof(T));
        MUMPS_INT ierr;
        ierr = MPI_Init(NULL, NULL);
        ierr = MPI_Comm_rank(MPI_COMM_WORLD, &pt->myid);

        // Initialize a MUMPS instance. Use MPI_COMM_WORLD
        pt->id.comm_fortran = MUMPS_USE_COMM_WORLD;
        pt->id.par = 1; pt->id.sym = 0;
        pt->id.job = MUMPS_JOB_INIT;
        pt->mumps_c(&pt->id);
        // Define the problem on the host
        if (pt->myid == 0) {
            pt->id.n = n; pt->id.nnz = nnz;
            pt->id.irn = irn; pt->id.jcn = jcn;
            pt->id.a = reinterpret_cast<typename MUMPS_Handler<T>::mumps_t*>(A->val);
            pt->id.rhs = reinterpret_cast<typename MUMPS_Handler<T>::mumps_t*>(pt->rhs);
        }
        if (!pt->verbose) {
            // No outputs
            pt->id.ICNTL(1)=-1; pt->id.ICNTL(2)=-1; pt->id.ICNTL(3)=-1; pt->id.ICNTL(4)=0;
        }

        // Call the MUMPS package (analyse, factorization and solve).
        pt->id.job = 6;
        pt->mumps_c(&pt->id);
        if (pt->id.infog[0] < 0) {
            pt->ssExceptMsg << "(PROC " << pt->myid << ") MUMPS ERROR: INFOG(1)=" << pt->id.infog[0]
                << ", INFOG(2)=" << pt->id.infog[1];
        } else {
            std::memcpy(out, reinterpret_cast<T*>(pt->rhs), n*sizeof(T));
            if (pt->id.infog[0] > 0) {
                std::cout << "(PROC " << pt->myid << ") MUMPS WARNING: INFOG(1)=" << pt->id.infog[0]
                    << ", INFOG(2)=" << pt->id.infog[1];
            }
            if (pt->verbose) {
                std::cout << "MUMPS out ===>" << std::endl;
                MUMPS_print_list("out", out, n);
                std::cout << "MUMPS: factorization and solve completed (time = "
                    << escript::gettime()-time0 << ")." << std::endl;
            }
        }
    }
#else // ESYS_HAVE_MUMPS
    throw PasoException("Paso: Not compiled with MUMPS.");
#endif
}

// output array data for debugging solver
// array length limit is 100 by default, use 0 for no limit
template <typename T>
void MUMPS_print_list(const char* name, const T* vals, const int n, const int max_n)
{
    std::cout << name << " = [ ";
    for (int i=0; i<n; i++) {
        if (i > 0) {
            std::cout << ", ";
        }
        std::cout << vals[i];
        if (max_n > 0) {
            if (i > max_n) {
                std::cout << ", ...";
                break;
            }
        }
    }
    std::cout << " ]" << std::endl;
}

} // namespace paso

#endif // __PASO_MUMPS_H__

