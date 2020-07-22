
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

#include "Paso.h"
#include "MUMPS.h"
#include "Options.h"
#include "PasoException.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>

namespace paso {

/// frees any MUMPS related data from the matrix
void MUMPS_free(SparseMatrix* A)
{
    if (A && A->solver_p) {
        MUMPS_Handler* pt = reinterpret_cast<MUMPS_Handler*>(A->solver_p);
#ifdef ESYS_HAVE_MUMPS
        // Terminate instance.
        if (pt->isComplex) {
            delete[] pt->crhs;
            pt->zid.job = MUMPS_JOB_END;
#ifdef _WIN32
            pt->zmumps_c(&pt->zid);
            FreeLibrary(pt->h_zmumps_c_dll);
#else
            zmumps_c(&pt->zid);
#endif
        } else {
            delete[] pt->rhs;
            pt->id.job = MUMPS_JOB_END;
#ifdef _WIN32
            pt->dmumps_c(&pt->id);
            FreeLibrary(pt->h_dmumps_c_dll);
#else
            dmumps_c(&pt->id);
#endif
        }
        if (pt->myid == 0) {
            std::string message = pt->ssExceptMsg.str();
            if (!message.empty()) {
                // terminating with solve error message
                throw PasoException(message);
            }
        }
// #ifdef ESYS_MPI
//         MUMPS_INT ierr = MPI_Finalize();
// #endif
        if (pt->verbose) {
            std::cout << "MUMPS: instance terminated." << std::endl;
        }
#endif
        delete pt;
        A->solver_p = NULL;
        //throw PasoException("Paso: ME done.");
    }
}


/// calls the solver
void MUMPS_solve(SparseMatrix_ptr A, double* out, double* in,
                   dim_t numRefinements, bool verbose)
{
#ifdef ESYS_HAVE_MUMPS
    if (! (A->type & (MATRIX_FORMAT_OFFSET1 + MATRIX_FORMAT_BLK1)) ) {
        throw PasoException("Paso: MUMPS requires CSR format with index offset 1 and block size 1.");
    }

    MUMPS_Handler* pt = reinterpret_cast<MUMPS_Handler*>(A->solver_p);
    if (pt == NULL) {
        pt = new MUMPS_Handler;
#ifdef _WIN32
        pt->h_dmumps_c_dll = LoadLibrary("libdmumps");
        if (pt->h_dmumps_c_dll == NULL)
            throw PasoException("Paso: MUMPS LoadLibrary failed - \"libdmumps\".");
        pt->dmumps_c = (DMUMPS_C_FUNC_PTR)GetProcAddress(pt->h_dmumps_c_dll, "dmumps_c");
        if (pt->dmumps_c == NULL)
            throw PasoException("Paso: MUMPS GetProcAddress failed - \"dmumps_c\".");
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
        pt->isComplex = false;
        pt->rhs = new double[n];
        std::memcpy(pt->rhs, in, n*sizeof(double));
// #ifdef ESYS_MPI
        // MUMPS_INT ierr;
        // ierr = MPI_Init(NULL, NULL);
        // ierr = MPI_Comm_rank(MPI_COMM_WORLD, &pt->myid);
// #endif

        // Initialize a MUMPS instance. Use MPI_COMM_WORLD
        pt->id.comm_fortran = MUMPS_USE_COMM_WORLD;
        pt->id.par = 1; 
        pt->id.sym = 0;
        pt->id.job = MUMPS_JOB_INIT;
#ifdef _WIN32
        pt->dmumps_c(&pt->id);
#else
        dmumps_c(&pt->id);
#endif
        // Define the problem on the host
        if (pt->myid == 0) {
            pt->id.n = n; pt->id.nnz = nnz;
            pt->id.irn = irn; pt->id.jcn = jcn;
            pt->id.a = A->val; pt->id.rhs = pt->rhs;
        }
        if (!pt->verbose) {
            // No outputs
            pt->id.ICNTL(1)=-1; pt->id.ICNTL(2)=-1; pt->id.ICNTL(3)=-1; pt->id.ICNTL(4)=0;
        }

        // Call the MUMPS package (analyse, factorization and solve).
        pt->id.job = 6;
#ifdef _WIN32
        pt->dmumps_c(&pt->id);
#else
        dmumps_c(&pt->id);
#endif
        if (pt->id.infog[0] < 0) {
            pt->ssExceptMsg << "(PROC " << pt->myid << ") MUMPS ERROR: INFOG(1)=" << pt->id.infog[0]
                << ", INFOG(2)=" << pt->id.infog[1];
        } else {
            std::memcpy(out, pt->rhs, n*sizeof(double));
            if (pt->id.infog[0] > 0) {
                std::cout << "(PROC " << pt->myid << ") MUMPS WARNING: INFOG(1)=" << pt->id.infog[0]
                    << ", INFOG(2)=" << pt->id.infog[1];
            }
            if (pt->verbose) {
                std::cout << "MUMPS: factorization and solve completed (time = "
                    << escript::gettime()-time0 << ")." << std::endl;
            }
        }
    }
#else // ESYS_HAVE_MUMPS
    throw PasoException("Paso: Not compiled with MUMPS.");
#endif
}

void MUMPS_solve(SparseMatrix_ptr A, cplx_t* out, cplx_t* in,
                   dim_t numRefinements, bool verbose)
{
#ifdef ESYS_HAVE_MUMPS
    if (! (A->type & (MATRIX_FORMAT_OFFSET1 + MATRIX_FORMAT_BLK1)) ) {
        throw PasoException("Paso: MUMPS requires CSR format with index offset 1 and block size 1.");
    }

    MUMPS_Handler* pt = reinterpret_cast<MUMPS_Handler*>(A->solver_p);
    if (pt == NULL) {
        pt = new MUMPS_Handler;
#ifdef _WIN32
        pt->h_zmumps_c_dll = LoadLibrary("libzmumps");
        if (pt->h_zmumps_c_dll == NULL)
            throw PasoException("Paso: MUMPS LoadLibrary failed - \"libzmumps\".");
        pt->zmumps_c = (ZMUMPS_C_FUNC_PTR)GetProcAddress(pt->h_zmumps_c_dll, "zmumps_c");
        if (pt->zmumps_c == NULL)
            throw PasoException("Paso: MUMPS GetProcAddress failed - \"zmumps_c\".");
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
        pt->isComplex = true;
        if (pt->verbose) {
            std::cout << "MUMPS in  ===>" << std::endl;
            std::cout << "isComplex = " << pt->isComplex << std::endl;
            std::cout << "n = " << n << std::endl;
            std::cout << "nnz = " << nnz << std::endl;
            std::cout << "cval = [";
            for (int i=0; i<nnz; i++) std::cout << A->cval[i] << ",";
            std::cout << "]" << std::endl;
            std::cout << "in = [";
            for (int i=0; i<n; i++) std::cout << "(" << in[i].real() << "," << in[i].imag() << "),";
            std::cout << "]" << std::endl;
            std::cout << "ptr = [";
            for (int i=0; i<n+1; i++) std::cout << A->pattern->ptr[i] << ",";
            std::cout << "]" << std::endl;
            std::cout << "index = [";
            for (int i=0; i<nnz; i++) std::cout << A->pattern->index[i] << ",";
            std::cout << "]" << std::endl;
            std::cout << "hb_row = [";
            for (int i=0; i<nnz; i++) std::cout << A->pattern->hb_row[i] << ",";
            std::cout << "]" << std::endl;
            std::cout << "hb_col = [";
            for (int i=0; i<nnz; i++) std::cout << A->pattern->hb_col[i] << ",";
            std::cout << "]" << std::endl;
        }
        pt->crhs = new cplx_t[n];
        std::memcpy(pt->crhs, in, n*sizeof(cplx_t));
// #ifdef ESYS_MPI
        // MUMPS_INT ierr;
        // ierr = MPI_Init(NULL, NULL);
        // ierr = MPI_Comm_rank(MPI_COMM_WORLD, &pt->myid);
// #endif

        // Initialize a MUMPS instance. Use MPI_COMM_WORLD
        pt->zid.comm_fortran = MUMPS_USE_COMM_WORLD;
        pt->zid.par = 1; 
        pt->zid.sym = 0;
        pt->zid.job = MUMPS_JOB_INIT;
#ifdef _WIN32
        pt->zmumps_c(&pt->zid);
#else
        zmumps_c(&pt->zid);
#endif
        // Define the problem on the host
        if (pt->myid == 0) {
            pt->zid.n = n; pt->zid.nnz = nnz;
            pt->zid.irn = irn; pt->zid.jcn = jcn;
            pt->zid.a = reinterpret_cast<ZMUMPS_COMPLEX*>(A->cval);
            pt->zid.rhs = reinterpret_cast<ZMUMPS_COMPLEX*>(pt->crhs);
        }
        if (!pt->verbose) {
            // No outputs
            pt->zid.ICNTL(1)=-1; pt->zid.ICNTL(2)=-1; pt->zid.ICNTL(3)=-1; pt->zid.ICNTL(4)=0;
        }

        // Call the MUMPS package (analyse, factorization and solve).
        pt->zid.job = 6;
#ifdef _WIN32
        pt->zmumps_c(&pt->zid);
#else
        zmumps_c(&pt->zid);
#endif
        if (pt->zid.infog[0] < 0) {
            pt->ssExceptMsg << "(PROC " << pt->myid << ") MUMPS ERROR: INFOG(1)=" << pt->zid.infog[0]
                << ", INFOG(2)=" << pt->zid.infog[1];
        } else {
            std::memcpy(out, reinterpret_cast<cplx_t*>(pt->crhs), n*sizeof(cplx_t));
            if (pt->zid.infog[0] > 0) {
                std::cout << "(PROC " << pt->myid << ") MUMPS WARNING: INFOG(1)=" << pt->zid.infog[0]
                    << ", INFOG(2)=" << pt->zid.infog[1];
            }
            if (pt->verbose) {
                std::cout << "MUMPS out ===>" << std::endl;
                std::cout << "out = [";
                for (int i=0; i<n; i++) std::cout << "(" << out[i].real() << "," << out[i].imag() << "),";
                std::cout << "]" << std::endl;
                std::cout << "MUMPS: factorization and solve completed (time = "
                    << escript::gettime()-time0 << ")." << std::endl;
            }
        }
    }
#else // ESYS_HAVE_MUMPS
    throw PasoException("Paso: Not compiled with MUMPS.");
#endif
}

} // namespace paso

