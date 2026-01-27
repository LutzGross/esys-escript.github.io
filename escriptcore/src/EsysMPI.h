
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __ESCRIPT_ESYSMPI_H__
#define __ESCRIPT_ESYSMPI_H__

#include "system_dep.h"

#include <escript/DataTypes.h>
#include "EsysException.h"
#include <ctime>
#include <sstream>

#include <boost/shared_ptr.hpp>
#include <boost/python.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ESYS_MPI
#include <mpi.h>

#ifdef ESYS_INDEXTYPE_LONG
#define MPI_DIM_T MPI_LONG
#else
#define MPI_DIM_T MPI_INT
#endif

#else
   typedef int MPI_Comm;
   typedef int MPI_Request;
   typedef int MPI_Op;
   typedef int MPI_Status;

   #define MPI_INT 6
   #define MPI_DOUBLE 11
   #define MPI_COMM_WORLD 91
   #define MPI_COMM_SELF 92
   #define MPI_COMM_NULL 0

   // MPI_Op replacements for non-MPI - these values are arbitrary
   #define MPI_SUM 100
   #define MPI_MIN 101
   #define MPI_MAX 102

   #define MPI_OP_NULL 17
   // end MPI_op

#endif // ESYS_MPI

#ifdef ESYS_HAVE_MPI4PY
#include <mpi4py/mpi4py.h>
#endif

namespace escript {

inline static MPI_Comm * extractMPICommunicator(boost::python::object py_comm)
{
    MPI_Comm *comm_p = 0;
    if ( py_comm.is_none() )  {
        *comm_p = MPI_COMM_WORLD;
    } else {
        #ifdef ESYS_HAVE_MPI4PY
        PyObject* py_obj = py_comm.ptr();
        comm_p = PyMPIComm_Get(py_obj);
        if (comm_p == NULL)
            throw EsysException("Null communicator.");
        #else
        *comm_p = MPI_COMM_WORLD;
        #endif
    }
    return comm_p;
}

class JMPI_;

typedef boost::shared_ptr<JMPI_> JMPI;

/// ensures MPI has been initialized, throws exception if not
/// This should be called before using MPI_COMM_WORLD or any MPI operations
/// Does not call MPI_Init - initialization must happen at a higher level
ESCRIPT_DLL_API
void ensureMPIInitialized();

/// creates a JMPI shared pointer from MPI communicator
ESCRIPT_DLL_API
JMPI makeInfo(MPI_Comm comm);

/// creates a JMPI shared pointer from optional Python mpi4py communicator
/// if py_comm is None or not provided, uses MPI_COMM_WORLD
/// requires ESYS_HAVE_MPI4PY to be enabled for custom communicators
ESCRIPT_DLL_API
JMPI makeInfoFromPyComm(boost::python::object py_comm);

/// converts an MPI_Comm to a Python mpi4py.MPI.Comm object
/// Returns None if MPI is not enabled or mpi4py is not available
ESCRIPT_DLL_API
boost::python::object makePyCommFromMPI(MPI_Comm comm);

class ESCRIPT_DLL_API JMPI_
{
public:
    ~JMPI_();

    ///
    DataTypes::dim_t setDistribution(DataTypes::index_t min_id,
                                     DataTypes::index_t max_id,
                                     DataTypes::index_t* distribution);

    ///
    void split(DataTypes::dim_t N, DataTypes::dim_t* local_N,
               DataTypes::index_t* offset);

    /// N = #CPUs, k is a CPU number but out of range or even negative.
    /// Return a CPU number in 0...N-1.
    inline int mod_rank(int k) const
    {
        int out=0;
#ifdef ESYS_MPI
        if (size > 1) {
            const int q = k/size;
            if (k > 0) {
               out=k-size*q;
            } else if (k < 0) {
               out=k-size*(q-1);
            }
        }
#endif
        return out;
    }

    /// appends MPI rank to a file name if MPI size > 1
    inline std::string appendRankToFileName(const std::string& fileName) const
    {
#ifdef ESYS_MPI
        if (size > 1) {
            std::stringstream ss;
            ss << fileName << '.';
            ss.fill('0');
            ss.width(4);
            ss << rank;
            return ss.str();
        }
#endif
        return fileName;
    }

    /// returns the current value of the message tag counter
    inline int counter() const
    {
        return msg_tag_counter;
    }

    /// increments the message tag counter by `i`
    inline void incCounter(int i=1)
    {
        msg_tag_counter+=i;
        // there is no particular significance here other than being 7 digits
        // and prime (because why not). It just needs to be big.
        msg_tag_counter %= 1010201;
    }

    /// sets the message tag counter to `value`
    inline void setCounter(int value)
    {
        msg_tag_counter = value%1010201;
    }

    /// returns true if this has a valid MPI communicator
    inline bool isValid() const
    {
        return comm!=MPI_COMM_NULL;
    }

    int size;
    int rank;
    MPI_Comm comm;

private:
    JMPI_(MPI_Comm comm);
    friend ESCRIPT_DLL_API JMPI makeInfo(MPI_Comm comm);
    int msg_tag_counter;
};

/// Everyone puts in their error code and everyone gets the largest one
ESCRIPT_DLL_API
bool checkResult(int input, int& output, const JMPI& comm);

/// ensure that the any ranks with an empty src argument end up with the
/// string from one of the other ranks.
/// With no MPI, it makes dest point at a copy of src.
ESCRIPT_DLL_API
bool shipString(const char* src, char** dest, MPI_Comm& comm);

/// returns the current ticks for timing
inline double gettime()
{
    double out;
#ifdef ESYS_MPI
    out = MPI_Wtime();
#else
#ifdef _OPENMP 
    out=omp_get_wtime();
#else
    out=((double) clock())/CLOCKS_PER_SEC;
#endif
#endif
    return out;
}




} // namespace escript

#endif // __ESCRIPT_ESYSMPI_H__

