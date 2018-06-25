
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
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

#ifndef __ESCRIPT_ESYSMPI_H__
#define __ESCRIPT_ESYSMPI_H__

#include <escript/DataTypes.h>

#include <ctime>
#include <sstream>

#include <boost/shared_ptr.hpp>

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
   #define MPI_COMM_NULL 0

   // MPI_Op replacements for non-MPI - these values are arbitrary
   #define MPI_SUM 100
   #define MPI_MIN 101
   #define MPI_MAX 102

   #define MPI_OP_NULL 17
   // end MPI_op

#endif // ESYS_MPI

namespace escript {

/** \brief tag reserved for use by SubWorld code
    This value should be higher than the modulus used in JMPI_::setCounter.
    Apart from that, its value is not particularly significant.
*/
inline int getSubWorldTag()
{
    return (('S'<< 24) + ('u' << 16) + ('b' << 8) + 'W')%1010201;
}

class JMPI_;

typedef boost::shared_ptr<JMPI_> JMPI;

/// creates a JMPI shared pointer from MPI communicator
/// if owncom is true, the communicator is freed when mpi info is destroyed.
JMPI makeInfo(MPI_Comm comm, bool owncom=false);

class JMPI_
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
    JMPI_(MPI_Comm comm, bool owncomm);
    friend JMPI makeInfo(MPI_Comm comm, bool owncom);

    bool ownscomm;
    int msg_tag_counter;
};

// Does not cope with nested calls
class NoCOMM_WORLD
{
public:
    NoCOMM_WORLD();
    ~NoCOMM_WORLD();
    static bool active();
};

/// Everyone puts in their error code and everyone gets the largest one
bool checkResult(int input, int& output, const JMPI& comm);

/// ensure that the any ranks with an empty src argument end up with the
/// string from one of the other ranks.
/// With no MPI, it makes dest point at a copy of src.
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

