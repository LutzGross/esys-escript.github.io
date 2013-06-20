
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


#ifndef INC_ESYS_MPI
#define INC_ESYS_MPI

#include "system_dep.h"
#include "types.h"

#include <sstream>

#ifdef ESYS_MPI
   #include "mpi_C.h"
#else
   typedef int MPI_Comm;
   typedef int MPI_Request;
   #define MPI_INT 6
   #define MPI_DOUBLE 11
   #define MPI_COMM_WORLD 91
#endif

typedef int Esys_MPI_rank;

#define ESYS_MPI_TODO 	{ fprintf( stdout, "\nTODO : %s:%d\n", __FILE__, __LINE__);	MPI_Finalize(); exit(1); }

/* Datatypes */

/** \brief Struct that holds MPI communicator, rank, size and a tag counter */
struct Esys_MPIInfo {
  dim_t reference_counter;
  int size;
  Esys_MPI_rank rank;
  MPI_Comm comm;
  int msg_tag_counter;
};

typedef struct Esys_MPIInfo Esys_MPIInfo;

/* Function prototypes */

ESYSUTILS_DLL_API
Esys_MPIInfo* Esys_MPIInfo_alloc( MPI_Comm comm );

ESYSUTILS_DLL_API
void          Esys_MPIInfo_free( Esys_MPIInfo* );

ESYSUTILS_DLL_API
Esys_MPIInfo *Esys_MPIInfo_getReference( Esys_MPIInfo* in );

ESYSUTILS_DLL_API
int           Esys_MPIInfo_initialized( void );

ESYSUTILS_DLL_API
index_t Esys_MPIInfo_mod(index_t n, index_t k);

ESYSUTILS_DLL_API
dim_t Esys_MPIInfo_setDistribution(Esys_MPIInfo* in ,index_t min_id,index_t max_id,index_t* distribution);

ESYSUTILS_DLL_API
void Esys_MPIInfo_Split( Esys_MPIInfo *mpi_info, dim_t n, dim_t* local_N,index_t* offset); 

ESYSUTILS_DLL_API
bool_t Esys_MPIInfo_noError( Esys_MPIInfo *mpi_info);

namespace esysUtils {

/// Appends MPI rank to a file name if MPI size > 1
ESYSUTILS_DLL_API
inline std::string appendRankToFileName(const std::string &fileName,
                                        int mpiSize, int mpiRank)
{
    std::stringstream ss(fileName);
    if (mpiSize > 1) {
        ss << '.';
        ss.fill('0');
        ss.width(4);
        ss << mpiRank;
    }
    std::string result(ss.str());
    return result;
}

} // namespace esysUtils

#endif /* INC_ESYS_MPI */

