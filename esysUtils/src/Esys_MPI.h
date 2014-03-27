
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


#ifndef INC_ESYS_MPI
#define INC_ESYS_MPI

#include "system_dep.h"
#include "types.h"

#include <sstream>
#include <boost/shared_ptr.hpp>

#ifdef ESYS_MPI
   #include "mpi_C.h"
#else
   typedef int MPI_Comm;
   typedef int MPI_Request;
   typedef int MPI_Op;
   #define MPI_INT 6
   #define MPI_DOUBLE 11
   #define MPI_COMM_WORLD 91
   #define MPI_COMM_NULL 0
   
// MPI_Op replacements for non-MPI - these values are arbitrary

   #define MPI_SUM 100

// end MPI_op

   
#endif

typedef int Esys_MPI_rank;

#define ESYS_MPI_TODO 	{ fprintf( stdout, "\nTODO : %s:%d\n", __FILE__, __LINE__);	MPI_Finalize(); exit(1); }

// Modding by 7 digit prime to avoid overflow
#define ESYS_MPI_INC_COUNTER(V,I) {(V).msg_tag_counter=((V).msg_tag_counter+(I))%1010201;}
#define ESYS_MPI_SET_COUNTER(V,I) {(V).msg_tag_counter=(I)%1010201;}

namespace esysUtils {

/** \brief tag reserved for use by SubWorld code
    this value should be higher than the modulus used in JMPI_::setCounter, apart from that, its value
    is not particularly significant.
*/
ESYSUTILS_DLL_API
inline const int getSubWorldTag()	
{
    return ('S'<< 24) + ('u' << 16) + ('b' << 8) + 'W';  
}
  
class JMPI_;

typedef boost::shared_ptr<JMPI_> JMPI;

class JMPI_
{
public:
    ~JMPI_();
    int size;
    Esys_MPI_rank rank;
    MPI_Comm comm;
    int msg_tag_counter;
    bool ownscomm;	// if true, destroy comm on destruct    
    
    dim_t setDistribution(index_t min_id,index_t max_id,index_t* distribution);
    void split(dim_t N, dim_t* local_N,index_t* offset);     
    
    void incCounter(int i)
    {
	msg_tag_counter+=i;
	msg_tag_counter%=1010201;		// there is no particular significance here other than being 7 digits 
    }					// and prime (because why not). It just needs to be big.
    
    void setCounter(int i)
    {
	msg_tag_counter%=1010201;
    }
private:
    JMPI_(MPI_Comm comm, bool ocomm);
    friend JMPI makeInfo(MPI_Comm comm, bool owncom);
};

JMPI makeInfo(MPI_Comm comm, bool owncom=false);

ESYSUTILS_DLL_API
bool Esys_MPIInfo_noError( const JMPI& mpi_info);

ESYSUTILS_DLL_API
index_t mod_rank(index_t n, index_t k);


/// Appends MPI rank to a file name if MPI size > 1
ESYSUTILS_DLL_API
inline std::string appendRankToFileName(const std::string &fileName,
                                        int mpiSize, int mpiRank)
{
    std::stringstream ss;
    ss << fileName;
    if (mpiSize > 1) {
        ss << '.';
        ss.fill('0');
        ss.width(4);
        ss << mpiRank;
    }
    std::string result(ss.str());
    return result;
}

/* has the have sub-communicators been created? */
bool getSplitWorld();
/* record that a sub-communicator has been created or used */
void splitWorld();

/* returns the max of inputs on all ranks -- or just sends the input back on nompi */
bool checkResult(int& input, int& output, MPI_Comm& comm);

// ensure that the any ranks with an empty src argument end up with the string from
// one of the other ranks
// with no-mpi, it makes dest point at a copy of src
bool shipString(const char* src, char** dest, MPI_Comm& comm);

} // namespace esysUtils

#endif /* INC_ESYS_MPI */

