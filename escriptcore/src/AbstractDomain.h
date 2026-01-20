
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


#ifndef __ESCRIPT_ABSTRACTDOMAIN_H__
#define __ESCRIPT_ABSTRACTDOMAIN_H__

#include "system_dep.h"
#include "DataTypes.h"
#include "EsysMPI.h"
#include "Pointers.h"

#include <boost/python/tuple.hpp>
#ifdef ESYS_HAVE_BOOST_NUMPY
#include <boost/python/numpy.hpp>
#endif

#ifdef ESYS_HAVE_MPI4PY
#include <mpi4py/mpi4py.h>
#endif

#include <vector>
#include <string>

namespace escript {

// class forward declarations
class Data;
class FunctionSpace;
class AbstractDomain;

typedef POINTER_WRAPPER_CLASS(AbstractDomain) Domain_ptr;
typedef POINTER_WRAPPER_CLASS(const AbstractDomain) const_Domain_ptr;

/**
   \brief
   Base class for all escript domains.
*/
class ESCRIPT_DLL_API AbstractDomain: public REFCOUNT_BASE_CLASS(AbstractDomain)
{
public:
    typedef int StatusType;

    /**
    \brief Returns smart pointer which is managing this object.
    If one does not exist yet it creates one.

    Note: This is _not_ equivalent to weak_ptr::lock.
    */
    Domain_ptr getPtr();

    const_Domain_ptr getPtr() const; 

    /**
     \brief
     Destructor for AbstractDomain.
    */
    virtual ~AbstractDomain() {}

    /**
     \brief
     returns a shared pointer to the MPI information wrapper for this domain
    */
    virtual JMPI getMPI() const { return m_mpiInfo; }

    /**
     \brief
     return the number of processors used for this domain
    */
    virtual int getMPISize() const { return m_mpiInfo->size; }

    /**
     \brief
     return the number MPI rank of this processor
    */
    virtual int getMPIRank() const { return m_mpiInfo->rank; }

    /**
     \brief
     If compiled for MPI then execute an MPI_Barrier, else do nothing
    */
    virtual void MPIBarrier() const {
#ifdef ESYS_MPI
        MPI_Barrier(m_mpiInfo->comm);
#endif
    }

    /**
     \brief
     Return true if this is the root MPI rank (rank 0), else false
    */
    virtual bool isRootRank() const { return m_mpiInfo->rank == 0; }

    /**
       \brief
       get the communicator for this domain.
       Returns an integer on non-MPI builds
       Routine must be implemented by the DomainAdapter.
    */
    virtual MPI_Comm getMPIComm() const { return m_mpiInfo->comm; }

    /**
       \brief
       get the communicator for this domain as a Python mpi4py.MPI.Comm object.
       Returns None if MPI or mpi4py support is not enabled.
       This allows Python code to retrieve the domain's communicator.
    */
    virtual boost::python::object getMPICommPython() const {
        return makePyCommFromMPI(m_mpiInfo->comm);
    }

    /**
       \brief
       Returns true if the given integer is a valid function space type
       for this domain.
    */
    virtual bool isValidFunctionSpaceType(int functionSpaceType) const = 0;

    /**
     \brief
     Return a description for this domain.
    */
    virtual std::string getDescription() const = 0;

    /**
     \brief
     Return a description for the given function space type code.
    */
    virtual std::string functionSpaceTypeAsString(int functionSpaceType) const = 0;

    /**
     \brief
      Returns the spatial dimension of the domain.

      This has to be implemented by the actual Domain adapter.
    */
    virtual int getDim() const = 0;

    /**
     \brief
     Return true if given domains are equal.
    */
    virtual bool operator==(const AbstractDomain& other) const = 0;

    /**
     \brief
     Return true if given domains are not equal.
    */
    virtual bool operator!=(const AbstractDomain& other) const = 0;

    /**
     \brief
     Writes the domain to an external file filename.

     This has to be implemented by the actual Domain adapter.
    */
    virtual void write(const std::string& filename) const = 0;

    /**
     \brief
     dumps the domain to an external file filename.

     This has to be implemented by the actual Domain adapter.
    */
    virtual void dump(const std::string& filename) const = 0;

    /**
       \brief
       Returns the number of data points per sample, and the number of samples
       as a pair.

        This has to be implemented by the actual Domain adapter.

        \param functionSpaceCode Input - Code for the function space type.
        \return pair, first - number of data points per sample,
                second - number of samples
    */
    virtual std::pair<int,DataTypes::dim_t> getDataShape(int functionSpaceCode) const = 0;

    /**
       \brief
       Return the tag key for the given sample number.
       \param functionSpaceType Input - The function space type.
       \param sampleNo Input - The sample number.
    */
    virtual int getTagFromSampleNo(int functionSpaceType, DataTypes::index_t sampleNo) const = 0;

    /**
       \brief
       sets a map from a clear tag name to a tag key
       \param name Input - tag name.
       \param tag Input - tag key.
    */
    virtual void setTagMap(const std::string& name,  int tag) = 0;

    /**
       \brief
       Return the tag key for tag name.
       \param name Input - tag name
    */
    virtual int getTag(const std::string& name) const = 0;

    /**
       \brief
       Returns True if name is a defined tag name
       \param name Input - tag name
    */
    virtual bool isValidTagName(const std::string& name) const;

    /**
       \brief
       Returns all tag names in a single string sperated by commas
    */
    virtual std::string showTagNames() const = 0;

    /**
       \brief
       Returns a borrowed pointer to the sample reference number id list
       \param functionSpaceType Input - The function space type.
    */
    // ESCRIPT_DLL_API class already exported
    virtual const DataTypes::dim_t* borrowSampleReferenceIDs(int functionSpaceType) const = 0;

    /**
       \brief
       Assigns new location to the domain.

       This has to be implemented by the actual Domain adapter.
    */
    virtual void setNewX(const escript::Data& arg) = 0;

    /**
       \brief
       Interpolates data given on source onto target where source and target
       have to be given on the same domain.

       This has to be implemented by the actual Domain adapter.
    */
    virtual void interpolateOnDomain(escript::Data& target, const escript::Data& source) const = 0;

    /**
       \brief True if interpolation is possible from source to target 
    */
    virtual bool probeInterpolationOnDomain(int functionSpaceType_source,
                                      int functionSpaceType_target) const = 0;

    /**
       \brief
       Preferred direction of interpolation.
       If you really need to test for a particular direction, then use
       probeInterpolation.

       \return 0 for not possible,  1 for possible and preferred, -1 other
       direction preferred (does not mean this direction is possible)
    */ 
    virtual signed char preferredInterpolationOnDomain(
                                      int functionSpaceType_source,
                                      int functionSpaceType_target) const = 0;

    /**
       \brief given a vector of FunctionSpace type codes, pass back a code
              which then can all be interpolated to.
       \note This method must be called on the domain which the FunctionSpaces
             point to
       \return true is result is valid, false if not
    */
    virtual bool
    commonFunctionSpace(const std::vector<int>& fs, int& resultcode) const = 0;

    /**
     \brief
     Interpolates data given on source onto target where source and target are given on different domains.

     This has to be implemented by the actual Domain adapter.
    */
    virtual void interpolateAcross(escript::Data& target, const escript::Data& source) const = 0;

    virtual bool probeInterpolationAcross(int functionSpaceType_source,
                                      const AbstractDomain& targetDomain,
                                      int functionSpaceType_target) const = 0;

    /**
     \brief
     Returns locations in the domain. The function space is chosen appropriately.
    */
    virtual escript::Data getX() const = 0;

#ifdef ESYS_HAVE_BOOST_NUMPY
    /**
     \brief
     Returns locations in the domain as a numpy ndarray. The function space is chosen appropriately.
    */
    virtual boost::python::numpy::ndarray getNumpyX() const = 0;
#endif

#if defined(ESYS_MPI) && defined(ESYS_HAVE_MPI4PY)
    // Changes the communicator
   // void setMPIComm(boost::python::object py_comm);
#endif

    /**
     \brief
     Return boundary normals. The function space is chosen appropriately.
    */
    virtual escript::Data getNormal() const = 0;

    /**
     \brief
     Returns the local size of samples. The function space is chosen appropriately.
    */
    virtual escript::Data getSize() const = 0;
  
    /**
     \brief
     Copies the location of data points on the domain into out.
     The actual function space to be considered
     is defined by out. out has to be defined on this.

     This has to be implemented by the actual Domain adapter.
    */
    virtual void setToX(escript::Data& out) const = 0;

    /**
     \brief
     Copies the surface normals at data points into out.
     The actual function space to be considered
     is defined by out. out has to be defined on this.

     This has to be implemented by the actual Domain adapter.
    */
    virtual void setToNormal(escript::Data& out) const = 0;

    /**
     \brief
     Copies the size of samples into out. The actual
     function space to be considered
     is defined by out. out has to be defined on this.

     This has to be implemented by the actual Domain adapter.
    */
    virtual void setToSize(escript::Data& out) const = 0;

    /**
     \brief
     Copies the gradient of arg into grad. The actual function space to be considered
     for the gradient is defined by grad. arg and grad have to be defined on this.

     This has to be implemented by the actual Domain adapter.
    */
    virtual void setToGradient(escript::Data& grad, const escript::Data& arg) const = 0;

    /**
    \brief True if this rank owns the sample(id)
    Must be implemented by the Domain adapter
    */
    virtual bool ownSample(int fs_code, DataTypes::index_t id) const = 0;

    /**
       \brief
       assigns new tag newTag to all samples of functionspace with a positive
       value of mask for any its sample point.
    */
    virtual void setTags(int functionSpaceType, int newTag, const escript::Data& mask) const = 0;

    /**
     \brief
     returns true if data on this domain and a function space of type functionSpaceCode has to
     considered as cell centered data.

     This has to be implemented by the actual Domain adapter.
    */
    virtual bool isCellOriented(int functionSpaceCode) const = 0;

    /**
     \brief
      Returns a status indicator of the domain. The status identifier should be unique over 
      the live time if the object but may be updated if changes to the domain happen, e.g. 
      modifications to its geometry. 

     This has to be implemented by the actual Domain adapter.
    */
    virtual StatusType getStatus() const;

    /**
     \brief
     Throw a standard exception. This function is called if any attempt 
     is made to use a base class function.
    */
    void throwStandardException(const std::string& functionName) const;

    /**
       \brief
       returns the number of tags in use and a pointer to an array with the
       number of tags in use
    */
    virtual int getNumberOfTagsInUse(int functionSpaceCode) const = 0;

    virtual const int* borrowListOfTagsInUse(int functionSpaceCode) const = 0;

    /**
    \brief Checks if this domain allows tags for the specified functionSpaceCode.
    */
    virtual bool canTag(int functionspacecode) const = 0;

    /**
    \brief returns the approximation order used for a function space functionSpaceCode
    */
    virtual int getApproximationOrder(const int functionSpaceCode) const = 0;

    virtual bool supportsContactElements() const = 0;

    /**
    * \brief true if this domain can handle to specified tuple of filter options.
    */
    virtual bool supportsFilter(const boost::python::tuple& t) const;

    /**
    * \brief Fills the data object with filtered random values
    */
    virtual escript::Data randomFill(const DataTypes::ShapeType& shape,
                                const FunctionSpace& what, long seed,
                                const boost::python::tuple& filter) const = 0;

protected:
    /**
     \brief
     MPI information wrapper for this domain.
     Contains communicator, rank, and size information.
     All domains should initialize this in their constructors.
    */
    JMPI m_mpiInfo;

    /**
     \brief
     Default constructor for AbstractDomain.
     Domains should call this and initialize m_mpiInfo in their constructors.
    */
    AbstractDomain() {}

    /**
     \brief
     Constructor that initializes MPI information.
     \param mpiInfo MPI information wrapper
    */
    explicit AbstractDomain(JMPI mpiInfo) : m_mpiInfo(mpiInfo) {}
};

} // end of namespace

#endif // __ESCRIPT_ABSTRACTDOMAIN_H__

