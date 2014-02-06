
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


#if !defined escript_AbstractDomain_20040609_H
#define escript_AbstractDomain_20040609_H

#ifdef BADPYTHONMACROS
// This hack is required for BSD/OSX builds with python 2.7
// (and possibly others).  It must be the first include.
// From bug reports online it seems that python redefines
// some c macros that are functions in c++.
// c++ doesn't like that!
#include <Python.h>
#undef BADPYTHONMACROS
#endif


#include "system_dep.h"

#include <vector>
#include <string>
#include <map>
#include <boost/python/dict.hpp>
#include <boost/python/list.hpp>
#include "esysUtils/Esys_MPI.h"


#include "Pointers.h"

namespace escript {
// class forward declarations
class Data;
/**
   \brief
   Base class for all escript domains.

   Description:
   Base class for all escript domains.
*/

class AbstractDomain;

typedef POINTER_WRAPPER_CLASS(AbstractDomain) Domain_ptr;
typedef POINTER_WRAPPER_CLASS(const AbstractDomain) const_Domain_ptr;

class AbstractDomain : public REFCOUNT_BASE_CLASS(AbstractDomain){

 public:

/**
\brief Returns smart pointer which is managing this object.
If one does not exist yet it creates one.

Note: This is _not_ equivalent to weak_ptr::lock.
*/
   ESCRIPT_DLL_API
   Domain_ptr getPtr();
   ESCRIPT_DLL_API
   const_Domain_ptr getPtr() const; 

   // structure holding values for X, size and normal
   typedef int StatusType;
   struct ValueBuffer
   {
       StatusType m_status;
       boost::shared_ptr<Data> m_data;
   };
   typedef struct ValueBuffer ValueBuffer;

   // 
   // map from function space type code to value buffer
   typedef std::map<int, ValueBuffer> BufferMapType;


  /**
     \brief
     Default constructor for AbstractDomain.

     Description:
     Default constructor for AbstractDomain. As the name suggests
     this is intended to be an abstract base class but by making it
     constructable we avoid a boost.python wrapper class. A call to 
     almost any of the base class functions will throw an exception
     as they are not intended to be used directly, but are overridden
     by the underlying solver package which escript is linked to.

     By default, this class is overridden by the class NullDomain.

     Preconditions:
     Describe any preconditions.

     Throws:
     Describe any exceptions thrown.
  */
  ESCRIPT_DLL_API
  AbstractDomain();

  /**
     \brief
     Destructor for AbstractDomain.

     Description:
     Destructor for AbstractDomain.
  */
  ESCRIPT_DLL_API
  virtual ~AbstractDomain();

  /**
     \brief
     return the number of processors used for this domain
  */
  ESCRIPT_DLL_API
  virtual int getMPISize() const;
  /**
     \brief
     return the number MPI rank of this processor
  */

  ESCRIPT_DLL_API
  virtual int getMPIRank() const;

  /**
     \brief
     If compiled for MPI then execute an MPI_Barrier, else do nothing
  */

  ESCRIPT_DLL_API
  virtual void MPIBarrier() const;
  /**
     \brief
     Return true if on MPI master, else false
  */

  ESCRIPT_DLL_API
  virtual bool onMasterProcessor() const;

  /**
  \brief get the communicator for this domain.
  Returns 0 on non-MPI builds
  Routine must be implemented by the DomainAdapter
  */
  ESCRIPT_DLL_API
  virtual
#ifdef ESYS_MPI
  MPI_Comm
#else
  unsigned int
#endif
  getMPIComm() const;

  /**
     \brief
     Returns true if the given integer is a valid function space type
     for this domain.
  */
  ESCRIPT_DLL_API
  virtual bool isValidFunctionSpaceType(int functionSpaceType) const;

  /**
     \brief
     Return a description for this domain.
  */
  ESCRIPT_DLL_API
  virtual std::string getDescription() const;

  /**
     \brief
     Return a description for the given function space type code.
  */
  ESCRIPT_DLL_API
  virtual std::string functionSpaceTypeAsString(int functionSpaceType) const;

  /**
     \brief
      Returns the spatial dimension of the domain.

      This has to be implemented by the actual Domain adapter.
  */
  ESCRIPT_DLL_API
  virtual int getDim() const;

  /**
     \brief
     Return true if given domains are equal.
  */
  ESCRIPT_DLL_API
  virtual bool operator==(const AbstractDomain& other) const;
  ESCRIPT_DLL_API
  virtual bool operator!=(const AbstractDomain& other) const;

  /**
     \brief
     Writes the domain to an external file filename.

     This has to be implemented by the actual Domain adapter.
  */
  ESCRIPT_DLL_API
  virtual void write(const std::string& filename) const;

  /**
     \brief
     dumps the domain to an external file filename.

     This has to be implemented by the actual Domain adapter.
  */
  ESCRIPT_DLL_API
  virtual void dump(const std::string& filename) const;

  /**
     \brief
     Return the number of data points per sample, and the number of samples as a pair.

     This has to be implemented by the actual Domain adapter.

     \param functionSpaceCode Input - Code for the function space type.
     \return pair, first - number of data points per sample, second - number of samples
  */
  ESCRIPT_DLL_API
  virtual std::pair<int,int> getDataShape(int functionSpaceCode) const;

  /**
     \brief
     Return the tag key for the given sample number.
     \param functionSpaceType Input - The function space type.
     \param sampleNo Input - The sample number.
  */
  ESCRIPT_DLL_API
  virtual int getTagFromSampleNo(int functionSpaceType, int sampleNo) const;

  /**
     \brief
     sets a map from a clear tag name to a tag key
     \param name Input - tag name.
     \param tag Input - tag key.
  */
  ESCRIPT_DLL_API
  virtual void setTagMap(const std::string& name,  int tag);

  /**
     \brief
     Return the tag key for tag name.
     \param name Input - tag name
  */
  ESCRIPT_DLL_API
  virtual int getTag(const std::string& name) const;

  /**
     \brief
     Returns True if name is a defined tag name
     \param name Input - tag name
  */
  ESCRIPT_DLL_API
  virtual bool isValidTagName(const std::string& name) const;

  /**
     \brief
     Returns all tag names in a single string sperated by commas
  */
  ESCRIPT_DLL_API
  virtual std::string showTagNames() const;

  /**
     \brief
     Return a borrowed pointer to the sample reference number id list
     \param functionSpaceType Input - The function space type.
  */
  ESCRIPT_DLL_API
  virtual const int* borrowSampleReferenceIDs(int functionSpaceType) const;

  /**
     \brief
     Assigns new location to the domain.

     This has to be implemented by the actual Domain adapter.
  */
  ESCRIPT_DLL_API
  virtual void setNewX(const escript::Data& arg);

  /**
     \brief
     Interpolates data given on source onto target where source and target have to be given on the same domain.

     This has to be implemented by the actual Domain adapter.
  */
  ESCRIPT_DLL_API
  virtual void interpolateOnDomain(escript::Data& target,const escript::Data& source) const;

  /**
  \brief True if interpolation is possible from source to target 
  */
  ESCRIPT_DLL_API
  virtual bool probeInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const;

  /**
   \brief Preferred direction of interpolation.

   If you really need to test for a particular direction, then use probeInterpolation.
   
   \return 0 for not possible,  1 for possible and preferred, -1 other direction preferred (does not mean this direction is possible)
  */ 
  virtual signed char preferredInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const;

  /**
    \brief given a vector of FunctionSpace type codes, pass back a code which then can all be interpolated to.
    \note This method must be called on the domain which the FunctionSpaces point to
    \return true is result is valid, false if not
  */
  ESCRIPT_DLL_API
  virtual
  bool
  commonFunctionSpace(const std::vector<int>& fs, int& resultcode) const;

  /**
     \brief
     Interpolates data given on source onto target where source and target are given on different domains.

     This has to be implemented by the actual Domain adapter.
  */
  ESCRIPT_DLL_API
  virtual void interpolateACross(escript::Data& target, const escript::Data& source) const;
  ESCRIPT_DLL_API
  virtual bool probeInterpolationACross(int functionSpaceType_source,const AbstractDomain& targetDomain, int functionSpaceType_target) const;

  /**
     \brief
     Returns locations in the domain. The function space is chosen appropriately.
  */
  ESCRIPT_DLL_API
  virtual escript::Data getX() const;

  /**
     \brief
     Return boundary normals. The function space is chosen appropriately.
  */
  ESCRIPT_DLL_API
  virtual escript::Data getNormal() const;

  /**
     \brief
     Returns the local size of samples. The function space is chosen appropriately.
  */
  ESCRIPT_DLL_API
  virtual escript::Data getSize() const;
  
  /**
     \brief
     Copies the location of data points on the domain into out.
     The actual function space to be considered
     is defined by out. out has to be defined on this.

     This has to be implemented by the actual Domain adapter.
  */
  ESCRIPT_DLL_API
  virtual void setToX(escript::Data& out) const;

  /**
     \brief
     Copies the surface normals at data points into out.
     The actual function space to be considered
     is defined by out. out has to be defined on this.

     This has to be implemented by the actual Domain adapter.
  */
  ESCRIPT_DLL_API
  virtual void setToNormal(escript::Data& out) const;

  /**
     \brief
     Copies the size of samples into out. The actual
     function space to be considered
     is defined by out. out has to be defined on this.

     This has to be implemented by the actual Domain adapter.
  */
  ESCRIPT_DLL_API
  virtual void setToSize(escript::Data& out) const;

  /**
     \brief
     Copies the gradient of arg into grad. The actual function space to be considered
     for the gradient is defined by grad. arg and grad have to be defined on this.

     This has to be implemented by the actual Domain adapter.
  */
  ESCRIPT_DLL_API
  virtual void setToGradient(escript::Data& grad, const escript::Data& arg) const;

  /**
  \brief True if this rank owns the sample(id)
  Must be implemented by the Domain adapter
  */
  ESCRIPT_DLL_API
  virtual bool ownSample(int fs_code, index_t id) const;

  /**
     \brief assigns new tag newTag to all samples of functionspace with a positive
     value of mask for any its sample point.

  */
  ESCRIPT_DLL_API
  virtual void setTags(const int functionSpaceType, const int newTag, const escript::Data& mask) const;

  /**
     \brief
     returns true if data on this domain and a function space of type functionSpaceCode has to
     considered as cell centered data.

     This has to be implemented by the actual Domain adapter.
  */
  ESCRIPT_DLL_API
  virtual bool isCellOriented(int functionSpaceCode) const;

  /**
     \brief
      Returns a status indicator of the domain. The status identifier should be unique over 
      the live time if the object but may be updated if changes to the domain happen, e.g. 
      modifications to its geometry. 

     This has to be implemented by the actual Domain adapter.
  */
  ESCRIPT_DLL_API
  virtual StatusType getStatus() const;

  /**
     \brief
     Throw a standard exception. This function is called if any attempt 
     is made to use a base class function.
  */
  ESCRIPT_DLL_API
  void throwStandardException(const std::string& functionName) const;

  /**
        \brief
                  return the number of tags in use and a pointer to an array with the number of tags in use
  */
  ESCRIPT_DLL_API
  virtual int getNumberOfTagsInUse(int functionSpaceCode) const;

  ESCRIPT_DLL_API
  virtual const int* borrowListOfTagsInUse(int functionSpaceCode) const;

  /**
    \brief Checks if this domain allows tags for the specified functionSpaceCode.
  */
  ESCRIPT_DLL_API
  virtual bool canTag(int functionspacecode) const;

   /**
   \brief returns the approximation order used for a function space functionSpaceCode
   */

   ESCRIPT_DLL_API
   virtual int getApproximationOrder(const int functionSpaceCode) const;

   ESCRIPT_DLL_API
   virtual bool supportsContactElements() const;

   

  /**
   * \brief true if this domain can handle to specified tuple of filter options.
  */
  ESCRIPT_DLL_API
  virtual bool supportsFilter(const boost::python::tuple& t) const;
  
  /**
   * \brief Fills the data object with filtered random values 
  */ 
  ESCRIPT_DLL_API
  virtual escript::Data randomFill(long seed, const boost::python::tuple& filter) const;   
   
 protected:

 private:

   // buffer for coordinates used by function spaces
   BufferMapType m_x_buffer;

   // buffer for normal vectors used by function spaces
   BufferMapType m_normal_buffer;

   // buffer for normal element size used by function spaces
   BufferMapType m_size_buffer;

};

} // end of namespace

#endif
