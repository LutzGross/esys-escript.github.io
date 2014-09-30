
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


#ifndef __ESCRIPT_FUNCTIONSPACE_H__
#define __ESCRIPT_FUNCTIONSPACE_H__

#include "system_dep.h"

#include "AbstractDomain.h"

#include <boost/python/list.hpp>
#include <list>

namespace escript {

//
// Forward declaration for class Data.
class Data;

class ESCRIPT_DLL_API FunctionSpace 
{
public:
    FunctionSpace();

    FunctionSpace(const_Domain_ptr domain, int functionSpaceType);

    FunctionSpace(const FunctionSpace& other);

    /**
      \brief Returns the function space type code.
      \note The meaning of the code depends on the domain object the
            FunctionSpace is built on.
    */
    int getTypeCode() const;

    /**
      \brief Returns the function space domain.
    */
    const_Domain_ptr getDomain() const;

    /**
      \brief Return the function space domain. Internal use only! This gets
             around some python difficulties by casting away the const.
             Do not use this in C++. 
    */
    Domain_ptr getDomainPython() const;

    /**
      \brief Returns true if this function space support tags
    */
    bool canTag() const;

    /**
      \brief Returns the approximation order used for this function space
    */
    int getApproximationOrder() const;

    /**
      \brief Assigns new tag newTag to all samples with a positive value of
             mask for any of its sample points.
    */
    void setTags(const int newTag, const escript::Data& mask) const;

    void setTagsByString(const std::string& name, const escript::Data& mask) const;
  
    /**
      \brief Returns the shape of the data needed to represent the function
             space.
    */
    std::pair<int,dim_t> getDataShape() const;

    /**
      \brief Comparison operator.  Returns true if function spaces are equal.
             (i.e. same domain and same function space type)
    */
    bool operator==(const FunctionSpace& other) const;

    bool operator!=(const FunctionSpace& other) const;

    /**
      \brief Returns a text description of the function space.
    */
    std::string toString() const;

    /**
      \brief Returns the tag associated with the given sample number.
    */
    int getTagFromSampleNo(dim_t sampleNo) const;

    /**
      \brief Returns the tag associated with the given data-point number.
    */
    int getTagFromDataPointNo(dim_t dataPointNo) const;

    /**
      \brief Returns the reference number associated with the given data-point
             number.
    */
    dim_t getReferenceIDFromDataPointNo(dim_t dataPointNo) const;

    /**
      \brief Returns the reference number associated with the given sample
             number. This function is not efficient. It is better to first call
             `borrowSampleReferenceIDs` and then when iterating over sampleNo
             to use sampleNo as an offset.
    */
    inline
    dim_t getReferenceIDOfSample(dim_t sampleNo) const {
        return borrowSampleReferenceIDs()[sampleNo];
    }

    /**
      \brief Does this process own the sample? For non-MPI builds will always
             return true
    */
    inline
    bool ownSample(dim_t sampleNo) const {
        return m_domain->ownSample(m_functionSpaceType, sampleNo);
    }

    /**
      \brief Returns a borrowed reference to the list of sample reference IDs
    */
    const dim_t* borrowSampleReferenceIDs() const;

    /**
      \brief Returns the spatial locations of the data points.
    */
    escript::Data getX() const;
 
    /**
      \brief Returns the surface normal field.
    */
    escript::Data getNormal() const;

    /**
      \brief Returns the sample size (e.g. the diameter of elements, radius
             of particles).
    */
    escript::Data getSize() const;

    /**
      \brief Returns the number of samples.
    */
    inline
    dim_t getNumSamples() const { return getDataShape().second; }

    /**
      \brief Returns the number of data points per sample.
    */
    inline
    int getNumDPPSample() const { return getNumDataPointsPerSample(); }

    inline
    int getNumDataPointsPerSample() const { return getDataShape().first; }

    /**
      \brief Return the number of spatial dimensions of the underlying domain.
    */
    inline
    int getDim() const { return getDomain()->getDim(); }

    /**
      \brief Returns a list of the tags used in this function space
    */
    boost::python::list getListOfTags() const;

    /**
      \brief Returns an stl list of the tags used in this function space
    */
    std::list<int> getListOfTagsSTL() const;

    /**
      \brief Returns the number of tags in use.
    */
    int getNumberOfTagsInUse() const;

    const int* borrowListOfTagsInUse() const;

    inline
    bool probeInterpolation(const FunctionSpace& other) const
    {
        if (*this == other)
            return true;
        const_Domain_ptr domain(getDomain());
        if (*domain == *other.getDomain()) {
            return domain->probeInterpolationOnDomain(
                        getTypeCode(), other.getTypeCode());
        }
        return domain->probeInterpolationACross(
                getTypeCode(), *(other.getDomain()), other.getTypeCode());
    }

private:
    /**
    \brief Assignment operator. 
           This method is only defined (private) to prevent people from
           using it.
    */
    FunctionSpace& operator=(const FunctionSpace& other);

    const_Domain_ptr m_domain;

    int m_functionSpaceType;
};

bool canInterpolate(FunctionSpace src, FunctionSpace dest);

} // end of namespace

#endif // __ESCRIPT_FUNCTIONSPACE_H__

