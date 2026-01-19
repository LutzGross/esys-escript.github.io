
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __ESCRIPT_NULLDOMAIN_H__
#define __ESCRIPT_NULLDOMAIN_H__

#include "system_dep.h"
#include "AbstractDomain.h"

namespace escript {

/**
   \brief
   NullDomain provides a null value for domain. Needed for the construction
   of a default FunctionSpace.

   Description:
   NullDomain provides a null value for domain. Needed for the construction
   of a default FunctionSpace. Inherits from AbstractDomain and overrides its
   methods.
   This domain supports a single type of FunctionSpace for which canTag is true.
   This compromise is needed to allow the default constructor of DataTagged to
   have a FunctionSpace which supports tagging.
   See notes on the borrowListOfTagsInUse() method.
*/
class ESCRIPT_DLL_API NullDomain : public AbstractDomain
{
private:
    static int NullDomainFS;
    static DataTypes::dim_t referenceID;
        
public:
    NullDomain() : AbstractDomain(makeInfo(MPI_COMM_NULL)) {}

    virtual bool isValidFunctionSpaceType(int fsCode) const {
        return fsCode==NullDomainFS;
    }

    virtual std::string getDescription() const;

    virtual std::string functionSpaceTypeAsString(int functionSpaceType) const;

    virtual void interpolateOnDomain(escript::Data& target,const escript::Data& source) const;

    virtual bool probeInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const;

    virtual void interpolateAcross(escript::Data& target, const escript::Data& source) const;

    virtual bool probeInterpolationAcross(int, const AbstractDomain&, int) const {
        return false;
    }

    virtual int getContinuousFunctionCode() const { return NullDomainFS; }
    virtual int getFunctionCode() const { return NullDomainFS; }
    virtual int getFunctionOnBoundaryCode() const { return NullDomainFS; }
    virtual int getFunctionOnContactZeroCode() const { return NullDomainFS; }
    virtual int getFunctionOnContactOneCode() const { return NullDomainFS; }
    virtual int getSolutionCode() const { return NullDomainFS; }
    virtual int getReducedSolutionCode() const { return NullDomainFS; }
    virtual int getDiracDeltaFunctionsCode() const { return NullDomainFS; }

    virtual std::pair<int,DataTypes::dim_t> getDataShape(int functionSpaceCode) const;

    virtual int getTagFromSampleNo(int, DataTypes::index_t) const { return 1; }

    virtual const DataTypes::dim_t* borrowSampleReferenceIDs(int) const { return &referenceID; }

    virtual int getDim() const { return 1; }

    virtual bool operator==(const AbstractDomain& other) const;

    virtual bool operator!=(const AbstractDomain& other) const {
        return !(*this==other);
    }

    virtual void write(const std::string& filename) const;

    virtual void dump(const std::string& filename) const;

    virtual void setTagMap(const std::string& name,  int tag);

    virtual int getTag(const std::string& name) const;

    virtual bool canTag(int) const { return true; }

    virtual std::string showTagNames() const;

    virtual int getNumberOfTagsInUse(int) const { return 1; }

    virtual const int* borrowListOfTagsInUse(int functionSpaceCode) const;

    virtual void setTags(int functionSpaceType, int newTag, const escript::Data& mask) const;

    bool supportsContactElements() const { return false; }

    virtual void setNewX(const escript::Data& arg);

    virtual signed char preferredInterpolationOnDomain(
                                      int functionSpaceType_source,
                                      int functionSpaceType_target) const;

    virtual bool commonFunctionSpace(const std::vector<int>& fs,
                                     int& resultcode) const;
    virtual bool isCellOriented(int functionSpaceCode) const;
    virtual int getApproximationOrder(const int functionSpaceCode) const;

    virtual escript::Data getX() const;
#ifdef ESYS_HAVE_BOOST_NUMPY
    virtual boost::python::numpy::ndarray getNumpyX() const;
#endif
    virtual escript::Data getNormal() const;
    virtual escript::Data getSize() const;
    virtual void setToX(escript::Data& out) const;
    virtual void setToNormal(escript::Data& out) const;
    virtual void setToSize(escript::Data& out) const;
    virtual void setToGradient(escript::Data& grad, const escript::Data& arg) const;
    virtual bool ownSample(int fs_code, DataTypes::index_t id) const;
    virtual escript::Data randomFill(const DataTypes::ShapeType& shape,
                                     const FunctionSpace& what, long seed,
                                     const boost::python::tuple& filter) const;
};

} // end of namespace

#endif // __ESCRIPT_NULLDOMAIN_H__

