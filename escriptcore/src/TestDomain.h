
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

#ifndef __ESCRIPT_TESTDOMAIN_H__
#define __ESCRIPT_TESTDOMAIN_H__

#include "system_dep.h"
#include "NullDomain.h"
#include <vector>

namespace escript {

/**
   \brief
   (Testing use only) Provides a domain to wrap a collection of values.

   This domain provides more functionality than NullDomain in that it can
   store varying numbers of samples and points per sample.

   It currently supports a single function space which does not support tagging.
   No effort has been made to make this work with MPI

   \warning This class exists to support testing and should not be used
   as a general domain without ensuring that it works the way you expect.
   Also, other doco says that this class can be removed without notice.

*/
class ESCRIPT_DLL_API TestDomain : public NullDomain
{
public:
    TestDomain(int pointspersample, int numsamples, int dpsize=1);

    virtual ~TestDomain();

    virtual int getMPISize() const;
    virtual int getMPIRank() const;

    virtual bool isValidFunctionSpaceType(int functionSpaceType) const;

    virtual std::string getDescription() const;

    virtual std::string functionSpaceTypeAsString(int functionSpaceType) const;

    virtual void interpolateOnDomain(escript::Data& target,const escript::Data& source) const;

    virtual bool probeInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const;

    bool commonFunctionSpace(const std::vector<int>& fs, int& resultcode) const;

    virtual escript::Data getX() const;

    virtual void interpolateAcross(escript::Data& target, const escript::Data& source) const;

    virtual bool probeInterpolationAcross(int functionSpaceType_source,const AbstractDomain& targetDomain, int functionSpaceType_target) const;

    virtual int getDefaultCode() const;
    virtual int getContinuousFunctionCode() const;
    virtual int getFunctionCode() const;
    virtual int getFunctionOnBoundaryCode() const;
    virtual int getFunctionOnContactZeroCode() const;
    virtual int getFunctionOnContactOneCode() const;
    virtual int getSolutionCode() const;
    virtual int getReducedSolutionCode() const;
    virtual int getDiracDeltaFunctionsCode() const;

    virtual std::pair<int,DataTypes::dim_t> getDataShape(int functionSpaceCode) const;

    virtual int getTagFromSampleNo(int functionSpaceType, DataTypes::index_t sampleNo) const;

    virtual const DataTypes::dim_t* borrowSampleReferenceIDs(int functionSpaceType) const;

    virtual int getDim() const;

    virtual bool operator==(const AbstractDomain& other) const;

    virtual bool operator!=(const AbstractDomain& other) const;

    virtual bool canTag(int functionSpaceCode) const;

    virtual int getNumberOfTagsInUse(int functionSpaceCode) const;

    virtual const int* borrowListOfTagsInUse(int functionSpaceCode) const;

    virtual escript::Data randomFill(const DataTypes::ShapeType& shape,
                                     const FunctionSpace& what, long seed,
                                     const boost::python::tuple& filter) const;
				     
    void addUsedTag(int t);
    void clearUsedTags();
    void assignTags(std::vector<int> t);
    void resetTagAssignments();

private:
    DataTypes::dim_t m_totalsamples;	// samples in all worlds  
    DataTypes::dim_t m_samples;       // number of samples
    DataTypes::dim_t m_originsample;
    DataTypes::dim_t m_endsample;

    
    int m_dpps;            // data points per sample
    int m_dpsize;          // how big are the datapoints?
    DataTypes::dim_t* m_samplerefids; // sample reference ids
    
    std::vector<int> mytags;
    std::vector<int> tag_assignment; 	// which tag is assigned to each sample
				// to make testing easier, the tags in use list is
				// controlled separately
};

ESCRIPT_DLL_API
FunctionSpace getTestDomainFunctionSpace(int dpps, DataTypes::dim_t samples, int dpsize);

} // end of namespace

#endif // __ESCRIPT_TESTDOMAIN_H__

