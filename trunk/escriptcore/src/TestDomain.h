
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

#ifndef __ESCRIPT_TESTDOMAIN_H__
#define __ESCRIPT_TESTDOMAIN_H__

#include "system_dep.h"

#include "NullDomain.h"

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
    virtual void MPIBarrier() const;
    virtual bool onMasterProcessor() const;
    virtual MPI_Comm getMPIComm() const;

    virtual bool isValidFunctionSpaceType(int functionSpaceType) const;

    virtual std::string getDescription() const;

    virtual std::string functionSpaceTypeAsString(int functionSpaceType) const;

    virtual void interpolateOnDomain(escript::Data& target,const escript::Data& source) const;

    virtual bool probeInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const;

    bool commonFunctionSpace(const std::vector<int>& fs, int& resultcode) const;

    virtual escript::Data getX() const;

    virtual void interpolateACross(escript::Data& target, const escript::Data& source) const;

    virtual bool probeInterpolationACross(int functionSpaceType_source,const AbstractDomain& targetDomain, int functionSpaceType_target) const;

    virtual int getDefaultCode() const;
    virtual int getContinuousFunctionCode() const;
    virtual int getFunctionCode() const;
    virtual int getFunctionOnBoundaryCode() const;
    virtual int getFunctionOnContactZeroCode() const;
    virtual int getFunctionOnContactOneCode() const;
    virtual int getSolutionCode() const;
    virtual int getReducedSolutionCode() const;
    virtual int getDiracDeltaFunctionsCode() const;

    virtual std::pair<int,dim_t> getDataShape(int functionSpaceCode) const;

    virtual int getTagFromSampleNo(int functionSpaceType, dim_t sampleNo) const;

    virtual const dim_t* borrowSampleReferenceIDs(int functionSpaceType) const;

    virtual int getDim() const;

    virtual bool operator==(const AbstractDomain& other) const;

    virtual bool operator!=(const AbstractDomain& other) const;

    virtual bool canTag(int functionSpaceCode) const;

    virtual int getNumberOfTagsInUse(int functionSpaceCode) const;

    virtual const int* borrowListOfTagsInUse(int functionSpaceCode) const;

    virtual escript::Data randomFill(const DataTypes::ShapeType& shape,
                                     const FunctionSpace& what, long seed,
                                     const boost::python::tuple& filter) const;

private:
    dim_t m_samples;       // number of samples
    int m_dpps;            // data points per sample
    int m_dpsize;          // how big are the datapoints?
    dim_t* m_samplerefids; // sample reference ids
};

ESCRIPT_DLL_API
FunctionSpace getTestDomainFunctionSpace(int dpps, dim_t samples, int dpsize);

} // end of namespace

#endif // __ESCRIPT_TESTDOMAIN_H__

