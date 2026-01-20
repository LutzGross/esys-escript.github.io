
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <sstream>

#include "TestDomain.h"
#include "Data.h"
#include "DomainException.h"
#include "Random.h"
#include "Utils.h" // for MPI functions

namespace escript {

namespace {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-const-variable"
const int defaultList[1]={0}; // an array to return in borrowListOfTagsInUse();
#pragma clang diagnostic pop
const int TestDomainFS=1;     // Null domains only support 1 functionspace type.
                              // The choice of =1 as the value is arbitrary
}

TestDomain::TestDomain(int pointspersample, int numsamples, int dpsize)
        : m_totalsamples(numsamples), m_samples(numsamples), m_dpps(pointspersample), m_dpsize(dpsize)
{
    // Override NullDomain's MPI_COMM_NULL with MPI_COMM_WORLD for testing
    m_mpiInfo = makeInfo(MPI_COMM_WORLD);
    int world=getMPISizeWorld();
    int rank=getMPIRankWorld();
    m_samples/=world;
    m_originsample=rank*m_samples;	// how many samples have gone to earlier ranks before we start    
    if (world > 1 && rank < numsamples%world) {
        m_samples++;
    }
    if ((world > 1) && (numsamples%world))
    {
	m_originsample+=(rank<numsamples%world?rank:numsamples%world);      
    }
    m_endsample=m_originsample+numsamples/world;
    if ((world > 1) && (rank < numsamples%world))
    {
	m_endsample++;      
    }
    m_endsample--;		// so that the end_sample is inclusive
    m_samplerefids=new DataTypes::dim_t[numsamples];
    for (DataTypes::dim_t i=0; i<numsamples; ++i) {
        m_samplerefids[i]=i+10; // the +10 is arbitrary.
    }                           // so these ids look different from others
    mytags.push_back(0);
    resetTagAssignments();
}

TestDomain::~TestDomain()
{
    delete[] m_samplerefids;
}

int TestDomain::getMPISize() const
{
    return getMPISizeWorld();
}

int TestDomain::getMPIRank() const
{
    return getMPIRankWorld();
}


bool TestDomain::isValidFunctionSpaceType(int functionSpaceType) const
{
    return (functionSpaceType==TestDomainFS);
}

std::string TestDomain::getDescription() const
{
    return "TestDomain";
}

std::string TestDomain::functionSpaceTypeAsString(int functionSpaceType) const
{
    return "Default_FunctionSpace";
}

void TestDomain::interpolateOnDomain(Data& target,const Data& source) const
{
    if (source.getFunctionSpace().getDomain().get()!=this)
        throw DomainException("Error - Illegal domain of interpolant.");
    if (target.getFunctionSpace().getDomain().get()!=this)
        throw DomainException("Error - Illegal domain of interpolation target.");
    target=source;
}

bool TestDomain::probeInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const
{
    if ((functionSpaceType_source!=functionSpaceType_target) || (functionSpaceType_target!=TestDomainFS))
    {
        throw DomainException("Error - Illegal function type for TestDomain.");
    }
    return true;
}

void TestDomain::interpolateAcross(Data& target, const Data& source) const
{
    throw DomainException("Error - interpolation to the TestDomain not supported.");
}

bool TestDomain::probeInterpolationAcross(int functionSpaceType_source,const AbstractDomain& targetDomain, int functionSpaceType_target) const
{
    return false;
}

ESCRIPT_DLL_API
bool TestDomain::commonFunctionSpace(const std::vector<int>& fs, int& resultcode) const
{
    for (int i=0;i<fs.size();++i) {
        if (fs[i] != TestDomainFS)
            return false;
    }
    resultcode=TestDomainFS;
    return true;
}

int TestDomain::getDefaultCode() const
{
    return TestDomainFS;
}

int TestDomain::getContinuousFunctionCode() const
{
    return TestDomainFS;
}

int TestDomain::getFunctionCode() const
{
    return TestDomainFS;
}

int TestDomain::getFunctionOnBoundaryCode() const
{
    return TestDomainFS;
}

int TestDomain::getFunctionOnContactZeroCode() const
{
    return TestDomainFS;
}

int TestDomain::getFunctionOnContactOneCode() const
{
    return TestDomainFS;
}

int TestDomain::getSolutionCode() const
{
    return TestDomainFS;
}

int TestDomain::getReducedSolutionCode() const
{
    return TestDomainFS;
}

int TestDomain::getDiracDeltaFunctionsCode() const
{
    return TestDomainFS;
}

std::pair<int,DataTypes::dim_t> TestDomain::getDataShape(int functionSpaceCode) const
{
    return std::pair<int,DataTypes::dim_t>(m_dpps,m_samples);
}

int TestDomain::getTagFromSampleNo(int functionSpaceType, DataTypes::index_t sampleNo) const
{
    if (sampleNo>=tag_assignment.size())
    {
	std::ostringstream oss;
	oss << "invalid sample number " << sampleNo << " of " << tag_assignment.size();
	throw DataException(oss.str());
    }
    return tag_assignment[sampleNo];
}

const DataTypes::dim_t* TestDomain::borrowSampleReferenceIDs(int functionSpaceType) const
{
    return m_samplerefids;
}

int TestDomain::getDim() const
{
    return 1;
}

bool TestDomain::operator==(const AbstractDomain& other) const
{
    const TestDomain* temp=dynamic_cast<const TestDomain*>(&other);
    return (temp != NULL);
}

bool TestDomain::operator!=(const AbstractDomain& other) const
{
  return !(*this==other);
}

bool TestDomain::canTag(int functionSpaceCode) const
{
    return true;
}

int TestDomain::getNumberOfTagsInUse(int functionSpaceCode) const
{
    return mytags.size();
}

escript::Data TestDomain::getX() const
{
    if (m_dpsize<2) {
        Data res(0,DataTypes::scalarShape,FunctionSpace( getPtr(), getDefaultCode()),true);
        DataTypes::RealVectorType& vec=res.getReady()->getVectorRW();
        for (DataTypes::dim_t i=0; i<m_samples; ++i) {
            for (int j=0; j<m_dpps; ++j) {
                vec[i*m_dpps+j]=m_originsample+i+(1.0*j)/m_dpps;
            }
        }
        return res;
    }
    DataTypes::ShapeType p;
    p.push_back(m_dpsize);
    Data res(0,p,FunctionSpace( getPtr(), getDefaultCode()),true);
    DataTypes::RealVectorType& vec=res.getReady()->getVectorRW();
    double majorstep=double(1)/m_dpps;
    double minorstep=majorstep*0.9/m_dpsize;
    for (DataTypes::dim_t i=0; i<m_samples; ++i) {
        for (int j=0; j<m_dpps; ++j) {
            for (int k=0; k<m_dpsize; ++k) {
                vec[i*m_dpsize*m_dpps+j*m_dpsize+k]=(i+m_originsample)+j*majorstep+k*minorstep;
            }
        }
    }
    return res;
}

escript::Data TestDomain::randomFill(const DataTypes::ShapeType& shape,
       const FunctionSpace& what, long seed, const boost::python::tuple& filter) const
{
    escript::Data towipe(0, shape, what, true);
    // since we just made this object, no sharing is possible and we don't
    // need to check for exclusive write
    escript::DataTypes::RealVectorType& dv=towipe.getExpandedVectorReference();
    const size_t dvsize=dv.size();
    escript::randomFillArray(seed, &(dv[0]), dvsize, getMPI());
    return towipe;
}

void TestDomain::addUsedTag(int t)
{
    for (auto i=mytags.begin();i!=mytags.end();++i)
    {
	if (*i==t)
	{
	    return;
	}
    }
    mytags.push_back(t);
}

void TestDomain::clearUsedTags()
{
    mytags.clear();
    mytags.push_back(0);
}

const int* TestDomain::borrowListOfTagsInUse(int functionSpaceCode) const
{
    return &mytags[0];
}

void TestDomain::assignTags(std::vector<int> t)
{
    if (t.size()!=m_totalsamples)
    {
	throw DataException("Programming error - Tag vector must be the same size as the number of samples.");
    }
    tag_assignment=std::vector<int>(m_samples);
    for (int i=m_originsample;i<=m_endsample;++i)
    {
	tag_assignment[i-m_originsample]=t[i];
    }
}


void TestDomain::resetTagAssignments()
{
    tag_assignment=std::vector<int>(m_samples);
    for (size_t i=0;i<m_samples;++i)
    {
	tag_assignment[i]=0;
    }
}


FunctionSpace getTestDomainFunctionSpace(int dpps, DataTypes::dim_t samples, int dpsize)
{
    TestDomain* td=new TestDomain(dpps, samples, dpsize);
    Domain_ptr p=Domain_ptr(td);
    return FunctionSpace(p, td->getDefaultCode());
}





}  // end of namespace

