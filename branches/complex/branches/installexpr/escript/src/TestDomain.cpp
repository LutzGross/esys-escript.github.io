
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


#include "DomainException.h"
#include "TestDomain.h" 
#include "Data.h"
#include "Utils.h"	// for MPI functions

namespace escript {

namespace {
const int defaultList[1]={0};		// an array to return in borrowListOfTagsInUse();
const int TestDomainFS=1;		// Null domains only support 1 functionspace type.
			// The choice of =1 as the value is arbitrary
}


TestDomain::TestDomain(int pointspersample, int numsamples, int dpsize)
	: m_samples(numsamples), m_dpps(pointspersample), m_dpsize(dpsize)
{
#ifdef ESYS_MPI
    int world=getMPISizeWorld();
    int rank=getMPIRankWorld();
    m_samples/=world;
    if (rank<(numsamples%world))
    {
	m_samples++;
    }
#endif
    m_samplerefids=new int[numsamples];
    for (int i=0;i<numsamples;++i)
    {
        m_samplerefids[i]=i+10;		// the +10 is arbitrary. 
    }					// so these ids look different from others
}
 

TestDomain::~TestDomain()
{    
    delete[] m_samplerefids;
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

void TestDomain::interpolateACross(Data& target, const Data& source) const
{
   throw DomainException("Error - interpolation to the TestDomain not supported.");
}

bool TestDomain::probeInterpolationACross(int functionSpaceType_source,const AbstractDomain& targetDomain, int functionSpaceType_target) const
{
   return false;
}

ESCRIPT_DLL_API
bool TestDomain::commonFunctionSpace(const std::vector<int>& fs, int& resultcode) const
{
    for (int i=0;i<fs.size();++i)
    {
	if (fs[i]!=TestDomainFS)
	{
		return false;
	}
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

std::pair<int,int> TestDomain::getDataShape(int functionSpaceCode) const
{
  return std::pair<int,int>(m_dpps,m_samples);
}

int TestDomain::getTagFromSampleNo(int functionSpaceType, int sampleNo) const
{
  //
  // return an arbitrary value
  // - In this case I have chosen to return the default tag
  return 0; 
}

const int* TestDomain::borrowSampleReferenceIDs(int functionSpaceType) const
{
  //
  // return an arbitrary value
  return m_samplerefids;
}

int TestDomain::getDim() const
{
  //
  // return an arbitrary value
  // Since this domain doesn't really have structure I guess 1 seems sensible
  return 1; 
}

bool TestDomain::operator==(const AbstractDomain& other) const
{
  const TestDomain* temp=dynamic_cast<const TestDomain*>(&other);
  if (temp!=0) {
    return true;
  } else {
    return false;
  }
}

bool TestDomain::operator!=(const AbstractDomain& other) const
{
  return(!(*this==other));
}



bool TestDomain::canTag(int functionSpaceCode) const
{
  return true;
}

int TestDomain::getNumberOfTagsInUse(int functionSpaceCode) const
{
  return 1;	// this is not arbitrary. It allows us to report that the default tag is in use
}

const int* TestDomain::borrowListOfTagsInUse(int functionSpaceCode) const
{
  return defaultList;
}

escript::Data TestDomain::getX() const
{
  if (m_dpsize<2)
  {  
  Data res(0,DataTypes::scalarShape,FunctionSpace( getPtr(), getDefaultCode()),true);
  DataTypes::ValueType& vec=res.getReady()->getVectorRW();
  for (int i=0;i<m_samples;++i)
  {
    for (int j=0;j<m_dpps;++j)
    {
	vec[i*m_dpps+j]=i+(1.0*j)/m_dpps;
    }
  }
  return res;
  }
  DataTypes::ShapeType p;
  p.push_back(m_dpsize);
  Data res(0,p,FunctionSpace( getPtr(), getDefaultCode()),true);
  DataTypes::ValueType& vec=res.getReady()->getVectorRW();
  double majorstep=double(1)/m_dpps;
  double minorstep=majorstep*0.9/m_dpsize;
  for (int i=0;i<m_samples;++i)
  {
    for (int j=0;j<m_dpps;++j)
    {
        for (int k=0;k<m_dpsize;++k)
	{
	    vec[i*m_dpsize*m_dpps+j*m_dpsize+k]=i+j*majorstep+k*minorstep;
	}
    }
  }
  return res;
  
  
}

FunctionSpace
getTestDomainFunctionSpace(int dpps, int samples, int dpsize)
{
    TestDomain* td=new TestDomain(dpps, samples, dpsize);
    Domain_ptr p=Domain_ptr(td);
    return FunctionSpace(p, td->getDefaultCode());
}

}  // end of namespace
