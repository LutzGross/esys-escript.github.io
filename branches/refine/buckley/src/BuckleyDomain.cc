
/*******************************************************
*
* Copyright (c) 2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include <iostream>
#include "BuckleyDomain.h"
#include "BuckleyException.h"

#include <escript/FunctionSpace.h>
#include <escript/FunctionSpaceFactory.h>
using namespace buckley;



namespace
{
  // Need to put some thought into how many spaces we actually need
    const int initialspace=1;   
    const int invalidspace=-1;
    
    const int notIMPLEMENTED=-1;
    
    const int ctsfn=initialspace;
}


BuckleyDomain::BuckleyDomain(double x, double y, double z)
:ot(x,y,z),modified(true),leaves(0),numpts(0),samplerefids(0)
{
}

BuckleyDomain::~BuckleyDomain()
{
    if (leaves)
    {
       delete[] leaves;
    }
    if (samplerefids)
    {
       delete[] samplerefids;
    }
}

bool BuckleyDomain::isValidFunctionSpaceType(int functionSpaceType) const
{
    return functionSpaceType==initialspace;
}

escript::Data BuckleyDomain::getX() const
{
   if (modified)	// is the cached data we have about this domain stale?
   {
      processMods();
   }   
   return escript::continuousFunction(*this).getX();
}


// typedef struct 
// {
//     escript::Vector* vec;  
//     size_t chunksize;
//     unkid id;
//   
// } idstruct;
// 
// 
// void copycoords(const OctCell& c, void* v)
// {
//     idstruct* s=reinterpret_cast<idstruct*>(v);
// #ifdef OMP    
//     if (id / s->chunksize==omp_get_num_thread())
//     {
// #endif      
//         s->vec[id++]
// 
// 
// #ifdef OMP    
//     }
// #endif
// }

bool BuckleyDomain::operator==(const AbstractDomain& other) const
{
  if (dynamic_cast<const BuckleyDomain*>(&other)==0)
  {
      return false;
  }
  return this==&(other);
}

bool BuckleyDomain::operator!=(const AbstractDomain& other) const
{
  if (dynamic_cast<const BuckleyDomain*>(&other)==0)
  {
      return true;
  }
  return this!=&(other);
}


const int* BuckleyDomain::borrowSampleReferenceIDs(int functionSpaceType) const
{
  return samplerefids;
}


void BuckleyDomain::processMods() const
{
      if (leaves!=0)
      {
	  delete [] leaves;
	  delete [] samplerefids;
      }
      leaves=ot.process(numpts);
      samplerefids=new int[numpts];
      for (int i=0;i<numpts;++i)
      {
	  samplerefids[i]=i;
      }
      modified=false;
}

void BuckleyDomain::setToX(escript::Data& arg) const
{
   const BuckleyDomain& argDomain=dynamic_cast<const BuckleyDomain&>(*(arg.getFunctionSpace().getDomain()));
   if (argDomain!=*this) 
      throw BuckleyException("Error - Illegal domain of data point locations");
   
   if (arg.getFunctionSpace().getTypeCode()==invalidspace)
   {
       throw BuckleyException("Error - Invalid functionspace for coordinate loading");
   }
   if (modified)	// is the cached data we have about this domain stale?
   {
      processMods();
   }
   
   if (arg.getFunctionSpace().getTypeCode()==getContinuousFunctionCode())	// values on nodes
   {
      escript::DataTypes::ValueType::ValueType vec=(&arg.getDataPointRW(0, 0));     
      int i;	// one day OMP-3 will be default
// Why isn't this parallelised??
// well it isn't threadsafe the way it is written
// up to 8 cells could touch a point and try to write its coords
// This can be fixed but I haven't done it yet
//       #pragma omp parallel for private(i) schedule(static)
      for (i=0;i<ot.leafCount();++i)
      {
	  const OctCell& oc=*leaves[i];
	  const double dx=oc.sides[0]/2;
	  const double dy=oc.sides[1]/2;
	  const double dz=oc.sides[2]/2;
	  const double cx=oc.centre[0];
	  const double cy=oc.centre[1];
	  const double cz=oc.centre[2];
	  
	  // So this section is likely to be horribly slow
	  // Once the system is working unroll this
	  for (unsigned int j=0;j<8;++j)
	  {
	      unkid destpoint=oc.leafinfo->pmap[j];
	      if (destpoint<2)
	      {
		  continue;	// it's a hanging node so don't produce
	      }
	      destpoint-=2;
	      double x=cx;
	      double y=cy;
	      double z=cz;
	      if (j>3)
	      {
		  z+=dz;
	      }
	      else
	      {
		  z-=dz;
	      }
	      switch (j)
	      {
		case 0:
		case 3:
		case 4:
		case 7: x-=dx; break;
		  
		default:  x+=dx;
		
	      }
	      
	      switch (j)
	      {
		case 0:
		case 1:
		case 4:
		case 5: y-=dy; break;
		
		default:  y+=dy;
	      }	      
	      
	      vec[destpoint*3]=x;
	      vec[destpoint*3+1]=y;
	      vec[destpoint*3+2]=z;
	  }
	
      }
   }
   else		// so it's not on the nodes and we need to do some interpolation
   {
      throw BuckleyException("Please don't use other functionspaces on this yet.");
     
   }
     
//    Dudley_Mesh* mesh=m_dudleyMesh.get();
//    // in case of values node coordinates we can do the job directly:
//    if (arg.getFunctionSpace().getTypeCode()==Nodes) {
//       escriptDataC _arg=arg.getDataC();
//       Dudley_Assemble_NodeCoordinates(mesh->Nodes,&_arg);
//    } else {
//       escript::Data tmp_data=Vector(0.0,continuousFunction(*this),true);
//       escriptDataC _tmp_data=tmp_data.getDataC();
//       Dudley_Assemble_NodeCoordinates(mesh->Nodes,&_tmp_data);
//       // this is then interpolated onto arg:
//       interpolateOnDomain(arg,tmp_data);
//    }
//    checkDudleyError();
}


void BuckleyDomain::refineAll(unsigned min_depth)
{
    modified=true;
    ot.allSplit(min_depth);  
}



std::string BuckleyDomain::getDescription() const
{
    return "Buckley domain";
}

int BuckleyDomain::getContinuousFunctionCode() const
{
    return ctsfn;    
}

int BuckleyDomain::getReducedContinuousFunctionCode() const
{
    return initialspace;
}

int BuckleyDomain::getFunctionCode() const
{
    return invalidspace;  
}

int BuckleyDomain::getReducedFunctionCode() const
{
    return invalidspace;  
}

int BuckleyDomain::getFunctionOnBoundaryCode() const
{
    return invalidspace;  
}

int BuckleyDomain::getReducedFunctionOnBoundaryCode() const
{
    return invalidspace;  
}

int BuckleyDomain::getFunctionOnContactZeroCode() const
{
    return invalidspace;  
}

int BuckleyDomain::getReducedFunctionOnContactZeroCode() const
{
    return invalidspace;
}

int BuckleyDomain::getFunctionOnContactOneCode() const
{
    return invalidspace;
}

int BuckleyDomain::getReducedFunctionOnContactOneCode() const
{
    return invalidspace;  
}

int BuckleyDomain::getSolutionCode() const
{
    return initialspace;  
}

int BuckleyDomain::getReducedSolutionCode() const
{
    return initialspace;  
}

// hahaha - no 
int BuckleyDomain::getDiracDeltaFunctionsCode() const
{
    return invalidspace  ;
}

int BuckleyDomain::getSystemMatrixTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const
{
    return notIMPLEMENTED;  
}

int BuckleyDomain::getTransportTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const
{
    return notIMPLEMENTED;  
}

void BuckleyDomain::setToIntegrals(std::vector<double>& integrals,const escript::Data& arg) const
{
    throw BuckleyException("Not Implemented");
}

void BuckleyDomain::addPDEToSystem(
                     escript::AbstractSystemMatrix& mat, escript::Data& rhs,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C, 
                     const escript::Data& D, const escript::Data& X, const escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact, const escript::Data& y_contact, 
                     const escript::Data& d_dirac, const escript::Data& y_dirac) const
{
    throw BuckleyException("Not Implemented");
}

void BuckleyDomain::addPDEToRHS(escript::Data& rhs,
                     const escript::Data& X, const escript::Data& Y,
                     const escript::Data& y, const escript::Data& y_contact, const escript::Data& y_dirac) const
{
    throw BuckleyException("Not Implemented");
}

void BuckleyDomain::addPDEToTransportProblem(
                     escript::AbstractTransportProblem& tp, escript::Data& source, 
                     const escript::Data& M,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C,const  escript::Data& D,
                     const  escript::Data& X,const  escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact,const escript::Data& y_contact,
                     const escript::Data& d_dirac,const escript::Data& y_dirac) const
{
    throw BuckleyException("Not Implemented");
}

escript::ASM_ptr BuckleyDomain::newSystemMatrix(
                      const int row_blocksize,
                      const escript::FunctionSpace& row_functionspace,
                      const int column_blocksize,
                      const escript::FunctionSpace& column_functionspace,
                      const int type) const
{
    throw BuckleyException("Not Implemented");
}

escript::ATP_ptr BuckleyDomain::newTransportProblem(
                      const bool useBackwardEuler,
                      const int blocksize,
                      const escript::FunctionSpace& functionspace,
                      const int type) const
{
    throw BuckleyException("Not Implemented");
}

int BuckleyDomain::getNumDataPointsGlobal() const
{
    throw BuckleyException("Not Implemented");
}


  BUCKLEY_DLL_API
std::pair<int,int> BuckleyDomain::getDataShape(int functionSpaceCode) const
{
   switch (functionSpaceCode)
   {
     case ctsfn: return std::pair<int,int>(1,numpts);  
     
   // for element instead of node based reps we should return 8? but we don't have any of those yet  
   }
   throw BuckleyException("Not Implemented");  
  
}

BUCKLEY_DLL_API
void BuckleyDomain::setNewX(const escript::Data& arg)
{
    throw BuckleyException("This domain does not support changing coordinates");  
}


BUCKLEY_DLL_API
void BuckleyDomain::Print_Mesh_Info(const bool full) const
{
    if (modified)
    {
        processMods();
    }
    std::cout << "Buckley tree with " << ot.leafCount();
    std::cout << ((ot.leafCount()>1)?" leaves ":" leaf ") << "and " << numpts << " unknowns\n";
}

void BuckleyDomain::refinePoint(double x, double y, double z, unsigned int desdepth)
{
    modified=true;
    ot.splitPoint(x, y, z, desdepth);  
}



void BuckleyDomain::collapseAll(unsigned desdepth)
{
    modified=true;
    ot.allCollapse(desdepth);
}

void BuckleyDomain::collapsePoint(double x, double y, double z, unsigned int desdepth)
{
    modified=true;
    ot.collapsePoint(x, y, z, desdepth);  
}



