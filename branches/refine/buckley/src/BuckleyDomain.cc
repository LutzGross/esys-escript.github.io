
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

#include <pasowrap/PatternBuilder.h>
#include <pasowrap/SystemMatrixAdapter.h>

using namespace buckley;
using namespace std;
using namespace paso;


namespace
{
  // Need to put some thought into how many spaces we actually need
    const int initialspace=1;   
    const int invalidspace=-1;
    
    const int notIMPLEMENTED=-1;
    
    const int ctsfn=initialspace;
    const int discfn=2;
    const int red_discfn=-3;	// reduced discontinuous function  (reduced elements in finley terms)
    const int disc_faces=-4;
    const int red_disc_faces=-6;
    
    
    
    const unsigned int NUM_QUAD=8;	// number of quadrature points required


escript::DataTypes::ShapeType make3()
{
    escript::DataTypes::ShapeType s;
    s.push_back(3);
    return s;
}

escript::DataTypes::ShapeType make3x3()
{
    escript::DataTypes::ShapeType s;
    s.push_back(3);
    s.push_back(3);
    return s;
}


    const escript::DataTypes::ShapeType VECTOR3=make3();
    const escript::DataTypes::ShapeType MAT3X3=make3x3();

}


BuckleyDomain::BuckleyDomain(double x, double y, double z)
:ot(x,y,z),modified(true),leaves(0),numpts(0),samplerefids(0), m_generation(1)
{
    m_mpiInfo = Esys_MPIInfo_alloc(MPI_COMM_WORLD);  
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
    Esys_MPIInfo_free(m_mpiInfo);
}

bool BuckleyDomain::isValidFunctionSpaceType(int functionSpaceType) const
{
    return (functionSpaceType==ctsfn || functionSpaceType==discfn);
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
   else if (arg.getFunctionSpace().getTypeCode()==getFunctionCode())	// not on the nodes
   {
      // This is actually simpler because there is no overlap between different leaves
      escript::DataTypes::ValueType::ValueType vec=(&arg.getDataPointRW(0, 0));     
      int i;	// one day OMP-3 will be default
      #pragma omp parallel for private(i) schedule(static)
      for (i=0;i<ot.leafCount();++i)
      {
	  const OctCell& oc=*leaves[i];
 
	  // So this section is likely to be horribly slow
	  // Once the system is working unroll this
	  for (unsigned int j=0;j<NUM_QUAD;++j)
	  {
	      double x, y, z;
	      oc.quadCoords(j,x,y,z);
	      vec[i*NUM_QUAD*3+3*j]=x;
	      vec[i*NUM_QUAD*3+3*j+1]=y;
	      vec[i*NUM_QUAD*3+3*j+2]=z;
	  }
	
      }      
   }
   else		
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
    if (!modified)
    {
        modified=true;
	m_generation++;
    }  
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
    return discfn;  
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
   // This setup chunk copied from finley
   int reduceRowOrder=0;
   int reduceColOrder=0;
   // is the domain right?
   const BuckleyDomain& row_domain=dynamic_cast<const BuckleyDomain&>(*(row_functionspace.getDomain()));
   if (row_domain!=*this) 
      throw BuckleyException("Error - domain of row function space does not match the domain of matrix generator.");
   const BuckleyDomain& col_domain=dynamic_cast<const BuckleyDomain&>(*(column_functionspace.getDomain()));
   if (col_domain!=*this) 
      throw BuckleyException("Error - domain of columnn function space does not match the domain of matrix generator.");
   // is the function space type right 
   
//    if (row_functionspace.getTypeCode()==DegreesOfFreedom) {
      reduceRowOrder=0;
//    } else if (row_functionspace.getTypeCode()==ReducedDegreesOfFreedom) {
//       reduceRowOrder=1;
//    } else {
//       throw FinleyAdapterException("Error - illegal function space type for system matrix rows.");
//    }
//    if (column_functionspace.getTypeCode()==DegreesOfFreedom) {
      reduceColOrder=0;
/*   } else if (column_functionspace.getTypeCode()==ReducedDegreesOfFreedom) {
      reduceColOrder=1;*/
//    } else {
//       throw FinleyAdapterException("Error - illegal function space type for system matrix columns.");
//    }
   // generate matrix:
   
   // This is just here to remind me that I need to consider multiple ranks
   const unsigned int numranks=1;
   PatternBuilder* pb=makePB( m_mpiInfo, numpts/numranks,26);
   
   index_t dummy[2];
   dummy[0]=0;
   dummy[1]=numpts;
   
   pb->setDistribution(dummy, 1, 0, false);
   pb->setDistribution(dummy, 1, 0, true);
   
   Paso_SharedComponents* send=Paso_SharedComponents_alloc( 0, 0, 0, 0, 0, 0, 0, m_mpiInfo);
   Paso_SharedComponents* recv=Paso_SharedComponents_alloc(0, 0, 0, 0, 0, 0, 0, m_mpiInfo);
   
   Paso_Connector* conn=Paso_Connector_alloc(send, recv);
   
   // again, dummy values for a sole rank
   Paso_SystemMatrixPattern* psystemMatrixPattern=pb->generatePattern(0, numpts, conn ); 
   Paso_SystemMatrix* sm=Paso_SystemMatrix_alloc(MATRIX_FORMAT_DEFAULT, psystemMatrixPattern, numpts, numpts, true);
   Paso_SystemMatrixPattern_free(psystemMatrixPattern);

   SystemMatrixAdapter* sma=paso::makeSystemMatrixAdapter(sm, row_blocksize, row_functionspace, column_blocksize, column_functionspace);
   
   return escript::ASM_ptr(sma);
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
     case discfn: return std::pair<int, int>(NUM_QUAD,ot.leafCount());
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
    std::cout << " Generation=" << m_generation << std::endl;
}

void BuckleyDomain::refinePoint(double x, double y, double z, unsigned int desdepth)
{
    if (!modified)
    {
        modified=true;
	m_generation++;
    }
    ot.splitPoint(x, y, z, desdepth);  
}



void BuckleyDomain::collapseAll(unsigned desdepth)
{
    if (!modified)
    {
        modified=true;
	m_generation++;
    }    
    ot.allCollapse(desdepth);
}

void BuckleyDomain::collapsePoint(double x, double y, double z, unsigned int desdepth)
{
    if (!modified)
    {
        modified=true;
	m_generation++;
    }
    ot.collapsePoint(x, y, z, desdepth);  
}

unsigned BuckleyDomain::getGeneration() const
{
    return m_generation;
}


std::string BuckleyDomain::functionSpaceTypeAsString(int functionSpaceType) const
{
    switch (functionSpaceType)
    {
    case ctsfn: return "ContinuousFunction";
    case discfn: return "DiscontinuousFunction";
    default: return "Invalid";  
    };
}

/*
// interpolate an element from ctsfn to discfn
// qvalues and values2 must be big enough to hold a datapoint for each quadrature point
// the results end up in the values2 array
void BuckleyDomain::interpolateElementFromCtsToDisc(const LeafInfo* li, size_t ptsize, 
						    double* qvalues, double* values2) const
{
    // qvalues does double duty.
    // first it is used to hold the values of the corner points

    // we will do this with a sequence of interpolations    
    // first allong the positive x-lines [0->1], [3->2], [4->5], [6->7]
    
    for (size_t k=0;k<ptsize;++k)
    {
        values2[k]=qvalues[k]+0.25*(qvalues[ptsize+k]-qvalues[k]);	//0
        values2[ptsize+k]=qvalues[k]+0.75*(qvalues[ptsize+k]-qvalues[k]); //1
      
        values2[3*ptsize+k]=qvalues[3*ptsize+k]+0.25*(qvalues[2*ptsize+k]-qvalues[3*ptsize+k]);	//3
        values2[2*ptsize+k]=qvalues[3*ptsize+k]+0.75*(qvalues[2*ptsize+k]-qvalues[3*ptsize+k]); //2 
	
        values2[4*ptsize+k]=qvalues[4*ptsize+k]+0.25*(qvalues[5*ptsize+k]-qvalues[4*ptsize+k]);	//4
        values2[5*ptsize+k]=qvalues[4*ptsize+k]+0.75*(qvalues[5*ptsize+k]-qvalues[4*ptsize+k]); //5
      
        values2[7*ptsize+k]=qvalues[7*ptsize+k]+0.25*(qvalues[6*ptsize+k]-qvalues[7*ptsize+k]);	//7
        values2[6*ptsize+k]=qvalues[7*ptsize+k]+0.75*(qvalues[6*ptsize+k]-qvalues[7*ptsize+k]); //6     		
    }
    
    // now along the postive y-lines from values2 to qvalues
    for (size_t k=0;k<ptsize;++k)
    {
        qvalues[k]=values2[k]+0.25*(values2[3*ptsize+k]-values2[k]);	//0
        qvalues[3*ptsize+k]=values2[k]+0.75*(values2[3*ptsize+k]-values2[k]); //3
      
        qvalues[1*ptsize+k]=values2[1*ptsize+k]+0.25*(values2[2*ptsize+k]-values2[1*ptsize+k]);	//1
        qvalues[2*ptsize+k]=values2[1*ptsize+k]+0.75*(values2[2*ptsize+k]-values2[1*ptsize+k]); //2 
	
        qvalues[4*ptsize+k]=values2[4*ptsize+k]+0.25*(values2[7*ptsize+k]-values2[4*ptsize+k]);	//4
        qvalues[7*ptsize+k]=values2[4*ptsize+k]+0.75*(values2[7*ptsize+k]-values2[4*ptsize+k]); //7
      
        qvalues[5*ptsize+k]=values2[5*ptsize+k]+0.25*(values2[6*ptsize+k]-values2[5*ptsize+k]);	//5
        qvalues[6*ptsize+k]=values2[5*ptsize+k]+0.75*(values2[6*ptsize+k]-values2[5*ptsize+k]); //6     		
    }    
    
    // now the z-order
    for (size_t k=0;k<ptsize;++k)
    {
        values2[k]=qvalues[k]+0.25*(qvalues[4*ptsize+k]-qvalues[k]);	//0
        values2[4*ptsize+k]=qvalues[k]+0.75*(qvalues[4*ptsize+k]-qvalues[k]); //4
      
        values2[1*ptsize+k]=qvalues[1*ptsize+k]+0.25*(qvalues[5*ptsize+k]-qvalues[1*ptsize+k]);	//1
        values2[5*ptsize+k]=qvalues[1*ptsize+k]+0.75*(qvalues[5*ptsize+k]-qvalues[1*ptsize+k]); //5 
	
        values2[2*ptsize+k]=qvalues[2*ptsize+k]+0.25*(qvalues[6*ptsize+k]-qvalues[2*ptsize+k]);	//2
        values2[6*ptsize+k]=qvalues[2*ptsize+k]+0.75*(qvalues[6*ptsize+k]-qvalues[2*ptsize+k]); //6
      
        values2[3*ptsize+k]=qvalues[3*ptsize+k]+0.25*(qvalues[7*ptsize+k]-qvalues[3*ptsize+k]);	//3
        values2[7*ptsize+k]=qvalues[3*ptsize+k]+0.75*(qvalues[7*ptsize+k]-qvalues[3*ptsize+k]); //7     		
    }        
       
}*/



#define GETHANGSAMPLE(LEAF, KID, VNAME) const register double* VNAME; if (li->pmap[KID]<2) {\
VNAME=buffer+buffcounter;\
if (whichchild<0) {whichchild=LEAF->whichChild();\
src1=const_cast<escript::Data&>(in).getSampleDataRO(li->pmap[whichchild]-2);}\
const double* src2=const_cast<escript::Data&>(in).getSampleDataRO(LEAF->parent->kids[KID]->leafinfo->pmap[KID]-2);\
for (int k=0;k<numComp;++k)\
{\
    buffer[buffcounter++]=(src1[k]+src2[k])/2;\
}\
} else {VNAME=in.getSampleDataRO(li->pmap[KID]-2);}

// Code from Lutz' magic generator
void BuckleyDomain::setToGradient(escript::Data& out, const escript::Data& cIn) const
{
    escript::Data& in = *const_cast<escript::Data*>(&cIn);
    const dim_t numComp = in.getDataPointSize();
    if (modified)	// is the cached data we have about this domain stale?
    {
        processMods();
    }
    bool err=false;
#pragma omp parallel
    {
      double* buffer=new double[numComp*8];	// we will never have 8 hanging nodes
      if (out.getFunctionSpace().getTypeCode() == discfn) {


#pragma omp for
        for (index_t leaf=0; leaf < ot.leafCount(); ++leaf) {
	    const LeafInfo* li=leaves[leaf]->leafinfo;
	  
	    const double h0 = leaves[leaf]->sides[0];
	    const double h1 = leaves[leaf]->sides[1];
	    const double h2 = leaves[leaf]->sides[2];

	    /* GENERATOR SNIP_GRAD_ELEMENTS TOP */
	    const double tmp0_22 = -0.044658198738520451079/h1;
	    const double tmp0_16 = 0.16666666666666666667/h0;
	    const double tmp0_33 = -0.62200846792814621559/h1;
	    const double tmp0_0 = -0.62200846792814621559/h0;
	    const double tmp0_21 = -0.16666666666666666667/h1;
	    const double tmp0_17 = 0.62200846792814621559/h0;
	    const double tmp0_52 = -0.044658198738520451079/h2;
	    const double tmp0_1 = -0.16666666666666666667/h0;
	    const double tmp0_20 = -0.62200846792814621559/h1;
	    const double tmp0_14 = -0.044658198738520451079/h0;
	    const double tmp0_53 = -0.62200846792814621559/h2;
	    const double tmp0_49 = 0.16666666666666666667/h2;
	    const double tmp0_2 = 0.16666666666666666667/h0;
	    const double tmp0_27 = -0.044658198738520451079/h1;
	    const double tmp0_15 = -0.16666666666666666667/h0;
	    const double tmp0_50 = -0.16666666666666666667/h2;
	    const double tmp0_48 = 0.62200846792814621559/h2;
	    const double tmp0_3 = 0.044658198738520451079/h0;
	    const double tmp0_26 = -0.16666666666666666667/h1;
	    const double tmp0_12 = -0.62200846792814621559/h0;
	    const double tmp0_51 = 0.044658198738520451079/h2;
	    const double tmp0_25 = 0.62200846792814621559/h1;
	    const double tmp0_13 = 0.16666666666666666667/h0;
	    const double tmp0_56 = 0.16666666666666666667/h2;
	    const double tmp0_24 = 0.16666666666666666667/h1;
	    const double tmp0_10 = 0.62200846792814621559/h0;
	    const double tmp0_57 = 0.62200846792814621559/h2;
	    const double tmp0_11 = -0.16666666666666666667/h0;
	    const double tmp0_54 = -0.044658198738520451079/h2;
	    const double tmp0_38 = 0.16666666666666666667/h1;
	    const double tmp0_34 = -0.044658198738520451079/h1;
	    const double tmp0_42 = 0.16666666666666666667/h2;
	    const double tmp0_35 = -0.16666666666666666667/h1;
	    const double tmp0_36 = -0.62200846792814621559/h1;
	    const double tmp0_41 = 0.62200846792814621559/h2;
	    const double tmp0_8 = 0.044658198738520451079/h0;
	    const double tmp0_37 = 0.62200846792814621559/h1;
	    const double tmp0_29 = 0.16666666666666666667/h1;
	    const double tmp0_40 = -0.62200846792814621559/h2;
	    const double tmp0_9 = 0.16666666666666666667/h0;
	    const double tmp0_30 = 0.62200846792814621559/h1;
	    const double tmp0_28 = -0.16666666666666666667/h1;
	    const double tmp0_43 = 0.044658198738520451079/h2;
	    const double tmp0_32 = 0.16666666666666666667/h1;
	    const double tmp0_31 = 0.044658198738520451079/h1;
	    const double tmp0_39 = 0.044658198738520451079/h1;
	    const double tmp0_58 = -0.62200846792814621559/h2;
	    const double tmp0_55 = 0.044658198738520451079/h2;
	    const double tmp0_18 = -0.62200846792814621559/h0;
	    const double tmp0_45 = -0.16666666666666666667/h2;
	    const double tmp0_59 = -0.16666666666666666667/h2;
	    const double tmp0_4 = -0.044658198738520451079/h0;
	    const double tmp0_19 = 0.044658198738520451079/h0;
	    const double tmp0_44 = -0.044658198738520451079/h2;
	    const double tmp0_5 = 0.62200846792814621559/h0;
	    const double tmp0_47 = 0.16666666666666666667/h2;
	    const double tmp0_6 = -0.16666666666666666667/h0;
	    const double tmp0_23 = 0.044658198738520451079/h1;
	    const double tmp0_46 = -0.16666666666666666667/h2;
	    const double tmp0_7 = -0.044658198738520451079/h0;
	  
	    int buffcounter=0;
	    const double* src1=0;
	    int whichchild=-1;

	    GETHANGSAMPLE(leaves[leaf], 0, f_000)
	    GETHANGSAMPLE(leaves[leaf], 4, f_001)
	    GETHANGSAMPLE(leaves[leaf], 5, f_101)
	    GETHANGSAMPLE(leaves[leaf], 6, f_111)
	    GETHANGSAMPLE(leaves[leaf], 2, f_110)
	    GETHANGSAMPLE(leaves[leaf], 7, f_011)
	    GETHANGSAMPLE(leaves[leaf], 3, f_010)
	    GETHANGSAMPLE(leaves[leaf], 1, f_100)
//	    const register double* f_000 = (li->pmap[0]>1)?in.getSampleDataRO(li->pmap[0]-2):HANG_INTERPOLATE(leaves[leaf], 0);	    
// 	    const register double* f_001 = (li->pmap[1]>1)?in.getSampleDataRO(li->pmap[1]-2):HANG_INTERPOLATE(leaves[leaf], 1);
// 	    const register double* f_101 = (li->pmap[5]>1)?in.getSampleDataRO(li->pmap[5]-2):HANG_INTERPOLATE(leaves[leaf], 5);
// 	    const register double* f_111 = (li->pmap[6]>1)?in.getSampleDataRO(li->pmap[6]-2):HANG_INTERPOLATE(leaves[leaf], 6);
// 	    const register double* f_110 = (li->pmap[7]>1)?in.getSampleDataRO(li->pmap[7]-2):HANG_INTERPOLATE(leaves[leaf], 7);
// 	    const register double* f_011 = (li->pmap[2]>1)?in.getSampleDataRO(li->pmap[2]-2):HANG_INTERPOLATE(leaves[leaf], 2);
// 	    const register double* f_010 = (li->pmap[3]>1)?in.getSampleDataRO(li->pmap[3]-2):HANG_INTERPOLATE(leaves[leaf], 3);
// 	    const register double* f_100 = (li->pmap[4]>1)?in.getSampleDataRO(li->pmap[4]-2):HANG_INTERPOLATE(leaves[leaf], 4);

	    
	    double* o = out.getSampleDataRW(leaf);

 
	    for (index_t i=0; i < numComp; ++i) {
		o[INDEX3(i,0,0,numComp,3)] = f_000[i]*tmp0_0 + f_011[i]*tmp0_4 + f_100[i]*tmp0_5 + f_111[i]*tmp0_3 + tmp0_1*(f_001[i] + f_010[i]) + tmp0_2*(f_101[i] + f_110[i]);
		o[INDEX3(i,1,0,numComp,3)] = f_000[i]*tmp0_20 + f_010[i]*tmp0_25 + f_101[i]*tmp0_22 + f_111[i]*tmp0_23 + tmp0_21*(f_001[i] + f_100[i]) + tmp0_24*(f_011[i] + f_110[i]);
		o[INDEX3(i,2,0,numComp,3)] = f_000[i]*tmp0_40 + f_001[i]*tmp0_41 + f_110[i]*tmp0_44 + f_111[i]*tmp0_43 + tmp0_42*(f_011[i] + f_101[i]) + tmp0_45*(f_010[i] + f_100[i]);
		o[INDEX3(i,0,1,numComp,3)] = f_000[i]*tmp0_0 + f_011[i]*tmp0_4 + f_100[i]*tmp0_5 + f_111[i]*tmp0_3 + tmp0_1*(f_001[i] + f_010[i]) + tmp0_2*(f_101[i] + f_110[i]);
		o[INDEX3(i,1,1,numComp,3)] = f_000[i]*tmp0_26 + f_001[i]*tmp0_27 + f_010[i]*tmp0_32 + f_011[i]*tmp0_31 + f_100[i]*tmp0_33 + f_101[i]*tmp0_28 + f_110[i]*tmp0_30 + f_111[i]*tmp0_29;
		o[INDEX3(i,2,1,numComp,3)] = f_000[i]*tmp0_46 + f_001[i]*tmp0_47 + f_010[i]*tmp0_52 + f_011[i]*tmp0_51 + f_100[i]*tmp0_53 + f_101[i]*tmp0_48 + f_110[i]*tmp0_50 + f_111[i]*tmp0_49;
		o[INDEX3(i,0,2,numComp,3)] = f_000[i]*tmp0_6 + f_001[i]*tmp0_7 + f_010[i]*tmp0_12 + f_011[i]*tmp0_11 + f_100[i]*tmp0_13 + f_101[i]*tmp0_8 + f_110[i]*tmp0_10 + f_111[i]*tmp0_9;
		o[INDEX3(i,1,2,numComp,3)] = f_000[i]*tmp0_20 + f_010[i]*tmp0_25 + f_101[i]*tmp0_22 + f_111[i]*tmp0_23 + tmp0_21*(f_001[i] + f_100[i]) + tmp0_24*(f_011[i] + f_110[i]);
		o[INDEX3(i,2,2,numComp,3)] = f_000[i]*tmp0_46 + f_001[i]*tmp0_47 + f_010[i]*tmp0_53 + f_011[i]*tmp0_48 + f_100[i]*tmp0_52 + f_101[i]*tmp0_51 + f_110[i]*tmp0_50 + f_111[i]*tmp0_49;
		o[INDEX3(i,0,3,numComp,3)] = f_000[i]*tmp0_6 + f_001[i]*tmp0_7 + f_010[i]*tmp0_12 + f_011[i]*tmp0_11 + f_100[i]*tmp0_13 + f_101[i]*tmp0_8 + f_110[i]*tmp0_10 + f_111[i]*tmp0_9;
		o[INDEX3(i,1,3,numComp,3)] = f_000[i]*tmp0_26 + f_001[i]*tmp0_27 + f_010[i]*tmp0_32 + f_011[i]*tmp0_31 + f_100[i]*tmp0_33 + f_101[i]*tmp0_28 + f_110[i]*tmp0_30 + f_111[i]*tmp0_29;
		o[INDEX3(i,2,3,numComp,3)] = f_000[i]*tmp0_54 + f_001[i]*tmp0_55 + f_110[i]*tmp0_58 + f_111[i]*tmp0_57 + tmp0_56*(f_011[i] + f_101[i]) + tmp0_59*(f_010[i] + f_100[i]);
		o[INDEX3(i,0,4,numComp,3)] = f_000[i]*tmp0_6 + f_001[i]*tmp0_12 + f_010[i]*tmp0_7 + f_011[i]*tmp0_11 + f_100[i]*tmp0_13 + f_101[i]*tmp0_10 + f_110[i]*tmp0_8 + f_111[i]*tmp0_9;
		o[INDEX3(i,1,4,numComp,3)] = f_000[i]*tmp0_26 + f_001[i]*tmp0_33 + f_010[i]*tmp0_32 + f_011[i]*tmp0_30 + f_100[i]*tmp0_27 + f_101[i]*tmp0_28 + f_110[i]*tmp0_31 + f_111[i]*tmp0_29;
		o[INDEX3(i,2,4,numComp,3)] = f_000[i]*tmp0_40 + f_001[i]*tmp0_41 + f_110[i]*tmp0_44 + f_111[i]*tmp0_43 + tmp0_42*(f_011[i] + f_101[i]) + tmp0_45*(f_010[i] + f_100[i]);
		o[INDEX3(i,0,5,numComp,3)] = f_000[i]*tmp0_6 + f_001[i]*tmp0_12 + f_010[i]*tmp0_7 + f_011[i]*tmp0_11 + f_100[i]*tmp0_13 + f_101[i]*tmp0_10 + f_110[i]*tmp0_8 + f_111[i]*tmp0_9;
		o[INDEX3(i,1,5,numComp,3)] = f_000[i]*tmp0_34 + f_010[i]*tmp0_39 + f_101[i]*tmp0_36 + f_111[i]*tmp0_37 + tmp0_35*(f_001[i] + f_100[i]) + tmp0_38*(f_011[i] + f_110[i]);
		o[INDEX3(i,2,5,numComp,3)] = f_000[i]*tmp0_46 + f_001[i]*tmp0_47 + f_010[i]*tmp0_52 + f_011[i]*tmp0_51 + f_100[i]*tmp0_53 + f_101[i]*tmp0_48 + f_110[i]*tmp0_50 + f_111[i]*tmp0_49;
		o[INDEX3(i,0,6,numComp,3)] = f_000[i]*tmp0_14 + f_011[i]*tmp0_18 + f_100[i]*tmp0_19 + f_111[i]*tmp0_17 + tmp0_15*(f_001[i] + f_010[i]) + tmp0_16*(f_101[i] + f_110[i]);
		o[INDEX3(i,1,6,numComp,3)] = f_000[i]*tmp0_26 + f_001[i]*tmp0_33 + f_010[i]*tmp0_32 + f_011[i]*tmp0_30 + f_100[i]*tmp0_27 + f_101[i]*tmp0_28 + f_110[i]*tmp0_31 + f_111[i]*tmp0_29;
		o[INDEX3(i,2,6,numComp,3)] = f_000[i]*tmp0_46 + f_001[i]*tmp0_47 + f_010[i]*tmp0_53 + f_011[i]*tmp0_48 + f_100[i]*tmp0_52 + f_101[i]*tmp0_51 + f_110[i]*tmp0_50 + f_111[i]*tmp0_49;
		o[INDEX3(i,0,7,numComp,3)] = f_000[i]*tmp0_14 + f_011[i]*tmp0_18 + f_100[i]*tmp0_19 + f_111[i]*tmp0_17 + tmp0_15*(f_001[i] + f_010[i]) + tmp0_16*(f_101[i] + f_110[i]);
		o[INDEX3(i,1,7,numComp,3)] = f_000[i]*tmp0_34 + f_010[i]*tmp0_39 + f_101[i]*tmp0_36 + f_111[i]*tmp0_37 + tmp0_35*(f_001[i] + f_100[i]) + tmp0_38*(f_011[i] + f_110[i]);
		o[INDEX3(i,2,7,numComp,3)] = f_000[i]*tmp0_54 + f_001[i]*tmp0_55 + f_110[i]*tmp0_58 + f_111[i]*tmp0_57 + tmp0_56*(f_011[i] + f_101[i]) + tmp0_59*(f_010[i] + f_100[i]);
	    } /* end of component loop i */
        } /* end of leaf loop */
        /* GENERATOR SNIP_GRAD_ELEMENTS BOTTOM */
    } else if (out.getFunctionSpace().getTypeCode() == red_discfn) {
        /* GENERATOR SNIP_GRAD_REDUCED_ELEMENTS TOP */
#pragma omp for
        for (index_t leaf=0; leaf < ot.leafCount(); ++leaf) {
	    const LeafInfo* li=leaves[leaf]->leafinfo;
	  
	    const double h0 = leaves[leaf]->sides[0];
	    const double h1 = leaves[leaf]->sides[1];
	    const double h2 = leaves[leaf]->sides[2];
	    
	    const double tmp0_0 = -0.25/h0;
	    const double tmp0_4 = -0.25/h2;
	    const double tmp0_1 = 0.25/h0;
	    const double tmp0_5 = 0.25/h2;
	    const double tmp0_2 = -0.25/h1;
	    const double tmp0_3 = 0.25/h1;

	    int buffcounter=0;
	    const double* src1=0;
	    int whichchild=-1;	    
	    GETHANGSAMPLE(leaves[leaf], 0, f_000)
	    GETHANGSAMPLE(leaves[leaf], 4, f_001)
	    GETHANGSAMPLE(leaves[leaf], 5, f_101)
	    GETHANGSAMPLE(leaves[leaf], 6, f_111)
	    GETHANGSAMPLE(leaves[leaf], 2, f_110)
	    GETHANGSAMPLE(leaves[leaf], 7, f_011)
	    GETHANGSAMPLE(leaves[leaf], 3, f_010)
	    GETHANGSAMPLE(leaves[leaf], 1, f_100)

	    
// 	    const register double* f_000 = (li->pmap[0]>1)?in.getSampleDataRO(li->pmap[0]-2):HANG_INTERPOLATE(leaves[leaf], 0);
// 	    const register double* f_001 = (li->pmap[1]>1)?in.getSampleDataRO(li->pmap[1]-2):HANG_INTERPOLATE(leaves[leaf], 1);
// 	    const register double* f_101 = (li->pmap[5]>1)?in.getSampleDataRO(li->pmap[5]-2):HANG_INTERPOLATE(leaves[leaf], 5);
// 	    const register double* f_111 = (li->pmap[6]>1)?in.getSampleDataRO(li->pmap[6]-2):HANG_INTERPOLATE(leaves[leaf], 6);
// 	    const register double* f_110 = (li->pmap[7]>1)?in.getSampleDataRO(li->pmap[7]-2):HANG_INTERPOLATE(leaves[leaf], 7);
// 	    const register double* f_011 = (li->pmap[2]>1)?in.getSampleDataRO(li->pmap[2]-2):HANG_INTERPOLATE(leaves[leaf], 2);
// 	    const register double* f_010 = (li->pmap[3]>1)?in.getSampleDataRO(li->pmap[3]-2):HANG_INTERPOLATE(leaves[leaf], 3);
// 	    const register double* f_100 = (li->pmap[4]>1)?in.getSampleDataRO(li->pmap[4]-2):HANG_INTERPOLATE(leaves[leaf], 4);	    

            double* o = out.getSampleDataRW(leaf);
	    for (index_t i=0; i < numComp; ++i) {
		o[INDEX3(i,0,0,numComp,3)] = tmp0_0*(f_000[i] + f_001[i] + f_010[i] + f_011[i]) + tmp0_1*(f_100[i] + f_101[i] + f_110[i] + f_111[i]);
		o[INDEX3(i,1,0,numComp,3)] = tmp0_2*(f_000[i] + f_001[i] + f_100[i] + f_101[i]) + tmp0_3*(f_010[i] + f_011[i] + f_110[i] + f_111[i]);
		o[INDEX3(i,2,0,numComp,3)] = tmp0_4*(f_000[i] + f_010[i] + f_100[i] + f_110[i]) + tmp0_5*(f_001[i] + f_011[i] + f_101[i] + f_111[i]);
	    } /* end of component loop i */
        } /* end of leaf loop */
        /* GENERATOR SNIP_GRAD_REDUCED_ELEMENTS BOTTOM */
    } else if (out.getFunctionSpace().getTypeCode() == disc_faces) {
        /* GENERATOR SNIP_GRAD_FACES TOP */
        if (face_cells[0].size() > 0) {		// left face

#pragma omp for
          for (index_t leaf=0; leaf < face_cells[0].size(); ++leaf) {
	    const LeafInfo* li=face_cells[0][leaf]->leafinfo;		// we aren't iterating over all leaves this time
	  
	    const double h0 = face_cells[0][leaf]->sides[0];
	    const double h1 = face_cells[0][leaf]->sides[1];
	    const double h2 = face_cells[0][leaf]->sides[2];

            const double tmp0_22 = 0.21132486540518711775/h1;
            const double tmp0_16 = 0.16666666666666666667/h0;
            const double tmp0_33 = 0.21132486540518711775/h2;
            const double tmp0_0 = -0.62200846792814621559/h0;
            const double tmp0_21 = -0.21132486540518711775/h1;
            const double tmp0_17 = 0.62200846792814621559/h0;
            const double tmp0_1 = -0.16666666666666666667/h0;
            const double tmp0_20 = -0.78867513459481288225/h1;
            const double tmp0_14 = -0.044658198738520451079/h0;
            const double tmp0_2 = 0.16666666666666666667/h0;
            const double tmp0_27 = 0.21132486540518711775/h1;
            const double tmp0_15 = -0.16666666666666666667/h0;
            const double tmp0_3 = 0.044658198738520451079/h0;
            const double tmp0_26 = 0.78867513459481288225/h1;
            const double tmp0_12 = -0.62200846792814621559/h0;
            const double tmp0_25 = -0.78867513459481288225/h1;
            const double tmp0_13 = 0.16666666666666666667/h0;
            const double tmp0_24 = -0.21132486540518711775/h1;
            const double tmp0_10 = 0.62200846792814621559/h0;
            const double tmp0_11 = -0.16666666666666666667/h0;
            const double tmp0_34 = 0.78867513459481288225/h2;
            const double tmp0_35 = -0.78867513459481288225/h2;
            const double tmp0_8 = 0.044658198738520451079/h0;
            const double tmp0_29 = 0.78867513459481288225/h2;
            const double tmp0_9 = 0.16666666666666666667/h0;
            const double tmp0_30 = 0.21132486540518711775/h2;
            const double tmp0_28 = -0.78867513459481288225/h2;
            const double tmp0_32 = -0.21132486540518711775/h2;
            const double tmp0_31 = -0.21132486540518711775/h2;
            const double tmp0_18 = -0.62200846792814621559/h0;
            const double tmp0_4 = -0.044658198738520451079/h0;
            const double tmp0_19 = 0.044658198738520451079/h0;
            const double tmp0_5 = 0.62200846792814621559/h0;
            const double tmp0_6 = -0.16666666666666666667/h0;
            const double tmp0_23 = 0.78867513459481288225/h1;
            const double tmp0_7 = -0.044658198738520451079/h0;	  
	    
	    int buffcounter=0;
	    const double* src1=0;
	    int whichchild=-1;	    
	    
	    GETHANGSAMPLE(face_cells[0][leaf], 0, f_000)
	    GETHANGSAMPLE(face_cells[0][leaf], 4, f_001)
	    GETHANGSAMPLE(face_cells[0][leaf], 5, f_101)
	    GETHANGSAMPLE(face_cells[0][leaf], 6, f_111)
	    GETHANGSAMPLE(face_cells[0][leaf], 2, f_110)
	    GETHANGSAMPLE(face_cells[0][leaf], 7, f_011)
	    GETHANGSAMPLE(face_cells[0][leaf], 3, f_010)
	    GETHANGSAMPLE(face_cells[0][leaf], 1, f_100)

/*	    const register double* f_000 = (li->pmap[0]>1)?in.getSampleDataRO(li->pmap[0]-2):HANG_INTERPOLATE(face_cells[0][leaf], 0);
	    const register double* f_001 = (li->pmap[1]>1)?in.getSampleDataRO(li->pmap[1]-2):HANG_INTERPOLATE(face_cells[0][leaf], 1);
	    const register double* f_101 = (li->pmap[5]>1)?in.getSampleDataRO(li->pmap[5]-2):HANG_INTERPOLATE(face_cells[0][leaf], 5);
	    const register double* f_111 = (li->pmap[6]>1)?in.getSampleDataRO(li->pmap[6]-2):HANG_INTERPOLATE(face_cells[0][leaf], 6);
	    const register double* f_110 = (li->pmap[7]>1)?in.getSampleDataRO(li->pmap[7]-2):HANG_INTERPOLATE(face_cells[0][leaf], 7);
	    const register double* f_011 = (li->pmap[2]>1)?in.getSampleDataRO(li->pmap[2]-2):HANG_INTERPOLATE(face_cells[0][leaf], 2);
	    const register double* f_010 = (li->pmap[3]>1)?in.getSampleDataRO(li->pmap[3]-2):HANG_INTERPOLATE(face_cells[0][leaf], 3);
	    const register double* f_100 = (li->pmap[4]>1)?in.getSampleDataRO(li->pmap[4]-2):HANG_INTERPOLATE(face_cells[0][leaf], 4);	  	*/    
	    

            double* o = out.getSampleDataRW(leaf);	// face 0 start at the beginning of the object
	    for (index_t i=0; i < numComp; ++i) {
		o[INDEX3(i,0,0,numComp,3)] = f_000[i]*tmp0_0 + f_011[i]*tmp0_4 + f_100[i]*tmp0_5 + f_111[i]*tmp0_3 + tmp0_1*(f_001[i] + f_010[i]) + tmp0_2*(f_101[i] + f_110[i]);
		o[INDEX3(i,1,0,numComp,3)] = f_000[i]*tmp0_20 + f_001[i]*tmp0_21 + f_010[i]*tmp0_23 + f_011[i]*tmp0_22;
		o[INDEX3(i,2,0,numComp,3)] = f_000[i]*tmp0_28 + f_001[i]*tmp0_29 + f_010[i]*tmp0_31 + f_011[i]*tmp0_30;
		o[INDEX3(i,0,1,numComp,3)] = f_000[i]*tmp0_6 + f_001[i]*tmp0_7 + f_010[i]*tmp0_12 + f_011[i]*tmp0_11 + f_100[i]*tmp0_13 + f_101[i]*tmp0_8 + f_110[i]*tmp0_10 + f_111[i]*tmp0_9;
		o[INDEX3(i,1,1,numComp,3)] = f_000[i]*tmp0_20 + f_001[i]*tmp0_21 + f_010[i]*tmp0_23 + f_011[i]*tmp0_22;
		o[INDEX3(i,2,1,numComp,3)] = f_000[i]*tmp0_32 + f_001[i]*tmp0_33 + f_010[i]*tmp0_35 + f_011[i]*tmp0_34;
		o[INDEX3(i,0,2,numComp,3)] = f_000[i]*tmp0_6 + f_001[i]*tmp0_12 + f_010[i]*tmp0_7 + f_011[i]*tmp0_11 + f_100[i]*tmp0_13 + f_101[i]*tmp0_10 + f_110[i]*tmp0_8 + f_111[i]*tmp0_9;
		o[INDEX3(i,1,2,numComp,3)] = f_000[i]*tmp0_24 + f_001[i]*tmp0_25 + f_010[i]*tmp0_27 + f_011[i]*tmp0_26;
		o[INDEX3(i,2,2,numComp,3)] = f_000[i]*tmp0_28 + f_001[i]*tmp0_29 + f_010[i]*tmp0_31 + f_011[i]*tmp0_30;
		o[INDEX3(i,0,3,numComp,3)] = f_000[i]*tmp0_14 + f_011[i]*tmp0_18 + f_100[i]*tmp0_19 + f_111[i]*tmp0_17 + tmp0_15*(f_001[i] + f_010[i]) + tmp0_16*(f_101[i] + f_110[i]);
		o[INDEX3(i,1,3,numComp,3)] = f_000[i]*tmp0_24 + f_001[i]*tmp0_25 + f_010[i]*tmp0_27 + f_011[i]*tmp0_26;
		o[INDEX3(i,2,3,numComp,3)] = f_000[i]*tmp0_32 + f_001[i]*tmp0_33 + f_010[i]*tmp0_35 + f_011[i]*tmp0_34;
	    } /* end of component loop i */
	  }    
        } /* end of face 0 */
        int baseoffset=face_cells[0].size();
        if (face_cells[1].size() > 0) {
#pragma omp for
          for (index_t leaf=0; leaf < face_cells[1].size(); ++leaf) {	  
	    const LeafInfo* li=face_cells[1][leaf]->leafinfo;		// we aren't iterating over all leaves this time
	  
	    const double h0 = face_cells[1][leaf]->sides[0];
	    const double h1 = face_cells[1][leaf]->sides[1];
	    const double h2 = face_cells[1][leaf]->sides[2];
	    
            const double tmp0_22 = 0.78867513459481288225/h1;
            const double tmp0_16 = 0.16666666666666666667/h0;
            const double tmp0_33 = 0.78867513459481288225/h2;
            const double tmp0_0 = -0.62200846792814621559/h0;
            const double tmp0_21 = 0.21132486540518711775/h1;
            const double tmp0_17 = 0.62200846792814621559/h0;
            const double tmp0_1 = -0.16666666666666666667/h0;
            const double tmp0_20 = -0.21132486540518711775/h1;
            const double tmp0_14 = -0.044658198738520451079/h0;
            const double tmp0_2 = 0.16666666666666666667/h0;
            const double tmp0_27 = -0.21132486540518711775/h1;
            const double tmp0_15 = -0.16666666666666666667/h0;
            const double tmp0_3 = 0.044658198738520451079/h0;
            const double tmp0_26 = 0.21132486540518711775/h1;
            const double tmp0_12 = -0.62200846792814621559/h0;
            const double tmp0_25 = 0.78867513459481288225/h1;
            const double tmp0_13 = 0.16666666666666666667/h0;
            const double tmp0_24 = -0.78867513459481288225/h1;
            const double tmp0_10 = 0.62200846792814621559/h0;
            const double tmp0_11 = -0.16666666666666666667/h0;
            const double tmp0_34 = -0.78867513459481288225/h2;
            const double tmp0_35 = -0.21132486540518711775/h2;
            const double tmp0_8 = 0.044658198738520451079/h0;
            const double tmp0_29 = 0.21132486540518711775/h2;
            const double tmp0_9 = 0.16666666666666666667/h0;
            const double tmp0_30 = -0.21132486540518711775/h2;
            const double tmp0_28 = 0.78867513459481288225/h2;
            const double tmp0_32 = 0.21132486540518711775/h2;
            const double tmp0_31 = -0.78867513459481288225/h2;
            const double tmp0_18 = -0.62200846792814621559/h0;
            const double tmp0_4 = -0.044658198738520451079/h0;
            const double tmp0_19 = 0.044658198738520451079/h0;
            const double tmp0_5 = 0.62200846792814621559/h0;
            const double tmp0_6 = -0.16666666666666666667/h0;
            const double tmp0_23 = -0.78867513459481288225/h1;
            const double tmp0_7 = -0.044658198738520451079/h0;

	    int buffcounter=0;
	    const double* src1=0;
	    int whichchild=-1;	    

	    GETHANGSAMPLE(face_cells[1][leaf], 0, f_000)
	    GETHANGSAMPLE(face_cells[1][leaf], 4, f_001)
	    GETHANGSAMPLE(face_cells[1][leaf], 5, f_101)
	    GETHANGSAMPLE(face_cells[1][leaf], 6, f_111)
	    GETHANGSAMPLE(face_cells[1][leaf], 2, f_110)
	    GETHANGSAMPLE(face_cells[1][leaf], 7, f_011)
	    GETHANGSAMPLE(face_cells[1][leaf], 3, f_010)
	    GETHANGSAMPLE(face_cells[1][leaf], 1, f_100)	    
	    
	    
/*	    const register double* f_000 = (li->pmap[0]>1)?in.getSampleDataRO(li->pmap[0]-2):HANG_INTERPOLATE(face_cells[1][leaf], 0);
	    const register double* f_001 = (li->pmap[1]>1)?in.getSampleDataRO(li->pmap[1]-2):HANG_INTERPOLATE(face_cells[1][leaf], 1);
	    const register double* f_101 = (li->pmap[5]>1)?in.getSampleDataRO(li->pmap[5]-2):HANG_INTERPOLATE(face_cells[1][leaf], 5);
	    const register double* f_111 = (li->pmap[6]>1)?in.getSampleDataRO(li->pmap[6]-2):HANG_INTERPOLATE(face_cells[1][leaf], 6);
	    const register double* f_110 = (li->pmap[7]>1)?in.getSampleDataRO(li->pmap[7]-2):HANG_INTERPOLATE(face_cells[1][leaf], 7);
	    const register double* f_011 = (li->pmap[2]>1)?in.getSampleDataRO(li->pmap[2]-2):HANG_INTERPOLATE(face_cells[1][leaf], 2);
	    const register double* f_010 = (li->pmap[3]>1)?in.getSampleDataRO(li->pmap[3]-2):HANG_INTERPOLATE(face_cells[1][leaf], 3);
	    const register double* f_100 = (li->pmap[4]>1)?in.getSampleDataRO(li->pmap[4]-2):HANG_INTERPOLATE(face_cells[1][leaf], 4);	*/    

            double* o = out.getSampleDataRW(baseoffset+leaf);
	    for (index_t i=0; i < numComp; ++i) {
		o[INDEX3(i,0,0,numComp,3)] = f_000[i]*tmp0_0 + f_011[i]*tmp0_4 + f_100[i]*tmp0_5 + f_111[i]*tmp0_3 + tmp0_1*(f_001[i] + f_010[i]) + tmp0_2*(f_101[i] + f_110[i]);
		o[INDEX3(i,1,0,numComp,3)] = f_100[i]*tmp0_23 + f_101[i]*tmp0_20 + f_110[i]*tmp0_22 + f_111[i]*tmp0_21;
		o[INDEX3(i,2,0,numComp,3)] = f_100[i]*tmp0_31 + f_101[i]*tmp0_28 + f_110[i]*tmp0_30 + f_111[i]*tmp0_29;
		o[INDEX3(i,0,1,numComp,3)] = f_000[i]*tmp0_6 + f_001[i]*tmp0_7 + f_010[i]*tmp0_12 + f_011[i]*tmp0_11 + f_100[i]*tmp0_13 + f_101[i]*tmp0_8 + f_110[i]*tmp0_10 + f_111[i]*tmp0_9;
		o[INDEX3(i,1,1,numComp,3)] = f_100[i]*tmp0_23 + f_101[i]*tmp0_20 + f_110[i]*tmp0_22 + f_111[i]*tmp0_21;
		o[INDEX3(i,2,1,numComp,3)] = f_100[i]*tmp0_35 + f_101[i]*tmp0_32 + f_110[i]*tmp0_34 + f_111[i]*tmp0_33;
		o[INDEX3(i,0,2,numComp,3)] = f_000[i]*tmp0_6 + f_001[i]*tmp0_12 + f_010[i]*tmp0_7 + f_011[i]*tmp0_11 + f_100[i]*tmp0_13 + f_101[i]*tmp0_10 + f_110[i]*tmp0_8 + f_111[i]*tmp0_9;
		o[INDEX3(i,1,2,numComp,3)] = f_100[i]*tmp0_27 + f_101[i]*tmp0_24 + f_110[i]*tmp0_26 + f_111[i]*tmp0_25;
		o[INDEX3(i,2,2,numComp,3)] = f_100[i]*tmp0_31 + f_101[i]*tmp0_28 + f_110[i]*tmp0_30 + f_111[i]*tmp0_29;
		o[INDEX3(i,0,3,numComp,3)] = f_000[i]*tmp0_14 + f_011[i]*tmp0_18 + f_100[i]*tmp0_19 + f_111[i]*tmp0_17 + tmp0_15*(f_001[i] + f_010[i]) + tmp0_16*(f_101[i] + f_110[i]);
		o[INDEX3(i,1,3,numComp,3)] = f_100[i]*tmp0_27 + f_101[i]*tmp0_24 + f_110[i]*tmp0_26 + f_111[i]*tmp0_25;
		o[INDEX3(i,2,3,numComp,3)] = f_100[i]*tmp0_35 + f_101[i]*tmp0_32 + f_110[i]*tmp0_34 + f_111[i]*tmp0_33;
	    } /* end of component loop i */
	  }    
        } /* end of face 1 */
        baseoffset+=face_cells[1].size();
        if (face_cells[2].size() > 0) {
#pragma omp for
          for (index_t leaf=0; leaf < face_cells[2].size(); ++leaf) {	  
	    const LeafInfo* li=face_cells[2][leaf]->leafinfo;		// we aren't iterating over all leaves this time
	    
	    const double h0 = face_cells[2][leaf]->sides[0];
	    const double h1 = face_cells[2][leaf]->sides[1];
	    const double h2 = face_cells[2][leaf]->sides[2];	    
	    
            const double tmp0_22 = -0.044658198738520451079/h1;
            const double tmp0_16 = -0.16666666666666666667/h1;
            const double tmp0_33 = 0.21132486540518711775/h2;
            const double tmp0_0 = -0.78867513459481288225/h0;
            const double tmp0_21 = 0.16666666666666666667/h1;
            const double tmp0_17 = -0.62200846792814621559/h1;
            const double tmp0_1 = -0.21132486540518711775/h0;
            const double tmp0_20 = 0.044658198738520451079/h1;
            const double tmp0_14 = -0.16666666666666666667/h1;
            const double tmp0_2 = 0.21132486540518711775/h0;
            const double tmp0_27 = 0.044658198738520451079/h1;
            const double tmp0_15 = -0.044658198738520451079/h1;
            const double tmp0_3 = 0.78867513459481288225/h0;
            const double tmp0_26 = 0.16666666666666666667/h1;
            const double tmp0_12 = 0.16666666666666666667/h1;
            const double tmp0_25 = 0.62200846792814621559/h1;
            const double tmp0_13 = 0.62200846792814621559/h1;
            const double tmp0_24 = -0.62200846792814621559/h1;
            const double tmp0_10 = -0.044658198738520451079/h1;
            const double tmp0_11 = 0.044658198738520451079/h1;
            const double tmp0_34 = 0.78867513459481288225/h2;
            const double tmp0_35 = -0.78867513459481288225/h2;
            const double tmp0_8 = -0.62200846792814621559/h1;
            const double tmp0_29 = 0.78867513459481288225/h2;
            const double tmp0_9 = -0.16666666666666666667/h1;
            const double tmp0_30 = 0.21132486540518711775/h2;
            const double tmp0_28 = -0.78867513459481288225/h2;
            const double tmp0_32 = -0.21132486540518711775/h2;
            const double tmp0_31 = -0.21132486540518711775/h2;
            const double tmp0_18 = 0.16666666666666666667/h1;
            const double tmp0_4 = -0.21132486540518711775/h0;
            const double tmp0_19 = 0.62200846792814621559/h1;
            const double tmp0_5 = -0.78867513459481288225/h0;
            const double tmp0_6 = 0.78867513459481288225/h0;
            const double tmp0_23 = -0.16666666666666666667/h1;
            const double tmp0_7 = 0.21132486540518711775/h0;
	    
	    int buffcounter=0;
	    const double* src1=0;
	    int whichchild=-1;	    

	    GETHANGSAMPLE(face_cells[2][leaf], 0, f_000)
	    GETHANGSAMPLE(face_cells[2][leaf], 4, f_001)
	    GETHANGSAMPLE(face_cells[2][leaf], 5, f_101)
	    GETHANGSAMPLE(face_cells[2][leaf], 6, f_111)
	    GETHANGSAMPLE(face_cells[2][leaf], 2, f_110)
	    GETHANGSAMPLE(face_cells[2][leaf], 7, f_011)
	    GETHANGSAMPLE(face_cells[2][leaf], 3, f_010)
	    GETHANGSAMPLE(face_cells[2][leaf], 1, f_100)	    
	    
	    
/*	    const register double* f_000 = (li->pmap[0]>1)?in.getSampleDataRO(li->pmap[0]-2):HANG_INTERPOLATE(face_cells[2][leaf], 0);
	    const register double* f_001 = (li->pmap[1]>1)?in.getSampleDataRO(li->pmap[1]-2):HANG_INTERPOLATE(face_cells[2][leaf], 1);
	    const register double* f_101 = (li->pmap[5]>1)?in.getSampleDataRO(li->pmap[5]-2):HANG_INTERPOLATE(face_cells[2][leaf], 5);
	    const register double* f_111 = (li->pmap[6]>1)?in.getSampleDataRO(li->pmap[6]-2):HANG_INTERPOLATE(face_cells[2][leaf], 6);
	    const register double* f_110 = (li->pmap[7]>1)?in.getSampleDataRO(li->pmap[7]-2):HANG_INTERPOLATE(face_cells[2][leaf], 7);
	    const register double* f_011 = (li->pmap[2]>1)?in.getSampleDataRO(li->pmap[2]-2):HANG_INTERPOLATE(face_cells[2][leaf], 2);
	    const register double* f_010 = (li->pmap[3]>1)?in.getSampleDataRO(li->pmap[3]-2):HANG_INTERPOLATE(face_cells[2][leaf], 3);
	    const register double* f_100 = (li->pmap[4]>1)?in.getSampleDataRO(li->pmap[4]-2):HANG_INTERPOLATE(face_cells[2][leaf], 4);	*/   	    
	    
	    double* o = out.getSampleDataRW(baseoffset+leaf);

	    for (index_t i=0; i < numComp; ++i) {
		o[INDEX3(i,0,0,numComp,3)] = f_000[i]*tmp0_0 + f_001[i]*tmp0_1 + f_100[i]*tmp0_3 + f_101[i]*tmp0_2;
		o[INDEX3(i,1,0,numComp,3)] = f_000[i]*tmp0_8 + f_010[i]*tmp0_13 + f_101[i]*tmp0_10 + f_111[i]*tmp0_11 + tmp0_12*(f_011[i] + f_110[i]) + tmp0_9*(f_001[i] + f_100[i]);
		o[INDEX3(i,2,0,numComp,3)] = f_000[i]*tmp0_28 + f_001[i]*tmp0_29 + f_100[i]*tmp0_31 + f_101[i]*tmp0_30;
		o[INDEX3(i,0,1,numComp,3)] = f_000[i]*tmp0_0 + f_001[i]*tmp0_1 + f_100[i]*tmp0_3 + f_101[i]*tmp0_2;
		o[INDEX3(i,1,1,numComp,3)] = f_000[i]*tmp0_14 + f_001[i]*tmp0_15 + f_010[i]*tmp0_21 + f_011[i]*tmp0_20 + f_100[i]*tmp0_17 + f_101[i]*tmp0_16 + f_110[i]*tmp0_19 + f_111[i]*tmp0_18;
		o[INDEX3(i,2,1,numComp,3)] = f_000[i]*tmp0_32 + f_001[i]*tmp0_33 + f_100[i]*tmp0_35 + f_101[i]*tmp0_34;
		o[INDEX3(i,0,2,numComp,3)] = f_000[i]*tmp0_4 + f_001[i]*tmp0_5 + f_100[i]*tmp0_7 + f_101[i]*tmp0_6;
		o[INDEX3(i,1,2,numComp,3)] = f_000[i]*tmp0_14 + f_001[i]*tmp0_17 + f_010[i]*tmp0_21 + f_011[i]*tmp0_19 + f_100[i]*tmp0_15 + f_101[i]*tmp0_16 + f_110[i]*tmp0_20 + f_111[i]*tmp0_18;
		o[INDEX3(i,2,2,numComp,3)] = f_000[i]*tmp0_28 + f_001[i]*tmp0_29 + f_100[i]*tmp0_31 + f_101[i]*tmp0_30;
		o[INDEX3(i,0,3,numComp,3)] = f_000[i]*tmp0_4 + f_001[i]*tmp0_5 + f_100[i]*tmp0_7 + f_101[i]*tmp0_6;
		o[INDEX3(i,1,3,numComp,3)] = f_000[i]*tmp0_22 + f_010[i]*tmp0_27 + f_101[i]*tmp0_24 + f_111[i]*tmp0_25 + tmp0_23*(f_001[i] + f_100[i]) + tmp0_26*(f_011[i] + f_110[i]);
		o[INDEX3(i,2,3,numComp,3)] = f_000[i]*tmp0_32 + f_001[i]*tmp0_33 + f_100[i]*tmp0_35 + f_101[i]*tmp0_34;
	    } /* end of component loop i */
	  }
        } /* end of face 2 */
        baseoffset+=face_cells[2].size();
        if (face_cells[3].size() > 0) {
#pragma omp for
          for (index_t leaf=0; leaf < face_cells[3].size(); ++leaf) {	  
	    const LeafInfo* li=face_cells[3][leaf]->leafinfo;		// we aren't iterating over all leaves this time	
	    
	    const double h0 = face_cells[3][leaf]->sides[0];
	    const double h1 = face_cells[3][leaf]->sides[1];
	    const double h2 = face_cells[3][leaf]->sides[2];
	    
            const double tmp0_22 = 0.16666666666666666667/h1;
            const double tmp0_16 = 0.16666666666666666667/h1;
            const double tmp0_33 = -0.78867513459481288225/h2;
            const double tmp0_0 = -0.21132486540518711775/h0;
            const double tmp0_21 = -0.62200846792814621559/h1;
            const double tmp0_17 = 0.16666666666666666667/h1;
            const double tmp0_1 = 0.78867513459481288225/h0;
            const double tmp0_20 = -0.16666666666666666667/h1;
            const double tmp0_14 = 0.044658198738520451079/h1;
            const double tmp0_2 = -0.78867513459481288225/h0;
            const double tmp0_27 = -0.62200846792814621559/h1;
            const double tmp0_15 = 0.62200846792814621559/h1;
            const double tmp0_3 = 0.21132486540518711775/h0;
            const double tmp0_26 = -0.16666666666666666667/h1;
            const double tmp0_12 = -0.16666666666666666667/h1;
            const double tmp0_25 = -0.044658198738520451079/h1;
            const double tmp0_13 = -0.044658198738520451079/h1;
            const double tmp0_24 = 0.62200846792814621559/h1;
            const double tmp0_10 = 0.044658198738520451079/h1;
            const double tmp0_11 = -0.62200846792814621559/h1;
            const double tmp0_34 = -0.21132486540518711775/h2;
            const double tmp0_35 = 0.78867513459481288225/h2;
            const double tmp0_8 = 0.16666666666666666667/h1;
            const double tmp0_29 = -0.21132486540518711775/h2;
            const double tmp0_9 = 0.62200846792814621559/h1;
            const double tmp0_30 = -0.78867513459481288225/h2;
            const double tmp0_28 = 0.78867513459481288225/h2;
            const double tmp0_32 = 0.21132486540518711775/h2;
            const double tmp0_31 = 0.21132486540518711775/h2;
            const double tmp0_18 = -0.16666666666666666667/h1;
            const double tmp0_4 = -0.78867513459481288225/h0;
            const double tmp0_19 = -0.044658198738520451079/h1;
            const double tmp0_5 = 0.21132486540518711775/h0;
            const double tmp0_6 = -0.21132486540518711775/h0;
            const double tmp0_23 = 0.044658198738520451079/h1;
            const double tmp0_7 = 0.78867513459481288225/h0;
	    
	    int buffcounter=0;
	    const double* src1=0;
	    int whichchild=-1;	    

	    GETHANGSAMPLE(face_cells[3][leaf], 0, f_000)
	    GETHANGSAMPLE(face_cells[3][leaf], 4, f_001)
	    GETHANGSAMPLE(face_cells[3][leaf], 5, f_101)
	    GETHANGSAMPLE(face_cells[3][leaf], 6, f_111)
	    GETHANGSAMPLE(face_cells[3][leaf], 2, f_110)
	    GETHANGSAMPLE(face_cells[3][leaf], 7, f_011)
	    GETHANGSAMPLE(face_cells[3][leaf], 3, f_010)
	    GETHANGSAMPLE(face_cells[3][leaf], 1, f_100)	    
	    
	    
/*	    const register double* f_000 = (li->pmap[0]>1)?in.getSampleDataRO(li->pmap[0]-2):HANG_INTERPOLATE(face_cells[3][leaf], 0);
	    const register double* f_001 = (li->pmap[1]>1)?in.getSampleDataRO(li->pmap[1]-2):HANG_INTERPOLATE(face_cells[3][leaf], 1);
	    const register double* f_101 = (li->pmap[5]>1)?in.getSampleDataRO(li->pmap[5]-2):HANG_INTERPOLATE(face_cells[3][leaf], 5);
	    const register double* f_111 = (li->pmap[6]>1)?in.getSampleDataRO(li->pmap[6]-2):HANG_INTERPOLATE(face_cells[3][leaf], 6);
	    const register double* f_110 = (li->pmap[7]>1)?in.getSampleDataRO(li->pmap[7]-2):HANG_INTERPOLATE(face_cells[3][leaf], 7);
	    const register double* f_011 = (li->pmap[2]>1)?in.getSampleDataRO(li->pmap[2]-2):HANG_INTERPOLATE(face_cells[3][leaf], 2);
	    const register double* f_010 = (li->pmap[3]>1)?in.getSampleDataRO(li->pmap[3]-2):HANG_INTERPOLATE(face_cells[3][leaf], 3);
	    const register double* f_100 = (li->pmap[4]>1)?in.getSampleDataRO(li->pmap[4]-2):HANG_INTERPOLATE(face_cells[3][leaf], 4);	   	*/    
	    
	    double* o = out.getSampleDataRW(baseoffset+leaf);	    
	    for (index_t i=0; i < numComp; ++i) {
		o[INDEX3(i,0,0,numComp,3)] = f_010[i]*tmp0_2 + f_011[i]*tmp0_0 + f_110[i]*tmp0_1 + f_111[i]*tmp0_3;
		o[INDEX3(i,1,0,numComp,3)] = f_000[i]*tmp0_11 + f_010[i]*tmp0_9 + f_101[i]*tmp0_13 + f_111[i]*tmp0_10 + tmp0_12*(f_001[i] + f_100[i]) + tmp0_8*(f_011[i] + f_110[i]);
		o[INDEX3(i,2,0,numComp,3)] = f_010[i]*tmp0_30 + f_011[i]*tmp0_28 + f_110[i]*tmp0_29 + f_111[i]*tmp0_31;
		o[INDEX3(i,0,1,numComp,3)] = f_010[i]*tmp0_2 + f_011[i]*tmp0_0 + f_110[i]*tmp0_1 + f_111[i]*tmp0_3;
		o[INDEX3(i,1,1,numComp,3)] = f_000[i]*tmp0_18 + f_001[i]*tmp0_19 + f_010[i]*tmp0_16 + f_011[i]*tmp0_14 + f_100[i]*tmp0_21 + f_101[i]*tmp0_20 + f_110[i]*tmp0_15 + f_111[i]*tmp0_17;
		o[INDEX3(i,2,1,numComp,3)] = f_010[i]*tmp0_34 + f_011[i]*tmp0_32 + f_110[i]*tmp0_33 + f_111[i]*tmp0_35;
		o[INDEX3(i,0,2,numComp,3)] = f_010[i]*tmp0_6 + f_011[i]*tmp0_4 + f_110[i]*tmp0_5 + f_111[i]*tmp0_7;
		o[INDEX3(i,1,2,numComp,3)] = f_000[i]*tmp0_18 + f_001[i]*tmp0_21 + f_010[i]*tmp0_16 + f_011[i]*tmp0_15 + f_100[i]*tmp0_19 + f_101[i]*tmp0_20 + f_110[i]*tmp0_14 + f_111[i]*tmp0_17;
		o[INDEX3(i,2,2,numComp,3)] = f_010[i]*tmp0_30 + f_011[i]*tmp0_28 + f_110[i]*tmp0_29 + f_111[i]*tmp0_31;
		o[INDEX3(i,0,3,numComp,3)] = f_010[i]*tmp0_6 + f_011[i]*tmp0_4 + f_110[i]*tmp0_5 + f_111[i]*tmp0_7;
		o[INDEX3(i,1,3,numComp,3)] = f_000[i]*tmp0_25 + f_010[i]*tmp0_23 + f_101[i]*tmp0_27 + f_111[i]*tmp0_24 + tmp0_22*(f_011[i] + f_110[i]) + tmp0_26*(f_001[i] + f_100[i]);
		o[INDEX3(i,2,3,numComp,3)] = f_010[i]*tmp0_34 + f_011[i]*tmp0_32 + f_110[i]*tmp0_33 + f_111[i]*tmp0_35;
	    } /* end of component loop i */
	  }
        } /* end of face 3 */
        baseoffset+=face_cells[3].size();
        if (face_cells[4].size() > 0) {
#pragma omp for
          for (index_t leaf=0; leaf < face_cells[4].size(); ++leaf) {	  
	    const LeafInfo* li=face_cells[4][leaf]->leafinfo;		// we aren't iterating over all leaves this time
	    
	    const double h0 = face_cells[4][leaf]->sides[0];
	    const double h1 = face_cells[4][leaf]->sides[1];
	    const double h2 = face_cells[4][leaf]->sides[2];        

	    const double tmp0_22 = -0.16666666666666666667/h2;
            const double tmp0_16 = -0.62200846792814621559/h2;
            const double tmp0_33 = 0.044658198738520451079/h2;
            const double tmp0_0 = -0.78867513459481288225/h0;
            const double tmp0_21 = 0.044658198738520451079/h2;
            const double tmp0_17 = -0.16666666666666666667/h2;
            const double tmp0_1 = 0.78867513459481288225/h0;
            const double tmp0_20 = 0.16666666666666666667/h2;
            const double tmp0_14 = 0.78867513459481288225/h1;
            const double tmp0_2 = 0.21132486540518711775/h0;
            const double tmp0_27 = 0.62200846792814621559/h2;
            const double tmp0_15 = 0.21132486540518711775/h1;
            const double tmp0_3 = -0.21132486540518711775/h0;
            const double tmp0_26 = 0.16666666666666666667/h2;
            const double tmp0_12 = -0.21132486540518711775/h1;
            const double tmp0_25 = -0.044658198738520451079/h2;
            const double tmp0_13 = -0.78867513459481288225/h1;
            const double tmp0_24 = -0.16666666666666666667/h2;
            const double tmp0_10 = 0.21132486540518711775/h1;
            const double tmp0_11 = 0.78867513459481288225/h1;
            const double tmp0_34 = 0.16666666666666666667/h2;
            const double tmp0_35 = 0.62200846792814621559/h2;
            const double tmp0_8 = -0.78867513459481288225/h1;
            const double tmp0_29 = 0.16666666666666666667/h2;
            const double tmp0_9 = -0.21132486540518711775/h1;
            const double tmp0_30 = -0.044658198738520451079/h2;
            const double tmp0_28 = 0.044658198738520451079/h2;
            const double tmp0_32 = -0.62200846792814621559/h2;
            const double tmp0_31 = -0.16666666666666666667/h2;
            const double tmp0_18 = -0.044658198738520451079/h2;
            const double tmp0_4 = -0.21132486540518711775/h0;
            const double tmp0_19 = 0.62200846792814621559/h2;
            const double tmp0_5 = 0.21132486540518711775/h0;
            const double tmp0_6 = 0.78867513459481288225/h0;
            const double tmp0_23 = -0.62200846792814621559/h2;
            const double tmp0_7 = -0.78867513459481288225/h0;
	    
	    int buffcounter=0;
	    const double* src1=0;
	    int whichchild=-1;	    
	    
	    GETHANGSAMPLE(face_cells[4][leaf], 0, f_000)
	    GETHANGSAMPLE(face_cells[4][leaf], 4, f_001)
	    GETHANGSAMPLE(face_cells[4][leaf], 5, f_101)
	    GETHANGSAMPLE(face_cells[4][leaf], 6, f_111)
	    GETHANGSAMPLE(face_cells[4][leaf], 2, f_110)
	    GETHANGSAMPLE(face_cells[4][leaf], 7, f_011)
	    GETHANGSAMPLE(face_cells[4][leaf], 3, f_010)
	    GETHANGSAMPLE(face_cells[4][leaf], 1, f_100)	    
	    
	    
/*	    const register double* f_000 = (li->pmap[0]>1)?in.getSampleDataRO(li->pmap[0]-2):HANG_INTERPOLATE(face_cells[4][leaf], 0);
	    const register double* f_001 = (li->pmap[1]>1)?in.getSampleDataRO(li->pmap[1]-2):HANG_INTERPOLATE(face_cells[4][leaf], 1);
	    const register double* f_101 = (li->pmap[5]>1)?in.getSampleDataRO(li->pmap[5]-2):HANG_INTERPOLATE(face_cells[4][leaf], 5);
	    const register double* f_111 = (li->pmap[6]>1)?in.getSampleDataRO(li->pmap[6]-2):HANG_INTERPOLATE(face_cells[4][leaf], 6);
	    const register double* f_110 = (li->pmap[7]>1)?in.getSampleDataRO(li->pmap[7]-2):HANG_INTERPOLATE(face_cells[4][leaf], 7);
	    const register double* f_011 = (li->pmap[2]>1)?in.getSampleDataRO(li->pmap[2]-2):HANG_INTERPOLATE(face_cells[4][leaf], 2);
	    const register double* f_010 = (li->pmap[3]>1)?in.getSampleDataRO(li->pmap[3]-2):HANG_INTERPOLATE(face_cells[4][leaf], 3);
	    const register double* f_100 = (li->pmap[4]>1)?in.getSampleDataRO(li->pmap[4]-2):HANG_INTERPOLATE(face_cells[4][leaf], 4);	   */	    
	    
	    double* o = out.getSampleDataRW(baseoffset+leaf);	    
	    for (index_t i=0; i < numComp; ++i) {
		o[INDEX3(i,0,0,numComp,3)] = f_000[i]*tmp0_0 + f_010[i]*tmp0_3 + f_100[i]*tmp0_1 + f_110[i]*tmp0_2;
		o[INDEX3(i,1,0,numComp,3)] = f_000[i]*tmp0_8 + f_010[i]*tmp0_11 + f_100[i]*tmp0_9 + f_110[i]*tmp0_10;
		o[INDEX3(i,2,0,numComp,3)] = f_000[i]*tmp0_16 + f_001[i]*tmp0_19 + f_110[i]*tmp0_18 + f_111[i]*tmp0_21 + tmp0_17*(f_010[i] + f_100[i]) + tmp0_20*(f_011[i] + f_101[i]);
		o[INDEX3(i,0,1,numComp,3)] = f_000[i]*tmp0_0 + f_010[i]*tmp0_3 + f_100[i]*tmp0_1 + f_110[i]*tmp0_2;
		o[INDEX3(i,1,1,numComp,3)] = f_000[i]*tmp0_12 + f_010[i]*tmp0_15 + f_100[i]*tmp0_13 + f_110[i]*tmp0_14;
		o[INDEX3(i,2,1,numComp,3)] = f_000[i]*tmp0_22 + f_001[i]*tmp0_26 + f_010[i]*tmp0_25 + f_011[i]*tmp0_28 + f_100[i]*tmp0_23 + f_101[i]*tmp0_27 + f_110[i]*tmp0_24 + f_111[i]*tmp0_29;
		o[INDEX3(i,0,2,numComp,3)] = f_000[i]*tmp0_4 + f_010[i]*tmp0_7 + f_100[i]*tmp0_5 + f_110[i]*tmp0_6;
		o[INDEX3(i,1,2,numComp,3)] = f_000[i]*tmp0_8 + f_010[i]*tmp0_11 + f_100[i]*tmp0_9 + f_110[i]*tmp0_10;
		o[INDEX3(i,2,2,numComp,3)] = f_000[i]*tmp0_22 + f_001[i]*tmp0_26 + f_010[i]*tmp0_23 + f_011[i]*tmp0_27 + f_100[i]*tmp0_25 + f_101[i]*tmp0_28 + f_110[i]*tmp0_24 + f_111[i]*tmp0_29;
		o[INDEX3(i,0,3,numComp,3)] = f_000[i]*tmp0_4 + f_010[i]*tmp0_7 + f_100[i]*tmp0_5 + f_110[i]*tmp0_6;
		o[INDEX3(i,1,3,numComp,3)] = f_000[i]*tmp0_12 + f_010[i]*tmp0_15 + f_100[i]*tmp0_13 + f_110[i]*tmp0_14;
		o[INDEX3(i,2,3,numComp,3)] = f_000[i]*tmp0_30 + f_001[i]*tmp0_33 + f_110[i]*tmp0_32 + f_111[i]*tmp0_35 + tmp0_31*(f_010[i] + f_100[i]) + tmp0_34*(f_011[i] + f_101[i]);
	    } /* end of component loop i */
	  }  
        } /* end of face 4 */
        baseoffset+=face_cells[4].size();
        if (face_cells[5].size() > 0) {
#pragma omp for
        for (index_t leaf=0; leaf < face_cells[5].size(); ++leaf) {	  
	    const LeafInfo* li=face_cells[5][leaf]->leafinfo;		// we aren't iterating over all leaves this time
	    
	    const double h0 = face_cells[5][leaf]->sides[0];
	    const double h1 = face_cells[5][leaf]->sides[1];
	    const double h2 = face_cells[5][leaf]->sides[2];        

            const double tmp0_22 = 0.16666666666666666667/h2;
            const double tmp0_16 = 0.62200846792814621559/h2;
            const double tmp0_33 = -0.044658198738520451079/h2;
            const double tmp0_0 = -0.78867513459481288225/h0;
            const double tmp0_21 = -0.16666666666666666667/h2;
            const double tmp0_17 = 0.16666666666666666667/h2;
            const double tmp0_1 = 0.78867513459481288225/h0;
            const double tmp0_20 = -0.044658198738520451079/h2;
            const double tmp0_14 = 0.21132486540518711775/h1;
            const double tmp0_2 = -0.21132486540518711775/h0;
            const double tmp0_27 = -0.16666666666666666667/h2;
            const double tmp0_15 = 0.78867513459481288225/h1;
            const double tmp0_3 = 0.21132486540518711775/h0;
            const double tmp0_26 = -0.16666666666666666667/h2;
            const double tmp0_12 = -0.21132486540518711775/h1;
            const double tmp0_25 = 0.16666666666666666667/h2;
            const double tmp0_13 = -0.78867513459481288225/h1;
            const double tmp0_24 = 0.044658198738520451079/h2;
            const double tmp0_10 = 0.78867513459481288225/h1;
            const double tmp0_11 = 0.21132486540518711775/h1;
            const double tmp0_34 = -0.62200846792814621559/h2;
            const double tmp0_35 = -0.16666666666666666667/h2;
            const double tmp0_8 = -0.78867513459481288225/h1;
            const double tmp0_29 = -0.62200846792814621559/h2;
            const double tmp0_9 = -0.21132486540518711775/h1;
            const double tmp0_30 = 0.044658198738520451079/h2;
            const double tmp0_28 = -0.044658198738520451079/h2;
            const double tmp0_32 = 0.62200846792814621559/h2;
            const double tmp0_31 = 0.16666666666666666667/h2;
            const double tmp0_18 = 0.044658198738520451079/h2;
            const double tmp0_4 = -0.21132486540518711775/h0;
            const double tmp0_19 = -0.62200846792814621559/h2;
            const double tmp0_5 = 0.21132486540518711775/h0;
            const double tmp0_6 = -0.78867513459481288225/h0;
            const double tmp0_23 = 0.62200846792814621559/h2;
            const double tmp0_7 = 0.78867513459481288225/h0;
	    
	    int buffcounter=0;
	    const double* src1=0;
	    int whichchild=-1;	    

	    GETHANGSAMPLE(face_cells[5][leaf], 0, f_000)
	    GETHANGSAMPLE(face_cells[5][leaf], 4, f_001)
	    GETHANGSAMPLE(face_cells[5][leaf], 5, f_101)
	    GETHANGSAMPLE(face_cells[5][leaf], 6, f_111)
	    GETHANGSAMPLE(face_cells[5][leaf], 2, f_110)
	    GETHANGSAMPLE(face_cells[5][leaf], 7, f_011)
	    GETHANGSAMPLE(face_cells[5][leaf], 3, f_010)
	    GETHANGSAMPLE(face_cells[5][leaf], 1, f_100)	    
	    
/*	    const register double* f_000 = (li->pmap[0]>1)?in.getSampleDataRO(li->pmap[0]-2):HANG_INTERPOLATE(face_cells[5][leaf], 0);
	    const register double* f_001 = (li->pmap[1]>1)?in.getSampleDataRO(li->pmap[1]-2):HANG_INTERPOLATE(face_cells[5][leaf], 1);
	    const register double* f_101 = (li->pmap[5]>1)?in.getSampleDataRO(li->pmap[5]-2):HANG_INTERPOLATE(face_cells[5][leaf], 5);
	    const register double* f_111 = (li->pmap[6]>1)?in.getSampleDataRO(li->pmap[6]-2):HANG_INTERPOLATE(face_cells[5][leaf], 6);
	    const register double* f_110 = (li->pmap[7]>1)?in.getSampleDataRO(li->pmap[7]-2):HANG_INTERPOLATE(face_cells[5][leaf], 7);
	    const register double* f_011 = (li->pmap[2]>1)?in.getSampleDataRO(li->pmap[2]-2):HANG_INTERPOLATE(face_cells[5][leaf], 2);
	    const register double* f_010 = (li->pmap[3]>1)?in.getSampleDataRO(li->pmap[3]-2):HANG_INTERPOLATE(face_cells[5][leaf], 3);
	    const register double* f_100 = (li->pmap[4]>1)?in.getSampleDataRO(li->pmap[4]-2):HANG_INTERPOLATE(face_cells[5][leaf], 4);	  */ 	    
	    
	    double* o = out.getSampleDataRW(baseoffset+leaf);	    
	    for (index_t i=0; i < numComp; ++i) {
		o[INDEX3(i,0,0,numComp,3)] = f_001[i]*tmp0_0 + f_011[i]*tmp0_2 + f_101[i]*tmp0_1 + f_111[i]*tmp0_3;
		o[INDEX3(i,1,0,numComp,3)] = f_001[i]*tmp0_8 + f_011[i]*tmp0_10 + f_101[i]*tmp0_9 + f_111[i]*tmp0_11;
		o[INDEX3(i,2,0,numComp,3)] = f_000[i]*tmp0_19 + f_001[i]*tmp0_16 + f_110[i]*tmp0_20 + f_111[i]*tmp0_18 + tmp0_17*(f_011[i] + f_101[i]) + tmp0_21*(f_010[i] + f_100[i]);
		o[INDEX3(i,0,1,numComp,3)] = f_001[i]*tmp0_0 + f_011[i]*tmp0_2 + f_101[i]*tmp0_1 + f_111[i]*tmp0_3;
		o[INDEX3(i,1,1,numComp,3)] = f_001[i]*tmp0_12 + f_011[i]*tmp0_14 + f_101[i]*tmp0_13 + f_111[i]*tmp0_15;
		o[INDEX3(i,2,1,numComp,3)] = f_000[i]*tmp0_26 + f_001[i]*tmp0_22 + f_010[i]*tmp0_28 + f_011[i]*tmp0_24 + f_100[i]*tmp0_29 + f_101[i]*tmp0_23 + f_110[i]*tmp0_27 + f_111[i]*tmp0_25;
		o[INDEX3(i,0,2,numComp,3)] = f_001[i]*tmp0_4 + f_011[i]*tmp0_6 + f_101[i]*tmp0_5 + f_111[i]*tmp0_7;
		o[INDEX3(i,1,2,numComp,3)] = f_001[i]*tmp0_8 + f_011[i]*tmp0_10 + f_101[i]*tmp0_9 + f_111[i]*tmp0_11;
		o[INDEX3(i,2,2,numComp,3)] = f_000[i]*tmp0_26 + f_001[i]*tmp0_22 + f_010[i]*tmp0_29 + f_011[i]*tmp0_23 + f_100[i]*tmp0_28 + f_101[i]*tmp0_24 + f_110[i]*tmp0_27 + f_111[i]*tmp0_25;
		o[INDEX3(i,0,3,numComp,3)] = f_001[i]*tmp0_4 + f_011[i]*tmp0_6 + f_101[i]*tmp0_5 + f_111[i]*tmp0_7;
		o[INDEX3(i,1,3,numComp,3)] = f_001[i]*tmp0_12 + f_011[i]*tmp0_14 + f_101[i]*tmp0_13 + f_111[i]*tmp0_15;
		o[INDEX3(i,2,3,numComp,3)] = f_000[i]*tmp0_33 + f_001[i]*tmp0_30 + f_110[i]*tmp0_34 + f_111[i]*tmp0_32 + tmp0_31*(f_011[i] + f_101[i]) + tmp0_35*(f_010[i] + f_100[i]);
	    } /* end of component loop i */
          } /* end of face 5 */
	}  
        /* GENERATOR SNIP_GRAD_FACES BOTTOM */
    } else if (out.getFunctionSpace().getTypeCode() == red_disc_faces) {
      throw BuckleyException("Haven't done reduced face function yet.");
#if 0      
        /* GENERATOR SNIP_GRAD_REDUCED_FACES TOP */
        if (m_faceOffset[0] > -1) {
            const double tmp0_0 = -0.25/h0;
            const double tmp0_4 = -0.5/h2;
            const double tmp0_1 = 0.25/h0;
            const double tmp0_5 = 0.5/h2;
            const double tmp0_2 = -0.5/h1;
            const double tmp0_3 = 0.5/h1;
#pragma omp parallel for
            for (index_t k2=0; k2 < m_NE2; ++k2) {
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    const register double* f_000 = in.getSampleDataRO(INDEX3(0,k1,k2, m_N0,m_N1));
                    const register double* f_001 = in.getSampleDataRO(INDEX3(0,k1,k2+1, m_N0,m_N1));
                    const register double* f_101 = in.getSampleDataRO(INDEX3(1,k1,k2+1, m_N0,m_N1));
                    const register double* f_011 = in.getSampleDataRO(INDEX3(0,k1+1,k2+1, m_N0,m_N1));
                    const register double* f_100 = in.getSampleDataRO(INDEX3(1,k1,k2, m_N0,m_N1));
                    const register double* f_110 = in.getSampleDataRO(INDEX3(1,k1+1,k2, m_N0,m_N1));
                    const register double* f_010 = in.getSampleDataRO(INDEX3(0,k1+1,k2, m_N0,m_N1));
                    const register double* f_111 = in.getSampleDataRO(INDEX3(1,k1+1,k2+1, m_N0,m_N1));
                    double* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE1));
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,3)] = tmp0_0*(f_000[i] + f_001[i] + f_010[i] + f_011[i]) + tmp0_1*(f_100[i] + f_101[i] + f_110[i] + f_111[i]);
                        o[INDEX3(i,1,0,numComp,3)] = tmp0_2*(f_000[i] + f_001[i]) + tmp0_3*(f_010[i] + f_011[i]);
                        o[INDEX3(i,2,0,numComp,3)] = tmp0_4*(f_000[i] + f_010[i]) + tmp0_5*(f_001[i] + f_011[i]);
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of k2 loop */
        } /* end of face 0 */
        if (m_faceOffset[1] > -1) {
            const double tmp0_0 = -0.25/h0;
            const double tmp0_4 = 0.5/h2;
            const double tmp0_1 = 0.25/h0;
            const double tmp0_5 = -0.5/h2;
            const double tmp0_2 = -0.5/h1;
            const double tmp0_3 = 0.5/h1;
#pragma omp parallel for
            for (index_t k2=0; k2 < m_NE2; ++k2) {
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    const register double* f_000 = in.getSampleDataRO(INDEX3(m_N0-2,k1,k2, m_N0,m_N1));
                    const register double* f_001 = in.getSampleDataRO(INDEX3(m_N0-2,k1,k2+1, m_N0,m_N1));
                    const register double* f_101 = in.getSampleDataRO(INDEX3(m_N0-1,k1,k2+1, m_N0,m_N1));
                    const register double* f_011 = in.getSampleDataRO(INDEX3(m_N0-2,k1+1,k2+1, m_N0,m_N1));
                    const register double* f_100 = in.getSampleDataRO(INDEX3(m_N0-1,k1,k2, m_N0,m_N1));
                    const register double* f_110 = in.getSampleDataRO(INDEX3(m_N0-1,k1+1,k2, m_N0,m_N1));
                    const register double* f_010 = in.getSampleDataRO(INDEX3(m_N0-2,k1+1,k2, m_N0,m_N1));
                    const register double* f_111 = in.getSampleDataRO(INDEX3(m_N0-1,k1+1,k2+1, m_N0,m_N1));
                    double* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE1));
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,3)] = tmp0_0*(f_000[i] + f_001[i] + f_010[i] + f_011[i]) + tmp0_1*(f_100[i] + f_101[i] + f_110[i] + f_111[i]);
                        o[INDEX3(i,1,0,numComp,3)] = tmp0_2*(f_100[i] + f_101[i]) + tmp0_3*(f_110[i] + f_111[i]);
                        o[INDEX3(i,2,0,numComp,3)] = tmp0_4*(f_101[i] + f_111[i]) + tmp0_5*(f_100[i] + f_110[i]);
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of k2 loop */
        } /* end of face 1 */
        if (m_faceOffset[2] > -1) {
            const double tmp0_0 = -0.5/h0;
            const double tmp0_4 = -0.5/h2;
            const double tmp0_1 = 0.5/h0;
            const double tmp0_5 = 0.5/h2;
            const double tmp0_2 = -0.25/h1;
            const double tmp0_3 = 0.25/h1;
#pragma omp parallel for
            for (index_t k2=0; k2 < m_NE2; ++k2) {
                for (index_t k0=0; k0 < m_NE0; ++k0) {
                    const register double* f_000 = in.getSampleDataRO(INDEX3(k0,0,k2, m_N0,m_N1));
                    const register double* f_001 = in.getSampleDataRO(INDEX3(k0,0,k2+1, m_N0,m_N1));
                    const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,0,k2+1, m_N0,m_N1));
                    const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,0,k2, m_N0,m_N1));
                    const register double* f_011 = in.getSampleDataRO(INDEX3(k0,1,k2+1, m_N0,m_N1));
                    const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,1,k2, m_N0,m_N1));
                    const register double* f_010 = in.getSampleDataRO(INDEX3(k0,1,k2, m_N0,m_N1));
                    const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,1,k2+1, m_N0,m_N1));
                    double* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE0));
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,3)] = tmp0_0*(f_000[i] + f_001[i]) + tmp0_1*(f_100[i] + f_101[i]);
                        o[INDEX3(i,1,0,numComp,3)] = tmp0_2*(f_000[i] + f_001[i] + f_100[i] + f_101[i]) + tmp0_3*(f_010[i] + f_011[i] + f_110[i] + f_111[i]);
                        o[INDEX3(i,2,0,numComp,3)] = tmp0_4*(f_000[i] + f_100[i]) + tmp0_5*(f_001[i] + f_101[i]);
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of k2 loop */
        } /* end of face 2 */
        if (m_faceOffset[3] > -1) {
            const double tmp0_0 = -0.5/h0;
            const double tmp0_4 = 0.5/h2;
            const double tmp0_1 = 0.5/h0;
            const double tmp0_5 = -0.5/h2;
            const double tmp0_2 = 0.25/h1;
            const double tmp0_3 = -0.25/h1;
#pragma omp parallel for
            for (index_t k2=0; k2 < m_NE2; ++k2) {
                for (index_t k0=0; k0 < m_NE0; ++k0) {
                    const register double* f_011 = in.getSampleDataRO(INDEX3(k0,m_N1-1,k2+1, m_N0,m_N1));
                    const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,m_N1-1,k2, m_N0,m_N1));
                    const register double* f_010 = in.getSampleDataRO(INDEX3(k0,m_N1-1,k2, m_N0,m_N1));
                    const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,m_N1-1,k2+1, m_N0,m_N1));
                    const register double* f_000 = in.getSampleDataRO(INDEX3(k0,m_N1-2,k2, m_N0,m_N1));
                    const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,m_N1-2,k2+1, m_N0,m_N1));
                    const register double* f_001 = in.getSampleDataRO(INDEX3(k0,m_N1-2,k2+1, m_N0,m_N1));
                    const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,m_N1-2,k2, m_N0,m_N1));
                    double* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE0));
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,3)] = tmp0_0*(f_010[i] + f_011[i]) + tmp0_1*(f_110[i] + f_111[i]);
                        o[INDEX3(i,1,0,numComp,3)] = tmp0_2*(f_010[i] + f_011[i] + f_110[i] + f_111[i]) + tmp0_3*(f_000[i] + f_001[i] + f_100[i] + f_101[i]);
                        o[INDEX3(i,2,0,numComp,3)] = tmp0_4*(f_011[i] + f_111[i]) + tmp0_5*(f_010[i] + f_110[i]);
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of k2 loop */
        } /* end of face 3 */
        if (m_faceOffset[4] > -1) {
            const double tmp0_0 = -0.5/h0;
            const double tmp0_4 = -0.25/h2;
            const double tmp0_1 = 0.5/h0;
            const double tmp0_5 = 0.25/h2;
            const double tmp0_2 = -0.5/h1;
            const double tmp0_3 = 0.5/h1;
#pragma omp parallel for
            for (index_t k1=0; k1 < m_NE1; ++k1) {
                for (index_t k0=0; k0 < m_NE0; ++k0) {
                    const register double* f_000 = in.getSampleDataRO(INDEX3(k0,k1,0, m_N0,m_N1));
                    const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,k1,0, m_N0,m_N1));
                    const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,k1+1,0, m_N0,m_N1));
                    const register double* f_010 = in.getSampleDataRO(INDEX3(k0,k1+1,0, m_N0,m_N1));
                    const register double* f_001 = in.getSampleDataRO(INDEX3(k0,k1,1, m_N0,m_N1));
                    const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,k1,1, m_N0,m_N1));
                    const register double* f_011 = in.getSampleDataRO(INDEX3(k0,k1+1,1, m_N0,m_N1));
                    const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,k1+1,1, m_N0,m_N1));
                    double* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE0));
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,3)] = tmp0_0*(f_000[i] + f_010[i]) + tmp0_1*(f_100[i] + f_110[i]);
                        o[INDEX3(i,1,0,numComp,3)] = tmp0_2*(f_000[i] + f_100[i]) + tmp0_3*(f_010[i] + f_110[i]);
                        o[INDEX3(i,2,0,numComp,3)] = tmp0_4*(f_000[i] + f_010[i] + f_100[i] + f_110[i]) + tmp0_5*(f_001[i] + f_011[i] + f_101[i] + f_111[i]);
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of k1 loop */
        } /* end of face 4 */
        if (m_faceOffset[5] > -1) {
            const double tmp0_0 = -0.5/h0;
            const double tmp0_4 = 0.25/h2;
            const double tmp0_1 = 0.5/h0;
            const double tmp0_5 = -0.25/h2;
            const double tmp0_2 = -0.5/h1;
            const double tmp0_3 = 0.5/h1;
#pragma omp parallel for
            for (index_t k1=0; k1 < m_NE1; ++k1) {
                for (index_t k0=0; k0 < m_NE0; ++k0) {
                    const register double* f_001 = in.getSampleDataRO(INDEX3(k0,k1,m_N2-1, m_N0,m_N1));
                    const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,k1,m_N2-1, m_N0,m_N1));
                    const register double* f_011 = in.getSampleDataRO(INDEX3(k0,k1+1,m_N2-1, m_N0,m_N1));
                    const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,k1+1,m_N2-1, m_N0,m_N1));
                    const register double* f_000 = in.getSampleDataRO(INDEX3(k0,k1,m_N2-2, m_N0,m_N1));
                    const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,k1,m_N2-2, m_N0,m_N1));
                    const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,k1+1,m_N2-2, m_N0,m_N1));
                    const register double* f_010 = in.getSampleDataRO(INDEX3(k0,k1+1,m_N2-2, m_N0,m_N1));
                    double* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE0));
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,3)] = tmp0_0*(f_001[i] + f_011[i]) + tmp0_1*(f_101[i] + f_111[i]);
                        o[INDEX3(i,1,0,numComp,3)] = tmp0_2*(f_001[i] + f_101[i]) + tmp0_3*(f_011[i] + f_111[i]);
                        o[INDEX3(i,2,0,numComp,3)] = tmp0_4*(f_001[i] + f_011[i] + f_101[i] + f_111[i]) + tmp0_5*(f_000[i] + f_010[i] + f_100[i] + f_110[i]);
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of k1 loop */
        } /* end of face 5 */
        /* GENERATOR SNIP_GRAD_REDUCED_FACES BOTTOM */
#endif	
      } else {
        err=true;
      }	
    }  // end omp parallel
    if (err)
    {
        stringstream msg;
        msg << "setToGradient() not implemented for "
            << functionSpaceTypeAsString(out.getFunctionSpace().getTypeCode());
        throw BuckleyException(msg.str());
    }
}

	    
#undef GETHANGSAMPLE



#if 0
void BuckleyDomain::setToGradient(escript::Data& grad, const escript::Data& arg) const
{
    // sanity checks
    if (!grad.getFunctionSpace().upToDate() ||  !arg.getFunctionSpace().upToDate())
    {
        throw BuckleyException("Stale functionspaces passed to setToGradient");
    }
    if ((grad.getFunctionSpace().getDomain().get()!=this) || (arg.getFunctionSpace().getDomain().get()!=this))
    {
        throw BuckleyException("Domain mismatch in setToGradient");
    }
    if ((grad.getFunctionSpace().getTypeCode()!=discfn) || (arg.getFunctionSpace().getTypeCode()!=ctsfn))
    {
        throw BuckleyException("Invalid functionspaces passed to setToGradient");  
    }
    if ((grad.getDataPointShape()!=MAT3X3) || (arg.getDataPointShape()!=VECTOR3))
    {
        throw BuckleyException("Expected shapes for grad inputs are (3,3) and (3,)");
    }
    if (modified)
    {
        processMods();
    }
    
       
    double* dest=grad.getSampleDataRW(0);
    for (int i=0;i<ot.leafCount();++i)
    {
        double quadpts[NUM_QUAD][3];
	for (int j=0;j<NUM_QUAD;++j)
	{
	    leaves[i]->quadCoords(j, quadpts[j][0], quadpts[j][1], quadpts[j][2]);
	}
        double values[8*3];	// the values of the corner points of this cell
        double values2[8*3];
	double valuetemp[8];
	double safepoint[3];
	const LeafInfo* li=leaves[i]->leafinfo;
	int whichchild=-1;	// which child of our parent are we
	for (int j=0;j<8;++j)
	{
	    // now we need to collate the values at the corner points
	    if (li->pmap[j]<2)	// hanging node so we will need to interpolate
	    {
	        // We know that if a node is hanging, then there must be at least one face
		// touching that edge which is less refined.
		// This means that nothing else touching this edge can be further refined.
		// based on this, if we know which child we are, then we can find the opposite
		// corner which we can use to interpolate with
		
		// to understand this code please see the hangskip table in FaceConsts.h
		if (whichchild==-1)   // we don't know which child we are yet
		{
		    whichchild=leaves[i]->whichChild();
		    const double* src=const_cast<escript::Data&>(arg).getSampleDataRO(li->pmap[whichchild]-2);
		    safepoint[0]=src[0];	// any hanging nodes in this cell can use the same
		    safepoint[1]=src[1];	// safepoint as the start of interpolation
		    safepoint[2]=src[2];
		}
		// The other end of the interpolation line matches the number of the hanging node
		const double* src=const_cast<escript::Data&>(arg).getSampleDataRO(leaves[i]->parent->kids[j]->leafinfo->pmap[j]-2);
		values[j*3+0]=(safepoint[0]+src[0])/2;
		values[j*3+1]=(safepoint[1]+src[1])/2;
		values[j*3+2]=(safepoint[2]+src[2])/2;
	    }
	    else
	    {
	        // for ctsfn each point is a sample
		const double* src=const_cast<escript::Data&>(arg).getSampleDataRO(li->pmap[j]-2);
	        values[j*3+0]=src[0];
		values[j*3+1]=src[1];
	        values[j*3+2]=src[2];
		// sort out that const_cast situation.
		// Why isn't getSampleDataRO  const?
		// Do I need to be more agressive with mutables?
		// Finley dodges this issue by asking for DataC and using that --- which is a constant operation
	    }
	}
	
	valuetemp[0]=values[0];
	valuetemp[1]=values[1];
	valuetemp[2]=values[2];
	
	// now we have all the values we can interpolate
        
	interpolateElementFromCtsToDisc(li, 3, values, values2);
        // ok so now values2 has the vectors for each of the points in it
	// first we want to translate so that the value from point_0 is at the origin	
	for (int j=0;j<8;++j)
	{
	    values2[j*3]-=valuetemp[0];
	    values2[j*3+1]-=valuetemp[1];
	    values2[j*3+2]-=valuetemp[2];  
	}
	// still need to compute each matrix

        double dx=leaves[i]->sides[0]/4;
	double dy=leaves[i]->sides[1]/4;
	double dz=leaves[i]->sides[2]/4;
	using namespace escript::DataTypes;
	
	for (int j=0; j<8;++j)
	{
	    double lx=dx, ly=dy, lz=dz;
	    if (j>3){lz*=3;}
	    if (j==1 || j==2 || j==5 || j==6) {lx*=3;}
	    if (j==2 || j==3 || j==6 || j==7) {ly*=3;}
	    
	    // is an element with a kink in it a problem here?
	    
	    dest[getRelIndex(MAT3X3, 0, 0)]=values2[j*3+0]/lx;
	    dest[getRelIndex(MAT3X3, 0, 1)]=values2[j*3+0]/ly;
	    dest[getRelIndex(MAT3X3, 0, 2)]=values2[j*3+0]/lz;
	    dest[getRelIndex(MAT3X3, 1, 0)]=values2[j*3+1]/lx;
	    dest[getRelIndex(MAT3X3, 1, 1)]=values2[j*3+1]/ly;
	    dest[getRelIndex(MAT3X3, 1, 2)]=values2[j*3+1]/lz;
	    dest[getRelIndex(MAT3X3, 2, 0)]=values2[j*3+2]/lx;
	    dest[getRelIndex(MAT3X3, 2, 1)]=values2[j*3+2]/ly;
	    dest[getRelIndex(MAT3X3, 2, 2)]=values2[j*3+2]/lz;
	    dest+=9;
	}

//	throw BuckleyException("This still does not poopulate the arrays corresponding to each point.");
	

//	dest+=3*3*8;		// move on to next element
    }
    
    
}

#endif
