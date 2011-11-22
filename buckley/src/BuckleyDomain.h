
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


#ifndef rdomain_20111025_H
#define rdomain_20111025_H

#include <escript/Data.h>
#include <escript/AbstractContinuousDomain.h>
#include <escript/AbstractSystemMatrix.h>

#include "OctTree.h"
#include "system_dep.h"
#include "BuckleyException.h"

namespace buckley {

 

class BuckleyDomain : public escript::AbstractContinuousDomain 
{

 public:

  BUCKLEY_DLL_API 
  BuckleyDomain(double xside=1, double yside=1, double zside=1);

  BUCKLEY_DLL_API 
  virtual ~BuckleyDomain();

  /**
     \brief
     Returns true if the given integer is a valid function space type
     for this domain.
  */
  BUCKLEY_DLL_API 
  virtual bool isValidFunctionSpaceType(int functionSpaceType) const;

  /**
     \brief
     Return a description for this domain
  */
  BUCKLEY_DLL_API 
  virtual std::string getDescription() const;

  /**
     \brief
     Return a continuous FunctionSpace code
  */
  BUCKLEY_DLL_API 
  virtual int getContinuousFunctionCode() const;

  /**
     \brief
     Return a continuous on reduced order FunctionSpace code
  */
  BUCKLEY_DLL_API 
  virtual int getReducedContinuousFunctionCode() const;

  /**
     \brief
     Return a function FunctionSpace code
  */
  BUCKLEY_DLL_API 
  virtual int getFunctionCode() const;

  /**
     \brief
     Return a function FunctionSpace code with reduced integration order
  */
  BUCKLEY_DLL_API 
  virtual int getReducedFunctionCode() const;

  /**
     \brief
     Return a function on boundary FunctionSpace code
  */
  BUCKLEY_DLL_API 
  virtual int getFunctionOnBoundaryCode() const;

  /**
     \brief
     Return a code for a function on boundary FunctionSpace with reduced integration order
  */
  BUCKLEY_DLL_API 
  virtual int getReducedFunctionOnBoundaryCode() const;


  /**
     \brief
     Return a FunctionOnContactZero code
  */
  BUCKLEY_DLL_API 
  virtual int getFunctionOnContactZeroCode() const;

  /**
     \brief
     Return a FunctionOnContactZero for reduced integration order code 
  */
  BUCKLEY_DLL_API 
  virtual int getReducedFunctionOnContactZeroCode() const;

  /**
     \brief
     Return a FunctionOnContactOne code
  */
  BUCKLEY_DLL_API 
  virtual int getFunctionOnContactOneCode() const;

  /**
     \brief
     Return a FunctionOnContactOne for reduced integration order code
  */
  BUCKLEY_DLL_API 
  virtual int getReducedFunctionOnContactOneCode() const;

  /**
     \brief
     Return a Solution code
  */
  BUCKLEY_DLL_API 
  virtual int getSolutionCode() const;

  /**
     \brief
     Return a ReducedSolution code
  */
  BUCKLEY_DLL_API 
  virtual int getReducedSolutionCode() const;

  /**
     \brief
     Return a DiracDeltaFunctions code
  */
  BUCKLEY_DLL_API 
  virtual int getDiracDeltaFunctionsCode() const;

  BUCKLEY_DLL_API
  virtual escript::Data getX() const;
  
  BUCKLEY_DLL_API
  virtual void setToX(escript::Data& arg) const  ;
  
  
  BUCKLEY_DLL_API
  virtual int getDim() const
  {
      return 3;
  }
 
  
  /**
     \brief
     return the identifier of the matrix type to be used for the global
     stiffness matrix when a particular solver package
     and symmetric matrix is used.
  */
  BUCKLEY_DLL_API 
  virtual int getSystemMatrixTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const;

  /**
     \brief
     return the identifier of the matrix type to be used for the global
     stiffness matrix when a particular solver package
     and symmetric matrix is used.
  */
  BUCKLEY_DLL_API 
  virtual int getTransportTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const;

  /**
     \brief
     copies the integrals of the function defined by arg into integrals.
     arg has to be defined on this.
     has to be implemented by the Domain Adapter.
  */
  BUCKLEY_DLL_API 
  virtual void setToIntegrals(std::vector<double>& integrals,const escript::Data& arg) const;

//  /**
//     \brief
//     Return the domain as const AbstractContinuousDomain&
//  */
//   BUCKLEY_DLL_API 
//   static const  AbstractContinuousDomain& asAbstractContinuousDomain(const AbstractDomain& domain);




  /**
     \brief
     adds a PDE onto the stiffness matrix mat and a rhs 
  */
  BUCKLEY_DLL_API
  virtual void addPDEToSystem(
                     escript::AbstractSystemMatrix& mat, escript::Data& rhs,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C, 
                     const escript::Data& D, const escript::Data& X, const escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact, const escript::Data& y_contact, 
                     const escript::Data& d_dirac, const escript::Data& y_dirac) const;

// We do not require this method at this level since the python side checks to ensure it exists
// before calling it.

//  /**
//     \brief
//     adds a PDE onto the lumped stiffness matrix matrix
//  */
//  BUCKLEY_DLL_API
//  virtual void addPDEToLumpedSystem(
//                     escript::Data& mat,
//                     const escript::Data& D, 
//                     const escript::Data& d) const;

  /**
     \brief
     adds a PDE onto the stiffness matrix mat and a rhs 
  */
  BUCKLEY_DLL_API
  virtual void addPDEToRHS(escript::Data& rhs,
                     const escript::Data& X, const escript::Data& Y,
                     const escript::Data& y, const escript::Data& y_contact, const escript::Data& y_dirac) const;
  /**
     \brief
     adds a PDE onto a transport problem
  */

  BUCKLEY_DLL_API
  virtual void addPDEToTransportProblem(
                     escript::AbstractTransportProblem& tp, escript::Data& source, 
                     const escript::Data& M,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C,const  escript::Data& D,
                     const  escript::Data& X,const  escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact,const escript::Data& y_contact,
                     const escript::Data& d_dirac,const escript::Data& y_dirac) const;

  /**
     \brief
    creates a SystemMatrixAdapter stiffness matrix and initializes it with zeros:
  */
  BUCKLEY_DLL_API
  virtual escript::ASM_ptr newSystemMatrix(
                      const int row_blocksize,
                      const escript::FunctionSpace& row_functionspace,
                      const int column_blocksize,
                      const escript::FunctionSpace& column_functionspace,
                      const int type) const;
  /**
   \brief 
    creates a TransportProblemAdapter 

  */

  BUCKLEY_DLL_API
  virtual escript::ATP_ptr newTransportProblem(
                      const bool useBackwardEuler,
                      const int blocksize,
                      const escript::FunctionSpace& functionspace,
                      const int type) const;

  /**
     \brief
     Return the number of data points summed across all MPI processes
  */
  BUCKLEY_DLL_API
  virtual int getNumDataPointsGlobal() const;

  /**
     \brief
     Return the number of data points per sample, and the number of samples as a pair.
     \param functionSpaceCode Input -
  */
  BUCKLEY_DLL_API
  virtual std::pair<int,int> getDataShape(int functionSpaceCode) const;

  /**
     \brief
     assigns new location to the domain
  */
  BUCKLEY_DLL_API
  virtual void setNewX(const escript::Data& arg);

  
  BUCKLEY_DLL_API
  virtual void refineAll(unsigned min_depth);
  
  /**
     \brief
     \param full
  */
  BUCKLEY_DLL_API
  virtual void Print_Mesh_Info(const bool full=false) const;
  
  BUCKLEY_DLL_API
  bool operator==(const AbstractDomain& other) const;
  
  BUCKLEY_DLL_API
  bool operator!=(const AbstractDomain& other) const;
  
  BUCKLEY_DLL_API
  const int* borrowSampleReferenceIDs(int functionSpaceType) const;
  
  BUCKLEY_DLL_API
  void refinePoint(double x, double y, double z, unsigned int desdepth);
  
  BUCKLEY_DLL_API
  void collapseAll(unsigned max_depth);
  
  BUCKLEY_DLL_API
  void collapsePoint(double x, double y, double z, unsigned int desdepth); 

  BUCKLEY_DLL_API
  unsigned  getGeneration() const; 
  
  BUCKLEY_DLL_API
  virtual std::string functionSpaceTypeAsString(int functionSpaceType) const;  
  

 protected:

 private:
   
    void processMods() const;
   
    OctTree ot;
    mutable bool modified;
    mutable const OctCell** leaves;
    mutable unkid numpts;	// number of independent (non-hanging) verticies
    mutable int* samplerefids;
    unsigned m_generation;

    
};

} // end of namespace

#endif
