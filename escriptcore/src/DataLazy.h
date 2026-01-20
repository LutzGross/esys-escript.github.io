
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

#ifndef __ESCRIPT_DATALAZY_H__
#define __ESCRIPT_DATALAZY_H__

#include "system_dep.h"
#include "DataAbstract.h"
#include "ArrayOps.h"		// for tensor_binary_op
#include "DataVector.h"		// for ElementType
#include "ES_optype.h"

#include <string>

//#define LAZY_NODE_STORAGE

namespace escript {


/**
\class escript::DataLazy
\brief Wraps an expression tree of other DataObjects.
The data will be evaluated when required.


NOTE: This class assumes that the Data being pointed at are immutable.
*/

class DataLazy;

typedef POINTER_WRAPPER_CLASS(DataLazy) DataLazy_ptr;
typedef POINTER_WRAPPER_CLASS(const DataLazy) const_DataLazy_ptr;

class DataLazy : public DataAbstract
{

typedef DataAbstract parent;
typedef DataTypes::ShapeType ShapeType;

public:
  /**
  \brief Create an IDENTITY DataLazy for the given DataAbstract.
  \param p DataAbstract to be wrapped.
  \throws DataException if p is lazy data or it is not constant, tagged or expanded.
  */
  ESCRIPT_DLL_API
  DataLazy(DataAbstract_ptr p);


  /**
  \brief Produce a DataLazy for a unary operation.
  \param left DataAbstract to be operated on.
  \param op unary operation to perform.
  \throws DataException if op is not a unary operation or if p cannot be converted to a DataLazy.
  Note that IDENTITY is not considered a unary operation.
  */
  ESCRIPT_DLL_API
  DataLazy(DataAbstract_ptr left, ES_optype op);

  /**
  \brief Produce a DataLazy for a unary operation.
  \param left DataAbstract to be operated on.
  \param op unary operation to perform.
  \param tol tolerance for operation
  \throws DataException if op is not a unary operation or if p cannot be converted to a DataLazy.
  Note that IDENTITY is not considered a unary operation.
  */
  ESCRIPT_DLL_API
  DataLazy(DataAbstract_ptr left, ES_optype op, double tol);

  /**
  \brief Produce a DataLazy for a unary operation which requires a parameter.
  \param left DataAbstract to be operated on.
  \param op unary operation to perform.
  \param axis_offset the parameter for the operation
  \throws DataException if op is not a unary operation or if p cannot be converted to a DataLazy.
  Note that IDENTITY is not considered a unary operation.
  */
  ESCRIPT_DLL_API  
  DataLazy(DataAbstract_ptr left, ES_optype op, int axis_offset);


  /**
  \brief Produce a DataLazy for a binary operation.
  \param left left operand
  \param right right operand
  \param op unary operation to perform.
  \throws DataException if op is not a binary operation or if left or right cannot be converted to a DataLazy.
  */
  ESCRIPT_DLL_API
  DataLazy(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op);

  /**
  \brief Produce a DataLazy for a binary operation with additional paramters.
  \param left left operand
  \param right right operand
  \param op unary operation to perform.
  \param axis_offset 
  \param transpose  
  \throws DataException if op is not a binary operation requiring parameters or if left or right cannot be converted to a DataLazy.
  */
  ESCRIPT_DLL_API
  DataLazy(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op, int axis_offset, int transpose);

  /**
  \brief Produce a DataLazy for a unary operation which requires two integer parameters.
  \param left DataAbstract to be operated on.
  \param op unary operation to perform.
  \param axis0 the first parameter for the operation
  \param axis1 the second parameter for the operation
  \throws DataException if op is not a unary operation or if p cannot be converted to a DataLazy.
  Note that IDENTITY is not considered a unary operation.
  */
  ESCRIPT_DLL_API
  DataLazy(DataAbstract_ptr left, ES_optype op, const int axis0, const int axis1);

  /**
  \brief Produce a DataLazy for a unary operation which requires two integer parameters.
  \param mask scalar mask to select values.
  \param left DataAbstract to use for true mask.
  \param right DataAbstract to use for false mask.
  */
  ESCRIPT_DLL_API
  DataLazy(DataAbstract_ptr mask, DataAbstract_ptr left, DataAbstract_ptr right/*, double tol*/);

  ESCRIPT_DLL_API
  ~DataLazy();

  /**
  \brief Evaluate the lazy expression.
  \return A DataReady with the value of the lazy expresion.
  */
  ESCRIPT_DLL_API
  DataReady_ptr 
  resolve();

  ESCRIPT_DLL_API
  std::string
  toString() const;

  ESCRIPT_DLL_API
  DataAbstract* 
  deepCopy() const;

  ESCRIPT_DLL_API
  DataAbstract* 
  zeroedCopy() const;  

  /**
     \brief
     This method throws an exception. It does not really make sense to ask this question of lazy data.
  */
  ESCRIPT_DLL_API
  DataTypes::RealVectorType::size_type
  getLength() const;


  ESCRIPT_DLL_API
  DataAbstract*
  getSlice(const DataTypes::RegionType& region) const;


  DataTypes::RealVectorType::size_type 
  getPointOffset(int sampleNo,
                 int dataPointNo) const;

  DataTypes::RealVectorType::size_type 
  getPointOffset(int sampleNo,
                 int dataPointNo);

  /**
    \return the largest samplesize required to evaluate the expression.
  */
  ESCRIPT_DLL_API
  size_t
  getMaxSampleSize() const;

   /**
  \brief Compute the value of the expression for the given sample.
  \return Vector which stores the value of the subexpression for the given sample.
  \param sampleNo Sample number to evaluate.
  \param roffset (output parameter) the offset in the return vector where the result begins.

  The return value will be an existing vector so do not deallocate it.
  */
  ESCRIPT_DLL_API
  const DataTypes::RealVectorType*
  resolveSample(int sampleNo, size_t& roffset) const; 
  
  ESCRIPT_DLL_API
  const DataTypes::CplxVectorType*
  resolveTypedSample(int sampleNo, size_t& roffset, DataTypes::cplx_t dummy) const; 
  
  ESCRIPT_DLL_API
  const DataTypes::RealVectorType*
  resolveTypedSample(int sampleNo, size_t& roffset, DataTypes::real_t dummy) const; 
  

  /**
  \brief if resolve() was called would it produce expanded data.
  */
  ESCRIPT_DLL_API
  bool
  actsExpanded() const;

  /**
     \brief Produces an IDENTITY DataLazy containing zero.
     The result will have the same shape and functionspace as before.
  */
  ESCRIPT_DLL_API
  virtual void
  setToZero();

  
  ESCRIPT_DLL_API
  void
  resolveGroupWorker(std::vector<DataLazy*>& dats);


private:
  int* m_sampleids;		// may be NULL
  mutable DataTypes::RealVectorType m_samples_r;
  mutable DataTypes::CplxVectorType m_samples_c;     
    
  mutable DataReady_ptr m_id;	//  For IDENTITY nodes, stores a wrapped value.
  mutable DataLazy_ptr m_left, m_right, m_mask;	// operands for operation.
  mutable ES_optype m_op;	// operation to perform.
  mutable ES_opgroup m_opgroup; // type of operation to perform

  size_t m_samplesize;	// number of values required to store a sample

  char m_readytype;	// E for expanded, T for tagged, C for constant

  int m_axis_offset;	// required extra info for general tensor product
  int m_transpose;	// offset and transpose are used for swapaxes as well
  int m_SL, m_SM, m_SR;	// computed properties used in general tensor product


  double m_tol;		// required extra info for <>0 and ==0

  mutable size_t m_children;
  mutable size_t m_height;

 

  /**
  Allocates sample storage at each node
  */
  void LazyNodeSetup();


  const DataTypes::RealVectorType*
  resolveNodeUnary(int tid, int sampleNo, size_t& roffset) const;
  
  const DataTypes::CplxVectorType*
  resolveNodeUnaryCplx(int tid, int sampleNo, size_t& roffset) const;  


  const DataTypes::RealVectorType*
  resolveNodeReduction(int tid, int sampleNo, size_t& roffset) const;  

  const DataTypes::CplxVectorType*
  resolveNodeReductionCplx(int tid, int sampleNo, size_t& roffset) const;  
  
  const DataTypes::RealVectorType*
  resolveNodeSample(int tid, int sampleNo, size_t& roffset) const;
  
  const DataTypes::CplxVectorType*
  resolveNodeSampleCplx(int tid, int sampleNo, size_t& roffset) const;  

  const DataTypes::RealVectorType*
  resolveNodeBinary(int tid, int sampleNo, size_t& roffset) const;
  
  const DataTypes::CplxVectorType*
  resolveNodeBinaryCplx(int tid, int sampleNo, size_t& roffset) const;
  

  const DataTypes::RealVectorType*
  resolveNodeNP1OUT(int tid, int sampleNo, size_t& roffset) const;

  const DataTypes::CplxVectorType*
  resolveNodeNP1OUTCplx(int tid, int sampleNo, size_t& roffset) const;
  
  const DataTypes::RealVectorType*
  resolveNodeNP1OUT_P(int tid, int sampleNo, size_t& roffset) const;

  const DataTypes::CplxVectorType*
  resolveNodeNP1OUT_PCplx(int tid, int sampleNo, size_t& roffset) const;
  
  
  const DataTypes::RealVectorType*
  resolveNodeTProd(int tid, int sampleNo, size_t& roffset) const;

  const DataTypes::CplxVectorType*
  resolveNodeTProdCplx(int tid, int sampleNo, size_t& roffset) const;

  
  const DataTypes::RealVectorType*
  resolveNodeNP1OUT_2P(int tid, int sampleNo, size_t& roffset) const;

  const DataTypes::CplxVectorType*
  resolveNodeNP1OUT_2PCplx(int tid, int sampleNo, size_t& roffset) const;
  
  const DataTypes::RealVectorType*
  resolveNodeCondEval(int tid, int sampleNo, size_t& roffset) const;

  const DataTypes::CplxVectorType*
  resolveNodeCondEvalCplx(int tid, int sampleNo, size_t& roffset) const;
  
  
  const DataTypes::CplxVectorType*
  resolveNodeUnary_C(int tid, int sampleNo, size_t& roffset) const;
  
  /**
  Does the work for toString. 
  */
  void
  intoString(std::ostringstream& oss) const;

  /**
  \warning internal use only!!
  */
  void
  intoTreeString(std::ostringstream& oss,std::string indent) const;

  /**
   \brief Converts the DataLazy into an IDENTITY storing the value of the expression.
   This method uses the original methods on the Data class to evaluate the expressions.
   For this reason, it should not be used on DataExpanded instances. (To do so would defeat
   the purpose of using DataLazy in the first place).
  */
  void
  collapse() const;		// converts the node into an IDENTITY node


  /**
  \brief Evaluates the expression using methods on Data.
  This does the work for the collapse method.
  For reasons of efficiency do not call this method on DataExpanded nodes.
  */
  DataReady_ptr
  collapseToReady() const;

  /**
  \brief resolve the expression can store it in the current node
  The current node will be converted to an identity node.
  */
  void
  resolveToIdentity();

  /**
  \brief helper method for resolveToIdentity and the identity constructor
  */
  void 
  makeIdentity(const DataReady_ptr& p);


  /**
  \brief resolve to a ReadyData object using storage at nodes
  */
  DataReady_ptr
  resolveNodeWorker();
  
  DataReady_ptr
  resolveNodeWorkerCplx();  

};

// If an expression is already complex, return the same expression.
// Otherwise, return the old expression with a promote operation
// above it
DataLazy_ptr makePromote(DataLazy_ptr p);

}

#endif // __ESCRIPT_DATALAZY_H__

