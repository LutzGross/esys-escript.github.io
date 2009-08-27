
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#if !defined escript_DataLazy_20081008_H
#define escript_DataLazy_20081008_H
#include "system_dep.h"

#include "DataAbstract.h"

#include <string>
#include <functional>

#include "LocalOps.h"		// for tensor_binary_op
#include "BufferGroup.h"
#include "DataVector.h"		// for ElementType


#define LAZY_NODE_STORAGE

namespace escript {

// For the purposes of unit testing and maintaining sanity, it is important that this enum be contiguous
enum ES_optype
{
	UNKNOWNOP=0,
	IDENTITY=1,
	ADD=2,
	SUB=3,
	MUL=4,
	DIV=5,
	POW=6,
	SIN=POW+1,
	COS=SIN+1,
	TAN=SIN+2,
	ASIN=SIN+3,
	ACOS=SIN+4,
	ATAN=SIN+5,
	SINH=SIN+6,
	COSH=SIN+7,
	TANH=SIN+8,
	ERF=SIN+9,
	ASINH=SIN+10,
	ACOSH=SIN+11,
	ATANH=SIN+12,
	LOG10=ATANH+1,
	LOG=LOG10+1,
	SIGN=LOG10+2,
	ABS=LOG10+3,
	NEG=LOG10+4,
	POS=LOG10+5,
	EXP=LOG10+6,
	SQRT=LOG10+7,
	RECIP=LOG10+8,
	GZ=RECIP+1,
	LZ=GZ+1,
	GEZ=GZ+2,
	LEZ=GZ+3,
	NEZ=GZ+4,
	EZ=GZ+5,
	SYM=EZ+1,
	NSYM=SYM+1,
	PROD=NSYM+1,
	TRANS=PROD+1,
	TRACE=TRANS+1,
	SWAP=TRACE+1
};

ESCRIPT_DLL_API
const std::string&
opToString(ES_optype op);

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
typedef DataTypes::ValueType ValueType;
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
  deepCopy();


  /**
     \brief
     This method throws an exception. It does not really make sense to ask this question of lazy data.
  */
  ESCRIPT_DLL_API
  ValueType::size_type
  getLength() const;


  ESCRIPT_DLL_API
  DataAbstract*
  getSlice(const DataTypes::RegionType& region) const;


  DataTypes::ValueType::size_type 
  getPointOffset(int sampleNo,
                 int dataPointNo) const;

  DataTypes::ValueType::size_type 
  getPointOffset(int sampleNo,
                 int dataPointNo);


  /**
    \return the number of samples which need to be stored to evaluate the expression.
  */
  ESCRIPT_DLL_API
  int
  getBuffsRequired() const;

  /**
    \return the largest samplesize required to evaluate the expression.
  */
  ESCRIPT_DLL_API
  size_t
  getMaxSampleSize() const;

  /**
    \return the size of the buffer required to evaulate a sample for this object
  */
  ESCRIPT_DLL_API
  size_t
  getSampleBufferSize() const;

   /**
  \brief Compute the value of the expression for the given sample.
  \return Vector which stores the value of the subexpression for the given sample.
  \param bg A BufferGroup to store intermediate results.
  \param sampleNo Sample number to evaluate.
  \param roffset (output parameter) the offset in the return vector where the result begins.

  The return value will be an existing vector so do not deallocate it.
  */
  ESCRIPT_DLL_API
  const ValueType*
  resolveSample(BufferGroup& bg, int sampleNo, size_t& roffset); 

  /**
  \brief if resolve() was called would it produce expanded data.
  */
  ESCRIPT_DLL_API
  bool
  actsExpanded() const;

  /**
  \brief if resolve() was called would it produce constant data.
  */
  ESCRIPT_DLL_API
  bool
  actsConstant() const;

  /**
     \brief Produces an IDENTITY DataLazy containing zero.
     The result will have the same shape and functionspace as before.
  */
  ESCRIPT_DLL_API
  virtual void
  setToZero();

private:
  DataReady_ptr m_id;	//  For IDENTITY nodes, stores a wrapped value.
  DataLazy_ptr m_left, m_right;	// operands for operation.
  ES_optype m_op;	// operation to perform.

  int m_buffsRequired;	// how many samples are required to evaluate this expression
  size_t m_samplesize;	// number of values required to store a sample

  char m_readytype;	// E for expanded, T for tagged, C for constant

  int m_axis_offset;	// required extra info for general tensor product
  int m_transpose;	// offset and transpose are used for swapaxes as well
  int m_SL, m_SM, m_SR;	// computed properties used in general tensor product


  double m_tol;		// required extra info for <>0 and ==0

  size_t m_maxsamplesize;	// largest samplesize required by any node in the expression
  size_t m_children;
  size_t m_height;

#ifdef LAZY_NODE_STORAGE

  int* m_sampleids;		// may be NULL
  DataVector m_samples;  

#endif // LAZY_NODE_STORAGE


#ifdef LAZY_NODE_STORAGE
  /**
  Allocates sample storage at each node
  */
  void LazyNodeSetup();


  const DataTypes::ValueType*
  resolveNodeUnary(int tid, int sampleNo, size_t& roffset);

  const DataTypes::ValueType*
  resolveNodeSample(int tid, int sampleNo, size_t& roffset);

  const DataTypes::ValueType*
  resolveNodeBinary(int tid, int sampleNo, size_t& roffset);

  const DataTypes::ValueType*
  resolveNodeNP1OUT(int tid, int sampleNo, size_t& roffset);

  const DataTypes::ValueType*
  resolveNodeNP1OUT_P(int tid, int sampleNo, size_t& roffset);

  const DataTypes::ValueType*
  resolveNodeTProd(int tid, int sampleNo, size_t& roffset);

  const DataTypes::ValueType*
  resolveNodeNP1OUT_2P(int tid, int sampleNo, size_t& roffset);
  
#endif

  /**
  Does the work for toString. 
  */
  void
  intoString(std::ostringstream& oss) const;

  /**
   \brief Converts the DataLazy into an IDENTITY storing the value of the expression.
   This method uses the original methods on the Data class to evaluate the expressions.
   For this reason, it should not be used on DataExpanded instances. (To do so would defeat
   the purpose of using DataLazy in the first place).
  */
  void
  collapse();		// converts the node into an IDENTITY node


  /**
  \brief Evaluates the expression using methods on Data.
  This does the work for the collapse method.
  For reasons of efficiency do not call this method on DataExpanded nodes.
  */
  DataReady_ptr
  collapseToReady();

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
  \brief resolve to a ReadyData object using a vector buffer.
  */
  DataReady_ptr
  resolveVectorWorker();

#ifdef LAZY_NODE_STORAGE
  /**
  \brief resolve to a ReadyData object using storage at nodes
  */
  DataReady_ptr
  resolveNodeWorker();
#endif

  /**
  \brief Compute the value of the expression for the given sample - using the vector buffer approach.
  \return Vector which stores the value of the subexpression for the given sample.
  \param v A vector to store intermediate results.
  \param offset Index in v to begin storing results.
  \param sampleNo Sample number to evaluate.
  \param roffset (output parameter) the offset in the return vector where the result begins.

  The return value will be an existing vector so do not deallocate it.
  */
  ESCRIPT_DLL_API
  const ValueType*
  resolveVectorSample(ValueType& v,  size_t offset, int sampleNo, size_t& roffset);


  /**
  \brief Compute the value of the expression (unary operation) for the given sample.
  \return Vector which stores the value of the subexpression for the given sample.
  \param v A vector to store intermediate results.
  \param offset Index in v to begin storing results.
  \param sampleNo Sample number to evaluate.
  \param roffset (output parameter) the offset in the return vector where the result begins.

  The return value will be an existing vector so do not deallocate it.
  If the result is stored in v it should be stored at the offset given.
  Everything from offset to the end of v should be considered available for this method to use.
  */
  ValueType*
  resolveUnary(ValueType& v,  size_t offset,int sampleNo,  size_t& roffset) const;

  /**
  \brief Compute the value of the expression (unary non-pointwise operation) for the given sample.
  \return Vector which stores the value of the subexpression for the given sample.
  \param v A vector to store intermediate results.
  \param offset Index in v to begin storing results.
  \param sampleNo Sample number to evaluate.
  \param roffset (output parameter) the offset in the return vector where the result begins.

  The return value will be an existing vector so do not deallocate it.
  If the result is stored in v it should be stored at the offset given.
  Everything from offset to the end of v should be considered available for this method to use.

  This method differs from the one above in that deals with operations that are not
  point-wise. That is, the answer cannot just be written on top of the input.
  Extra buffers are required for these operations.
  */

  ValueType*
  resolveNP1OUT(ValueType& v, size_t offset, int sampleNo, size_t& roffset) const;

/**
  \brief Compute the value of the expression (unary operation) for the given sample.
  \return Vector which stores the value of the subexpression for the given sample.
  \param v A vector to store intermediate results.
  \param offset Index in v to begin storing results.
  \param sampleNo Sample number to evaluate.
  \param roffset (output parameter) the offset in the return vector where the result begins.

  The return value will be an existing vector so do not deallocate it.
  If the result is stored in v it should be stored at the offset given.
  Everything from offset to the end of v should be considered available for this method to use.
*/
DataTypes::ValueType*
resolveNP1OUT_P(ValueType& v, size_t offset, int sampleNo, size_t& roffset) const;

/**
  \brief Compute the value of the expression (unary operation with int params) for the given sample.
  \return Vector which stores the value of the subexpression for the given sample.
  \param v A vector to store intermediate results.
  \param offset Index in v to begin storing results.
  \param sampleNo Sample number to evaluate.
  \param roffset (output parameter) the offset in the return vector where the result begins.

  The return value will be an existing vector so do not deallocate it.
  If the result is stored in v it should be stored at the offset given.
  Everything from offset to the end of v should be considered available for this method to use.
*/
DataTypes::ValueType*
resolveNP1OUT_2P(ValueType& v, size_t offset, int sampleNo, size_t& roffset) const;


  /**
  \brief Compute the value of the expression (binary operation) for the given sample.
  \return Vector which stores the value of the subexpression for the given sample.
  \param v A vector to store intermediate results.
  \param offset Index in v to begin storing results.
  \param sampleNo Sample number to evaluate.
  \param roffset (output parameter) the offset in the return vector where the result begins.

  The return value will be an existing vector so do not deallocate it.
  If the result is stored in v it should be stored at the offset given.
  Everything from offset to the end of v should be considered available for this method to use.
  */
  ValueType*
  resolveBinary(ValueType& v,  size_t offset,int sampleNo,  size_t& roffset) const;

  /**
  \brief Compute the value of the expression (tensor product) for the given sample.
  \return Vector which stores the value of the subexpression for the given sample.
  \param v A vector to store intermediate results.
  \param offset Index in v to begin storing results.
  \param sampleNo Sample number to evaluate.
  \param roffset (output parameter) the offset in the return vector where the result begins.

  The return value will be an existing vector so do not deallocate it.
  If the result is stored in v it should be stored at the offset given.
  Everything from offset to the end of v should be considered available for this method to use.
  */
  DataTypes::ValueType*
  resolveTProd(ValueType& v,  size_t offset, int sampleNo, size_t& roffset) const;

};

}
#endif
