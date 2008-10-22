
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
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

namespace escript {

enum ES_optype
{
	UNKNOWNOP=0,
	IDENTITY=1,
	ADD=2,
	SUB=3,
	MUL=4,
	DIV=5,
	SIN=6,
	COS=7,
	TAN=8,
	ASIN=9,
	ACOS=10,
	ATAN=11,
	SINH=12,
	COSH=13,
	TANH=14,
	ERF=15,
	ASINH=16,
	ACOSH=17,
	ATANH=18,
	LOG10=19,
	LOG=20,
	SIGN=21,
	ABS=22,
	NEG=23,
	POS=24,
	EXP=25,
	SQRT=26,
	RECIP=27,
	GZ=28,
	LZ=29,
	GEZ=30,
	LEZ=31
};

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
  \brief Produce a DataLazy for a binary operation.
  \param left left operand
  \param right right operand
  \param op unary operation to perform.
  \throws DataException if op is not a binary operation or if left or right cannot be converted to a DataLazy.
  */
  ESCRIPT_DLL_API
  DataLazy(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op);

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
     Return the number of doubles that would be stored for this Data object if it were resolved.
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


  /**
    \return the number of samples which need to be stored to evaluate the expression.
  */
  ESCRIPT_DLL_API
  int
  getBuffsRequired() const;

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
  size_t m_length;	// number of values represented by the operation

  int m_buffsRequired;	// how many samples are required to evaluate this expression
  size_t m_samplesize;	// number of values required to store a sample

  char m_readytype;	// E for expanded, T for tagged, C for constant


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
  \brief Compute the value of the expression for the given sample.
  \return Vector which stores the value of the subexpression for the given sample.
  \param v A vector to store intermediate results.
  \param offset Index in v to begin storing results.
  \param sampleNo Sample number to evaluate.
  \param roffset (output parameter) the offset in the return vector where the result begins.

  The return value will be an existing vector so do not deallocate it.
  */
  const ValueType*
  resolveSample(ValueType& v,  size_t offset ,int sampleNo, size_t& roffset);

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
  resolveUnary(ValueType& v,  size_t offset,int sampleNo,  size_t& roffset) const;

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

};

}
#endif
