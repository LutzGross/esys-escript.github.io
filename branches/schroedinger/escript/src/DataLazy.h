
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
//#include "DataTypes.h"
//#include "FunctionSpace.h"

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
The values of DataPoints are computed when requested rather than all at once.

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
  ESCRIPT_DLL_API
  DataLazy(DataAbstract_ptr p);

  ESCRIPT_DLL_API
  DataLazy(DataAbstract_ptr left, ES_optype op);


/*
  ESCRIPT_DLL_API
  DataLazy(DataLazy_ptr left, DataLazy_ptr right, ES_optype op);*/

  ESCRIPT_DLL_API
  DataLazy(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op);

  ESCRIPT_DLL_API
  ~DataLazy();



  /**
  \brief Compute all data points in the expression tree
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


  ESCRIPT_DLL_API
  int
  getBuffsRequired() const;



private:
  DataReady_ptr m_id;
  DataLazy_ptr m_left, m_right;
  ES_optype m_op;
  size_t m_length;	// number of values represented by the operation

  int m_buffsRequired;	// how many buffers are required to evaluate this expression
  size_t m_samplesize;	// number of values required to store a sample

  const double*
  resolveSample(ValueType& v,int sampleNo,  size_t offset );

  const double*
  resolveSample2(ValueType& v,int sampleNo,  size_t offset );


  void
  intoString(std::ostringstream& oss) const;

  char m_readytype;

  void
  collapse();		// converts the node into an IDENTITY node

  DataReady_ptr
  collapseToReady();

  const double*
  resolveUnary(ValueType& v,int sampleNo,  size_t offset) const;

  const double*
  resolveBinary(ValueType& v,int sampleNo,  size_t offset) const;


};

}
#endif
