
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


#include "DataLazy.h"
#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif
#ifdef PASO_MPI
#include <mpi.h>
#endif
#include "FunctionSpace.h"
#include "DataTypes.h"
#include "Data.h"

using namespace std;
using namespace boost;

namespace escript
{

const std::string&
opToString(ES_optype op);

namespace
{

// return the FunctionSpace of the result of "left op right"
FunctionSpace
resultFS(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op)
{
	// perhaps this should call interpolate and throw or something?
	// maybe we need an interpolate node -
	// that way, if interpolate is required in any other op we can just throw a 
	// programming error exception.


	if (left->getFunctionSpace()!=right->getFunctionSpace())
	{
		throw DataException("FunctionSpaces not equal - interpolation not supported on lazy data.");
	}
	return left->getFunctionSpace();
}

// return the shape of the result of "left op right"
DataTypes::ShapeType
resultShape(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op)
{
	if (left->getShape()!=right->getShape())
	{
		throw DataException("Shapes not the same - shapes must match for lazy data.");
	}
	return left->getShape();
}

size_t
resultLength(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op)
{
   switch(op)
   {
//   case IDENTITY: return left->getLength();
   case ADD:	// the length is preserved in these ops
   case SUB:
   case MUL:
   case DIV: return left->getLength();
   default: 
	throw DataException("Programmer Error - attempt to getLength() for operator "+opToString(op)+".");

   }
}

int
calcBuffs(const DataLazy_ptr& left, const DataLazy_ptr& right, ES_optype op)
{
   switch(op)
   {
   case IDENTITY: return 0;
   case ADD:	// the length is preserved in these ops
   case SUB:
   case MUL:
   case DIV: return max(left->getBuffsRequired(),right->getBuffsRequired());
   default: 
	throw DataException("Programmer Error - attempt to calcBuffs() for operator "+opToString(op)+".");
   }
}

string ES_opstrings[]={"UNKNOWN","IDENTITY","+","-","*","/"};
int ES_opcount=5;

// void performOp(ValueType& v, int startOffset, ES_optype op, int m_samplesize)
// {
//    switch(op)
//    {
//    case ADD:  DataMaths::binaryOp(v,getShape(),startOffset,v,getShape(),
// 		startOffset+m_samplesize,plus<double>());
// 	      break;	
//    case SUB:  DataMaths::binaryOp(v,getShape(),startOffset,v,getShape(),
// 		startOffset+m_samplesize,minus<double>());
// 	      break;
//    case MUL:  DataMaths::binaryOp(v,getShape(),startOffset,v,getShape(),
// 		startOffset+m_samplesize,multiplies<double>());
// 	      break;
//    case DIV:  DataMaths::binaryOp(v,getShape(),startOffset,v,getShape(),
// 		startOffset+m_samplesize,divides<double>());
// 	      break;
//    default: 
// 	throw DataException("Programmer Error - attempt to performOp() for operator "+opToString(op)+".");
//    }
// 
// }

}	// end anonymous namespace


const std::string&
opToString(ES_optype op)
{
  if (op<0 || op>=ES_opcount) 
  {
    op=UNKNOWNOP;
  }
  return ES_opstrings[op];
}


DataLazy::DataLazy(DataAbstract_ptr p)
	: parent(p->getFunctionSpace(),p->getShape()),
	m_op(IDENTITY)
{
   if (p->isLazy())
   {
	// TODO: fix this.   We could make the new node a copy of p?
	// I don't want identity of Lazy.
	// Question: Why would that be so bad?
	// Answer: We assume that the child of ID is something we can call getVector on
	throw DataException("Programmer error - attempt to create identity from a DataLazy.");
   }
   else
   {
	m_id=dynamic_pointer_cast<DataReady>(p);
   }
   m_length=p->getLength();
   m_buffsRequired=0;
   m_samplesize=getNumDPPSample()*getNoValues();
cout << "(1)Lazy created with " << m_samplesize << endl;
}

DataLazy::DataLazy(DataLazy_ptr left, DataLazy_ptr right, ES_optype op)
	: parent(resultFS(left,right,op), resultShape(left,right,op)),
	m_left(left),
	m_right(right),
	m_op(op)
{
   m_length=resultLength(m_left,m_right,m_op);
   m_samplesize=getNumDPPSample()*getNoValues();
   m_buffsRequired=calcBuffs(m_left, m_right, m_op);
cout << "(2)Lazy created with " << m_samplesize << endl;
}

DataLazy::DataLazy(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op)
	: parent(resultFS(left,right,op), resultShape(left,right,op)),
	m_op(op)
{
   if (left->isLazy())
   {
	m_left=dynamic_pointer_cast<DataLazy>(left);
   }
   else
   {
	m_left=DataLazy_ptr(new DataLazy(left));
   }
   if (right->isLazy())
   {
	m_right=dynamic_pointer_cast<DataLazy>(right);
   }
   else
   {
	m_right=DataLazy_ptr(new DataLazy(right));
   }

   m_length=resultLength(m_left,m_right,m_op);
   m_samplesize=getNumDPPSample()*getNoValues();
   m_buffsRequired=calcBuffs(m_left, m_right,m_op);
cout << "(3)Lazy created with " << m_samplesize << endl;
}


DataLazy::~DataLazy()
{
}


int
DataLazy::getBuffsRequired() const
{
	return m_buffsRequired;
}


// the vector and the offset are a place where the method could write its data if it wishes
// it is not obligated to do so. For example, if it has its own storage already, it can use that.
// hence the return value to indicate where the data is actually stored.
// regardless, the storage should be assumed to be used, even if it isn't.
const double*
DataLazy::resolveSample(ValueType& v,int sampleNo,  size_t offset ) const
{
  if (m_op==IDENTITY)	// copy the contents into the vector
  {
cout << "Begin ID" << endl;
cout << "dpps=" << getNumDPPSample() << " novals=" << getNoValues() << endl;
    const ValueType& vec=m_id->getVector();
//     size_t srcOffset=m_id->getPointOffset(sampleNo, 0);
// cout << "v.size()=" << v.size() << " vec=" << vec.size() << endl;
//     for (size_t i=0;i<m_samplesize;++i,++srcOffset,++offset)
//     {
// cout << "Trying offset=" << offset << " srcOffset=" << srcOffset << endl;
// 	v[offset]=vec[srcOffset];	
//     }
cout << "End ID - returning offset " << m_id->getPointOffset(sampleNo, 0) << " of vector@" << &vec<<endl;
    return &(vec[m_id->getPointOffset(sampleNo, 0)]);
//     return;
  }
cout << "Begin op";
  size_t rightoffset=offset+m_samplesize;
  const double* left=m_left->resolveSample(v,sampleNo,offset);
  const double* right=m_right->resolveSample(v,sampleNo,rightoffset);
  double* result=&(v[offset]);
cout << "left=" << left << " right=" << right << " result=" << result << endl;
//  for (int i=0;i<getNumDPPSample();++i)
  {
    switch(m_op)
    {
    case ADD:		// since these are pointwise ops, pretend each sample is one point
	tensor_binary_operation(m_samplesize, left, right, result, plus<double>());
	break;
    case SUB:		
	tensor_binary_operation(m_samplesize, left, right, result, minus<double>());
	break;
    case MUL:		
	tensor_binary_operation(m_samplesize, left, right, result, multiplies<double>());
	break;
    case DIV:		
	tensor_binary_operation(m_samplesize, left, right, result, divides<double>());
	break;
    default:
	throw DataException("Programmer error - do not know how to resolve operator "+opToString(m_op)+".");
    }
  }
  return result;
}

DataReady_ptr
DataLazy::resolve()
{
  // This is broken!     We need to have a buffer per thread!
  // so the allocation of v needs to move inside the loop somehow

cout << "Sample size=" << m_samplesize << endl;
cout << "Buffers=" << m_buffsRequired << endl;

  ValueType v(m_samplesize*(max(1,m_buffsRequired)+1));	// the +1 comes from the fact that I want to have a safe
							// space for the RHS of ops to write to even if they don't
							// need it.
cout << "Buffer created with size=" << v.size() << endl;
  ValueType dummy(getNoValues());
  DataExpanded* result=new DataExpanded(getFunctionSpace(),getShape(),dummy);
  ValueType& resvec=result->getVector();
  DataReady_ptr resptr=DataReady_ptr(result);
  int sample;
  #pragma omp parallel for private(sample) schedule(static)
  for (sample=0;sample<getNumSamples();++sample)
  {
cout << "Processing sample#" << sample << endl;
    resolveSample(v,sample,0);
cout << "Copying#" << sample << endl;
    for (int i=0;i<m_samplesize;++i)	// copy values into the output vector
    {
	resvec[i]=v[i];
    }
  }
  return resptr;
}

std::string
DataLazy::toString() const
{
  return "Lazy evaluation object. No details available.";
}

// Note that in this case, deepCopy does not make copies of the leaves.
// Hopefully copy on write (or whatever we end up using) will take care of this.
DataAbstract* 
DataLazy::deepCopy()
{
  if (m_op==IDENTITY)
  {
	return new DataLazy(m_left);	// we don't need to copy the child here
  }
  return new DataLazy(m_left->deepCopy()->getPtr(),m_right->deepCopy()->getPtr(),m_op); 
}


DataTypes::ValueType::size_type
DataLazy::getLength() const
{
  return m_length;
}


DataAbstract*
DataLazy::getSlice(const DataTypes::RegionType& region) const
{
  throw DataException("getSlice - not implemented for Lazy objects.");
}

DataTypes::ValueType::size_type 
DataLazy::getPointOffset(int sampleNo,
                 int dataPointNo) const
{
  throw DataException("getPointOffset - not implemented for Lazy objects - yet.");
}

}	// end namespace
