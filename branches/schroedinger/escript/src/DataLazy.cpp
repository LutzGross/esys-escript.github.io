
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
#ifdef _OPENMP
#include <omp.h>
#endif
#include "FunctionSpace.h"
#include "DataTypes.h"
#include "Data.h"
#include "UnaryFuncs.h"		// for escript::fsign

using namespace std;
using namespace boost;

namespace escript
{

const std::string&
opToString(ES_optype op);

namespace
{



enum ES_opgroup
{
   G_UNKNOWN,
   G_IDENTITY,
   G_BINARY,
   G_UNARY
};




string ES_opstrings[]={"UNKNOWN","IDENTITY","+","-","*","/","sin","cos","tan",
			"asin","acos","atan","sinh","cosh","tanh","erf",
			"asinh","acosh","atanh",
			"log10","log","sign","abs","neg","pos"};
int ES_opcount=25;
ES_opgroup opgroups[]={G_UNKNOWN,G_IDENTITY,G_BINARY,G_BINARY,G_BINARY,G_BINARY,G_UNARY,G_UNARY,G_UNARY, //9
			G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,	// 16
			G_UNARY,G_UNARY,G_UNARY,					// 19
			G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY};		// 25

inline
ES_opgroup
getOpgroup(ES_optype op)
{
  return opgroups[op];
}

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
   switch (getOpgroup(op))
   {
   case G_BINARY: return left->getLength();
   case G_UNARY: return left->getLength();
   default: 
	throw DataException("Programmer Error - attempt to getLength() for operator "+opToString(op)+".");
   }
}

int
calcBuffs(const DataLazy_ptr& left, const DataLazy_ptr& right, ES_optype op)
{
   switch(getOpgroup(op))
   {
   case G_IDENTITY: return 1;
   case G_BINARY: return max(left->getBuffsRequired(),right->getBuffsRequired()+1);
   case G_UNARY: return max(left->getBuffsRequired(),1);
   default: 
	throw DataException("Programmer Error - attempt to calcBuffs() for operator "+opToString(op)+".");
   }
}



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
   m_buffsRequired=1;
   m_samplesize=getNumDPPSample()*getNoValues();
cout << "(1)Lazy created with " << m_samplesize << endl;
}

DataLazy::DataLazy(DataAbstract_ptr left, ES_optype op)
	: parent(left->getFunctionSpace(),left->getShape()),
	m_op(op)
{
   if (getOpgroup(op)!=G_UNARY)
   {
	throw DataException("Programmer error - constructor DataLazy(left, op) will only process UNARY operations.");
   }
   DataLazy_ptr lleft;
   if (!left->isLazy())
   {
	lleft=DataLazy_ptr(new DataLazy(left));
   }
   else
   {
	lleft=dynamic_pointer_cast<DataLazy>(left);
   }
   m_length=left->getLength();
   m_left=lleft;
   m_buffsRequired=1;
   m_samplesize=getNumDPPSample()*getNoValues();
}


DataLazy::DataLazy(DataLazy_ptr left, DataLazy_ptr right, ES_optype op)
	: parent(resultFS(left,right,op), resultShape(left,right,op)),
	m_left(left),
	m_right(right),
	m_op(op)
{
   if (getOpgroup(op)!=G_BINARY)
   {
	throw DataException("Programmer error - constructor DataLazy(left, right, op) will only process BINARY operations.");
   }
   m_length=resultLength(m_left,m_right,m_op);
   m_samplesize=getNumDPPSample()*getNoValues();
   m_buffsRequired=calcBuffs(m_left, m_right, m_op);
cout << "(2)Lazy created with " << m_samplesize << endl;
}

DataLazy::DataLazy(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op)
	: parent(resultFS(left,right,op), resultShape(left,right,op)),
	m_op(op)
{
   if (getOpgroup(op)!=G_BINARY)
   {
	throw DataException("Programmer error - constructor DataLazy(left, op) will only process BINARY operations.");
   }
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
// Hence the return value to indicate where the data is actually stored.
// Regardless, the storage should be assumed to be used, even if it isn't.
const double*
DataLazy::resolveSample(ValueType& v,int sampleNo,  size_t offset ) const
{
  if (m_op==IDENTITY)	
  {
    const ValueType& vec=m_id->getVector();
    return &(vec[m_id->getPointOffset(sampleNo, 0)]);
  }
  size_t rightoffset=offset+m_samplesize;
  const double* left=m_left->resolveSample(v,sampleNo,offset);
  const double* right=0;
  if (getOpgroup(m_op)==G_BINARY)
  {
    right=m_right->resolveSample(v,sampleNo,rightoffset);
  }
  double* result=&(v[offset]);
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
// unary ops
    case SIN:
	tensor_unary_operation(m_samplesize, left, result, ::sin);
	break;
    case COS:
	tensor_unary_operation(m_samplesize, left, result, ::cos);
	break;
    case TAN:
	tensor_unary_operation(m_samplesize, left, result, ::tan);
	break;
    case ASIN:
	tensor_unary_operation(m_samplesize, left, result, ::asin);
	break;
    case ACOS:
	tensor_unary_operation(m_samplesize, left, result, ::acos);
	break;
    case ATAN:
	tensor_unary_operation(m_samplesize, left, result, ::atan);
	break;
    case SINH:
	tensor_unary_operation(m_samplesize, left, result, ::sinh);
	break;
    case COSH:
	tensor_unary_operation(m_samplesize, left, result, ::cosh);
	break;
    case TANH:
	tensor_unary_operation(m_samplesize, left, result, ::tanh);
	break;
    case ERF:
#ifdef _WIN32
	throw DataException("Error - Data:: erf function is not supported on _WIN32 platforms.");
#else
	tensor_unary_operation(m_samplesize, left, result, ::erf);
	break;
#endif
   case ASINH:
#ifdef _WIN32
	tensor_unary_operation(m_samplesize, left, result, escript::asinh_substitute);
#else
	tensor_unary_operation(m_samplesize, left, result, ::asinh);
#endif   
	break;
   case ACOSH:
#ifdef _WIN32
	tensor_unary_operation(m_samplesize, left, result, escript::acosh_substitute);
#else
	tensor_unary_operation(m_samplesize, left, result, ::acosh);
#endif   
	break;
   case ATANH:
#ifdef _WIN32
	tensor_unary_operation(m_samplesize, left, result, escript::atanh_substitute);
#else
	tensor_unary_operation(m_samplesize, left, result, ::atanh);
#endif   
	break;
    case LOG10:
	tensor_unary_operation(m_samplesize, left, result, ::log10);
	break;
    case LOG:
	tensor_unary_operation(m_samplesize, left, result, ::log);
	break;
    case SIGN:
	tensor_unary_operation(m_samplesize, left, result, escript::fsign);
	break;
    case ABS:
	tensor_unary_operation(m_samplesize, left, result, ::fabs);
	break;
    case NEG:
	tensor_unary_operation(m_samplesize, left, result, negate<double>());
	break;
    case POS:
	// it doesn't mean anything for delayed.
	// it will just trigger a deep copy of the lazy object
	throw DataException("Programmer error - POS not supported for lazy data.");
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

  size_t threadbuffersize=m_samplesize*(max(1,m_buffsRequired)+1);
  int numthreads=1;
#ifdef _OPENMP
  numthreads=omp_get_max_threads();
  int threadnum=0;
#endif 
  ValueType v(numthreads*threadbuffersize);	
cout << "Buffer created with size=" << v.size() << endl;
  ValueType dummy(getNoValues());
  DataExpanded* result=new DataExpanded(getFunctionSpace(),getShape(),dummy);
  ValueType& resvec=result->getVector();
  DataReady_ptr resptr=DataReady_ptr(result);
  int sample;
  int resoffset;
  int totalsamples=getNumSamples();
  #pragma omp parallel for private(sample,resoffset,threadnum) schedule(static)
  for (sample=0;sample<totalsamples;++sample)
  {
#ifdef _OPENMP
    const double* res=resolveSample(v,sample,threadbuffersize*omp_get_thread_num());
#else
    const double* res=resolveSample(v,sample,0);   // this would normally be v, but not if its a single IDENTITY op.
#endif
    resoffset=result->getPointOffset(sample,0);
    for (int i=0;i<m_samplesize;++i,++resoffset)	// copy values into the output vector
    {
	resvec[resoffset]=res[i];
    }
  }
  return resptr;
}

std::string
DataLazy::toString() const
{
  ostringstream oss;
  oss << "Lazy Data:";
  intoString(oss);
  return oss.str();
}

void
DataLazy::intoString(ostringstream& oss) const
{
  switch (getOpgroup(m_op))
  {
  case G_IDENTITY: 
	oss << '@' << m_id.get();
	break;
  case G_BINARY:
	oss << '(';
	m_left->intoString(oss);
	oss << ' ' << opToString(m_op) << ' ';
	m_right->intoString(oss);
	oss << ')';
	break;
  case G_UNARY:
	oss << opToString(m_op) << '(';
	m_left->intoString(oss);
	oss << ')';
	break;
  default:
	oss << "UNKNOWN";
  }
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
