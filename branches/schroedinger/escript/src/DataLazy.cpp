
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
			"log10","log","sign","abs","neg","pos","exp","sqrt",
			"1/","where>0","where<0","where>=0","where<=0"};
int ES_opcount=32;
ES_opgroup opgroups[]={G_UNKNOWN,G_IDENTITY,G_BINARY,G_BINARY,G_BINARY,G_BINARY,G_UNARY,G_UNARY,G_UNARY, //9
			G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,	// 16
			G_UNARY,G_UNARY,G_UNARY,					// 19
			G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,		// 27
			G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY};
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
  	if(p->isConstant()) {m_readytype='C';}
  	else if(p->isExpanded()) {m_readytype='E';}
  	else if (p->isTagged()) {m_readytype='T';}
	else {throw DataException("Unknown DataReady instance in DataLazy constructor.");}
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
   m_readytype=lleft->m_readytype;
   m_length=left->getLength();
   m_left=lleft;
   m_buffsRequired=1;
   m_samplesize=getNumDPPSample()*getNoValues();
}


// DataLazy::DataLazy(DataLazy_ptr left, DataLazy_ptr right, ES_optype op)
// 	: parent(resultFS(left,right,op), resultShape(left,right,op)),
// 	m_left(left),
// 	m_right(right),
// 	m_op(op)
// {
//    if (getOpgroup(op)!=G_BINARY)
//    {
// 	throw DataException("Programmer error - constructor DataLazy(left, right, op) will only process BINARY operations.");
//    }
//    m_length=resultLength(m_left,m_right,m_op);
//    m_samplesize=getNumDPPSample()*getNoValues();
//    m_buffsRequired=calcBuffs(m_left, m_right, m_op);
// cout << "(2)Lazy created with " << m_samplesize << endl;
// }

DataLazy::DataLazy(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op)
	: parent(resultFS(left,right,op), resultShape(left,right,op)),
	m_op(op)
{
   if (getOpgroup(op)!=G_BINARY)
   {
	throw DataException("Programmer error - constructor DataLazy(left, right, op) will only process BINARY operations.");
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
   char lt=m_left->m_readytype;
   char rt=m_right->m_readytype;
   if (lt=='E' || rt=='E')
   {
	m_readytype='E';
   }
   else if (lt=='T' || rt=='T')
   {
	m_readytype='T';
   }
   else
   {
	m_readytype='C';
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


DataReady_ptr
DataLazy::collapseToReady()
{
  if (m_readytype=='E')
  {	// this is more an efficiency concern than anything else
    throw DataException("Programmer Error - do not use collapse on Expanded data.");
  }
  if (m_op==IDENTITY)
  {
    return m_id;
  }
  DataReady_ptr pleft=m_left->collapseToReady();
  Data left(pleft);
  Data right;
  if (getOpgroup(m_op)==G_BINARY)
  {
    right=Data(m_right->collapseToReady());
  }
  Data result;
  switch(m_op)
  {
    case ADD:
	result=left+right;
	break;
    case SUB:		
	result=left-right;
	break;
    case MUL:		
	result=left*right;
	break;
    case DIV:		
	result=left/right;
	break;
    case SIN:
	result=left.sin();	
	break;
    case COS:
	result=left.cos();
	break;
    case TAN:
	result=left.tan();
	break;
    case ASIN:
	result=left.asin();
	break;
    case ACOS:
	result=left.acos();
	break;
    case ATAN:
	result=left.atan();
	break;
    case SINH:
	result=left.sinh();
	break;
    case COSH:
	result=left.cosh();
	break;
    case TANH:
	result=left.tanh();
	break;
    case ERF:
	result=left.erf();
	break;
   case ASINH:
	result=left.asinh();
	break;
   case ACOSH:
	result=left.acosh();
	break;
   case ATANH:
	result=left.atanh();
	break;
    case LOG10:
	result=left.log10();
	break;
    case LOG:
	result=left.log();
	break;
    case SIGN:
	result=left.sign();
	break;
    case ABS:
	result=left.abs();
	break;
    case NEG:
	result=left.neg();
	break;
    case POS:
	// it doesn't mean anything for delayed.
	// it will just trigger a deep copy of the lazy object
	throw DataException("Programmer error - POS not supported for lazy data.");
	break;
    case EXP:
	result=left.exp();
	break;
    case SQRT:
	result=left.sqrt();
	break;
    case RECIP:
	result=left.oneOver();
	break;
    case GZ:
	result=left.wherePositive();
	break;
    case LZ:
	result=left.whereNegative();
	break;
    case GEZ:
	result=left.whereNonNegative();
	break;
    case LEZ:
	result=left.whereNonPositive();
	break;
    default:
	throw DataException("Programmer error - do not know how to resolve operator "+opToString(m_op)+".");
  }
  return result.borrowReadyPtr();
}

void
DataLazy::collapse()
{
  if (m_op==IDENTITY)
  {
	return;
  }
  if (m_readytype=='E')
  {	// this is more an efficiency concern than anything else
    throw DataException("Programmer Error - do not use collapse on Expanded data.");
  }
  m_id=collapseToReady();
  m_op=IDENTITY;
}

DataTypes::ValueType*
DataLazy::resolveUnary(ValueType& v, size_t offset, int sampleNo, size_t& roffset) const
{
	// we assume that any collapsing has been done before we get here
	// since we only have one argument we don't need to think about only
	// processing single points.
  if (m_readytype!='E')
  {
    throw DataException("Programmer error - resolveUnary should only be called on expanded Data.");
  }
  const ValueType* vleft=m_left->resolveSample(v,offset,sampleNo,roffset);
  const double* left=&((*vleft)[roffset]);
  double* result=&(v[offset]);
  roffset=offset;
  switch (m_op)
  {
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
    case EXP:
	tensor_unary_operation(m_samplesize, left, result, ::exp);
	break;
    case SQRT:
	tensor_unary_operation(m_samplesize, left, result, ::sqrt);
	break;
    case RECIP:
	tensor_unary_operation(m_samplesize, left, result, bind1st(divides<double>(),1.));
	break;
    case GZ:
	tensor_unary_operation(m_samplesize, left, result, bind2nd(greater<double>(),0.0));
	break;
    case LZ:
	tensor_unary_operation(m_samplesize, left, result, bind2nd(less<double>(),0.0));
	break;
    case GEZ:
	tensor_unary_operation(m_samplesize, left, result, bind2nd(greater_equal<double>(),0.0));
	break;
    case LEZ:
	tensor_unary_operation(m_samplesize, left, result, bind2nd(less_equal<double>(),0.0));
	break;

    default:
	throw DataException("Programmer error - resolveUnary can not resolve operator "+opToString(m_op)+".");
  }
  return &v;
}



// const double*
// DataLazy::resolveUnary(ValueType& v,int sampleNo,  size_t offset) const
// {
// 	// we assume that any collapsing has been done before we get here
// 	// since we only have one argument we don't need to think about only
// 	// processing single points.
//   if (m_readytype!='E')
//   {
//     throw DataException("Programmer error - resolveUnary should only be called on expanded Data.");
//   }
//   const double* left=m_left->resolveSample(v,sampleNo,offset);
//   double* result=&(v[offset]);
//   switch (m_op)
//   {
//     case SIN:	
// 	tensor_unary_operation(m_samplesize, left, result, ::sin);
// 	break;
//     case COS:
// 	tensor_unary_operation(m_samplesize, left, result, ::cos);
// 	break;
//     case TAN:
// 	tensor_unary_operation(m_samplesize, left, result, ::tan);
// 	break;
//     case ASIN:
// 	tensor_unary_operation(m_samplesize, left, result, ::asin);
// 	break;
//     case ACOS:
// 	tensor_unary_operation(m_samplesize, left, result, ::acos);
// 	break;
//     case ATAN:
// 	tensor_unary_operation(m_samplesize, left, result, ::atan);
// 	break;
//     case SINH:
// 	tensor_unary_operation(m_samplesize, left, result, ::sinh);
// 	break;
//     case COSH:
// 	tensor_unary_operation(m_samplesize, left, result, ::cosh);
// 	break;
//     case TANH:
// 	tensor_unary_operation(m_samplesize, left, result, ::tanh);
// 	break;
//     case ERF:
// #ifdef _WIN32
// 	throw DataException("Error - Data:: erf function is not supported on _WIN32 platforms.");
// #else
// 	tensor_unary_operation(m_samplesize, left, result, ::erf);
// 	break;
// #endif
//    case ASINH:
// #ifdef _WIN32
// 	tensor_unary_operation(m_samplesize, left, result, escript::asinh_substitute);
// #else
// 	tensor_unary_operation(m_samplesize, left, result, ::asinh);
// #endif   
// 	break;
//    case ACOSH:
// #ifdef _WIN32
// 	tensor_unary_operation(m_samplesize, left, result, escript::acosh_substitute);
// #else
// 	tensor_unary_operation(m_samplesize, left, result, ::acosh);
// #endif   
// 	break;
//    case ATANH:
// #ifdef _WIN32
// 	tensor_unary_operation(m_samplesize, left, result, escript::atanh_substitute);
// #else
// 	tensor_unary_operation(m_samplesize, left, result, ::atanh);
// #endif   
// 	break;
//     case LOG10:
// 	tensor_unary_operation(m_samplesize, left, result, ::log10);
// 	break;
//     case LOG:
// 	tensor_unary_operation(m_samplesize, left, result, ::log);
// 	break;
//     case SIGN:
// 	tensor_unary_operation(m_samplesize, left, result, escript::fsign);
// 	break;
//     case ABS:
// 	tensor_unary_operation(m_samplesize, left, result, ::fabs);
// 	break;
//     case NEG:
// 	tensor_unary_operation(m_samplesize, left, result, negate<double>());
// 	break;
//     case POS:
// 	// it doesn't mean anything for delayed.
// 	// it will just trigger a deep copy of the lazy object
// 	throw DataException("Programmer error - POS not supported for lazy data.");
// 	break;
//     case EXP:
// 	tensor_unary_operation(m_samplesize, left, result, ::exp);
// 	break;
//     case SQRT:
// 	tensor_unary_operation(m_samplesize, left, result, ::sqrt);
// 	break;
//     case RECIP:
// 	tensor_unary_operation(m_samplesize, left, result, bind1st(divides<double>(),1.));
// 	break;
//     case GZ:
// 	tensor_unary_operation(m_samplesize, left, result, bind2nd(greater<double>(),0.0));
// 	break;
//     case LZ:
// 	tensor_unary_operation(m_samplesize, left, result, bind2nd(less<double>(),0.0));
// 	break;
//     case GEZ:
// 	tensor_unary_operation(m_samplesize, left, result, bind2nd(greater_equal<double>(),0.0));
// 	break;
//     case LEZ:
// 	tensor_unary_operation(m_samplesize, left, result, bind2nd(less_equal<double>(),0.0));
// 	break;
// 
//     default:
// 	throw DataException("Programmer error - resolveUnary can not resolve operator "+opToString(m_op)+".");
//   }
//   return result;
// }

#define PROC_OP(X) \
	for (int i=0;i<steps;++i,resultp+=getNoValues()) \
	{ \
cout << "Step#" << i << " chunk=" << chunksize << endl; \
cout << left[0] << left[1] << left[2] << endl; \
cout << right[0] << right[1] << right[2] << endl; \
	   tensor_binary_operation(chunksize, left, right, resultp, X); \
	   left+=leftStep; \
	   right+=rightStep; \
cout << "Result=" << result << " " << result[0] << result[1] << result[2] << endl; \
	}

DataTypes::ValueType*
DataLazy::resolveBinary(ValueType& v,  size_t offset ,int sampleNo, size_t& roffset) const
{
	// again we assume that all collapsing has already been done
	// so we have at least one expanded child.
	// however, we could still have one of the children being not expanded.

cout << "Resolve binary: " << toString() << endl;

  size_t lroffset=0, rroffset=0;

  bool leftExp=(m_left->m_readytype=='E');
  bool rightExp=(m_right->m_readytype=='E');
  bool bigloops=((leftExp && rightExp) || (!leftExp && !rightExp));	// is processing in single step
  int steps=(bigloops?1:getNumDPPSample());
  size_t chunksize=(bigloops? m_samplesize : getNoValues());
  int leftStep=((leftExp && !rightExp)? getNoValues() : 0);
  int rightStep=((rightExp && !leftExp)? getNoValues() : 0);

  const ValueType* left=m_left->resolveSample(v,offset,sampleNo,lroffset);
  const ValueType* right=m_right->resolveSample(v,offset,sampleNo,rroffset);	
  	// now we need to know which args are expanded
cout << "left=" << left << " right=" << right << endl;
cout << "(Length) l=" << left->size() << " r=" << right->size() << " res=" << v.size() << endl;
  double* resultp=&(v[offset]);
  switch(m_op)
  {
    case ADD:
	for (int i=0;i<steps;++i,resultp+=getNoValues()) 
	{ 
cerr << "Step#" << i << " chunk=" << chunksize << endl; 
cerr << left << "[" << lroffset << "] " << right << "[" << rroffset << "]" << endl;
	   tensor_binary_operation(chunksize, &((*left)[lroffset]), &((*right)[rroffset]), resultp, plus<double>()); 
	   lroffset+=leftStep; 
	   rroffset+=rightStep; 
cerr << "left=" << lroffset << " right=" << rroffset << endl;
	}
	break;
// need to fill in the rest
    default:
	throw DataException("Programmer error - resolveBinary can not resolve operator "+opToString(m_op)+".");
  }
  roffset=offset;
  return &v;
}



// #define PROC_OP(X) \
// 	for (int i=0;i<steps;++i,resultp+=getNoValues()) \
// 	{ \
// cout << "Step#" << i << " chunk=" << chunksize << endl; \
// cout << left[0] << left[1] << left[2] << endl; \
// cout << right[0] << right[1] << right[2] << endl; \
// 	   tensor_binary_operation(chunksize, left, right, resultp, X); \
// 	   left+=leftStep; \
// 	   right+=rightStep; \
// cout << "Result=" << result << " " << result[0] << result[1] << result[2] << endl; \
// 	}
// 
// const double*
// DataLazy::resolveBinary(ValueType& v,int sampleNo,  size_t offset) const
// {
// 	// again we assume that all collapsing has already been done
// 	// so we have at least one expanded child.
// 	// however, we could still have one of the children being not expanded.
// 
// cout << "Resolve binary: " << toString() << endl;
// 
//   const double* left=m_left->resolveSample(v,sampleNo,offset);
// // cout << "Done Left " << /*left[0] << left[1] << left[2] << */endl;
//   const double* right=m_right->resolveSample(v,sampleNo,offset);	
// // cout << "Done Right"  << /*right[0] << right[1] << right[2] <<*/ endl;
//   	// now we need to know which args are expanded
//   bool leftExp=(m_left->m_readytype=='E');
//   bool rightExp=(m_right->m_readytype=='E');
//   bool bigloops=((leftExp && rightExp) || (!leftExp && !rightExp));	// is processing in single step
//   int steps=(bigloops?1:getNumSamples());
//   size_t chunksize=(bigloops? m_samplesize : getNoValues());
//   int leftStep=((leftExp && !rightExp)? getNoValues() : 0);
//   int rightStep=((rightExp && !leftExp)? getNoValues() : 0);
// cout << "left=" << left << " right=" << right << endl;
//   double* result=&(v[offset]);
//   double* resultp=result;
//   switch(m_op)
//   {
//     case ADD:
// 	for (int i=0;i<steps;++i,resultp+=getNoValues()) 
// 	{ 
// cout << "Step#" << i << " chunk=" << chunksize << endl; \
// // cout << left[0] << left[1] << left[2] << endl; 
// // cout << right[0] << right[1] << right[2] << endl; 
// 	   tensor_binary_operation(chunksize, left, right, resultp, plus<double>()); 
// cout << "left=" << left << " right=" << right << " resp=" << resultp << endl;
// 	   left+=leftStep; 
// 	   right+=rightStep; 
// cout << "left=" << left << " right=" << right << endl;
// // cout << "Result=" << result << " " << result[0] << result[1] << result[2] << endl; 
// 	}
// 	break;
// // need to fill in the rest
//     default:
// 	throw DataException("Programmer error - resolveBinay can not resolve operator "+opToString(m_op)+".");
//   }
// // cout << "About to return "  << result[0] << result[1] << result[2] << endl;;
//   return result;
// }

// // the vector and the offset are a place where the method could write its data if it wishes
// // it is not obligated to do so. For example, if it has its own storage already, it can use that.
// // Hence the return value to indicate where the data is actually stored.
// // Regardless, the storage should be assumed to be used, even if it isn't.
// const double*
// DataLazy::resolveSample(ValueType& v,int sampleNo,  size_t offset )
// {
// cout << "Resolve sample " << toString() << endl;
// 	// collapse so we have a 'E' node or an IDENTITY for some other type
//   if (m_readytype!='E' && m_op!=IDENTITY)
//   {
// 	collapse();
//   }
//   if (m_op==IDENTITY)	
//   {
//     const ValueType& vec=m_id->getVector();
//     if (m_readytype=='C')
//     {
// 	return &(vec[0]);
//     }
//     return &(vec[m_id->getPointOffset(sampleNo, 0)]);
//   }
//   if (m_readytype!='E')
//   {
//     throw DataException("Programmer Error - Collapse did not produce an expanded node.");
//   }
//   switch (getOpgroup(m_op))
//   {
//   case G_UNARY: return resolveUnary(v,sampleNo,offset);
//   case G_BINARY: return resolveBinary(v,sampleNo,offset);
//   default:
//     throw DataException("Programmer Error - resolveSample does not know how to process "+opToString(m_op)+".");
//   }
// }



// the vector and the offset are a place where the method could write its data if it wishes
// it is not obligated to do so. For example, if it has its own storage already, it can use that.
// Hence the return value to indicate where the data is actually stored.
// Regardless, the storage should be assumed to be used, even if it isn't.

// the roffset is the offset within the returned vector where the data begins
const DataTypes::ValueType*
DataLazy::resolveSample(ValueType& v, size_t offset, int sampleNo, size_t& roffset)
{
cout << "Resolve sample " << toString() << endl;
	// collapse so we have a 'E' node or an IDENTITY for some other type
  if (m_readytype!='E' && m_op!=IDENTITY)
  {
	collapse();
  }
  if (m_op==IDENTITY)	
  {
    const ValueType& vec=m_id->getVector();
    if (m_readytype=='C')
    {
	roffset=0;
	return &(vec);
    }
    roffset=m_id->getPointOffset(sampleNo, 0);
    return &(vec);
  }
  if (m_readytype!='E')
  {
    throw DataException("Programmer Error - Collapse did not produce an expanded node.");
  }
  switch (getOpgroup(m_op))
  {
  case G_UNARY: return resolveUnary(v, offset,sampleNo,roffset);
  case G_BINARY: return resolveBinary(v, offset,sampleNo,roffset);
  default:
    throw DataException("Programmer Error - resolveSample does not know how to process "+opToString(m_op)+".");
  }
}




// This version uses double* trying again with vectors
// DataReady_ptr
// DataLazy::resolve()
// {
// 
// cout << "Sample size=" << m_samplesize << endl;
// cout << "Buffers=" << m_buffsRequired << endl;
// 
//   if (m_readytype!='E')
//   {
//     collapse();
//   }
//   if (m_op==IDENTITY)
//   {
//     return m_id;
//   }
//   	// from this point on we must have m_op!=IDENTITY and m_readytype=='E'
//   size_t threadbuffersize=m_samplesize*(max(1,m_buffsRequired)+1);
//   int numthreads=1;
// #ifdef _OPENMP
//   numthreads=getNumberOfThreads();
//   int threadnum=0;
// #endif 
//   ValueType v(numthreads*threadbuffersize);	
// cout << "Buffer created with size=" << v.size() << endl;
//   DataExpanded* result=new DataExpanded(getFunctionSpace(),getShape(),  ValueType(getNoValues()));
//   ValueType& resvec=result->getVector();
//   DataReady_ptr resptr=DataReady_ptr(result);
//   int sample;
//   int resoffset;
//   int totalsamples=getNumSamples();
//   const double* res=0;
//   #pragma omp parallel for private(sample,resoffset,threadnum,res) schedule(static)
//   for (sample=0;sample<totalsamples;++sample)
//   {
// cout << "################################# " << sample << endl;
// #ifdef _OPENMP
//     res=resolveSample(v,sample,threadbuffersize*omp_get_thread_num());
// #else
//     res=resolveSample(v,sample,0);   // this would normally be v, but not if its a single IDENTITY op.
// #endif
// cerr << "-------------------------------- " << endl;
//     resoffset=result->getPointOffset(sample,0);
// cerr << "offset=" << resoffset << endl;
//     for (unsigned int i=0;i<m_samplesize;++i,++resoffset)	// copy values into the output vector
//     {
// 	resvec[resoffset]=res[i];
//     }
// cerr << "*********************************" << endl;
//   }
//   return resptr;
// }


DataReady_ptr
DataLazy::resolve()
{

cout << "Sample size=" << m_samplesize << endl;
cout << "Buffers=" << m_buffsRequired << endl;

  if (m_readytype!='E')
  {
    collapse();
  }
  if (m_op==IDENTITY)
  {
    return m_id;
  }
  	// from this point on we must have m_op!=IDENTITY and m_readytype=='E'
  size_t threadbuffersize=m_samplesize*(max(1,m_buffsRequired)+1);
  int numthreads=1;
#ifdef _OPENMP
  numthreads=getNumberOfThreads();
  int threadnum=0;
#endif 
  ValueType v(numthreads*threadbuffersize);	
cout << "Buffer created with size=" << v.size() << endl;
  DataExpanded* result=new DataExpanded(getFunctionSpace(),getShape(),  ValueType(getNoValues()));
  ValueType& resvec=result->getVector();
  DataReady_ptr resptr=DataReady_ptr(result);
  int sample;
  size_t outoffset;		// offset in the output data
  int totalsamples=getNumSamples();
  const ValueType* res=0;
  size_t resoffset=0;
  #pragma omp parallel for private(sample,resoffset,outoffset,threadnum,res) schedule(static)
  for (sample=0;sample<totalsamples;++sample)
  {
cout << "################################# " << sample << endl;
#ifdef _OPENMP
    res=resolveSample(v,threadbuffersize*omp_get_thread_num(),sample,resoffset);
#else
    res=resolveSample(v,0,sample,resoffset);   // this would normally be v, but not if its a single IDENTITY op.
#endif
cerr << "-------------------------------- " << endl;
    outoffset=result->getPointOffset(sample,0);
cerr << "offset=" << outoffset << endl;
    for (unsigned int i=0;i<m_samplesize;++i,++outoffset,++resoffset)	// copy values into the output vector
    {
	resvec[outoffset]=(*res)[resoffset];
    }
cerr << "*********************************" << endl;
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
	if (m_id->isExpanded())
	{
	   oss << "E";
	}
	else if (m_id->isTagged())
	{
	  oss << "T";
	}
	else if (m_id->isConstant())
	{
	  oss << "C";
	}
	else
	{
	  oss << "?";
	}
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
