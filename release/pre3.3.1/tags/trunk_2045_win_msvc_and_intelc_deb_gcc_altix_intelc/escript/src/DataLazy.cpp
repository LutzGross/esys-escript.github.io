
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
#include "Utils.h"

/*
How does DataLazy work?
~~~~~~~~~~~~~~~~~~~~~~~

Each instance represents a single operation on one or two other DataLazy instances. These arguments are normally
denoted left and right.

A special operation, IDENTITY, stores an instance of DataReady in the m_id member.
This means that all "internal" nodes in the structure are instances of DataLazy.

Each operation has a string representation as well as an opgroup - eg G_IDENTITY, G_BINARY, ...
Note that IDENITY is not considered a unary operation.

I am avoiding calling the structure formed a tree because it is not guaranteed to be one (eg c=a+a). 
It must however form a DAG (directed acyclic graph).
I will refer to individual DataLazy objects with the structure as nodes.

Each node also stores:
- m_readytype \in {'E','T','C','?'} ~ indicates what sort of DataReady would be produced if the expression was
	evaluated.
- m_length ~ how many values would be stored in the answer if the expression was evaluated.
- m_buffsrequired ~ the larged number of samples which would need to be kept simultaneously in order to
	evaluate the expression.
- m_samplesize ~ the number of doubles stored in a sample.

When a new node is created, the above values are computed based on the values in the child nodes.
Eg: if left requires 4 samples and right requires 6 then left+right requires 7 samples.

The resolve method, which produces a DataReady from a DataLazy, does the following:
1) Create a DataReady to hold the new result.
2) Allocate a vector (v) big enough to hold m_buffsrequired samples.
3) For each sample, call resolveSample with v, to get its values and copy them into the result object.

(In the case of OMP, multiple samples are resolved in parallel so the vector needs to be larger.)

resolveSample returns a Vector* and an offset within that vector where the result is stored.
Normally, this would be v, but for identity nodes their internal vector is returned instead.

The convention that I use, is that the resolve methods should store their results starting at the offset they are passed.

For expressions which evaluate to Constant or Tagged, there is a different evaluation method.
The collapse method invokes the (non-lazy) operations on the Data class to evaluate the expression.

To add a new operator you need to do the following (plus anything I might have forgotten):
1) Add to the ES_optype.
2) determine what opgroup your operation belongs to (X)
3) add a string for the op to the end of ES_opstrings
4) increase ES_opcount
5) add an entry (X) to opgroups
6) add an entry to the switch in collapseToReady
7) add an entry to resolveX
*/


using namespace std;
using namespace boost;

namespace escript
{

namespace
{

enum ES_opgroup
{
   G_UNKNOWN,
   G_IDENTITY,
   G_BINARY,		// pointwise operations with two arguments
   G_UNARY,		// pointwise operations with one argument
   G_NP1OUT		// non-pointwise op with one output
};




string ES_opstrings[]={"UNKNOWN","IDENTITY","+","-","*","/","^",
			"sin","cos","tan",
			"asin","acos","atan","sinh","cosh","tanh","erf",
			"asinh","acosh","atanh",
			"log10","log","sign","abs","neg","pos","exp","sqrt",
			"1/","where>0","where<0","where>=0","where<=0",
			"symmetric","nonsymmetric"};
int ES_opcount=35;
ES_opgroup opgroups[]={G_UNKNOWN,G_IDENTITY,G_BINARY,G_BINARY,G_BINARY,G_BINARY, G_BINARY,
			G_UNARY,G_UNARY,G_UNARY, //10
			G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,	// 17
			G_UNARY,G_UNARY,G_UNARY,					// 20
			G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,		// 28
			G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,			// 33
			G_NP1OUT,G_NP1OUT};
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

  FunctionSpace l=left->getFunctionSpace();
  FunctionSpace r=right->getFunctionSpace();
  if (l!=r)
  {
    if (r.probeInterpolation(l))
    {
	return l;
    }
    if (l.probeInterpolation(r))
    {
	return r;
    }
    throw DataException("Cannot interpolate between the FunctionSpaces given for operation "+opToString(op)+".");
  }
  return l;
}

// return the shape of the result of "left op right"
DataTypes::ShapeType
resultShape(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op)
{
	if (left->getShape()!=right->getShape())
	{
	  if ((getOpgroup(op)!=G_BINARY) && (getOpgroup(op)!=G_NP1OUT))
	  {
		throw DataException("Shapes not the name - shapes must match for (point)binary operations.");
	  }
	  if (left->getRank()==0)	// we need to allow scalar * anything
	  {
		return right->getShape();
	  }
	  if (right->getRank()==0)
	  {
		return left->getShape();
	  }
	  throw DataException("Shapes not the same - arguments must have matching shapes (or be scalars) for (point)binary operations on lazy data.");
	}
	return left->getShape();
}

// determine the number of points in the result of "left op right"
size_t
resultLength(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op)
{
   switch (getOpgroup(op))
   {
   case G_BINARY: return left->getLength();
   case G_UNARY: return left->getLength();
   case G_NP1OUT: return left->getLength();
   default: 
	throw DataException("Programmer Error - attempt to getLength() for operator "+opToString(op)+".");
   }
}

// determine the number of samples requires to evaluate an expression combining left and right
// NP1OUT needs an extra buffer because we can't write the answers over the top of the input.
int
calcBuffs(const DataLazy_ptr& left, const DataLazy_ptr& right, ES_optype op)
{
   switch(getOpgroup(op))
   {
   case G_IDENTITY: return 1;
   case G_BINARY: return max(left->getBuffsRequired(),right->getBuffsRequired()+1);
   case G_UNARY: return max(left->getBuffsRequired(),1);
   case G_NP1OUT: return 1+max(left->getBuffsRequired(),1);
   default: 
	throw DataException("Programmer Error - attempt to calcBuffs() for operator "+opToString(op)+".");
   }
}


}	// end anonymous namespace



// Return a string representing the operation
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
   if ((getOpgroup(op)!=G_UNARY) && (getOpgroup(op)!=G_NP1OUT))
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
   m_buffsRequired=calcBuffs(m_left, m_right,m_op);	// yeah m_right will be null at this point
   m_samplesize=getNumDPPSample()*getNoValues();
}


// In this constructor we need to consider interpolation
DataLazy::DataLazy(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op)
	: parent(resultFS(left,right,op), resultShape(left,right,op)),
	m_op(op)
{
   if ((getOpgroup(op)!=G_BINARY))
   {
cout << opToString(op) << endl;
	throw DataException("Programmer error - constructor DataLazy(left, right, op) will only process BINARY operations.");
   }

   if (getFunctionSpace()!=left->getFunctionSpace())	// left needs to be interpolated
   {
	FunctionSpace fs=getFunctionSpace();
	Data ltemp(left);
	Data tmp(ltemp,fs);
	left=tmp.borrowDataPtr();
   }
   if (getFunctionSpace()!=right->getFunctionSpace())	// left needs to be interpolated
   {
	Data tmp(Data(right),getFunctionSpace());
	right=tmp.borrowDataPtr();
   }
   left->operandCheck(*right);

   if (left->isLazy())			// the children need to be DataLazy. Wrap them in IDENTITY if required
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


/*
  \brief Evaluates the expression using methods on Data.
  This does the work for the collapse method.
  For reasons of efficiency do not call this method on DataExpanded nodes.
*/
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
    case SYM:
	result=left.symmetric();
	break;
    case NSYM:
	result=left.nonsymmetric();
	break;
    default:
	throw DataException("Programmer error - collapseToReady does not know how to resolve operator "+opToString(m_op)+".");
  }
  return result.borrowReadyPtr();
}

/*
   \brief Converts the DataLazy into an IDENTITY storing the value of the expression.
   This method uses the original methods on the Data class to evaluate the expressions.
   For this reason, it should not be used on DataExpanded instances. (To do so would defeat
   the purpose of using DataLazy in the first place).
*/
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

/*
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
	tensor_unary_operation<double (*)(double)>(m_samplesize, left, result, ::sin);
	break;
    case COS:
	tensor_unary_operation<double (*)(double)>(m_samplesize, left, result, ::cos);
	break;
    case TAN:
	tensor_unary_operation<double (*)(double)>(m_samplesize, left, result, ::tan);
	break;
    case ASIN:
	tensor_unary_operation<double (*)(double)>(m_samplesize, left, result, ::asin);
	break;
    case ACOS:
	tensor_unary_operation<double (*)(double)>(m_samplesize, left, result, ::acos);
	break;
    case ATAN:
	tensor_unary_operation<double (*)(double)>(m_samplesize, left, result, ::atan);
	break;
    case SINH:
	tensor_unary_operation<double (*)(double)>(m_samplesize, left, result, ::sinh);
	break;
    case COSH:
	tensor_unary_operation<double (*)(double)>(m_samplesize, left, result, ::cosh);
	break;
    case TANH:
	tensor_unary_operation<double (*)(double)>(m_samplesize, left, result, ::tanh);
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
	tensor_unary_operation<double (*)(double)>(m_samplesize, left, result, ::log10);
	break;
    case LOG:
	tensor_unary_operation<double (*)(double)>(m_samplesize, left, result, ::log);
	break;
    case SIGN:
	tensor_unary_operation(m_samplesize, left, result, escript::fsign);
	break;
    case ABS:
	tensor_unary_operation<double (*)(double)>(m_samplesize, left, result, ::fabs);
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
	tensor_unary_operation<double (*)(double)>(m_samplesize, left, result, ::exp);
	break;
    case SQRT:
	tensor_unary_operation<double (*)(double)>(m_samplesize, left, result, ::sqrt);
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


/*
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
DataLazy::resolveNP1OUT(ValueType& v, size_t offset, int sampleNo, size_t& roffset) const
{
	// we assume that any collapsing has been done before we get here
	// since we only have one argument we don't need to think about only
	// processing single points.
  if (m_readytype!='E')
  {
    throw DataException("Programmer error - resolveNP1OUT should only be called on expanded Data.");
  }
	// since we can't write the result over the input, we need a result offset further along
  size_t subroffset=roffset+m_samplesize;
  const ValueType* vleft=m_left->resolveSample(v,offset,sampleNo,subroffset);
  roffset=offset;
  switch (m_op)
  {
    case SYM:
	DataMaths::symmetric(*vleft,m_left->getShape(),subroffset, v, getShape(), offset);
	break;
    case NSYM:
	DataMaths::nonsymmetric(*vleft,m_left->getShape(),subroffset, v, getShape(), offset);
	break;
    default:
	throw DataException("Programmer error - resolveNP1OUT can not resolve operator "+opToString(m_op)+".");
  }
  return &v;
}




#define PROC_OP(TYPE,X)                               \
	for (int i=0;i<steps;++i,resultp+=resultStep) \
	{ \
	   tensor_binary_operation< TYPE >(chunksize, &((*left)[lroffset]), &((*right)[rroffset]), resultp, X); \
	   lroffset+=leftStep; \
	   rroffset+=rightStep; \
	}

/*
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
// This method assumes that any subexpressions which evaluate to Constant or Tagged Data
// have already been collapsed to IDENTITY. So we must have at least one expanded child.
// If both children are expanded, then we can process them in a single operation (we treat
// the whole sample as one big datapoint.
// If one of the children is not expanded, then we need to treat each point in the sample
// individually.
// There is an additional complication when scalar operations are considered.
// For example, 2+Vector.
// In this case each double within the point is treated individually
DataTypes::ValueType*
DataLazy::resolveBinary(ValueType& v,  size_t offset, int sampleNo, size_t& roffset) const
{
cout << "Resolve binary: " << toString() << endl;

  size_t lroffset=0, rroffset=0;	// offsets in the left and right result vectors
	// first work out which of the children are expanded
  bool leftExp=(m_left->m_readytype=='E');
  bool rightExp=(m_right->m_readytype=='E');
  bool bigloops=((leftExp && rightExp) || (!leftExp && !rightExp));	// is processing in single step?
  int steps=(bigloops?1:getNumDPPSample());
  size_t chunksize=(bigloops? m_samplesize : getNoValues());	// if bigloops, pretend the whole sample is a datapoint
  if (m_left->getRank()!=m_right->getRank())	// need to deal with scalar * ? ops
  {
 	EsysAssert((m_left->getRank()==0) || (m_right->getRank()==0), "Error - Ranks must match unless one is 0.");
	steps=getNumDPPSample()*max(m_left->getNoValues(),m_right->getNoValues());
	chunksize=1;	// for scalar
  }		
  int leftStep=((leftExp && !rightExp)? m_right->getNoValues() : 0);
  int rightStep=((rightExp && !leftExp)? m_left->getNoValues() : 0);
  int resultStep=max(leftStep,rightStep);	// only one (at most) should be !=0
	// Get the values of sub-expressions
  const ValueType* left=m_left->resolveSample(v,offset,sampleNo,lroffset);
  const ValueType* right=m_right->resolveSample(v,offset+m_samplesize,sampleNo,rroffset); // Note
	// the right child starts further along.
  double* resultp=&(v[offset]);		// results are stored at the vector offset we recieved
  switch(m_op)
  {
    case ADD:
        PROC_OP(NO_ARG,plus<double>());
	break;
    case SUB:
	PROC_OP(NO_ARG,minus<double>());
	break;
    case MUL:
	PROC_OP(NO_ARG,multiplies<double>());
	break;
    case DIV:
	PROC_OP(NO_ARG,divides<double>());
	break;
    case POW:
       PROC_OP(double (double,double),::pow);
	break;
    default:
	throw DataException("Programmer error - resolveBinary can not resolve operator "+opToString(m_op)+".");
  }
  roffset=offset;	
  return &v;
}



/*
  \brief Compute the value of the expression for the given sample.
  \return Vector which stores the value of the subexpression for the given sample.
  \param v A vector to store intermediate results.
  \param offset Index in v to begin storing results.
  \param sampleNo Sample number to evaluate.
  \param roffset (output parameter) the offset in the return vector where the result begins.

  The return value will be an existing vector so do not deallocate it.
*/
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


// To simplify the memory management, all threads operate on one large vector, rather than one each.
// Each sample is evaluated independently and copied into the result DataExpanded.
DataReady_ptr
DataLazy::resolve()
{

cout << "Sample size=" << m_samplesize << endl;
cout << "Buffers=" << m_buffsRequired << endl;

  if (m_readytype!='E')		// if the whole sub-expression is Constant or Tagged, then evaluate it normally
  {
    collapse();
  }
  if (m_op==IDENTITY)		// So a lazy expression of Constant or Tagged data will be returned here. 
  {
    return m_id;
  }
  	// from this point on we must have m_op!=IDENTITY and m_readytype=='E'
  size_t threadbuffersize=m_samplesize*(max(1,m_buffsRequired));	// Each thread needs to have enough
	// storage to evaluate its expression
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
  const ValueType* res=0;	// Vector storing the answer
  size_t resoffset=0;		// where in the vector to find the answer
  #pragma omp parallel for private(sample,resoffset,outoffset,threadnum,res) schedule(static)
  for (sample=0;sample<totalsamples;++sample)
  {
cout << "################################# " << sample << endl;
#ifdef _OPENMP
    res=resolveSample(v,threadbuffersize*omp_get_thread_num(),sample,resoffset);
#else
    res=resolveSample(v,0,sample,resoffset);   // res would normally be v, but not if its a single IDENTITY op.
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
  case G_NP1OUT:
	oss << opToString(m_op) << '(';
	m_left->intoString(oss);
	oss << ')';
	break;
  default:
	oss << "UNKNOWN";
  }
}

DataAbstract* 
DataLazy::deepCopy()
{
  switch (getOpgroup(m_op))
  {
  case G_IDENTITY:  return new DataLazy(m_id->deepCopy()->getPtr());
  case G_UNARY:	return new DataLazy(m_left->deepCopy()->getPtr(),m_op);
  case G_BINARY:	return new DataLazy(m_left->deepCopy()->getPtr(),m_right->deepCopy()->getPtr(),m_op);
  default:
	throw DataException("Programmer error - do not know how to deepcopy operator "+opToString(m_op)+".");
  }
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


// To do this we need to rely on our child nodes
DataTypes::ValueType::size_type 
DataLazy::getPointOffset(int sampleNo,
                 int dataPointNo)
{
  if (m_op==IDENTITY)
  {
	return m_id->getPointOffset(sampleNo,dataPointNo);
  }
  if (m_readytype!='E')
  {
	collapse();
	return m_id->getPointOffset(sampleNo,dataPointNo);
  }
  // at this point we do not have an identity node and the expression will be Expanded
  // so we only need to know which child to ask
  if (m_left->m_readytype=='E')
  {
	return m_left->getPointOffset(sampleNo,dataPointNo);
  }
  else
  {
	return m_right->getPointOffset(sampleNo,dataPointNo);
  }
}

// To do this we need to rely on our child nodes
DataTypes::ValueType::size_type 
DataLazy::getPointOffset(int sampleNo,
                 int dataPointNo) const
{
  if (m_op==IDENTITY)
  {
	return m_id->getPointOffset(sampleNo,dataPointNo);
  }
  if (m_readytype=='E')
  {
    // at this point we do not have an identity node and the expression will be Expanded
    // so we only need to know which child to ask
    if (m_left->m_readytype=='E')
    {
	return m_left->getPointOffset(sampleNo,dataPointNo);
    }
    else
    {
	return m_right->getPointOffset(sampleNo,dataPointNo);
    }
  }
  if (m_readytype=='C')
  {
	return m_left->getPointOffset(sampleNo,dataPointNo); // which child doesn't matter
  }
  throw DataException("Programmer error - getPointOffset on lazy data may require collapsing (but this object is marked const).");
}

// It would seem that DataTagged will need to be treated differently since even after setting all tags
// to zero, all the tags from all the DataTags would be in the result.
// However since they all have the same value (0) whether they are there or not should not matter.
// So I have decided that for all types this method will create a constant 0.
// It can be promoted up as required.
// A possible efficiency concern might be expanded->constant->expanded which has an extra memory management
// but we can deal with that if it arrises.
void
DataLazy::setToZero()
{
  DataTypes::ValueType v(getNoValues(),0);
  m_id=DataReady_ptr(new DataConstant(getFunctionSpace(),getShape(),v));
  m_op=IDENTITY;
  m_right.reset();   
  m_left.reset();
  m_readytype='C';
  m_buffsRequired=1;
}

}	// end namespace
