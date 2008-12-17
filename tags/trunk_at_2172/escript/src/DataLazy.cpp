
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

// #define LAZYDEBUG(X) if (privdebug){X;} 
#define LAZYDEBUG(X)
namespace
{
bool privdebug=false;

#define ENABLEDEBUG privdebug=true;
#define DISABLEDEBUG privdebug=false;
}

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

To add a new operator you need to do the following (plus anything I might have forgotten - adding a new group for example):
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
   G_UNARY_P,		// pointwise operations with one argument, requiring a parameter
   G_NP1OUT,		// non-pointwise op with one output
   G_NP1OUT_P,		// non-pointwise op with one output requiring a parameter
   G_TENSORPROD		// general tensor product
};




string ES_opstrings[]={"UNKNOWN","IDENTITY","+","-","*","/","^",
			"sin","cos","tan",
			"asin","acos","atan","sinh","cosh","tanh","erf",
			"asinh","acosh","atanh",
			"log10","log","sign","abs","neg","pos","exp","sqrt",
			"1/","where>0","where<0","where>=0","where<=0", "where<>0","where=0",
			"symmetric","nonsymmetric",
			"prod",
			"transpose", "trace"};
int ES_opcount=40;
ES_opgroup opgroups[]={G_UNKNOWN,G_IDENTITY,G_BINARY,G_BINARY,G_BINARY,G_BINARY, G_BINARY,
			G_UNARY,G_UNARY,G_UNARY, //10
			G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,	// 17
			G_UNARY,G_UNARY,G_UNARY,					// 20
			G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,	// 28
			G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY, G_UNARY_P, G_UNARY_P,		// 35
			G_NP1OUT,G_NP1OUT,
			G_TENSORPROD,
			G_NP1OUT_P, G_NP1OUT_P};
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
// the shapes resulting from tensor product are more complex to compute so are worked out elsewhere
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

// return the shape for "op left"

DataTypes::ShapeType
resultShape(DataAbstract_ptr left, ES_optype op)
{
	switch(op)
	{
    	case TRANS:
		return left->getShape();
	break;
   	case TRACE:
		return DataTypes::scalarShape;
	break;
    	default:
	throw DataException("Programmer error - resultShape(left,op) can't compute shapes for operator "+opToString(op)+".");
	}
}

// determine the output shape for the general tensor product operation
// the additional parameters return information required later for the product
// the majority of this code is copy pasted from C_General_Tensor_Product
DataTypes::ShapeType
GTPShape(DataAbstract_ptr left, DataAbstract_ptr right, int axis_offset, int transpose, int& SL, int& SM, int& SR)
{
	
  // Get rank and shape of inputs
  int rank0 = left->getRank();
  int rank1 = right->getRank();
  const DataTypes::ShapeType& shape0 = left->getShape();
  const DataTypes::ShapeType& shape1 = right->getShape();

  // Prepare for the loops of the product and verify compatibility of shapes
  int start0=0, start1=0;
  if (transpose == 0)		{}
  else if (transpose == 1)	{ start0 = axis_offset; }
  else if (transpose == 2)	{ start1 = rank1-axis_offset; }
  else				{ throw DataException("DataLazy GeneralTensorProduct Constructor: Error - transpose should be 0, 1 or 2"); }

  if (rank0<axis_offset)
  {
	throw DataException("DataLazy GeneralTensorProduct Constructor: Error - rank of left < axisoffset");
  }

  // Adjust the shapes for transpose
  DataTypes::ShapeType tmpShape0(rank0);	// pre-sizing the vectors rather
  DataTypes::ShapeType tmpShape1(rank1);	// than using push_back
  for (int i=0; i<rank0; i++)	{ tmpShape0[i]=shape0[(i+start0)%rank0]; }
  for (int i=0; i<rank1; i++)	{ tmpShape1[i]=shape1[(i+start1)%rank1]; }

  // Prepare for the loops of the product
  SL=1, SM=1, SR=1;
  for (int i=0; i<rank0-axis_offset; i++)	{
    SL *= tmpShape0[i];
  }
  for (int i=rank0-axis_offset; i<rank0; i++)	{
    if (tmpShape0[i] != tmpShape1[i-(rank0-axis_offset)]) {
      throw DataException("C_GeneralTensorProduct: Error - incompatible shapes");
    }
    SM *= tmpShape0[i];
  }
  for (int i=axis_offset; i<rank1; i++)		{
    SR *= tmpShape1[i];
  }

  // Define the shape of the output (rank of shape is the sum of the loop ranges below)
  DataTypes::ShapeType shape2(rank0+rank1-2*axis_offset);	
  {			// block to limit the scope of out_index
     int out_index=0;
     for (int i=0; i<rank0-axis_offset; i++, ++out_index) { shape2[out_index]=tmpShape0[i]; } // First part of arg_0_Z
     for (int i=axis_offset; i<rank1; i++, ++out_index)   { shape2[out_index]=tmpShape1[i]; } // Last part of arg_1_Z
  }

  if (shape2.size()>ESCRIPT_MAX_DATA_RANK)
  {
     ostringstream os;
     os << "C_GeneralTensorProduct: Error - Attempt to create a rank " << shape2.size() << " object. The maximum rank is " << ESCRIPT_MAX_DATA_RANK << ".";
     throw DataException(os.str());
  }

  return shape2;
}


// determine the number of points in the result of "left op right"
// note that determining the resultLength for G_TENSORPROD is more complex and will not be processed here
// size_t
// resultLength(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op)
// {
//    switch (getOpgroup(op))
//    {
//    case G_BINARY: return left->getLength();
//    case G_UNARY: return left->getLength();
//    case G_NP1OUT: return left->getLength();
//    default: 
// 	throw DataException("Programmer Error - attempt to getLength() for operator "+opToString(op)+".");
//    }
// }

// determine the number of samples requires to evaluate an expression combining left and right
// NP1OUT needs an extra buffer because we can't write the answers over the top of the input.
// The same goes for G_TENSORPROD
// It might seem that pointwise binary ops (G_BINARY) could be written over the top of the lefths.
// This would be true were it not for the possibility that the LHS could be a scalar which needs to be examined
// multiple times
int
calcBuffs(const DataLazy_ptr& left, const DataLazy_ptr& right, ES_optype op)
{
   switch(getOpgroup(op))
   {
   case G_IDENTITY: return 1;
   case G_BINARY: return 1+max(left->getBuffsRequired(),right->getBuffsRequired()+1);
   case G_UNARY: 
   case G_UNARY_P: return max(left->getBuffsRequired(),1);
   case G_NP1OUT: return 1+max(left->getBuffsRequired(),1);
   case G_NP1OUT_P: return 1+max(left->getBuffsRequired(),1);
   case G_TENSORPROD: return 1+max(left->getBuffsRequired(),right->getBuffsRequired()+1);
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
	m_op(IDENTITY),
	m_axis_offset(0),
	m_transpose(0),
	m_SL(0), m_SM(0), m_SR(0)
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
   m_buffsRequired=1;
   m_samplesize=getNumDPPSample()*getNoValues();
   m_maxsamplesize=m_samplesize;
LAZYDEBUG(cout << "(1)Lazy created with " << m_samplesize << endl;)
}




DataLazy::DataLazy(DataAbstract_ptr left, ES_optype op)
	: parent(left->getFunctionSpace(),left->getShape()),
	m_op(op),
	m_axis_offset(0),
	m_transpose(0),
	m_SL(0), m_SM(0), m_SR(0)
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
   m_left=lleft;
   m_buffsRequired=calcBuffs(m_left, m_right,m_op);	// yeah m_right will be null at this point
   m_samplesize=getNumDPPSample()*getNoValues();
   m_maxsamplesize=max(m_samplesize,m_left->getMaxSampleSize());
}


// In this constructor we need to consider interpolation
DataLazy::DataLazy(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op)
	: parent(resultFS(left,right,op), resultShape(left,right,op)),
	m_op(op),
	m_SL(0), m_SM(0), m_SR(0)
{
   if ((getOpgroup(op)!=G_BINARY))
   {
	throw DataException("Programmer error - constructor DataLazy(left, right, op) will only process BINARY operations.");
   }

   if (getFunctionSpace()!=left->getFunctionSpace())	// left needs to be interpolated
   {
	FunctionSpace fs=getFunctionSpace();
	Data ltemp(left);
	Data tmp(ltemp,fs);
	left=tmp.borrowDataPtr();
   }
   if (getFunctionSpace()!=right->getFunctionSpace())	// right needs to be interpolated
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
   m_samplesize=getNumDPPSample()*getNoValues();
   m_maxsamplesize=max(max(m_samplesize,m_right->getMaxSampleSize()),m_left->getMaxSampleSize());	
   m_buffsRequired=calcBuffs(m_left, m_right,m_op);
LAZYDEBUG(cout << "(3)Lazy created with " << m_samplesize << endl;)
}

DataLazy::DataLazy(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op, int axis_offset, int transpose)
	: parent(resultFS(left,right,op), GTPShape(left,right, axis_offset, transpose, m_SL,m_SM, m_SR)),
	m_op(op),
	m_axis_offset(axis_offset),
	m_transpose(transpose)
{
   if ((getOpgroup(op)!=G_TENSORPROD))
   {
	throw DataException("Programmer error - constructor DataLazy(left, right, op, ax, tr) will only process BINARY operations which require parameters.");
   }
   if ((transpose>2) || (transpose<0))
   {
	throw DataException("DataLazy GeneralTensorProduct constructor: Error - transpose should be 0, 1 or 2");
   }
   if (getFunctionSpace()!=left->getFunctionSpace())	// left needs to be interpolated
   {
	FunctionSpace fs=getFunctionSpace();
	Data ltemp(left);
	Data tmp(ltemp,fs);
	left=tmp.borrowDataPtr();
   }
   if (getFunctionSpace()!=right->getFunctionSpace())	// right needs to be interpolated
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
   m_samplesize=getNumDPPSample()*getNoValues();
   m_maxsamplesize=max(max(m_samplesize,m_right->getMaxSampleSize()),m_left->getMaxSampleSize());	
   m_buffsRequired=calcBuffs(m_left, m_right,m_op);
LAZYDEBUG(cout << "(4)Lazy created with " << m_samplesize << endl;)
}


DataLazy::DataLazy(DataAbstract_ptr left, ES_optype op, int axis_offset)
	: parent(left->getFunctionSpace(), resultShape(left,op)),
	m_op(op),
	m_axis_offset(axis_offset),
	m_transpose(0),
	m_tol(0)
{
   if ((getOpgroup(op)!=G_NP1OUT_P))
   {
	throw DataException("Programmer error - constructor DataLazy(left, op, ax) will only process UNARY operations which require parameters.");
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
   m_left=lleft;
   m_buffsRequired=calcBuffs(m_left, m_right,m_op);	// yeah m_right will be null at this point
   m_samplesize=getNumDPPSample()*getNoValues();
   m_maxsamplesize=max(m_samplesize,m_left->getMaxSampleSize());
LAZYDEBUG(cout << "(5)Lazy created with " << m_samplesize << endl;)
}

DataLazy::DataLazy(DataAbstract_ptr left, ES_optype op, double tol)
	: parent(left->getFunctionSpace(), left->getShape()),
	m_op(op),
	m_axis_offset(0),
	m_transpose(0),
	m_tol(tol)
{
   if ((getOpgroup(op)!=G_UNARY_P))
   {
	throw DataException("Programmer error - constructor DataLazy(left, op, tol) will only process UNARY operations which require parameters.");
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
   m_left=lleft;
   m_buffsRequired=calcBuffs(m_left, m_right,m_op);	// yeah m_right will be null at this point
   m_samplesize=getNumDPPSample()*getNoValues();
   m_maxsamplesize=max(m_samplesize,m_left->getMaxSampleSize());
LAZYDEBUG(cout << "(6)Lazy created with " << m_samplesize << endl;)
}

DataLazy::~DataLazy()
{
}


int
DataLazy::getBuffsRequired() const
{
	return m_buffsRequired;
}


size_t
DataLazy::getMaxSampleSize() const
{
	return m_maxsamplesize;
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
  if ((getOpgroup(m_op)==G_BINARY) || (getOpgroup(m_op)==G_TENSORPROD))
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
    case NEZ:
	result=left.whereNonZero(m_tol);
	break;
    case EZ:
	result=left.whereZero(m_tol);
	break;
    case SYM:
	result=left.symmetric();
	break;
    case NSYM:
	result=left.nonsymmetric();
	break;
    case PROD:
	result=C_GeneralTensorProduct(left,right,m_axis_offset, m_transpose);
	break;
    case TRANS:
	result=left.transpose(m_axis_offset);
	break;
    case TRACE:
	result=left.trace(m_axis_offset);
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
#if defined (_WIN32) && !defined(__INTEL_COMPILER)
	throw DataException("Error - Data:: erf function is not supported on _WIN32 platforms.");
#else
	tensor_unary_operation(m_samplesize, left, result, ::erf);
	break;
#endif
   case ASINH:
#if defined (_WIN32) && !defined(__INTEL_COMPILER)
	tensor_unary_operation(m_samplesize, left, result, escript::asinh_substitute);
#else
	tensor_unary_operation(m_samplesize, left, result, ::asinh);
#endif   
	break;
   case ACOSH:
#if defined (_WIN32) && !defined(__INTEL_COMPILER)
	tensor_unary_operation(m_samplesize, left, result, escript::acosh_substitute);
#else
	tensor_unary_operation(m_samplesize, left, result, ::acosh);
#endif   
	break;
   case ATANH:
#if defined (_WIN32) && !defined(__INTEL_COMPILER)
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
// There are actually G_UNARY_P but I don't see a compelling reason to treat them differently
    case NEZ:
	tensor_unary_operation(m_samplesize, left, result, bind2nd(AbsGT(),m_tol));
	break;
    case EZ:
	tensor_unary_operation(m_samplesize, left, result, bind2nd(AbsLTE(),m_tol));
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
LAZYDEBUG(cerr << "subroffset=" << subroffset << endl;)
  const ValueType* vleft=m_left->resolveSample(v,offset+m_left->m_samplesize,sampleNo,subroffset);
  roffset=offset;
  size_t loop=0;
  size_t numsteps=(m_readytype=='E')?getNumDPPSample():1;
  size_t step=getNoValues();
  switch (m_op)
  {
    case SYM:
	for (loop=0;loop<numsteps;++loop)
	{
	    DataMaths::symmetric(*vleft,m_left->getShape(),subroffset, v, getShape(), offset);
	    subroffset+=step;
	    offset+=step;
	}
	break;
    case NSYM:
	for (loop=0;loop<numsteps;++loop)
	{
	    DataMaths::nonsymmetric(*vleft,m_left->getShape(),subroffset, v, getShape(), offset);
	    subroffset+=step;
	    offset+=step;
	}
	break;
    default:
	throw DataException("Programmer error - resolveNP1OUT can not resolve operator "+opToString(m_op)+".");
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
DataLazy::resolveNP1OUT_P(ValueType& v, size_t offset, int sampleNo, size_t& roffset) const
{
	// we assume that any collapsing has been done before we get here
	// since we only have one argument we don't need to think about only
	// processing single points.
  if (m_readytype!='E')
  {
    throw DataException("Programmer error - resolveNP1OUT_P should only be called on expanded Data.");
  }
	// since we can't write the result over the input, we need a result offset further along
  size_t subroffset;
  const ValueType* vleft=m_left->resolveSample(v,offset+m_left->m_samplesize,sampleNo,subroffset);
LAZYDEBUG(cerr << "srcsamplesize=" << offset+m_left->m_samplesize << " beg=" << subroffset << endl;)
LAZYDEBUG(cerr << "Offset for 5800=" << getPointOffset(5800/getNumDPPSample(),5800%getNumDPPSample()) << endl;)
  roffset=offset;
  size_t loop=0;
  size_t numsteps=(m_readytype=='E')?getNumDPPSample():1;
  size_t outstep=getNoValues();
  size_t instep=m_left->getNoValues();
LAZYDEBUG(cerr << "instep=" << instep << " outstep=" << outstep<< " numsteps=" << numsteps << endl;)
  switch (m_op)
  {
    case TRACE:
	for (loop=0;loop<numsteps;++loop)
	{
size_t zz=sampleNo*getNumDPPSample()+loop;
if (zz==5800)
{
LAZYDEBUG(cerr << "point=" <<  zz<< endl;)
LAZYDEBUG(cerr << "Input to  trace=" << DataTypes::pointToString(*vleft,m_left->getShape(),subroffset,"") << endl;)
LAZYDEBUG(cerr << "Offset for point=" << getPointOffset(5800/getNumDPPSample(),5800%getNumDPPSample()) << " vs ";)
LAZYDEBUG(cerr << subroffset << endl;)
LAZYDEBUG(cerr << "output=" << offset << endl;)
}
            DataMaths::trace(*vleft,m_left->getShape(),subroffset, v ,getShape(),offset,m_axis_offset);
if (zz==5800)
{
LAZYDEBUG(cerr << "Result of trace=" << DataTypes::pointToString(v,getShape(),offset,"") << endl;)
}
	    subroffset+=instep;
	    offset+=outstep;
	}
	break;
    case TRANS:
	for (loop=0;loop<numsteps;++loop)
	{
            DataMaths::transpose(*vleft,m_left->getShape(),subroffset, v,getShape(),offset,m_axis_offset);
	    subroffset+=instep;
	    offset+=outstep;
	}
	break;
    default:
	throw DataException("Programmer error - resolveNP1OUTP can not resolve operator "+opToString(m_op)+".");
  }
  return &v;
}


#define PROC_OP(TYPE,X)                               \
	for (int j=0;j<onumsteps;++j)\
	{\
	  for (int i=0;i<numsteps;++i,resultp+=resultStep) \
	  { \
LAZYDEBUG(cout << "[left,right]=[" << lroffset << "," << rroffset << "]" << endl;)\
LAZYDEBUG(cout << "{left,right}={" << (*left)[lroffset] << "," << (*right)[rroffset] << "}\n";)\
	     tensor_binary_operation< TYPE >(chunksize, &((*left)[lroffset]), &((*right)[rroffset]), resultp, X); \
LAZYDEBUG(cout << " result=      " << resultp[0] << endl;) \
	     lroffset+=leftstep; \
	     rroffset+=rightstep; \
	  }\
	  lroffset+=oleftstep;\
	  rroffset+=orightstep;\
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
LAZYDEBUG(cout << "Resolve binary: " << toString() << endl;)

  size_t lroffset=0, rroffset=0;	// offsets in the left and right result vectors
	// first work out which of the children are expanded
  bool leftExp=(m_left->m_readytype=='E');
  bool rightExp=(m_right->m_readytype=='E');
  if (!leftExp && !rightExp)
  {
	throw DataException("Programmer Error - please use collapse if neither argument has type 'E'.");
  }
  bool leftScalar=(m_left->getRank()==0);
  bool rightScalar=(m_right->getRank()==0);
  if ((m_left->getRank()!=m_right->getRank()) && (!leftScalar && !rightScalar))
  {
	throw DataException("resolveBinary - ranks of arguments must match unless one of them is scalar."); 
  }
  size_t leftsize=m_left->getNoValues();
  size_t rightsize=m_right->getNoValues();
  size_t chunksize=1;			// how many doubles will be processed in one go
  int leftstep=0;		// how far should the left offset advance after each step
  int rightstep=0;
  int numsteps=0;		// total number of steps for the inner loop
  int oleftstep=0;	// the o variables refer to the outer loop
  int orightstep=0;	// The outer loop is only required in cases where there is an extended scalar
  int onumsteps=1;
  
  bool LES=(leftExp && leftScalar);	// Left is an expanded scalar
  bool RES=(rightExp && rightScalar);
  bool LS=(!leftExp && leftScalar);	// left is a single scalar
  bool RS=(!rightExp && rightScalar);
  bool LN=(!leftExp && !leftScalar);	// left is a single non-scalar
  bool RN=(!rightExp && !rightScalar);
  bool LEN=(leftExp && !leftScalar);	// left is an expanded non-scalar
  bool REN=(rightExp && !rightScalar);

  if ((LES && RES) || (LEN && REN))	// both are Expanded scalars or both are expanded non-scalars
  {
	chunksize=m_left->getNumDPPSample()*leftsize;
	leftstep=0;
	rightstep=0;
	numsteps=1;
  }
  else if (LES || RES)
  {
	chunksize=1;
	if (LES)		// left is an expanded scalar
	{
		if (RS)
		{
		   leftstep=1;
		   rightstep=0;
		   numsteps=m_left->getNumDPPSample();
		}
		else		// RN or REN
		{
		   leftstep=0;
		   oleftstep=1;
		   rightstep=1;
		   orightstep=(RN ? -(int)rightsize : 0);
		   numsteps=rightsize;
		   onumsteps=m_left->getNumDPPSample();
		}
	}
	else		// right is an expanded scalar
	{
		if (LS)
		{
		   rightstep=1;
		   leftstep=0;
		   numsteps=m_right->getNumDPPSample();
		}
		else
		{
		   rightstep=0;
		   orightstep=1;
		   leftstep=1;
		   oleftstep=(RN ? -(int)leftsize : 0);
		   numsteps=leftsize;
		   onumsteps=m_right->getNumDPPSample();
		}
	}
  }
  else 	// this leaves (LEN, RS), (LEN, RN) and their transposes
  {
	if (LEN)	// and Right will be a single value 
	{
		chunksize=rightsize;
		leftstep=rightsize;
	   	rightstep=0;
		numsteps=m_left->getNumDPPSample();
		if (RS)
		{
		   numsteps*=leftsize;
		}
	}
	else	// REN
	{
		chunksize=leftsize;
		rightstep=leftsize;
		leftstep=0;
		numsteps=m_right->getNumDPPSample();
		if (LS)
		{
		   numsteps*=rightsize;
		}
	}
  }

  int resultStep=max(leftstep,rightstep);	// only one (at most) should be !=0
	// Get the values of sub-expressions
  const ValueType* left=m_left->resolveSample(v,offset+getMaxSampleSize(),sampleNo,lroffset);	// see note on 
	// calcBufss for why we can't put offset as the 2nd param above
  const ValueType* right=m_right->resolveSample(v,offset+2*getMaxSampleSize(),sampleNo,rroffset); // Note
	// the right child starts further along.
LAZYDEBUG(cout << "Post sub calls in " << toString() << endl;)
LAZYDEBUG(cout << "shapes=" << DataTypes::shapeToString(m_left->getShape()) << "," << DataTypes::shapeToString(m_right->getShape()) << endl;)
LAZYDEBUG(cout << "chunksize=" << chunksize << endl << "leftstep=" << leftstep << " rightstep=" << rightstep;)
LAZYDEBUG(cout << " numsteps=" << numsteps << endl << "oleftstep=" << oleftstep << " orightstep=" << orightstep;)
LAZYDEBUG(cout << "onumsteps=" << onumsteps << endl;)
LAZYDEBUG(cout << " DPPS=" << m_left->getNumDPPSample() << "," <<m_right->getNumDPPSample() << endl;)
LAZYDEBUG(cout << "" << LS << RS << LN << RN << LES << RES <<LEN << REN <<   endl;)
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
// This method assumes that any subexpressions which evaluate to Constant or Tagged Data
// have already been collapsed to IDENTITY. So we must have at least one expanded child.
// unlike the other resolve helpers, we must treat these datapoints separately.
DataTypes::ValueType*
DataLazy::resolveTProd(ValueType& v,  size_t offset, int sampleNo, size_t& roffset) const
{
LAZYDEBUG(cout << "Resolve TensorProduct: " << toString() << endl;)

  size_t lroffset=0, rroffset=0;	// offsets in the left and right result vectors
	// first work out which of the children are expanded
  bool leftExp=(m_left->m_readytype=='E');
  bool rightExp=(m_right->m_readytype=='E');
  int steps=getNumDPPSample();
  int leftStep=((leftExp && !rightExp)? m_right->getNoValues() : 0);
  int rightStep=((rightExp && !leftExp)? m_left->getNoValues() : 0);
  int resultStep=max(leftStep,rightStep);	// only one (at most) should be !=0
	// Get the values of sub-expressions (leave a gap of one sample for the result).
  int gap=offset+m_left->getMaxSampleSize();	// actually only needs to be m_left->m_samplesize
  const ValueType* left=m_left->resolveSample(v,gap,sampleNo,lroffset);
  gap+=m_right->getMaxSampleSize();
  const ValueType* right=m_right->resolveSample(v,gap,sampleNo,rroffset); 
LAZYDEBUG(cout << "Post sub calls: " << toString() << endl;)
LAZYDEBUG(cout << "LeftExp=" << leftExp << " rightExp=" << rightExp << endl;)
LAZYDEBUG(cout << "LeftR=" << m_left->getRank() << " rightExp=" << m_right->getRank() << endl;)
LAZYDEBUG(cout << "LeftSize=" << m_left->getNoValues() << " RightSize=" << m_right->getNoValues() << endl;)
LAZYDEBUG(cout << "m_samplesize=" << m_samplesize << endl;)
LAZYDEBUG(cout << "outputshape=" << DataTypes::shapeToString(getShape()) << endl;)
  double* resultp=&(v[offset]);		// results are stored at the vector offset we recieved
  switch(m_op)
  {
    case PROD:
	for (int i=0;i<steps;++i,resultp+=resultStep)
	{
LAZYDEBUG(cout << "lroffset=" << lroffset << "rroffset=" << rroffset << endl;)
LAZYDEBUG(cout << "l*=" << left << " r*=" << right << endl;)
LAZYDEBUG(cout << "m_SL=" << m_SL << " m_SM=" << m_SM << " m_SR=" << m_SR << endl;) 
    	  const double *ptr_0 = &((*left)[lroffset]);
    	  const double *ptr_1 = &((*right)[rroffset]);
LAZYDEBUG(cout << DataTypes::pointToString(*left, m_left->getShape(),lroffset,"LEFT") << endl;)
LAZYDEBUG(cout << DataTypes::pointToString(*right,m_right->getShape(),rroffset, "RIGHT") << endl;)
    	  matrix_matrix_product(m_SL, m_SM, m_SR, ptr_0, ptr_1, resultp, m_transpose);
	  lroffset+=leftStep;
	  rroffset+=rightStep;
	}
	break;
    default:
	throw DataException("Programmer error - resolveTProduct can not resolve operator "+opToString(m_op)+".");
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
LAZYDEBUG(cout << "Resolve sample " << toString() << endl;)
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
LAZYDEBUG(cout << "Finish  sample " << toString() << endl;)
	return &(vec);
    }
    roffset=m_id->getPointOffset(sampleNo, 0);
LAZYDEBUG(cout << "Finish  sample " << toString() << endl;)
    return &(vec);
  }
  if (m_readytype!='E')
  {
    throw DataException("Programmer Error - Collapse did not produce an expanded node.");
  }
  switch (getOpgroup(m_op))
  {
  case G_UNARY:
  case G_UNARY_P: return resolveUnary(v, offset,sampleNo,roffset);
  case G_BINARY: return resolveBinary(v, offset,sampleNo,roffset);
  case G_NP1OUT: return resolveNP1OUT(v, offset, sampleNo,roffset);
  case G_NP1OUT_P: return resolveNP1OUT_P(v, offset, sampleNo,roffset);
  case G_TENSORPROD: return resolveTProd(v,offset, sampleNo,roffset);
  default:
    throw DataException("Programmer Error - resolveSample does not know how to process "+opToString(m_op)+".");
  }

}


// To simplify the memory management, all threads operate on one large vector, rather than one each.
// Each sample is evaluated independently and copied into the result DataExpanded.
DataReady_ptr
DataLazy::resolve()
{

LAZYDEBUG(cout << "Sample size=" << m_samplesize << endl;)
LAZYDEBUG(cout << "Buffers=" << m_buffsRequired << endl;)

  if (m_readytype!='E')		// if the whole sub-expression is Constant or Tagged, then evaluate it normally
  {
    collapse();
  }
  if (m_op==IDENTITY)		// So a lazy expression of Constant or Tagged data will be returned here. 
  {
    return m_id;
  }
  	// from this point on we must have m_op!=IDENTITY and m_readytype=='E'
  size_t threadbuffersize=m_maxsamplesize*(max(1,m_buffsRequired));	// Each thread needs to have enough
	// storage to evaluate its expression
  int numthreads=1;
#ifdef _OPENMP
  numthreads=getNumberOfThreads();
#endif 
  ValueType v(numthreads*threadbuffersize);	
LAZYDEBUG(cout << "Buffer created with size=" << v.size() << endl;)
  DataExpanded* result=new DataExpanded(getFunctionSpace(),getShape(),  ValueType(getNoValues()));
  ValueType& resvec=result->getVector();
  DataReady_ptr resptr=DataReady_ptr(result);
  int sample;
  size_t outoffset;		// offset in the output data
  int totalsamples=getNumSamples();
  const ValueType* res=0;	// Vector storing the answer
  size_t resoffset=0;		// where in the vector to find the answer
LAZYDEBUG(cout << "Total number of samples=" <<totalsamples << endl;)
  #pragma omp parallel for private(sample,resoffset,outoffset,res) schedule(static)
  for (sample=0;sample<totalsamples;++sample)
  {
      if (sample==0)  {ENABLEDEBUG}

//       if (sample==5800/getNumDPPSample())  {ENABLEDEBUG}
LAZYDEBUG(cout << "################################# " << sample << endl;)
#ifdef _OPENMP
    res=resolveSample(v,threadbuffersize*omp_get_thread_num(),sample,resoffset);
#else
    res=resolveSample(v,0,sample,resoffset);   // res would normally be v, but not if its a single IDENTITY op.
#endif
LAZYDEBUG(cerr << "-------------------------------- " << endl;)
LAZYDEBUG(cerr<< "Copying sample#" << sample << endl;)
    outoffset=result->getPointOffset(sample,0);
LAZYDEBUG(cerr << "offset=" << outoffset << " from offset=" << resoffset << " " << m_samplesize << " doubles" << endl;)
    for (unsigned int i=0;i<m_samplesize;++i,++outoffset,++resoffset)	// copy values into the output vector
    {
LAZYDEBUG(cerr << "outoffset=" << outoffset << " resoffset=" << resoffset << " " << (*res)[resoffset]<< endl;)
	resvec[outoffset]=(*res)[resoffset];
    }
LAZYDEBUG(cerr << DataTypes::pointToString(resvec,getShape(),outoffset-m_samplesize+DataTypes::noValues(getShape()),"Final result:") << endl;)
LAZYDEBUG(cerr << "*********************************" << endl;)
    DISABLEDEBUG
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
  case G_UNARY_P:
  case G_NP1OUT:
  case G_NP1OUT_P:
	oss << opToString(m_op) << '(';
	m_left->intoString(oss);
	oss << ')';
	break;
  case G_TENSORPROD:
	oss << opToString(m_op) << '(';
	m_left->intoString(oss);
	oss << ", ";
	m_right->intoString(oss);
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
  case G_NP1OUT: return new DataLazy(m_left->deepCopy()->getPtr(), m_right->deepCopy()->getPtr(),m_op);
  case G_TENSORPROD: return new DataLazy(m_left->deepCopy()->getPtr(), m_right->deepCopy()->getPtr(), m_op, m_axis_offset, m_transpose);
  default:
	throw DataException("Programmer error - do not know how to deepcopy operator "+opToString(m_op)+".");
  }
}


// There is no single, natural interpretation of getLength on DataLazy.
// Instances of DataReady can look at the size of their vectors.
// For lazy though, it could be the size the data would be if it were resolved;
// or it could be some function of the lengths of the DataReady instances which 
// form part of the expression.
// Rather than have people making assumptions, I have disabled the method.
DataTypes::ValueType::size_type
DataLazy::getLength() const
{
  throw DataException("getLength() does not make sense for lazy data.");
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
