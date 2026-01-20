
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "DataLazy.h"
#include "Data.h"
#include "DataTypes.h"
#include "EscriptParams.h"
#include "FunctionSpace.h"
#include "Utils.h"
#include "DataVectorOps.h"

#include <iomanip> // for some fancy formatting in debug

using namespace escript::DataTypes;

#define NO_ARG

// #define LAZYDEBUG(X) if (privdebug){X;} 
#define LAZYDEBUG(X)
namespace
{
bool privdebug=false;

#define ENABLEDEBUG privdebug=true;
#define DISABLEDEBUG privdebug=false;
}

//#define SIZELIMIT if ((m_height>escript::escriptParams.getInt("TOO_MANY_LEVELS")) || (m_children>escript::escriptParams.getInt("TOO_MANY_NODES"))) {cerr << "\n!!!!!!! SIZE LIMIT EXCEEDED " << m_children << ";" << m_height << endl << toString() << endl;resolveToIdentity();}

//#define SIZELIMIT if ((m_height>escript::escriptParams.getInt("TOO_MANY_LEVELS")) || (m_children>escript::escriptParams.getInt("TOO_MANY_NODES"))) {cerr << "SIZE LIMIT EXCEEDED " << m_height << endl;resolveToIdentity();}


#define SIZELIMIT \
    if (m_height > escript::escriptParams.getTooManyLevels()) {\
        if (escript::escriptParams.getLazyVerbose()) {\
            cerr << "SIZE LIMIT EXCEEDED height=" << m_height << endl;\
        }\
        resolveToIdentity();\
    }

/*
How does DataLazy work?
~~~~~~~~~~~~~~~~~~~~~~~

Each instance represents a single operation on one or two other DataLazy instances. These arguments are normally
denoted left and right.

A special operation, IDENTITY, stores an instance of DataReady in the m_id member.
This means that all "internal" nodes in the structure are instances of DataLazy.

Each operation has a string representation as well as an opgroup - eg G_IDENTITY, G_BINARY, ...
Note that IDENTITY is not considered a unary operation.

I am avoiding calling the structure formed a tree because it is not guaranteed to be one (eg c=a+a). 
It must however form a DAG (directed acyclic graph).
I will refer to individual DataLazy objects with the structure as nodes.

Each node also stores:
- m_readytype \in {'E','T','C','?'} ~ indicates what sort of DataReady would be produced if the expression was
        evaluated.
- m_buffsrequired ~ the large number of samples which would need to be kept simultaneously in order to
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
    
    
DataLazy_ptr makePromote(DataLazy_ptr p)
{
    if (p->isComplex())
    {
        return p;
    }
    DataLazy* temp=new DataLazy(p, PROM);
    return DataLazy_ptr(temp);
}
    

namespace
{


// enabling this will print out when ever the maximum stacksize used by resolve increases
// it assumes _OPENMP is also in use
//#define LAZY_STACK_PROF



#ifndef _OPENMP
  #ifdef LAZY_STACK_PROF
  #undef LAZY_STACK_PROF
  #endif
#endif


#ifdef LAZY_STACK_PROF
std::vector<void*> stackstart(getNumberOfThreads());
std::vector<void*> stackend(getNumberOfThreads());
size_t maxstackuse=0;
#endif


inline int max3(int a, int b, int c)
{
    int t=(a>b?a:b);
    return (t>c?t:c);

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
    signed char res=r.getDomain()->preferredInterpolationOnDomain(r.getTypeCode(), l.getTypeCode());
    if (res==1)
    {
        return l;
    }
    if (res==-1)
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

          if (left->getRank()==0)       // we need to allow scalar * anything
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
resultShape(DataAbstract_ptr left, ES_optype op, int axis_offset)
{
        switch(op)
        {
        case TRANS:
           {                    // for the scoping of variables
                const DataTypes::ShapeType& s=left->getShape();
                DataTypes::ShapeType sh;
                int rank=left->getRank();
                if (axis_offset<0 || axis_offset>rank)
                {
            stringstream e;
            e << "Error - Data::transpose must have 0 <= axis_offset <= rank=" << rank;
            throw DataException(e.str());
        }
        for (int i=0; i<rank; i++)
                {
                   int index = (axis_offset+i)%rank;
           sh.push_back(s[index]); // Append to new shape
        }
                return sh;
           }
        break;
        case TRACE:
           {
                int rank=left->getRank();
                if (rank<2)
                {
                   throw DataException("Trace can only be computed for objects with rank 2 or greater.");
                }
                if ((axis_offset>rank-2) || (axis_offset<0))
                {
                   throw DataException("Trace: axis offset must lie between 0 and rank-2 inclusive.");
                }
                if (rank==2)
                {
                   return DataTypes::scalarShape;
                }
                else if (rank==3)
                {
                   DataTypes::ShapeType sh;
                   if (axis_offset==0)
                   {
                        sh.push_back(left->getShape()[2]);
                   }
                   else         // offset==1
                   {
                        sh.push_back(left->getShape()[0]);
                   }
                   return sh;
                }
                else if (rank==4)
                {
                   DataTypes::ShapeType sh;
                   const DataTypes::ShapeType& s=left->getShape();
                   if (axis_offset==0)
                   {
                        sh.push_back(s[2]);
                        sh.push_back(s[3]);
                   }
                   else if (axis_offset==1)
                   {
                        sh.push_back(s[0]);
                        sh.push_back(s[3]);
                   }
                   else         // offset==2
                   {
                        sh.push_back(s[0]);
                        sh.push_back(s[1]);
                   }
                   return sh;
                }
                else            // unknown rank
                {
                   throw DataException("Error - Data::trace can only be calculated for rank 2, 3 or 4 object.");
                }
           }
        break;
        default:
        throw DataException("Programmer error - resultShape(left,op) can't compute shapes for operator "+opToString(op)+".");
        }
}

DataTypes::ShapeType
SwapShape(DataAbstract_ptr left, const int axis0, const int axis1)
{
     // This code taken from the Data.cpp swapaxes() method
     // Some of the checks are probably redundant here
     int axis0_tmp,axis1_tmp;
     const DataTypes::ShapeType& s=left->getShape();
     DataTypes::ShapeType out_shape;
     // Here's the equivalent of python s_out=s[axis_offset:]+s[:axis_offset]
     // which goes thru all shape vector elements starting with axis_offset (at index=rank wrap around to 0)
     int rank=left->getRank();
     if (rank<2) {
        throw DataException("Error - Data::swapaxes argument must have at least rank 2.");
     }
     if (axis0<0 || axis0>rank-1) {
        stringstream e;
        e << "Error - Data::swapaxes: axis0 must be between 0 and rank-1=" << (rank-1);
        throw DataException(e.str());
     }
     if (axis1<0 || axis1>rank-1) {
        stringstream e;
        e << "Error - Data::swapaxes: axis1 must be between 0 and rank-1=" << (rank-1);
        throw DataException(e.str());
     }
     if (axis0 == axis1) {
         throw DataException("Error - Data::swapaxes: axis indices must be different.");
     }
     if (axis0 > axis1) {
         axis0_tmp=axis1;
         axis1_tmp=axis0;
     } else {
         axis0_tmp=axis0;
         axis1_tmp=axis1;
     }
     for (int i=0; i<rank; i++) {
       if (i == axis0_tmp) {
          out_shape.push_back(s[axis1_tmp]);
       } else if (i == axis1_tmp) {
          out_shape.push_back(s[axis0_tmp]);
       } else {
          out_shape.push_back(s[i]);
       }
     }
    return out_shape;
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
  if (transpose == 0)           {}
  else if (transpose == 1)      { start0 = axis_offset; }
  else if (transpose == 2)      { start1 = rank1-axis_offset; }
  else                          { throw DataException("DataLazy GeneralTensorProduct Constructor: Error - transpose should be 0, 1 or 2"); }

  if (rank0<axis_offset)
  {
        throw DataException("DataLazy GeneralTensorProduct Constructor: Error - rank of left < axisoffset");
  }

  // Adjust the shapes for transpose
  DataTypes::ShapeType tmpShape0(rank0);        // pre-sizing the vectors rather
  DataTypes::ShapeType tmpShape1(rank1);        // than using push_back
  for (int i=0; i<rank0; i++)   { tmpShape0[i]=shape0[(i+start0)%rank0]; }
  for (int i=0; i<rank1; i++)   { tmpShape1[i]=shape1[(i+start1)%rank1]; }

  // Prepare for the loops of the product
  SL=1, SM=1, SR=1;
  for (int i=0; i<rank0-axis_offset; i++)       {
    SL *= tmpShape0[i];
  }
  for (int i=rank0-axis_offset; i<rank0; i++)   {
    if (tmpShape0[i] != tmpShape1[i-(rank0-axis_offset)]) {
      throw DataException("C_GeneralTensorProduct: Error - incompatible shapes");
    }
    SM *= tmpShape0[i];
  }
  for (int i=axis_offset; i<rank1; i++)         {
    SR *= tmpShape1[i];
  }

  // Define the shape of the output (rank of shape is the sum of the loop ranges below)
  DataTypes::ShapeType shape2(rank0+rank1-2*axis_offset);       
  {                     // block to limit the scope of out_index
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

}       // end anonymous namespace

void DataLazy::LazyNodeSetup()
{
#ifdef _OPENMP
    int numthreads=omp_get_max_threads();
    if (m_iscompl)
    {
        m_samples_c.resize(numthreads*m_samplesize);
    }
    else
    {
        m_samples_r.resize(numthreads*m_samplesize);
    }    
    m_sampleids=new int[numthreads];
    for (int i=0;i<numthreads;++i) 
    { 
        m_sampleids[i]=-1;  
    }
#else
    if (m_iscompl)
    {
        m_samples_c.resize(m_samplesize);
    }
    else
    {
        m_samples_r.resize(m_samplesize);
    }
    m_sampleids=new int[1];
    m_sampleids[0]=-1;
#endif  // _OPENMP
}


// Creates an identity node
DataLazy::DataLazy(DataAbstract_ptr p)
        : parent(p->getFunctionSpace(),p->getShape())
        ,m_sampleids(0),
        m_samples_r(1),
        m_op(IDENTITY),
        m_opgroup(getOpgroup(m_op))
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
        DataReady_ptr dr=dynamic_pointer_cast<DataReady>(p);
        makeIdentity(dr);
LAZYDEBUG(cout << "Wrapping " << dr.get() << " id=" << m_id.get() << endl;)
   }
LAZYDEBUG(cout << "(1)Lazy created with " << m_samplesize << endl;)
}

// Wrapping a unary operation
DataLazy::DataLazy(DataAbstract_ptr left, ES_optype op)
        : parent(left->getFunctionSpace(),(getOpgroup(op)!=G_REDUCTION)?left->getShape():DataTypes::scalarShape),
        m_op(op),
        m_opgroup(getOpgroup(m_op)),                  
        m_axis_offset(0),
        m_transpose(0),
        m_SL(0), m_SM(0), m_SR(0)
{
   ES_opgroup gop=getOpgroup(op);
   if ((gop!=G_UNARY) && (gop!=G_NP1OUT) && (gop!=G_REDUCTION) && (gop!=G_UNARY_C) && (gop!=G_UNARY_R))
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
   m_samplesize=getNumDPPSample()*getNoValues();
   m_children=m_left->m_children+1;
   m_height=m_left->m_height+1;
   if (gop==G_UNARY_R)
   {
       m_iscompl=false;
   }
   else if (gop==G_UNARY_C)
   {
       m_iscompl=true;
   }
   else
   {
       m_iscompl=left->isComplex();
   }
   LazyNodeSetup();
   if ((lleft->m_readytype!='E') && (op!=IDENTITY))
   {
       collapse();
   }
   SIZELIMIT
}


// In this constructor we need to consider interpolation and promotion
DataLazy::DataLazy(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op)
        : parent(resultFS(left,right,op), resultShape(left,right,op)),
        m_op(op),
        m_opgroup(getOpgroup(m_op)),                  
        m_SL(0), m_SM(0), m_SR(0)
{
LAZYDEBUG(cout << "Forming operator with " << left.get() << " " << right.get() << endl;)
   if ((getOpgroup(op)!=G_BINARY))
   {
        throw DataException("Programmer error - constructor DataLazy(left, right, op) will only process BINARY operations.");
   }

   if (getFunctionSpace()!=left->getFunctionSpace())    // left needs to be interpolated
   {
        FunctionSpace fs=getFunctionSpace();
        Data ltemp(left);
        Data tmp(ltemp,fs);
        left=tmp.borrowDataPtr();
   }
   if (getFunctionSpace()!=right->getFunctionSpace())   // right needs to be interpolated
   {
        Data tmp(Data(right),getFunctionSpace());
        right=tmp.borrowDataPtr();
LAZYDEBUG(cout << "Right interpolation required " << right.get() << endl;)
   }
   left->operandCheck(*right);

   if (left->isLazy())                  // the children need to be DataLazy. Wrap them in IDENTITY if required
   {
        m_left=dynamic_pointer_cast<DataLazy>(left);
LAZYDEBUG(cout << "Left is " << m_left->toString() << endl;)
   }
   else
   {
        m_left=DataLazy_ptr(new DataLazy(left));
LAZYDEBUG(cout << "Left " << left.get() << " wrapped " << m_left->m_id.get() << endl;)
   }
   if (right->isLazy())
   {
        m_right=dynamic_pointer_cast<DataLazy>(right);
LAZYDEBUG(cout << "Right is " << m_right->toString() << endl;)
   }
   else
   {
        m_right=DataLazy_ptr(new DataLazy(right));
LAZYDEBUG(cout << "Right " << right.get() << " wrapped " << m_right->m_id.get() << endl;)
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
   if ((lt!='E') && (m_left->m_op!=IDENTITY))
   {
       m_left->collapse();
   }
   if ((rt!='E') && (m_right->m_op!=IDENTITY))
   {
       m_right->collapse();
   }
   m_samplesize=getNumDPPSample()*getNoValues();
   m_children=m_left->m_children+m_right->m_children+2;
   m_height=max(m_left->m_height,m_right->m_height)+1;
   
   // now we need to work out if we need to promote anything
   if (left->isComplex()!=right->isComplex())
   {
       if (left->isComplex())
       {
           m_right=makePromote(m_right);
       }
       else
       {
           m_left=makePromote(m_left);
       }
   }
   m_iscompl=m_left->isComplex();
   LazyNodeSetup();
   if ((m_readytype!='E') && (m_op!=IDENTITY))
   {
	collapse();
   }
   SIZELIMIT
LAZYDEBUG(cout << "(3)Lazy created with " << m_samplesize << endl;)
}

// need to consider promotion
DataLazy::DataLazy(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op, int axis_offset, int transpose)
        : parent(resultFS(left,right,op), GTPShape(left,right, axis_offset, transpose, m_SL,m_SM, m_SR)),
        m_op(op),
        m_opgroup(getOpgroup(m_op)),                  
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
   if (getFunctionSpace()!=left->getFunctionSpace())    // left needs to be interpolated
   {
        FunctionSpace fs=getFunctionSpace();
        Data ltemp(left);
        Data tmp(ltemp,fs);
        left=tmp.borrowDataPtr();
   }
   if (getFunctionSpace()!=right->getFunctionSpace())   // right needs to be interpolated
   {
        Data tmp(Data(right),getFunctionSpace());
        right=tmp.borrowDataPtr();
   }
//    left->operandCheck(*right);

   if (left->isLazy())                  // the children need to be DataLazy. Wrap them in IDENTITY if required
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
   if ((lt!='E') && (m_left->m_op!=IDENTITY))
   {
       m_left->collapse();
   }
   if ((rt!='E') && (m_right->m_op!=IDENTITY))
   {
       m_right->collapse();
   }   
   m_samplesize=getNumDPPSample()*getNoValues();
   m_children=m_left->m_children+m_right->m_children+2;
   m_height=max(m_left->m_height,m_right->m_height)+1;

   // now we need to work out if we need to promote anything
   if (left->isComplex()!=right->isComplex())
   {
       if (left->isComplex())
       {
           m_right=makePromote(m_right);
       }
       else
       {
           m_left=makePromote(m_left);
       }
   }
   m_iscompl=m_left->isComplex();      
   LazyNodeSetup();
   if ((m_readytype!='E') && (m_op!=IDENTITY))
   {
	collapse();
   }
   SIZELIMIT
LAZYDEBUG(cout << "(4)Lazy created with " << m_samplesize << endl;)
}


DataLazy::DataLazy(DataAbstract_ptr left, ES_optype op, int axis_offset)
        : parent(left->getFunctionSpace(), resultShape(left,op, axis_offset)),
        m_op(op),
        m_opgroup(getOpgroup(m_op)),                  
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
   if ((m_left->m_readytype!='E') && (m_left->m_op!=IDENTITY))
   {
       m_left->collapse();
   }
   m_samplesize=getNumDPPSample()*getNoValues();
   m_children=m_left->m_children+1;
   m_height=m_left->m_height+1;
   m_iscompl=left->isComplex();
   LazyNodeSetup();
   if ((m_readytype!='E') && (m_op!=IDENTITY))
   {
	collapse();
   }
   SIZELIMIT
LAZYDEBUG(cout << "(5)Lazy created with " << m_samplesize << endl;)
}

// G_UNARY_P or G_UNARY_PR
DataLazy::DataLazy(DataAbstract_ptr left, ES_optype op, double tol)
        : parent(left->getFunctionSpace(), left->getShape()),
        m_op(op),
        m_opgroup(getOpgroup(m_op)),                  
        m_axis_offset(0),
        m_transpose(0),
        m_tol(tol)
{
   if ((m_opgroup!=G_UNARY_P) && (m_opgroup!=G_UNARY_PR))
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
   if ((m_left->m_readytype!='E') && (m_left->m_op!=IDENTITY))
   {
       m_left->collapse();
   }   
   m_samplesize=getNumDPPSample()*getNoValues();
   m_children=m_left->m_children+1;
   m_height=m_left->m_height+1;
   if (m_opgroup==G_UNARY_PR)
   {
       m_iscompl=false;
   }
   else
   {
       m_iscompl=left->isComplex();
   }
   LazyNodeSetup();
   if ((m_readytype!='E') && (m_op!=IDENTITY))
   {
	collapse();
   }
   SIZELIMIT
LAZYDEBUG(cout << "(6)Lazy created with " << m_samplesize << endl;)
}


DataLazy::DataLazy(DataAbstract_ptr left, ES_optype op, const int axis0, const int axis1)
        : parent(left->getFunctionSpace(), SwapShape(left,axis0,axis1)),
        m_op(op),
        m_opgroup(getOpgroup(m_op)),                  
        m_axis_offset(axis0),
        m_transpose(axis1),
        m_tol(0)
{
   if ((getOpgroup(op)!=G_NP1OUT_2P))
   {
        throw DataException("Programmer error - constructor DataLazy(left, op, tol) will only process UNARY operations which require two integer parameters.");
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
   if ((m_left->m_readytype!='E') && (m_left->m_op!=IDENTITY))
   {
       m_left->collapse();
   }   
   m_samplesize=getNumDPPSample()*getNoValues();
   m_children=m_left->m_children+1;
   m_height=m_left->m_height+1;
   m_iscompl=left->isComplex();
   LazyNodeSetup();
   if ((m_readytype!='E') && (m_op!=IDENTITY))
   {
	collapse();
   }
   SIZELIMIT
LAZYDEBUG(cout << "(7)Lazy created with " << m_samplesize << endl;)
}




DataLazy::DataLazy(DataAbstract_ptr mask, DataAbstract_ptr left, DataAbstract_ptr right/*, double tol*/)
        : parent(left->getFunctionSpace(), left->getShape()),
        m_op(CONDEVAL),
        m_opgroup(getOpgroup(m_op)),                  
        m_axis_offset(0),
        m_transpose(0),
        m_tol(0)
{
   DataLazy_ptr lmask;
   DataLazy_ptr lleft;
   DataLazy_ptr lright;
   if (!mask->isLazy())
   {
        lmask=DataLazy_ptr(new DataLazy(mask));
   }
   else
   {
        lmask=dynamic_pointer_cast<DataLazy>(mask);
   }
   if (!left->isLazy())
   {
        lleft=DataLazy_ptr(new DataLazy(left));
   }
   else
   {
        lleft=dynamic_pointer_cast<DataLazy>(left);
   }
   if (!right->isLazy())
   {
        lright=DataLazy_ptr(new DataLazy(right));
   }
   else
   {
        lright=dynamic_pointer_cast<DataLazy>(right);
   }
   m_readytype=lmask->m_readytype;
   if ((lleft->m_readytype!=lright->m_readytype) || (lmask->m_readytype!=lleft->m_readytype))
   {
        throw DataException("Programmer Error - condEval arguments must have the same readytype");
   }
   m_left=lleft;
   m_right=lright;
   m_mask=lmask;
   m_samplesize=getNumDPPSample()*getNoValues();
   m_children=m_left->m_children+m_right->m_children+m_mask->m_children+1;
   m_height=max3(m_left->m_height,m_right->m_height,m_mask->m_height)+1;
   
   // now we need to work out if we need to promote anything
   if (left->isComplex()!=right->isComplex())
   {
       if (left->isComplex())
       {
           m_right=makePromote(m_right);
       }
       else
       {
           m_left=makePromote(m_left);
       }
   }
   if ((m_left->m_readytype!='E') && (m_left->m_op!=IDENTITY))
   {
       m_left->collapse();
   }
   if ((m_right->m_readytype!='E') && (m_right->m_op!=IDENTITY))
   {
       m_right->collapse();
   }   
   if ((m_mask->m_readytype!='E') && (m_mask->m_op!=IDENTITY))
   {
       m_mask->collapse();
   }   
   m_iscompl=left->isComplex();   
   LazyNodeSetup();
   if ((m_readytype!='E') && (m_op!=IDENTITY))
   {
	collapse();
   }
   SIZELIMIT
LAZYDEBUG(cout << "(8)Lazy created with " << m_samplesize << endl;)
}



DataLazy::~DataLazy()
{
   delete[] m_sampleids;
}


/*
  \brief Evaluates the expression using methods on Data.
  This does the work for the collapse method.
  For reasons of efficiency do not call this method on DataExpanded nodes.
*/
DataReady_ptr
DataLazy::collapseToReady() const
{
  if (m_readytype=='E')
  {     // this is more an efficiency concern than anything else
    throw DataException("Programmer Error - do not use collapse on Expanded data.");
  }
  if (m_op==IDENTITY)
  {
    return m_id;
  }
  DataReady_ptr pleft=m_left->collapseToReady();
  Data left(pleft);
  Data right;
  if ((m_opgroup==G_BINARY) || (m_opgroup==G_TENSORPROD))
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
    case POW:
        result=left.powD(right);
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
        result=left.antisymmetric();
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
    case SWAP:
        result=left.swapaxes(m_axis_offset, m_transpose);
        break;
    case MINVAL:
        result=left.minval();
        break;
    case MAXVAL:
        result=left.minval();
        break;
    case HER:
        result=left.hermitian();
        break;
    case NHER:
        result=left.antihermitian();
        break;
    case PROM:
        result.copy(left);
        result.complicate();
        break;
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
DataLazy::collapse() const
{
  if (m_op==IDENTITY)
  {
        return;
  }
  if (m_readytype=='E')
  {     // this is more an efficiency concern than anything else
    throw DataException("Programmer Error - do not use collapse on Expanded data.");
  }
  m_id=collapseToReady();
  m_op=IDENTITY;
  m_opgroup=getOpgroup(m_op);
  m_children=m_height=0;
}

// The result will be stored in m_samples
// The return value is a pointer to the DataVector, offset is the offset within the return value
const DataTypes::RealVectorType*
DataLazy::resolveNodeSample(int tid, int sampleNo, size_t& roffset) const
{
LAZYDEBUG(cout << "Resolve sample " << toString() << endl;)
        // collapse so we have a 'E' node or an IDENTITY for some other type
  if (m_readytype!='E' && m_op!=IDENTITY)
  {
        collapse();
  }
  if (m_op==IDENTITY)   
  {
    const RealVectorType& vec=m_id->getVectorRO();
    roffset=m_id->getPointOffset(sampleNo, 0);
#ifdef LAZY_STACK_PROF
int x;
if (&x<stackend[omp_get_thread_num()])
{
       stackend[omp_get_thread_num()]=&x;
}
#endif
    return &(vec);
  }
  if (m_readytype!='E')
  {
    throw DataException("Programmer Error - Collapse did not produce an expanded node.");
  }
  if (m_sampleids[tid]==sampleNo)
  {
        roffset=tid*m_samplesize;
        return &(m_samples_r);            // sample is already resolved
  }
  m_sampleids[tid]=sampleNo;

  switch (m_opgroup)
  {
        // note that G_UNARY_C is not listed here - this is deliberate
        // by definition this type can not produce a real_t answer
        // so processing it here would be a bug
  case G_UNARY:
  case G_UNARY_R:
  case G_UNARY_PR:
  case G_UNARY_P: return resolveNodeUnary(tid, sampleNo, roffset);
  case G_BINARY: return resolveNodeBinary(tid, sampleNo, roffset);
  case G_NP1OUT: return resolveNodeNP1OUT(tid, sampleNo, roffset);
  case G_NP1OUT_P: return resolveNodeNP1OUT_P(tid, sampleNo, roffset);
  case G_TENSORPROD: return resolveNodeTProd(tid, sampleNo, roffset);
  case G_NP1OUT_2P: return resolveNodeNP1OUT_2P(tid, sampleNo, roffset);
  case G_REDUCTION: return resolveNodeReduction(tid, sampleNo, roffset);
  case G_CONDEVAL: return resolveNodeCondEval(tid, sampleNo, roffset);
  default:
    throw DataException("Programmer Error - resolveNodeSample does not know how to process "+opToString(m_op)+" in group "+groupToString(m_opgroup)+".");
  }
}

// The result will be stored in m_samples
// The return value is a pointer to the DataVector, offset is the offset within the return value
const DataTypes::CplxVectorType*
DataLazy::resolveNodeSampleCplx(int tid, int sampleNo, size_t& roffset) const
{
LAZYDEBUG(cout << "Resolve sample " << toString() << endl;)
        // collapse so we have a 'E' node or an IDENTITY for some other type
  if (m_readytype!='E' && m_op!=IDENTITY)
  {
      throw DataException("Programmer Error - attempt to collapse inside resolveNodeSampleCplx. This should not happen.");
  }
  if (m_op==IDENTITY)   
  {
    const CplxVectorType& vec=m_id->getVectorROC();
    roffset=m_id->getPointOffset(sampleNo, 0);
#ifdef LAZY_STACK_PROF
int x;
if (&x<stackend[omp_get_thread_num()])
{
       stackend[omp_get_thread_num()]=&x;
}
#endif
    return &(vec);
  }
  if (m_readytype!='E')
  {
    throw DataException("Programmer Error - Collapse did not produce an expanded node.");
  }
  if (m_sampleids[tid]==sampleNo)
  {
        roffset=tid*m_samplesize;
        return &(m_samples_c);            // sample is already resolved
  }
  m_sampleids[tid]=sampleNo;

  switch (m_opgroup)
  {
  case G_UNARY:
  case G_UNARY_P: return resolveNodeUnaryCplx(tid, sampleNo, roffset);
  case G_UNARY_C: return resolveNodeUnary_C(tid, sampleNo, roffset);  
  case G_BINARY: return resolveNodeBinaryCplx(tid, sampleNo, roffset);
  case G_NP1OUT: return resolveNodeNP1OUTCplx(tid, sampleNo, roffset);
  case G_NP1OUT_P: return resolveNodeNP1OUT_PCplx(tid, sampleNo, roffset);
  case G_TENSORPROD: return resolveNodeTProdCplx(tid, sampleNo, roffset);
  case G_NP1OUT_2P: return resolveNodeNP1OUT_2PCplx(tid, sampleNo, roffset);
  case G_REDUCTION: return resolveNodeReductionCplx(tid, sampleNo, roffset);
  case G_CONDEVAL: return resolveNodeCondEvalCplx(tid, sampleNo, roffset);
  default:
    throw DataException("Programmer Error - resolveNodeSampleCplx does not know how to process "+opToString(m_op)+".");
  }
}

const DataTypes::RealVectorType*
DataLazy::resolveNodeUnary(int tid, int sampleNo, size_t& roffset) const
{
        // we assume that any collapsing has been done before we get here
        // since we only have one argument we don't need to think about only
        // processing single points.
        // we will also know we won't get identity nodes
  if (m_readytype!='E')
  {
    throw DataException("Programmer error - resolveUnary should only be called on expanded Data.");
  }
  if (m_op==IDENTITY)
  {
    throw DataException("Programmer error - resolveNodeUnary should not be called on identity nodes.");
  }
  if (m_op==POS)
  {
        // this should be prevented earlier
        // operation is meaningless for lazy
            throw DataException("Programmer error - POS not supported for lazy data.");    
  }
  // The tricky bit here is that we expect a real output, but for G_UNARY_R, the child nodes could be complex
  if (((m_opgroup==G_UNARY_R) || (m_opgroup==G_UNARY_PR)) && (m_left->isComplex()))
  {
    const DataTypes::CplxVectorType* leftres=m_left->resolveNodeSampleCplx(tid, sampleNo, roffset);
    const cplx_t* left=&((*leftres)[roffset]);
    roffset=m_samplesize*tid;
    real_t* result=&(m_samples_r[roffset]);
    tensor_unary_array_operation_real(m_samplesize,
                                left,
                                result,
                                m_op,
                                m_tol);
  }
  else
  {
    const DataTypes::RealVectorType* leftres=m_left->resolveNodeSample(tid, sampleNo, roffset);
    const real_t* left=&((*leftres)[roffset]);
    roffset=m_samplesize*tid;
    real_t* result=&(m_samples_r[roffset]);
    tensor_unary_array_operation(m_samplesize,
                                left,
                                result,
                                m_op,
                                m_tol);
  }  
  return &(m_samples_r);
}

const DataTypes::CplxVectorType*
DataLazy::resolveNodeUnaryCplx(int tid, int sampleNo, size_t& roffset) const
{
        // we assume that any collapsing has been done before we get here
        // since we only have one argument we don't need to think about only
        // processing single points.
        // we will also know we won't get identity nodes
  if (m_readytype!='E')
  {
    throw DataException("Programmer error - resolveUnary should only be called on expanded Data.");
  }
  if (m_op==IDENTITY)
  {
    throw DataException("Programmer error - resolveNodeUnary should not be called on identity nodes.");
  }
  if (m_op==POS)
  {
	// this should be prevented earlier
	// operation is meaningless for lazy
        throw DataException("Programmer error - POS not supported for lazy data.");    
  }
  
  roffset=m_samplesize*tid;
  DataTypes::cplx_t* result=&(m_samples_c[roffset]); 
  if (m_op==PROM)
  {
    size_t subroffset=0;
    const DataTypes::RealVectorType* leftres=m_left->resolveNodeSample(tid, sampleNo, subroffset);
    const DataTypes::real_t* left=&((*leftres)[subroffset]);      
      
    tensor_unary_promote(m_samplesize, left, result);      
  }
  else
  {
    size_t subroffset=0;
    const DataTypes::CplxVectorType* leftres=m_left->resolveNodeSampleCplx(tid, sampleNo, subroffset);
    const DataTypes::cplx_t* left=&((*leftres)[subroffset]);      
    tensor_unary_array_operation(m_samplesize,
                                left,
                                result,
                                m_op,
                                m_tol);  
  }
  return &(m_samples_c);
}

const DataTypes::RealVectorType*
DataLazy::resolveNodeReduction(int tid, int sampleNo, size_t& roffset) const
{
        // we assume that any collapsing has been done before we get here
        // since we only have one argument we don't need to think about only
        // processing single points.
        // we will also know we won't get identity nodes
  if (m_readytype!='E')
  {
    throw DataException("Programmer error - resolveUnary should only be called on expanded Data.");
  }
  if (m_op==IDENTITY)
  {
    throw DataException("Programmer error - resolveNodeUnary should not be called on identity nodes.");
  }
  size_t loffset=0;
  const DataTypes::RealVectorType* leftres=m_left->resolveNodeSample(tid, sampleNo, loffset);

  roffset=m_samplesize*tid;
  unsigned int ndpps=getNumDPPSample();
  unsigned int psize=DataTypes::noValues(m_left->getShape());
  double* result=&(m_samples_r[roffset]);
  switch (m_op)
  {
    case MINVAL:
        {
          for (unsigned int z=0;z<ndpps;++z)
          {
            FMin op;
            *result=escript::reductionOpVector(*leftres, m_left->getShape(), loffset, op, numeric_limits<double>::max());
            loffset+=psize;
            result++;
          }
        }
        break;
    case MAXVAL:
        {
          for (unsigned int z=0;z<ndpps;++z)
          {
          FMax op;
          *result=escript::reductionOpVector(*leftres, m_left->getShape(), loffset, op, numeric_limits<double>::max()*-1);
          loffset+=psize;
          result++;
          }
        }
        break;
    default:
        throw DataException("Programmer error - resolveUnary can not resolve operator "+opToString(m_op)+".");
  }
  return &(m_samples_r);
}

const DataTypes::CplxVectorType*
DataLazy::resolveNodeReductionCplx(int tid, int sampleNo, size_t& roffset) const
{
        // we assume that any collapsing has been done before we get here
        // since we only have one argument we don't need to think about only
        // processing single points.
        // we will also know we won't get identity nodes
  if (m_readytype!='E')
  {
    throw DataException("Programmer error - resolveReductionCplx should only be called on expanded Data.");
  }
  if (m_op==IDENTITY)
  {
    throw DataException("Programmer error - resolveNodeReductionCplx should not be called on identity nodes.");
  }
  throw DataException("Programmer error - reduction operations MIN and MAX not supported for complex values.");
}


const DataTypes::RealVectorType*
DataLazy::resolveNodeNP1OUT(int tid, int sampleNo, size_t& roffset) const
{
        // we assume that any collapsing has been done before we get here
        // since we only have one argument we don't need to think about only
        // processing single points.
  if (m_readytype!='E')
  {
    throw DataException("Programmer error - resolveNodeNP1OUT should only be called on expanded Data.");
  }
  if (m_op==IDENTITY)
  {
    throw DataException("Programmer error - resolveNodeNP1OUT should not be called on identity nodes.");
  }
  size_t subroffset;
  const RealVectorType* leftres=m_left->resolveNodeSample(tid, sampleNo, subroffset);
  roffset=m_samplesize*tid;
  size_t loop=0;
  size_t numsteps=(m_readytype=='E')?getNumDPPSample():1;
  size_t step=getNoValues();
  size_t offset=roffset;
  switch (m_op)
  {
    case SYM:
        for (loop=0;loop<numsteps;++loop)
        {
            escript::symmetric(*leftres,m_left->getShape(),subroffset, m_samples_r, getShape(), offset);
            subroffset+=step;
            offset+=step;
        }
        break;
    case NSYM:
        for (loop=0;loop<numsteps;++loop)
        {
            escript::antisymmetric(*leftres,m_left->getShape(),subroffset, m_samples_r, getShape(), offset);
            subroffset+=step;
            offset+=step;
        }
        break;
    default:
        throw DataException("Programmer error - resolveNP1OUT can not resolve operator "+opToString(m_op)+".");
  }
  return &m_samples_r;
}


const DataTypes::CplxVectorType*
DataLazy::resolveNodeUnary_C(int tid, int sampleNo, size_t& roffset) const
{
        // we assume that any collapsing has been done before we get here
        // since we only have one argument we don't need to think about only
        // processing single points.
  if (m_readytype!='E')
  {
    throw DataException("Programmer error - resolveNodeUnaryC should only be called on expanded Data.");
  }
  if (m_op==IDENTITY)
  {
    throw DataException("Programmer error - resolveNodeUnaryC should not be called on identity nodes.");
  }
  if (m_op!=PROM)
  {
    throw DataException("Programmer error - resolveNodeUnaryC can not resolve operator "+opToString(m_op)+".");
  }
  const DataTypes::RealVectorType* leftres=m_left->resolveNodeSample(tid, sampleNo, roffset);
  const real_t* left=&((*leftres)[roffset]);
  roffset=m_samplesize*tid;
  cplx_t* result=&(m_samples_c[roffset]);
  // rather than calling the tensor_unary_array_operation, I'm just coding this one here
  for (int i = 0; i < m_samplesize; ++i) {
    result[i] = static_cast<cplx_t>(left[i]);
  }
  return &(m_samples_c);  
}

const DataTypes::CplxVectorType*
DataLazy::resolveNodeNP1OUTCplx(int tid, int sampleNo, size_t& roffset) const
{
        // we assume that any collapsing has been done before we get here
        // since we only have one argument we don't need to think about only
        // processing single points.
  if (m_readytype!='E')
  {
    throw DataException("Programmer error - resolveNodeNP1OUT should only be called on expanded Data.");
  }
  if (m_op==IDENTITY)
  {
    throw DataException("Programmer error - resolveNodeNP1OUT should not be called on identity nodes.");
  }
  size_t subroffset;
  const CplxVectorType* leftres=m_left->resolveNodeSampleCplx(tid, sampleNo, subroffset);
  roffset=m_samplesize*tid;
  size_t loop=0;
  size_t numsteps=(m_readytype=='E')?getNumDPPSample():1;
  size_t step=getNoValues();
  size_t offset=roffset;
  switch (m_op)
  {
    case SYM:
        for (loop=0;loop<numsteps;++loop)
        {
            escript::symmetric(*leftres,m_left->getShape(),subroffset, m_samples_c, getShape(), offset);
            subroffset+=step;
            offset+=step;
        }
        break;
    case NSYM:
        for (loop=0;loop<numsteps;++loop)
        {
            escript::antisymmetric(*leftres,m_left->getShape(),subroffset, m_samples_c, getShape(), offset);
            subroffset+=step;
            offset+=step;
        }
        break;
    case HER:
        for (loop=0;loop<numsteps;++loop)
        {
            escript::hermitian(*leftres,m_left->getShape(),subroffset, m_samples_c, getShape(), offset);
            subroffset+=step;
            offset+=step;
        }
        break;
    case NHER:
        for (loop=0;loop<numsteps;++loop)
        {
            escript::antihermitian(*leftres,m_left->getShape(),subroffset, m_samples_c, getShape(), offset);
            subroffset+=step;
            offset+=step;
        }
        break;
    default:
        throw DataException("Programmer error - resolveNP1OUT can not resolve operator "+opToString(m_op)+".");
  }
  return &m_samples_c;
}

const DataTypes::RealVectorType*
DataLazy::resolveNodeNP1OUT_P(int tid, int sampleNo, size_t& roffset) const
{
        // we assume that any collapsing has been done before we get here
        // since we only have one argument we don't need to think about only
        // processing single points.
  if (m_readytype!='E')
  {
    throw DataException("Programmer error - resolveNodeNP1OUT_P should only be called on expanded Data.");
  }
  if (m_op==IDENTITY)
  {
    throw DataException("Programmer error - resolveNodeNP1OUT_P should not be called on identity nodes.");
  }
  size_t subroffset;
  size_t offset;
  const RealVectorType* leftres=m_left->resolveNodeSample(tid, sampleNo, subroffset);
  roffset=m_samplesize*tid;
  offset=roffset;
  size_t loop=0;
  size_t numsteps=(m_readytype=='E')?getNumDPPSample():1;
  size_t outstep=getNoValues();
  size_t instep=m_left->getNoValues();
  switch (m_op)
  {
    case TRACE:
        for (loop=0;loop<numsteps;++loop)
        {
            escript::trace(*leftres,m_left->getShape(),subroffset, m_samples_r ,getShape(),offset,m_axis_offset);
            subroffset+=instep;
            offset+=outstep;
        }
        break;
    case TRANS:
        for (loop=0;loop<numsteps;++loop)
        {
            escript::transpose(*leftres,m_left->getShape(),subroffset, m_samples_r, getShape(),offset,m_axis_offset);
            subroffset+=instep;
            offset+=outstep;
        }
        break;
    default:
        throw DataException("Programmer error - resolveNP1OUTP can not resolve operator "+opToString(m_op)+".");
  }
  return &m_samples_r;
}

const DataTypes::CplxVectorType*
DataLazy::resolveNodeNP1OUT_PCplx(int tid, int sampleNo, size_t& roffset) const
{
        // we assume that any collapsing has been done before we get here
        // since we only have one argument we don't need to think about only
        // processing single points.
  if (m_readytype!='E')
  {
    throw DataException("Programmer error - resolveNodeNP1OUT_P should only be called on expanded Data.");
  }
  if (m_op==IDENTITY)
  {
    throw DataException("Programmer error - resolveNodeNP1OUT_P should not be called on identity nodes.");
  }
  size_t subroffset;
  size_t offset;
  const CplxVectorType* leftres=m_left->resolveNodeSampleCplx(tid, sampleNo, subroffset);
  roffset=m_samplesize*tid;
  offset=roffset;
  size_t loop=0;
  size_t numsteps=(m_readytype=='E')?getNumDPPSample():1;
  size_t outstep=getNoValues();
  size_t instep=m_left->getNoValues();
  switch (m_op)
  {
    case TRACE:
        for (loop=0;loop<numsteps;++loop)
        {
            escript::trace(*leftres,m_left->getShape(),subroffset, m_samples_c ,getShape(),offset,m_axis_offset);
            subroffset+=instep;
            offset+=outstep;
        }
        break;
    case TRANS:
        for (loop=0;loop<numsteps;++loop)
        {
            escript::transpose(*leftres,m_left->getShape(),subroffset, m_samples_c, getShape(),offset,m_axis_offset);
            subroffset+=instep;
            offset+=outstep;
        }
        break;
    default:
        throw DataException("Programmer error - resolveNP1OUTP can not resolve operator "+opToString(m_op)+".");
  }
  return &m_samples_c;
}

const DataTypes::RealVectorType*
DataLazy::resolveNodeNP1OUT_2P(int tid, int sampleNo, size_t& roffset) const
{
  if (m_readytype!='E')
  {
    throw DataException("Programmer error - resolveNodeNP1OUT_2P should only be called on expanded Data.");
  }
  if (m_op==IDENTITY)
  {
    throw DataException("Programmer error - resolveNodeNP1OUT_2P should not be called on identity nodes.");
  }
  size_t subroffset;
  size_t offset;
  const RealVectorType* leftres=m_left->resolveNodeSample(tid, sampleNo, subroffset);
  roffset=m_samplesize*tid;
  offset=roffset;
  size_t loop=0;
  size_t numsteps=(m_readytype=='E')?getNumDPPSample():1;
  size_t outstep=getNoValues();
  size_t instep=m_left->getNoValues();
  switch (m_op)
  {
    case SWAP:
        for (loop=0;loop<numsteps;++loop)
        {
            escript::swapaxes(*leftres,m_left->getShape(),subroffset, m_samples_r, getShape(),offset, m_axis_offset, m_transpose);
            subroffset+=instep;
            offset+=outstep;
        }
        break;
    default:
        throw DataException("Programmer error - resolveNodeNP1OUT2P can not resolve operator "+opToString(m_op)+".");
  }
  return &m_samples_r;
}


const DataTypes::CplxVectorType*
DataLazy::resolveNodeNP1OUT_2PCplx(int tid, int sampleNo, size_t& roffset) const
{
  if (m_readytype!='E')
  {
    throw DataException("Programmer error - resolveNodeNP1OUT_2P should only be called on expanded Data.");
  }
  if (m_op==IDENTITY)
  {
    throw DataException("Programmer error - resolveNodeNP1OUT_2P should not be called on identity nodes.");
  }
  size_t subroffset;
  size_t offset;
  const CplxVectorType* leftres=m_left->resolveNodeSampleCplx(tid, sampleNo, subroffset);
  roffset=m_samplesize*tid;
  offset=roffset;
  size_t loop=0;
  size_t numsteps=(m_readytype=='E')?getNumDPPSample():1;
  size_t outstep=getNoValues();
  size_t instep=m_left->getNoValues();
  switch (m_op)
  {
    case SWAP:
        for (loop=0;loop<numsteps;++loop)
        {
            escript::swapaxes(*leftres,m_left->getShape(),subroffset, m_samples_c, getShape(),offset, m_axis_offset, m_transpose);
            subroffset+=instep;
            offset+=outstep;
        }
        break;
    default:
        throw DataException("Programmer error - resolveNodeNP1OUT2P can not resolve operator "+opToString(m_op)+".");
  }
  return &m_samples_c;
}

const DataTypes::RealVectorType*
DataLazy::resolveNodeCondEval(int tid, int sampleNo, size_t& roffset) const
{
  if (m_readytype!='E')
  {
    throw DataException("Programmer error - resolveNodeCondEval should only be called on expanded Data.");
  }
  if (m_op!=CONDEVAL)
  {
    throw DataException("Programmer error - resolveNodeCondEval should only be called on CONDEVAL nodes.");
  }
  size_t subroffset;

  const RealVectorType* maskres=m_mask->resolveNodeSample(tid, sampleNo, subroffset);
  const RealVectorType* srcres=0;
  if ((*maskres)[subroffset]>0)
  {
        srcres=m_left->resolveNodeSample(tid, sampleNo, subroffset);
  }
  else
  {
        srcres=m_right->resolveNodeSample(tid, sampleNo, subroffset);
  }

  // Now we need to copy the result

  roffset=m_samplesize*tid;
  for (int i=0;i<m_samplesize;++i)
  {
        m_samples_r[roffset+i]=(*srcres)[subroffset+i];   
  }

  return &m_samples_r;
}

const DataTypes::CplxVectorType*
DataLazy::resolveNodeCondEvalCplx(int tid, int sampleNo, size_t& roffset) const
{
  if (m_readytype!='E')
  {
    throw DataException("Programmer error - resolveNodeCondEval should only be called on expanded Data.");
  }
  if (m_op!=CONDEVAL)
  {
    throw DataException("Programmer error - resolveNodeCondEval should only be called on CONDEVAL nodes.");
  }
  size_t subroffset;

  const RealVectorType* maskres=m_mask->resolveNodeSample(tid, sampleNo, subroffset);
  const CplxVectorType* srcres=0;
  if ((*maskres)[subroffset]>0)
  {
        srcres=m_left->resolveNodeSampleCplx(tid, sampleNo, subroffset);
  }
  else
  {
        srcres=m_right->resolveNodeSampleCplx(tid, sampleNo, subroffset);
  }

  // Now we need to copy the result

  roffset=m_samplesize*tid;
  for (int i=0;i<m_samplesize;++i)
  {
        m_samples_c[roffset+i]=(*srcres)[subroffset+i];   
  }

  return &m_samples_c;
}


// This method assumes that any subexpressions which evaluate to Constant or Tagged Data
// have already been collapsed to IDENTITY. So we must have at least one expanded child.
// If both children are expanded, then we can process them in a single operation (we treat
// the whole sample as one big datapoint.
// If one of the children is not expanded, then we need to treat each point in the sample
// individually.
// There is an additional complication when scalar operations are considered.
// For example, 2+Vector.
// In this case each double within the point is treated individually
const DataTypes::RealVectorType*
DataLazy::resolveNodeBinary(int tid, int sampleNo, size_t& roffset) const
{
LAZYDEBUG(cout << "Resolve binary: " << toString() << endl;)

  size_t lroffset=0, rroffset=0;        // offsets in the left and right result vectors
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
  size_t chunksize=1;                   // how many doubles will be processed in one go
  int leftstep=0;               // how far should the left offset advance after each step
  int rightstep=0;
  int numsteps=0;               // total number of steps for the inner loop
  int oleftstep=0;      // the o variables refer to the outer loop
  int orightstep=0;     // The outer loop is only required in cases where there is an extended scalar
  int onumsteps=1;
  
  bool LES=(leftExp && leftScalar);     // Left is an expanded scalar
  bool RES=(rightExp && rightScalar);
  bool LS=(!leftExp && leftScalar);     // left is a single scalar
  bool RS=(!rightExp && rightScalar);
  bool LN=(!leftExp && !leftScalar);    // left is a single non-scalar
  bool RN=(!rightExp && !rightScalar);
  bool LEN=(leftExp && !leftScalar);    // left is an expanded non-scalar
  bool REN=(rightExp && !rightScalar);

  if ((LES && RES) || (LEN && REN))     // both are Expanded scalars or both are expanded non-scalars
  {
        chunksize=m_left->getNumDPPSample()*leftsize;
        leftstep=0;
        rightstep=0;
        numsteps=1;
  }
  else if (LES || RES)
  {
        chunksize=1;
        if (LES)                // left is an expanded scalar
        {
                if (RS)
                {
                   leftstep=1;
                   rightstep=0;
                   numsteps=m_left->getNumDPPSample();
                }
                else            // RN or REN
                {
                   leftstep=0;
                   oleftstep=1;
                   rightstep=1;
                   orightstep=(RN ? -(int)rightsize : 0);
                   numsteps=rightsize;
                   onumsteps=m_left->getNumDPPSample();
                }
        }
        else            // right is an expanded scalar
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
                   oleftstep=(LN ? -(int)leftsize : 0);
                   numsteps=leftsize;
                   onumsteps=m_right->getNumDPPSample();
                }
        }
  }
  else  // this leaves (LEN, RS), (LEN, RN) and their transposes
  {
        if (LEN)        // and Right will be a single value 
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
        else    // REN
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

  int resultStep=max(leftstep,rightstep);       // only one (at most) should be !=0
        // Get the values of sub-expressions
  const RealVectorType* left=m_left->resolveNodeSample(tid,sampleNo,lroffset);       
  const RealVectorType* right=m_right->resolveNodeSample(tid,sampleNo,rroffset);
LAZYDEBUG(cout << "Post sub calls in " << toString() << endl;)
LAZYDEBUG(cout << "shapes=" << DataTypes::shapeToString(m_left->getShape()) << "," << DataTypes::shapeToString(m_right->getShape()) << endl;)
LAZYDEBUG(cout << "chunksize=" << chunksize << endl << "leftstep=" << leftstep << " rightstep=" << rightstep;)
LAZYDEBUG(cout << " numsteps=" << numsteps << endl << "oleftstep=" << oleftstep << " orightstep=" << orightstep;)
LAZYDEBUG(cout << "onumsteps=" << onumsteps << endl;)
LAZYDEBUG(cout << " DPPS=" << m_left->getNumDPPSample() << "," <<m_right->getNumDPPSample() << endl;)
LAZYDEBUG(cout << "" << LS << RS << LN << RN << LES << RES <<LEN << REN <<   endl;)

LAZYDEBUG(cout << "Left res["<< lroffset<< "]=" << (*left)[lroffset] << endl;)
LAZYDEBUG(cout << "Right res["<< rroffset<< "]=" << (*right)[rroffset] << endl;)


  roffset=m_samplesize*tid;
  double* resultp=&(m_samples_r[roffset]);                // results are stored at the vector offset we received
  switch(m_op)
  {
    case ADD:
        //PROC_OP(NO_ARG,plus<double>());
      escript::binaryOpVectorLazyArithmeticHelper<real_t, real_t, real_t>(resultp, 
			 &(*left)[0],
			 &(*right)[0],
			 chunksize,
			 onumsteps,
			 numsteps,
			 resultStep,
			 leftstep,
			 rightstep,
			 oleftstep,
			 orightstep,
			 lroffset,
			 rroffset,
			 escript::ES_optype::ADD);	
        break;
    case SUB:
      escript::binaryOpVectorLazyArithmeticHelper<real_t, real_t, real_t>(resultp, 
			 &(*left)[0],
			 &(*right)[0],
			 chunksize,
			 onumsteps,
			 numsteps,
			 resultStep,
			 leftstep,
			 rightstep,
			 oleftstep,
			 orightstep,
			 lroffset,
			 rroffset,
			 escript::ES_optype::SUB);	      
        //PROC_OP(NO_ARG,minus<double>());
        break;
    case MUL:
        //PROC_OP(NO_ARG,multiplies<double>());
      escript::binaryOpVectorLazyArithmeticHelper<real_t, real_t, real_t>(resultp, 
			 &(*left)[0],
			 &(*right)[0],
			 chunksize,
			 onumsteps,
			 numsteps,
			 resultStep,
			 leftstep,
			 rightstep,
			 oleftstep,
			 orightstep,
			 lroffset,
			 rroffset,
			 escript::ES_optype::MUL);	      
        break;
    case DIV:
        //PROC_OP(NO_ARG,divides<double>());
      escript::binaryOpVectorLazyArithmeticHelper<real_t, real_t, real_t>(resultp, 
			 &(*left)[0],
			 &(*right)[0],
			 chunksize,
			 onumsteps,
			 numsteps,
			 resultStep,
			 leftstep,
			 rightstep,
			 oleftstep,
			 orightstep,
			 lroffset,
			 rroffset,
			 escript::ES_optype::DIV);	      
        break;
    case POW:
       //PROC_OP(double (double,double),::pow);
      escript::binaryOpVectorLazyArithmeticHelper<real_t, real_t, real_t>(resultp, 
			 &(*left)[0],
			 &(*right)[0],
			 chunksize,
			 onumsteps,
			 numsteps,
			 resultStep,
			 leftstep,
			 rightstep,
			 oleftstep,
			 orightstep,
			 lroffset,
			 rroffset,
			 escript::ES_optype::POW);	      
        break;
    default:
        throw DataException("Programmer error - resolveBinary can not resolve operator "+opToString(m_op)+".");
  }
LAZYDEBUG(cout << "Result res[" << roffset<< "]" << m_samples_r[roffset] << endl;)
  return &m_samples_r;
}

const DataTypes::CplxVectorType*
DataLazy::resolveNodeBinaryCplx(int tid, int sampleNo, size_t& roffset) const
{
LAZYDEBUG(cout << "Resolve binary: " << toString() << endl;)

  size_t lroffset=0, rroffset=0;        // offsets in the left and right result vectors
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
  size_t chunksize=1;                   // how many doubles will be processed in one go
  int leftstep=0;               // how far should the left offset advance after each step
  int rightstep=0;
  int numsteps=0;               // total number of steps for the inner loop
  int oleftstep=0;      // the o variables refer to the outer loop
  int orightstep=0;     // The outer loop is only required in cases where there is an extended scalar
  int onumsteps=1;
  
  bool LES=(leftExp && leftScalar);     // Left is an expanded scalar
  bool RES=(rightExp && rightScalar);
  bool LS=(!leftExp && leftScalar);     // left is a single scalar
  bool RS=(!rightExp && rightScalar);
  bool LN=(!leftExp && !leftScalar);    // left is a single non-scalar
  bool RN=(!rightExp && !rightScalar);
  bool LEN=(leftExp && !leftScalar);    // left is an expanded non-scalar
  bool REN=(rightExp && !rightScalar);

  if ((LES && RES) || (LEN && REN))     // both are Expanded scalars or both are expanded non-scalars
  {
        chunksize=m_left->getNumDPPSample()*leftsize;
        leftstep=0;
        rightstep=0;
        numsteps=1;
  }
  else if (LES || RES)
  {
        chunksize=1;
        if (LES)                // left is an expanded scalar
        {
                if (RS)
                {
                   leftstep=1;
                   rightstep=0;
                   numsteps=m_left->getNumDPPSample();
                }
                else            // RN or REN
                {
                   leftstep=0;
                   oleftstep=1;
                   rightstep=1;
                   orightstep=(RN ? -(int)rightsize : 0);
                   numsteps=rightsize;
                   onumsteps=m_left->getNumDPPSample();
                }
        }
        else            // right is an expanded scalar
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
                   oleftstep=(LN ? -(int)leftsize : 0);
                   numsteps=leftsize;
                   onumsteps=m_right->getNumDPPSample();
                }
        }
  }
  else  // this leaves (LEN, RS), (LEN, RN) and their transposes
  {
        if (LEN)        // and Right will be a single value 
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
        else    // REN
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

  int resultStep=max(leftstep,rightstep);       // only one (at most) should be !=0
        // Get the values of sub-expressions
  const CplxVectorType* left=m_left->resolveNodeSampleCplx(tid,sampleNo,lroffset);       
  const CplxVectorType* right=m_right->resolveNodeSampleCplx(tid,sampleNo,rroffset);
LAZYDEBUG(cout << "Post sub calls in " << toString() << endl;)
LAZYDEBUG(cout << "shapes=" << DataTypes::shapeToString(m_left->getShape()) << "," << DataTypes::shapeToString(m_right->getShape()) << endl;)
LAZYDEBUG(cout << "chunksize=" << chunksize << endl << "leftstep=" << leftstep << " rightstep=" << rightstep;)
LAZYDEBUG(cout << " numsteps=" << numsteps << endl << "oleftstep=" << oleftstep << " orightstep=" << orightstep;)
LAZYDEBUG(cout << "onumsteps=" << onumsteps << endl;)
LAZYDEBUG(cout << " DPPS=" << m_left->getNumDPPSample() << "," <<m_right->getNumDPPSample() << endl;)
LAZYDEBUG(cout << "" << LS << RS << LN << RN << LES << RES <<LEN << REN <<   endl;)

LAZYDEBUG(cout << "Left res["<< lroffset<< "]=" << (*left)[lroffset] << endl;)
LAZYDEBUG(cout << "Right res["<< rroffset<< "]=" << (*right)[rroffset] << endl;)


  roffset=m_samplesize*tid;
  cplx_t* resultp=&(m_samples_c[roffset]);                // results are stored at the vector offset we received
  switch(m_op)
  {
    case ADD:
        //PROC_OP(NO_ARG,plus<double>());
      escript::binaryOpVectorLazyArithmeticHelper<cplx_t, cplx_t, cplx_t>(resultp, 
			 &(*left)[0],
			 &(*right)[0],
			 chunksize,
			 onumsteps,
			 numsteps,
			 resultStep,
			 leftstep,
			 rightstep,
			 oleftstep,
			 orightstep,
			 lroffset,
			 rroffset,
			 escript::ES_optype::ADD);	
        break;
    case SUB:
      escript::binaryOpVectorLazyArithmeticHelper<cplx_t, cplx_t, cplx_t>(resultp, 
			 &(*left)[0],
			 &(*right)[0],
			 chunksize,
			 onumsteps,
			 numsteps,
			 resultStep,
			 leftstep,
			 rightstep,
			 oleftstep,
			 orightstep,
			 lroffset,
			 rroffset,
			 escript::ES_optype::SUB);	      
        //PROC_OP(NO_ARG,minus<double>());
        break;
    case MUL:
        //PROC_OP(NO_ARG,multiplies<double>());
      escript::binaryOpVectorLazyArithmeticHelper<cplx_t, cplx_t, cplx_t>(resultp, 
			 &(*left)[0],
			 &(*right)[0],
			 chunksize,
			 onumsteps,
			 numsteps,
			 resultStep,
			 leftstep,
			 rightstep,
			 oleftstep,
			 orightstep,
			 lroffset,
			 rroffset,
			 escript::ES_optype::MUL);	      
        break;
    case DIV:
        //PROC_OP(NO_ARG,divides<double>());
      escript::binaryOpVectorLazyArithmeticHelper<cplx_t, cplx_t, cplx_t>(resultp, 
			 &(*left)[0],
			 &(*right)[0],
			 chunksize,
			 onumsteps,
			 numsteps,
			 resultStep,
			 leftstep,
			 rightstep,
			 oleftstep,
			 orightstep,
			 lroffset,
			 rroffset,
			 escript::ES_optype::DIV);	      
        break;
    case POW:
       //PROC_OP(double (double,double),::pow);
      escript::binaryOpVectorLazyArithmeticHelper<cplx_t, cplx_t, cplx_t>(resultp, 
			 &(*left)[0],
			 &(*right)[0],
			 chunksize,
			 onumsteps,
			 numsteps,
			 resultStep,
			 leftstep,
			 rightstep,
			 oleftstep,
			 orightstep,
			 lroffset,
			 rroffset,
			 escript::ES_optype::POW);	      
        break;
    default:
        throw DataException("Programmer error - resolveBinary can not resolve operator "+opToString(m_op)+".");
  }
LAZYDEBUG(cout << "Result res[" << roffset<< "]" << m_samples_c[roffset] << endl;)
  return &m_samples_c;
}


// This method assumes that any subexpressions which evaluate to Constant or Tagged Data
// have already been collapsed to IDENTITY. So we must have at least one expanded child.
// unlike the other resolve helpers, we must treat these datapoints separately.
const DataTypes::RealVectorType*
DataLazy::resolveNodeTProd(int tid, int sampleNo, size_t& roffset) const
{
LAZYDEBUG(cout << "Resolve TensorProduct: " << toString() << endl;)
  size_t lroffset=0, rroffset=0;        // offsets in the left and right result vectors
        // first work out which of the children are expanded
  bool leftExp=(m_left->m_readytype=='E');
  bool rightExp=(m_right->m_readytype=='E');
  int steps=getNumDPPSample();
  int leftStep=(leftExp? m_left->getNoValues() : 0);            // do not have scalars as input to this method
  int rightStep=(rightExp?m_right->getNoValues() : 0);

  int resultStep=getNoValues();
  roffset=m_samplesize*tid;
  size_t offset=roffset;
  const RealVectorType* left=m_left->resolveNodeSample(tid, sampleNo, lroffset);
  const RealVectorType* right=m_right->resolveNodeSample(tid, sampleNo, rroffset);
LAZYDEBUG(cerr << "[Left shape]=" << DataTypes::shapeToString(m_left->getShape()) << "\n[Right shape]=" << DataTypes::shapeToString(m_right->getShape()) << " result=" <<DataTypes::shapeToString(getShape()) <<  endl;
cout << getNoValues() << endl;)


LAZYDEBUG(cerr << "Post sub calls: " << toString() << endl;)
LAZYDEBUG(cout << "LeftExp=" << leftExp << " rightExp=" << rightExp << endl;)
LAZYDEBUG(cout << "LeftR=" << m_left->getRank() << " rightExp=" << m_right->getRank() << endl;)
LAZYDEBUG(cout << "LeftSize=" << m_left->getNoValues() << " RightSize=" << m_right->getNoValues() << endl;)
LAZYDEBUG(cout << "m_samplesize=" << m_samplesize << endl;)
LAZYDEBUG(cout << "outputshape=" << DataTypes::shapeToString(getShape()) << endl;)
LAZYDEBUG(cout << "DPPS=" << m_right->getNumDPPSample() <<"."<<endl;)

  double* resultp=&(m_samples_r[offset]);         // results are stored at the vector offset we received
  switch(m_op)
  {
    case PROD:
        for (int i=0;i<steps;++i,resultp+=resultStep)
        {
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
  return &m_samples_r;
}

const DataTypes::CplxVectorType*
DataLazy::resolveNodeTProdCplx(int tid, int sampleNo, size_t& roffset) const
{
LAZYDEBUG(cout << "Resolve TensorProduct: " << toString() << endl;)

  size_t lroffset=0, rroffset=0;        // offsets in the left and right result vectors
        // first work out which of the children are expanded
  bool leftExp=(m_left->m_readytype=='E');
  bool rightExp=(m_right->m_readytype=='E');
  int steps=getNumDPPSample();
  int leftStep=(leftExp? m_left->getNoValues() : 0);            // do not have scalars as input to this method
  int rightStep=(rightExp?m_right->getNoValues() : 0);

  int resultStep=getNoValues();
  roffset=m_samplesize*tid;
  size_t offset=roffset;

  const CplxVectorType* left=m_left->resolveNodeSampleCplx(tid, sampleNo, lroffset);

  const CplxVectorType* right=m_right->resolveNodeSampleCplx(tid, sampleNo, rroffset);

LAZYDEBUG(cerr << "[Left shape]=" << DataTypes::shapeToString(m_left->getShape()) << "\n[Right shape]=" << DataTypes::shapeToString(m_right->getShape()) << " result=" <<DataTypes::shapeToString(getShape()) <<  endl;
cout << getNoValues() << endl;)


LAZYDEBUG(cerr << "Post sub calls: " << toString() << endl;)
LAZYDEBUG(cout << "LeftExp=" << leftExp << " rightExp=" << rightExp << endl;)
LAZYDEBUG(cout << "LeftR=" << m_left->getRank() << " rightExp=" << m_right->getRank() << endl;)
LAZYDEBUG(cout << "LeftSize=" << m_left->getNoValues() << " RightSize=" << m_right->getNoValues() << endl;)
LAZYDEBUG(cout << "m_samplesize=" << m_samplesize << endl;)
LAZYDEBUG(cout << "outputshape=" << DataTypes::shapeToString(getShape()) << endl;)
LAZYDEBUG(cout << "DPPS=" << m_right->getNumDPPSample() <<"."<<endl;)

  cplx_t* resultp=&(m_samples_c[offset]);         // results are stored at the vector offset we received
  switch(m_op)
  {
    case PROD:
        for (int i=0;i<steps;++i,resultp+=resultStep)
        {
          const cplx_t *ptr_0 = &((*left)[lroffset]);
          const cplx_t *ptr_1 = &((*right)[rroffset]);

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
   return &m_samples_c;
}

const DataTypes::RealVectorType*
DataLazy::resolveSample(int sampleNo, size_t& roffset) const
{
#ifdef _OPENMP
        int tid=omp_get_thread_num();
#else
        int tid=0;
#endif 

#ifdef LAZY_STACK_PROF
        stackstart[tid]=&tid;
        stackend[tid]=&tid;
        const DataTypes::RealVectorType* r=resolveNodeSample(tid, sampleNo, roffset);
        size_t d=(size_t)stackstart[tid]-(size_t)stackend[tid];
        #pragma omp critical
        if (d>maxstackuse)
        {
cout << "Max resolve Stack use " << d << endl;
                maxstackuse=d;
        }
        return r;
#else
        return resolveNodeSample(tid, sampleNo, roffset);
#endif
}


const DataTypes::RealVectorType*
DataLazy::resolveTypedSample(int sampleNo, size_t& roffset, DataTypes::real_t dummy) const
{
#ifdef _OPENMP
        int tid=omp_get_thread_num();
#else
        int tid=0;
#endif 

#ifdef LAZY_STACK_PROF
        stackstart[tid]=&tid;
        stackend[tid]=&tid;
        const DataTypes::RealVectorType* r=resolveNodeSample(tid, sampleNo, roffset);
        size_t d=(size_t)stackstart[tid]-(size_t)stackend[tid];
        #pragma omp critical
        if (d>maxstackuse)
        {
cout << "Max resolve Stack use " << d << endl;
                maxstackuse=d;
        }
        return r;
#else
        return resolveNodeSample(tid, sampleNo, roffset);
#endif
}

const DataTypes::CplxVectorType*
DataLazy::resolveTypedSample(int sampleNo, size_t& roffset, DataTypes::cplx_t dummy) const
{
#ifdef _OPENMP
        int tid=omp_get_thread_num();
#else
        int tid=0;
#endif 

#ifdef LAZY_STACK_PROF
        stackstart[tid]=&tid;
        stackend[tid]=&tid;
        const DataTypes::CplxVectorType* r=resolveNodeSampleCplx(tid, sampleNo, roffset);
        size_t d=(size_t)stackstart[tid]-(size_t)stackend[tid];
        #pragma omp critical
        if (d>maxstackuse)
        {
cout << "Max resolve Stack use " << d << endl;
                maxstackuse=d;
        }
        return r;
#else
        return resolveNodeSampleCplx(tid, sampleNo, roffset);
#endif
}



// This needs to do the work of the identity constructor
void
DataLazy::resolveToIdentity()
{
   if (m_op==IDENTITY)
        return;
   if (isComplex())
   {
        DataReady_ptr p=resolveNodeWorkerCplx();
        makeIdentity(p);
   }
   else
   {
        DataReady_ptr p=resolveNodeWorker();
        makeIdentity(p);
   }
}

void DataLazy::makeIdentity(const DataReady_ptr& p)
{
   m_axis_offset=0;
   m_transpose=0;
   m_SL=m_SM=m_SR=0;
   m_children=m_height=0;
   m_id=p;
   if(p->isConstant()) {m_readytype='C';}
   else if(p->isExpanded()) {m_readytype='E';}
   else if (p->isTagged()) {m_readytype='T';}
   else {throw DataException("Unknown DataReady instance in convertToIdentity constructor.");}
   m_samplesize=p->getNumDPPSample()*p->getNoValues();
   m_left.reset();
   m_right.reset();
   m_iscompl=p->isComplex();
   m_op=IDENTITY;
   m_opgroup=getOpgroup(m_op);
}


DataReady_ptr
DataLazy::resolve()
{
    resolveToIdentity();
    return m_id;
}


/* This is really a static method but I think that caused problems in windows */
void
DataLazy::resolveGroupWorker(std::vector<DataLazy*>& dats)
{
  if (dats.empty())
  {
        return;
  }
  vector<DataLazy*> work;
  FunctionSpace fs=dats[0]->getFunctionSpace();
  bool match=true;
  for (int i=dats.size()-1;i>=0;--i)
  {
        if (dats[i]->m_readytype!='E')
        {
                dats[i]->collapse();
        }
        if (dats[i]->m_op!=IDENTITY)
        {
                work.push_back(dats[i]);
                if (fs!=dats[i]->getFunctionSpace())
                {
                        match=false;
                }
        }
  }
  if (work.empty())
  {
        return;         // no work to do
  }
  if (match)    // all functionspaces match.  Yes I realise this is overly strict
  {             // it is possible that dats[0] is one of the objects which we discarded and
                // all the other functionspaces match.
        vector<DataExpanded*> dep;
        vector<RealVectorType*> vecs;
        for (int i=0;i<work.size();++i)
        {
                dep.push_back(new DataExpanded(fs,work[i]->getShape(), RealVectorType(work[i]->getNoValues())));
                vecs.push_back(&(dep[i]->getVectorRW()));
        }
        int totalsamples=work[0]->getNumSamples();
        const RealVectorType* res=0; // Storage for answer
        int sample;
        #pragma omp parallel private(sample, res)
        {
            size_t roffset=0;
            #pragma omp for schedule(static)
            for (sample=0;sample<totalsamples;++sample)
            {
                roffset=0;
                int j;
                for (j=work.size()-1;j>=0;--j)
                {
#ifdef _OPENMP
                    res=work[j]->resolveNodeSample(omp_get_thread_num(),sample,roffset);
#else
                    res=work[j]->resolveNodeSample(0,sample,roffset);
#endif
                    RealVectorType::size_type outoffset=dep[j]->getPointOffset(sample,0);
                    memcpy(&((*vecs[j])[outoffset]),&((*res)[roffset]),work[j]->m_samplesize*sizeof(RealVectorType::ElementType));
                }
            }
        }
        // Now we need to load the new results as identity ops into the lazy nodes
        for (int i=work.size()-1;i>=0;--i)
        {
            work[i]->makeIdentity(REFCOUNTNS::dynamic_pointer_cast<DataReady>(dep[i]->getPtr()));
        }
  }
  else  // functionspaces do not match
  {
        for (int i=0;i<work.size();++i)
        {
                work[i]->resolveToIdentity();
        }
  }
}



// This version of resolve uses storage in each node to hold results
DataReady_ptr
DataLazy::resolveNodeWorker()
{
  if (m_readytype!='E')         // if the whole sub-expression is Constant or Tagged, then evaluate it normally
  {
    collapse();
  }
  if (m_op==IDENTITY)           // So a lazy expression of Constant or Tagged data will be returned here. 
  {
    return m_id;
  }
        // from this point on we must have m_op!=IDENTITY and m_readytype=='E'
  DataExpanded* result=new DataExpanded(getFunctionSpace(),getShape(),  RealVectorType(getNoValues()));
  RealVectorType& resvec=result->getVectorRW();
  DataReady_ptr resptr=DataReady_ptr(result);

  int sample;
  int totalsamples=getNumSamples();
  const RealVectorType* res=0;       // Storage for answer
LAZYDEBUG(cout << "Total number of samples=" <<totalsamples << endl;)
  #pragma omp parallel private(sample,res)
  {
        size_t roffset=0;
#ifdef LAZY_STACK_PROF
        stackstart[omp_get_thread_num()]=&roffset;
        stackend[omp_get_thread_num()]=&roffset;
#endif
        #pragma omp for schedule(static)
        for (sample=0;sample<totalsamples;++sample)
        {
                roffset=0;
#ifdef _OPENMP
                res=resolveNodeSample(omp_get_thread_num(),sample,roffset);
#else
                res=resolveNodeSample(0,sample,roffset);
#endif
LAZYDEBUG(cout << "Sample #" << sample << endl;)
LAZYDEBUG(cout << "Final res[" << roffset<< "]=" << (*res)[roffset] << (*res)[roffset]<< endl; )
                RealVectorType::size_type outoffset=result->getPointOffset(sample,0);
                memcpy(&(resvec[outoffset]),&((*res)[roffset]),m_samplesize*sizeof(RealVectorType::ElementType));
        }
  }
#ifdef LAZY_STACK_PROF
  for (int i=0;i<getNumberOfThreads();++i)
  {
        size_t r=((size_t)stackstart[i] - (size_t)stackend[i]);
//      cout << i << " " << stackstart[i] << " .. " << stackend[i] << " = " <<  r << endl;
        if (r>maxstackuse)
        {
                maxstackuse=r;
        }
  }
  cout << "Max resolve Stack use=" << maxstackuse << endl;
#endif
  return resptr;
}

// This version should only be called on complex lazy nodes
DataReady_ptr
DataLazy::resolveNodeWorkerCplx()
{
  if (m_readytype!='E')         // if the whole sub-expression is Constant or Tagged, then evaluate it normally
  {
    collapse();
  }
  if (m_op==IDENTITY)           // So a lazy expression of Constant or Tagged data will be returned here. 
  {
    return m_id;
  }
        // from this point on we must have m_op!=IDENTITY and m_readytype=='E'
  DataExpanded* result=new DataExpanded(getFunctionSpace(),getShape(),  CplxVectorType(getNoValues()));
  CplxVectorType& resvec=result->getVectorRWC();
  DataReady_ptr resptr=DataReady_ptr(result);

  int sample;
  int totalsamples=getNumSamples();
  const CplxVectorType* res=0;       // Storage for answer
LAZYDEBUG(cout << "Total number of samples=" <<totalsamples << endl;)
  #pragma omp parallel private(sample,res)
  {
        size_t roffset=0;
#ifdef LAZY_STACK_PROF
        stackstart[omp_get_thread_num()]=&roffset;
        stackend[omp_get_thread_num()]=&roffset;
#endif
        #pragma omp for schedule(static)
        for (sample=0;sample<totalsamples;++sample)
        {
                roffset=0;
#ifdef _OPENMP
                res=resolveNodeSampleCplx(omp_get_thread_num(),sample,roffset);
#else
                res=resolveNodeSampleCplx(0,sample,roffset);
#endif
LAZYDEBUG(cout << "Sample #" << sample << endl;)
LAZYDEBUG(cout << "Final res[" << roffset<< "]=" << (*res)[roffset] << (*res)[roffset]<< endl; )
                CplxVectorType::size_type outoffset=result->getPointOffset(sample,0);
                memcpy(&(resvec[outoffset]),&((*res)[roffset]),m_samplesize*sizeof(CplxVectorType::ElementType));
        }
  }
#ifdef LAZY_STACK_PROF
  for (int i=0;i<getNumberOfThreads();++i)
  {
        size_t r=((size_t)stackstart[i] - (size_t)stackend[i]);
//      cout << i << " " << stackstart[i] << " .. " << stackend[i] << " = " <<  r << endl;
        if (r>maxstackuse)
        {
                maxstackuse=r;
        }
  }
  cout << "Max resolve Stack use=" << maxstackuse << endl;
#endif
  return resptr;
}


std::string
DataLazy::toString() const
{
  ostringstream oss;
  oss << "Lazy Data: [depth=" << m_height<< "] ";
  switch (escriptParams.getLazyStrFmt())
  {
  case 1:       // tree format
        oss << endl;
        intoTreeString(oss,""); 
        break;
  case 2:       // just the depth
        break;
  default:
        intoString(oss);
        break;
  }
  return oss.str();
}


void
DataLazy::intoString(ostringstream& oss) const
{
//    oss << "[" << m_children <<";"<<m_height <<"]";
  switch (m_opgroup)
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
        if (m_id->isComplex())
        {
          oss << "j";  
        }
        oss << '@' << m_id.get();
        break;
  case G_BINARY:
        oss << '(';
        m_left->intoString(oss);
        oss << ' ' << opToString(m_op) << ' ';
        m_right->intoString(oss);
        oss << ')';
        if (isComplex())
        {
            oss << "j";
        }            
        break;
  case G_UNARY:
  case G_UNARY_P:
  case G_UNARY_PR:
  case G_UNARY_R:
  case G_NP1OUT:
  case G_NP1OUT_P:
  case G_REDUCTION:
        oss << opToString(m_op) << '(';
        m_left->intoString(oss);
        oss << ')';
        if (isComplex())
        {
            oss << "j";
        }                    
        break;
  case G_TENSORPROD:
        oss << opToString(m_op) << '(';
        m_left->intoString(oss);
        oss << ", ";
        m_right->intoString(oss);
        oss << ')'; 
        if (isComplex())
        {
            oss << "j";
        }                    
        break;
  case G_NP1OUT_2P:
        oss << opToString(m_op) << '(';
        m_left->intoString(oss);
        oss << ", " << m_axis_offset << ", " << m_transpose;
        oss << ')';
        if (isComplex())
        {
            oss << "j";
        }                    
        break;
  case G_CONDEVAL:
        oss << opToString(m_op)<< '(' ;
        m_mask->intoString(oss);
        oss << " ? ";
        m_left->intoString(oss);
        oss << " : ";
        m_right->intoString(oss); 
        oss << ')';
        if (isComplex())
        {
            oss << "j";
        }                    
        break;
  case G_UNARY_C:
        oss << opToString(m_op) <<'(';
        m_left->intoString(oss);
        oss << ')';
        if (isComplex())
        {
            oss << "j";
        }                    
        break;
  default:
        oss << "UNKNOWN";
        if (isComplex())
        {
            oss << "j";
        }                    
  }
}


void
DataLazy::intoTreeString(ostringstream& oss, string indent) const
{
  oss << '[' << m_rank << ':' << setw(3) << m_samplesize << "] " << indent;
  ES_opgroup gop=m_opgroup;
  switch (gop)
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
        if (isComplex())
        {
            oss << "j";
        }
        oss << '@' << m_id.get() << endl;
        break;
  case G_BINARY:
        oss << opToString(m_op);
        if (isComplex())
        {
            oss << 'j';
        }
        oss << endl;
        indent+='.';
        m_left->intoTreeString(oss, indent);
        m_right->intoTreeString(oss, indent);
        break;
  case G_UNARY:
  case G_UNARY_P:
  case G_UNARY_PR:      
  case G_UNARY_R:      
  case G_NP1OUT:
  case G_NP1OUT_P:
  case G_REDUCTION:
  case G_UNARY_C:
        oss << opToString(m_op);
        if (isComplex())
        {
            oss << 'j';
        }
        oss << endl;
        indent+='.';
        m_left->intoTreeString(oss, indent);
        break;
  case G_TENSORPROD:
        oss << opToString(m_op);
        if (isComplex())
        {
            oss << 'j';
        }
        oss << endl;
        indent+='.';
        m_left->intoTreeString(oss, indent);
        m_right->intoTreeString(oss, indent);
        break;
  case G_NP1OUT_2P:
        oss << opToString(m_op);
        if (isComplex())
        {
            oss << 'j';
        }
        oss << ", " << m_axis_offset << ", " << m_transpose<< endl;
        indent+='.';
        m_left->intoTreeString(oss, indent);
        break;
  default:
        oss << "UNKNOWN";
        if (isComplex())
        {
            oss << 'j';
        }
        oss << endl;                
  }
}


DataAbstract* 
DataLazy::deepCopy() const
{
  switch (m_opgroup)
  {
  case G_IDENTITY:  return new DataLazy(m_id->deepCopy()->getPtr());
  case G_UNARY: 
  case G_UNARY_R:
  case G_UNARY_C:
  case G_REDUCTION:      return new DataLazy(m_left->deepCopy()->getPtr(),m_op);
  case G_UNARY_PR:
  case G_UNARY_P:       return new DataLazy(m_left->deepCopy()->getPtr(), m_op, m_tol);
  case G_BINARY:        return new DataLazy(m_left->deepCopy()->getPtr(),m_right->deepCopy()->getPtr(),m_op);
  case G_NP1OUT: return new DataLazy(m_left->deepCopy()->getPtr(), m_right->deepCopy()->getPtr(),m_op);
  case G_TENSORPROD: return new DataLazy(m_left->deepCopy()->getPtr(), m_right->deepCopy()->getPtr(), m_op, m_axis_offset, m_transpose);
  case G_NP1OUT_P:   return new DataLazy(m_left->deepCopy()->getPtr(),m_op,  m_axis_offset);
  case G_NP1OUT_2P:  return new DataLazy(m_left->deepCopy()->getPtr(), m_op, m_axis_offset, m_transpose);
  case G_CONDEVAL:  return new DataLazy(m_mask->deepCopy()->getPtr(), m_left->deepCopy()->getPtr(), m_right->deepCopy()->getPtr());
  default:
        throw DataException("Programmer error - do not know how to deepcopy operator "+opToString(m_op)+".");
  }
}

// For this, we don't care what op you were doing because the answer is now zero
DataAbstract* 
DataLazy::zeroedCopy() const
{
  return new DataLazy(m_id->zeroedCopy()->getPtr());
}

// There is no single, natural interpretation of getLength on DataLazy.
// Instances of DataReady can look at the size of their vectors.
// For lazy though, it could be the size the data would be if it were resolved;
// or it could be some function of the lengths of the DataReady instances which 
// form part of the expression.
// Rather than have people making assumptions, I have disabled the method.
DataTypes::RealVectorType::size_type
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
DataTypes::RealVectorType::size_type 
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
DataTypes::RealVectorType::size_type 
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


// I have decided to let Data:: handle this issue.
void
DataLazy::setToZero()
{
//   DataTypes::RealVectorType v(getNoValues(),0);
//   m_id=DataReady_ptr(new DataConstant(getFunctionSpace(),getShape(),v));
//   m_op=IDENTITY;
//   m_right.reset();   
//   m_left.reset();
//   m_readytype='C';
//   m_buffsRequired=1;

  (void)privdebug;  // to stop the compiler complaining about unused privdebug
  throw DataException("Programmer error - setToZero not supported for DataLazy (DataLazy objects should be read only).");
}

bool
DataLazy::actsExpanded() const
{
        return (m_readytype=='E');
}

} // end namespace

