
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


#include "Data.h"

#include "DataExpanded.h"
#include "DataConstant.h"
#include "DataTagged.h"
#include "DataEmpty.h"
#include "DataLazy.h"
#include "FunctionSpaceFactory.h"
#include "AbstractContinuousDomain.h"
#include "UnaryFuncs.h"
#include "FunctionSpaceException.h"
#include "EscriptParams.h"

extern "C" {
#include "esysUtils/blocktimer.h"
}

#include <fstream>
#include <algorithm>
#include <vector>
#include <functional>
#include <sstream>	// so we can throw messages about ranks

#include <boost/python/dict.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/long.hpp>
#include "WrappedArray.h"

using namespace std;
using namespace boost::python;
using namespace boost;
using namespace escript;

// ensure the current object is not a DataLazy
// The idea was that we could add an optional warning whenever a resolve is forced
// #define forceResolve() if (isLazy()) {#resolve();}

#define AUTOLAZYON escriptParams.getAUTOLAZY()
#define MAKELAZYOP(X)   if (isLazy() || (AUTOLAZYON && m_data->isExpanded())) \
  {\
	DataLazy* c=new DataLazy(borrowDataPtr(),X);\
	return Data(c);\
  }
#define MAKELAZYOPOFF(X,Y) if (isLazy() || (AUTOLAZYON && m_data->isExpanded())) \
  {\
	DataLazy* c=new DataLazy(borrowDataPtr(),X,Y);\
	return Data(c);\
  }

#define MAKELAZYOP2(X,Y,Z) if (isLazy() || (AUTOLAZYON && m_data->isExpanded())) \
  {\
	DataLazy* c=new DataLazy(borrowDataPtr(),X,Y,Z);\
	return Data(c);\
  }

#define MAKELAZYBINSELF(R,X)   if (isLazy() || R.isLazy() || (AUTOLAZYON && (isExpanded() || R.isExpanded()))) \
  {\
	DataLazy* c=new DataLazy(m_data,R.borrowDataPtr(),X);\
/*         m_data=c->getPtr();*/     set_m_data(c->getPtr());\
	return (*this);\
  }

// like the above but returns a new data rather than *this
#define MAKELAZYBIN(R,X)   if (isLazy() || R.isLazy() || (AUTOLAZYON && (isExpanded() || R.isExpanded()))) \
  {\
	DataLazy* c=new DataLazy(m_data,R.borrowDataPtr(),X);\
	return Data(c);\
  }

#define MAKELAZYBIN2(L,R,X)   if (L.isLazy() || R.isLazy() || (AUTOLAZYON && (L.isExpanded() || R.isExpanded()))) \
  {\
	DataLazy* c=new DataLazy(L.borrowDataPtr(),R.borrowDataPtr(),X);\
	return Data(c);\
  }

#define CHECK_DO_CRES escriptParams.getRESOLVE_COLLECTIVE()

namespace
{

template <class ARR>
inline
boost::python::tuple
pointToTuple1(const DataTypes::ShapeType& shape, ARR v, unsigned long offset)
{
    using namespace boost::python;
    using boost::python::tuple;
    using boost::python::list;

    list l;
    unsigned int dim=shape[0];
    for (size_t i=0;i<dim;++i)
    {
	l.append(v[i+offset]);
    }
    return tuple(l);
}

template <class ARR>
inline
boost::python::tuple
pointToTuple2(const DataTypes::ShapeType& shape, ARR v, unsigned long offset)
{
    using namespace boost::python;
    using boost::python::tuple;
    using boost::python::list;

    unsigned int shape0=shape[0];
    unsigned int shape1=shape[1];
    list lj;
    for (size_t j=0;j<shape0;++j)
    {
        list li;
	for (size_t i=0;i<shape1;++i)
	{
	    li.append(v[offset+DataTypes::getRelIndex(shape,j,i)]);
	}
	lj.append(tuple(li));
    }
    return tuple(lj);
}

template <class ARR>
inline
boost::python::tuple
pointToTuple3(const DataTypes::ShapeType& shape, ARR v, unsigned long offset)
{
    using namespace boost::python;
    using boost::python::tuple;
    using boost::python::list;

    unsigned int shape0=shape[0];
    unsigned int shape1=shape[1];
    unsigned int shape2=shape[2];

    list lk;
    for (size_t k=0;k<shape0;++k)
    {
        list lj;
	for (size_t j=0;j<shape1;++j)
	{
	    list li;
	    for (size_t i=0;i<shape2;++i)
	    {
                li.append(v[offset+DataTypes::getRelIndex(shape,k,j,i)]);
            }
	    lj.append(tuple(li));
        }
        lk.append(tuple(lj));
    }
    return tuple(lk);
}

template <class ARR>
inline
boost::python::tuple
pointToTuple4(const DataTypes::ShapeType& shape, ARR v, unsigned long offset)
{
    using namespace boost::python;
    using boost::python::tuple;
    using boost::python::list;

    unsigned int shape0=shape[0];
    unsigned int shape1=shape[1];
    unsigned int shape2=shape[2];
    unsigned int shape3=shape[3];

    list ll;
    for (size_t l=0;l<shape0;++l)
    {
        list lk;
	for (size_t k=0;k<shape1;++k)
	{
            list lj;
            for (size_t j=0;j<shape2;++j)
            {
                list li;
                for (size_t i=0;i<shape3;++i)
                {
                    li.append(v[offset+DataTypes::getRelIndex(shape,l,k,j,i)]);
                }
                lj.append(tuple(li));
            }
            lk.append(tuple(lj));
	}
        ll.append(tuple(lk));
    }
    return tuple(ll);
}


// This should be safer once the DataC RO changes have been brought in
template <class ARR>
boost::python::tuple
pointToTuple( const DataTypes::ShapeType& shape,ARR v)
{
   using namespace boost::python;
   using boost::python::list;
   int rank=shape.size();
   if (rank==0)
   {
	return make_tuple(v[0]);
   }
   else if (rank==1)
   {
        return pointToTuple1(shape,v,0);
   }
   else if (rank==2)
   {
	return pointToTuple2(shape,v,0);
   }
   else if (rank==3)
   {
	return pointToTuple3(shape,v,0);
   }
   else if (rank==4)
   {
	return pointToTuple4(shape,v,0);
   }
   else
     throw DataException("Unknown rank in pointToTuple.");
}

}  // anonymous namespace

Data::Data()
	: m_shared(false), m_lazy(false)
{
  //
  // Default data is type DataEmpty
  DataAbstract* temp=new DataEmpty();
//   m_data=temp->getPtr();
  set_m_data(temp->getPtr());
  m_protected=false;
}

Data::Data(double value,
           const tuple& shape,
           const FunctionSpace& what,
           bool expanded)
	: m_shared(false), m_lazy(false)
{
  DataTypes::ShapeType dataPointShape;
  for (int i = 0; i < shape.attr("__len__")(); ++i) {
    dataPointShape.push_back(extract<const int>(shape[i]));
  }

  int len = DataTypes::noValues(dataPointShape);
  DataVector temp_data(len,value,len);
  initialise(temp_data, dataPointShape, what, expanded);
  m_protected=false;
}

Data::Data(double value,
	   const DataTypes::ShapeType& dataPointShape,
	   const FunctionSpace& what,
           bool expanded)
	: m_shared(false), m_lazy(false)
{
  int len = DataTypes::noValues(dataPointShape);
  DataVector temp_data(len,value,len);
  initialise(temp_data, dataPointShape, what, expanded);
  m_protected=false;
}

Data::Data(const Data& inData)
	: m_shared(false), m_lazy(false)
{
  set_m_data(inData.m_data);
  m_protected=inData.isProtected();
}


Data::Data(const Data& inData,
           const DataTypes::RegionType& region)
	: m_shared(false), m_lazy(false)
{
  DataAbstract_ptr dat=inData.m_data;
  if (inData.isLazy())
  {
	dat=inData.m_data->resolve();
  }
  else
  {
	dat=inData.m_data;
  }
  //
  // Create Data which is a slice of another Data
  DataAbstract* tmp = dat->getSlice(region);
  set_m_data(DataAbstract_ptr(tmp));
  m_protected=false;

}

Data::Data(const Data& inData,
           const FunctionSpace& functionspace)
	: m_shared(false), m_lazy(false)
{
  if (inData.isEmpty())
  {
    throw DataException("Error - will not interpolate for instances of DataEmpty.");
  }
  if (inData.getFunctionSpace()==functionspace) {
    set_m_data(inData.m_data);
  } 
  else 
  {

    if (inData.isConstant()) {	// for a constant function, we just need to use the new function space
      if (!inData.probeInterpolation(functionspace))
      {           // Even though this is constant, we still need to check whether interpolation is allowed
	throw FunctionSpaceException("Cannot interpolate across to the domain of the specified FunctionSpace. (DataConstant)");
      }
      // if the data is not lazy, this will just be a cast to DataReady
      DataReady_ptr dr=inData.m_data->resolve();
      DataConstant* dc=new DataConstant(functionspace,inData.m_data->getShape(),dr->getVectorRO());	
//       m_data=DataAbstract_ptr(dc); 
      set_m_data(DataAbstract_ptr(dc));
    } else {
      Data tmp(0,inData.getDataPointShape(),functionspace,true);
      // Note: Must use a reference or pointer to a derived object
      // in order to get polymorphic behaviour. Shouldn't really
      // be able to create an instance of AbstractDomain but that was done
      // as a boost:python work around which may no longer be required.
      /*const AbstractDomain& inDataDomain=inData.getDomain();*/
      const_Domain_ptr inDataDomain=inData.getDomain();
      if  (inDataDomain==functionspace.getDomain()) {
        inDataDomain->interpolateOnDomain(tmp,inData);
      } else {
        inDataDomain->interpolateACross(tmp,inData);
      }
//       m_data=tmp.m_data;
      set_m_data(tmp.m_data);
    }
  }
  m_protected=false;
}

Data::Data(DataAbstract* underlyingdata)
	: m_shared(false), m_lazy(false)
{
	set_m_data(underlyingdata->getPtr());
	m_protected=false;
}

Data::Data(DataAbstract_ptr underlyingdata)
	: m_shared(false), m_lazy(false)
{
	set_m_data(underlyingdata);
	m_protected=false;
}

Data::Data(const DataTypes::ValueType& value,
		 const DataTypes::ShapeType& shape,
                 const FunctionSpace& what,
                 bool expanded)
	: m_shared(false), m_lazy(false)
{
   initialise(value,shape,what,expanded);
   m_protected=false;
}


Data::Data(const object& value,
	   const FunctionSpace& what,
           bool expanded)
	: m_shared(false), m_lazy(false)
{
  WrappedArray w(value);
  initialise(w,what,expanded);
  m_protected=false;
}


Data::Data(const object& value,
           const Data& other)
	: m_shared(false), m_lazy(false)
{
  WrappedArray w(value);

  // extract the shape of the array
  const DataTypes::ShapeType& tempShape=w.getShape();
  if (w.getRank()==0) {


    // get the space for the data vector
    int len1 = DataTypes::noValues(tempShape);
    DataVector temp_data(len1, 0.0, len1);
    temp_data.copyFromArray(w,1);

    int len = DataTypes::noValues(other.getDataPointShape());

    DataVector temp2_data(len, temp_data[0], len);
    DataConstant* t=new DataConstant(other.getFunctionSpace(),other.getDataPointShape(),temp2_data);
//     m_data=DataAbstract_ptr(t);
    set_m_data(DataAbstract_ptr(t));

  } else {
    //
    // Create a DataConstant with the same sample shape as other
    DataConstant* t=new DataConstant(w,other.getFunctionSpace());
//     m_data=DataAbstract_ptr(t);
    set_m_data(DataAbstract_ptr(t));
  }
  m_protected=false;
}

Data::~Data()
{
  set_m_data(DataAbstract_ptr());
}


// only call in thread safe contexts.
// This method should be atomic
void Data::set_m_data(DataAbstract_ptr p)
{
  if (m_data.get()!=0)	// release old ownership
  {
	m_data->removeOwner(this);
  }
  if (p.get()!=0)
  {
	m_data=p;
	m_data->addOwner(this);
	m_shared=m_data->isShared();
	m_lazy=m_data->isLazy();
  }
}

void Data::initialise(const WrappedArray& value,
                 const FunctionSpace& what,
                 bool expanded)
{
  //
  // Construct a Data object of the appropriate type.
  // Construct the object first as there seems to be a bug which causes
  // undefined behaviour if an exception is thrown during construction
  // within the shared_ptr constructor.
  if (expanded) {
    DataAbstract* temp=new DataExpanded(value, what);
//     m_data=temp->getPtr();
    set_m_data(temp->getPtr());
  } else {
    DataAbstract* temp=new DataConstant(value, what);
//     m_data=temp->getPtr();
    set_m_data(temp->getPtr());
  }
}


void
Data::initialise(const DataTypes::ValueType& value,
		 const DataTypes::ShapeType& shape,
                 const FunctionSpace& what,
                 bool expanded)
{
  //
  // Construct a Data object of the appropriate type.
  // Construct the object first as there seems to be a bug which causes
  // undefined behaviour if an exception is thrown during construction
  // within the shared_ptr constructor.
  if (expanded) {
    DataAbstract* temp=new DataExpanded(what, shape, value);
//     m_data=temp->getPtr();
    set_m_data(temp->getPtr());
  } else {
    DataAbstract* temp=new DataConstant(what, shape, value);
//     m_data=temp->getPtr();
    set_m_data(temp->getPtr());
  }
}


escriptDataC
Data::getDataC()
{
  escriptDataC temp;
  temp.m_dataPtr=(void*)this;
  return temp;
}

escriptDataC
Data::getDataC() const
{
  escriptDataC temp;
  temp.m_dataPtr=(void*)this;
  return temp;
}

size_t
Data::getSampleBufferSize() const
{
  return m_data->getSampleBufferSize();
}


const boost::python::tuple
Data::getShapeTuple() const
{
  const DataTypes::ShapeType& shape=getDataPointShape();
  switch(getDataPointRank()) {
     case 0:
        return make_tuple();
     case 1:
        return make_tuple(long_(shape[0]));
     case 2:
        return make_tuple(long_(shape[0]),long_(shape[1]));
     case 3:
        return make_tuple(long_(shape[0]),long_(shape[1]),long_(shape[2]));
     case 4:
        return make_tuple(long_(shape[0]),long_(shape[1]),long_(shape[2]),long_(shape[3]));
     default:
        throw DataException("Error - illegal Data rank.");
  }
}


// The different name is needed because boost has trouble with overloaded functions.
// It can't work out what type the function is based soley on its name.
// There are ways to fix this involving creating function pointer variables for each form
// but there doesn't seem to be a need given that the methods have the same name from the python point of view
Data
Data::copySelf()
{
   DataAbstract* temp=m_data->deepCopy();
   return Data(temp);
}

void
Data::copy(const Data& other)
{
  DataAbstract* temp=other.m_data->deepCopy();
  DataAbstract_ptr p=temp->getPtr();
//   m_data=p;
  set_m_data(p);
}


Data
Data::delay()
{
  if (!isLazy())
  {
      DataLazy* dl=new DataLazy(m_data);
      return Data(dl);
  }
  return *this;
}

void
Data::delaySelf()
{
  if (!isLazy())
  {
// 	m_data=(new DataLazy(m_data))->getPtr();
	set_m_data((new DataLazy(m_data))->getPtr());
  }
}


// For lazy data, it would seem that DataTagged will need to be treated differently since even after setting all tags
// to zero, all the tags from all the DataTags would be in the result.
// However since they all have the same value (0) whether they are there or not should not matter.
// So I have decided that for all types this method will create a constant 0.
// It can be promoted up as required.
// A possible efficiency concern might be expanded->constant->expanded which has an extra memory management
// but we can deal with that if it arises.
//
void
Data::setToZero()
{
  if (isEmpty())
  {
     throw DataException("Error - Operations not permitted on instances of DataEmpty.");
  }
  if (isLazy())
  {
     DataTypes::ValueType v(getNoValues(),0);
     DataConstant* dc=new DataConstant(getFunctionSpace(),getDataPointShape(),v);
     DataLazy* dl=new DataLazy(dc->getPtr());
     set_m_data(dl->getPtr());
  }
  else
  {
     exclusiveWrite();
     m_data->setToZero();
  }
}


void
Data::copyWithMask(const Data& other,
                   const Data& mask)
{
  // 1. Interpolate if required so all Datas use the same FS as this
  // 2. Tag or Expand so that all Data's are the same type
  // 3. Iterate over the data vectors copying values where mask is >0
  if (other.isEmpty() || mask.isEmpty())
  {
	throw DataException("Error - copyWithMask not permitted using instances of DataEmpty.");
  }
  Data other2(other);
  Data mask2(mask);
  other2.resolve();
  mask2.resolve();
  this->resolve();
  FunctionSpace myFS=getFunctionSpace();
  FunctionSpace oFS=other2.getFunctionSpace();
  FunctionSpace mFS=mask2.getFunctionSpace();
  if (oFS!=myFS)
  {
     if (other2.probeInterpolation(myFS))
     {
	other2=other2.interpolate(myFS);
     }
     else
     {
	throw DataException("Error - copyWithMask: other FunctionSpace is not compatible with this one.");
     }
  }
  if (mFS!=myFS)
  {
     if (mask2.probeInterpolation(myFS))
     {
	mask2=mask2.interpolate(myFS);
     }
     else
     {
	throw DataException("Error - copyWithMask: mask FunctionSpace is not compatible with this one.");
     }
  }
  			// Ensure that all args have the same type
  if (this->isExpanded() || mask2.isExpanded() || other2.isExpanded())
  {
	this->expand();
	other2.expand();
	mask2.expand();
  }
  else if (this->isTagged() || mask2.isTagged() || other2.isTagged())
  {
	this->tag();
	other2.tag();
	mask2.tag();
  }
  else if (this->isConstant() && mask2.isConstant() && other2.isConstant())
  {
  }
  else
  {
	throw DataException("Error - Unknown DataAbstract passed to copyWithMask.");
  }
  unsigned int selfrank=getDataPointRank();
  unsigned int otherrank=other2.getDataPointRank();
  unsigned int maskrank=mask2.getDataPointRank();
  if ((selfrank==0) && (otherrank>0 || maskrank>0))
  {
	// to get here we must be copying from a large object into a scalar
	// I am not allowing this.
	// If you are calling copyWithMask then you are considering keeping some existing values
	// and so I'm going to assume that you don't want your data objects getting a new shape.
	throw DataException("Attempt to copyWithMask into a scalar from an object or mask with rank>0.");
  }
  exclusiveWrite();
  // Now we iterate over the elements
  DataVector& self=getReady()->getVectorRW();;
  const DataVector& ovec=other2.getReadyPtr()->getVectorRO();
  const DataVector& mvec=mask2.getReadyPtr()->getVectorRO();

  if ((selfrank>0) && (otherrank==0) &&(maskrank==0))
  {
	// Not allowing this combination.
	// it is not clear what the rank of the target should be.
	// Should it be filled with the scalar (rank stays the same); 
	// or should the target object be reshaped to be a scalar as well.
	throw DataException("Attempt to copyWithMask from scalar mask and data into non-scalar target.");
  }
  if ((selfrank>0) && (otherrank>0) &&(maskrank==0))
  {
	if (mvec[0]>0)		// copy whole object if scalar is >0
	{
	    copy(other);
	}
	return;
  }
  if (isTagged())		// so all objects involved will also be tagged
  {
	// note the !
	if (!((getDataPointShape()==mask2.getDataPointShape()) && 
		((other2.getDataPointShape()==mask2.getDataPointShape()) || (otherrank==0))))
	{
		throw DataException("copyWithMask, shape mismatch.");
	}

	// We need to consider the possibility that tags are missing or in the wrong order
	// My guiding assumption here is: All tagged Datas are assumed to have the default value for
	// all tags which are not explicitly defined

	const DataTagged* mptr=dynamic_cast<const DataTagged*>(mask2.m_data.get());
	const DataTagged* optr=dynamic_cast<const DataTagged*>(other2.m_data.get());
	DataTagged* tptr=dynamic_cast<DataTagged*>(m_data.get());

	// first, add any tags missing from other or mask
	const DataTagged::DataMapType& olookup=optr->getTagLookup();
        const DataTagged::DataMapType& mlookup=mptr->getTagLookup();
	const DataTagged::DataMapType& tlookup=tptr->getTagLookup();
 	DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
	for (i=olookup.begin();i!=olookup.end();i++)
	{
           tptr->addTag(i->first); 
        }
        for (i=mlookup.begin();i!=mlookup.end();i++) {
           tptr->addTag(i->first);
        }
	// now we know that *this has all the required tags but they aren't guaranteed to be in
	// the same order

	// There are two possibilities: 1. all objects have the same rank. 2. other is a scalar
	if ((selfrank==otherrank) && (otherrank==maskrank))
	{
		for (i=tlookup.begin();i!=tlookup.end();i++)
		{
			// get the target offset
			DataTypes::ValueType::size_type toff=tptr->getOffsetForTag(i->first);
           		DataTypes::ValueType::size_type moff=mptr->getOffsetForTag(i->first);
			DataTypes::ValueType::size_type ooff=optr->getOffsetForTag(i->first);
			for (int j=0;j<getDataPointSize();++j)
			{
				if (mvec[j+moff]>0)
				{
					self[j+toff]=ovec[j+ooff];
				}
			}
        	}
		// now for the default value
		for (int j=0;j<getDataPointSize();++j)
		{
			if (mvec[j+mptr->getDefaultOffset()]>0)
			{
				self[j+tptr->getDefaultOffset()]=ovec[j+optr->getDefaultOffset()];
			}
		}
	}
	else	// other is a scalar
	{
		for (i=tlookup.begin();i!=tlookup.end();i++)
		{
			// get the target offset
			DataTypes::ValueType::size_type toff=tptr->getOffsetForTag(i->first);
           		DataTypes::ValueType::size_type moff=mptr->getOffsetForTag(i->first);
			DataTypes::ValueType::size_type ooff=optr->getOffsetForTag(i->first);
			for (int j=0;j<getDataPointSize();++j)
			{
				if (mvec[j+moff]>0)
				{
					self[j+toff]=ovec[ooff];
				}
			}
        	}
		// now for the default value
		for (int j=0;j<getDataPointSize();++j)
		{
			if (mvec[j+mptr->getDefaultOffset()]>0)
			{
				self[j+tptr->getDefaultOffset()]=ovec[0];
			}
		}
	}

	return;			// ugly
  }
  // mixed scalar and non-scalar operation
  if ((selfrank>0) && (otherrank==0) && (mask2.getDataPointShape()==getDataPointShape()))
  {
    	size_t num_points=self.size();
  	// OPENMP 3.0 allows unsigned loop vars.
	#if defined(_OPENMP) && (_OPENMP < 200805)
  	long i;
	#else
  	size_t i;
	#endif	
	size_t psize=getDataPointSize();	
	#pragma omp parallel for private(i) schedule(static)
  	for (i=0;i<num_points;++i)
  	{
		if (mvec[i]>0)
		{
	   	    self[i]=ovec[i/psize];		// since this is expanded there is one scalar 
		}					// dest point
  	} 
	return;			// ugly!
  }
  // tagged data is already taken care of so we only need to worry about shapes
  // special cases with scalars are already dealt with so all we need to worry about is shape
  if ((getDataPointShape()!=other2.getDataPointShape()) || getDataPointShape()!=mask2.getDataPointShape())
  {
	ostringstream oss;
	oss <<"Error - size mismatch in arguments to copyWithMask.";
	oss << "\nself_shape=" << DataTypes::shapeToString(getDataPointShape());
	oss << " other2_shape=" << DataTypes::shapeToString(other2.getDataPointShape());
	oss << " mask2_shape=" << DataTypes::shapeToString(mask2.getDataPointShape());
	throw DataException(oss.str());
  }
  size_t num_points=self.size();

  // OPENMP 3.0 allows unsigned loop vars.
#if defined(_OPENMP) && (_OPENMP < 200805)
  long i;
#else
  size_t i;
#endif
  #pragma omp parallel for private(i) schedule(static)
  for (i=0;i<num_points;++i)
  {
	if (mvec[i]>0)
	{
	   self[i]=ovec[i];
	}
  }
}

bool
Data::isExpanded() const
{
  DataExpanded* temp=dynamic_cast<DataExpanded*>(m_data.get());
  return (temp!=0);
}

bool
Data::actsExpanded() const
{
  return m_data->actsExpanded();

}


bool
Data::isTagged() const
{
  DataTagged* temp=dynamic_cast<DataTagged*>(m_data.get());
  return (temp!=0);
}

bool
Data::isEmpty() const
{
  DataEmpty* temp=dynamic_cast<DataEmpty*>(m_data.get());
  return (temp!=0);
}

bool
Data::isConstant() const
{
  DataConstant* temp=dynamic_cast<DataConstant*>(m_data.get());
  return (temp!=0);
}

bool
Data::isLazy() const
{
  return m_lazy;	// not asking m_data because we need to be able to ask this while m_data is changing
}

// at the moment this is synonymous with !isLazy() but that could change
bool
Data::isReady() const
{
  return (dynamic_cast<DataReady*>(m_data.get())!=0);
}


void
Data::setProtection()
{
   m_protected=true;
}

bool
Data::isProtected() const
{
   return m_protected;
}



void
Data::expand()
{
  if (isConstant()) {
    DataConstant* tempDataConst=dynamic_cast<DataConstant*>(m_data.get());
    DataAbstract* temp=new DataExpanded(*tempDataConst);
//     m_data=temp->getPtr();
    set_m_data(temp->getPtr());
  } else if (isTagged()) {
    DataTagged* tempDataTag=dynamic_cast<DataTagged*>(m_data.get());
    DataAbstract* temp=new DataExpanded(*tempDataTag);
//     m_data=temp->getPtr();
    set_m_data(temp->getPtr());
  } else if (isExpanded()) {
    //
    // do nothing
  } else if (isEmpty()) {
    throw DataException("Error - Expansion of DataEmpty not possible.");
  } else if (isLazy()) {
    resolve();
    expand();		// resolve might not give us expanded data
  } else {
    throw DataException("Error - Expansion not implemented for this Data type.");
  }
}

void
Data::tag()
{
  if (isConstant()) {
    DataConstant* tempDataConst=dynamic_cast<DataConstant*>(m_data.get());
    DataAbstract* temp=new DataTagged(*tempDataConst);
//     m_data=temp->getPtr();
    set_m_data(temp->getPtr());
  } else if (isTagged()) {
    // do nothing
  } else if (isExpanded()) {
    throw DataException("Error - Creating tag data from DataExpanded not possible.");
  } else if (isEmpty()) {
    throw DataException("Error - Creating tag data from DataEmpty not possible.");
  } else if (isLazy()) {
     DataAbstract_ptr res=m_data->resolve();
     if (m_data->isExpanded())
     {
	throw DataException("Error - data would resolve to DataExpanded, tagging is not possible.");
     }
//      m_data=res;	
     set_m_data(res);
     tag();
  } else {
    throw DataException("Error - Tagging not implemented for this Data type.");
  }
}

void
Data::resolve()
{
  if (isLazy())
  {
//      m_data=m_data->resolve();
    set_m_data(m_data->resolve());
  }
}

void 
Data::requireWrite()
{
  resolve();
  exclusiveWrite();
}

Data
Data::oneOver() const
{
  MAKELAZYOP(RECIP)
  return C_TensorUnaryOperation(*this, bind1st(divides<double>(),1.));
}

Data
Data::wherePositive() const
{
  MAKELAZYOP(GZ)
  return C_TensorUnaryOperation(*this, bind2nd(greater<double>(),0.0));
}

Data
Data::whereNegative() const
{
  MAKELAZYOP(LZ)
  return C_TensorUnaryOperation(*this, bind2nd(less<double>(),0.0));
}

Data
Data::whereNonNegative() const
{
  MAKELAZYOP(GEZ)
  return C_TensorUnaryOperation(*this, bind2nd(greater_equal<double>(),0.0));
}

Data
Data::whereNonPositive() const
{
  MAKELAZYOP(LEZ)
  return C_TensorUnaryOperation(*this, bind2nd(less_equal<double>(),0.0));
}

Data
Data::whereZero(double tol) const
{
//   Data dataAbs=abs();
//   return C_TensorUnaryOperation(dataAbs, bind2nd(less_equal<double>(),tol));
   MAKELAZYOPOFF(EZ,tol)
   return C_TensorUnaryOperation(*this, bind2nd(AbsLTE(),tol));

}

Data
Data::whereNonZero(double tol) const
{
//   Data dataAbs=abs();
//   return C_TensorUnaryOperation(dataAbs, bind2nd(greater<double>(),tol));
  MAKELAZYOPOFF(NEZ,tol)
  return C_TensorUnaryOperation(*this, bind2nd(AbsGT(),tol));

}

Data
Data::interpolate(const FunctionSpace& functionspace) const
{
  return Data(*this,functionspace);
}

bool
Data::probeInterpolation(const FunctionSpace& functionspace) const
{
  return getFunctionSpace().probeInterpolation(functionspace);
}

Data
Data::gradOn(const FunctionSpace& functionspace) const
{
  if (isEmpty())
  {
	throw DataException("Error - operation not permitted on instances of DataEmpty.");
  }
  double blocktimer_start = blocktimer_time();
  if (functionspace.getDomain()!=getDomain())
    throw DataException("Error - gradient cannot be calculated on different domains.");
  DataTypes::ShapeType grad_shape=getDataPointShape();
  grad_shape.push_back(functionspace.getDim());
  Data out(0.0,grad_shape,functionspace,true);
  getDomain()->setToGradient(out,*this);
  blocktimer_increment("grad()", blocktimer_start);
  return out;
}

Data
Data::grad() const
{
  if (isEmpty())
  {
	throw DataException("Error - operation not permitted on instances of DataEmpty.");
  }
  return gradOn(escript::function(*getDomain()));
}

int
Data::getDataPointSize() const
{
  return m_data->getNoValues();
}


DataTypes::ValueType::size_type
Data::getLength() const
{
  return m_data->getLength();
}


// There is no parallelism here ... elements need to be added in the correct order.
//   If we could presize the list and then fill in the elements it might work
//   This would need setting elements to be threadsafe.
//   Having mulitple C threads calling into one interpreter is aparently a no-no.
const boost::python::object
Data::toListOfTuples(bool scalarastuple)
{
    using namespace boost::python;
    using boost::python::list;
    if (get_MPISize()>1)
    {
        throw DataException("::toListOfTuples is not available for MPI with more than one process.");
    }
    unsigned int rank=getDataPointRank();
    unsigned int size=getDataPointSize();

    int npoints=getNumDataPoints();
    expand();			// This will also resolve if required
    const DataTypes::ValueType& vec=getReady()->getVectorRO();
    boost::python::list temp;
    temp.append(object());
    boost::python::list res(temp*npoints);// presize the list by the "[None] * npoints"  trick
    if (rank==0)
    {
        long count;
        if (scalarastuple)
        {
            for (count=0;count<npoints;++count)
            {
		res[count]=make_tuple(vec[count]);
            }
        }
        else
        {
            for (count=0;count<npoints;++count)
            {
                res[count]=vec[count];
            }
        }
    }
    else if (rank==1)
    {
        size_t count;
        size_t offset=0;
        for (count=0;count<npoints;++count,offset+=size)
        {
            res[count]=pointToTuple1(getDataPointShape(), vec, offset);
        }
    }
    else if (rank==2)
    {
        size_t count;
        size_t offset=0;
        for (count=0;count<npoints;++count,offset+=size)
        {
	    res[count]=pointToTuple2(getDataPointShape(), vec, offset);
        }
    }
    else if (rank==3)
    {
        size_t count;
        size_t offset=0;
        for (count=0;count<npoints;++count,offset+=size)
        {
            res[count]=pointToTuple3(getDataPointShape(), vec, offset);
        }
    }
    else if (rank==4)
    {
        size_t count;
        size_t offset=0;
        for (count=0;count<npoints;++count,offset+=size)
        {
            res[count]=pointToTuple4(getDataPointShape(), vec, offset);
        }
    }
    else
    {
        throw DataException("Unknown rank in ::toListOfTuples()");
    }
    return res;
}

const boost::python::object
Data::getValueOfDataPointAsTuple(int dataPointNo)
{
   forceResolve();
   if (getNumDataPointsPerSample()>0) {
       int sampleNo = dataPointNo/getNumDataPointsPerSample();
       int dataPointNoInSample = dataPointNo - sampleNo * getNumDataPointsPerSample();
       //
       // Check a valid sample number has been supplied
       if ((sampleNo >= getNumSamples()) || (sampleNo < 0 )) {
           throw DataException("Error - Data::getValueOfDataPointAsTuple: invalid sampleNo.");
       }

       //
       // Check a valid data point number has been supplied
       if ((dataPointNoInSample >= getNumDataPointsPerSample()) || (dataPointNoInSample < 0)) {
           throw DataException("Error - Data::getValueOfDataPointAsTuple: invalid dataPointNoInSample.");
       }
       // TODO: global error handling
       DataTypes::ValueType::size_type offset=getDataOffset(sampleNo, dataPointNoInSample);
       return pointToTuple(getDataPointShape(),&(getDataAtOffsetRO(offset)));
  }
  else
  {
	// The pre-numpy method would return an empty array of the given shape
	// I'm going to throw an exception because if we have zero points per sample we have problems
	throw DataException("Error - need at least 1 datapoint per sample.");
  }

}


void
Data::setValueOfDataPointToPyObject(int dataPointNo, const boost::python::object& py_object)
{
    // this will throw if the value cannot be represented
    setValueOfDataPointToArray(dataPointNo,py_object);
}

void
Data::setValueOfDataPointToArray(int dataPointNo, const boost::python::object& obj)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }

  WrappedArray w(obj);
  //
  // check rank
  if (static_cast<unsigned int>(w.getRank())<getDataPointRank())
      throw DataException("Rank of array does not match Data object rank");

  //
  // check shape of array
  for (unsigned int i=0; i<getDataPointRank(); i++) {
    if (w.getShape()[i]!=getDataPointShape()[i])
       throw DataException("Shape of array does not match Data object rank");
  }

  exclusiveWrite();

  //
  // make sure data is expanded:
  //
  if (!isExpanded()) {
    expand();
  }
  if (getNumDataPointsPerSample()>0) {
       int sampleNo = dataPointNo/getNumDataPointsPerSample();
       int dataPointNoInSample = dataPointNo - sampleNo * getNumDataPointsPerSample();
       m_data->copyToDataPoint(sampleNo, dataPointNoInSample,w);
  } else {
       m_data->copyToDataPoint(-1, 0,w);
  }
}

void
Data::setValueOfDataPoint(int dataPointNo, const double value)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  //
  // make sure data is expanded:
  exclusiveWrite();
  if (!isExpanded()) {
    expand();
  }
  if (getNumDataPointsPerSample()>0) {
       int sampleNo = dataPointNo/getNumDataPointsPerSample();
       int dataPointNoInSample = dataPointNo - sampleNo * getNumDataPointsPerSample();
       m_data->copyToDataPoint(sampleNo, dataPointNoInSample,value);
  } else {
       m_data->copyToDataPoint(-1, 0,value);
  }
}

const
boost::python::object
Data::getValueOfGlobalDataPointAsTuple(int procNo, int dataPointNo)
{
  // This could be lazier than it is now
  forceResolve();

  // copy datapoint into a buffer
  // broadcast buffer to all nodes
  // convert buffer to tuple
  // return tuple

  const DataTypes::ShapeType& dataPointShape = getDataPointShape();
  size_t length=DataTypes::noValues(dataPointShape);

  // added for the MPI communication
  double *tmpData = new double[length];

  // updated for the MPI case
  if( get_MPIRank()==procNo ){
      if (getNumDataPointsPerSample()>0) {
	int sampleNo = dataPointNo/getNumDataPointsPerSample();
	int dataPointNoInSample = dataPointNo - sampleNo * getNumDataPointsPerSample();
	//
	// Check a valid sample number has been supplied
	if ((sampleNo >= getNumSamples()) || (sampleNo < 0 )) {
		throw DataException("Error - Data::getValueOfGlobalDataPointAsTuple: invalid sampleNo.");
	}

	//
	// Check a valid data point number has been supplied
	if ((dataPointNoInSample >= getNumDataPointsPerSample()) || (dataPointNoInSample < 0)) {
		throw DataException("Error - Data::getValueOfGlobalDataPointAsTuple: invalid dataPointNoInSample.");
	}
	// TODO: global error handling
	DataTypes::ValueType::size_type offset=getDataOffset(sampleNo, dataPointNoInSample);

	memcpy(tmpData,&(getDataAtOffsetRO(offset)),length*sizeof(double));
     }
  }
#ifdef PASO_MPI
  // broadcast the data to all other processes
  MPI_Bcast( tmpData, length, MPI_DOUBLE, procNo, get_MPIComm() );
#endif

  boost::python::tuple t=pointToTuple(dataPointShape,tmpData);
  delete [] tmpData;
  //
  // return the loaded array
  return t;

}


boost::python::object
Data::integrateToTuple_const() const
{
  if (isLazy())
  {
	throw DataException("Error - cannot integrate for constant lazy data.");
  }
  return integrateWorker();
}

boost::python::object
Data::integrateToTuple()
{
  if (isLazy())
  {
	expand(); 	// Can't do a non-resolving version of this without changing the domain object
  }			// see the dom->setToIntegrals call. Not saying it can't be done, just not doing it yet.
  return integrateWorker();

}

boost::python::object
Data::integrateWorker() const
{
  DataTypes::ShapeType shape = getDataPointShape();
  int dataPointSize = getDataPointSize();

  //
  // calculate the integral values
  vector<double> integrals(dataPointSize);
  vector<double> integrals_local(dataPointSize);
  const AbstractContinuousDomain* dom=dynamic_cast<const AbstractContinuousDomain*>(getDomain().get());
  if (dom==0)
  {				
    throw DataException("Can not integrate over non-continuous domains.");
  }
#ifdef PASO_MPI
  dom->setToIntegrals(integrals_local,*this);
  // Global sum: use an array instead of a vector because elements of array are guaranteed to be contiguous in memory
  double *tmp = new double[dataPointSize];
  double *tmp_local = new double[dataPointSize];
  for (int i=0; i<dataPointSize; i++) { tmp_local[i] = integrals_local[i]; }
  MPI_Allreduce( &tmp_local[0], &tmp[0], dataPointSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  for (int i=0; i<dataPointSize; i++) { integrals[i] = tmp[i]; }
  tuple result=pointToTuple(shape,tmp);
  delete[] tmp;
  delete[] tmp_local;
#else
  dom->setToIntegrals(integrals,*this);
/*  double *tmp = new double[dataPointSize];
  for (int i=0; i<dataPointSize; i++) { tmp[i]=integrals[i]; }*/
  tuple result=pointToTuple(shape,integrals);
//   delete tmp;
#endif


  return result;
}

Data
Data::sin() const
{
  MAKELAZYOP(SIN)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::sin);
}

Data
Data::cos() const
{
  MAKELAZYOP(COS)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::cos);
}

Data
Data::tan() const
{
  MAKELAZYOP(TAN)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::tan);
}

Data
Data::asin() const
{
  MAKELAZYOP(ASIN)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::asin);
}

Data
Data::acos() const
{
  MAKELAZYOP(ACOS)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::acos);
}


Data
Data::atan() const
{
  MAKELAZYOP(ATAN)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::atan);
}

Data
Data::sinh() const
{
  MAKELAZYOP(SINH)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::sinh);
}

Data
Data::cosh() const
{
  MAKELAZYOP(COSH)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::cosh);
}

Data
Data::tanh() const
{
  MAKELAZYOP(TANH)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::tanh);
}


Data
Data::erf() const
{
#if defined (_WIN32) && !defined(__INTEL_COMPILER)
  throw DataException("Error - Data:: erf function is not supported on _WIN32 platforms.");
#else
  MAKELAZYOP(ERF)
  return C_TensorUnaryOperation(*this, ::erf);
#endif
}

Data
Data::asinh() const
{
  MAKELAZYOP(ASINH)
#if defined (_WIN32) && !defined(__INTEL_COMPILER)
  return C_TensorUnaryOperation(*this, escript::asinh_substitute);
#else
  return C_TensorUnaryOperation(*this, ::asinh);
#endif
}

Data
Data::acosh() const
{
  MAKELAZYOP(ACOSH)
#if defined (_WIN32) && !defined(__INTEL_COMPILER)
  return C_TensorUnaryOperation(*this, escript::acosh_substitute);
#else
  return C_TensorUnaryOperation(*this, ::acosh);
#endif
}

Data
Data::atanh() const
{
  MAKELAZYOP(ATANH)
#if defined (_WIN32) && !defined(__INTEL_COMPILER)
  return C_TensorUnaryOperation(*this, escript::atanh_substitute);
#else
  return C_TensorUnaryOperation(*this, ::atanh);
#endif
}

Data
Data::log10() const
{
  MAKELAZYOP(LOG10)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::log10);
}

Data
Data::log() const
{
  MAKELAZYOP(LOG)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::log);
}

Data
Data::sign() const
{
  MAKELAZYOP(SIGN)
  return C_TensorUnaryOperation(*this, escript::fsign);
}

Data
Data::abs() const
{
  MAKELAZYOP(ABS)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::fabs);
}

Data
Data::neg() const
{
  MAKELAZYOP(NEG)
  return C_TensorUnaryOperation(*this, negate<double>());
}

Data
Data::pos() const
{
	// not doing lazy check here is deliberate.
	// since a deep copy of lazy data should be cheap, I'll just let it happen now
  Data result;
  // perform a deep copy
  result.copy(*this);
  return result;
}

Data
Data::exp() const
{
  MAKELAZYOP(EXP)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::exp);
}

Data
Data::sqrt() const
{
  MAKELAZYOP(SQRT)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::sqrt);
}

double
Data::Lsup_const() const
{
   if (isLazy())
   {
	throw DataException("Error - cannot compute Lsup for constant lazy data.");
   }
   return LsupWorker();
}

double
Data::Lsup() 
{
   if (isLazy())
   {
	if (!actsExpanded() || CHECK_DO_CRES)
	{
	    resolve();
	}
	else
	{
#ifdef PASO_MPI
	    return lazyAlgWorker<AbsMax>(0,MPI_MAX);
#else
	    return lazyAlgWorker<AbsMax>(0);
#endif
	}
   }
   return LsupWorker();
}

double
Data::sup_const() const
{
   if (isLazy())
   {
	throw DataException("Error - cannot compute sup for constant lazy data.");
   }
   return supWorker();
}

double
Data::sup() 
{
   if (isLazy())
   {
	if (!actsExpanded() || CHECK_DO_CRES)
	{
	    resolve();
	}
	else
	{
#ifdef PASO_MPI
	    return lazyAlgWorker<FMax>(numeric_limits<double>::max()*-1, MPI_MAX);
#else
	    return lazyAlgWorker<FMax>(numeric_limits<double>::max()*-1);
#endif
	}
   }
   return supWorker();
}

double
Data::inf_const() const
{
   if (isLazy())
   {
	throw DataException("Error - cannot compute inf for constant lazy data.");
   }
   return infWorker();
}

double
Data::inf() 
{
   if (isLazy())
   {
	if (!actsExpanded() || CHECK_DO_CRES)
	{
	    resolve();
	}
	else
	{
#ifdef PASO_MPI
	    return lazyAlgWorker<FMin>(numeric_limits<double>::max(), MPI_MIN);
#else
	    return lazyAlgWorker<FMin>(numeric_limits<double>::max());
#endif
	}
   }
   return infWorker();
}

template <class BinaryOp>
double
#ifdef PASO_MPI
Data::lazyAlgWorker(double init, MPI_Op mpiop_type)
#else
Data::lazyAlgWorker(double init)
#endif
{
   if (!isLazy() || !m_data->actsExpanded())
   {
	throw DataException("Error - lazyAlgWorker can only be called on lazy(expanded) data.");
   }
   DataLazy* dl=dynamic_cast<DataLazy*>(m_data.get());
   EsysAssert((dl!=0), "Programming error - casting to DataLazy.");
   BufferGroup* bg=allocSampleBuffer();
   double val=init;
   int i=0;
   const size_t numsamples=getNumSamples();
   const size_t samplesize=getNoValues()*getNumDataPointsPerSample();
   BinaryOp operation;
   #pragma omp parallel private(i)
   {
	double localtot=init;
	#pragma omp for schedule(static) private(i)
	for (i=0;i<numsamples;++i)
	{
	    size_t roffset=0;
	    const DataTypes::ValueType* v=dl->resolveSample(*bg, i, roffset);
	    // Now we have the sample, run operation on all points
	    for (size_t j=0;j<samplesize;++j)
	    {
		localtot=operation(localtot,(*v)[j+roffset]);
	    }
	}
	#pragma omp critical
	val=operation(val,localtot);
   }
   freeSampleBuffer(bg);
#ifdef PASO_MPI
   double globalValue;
   MPI_Allreduce( &val, &globalValue, 1, MPI_DOUBLE, mpiop_type, MPI_COMM_WORLD );
   return globalValue;
#else
   return val;
#endif
}

double
Data::LsupWorker() const
{
  double localValue;
  //
  // set the initial absolute maximum value to zero

  AbsMax abs_max_func;
  localValue = algorithm(abs_max_func,0);
#ifdef PASO_MPI
  double globalValue;
  MPI_Allreduce( &localValue, &globalValue, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
  return globalValue;
#else
  return localValue;
#endif
}

double
Data::supWorker() const
{
  double localValue;
  //
  // set the initial maximum value to min possible double
  FMax fmax_func;
  localValue = algorithm(fmax_func,numeric_limits<double>::max()*-1);
#ifdef PASO_MPI
  double globalValue;
  MPI_Allreduce( &localValue, &globalValue, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
  return globalValue;
#else
  return localValue;
#endif
}

double
Data::infWorker() const
{
  double localValue;
  //
  // set the initial minimum value to max possible double
  FMin fmin_func;
  localValue = algorithm(fmin_func,numeric_limits<double>::max());
#ifdef PASO_MPI
  double globalValue;
  MPI_Allreduce( &localValue, &globalValue, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
  return globalValue;
#else
  return localValue;
#endif
}

/* global reduction */


inline Data
Data::minval_nonlazy() const
{
  //
  // set the initial minimum value to max possible double
  FMin fmin_func;
  return dp_algorithm(fmin_func,numeric_limits<double>::max());
}


inline Data
Data::maxval_nonlazy() const
{
  //
  // set the initial maximum value to min possible double
  FMax fmax_func;
  return dp_algorithm(fmax_func,numeric_limits<double>::max()*-1);
}



Data
Data::maxval() const
{
   MAKELAZYOP(MAXVAL)
   return maxval_nonlazy();
}


Data
Data::minval() const
{
  MAKELAZYOP(MINVAL)
  return minval_nonlazy();
}


Data
Data::swapaxes(const int axis0, const int axis1) const
{
     int axis0_tmp,axis1_tmp;
     DataTypes::ShapeType s=getDataPointShape();
     DataTypes::ShapeType ev_shape;
     // Here's the equivalent of python s_out=s[axis_offset:]+s[:axis_offset]
     // which goes thru all shape vector elements starting with axis_offset (at index=rank wrap around to 0)
     int rank=getDataPointRank();
     if (rank<2) {
        throw DataException("Error - Data::swapaxes argument must have at least rank 2.");
     }
     if (axis0<0 || axis0>rank-1) {
        throw DataException("Error - Data::swapaxes: axis0 must be between 0 and rank-1=" + rank-1);
     }
     if (axis1<0 || axis1>rank-1) {
         throw DataException("Error - Data::swapaxes: axis1 must be between 0 and rank-1=" + rank-1);
     }
     if (axis0 == axis1) {
         throw DataException("Error - Data::swapaxes: axis indices must be different.");
     }
     MAKELAZYOP2(SWAP,axis0,axis1)
     if (axis0 > axis1)
     {
	axis0_tmp=axis1;
	axis1_tmp=axis0;
     }
     else
     {
	axis0_tmp=axis0;
	axis1_tmp=axis1;
     }
     for (int i=0; i<rank; i++)
     {
	if (i == axis0_tmp)
	{
		ev_shape.push_back(s[axis1_tmp]);
	} 
	else if (i == axis1_tmp)
	{
		ev_shape.push_back(s[axis0_tmp]);
	}
	else
	{
		ev_shape.push_back(s[i]);
	}
     }
     Data ev(0.,ev_shape,getFunctionSpace());
     ev.typeMatchRight(*this);
     m_data->swapaxes(ev.m_data.get(), axis0_tmp, axis1_tmp);
     return ev;
}

Data
Data::symmetric() const
{
     // check input
     DataTypes::ShapeType s=getDataPointShape();
     if (getDataPointRank()==2) {
        if(s[0] != s[1])
           throw DataException("Error - Data::symmetric can only be calculated for rank 2 object with equal first and second dimension.");
     }
     else if (getDataPointRank()==4) {
        if(!(s[0] == s[2] && s[1] == s[3]))
           throw DataException("Error - Data::symmetric can only be calculated for rank 4 object with dim0==dim2 and dim1==dim3.");
     }
     else {
        throw DataException("Error - Data::symmetric can only be calculated for rank 2 or 4 object.");
     }
     MAKELAZYOP(SYM)
     Data ev(0.,getDataPointShape(),getFunctionSpace());
     ev.typeMatchRight(*this);
     m_data->symmetric(ev.m_data.get());
     return ev;
}

Data
Data::nonsymmetric() const
{
     MAKELAZYOP(NSYM)
     // check input
     DataTypes::ShapeType s=getDataPointShape();
     if (getDataPointRank()==2) {
        if(s[0] != s[1])
           throw DataException("Error - Data::nonsymmetric can only be calculated for rank 2 object with equal first and second dimension.");
        DataTypes::ShapeType ev_shape;
        ev_shape.push_back(s[0]);
        ev_shape.push_back(s[1]);
        Data ev(0.,ev_shape,getFunctionSpace());
        ev.typeMatchRight(*this);
        m_data->nonsymmetric(ev.m_data.get());
        return ev;
     }
     else if (getDataPointRank()==4) {
        if(!(s[0] == s[2] && s[1] == s[3]))
           throw DataException("Error - Data::nonsymmetric can only be calculated for rank 4 object with dim0==dim2 and dim1==dim3.");
        DataTypes::ShapeType ev_shape;
        ev_shape.push_back(s[0]);
        ev_shape.push_back(s[1]);
        ev_shape.push_back(s[2]);
        ev_shape.push_back(s[3]);
        Data ev(0.,ev_shape,getFunctionSpace());
        ev.typeMatchRight(*this);
        m_data->nonsymmetric(ev.m_data.get());
        return ev;
     }
     else {
        throw DataException("Error - Data::nonsymmetric can only be calculated for rank 2 or 4 object.");
     }
}

Data
Data::trace(int axis_offset) const
{     
     MAKELAZYOPOFF(TRACE,axis_offset)
     if ((axis_offset<0) || (axis_offset>getDataPointRank()))
     {
	throw DataException("Error - Data::trace, axis_offset must be between 0 and rank-2 inclusive.");
     }
     DataTypes::ShapeType s=getDataPointShape();
     if (getDataPointRank()==2) {
        DataTypes::ShapeType ev_shape;
        Data ev(0.,ev_shape,getFunctionSpace());
        ev.typeMatchRight(*this);
        m_data->trace(ev.m_data.get(), axis_offset);
        return ev;
     }
     if (getDataPointRank()==3) {
        DataTypes::ShapeType ev_shape;
        if (axis_offset==0) {
          int s2=s[2];
          ev_shape.push_back(s2);
        }
        else if (axis_offset==1) {
          int s0=s[0];
          ev_shape.push_back(s0);
        }
        Data ev(0.,ev_shape,getFunctionSpace());
        ev.typeMatchRight(*this);
        m_data->trace(ev.m_data.get(), axis_offset);
        return ev;
     }
     if (getDataPointRank()==4) {
        DataTypes::ShapeType ev_shape;
        if (axis_offset==0) {
          ev_shape.push_back(s[2]);
          ev_shape.push_back(s[3]);
        }
        else if (axis_offset==1) {
          ev_shape.push_back(s[0]);
          ev_shape.push_back(s[3]);
        }
	else if (axis_offset==2) {
	  ev_shape.push_back(s[0]);
	  ev_shape.push_back(s[1]);
	}
        Data ev(0.,ev_shape,getFunctionSpace());
        ev.typeMatchRight(*this);
	m_data->trace(ev.m_data.get(), axis_offset);
        return ev;
     }
     else {
        throw DataException("Error - Data::trace can only be calculated for rank 2, 3 or 4 object.");
     }
}

Data
Data::transpose(int axis_offset) const
{     
     MAKELAZYOPOFF(TRANS,axis_offset)
     DataTypes::ShapeType s=getDataPointShape();
     DataTypes::ShapeType ev_shape;
     // Here's the equivalent of python s_out=s[axis_offset:]+s[:axis_offset]
     // which goes thru all shape vector elements starting with axis_offset (at index=rank wrap around to 0)
     int rank=getDataPointRank();
     if (axis_offset<0 || axis_offset>rank) {
        throw DataException("Error - Data::transpose must have 0 <= axis_offset <= rank=" + rank);
     }
     for (int i=0; i<rank; i++) {

       int index = (axis_offset+i)%rank;
       ev_shape.push_back(s[index]); // Append to new shape
     }
     Data ev(0.,ev_shape,getFunctionSpace());
     ev.typeMatchRight(*this);
     m_data->transpose(ev.m_data.get(), axis_offset);
     return ev;
}

Data
Data::eigenvalues() const
{
     if (isLazy())
     {
	Data temp(*this);	// to get around the fact that you can't resolve a const Data
	temp.resolve();
	return temp.eigenvalues();
     }
     // check input
     DataTypes::ShapeType s=getDataPointShape();
     if (getDataPointRank()!=2)
        throw DataException("Error - Data::eigenvalues can only be calculated for rank 2 object.");
     if(s[0] != s[1])
        throw DataException("Error - Data::eigenvalues can only be calculated for object with equal first and second dimension.");
     // create return
     DataTypes::ShapeType ev_shape(1,s[0]);
     Data ev(0.,ev_shape,getFunctionSpace());
     ev.typeMatchRight(*this);
     m_data->eigenvalues(ev.m_data.get());
     return ev;
}

const boost::python::tuple
Data::eigenvalues_and_eigenvectors(const double tol) const
{
     if (isLazy())
     {
	Data temp(*this);	// to get around the fact that you can't resolve a const Data
	temp.resolve();
	return temp.eigenvalues_and_eigenvectors(tol);
     }
     DataTypes::ShapeType s=getDataPointShape();
     if (getDataPointRank()!=2)
        throw DataException("Error - Data::eigenvalues and eigenvectors can only be calculated for rank 2 object.");
     if(s[0] != s[1])
        throw DataException("Error - Data::eigenvalues and eigenvectors can only be calculated for object with equal first and second dimension.");
     // create return
     DataTypes::ShapeType ev_shape(1,s[0]);
     Data ev(0.,ev_shape,getFunctionSpace());
     ev.typeMatchRight(*this);
     DataTypes::ShapeType V_shape(2,s[0]);
     Data V(0.,V_shape,getFunctionSpace());
     V.typeMatchRight(*this);
     m_data->eigenvalues_and_eigenvectors(ev.m_data.get(),V.m_data.get(),tol);
     return make_tuple(boost::python::object(ev),boost::python::object(V));
}

const boost::python::tuple
Data::minGlobalDataPoint() const
{
  // NB: calc_minGlobalDataPoint( had to be split off from minGlobalDataPoint( as boost::make_tuple causes an
  // abort (for unknown reasons) if there are openmp directives with it in the
  // surrounding function

  int DataPointNo;
  int ProcNo;
  calc_minGlobalDataPoint(ProcNo,DataPointNo);
  return make_tuple(ProcNo,DataPointNo);
}

void
Data::calc_minGlobalDataPoint(int& ProcNo,
       	                int& DataPointNo) const
{
  if (isLazy())
  {
    Data temp(*this);	// to get around the fact that you can't resolve a const Data
    temp.resolve();
    return temp.calc_minGlobalDataPoint(ProcNo,DataPointNo);
  }
  int i,j;
  int lowi=0,lowj=0;
  double min=numeric_limits<double>::max();

  Data temp=minval_nonlazy();	// need to do this to prevent autolazy from reintroducing laziness

  int numSamples=temp.getNumSamples();
  int numDPPSample=temp.getNumDataPointsPerSample();

  double next,local_min;
  int local_lowi=0,local_lowj=0;	

  #pragma omp parallel firstprivate(local_lowi,local_lowj) private(next,local_min)
  {
    local_min=min;
    #pragma omp for private(i,j) schedule(static)
    for (i=0; i<numSamples; i++) {
      for (j=0; j<numDPPSample; j++) {
        next=temp.getDataAtOffsetRO(temp.getDataOffset(i,j));
        if (next<local_min) {
          local_min=next;
          local_lowi=i;
          local_lowj=j;
        }
      }
    }
    #pragma omp critical
    if (local_min<min) {	// If we found a smaller value than our sentinel
      min=local_min;
      lowi=local_lowi;
      lowj=local_lowj;
    }
  }

#ifdef PASO_MPI
  // determine the processor on which the minimum occurs
  next = temp.getDataPointRO(lowi,lowj);
  int lowProc = 0;
  double *globalMins = new double[get_MPISize()+1];
  int error;
  error = MPI_Gather ( &next, 1, MPI_DOUBLE, globalMins, 1, MPI_DOUBLE, 0, get_MPIComm() );

  if( get_MPIRank()==0 ){
	next = globalMins[lowProc];
	for( i=1; i<get_MPISize(); i++ )
		if( next>globalMins[i] ){
			lowProc = i;
			next = globalMins[i];
		}
  }
  MPI_Bcast( &lowProc, 1, MPI_INT, 0, get_MPIComm() );
  DataPointNo = lowj + lowi * numDPPSample;
  MPI_Bcast(&DataPointNo, 1, MPI_INT, lowProc, get_MPIComm() );
  delete [] globalMins;
  ProcNo = lowProc;
#else
  ProcNo = 0;
  DataPointNo = lowj + lowi * numDPPSample;
#endif
}


const boost::python::tuple
Data::maxGlobalDataPoint() const
{
  int DataPointNo;
  int ProcNo;
  calc_maxGlobalDataPoint(ProcNo,DataPointNo);
  return make_tuple(ProcNo,DataPointNo);
}

void
Data::calc_maxGlobalDataPoint(int& ProcNo,
       	                int& DataPointNo) const
{
  if (isLazy())
  {
    Data temp(*this);	// to get around the fact that you can't resolve a const Data
    temp.resolve();
    return temp.calc_maxGlobalDataPoint(ProcNo,DataPointNo);
  }
  int i,j;
  int highi=0,highj=0;
//-------------
  double max= -numeric_limits<double>::max();

  Data temp=maxval_nonlazy();   // need to do this to prevent autolazy from reintroducing laziness

  int numSamples=temp.getNumSamples();
  int numDPPSample=temp.getNumDataPointsPerSample();

  double next,local_max;
  int local_highi=0,local_highj=0;	

  #pragma omp parallel firstprivate(local_highi,local_highj) private(next,local_max)
  {
    local_max=max;
    #pragma omp for private(i,j) schedule(static)
    for (i=0; i<numSamples; i++) {
      for (j=0; j<numDPPSample; j++) {
        next=temp.getDataAtOffsetRO(temp.getDataOffset(i,j));
        if (next>local_max) {
          local_max=next;
          local_highi=i;
          local_highj=j;
        }
      }
    }
    #pragma omp critical
    if (local_max>max) {	// If we found a larger value than our sentinel
      max=local_max;
      highi=local_highi;
      highj=local_highj;
    }
  }
#ifdef PASO_MPI
  // determine the processor on which the maximum occurs
  next = temp.getDataPointRO(highi,highj);
  int highProc = 0;
  double *globalMaxs = new double[get_MPISize()+1];
  int error;
  error = MPI_Gather ( &next, 1, MPI_DOUBLE, globalMaxs, 1, MPI_DOUBLE, 0, get_MPIComm() );
  if( get_MPIRank()==0 ){
    next = globalMaxs[highProc];
    for( i=1; i<get_MPISize(); i++ )
    {
	if( next<globalMaxs[i] )
	{
		highProc = i;
		next = globalMaxs[i];
	}
    }
  }
  MPI_Bcast( &highProc, 1, MPI_INT, 0, get_MPIComm() );
  DataPointNo = highj + highi * numDPPSample;  
  MPI_Bcast(&DataPointNo, 1, MPI_INT, highProc, get_MPIComm() );

  delete [] globalMaxs;
  ProcNo = highProc;
#else
  ProcNo = 0;
  DataPointNo = highj + highi * numDPPSample;
#endif
}

void
Data::saveDX(std::string fileName) const
{
  if (isEmpty())
  {
    throw DataException("Error - Operations not permitted on instances of DataEmpty.");
  }
  if (isLazy())
  {
     Data temp(*this);	// to get around the fact that you can't resolve a const Data
     temp.resolve();
     temp.saveDX(fileName);
     return;
  }
  boost::python::dict args;
  args["data"]=boost::python::object(this);
  getDomain()->saveDX(fileName,args);
  return;
}

void
Data::saveVTK(std::string fileName) const
{
  if (isEmpty())
  {
    throw DataException("Error - Operations not permitted on instances of DataEmpty.");
  }
  if (isLazy())
  {
     Data temp(*this);	// to get around the fact that you can't resolve a const Data
     temp.resolve();
     temp.saveVTK(fileName);
     return;
  }
  boost::python::dict args;
  args["data"]=boost::python::object(this);
  getDomain()->saveVTK(fileName,args,"","");
  return;
}



Data&
Data::operator+=(const Data& right)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  MAKELAZYBINSELF(right,ADD)	// for lazy + is equivalent to +=
  exclusiveWrite();			// Since Lazy data does not modify its leaves we only need to worry here
  binaryOp(right,plus<double>());
  return (*this);
}

Data&
Data::operator+=(const boost::python::object& right)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  Data tmp(right,getFunctionSpace(),false);
  (*this)+=tmp;
  return *this;
}

// Hmmm, operator= makes a deep copy but the copy constructor does not?
Data&
Data::operator=(const Data& other)
{
  m_protected=false;		// since any changes should be caught by exclusiveWrite();
//   m_data=other.m_data;
  set_m_data(other.m_data);
  return (*this);
}

Data&
Data::operator-=(const Data& right)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  MAKELAZYBINSELF(right,SUB)
  exclusiveWrite();
  binaryOp(right,minus<double>());
  return (*this);
}

Data&
Data::operator-=(const boost::python::object& right)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  Data tmp(right,getFunctionSpace(),false);
  (*this)-=tmp;
  return (*this);
}

Data&
Data::operator*=(const Data& right)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  MAKELAZYBINSELF(right,MUL)
  exclusiveWrite();
  binaryOp(right,multiplies<double>());
  return (*this);
}

Data&
Data::operator*=(const boost::python::object& right)
{  
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  Data tmp(right,getFunctionSpace(),false);
  (*this)*=tmp;
  return (*this);
}

Data&
Data::operator/=(const Data& right)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  MAKELAZYBINSELF(right,DIV)
  exclusiveWrite();
  binaryOp(right,divides<double>());
  return (*this);
}

Data&
Data::operator/=(const boost::python::object& right)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  Data tmp(right,getFunctionSpace(),false);
  (*this)/=tmp;
  return (*this);
}

Data
Data::rpowO(const boost::python::object& left) const
{
  Data left_d(left,*this);
  return left_d.powD(*this);
}

Data
Data::powO(const boost::python::object& right) const
{
  Data tmp(right,getFunctionSpace(),false);
  return powD(tmp);
}

Data
Data::powD(const Data& right) const
{
  MAKELAZYBIN(right,POW)
  return C_TensorBinaryOperation<double (*)(double, double)>(*this, right, ::pow);
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator+(const Data& left, const Data& right)
{
  MAKELAZYBIN2(left,right,ADD)
  return C_TensorBinaryOperation(left, right, plus<double>());
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator-(const Data& left, const Data& right)
{
  MAKELAZYBIN2(left,right,SUB)
  return C_TensorBinaryOperation(left, right, minus<double>());
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator*(const Data& left, const Data& right)
{
  MAKELAZYBIN2(left,right,MUL)
  return C_TensorBinaryOperation(left, right, multiplies<double>());
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator/(const Data& left, const Data& right)
{
  MAKELAZYBIN2(left,right,DIV)
  return C_TensorBinaryOperation(left, right, divides<double>());
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator+(const Data& left, const boost::python::object& right)
{
  Data tmp(right,left.getFunctionSpace(),false);
  MAKELAZYBIN2(left,tmp,ADD)
  return left+tmp;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator-(const Data& left, const boost::python::object& right)
{
  Data tmp(right,left.getFunctionSpace(),false);
  MAKELAZYBIN2(left,tmp,SUB)
  return left-tmp;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator*(const Data& left, const boost::python::object& right)
{
  Data tmp(right,left.getFunctionSpace(),false);
  MAKELAZYBIN2(left,tmp,MUL)
  return left*tmp;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator/(const Data& left, const boost::python::object& right)
{
  Data tmp(right,left.getFunctionSpace(),false);
  MAKELAZYBIN2(left,tmp,DIV)
  return left/tmp;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator+(const boost::python::object& left, const Data& right)
{
  Data tmp(left,right.getFunctionSpace(),false);
  MAKELAZYBIN2(tmp,right,ADD)
  return tmp+right;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator-(const boost::python::object& left, const Data& right)
{
  Data tmp(left,right.getFunctionSpace(),false);
  MAKELAZYBIN2(tmp,right,SUB)
  return tmp-right;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator*(const boost::python::object& left, const Data& right)
{
  Data tmp(left,right.getFunctionSpace(),false);
  MAKELAZYBIN2(tmp,right,MUL)
  return tmp*right;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator/(const boost::python::object& left, const Data& right)
{
  Data tmp(left,right.getFunctionSpace(),false);
  MAKELAZYBIN2(tmp,right,DIV)
  return tmp/right;
}


/* TODO */
/* global reduction */
Data
Data::getItem(const boost::python::object& key) const
{

  DataTypes::RegionType slice_region=DataTypes::getSliceRegion(getDataPointShape(),key);

  if (slice_region.size()!=getDataPointRank()) {
    throw DataException("Error - slice size does not match Data rank.");
  }

  return getSlice(slice_region);
}

/* TODO */
/* global reduction */
Data
Data::getSlice(const DataTypes::RegionType& region) const
{
  return Data(*this,region);
}

/* TODO */
/* global reduction */
void
Data::setItemO(const boost::python::object& key,
               const boost::python::object& value)
{
  Data tempData(value,getFunctionSpace());
  setItemD(key,tempData);
}

void
Data::setItemD(const boost::python::object& key,
               const Data& value)
{
  DataTypes::RegionType slice_region=DataTypes::getSliceRegion(getDataPointShape(),key);
  if (slice_region.size()!=getDataPointRank()) {
    throw DataException("Error - slice size does not match Data rank.");
  }
  exclusiveWrite();
  if (getFunctionSpace()!=value.getFunctionSpace()) {
     setSlice(Data(value,getFunctionSpace()),slice_region);
  } else {
     setSlice(value,slice_region);
  }
}

void
Data::setSlice(const Data& value,
               const DataTypes::RegionType& region)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  forceResolve();
  exclusiveWrite();		// In case someone finds a way to call this without going through setItemD
  Data tempValue(value);
  typeMatchLeft(tempValue);
  typeMatchRight(tempValue);
  getReady()->setSlice(tempValue.m_data.get(),region);
}

void
Data::typeMatchLeft(Data& right) const
{
  if (right.isLazy() && !isLazy())
  {
    right.resolve();
  }
  if (isExpanded()){
    right.expand();
  } else if (isTagged()) {
    if (right.isConstant()) {
      right.tag();
    }
  }
}

void
Data::typeMatchRight(const Data& right)
{
  if (isLazy() && !right.isLazy())
  {
    resolve();
  }
  if (isTagged()) {
    if (right.isExpanded()) {
      expand();
    }
  } else if (isConstant()) {
    if (right.isExpanded()) {
      expand();
    } else if (right.isTagged()) {
      tag();
    }
  }
}

// The normal TaggedValue adds the tag if it is not already present
// This form does not. It throws instead.
// This is because the names are maintained by the domain and cannot be added
// without knowing the tag number to map it to.
void
Data::setTaggedValueByName(std::string name,
                           const boost::python::object& value)
{
     if (getFunctionSpace().getDomain()->isValidTagName(name)) {
	forceResolve();
	exclusiveWrite();
        int tagKey=getFunctionSpace().getDomain()->getTag(name);
        setTaggedValue(tagKey,value);
     }
     else
     {					// The 
	throw DataException("Error - unknown tag in setTaggedValueByName.");
     }
}

void
Data::setTaggedValue(int tagKey,
                     const boost::python::object& value)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  //
  // Ensure underlying data object is of type DataTagged
  forceResolve();
  exclusiveWrite();
  if (isConstant()) tag();
  WrappedArray w(value);

  DataVector temp_data2;
  temp_data2.copyFromArray(w,1);

  m_data->setTaggedValue(tagKey,w.getShape(), temp_data2);
}


void
Data::setTaggedValueFromCPP(int tagKey,
			    const DataTypes::ShapeType& pointshape,
                            const DataTypes::ValueType& value,
			    int dataOffset)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  //
  // Ensure underlying data object is of type DataTagged
  forceResolve();
  if (isConstant()) tag();
  exclusiveWrite();
  //
  // Call DataAbstract::setTaggedValue
  m_data->setTaggedValue(tagKey,pointshape, value, dataOffset);
}

int
Data::getTagNumber(int dpno)
{
  if (isEmpty())
  {
	throw DataException("Error - operation not permitted on instances of DataEmpty.");
  }
  return getFunctionSpace().getTagFromDataPointNo(dpno);
}


ostream& escript::operator<<(ostream& o, const Data& data)
{
  o << data.toString();
  return o;
}

Data
escript::C_GeneralTensorProduct(Data& arg_0,
                     Data& arg_1,
                     int axis_offset,
                     int transpose)
{
  // General tensor product: res(SL x SR) = arg_0(SL x SM) * arg_1(SM x SR)
  // SM is the product of the last axis_offset entries in arg_0.getShape().

  // deal with any lazy data
//   if (arg_0.isLazy()) {arg_0.resolve();}
//   if (arg_1.isLazy()) {arg_1.resolve();}
  if (arg_0.isLazy() || arg_1.isLazy() || (AUTOLAZYON && (arg_0.isExpanded() || arg_1.isExpanded())))
  {
	DataLazy* c=new DataLazy(arg_0.borrowDataPtr(), arg_1.borrowDataPtr(), PROD, axis_offset,transpose);
	return Data(c);
  }

  // Interpolate if necessary and find an appropriate function space
  Data arg_0_Z, arg_1_Z;
  if (arg_0.getFunctionSpace()!=arg_1.getFunctionSpace()) {
    if (arg_0.probeInterpolation(arg_1.getFunctionSpace())) {
      arg_0_Z = arg_0.interpolate(arg_1.getFunctionSpace());
      arg_1_Z = Data(arg_1);
    }
    else if (arg_1.probeInterpolation(arg_0.getFunctionSpace())) {
      arg_1_Z=arg_1.interpolate(arg_0.getFunctionSpace());
      arg_0_Z =Data(arg_0);
    }
    else {
      throw DataException("Error - C_GeneralTensorProduct: arguments have incompatible function spaces.");
    }
  } else {
      arg_0_Z = Data(arg_0);
      arg_1_Z = Data(arg_1);
  }
  // Get rank and shape of inputs
  int rank0 = arg_0_Z.getDataPointRank();
  int rank1 = arg_1_Z.getDataPointRank();
  const DataTypes::ShapeType& shape0 = arg_0_Z.getDataPointShape();
  const DataTypes::ShapeType& shape1 = arg_1_Z.getDataPointShape();

  // Prepare for the loops of the product and verify compatibility of shapes
  int start0=0, start1=0;
  if (transpose == 0)		{}
  else if (transpose == 1)	{ start0 = axis_offset; }
  else if (transpose == 2)	{ start1 = rank1-axis_offset; }
  else				{ throw DataException("C_GeneralTensorProduct: Error - transpose should be 0, 1 or 2"); }


  // Adjust the shapes for transpose
  DataTypes::ShapeType tmpShape0(rank0);	// pre-sizing the vectors rather
  DataTypes::ShapeType tmpShape1(rank1);	// than using push_back
  for (int i=0; i<rank0; i++)	{ tmpShape0[i]=shape0[(i+start0)%rank0]; }
  for (int i=0; i<rank1; i++)	{ tmpShape1[i]=shape1[(i+start1)%rank1]; }

#if 0
  // For debugging: show shape after transpose
  char tmp[100];
  std::string shapeStr;
  shapeStr = "(";
  for (int i=0; i<rank0; i++)	{ sprintf(tmp, "%d,", tmpShape0[i]); shapeStr += tmp; }
  shapeStr += ")";
  cout << "C_GeneralTensorProduct: Shape of arg0 is " << shapeStr << endl;
  shapeStr = "(";
  for (int i=0; i<rank1; i++)	{ sprintf(tmp, "%d,", tmpShape1[i]); shapeStr += tmp; }
  shapeStr += ")";
  cout << "C_GeneralTensorProduct: Shape of arg1 is " << shapeStr << endl;
#endif

  // Prepare for the loops of the product
  int SL=1, SM=1, SR=1;
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

  // Declare output Data object
  Data res;

  if      (arg_0_Z.isConstant()   && arg_1_Z.isConstant()) {
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace());	// DataConstant output
    const double *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(0));
    const double *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(0));
    double *ptr_2 = &(res.getDataAtOffsetRW(0));
    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
  }
  else if (arg_0_Z.isConstant()   && arg_1_Z.isTagged()) {

    // Prepare the DataConstant input
    DataConstant* tmp_0=dynamic_cast<DataConstant*>(arg_0_Z.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataConstant."); }

    // Borrow DataTagged input from Data object
    DataTagged* tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());
    if (tmp_1==0) { throw DataException("GTP_1 Programming error - casting to DataTagged."); }

    // Prepare a DataTagged output 2
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace());	// DataTagged output
    res.tag();
    DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

    // Prepare offset into DataConstant
    int offset_0 = tmp_0->getPointOffset(0,0);
    const double *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0));

    const double *ptr_1 = &(tmp_1->getDefaultValueRO(0));
    double *ptr_2 = &(tmp_2->getDefaultValueRW(0));

    // Compute an MVP for the default
    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    // Compute an MVP for each tag
    const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
    DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
    for (i=lookup_1.begin();i!=lookup_1.end();i++) {
      tmp_2->addTag(i->first);

      const double *ptr_1 = &(tmp_1->getDataByTagRO(i->first,0));
      double *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0));
	
      matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    }

  }
  else if (arg_0_Z.isConstant()   && arg_1_Z.isExpanded()) {

    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
    DataConstant* tmp_0=dynamic_cast<DataConstant*>(arg_0_Z.borrowData());
    DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
    DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataConstant."); }
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    int sampleNo_1,dataPointNo_1;
    int numSamples_1 = arg_1_Z.getNumSamples();
    int numDataPointsPerSample_1 = arg_1_Z.getNumDataPointsPerSample();
    int offset_0 = tmp_0->getPointOffset(0,0);
    #pragma omp parallel for private(sampleNo_1,dataPointNo_1) schedule(static)
    for (sampleNo_1 = 0; sampleNo_1 < numSamples_1; sampleNo_1++) {
      for (dataPointNo_1 = 0; dataPointNo_1 < numDataPointsPerSample_1; dataPointNo_1++) {
        int offset_1 = tmp_1->getPointOffset(sampleNo_1,dataPointNo_1);
        int offset_2 = tmp_2->getPointOffset(sampleNo_1,dataPointNo_1);
        const double *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0));
        const double *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1));
        double *ptr_2 = &(res.getDataAtOffsetRW(offset_2));
        matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
      }
    }

  }
  else if (arg_0_Z.isTagged()     && arg_1_Z.isConstant()) {

    // Borrow DataTagged input from Data object
    DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());
    if (tmp_0==0) { throw DataException("GTP_0 Programming error - casting to DataTagged."); }

    // Prepare the DataConstant input
    DataConstant* tmp_1=dynamic_cast<DataConstant*>(arg_1_Z.borrowData());
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataConstant."); }

    // Prepare a DataTagged output 2
    res = Data(0.0, shape2, arg_0_Z.getFunctionSpace());	// DataTagged output
    res.tag();
    DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

    // Prepare offset into DataConstant
    int offset_1 = tmp_1->getPointOffset(0,0);
    const double *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1));
    const double *ptr_0 = &(tmp_0->getDefaultValueRO(0));
    double *ptr_2 = &(tmp_2->getDefaultValueRW(0));

    // Compute an MVP for the default
    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    // Compute an MVP for each tag
    const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
    DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
    for (i=lookup_0.begin();i!=lookup_0.end();i++) {

      tmp_2->addTag(i->first);
      const double *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0));
      double *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0));
      matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    }

  }
  else if (arg_0_Z.isTagged()     && arg_1_Z.isTagged()) {

    // Borrow DataTagged input from Data object
    DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

    // Borrow DataTagged input from Data object
    DataTagged* tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

    // Prepare a DataTagged output 2
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace());
    res.tag();	// DataTagged output
    DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

    const double *ptr_0 = &(tmp_0->getDefaultValueRO(0));
    const double *ptr_1 = &(tmp_1->getDefaultValueRO(0));
    double *ptr_2 = &(tmp_2->getDefaultValueRW(0));

    // Compute an MVP for the default
    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    // Merge the tags
    DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
    const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
    const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
    for (i=lookup_0.begin();i!=lookup_0.end();i++) {
      tmp_2->addTag(i->first); // use tmp_2 to get correct shape
    }
    for (i=lookup_1.begin();i!=lookup_1.end();i++) {
      tmp_2->addTag(i->first);
    }
    // Compute an MVP for each tag
    const DataTagged::DataMapType& lookup_2=tmp_2->getTagLookup();
    for (i=lookup_2.begin();i!=lookup_2.end();i++) {
      const double *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0));
      const double *ptr_1 = &(tmp_1->getDataByTagRO(i->first,0));
      double *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0));

      matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    }

  }
  else if (arg_0_Z.isTagged()     && arg_1_Z.isExpanded()) {

    // After finding a common function space above the two inputs have the same numSamples and num DPPS
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
    DataTagged*   tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());
    DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
    DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataTagged."); }
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    int sampleNo_0,dataPointNo_0;
    int numSamples_0 = arg_0_Z.getNumSamples();
    int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
    #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
    for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
      int offset_0 = tmp_0->getPointOffset(sampleNo_0,0); // They're all the same, so just use #0
      const double *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0));
      for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
        int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
        int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
        const double *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1));
        double *ptr_2 = &(res.getDataAtOffsetRW(offset_2));
        matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
      }
    }

  }
  else if (arg_0_Z.isExpanded()   && arg_1_Z.isConstant()) {

    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
    DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
    DataConstant* tmp_1=dynamic_cast<DataConstant*>(arg_1_Z.borrowData());
    DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataConstant."); }
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    int sampleNo_0,dataPointNo_0;
    int numSamples_0 = arg_0_Z.getNumSamples();
    int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
    int offset_1 = tmp_1->getPointOffset(0,0);
    #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
    for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
      for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
        int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
        int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
        const double *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0));
        const double *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1));
        double *ptr_2 = &(res.getDataAtOffsetRW(offset_2));
        matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
      }
    }


  }
  else if (arg_0_Z.isExpanded()   && arg_1_Z.isTagged()) {

    // After finding a common function space above the two inputs have the same numSamples and num DPPS
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
    DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
    DataTagged*   tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());
    DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataTagged."); }
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    int sampleNo_0,dataPointNo_0;
    int numSamples_0 = arg_0_Z.getNumSamples();
    int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
    #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
    for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
      int offset_1 = tmp_1->getPointOffset(sampleNo_0,0);
      const double *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1));
      for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
        int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
        int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
        const double *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0));
        double *ptr_2 = &(res.getDataAtOffsetRW(offset_2));
        matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
      }
    }

  }
  else if (arg_0_Z.isExpanded()   && arg_1_Z.isExpanded()) {

    // After finding a common function space above the two inputs have the same numSamples and num DPPS
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
    DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
    DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
    DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    int sampleNo_0,dataPointNo_0;
    int numSamples_0 = arg_0_Z.getNumSamples();
    int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
    #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
    for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
      for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
        int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
        int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
        int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
        const double *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0));
        const double *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1));
        double *ptr_2 = &(res.getDataAtOffsetRW(offset_2));
        matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
      }
    }

  }
  else {
    throw DataException("Error - C_GeneralTensorProduct: unknown combination of inputs");
  }

  return res;
}

DataAbstract*
Data::borrowData() const
{
  return m_data.get();
}

// Not all that happy about returning a non-const from a const
DataAbstract_ptr
Data::borrowDataPtr() const
{
  return m_data;
}

// Not all that happy about returning a non-const from a const
DataReady_ptr
Data::borrowReadyPtr() const
{
   DataReady_ptr dr=dynamic_pointer_cast<DataReady>(m_data);
   EsysAssert((dr!=0), "Error - casting to DataReady.");
   return dr;
}

std::string
Data::toString() const
{
    if (!m_data->isEmpty() &&
	!m_data->isLazy() && 
	getLength()>escriptParams.getInt("TOO_MANY_LINES"))
    {
	stringstream temp;
	temp << "Summary: inf="<< inf_const() << " sup=" << sup_const() << " data points=" << getNumDataPoints();
	return  temp.str();
    }
    return m_data->toString();
}


// This method is not thread-safe
DataTypes::ValueType::reference
Data::getDataAtOffsetRW(DataTypes::ValueType::size_type i)
{
    checkExclusiveWrite();
    return getReady()->getDataAtOffsetRW(i);
}

// This method is not thread-safe
DataTypes::ValueType::const_reference
Data::getDataAtOffsetRO(DataTypes::ValueType::size_type i)
{
    forceResolve();
    return getReady()->getDataAtOffsetRO(i);
}


// DataTypes::ValueType::const_reference
// Data::getDataAtOffsetRO(DataTypes::ValueType::size_type i) const
// {
//     if (isLazy())
//     {
// 	throw DataException("Programmer error - getDataAtOffsetRO() not permitted on Lazy Data (object is const which prevents resolving).");
//     }
//     return getReady()->getDataAtOffsetRO(i);
// }


DataTypes::ValueType::const_reference
Data::getDataPointRO(int sampleNo, int dataPointNo)
{
  forceResolve();
  if (!isReady())
  {
	throw DataException("Programmer error -getDataPointRO() not permitted on Lazy Data.");
  }
  else
  {
	const DataReady* dr=getReady();
	return dr->getDataAtOffsetRO(dr->getPointOffset(sampleNo, dataPointNo));
  }
}




DataTypes::ValueType::reference
Data::getDataPointRW(int sampleNo, int dataPointNo)
{
  checkExclusiveWrite();
  DataReady* dr=getReady();
  return dr->getDataAtOffsetRW(dr->getPointOffset(sampleNo, dataPointNo));
}

BufferGroup* 
Data::allocSampleBuffer() const
{
     if (isLazy())
     {
	#ifdef _OPENMP
	int tnum=omp_get_max_threads();
	#else
	int tnum=1;
	#endif
	return new BufferGroup(getSampleBufferSize(),tnum);
     }
     else
     {
	return NULL;
     }
}

void
Data::freeSampleBuffer(BufferGroup* bufferg)
{
     if (bufferg!=0)
     {
	delete bufferg;
     }
}


Data
Data::interpolateFromTable2DP(boost::python::object table, double Amin, double Astep,
		Data& B, double Bmin, double Bstep, double undef, bool check_boundaries)
{
    WrappedArray t(table);
    return interpolateFromTable2D(t, Amin, Astep, undef, B, Bmin, Bstep,check_boundaries);
}

Data
Data::interpolateFromTable1DP(boost::python::object table, double Amin, double Astep,
					      double undef,bool check_boundaries)
{
	WrappedArray t(table);
	return interpolateFromTable1D(t, Amin, Astep, undef, check_boundaries);
}


Data
Data::interpolateFromTable1D(const WrappedArray& table, double Amin, double Astep,
			     double undef, bool check_boundaries)
{
	table.convertArray();		// critical!  Calling getElt on an unconverted array is not thread safe
	int error=0;
	if ((getDataPointRank()!=0))
	{
		throw DataException("Input to 1D interpolation must be scalar");
	}
	if (table.getRank()!=1)
	{
		throw DataException("Table for 1D interpolation must be 1D");
	}
	if (Astep<=0)
	{
		throw DataException("Astep must be positive");
	}
	if (!isExpanded())
	{
		expand();
	}
	Data res(0, DataTypes::scalarShape, getFunctionSpace(), true);
	try
	{
		int numpts=getNumDataPoints();
		const DataVector& adat=getReady()->getVectorRO();
		DataVector& rdat=res.getReady()->getVectorRW();
		int twidth=table.getShape()[0]-1;
		bool haserror=false;
		int l=0;
		#pragma omp parallel for private(l) schedule(static)
		for (l=0;l<numpts; ++l)
		{
		  #pragma omp flush(haserror)		// In case haserror was in register
		  if (!haserror)		
		  {
		   int lerror=0;
		   try
		   {
		     do			// so we can use break
		     {
			double a=adat[l];
			int x=static_cast<int>(((a-Amin)/Astep));
			if (check_boundaries) {
				if ((a<Amin) || (x<0))
				{
					lerror=1;
					break;	
				}
				if (a>Amin+Astep*twidth) 
				{
				        lerror=4;
					break;
				}
			} 
			if (x<0) x=0;
			if (x>twidth) x=twidth;

			if (x==twidth)			// value is on the far end of the table
			{
				double e=table.getElt(x);
				if (e>undef)
				{
					lerror=2;
					break; 
				}
				rdat[l]=e;
			}
			else		// x and y are in bounds
			{
				double e=table.getElt(x);
				double w=table.getElt(x+1);
				if ((e>undef) || (w>undef))
				{
					lerror=2;
					break; 
				}
		// map x*Astep <= a << (x+1)*Astep to [-1,1] 
				double la = 2.0*(a-Amin-(x*Astep))/Astep-1;
				rdat[l]=((1-la)*e + (1+la)*w)/2;
			}	
		      } while (false);
		    } catch (DataException d)
		    {
			    lerror=3;
		    }
		    if (lerror!=0)
		    {
			#pragma omp critical	// Doco says there is a flush associated with critical
			{
			    haserror=true;	// We only care that one error is recorded. We don't care which 
			    error=lerror;	// one
			}
		    }
		  } // if (!error)
		}	// parallelised for
	} catch (DataException d)
	{
		error=3;		// this is outside the parallel region so assign directly
	}
#ifdef PASO_MPI
	int rerror=0;
	MPI_Allreduce( &error, &rerror, 1, MPI_INT, MPI_MAX, get_MPIComm() );
	error=rerror;
#endif
	if (error)
	{
		switch (error)
		{
			case 1: throw DataException("Value below lower table range.");
			case 2: throw DataException("Interpolated value too large");
			case 4: throw DataException("Value greater than upper table range.");
			default:
				throw DataException("Unknown error in interpolation");		
		}
	}
	return res;
}

		
Data
Data::interpolateFromTable2D(const WrappedArray& table, double Amin, double Astep,
                       double undef, Data& B, double Bmin, double Bstep, bool check_boundaries)
{
    table.convertArray();		// critical!   Calling getElt on an unconverted array is not thread safe
    int error=0;
    if ((getDataPointRank()!=0) || (B.getDataPointRank()!=0))
    {
        throw DataException("Inputs to 2D interpolation must be scalar");
    }
    if (table.getRank()!=2)
    {
	throw DataException("Table for 2D interpolation must be 2D");
    }
    if ((Astep<=0) || (Bstep<=0))
    {
	throw DataException("Astep and Bstep must be postive");
    }
    if (getFunctionSpace()!=B.getFunctionSpace())
    {
	Data n=B.interpolate(getFunctionSpace());
	return interpolateFromTable2D(table, Amin, Astep, undef, 
		n , Bmin, Bstep, check_boundaries);
    }
    if (!isExpanded())
    {
	expand();
    }
    if (!B.isExpanded())
    {
	B.expand();
    }
    Data res(0, DataTypes::scalarShape, getFunctionSpace(), true);
    try
    {
	int numpts=getNumDataPoints();
	const DataVector& adat=getReady()->getVectorRO();
	const DataVector& bdat=B.getReady()->getVectorRO();
	DataVector& rdat=res.getReady()->getVectorRW();
	const DataTypes::ShapeType& ts=table.getShape();
	int twx=ts[0]-1;	// table width x
	int twy=ts[1]-1;	// table width y
	bool haserror=false;
	int l=0;
	#pragma omp parallel for private(l) schedule(static) 
	for (l=0; l<numpts; ++l)
	{
	   #pragma omp flush(haserror)		// In case haserror was in register
	   if (!haserror)		
	   {
	     int lerror=0;
	     try
	     {
	       do
	       {
		double a=adat[l];
		double b=bdat[l];
		int x=static_cast<int>(((a-Amin)/Astep));
		int y=static_cast<int>(((b-Bmin)/Bstep));
                if (check_boundaries) {
			    if ( (a<Amin) || (b<Bmin) || (x<0) || (y<0) )
			    {
				#pragma omp critical
				{
				    lerror=1;
				}
				break;	
			    }
			    if ( (a>Amin+Astep*twx) || (b>Bmin+Bstep*twy) )
			    {
				#pragma omp critical
				{
				    lerror=4;
				}
				break;
			    }
                } 
                if (x<0) x=0;
                if (y<0) y=0;
                if (x>twx) x=twx;
                if (y>twx) y=twy;

                if (x == twx ) {
                     if (y == twy ) {
		         double sw=table.getElt(x,y);
		         if ((sw>undef))
		         {
			     #pragma omp critical
			     {
			         lerror=2;
			     }
			     break; 
		         }
		         rdat[l]=sw;

                     } else {
		         double sw=table.getElt(x,y);
		         double nw=table.getElt(x,y+1);
		         if ((sw>undef) || (nw>undef))
		         {
			     #pragma omp critical
			     {
			        lerror=2;
			     }
			     break; 
		         }
		         double lb = 2.0*(b-Bmin-(y*Bstep))/Bstep-1;
		         rdat[l]=((1-lb)*sw + (1+lb)*nw )/2.;

                     }
                } else {
                     if (y == twy ) {
		         double sw=table.getElt(x,y);
		         double se=table.getElt(x+1,y);
		         if ((sw>undef) || (se>undef) )
		         {
			     #pragma omp critical
			     {
			     	lerror=2;
			     }
			     break; 
		         }
		         double la = 2.0*(a-Amin-(x*Astep))/Astep-1;
		         rdat[l]=((1-la)*sw + (1+la)*se )/2;

                     } else {
		         double sw=table.getElt(x,y);
		         double nw=table.getElt(x,y+1);
		         double se=table.getElt(x+1,y);
		         double ne=table.getElt(x+1,y+1);
		         if ((sw>undef) || (nw>undef) || (se>undef) || (ne>undef))
		         {
			    #pragma omp critical
			    {
			         lerror=2;
			    }
			    break; 
		         }
		         // map x*Astep <= a << (x+1)*Astep to [-1,1] 
		         // same with b
		         double la = 2.0*(a-Amin-(x*Astep))/Astep-1;
		         double lb = 2.0*(b-Bmin-(y*Bstep))/Bstep-1;
		         rdat[l]=((1-la)*(1-lb)*sw + (1-la)*(1+lb)*nw +
			          (1+la)*(1-lb)*se + (1+la)*(1+lb)*ne)/4;
                     }
                }
	       } while (false);
	      } catch (DataException d)
	      {
		lerror=3;
	      }
	      if (lerror!=0)
	      {
			#pragma omp critical	// Doco says there is a flush associated with critical
			{
			    haserror=true;	// We only care that one error is recorded. We don't care which 
			    error=lerror;	// one
			}		
	      }
	    }  // if (!error)
	}	// parallel for
    } catch (DataException d)
    {
	    error=3;
    }
#ifdef PASO_MPI
    int rerror=0;
    MPI_Allreduce( &error, &rerror, 1, MPI_INT, MPI_MAX, get_MPIComm() );
    error=rerror;
#endif
    if (error)
    {
	switch (error)
	{
			case 1: throw DataException("Value below lower table range.");
			case 2: throw DataException("Interpolated value too large");
			case 4: throw DataException("Value greater than upper table range.");
	                default:
		              throw DataException("Unknown error in interpolation");		
	}
    }
    return res;
}


/* Member functions specific to the MPI implementation */

void
Data::print()
{
  int i,j;

  printf( "Data is %dX%d\n", getNumSamples(), getNumDataPointsPerSample() );
  for( i=0; i<getNumSamples(); i++ )
  {
    printf( "[%6d]", i );
    for( j=0; j<getNumDataPointsPerSample(); j++ )
      printf( "\t%10.7g", (getSampleDataRW(i))[j] );	// doesn't really need RW access
    printf( "\n" );
  }
}
void
Data::dump(const std::string fileName) const
{
  try
  {
	if (isLazy())
	{
	  Data temp(*this);	// this is to get a non-const object which we can resolve
	  temp.resolve();
	  temp.dump(fileName);
	}
	else
	{
          return m_data->dump(fileName);
	}
  }
  catch (std::exception& e)
  {
        cout << e.what() << endl;
  }
}

int
Data::get_MPISize() const
{
	int size;
#ifdef PASO_MPI
	int error;
	error = MPI_Comm_size( get_MPIComm(), &size );
#else
	size = 1;
#endif
	return size;
}

int
Data::get_MPIRank() const
{
	int rank;
#ifdef PASO_MPI
	int error;
	error = MPI_Comm_rank( get_MPIComm(), &rank );
#else
	rank = 0;
#endif
	return rank;
}

MPI_Comm
Data::get_MPIComm() const
{
#ifdef PASO_MPI
	return MPI_COMM_WORLD;
#else
	return -1;
#endif
}

