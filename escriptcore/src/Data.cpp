
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "Data.h"

#include "AbstractContinuousDomain.h"
#include "DataConstant.h"
#include "DataEmpty.h"
#include "DataExpanded.h"
#include "DataLazy.h"
#include "DataTagged.h"
#include "EscriptParams.h"
#include "FunctionSpaceException.h"
#include "FunctionSpaceFactory.h"
#include "BinaryDataReadyOps.h"

#ifdef IKNOWWHATIMDOING
#include "Dodgy.h"
#endif

#include <algorithm>
#include <fstream>
#include <functional>
#include <sstream>      // so we can throw messages about ranks
#include <vector>
#include <iostream>

#include <boost/python/dict.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/long.hpp>
#include "WrappedArray.h"

namespace bp = boost::python;
using namespace escript;
using namespace escript::DataTypes;
using namespace std;
using DataTypes::real_t;
using DataTypes::cplx_t;


#define THROWONCOMPLEX if (m_data->isComplex()){throw DataException("Operation does not support complex objects");}
#define THROWONCOMPLEXA(Z) if (Z.isComplex()){throw DataException("Operation does not support complex objects");}

// ensure the current object is not a DataLazy
// The idea was that we could add an optional warning whenever a resolve is forced
// #define forceResolve() if (isLazy()) {#resolve();}

#define AUTOLAZYON escriptParams.getAutoLazy()
#define MAKELAZYOP(X) do {\
  if (isLazy() || (AUTOLAZYON && m_data->isExpanded())) \
  {\
        DataLazy* c=new DataLazy(borrowDataPtr(),X);\
        return Data(c);\
  }\
}while(0)

#define MAKELAZYOPOFF(X,Y) do {\
  if (isLazy() || (AUTOLAZYON && m_data->isExpanded())) \
  {\
        DataLazy* c=new DataLazy(borrowDataPtr(),X,Y);\
        return Data(c);\
  }\
}while(0)

#define MAKELAZYOP2(X,Y,Z) do {\
  if (isLazy() || (AUTOLAZYON && m_data->isExpanded())) \
  {\
        DataLazy* c=new DataLazy(borrowDataPtr(),X,Y,Z);\
        return Data(c);\
  }\
}while(0)

#define MAKELAZYBINSELF(R,X) do {\
  if (isLazy() || R.isLazy() || (AUTOLAZYON && (isExpanded() || R.isExpanded()))) \
  {\
        DataLazy* c=new DataLazy(m_data,R.borrowDataPtr(),X);\
/*         m_data=c->getPtr();*/     set_m_data(c->getPtr());\
        return (*this);\
  }\
}while(0)

// like the above but returns a new data rather than *this
#define MAKELAZYBIN(R,X) do {\
  if (isLazy() || R.isLazy() || (AUTOLAZYON && (isExpanded() || R.isExpanded()))) \
  {\
        DataLazy* c=new DataLazy(m_data,R.borrowDataPtr(),X);\
        return Data(c);\
  }\
}while(0)

#define MAKELAZYBIN2(L,R,X) do {\
  if (L.isLazy() || R.isLazy() || (AUTOLAZYON && (L.isExpanded() || R.isExpanded()))) \
  {\
/*  if (L.isComplex() || R.isComplex()) \
  {\
      throw DataException("Lazy operations on complex not supported yet");\
  }\
*/        DataLazy* c=new DataLazy(L.borrowDataPtr(),R.borrowDataPtr(),X);\
        return Data(c);\
  }\
}while(0)

#define CHECK_DO_CRES escriptParams.getResolveCollective()


namespace
{

template <class ARR>
inline
bp::tuple
pointToTuple1(const DataTypes::ShapeType& shape, ARR v, unsigned long offset)
{
    bp::list l;
    unsigned int dim=shape[0];
    for (size_t i=0;i<dim;++i)
    {
        l.append(v[i+offset]);
    }
    return bp::tuple(l);
}

template <class ARR>
inline
bp::tuple
pointToTuple2(const DataTypes::ShapeType& shape, ARR v, unsigned long offset)
{
    unsigned int shape0=shape[0];
    unsigned int shape1=shape[1];
    bp::list lj;
    for (size_t j=0;j<shape0;++j)
    {
        bp::list li;
        for (size_t i=0;i<shape1;++i)
        {
            li.append(v[offset+DataTypes::getRelIndex(shape,j,i)]);
        }
        lj.append(bp::tuple(li));
    }
    return bp::tuple(lj);
}

template <class ARR>
inline
bp::tuple
pointToTuple3(const DataTypes::ShapeType& shape, ARR v, unsigned long offset)
{
    unsigned int shape0=shape[0];
    unsigned int shape1=shape[1];
    unsigned int shape2=shape[2];

    bp::list lk;
    for (size_t k=0;k<shape0;++k)
    {
        bp::list lj;
        for (size_t j=0;j<shape1;++j)
        {
            bp::list li;
            for (size_t i=0;i<shape2;++i)
            {
                li.append(v[offset+DataTypes::getRelIndex(shape,k,j,i)]);
            }
            lj.append(bp::tuple(li));
        }
        lk.append(bp::tuple(lj));
    }
    return bp::tuple(lk);
}

template <class ARR>
inline
bp::tuple
pointToTuple4(const DataTypes::ShapeType& shape, ARR v, unsigned long offset)
{
    unsigned int shape0=shape[0];
    unsigned int shape1=shape[1];
    unsigned int shape2=shape[2];
    unsigned int shape3=shape[3];

    bp::list ll;
    for (size_t l=0;l<shape0;++l)
    {
        bp::list lk;
        for (size_t k=0;k<shape1;++k)
        {
            bp::list lj;
            for (size_t j=0;j<shape2;++j)
            {
                bp::list li;
                for (size_t i=0;i<shape3;++i)
                {
                    li.append(v[offset+DataTypes::getRelIndex(shape,l,k,j,i)]);
                }
                lj.append(bp::tuple(li));
            }
            lk.append(bp::tuple(lj));
        }
        ll.append(bp::tuple(lk));
    }
    return bp::tuple(ll);
}


// This should be safer once the DataC RO changes have been brought in
template <class ARR>
bp::tuple
pointToTuple( const DataTypes::ShapeType& shape,ARR v)
{
    int rank=shape.size();
    if (rank==0)
    {
        return bp::make_tuple(v[0]);
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
    :  m_lazy(false)
{
    //
    // Default data is type DataEmpty
    DataAbstract* temp=new DataEmpty();
    set_m_data(temp->getPtr());
    m_protected=false;
}



Data::Data(real_t value,
           const DataTypes::ShapeType& dataPointShape,
           const FunctionSpace& what,
           bool expanded)
        :  m_lazy(false)
{
    initialise(value, dataPointShape, what, expanded);
    m_protected=false;
}

Data::Data(cplx_t value,
           const DataTypes::ShapeType& dataPointShape,
           const FunctionSpace& what,
           bool expanded)
        :  m_lazy(false)
{
    initialise(value, dataPointShape, what, expanded);
    m_protected=false;
}

// Ordering is: shape, functionspace, expanded
Data::Data(boost::python::object value,
       boost::python::object par2,
       boost::python::object par3,
       boost::python::object par4)
{
    if (value.is_none())
    {
	throw DataException("Data constructor from python - first argument must not be None.");
    }
      // now to enforce contiguous Nones
    if ((par2.is_none() && (!par3.is_none() || !par4.is_none())) ||
        (par3.is_none() && !par4.is_none()))
    {
	throw DataException("Data constructor from python - arguments must be omitted from the right.");
    }
  
    // what is the first arg
    boost::python::extract<DataTypes::cplx_t> exc(value);
    boost::python::extract<DataTypes::real_t> exr(value);
    boost::python::extract<Data> exd(value);
    if (exc.check() || exr.check())
    {
	//(value, shape, functionspace, expanded), but shape could be missing in which case => ()
        DataTypes::ShapeType dataPointShape;	// default to scalar
	if (!par2.is_none())
	{
	    boost::python::extract<boost::python::tuple> ex2t(par2);
	    boost::python::extract<FunctionSpace> ex2fs(par2);
	    if (ex2t.check())
	    {
		for (int i = 0; i < par2.attr("__len__")(); ++i) {
		    dataPointShape.push_back(bp::extract<const int>(par2[i]));
		}	    
	    }
	    else if (ex2fs.check())
	    {
		// shape will default to ()
	        // now we copy the params up to make it look like they did specify one
	        par4=par3;
		par3=par2;
	    }
	    else
	    {
		throw DataException("Data constructor from python - expected a tuple or None as second argument.");
	    }
	}
	boost::python::extract<FunctionSpace> ex3fs(par3);	
	if (!par3.is_none())
	{
	    if (!ex3fs.check())
	    {
		throw DataException("Data constructor from python - expected a FunctionSpace or None as third argument.");
	    }
	}
	bool expa=false;
	if (!par4.is_none())
	{
	    boost::python::extract<bool> ex4b(par4);
	    if (!ex4b.check())
	    {
		throw DataException("Data constructor from python - expected a boolean or None as fourth argument.");	      
	    }
	    expa=ex4b();
	}
	if (exr.check())
	{
	    int len = DataTypes::noValues(dataPointShape);
	    RealVectorType temp_data(len,exr(),len);
	    initialise(temp_data, dataPointShape, par3.is_none()?FunctionSpace():ex3fs(), expa);
	    m_protected=false;	  	  
	}
	else
	{
	    int len = DataTypes::noValues(dataPointShape);
	    CplxVectorType temp_data(len,exc(),len);
	    initialise(temp_data, dataPointShape, par3.is_none()?FunctionSpace():ex3fs(), expa);
	    m_protected=false;	  
	}
    }
    else if (exd.check())	// Construct from (Data, [FunctionSpace]
    {
	boost::python::extract<FunctionSpace> ex2fs(par2);	      
	if (!par2.is_none())
	{
	    if (!ex2fs.check())
	    {
		throw DataException("Data constructor from python - expected a FunctionSpace or None as second argument.");	      
	    }
	}        
        init_from_data_and_fs(exd(), par2.is_none()?FunctionSpace():ex2fs());
    }
    else	//non-data, non-scalar first argument
    {   //               2         3
	//(value, functionspace, expanded)      
        if (!par4.is_none())
	{
	    throw DataException("Data constructor from python - unexpected fourth argument.");
	}
	bool expa=false;
	if (!par3.is_none())
	{
	    boost::python::extract<bool> ex3b(par3);
	    if (!ex3b.check())
	    {
		throw DataException("Data constructor from python - expected a boolean or None as third argument.");	      
	    }
	    expa=ex3b();
	}
	boost::python::extract<FunctionSpace> ex2fs(par2);	 
	if (!par2.is_none())
	{
	    if (!ex2fs.check())
	    {
		throw DataException("Data constructor from python - expected a FunctionSpace or None as second argument.");	      
	    }
	}
	WrappedArray w(value);
	initialise(w,par2.is_none()?FunctionSpace():ex2fs(),expa);
	m_protected=false;
    }
}

Data::Data(const Data& inData)
        :  m_lazy(false)
{
    set_m_data(inData.m_data);
    m_protected=inData.isProtected();
}


Data::Data(const Data& inData,
           const DataTypes::RegionType& region)
        :  m_lazy(false)
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

void Data::init_from_data_and_fs(const Data& inData,
           const FunctionSpace& functionspace)
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
        if (inData.isConstant()) 
        {  // for a constant function, we just need to use the new function space
            if (!inData.probeInterpolation(functionspace))
            { // Even though this is constant, we still need to check whether interpolation is allowed
                throw FunctionSpaceException("Cannot interpolate across to the domain of the specified FunctionSpace. (DataConstant)");
            }
            // if the data is not lazy, this will just be a cast to DataReady
            DataReady_ptr dr=inData.m_data->resolve();
            DataConstant* dc=0;
            if (inData.isComplex())
            {
                dc=new DataConstant(functionspace,inData.m_data->getShape(),dr->getTypedVectorRO((DataTypes::cplx_t)0));     	      
            }
            else
            {
                dc=new DataConstant(functionspace,inData.m_data->getShape(),dr->getTypedVectorRO((DataTypes::real_t)0));     	      
            }
            set_m_data(DataAbstract_ptr(dc));
        } 
        else 
        {
            if (inData.isComplex())
            {
                Data tmp((DataTypes::cplx_t)0,inData.getDataPointShape(),functionspace,true);
                const_Domain_ptr inDataDomain=inData.getDomain();
                if  (inDataDomain==functionspace.getDomain()) {
                    if (inData.isLazy())
                    {
                        Data temp(inData);
                        temp.resolve();     // since complex lazy is not a thing
                        inDataDomain->interpolateOnDomain(tmp,temp);
                    }
                    else
                    {
                        inDataDomain->interpolateOnDomain(tmp,inData);
                    }
                } 
                else 
                {
                    if (inData.isLazy())
                    {
                        Data temp(inData);
                        temp.resolve();
                        inDataDomain->interpolateAcross(tmp,temp);
                    }
                    else
                    {
                        inDataDomain->interpolateAcross(tmp,inData);
                    }
                }
                set_m_data(tmp.m_data);	      	      
            }
            else
            {
                Data tmp(0,inData.getDataPointShape(),functionspace,true);
                const_Domain_ptr inDataDomain=inData.getDomain();
                if  (inDataDomain==functionspace.getDomain()) {
                    inDataDomain->interpolateOnDomain(tmp,inData);
                } else {
                    inDataDomain->interpolateAcross(tmp,inData);
                }
                set_m_data(tmp.m_data);
            }
        }
    }
    m_protected=false;
}


Data::Data(const Data& inData,
           const FunctionSpace& functionspace)
        :  m_lazy(false)
{
    init_from_data_and_fs(inData, functionspace);  
}

Data::Data(DataAbstract* underlyingdata)
        :  m_lazy(false)
{
    set_m_data(underlyingdata->getPtr());
    m_protected=false;
}

Data::Data(DataAbstract_ptr underlyingdata)
        :  m_lazy(false)
{
    set_m_data(underlyingdata);
    m_protected=false;
}

Data::Data(const DataTypes::RealVectorType& value,
           const DataTypes::ShapeType& shape,
           const FunctionSpace& what,
           bool expanded)
        :  m_lazy(false)
{
    initialise(value,shape,what,expanded);
    m_protected=false;
}




Data::Data(const WrappedArray& w, const FunctionSpace& what,
           bool expanded)
           : m_lazy(false)
{
    initialise(w,what,expanded);  
    m_protected=false; 
}


Data::Data(const boost::python::object& value,
           const Data& other)
        :  m_lazy(false)
{
    WrappedArray w(value);

    // extract the shape of the array
    const DataTypes::ShapeType& tempShape=w.getShape();
    if (w.getRank()==0) {

        if (w.isComplex()==true)
        {
            // get the space for the data vector
            int len1 = DataTypes::noValues(tempShape);
            CplxVectorType temp_data(len1, 0.0, len1);
            temp_data.copyFromArray(w,1);

            int len = DataTypes::noValues(other.getDataPointShape());

            CplxVectorType temp2_data(len, temp_data[0], len);
            DataConstant* t=new DataConstant(other.getFunctionSpace(),other.getDataPointShape(),temp2_data);
            set_m_data(DataAbstract_ptr(t));
        }
        else
        {
            // get the space for the data vector
            int len1 = DataTypes::noValues(tempShape);
            RealVectorType temp_data(len1, 0.0, len1);
            temp_data.copyFromArray(w,1);

            int len = DataTypes::noValues(other.getDataPointShape());

            RealVectorType temp2_data(len, temp_data[0], len);
            DataConstant* t=new DataConstant(other.getFunctionSpace(),other.getDataPointShape(),temp2_data);
            set_m_data(DataAbstract_ptr(t));
        }
    } else {
        //
        // Create a DataConstant with the same sample shape as other
        DataConstant* t=new DataConstant(w,other.getFunctionSpace());
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
    if (p.get()!=0)
    {
        m_data=p;
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
        set_m_data(temp->getPtr());
    } else {
        DataAbstract* temp=new DataConstant(value, what);
        set_m_data(temp->getPtr());
    }
}


void
Data::initialise(const DataTypes::RealVectorType& value,
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
        set_m_data(temp->getPtr());
    } else {
        DataAbstract* temp=new DataConstant(what, shape, value);
        set_m_data(temp->getPtr());
    }
}

void
Data::initialise(const DataTypes::CplxVectorType& value,
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
        set_m_data(temp->getPtr());
    } else {
        DataAbstract* temp=new DataConstant(what, shape, value);
        set_m_data(temp->getPtr());
    }
}


void
Data::initialise(const real_t value,
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
        DataAbstract_ptr p(temp);
        set_m_data(p);
    } else {
        DataAbstract* temp=new DataConstant(what, shape, value);
        DataAbstract_ptr p(temp);
        set_m_data(p);
    }
}


void
Data::initialise(const cplx_t value,
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
        DataAbstract_ptr p(temp);
        set_m_data(p);
    } else {
        DataAbstract* temp=new DataConstant(what, shape, value);
        DataAbstract_ptr p(temp);
        set_m_data(p);
    }
}

const bp::tuple
Data::getShapeTuple() const
{
    const DataTypes::ShapeType& shape=getDataPointShape();
    switch(getDataPointRank()) {
        case 0:
            return bp::make_tuple();
        case 1:
            return bp::make_tuple(bp::long_(shape[0]));
        case 2:
            return bp::make_tuple(bp::long_(shape[0]),bp::long_(shape[1]));
        case 3:
            return bp::make_tuple(bp::long_(shape[0]),bp::long_(shape[1]),bp::long_(shape[2]));
        case 4:
            return bp::make_tuple(bp::long_(shape[0]),bp::long_(shape[1]),bp::long_(shape[2]),bp::long_(shape[3]));
        default:
            throw DataException("Error - illegal Data rank.");
    }
}


long
Data::getShapeProduct() const
{
    const DataTypes::ShapeType& shape=getDataPointShape();
    switch(getDataPointRank()) {
        case 0:
            return 1;
        case 1:
            return shape[0];
        case 2:
            return shape[0]*shape[1];
        case 3:
            return shape[0]*shape[1]*shape[2];
        case 4:
            return shape[0]*shape[1]*shape[2]*shape[3];
        default:
            throw DataException("Error - illegal Data rank.");
    }
}


// The different name is needed because boost has trouble with overloaded functions.
// It can't work out what type the function is based solely on its name.
// There are ways to fix this involving creating function pointer variables for each form
// but there doesn't seem to be a need given that the methods have the same name from the python point of view
Data
Data::copySelf() const
{
    DataAbstract* temp=m_data->deepCopy();
    return Data(temp);
}

void
Data::copy(const Data& other)
{
    DataAbstract* temp=other.m_data->deepCopy();
    DataAbstract_ptr p=temp->getPtr();
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
        throw DataException("Error - Operations (setToZero)  permitted on instances of DataEmpty.");
    }
    if (isLazy())
    {
	// This can't actually be complex yet but Just putting in this for later
        if (isComplex())
	{
	    throw DataException("Programmer Error - setToZero is not supported on lazy complex values.");
	}
        DataTypes::RealVectorType v(getNoValues(),0);
        DataConstant* dc=new DataConstant(getFunctionSpace(),getDataPointShape(),v);
        DataLazy* dl=new DataLazy(dc->getPtr());
        set_m_data(dl->getPtr());
    }
    else
    {
	// we don't want to call exclusiveWrite() here because
        // as soon as we get the copy we'd overwrite it
        if (isShared())
        {	  
                DataAbstract* t=m_data->zeroedCopy();
                set_m_data(DataAbstract_ptr(t));
        }	
	else
	{
	    m_data->setToZero();
	}
    }
}


void
Data::copyWithMask(const Data& other,
                   const Data& mask)
{
    // 1. Interpolate if required so all Data use the same FS as this
    // 2. Tag or Expand so that all Data's are the same type
    // 3. Iterate over the data vectors copying values where mask is >0
    if (other.isEmpty() || mask.isEmpty())
    {
        throw DataException("Error - copyWithMask not permitted using instances of DataEmpty.");
    }
    if (mask.isComplex())
    {
        throw DataException("Error - copyWithMask not permitted using a complex mask.");
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
    if ((selfrank>0) && (otherrank==0) &&(maskrank==0))
    {
        // Not allowing this combination.
        // it is not clear what the rank of the target should be.
        // Should it be filled with the scalar (rank stays the same); 
        // or should the target object be reshaped to be a scalar as well.
        throw DataException("Attempt to copyWithMask from scalar mask and data into non-scalar target.");
    }
    
    if (isComplex()!=other2.isComplex())
    {
        complicate();
        other2.complicate();
    }
    
    exclusiveWrite();
    if (!isComplex())
    {
	maskWorker(other2, mask2, DataTypes::real_t(0));
    }
    else
    {
	maskWorker(other2, mask2, DataTypes::cplx_t(0));
    }
}

template <typename S>
void 
Data::maskWorker(Data& other2, Data& mask2, S sentinel)
{
     // Now we iterate over the elements
    auto& self=getReady()->getTypedVectorRW(sentinel);
    auto& ovec=other2.getReadyPtr()->getTypedVectorRO(sentinel);
    auto& mvec=mask2.getReadyPtr()->getTypedVectorRO(DataTypes::real_t(0));

    unsigned int selfrank=getDataPointRank();
    unsigned int otherrank=other2.getDataPointRank();
    unsigned int maskrank=mask2.getDataPointRank();    

    if ((selfrank>0) && (otherrank>0) &&(maskrank==0))
    {
        if (mvec[0]>0)          // copy whole object if scalar is >0
        {
            copy(other2);
        }
        return;
    }
    if (isTagged())               // so all objects involved will also be tagged
    {
        // note the !
        if (!((getDataPointShape()==mask2.getDataPointShape()) && 
                ((other2.getDataPointShape()==mask2.getDataPointShape()) || (otherrank==0))))
        {
            throw DataException("copyWithMask, shape mismatch.");
        }

        // We need to consider the possibility that tags are missing or in the wrong order
        // My guiding assumption here is: All tagged Data are assumed to have the default value for
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
	        // yes this is explicitly RealType but size_t should not vary
	        DataTypes::RealVectorType::size_type toff=tptr->getOffsetForTag(i->first);
                DataTypes::RealVectorType::size_type moff=mptr->getOffsetForTag(i->first);
                DataTypes::RealVectorType::size_type ooff=optr->getOffsetForTag(i->first);
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
        else    // other is a scalar
        {
            for (i=tlookup.begin();i!=tlookup.end();i++)
            {
                // get the target offset
	        DataTypes::RealVectorType::size_type toff=tptr->getOffsetForTag(i->first);
                DataTypes::RealVectorType::size_type moff=mptr->getOffsetForTag(i->first);
                DataTypes::RealVectorType::size_type ooff=optr->getOffsetForTag(i->first);	      
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

        return;  // ugly
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
                self[i]=ovec[i/psize]; // since this is expanded there is one scalar 
            }                          // dest point
        } 
        return; // ugly!
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
    return m_lazy;  // not asking m_data because we need to be able to ask this while m_data is changing
}

// at the moment this is synonymous with !isLazy() but that could change
bool
Data::isReady() const
{
    return (dynamic_cast<DataReady*>(m_data.get())!=0);
}


bool
Data::isComplex() const
{
    return m_data->isComplex();
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
        set_m_data(temp->getPtr());
    } else if (isTagged()) {
        DataTagged* tempDataTag=dynamic_cast<DataTagged*>(m_data.get());
        DataAbstract* temp=new DataExpanded(*tempDataTag);
        set_m_data(temp->getPtr());
    } else if (isExpanded()) {
        //
        // do nothing
    } else if (isEmpty()) {
        throw DataException("Error - Expansion of DataEmpty not possible.");
    } else if (isLazy()) {
        resolve();
        expand();   // resolve might not give us expanded data
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
    MAKELAZYOP(RECIP);
    return C_TensorUnaryOperation(*this, escript::ES_optype::RECIP);    
}

Data
Data::wherePositive() const
{
    if (isComplex())
    {
        throw DataException("The wherePositive operation is not supported for complex data.");
    }
    MAKELAZYOP(GZ);
    return C_TensorUnaryOperation(*this, escript::ES_optype::GZ);    
}

Data
Data::whereNegative() const
{
    if (isComplex())
    {
        throw DataException("The whereNegative operation is not supported for complex data.");
    }
    MAKELAZYOP(LZ);
    return C_TensorUnaryOperation(*this, escript::ES_optype::LZ);    
}

Data
Data::whereNonNegative() const
{
    if (isComplex())
    {
        throw DataException("The whereNonNegative operation is not supported for complex data.");
    }
    MAKELAZYOP(GEZ);
    return C_TensorUnaryOperation(*this, escript::ES_optype::GEZ);    
}

Data
Data::whereNonPositive() const
{
    if (isComplex())
    {
        throw DataException("The whereNonPositive operation is not supported for complex data.");
    }
    MAKELAZYOP(LEZ);
    return C_TensorUnaryOperation(*this, escript::ES_optype::LEZ);
}

Data
Data::whereZero(real_t tol) const
{
    MAKELAZYOPOFF(EZ,tol);
    return C_TensorUnaryOperation(*this, escript::ES_optype::EZ, tol);
}

Data
Data::whereNonZero(real_t tol) const
{
    MAKELAZYOPOFF(NEZ,tol);
    return C_TensorUnaryOperation(*this, escript::ES_optype::NEZ, tol);
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
    if (functionspace.getDomain()!=getDomain())
        throw DataException("Error - gradient cannot be calculated on different domains.");
    DataTypes::ShapeType grad_shape=getDataPointShape();
    grad_shape.push_back(functionspace.getDim());
    Data out(0.0,grad_shape,functionspace,true);
    if (isComplex())
        out.complicate();    
    if (isLazy() && isComplex())
    {
        Data temp(*this);
        temp.resolve();     // since lazy complex is not a thing       
        getDomain()->setToGradient(out, temp);
    } 
    else
    {
        getDomain()->setToGradient(out,*this);
    }
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


DataTypes::RealVectorType::size_type
Data::getLength() const
{
    return m_data->getLength();
}


// There is no parallelism here ... elements need to be added in the correct order.
//   If we could pre-size the list and then fill in the elements it might work
//   This would need setting elements to be threadsafe.
//   Having multiple C threads calling into one interpreter is apparently a no-no.
const bp::object
Data::toListOfTuples(bool scalarastuple)
{
    if (get_MPISize()>1)
    {
        throw DataException("::toListOfTuples is not available for MPI with more than one process.");
    }
    unsigned int rank=getDataPointRank();
    unsigned int size=getDataPointSize();

    int npoints=getNumDataPoints();
    expand();                   // This will also resolve if required
    bp::list temp;
    temp.append(bp::object());
    bp::list res(temp*npoints);// pre-size the list by the "[None] * npoints"  trick
    if (isComplex())
    {
        const DataTypes::CplxVectorType& vec=getReady()->getTypedVectorRO(cplx_t(0));
        if (rank==0)
        {
            long count;
            if (scalarastuple)
            {
                for (count=0;count<npoints;++count)
                {
                    res[count]=bp::make_tuple(vec[count]);
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
    }
    else
    {
        const DataTypes::RealVectorType& vec=getReady()->getVectorRO();
        if (rank==0)
        {
            long count;
            if (scalarastuple)
            {
                for (count=0;count<npoints;++count)
                {
                    res[count]=bp::make_tuple(vec[count]);
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
    }
    return std::move(res);
}

const bp::object
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
        if (isComplex())
	{
	    DataTypes::CplxVectorType::size_type offset=getDataOffset(sampleNo, dataPointNoInSample);
	    return pointToTuple(getDataPointShape(),&(getDataAtOffsetRO(offset, cplx_t(0))));
	}
	else
	{
	    DataTypes::RealVectorType::size_type offset=getDataOffset(sampleNo, dataPointNoInSample);
	    return pointToTuple(getDataPointShape(),&(getDataAtOffsetRO(offset, real_t(0))));
	}
    }
    else
    {
        // The pre-numpy method would return an empty array of the given shape
        // I'm going to throw an exception because if we have zero points per sample we have problems
        throw DataException("Error - need at least 1 datapoint per sample.");
    }
}


void
Data::setValueOfDataPointToPyObject(int dataPointNo, const bp::object& py_object)
{
    // this will throw if the value cannot be represented
    setValueOfDataPointToArray(dataPointNo,py_object);
}


void
Data::setTupleForGlobalDataPoint(int id, int proc, bp::object v)
{
    THROWONCOMPLEX
#ifdef ESYS_MPI 
    int error=0;
#endif
    if( get_MPIRank()==proc )
    {
        try
        {
            bp::extract<real_t> dex(v);
            if (dex.check())
            {
                setValueOfDataPoint(id, dex());
            }
            else
            {
                setValueOfDataPointToArray(id, v);          
            }
        } 
        catch (...)
        {
#ifdef ESYS_MPI 
            error=1;
            int e2;
            MPI_Allreduce( &error, &e2, 1, MPI_INT, MPI_SUM, get_MPIComm() );      
#endif      
            // participate in gather
            throw;
        }
    }
#ifdef ESYS_MPI
    int e2;
    // If we get here, then either we succeeded or it was on another rank
    MPI_Allreduce( &error, &e2, 1, MPI_INT, MPI_MAX, get_MPIComm() );       
    // participate in gather
    if (e2)
    {
        throw DataException("Error in another rank performing setTupleForGlobalDataPoint");      
    }
#endif    
}


void
Data::setValueOfDataPointToArray(int dataPointNo, const bp::object& obj)
{
    if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
    }

    WrappedArray w(obj);
    if (w.isComplex() && (static_cast<unsigned int>(w.getRank())==0))
    {
	cplx_t v=w.getEltC();
	setValueOfDataPointC(dataPointNo, v);
	return;
    }
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
Data::setValueOfDataPoint(int dataPointNo, const real_t value)
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

void
Data::setValueOfDataPointC(int dataPointNo, const cplx_t value)
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
bp::object
Data::getValueOfGlobalDataPointAsTuple(int procNo, int dataPointNo)
{

    // copy datapoint into a buffer
    // broadcast buffer to all nodes
    // convert buffer to tuple
    // return tuple

    bp::tuple t;

    // This could be lazier than it is now
    forceResolve();

    int ndpps = getNumDataPointsPerSample();
    int sampleNo = dataPointNo/ndpps;
    int dataPointNoInSample = dataPointNo - sampleNo * ndpps;

    const DataTypes::ShapeType& dataPointShape = getDataPointShape();
    size_t length=DataTypes::noValues(dataPointShape);

    // TODO: global error handling
    if (get_MPIRank() == procNo && ndpps > 0 && ((sampleNo >= getNumSamples()) || (sampleNo < 0 ))) {
        throw DataException("Error - Data::getValueOfGlobalDataPointAsTuple: invalid sampleNo.");
    }

    if (get_MPIRank() == procNo && ndpps > 0 && ((dataPointNoInSample >= ndpps) || (dataPointNoInSample < 0))) {
        throw DataException("Error - Data::getValueOfGlobalDataPointAsTuple: invalid dataPointNoInSample.");
    }

    if(isComplex()){

        cplx_t *tmpData = new cplx_t[length];

        if( get_MPIRank() == procNo && ndpps > 0){
            DataTypes::CplxVectorType::size_type offset = getDataOffset(sampleNo, dataPointNoInSample);
            memcpy(tmpData,&(getDataAtOffsetRO(offset, static_cast<DataTypes::cplx_t>(0))), length*sizeof(cplx_t));
        }

        #ifdef ESYS_MPI
            // broadcast the data to all other processes
            MPI_Bcast(tmpData, length, MPI_DOUBLE_COMPLEX, procNo, get_MPIComm());
        #endif

        t = pointToTuple(dataPointShape,tmpData);
        delete [] tmpData;

    } else {

        real_t *tmpData = new real_t[length];

        if( get_MPIRank() == procNo && ndpps > 0){
            DataTypes::RealVectorType::size_type offset = getDataOffset(sampleNo, dataPointNoInSample);
            memcpy(tmpData,&(getDataAtOffsetRO(offset, static_cast<DataTypes::real_t>(0))), length*sizeof(real_t));
        }

        #ifdef ESYS_MPI
            // broadcast the data to all other processes
            MPI_Bcast(tmpData, length, MPI_DOUBLE, procNo, get_MPIComm());
        #endif

        t = pointToTuple(dataPointShape,tmpData);
        delete [] tmpData;

    }
    
    // return the loaded array
    return std::move(t);
}


bp::object
Data::integrateToTuple_const() const
{
    if (isLazy())
    {
        throw DataException("Error - cannot integrate for constant lazy data.");
    }
    if (isComplex()) {
        return integrateWorker<cplx_t>();
    } else {
        return integrateWorker<real_t>();
    }
}

bp::object
Data::integrateToTuple()
{
    if (isLazy())
    {
        expand();       // Can't do a non-resolving version of this without changing the domain object
    }                     // see the dom->setToIntegrals call. Not saying it can't be done, just not doing it yet.
    if (isComplex()) {
        return integrateWorker<cplx_t>();
    } else {
        return integrateWorker<real_t>();
    }
}

template<typename Scalar>
bp::object
Data::integrateWorker() const
{
    DataTypes::ShapeType shape = getDataPointShape();
    int dataPointSize = getDataPointSize();

    //
    // calculate the integral values
    vector<Scalar> integrals(dataPointSize);
    vector<Scalar> integrals_local(dataPointSize);
    const AbstractContinuousDomain* dom=dynamic_cast<const AbstractContinuousDomain*>(getDomain().get());
    if (dom == 0)
    {                             
        throw DataException("Can not integrate over non-continuous domains.");
    }
   
#ifdef ESYS_MPI
    if (isLazy() && isComplex())
    {
        Data temp(*this);
        temp.resolve();
        dom->setToIntegrals(integrals_local, temp);
    }
    else
    {
        dom->setToIntegrals(integrals_local, *this);
    }
    
    // Global sum: use an array instead of a vector because elements of array are guaranteed to be contiguous in memory
    Scalar* tmp = new Scalar[dataPointSize];
    Scalar* tmp_local = new Scalar[dataPointSize];
    for (int i=0; i<dataPointSize; i++) { tmp_local[i] = integrals_local[i]; }
    if (sizeof(Scalar) == sizeof(double)) {
        MPI_Allreduce(&tmp_local[0], &tmp[0], dataPointSize, MPI_DOUBLE, MPI_SUM, get_MPIComm());
    } else {
        MPI_Allreduce(&tmp_local[0], &tmp[0], dataPointSize, MPI_DOUBLE_COMPLEX, MPI_SUM, get_MPIComm());
    }
    for (int i=0; i<dataPointSize; i++) { integrals[i] = tmp[i]; }
    bp::tuple result = pointToTuple(shape, tmp);
    delete[] tmp;
    delete[] tmp_local;
#else
    if (isLazy() && isComplex())
    {
        Data temp(*this);
        temp.resolve();
        dom->setToIntegrals(integrals, temp);
    }
    else
    {
        dom->setToIntegrals(integrals, *this);
    }
    bp::tuple result = pointToTuple(shape, integrals);
#endif

    return std::move(result);
}

Data
Data::besselFirstKind(int order)
{
    THROWONCOMPLEX
    return bessel(order,boost::math::cyl_bessel_j);
}

Data
Data::besselSecondKind(int order)
{
    THROWONCOMPLEX
    return bessel(order,boost::math::cyl_neumann);
}

Data
Data::bessel(int order, real_t (*besselfunc) (int,real_t) )
{
    THROWONCOMPLEX
    DataTypes::real_t wantreal=0;
    if (isEmpty())  // do this before we attempt to interpolate
    {
     throw DataException("Error - Operations (bessel) not permitted on instances of DataEmpty.");
    }
    if (isLazy())
    {
        resolve();
    }
    // Interpolate if necessary and find an appropriate function space
    Data arg_0_Z = Data(*this);

    // Get rank and shape of inputs
    const DataTypes::ShapeType& shape0 = arg_0_Z.getDataPointShape();
    int size0 = arg_0_Z.getDataPointSize();

    // Declare output Data object
    Data res;

    if (arg_0_Z.isConstant()) {
        res = Data(0.0, shape0, arg_0_Z.getFunctionSpace(),false);      // DataConstant output
        const real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(0, wantreal));
        real_t *ptr_2 = &(res.getDataAtOffsetRW(0, wantreal));
        for (int i = 0; i < size0; ++i) {
            ptr_2[i] = besselfunc(order, ptr_0[i]);
        }
    }
    else if (arg_0_Z.isTagged()) {

        // Borrow DataTagged input from Data object
        DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());

        // Prepare a DataTagged output 2
        res = Data(0.0, shape0, arg_0_Z.getFunctionSpace(), false);   // DataTagged output
        res.tag();
        DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());

        // Get the pointers to the actual data
        const real_t *ptr_0 = &(tmp_0->getDefaultValueRO(0, wantreal));
        real_t *ptr_2 = &(tmp_2->getDefaultValueRW(0, wantreal));
        // Compute a result for the default
        for (int i = 0; i < size0; ++i) {
            ptr_2[i] = besselfunc(order, ptr_0[i]);
        }
        // Compute a result for each tag
        const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
        DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
        for (i=lookup_0.begin();i!=lookup_0.end();i++) {
            tmp_2->addTag(i->first);
            const real_t *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0));
            real_t *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0));
            for (int i = 0; i < size0; ++i) {
                ptr_2[i] = besselfunc(order, ptr_0[i]);
            }

        }

    }
    else if (arg_0_Z.isExpanded()) {

        res = Data(0.0, shape0, arg_0_Z.getFunctionSpace(),true); // DataExpanded output
        DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
        DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

        int sampleNo_0,dataPointNo_0;
        int numSamples_0 = arg_0_Z.getNumSamples();
        int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
	DataTypes::real_t wantreal=0;
        #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
        for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
            dataPointNo_0=0;
        //      for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
            int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
            int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
            const real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, wantreal));
            real_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, wantreal));
            for (int i = 0; i < size0*numDataPointsPerSample_0; ++i) {
                ptr_2[i] = besselfunc(order, ptr_0[i]);
            }
        //      }
        }
  }
  else {
    throw DataException("Error - Bessel function: unknown combination of inputs");
  }

  return res;
}

Data
Data::conjugate() const
{
    if (isLazy())
    {
        Data temp(*this);
        temp.resolve();
        return temp.conjugate();
    }
    if (isComplex())
    {
        return C_TensorUnaryOperation(*this, escript::ES_optype::CONJ);      
    }
    else
    {
        return copySelf();
    }
}


Data
Data::real() const
{
    if (isLazy())
    {
        Data temp(*this);
        temp.resolve();
        return temp.real();
    }
    if (isComplex())
    {
        return C_TensorUnaryOperation(*this, escript::ES_optype::REAL);      
    }
    else
    {
        return copySelf();
    }
}

Data
Data::imag() const
{
    if (isLazy())
    {
        Data temp(*this);
        temp.resolve();
        return temp.imag();
    }
    if (isComplex())
    {
        return C_TensorUnaryOperation(*this, escript::ES_optype::IMAG);      
    }
    else
    {
        return copySelf()*Data(0, m_data->getShape(), getFunctionSpace(),false);      // return an object with same tags etc but all values 0
                                // This is not efficient, but why are you taking imag of R anyway?
    }
}

Data
Data::phase() const
{
    if (isLazy())
    {
        Data temp(*this);
        temp.resolve();
        return temp.phase();
    }
    if (isComplex())
    {
        return C_TensorUnaryOperation(*this, escript::ES_optype::PHS);      
    }
    else
    {
        return whereNegative()*Data(M_PI, DataTypes::scalarShape, getFunctionSpace(), false);
    }
}





Data
Data::sin() const
{
    MAKELAZYOP(SIN);
    return C_TensorUnaryOperation(*this, escript::ES_optype::SIN);
}

Data
Data::cos() const
{
    MAKELAZYOP(COS);
    return C_TensorUnaryOperation(*this, escript::ES_optype::COS);
}

Data
Data::tan() const
{
    MAKELAZYOP(TAN);
    return C_TensorUnaryOperation(*this, escript::ES_optype::TAN);
}

Data
Data::asin() const
{
    MAKELAZYOP(ASIN);
    return C_TensorUnaryOperation(*this, escript::ES_optype::ASIN);
}

Data
Data::acos() const
{
    MAKELAZYOP(ACOS);
    return C_TensorUnaryOperation(*this, escript::ES_optype::ACOS);
}


Data
Data::atan() const
{
    MAKELAZYOP(ATAN);
    return C_TensorUnaryOperation(*this, escript::ES_optype::ATAN);
}

Data
Data::sinh() const
{
    MAKELAZYOP(SINH);
    return C_TensorUnaryOperation(*this, escript::ES_optype::SINH);
}

Data
Data::cosh() const
{
    MAKELAZYOP(COSH);
    return C_TensorUnaryOperation(*this, escript::ES_optype::COSH);
}

Data
Data::tanh() const
{
    MAKELAZYOP(TANH);
    return C_TensorUnaryOperation(*this, escript::ES_optype::TANH);
}


Data
Data::erf() const
{
    MAKELAZYOP(ERF);
    return C_TensorUnaryOperation(*this, escript::ES_optype::ERF);
}

Data
Data::asinh() const
{
    MAKELAZYOP(ASINH);
    return C_TensorUnaryOperation(*this, escript::ES_optype::ASINH);
}

Data
Data::acosh() const
{
    MAKELAZYOP(ACOSH);
    return C_TensorUnaryOperation(*this, escript::ES_optype::ACOSH);    
}

Data
Data::atanh() const
{
    MAKELAZYOP(ATANH);
    return C_TensorUnaryOperation(*this, escript::ES_optype::ATANH); 
}

Data
Data::log10() const
{
    MAKELAZYOP(LOG10);
    return C_TensorUnaryOperation(*this, escript::ES_optype::LOG10);
}

Data
Data::log() const
{
    MAKELAZYOP(LOG);
    return C_TensorUnaryOperation(*this, escript::ES_optype::LOG);
}

Data
Data::sign() const
{
    THROWONCOMPLEX
    MAKELAZYOP(SIGN);
    return C_TensorUnaryOperation(*this, escript::ES_optype::SIGN);
}

Data
Data::abs() const
{
    MAKELAZYOP(ABS);
    return C_TensorUnaryOperation(*this, escript::ES_optype::ABS);
}

Data
Data::neg() const
{
    MAKELAZYOP(NEG);
    return C_TensorUnaryOperation(*this, escript::ES_optype::NEG);
}

Data
Data::pos() const
{
    THROWONCOMPLEX
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
    MAKELAZYOP(EXP);
    return C_TensorUnaryOperation(*this, escript::ES_optype::EXP);
}

Data
Data::sqrt() const
{
    MAKELAZYOP(SQRT);
    return C_TensorUnaryOperation(*this, escript::ES_optype::SQRT);
}

real_t
Data::Lsup_const() const
{
    if (isLazy())
    {
        throw DataException("Error - cannot compute Lsup for constant lazy data.");
    }
    return LsupWorker();
}

real_t
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
#ifdef ESYS_MPI
            if (isComplex())
            {
                return lazyAlgWorker<AbsMax<cplx_t> >(0,MPI_MAX);
            }
            else
            {
                return lazyAlgWorker<AbsMax<real_t> >(0,MPI_MAX);             
            }
#else
            if (isComplex())
            {
                return lazyAlgWorker<AbsMax<cplx_t> >(0);
            }
            else
            {
                return lazyAlgWorker<AbsMax<real_t> >(0);
            }
#endif
        }
    }
    return LsupWorker();
}

real_t
Data::sup_const() const
{
    if (isComplex())
    {
        throw DataException("Error Cannot compute sup() for complex data.");
    }  
    if (isLazy())
    {
        throw DataException("Error - cannot compute sup for constant lazy data.");
    }
    return supWorker();
}

real_t
Data::sup() 
{
    if (isComplex())
    {
        throw DataException("Error Cannot compute sup() for complex data.");
    }
    if (isLazy())
    {
        if (!actsExpanded() || CHECK_DO_CRES)
        {
            resolve();
        }
        else
        {
#ifdef ESYS_MPI
            return lazyAlgWorker<FMax>(numeric_limits<real_t>::max()*-1, MPI_MAX);
#else
            return lazyAlgWorker<FMax>(numeric_limits<real_t>::max()*-1);
#endif
        }
    }
    return supWorker();
}

real_t
Data::inf_const() const
{
    if (isComplex())
    {
        throw DataException("Error Cannot compute inf() for complex data.");
    }  
    if (isLazy())
    {
        throw DataException("Error - cannot compute inf for constant lazy data.");
    }
    return infWorker();
}

real_t
Data::inf() 
{
    if (isComplex())
    {
        throw DataException("Error Cannot compute inf() for complex data.");
    }  
    if (isLazy())
    {
        if (!actsExpanded() || CHECK_DO_CRES)
        {
            resolve();
        }
        else
        {
#ifdef ESYS_MPI
            return lazyAlgWorker<FMin>(numeric_limits<real_t>::max(), MPI_MIN);
#else
            return lazyAlgWorker<FMin>(numeric_limits<real_t>::max());
#endif
        }
    }
    return infWorker();
}

template <class BinaryOp>
real_t
#ifdef ESYS_MPI
Data::lazyAlgWorker(real_t init, MPI_Op mpiop_type)
#else
Data::lazyAlgWorker(real_t init)
#endif
{
    if (!isLazy() || !m_data->actsExpanded())
    {
        throw DataException("Error - lazyAlgWorker can only be called on lazy(expanded) data.");
    }
    DataLazy* dl=dynamic_cast<DataLazy*>(m_data.get());
    ESYS_ASSERT(dl!=NULL, "Programming error - casting to DataLazy.");
    real_t val=init;
    int i=0;
    const size_t numsamples=getNumSamples();
    const size_t samplesize=getNoValues()*getNumDataPointsPerSample();
    BinaryOp operation;
    real_t localValue=0, globalValue;    
    #pragma omp parallel private(i)
    {
        real_t localtot=init;
        #pragma omp for schedule(static) private(i)
        for (i=0;i<numsamples;++i)
        {
            size_t roffset=0;
            auto v=dl->resolveTypedSample(i, roffset, typename BinaryOp::second_argument_type(0));
            // Now we have the sample, run operation on all points
            for (size_t j=0;j<samplesize;++j)
            {
                localtot=operation(localtot,(*v)[j+roffset]);
            }
            if (escript::vectorHasNaN(*v,roffset, samplesize))
            {
                #pragma omp critical
                {
                    localValue=1.0;
                }
            }
        }
        #pragma omp critical
        val=operation(val,localtot);
    }
#ifdef ESYS_MPI
    MPI_Allreduce( &localValue, &globalValue, 1, MPI_DOUBLE, MPI_MAX, get_MPIComm() );
#else
    globalValue=localValue;
#endif
    if (globalValue!=0)
    {
        return makeNaN();
    }
#ifdef ESYS_MPI
    MPI_Allreduce( &val, &globalValue, 1, MPI_DOUBLE, mpiop_type, get_MPIComm() );
    return globalValue;
#else
    return val;
#endif
}

bool
Data::hasNaN()
{
    if (isLazy())
    {
        resolve();
    }
    bool haveNaN=getReady()->hasNaN(); 
    return haveNaN;
}

void
Data::replaceNaN(real_t value)
{
    if (isLazy())
    {
        resolve();
    }
    getReady()->replaceNaN(value); 
}

void
Data::replaceNaN(cplx_t value)
{
    if (isLazy())
    {
        resolve();
    }
    getReady()->replaceNaN(value); 
}

void
Data::replaceNaNPython(boost::python::object obj)
{
    boost::python::extract<DataTypes::real_t> exr(obj);
    if (exr.check())
    {
	replaceNaN(exr());
    }
    else
    {
	replaceNaN(boost::python::extract<DataTypes::cplx_t>(obj)());
    }
}

bool
Data::hasInf()
{
    if (isLazy())
    {
        resolve();
    }
    bool haveInf=getReady()->hasInf(); 
    return haveInf;
}

void
Data::replaceInf(real_t value)
{
    if (isLazy())
    {
        resolve();
    }
    getReady()->replaceInf(value); 
}

void
Data::replaceInf(cplx_t value)
{
    if (isLazy())
    {
        resolve();
    }
    getReady()->replaceInf(value); 
}

void
Data::replaceInfPython(boost::python::object obj)
{
    boost::python::extract<DataTypes::real_t> exr(obj);
    if (exr.check())
    {
	replaceInf(exr());
    }
    else
    {
	replaceInf(boost::python::extract<DataTypes::cplx_t>(obj)());
    }
}

// Do not call this on Lazy Data use the proper entry point
real_t
Data::LsupWorker() const
{
    bool haveNaN=getReady()->hasNaN();

  
#ifdef ESYS_MPI
    int nanchecker=0;    
    if (haveNaN)
    {
        nanchecker=1.0;
    }
    int globalnan;
    MPI_Allreduce( &nanchecker, &globalnan, 1, MPI_INT, MPI_MAX, get_MPIComm() );
    if (globalnan!=0)
    {
        return makeNaN();
    }
#else
    if (haveNaN) 
    {
        return makeNaN();
    }
#endif

    //
    // set the initial absolute maximum value to zero
    if (isComplex())
    {
        AbsMax<cplx_t> abs_max_func;
        real_t localValue=0;
        localValue = reduction(abs_max_func,0);

    #ifdef ESYS_MPI
        real_t globalValue=0;
        MPI_Allreduce( &localValue, &globalValue, 1, MPI_DOUBLE, MPI_MAX, get_MPIComm() );
        return globalValue;
    #else
        return localValue;
    #endif
    }
    else
    {  
        AbsMax<real_t> abs_max_func;
        real_t localValue=0;
        localValue = reduction(abs_max_func,0);

    #ifdef ESYS_MPI
        real_t globalValue=0;   
        MPI_Allreduce( &localValue, &globalValue, 1, MPI_DOUBLE, MPI_MAX, get_MPIComm() );
        return globalValue;
    #else
        return localValue;
    #endif
    }
}

real_t
Data::supWorker() const
{
    bool haveNaN=getReady()->hasNaN();
    real_t localValue=0;

#ifdef ESYS_MPI
    if (haveNaN)
    {
        localValue=1.0;
    }
    real_t globalValue;
    MPI_Allreduce( &localValue, &globalValue, 1, MPI_DOUBLE, MPI_MAX, get_MPIComm() );
    if (globalValue!=0)
    {
        return makeNaN();
    }
#else
    if (haveNaN) 
    {
        return makeNaN();
    }
#endif

    //
    // set the initial maximum value to min possible real_t
    FMax fmax_func;
    if (hasNoSamples())
    {
	localValue = numeric_limits<real_t>::infinity()*-1;
    }
    else
    {
	localValue = reduction(fmax_func,numeric_limits<real_t>::infinity()*-1);      
    }
    #ifdef ESYS_MPI
    MPI_Allreduce( &localValue, &globalValue, 1, MPI_DOUBLE, MPI_MAX, get_MPIComm() );
    return globalValue;
#else
    return localValue;
#endif
}

real_t
Data::infWorker() const
{
    bool haveNaN=getReady()->hasNaN();
    real_t localValue=0;

#ifdef ESYS_MPI
    if (haveNaN)
    {
        localValue=1.0;
    }
    real_t globalValue;
    MPI_Allreduce( &localValue, &globalValue, 1, MPI_DOUBLE, MPI_MAX, get_MPIComm() );
    if (globalValue!=0)
    {
        return makeNaN();
    }
#else
    if (haveNaN) 
    {
        return makeNaN();
    }
#endif
    //
    // set the initial minimum value to max possible real_t
    FMin fmin_func;
    if (hasNoSamples())
    {
	localValue = numeric_limits<real_t>::infinity();
    }
    else
    {
	localValue = reduction(fmin_func,numeric_limits<real_t>::infinity());
    }
#ifdef ESYS_MPI
    MPI_Allreduce( &localValue, &globalValue, 1, MPI_DOUBLE, MPI_MIN, get_MPIComm() );
    return globalValue;
#else
    return localValue;
#endif
}

/* global reduction */


inline Data
Data::minval_nonlazy() const
{
    THROWONCOMPLEX
    //
    // set the initial minimum value to max possible real_t
    FMin fmin_func;
    return dp_algorithm(fmin_func,numeric_limits<real_t>::max());
}


inline Data
Data::maxval_nonlazy() const
{
    THROWONCOMPLEX
    //
    // set the initial maximum value to min possible real_t
    FMax fmax_func;
    return dp_algorithm(fmax_func,numeric_limits<real_t>::max()*-1);
}


Data
Data::maxval() const
{
    THROWONCOMPLEX
    MAKELAZYOP(MAXVAL);
    return maxval_nonlazy();
}


Data
Data::minval() const
{
    THROWONCOMPLEX
    MAKELAZYOP(MINVAL);
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
    MAKELAZYOP2(SWAP,axis0,axis1);
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
    Data ev(0.,ev_shape,getFunctionSpace(), false);
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
    MAKELAZYOP(SYM);
    Data ev(0.,getDataPointShape(),getFunctionSpace(), false);
    ev.typeMatchRight(*this);
    m_data->symmetric(ev.m_data.get());
    return ev;
}

Data
Data::antisymmetric() const
{
    // check input
    DataTypes::ShapeType s=getDataPointShape();
    if (getDataPointRank()==2) {
        if(s[0] != s[1])
            throw DataException("Error - Data::antisymmetric can only be calculated for rank 2 object with equal first and second dimension.");
	MAKELAZYOP(NSYM);
        DataTypes::ShapeType ev_shape;
        ev_shape.push_back(s[0]);
        ev_shape.push_back(s[1]);
        Data ev(0.,ev_shape,getFunctionSpace(),false);
        ev.typeMatchRight(*this);
        m_data->antisymmetric(ev.m_data.get());
        return ev;
    }
    else if (getDataPointRank()==4) {
        if(!(s[0] == s[2] && s[1] == s[3]))
            throw DataException("Error - Data::antisymmetric can only be calculated for rank 4 object with dim0==dim2 and dim1==dim3.");
	MAKELAZYOP(NSYM);
        DataTypes::ShapeType ev_shape;
        ev_shape.push_back(s[0]);
        ev_shape.push_back(s[1]);
        ev_shape.push_back(s[2]);
        ev_shape.push_back(s[3]);
        Data ev(0.,ev_shape,getFunctionSpace(),false);
        ev.typeMatchRight(*this);
        m_data->antisymmetric(ev.m_data.get());
        return ev;
    }
    else {
        throw DataException("Error - Data::antisymmetric can only be calculated for rank 2 or 4 object.");
    }
}

Data
Data::hermitian() const
{
    if (!isComplex())
    {
        return symmetric();
    }
    // check input
    DataTypes::ShapeType s=getDataPointShape();
    if (getDataPointRank()==2) {
        if(s[0] != s[1])
            throw DataException("Error - Data::hermitian can only be calculated for rank 2 object with equal first and second dimension.");
    }
    else if (getDataPointRank()==4) {
        if(!(s[0] == s[2] && s[1] == s[3]))
            throw DataException("Error - Data::hermitian can only be calculated for rank 4 object with dim0==dim2 and dim1==dim3.");
    }
    else {
        throw DataException("Error - Data::hermitian can only be calculated for rank 2 or 4 object.");
    }
    MAKELAZYOP(HER);
    Data ev(0.,getDataPointShape(),getFunctionSpace(), false);
    ev.typeMatchRight(*this);
    m_data->hermitian(ev.m_data.get());
    return ev;
}

Data
Data::antihermitian() const
{
    if (!isComplex())
    {
        return antisymmetric();
    }  
    // check input
    DataTypes::ShapeType s=getDataPointShape();
    if (getDataPointRank()==2) {
        if(s[0] != s[1])
            throw DataException("Error - Data::antihermitian can only be calculated for rank 2 object with equal first and second dimension.");
	MAKELAZYOP(NHER);
        DataTypes::ShapeType ev_shape;
        ev_shape.push_back(s[0]);
        ev_shape.push_back(s[1]);
        Data ev(0.,ev_shape,getFunctionSpace(), false);
        ev.typeMatchRight(*this);
        m_data->antihermitian(ev.m_data.get());
        return ev;
    }
    else if (getDataPointRank()==4) {
        if(!(s[0] == s[2] && s[1] == s[3]))
            throw DataException("Error - Data::antihermitian can only be calculated for rank 4 object with dim0==dim2 and dim1==dim3.");
	MAKELAZYOP(NHER);
        DataTypes::ShapeType ev_shape;
        ev_shape.push_back(s[0]);
        ev_shape.push_back(s[1]);
        ev_shape.push_back(s[2]);
        ev_shape.push_back(s[3]);
        Data ev(0.,ev_shape,getFunctionSpace(), false);
        ev.typeMatchRight(*this);
        m_data->antihermitian(ev.m_data.get());
        return ev;
    }
    else {
        throw DataException("Error - Data::antihermitian can only be calculated for rank 2 or 4 object.");
    }
}

Data
Data::trace(int axis_offset) const
{
    MAKELAZYOPOFF(TRACE,axis_offset);
    if ((axis_offset<0) || (axis_offset>getDataPointRank()))
    {
        throw DataException("Error - Data::trace, axis_offset must be between 0 and rank-2 inclusive.");
    }
    DataTypes::ShapeType s=getDataPointShape();
    if (getDataPointRank()==2) {
        DataTypes::ShapeType ev_shape;
        Data ev(0.,ev_shape,getFunctionSpace(),false);
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
        Data ev(0.,ev_shape,getFunctionSpace(),false);
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
        Data ev(0.,ev_shape,getFunctionSpace(),false);
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
    MAKELAZYOPOFF(TRANS,axis_offset);
    DataTypes::ShapeType s=getDataPointShape();
    DataTypes::ShapeType ev_shape;
    // Here's the equivalent of python s_out=s[axis_offset:]+s[:axis_offset]
    // which goes thru all shape vector elements starting with axis_offset (at index=rank wrap around to 0)
    int rank=getDataPointRank();
    if (axis_offset<0 || axis_offset>rank) {
        stringstream e;
        e << "Error - Data::transpose must have 0 <= axis_offset <= rank=" << rank;
        throw DataException(e.str());
    }
    for (int i=0; i<rank; i++) {
        int index = (axis_offset+i)%rank;
        ev_shape.push_back(s[index]); // Append to new shape
    }
    Data ev(0.,ev_shape,getFunctionSpace(), false);
    ev.typeMatchRight(*this);
    m_data->transpose(ev.m_data.get(), axis_offset);
    return ev;
}

Data
Data::eigenvalues() const
{
    if (isLazy())
    {
        Data temp(*this);       // to get around the fact that you can't resolve a const Data
        temp.resolve();
        return temp.eigenvalues();
    }
    // check input
    DataTypes::ShapeType s=getDataPointShape();
    if (getDataPointRank()!=2)
        throw DataException("Error - Data::eigenvalues can only be calculated for rank 2 object.");
    if(s[0] != s[1])
        throw DataException("Error - Data::eigenvalues can only be calculated for object with equal first and second dimension.");
    if (isComplex() && (s[0]>2))
    {
        throw DataException("Error - Data::eigenvalues not supported for complex 3x3.");
    }
    // create return
    DataTypes::ShapeType ev_shape(1,s[0]);
    Data ev(0.,ev_shape,getFunctionSpace(),false);
    ev.typeMatchRight(*this);
    m_data->eigenvalues(ev.m_data.get());
    return ev;
}

const bp::tuple
Data::eigenvalues_and_eigenvectors(const real_t tol) const
{
    THROWONCOMPLEX
    if (isLazy())
    {
        Data temp(*this);       // to get around the fact that you can't resolve a const Data
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
    Data ev(0.,ev_shape,getFunctionSpace(), false);
    ev.typeMatchRight(*this);
    DataTypes::ShapeType V_shape(2,s[0]);
    Data V(0.,V_shape,getFunctionSpace(), false);
    V.typeMatchRight(*this);
    m_data->eigenvalues_and_eigenvectors(ev.m_data.get(),V.m_data.get(),tol);
    return bp::make_tuple(bp::object(ev),bp::object(V));
}

const bp::tuple
Data::minGlobalDataPoint() const
{
    // NB: calc_minGlobalDataPoint( had to be split off from minGlobalDataPoint( as boost::make_tuple causes an
    // abort (for unknown reasons) if there are openmp directives with it in the
    // surrounding function

    THROWONCOMPLEX
    int DataPointNo;
    int ProcNo;
    calc_minGlobalDataPoint(ProcNo,DataPointNo);
    if (ProcNo==-1) {
	throw DataException("There are no values to find minimum of.");
    }
    return bp::make_tuple(ProcNo,DataPointNo);
}

void
Data::calc_minGlobalDataPoint(int& ProcNo,
                        int& DataPointNo) const
{
    THROWONCOMPLEX
    if (isLazy())
    {
        Data temp(*this);   // to get around the fact that you can't resolve a const Data
        temp.resolve();
        return temp.calc_minGlobalDataPoint(ProcNo,DataPointNo);
    }
    int i,j;
    int lowi=0,lowj=0;
    real_t min=numeric_limits<real_t>::max();

    Data temp=minval_nonlazy();   // need to do this to prevent autolazy from reintroducing laziness

    int numSamples=temp.getNumSamples();
    int numDPPSample=temp.getNumDataPointsPerSample();

    real_t local_val, local_min;
#ifdef ESYS_MPI
    real_t next[2];
#endif
    int local_lowi=0,local_lowj=0;        

#pragma omp parallel firstprivate(local_lowi,local_lowj) private(local_val,local_min)
    {
        local_min=min;
	DataTypes::real_t wantreal=0;
#pragma omp for private(i,j) schedule(static)
        for (i=0; i<numSamples; i++) {
            for (j=0; j<numDPPSample; j++) {
                local_val=temp.getDataAtOffsetRO(temp.getDataOffset(i,j), wantreal);
                if (local_val<local_min) {
                    local_min=local_val;
                    local_lowi=i;
                    local_lowj=j;
                }
            }
        }
#pragma omp critical
        if (local_min<min) { // If we found a smaller value than our sentinel
            min=local_min;
            lowi=local_lowi;
            lowj=local_lowj;
        }
    } // parallel section

#ifdef ESYS_MPI
    // determine the processor on which the minimum occurs
    next[0] = min;
    next[1] = numSamples;
    int lowProc = 0;
    real_t *globalMins = new real_t[get_MPISize()*2+1];
    /*int error =*/ MPI_Gather (next, 2, MPI_DOUBLE, globalMins, 2, MPI_DOUBLE, 0, get_MPIComm() );

    if ( get_MPIRank()==0 ) {
        for (lowProc=0; lowProc<get_MPISize(); lowProc++) {
            if (globalMins[lowProc*2+1] > 0) break;
        }
	if (lowProc<get_MPISize()) {
            min = globalMins[lowProc*2];
            for( i=lowProc+1; i<get_MPISize(); i++ )
                if( globalMins[i*2+1]>0 && min>globalMins[i*2] ) {
                    lowProc = i;
                    min = globalMins[i*2];
                }
        }
    }
    MPI_Bcast( &lowProc, 1, MPI_INT, 0, get_MPIComm() );
    DataPointNo = lowj + lowi * numDPPSample;
    if (lowProc>=get_MPISize()) {
        ProcNo = -1;
    } else {
        MPI_Bcast(&DataPointNo, 1, MPI_INT, lowProc, get_MPIComm() );
        ProcNo = lowProc;
    }
    delete [] globalMins;
#else
    ProcNo = 0;
    DataPointNo = lowj + lowi * numDPPSample;
#endif
}


const bp::tuple
Data::maxGlobalDataPoint() const
{
    THROWONCOMPLEX
    int DataPointNo;
    int ProcNo;
    calc_maxGlobalDataPoint(ProcNo,DataPointNo);
    return bp::make_tuple(ProcNo,DataPointNo);
}

void
Data::calc_maxGlobalDataPoint(int& ProcNo,
                        int& DataPointNo) const
{
    if (isLazy())
    {
        Data temp(*this);   // to get around the fact that you can't resolve a const Data
        temp.resolve();
        return temp.calc_maxGlobalDataPoint(ProcNo,DataPointNo);
    }
    THROWONCOMPLEX
    int i,j;
    int highi=0,highj=0;
    //-------------
    real_t max= -numeric_limits<real_t>::max();

    Data temp=maxval_nonlazy();   // need to do this to prevent autolazy from reintroducing laziness

    int numSamples=temp.getNumSamples();
    int numDPPSample=temp.getNumDataPointsPerSample();

    real_t local_val, local_max;
#ifdef ESYS_MPI
    real_t next[2];
#endif
    int local_highi=0,local_highj=0;      

#pragma omp parallel firstprivate(local_highi,local_highj) private(local_val,local_max)
    {
        local_max=max;
	DataTypes::real_t wantreal=0;
#pragma omp for private(i,j) schedule(static)
        for (i=0; i<numSamples; i++) {
            for (j=0; j<numDPPSample; j++) {
                local_val=temp.getDataAtOffsetRO(temp.getDataOffset(i,j), wantreal);
                if (local_val>local_max) {
                    local_max=local_val;
                    local_highi=i;
                    local_highj=j;
                }
            }
        }
#pragma omp critical
        if (local_max>max) { // If we found a larger value than our sentinel
            max=local_max;
            highi=local_highi;
            highj=local_highj;
        }
    }
#ifdef ESYS_MPI
    // determine the processor on which the maximum occurs
    next[0] = max;
    next[1] = numSamples;
    int highProc = 0;
    real_t *globalMaxs = new real_t[get_MPISize()*2+1];
    /*int error =*/ MPI_Gather ( next, 2, MPI_DOUBLE, globalMaxs, 2, MPI_DOUBLE, 0, get_MPIComm() );
    if( get_MPIRank()==0 ){
        for (highProc=0; highProc<get_MPISize(); highProc++)
            if (globalMaxs[highProc*2+1] > 0) break;
        max = globalMaxs[highProc*2];
        for( i=highProc+1; i<get_MPISize(); i++ )
        {
            if( globalMaxs[i*2+1]>0 && max<globalMaxs[i*2] )
            {
                highProc = i;
                max = globalMaxs[i*2];
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

Data&
Data::operator+=(const Data& right)
{
    if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
    }
    MAKELAZYBINSELF(right,ADD);    // for lazy + is equivalent to +=
    exclusiveWrite();                     // Since Lazy data does not modify its leaves we only need to worry here
    if (!isComplex() && right.isComplex())
    {
        complicate();
    }
    TensorSelfUpdateBinaryOperation(right, escript::ES_optype::ADD);  
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
    m_protected=false; // since any changes should be caught by exclusiveWrite()
    set_m_data(other.m_data);
    return (*this);
}

Data&
Data::operator-=(const Data& right)
{
    if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
    }
    MAKELAZYBINSELF(right,SUB);
    exclusiveWrite();
    if (!isComplex() && right.isComplex())
    {
        complicate();
    }
    TensorSelfUpdateBinaryOperation(right, escript::ES_optype::SUB);
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
    MAKELAZYBINSELF(right,MUL);
    exclusiveWrite();
    if (!isComplex() && right.isComplex())
    {
        complicate();
    }
    TensorSelfUpdateBinaryOperation(right, escript::ES_optype::MUL);
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
    MAKELAZYBINSELF(right,DIV);
    exclusiveWrite();
    if (!isComplex() && right.isComplex())
    {
        complicate();
    }
    TensorSelfUpdateBinaryOperation(right, escript::ES_optype::DIV);
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

/* Be careful trying to make this operation lazy.
At time of writing, resolve() and resolveSample() do not throw.
Changing this would mean that any resolve call would need to use MPI (to check for global errors)
*/
Data
Data::matrixInverse() const
{
    if (isLazy())       // Cannot use lazy for this because individual inversions could throw.
    {
        Data d(*this);
        d.resolve();
        return d.matrixInverse();
    }
    THROWONCOMPLEX
    Data out(0.,getDataPointShape(),getFunctionSpace(),false);
    out.typeMatchRight(*this);

    DataReady* drp=out.getReadyPtr().get();
    int errcode=m_data->matrixInverse(drp);
#ifdef ESYS_MPI
    int globalval=0;
    MPI_Allreduce( &errcode, &globalval, 1, MPI_INT, MPI_MAX, get_MPIComm() );
    errcode=globalval;
#endif
    if (errcode)
    {
        escript::matrixInverseError(errcode); // throws exceptions
    }
    return out;
}


Data
Data::rpowO(const bp::object& left) const
{
    Data left_d(left,*this);
    return left_d.powD(*this);
}

Data
Data::powO(const bp::object& right) const
{
    Data tmp(right,getFunctionSpace(),false);
    return powD(tmp);
}

Data
Data::powD(const Data& right) const
{
    MAKELAZYBIN(right,POW);
    
    return C_TensorBinaryOperation(*this, right, ES_optype::POW);    
    
}


//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator+(const Data& left, const Data& right)
{
    MAKELAZYBIN2(left,right,ADD);
    
    return C_TensorBinaryOperation(left, right, ES_optype::ADD);
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator-(const Data& left, const Data& right)
{
    MAKELAZYBIN2(left,right,SUB);
    return C_TensorBinaryOperation(left, right, ES_optype::SUB);    
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator*(const Data& left, const Data& right)
{
    MAKELAZYBIN2(left,right,MUL);    
    
    return C_TensorBinaryOperation(left, right, ES_optype::MUL);        
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator/(const Data& left, const Data& right)
{
    MAKELAZYBIN2(left,right,DIV);
    return C_TensorBinaryOperation(left, right, ES_optype::DIV);        
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator+(const Data& left, const boost::python::object& right)
{
    Data tmp(right,left.getFunctionSpace(),false);
    MAKELAZYBIN2(left,tmp,ADD);
    return left+tmp;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator-(const Data& left, const boost::python::object& right)
{
    Data tmp(right,left.getFunctionSpace(),false);
    MAKELAZYBIN2(left,tmp,SUB);
    return left-tmp;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator*(const Data& left, const boost::python::object& right)
{
    Data tmp(right,left.getFunctionSpace(),false);
    MAKELAZYBIN2(left,tmp,MUL);
    return left*tmp;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator/(const Data& left, const boost::python::object& right)
{
    Data tmp(right,left.getFunctionSpace(),false);
    MAKELAZYBIN2(left,tmp,DIV);
    return left/tmp;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator+(const boost::python::object& left, const Data& right)
{
    Data tmp(left,right.getFunctionSpace(),false);
    MAKELAZYBIN2(tmp,right,ADD);
    return tmp+right;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator-(const boost::python::object& left, const Data& right)
{
    Data tmp(left,right.getFunctionSpace(),false);
    MAKELAZYBIN2(tmp,right,SUB);
    return tmp-right;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator*(const boost::python::object& left, const Data& right)
{
    Data tmp(left,right.getFunctionSpace(),false);
    MAKELAZYBIN2(tmp,right,MUL);
    return tmp*right;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator/(const boost::python::object& left, const Data& right)
{
    Data tmp(left,right.getFunctionSpace(),false);
    MAKELAZYBIN2(tmp,right,DIV);
    return tmp/right;
}


/* TODO */
/* global reduction */
Data
Data::getItem(const bp::object& key) const
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
Data::setItemO(const bp::object& key,
               const bp::object& value)
{
    Data tempData(value,getFunctionSpace(),false);
    setItemD(key,tempData);
}

void
Data::setItemD(const bp::object& key,
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
    exclusiveWrite(); // In case someone finds a way to call this without going through setItemD
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
    if (isComplex())
    {
        right.complicate();
    }
    if (isExpanded()) {
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
    if (right.isComplex())
    {
	complicate();
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
                           const bp::object& value)
{
    if (getFunctionSpace().getDomain()->isValidTagName(name)) {
        forceResolve();
        exclusiveWrite();
        int tagKey=getFunctionSpace().getDomain()->getTag(name);
        setTaggedValue(tagKey,value);
    }
    else
    {
        std::string msg="Error - unknown tag ("+name+") in setTaggedValueByName.";
        throw DataException(msg);
    }
}

void
Data::setTaggedValue(int tagKey,
                     const bp::object& value)
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

    if (w.isComplex())
    {
        CplxVectorType temp_data2;
        temp_data2.copyFromArray(w,1);

        m_data->setTaggedValue(tagKey,w.getShape(), temp_data2);

    }
    else
    {
        RealVectorType temp_data2;
        temp_data2.copyFromArray(w,1);
        if (isComplex())	// set real value in complex
        {
	    CplxVectorType temp_data3;
            fillComplexFromReal(temp_data2,temp_data3);
            m_data->setTaggedValue(tagKey,w.getShape(), temp_data3);
        }
        else
        {
            m_data->setTaggedValue(tagKey,w.getShape(), temp_data2);
        }
    }
}


void
Data::setTaggedValueFromCPP(int tagKey,
                            const DataTypes::ShapeType& pointshape,
                            const DataTypes::RealVectorType& value,
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

void
Data::setTaggedValueFromCPP(int tagKey,
                            const DataTypes::ShapeType& pointshape,
                            const DataTypes::CplxVectorType& value,
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
    if (transpose == 0)           {}
    else if (transpose == 1)      { start0 = axis_offset; }
    else if (transpose == 2)      { start1 = rank1-axis_offset; }
    else                          { throw DataException("C_GeneralTensorProduct: Error - transpose should be 0, 1 or 2"); }


    // Adjust the shapes for transpose
    DataTypes::ShapeType tmpShape0(rank0); // pre-sizing the vectors rather
    DataTypes::ShapeType tmpShape1(rank1); // than using push_back
    for (int i=0; i<rank0; i++)   { tmpShape0[i]=shape0[(i+start0)%rank0]; }
    for (int i=0; i<rank1; i++)   { tmpShape1[i]=shape1[(i+start1)%rank1]; }

    // Prepare for the loops of the product
    int SL=1, SM=1, SR=1;
    for (int i=0; i<rank0-axis_offset; i++) {
        SL *= tmpShape0[i];
    }
    for (int i=rank0-axis_offset; i<rank0; i++) {
        if (tmpShape0[i] != tmpShape1[i-(rank0-axis_offset)]) {
            throw DataException("C_GeneralTensorProduct: Error - incompatible shapes");
        }
        SM *= tmpShape0[i];
    }
    for (int i=axis_offset; i<rank1; i++) {
        SR *= tmpShape1[i];
    }

    // Define the shape of the output (rank of shape is the sum of the loop ranges below)
    DataTypes::ShapeType shape2(rank0+rank1-2*axis_offset);       
    { // block to limit the scope of out_index
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

    bool complexresult=arg_0_Z.isComplex() || arg_1_Z.isComplex();
    cplx_t dummyc=0;
    real_t dummyr=0;    
    if (arg_0_Z.isConstant() && arg_1_Z.isConstant()) {
        res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(), false);        // DataConstant output
	if (complexresult)
	{
	    res.complicate();
	}
	if (arg_0_Z.isComplex())
	{
	    if (arg_1_Z.isComplex())
	    {
		cplx_t dummyc=0;
	        res.complicate();
		const cplx_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(0,dummyc));
		const cplx_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(0,dummyc));
		cplx_t *ptr_2 = &(res.getDataAtOffsetRW(0,dummyc));
		matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);		
	    }
	    else	// right is real
	    {
		cplx_t dummyc=0;
		real_t dummyr=0;
	        res.complicate();
		const cplx_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(0,dummyc));
		const real_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(0,dummyr));
		cplx_t *ptr_2 = &(res.getDataAtOffsetRW(0,dummyc));
		matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);		
	    }
	}
	else	// arg_0_Z is real
	{
	    if (arg_1_Z.isComplex())
	    {
		cplx_t dummyc=0;
		real_t dummyr=0;
	        res.complicate();
		const real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(0,dummyr));
		const cplx_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(0,dummyc));
		cplx_t *ptr_2 = &(res.getDataAtOffsetRW(0,dummyc));
		matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);	      
	    }
	    else
	    {
	        real_t dummyr=0;
		const real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(0,0));
		const real_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(0,0));
		real_t *ptr_2 = &(res.getDataAtOffsetRW(0, dummyr));
		matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);	      
	    }
	}

    }
    else if (arg_0_Z.isConstant()   && arg_1_Z.isTagged()) {

        // Prepare the DataConstant input
        DataConstant* tmp_0=dynamic_cast<DataConstant*>(arg_0_Z.borrowData());
        if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataConstant."); }

        // Borrow DataTagged input from Data object
        DataTagged* tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());
        if (tmp_1==0) { throw DataException("GTP_1 Programming error - casting to DataTagged."); }

        // Prepare a DataTagged output 2
        res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(), false);        // DataTagged output
        res.tag();
	if (complexresult)
	{
	    res.complicate();
	}	
        DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());
        if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

	cplx_t dummyc=0;
	real_t dummyr=0;
        
	if (arg_0_Z.isComplex())
	{
	    if (arg_1_Z.isComplex())
	    {        
		// Prepare offset into DataConstant
		int offset_0 = tmp_0->getPointOffset(0,0);
		const cplx_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyc));

		const cplx_t *ptr_1 = &(tmp_1->getDefaultValueRO(0, dummyc));
		cplx_t *ptr_2 = &(tmp_2->getDefaultValueRW(0, dummyc));

		// Compute an MVP for the default
		matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		// Compute an MVP for each tag
		const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
		DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
		for (i=lookup_1.begin();i!=lookup_1.end();i++) {
		    tmp_2->addTag(i->first);

		    const cplx_t *ptr_1 = &(tmp_1->getDataByTagRO(i->first,0, dummyc));
		    cplx_t *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0, dummyc));
		
		    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		}		
	    }
	    else
	    {
		// Prepare offset into DataConstant
		int offset_0 = tmp_0->getPointOffset(0,0);
		const cplx_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyc));

		const real_t *ptr_1 = &(tmp_1->getDefaultValueRO(0, dummyr));
		cplx_t *ptr_2 = &(tmp_2->getDefaultValueRW(0, dummyc));

		// Compute an MVP for the default
		matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		// Compute an MVP for each tag
		const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
		DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
		for (i=lookup_1.begin();i!=lookup_1.end();i++) {
		    tmp_2->addTag(i->first);

		    const real_t *ptr_1 = &(tmp_1->getDataByTagRO(i->first,0, dummyr));
		    cplx_t *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0, dummyc));
		
		    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		}		      
	    }
	}
	else	// arg_0 is real
	{
	    if (arg_1_Z.isComplex())
	    {        
		// Prepare offset into DataConstant
		int offset_0 = tmp_0->getPointOffset(0,0);
		const real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyr));

		const cplx_t *ptr_1 = &(tmp_1->getDefaultValueRO(0, dummyc));
		cplx_t *ptr_2 = &(tmp_2->getDefaultValueRW(0,dummyc));	// the result

		// Compute an MVP for the default
		matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		// Compute an MVP for each tag
		const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
		DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
		for (i=lookup_1.begin();i!=lookup_1.end();i++) {
		    tmp_2->addTag(i->first);

		    const cplx_t *ptr_1 = &(tmp_1->getDataByTagRO(i->first,0,dummyc));
		    cplx_t *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0,dummyc));
		
		    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		}	      
	    }
	    else
	    {
		// Prepare offset into DataConstant
	        real_t dummyr=0;
		int offset_0 = tmp_0->getPointOffset(0,0);
		const real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyr));

		const real_t *ptr_1 = &(tmp_1->getDefaultValueRO(0));
		real_t *ptr_2 = &(tmp_2->getDefaultValueRW(0));

		// Compute an MVP for the default
		matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		// Compute an MVP for each tag
		const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
		DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
		for (i=lookup_1.begin();i!=lookup_1.end();i++) {
		    tmp_2->addTag(i->first);

		    const real_t *ptr_1 = &(tmp_1->getDataByTagRO(i->first,0));
		    real_t *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0));
		
		    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		}	      
	    }
	}
    }
    else if (arg_0_Z.isConstant() && arg_1_Z.isExpanded()) {

	res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
	if (complexresult)
	{
	    res.complicate();
	}	
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
        
	if (arg_0_Z.isComplex())
	{
	    if (arg_1_Z.isComplex())
	    {   
		#pragma omp parallel for private(sampleNo_1,dataPointNo_1) schedule(static)
		for (sampleNo_1 = 0; sampleNo_1 < numSamples_1; sampleNo_1++) {
		    for (dataPointNo_1 = 0; dataPointNo_1 < numDataPointsPerSample_1; dataPointNo_1++) {
			int offset_1 = tmp_1->getPointOffset(sampleNo_1,dataPointNo_1);
			int offset_2 = tmp_2->getPointOffset(sampleNo_1,dataPointNo_1);
			const cplx_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyc));
			const cplx_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyc));
			cplx_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyc));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}
	    }
	    else
	    {
		#pragma omp parallel for private(sampleNo_1,dataPointNo_1) schedule(static)
		for (sampleNo_1 = 0; sampleNo_1 < numSamples_1; sampleNo_1++) {
		    for (dataPointNo_1 = 0; dataPointNo_1 < numDataPointsPerSample_1; dataPointNo_1++) {
			int offset_1 = tmp_1->getPointOffset(sampleNo_1,dataPointNo_1);
			int offset_2 = tmp_2->getPointOffset(sampleNo_1,dataPointNo_1);
			const cplx_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyc));
			const real_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyr));
			cplx_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyc));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}	      
	    }
	}
	else
	{
	    if (arg_1_Z.isComplex())
	    {  
		#pragma omp parallel for private(sampleNo_1,dataPointNo_1) schedule(static)
		for (sampleNo_1 = 0; sampleNo_1 < numSamples_1; sampleNo_1++) {
		    for (dataPointNo_1 = 0; dataPointNo_1 < numDataPointsPerSample_1; dataPointNo_1++) {
			int offset_1 = tmp_1->getPointOffset(sampleNo_1,dataPointNo_1);
			int offset_2 = tmp_2->getPointOffset(sampleNo_1,dataPointNo_1);
			const real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyr));
			const cplx_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyc));
			cplx_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyc));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}	      
	    }
	    else	// both real
	    {
		#pragma omp parallel for private(sampleNo_1,dataPointNo_1) schedule(static)
		for (sampleNo_1 = 0; sampleNo_1 < numSamples_1; sampleNo_1++) {
		    for (dataPointNo_1 = 0; dataPointNo_1 < numDataPointsPerSample_1; dataPointNo_1++) {
			int offset_1 = tmp_1->getPointOffset(sampleNo_1,dataPointNo_1);
			int offset_2 = tmp_2->getPointOffset(sampleNo_1,dataPointNo_1);
			const real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyr));
			const real_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyr));
			real_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}	      
	    }	  
	}
    }
    else if (arg_0_Z.isTagged() && arg_1_Z.isConstant()) {
	// Borrow DataTagged input from Data object
	DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());
	if (tmp_0==0) { throw DataException("GTP_0 Programming error - casting to DataTagged."); }

	// Prepare the DataConstant input
	DataConstant* tmp_1=dynamic_cast<DataConstant*>(arg_1_Z.borrowData());
	if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataConstant."); }

	// Prepare a DataTagged output 2
	res = Data(0.0, shape2, arg_0_Z.getFunctionSpace(), false);        // DataTagged output
	res.tag();
	if (complexresult)
	{
	    res.complicate();
	}	
	DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());
	if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

	// Prepare offset into DataConstant
	int offset_1 = tmp_1->getPointOffset(0,0);      
	if (arg_0_Z.isComplex())
	{
	    if (arg_1_Z.isComplex())
	    {   

		const cplx_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyc));
		const cplx_t *ptr_0 = &(tmp_0->getDefaultValueRO(0, dummyc));
		cplx_t *ptr_2 = &(tmp_2->getDefaultValueRW(0, dummyc));

		// Compute an MVP for the default
		matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		// Compute an MVP for each tag
		const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
		DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
		for (i=lookup_0.begin();i!=lookup_0.end();i++) {

		    tmp_2->addTag(i->first);
		    const cplx_t *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0, dummyc));
		    cplx_t *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0, dummyc));
		    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		}	      
	    }
	    else	// arg_1_Z is real
	    {
		const real_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyr));
		const cplx_t *ptr_0 = &(tmp_0->getDefaultValueRO(0, dummyc));
		cplx_t *ptr_2 = &(tmp_2->getDefaultValueRW(0, dummyc));

		// Compute an MVP for the default
		matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		// Compute an MVP for each tag
		const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
		DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
		for (i=lookup_0.begin();i!=lookup_0.end();i++) {

		    tmp_2->addTag(i->first);
		    const cplx_t *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0, dummyc));
		    cplx_t *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0, dummyc));
		    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		}	      
	    }
	}
	else
	{
	    if (arg_1_Z.isComplex())
	    {   
		const cplx_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyc));
		const real_t *ptr_0 = &(tmp_0->getDefaultValueRO(0, dummyr));
		cplx_t *ptr_2 = &(tmp_2->getDefaultValueRW(0, dummyc));

		// Compute an MVP for the default
		matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		// Compute an MVP for each tag
		const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
		DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
		for (i=lookup_0.begin();i!=lookup_0.end();i++) {

		    tmp_2->addTag(i->first);
		    const real_t *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0, dummyr));
		    cplx_t *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0, dummyc));
		    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		}	      
	    }
	    else
	    {
		real_t dummyr=0;
		const real_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyr));
		const real_t *ptr_0 = &(tmp_0->getDefaultValueRO(0));
		real_t *ptr_2 = &(tmp_2->getDefaultValueRW(0));

		// Compute an MVP for the default
		matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		// Compute an MVP for each tag
		const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
		DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
		for (i=lookup_0.begin();i!=lookup_0.end();i++) {

		    tmp_2->addTag(i->first);
		    const real_t *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0));
		    real_t *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0));
		    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		}	      
	    }
	}
    }
    else if (arg_0_Z.isTagged() && arg_1_Z.isTagged()) {
        // Borrow DataTagged input from Data object
        DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());
        if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

        // Borrow DataTagged input from Data object
        DataTagged* tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());
        if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

        // Prepare a DataTagged output 2
        res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(), false);
        res.tag();  // DataTagged output
	if (complexresult)
	{
	    res.complicate();
	}
        DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());
        if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataTagged."); }
      
	if (arg_0_Z.isComplex())
	{
	    if (arg_1_Z.isComplex())
	    {   
		const cplx_t *ptr_0 = &(tmp_0->getDefaultValueRO(0, dummyc));
		const cplx_t *ptr_1 = &(tmp_1->getDefaultValueRO(0, dummyc));
		cplx_t *ptr_2 = &(tmp_2->getDefaultValueRW(0, dummyc));

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
		    const cplx_t *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0, dummyc));
		    const cplx_t *ptr_1 = &(tmp_1->getDataByTagRO(i->first,0, dummyc));
		    cplx_t *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0, dummyc));

		    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		}	      
	    }
	    else	// arg_1_Z is real
	    {
		const cplx_t *ptr_0 = &(tmp_0->getDefaultValueRO(0, dummyc));
		const real_t *ptr_1 = &(tmp_1->getDefaultValueRO(0, dummyr));
		cplx_t *ptr_2 = &(tmp_2->getDefaultValueRW(0, dummyc));

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
		    const cplx_t *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0, dummyc));
		    const real_t *ptr_1 = &(tmp_1->getDataByTagRO(i->first,0, dummyr));
		    cplx_t *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0, dummyc));

		    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		}	      
	    }
	}
	else	// arg_0_Z is real
	{
	    if (arg_1_Z.isComplex())
	    {   
		const real_t *ptr_0 = &(tmp_0->getDefaultValueRO(0, dummyr));
		const cplx_t *ptr_1 = &(tmp_1->getDefaultValueRO(0, dummyc));
		cplx_t *ptr_2 = &(tmp_2->getDefaultValueRW(0, dummyc));

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
		    const real_t *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0, dummyr));
		    const cplx_t *ptr_1 = &(tmp_1->getDataByTagRO(i->first,0, dummyc));
		    cplx_t *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0,dummyc));

		    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		}		      
	    }
	    else
	    {
		const real_t *ptr_0 = &(tmp_0->getDefaultValueRO(0));
		const real_t *ptr_1 = &(tmp_1->getDefaultValueRO(0));
		real_t *ptr_2 = &(tmp_2->getDefaultValueRW(0));

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
		    const real_t *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0));
		    const real_t *ptr_1 = &(tmp_1->getDataByTagRO(i->first,0));
		    real_t *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0));

		    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		}	      
	    }
	}
    }
    else if (arg_0_Z.isTagged() && arg_1_Z.isExpanded()) {
      
	// After finding a common function space above the two inputs have the same numSamples and num DPPS
	res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
	if (complexresult)
	{
	    res.complicate();
	}
	DataTagged*   tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());
	DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
	DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());    
	if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataTagged."); }
	if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
	if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
	int sampleNo_0,dataPointNo_0;
	int numSamples_0 = arg_0_Z.getNumSamples();
	int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();	
	if (arg_0_Z.isComplex())
	{
	    if (arg_1_Z.isComplex())
	    {   
		#pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
		for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
		    int offset_0 = tmp_0->getPointOffset(sampleNo_0,0); // They're all the same, so just use #0
		    const cplx_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyc));
		    for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
			int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
			const cplx_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyc));
			cplx_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyc));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}	      
	    }
	    else	// arg_1_Z is real
	    {
		#pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
		for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
		    int offset_0 = tmp_0->getPointOffset(sampleNo_0,0); // They're all the same, so just use #0
		    const cplx_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyc));
		    for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
			int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
			const real_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyr));
			cplx_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyc));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}	      
	    }
	}
	else	// arg_0_Z is real
	{
	    if (arg_1_Z.isComplex())
	    {   
		#pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
		for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
		    int offset_0 = tmp_0->getPointOffset(sampleNo_0,0); // They're all the same, so just use #0
		    const real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyr));
		    for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
			int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
			const cplx_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyc));
			cplx_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyc));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}	      
	    }
	    else
	    {
		#pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
		for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
		    int offset_0 = tmp_0->getPointOffset(sampleNo_0,0); // They're all the same, so just use #0
		    const real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyr));
		    for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
			int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
			const real_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyr));
			real_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}	      
	    }
	}
    }
    else if (arg_0_Z.isExpanded() && arg_1_Z.isConstant()) {
        res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
	if (complexresult)
	{
	    res.complicate();
	}
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
	if (arg_0_Z.isComplex())
	{
	    if (arg_1_Z.isComplex())
	    {   
		#pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
		for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
		    for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
			int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
			const cplx_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyc));
			const cplx_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyc));
			cplx_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyc));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}		      
	    }
	    else	// arg_1_Z is real
	    {
		#pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
		for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
		    for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
			int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
			const cplx_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyc));
			const real_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyr));
			cplx_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyc));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}	      
	    }
	}
	else	// arg_0_Z is real
	{
	    if (arg_1_Z.isComplex())
	    {   
		#pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
		for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
		    for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
			int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
			const real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyr));
			const cplx_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyc));
			cplx_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyc));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }	      
		}
	    }
	    else
	    {
		#pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
		for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
		    for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
			int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
			const real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyr));
			const real_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyr));
			real_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}	      
	    }
	}


    }
    else if (arg_0_Z.isExpanded() && arg_1_Z.isTagged()) {
        // After finding a common function space above the two inputs have the same numSamples and num DPPS
        res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
	if (complexresult)
	{
	    res.complicate();
	}
        DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
        DataTagged*   tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());
        DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());
        if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
        if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataTagged."); }
        if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
        int sampleNo_0,dataPointNo_0;
        int numSamples_0 = arg_0_Z.getNumSamples();
        int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();      
	if (arg_0_Z.isComplex())
	{
	    if (arg_1_Z.isComplex())
	    {   
		#pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
		for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
		    int offset_1 = tmp_1->getPointOffset(sampleNo_0,0);
		    const cplx_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyc));
		    for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
			int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
			const cplx_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyc));
			cplx_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyc));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}	      
	    }
	    else
	    {
		#pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
		for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
		    int offset_1 = tmp_1->getPointOffset(sampleNo_0,0);
		    const real_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyr));
		    for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
			int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
			const cplx_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyc));
			cplx_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyc));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}	      
	    }
	}
	else
	{
	    if (arg_1_Z.isComplex())
	    {   
		#pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
		for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
		    int offset_1 = tmp_1->getPointOffset(sampleNo_0,0);
		    const cplx_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyc));
		    for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
			int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
			const real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyr));
			cplx_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyc));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}	      
	    }
	    else
	    {
		#pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
		for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
		    int offset_1 = tmp_1->getPointOffset(sampleNo_0,0);
		    const real_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyr));
		    for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
			int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
			const real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyr));
			real_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}	      
	    }
	}


    }
    else if (arg_0_Z.isExpanded() && arg_1_Z.isExpanded()) {
        // After finding a common function space above the two inputs have the same numSamples and num DPPS
        res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
	if (complexresult)
	{
	    res.complicate();
	}
        DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
        DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
        DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());
        if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
        if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
        if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
        int sampleNo_0,dataPointNo_0;
        int numSamples_0 = arg_0_Z.getNumSamples();
        int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();      
	if (arg_0_Z.isComplex())
	{
	    if (arg_1_Z.isComplex())
	    {   
		#pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
		for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
		    for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
			int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
			const cplx_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyc));
			const cplx_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyc));
			cplx_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyc));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}	      
	    }
	    else
	    {
		#pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
		for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
		    for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
			int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
			const cplx_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyc));
			const real_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyr));
			cplx_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyc));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}	      
	    }
	}
	else
	{
	    if (arg_1_Z.isComplex())
	    {   
		#pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
		for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
		    for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
			int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
			const real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyr));
			const cplx_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyc));
			cplx_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyc));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}	      
	    }
	    else
	    {
		#pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
		for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
		    for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
			int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
			int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
			const real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyr));
			const real_t *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummyr));
			real_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
			matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
		    }
		}	      
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
    DataReady_ptr dr=REFCOUNTNS::dynamic_pointer_cast<DataReady>(m_data);
    ESYS_ASSERT(dr!=NULL, "Casting to DataReady.");
    return dr;
}

std::string
Data::toString() const
{
    int localNeedSummary=0;
#ifdef ESYS_MPI
    int globalNeedSummary=0;
#endif
    if (!m_data->isEmpty() &&
        m_data->isExpanded() &&
        !m_data->isLazy() && 
        getLength() > escriptParams.getTooManyLines())
    {
        localNeedSummary=1;
    }

#ifdef ESYS_MPI
    MPI_Allreduce( &localNeedSummary, &globalNeedSummary, 1, MPI_INT, MPI_MAX, get_MPIComm() );
    localNeedSummary=globalNeedSummary;
#endif

    if (localNeedSummary){
        if (isComplex())
    	{
    	    stringstream temp;
    	    temp << "Summary: Lsup="<< Lsup_const() << " data points=" << getNumDataPoints();
    	    return  temp.str();
    	}
    	else
    	{
    	    stringstream temp;
    	    temp << "Summary: inf="<< inf_const() << " sup=" << sup_const() << " data points=" << getNumDataPoints();
    	    return  temp.str();
    	}
    }
    return m_data->toString();
}


// This method is not thread-safe
DataTypes::RealVectorType::reference
Data::getDataAtOffsetRW(DataTypes::RealVectorType::size_type i, DataTypes::real_t dummy)
{
    checkExclusiveWrite();
    return getReady()->getDataAtOffsetRW(i);
}

// This method is not thread-safe
DataTypes::RealVectorType::const_reference
Data::getDataAtOffsetRO(DataTypes::RealVectorType::size_type i, DataTypes::real_t dummy)
{
    forceResolve();
    return getReady()->getDataAtOffsetRO(i);
}

// This method is not thread-safe
DataTypes::CplxVectorType::reference
Data::getDataAtOffsetRW(DataTypes::RealVectorType::size_type i, DataTypes::cplx_t dummy)
{
    checkExclusiveWrite();
    return getReady()->getDataAtOffsetRWC(i);
}

// This method is not thread-safe
DataTypes::CplxVectorType::const_reference
Data::getDataAtOffsetRO(DataTypes::RealVectorType::size_type i, DataTypes::cplx_t dummy)
{
    forceResolve();
    return getReady()->getDataAtOffsetROC(i);
}

DataTypes::RealVectorType::const_reference
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


DataTypes::RealVectorType::reference
Data::getDataPointRW(int sampleNo, int dataPointNo)
{
    checkExclusiveWrite();
    DataReady* dr=getReady();
    return dr->getDataAtOffsetRW(dr->getPointOffset(sampleNo, dataPointNo));
}

Data
Data::interpolateFromTable3DP(bp::object table, real_t Amin, real_t Astep,
                Data& B, real_t Bmin, real_t Bstep,
                Data& C, real_t Cmin, real_t Cstep,
                real_t undef, bool check_boundaries)
{
    WrappedArray t(table);
    return interpolateFromTable3D(t, Amin, Astep, undef, B, Bmin, Bstep, C, Cmin, Cstep, check_boundaries);
}

Data
Data::interpolateFromTable2DP(bp::object table, real_t Amin, real_t Astep,
                Data& B, real_t Bmin, real_t Bstep,
                real_t undef, bool check_boundaries)
{
    WrappedArray t(table);
    return interpolateFromTable2D(t, Amin, Astep, undef, B, Bmin, Bstep,check_boundaries);
}

Data
Data::interpolateFromTable1DP(bp::object table, real_t Amin, real_t Astep,
                              real_t undef,bool check_boundaries)
{
    WrappedArray t(table);
    return interpolateFromTable1D(t, Amin, Astep, undef, check_boundaries);
}


Data
Data::interpolateFromTable1D(const WrappedArray& table, real_t Amin,
                             real_t Astep, real_t undef, bool check_boundaries)
{
    table.convertArray(); // critical! Calling getElt on an unconverted array is not thread safe
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
    int numpts=getNumDataPoints();
    int twidth=table.getShape()[0]-1;       
    bool haserror=false;
    const RealVectorType* adat=0;
    RealVectorType* rdat=0;
    try
    {
        adat=&(getReady()->getVectorRO());
        rdat=&(res.getReady()->getVectorRW());
    } catch (DataException& d)
    {
        haserror=true;
        error=3;
    }
    if (!haserror)
    {
        int l=0;
#pragma omp parallel for private(l) schedule(static)
        for (l=0;l<numpts; ++l)
        {
            int lerror=0;
#pragma omp flush(haserror) // In case haserror was in register
            if (!haserror)                
            {
                real_t a=(*adat)[l];
                int x=static_cast<int>(((a-Amin)/Astep));
                if (check_boundaries)
                {
                    if ((a<Amin) || (x<0))
                    {
                        lerror=1;
                    } 
                    else if (a>Amin+Astep*twidth) 
                    {
                        lerror=4;
                    }
                }
                if (!lerror)
                {
                    if (x<0) x=0;
                    if (x>twidth) x=twidth;
                    try {
                        if (x==twidth) // value is on the far end of the table
                        {
                            real_t e=table.getElt(x);
                            if (e>undef)
                            {
                                lerror=2;
                            }
                            else
                            {
                                (*rdat)[l]=e;
                            }
                        }
                        else            // x and y are in bounds
                        {
                            real_t e=table.getElt(x);
                            real_t w=table.getElt(x+1);
                            if ((e>undef) || (w>undef))
                            {
                                lerror=2;
                            }
                            else
                            {
                                // map x*Astep <= a << (x+1)*Astep to [-1,1] 
                                real_t la = 2.0*(a-Amin-(x*Astep))/Astep-1;
                                (*rdat)[l]=((1-la)*e + (1+la)*w)/2;
                            }
                        }
                    }
                    catch (DataException& d)
                    {
                        lerror=3;
                    }   
                } // if !lerror
                if (lerror!=0)
                {
#pragma omp critical // Doco says there is a flush associated with critical
                    {
                        haserror=true; // We only care that error is recorded.
                        error=lerror;  // We don't care which one.
                    }
                }
            } // if (!error)
        }   // parallelised for
    } // if !haserror
#ifdef ESYS_MPI
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
Data::interpolateFromTable2D(const WrappedArray& table, real_t Amin,
                             real_t Astep, real_t undef, Data& B, real_t Bmin,
                             real_t Bstep, bool check_boundaries)
{
    table.convertArray(); // critical! Calling getElt on an unconverted array is not thread safe
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
        throw DataException("All step components must be strictly positive.");
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

    int numpts=getNumDataPoints();
    const RealVectorType* adat=0;
    const RealVectorType* bdat=0;
    RealVectorType* rdat=0;
    const DataTypes::ShapeType& ts=table.getShape();
    try
    {
        adat=&(getReady()->getVectorRO());
        bdat=&(B.getReady()->getVectorRO());
        rdat=&(res.getReady()->getVectorRW());
    }
    catch (DataException& e)
    {
        error=3;
    }
    if (!error)
    {
        int twx=ts[1]-1;        // table width x
        int twy=ts[0]-1;        // table width y

        bool haserror=false;
        int l=0;
#pragma omp parallel for private(l) shared(res,rdat, adat, bdat) schedule(static) 
        for (l=0; l<numpts; ++l)
        {
#pragma omp flush(haserror) // In case haserror was in register
           if (!haserror)               
           {
                int lerror=0;
                real_t a=(*adat)[l];
                real_t b=(*bdat)[l];
                int x=static_cast<int>(((a-Amin)/Astep));
                int y=static_cast<int>(((b-Bmin)/Bstep));
                if (check_boundaries)
                {
                    if ( (a<Amin) || (b<Bmin) || (x<0) || (y<0))
                    {
                        lerror=1;
                    }
                    else if ( (a>Amin+Astep*twx) || (b>Bmin+Bstep*twy))
                    {
                        lerror=4;
                    }
                } 
                if (lerror==0)
                {
                    if (x<0) x=0;
                    if (y<0) y=0;

                    if (x>twx) x=twx;
                    if (y>twy) y=twy;
                    try
                    {
                        int nx=x+1;
                        int ny=y+1;

                        real_t la=0; // map position of a between x and nx to [-1,1]
                        real_t lb=0;
                        real_t weight=4;

                        // now we work out which terms we should be considering
                        bool usex=(x!=twx);
                        bool usey=(y!=twy);

                        la = 2.0*(a-Amin-(x*Astep))/Astep-1;
                        lb = 2.0*(b-Bmin-(y*Bstep))/Bstep-1;

                        real_t sw=table.getElt(y,x);
                        real_t nw=usey?table.getElt(ny,x):0; // 0 because if !usey ny does not actually exist
                        real_t se=usex?table.getElt(y,nx):0;
                        real_t ne=(usex&&usey)?table.getElt(ny,nx):0;                

                        real_t ans=(1-la)*(1-lb)*sw +
                                   (1-la)*(1+lb)*nw +
                                   (1+la)*(1-lb)*se +
                                   (1+la)*(1+lb)*ne;
                        ans/=weight;
                        (*rdat)[l]=ans;
                        // this code does not check to see if any of the points used in the interpolation are undef
                        if (ans>undef)
                        {
                            lerror=2;
                        }
                    }
                    catch (DataException& d)
                    {
                        lerror=3;
                    }
                }
                if (lerror!=0)
                {
#pragma omp critical  // Doco says there is a flush associated with critical
                    {
                        error=lerror;
                    }               
                }
            }  // if (!haserror)
        } // parallel for
    } // !error
#ifdef ESYS_MPI
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
Data::interpolateFromTable3D(const WrappedArray& table, real_t Amin,
                             real_t Astep, real_t undef, Data& B, real_t Bmin,
                             real_t Bstep, Data& C, real_t Cmin, real_t Cstep,
                             bool check_boundaries)
{
    table.convertArray(); // critical! Calling getElt on an unconverted array is not thread safe
    int error=0;
    if ((getDataPointRank()!=0) || (B.getDataPointRank()!=0) || (C.getDataPointRank()!=0))
    {
        throw DataException("Inputs to 3D interpolation must be scalar");
    }
    if (table.getRank()!=3)
    {
        throw DataException("Table for 3D interpolation must be 3D");
    }
    if ((Astep<=0) || (Bstep<=0) || (Cstep<=0))
    {
        throw DataException("All step components must be strictly positive.");
    }
    if (getFunctionSpace()!=B.getFunctionSpace())
    {
        Data n=B.interpolate(getFunctionSpace());
        return interpolateFromTable3D(table, Amin, Astep, undef, 
                n , Bmin, Bstep, C, Cmin, Cstep, check_boundaries);
    }
    if (getFunctionSpace()!=C.getFunctionSpace())
    {
        Data n=C.interpolate(getFunctionSpace());
        return interpolateFromTable3D(table, Amin, Astep, undef, 
                B , Bmin, Bstep, n, Cmin, Cstep, check_boundaries);
    }

    if (!isExpanded())
    {
        expand();
    }
    if (!B.isExpanded())
    {
        B.expand();
    }
    if (!C.isExpanded())
    {
        C.expand();
    }

    Data res(0, DataTypes::scalarShape, getFunctionSpace(), true);

    int numpts=getNumDataPoints();
    const RealVectorType* adat=0;
    const RealVectorType* bdat=0;
    const RealVectorType* cdat=0;
    RealVectorType* rdat=0;
    const DataTypes::ShapeType& ts=table.getShape();
    try
    {
        adat=&(getReady()->getVectorRO());
        bdat=&(B.getReady()->getVectorRO());
        cdat=&(C.getReady()->getVectorRO());
        rdat=&(res.getReady()->getVectorRW());
    }
    catch (DataException& e)
    {
        error=3;
    }
    if (!error)
    {
        int twx=ts[2]-1;        // table width x
        int twy=ts[1]-1;        // table width y
        int twz=ts[0]-1;        // table width z

        bool haserror=false;
        int l=0;
#pragma omp parallel for private(l) shared(res,rdat, adat, bdat) schedule(static) 
        for (l=0; l<numpts; ++l)
        {
#pragma omp flush(haserror) // In case haserror was in register
           if (!haserror)               
           {
                int lerror=0;
                real_t a=(*adat)[l];
                real_t b=(*bdat)[l];
                real_t c=(*cdat)[l];
                int x=static_cast<int>(((a-Amin)/Astep));
                int y=static_cast<int>(((b-Bmin)/Bstep));
                int z=static_cast<int>(((c-Cmin)/Cstep));
                if (check_boundaries)
                {
                    if ( (a<Amin) || (b<Bmin) || (c<Cmin)|| (x<0) || (y<0) || (z<0))
                    {
                        lerror=1;
                    }
                    else if ( (a>Amin+Astep*twx) || (b>Bmin+Bstep*twy) || (c>Cmin+Cstep*twz))
                    {
                        lerror=4;
                    }
                } 
                if (lerror==0)
                {
                    if (x<0) x=0;
                    if (y<0) y=0;
                    if (z<0) z=0;

                    if (x>twx) x=twx;
                    if (y>twy) y=twy;
                    if (z>twz) z=twz;
                    try
                    {
                        int nx=x+1;
                        int ny=y+1;
                        int nz=z+1;
                        real_t la=0; // map position of a between x and nx to [-1,1]
                        real_t lb=0;
                        real_t lc=0;
                        real_t weight=8;

                        // now we work out which terms we should be considering
                        bool usex=(x!=twx);
                        bool usey=(y!=twy);
                        bool usez=(z!=twz);

                        la = 2.0*(a-Amin-(x*Astep))/Astep-1;
                        lb = 2.0*(b-Bmin-(y*Bstep))/Bstep-1;
                        lc = 2.0*(c-Cmin-(z*Cstep))/Cstep-1;

                        real_t swb=table.getElt(z,y,x);
                        real_t swt=usez?table.getElt(nz,y,x):0;
                        real_t nwb=usey?table.getElt(z,ny,x):0;
                        real_t nwt=(usey&&usez)?table.getElt(nz,ny,x):0;
                        real_t seb=usex?table.getElt(z,y,nx):0;
                        real_t set=(usex&&usez)?table.getElt(nz,y,nx):0;
                        real_t neb=(usex&&usey)?table.getElt(z,ny,nx):0;
                        real_t net=(usex&&usey&&usez)?table.getElt(nz,ny,nx):0;

                        real_t ans=(1-la)*(1-lb)*(1-lc)*swb +
                                   (1-la)*(1-lb)*(1+lc)*swt +
                                   (1-la)*(1+lb)*(1-lc)*nwb +
                                   (1-la)*(1+lb)*(1+lc)*nwt +
                                   (1+la)*(1-lb)*(1-lc)*seb +
                                   (1+la)*(1-lb)*(1+lc)*set +
                                   (1+la)*(1+lb)*(1-lc)*neb +
                                   (1+la)*(1+lb)*(1+lc)*net;
                        ans/=weight;
                        (*rdat)[l]=ans;
                        // this code does not check to see if any of the points used in the interpolation are undef
                        if (ans>undef)
                        {
                            lerror=2;
                        }
                    }
                    catch (DataException& d)
                    {
                        lerror=3;
                    }
                }
                if (lerror!=0)
                {
#pragma omp critical // Doco says there is a flush associated with critical
                    {
                        error=lerror;
                    }               
                }
            }  // if (!haserror)
        } // parallel for
    } // !error
#ifdef ESYS_MPI
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

Data Data::nonuniforminterp(boost::python::object in, boost::python::object out, bool check_boundaries)
{
    WrappedArray win(in);
    win.convertArray();
    WrappedArray wout(out);
    wout.convertArray();
    if ((win.getRank()!=1) || (wout.getRank()!=1) || (win.getShape()[0]<1))
    {
        throw DataException("Input and output must be arrays/lists of scalars");
    }
    if (win.getShape()!=win.getShape())
    {
        throw DataException("Input and output must contain the same number of items.");
    }
    if (getDataPointRank()!=0)
    {
        throw DataException("The data being interpolated must be scalar.");
    }
    // now create an object the same size as this one
    // We'll expand it later if we need to
    expand();
    Data result(0, DataTypes::scalarShape, getFunctionSpace(), true);  
    int numpts=getNumDataPoints();
    const RealVectorType& sdat=getReady()->getVectorRO();
    RealVectorType& rdat=result.getReady()->getVectorRW();
    real_t maxlimit=win.getElt(win.getShape()[0]-1);
    real_t maxout=wout.getElt(wout.getShape()[0]-1);
    int ipoints=win.getShape()[0];
    int l=0;
    bool error=false;
    #pragma omp parallel for private(l) schedule(static) 
    for (l=0; l<numpts; ++l)
    {
        if ((sdat)[l]<win.getElt(0))
        {
           if (check_boundaries)
           {
               error=true;              // Could have done an early exit but I'm not sure it's worth it
           }
           else
           {
               rdat[l]=wout.getElt(0);
           }
        }
        else if (sdat[l]>maxlimit)
        {
           if (check_boundaries)
           {
               error=true;              // Could have done an early exit but I'm not sure it's worth it
           }
           else
           {
               rdat[l]=maxout;
           }
        }
        else
        {
            int i=0;
            for (;i<ipoints-2;++i)
            {
                if (sdat[l]<win.getElt(i+1))
                {
                    break;
                }
            }
            // we must have found one by this point or we would have triggered earlier branches
            rdat[l]=(wout.getElt(i+1)-wout.getElt(i))/(win.getElt(i+1)-win.getElt(i)) * (sdat[l]-win.getElt(i)) + wout.getElt(i);
        }
    }
    if (error)  // we had an illegal value (below the start threshold)
    {
        throw DataException("Data being interpolated contains a value outside the range given.");
    }
    return result;
}

Data Data::nonuniformslope(boost::python::object in, boost::python::object out, bool check_boundaries)
{
    WrappedArray win(in);
    win.convertArray();
    WrappedArray wout(out);
    wout.convertArray();
    if ((win.getRank()!=1) || (wout.getRank()!=1) || (win.getShape()[0]<1))
    {
        throw DataException("Input and output must be arrays/lists of scalars");
    }
    if (win.getShape()!=win.getShape())
    {
        throw DataException("Input and output must contain the same number of items.");
    }
    if (getDataPointRank()!=0)
    {
        throw DataException("The data being interpolated must be scalar.");
    }
    // now create an object the same size as this one
    // We'll expand it later if we need to
    expand();
    Data result(0, DataTypes::scalarShape, getFunctionSpace(), true);  
    int numpts=getNumDataPoints();
    const RealVectorType& sdat=getReady()->getVectorRO();
    RealVectorType& rdat=result.getReady()->getVectorRW();
    real_t maxlimit=win.getElt(win.getShape()[0]-1);
    int ipoints=win.getShape()[0];
    int l=0;
    bool error=false;
    #pragma omp parallel for private(l) schedule(static) 
    for (l=0; l<numpts; ++l)
    {
        if ((sdat)[l]<win.getElt(0))
        {
           if (check_boundaries)
           {
               error=true;              // Could have done an early exit but I'm not sure it's worth it
           }
           else
           {
               rdat[l]=0;
           }
        }
        else if (sdat[l]>maxlimit)
        {
           if (check_boundaries)
           {
               error=true;              // Could have done an early exit but I'm not sure it's worth it
           }
           else
           {
               rdat[l]=0;
           }
        }
        else
        {
            int i=0;
            for (;i<ipoints-2;++i)
            {
                if (sdat[l]<=win.getElt(i+1))
                {
                    break;
                }
            }
            // we must have found one by this point or we would have triggered earlier branches
            rdat[l]=(wout.getElt(i+1)-wout.getElt(i))/(win.getElt(i+1)-win.getElt(i));
        }
    }
    if (error)  // we had an illegal value (below the start threwshold)
    {
        throw DataException("Data being interpolated contains a value outside the range given.");
    }
    return result;
}


/* Member functions specific to the MPI implementation */

void
Data::print()
{
    int i,j;

    printf( "Data is %dX%d\n", getNumSamples(), getNumDataPointsPerSample() );
    if (isComplex())
    {
	for( i=0; i<getNumSamples(); i++ )
	{
	    printf( "[%6d]", i );
	    for( j=0; j<getNumDataPointsPerSample(); j++ )
	    {
	        DataTypes::cplx_t* v=getSampleDataRW(i, static_cast<DataTypes::cplx_t>(0));
		printf( "\t%10.7g,%10.7g", std::real(v[j]), std::imag(v[j]) );    // doesn't really need RW access
	    }
	    printf( "\n" );
	}      
    }
    else
    {
	for( i=0; i<getNumSamples(); i++ )
	{
	    printf( "[%6d]", i );
	    for( j=0; j<getNumDataPointsPerSample(); j++ )
		printf( "\t%10.7g", (getSampleDataRW(i, static_cast<DataTypes::real_t>(0)))[j] );    // doesn't really need RW access
	    printf( "\n" );
	}
    }
}

void
Data::dump(const std::string fileName) const
{
     if (isLazy())
     {
         Data temp(*this);     // this is to get a non-const object which we can resolve
         temp.resolve();
         temp.dump(fileName);
     }
     else
     {
     #ifdef ESYS_HAVE_HDF5
        try
        {
            JMPI mpiInfo(getFunctionSpace().getDomain()->getMPI());
            const std::string newFileName(mpiInfo->appendRankToFileName(fileName));
            H5::H5File h5_file(newFileName, H5F_ACC_TRUNC);

            if (isComplex()) {
                Data temp_re = (*this).real();
                Data temp_im = (*this).imag();
                temp_re.m_data->dump_hdf5(h5_file.createGroup("Data_Re"));
                temp_im.m_data->dump_hdf5(h5_file.createGroup("Data_Im"));
            } else {
                m_data->dump_hdf5(h5_file.createGroup("Data"));
            }
        }
        catch (std::exception& e)
        {
            std::cout << e.what() << std::endl;
        }
    #else
        throw DataException("Data::dump_hdf5: not configured with HDF5. Please contact your installation manager.");
    #endif
    }
}

int
Data::get_MPISize() const
{
    int size;
#ifdef ESYS_MPI
    MPI_Comm comm = get_MPIComm();
    if (comm == MPI_COMM_NULL) {
        size = 1;  // NullDomain uses MPI_COMM_NULL, treat as single process
    } else {
        /*int error =*/ MPI_Comm_size( comm, &size );
    }
#else
    size = 1;
#endif
    return size;
}

int
Data::get_MPIRank() const
{
    int rank;
#ifdef ESYS_MPI
    MPI_Comm comm = get_MPIComm();
    if (comm == MPI_COMM_NULL) {
        rank = 0;  // NullDomain uses MPI_COMM_NULL, treat as rank 0
    } else {
        /*int error =*/ MPI_Comm_rank( comm, &rank );
    }
#else
    rank = 0;
#endif
    return rank;
}

MPI_Comm
Data::get_MPIComm() const
{
#ifdef ESYS_MPI
        MPI_Comm comm = getDomain()->getMPIComm();
        // NullDomain uses MPI_COMM_NULL, return MPI_COMM_SELF for safe single-process operations
        if (comm == MPI_COMM_NULL) {
            return MPI_COMM_SELF;
        }
        return comm;
#else
        return -1;
#endif
}

#ifdef IKNOWWHATIMDOING

// Considered having a generic option argument for extra info for function.
// If you pass a python object in we may have threading issues, if it's a C object then it's not python friendly.
// It's better if the library supplying the C function has its own interface for doing configuration
ESCRIPT_DLL_API
Data
escript::applyBinaryCFunction(bp::object cfunc, bp::tuple shape, escript::Data& ind, escript::Data& ine)
{
    int err=255;
    PyObject* p=cfunc.ptr();
    if (!PyCObject_Check(p))
    {
        throw DataException("applyBinaryCFunction: function must be a PyCObject.");
    }
    void* v=PyCObject_AsVoidPtr(p);
    binOpFnPtr func=binOpFnPtrFromVoidPtr(v);
    Data d(ind),e(ine);
    if (ind.getFunctionSpace()!=ine.getFunctionSpace())
    {
        if (ind.getDomain()!=ine.getDomain())
        {
            throw DataException("applyBinaryCFunction: can't interpolate between domains.");
        }
        std::vector<int> fstypes(2);
        fstypes[0]=ind.getFunctionSpace().getTypeCode();
        fstypes[1]=ine.getFunctionSpace().getTypeCode();
        int bestfs;
        ind.getDomain()->commonFunctionSpace(fstypes, bestfs);
        FunctionSpace best(ind.getDomain(),bestfs);
        d=ind.interpolate(best);
        e=ine.interpolate(best);
    }
    // now we find the result shape
    DataTypes::ShapeType resultshape;
    for (int i = 0; i < shape.attr("__len__")(); ++i) {
        resultshape.push_back(bp::extract<const int>(shape[i]));
    }

    if (d.isLazy() && !d.actsExpanded())
    {
        d.resolve(); // If you aren't expanded you probably won't get benefit from lazy anyway
    }
    if (e.isLazy() && !e.actsExpanded())
    {
        e.resolve();
    }
    Data res(0,resultshape,d.getFunctionSpace());
    int dpointsize=d.getNoValues();
    int epointsize=e.getNoValues();
    int rpointsize=res.getNoValues();
    if (d.actsExpanded() && !e.actsExpanded())
    {
        e.expand();
    }
    else if (!d.actsExpanded() && e.actsExpanded())
    {
        d.expand();
    }
    else if (d.isTagged() && e.isConstant())
    {
        e.tag();
    }
    else if (e.isTagged() && d.isConstant())
    {
        d.tag();
    }
    if (d.isConstant() && e.isConstant())
    {
        const real_t* src=d.getSampleDataRO(0);
        const real_t* src2=e.getSampleDataRO(0);
        real_t* dest=res.getSampleDataRW(0);
        err=func(dest,src,src2,rpointsize, dpointsize, epointsize);
    }
    else if (d.isTagged() && e.isTagged())
    {
        res.tag();
        DataTagged& srcd=*dynamic_cast<DataTagged*>(d.m_data.get());
        DataTagged& srce=*dynamic_cast<DataTagged*>(e.m_data.get());
        DataTagged& destd=*dynamic_cast<DataTagged*>(res.m_data.get());
        std::list<int> alltags;
        const DataTagged::DataMapType& srcLookupd=srcd.getTagLookup();
        DataTagged::DataMapType::const_iterator i;
        DataTagged::DataMapType::const_iterator srcLookupEnd=srcLookupd.end();
        for (i=srcLookupd.begin();i!=srcLookupEnd;i++)
        {
           alltags.push_back(i->first);
        }               
        const DataTagged::DataMapType& srcLookupe=srce.getTagLookup();
        srcLookupEnd=srcLookupe.end();
        for (i=srcLookupe.begin();i!=srcLookupEnd;i++)
        {
            if (find(alltags.begin(), alltags.end(), i->first)==alltags.end())   // we have already seen the tag
            {
                alltags.push_back(i->first);
            }
        }
        err=0;
        // now all tags will be a complete list of tags from both inputs
        for (std::list<int>::iterator j=alltags.begin();(j!=alltags.end()) && (err==0);++j)
        {
            destd.addTag(*j);
            const real_t *ptr_0 = &(srcd.getDataByTagRO(*j,0));
            const real_t *ptr_1 = &(srce.getDataByTagRO(*j,0));
            real_t *ptr_2 = &(destd.getDataByTagRW(*j,0));
            err=func(ptr_2,ptr_0,ptr_1,rpointsize, dpointsize, epointsize);
        }
        if (err==0)
        {
            // now we do the default tag
            const real_t *ptr_0 = &(srcd.getDefaultValueRO(0));
            const real_t *ptr_1 = &(srce.getDefaultValueRO(0));
            real_t *ptr_2 = &(destd.getDefaultValueRW(0));
            err=func(ptr_2,ptr_0,ptr_1,rpointsize, dpointsize, epointsize);
        }
    }
    else if (e.actsExpanded() && d.actsExpanded())
    {
        res.expand();
        int numsamples=d.getNumSamples();
        
        int dpps=d.getNumDataPointsPerSample();
        err=0;
#pragma omp parallel shared(err)
        {
           int localerr=0;
           int sampleid;
#pragma omp for schedule(dynamic)
           for (sampleid=0;sampleid<numsamples;++sampleid)
           {
                if(!localerr)
                {
                    const real_t* src=d.getSampleDataRO(sampleid);
                    const real_t* src2=e.getSampleDataRO(sampleid);
                    real_t* dest=res.getSampleDataRW(sampleid);
                    for (int pointnum=0;pointnum<dpps;++pointnum)
                    {
                        localerr=func(dest,src,src2,rpointsize, dpointsize, epointsize);
                        if (localerr!=0)
                        {
                            break;
                        }
                        src+=dpointsize;
                        src2+=epointsize;
                        dest+=rpointsize;
                    }
                }
            }
            if (localerr)
            {
#pragma omp critical
                err=localerr;
            }
        }
    }
    else
    {
        throw DataException("applyBinaryCFunction: Unsupported combination of inputs.");
    }
#ifdef ESYS_MPI
    int global;
    MPI_Allreduce(&err, &global, 1, MPI_INT, MPI_MAX, get_MPIComm());
    err=global;
#endif
    if (err>0)
    {
        ostringstream oss;
        oss << "applyBinaryCFunction: error code " << err << " from C function.";
        throw DataException(oss.str());
    }
    return res;
}

#endif // IKNOWWHATIMDOING

Data
escript::condEval(escript::Data& mask, escript::Data& trueval, escript::Data& falseval)
{
    if (trueval.isComplex()!=falseval.isComplex())
    {
        trueval.complicate();
        falseval.complicate();
    }
    if (trueval.isComplex())
    {
        return condEvalWorker(mask, trueval, falseval, static_cast<DataTypes::cplx_t>(0));
    }
    else
    {
        return condEvalWorker(mask, trueval, falseval, static_cast<DataTypes::real_t>(0));
    }
}


template <typename S>
Data
escript::condEvalWorker(escript::Data& mask, escript::Data& trueval, escript::Data& falseval, S sentinel)
{
    // First, we need to make sure that trueval and falseval are compatible types.
    // also need to ensure that mask is a proper type for a mask
    // Need to choose a functionspace and shape for result
    // need to catch DataEmpty as well ?

    // only allowing scalar masks

    if (mask.getDataPointRank()!=0) 
    {
        throw DataException("Only supporting scalar masks");
        // Allowing people to slice two different objects together within a single datapoint - not allowing that

    }

    if (trueval.getDataPointShape()!=falseval.getDataPointShape()) 
    {
        throw DataException("condEval: shapes of true and false values must match.");
    }
    
    FunctionSpace fs=trueval.getFunctionSpace();        // should check this for compatibility as well
    if (trueval.getFunctionSpace()!=falseval.getFunctionSpace())
    {
        throw DataException("condEval: FunctionSpaces must match.");
    }
    // We aren't going to both with anything lazy except expanded data
    if (mask.isLazy() && !mask.actsExpanded())
    {
        mask.resolve();
    }
    if (trueval.isLazy() && !trueval.actsExpanded())
    {
        trueval.resolve();
    }
    if (falseval.isLazy() && !falseval.actsExpanded())
    {
        falseval.resolve();
    }

    if (mask.isConstant() && trueval.isConstant() && falseval.isConstant())
    {
        Data result(0,trueval.getDataPointShape(), fs , false);
        if (mask.getSampleDataRO(0,static_cast<DataTypes::real_t>(0))[0]>0)	// mask is always real
        {
            result.copy(trueval);
        }
        else
        {
            result.copy(falseval);
        }
        return result;
    }
    // Now we need to promote to correct ReadyData types
    // If they are lazy, they must be expanded
    if (mask.actsExpanded() || trueval.actsExpanded() || falseval.actsExpanded())
    {
        if (!mask.isLazy()) {mask.expand();}
        if (!trueval.isLazy()) {trueval.expand();}
        if (!falseval.isLazy()) {falseval.expand();}
    }
    else if (mask.isTagged() || trueval.isTagged() || falseval.isTagged())
    {
        mask.tag();
        trueval.tag();
        falseval.tag();
    }
    // by this point all data will be of the same ready type.
    if (mask.isTagged())
    {
        Data result(0,trueval.getDataPointShape(), fs , false);
        result.tag();
        DataTagged* rdat=dynamic_cast<DataTagged*>(result.getReady());
        const DataTagged* tdat=dynamic_cast<const DataTagged*>(trueval.getReady());
        const DataTagged* fdat=dynamic_cast<const DataTagged*>(falseval.getReady());
        const DataTagged* mdat=dynamic_cast<DataTagged*>(mask.getReady());
        //RealVectorType::const_pointer srcptr;
	const S* srcptr;

        // default value first
        if (mdat->getDefaultValueRO(0)>0)
        {
            srcptr=&(tdat->getDefaultValueRO(0, sentinel));
        } else {
            srcptr=&(fdat->getDefaultValueRO(0, sentinel));
        }
        for (int i=0;i<trueval.getDataPointSize();++i)
        {
            *(&(rdat->getDefaultValueRW(0, sentinel))+i)=*(srcptr+i);
        }

        // now we copy the tags from the mask - if the mask does not have it then it doesn't appear
        const DataTagged::DataMapType& maskLookup=mdat->getTagLookup();
        DataTagged::DataMapType::const_iterator it;
        DataTagged::DataMapType::const_iterator thisLookupEnd=maskLookup.end();
        for (it=maskLookup.begin();it!=thisLookupEnd;it++)
        {
            if (mdat->getDataByTagRO(it->first,0)>0)
            {
                rdat->addTaggedValue(it->first,trueval.getDataPointShape(), tdat->getVectorRO(), tdat->getOffsetForTag(it->first));
            }
            else
            {
                rdat->addTaggedValue(it->first,falseval.getDataPointShape(), fdat->getVectorRO(), fdat->getOffsetForTag(it->first));
            }
        }

        return result;
    }
    if (!trueval.actsExpanded() || !falseval.actsExpanded() || !mask.actsExpanded())
    {
        throw DataException("Programmer Error - Only actsExpanded Data should reach this point.");
    }
    else if (mask.actsExpanded() && trueval.actsExpanded() && falseval.actsExpanded())
    {
        // Here is the code for all expanded objects
        // this code will handle lazy data without expanding it just fine but lets allow people to lazify things

        if (mask.isLazy() || trueval.isLazy() || falseval.isLazy() || AUTOLAZYON)
        {
            DataAbstract_ptr pm=mask.borrowDataPtr();
            DataAbstract_ptr pt=trueval.borrowDataPtr();
            DataAbstract_ptr pf=falseval.borrowDataPtr();
            // now we create a lazy node for this
            DataLazy* p=new DataLazy(pm, pt, pf);
            return Data(p);
        }
        else
        {
            Data result(sentinel,trueval.getDataPointShape(), fs , true);  // Need to support non-expanded as well
            // OPENMP 3.0 allows unsigned loop vars.
#if defined(_OPENMP) && (_OPENMP < 200805)
            long i;
#else
            size_t i;
#endif
            auto& rvec=result.getReady()->getTypedVectorRW(sentinel);      // don't need to get acquireWrite since we made it
            unsigned int psize=result.getDataPointSize();
                
            size_t numsamples=result.getNumSamples();
            size_t dppsample=result.getNumDataPointsPerSample();
#pragma omp parallel for private(i) schedule(static)
            for (i=0;i<numsamples;++i)
            {
                // We are assuming that the first datapoint in the sample determines which side to use
                // for the whole sample.
                const decltype(sentinel)* src=0;
                const DataTypes::real_t* masksample=mask.getSampleDataRO(i, static_cast<DataTypes::real_t>(0));
                if (masksample[0]>0)    // first scalar determines whole sample
                {
                    src=trueval.getSampleDataRO(i, sentinel);
                }
                else
                {
                    src=falseval.getSampleDataRO(i, sentinel);
                }
                for (int j=0;j<dppsample;++j)
                {
                    size_t offset=j*psize;
                    for (long k=0;k<psize;++k)
                    {
                        rvec[i*dppsample*psize+offset+k]=(src)[offset+k];
                    }
                }
                
            }
            return result;
        }
    } else {
        throw DataException("condEval: Unsupported combination of DataAbstracts");
    }
}

DataTypes::RealVectorType& Data::getExpandedVectorReference(DataTypes::real_t dummy)
{
    if (!isExpanded())
    {
        expand();
    }
    return getReady()->getTypedVectorRW(dummy);
}

DataTypes::CplxVectorType& Data::getExpandedVectorReference(DataTypes::cplx_t dummy)
{
    if (!isExpanded())
    {
        expand();
    }
    return getReady()->getTypedVectorRW(dummy);
}

size_t Data::getNumberOfTaggedValues() const
{
    if (isTagged())
    {
        return m_data->getTagCount();
    }
    else
    {
        return 0;
    }
}


Data escript::randomData(const boost::python::tuple& shape,
       const FunctionSpace& what,
       long seed, const boost::python::tuple& filter)
{
  
    DataTypes::ShapeType dataPointShape;
    for (int i = 0; i < shape.attr("__len__")(); ++i) {
        dataPointShape.push_back(bp::extract<const int>(shape[i]));
    }  
  
  
    // first check what they have asked for in filter
    // does our domain support this?
    if (what.getDomain()->supportsFilter(filter))
    {
        return what.getDomain()->randomFill(dataPointShape, what, seed, filter);
    }
    else
    {
        throw DataException("The specified domain does not support those filter options.");
    }
}


namespace {
  
bp::object getNotImplemented()
{
    static bp::object notimpl=bp::object(bp::handle<>(
                bp::borrowed(PyImport_AddModule("__main__"))))
                        .attr("__builtins__").attr("NotImplemented");  
    return notimpl;
}
  
}

/* Implement part of pythons operator methods.
** We are doing this rather than just using boost's overloading shortcuts because we
** want to be able to return NotImplemented.
*/
bp::object Data::__add__(const bp::object& right)
{
    bp::extract<Data> data_extractor(right);
    if (data_extractor.check())   // if this is wrapping a Data
    {
        return bp::object(*this+data_extractor());
    }
    bool wrapok=false;  // If the object can't be wrapped we should return NotImplemented
    try                 // if the exception is due to something else we should rethrow it
    {
        WrappedArray w(right);
        wrapok=true;
        return bp::object(*this+Data(w, this->getFunctionSpace(), false));    
    }
    catch (DataException& e)
    {
        if (wrapok)
        {
            throw e;
        }
        return getNotImplemented();
    }
}

bp::object Data::__sub__(const bp::object& right)
{
    bp::extract<Data> data_extractor(right);
    if (data_extractor.check())   // if this is wrapping a Data
    {
        return bp::object(*this-data_extractor());
    }
    bool wrapok=false;  // If the object can't be wrapped we should return NotImplemented
    try                 // if the exception is due to something else we should rethrow it
    {
        WrappedArray w(right);
        wrapok=true;
        return bp::object(*this-Data(w, this->getFunctionSpace(), false));  
    }
    catch (DataException& e)
    {
        if (wrapok)
        {
            throw e;
        }
        return getNotImplemented();
    }    
}

bp::object Data::__rsub__(const bp::object& right)
{
    bp::extract<Data> data_extractor(right);
    if (data_extractor.check())   // if this is wrapping a Data
    {
        return bp::object(data_extractor()-*this);
    }
    bool wrapok=false;  // If the object can't be wrapped we should return NotImplemented
    try                 // if the exception is due to something else we should rethrow it
    {
        WrappedArray w(right);
        wrapok=true;
        return bp::object(Data(w, this->getFunctionSpace(),false)-*this); 
    }
    catch (DataException& e)
    {
        if (wrapok)
        {
            throw e;
        }
        return getNotImplemented();
    }      
 
}

bp::object Data::__mul__(const bp::object& right)
{
    bp::extract<Data> data_extractor(right);
    if (data_extractor.check())   // if this is wrapping a Data
    {
        return bp::object(*this*data_extractor());
    }
    bool wrapok=false;  // If the object can't be wrapped we should return NotImplemented
    try                 // if the exception is due to something else we should rethrow it
    {
        WrappedArray w(right);
        wrapok=true;
        return bp::object(*this*Data(w, this->getFunctionSpace(),false));  
    }
    catch (DataException& e)
    {
        if (wrapok)
        {
            throw e;
        }
        return getNotImplemented();
    }       
}

bp::object Data::__div__(const bp::object& right)
{
    bp::extract<Data> data_extractor(right);
    if (data_extractor.check())   // if this is wrapping a Data
    {
        return bp::object(*this/data_extractor());
    }
    bool wrapok=false;  // If the object can't be wrapped we should return NotImplemented
    try                 // if the exception is due to something else we should rethrow it
    {
        WrappedArray w(right);
        wrapok=true;
        return bp::object(*this/Data(w, this->getFunctionSpace(),false));  
    }
    catch (DataException& e)
    {
        if (wrapok)
        {
            throw e;
        }
        return getNotImplemented();
    }     
}

bp::object Data::__rdiv__(const bp::object& right)
{
    bp::extract<Data> data_extractor(right);
    if (data_extractor.check())   // if this is wrapping a Data
    {
        return bp::object(data_extractor()/(*this));
    }
    bool wrapok=false;  // If the object can't be wrapped we should return NotImplemented
    try                 // if the exception is due to something else we should rethrow it
    {
        WrappedArray w(right);
        wrapok=true;
        return bp::object(Data(w, this->getFunctionSpace(),false)/(*this));  
    }
    catch (DataException& e)
    {
        if (wrapok)
        {
            throw e;
        }
        return getNotImplemented();
    }         
}

void Data::complicate()
{
    if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
    }  
    
    if (m_data->isLazy())
    {
            // This is different to the other types because instead of switching from 
            // one internal storage vector to the other within the same node
            // m_data needs to be replaced with an new root (promote) node.
        DataLazy_ptr nn=dynamic_pointer_cast<DataLazy>(m_data);
        DataLazy_ptr res=makePromote(nn);
        set_m_data(res);
    }
    else
    {
        m_data->complicate();
    }
}

Data
escript::C_TensorUnaryOperation(Data const &arg_0,
                       escript::ES_optype operation,
                       DataTypes::real_t tol)
{
  if (arg_0.isEmpty())  // do this before we attempt to interpolate
  {
     throw DataException("Error - Operations (C_TensorUnaryOperation) not permitted on instances of DataEmpty.");
  }
  if (arg_0.isLazy())
  {
     throw DataException("Error - Operations not permitted on lazy data.");
  }
  
  if (arg_0.isComplex() && !supports_cplx(operation))
  {
      throw DataException("Error - the requested operation does not support complex values");
  }
  
  // Interpolate if necessary and find an appropriate function space
  Data arg_0_Z = Data(arg_0);

  // Get rank and shape of inputs
  const DataTypes::ShapeType& shape0 = arg_0_Z.getDataPointShape();
  int size0 = arg_0_Z.getDataPointSize();
  
  // Declare output Data object
  Data res;
  bool emptyResult=(arg_0_Z.getNumSamples()==0);
  if (arg_0_Z.isConstant()) {
    if (arg_0_Z.isComplex())                    // this is not taking into account cplx->real
    {
        DataTypes::cplx_t dummy=0;
        res = Data(0.0, shape0, arg_0_Z.getFunctionSpace(),0);      // DataConstant output
        const DataTypes::cplx_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(0, dummy));
        if (always_real(operation))
        {
	    if (emptyResult)
	    {
		return res;
	    }
            DataTypes::real_t *ptr_2 = &(res.getDataAtOffsetRW(0, (real_t)(0)));
            tensor_unary_array_operation_real(size0, ptr_0, ptr_2, operation, tol);       
        }
        else
        {
            res.complicate();
	    if (emptyResult)
	    {
		return res;
	    }	    
            DataTypes::cplx_t *ptr_2 = &(res.getDataAtOffsetRW(0, dummy));
            tensor_unary_array_operation(size0, ptr_0, ptr_2, operation, tol);
        }
    }
    else
    {
        // This currently does not call the tensor_unary_array_operation_real
        // functions like .real() and .imag() but they are caught in the Data interface
        DataTypes::real_t dummy=0;
        res = Data(0.0, shape0, arg_0_Z.getFunctionSpace(),false);      // DataConstant output
	if (emptyResult)
	{
	    return res;
	}
	
        const DataTypes::real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(0, dummy));
        DataTypes::real_t *ptr_2 = &(res.getDataAtOffsetRW(0, dummy));
        if (always_real(operation))
        {
            tensor_unary_array_operation_real(size0, ptr_0, ptr_2, operation, tol);
        }
        else
        {
            tensor_unary_array_operation(size0, ptr_0, ptr_2, operation, tol);
        }
    }
  }
  else if (arg_0_Z.isTagged()) {

    // Borrow DataTagged input from Data object
    DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());

    // Prepare a DataTagged output 2
    res = Data(0.0, shape0, arg_0_Z.getFunctionSpace(),false);   // DataTagged output


    if (arg_0_Z.isComplex())
    {
        if (always_real(operation))
        {
            res.tag();
	    if (emptyResult)
	    {
		return res;
	    }
	    
            DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());      
          
            DataTypes::cplx_t dummy=0;
            // Get the pointers to the actual data
            const DataTypes::cplx_t *ptr_0 = &(tmp_0->getDefaultValueRO(0,dummy));
            DataTypes::real_t *ptr_2 = &(tmp_2->getDefaultValueRW(0,real_t(0)));
            // Compute a result for the default
            tensor_unary_array_operation_real(size0, ptr_0, ptr_2, operation, tol);
            // Compute a result for each tag
            const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
            DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
            for (i=lookup_0.begin();i!=lookup_0.end();i++) {
              tmp_2->addTag(i->first);
              const DataTypes::cplx_t *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0, dummy));
              DataTypes::real_t *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0, real_t(0)));
              tensor_unary_array_operation_real(size0, ptr_0, ptr_2, operation, tol);
            }
        }
        else
        {
            res.complicate();
            res.tag();
	    if (emptyResult)
	    {
		return res;
	    }
	    
            DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());      
          
            DataTypes::cplx_t dummy=0;
            // Get the pointers to the actual data
            const DataTypes::cplx_t *ptr_0 = &(tmp_0->getDefaultValueRO(0,dummy));
            DataTypes::cplx_t *ptr_2 = &(tmp_2->getDefaultValueRW(0,dummy));
            // Compute a result for the default
            tensor_unary_array_operation(size0, ptr_0, ptr_2, operation, tol);
            // Compute a result for each tag
            const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
            DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
            for (i=lookup_0.begin();i!=lookup_0.end();i++) {
              tmp_2->addTag(i->first);
              const DataTypes::cplx_t *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0, dummy));
              DataTypes::cplx_t *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0, dummy));
              tensor_unary_array_operation(size0, ptr_0, ptr_2, operation, tol);
            }
        }
    }
    else
    {
      
        res.tag();
	if (emptyResult)
	{
	    return res;
	}
	
        DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());      
      
        // Get the pointers to the actual data
        const DataTypes::real_t *ptr_0 = &(tmp_0->getDefaultValueRO(0));
        DataTypes::real_t *ptr_2 = &(tmp_2->getDefaultValueRW(0));
        // Compute a result for the default
        if (always_real(operation))
        {
            tensor_unary_array_operation_real(size0, ptr_0, ptr_2, operation, tol);       
        }
        else
        {
            tensor_unary_array_operation(size0, ptr_0, ptr_2, operation, tol);
        }
        // Compute a result for each tag
        const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
        DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
        for (i=lookup_0.begin();i!=lookup_0.end();i++) {
          tmp_2->addTag(i->first);
          const DataTypes::real_t *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0));
          DataTypes::real_t *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0));
          if (always_real(operation))
          {
              tensor_unary_array_operation_real(size0, ptr_0, ptr_2, operation, tol);
          }
          else
          {
              tensor_unary_array_operation(size0, ptr_0, ptr_2, operation, tol);
          }
        }
    }
  }
  else if (arg_0_Z.isExpanded()) 
  {

    res = Data(0.0, shape0, arg_0_Z.getFunctionSpace(),true); // DataExpanded output
    if (arg_0_Z.isComplex() && !always_real(operation))
    {
        res.complicate();
    }
    if (emptyResult)
    {
	return res;
    }
    
    DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
    DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

    int sampleNo_0,dataPointNo_0;
    int numSamples_0 = arg_0_Z.getNumSamples();
    int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
    if (arg_0_Z.isComplex())
    {
        if (always_real(operation))
        {
            DataTypes::cplx_t dummy=0;
            #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
            for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
                dataPointNo_0=0;
                int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
                int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
                const DataTypes::cplx_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy));
                DataTypes::real_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, real_t(0)));
                tensor_unary_array_operation_real(size0*numDataPointsPerSample_0, ptr_0, ptr_2, operation, tol);
            }             
        }
        else
        {
            DataTypes::cplx_t dummy=0;
            #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
            for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
                dataPointNo_0=0;
                int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
                int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
                const DataTypes::cplx_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy));
                DataTypes::cplx_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummy));
                tensor_unary_array_operation(size0*numDataPointsPerSample_0, ptr_0, ptr_2, operation, tol);
            }     
        }
    }
    else
    {
        // we require storage to be contiguous so let's do it in one chunk
        #pragma omp parallel private(sampleNo_0,dataPointNo_0)
        {
#ifdef _OPENMP
            int tid=omp_get_thread_num();
            int mt=omp_get_num_threads();
            int rem=numSamples_0%mt;
            size_t samples_per=numSamples_0/mt;
            size_t startsample=samples_per*tid+((tid<rem)?tid:rem);
            size_t nextsample=samples_per*(tid+1)+(((tid+1)<rem)?(tid+1):rem);
            size_t ulimit=min<size_t>(nextsample, numSamples_0);
            size_t samples=ulimit-startsample;    
#else
            size_t startsample=0;
            size_t samples=numSamples_0;
#endif      
            if (startsample<numSamples_0)
            {
	        real_t dummyr=0;
                size_t offset_0 = tmp_0->getPointOffset(startsample,0);
                size_t offset_2 = tmp_2->getPointOffset(startsample,0);
                const DataTypes::real_t *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummyr));
                DataTypes::real_t *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
                if (always_real(operation))
                {
                    tensor_unary_array_operation_real(size0*samples*numDataPointsPerSample_0, ptr_0, ptr_2, operation, tol);
                }
                else
                {
                    tensor_unary_array_operation(size0*samples*numDataPointsPerSample_0, ptr_0, ptr_2, operation, tol);
                }
            }
        }           
    }
  }
  else {
    throw DataException("Error - C_TensorUnaryOperation: unknown combination of inputs");
  }

  return res;
}

Data
escript::C_TensorBinaryOperation(Data const &arg_0,
                        Data const &arg_1,
                        escript::ES_optype operation)
{
  if (arg_0.isEmpty() || arg_1.isEmpty())
  {
     throw DataException("Error - Operations (C_TensorBinaryOperation) not permitted on instances of DataEmpty.");
  }
  if (arg_0.isLazy() || arg_1.isLazy())
  {
     throw DataException("Error - Operations not permitted on lazy data.");
  }
  
  // Interpolate if necessary and find an appropriate function space
  Data arg_0_Z, arg_1_Z;
  FunctionSpace fsl=arg_0.getFunctionSpace();
  FunctionSpace fsr=arg_1.getFunctionSpace();
  if (fsl!=fsr) {
     signed char intres=fsl.getDomain()->preferredInterpolationOnDomain(fsr.getTypeCode(), fsl.getTypeCode());
     if (intres==0)
     {
         std::string msg="Error - C_TensorBinaryOperation: arguments have incompatible function spaces.";
         msg+=fsl.toString();
         msg+=" ";
         msg+=fsr.toString();
         throw DataException(msg.c_str());
     } 
     else if (intres==1)
     {
      arg_1_Z=arg_1.interpolate(arg_0.getFunctionSpace());
      arg_0_Z =Data(arg_0);      
     }
     else	// reverse interpolation preferred
     {
      arg_0_Z = arg_0.interpolate(arg_1.getFunctionSpace());
      arg_1_Z = Data(arg_1);
     }    
  } else {
      arg_0_Z = Data(arg_0);
      arg_1_Z = Data(arg_1);
  }
  DataTypes::ShapeType shape0 = arg_0_Z.getDataPointShape();
  DataTypes::ShapeType shape1 = arg_1_Z.getDataPointShape();
  
  DataTypes::ShapeType resultshape=((arg_0_Z.getDataPointRank()!=0)?shape0:shape1);

  bool emptyResult=((arg_0_Z.getNumSamples()==0) || (arg_1_Z.getNumSamples()==0));
  if ((shape0==shape1) || (arg_0_Z.getDataPointRank()==0) || (arg_1_Z.getDataPointRank()==0))
  {
    if (arg_0_Z.isConstant()   && arg_1_Z.isConstant())
    {
      Data res(0.0, resultshape, arg_1_Z.getFunctionSpace(),false);      // DataConstant output
      if (arg_0_Z.isComplex() || arg_1_Z.isComplex())
      {
        res.complicate();
      }
      if (!emptyResult)
      {      
          binaryOpDataCCC(*dynamic_cast<DataConstant*>(res.borrowData()), *dynamic_cast<const DataConstant*>(arg_0_Z.borrowData()), *dynamic_cast<const DataConstant*>(arg_1_Z.borrowData()), operation);
      }
      return res;
    }
    else if (arg_0_Z.isConstant()   && arg_1_Z.isTagged())
    {
      Data res(0.0, resultshape, arg_1_Z.getFunctionSpace(), false);      // DataTagged output
      if (arg_0_Z.isComplex() || arg_1_Z.isComplex())
      {
        res.complicate();
      }
      res.tag();
      if (!emptyResult)
      {
          binaryOpDataTCT(*dynamic_cast<DataTagged*>(res.borrowData()), *dynamic_cast<const DataConstant*>(arg_0_Z.borrowData()), *dynamic_cast<const DataTagged*>(arg_1_Z.borrowData()), operation);
      }
      return res;
    }
    else if (arg_0_Z.isConstant()   && arg_1_Z.isExpanded())
    {
      Data res(0.0, resultshape, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      if (arg_0_Z.isComplex() || arg_1_Z.isComplex())
      {
        res.complicate();
      }
      if (!emptyResult)
      {
          binaryOpDataECE(*dynamic_cast<DataExpanded*>(res.borrowData()), *dynamic_cast<const DataConstant*>(arg_0_Z.borrowData()), *dynamic_cast<const DataExpanded*>(arg_1_Z.borrowData()), operation);
      }
      return res;
    }
    else if (arg_0_Z.isTagged()     && arg_1_Z.isConstant())
    {
      Data res(0.0, resultshape, arg_0_Z.getFunctionSpace(),false);      // DataTagged output
      if (arg_0_Z.isComplex() || arg_1_Z.isComplex())
      {
        res.complicate();
      }
      res.tag();
      if (!emptyResult)
      {
          binaryOpDataTTC(*dynamic_cast<DataTagged*>(res.borrowData()), *dynamic_cast<const DataTagged*>(arg_0_Z.borrowData()), *dynamic_cast<const DataConstant*>(arg_1_Z.borrowData()), operation);
      }
      return res;
    }
    else if (arg_0_Z.isTagged()     && arg_1_Z.isTagged())
    {
      Data res(0.0, resultshape, arg_1_Z.getFunctionSpace(), false);
      if (arg_0_Z.isComplex() || arg_1_Z.isComplex())
      {
        res.complicate();
      }
      res.tag();        // DataTagged output
      if (!emptyResult)
      {
          binaryOpDataTTT(*dynamic_cast<DataTagged*>(res.borrowData()), *dynamic_cast<const DataTagged*>(arg_0_Z.borrowData()), *dynamic_cast<const DataTagged*>(arg_1_Z.borrowData()), operation);
      }
      return res;
    }
    else if (arg_0_Z.isTagged()     && arg_1_Z.isExpanded())
    {
      Data res(0.0, resultshape, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      if (arg_0_Z.isComplex() || arg_1_Z.isComplex())
      {
        res.complicate();
      }
      if (!emptyResult)
      {
          binaryOpDataETE(*dynamic_cast<DataExpanded*>(res.borrowData()), *dynamic_cast<const DataTagged*>(arg_0_Z.borrowData()), *dynamic_cast<const DataExpanded*>(arg_1_Z.borrowData()), operation);
      }
      return res;
    }
    else if (arg_0_Z.isExpanded()   && arg_1_Z.isConstant()) {
      Data res(0.0, resultshape, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      if (arg_0_Z.isComplex() || arg_1_Z.isComplex())
      {
        res.complicate();
      }
      if (!emptyResult)
      {
          binaryOpDataEEC(*dynamic_cast<DataExpanded*>(res.borrowData()), *dynamic_cast<const DataExpanded*>(arg_0_Z.borrowData()), *dynamic_cast<const DataConstant*>(arg_1_Z.borrowData()), operation);
      }
      return res;
    }
    else if (arg_0_Z.isExpanded()   && arg_1_Z.isTagged()) {
      Data res(0.0, resultshape, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      if (arg_0_Z.isComplex() || arg_1_Z.isComplex())
      {
        res.complicate();
      }
      if (!emptyResult)
      {
          binaryOpDataEET(*dynamic_cast<DataExpanded*>(res.borrowData()), *dynamic_cast<const DataExpanded*>(arg_0_Z.borrowData()), *dynamic_cast<const DataTagged*>(arg_1_Z.borrowData()), operation);
      }
      return res;
    }
    else if (arg_0_Z.isExpanded()   && arg_1_Z.isExpanded()) {
      Data res(0.0, resultshape, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      if (arg_0_Z.isComplex() || arg_1_Z.isComplex())
      {
        res.complicate();
      }
      if (!emptyResult)
      {
          binaryOpDataEEE(*dynamic_cast<DataExpanded*>(res.borrowData()), *dynamic_cast<const DataExpanded*>(arg_0_Z.borrowData()), *dynamic_cast<const DataExpanded*>(arg_1_Z.borrowData()), operation);
      }
      return res;
    }
    else {
      throw DataException("Error - C_TensorBinaryOperation: unknown combination of inputs");
    }
  } else {
    throw DataException("Error - C_TensorBinaryOperation: arguments have incompatible shapes");
  }
}


void
Data::TensorSelfUpdateBinaryOperation(const Data& right,
                   escript::ES_optype operation)
{
   //
   // if this has a rank of zero promote it to the rank of the RHS
   if (getDataPointRank()==0 && right.getDataPointRank()!=0) {
     throw DataException("Error - attempt to update rank zero object with object with rank bigger than zero.");
   }

   if (isLazy() || right.isLazy())
   {
     throw DataException("Programmer error - attempt to call binaryOp with Lazy Data.");
   }
   //
   // initially make the temporary a shallow copy
   Data tempRight(right);
   FunctionSpace fsl=getFunctionSpace();
   FunctionSpace fsr=right.getFunctionSpace();
   if (fsl!=fsr) {
     signed char intres=fsl.getDomain()->preferredInterpolationOnDomain(fsr.getTypeCode(), fsl.getTypeCode());
     if (intres==0)
     {
         std::string msg="Error - attempt to combine incompatible FunctionSpaces.";
         msg+=fsl.toString();
         msg+="  ";
         msg+=fsr.toString();
         throw DataException(msg.c_str());
     } 
     else if (intres==1)
     {
       // an interpolation is required so create a new Data
       tempRight=Data(right,fsl);
     }
     else       // reverse interpolation preferred
     {
        // interpolate onto the RHS function space
       Data tempLeft(*this,fsr);
       set_m_data(tempLeft.m_data);
     }
   }
   operandCheck(tempRight);
   //
   // ensure this has the right type for the RHS
   typeMatchRight(tempRight);
   //
   // Need to cast to the concrete types so that the correct binaryOp
   // is called.
   if (isExpanded()) {
     //
     // Expanded data will be done in parallel, the right hand side can be
     // of any data type
     DataExpanded* leftC=dynamic_cast<DataExpanded*>(m_data.get());
     ESYS_ASSERT(leftC!=NULL, "Programming error - casting to DataExpanded.");
     
     if (right.isExpanded())
     {
	binaryOpDataEEE(*leftC, *leftC, *dynamic_cast<const DataExpanded*>(tempRight.getReady()), operation);
     }
     else if (right.isTagged())
     {
	binaryOpDataEET(*leftC, *leftC, *dynamic_cast<const DataTagged*>(tempRight.getReady()), operation);
     }
     else	// it's constant
     {
	binaryOpDataEEC(*leftC, *leftC, *dynamic_cast<const DataConstant*>(tempRight.getReady()), operation);
     }
       
     //escript::binaryOpDataReady(*leftC,*(tempRight.getReady()),operation);
   } else if (isTagged()) {
     //
     // Tagged data is operated on serially, the right hand side can be
     // either DataConstant or DataTagged
     DataTagged* leftC=dynamic_cast<DataTagged*>(m_data.get());
     ESYS_ASSERT(leftC!=NULL, "Programming error - casting to DataTagged.");
     if (right.isTagged()) {
       DataTagged* rightC=dynamic_cast<DataTagged*>(tempRight.m_data.get());
       ESYS_ASSERT(rightC!=NULL, "Programming error - casting to DataTagged.");
       binaryOpDataTTT(*leftC, *leftC, *rightC, operation);
       //escript::binaryOpDataReady(*leftC,*rightC,operation);
     } else {
       DataConstant* rightC=dynamic_cast<DataConstant*>(tempRight.m_data.get());
       ESYS_ASSERT(rightC!=NULL, "Programming error - casting to DataConstant.");
       binaryOpDataTTC(*leftC, *leftC, *rightC, operation);
       //escript::binaryOpDataReady(*leftC,*rightC,operation);
     }
   } else if (isConstant()) {
     DataConstant* leftC=dynamic_cast<DataConstant*>(m_data.get());
     DataConstant* rightC=dynamic_cast<DataConstant*>(tempRight.m_data.get());
     ESYS_ASSERT(leftC!=NULL && rightC!=NULL, "Programming error - casting to DataConstant.");
     binaryOpDataCCC(*leftC, *leftC, *rightC, operation);
     //escript::binaryOpDataReady(*leftC,*rightC,operation);
   }  
}

#if 0
void
Data::binaryDataOp(const Data& right,
                   escript::ES_optype operation)
{
   //
   // if this has a rank of zero promote it to the rank of the RHS
   if (getDataPointRank()==0 && right.getDataPointRank()!=0) {
     throw DataException("Error - attempt to update rank zero object with object with rank bigger than zero.");
   }

   if (isLazy() || right.isLazy())
   {
     throw DataException("Programmer error - attempt to call binaryOp with Lazy Data.");
   }
   //
   // initially make the temporary a shallow copy
   Data tempRight(right);
   FunctionSpace fsl=getFunctionSpace();
   FunctionSpace fsr=right.getFunctionSpace();
   if (fsl!=fsr) {
     signed char intres=fsl.getDomain()->preferredInterpolationOnDomain(fsr.getTypeCode(), fsl.getTypeCode());
     if (intres==0)
     {
         std::string msg="Error - attempt to combine incompatible FunctionSpaces.";
         msg+=fsl.toString();
         msg+="  ";
         msg+=fsr.toString();
         throw DataException(msg.c_str());
     } 
     else if (intres==1)
     {
       // an interpolation is required so create a new Data
       tempRight=Data(right,fsl);
     }
     else       // reverse interpolation preferred
     {
        // interpolate onto the RHS function space
       Data tempLeft(*this,fsr);
       set_m_data(tempLeft.m_data);
     }
   }
   operandCheck(tempRight);
   //
   // ensure this has the right type for the RHS
   typeMatchRight(tempRight);
   //
   // Need to cast to the concrete types so that the correct binaryOp
   // is called.
   if (isExpanded()) {
     //
     // Expanded data will be done in parallel, the right hand side can be
     // of any data type
     DataExpanded* leftC=dynamic_cast<DataExpanded*>(m_data.get());
     ESYS_ASSERT(leftC!=NULL, "Programming error - casting to DataExpanded.");
     escript::binaryOpDataReady(*leftC,*(tempRight.getReady()),operation);
   } else if (isTagged()) {
     //
     // Tagged data is operated on serially, the right hand side can be
     // either DataConstant or DataTagged
     DataTagged* leftC=dynamic_cast<DataTagged*>(m_data.get());
     ESYS_ASSERT(leftC!=NULL, "Programming error - casting to DataTagged.");
     if (right.isTagged()) {
       DataTagged* rightC=dynamic_cast<DataTagged*>(tempRight.m_data.get());
       ESYS_ASSERT(rightC!=NULL, "Programming error - casting to DataTagged.");
       escript::binaryOpDataReady(*leftC,*rightC,operation);
     } else {
       DataConstant* rightC=dynamic_cast<DataConstant*>(tempRight.m_data.get());
       ESYS_ASSERT(rightC!=NULL, "Programming error - casting to DataConstant.");
       escript::binaryOpDataReady(*leftC,*rightC,operation);
     }
   } else if (isConstant()) {
     DataConstant* leftC=dynamic_cast<DataConstant*>(m_data.get());
     DataConstant* rightC=dynamic_cast<DataConstant*>(tempRight.m_data.get());
     ESYS_ASSERT(leftC!=NULL && rightC!=NULL,
             "Programming error - casting to DataConstant.");
     escript::binaryOpDataReady(*leftC,*rightC,operation);
   }  
}

#endif

