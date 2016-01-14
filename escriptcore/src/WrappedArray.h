
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/** \file WrappedArray.h */

#ifndef WrappedArray_20081202_H
#define WrappedArray_20081202_H
#include "system_dep.h"
#include "DataTypes.h"
#include "boost/python/extract.hpp"
#include <complex>

namespace escript
{

class WrappedArray
{
public:
	typedef  std::complex<double> complextype;
	WrappedArray(const boost::python::object& obj_in);
	~WrappedArray();
	unsigned int getRank() const;
	const DataTypes::ShapeType& getShape() const;
	bool isComplex() const;
	double getElt() const;
	double getElt(unsigned int i) const;
	double getElt(unsigned int i, unsigned int j) const;
	double getElt(unsigned int i, unsigned int j, unsigned int k) const;
	double getElt(unsigned int i, unsigned int j, unsigned int k, unsigned int m) const;
	
	complextype getEltC() const;
	complextype getEltC(unsigned int i) const;
	complextype getEltC(unsigned int i, unsigned int j) const;
	complextype getEltC(unsigned int i, unsigned int j, unsigned int k) const;
	complextype getEltC(unsigned int i, unsigned int j, unsigned int k, unsigned int m) const;
	
	
	void convertArray() const;
private:
	void convertArrayR() const;
	void convertArrayC() const;
	template<typename T> void convertNumpyArray(const T* array, const std::vector<int>& strides) const;
	template<typename T> void convertNumpyArrayC(const T* array, const std::vector<int>& strides) const;
	const boost::python::object& obj;
	int rank;	
	mutable bool converted;		// has the array been converted to a C array
	bool iscomplex;		// is the wrapped array storing complex values?	
	escript::DataTypes::ShapeType shape;
	double scalar_r;
	complextype scalar_c;
	mutable double* dat_r;			// real data
	mutable complextype* dat_c;	// complex data   - only one of these members should be used
};


inline bool WrappedArray::isComplex() const
{
    return iscomplex;
}

inline unsigned int 
WrappedArray::getRank() const
{
	return rank;
}

inline const DataTypes::ShapeType& 
WrappedArray::getShape() const
{
	return shape;
}

inline double
WrappedArray::getElt() const
{
	return scalar_r;
}


inline double
WrappedArray::getElt(unsigned int i) const
{  // __float__ added to deal with numpy. If this causes problems we may have to register a custom converter
    if (iscomplex)
    {
      return nan("");
    }
    return (dat_r!=0)?dat_r[i]:(boost::python::extract<double>(obj[i].attr("__float__")()));	
}

inline
double 
WrappedArray::getElt(unsigned int i, unsigned int j) const
{
    if (iscomplex)
    {
      return nan("");
    }  
    return (dat_r!=0)?dat_r[DataTypes::getRelIndex(shape,i,j)]:(boost::python::extract<double>(obj[i][j].attr("__float__")()));
}

inline
double 
WrappedArray::getElt(unsigned int i, unsigned int j, unsigned int k) const
{
    if (iscomplex)
    {
      return nan("");
    }    
    return (dat_r!=0)?dat_r[DataTypes::getRelIndex(shape,i,j,k)]:(boost::python::extract<double>(obj[i][j][k].attr("__float__")()));
}

inline
double 
WrappedArray::getElt(unsigned int i, unsigned int j, unsigned int k, unsigned int m) const
{
    if (iscomplex)
    {
      return nan("");
    }  
    return (dat_r!=0)?dat_r[DataTypes::getRelIndex(shape,i,j,k,m)]:(boost::python::extract<double>(obj[i][j][k][m].attr("__float__")()));
}





inline WrappedArray::complextype
WrappedArray::getEltC() const
{
	return scalar_c;
}


inline WrappedArray::complextype
WrappedArray::getEltC(unsigned int i) const
{
    if (!iscomplex)
    {
      return nan("");
    }
    return (dat_c!=0)?dat_c[i]:(boost::python::extract<complextype>(obj[i]));	// don't know if this will work with numpy	
}

inline
WrappedArray::complextype 
WrappedArray::getEltC(unsigned int i, unsigned int j) const
{
    if (!iscomplex)
    {
      return nan("");
    }  
    return (dat_c!=0)?dat_c[DataTypes::getRelIndex(shape,i,j)]:(boost::python::extract<complextype>(obj[i][j]));
}

inline
WrappedArray::complextype 
WrappedArray::getEltC(unsigned int i, unsigned int j, unsigned int k) const
{
    if (!iscomplex)
    {
      return nan("");
    }    
    return (dat_c!=0)?dat_c[DataTypes::getRelIndex(shape,i,j,k)]:(boost::python::extract<complextype>(obj[i][j][k]));
}

inline
WrappedArray::complextype 
WrappedArray::getEltC(unsigned int i, unsigned int j, unsigned int k, unsigned int m) const
{
    if (!iscomplex)
    {
      return nan("");
    }  
    return (dat_c!=0)?dat_c[DataTypes::getRelIndex(shape,i,j,k,m)]:(boost::python::extract<complextype>(obj[i][j][k][m]));
}





}

#endif

