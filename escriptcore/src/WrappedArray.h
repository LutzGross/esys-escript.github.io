
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


/** \file WrappedArray.h */

#ifndef WrappedArray_20081202_H
#define WrappedArray_20081202_H
#include "system_dep.h"
#include "DataTypes.h"
#include "boost/python/extract.hpp"
#include "boost/python/object.hpp"
#include <complex>

namespace escript
{

class ESCRIPT_DLL_API WrappedArray
{
public:
	WrappedArray(const boost::python::object& obj_in);
	~WrappedArray();
	unsigned int getRank() const;
	const DataTypes::ShapeType& getShape() const;
	bool isComplex() const;
	DataTypes::real_t getElt() const;
	DataTypes::real_t getElt(unsigned int i) const;
	DataTypes::real_t getElt(unsigned int i, unsigned int j) const;
	DataTypes::real_t getElt(unsigned int i, unsigned int j, unsigned int k) const;
	DataTypes::real_t getElt(unsigned int i, unsigned int j, unsigned int k, unsigned int m) const;
	
	DataTypes::cplx_t getEltC() const;
	DataTypes::cplx_t getEltC(unsigned int i) const;
	DataTypes::cplx_t getEltC(unsigned int i, unsigned int j) const;
	DataTypes::cplx_t getEltC(unsigned int i, unsigned int j, unsigned int k) const;
	DataTypes::cplx_t getEltC(unsigned int i, unsigned int j, unsigned int k, unsigned int m) const;
	
	
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
	DataTypes::real_t scalar_r;
	DataTypes::cplx_t scalar_c;
	mutable DataTypes::real_t* dat_r;			// real data
	mutable DataTypes::cplx_t* dat_c;	// complex data   - only one of these members should be used
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

inline DataTypes::real_t
WrappedArray::getElt() const
{
    if (iscomplex)
    {
      return nan("");
    }  
    return scalar_r;
}


inline DataTypes::real_t
WrappedArray::getElt(unsigned int i) const
{  // __float__ added to deal with numpy. If this causes problems we may have to register a custom converter
    if (iscomplex)
    {
      return nan("");
    }
    return (dat_r!=0)?dat_r[i]:(boost::python::extract<DataTypes::real_t>(obj[i].attr("__float__")()));	
}

inline
DataTypes::real_t 
WrappedArray::getElt(unsigned int i, unsigned int j) const
{
    if (iscomplex)
    {
      return nan("");
    }  
    return (dat_r!=0)?dat_r[DataTypes::getRelIndex(shape,i,j)]:(boost::python::extract<DataTypes::real_t>(obj[i][j].attr("__float__")()));
}

inline
DataTypes::real_t 
WrappedArray::getElt(unsigned int i, unsigned int j, unsigned int k) const
{
    if (iscomplex)
    {
      return nan("");
    }    
    return (dat_r!=0)?dat_r[DataTypes::getRelIndex(shape,i,j,k)]:(boost::python::extract<DataTypes::real_t>(obj[i][j][k].attr("__float__")()));
}

inline
DataTypes::real_t 
WrappedArray::getElt(unsigned int i, unsigned int j, unsigned int k, unsigned int m) const
{
    if (iscomplex)
    {
      return nan("");
    }  
    return (dat_r!=0)?dat_r[DataTypes::getRelIndex(shape,i,j,k,m)]:(boost::python::extract<DataTypes::real_t>(obj[i][j][k][m].attr("__float__")()));
}





inline DataTypes::cplx_t
WrappedArray::getEltC() const
{
    if (!iscomplex)
    {
      return scalar_r;
    }  
    return scalar_c;
}


inline DataTypes::cplx_t
WrappedArray::getEltC(unsigned int i) const
{
    if (!iscomplex)	// let's try to get a real value out instead
    {
      return (dat_r!=0)?dat_r[i]:(boost::python::extract<DataTypes::real_t>(obj[i]));
    }
    return (dat_c!=0)?dat_c[i]:(boost::python::extract<DataTypes::cplx_t>(obj[i]));	// don't know if this will work with numpy	
}

inline
DataTypes::cplx_t 
WrappedArray::getEltC(unsigned int i, unsigned int j) const
{
    if (!iscomplex)
    {
       return (dat_r!=0)?dat_r[DataTypes::getRelIndex(shape,i,j)]:(boost::python::extract<DataTypes::real_t>(obj[i][j]));
    }  
    return (dat_c!=0)?dat_c[DataTypes::getRelIndex(shape,i,j)]:(boost::python::extract<DataTypes::cplx_t>(obj[i][j]));
}

inline
DataTypes::cplx_t 
WrappedArray::getEltC(unsigned int i, unsigned int j, unsigned int k) const
{
    if (!iscomplex)
    {
      return (dat_r!=0)?dat_r[DataTypes::getRelIndex(shape,i,j,k)]:(boost::python::extract<DataTypes::real_t>(obj[i][j][k]));
    }    
    return (dat_c!=0)?dat_c[DataTypes::getRelIndex(shape,i,j,k)]:(boost::python::extract<DataTypes::cplx_t>(obj[i][j][k]));
}

inline
DataTypes::cplx_t 
WrappedArray::getEltC(unsigned int i, unsigned int j, unsigned int k, unsigned int m) const
{
    if (!iscomplex)
    {
      return (dat_r!=0)?dat_r[DataTypes::getRelIndex(shape,i,j,k,m)]:(boost::python::extract<DataTypes::real_t>(obj[i][j][k][m]));
    }  
    return (dat_c!=0)?dat_c[DataTypes::getRelIndex(shape,i,j,k,m)]:(boost::python::extract<DataTypes::cplx_t>(obj[i][j][k][m]));
}





}

#endif

