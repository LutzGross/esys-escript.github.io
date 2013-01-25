
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


/** \file WrappedArray.h */

#ifndef WrappedArray_20081202_H
#define WrappedArray_20081202_H
#include "system_dep.h"
#include "DataTypes.h"
#include "boost/python/extract.hpp"

namespace escript
{

class WrappedArray
{
public:
	WrappedArray(const boost::python::object& obj_in);
	unsigned int getRank() const;
	const DataTypes::ShapeType& getShape() const;
	double getElt() const;
	double getElt(unsigned int i) const;
	double getElt(unsigned int i, unsigned int j) const;
	double getElt(unsigned int i, unsigned int j, unsigned int k) const;
	double getElt(unsigned int i, unsigned int j, unsigned int k, unsigned int m) const;
private:
	const boost::python::object& obj;
	int rank;
	escript::DataTypes::ShapeType shape;
	double m_scalar;
};

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
	return m_scalar;
}


inline double
WrappedArray::getElt(unsigned int i) const
{
	return boost::python::extract<double>(obj[i]);
}

inline
double 
WrappedArray::getElt(unsigned int i, unsigned int j) const
{
	return boost::python::extract<double>(obj[i][j]);
}

inline
double 
WrappedArray::getElt(unsigned int i, unsigned int j, unsigned int k) const
{
	return boost::python::extract<double>(obj[i][j][k]);
}

inline
double 
WrappedArray::getElt(unsigned int i, unsigned int j, unsigned int k, unsigned int m) const
{
	return boost::python::extract<double>(obj[i][j][k][m]);
}

}

#endif