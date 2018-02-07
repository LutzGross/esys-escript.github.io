
/*****************************************************************************
*
* Copyright (c) 2009-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#ifndef LAPACKINVERSEHELPER_H
#define LAPACKINVERSEHELPER_H

namespace escript
{

/**
  Stores the memory required by different lapack implementations for matrix inverse
*/
class LapackInverseHelper
{
public:
	LapackInverseHelper(int N);
	~LapackInverseHelper();
	int invert(double* matrix);
private:
	int* piv;
	double* work;
	int N;
	int lwork;
};

}   // end of escript namespace

#endif

