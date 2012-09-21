
/*****************************************************************************
*
* Copyright (c) 2009-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
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