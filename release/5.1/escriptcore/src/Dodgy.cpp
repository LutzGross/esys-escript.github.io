
/*****************************************************************************
*
* Copyright (c) 2009-2017 by The University of Queensland
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

#include "Dodgy.h"

namespace escript
{

binOpFnPtr binOpFnPtrFromVoidPtr(void* v)
{
	return reinterpret_cast<binOpFnPtr>(v);

}




void* voidPtrFromBinOpFnPtr(binOpFnPtr f)
{

	return reinterpret_cast<void*>(f);
}





}