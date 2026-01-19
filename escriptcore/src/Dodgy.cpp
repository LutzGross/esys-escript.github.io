
/*****************************************************************************
*
* Copyright (c) 2009-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
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
