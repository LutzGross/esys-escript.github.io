
/*******************************************************
*
* Copyright (c) 2009-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

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