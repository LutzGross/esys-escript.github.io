
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef ES_DODGY_H
#define ES_DODGY_H

// This file contains functions which do bad things like violate the C standard
// Functions in here should not go into general use and should only be used on
// systems where the users are willing to accept responsibility for things working.

// Do not add functions here without explaining why they are problematic

namespace escript
{
// first param is output, next two are input data1 and2 the 3 ints are the size [in doubles] 
// of the first three params in order
typedef int (*binOpFnPtr)(double*, const double*, const double*, int, int, int);

/**
** Casts a void* to a function pointer taking 3 double* and 3 int
** \warning The C standard does not guarantee that void* and function pointer
** are the same size.
**
** Why do we have this? To allow us to pass function pointers between modules via python.
** To do this without using a custom object we use void*
*/
binOpFnPtr binOpFnPtrFromVoidPtr(void* v);

/**
Why do we have this? Only to act as an inverse to the above function
*/
void* voidPtrFromBinOpFnPtr(binOpFnPtr f);

}

#endif
