
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#include "ArrayOps.h"

using namespace escript;

namespace escript
{

bool supports_cplx(escript::ES_optype operation)
{
    switch (operation)
    {
    case NEG:
    case SIN: 
    case COS: 
    case TAN: 
    case ASIN: 
    case ACOS: 
    case ATAN: 
    case SINH: 
    case COSH: 
    case TANH: return true;
    case ERF: return false;
    case ASINH: 
    case ACOSH: 
    case ATANH: 
    case LOG10: 
    case LOG: return true;
    case SIGN: return false;
    case ABS: 
    case EXP: 
    case SQRT: return true;
    case EZ:
    case NEZ:return true;
    case GZ:
    case GEZ:
    case LZ:
    case LEZ: return false;   
    case CONJ: return true;
    case REAL: return true;
    case IMAG: return true;
    case RECIP: return true;
    case PHS: return true;
    default:
      return false;	// let's be conservative
  }  
}

bool always_real(escript::ES_optype operation)
{
    return ((operation==REAL) || (operation==IMAG) || (operation==EZ) || (operation==NEZ) || (operation==ABS) || (operation==PHS));
}


}
