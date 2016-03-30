
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


#include "ArrayOps.h"

using namespace escript;

namespace escript
{

bool supports_cplx(escript::ESFunction operation)
{
    switch (operation)
    {
    case NEGF:
    case SINF: 
    case COSF: 
    case TANF: 
    case ASINF: 
    case ACOSF: 
    case ATANF: 
    case SINHF: 
    case COSHF: 
    case TANHF: return true;
    case ERFF: return false;
    case ASINHF: 
    case ACOSHF: 
    case ATANHF: 
    case LOG10F: 
    case LOGF: return true;
    case SIGNF: return false;
    case ABSF: 
    case EXPF: 
    case SQRTF: return true;
    case EQZEROF:
    case NEQZEROF:return true;
    case GTZEROF:
    case GEZEROF:
    case LTZEROF:
    case LEZEROF: return false;   
    case CONJF: return true;
    case REALF: return true;
    case IMAGF: return true;
    case INVF: return true;
    default:
      return false;	// let's be conservative
  }  
}

bool always_real(escript::ESFunction operation)
{
    return ((operation==REALF) || (operation==IMAGF) || (operation==EQZEROF) || (operation==NEQZEROF) || (operation==ABSF));
}


}