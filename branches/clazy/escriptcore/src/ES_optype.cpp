
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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


#include "ES_optype.h"
#include <string>

using namespace escript;

namespace
{

std::string ES_opstrings[]={"UNKNOWN","IDENTITY","+","-","*","/","^",
                        "sin","cos","tan",
                        "asin","acos","atan","sinh","cosh","tanh","erf",
                        "asinh","acosh","atanh",
                        "log10","log","sign","abs","neg","pos","exp","sqrt",
                        "1/","where>0","where<0","where>=0","where<=0", "where<>0","where=0",
                        "symmetric","antisymmetric",
                        "prod",
                        "transpose", "trace",
                        "swapaxes",
                        "minval", "maxval",
                        "condEval",
                        "hermitian","antihermitian",
			"real","imaginary","conjugate",
			"<", ">", ">=", "<=",
			"phase",
            "promote"
};


ES_opgroup opgroups[]={G_UNKNOWN,G_IDENTITY,G_BINARY,G_BINARY,G_BINARY,G_BINARY, G_BINARY,
                        G_UNARY,G_UNARY,G_UNARY, //10
                        G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,        // 17
                        G_UNARY,G_UNARY,G_UNARY,                                        // 20
                        G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY,        // 28
                        G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY, G_UNARY_P, G_UNARY_P,          // 35
                        G_NP1OUT,G_NP1OUT,
                        G_TENSORPROD,
                        G_NP1OUT_P, G_NP1OUT_P,
                        G_NP1OUT_2P,
                        G_REDUCTION, G_REDUCTION,
                        G_CONDEVAL,
                        G_UNARY,G_UNARY,
                        G_UNARY_R, G_UNARY_R, G_UNARY,
			G_UNARY_R, G_UNARY_R, G_UNARY_R, G_UNARY_R,
			G_UNARY_R,
            G_UNARY_C
};


int ES_opcount=55;
}

// Return a string representing the operation
const std::string&
escript::opToString(ES_optype op)
{
  if (op<0 || op>=ES_opcount) 
  {
    op=UNKNOWNOP;
  }
  return ES_opstrings[op];
}

ES_opgroup
escript::getOpgroup(ES_optype op)
{
  return opgroups[op];
}
