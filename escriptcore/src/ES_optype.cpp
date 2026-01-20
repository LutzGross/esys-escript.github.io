
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
                        G_UNARY,G_UNARY,G_UNARY,G_UNARY_R,G_UNARY,G_UNARY,G_UNARY,G_UNARY,        // 28
                        G_UNARY,G_UNARY,G_UNARY,G_UNARY,G_UNARY, G_UNARY_PR, G_UNARY_PR,          // 35
                        G_NP1OUT,G_NP1OUT,
                        G_TENSORPROD,
                        G_NP1OUT_P, G_NP1OUT_P,
                        G_NP1OUT_2P,
                        G_REDUCTION, G_REDUCTION,
                        G_CONDEVAL,
                        G_NP1OUT,G_NP1OUT,
                        G_UNARY_R, G_UNARY_R, G_UNARY,
			G_UNARY_R, G_UNARY_R, G_UNARY_R, G_UNARY_R,
			G_UNARY_R,
            G_UNARY_C
};

std::string ES_groupstrings[]={
   "G_UNKNOWN",
   "G_IDENTITY",
   "G_BINARY",            // pointwise operations with two arguments
   "G_UNARY",             // pointwise operations with one argument
   "G_UNARY_P",           // pointwise operations with one argument, requiring a parameter
   "G_UNARY_R",		// pointwise operations with one argument, always real output
   "G_NP1OUT",            // non-pointwise op with one output
   "G_NP1OUT_P",          // non-pointwise op with one output requiring a parameter
   "G_TENSORPROD",        // general tensor product
   "G_NP1OUT_2P",         // non-pointwise op with one output requiring two params
   "G_REDUCTION",         // non-pointwise unary op with a scalar output
   "G_CONDEVAL",
   "G_UNARY_C"            // pointwise operations with one argument, always cplx output    
   "G_UNARY_PR"           // G_UNARY_P but always real output
};

int ES_opcount=55;
int ES_groupcount=14;
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

const std::string&
escript::groupToString(ES_opgroup g)
{
  if (g<0 || g>=ES_groupcount) 
  {
    g=G_UNKNOWN;
  }
  return ES_groupstrings[g];
}
