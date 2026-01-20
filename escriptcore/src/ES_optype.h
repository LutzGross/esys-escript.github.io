
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

#ifndef __ESCRIPT_ESOPTYPE_H__
#define __ESCRIPT_ESOPTYPE_H__

#include "system_dep.h"

#include <string>

namespace escript
{

// For the purposes of unit testing and maintaining sanity, it is important that this enum be contiguous
enum ES_optype
{
	UNKNOWNOP=0,
	IDENTITY=1,
	ADD=2,
	SUB=3,
	MUL=4,
	DIV=5,
	POW=6,
	SIN=POW+1,
	COS=SIN+1,
	TAN=SIN+2,
	ASIN=SIN+3,
	ACOS=SIN+4,
	ATAN=SIN+5,
	SINH=SIN+6,
	COSH=SIN+7,
	TANH=SIN+8,
	ERF=SIN+9,
	ASINH=SIN+10,
	ACOSH=SIN+11,
	ATANH=SIN+12,
	LOG10=ATANH+1,
	LOG=LOG10+1,
	SIGN=LOG10+2,
	ABS=LOG10+3,
	NEG=LOG10+4,
	POS=LOG10+5,
	EXP=LOG10+6,
	SQRT=LOG10+7,
	RECIP=LOG10+8,
	GZ=RECIP+1,
	LZ=GZ+1,	// <0
	GEZ=GZ+2,	// >=0
	LEZ=GZ+3,	// <=0
	NEZ=GZ+4,	// >=0
	EZ=GZ+5,
	SYM=EZ+1,
	NSYM=SYM+1,
	PROD=NSYM+1,
	TRANS=PROD+1,
	TRACE=TRANS+1,
	SWAP=TRACE+1,
	MINVAL=SWAP+1,
	MAXVAL=MINVAL+1,
	CONDEVAL=MAXVAL+1,
	HER=CONDEVAL+1,		// hermitian
	NHER=HER+1,              // antihermitian
	REAL=NHER+1,
	IMAG=REAL+1,
	CONJ=IMAG+1,
	LESS=CONJ+1,		// a<b
	GREATER=LESS+1,
	GREATER_EQUAL=GREATER+1,
	LESS_EQUAL=GREATER_EQUAL+1,
	PHS=LESS_EQUAL+1,	// phase
	PROM=PHS+1         // promote real to complex
};

ESCRIPT_DLL_API
const std::string&
opToString(ES_optype op);

enum ES_opgroup
{
   G_UNKNOWN,
   G_IDENTITY,
   G_BINARY,            // pointwise operations with two arguments
   G_UNARY,             // pointwise operations with one argument
   G_UNARY_P,           // pointwise operations with one argument, requiring a parameter
   G_UNARY_R,		// pointwise operations with one argument, always real output
   G_NP1OUT,            // non-pointwise op with one output
   G_NP1OUT_P,          // non-pointwise op with one output requiring a parameter
   G_TENSORPROD,        // general tensor product
   G_NP1OUT_2P,         // non-pointwise op with one output requiring two params
   G_REDUCTION,         // non-pointwise unary op with a scalar output
   G_CONDEVAL,
   G_UNARY_C,           // pointwise operations with one argument, always cplx output
   G_UNARY_PR           // G_UNARY_P but always real output   
};



ES_opgroup
getOpgroup(ES_optype op);

const std::string&
groupToString(ES_opgroup g);
}

#endif
