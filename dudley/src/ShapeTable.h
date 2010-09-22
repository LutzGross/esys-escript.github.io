
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/* Shape Function info
These tables are a much simplified version of content from finley's ShapeFunctions files

This file is not to be included in .h files - only .c files should have any use for it
*/

#ifndef SHAPETABLE_DUDLEY
#define SHAPETABLE_DUDLEY

#include "paso/Common.h"	// I just want the types not all the includes that get dragged in - fix that

// These are constructed from dsdv in ShapeFunction.c in finley
// The first two are just there for functions that want a pointer
static const double DTDV_0D[1][1]={{0}};
static const double DTDV_1D[2][2]={{-1.,1},{-1.,1.}};
static const double DTDV_2D[3][2]={{-1,-1}, {1,0}, {0,1}};
static const double DTDV_3D[4][3]={{-1, -1, -1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};


// If these appear to be in a different order to finley it is because finley uses macros to hide Fortran array ordering
static const double DTDV_2D_alt[3*3][2]={{-1,1}, {0,-1.}, {0,1},
{-1,1}, {0,-1.}, {0,1},
{-1,1}, {0,-1.}, {0,1}
};
// the repetition is a hack
// Why didn't I just reorder DTDV_2D?   Well some code apparently depends on the order as written.
// should probably fix that


// [0] is reduced quadrature, [1] is full quadrature
// in order the positions are POINT, LINE, TRI, TET
static const double QuadWeight[4][2]={{0, 0}, {1., 0.5}, {0.5, 1./6}, {1./6, 1./24}};

static const dim_t QuadNums[4][2] ={{0,0}, {1,2}, {1,3}, {1,4}};

//shape functions at quadrature nodes
bool_t getQuadShape(dim_t sim, bool_t reduced, const double** shapearr);







#endif 



