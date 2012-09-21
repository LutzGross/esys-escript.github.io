
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/

/* Shape Function info
These tables are a much simplified version of content from finley's ShapeFunctions files

This file is not to be included in .h files - only .c files should have any use for it
*/

#ifndef SHAPETABLE_DUDLEY
#define SHAPETABLE_DUDLEY

#include "esysUtils/types.h"	

#include "ElementType.h"

/* These are constructed from dsdv in ShapeFunction.c in finley
   The first two are just there for functions that want a pointer
*/
static const double DTDV_0D[1][1] = { {0} };
static const double DTDV_1D[2][2] = { {-1., 1}, {-1., 1.} };

/* The repetition here is a hack to prevent out of bounds access */
static const double DTDV_2D[3 * 3][2] = { {-1, 1}, {0, -1.}, {0, 1},
{-1, 1}, {0, -1.}, {0, 1},
{-1, 1}, {0, -1.}, {0, 1}
};
static const double DTDV_3D[4][3] = { {-1, -1, -1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };

/* Index by the following by Dudley_ElementTypeID
 * The number of local dimensions (as opposed to dimension of the embedding space) */
static const dim_t localDims[8] = { 0, 1, 2, 3, 0, 1, 2, 0 };
static const dim_t Dims[8] = { 0, 1, 2, 3, 1, 2, 3, 0 };

/* the following lists are only used for face elements defined by numNodesOnFace>0 */
static const dim_t numNodesOnFaceMap[8] = { 1, 2, 3, 4, 1, 2, 4, -1 };	/* if the element is allowed as a face element, numNodesOnFace defines the number of nodes defining the face */
static const dim_t shiftNodesMap[8][4] = { {0}, {1, 0}, {1, 2, 0}, {-1}, {0, 1, 2}, {1, 0, 2}, {1, 2, 0, 3}, {0} };	/* defines a permutation of the nodes which rotates the nodes on the face */
static const dim_t reverseNodesMap[8][4] = { {-1}, {-1}, {0, 2, 1}, {-1}, {-1}, {-1}, {0, 2, 1, 3}, {0} };	/* reverses the order of the nodes on a face. the permutation has keep 0 fixed. */
					      /* shiftNodes={-1} or reverseNodes={-1} are ignored. */

/* [0] is reduced quadrature, [1] is full quadrature */
/* in order the positions are POINT, LINE, TRI, TET */
static const double QuadWeight[4][2] = { {0, 0}, {1., 0.5}, {0.5, 1. / 6}, {1. / 6, 1. / 24} };

/* number of quadrature points per element */
static const dim_t QuadNums[4][2] = { {0, 0}, {1, 2}, {1, 3}, {1, 4} };

/*shape functions at quadrature nodes */
bool_t getQuadShape(dim_t sim, bool_t reduced, const double **shapearr);

const char *getElementName(Dudley_ElementTypeId id);

#endif
