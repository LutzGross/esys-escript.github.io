
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


/**************************************************************/

/*   Dudley: integration schemes for element shapes Tri, Quad, Hex, Tet, Line, Point */

/**************************************************************/

#ifndef INC_DUDLEY_QUADRATURE
#define INC_DUDLEY_QUADRATURE

/**************************************************************/

#include "Dudley.h"

/**************************************************************/

#define MAX_numQuadNodesLine 10

typedef enum {
  PointQuad,
  LineQuad,
  TriQuad,
  TetQuad,
  NoQuad   /* marks end of list */
} Dudley_QuadTypeId;

typedef void (Dudley_Quad_getNodes) (dim_t,double*,double*);
typedef dim_t (Dudley_Quad_getNumNodes) (dim_t);
typedef dim_t(Dudley_Quad_getMacro)(dim_t numSubElements, int numQuadNodes, double* quadNodes, double* quadWeights, 
                                        dim_t numF, double* dFdv, 
					dim_t new_len, double* new_quadNodes, double* new_quadWeights,
                                        double* new_dFfv );

typedef struct Dudley_QuadInfo {
  Dudley_QuadTypeId TypeId;                  /* the id */
  char* Name;                                /* the name in text form e.g. Line,Rec,... */
  dim_t numDim;                              /* spacial dimension */
  dim_t numVertices;                         /* number of vertices of the element */
  Dudley_Quad_getNodes* getQuadNodes;        /* function to set the quadrature points for a given order */
  Dudley_Quad_getNumNodes* getNumQuadNodes;  /* function selects the number of quadrature nodes for a given accuracy order */
}  Dudley_QuadInfo;

/**************************************************************/

/*     Interfaces: */

Dudley_Quad_getNodes Dudley_Quad_getNodesTri;
Dudley_Quad_getNodes Dudley_Quad_getNodesTet;
Dudley_Quad_getNodes Dudley_Quad_getNodesLine;
Dudley_Quad_getNodes Dudley_Quad_getNodesPoint;
Dudley_Quad_getNodes Dudley_Quad_getNodesTriOnFace;
Dudley_Quad_getNodes Dudley_Quad_getNodesLineOnFace;
Dudley_Quad_getNodes Dudley_Quad_getNodesPointOnFace;

Dudley_Quad_getNumNodes Dudley_Quad_getNumNodesPoint;
Dudley_Quad_getNumNodes Dudley_Quad_getNumNodesLine;
Dudley_Quad_getNumNodes Dudley_Quad_getNumNodesTri;
Dudley_Quad_getNumNodes Dudley_Quad_getNumNodesTet;

void Dudley_Quad_makeNodesOnFace(dim_t, dim_t,double*,double*, Dudley_Quad_getNodes);
Dudley_QuadInfo* Dudley_QuadInfo_getInfo(Dudley_QuadTypeId id);

#endif /* #ifndef INC_DUDLEY_QUADRATURE */

