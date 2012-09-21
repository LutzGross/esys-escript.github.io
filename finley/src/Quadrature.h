
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


/************************************************************************************/

/*   Finley: integration schemes for element shapes Tri, Quad, Hex, Tet, Line, Point */

/************************************************************************************/

#ifndef INC_FINLEY_QUADRATURE
#define INC_FINLEY_QUADRATURE

/************************************************************************************/

#include "Finley.h"

/************************************************************************************/

#define MAX_numQuadNodesLine 10

typedef enum {
  PointQuad,
  LineQuad,
  TriQuad,
  RecQuad,
  TetQuad,
  HexQuad,
  NoQuad   /* marks end of list */
} Finley_QuadTypeId;

typedef void (Finley_Quad_getNodes) (dim_t,double*,double*);
typedef dim_t (Finley_Quad_getNumNodes) (dim_t);
typedef dim_t(Finley_Quad_getMacro)(dim_t numSubElements, int numQuadNodes, double* quadNodes, double* quadWeights, 
                                        dim_t numF, double* dFdv, 
					dim_t new_len, double* new_quadNodes, double* new_quadWeights,
                                        double* new_dFfv );

typedef struct Finley_QuadInfo {
  Finley_QuadTypeId TypeId;                  /* the id */
  char* Name;                                /* the name in text form e.g. Line,Rec,... */
  dim_t numDim;                              /* spatial dimension */
  dim_t numVertices;                         /* number of vertices of the element */
  Finley_Quad_getNodes* getQuadNodes;        /* function to set the quadrature points for a given order */
  Finley_Quad_getNumNodes* getNumQuadNodes;  /* function selects the number of quadrature nodes for a given accuracy order */
  Finley_Quad_getMacro *getMacro;         		 /* transfers a given quadrature scheme to a macro element structure */
}  Finley_QuadInfo;

/************************************************************************************/

/*     Interfaces: */


Finley_Quad_getMacro Finley_Quad_MacroPoint;
Finley_Quad_getMacro Finley_Quad_MacroLine;
Finley_Quad_getMacro Finley_Quad_MacroTri;
Finley_Quad_getMacro Finley_Quad_MacroRec;
Finley_Quad_getMacro Finley_Quad_MacroTet;
Finley_Quad_getMacro Finley_Quad_MacroHex;


Finley_Quad_getNodes Finley_Quad_getNodesTri;
Finley_Quad_getNodes Finley_Quad_getNodesTet;
Finley_Quad_getNodes Finley_Quad_getNodesRec;
Finley_Quad_getNodes Finley_Quad_getNodesHex;
Finley_Quad_getNodes Finley_Quad_getNodesLine;
Finley_Quad_getNodes Finley_Quad_getNodesPoint;
Finley_Quad_getNodes Finley_Quad_getNodesTriOnFace;
Finley_Quad_getNodes Finley_Quad_getNodesRecOnFace;
Finley_Quad_getNodes Finley_Quad_getNodesLineOnFace;
Finley_Quad_getNodes Finley_Quad_getNodesPointOnFace;
Finley_Quad_getNodes Finley_Quad_getNodesTriMacro;
Finley_Quad_getNodes Finley_Quad_getNodesTetMacro;
Finley_Quad_getNodes Finley_Quad_getNodesRecMacro;
Finley_Quad_getNodes Finley_Quad_getNodesHexMacro;
Finley_Quad_getNodes Finley_Quad_getNodesLineMacro;



Finley_Quad_getNumNodes Finley_Quad_getNumNodesPoint;
Finley_Quad_getNumNodes Finley_Quad_getNumNodesLine;
Finley_Quad_getNumNodes Finley_Quad_getNumNodesTri;
Finley_Quad_getNumNodes Finley_Quad_getNumNodesRec;
Finley_Quad_getNumNodes Finley_Quad_getNumNodesTet;
Finley_Quad_getNumNodes Finley_Quad_getNumNodesHex;

void Finley_Quad_makeNodesOnFace(dim_t, dim_t,double*,double*, Finley_Quad_getNodes);
Finley_QuadInfo* Finley_QuadInfo_getInfo(Finley_QuadTypeId id);

#endif /* #ifndef INC_FINLEY_QUADRATURE */

