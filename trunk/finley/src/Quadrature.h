
/*****************************************************************************
*
* Copyright (c) 2003-2015 by University of Queensland
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


/****************************************************************************

  Finley: integration schemes for element shapes Tri, Quad, Hex, Tet, Line,
          Point

*****************************************************************************/

#ifndef __FINLEY_QUADRATURE_H__
#define __FINLEY_QUADRATURE_H__

#include "Finley.h"

#define MAX_numQuadNodesLine 10

namespace finley {

typedef enum {
    PointQuad,
    LineQuad,
    TriQuad,
    RecQuad,
    TetQuad,
    HexQuad,
    NoQuad   // marks end of list
} QuadTypeId;

typedef void (Quad_getNodes) (int, std::vector<double>&, std::vector<double>&);
typedef int (Quad_getNumNodes) (int);
typedef int (Quad_getMacro) (int numSubElements, int numQuadNodes,
                             const double* quadNodes,
                             const double* quadWeights,
                             int numF, const double* dFdv,
                             int new_len, double* new_quadNodes,
                             double* new_quadWeights, double* new_dFfv);

struct QuadInfo {
    /// quadrature type id
    QuadTypeId TypeId;
    /// the name in text form e.g. "Line", "Rec", ...
    const char* Name;
    /// number of spatial dimensions
    int numDim;
    /// number of vertices of the element
    int numVertices;
    /// function that returns the quadrature points for a given order
    Quad_getNodes* getQuadNodes;
    /// function that returns the number of quadrature nodes for a given
    /// accuracy order
    Quad_getNumNodes* getNumQuadNodes;
    /// transfers a given quadrature scheme to a macro element structure
    Quad_getMacro *getMacro;
};


/****** Interfaces ******/

Quad_getMacro Quad_MacroPoint;
Quad_getMacro Quad_MacroLine;
Quad_getMacro Quad_MacroTri;
Quad_getMacro Quad_MacroRec;
Quad_getMacro Quad_MacroTet;
Quad_getMacro Quad_MacroHex;

Quad_getNodes Quad_getNodesTri;
Quad_getNodes Quad_getNodesTet;
Quad_getNodes Quad_getNodesRec;
Quad_getNodes Quad_getNodesHex;
Quad_getNodes Quad_getNodesLine;
Quad_getNodes Quad_getNodesPoint;
Quad_getNodes Quad_getNodesTriOnFace;
Quad_getNodes Quad_getNodesRecOnFace;
Quad_getNodes Quad_getNodesLineOnFace;
Quad_getNodes Quad_getNodesPointOnFace;
Quad_getNodes Quad_getNodesTriMacro;
Quad_getNodes Quad_getNodesTetMacro;
Quad_getNodes Quad_getNodesRecMacro;
Quad_getNodes Quad_getNodesHexMacro;
Quad_getNodes Quad_getNodesLineMacro;

Quad_getNumNodes Quad_getNumNodesPoint;
Quad_getNumNodes Quad_getNumNodesLine;
Quad_getNumNodes Quad_getNumNodesTri;
Quad_getNumNodes Quad_getNumNodesRec;
Quad_getNumNodes Quad_getNumNodesTet;
Quad_getNumNodes Quad_getNumNodesHex;

void Quad_makeNodesOnFace(int, int, double*, double*, Quad_getNodes);
const QuadInfo* QuadInfo_getInfo(QuadTypeId id);

} // namespace finley

#endif // __FINLEY_QUADRATURE_H__

