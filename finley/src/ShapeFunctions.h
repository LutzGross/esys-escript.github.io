
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __FINLEY_SHAPEFUNCTIONS_H__
#define __FINLEY_SHAPEFUNCTIONS_H__

#include "system_dep.h"

#include "Finley.h"

#include <boost/shared_ptr.hpp>

#define S_INDEX(_J_,_I_,_NUMNODES_) INDEX2(_J_,_I_,_NUMNODES_)
#define DSDV_INDEX(_J_,_K_,_I_,_NUMNODES_,_DIM_) INDEX3(_J_,_K_,_I_,_NUMNODES_,_DIM_)

namespace finley {

typedef enum {
    Point1Shape,
    Line2Shape,
    Line3Shape,
    Line4Shape,
    Tri3Shape,
    Tri6Shape,
    Tri9Shape,
    Tri10Shape,
    Rec4Shape,
    Rec8Shape,
    Rec9Shape,
    Rec12Shape,
    Rec16Shape,
    Tet4Shape,
    Tet10Shape,
    Tet16Shape,
    Hex8Shape,
    Hex20Shape,
    Hex27Shape,
    Hex32Shape,
    NoShape  // marks end of list
} ShapeFunctionTypeId;


typedef void (ShapeFunction_Evaluation) (int, std::vector<double>&, std::vector<double>&, std::vector<double>&);

/// this struct holds the definition of the shape functions on an element
struct ShapeFunctionInfo {
    /// shape function type
    ShapeFunctionTypeId TypeId;
    /// the name in text form e.g. "Line2", "Rec12", ...
    const char* Name;
    /// number of spatial dimensions
    int numDim;
    /// number of shape functions
    int numShapes;
    /// order of the shape functions
    int numOrder;
    /// number of vertices of the element
    int numVertices;
    /// function to evaluate the shape functions at a set of points
    ShapeFunction_Evaluation* getValues;
};


/// this struct holds the evaluation of a shape function on a quadrature scheme
struct FINLEY_DLL_API ShapeFunction {
    ShapeFunction(ShapeFunctionTypeId id, int numQuadDim, int numQuadNodes,
                  const std::vector<double>& QuadNodes,
                  const std::vector<double>& QuadWeights);

    ShapeFunctionTypeId getTypeId(const char*);

    static const ShapeFunctionInfo* getInfo(ShapeFunctionTypeId id);

    /// shape function information
    const ShapeFunctionInfo* Type;
    /// number of quadrature points
    int numQuadNodes;
    /// coordinates of quadrature nodes
    std::vector<double> QuadNodes;
    /// weights of the quadrature scheme
    std::vector<double> QuadWeights;
    /// shape functions at quadrature nodes
    std::vector<double> S;
    /// derivative of the shape functions at quadrature nodes
    std::vector<double> dSdv;
};

typedef boost::shared_ptr<const ShapeFunction> const_ShapeFunction_ptr;

ShapeFunction_Evaluation Shape_Point1;
ShapeFunction_Evaluation Shape_Line2;
ShapeFunction_Evaluation Shape_Line3;
ShapeFunction_Evaluation Shape_Line4;
ShapeFunction_Evaluation Shape_Tri3;
ShapeFunction_Evaluation Shape_Tri6;
ShapeFunction_Evaluation Shape_Tri9;
ShapeFunction_Evaluation Shape_Tri10;
ShapeFunction_Evaluation Shape_Rec4;
ShapeFunction_Evaluation Shape_Rec8;
ShapeFunction_Evaluation Shape_Rec9;
ShapeFunction_Evaluation Shape_Rec12;
ShapeFunction_Evaluation Shape_Rec16;
ShapeFunction_Evaluation Shape_Tet4;
ShapeFunction_Evaluation Shape_Tet10;
ShapeFunction_Evaluation Shape_Tet16;
ShapeFunction_Evaluation Shape_Hex8;
ShapeFunction_Evaluation Shape_Hex20;
ShapeFunction_Evaluation Shape_Hex27;
ShapeFunction_Evaluation Shape_Hex32;

} // namespace finley

#endif // __FINLEY_SHAPEFUNCTIONS_H__

