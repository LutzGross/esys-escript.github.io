
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


/****************************************************************************

  Finley: Reference elements

*****************************************************************************/

#ifndef __FINLEY_REFERENCEELEMENTS_H__
#define __FINLEY_REFERENCEELEMENTS_H__

#include "system_dep.h"

#include "Finley.h"
#include "ShapeFunctions.h"
#include "Quadrature.h"

//  The ids of the allowed reference elements:
#define MAX_numNodes 64
#define MAX_numSubElements 8
#define MAX_numSides 2

namespace finley {

typedef enum {
    Point1,
    Line2,
    Line3,
    Line4,
    Tri3,
    Tri6,
    Tri9,
    Tri10,
    Rec4,
    Rec8,
    Rec9,
    Rec12,
    Rec16,
    Tet4,
    Tet10,
    Tet16,
    Hex8,
    Hex20,
    Hex27,
    Hex32,
    Line2Face,
    Line3Face,
    Line4Face,
    Tri3Face,
    Tri6Face,
    Tri9Face,
    Tri10Face,
    Rec4Face,
    Rec8Face,
    Rec9Face,
    Rec12Face,
    Rec16Face,
    Tet4Face,
    Tet10Face,
    Tet16Face,
    Hex8Face,
    Hex20Face,
    Hex27Face,
    Hex32Face,
    Point1_Contact,
    Line2_Contact,
    Line3_Contact,
    Line4_Contact,
    Tri3_Contact,
    Tri6_Contact,
    Tri9_Contact,
    Tri10_Contact,
    Rec4_Contact,
    Rec8_Contact,
    Rec9_Contact,
    Rec12_Contact,
    Rec16_Contact,
    Line2Face_Contact,
    Line3Face_Contact,
    Line4Face_Contact,
    Tri3Face_Contact,
    Tri6Face_Contact,
    Tri9Face_Contact,
    Tri10Face_Contact,
    Rec4Face_Contact,
    Rec8Face_Contact,
    Rec9Face_Contact,
    Rec12Face_Contact,
    Rec16Face_Contact,
    Tet4Face_Contact,
    Tet10Face_Contact,
    Tet16Face_Contact,
    Hex8Face_Contact,
    Hex20Face_Contact,
    Hex27Face_Contact,
    Hex32Face_Contact,
    Line3Macro,
    Tri6Macro,
    Rec9Macro,
    Tet10Macro,
    Hex27Macro,
    NoRef // marks end of list
} ElementTypeId;


/// this struct holds the definition of the reference element
struct FINLEY_DLL_API ReferenceElementInfo {
    /// the type
    ElementTypeId TypeId;
    /// the name in text form e.g. "Line1", "Rec12", ...
    const char* Name;
    /// number of nodes defining the element
    int numNodes;
    /// number of subelements (>1 if macro elements are used)
    int numSubElements;
    /// the number of sides the element supports (=2 if contact elements are
    /// used, otherwise =1)
    int numSides;
    /// offset to the side nodes: offsets[s]...offsets[s+1]-1 refer to the
    /// nodes to be used for side s
    int offsets[MAX_numSides+1];

    /// type id of the linear version of the element
    ElementTypeId LinearTypeId;
    /// stores the list of nodes defining the linear or macro element
    int linearNodes[MAX_numNodes*MAX_numSides];
    /// quadrature scheme
    QuadTypeId Quadrature;
    /// shape function for parametrization of the element
    ShapeFunctionTypeId Parametrization;
    /// shape function for the basis functions
    ShapeFunctionTypeId BasisFunctions;

    /// the list of nodes defining the subelements, i.e.
    /// subElementNodes[INDEX2(i,s,BasisFunctions->numShape*numSides)] is
    /// the i-th node in the s-th subelement
    int subElementNodes[MAX_numNodes*MAX_numSides*MAX_numSubElements];

    /// deprecated
    int numRelevantGeoNodes;
    int relevantGeoNodes[MAX_numNodes];

    /// if the element is allowed as a face element, numNodesOnFace defines
    /// the number of nodes defining the face
    int numNodesOnFace;

    // the following lists are only used for face elements defined by
    // numNodesOnFace>0:

    /// list of the nodes defining the face
    int faceNodes[MAX_numNodes];

    // shiftNodes={-1} or reverseNodes={-1} are ignored.
    /// defines a permutation of the nodes which rotates the nodes on the face
    int shiftNodes[MAX_numNodes];
    /// reverses the order of the nodes on a face. The permutation has to keep
    /// 0 fixed.
    int reverseNodes[MAX_numNodes];
};


/// this struct holds the realization of a reference element
struct FINLEY_DLL_API ReferenceElement {
    /// constructor with type ID and integration order
    ReferenceElement(ElementTypeId id, int order);

    /// destructor
    ~ReferenceElement();

    /// returns the element information structure for the given type id
    static const ReferenceElementInfo* getInfo(ElementTypeId id);

    /// returns the element type id from its textual representation
    static ElementTypeId getTypeId(const char*);

    ///
    int getNumNodes() const { return Type->numNodes; }

    /// type of the reference element
    const ReferenceElementInfo* Type;
    /// type of the linear reference element
    const ReferenceElementInfo* LinearType;
    /// used integration order
    int integrationOrder;
    int numNodes;
    int numLocalDim;
    int numLinearNodes;
    const_ShapeFunction_ptr Parametrization;
    const_ShapeFunction_ptr BasisFunctions;
    const_ShapeFunction_ptr LinearBasisFunctions;
    /// pointer to derivatives to basis function corresponding to the
    /// Parametrization of quad points
    double* DBasisFunctionDv;
    /// if true indicates that DBasisFunctionDv is shared with another object
    /// which is managing it
    bool DBasisFunctionDvShared;
};

typedef boost::shared_ptr<ReferenceElement> ReferenceElement_ptr;
typedef boost::shared_ptr<const ReferenceElement> const_ReferenceElement_ptr;

} // namespace finley

#endif // __FINLEY_REFERENCEELEMENTS_H__

