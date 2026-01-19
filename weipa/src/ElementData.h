
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

#ifndef __WEIPA_ELEMENTDATA_H__
#define __WEIPA_ELEMENTDATA_H__

#include <weipa/weipa.h>
#include <weipa/vtkCellType.h>
#include <ostream>

namespace weipa {
    
typedef enum {
    ZONETYPE_UNKNOWN=0,
    ZONETYPE_BEAM=VTK_LINE,
    ZONETYPE_HEX=VTK_HEXAHEDRON,
    ZONETYPE_POLYGON=VTK_POLYGON,
    ZONETYPE_QUAD=VTK_QUAD,
    ZONETYPE_TET=VTK_TETRA,
    ZONETYPE_TRIANGLE=VTK_TRIANGLE
} ZoneType;

/// \brief This struct holds a mask (0's and 1's) that indicates which quad
///        nodes contribute to a sub-element when full element order is used.
///        factor[i] contains the number of non-zeroes in mask[i].
struct QuadMaskInfo {
    std::vector<IntVec> mask;
    IntVec factor;
};

/// \brief Stores and manipulates one type of domain elements.
///
/// \note The corresponding mesh nodes are not part of this class but are
///       stored in a NodeData instance.
class ElementData
{
public:
    /// \brief Writes connectivity data to a stream in VTK text format.
    virtual void writeConnectivityVTK(std::ostream& os) = 0;

    /// \brief Returns the names of the meshes associated with the elements.
    virtual StringVec getMeshNames() const = 0;

    /// \brief Returns a vector with the mesh variable names.
    virtual StringVec getVarNames() const = 0;

    /// \brief Returns the number of elements.
    virtual int getNumElements() const = 0;

    /// \brief Returns the number of nodes per element.
    virtual int getNodesPerElement() const = 0;

    /// \brief Returns the number of "ghost" elements.
    virtual int getGhostCount() const = 0;

    /// \brief Returns the element type.
    virtual ZoneType getType() const = 0;

    /// \brief Returns a vector of the node IDs used by the elements.
    virtual const IntVec& getNodeList() const = 0;

    /// \brief Returns a vector of element IDs.
    virtual const IntVec& getIDs() const = 0;

    /// \brief Returns the NodeData instance used by the elements.
    virtual NodeData_ptr getNodes() const = 0;

    /// \brief Returns the reduced elements if available.
    virtual ElementData_ptr getReducedElements() const = 0;
 
    /// \brief Returns a QuadMaskInfo structure for given functionspace code.
    virtual const QuadMaskInfo& getQuadMask(int fsCode) const = 0;
 
    /// \brief If the original element type is not supported they are
    ///        subdivided into N smaller elements (e.g. one Rec9 -> four Rec4)
    ///        and this method returns the multiplication factor N.
    virtual int getElementFactor() const = 0;

protected:
    virtual ~ElementData() {}
};

} // namespace weipa

#endif // __WEIPA_ELEMENTDATA_H__

