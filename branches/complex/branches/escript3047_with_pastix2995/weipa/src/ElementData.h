
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

#ifndef __ELEMENTDATA_H__
#define __ELEMENTDATA_H__

#include <weipa/weipa.h>

extern "C" {
#include <finley/ReferenceElements.h> // for ElementTypeId
}

#include <finley/vtkCellType.h>
#include <ostream>

class DBfile;
class NcFile;

struct Finley_ElementFile;

namespace weipa {
    
class NodeData;

typedef enum {
    ZONETYPE_UNKNOWN=0,
    ZONETYPE_BEAM=VTK_LINE,
    ZONETYPE_HEX=VTK_HEXAHEDRON,
    ZONETYPE_POLYGON=VTK_POLYGON,
    ZONETYPE_QUAD=VTK_QUAD,
    ZONETYPE_TET=VTK_TETRA,
    ZONETYPE_TRIANGLE=VTK_TRIANGLE
} ZoneType;

/// \brief Holds information that is used to convert from finley element types
///        to elements supported by Silo and VTK.
struct FinleyElementInfo
{
    ZoneType elementType, reducedElementType;
    int elementFactor;
    int elementSize, reducedElementSize;
    const size_t* multiCellIndices;
    bool useQuadNodes;
    int quadDim;
};

/// \brief This struct holds a mask (0's and 1's) that indicates which quad
///        nodes contribute to a sub-element when full element order is used.
///        factor[i] contains the number of non-zeroes in mask[i].
struct QuadMaskInfo {
    std::vector<IntVec> mask;
    IntVec factor;
};

/// \brief Stores and manipulates one type of finley mesh elements (cells,
///        faces or contacts).
///
/// This class provides functionality to manipulate a finley element file.
/// It is able to load element data from NetCDF files or retrieve it from
/// a Finley_ElementFile instance.
///
/// \note The corresponding mesh nodes are not part of this class but are
///       stored in a NodeData instance.
class ElementData
{
public:

    /// \brief Constructor with name and accompanying NodeData object.
    WEIPA_DLL_API
    ElementData(const std::string& elementName, NodeData_ptr nodes);

    /// \brief Copy constructor
    WEIPA_DLL_API
    ElementData(const ElementData& e);

    /// \brief Initialises with data from a Finley_ElementFile instance.
    WEIPA_DLL_API
    bool initFromFinley(const Finley_ElementFile* finleyFile);

    /// \brief Destructor
    WEIPA_DLL_API
    ~ElementData() {}

    /// \brief Reads element data from escript/finley NetCDF file.
    WEIPA_DLL_API
    bool readFromNc(NcFile* ncfile);

    /// \brief Moves "ghost" elements (whose owner does not match ownIndex)
    ///        and the corresponding data to the end of the arrays.
    WEIPA_DLL_API
    void reorderGhostZones(int ownIndex);

    /// \brief Removes "ghost" elements.
    WEIPA_DLL_API
    void removeGhostZones(int ownIndex);

    /// \brief Writes connectivity data to a stream in VTK text format.
    WEIPA_DLL_API
    void writeConnectivityVTK(std::ostream& os);

    /// \brief Writes element data into given directory in given Silo file.
    ///
    /// Since the mesh depends on element information this method also writes
    /// the node mesh itself. If Silo was not available at compile time or if
    /// a Silo function fails this method returns false.
    WEIPA_DLL_API
    bool writeToSilo(DBfile* dbfile, const std::string& siloPath);

    /// \brief Returns the names of the meshes associated with the elements.
    WEIPA_DLL_API
    StringVec getMeshNames() const;

    /// \brief Returns a vector with the mesh variable names.
    WEIPA_DLL_API
    StringVec getVarNames() const;

    /// \brief Returns the number of elements.
    WEIPA_DLL_API
    int getNumElements() const { return numElements; }

    /// \brief Returns the number of nodes per element.
    WEIPA_DLL_API
    int getNodesPerElement() const { return nodesPerElement; }

    /// \brief Returns the number of "ghost" elements.
    WEIPA_DLL_API
    int getGhostCount() const { return numGhostElements; }

    /// \brief Returns the type of the elements.
    WEIPA_DLL_API
    ZoneType getType() const { return type; }

    /// \brief Returns the original type id of the Finley reference elements.
    WEIPA_DLL_API
    ElementTypeId getFinleyTypeId() const { return finleyTypeId; }

    /// \brief Returns a vector of the node IDs used by the elements.
    WEIPA_DLL_API
    const IntVec& getNodeList() const { return nodes; }

    /// \brief Returns a vector of element IDs.
    WEIPA_DLL_API
    const IntVec& getIDs() const { return ID; }

    /// \brief Returns an array of data values for the name provided.
    ///
    /// The name must be one of the names returned from getVarNames().
    WEIPA_DLL_API
    const IntVec& getVarDataByName(const std::string varName) const;

    /// \brief Returns the node mesh instance used by the elements.
    WEIPA_DLL_API
    NodeData_ptr getNodeMesh() const { return nodeMesh; }

    /// \brief Returns the reduced elements.
    WEIPA_DLL_API
    ElementData_ptr getReducedElements() const { return reducedElements; }
 
    /// \brief Returns a QuadMaskInfo structure for given functionspace code.
    const QuadMaskInfo& getQuadMask(int functionSpace) const;
 
    /// \brief If the original element type is not supported they are
    ///        subdivided into N smaller elements (e.g. one Rec9 -> four Rec4)
    ///        and this method returns the multiplication factor N.
    int getElementFactor() const { return elementFactor; }

private:
    ElementData() {}
    FinleyElementInfo getFinleyTypeInfo(ElementTypeId typeId);
    void buildMeshes();
    void buildReducedElements(const FinleyElementInfo& f);
    IntVec prepareGhostIndices(int ownIndex);
    void reorderArray(IntVec& v, const IntVec& idx, int elementsPerIndex);
    QuadMaskInfo buildQuadMask(const CoordArray& quadNodes, int numQNodes);

    ElementData_ptr reducedElements;
    NodeData_ptr nodeMesh;
    NodeData_ptr originalMesh;
    std::string name;
    int numElements;
    int numGhostElements;
    int nodesPerElement;
    ZoneType type;
    ElementTypeId finleyTypeId;
    IntVec nodes;
    IntVec color, ID, tag;
    IntVec owner;
    QuadMaskInfo quadMask, reducedQuadMask;
    int elementFactor;
};

} // namespace weipa

#endif // __ELEMENTDATA_H__

