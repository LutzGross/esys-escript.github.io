
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
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

#include <escriptexport/escriptexport.h>
#include <finley/ReferenceElements.h> // for ElementTypeId
#include <finley/vtkCellType.h>

class DBfile;
class NcFile;

struct Finley_ElementFile;

namespace escriptexport {
    
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
    ESCRIPTEXPORT_DLL_API
    ElementData(const std::string& elementName, NodeData_ptr nodes);

    /// \brief Copy constructor
    ESCRIPTEXPORT_DLL_API
    ElementData(const ElementData& e);

    /// \brief Initialises with data from a Finley_ElementFile instance.
    ESCRIPTEXPORT_DLL_API
    bool initFromFinley(const Finley_ElementFile* finleyFile);

    /// \brief Destructor
    ESCRIPTEXPORT_DLL_API
    ~ElementData();

    /// \brief Reads element data from escript/finley NetCDF file.
    ESCRIPTEXPORT_DLL_API
    bool readFromNc(NcFile* ncfile);

    /// \brief Moves "ghost" elements (whose owner does not match ownIndex)
    ///        and the corresponding data to the end of the arrays.
    ESCRIPTEXPORT_DLL_API
    void reorderGhostZones(int ownIndex);

    /// \brief Removes "ghost" elements.
    ESCRIPTEXPORT_DLL_API
    void removeGhostZones(int ownIndex);

    /// \brief Writes element data into given directory in given Silo file.
    ///
    /// Since the mesh depends on element information this method also writes
    /// the node mesh itself. If Silo was not available at compile time or if
    /// a Silo function fails this method returns false.
    ESCRIPTEXPORT_DLL_API
    bool writeToSilo(DBfile* dbfile, const std::string& siloPath);

    /// \brief Returns the names of the meshes associated with the elements.
    ESCRIPTEXPORT_DLL_API
    StringVec getMeshNames() const;

    /// \brief Returns a vector with the mesh variable names.
    ESCRIPTEXPORT_DLL_API
    StringVec getVarNames() const;

    /// \brief Returns the number of elements.
    ESCRIPTEXPORT_DLL_API
    int getNumElements() const { return numElements; }

    /// \brief Returns the number of nodes per element.
    ESCRIPTEXPORT_DLL_API
    int getNodesPerElement() const { return nodesPerElement; }

    /// \brief Returns the number of "ghost" elements.
    ESCRIPTEXPORT_DLL_API
    int getGhostCount() const { return numGhostElements; }

    /// \brief Returns the type of the elements.
    ESCRIPTEXPORT_DLL_API
    ZoneType getType() const { return type; }

    /// \brief Returns a vector of the node IDs used by the elements.
    ESCRIPTEXPORT_DLL_API
    const IntVec& getNodeList() const { return nodes; }

    /// \brief Returns a vector of element IDs.
    ESCRIPTEXPORT_DLL_API
    const IntVec& getIDs() const { return ID; }

    /// \brief Returns a mapping of element IDs to array index.
    ESCRIPTEXPORT_DLL_API
    IndexMap getIndexMap() const;

    /// \brief Returns an array of data values for the name provided.
    ///
    /// The name must be one of the names returned from getVarNames().
    ESCRIPTEXPORT_DLL_API
    const IntVec& getVarDataByName(const std::string varName) const;

    /// \brief Returns the node mesh instance used by the elements.
    ESCRIPTEXPORT_DLL_API
    NodeData_ptr getNodeMesh() const { return nodeMesh; }

    /// \brief Returns the reduced elements.
    ESCRIPTEXPORT_DLL_API
    ElementData_ptr getReducedElements() const { return reducedElements; }
 
private:
    ElementData() {}
    FinleyElementInfo getFinleyTypeInfo(ElementTypeId typeId);
    void buildMeshes();
    void buildReducedElements(const FinleyElementInfo& f);
    IntVec prepareGhostIndices(int ownIndex);
    void reorderArray(IntVec& v, const IntVec& idx, int elementsPerIndex);

    ElementData_ptr reducedElements;
    NodeData_ptr nodeMesh;
    NodeData_ptr originalMesh;
    std::string name;
    int numElements;
    int numGhostElements;
    int nodesPerElement;
    ZoneType type;
    IntVec nodes;
    IntVec color, ID, tag;
    IntVec owner;
};


inline IndexMap ElementData::getIndexMap() const
{
    IndexMap ID2idx;
    size_t idx = 0;
    IntVec::const_iterator idIt;
    for (idIt = ID.begin(); idIt != ID.end(); idIt++, idx++)
        ID2idx[*idIt] = idx;
    return ID2idx;
}

} // namespace escriptexport

#endif // __ELEMENTDATA_H__

