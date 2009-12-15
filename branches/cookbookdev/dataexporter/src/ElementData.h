
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

class DBfile;
class NcFile;

struct Finley_ElementFile;

namespace escriptexport {
    
class NodeData;

typedef enum {
    ZONETYPE_BEAM=1,
    ZONETYPE_HEX,
    ZONETYPE_POLYGON,
    ZONETYPE_QUAD,
    ZONETYPE_TET,
    ZONETYPE_TRIANGLE
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

    /// \brief Returns the number of reduced elements.
    ESCRIPTEXPORT_DLL_API
    int getReducedNumElements() const { return reducedNumElements; }

    /// \brief Returns the number of nodes per element.
    ESCRIPTEXPORT_DLL_API
    int getNodesPerElement() const { return nodesPerElement; }

    /// \brief Returns the number of nodes per reduced element.
    ESCRIPTEXPORT_DLL_API
    int getReducedNodesPerElement() const { return reducedNodesPerElement; }

    /// \brief Returns the number of "ghost" elements.
    ESCRIPTEXPORT_DLL_API
    int getGhostCount() const { return numGhostElements; }

    /// \brief Returns the number of "ghost" reduced elements.
    ESCRIPTEXPORT_DLL_API
    int getReducedGhostCount() const { return numReducedGhostElements; }

    /// \brief Returns the type of the elements.
    ESCRIPTEXPORT_DLL_API
    ZoneType getType() const { return type; }

    /// \brief Returns the type of reduced elements.
    ESCRIPTEXPORT_DLL_API
    ZoneType getReducedType() const { return reducedType; }

    /// \brief Returns a vector of the node IDs used by the elements.
    ESCRIPTEXPORT_DLL_API
    const IntVec& getNodeList() const { return nodes; }

    /// \brief Returns a vector of the node IDs used by the reduced elements.
    ESCRIPTEXPORT_DLL_API
    const IntVec& getReducedNodeList() const { return reducedNodes; }

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
    NodeData_ptr getNodeMesh() const { return fullMesh; }

    /// \brief Returns the node mesh instance used by the reduced elements.
    ESCRIPTEXPORT_DLL_API
    NodeData_ptr getReducedNodeMesh() const { return reducedMesh; }
 
private:
    ElementData() {}
    FinleyElementInfo getFinleyTypeInfo(ElementTypeId typeId);
    void buildMeshes();
    void buildReducedElements(const FinleyElementInfo& f);
    void prepareGhostIndices(int ownIndex, IntVec& indexArray,
                             IntVec& reducedIndexArray);
    void reorderArray(IntVec& v, const IntVec& idx, int elementsPerIndex);

    std::string name;
    int numElements, reducedNumElements;
    int numGhostElements, numReducedGhostElements;
    NodeData_ptr fullMesh;
    NodeData_ptr reducedMesh;
    NodeData_ptr originalNodes;
    bool fullMeshIsOriginalMesh;

    int numDims;
    ZoneType type, reducedType;
    int nodesPerElement, reducedNodesPerElement;
    IntVec nodes, reducedNodes;
    IntVec color, ID, tag;
    IntVec owner, reducedOwner;
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

