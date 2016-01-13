
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

//
// ElementData.h
//
#ifndef __ELEMENTDATA_H__
#define __ELEMENTDATA_H__

#include <finley/ReferenceElements.h> // for ElementTypeId
#include <escriptreader/Mesh.h>

class DBfile;
class NcFile;

namespace EscriptReader {
    
typedef enum {
    ZONETYPE_BEAM=1,
    ZONETYPE_HEX,
    ZONETYPE_POLYGON,
    ZONETYPE_QUAD,
    ZONETYPE_TET,
    ZONETYPE_TRIANGLE
} ZoneType;

/// \brief Holds information that is used to convert from escript element types
/// to Silo.
struct FinleyElementInfo
{
    ZoneType elementType, reducedElementType;
    int elementFactor;
    int elementSize, reducedElementSize;
    const size_t* multiCellIndices;
};

/// \brief A class that stores and manipulates one type of escript mesh elements
/// (elements, faces, contacts or points).
///
/// The corresponding mesh nodes are not part of this class but stored in a Mesh
/// object which has to be provided when constructing an ElementData object.
class ElementData
{
    friend class DataVar;
    friend class MeshWithElements;
public:

    /// \brief Constructor with name and accompanying Mesh object.
    ElementData(const std::string& elementName, const Mesh* mainMesh);

    /// \brief Copy constructor
    ElementData(const ElementData& e);

    /// \brief Destructor
    ~ElementData();

    /// \brief Reads element data from escript NetCDF file.
    bool readFromNc(NcFile* ncfile);

    /// \brief Moves "ghost" elements and corresponding data to the end of the
    /// arrays.
    void handleGhostZones(int ownIndex);

    /// \brief Removes "ghost" elements.
    void removeGhostZones();

    /// \brief Writes element data into given directory in given Silo file.
    ///
    /// Since the mesh depends on element information this method also writes
    /// the mesh itself. If Silo was not available at compile time or if a Silo
    /// function fails this method returns false.
    bool writeToSilo(DBfile* dbfile, const std::string& siloPath);

    /// \brief Returns the names of the meshes associated with the elements.
    StringVec getMeshNames() const;

    /// \brief Returns a vector with the mesh variable names.
    StringVec getVarNames() const;

    /// \brief Returns the number of elements.
    int getCount() const { return count; }

    /// \brief Returns the number of reduced elements.
    int getReducedCount() const { return reducedCount; }

    /// \brief Returns the number of nodes per element.
    int getNodesPerElement() const { return nodesPerElement; }

    /// \brief Returns the number of nodes per reduced element.
    int getReducedNodesPerElement() const { return reducedNodesPerElement; }

    /// \brief Returns the number of "ghost" elements.
    int getGhostCount() const { return numGhostElements; }

    /// \brief Returns the number of "ghost" reduced elements.
    int getReducedGhostCount() const { return numReducedGhostElements; }

    /// \brief Returns the type of the elements.
    ZoneType getType() const { return type; }

    /// \brief Returns the type of reduced elements.
    ZoneType getReducedType() const { return reducedType; }

    /// \brief Returns a vector of the node IDs used by the elements.
    const IntVec& getNodeList() const { return nodes; }

    /// \brief Returns a vector of the node IDs used by the reduced elements.
    const IntVec& getReducedNodeList() const { return reducedNodes; }

    /// \brief Returns a vector of element IDs.
    const IntVec& getIDs() const { return ID; }

private:
    ElementData() {}
    FinleyElementInfo getFinleyTypeInfo(ElementTypeId typeId);
    void buildIndexMap();
    void buildMeshes();
    void buildReducedElements(const FinleyElementInfo& f);
    void prepareGhostIndices(int ownIndex);
    void reorderArray(IntVec& v, const IntVec& idx, int elementsPerIndex);

    const IntVec& getVarDataByName(const std::string varName) const;

    std::string name;
    int count, reducedCount;
    int numGhostElements, numReducedGhostElements;
    IndexMap ID2idx;
    IntVec indexArray, reducedIndexArray;
    Mesh* fullMesh;
    Mesh* reducedMesh;
    const Mesh* originalMesh;
    bool fullMeshIsOriginalMesh;

    int numDims;
    ZoneType type, reducedType;
    int nodesPerElement, reducedNodesPerElement;
    IntVec nodes, reducedNodes;
    IntVec color, ID, tag;
    IntVec owner, reducedOwner;
};

//
//
//
inline void ElementData::buildIndexMap()
{
    ID2idx.clear();
    size_t idx = 0;
    IntVec::const_iterator idIt;
    for (idIt = ID.begin(); idIt != ID.end(); idIt++, idx++)
        ID2idx[*idIt] = idx;
}

} // namespace EscriptReader

#endif // __ELEMENTDATA_H__

