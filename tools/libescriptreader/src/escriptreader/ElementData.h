
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

struct FinleyElementInfo
{
    int elementType, reducedElementType;
    int elementFactor;
    int elementSize, reducedElementSize;
    const size_t* multiCellIndices;
};

class DBfile;
class NcFile;

//
//
//
class ElementData
{
    friend class DataVar;
    friend class MeshWithElements;
public:
    ElementData(const std::string& elementName, const Mesh* mainMesh);

    /// Copy constructor
    ElementData(const ElementData& e);
    
    /// Virtual destructor
    virtual ~ElementData();
    
    bool readFromNc(NcFile* ncfile);
    void handleGhostZones(int ownIndex);
    void removeGhostZones();
    bool writeToSilo(DBfile* dbfile, const std::string& siloPath);

    StringVec getMeshNames() const;
    StringVec getVarNames() const;
    int getCount() const { return count; }
    int getReducedCount() const { return reducedCount; }
    int getNodesPerElement() const { return nodesPerElement; }
    int getReducedNodesPerElement() const { return reducedNodesPerElement; }
    int getGhostCount() const { return numGhostElements; }
    int getReducedGhostCount() const { return numReducedGhostElements; }
    int getType() const { return type; }
    int getReducedType() const { return reducedType; }
    const IntVec& getNodeList() const { return nodes; }
    const IntVec& getReducedNodeList() const { return reducedNodes; }
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
    int type, reducedType;
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


#endif // __ELEMENTDATA_H__

