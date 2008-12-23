
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
// Mesh.h
//
#ifndef __MESH_H__
#define __MESH_H__

#include <escriptreader/common.h>

class DBfile;

namespace EscriptReader {

class ElementData;

//
//
//
class Mesh
{
    friend class ElementData;
public:
    /// Constructor with mesh coordinates, dimensionality and size
    Mesh(CoordArray c, int nDims, int nNodes);

    /// Copy constructor
    Mesh(const Mesh& m);

    /// Virtual destructor
    virtual ~Mesh();

    virtual bool readFromNc(const std::string& ncFile);
    virtual bool writeToSilo(DBfile* dbfile, const std::string& pathInSilo);

    void setName(const std::string& n) { name = n; }
    std::string getName() const { return name; }

    void setSiloPath(const std::string& p) { siloPath = p; }
    std::string getSiloPath() const { return siloPath; }

    /// Returns full Silo mesh name, e.g. "/block0000/Elements"
    std::string getFullSiloName() const;

    const IntVec& getNodeIDs() const { return nodeID; }
    void setNodeIDs(const IntVec& ids) { nodeID = ids; }

    const IndexMap& getIndexMap() const { return nodeID2idx; }
    void setIndexMap(const IndexMap& map) { nodeID2idx = map; }

    const CoordArray& getCoords() const { return coords; }
    int getNumDims() const { return numDims; }
    int getNumNodes() const { return numNodes; }

protected:
    /// Protected default constructor
    Mesh() {}
    
    void buildIndexMap();

    CoordArray coords;       // x, y[, z] coordinates of mesh nodes
    int numDims;             // dimensionality of mesh (2 or 3)
    int numNodes;            // number of mesh nodes
    IntVec nodeID;           // node IDs
    IndexMap nodeID2idx;     // mapping of node ID -> array index
    std::string name;        // the name of this mesh
    std::string siloPath;    // absolute path in Silo file
};


//
//
//
inline void Mesh::buildIndexMap()
{
    nodeID2idx.clear();
    int idx = 0;
    IntVec::const_iterator idIt;
    for (idIt = nodeID.begin(); idIt != nodeID.end(); idIt++, idx++)
        nodeID2idx[*idIt] = idx;
}

inline std::string Mesh::getFullSiloName() const
{
    if (siloPath == "/")
        return siloPath + name;
    else
        return siloPath + std::string("/") + name;
}

} // namespace EscriptReader

#endif // __MESH_H__

