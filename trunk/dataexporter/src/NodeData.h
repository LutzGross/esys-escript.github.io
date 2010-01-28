
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

#ifndef __NODEDATA_H__
#define __NODEDATA_H__

#include <escriptexport/escriptexport.h>

class DBfile;
class NcFile;
struct Finley_NodeFile;

namespace escriptexport {

/// \brief
class NodeData
{
public:
    /// \brief Constructor with mesh name.
    ESCRIPTEXPORT_DLL_API
    NodeData(const std::string& meshName);
    
    ESCRIPTEXPORT_DLL_API
    NodeData(NodeData_ptr fullNodes, IntVec& requiredNodes,
             const std::string& meshName);

    /// \brief Copy constructor.
    ESCRIPTEXPORT_DLL_API
    NodeData(const NodeData& m);

    /// \brief Virtual destructor
    ESCRIPTEXPORT_DLL_API
    virtual ~NodeData();

    /// \brief Initialises with finley node file.
    ESCRIPTEXPORT_DLL_API
    bool initFromFinley(const Finley_NodeFile* finleyFile);

    /// \brief Reads node data from a NetCDF file.
    ESCRIPTEXPORT_DLL_API
    bool readFromNc(NcFile* ncFile);

    /// \brief Writes node data to a Silo file.
    ESCRIPTEXPORT_DLL_API
    bool writeToSilo(DBfile* dbfile);

    /// \brief Sets the silo path to be used when saving.
    ESCRIPTEXPORT_DLL_API
    void setSiloPath(const std::string& path) { siloPath = path; }

    /// \brief Returns an array of nodal data by the given name.
    ///
    /// The name must be one of the names returned by getVarNames().
    ESCRIPTEXPORT_DLL_API
    const IntVec& getVarDataByName(const std::string& name) const;

    /// \brief Returns a vector with the mesh variable names.
    ESCRIPTEXPORT_DLL_API
    StringVec getVarNames() const;

    /// \brief Returns the name of this node mesh.
    ESCRIPTEXPORT_DLL_API
    std::string getName() const { return name; }

    /// \brief Returns full Silo mesh name, e.g. "/block0000/Nodes".
    ESCRIPTEXPORT_DLL_API
    std::string getFullSiloName() const;
    
    /// \brief Returns the node ID array.
    ESCRIPTEXPORT_DLL_API
    const IntVec& getNodeIDs() const { return nodeID; }

    /// \brief Returns the coordinates of the mesh nodes.
    ESCRIPTEXPORT_DLL_API
    const CoordArray& getCoords() const { return coords; }

    /// \brief Returns the dimensionality of this mesh (2 or 3).
    ESCRIPTEXPORT_DLL_API
    int getNumDims() const { return numDims; }

    /// \brief Returns the number of mesh nodes.
    ESCRIPTEXPORT_DLL_API
    int getNumNodes() const { return numNodes; }

protected:
    CoordArray coords;         /// x, y[, z] coordinates of nodes
    int numDims;               /// dimensionality (2 or 3)
    int numNodes;              /// number of nodes
    IntVec nodeID;             /// node IDs
    IntVec nodeTag, nodeGDOF, nodeGNI, nodeGRDFI, nodeGRNI;
    IntVec nodeDist;           /// node distribution
    std::string name;          /// the name of this node mesh
    std::string siloPath;      /// the path to this mesh within the SILO file
};


inline std::string NodeData::getFullSiloName() const
{
    std::string result(siloPath);
    if (result.length() == 0 || *result.rbegin() != '/')
        result += '/';
    result += name;
    return result;
}

} // namespace escriptexport

#endif // __NODEDATA_H__

