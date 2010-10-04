
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

#ifndef __WEIPA_FINLEYNODES_H__
#define __WEIPA_FINLEYNODES_H__

#include <weipa/NodeData.h>

class DBfile;
class NcFile;
struct Finley_NodeFile;

namespace weipa {

class FinleyNodes;
typedef boost::shared_ptr<FinleyNodes> FinleyNodes_ptr;

/// \brief Stores and manipulates finley mesh nodes.
///
/// This class provides functionality to manipulate a finley node file.
/// It is able to load node data from dump files or retrieve it from a
/// Finley_NodeFile instance.
class FinleyNodes : public NodeData
{
public:
    /// \brief Constructor with mesh name.
    FinleyNodes(const std::string& meshName);

    FinleyNodes(FinleyNodes_ptr fullNodes, IntVec& requiredNodes,
                const std::string& meshName);

    /// \brief Copy constructor.
    FinleyNodes(const FinleyNodes& m);

    /// \brief Virtual destructor
    virtual ~FinleyNodes();

    /// \brief Initialises with finley node file.
    bool initFromFinley(const Finley_NodeFile* finleyFile);

    /// \brief Reads node data from a NetCDF file.
    bool readFromNc(NcFile* ncFile);

    /// \brief Writes node data to a Silo file.
    bool writeToSilo(DBfile* dbfile);

    /// \brief Writes coordinates to a stream in VTK text format.
    virtual void writeCoordinatesVTK(std::ostream& os, int ownIndex);

    /// \brief Sets the silo path to be used when saving.
    void setSiloPath(const std::string& path) { siloPath = path; }

    /// \brief Returns an array of nodal data by the given name.
    ///
    /// The name must be one of the names returned by getVarNames().
    const IntVec& getVarDataByName(const std::string& name) const;

    /// \brief Returns a vector with the mesh variable names.
    virtual StringVec getVarNames() const;

    /// \brief Returns the name of this node mesh.
    virtual std::string getName() const { return name; }

    /// \brief Returns full Silo mesh name, e.g. "/block0000/Nodes".
    std::string getFullSiloName() const;

    /// \brief Returns the node ID array.
    virtual const IntVec& getNodeIDs() const { return nodeID; }

    /// \brief Returns the node distribution array
    virtual const IntVec& getNodeDistribution() const { return nodeDist; }

    /// \brief Returns the global node index array.
    virtual const IntVec& getGlobalNodeIndices() const { return nodeGNI; }

    /// \brief Returns the coordinates of the mesh nodes.
    virtual const CoordArray& getCoords() const { return coords; }

    /// \brief Returns the dimensionality of this mesh (2 or 3).
    virtual int getNumDims() const { return numDims; }

    /// \brief Returns the number of mesh nodes.
    virtual int getNumNodes() const { return numNodes; }

    /// \brief Returns the total number of mesh nodes for a distributed mesh.
    virtual int getGlobalNumNodes() const;

protected:
    CoordArray coords;     /// x, y[, z] coordinates of nodes
    int numDims;           /// dimensionality (2 or 3)
    int numNodes;          /// number of nodes
    IntVec nodeID;         /// node IDs
    IntVec nodeTag, nodeGDOF, nodeGNI, nodeGRDFI, nodeGRNI;
    IntVec nodeDist;       /// node distribution
    std::string name;      /// the name of this node mesh
    std::string siloPath;  /// the path to this mesh within the SILO file
};


inline std::string FinleyNodes::getFullSiloName() const
{
    std::string result(siloPath);
    if (result.length() == 0 || *result.rbegin() != '/')
        result += '/';
    result += name;
    return result;
}

} // namespace weipa

#endif // __WEIPA_FINLEYNODES_H__

