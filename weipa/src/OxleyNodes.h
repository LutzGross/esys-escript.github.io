
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

#ifndef __WEIPA_OXLEYNODES_H__
#define __WEIPA_OXLEYNODES_H__

#include <weipa/NodeData.h>

class DBfile;

namespace oxley {
class OxleyDomain;
}

namespace weipa {

class OxleyNodes;
typedef boost::shared_ptr<OxleyNodes> OxleyNodes_ptr;

/// \brief Stores and manipulates oxley mesh nodes.
///
/// This class provides functionality to manipulate the nodes of a oxley
/// domain. It is able to load node data from dump files or retrieve it from
/// a OxleyDomain instance.
class OxleyNodes : public NodeData
{
public:
    /// \brief Constructor with mesh name.
    OxleyNodes(const std::string& meshName);

    OxleyNodes(OxleyNodes_ptr fullNodes, IntVec& requiredNodes,
                const std::string& meshName);

    /// \brief Copy constructor.
    OxleyNodes(const OxleyNodes& m);

    /// \brief Virtual destructor
    virtual ~OxleyNodes();

    /// \brief Initialises with oxley domain.
    bool initFromOxley(const oxley::OxleyDomain* oxleyDomain);

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
    virtual const IntVec& getGlobalNodeIndices() const { return nodeID; }

    /// \brief Returns the coordinates of the mesh nodes.
    virtual const CoordArray& getCoords() const { return coords; }

    /// \brief Returns the dimensionality of this mesh (2 or 3).
    virtual int getNumDims() const { return numDims; }

    /// \brief Returns the number of mesh nodes.
    virtual int getNumNodes() const { return numNodes; }

    /// \brief Returns the total number of mesh nodes for a distributed mesh.
    virtual int getGlobalNumNodes() const { return globalNumNodes; }

protected:
    CoordArray coords;     /// x, y[, z] coordinates of nodes
    int numDims;           /// dimensionality (2 or 3)
    int numNodes;          /// number of nodes
    int globalNumNodes;    /// global number of nodes
    IntVec nodeID;         /// node IDs
    IntVec nodeTag;        /// node tags
    IntVec nodeDist;       /// node distribution
    std::string name;      /// the name of this node mesh
    std::string siloPath;  /// the path to this mesh within the SILO file
};


inline std::string OxleyNodes::getFullSiloName() const
{
    std::string result(siloPath);
    if (result.length() == 0 || *result.rbegin() != '/')
        result += '/';
    result += name;
    return result;
}

} // namespace weipa

#endif // __WEIPA_OXLEYNODES_H__

