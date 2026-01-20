
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __FINLEY_NODEFILE_H__
#define __FINLEY_NODEFILE_H__

#define MAX_numDim 3

#include "Finley.h"
#include "NodeMapping.h"

#include <escript/Distribution.h>

#ifdef ESYS_HAVE_PASO
#include <paso/Coupler.h>
#endif
#ifdef ESYS_HAVE_TRILINOS
#include <trilinoswrap/types.h>
#endif

namespace escript {
    struct IndexList;
}

namespace finley {

class NodeFile
{
public:

    /// constructor - creates empty node file.
    /// Use allocTable() to allocate the node table (Id,Coordinates).
    NodeFile(int nDim, escript::JMPI MPIInfo);

    /// destructor
    ~NodeFile();

    /// allocates the node table within this node file to hold numNodes nodes.
    void allocTable(dim_t numNodes);

    /// empties the node table and frees all memory
    void freeTable();

    void print() const;

    inline index_t getFirstNode() const;
    inline index_t getLastNode() const;
    inline dim_t getGlobalNumNodes() const;
    inline const index_t* borrowGlobalNodesIndex() const;

    inline index_t getFirstReducedNode() const;
    inline index_t getLastReducedNode() const;
    inline index_t getGlobalNumReducedNodes() const;
    inline const index_t* borrowGlobalReducedNodesIndex() const;

    /// returns the number of FEM nodes (on this rank)
    inline dim_t getNumNodes() const;

    /// returns the number of reduced order FEM nodes (on this rank)
    inline dim_t getNumReducedNodes() const;

    /// returns the number of degrees of freedom (on this rank)
    inline dim_t getNumDegreesOfFreedom() const;

    /// returns the number of reduced order degrees of freedom (on this rank)
    inline dim_t getNumReducedDegreesOfFreedom() const;

    /// returns the number of degrees of freedom targets (own and shared)
    inline dim_t getNumDegreesOfFreedomTargets() const;

    /// returns the number of reduced degrees of freedom targets (own and shared)
    inline dim_t getNumReducedDegreesOfFreedomTargets() const;

    inline const IndexVector& borrowReducedNodesTarget() const;
    inline const IndexVector& borrowDegreesOfFreedomTarget() const;
    inline const IndexVector& borrowNodesTarget() const;
    inline const IndexVector& borrowReducedDegreesOfFreedomTarget() const;

    inline const index_t* borrowTargetReducedNodes() const;
    inline const index_t* borrowTargetDegreesOfFreedom() const;

    /// returns the mapping from local nodes to a target
    inline const index_t* borrowTargetNodes() const;
    inline const index_t* borrowTargetReducedDegreesOfFreedom() const;

    inline void updateTagList();

    /// creates a dense labeling of the global degrees of freedom and returns
    /// the new number of global degrees of freedom
    dim_t createDenseDOFLabeling();

    dim_t createDenseNodeLabeling(IndexVector& nodeDistribution,
                                  const IndexVector& dofDistribution);

    dim_t createDenseReducedLabeling(const std::vector<short>& reducedMask,
                                     bool useNodes);

    void createNodeMappings(const IndexVector& indexReducedNodes,
                            const IndexVector& dofDistribution,
                            const IndexVector& nodeDistribution);


    void assignMPIRankToDOFs(std::vector<int>& mpiRankOfDOF,
                             const IndexVector& distribution);

    void copyTable(index_t offset, index_t idOffset, index_t dofOffset,
                   const NodeFile* in);

    /// gathers nodes from the NodeFile `in` using the entries in
    /// index[0:numNodes-1] which are between min_index and max_index
    /// (exclusive)
    void gather(const index_t* index, const NodeFile* in);

    void gather_global(const index_t* index, const NodeFile* in);
    void scatter(const index_t* index, const NodeFile* in);

    void setCoordinates(const escript::Data& newX);

    /// set tags to newTag where mask > 0
    void setTags(int newTag, const escript::Data& mask);

    std::pair<index_t,index_t> getDOFRange() const;

private:
    std::pair<index_t,index_t> getGlobalIdRange() const;
    std::pair<index_t,index_t> getGlobalDOFRange() const;
    std::pair<index_t,index_t> getGlobalNodeIDIndexRange() const;
    dim_t prepareLabeling(const std::vector<short>& mask, IndexVector& buffer,
                          IndexVector& distribution, bool useNodes);
    void createDOFMappingAndCoupling(bool reduced);

    NodeMapping nodesMapping;
    NodeMapping degreesOfFreedomMapping;
    NodeMapping reducedDegreesOfFreedomMapping;

    /// number of nodes
    dim_t numNodes;

public:
    ///////////////////////////////////////
    // these should be private as well.

    NodeMapping reducedNodesMapping;

    /// MPI information
    escript::JMPI MPIInfo;
    /// number of spatial dimensions
    int numDim;
    /// Id[i] is the unique ID number of FEM node i
    index_t* Id;
    /// Tag[i] is the tag of node i
    int* Tag;
    /// vector of tags which are actually used
    std::vector<int> tagsInUse;

    /// globalDegreesOfFreedom[i] is the global degree of freedom assigned
    /// to node i. This index is used to consider periodic boundary conditions
    /// by assigning the same degree of freedom to different nodes.
    index_t* globalDegreesOfFreedom;
    /// Coordinates[INDEX2(k,i,numDim)] is the k-th coordinate of node i
    double* Coordinates;
    /// assigns each local node a global unique ID in a dense labeling of
    /// reduced DOF. Value <0 indicates that the DOF is not used.
    index_t* globalReducedDOFIndex;
    /// assigns each local reduced node a global unique ID in a dense labeling
    /// Value <0 indicates that the DOF is not used
    index_t* globalReducedNodesIndex;
    /// assigns each local node a global unique ID in a dense labeling
    index_t* globalNodesIndex;

    /// MPI distribution of nodes
    escript::Distribution_ptr nodesDistribution;
    escript::Distribution_ptr reducedNodesDistribution;
    escript::Distribution_ptr degreesOfFreedomDistribution;
    escript::Distribution_ptr reducedDegreesOfFreedomDistribution;

#ifdef ESYS_HAVE_PASO
    paso::Connector_ptr degreesOfFreedomConnector;
    paso::Connector_ptr reducedDegreesOfFreedomConnector;
#endif
#ifdef ESYS_HAVE_TRILINOS
    esys_trilinos::const_TrilinosMap_ptr trilinosRowMap;
    esys_trilinos::const_TrilinosMap_ptr trilinosReducedRowMap;
    esys_trilinos::const_TrilinosMap_ptr trilinosColMap;
    esys_trilinos::const_TrilinosMap_ptr trilinosReducedColMap;
#endif
  
    // these are the packed versions of Id
    index_t* reducedNodesId;        
    index_t* degreesOfFreedomId;
    index_t* reducedDegreesOfFreedomId;

    /// the status counts the updates done on the node coordinates.
    /// The value is increased by 1 when the node coordinates are updated.
    int status;
};

//
// implementation of inline methods
//

inline index_t NodeFile::getFirstNode() const
{
    return nodesDistribution->getFirstComponent();
}

inline index_t NodeFile::getLastNode() const
{
    return nodesDistribution->getLastComponent();
}

inline dim_t NodeFile::getGlobalNumNodes() const
{
    return nodesDistribution->getGlobalNumComponents();
}

inline const index_t* NodeFile::borrowGlobalNodesIndex() const
{
    return globalNodesIndex;
}

inline index_t NodeFile::getFirstReducedNode() const
{
    return reducedNodesDistribution->getFirstComponent();
}

inline index_t NodeFile::getLastReducedNode() const
{
    return reducedNodesDistribution->getLastComponent();
}

inline dim_t NodeFile::getGlobalNumReducedNodes() const
{
    return reducedNodesDistribution->getGlobalNumComponents();
}

inline const index_t* NodeFile::borrowGlobalReducedNodesIndex() const
{
    return globalReducedNodesIndex;
}

inline dim_t NodeFile::getNumNodes() const
{
    return numNodes;
}

inline dim_t NodeFile::getNumReducedNodes() const
{
    return reducedNodesMapping.getNumTargets();
}

inline dim_t NodeFile::getNumDegreesOfFreedom() const
{
    return degreesOfFreedomDistribution->getMyNumComponents();
}

inline dim_t NodeFile::getNumReducedDegreesOfFreedom() const
{
    return reducedDegreesOfFreedomDistribution->getMyNumComponents();
}

inline dim_t NodeFile::getNumDegreesOfFreedomTargets() const
{
    return degreesOfFreedomMapping.getNumTargets();
}

inline dim_t NodeFile::getNumReducedDegreesOfFreedomTargets() const
{
    return reducedDegreesOfFreedomMapping.getNumTargets();
}

inline const IndexVector& NodeFile::borrowNodesTarget() const
{
    return nodesMapping.map;
}

inline const IndexVector& NodeFile::borrowReducedNodesTarget() const
{
    return reducedNodesMapping.map;
}

inline const IndexVector& NodeFile::borrowDegreesOfFreedomTarget() const
{
    return degreesOfFreedomMapping.map;
}

inline const IndexVector& NodeFile::borrowReducedDegreesOfFreedomTarget() const
{
    return reducedDegreesOfFreedomMapping.map;
}

inline const index_t* NodeFile::borrowTargetNodes() const
{
    return &nodesMapping.target[0];
}

inline const index_t* NodeFile::borrowTargetReducedNodes() const
{
    return &reducedNodesMapping.target[0];
}

inline const index_t* NodeFile::borrowTargetDegreesOfFreedom() const
{
    return &degreesOfFreedomMapping.target[0];
}

inline const index_t* NodeFile::borrowTargetReducedDegreesOfFreedom() const
{
    return &reducedDegreesOfFreedomMapping.target[0];
}

inline void NodeFile::updateTagList()
{
    util::setValuesInUse(Tag, numNodes, tagsInUse, MPIInfo);
}


} // namespace finley

#endif // __FINLEY_NODEFILE_H__

