
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#ifndef __DUDLEY_NODEFILE_H__
#define __DUDLEY_NODEFILE_H__

#include "Dudley.h"
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

namespace dudley {

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

    /// returns the number of FEM nodes (on this rank)
    inline dim_t getNumNodes() const;

    /// returns the number of degrees of freedom (on this rank)
    inline dim_t getNumDegreesOfFreedom() const;

    /// returns the number of degrees of freedom targets (own and shared)
    inline dim_t getNumDegreesOfFreedomTargets() const;

    /// returns the mapping from target to the local nodes
    inline const index_t* borrowNodesTarget() const;

    /// returns the mapping from target to the local degrees of freedom
    inline const index_t* borrowDegreesOfFreedomTarget() const;

    /// returns the mapping from local degrees of freedom to a target
    inline const index_t* borrowTargetDegreesOfFreedom() const;

    /// returns the mapping from local nodes to a target
    inline const index_t* borrowTargetNodes() const;

    inline void updateTagList();

    /// creates a dense labeling of the global degrees of freedom and returns
    /// the new number of global degrees of freedom
    dim_t createDenseDOFLabeling();

    dim_t createDenseNodeLabeling(IndexVector& nodeDistribution,
                                  const IndexVector& dofDistribution);

    void createNodeMappings(const IndexVector& dofDistribution,
                            const IndexVector& nodeDistribution);

    void assignMPIRankToDOFs(int* mpiRankOfDOF,
                             const IndexVector& distribution);

    void copyTable(index_t offset, index_t idOffset, index_t dofOffset,
                   const NodeFile* in);

    /// gathers nodes from the NodeFile `in` using the entries in
    /// index[0:numNodes-1] which are between min_index and max_index
    /// (exclusive)
    void gather(const index_t* index, const NodeFile* in);

    void gather_global(const index_t* index, const NodeFile* in);

    void setCoordinates(const escript::Data& newX);

    /// set tags to newTag where mask > 0
    void setTags(int newTag, const escript::Data& mask);

    std::pair<index_t,index_t> getDOFRange() const;

#ifdef ESYS_HAVE_TRILINOS
    void createTrilinosGraph(const escript::IndexList* indexList);
    esys_trilinos::const_TrilinosGraph_ptr getTrilinosGraph() const {
        return m_graph;
    }
#endif

private:
    std::pair<index_t,index_t> getGlobalIdRange() const;
    std::pair<index_t,index_t> getGlobalDOFRange() const;
    std::pair<index_t,index_t> getGlobalNodeIDIndexRange() const;
    dim_t prepareLabeling(const std::vector<short>& mask, IndexVector& buffer,
                          IndexVector& distribution, bool useNodes);
    void createDOFMappingAndCoupling();

    NodeMapping nodesMapping;
    NodeMapping degreesOfFreedomMapping;

    /// number of nodes
    dim_t numNodes;

#ifdef ESYS_HAVE_TRILINOS
    /// Trilinos graph structure, cached for efficiency
    esys_trilinos::const_TrilinosGraph_ptr m_graph;
#endif

public:
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
    /// assigns each local node a global unique ID in a dense labeling
    index_t* globalNodesIndex;

    /// MPI distribution of nodes
    escript::Distribution_ptr nodesDistribution;

    /// MPI distribution of degrees of freedom
    escript::Distribution_ptr dofDistribution;

#ifdef ESYS_HAVE_PASO
    paso::Connector_ptr degreesOfFreedomConnector;
#endif
    // these are the packed versions of Id
    index_t* degreesOfFreedomId;

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

inline dim_t NodeFile::getNumNodes() const
{
    return numNodes;
}

inline dim_t NodeFile::getNumDegreesOfFreedom() const
{
    return dofDistribution->getMyNumComponents();
}

inline dim_t NodeFile::getNumDegreesOfFreedomTargets() const
{
    return degreesOfFreedomMapping.numTargets;
}

inline const index_t* NodeFile::borrowNodesTarget() const
{
    return nodesMapping.map;
}

inline const index_t* NodeFile::borrowDegreesOfFreedomTarget() const
{
    return degreesOfFreedomMapping.map;
}

inline const index_t* NodeFile::borrowTargetNodes() const
{
    return nodesMapping.target;
}

inline const index_t* NodeFile::borrowTargetDegreesOfFreedom() const
{
    return degreesOfFreedomMapping.target;
}

inline void NodeFile::updateTagList()
{
    util::setValuesInUse(Tag, numNodes, tagsInUse, MPIInfo);
}


} // namespace dudley

#endif // __DUDLEY_NODEFILE_H__

