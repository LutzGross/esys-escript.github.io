
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

#ifndef __FINLEY_NODEFILE_H__
#define __FINLEY_NODEFILE_H__

#define MAX_numDim 3

#include "Finley.h"
#include "NodeMapping.h"

#include <paso/Coupler.h>
#include <paso/Distribution.h>


namespace finley {

class NodeFile
{
public:
    NodeFile(int nDim, escript::JMPI& mpiInfo);
    ~NodeFile();

    void allocTable(dim_t numNodes);
    void freeTable();

    void print() const;
    inline index_t getFirstNode() const;
    inline index_t getLastNode() const;
    inline index_t getGlobalNumNodes() const;
    inline index_t* borrowGlobalNodesIndex() const;

    inline index_t getFirstReducedNode() const;
    inline index_t getLastReducedNode() const;
    inline index_t getGlobalNumReducedNodes() const;
    inline index_t* borrowGlobalReducedNodesIndex() const;

    /// returns the number of FEM nodes
    inline dim_t getNumNodes() const;
    inline dim_t getNumReducedNodes() const;
    inline dim_t getNumDegreesOfFreedom() const;
    inline dim_t getNumReducedDegreesOfFreedom() const;

    inline const std::vector<index_t>& borrowReducedNodesTarget() const;
    inline const std::vector<index_t>& borrowDegreesOfFreedomTarget() const;
    inline const std::vector<index_t>& borrowNodesTarget() const;
    inline const std::vector<index_t>& borrowReducedDegreesOfFreedomTarget() const;

    inline const index_t* borrowTargetReducedNodes() const;
    inline const index_t* borrowTargetDegreesOfFreedom() const;
    inline const index_t* borrowTargetNodes() const;
    inline const index_t* borrowTargetReducedDegreesOfFreedom() const;

    void createNodeMappings(const std::vector<index_t>& indexReducedNodes,
                            const std::vector<index_t>& dofDistribution,
                            const std::vector<index_t>& nodeDistribution);
    dim_t createDenseDOFLabeling();
    dim_t createDenseNodeLabeling(std::vector<index_t>& nodeDistribution,
                                  const std::vector<index_t>& dofDistribution);
    dim_t createDenseReducedLabeling(const std::vector<short>& reducedMask,
                                     bool useNodes);
    void assignMPIRankToDOFs(std::vector<int>& mpiRankOfDOF, const std::vector<index_t>& distribution);

    void copyTable(index_t offset, index_t idOffset, index_t dofOffset,
                   const NodeFile* in);
    void gather(const index_t* index, const NodeFile* in);
    void gather_global(const index_t* index, const NodeFile* in);
    void scatter(const index_t* index, const NodeFile* in);

    void setCoordinates(const escript::Data& newX);
    void setTags(const int newTag, const escript::Data& mask);
    inline void updateTagList();

    std::pair<index_t,index_t> getDOFRange() const;

private:
    std::pair<index_t,index_t> getGlobalIdRange() const;
    std::pair<index_t,index_t> getGlobalDOFRange() const;
    std::pair<index_t,index_t> getGlobalNodeIDIndexRange() const;
    dim_t prepareLabeling(const std::vector<short>& mask,
                          std::vector<index_t>& buffer,
                          std::vector<index_t>& distribution, bool useNodes);
    void createDOFMappingAndCoupling(bool reduced);

    NodeMapping nodesMapping;
 
public:
    ///////////////////////////////////////
    // these should be private as well.

    NodeMapping reducedNodesMapping;
    NodeMapping degreesOfFreedomMapping;
    NodeMapping reducedDegreesOfFreedomMapping;

    /// MPI information
    escript::JMPI MPIInfo;
    /// number of nodes
    dim_t numNodes;
    /// number of spatial dimensions
    int numDim;
    /// Id[i] is the id number of node i. It needs to be unique.
    index_t *Id;
    /// Tag[i] is the tag of node i.
    int *Tag;
    /// vector of tags which are actually used
    std::vector<int> tagsInUse;
    /// globalDegreesOfFreedom[i] is the global degree of freedom assigned
    /// to node i. This index is used to consider periodic boundary conditions
    /// by assigning the same degreesOfFreedom to the same node.
    index_t* globalDegreesOfFreedom;
    /// Coordinates[INDEX2(k,i,numDim)] is the k-th coordinate of node i
    double *Coordinates;
    /// assigns each local node a global unique Id in a dense labeling of
    /// reduced DOF. Value <0 indicates that the DOF is not used.
    index_t *globalReducedDOFIndex;
    /// assigns each local node a global unique Id in a dense labeling.
    /// Value <0 indicates that the DOF is not used
    index_t *globalReducedNodesIndex;
    /// assigns each local reduced node a global unique Id in a dense labeling
    index_t *globalNodesIndex;

    paso::Distribution_ptr nodesDistribution;
    paso::Distribution_ptr reducedNodesDistribution;
    paso::Distribution_ptr degreesOfFreedomDistribution;
    paso::Distribution_ptr reducedDegreesOfFreedomDistribution;

    paso::Connector_ptr degreesOfFreedomConnector;
    paso::Connector_ptr reducedDegreesOfFreedomConnector;
  
    /// these are the packed versions of Id
    index_t *reducedNodesId;        
    index_t *degreesOfFreedomId;
    index_t *reducedDegreesOfFreedomId;

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

inline index_t NodeFile::getGlobalNumNodes() const
{
    return nodesDistribution->getGlobalNumComponents();
}

inline index_t* NodeFile::borrowGlobalNodesIndex() const
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

inline index_t* NodeFile::borrowGlobalReducedNodesIndex() const
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

inline const std::vector<index_t>& NodeFile::borrowNodesTarget() const
{
    return nodesMapping.map;
}

inline const std::vector<index_t>& NodeFile::borrowReducedNodesTarget() const
{
    return reducedNodesMapping.map;
}

inline const std::vector<index_t>& NodeFile::borrowDegreesOfFreedomTarget() const
{
    return degreesOfFreedomMapping.map;
}

inline const std::vector<index_t>& NodeFile::borrowReducedDegreesOfFreedomTarget() const
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

