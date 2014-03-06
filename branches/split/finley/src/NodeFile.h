
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
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
#include <paso/Distribution.h>

struct Esys_MPIInfo;
struct Paso_Connector;

namespace finley {

class NodeFile
{
public:
    NodeFile(int nDim, Esys_MPIInfo *mpiInfo);
    ~NodeFile();

    void allocTable(int numNodes);
    void freeTable();

    void print() const;
    inline int getFirstNode() const;
    inline int getLastNode() const;
    inline int getGlobalNumNodes() const;
    inline int* borrowGlobalNodesIndex() const;

    inline int getFirstReducedNode() const;
    inline int getLastReducedNode() const;
    inline int getGlobalNumReducedNodes() const;
    inline int* borrowGlobalReducedNodesIndex() const;

    /// returns the number of FEM nodes
    inline int getNumNodes() const;
    inline int getNumReducedNodes() const;
    inline int getNumDegreesOfFreedom() const;
    inline int getNumReducedDegreesOfFreedom() const;

    inline const std::vector<int>& borrowReducedNodesTarget() const;
    inline const std::vector<int>& borrowDegreesOfFreedomTarget() const;
    inline const std::vector<int>& borrowNodesTarget() const;
    inline const std::vector<int>& borrowReducedDegreesOfFreedomTarget() const;

    inline const int* borrowTargetReducedNodes() const;
    inline const int* borrowTargetDegreesOfFreedom() const;
    inline const int* borrowTargetNodes() const;
    inline const int* borrowTargetReducedDegreesOfFreedom() const;

    void createNodeMappings(const std::vector<int>& indexReducedNodes,
                            const std::vector<int>& dofDistribution,
                            const std::vector<int>& nodeDistribution);
    int createDenseDOFLabeling();
    int createDenseNodeLabeling(std::vector<int>& nodeDistribution,
                                const std::vector<int>& dofDistribution);
    int createDenseReducedLabeling(const std::vector<short>& reducedMask,
                                   bool useNodes);
    void assignMPIRankToDOFs(std::vector<int>& mpiRankOfDOF, const std::vector<int>& distribution);

    void copyTable(int offset, int idOffset, int dofOffset, const NodeFile* in);
    void gather(int* index, const NodeFile* in);
    void gather_global(const std::vector<int>& index, const NodeFile* in);
    void scatter(int* index, const NodeFile* in);

    void setCoordinates(const escript::Data& newX);
    void setTags(const int newTag, const escript::Data& mask);
    inline void updateTagList();

    std::pair<int,int> getDOFRange() const;

private:
    std::pair<int,int> getGlobalIdRange() const;
    std::pair<int,int> getGlobalDOFRange() const;
    std::pair<int,int> getGlobalNodeIDIndexRange() const;
    int prepareLabeling(const std::vector<short>& mask,
                        std::vector<int>& buffer,
                        std::vector<int>& distribution, bool useNodes);
    void createDOFMappingAndCoupling(bool reduced);

    NodeMapping nodesMapping;
 
public:
    ///////////////////////////////////////
    // these should be private as well.

    NodeMapping reducedNodesMapping;
    NodeMapping degreesOfFreedomMapping;
    NodeMapping reducedDegreesOfFreedomMapping;

    /// MPI information
    Esys_MPIInfo *MPIInfo;
    /// number of nodes
    int numNodes;
    /// number of spatial dimensions
    int numDim;
    /// Id[i] is the id number of node i. It needs to be unique.
    int *Id;
    /// Tag[i] is the tag of node i.
    int *Tag;
    /// vector of tags which are actually used
    std::vector<int> tagsInUse;
    /// globalDegreesOfFreedom[i] is the global degree of freedom assigned
    /// to node i. This index is used to consider periodic boundary conditions
    /// by assigning the same degreesOfFreedom to the same node.
    int* globalDegreesOfFreedom;
    /// Coordinates[INDEX2(k,i,numDim)] is the k-th coordinate of node i
    double *Coordinates;
    /// assigns each local node a global unique Id in a dense labeling of
    /// reduced DOF. Value <0 indicates that the DOF is not used.
    int *globalReducedDOFIndex;
    /// assigns each local node a global unique Id in a dense labeling.
    /// Value <0 indicates that the DOF is not used
    int *globalReducedNodesIndex;
    /// assigns each local reduced node a global unique Id in a dense labeling
    int *globalNodesIndex;

    Paso_Distribution *nodesDistribution;
    Paso_Distribution *reducedNodesDistribution;
    Paso_Distribution *degreesOfFreedomDistribution;
    Paso_Distribution *reducedDegreesOfFreedomDistribution;

    Paso_Connector* degreesOfFreedomConnector;
    Paso_Connector *reducedDegreesOfFreedomConnector;
  
    /// these are the packed versions of Id
    int *reducedNodesId;        
    int *degreesOfFreedomId;
    int *reducedDegreesOfFreedomId;

    /// the status counts the updates done on the node coordinates.
    /// The value is increased by 1 when the node coordinates are updated.
    int status;
};

//
// implementation of inline methods
//

inline int NodeFile::getFirstNode() const
{
    return Paso_Distribution_getFirstComponent(nodesDistribution);
}

inline int NodeFile::getLastNode() const
{
    return Paso_Distribution_getLastComponent(nodesDistribution);
}

inline int NodeFile::getGlobalNumNodes() const
{
    return Paso_Distribution_getGlobalNumComponents(nodesDistribution);
}

inline int* NodeFile::borrowGlobalNodesIndex() const
{
    return globalNodesIndex;
}

inline int NodeFile::getFirstReducedNode() const
{
    return Paso_Distribution_getFirstComponent(reducedNodesDistribution);
}

inline int NodeFile::getLastReducedNode() const
{
    return Paso_Distribution_getLastComponent(reducedNodesDistribution);
}

inline int NodeFile::getGlobalNumReducedNodes() const
{
    return Paso_Distribution_getGlobalNumComponents(reducedNodesDistribution);
}

inline int* NodeFile::borrowGlobalReducedNodesIndex() const
{
    return globalReducedNodesIndex;
}

inline int NodeFile::getNumNodes() const
{
    return numNodes;
}

inline int NodeFile::getNumReducedNodes() const
{
    return reducedNodesMapping.getNumTargets();
}

inline int NodeFile::getNumDegreesOfFreedom() const
{
    return Paso_Distribution_getMyNumComponents(degreesOfFreedomDistribution);
}

inline int NodeFile::getNumReducedDegreesOfFreedom() const
{
    return Paso_Distribution_getMyNumComponents(reducedDegreesOfFreedomDistribution);
}

inline const std::vector<int>& NodeFile::borrowNodesTarget() const
{
    return nodesMapping.map;
}

inline const std::vector<int>& NodeFile::borrowReducedNodesTarget() const
{
    return reducedNodesMapping.map;
}

inline const std::vector<int>& NodeFile::borrowDegreesOfFreedomTarget() const
{
    return degreesOfFreedomMapping.map;
}

inline const std::vector<int>& NodeFile::borrowReducedDegreesOfFreedomTarget() const
{
    return reducedDegreesOfFreedomMapping.map;
}

inline const int* NodeFile::borrowTargetNodes() const
{
    return &nodesMapping.target[0];
}

inline const int* NodeFile::borrowTargetReducedNodes() const
{
    return &reducedNodesMapping.target[0];
}

inline const int* NodeFile::borrowTargetDegreesOfFreedom() const
{
    return &degreesOfFreedomMapping.target[0];
}

inline const int* NodeFile::borrowTargetReducedDegreesOfFreedom() const
{
    return &reducedDegreesOfFreedomMapping.target[0];
}

inline void NodeFile::updateTagList()
{
    util::setValuesInUse(Tag, numNodes, tagsInUse, MPIInfo);
}


} // namespace finley

#endif // __FINLEY_NODEFILE_H__

