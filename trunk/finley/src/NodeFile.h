
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/

#ifndef __FINLEY_NODEFILE_H__
#define __FINLEY_NODEFILE_H__

#define MAX_numDim 3

#include "Finley.h"
#include "NodeMapping.h"
#include "paso/Distribution.h"
#include "paso/Coupler.h"
#include "esysUtils/Esys_MPI.h"

#include <vector>

namespace escript {
    class Data;
}

namespace finley {

class NodeFile
{
public:
    NodeFile(int nDim, Esys_MPIInfo *mpiInfo);
    ~NodeFile();

    void allocTable(int numNodes);
    void freeTable();

    inline int getFirstReducedNode();
    inline int getLastReducedNode();
    inline int getGlobalNumReducedNodes();
    inline int* borrowGlobalReducedNodesIndex();

    inline int getFirstNode();
    inline int getLastNode();
    inline int getGlobalNumNodes();
    inline int* borrowGlobalNodesIndex();

    inline int getNumReducedNodes();
    inline int getNumDegreesOfFreedom();
    inline int getNumNodes();
    inline int getNumReducedDegreesOfFreedom();

    inline int* borrowTargetReducedNodes();
    inline int* borrowTargetDegreesOfFreedom();
    inline int* borrowTargetNodes();
    inline int* borrowTargetReducedDegreesOfFreedom();

    inline int* borrowReducedNodesTarget();
    inline int* borrowDegreesOfFreedomTarget();
    inline int* borrowNodesTarget();
    inline int* borrowReducedDegreesOfFreedomTarget();

    int createDenseDOFLabeling();
    int createDenseNodeLabeling(int* node_distribution, const int* dof_distribution);
    int createDenseReducedLabeling(int* reducedMask, bool useNodes);
    void assignMPIRankToDOFs(int* mpiRankOfDOF, int *distribution);

    void copyTable(int offset, int idOffset, int dofOffset, const NodeFile* in);
    void gather(int* index, const NodeFile* in);
    void gather_global(int* index, const NodeFile* in);
    void scatter(int* index, const NodeFile* in);

    void setCoordinates(const escript::Data& newX);
    void setTags(const int newTag, const escript::Data& mask);
    void updateTagList();

    std::pair<int,int> getDOFRange() const;

private:
    std::pair<int,int> getGlobalIdRange() const;
    std::pair<int,int> getGlobalDOFRange() const;
    std::pair<int,int> getGlobalNodeIDIndexRange() const;
    int prepareLabeling(int* mask, std::vector<int>& buffer,
                        std::vector<int>& distribution, bool useNodes);

public:
    ///////////////////////////////////////

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

    Finley_NodeMapping *nodesMapping;
    Finley_NodeMapping *reducedNodesMapping;
    Finley_NodeMapping *degreesOfFreedomMapping;
    Finley_NodeMapping *reducedDegreesOfFreedomMapping;
 
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

inline int NodeFile::getFirstReducedNode()
{
    return Paso_Distribution_getFirstComponent(reducedNodesDistribution);
}

inline int NodeFile::getLastReducedNode()
{
    return Paso_Distribution_getLastComponent(reducedNodesDistribution);
}

inline int NodeFile::getGlobalNumReducedNodes()
{
    return Paso_Distribution_getGlobalNumComponents(reducedNodesDistribution);
}

inline int* NodeFile::borrowGlobalReducedNodesIndex()
{
    return globalReducedNodesIndex;
}

inline int NodeFile::getFirstNode()
{
    return Paso_Distribution_getFirstComponent(nodesDistribution);
}

inline int NodeFile::getLastNode()
{
    return Paso_Distribution_getLastComponent(nodesDistribution);
}

inline int NodeFile::getGlobalNumNodes()
{
    return Paso_Distribution_getGlobalNumComponents(nodesDistribution);
}

inline int* NodeFile::borrowGlobalNodesIndex()
{
    return globalNodesIndex;
}

inline int NodeFile::getNumReducedNodes()
{
    return reducedNodesMapping->numTargets;
}

inline int NodeFile::getNumDegreesOfFreedom()
{
    return Paso_Distribution_getMyNumComponents(degreesOfFreedomDistribution);
}

inline int NodeFile::getNumNodes()
{
    return nodesMapping->numNodes;
}

inline int NodeFile::getNumReducedDegreesOfFreedom()
{
    return Paso_Distribution_getMyNumComponents(reducedDegreesOfFreedomDistribution);
}

inline int* NodeFile::borrowTargetReducedNodes()
{
    return reducedNodesMapping->target;
}

inline int* NodeFile::borrowTargetDegreesOfFreedom()
{
    return degreesOfFreedomMapping->target;
}

inline int* NodeFile::borrowTargetNodes()
{
    return nodesMapping->target;
}

inline int* NodeFile::borrowTargetReducedDegreesOfFreedom()
{
    return reducedDegreesOfFreedomMapping->target;
}

inline int* NodeFile::borrowReducedNodesTarget()
{
    return reducedNodesMapping->map;
}

inline int* NodeFile::borrowDegreesOfFreedomTarget()
{
    return degreesOfFreedomMapping->map;
}

inline int* NodeFile::borrowNodesTarget()
{
    return nodesMapping->map;
}

inline int* NodeFile::borrowReducedDegreesOfFreedomTarget()
{
    return reducedDegreesOfFreedomMapping->map;
}


} // namespace finley

#endif // __FINLEY_NODEFILE_H__

