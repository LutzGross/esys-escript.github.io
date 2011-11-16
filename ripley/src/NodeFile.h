
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#ifndef __RIPLEY_NODEFILE_H__
#define __RIPLEY_NODEFILE_H__

#include <ripley/Ripley.h>
#include <ripley/RipleyException.h>
#include <ripley/Util.h>

#ifdef USE_NETCDF
class NcFile;
#endif

struct Paso_Connector;
struct Paso_Distribution;

namespace escript {
class Data;
}

namespace ripley {

class RipleyDomain;
class NodeFile;
typedef boost::shared_ptr<NodeFile> NodeFile_ptr;

struct NodeMapping
{
    index_t numTargets;
    IndexVector map;
    IndexVector target;

    NodeMapping(const IndexVector &tgt, index_t unused)
    {
        target = tgt;
        IndexPair minMaxTarget = getFlaggedMinMax(tgt, unused);
        if (minMaxTarget.first < 0)
            throw RipleyException("NodeMapping: Target has negative entry");
        numTargets = (minMaxTarget.first<minMaxTarget.second ? minMaxTarget.second+1 : 0);
        map.assign(numTargets, -1);
#pragma omp parallel for
        for (index_t i=0; i<tgt.size(); i++)
            if (target[i] != unused)
                map[target[i]] = i;
#pragma omp parallel for
        for (index_t i=0; i<numTargets; i++)
            if (map[i] == -1)
                throw RipleyException("NodeMapping: Target does not define a continuous labeling");
    }
};

typedef boost::shared_ptr<NodeMapping> NodeMapping_ptr;

/**
   \brief
   Ripley_NodeFile stores information about the nodes of a ripley domain.
*/

class NodeFile
{
public:
    /**
       \brief
       Constructor with dimensionality and MPI info structure.
    */
    RIPLEY_DLL_API
    NodeFile(dim_t numDim, Esys_MPIInfo *mpiInfo);

    /**
       \brief
       Destructor.
    */
    RIPLEY_DLL_API
    ~NodeFile();

#ifdef USE_NETCDF
    /**
       \brief
       Reads ID, coordinates, tags etc. from a netCDF file.
    */
    void readFromNetCDF(NcFile &dataFile, dim_t numNodes);

    /**
       \brief
       Stores ID, coordinates, tags etc. into a netCDF file.
    */
    void dumpToNetCDF(NcFile &dataFile);
#endif

    /**
       \brief
       Sets new node coordinates from an escript Data object.
    */
    void setCoordinates(const escript::Data &coords);

    /**
       \brief
       Returns a reference to the coordinates vector.
    */
    const std::vector<double> &getCoordinates() const { return m_coordinates; }

    /**
       \brief
       Returns the active MPI info structure.
    */
    Esys_MPIInfo* getMPIInfo() const { return m_mpiInfo; }

    /**
       \brief
       Returns a reference to the vector of node IDs.
    */
    const IndexVector &getIdVector() const { return m_id; }

    /**
       \brief
       Returns a reference to the vector of node tags.
    */
    const IndexVector &getTagVector() const { return m_tag; }

    /**
       \brief
       Returns a reference to the vector of global degrees of freedom.
    */
    const IndexVector &getGlobalDegreesOfFreedom() const { return m_globalDegreesOfFreedom; }

    /**
       \brief
       Returns a reference to the global nodes index vector.
    */
    const IndexVector &getGlobalNodesIndex() const { return m_globalNodesIndex; }

    /**
       \brief
       Returns a reference to the global reduced nodes index vector.
    */
    const IndexVector &getGlobalReducedNodesIndex() const { return m_globalReducedNodesIndex; }

    /**
       \brief
       Returns a reference to the global reduced DOF index vector.
    */
    const IndexVector &getGlobalReducedDOFIndex() const { return m_globalReducedDOFIndex; }

    /**
       \brief
       Returns a reference to the reduced nodes ID vector.
    */
    const IndexVector &getReducedNodesId() const { return m_reducedNodesId; }

    /**
       \brief
       Returns a reference to the degrees of freedom ID vector.
    */
    const IndexVector &getDegreesOfFreedomId() const { return m_degreesOfFreedomId; }

    /**
       \brief
       Returns a reference to the reduced DOF ID vector.
    */
    const IndexVector &getReducedDegreesOfFreedomId() const { return m_reducedDegreesOfFreedomId; }

    /**
       \brief
       Returns a reference to the vector of tags that are in use.
    */
    const IndexVector &getTagsInUse() const { return m_tagsInUse; }

    /**
       \brief
       Returns the reduced degrees of freedom mapping.
    */
    NodeMapping_ptr getReducedDegreesOfFreedomMapping() const { return m_reducedDegreesOfFreedomMapping; }

    /**
       \brief
       Returns the degrees of freedom mapping.
    */
    NodeMapping_ptr getDegreesOfFreedomMapping() const { return m_degreesOfFreedomMapping; }

    /**
       \brief
       Returns the reduced degrees of freedom connector.
    */
    Paso_Connector *getReducedDegreesOfFreedomConnector() const { return m_reducedDegreesOfFreedomConnector; }

    /**
       \brief
       Returns the degrees of freedom connector.
    */
    Paso_Connector *getDegreesOfFreedomConnector() const { return m_degreesOfFreedomConnector; }

    /**
       \brief
       Returns the node distribution.
    */
    Paso_Distribution *getNodesDistribution() const { return m_nodesDistribution; }

    /**
       \brief
       Returns the reduced degrees of freedom distribution.
    */
    Paso_Distribution *getReducedDegreesOfFreedomDistribution() const { return m_reducedDegreesOfFreedomDistribution; }

    /**
       \brief
       Returns the degrees of freedom distribution.
    */
    Paso_Distribution *getDegreesOfFreedomDistribution() const { return m_degreesOfFreedomDistribution; }

    /**
       \brief
       Returns the number of spatial dimensions.
    */
    dim_t getNumDim() const { return m_numDim; }

    /**
       \brief
       Returns the number of nodes in this node file.
    */
    dim_t getNumNodes() const { return m_id.size(); }

    /**
       \brief
       Returns the number of reduced nodes in this node file.
    */
    dim_t getNumReducedNodes() const { return m_reducedNodesMapping->numTargets; }

    /**
       \brief
       Returns the number of DOF in this node file.
    */
    dim_t getNumDegreesOfFreedom() const;

    /**
       \brief
       Returns the number of reduced DOF in this node file.
    */
    dim_t getNumReducedDegreesOfFreedom() const;

    /**
       \brief
       Returns the global number of nodes.
    */
    dim_t getGlobalNumNodes() const;

    /**
       \brief
       Returns the global number of reduced nodes.
    */
    dim_t getGlobalNumReducedNodes() const;

    /**
       \brief
       Returns the index of the first node for the current process.
    */
    index_t getFirstNode() const;

    /**
       \brief
       Returns the index of the last node for the current process.
    */
    index_t getLastNode() const;

    /**
       \brief
       Returns the index of the first reduced node for the current process.
    */
    index_t getFirstReducedNode() const;

    /**
       \brief
       Returns the index of the last reduced node for the current process.
    */
    index_t getLastReducedNode() const;

    /**
       \brief
       Copies node coordinates into an expanded Data object.
    */
    void assembleCoordinates(escript::Data &arg) const;

    /**
       \brief
       Returns the global index for node id.
    */
    index_t getIndexForGlobalNode(index_t id) const { return m_globalNodesIndex[id]; }

    /**
       \brief
       Returns the global index for reduced node id.
    */
    index_t getIndexForGlobalReducedNode(index_t id) const { return m_globalReducedNodesIndex[id]; }

    /**
       \brief
       Sets tags to newTag where mask>0.
    */
    void setTags(int newTag, const escript::Data &mask);

    /**
       \brief
       Refreshes the list of tags that are in use.
    */
    void updateTagsInUse();

    /**
       \brief
       Sets tags to newTag where mask>0.
    */
    dim_t getNumberOfTagsInUse() const { return m_tagsInUse.size(); }

    /**
       \brief
    */
    NodeFile_ptr gatherGlobal(IndexVector &index);

    /**
       \brief
       Gathers node data from the input vectors using entries in
       index[0:numNodes-1] which are between minIndex and maxIndex (exclusive).
    */
    void gatherEntries(const IndexVector &index, index_t minIndex,
                       index_t maxIndex, IndexVector &idIn, IndexVector &tagIn,
                       IndexVector &gDOFin, std::vector<double> &coordIn);

    /**
       \brief
       Scatters node data to the output vectors using entries in
       index[0:numNodes-1] which are between minIndex and maxIndex (exclusive).
    */
    void scatterEntries(const IndexVector &index, index_t minIndex,
                        index_t maxIndex, IndexVector &idOut,
                        IndexVector &tagOut, IndexVector &gDOFout,
                        std::vector<double> &coordOut);

    /**
       \brief
       Swaps the contents of the vectors with the node data vectors.
    */
    void swapEntries(IndexVector &idIn, IndexVector &tagIn,
                     IndexVector &gDOFin, std::vector<double> &coordIn);

    /**
       \brief
       Returns the minimum and maximum global degrees of freedom IDs as a pair.
    */
    IndexPair getGlobalDOFRange() const;

    /**
       \brief
       Returns the minimum and maximum global node indices as a pair.
    */
    IndexPair getGlobalNodeIndexRange() const;

    /**
       \brief
       Returns the minimum and maximum global node IDs as a pair.
    */
    IndexPair getGlobalIdRange() const;

    /**
       \brief
       Returns the minimum and maximum degrees of freedom IDs as a pair.
    */
    IndexPair getDOFRange() const;

    /**
       \brief
       Returns a vector which contains the processor id that owns each
       global degree of freedom.
    */
    RankVector getOwnerOfDOFs(const IndexVector &distribution);

    /**
       \brief
       Creates a dense labeling of the global degrees of freedom and
       returns the new number of global DoF.
    */
    dim_t createDenseDOFLabeling();

    /**
       \brief
       Creates a dense labeling of the reduced global degrees of freedom and
       returns the new number of reduced global DoF.
    */
    dim_t createDenseReducedDOFLabeling(const IndexVector &mask);

    /**
       \brief
       Creates a dense labeling of the global nodes and returns the new
       number of global nodes.
    */
    dim_t createDenseNodeLabeling(IndexVector &nodeDistribution,
                                  const IndexVector &dofDistribution);

    /**
       \brief
       Creates a dense labeling of the reduced global nodes and
       returns the new number of reduced global nodes.
    */
    dim_t createDenseReducedNodeLabeling(const IndexVector &mask);

    /**
       \brief
       Returns the status (=number of updates done on the coordinates).
    */
    int getStatus() const { return m_status; }

    /**
       \brief
       Assigns degrees of freedom with given distribution and IDs.
    */
    void resetGlobalDegreesOfFreedom(IndexVector &newGlobalDOFid,
                                     const IndexVector &distribution);
    /**
       \brief
       Helper
    */
    void createNodeFileMappings(const IndexVector &indexReducedNodes,
                                const IndexVector &dof_first_component,
                                const IndexVector &nodes_first_component);

private:
    /**
       \brief
       Helper
    */
    void createDOFMappingAndCoupling(bool useReduced);

    // MPI information
    Esys_MPIInfo *m_mpiInfo;

    // number of spatial dimensions
    dim_t m_numDim;

    // counts the updates done on the node coordinates. The value is
    // incremented when the node coordinates are updated.
    int m_status;

    // coordinates[INDEX2(k,i,numDim)] is the k-th coordinate of node i
    std::vector<double> m_coordinates;

    // m_id[i] is the unique ID of node i
    IndexVector m_id;

    // m_tag[i] is the tag of node i
    IndexVector m_tag;

    // array of tags which are actually used
    IndexVector m_tagsInUse;

    // gDOF[i] is the global degree of freedom assigned to node i.
    IndexVector m_globalDegreesOfFreedom;
    IndexVector m_globalReducedDOFIndex;
    IndexVector m_globalNodesIndex;
    IndexVector m_globalReducedNodesIndex;

    NodeMapping_ptr m_nodesMapping;
    NodeMapping_ptr m_reducedNodesMapping;
    NodeMapping_ptr m_degreesOfFreedomMapping;
    NodeMapping_ptr m_reducedDegreesOfFreedomMapping;

    Paso_Distribution *m_nodesDistribution;
    Paso_Distribution *m_reducedNodesDistribution;
    Paso_Distribution *m_degreesOfFreedomDistribution;
    Paso_Distribution *m_reducedDegreesOfFreedomDistribution;

    Paso_Connector *m_degreesOfFreedomConnector;
    Paso_Connector *m_reducedDegreesOfFreedomConnector;

    // packed versions of id's
    IndexVector m_reducedNodesId;
    IndexVector m_degreesOfFreedomId;
    IndexVector m_reducedDegreesOfFreedomId;

};

} // end of namespace ripley

#endif // __RIPLEY_NODEFILE_H__

