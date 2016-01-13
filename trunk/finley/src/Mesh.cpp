/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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


/****************************************************************************

  Finley: Mesh

*****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"


#include "Mesh.h"
#include "IndexList.h"
#include <boost/scoped_array.hpp>

#include "CPPAdapter/FinleyAdapterException.h"

namespace finley {

/// Constructor.
/// Allocates a Mesh with given name and dimensionality
Mesh::Mesh(const std::string name, int numDim, esysUtils::JMPI& mpi_info) :
    m_name(name),
    approximationOrder(-1),
    reducedApproximationOrder(-1),
    integrationOrder(-1),
    reducedIntegrationOrder(-1),
    Elements(NULL),
    FaceElements(NULL),
    ContactElements(NULL),
    Points(NULL)
{
    MPIInfo = mpi_info;

    // allocate node table
    Nodes = new NodeFile(numDim, mpi_info);
}

/// destructor
Mesh::~Mesh()
{
    delete Nodes;
    delete FaceElements;
    delete Elements;
    delete ContactElements;
    delete Points;
    tagMap.clear();
}

void Mesh::setElements(ElementFile *elements)
{
    delete Elements;
    Elements=elements;
}

void Mesh::setFaceElements(ElementFile *elements)
{
    delete FaceElements;
    FaceElements=elements;
}

void Mesh::setContactElements(ElementFile *elements)
{
    delete ContactElements;
    ContactElements=elements;
}

void Mesh::setPoints(ElementFile *elements)
{
    delete Points;
    Points=elements;
}

void Mesh::setOrders() 
{
    const int ORDER_MAX=9999999;
    int locals[4] = { ORDER_MAX, ORDER_MAX, ORDER_MAX, ORDER_MAX };

    if (Elements != NULL && Elements->numElements > 0) {
        locals[0]=std::min(locals[0], Elements->referenceElementSet->referenceElement->BasisFunctions->Type->numOrder);
        locals[1]=std::min(locals[1], Elements->referenceElementSet->referenceElement->LinearBasisFunctions->Type->numOrder);
        locals[2]=std::min(locals[2], Elements->referenceElementSet->referenceElement->integrationOrder);
        locals[3]=std::min(locals[3], Elements->referenceElementSet->referenceElementReducedQuadrature->integrationOrder);
    }
    if (FaceElements != NULL && FaceElements->numElements > 0) {
        locals[0]=std::min(locals[0], FaceElements->referenceElementSet->referenceElement->BasisFunctions->Type->numOrder);
        locals[1]=std::min(locals[1], FaceElements->referenceElementSet->referenceElement->LinearBasisFunctions->Type->numOrder);
        locals[2]=std::min(locals[2], FaceElements->referenceElementSet->referenceElement->integrationOrder);
        locals[3]=std::min(locals[3], FaceElements->referenceElementSet->referenceElementReducedQuadrature->integrationOrder);
    }
    if (ContactElements != NULL && ContactElements->numElements > 0) {
        locals[0]=std::min(locals[0], ContactElements->referenceElementSet->referenceElement->BasisFunctions->Type->numOrder);
        locals[1]=std::min(locals[1], ContactElements->referenceElementSet->referenceElement->LinearBasisFunctions->Type->numOrder);
        locals[2]=std::min(locals[2], ContactElements->referenceElementSet->referenceElement->integrationOrder);
        locals[3]=std::min(locals[3], ContactElements->referenceElementSet->referenceElementReducedQuadrature->integrationOrder);
    }

#ifdef ESYS_MPI
    int globals[4];
    MPI_Allreduce(locals, globals, 4, MPI_INT, MPI_MIN, MPIInfo->comm);
    approximationOrder=(globals[0] < ORDER_MAX ? globals[0] : -1);
    reducedApproximationOrder=(globals[1] < ORDER_MAX ? globals[1] : -1);
    integrationOrder=(globals[2] < ORDER_MAX ? globals[2] : -1);
    reducedIntegrationOrder=(globals[3] < ORDER_MAX ? globals[3] : -1);
#else
    approximationOrder=(locals[0] < ORDER_MAX ? locals[0] : -1);
    reducedApproximationOrder=(locals[1] < ORDER_MAX ? locals[1] : -1);
    integrationOrder=(locals[2] < ORDER_MAX ? locals[2] : -1);
    reducedIntegrationOrder=(locals[3] < ORDER_MAX ? locals[3] : -1);
#endif
}

/// creates node mappings without (re-)distributing anything
void Mesh::createMappings(const std::vector<index_t>& dofDistribution,
                          const std::vector<index_t>& nodeDistribution)
{
    std::vector<short> maskReducedNodes(Nodes->numNodes, -1);
    markNodes(maskReducedNodes, 0, true);
    std::vector<index_t> indexReducedNodes = util::packMask(maskReducedNodes);
    if (noError())
        Nodes->createNodeMappings(indexReducedNodes, dofDistribution,
                                  nodeDistribution);
}

/// redistributes the Nodes and Elements including overlap
/// according to the DOF distribution. It will create an element colouring
/// but will not create any mappings.
void Mesh::distributeByRankOfDOF(const std::vector<index_t>& dof_distribution)
{
    std::vector<int> mpiRankOfDOF(Nodes->numNodes);
    Nodes->assignMPIRankToDOFs(mpiRankOfDOF, dof_distribution);

    // first, the elements are redistributed according to mpiRankOfDOF
    // at the input the Node tables refer to the local labeling of the nodes
    // while at the output they refer to the global labeling which is rectified
    // in the next step
    if (noError())
        Elements->distributeByRankOfDOF(mpiRankOfDOF, Nodes->Id);
    if (noError())
        FaceElements->distributeByRankOfDOF(mpiRankOfDOF, Nodes->Id);
    if (noError())
        ContactElements->distributeByRankOfDOF(mpiRankOfDOF, Nodes->Id);
    if (noError())
        Points->distributeByRankOfDOF(mpiRankOfDOF, Nodes->Id);

    // resolve the node ids
    if (noError())
        resolveNodeIds();

    // create a local labeling of the DOFs
    const std::pair<index_t,index_t> dof_range(Nodes->getDOFRange());
    const index_t len=dof_range.second-dof_range.first+1;
    // local mask for used nodes
    std::vector<index_t> localDOF_mask(len, -1);
    std::vector<index_t> localDOF_map(Nodes->numNodes, -1);

#pragma omp parallel for
    for (index_t n=0; n<Nodes->numNodes; n++) {
#ifdef BOUNDS_CHECK
        if ((Nodes->globalDegreesOfFreedom[n]-dof_range.first) >= len ||
                (Nodes->globalDegreesOfFreedom[n]-dof_range.first) < 0) {
            printf("BOUNDS_CHECK %s %d\n", __FILE__, __LINE__);
            exit(1);
        }
#endif
        localDOF_mask[Nodes->globalDegreesOfFreedom[n]-dof_range.first]=n;
    }

    index_t numDOFs=0;
    for (int n=0; n<len; n++) {
        const index_t k=localDOF_mask[n];
        if (k>=0) {
             localDOF_mask[n]=numDOFs;
             numDOFs++;
          }
    }
#pragma omp parallel for
    for (index_t n=0; n<Nodes->numNodes; n++) {
        const index_t k=localDOF_mask[Nodes->globalDegreesOfFreedom[n]-dof_range.first];
        localDOF_map[n]=k;
    }
    // create element coloring
    if (noError())
        createColoring(localDOF_map);
}

/// prints the mesh details to standard output
void Mesh::print()
{
    // write header
    printf("Mesh name: %s\n", m_name.c_str());
  
    // write nodes
    Nodes->print();
  
    // write elements
    if (Elements) {
        std::cout << "=== "
                 << Elements->referenceElementSet->referenceElement->Type->Name
                 << ":\nnumber of elements=" << Elements->numElements
                 << "\ncolor range=[" << Elements->minColor << ","
                 << Elements->maxColor << "]\n";
        if (Elements->numElements > 0) {
            const int NN=Elements->referenceElementSet->referenceElement->Type->numNodes;
            const int NN2=Elements->numNodes;
            std::cout << "Id,Tag,Owner,Color,Nodes" << std::endl;
            for (index_t i=0; i<Elements->numElements; i++) {
                std::cout << Elements->Id[i] << "," << Elements->Tag[i] << ","
                    << Elements->Owner[i] << "," << Elements->Color[i] << ",";
                for (int j=0; j<NN; j++)
                    std::cout << " " << Nodes->Id[Elements->Nodes[INDEX2(j,i,NN2)]];
                std::cout << std::endl;
            }
        }
    }

    // write face elements
    if (FaceElements) {
        std::cout << "=== "
                 << FaceElements->referenceElementSet->referenceElement->Type->Name
                 << ":\nnumber of elements=" << FaceElements->numElements
                 << "\ncolor range=[" << FaceElements->minColor << ","
                 << FaceElements->maxColor << "]\n";
        if (FaceElements->numElements > 0) {
            const int NN=FaceElements->referenceElementSet->referenceElement->Type->numNodes;
            const int NN2=FaceElements->numNodes;
            std::cout << "Id,Tag,Owner,Color,Nodes" << std::endl;
            for (index_t i=0; i<FaceElements->numElements; i++) {
                std::cout << FaceElements->Id[i] << "," << FaceElements->Tag[i]
                    << "," << FaceElements->Owner[i] << ","
                    << FaceElements->Color[i] << ",";
                for (int j=0; j<NN; j++)
                    std::cout << " " << Nodes->Id[FaceElements->Nodes[INDEX2(j,i,NN2)]];
                std::cout << std::endl;
            }
        }
    }

    // write Contact elements
    if (ContactElements) {
        std::cout << "=== "
                 << ContactElements->referenceElementSet->referenceElement->Type->Name
                 << ":\nnumber of elements=" << ContactElements->numElements
                 << "\ncolor range=[" << ContactElements->minColor << ","
                 << ContactElements->maxColor << "]\n";
        if (ContactElements->numElements > 0) {
            const int NN=ContactElements->referenceElementSet->referenceElement->Type->numNodes;
            const int NN2=ContactElements->numNodes;
            std::cout << "Id,Tag,Owner,Color,Nodes" << std::endl;
            for (index_t i=0; i<ContactElements->numElements; i++) {
                std::cout << ContactElements->Id[i] << ","
                    << ContactElements->Tag[i] << ","
                    << ContactElements->Owner[i] << ","
                    << ContactElements->Color[i] << ",";
                for (int j=0; j<NN; j++)
                    std::cout << " " << Nodes->Id[ContactElements->Nodes[INDEX2(j,i,NN2)]];
                std::cout << std::endl;
            }
        }
    }
  
    // write points
    if (Points) {
        std::cout << "=== "
                 << Points->referenceElementSet->referenceElement->Type->Name
                 << ":\nnumber of elements=" << Points->numElements
                 << "\ncolor range=[" << Points->minColor << ","
                 << Points->maxColor << "]\n";
        if (Points->numElements > 0) {
            const int NN=Points->referenceElementSet->referenceElement->Type->numNodes;
            const int NN2=Points->numNodes;
            std::cout << "Id,Tag,Owner,Color,Nodes" << std::endl;
            for (index_t i=0; i<Points->numElements; i++) {
                std::cout << Points->Id[i] << "," << Points->Tag[i] << ","
                    << Points->Owner[i] << "," << Points->Color[i] << ",";
                for (int j=0; j<NN; j++)
                    std::cout << " " << Nodes->Id[Points->Nodes[INDEX2(j,i,NN2)]];
                std::cout << std::endl;
            }
        }
    }
}

void Mesh::markNodes(std::vector<short>& mask, int offset, bool useLinear)
{
    Elements->markNodes(mask, offset, useLinear);
    FaceElements->markNodes(mask, offset, useLinear);
    ContactElements->markNodes(mask, offset, useLinear);
    Points->markNodes(mask, offset, useLinear);
}

void Mesh::markDOFsConnectedToRange(int* mask, int offset, int marker, 
                                    index_t firstDOF, index_t lastDOF,
                                    bool useLinear)
{
    const index_t *dofIndex = (useLinear ? Nodes->globalReducedDOFIndex
                                     : Nodes->globalDegreesOfFreedom);
    Elements->markDOFsConnectedToRange(mask, offset, marker, firstDOF, lastDOF,
            dofIndex, useLinear);
    FaceElements->markDOFsConnectedToRange(mask, offset, marker, firstDOF,
            lastDOF, dofIndex, useLinear);
    ContactElements->markDOFsConnectedToRange(mask, offset, marker, firstDOF,
            lastDOF, dofIndex, useLinear);
    Points->markDOFsConnectedToRange(mask, offset, marker, firstDOF, lastDOF,
            dofIndex, useLinear);
}

/// optimizes the labeling of the DOFs on each processor
void Mesh::optimizeDOFLabeling(const std::vector<index_t>& distribution)
{
    const int myRank=MPIInfo->rank;
    const int mpiSize=MPIInfo->size;
    const index_t myFirstVertex=distribution[myRank];
    const index_t myLastVertex=distribution[myRank+1];
    const dim_t myNumVertices=myLastVertex-myFirstVertex;
    index_t len=0;
    for (int p=0; p<mpiSize; ++p)
        len=std::max(len, distribution[p+1]-distribution[p]);

    boost::scoped_array<IndexList> index_list(new IndexList[myNumVertices]);
    std::vector<index_t> newGlobalDOFID(len);
    // create the adjacency structure xadj and adjncy
#pragma omp parallel
    {
        // insert contributions from element matrices into columns index
        IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                myFirstVertex, myLastVertex, Elements,
                Nodes->globalDegreesOfFreedom,
                Nodes->globalDegreesOfFreedom);
        IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                myFirstVertex, myLastVertex, FaceElements,
                Nodes->globalDegreesOfFreedom,
                Nodes->globalDegreesOfFreedom);
        IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                myFirstVertex, myLastVertex, ContactElements,
                Nodes->globalDegreesOfFreedom,
                Nodes->globalDegreesOfFreedom);
        IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                myFirstVertex, myLastVertex, Points,
                Nodes->globalDegreesOfFreedom,
                Nodes->globalDegreesOfFreedom);
    }
    // create the local matrix pattern
    paso::Pattern_ptr pattern=paso::Pattern::fromIndexListArray(0,
            myNumVertices, index_list.get(), myFirstVertex, myLastVertex,
            -myFirstVertex);

    if (noError())
        pattern->reduceBandwidth(&newGlobalDOFID[0]); 

    esysUtils::Esys_MPIInfo_noError(MPIInfo);

    if (noError()) {
        // shift new labeling to create a global id
#pragma omp parallel for
        for (int i=0; i<myNumVertices; ++i)
            newGlobalDOFID[i]+=myFirstVertex;

        // distribute new labeling to other processors
#ifdef ESYS_MPI
        const int dest=esysUtils::mod_rank(mpiSize, myRank + 1);
        const int source=esysUtils::mod_rank(mpiSize, myRank - 1);
#endif
        int current_rank=myRank;
        for (int p=0; p<mpiSize; ++p) {
            const index_t firstVertex=distribution[current_rank];
            const index_t lastVertex=distribution[current_rank+1];
#pragma omp parallel for
            for (index_t i=0; i<Nodes->numNodes; ++i) {
                const index_t k=Nodes->globalDegreesOfFreedom[i];
                if (firstVertex<=k && k<lastVertex) {
                    Nodes->globalDegreesOfFreedom[i]=newGlobalDOFID[k-firstVertex];
                }
            }
   
            if (p<mpiSize-1) { // the final send can be skipped
#ifdef ESYS_MPI
                MPI_Status status;
                MPI_Sendrecv_replace(&newGlobalDOFID[0], len, MPI_DIM_T,
                                     dest, MPIInfo->msg_tag_counter,
                                     source, MPIInfo->msg_tag_counter,
                                     MPIInfo->comm, &status);
#endif
                MPIInfo->msg_tag_counter++;
                current_rank=esysUtils::mod_rank(mpiSize, current_rank-1);
            }
        }
    }
}

/// prepares the mesh for further use
void Mesh::prepare(bool optimize)
{
    setOrders();

    // first step is to distribute the elements according to a global
    // distribution of DOF
    std::vector<index_t> distribution(MPIInfo->size+1);

    // first we create dense labeling for the DOFs
    index_t newGlobalNumDOFs=Nodes->createDenseDOFLabeling();

    // create a distribution of the global DOFs and determine the MPI rank
    // controlling the DOFs on this processor
    MPIInfo->setDistribution(0, newGlobalNumDOFs-1, &distribution[0]);

    // now the mesh is re-distributed according to the distribution vector
    // this will redistribute the Nodes and Elements including overlap and
    // will create an element coloring but will not create any mappings
    // (see later in this function)
    if (noError())
        distributeByRankOfDOF(distribution);

    // at this stage we are able to start an optimization of the DOF
    // distribution using ParMetis. On return distribution is altered and
    // new DOF IDs have been assigned
    if (noError() && optimize && MPIInfo->size>1) {
        optimizeDOFDistribution(distribution); 
        if (noError())
            distributeByRankOfDOF(distribution);
    }
    // the local labelling of the degrees of freedom is optimized
    if (noError() && optimize) {
        optimizeDOFLabeling(distribution); 
    }
    // rearrange elements with the aim of bringing elements closer to memory
    // locations of the nodes (distributed shared memory!):
    optimizeElementOrdering();

    // create the global indices
    if (noError()) {
        std::vector<short> maskReducedNodes(Nodes->numNodes, -1);
        std::vector<index_t> nodeDistribution(MPIInfo->size+1);
        markNodes(maskReducedNodes, 0, true);
        std::vector<index_t> indexReducedNodes = util::packMask(maskReducedNodes);

        Nodes->createDenseNodeLabeling(nodeDistribution, distribution); 
        // created reduced DOF labeling
        Nodes->createDenseReducedLabeling(maskReducedNodes, false); 
        // created reduced node labeling
        Nodes->createDenseReducedLabeling(maskReducedNodes, true);

        // create the missing mappings
        if (noError())
            Nodes->createNodeMappings(indexReducedNodes, distribution, nodeDistribution);
    }

    updateTagList();
}

/// tries to reduce the number of colours for all element files
void Mesh::createColoring(const std::vector<index_t>& dofMap)
{
    if (noError())
        Elements->createColoring(dofMap);
    if (noError())
        FaceElements->createColoring(dofMap);
    if (noError())
        Points->createColoring(dofMap);
    if (noError())
        ContactElements->createColoring(dofMap);
}

/// redistributes elements to minimize communication during assemblage
void Mesh::optimizeElementOrdering()
{
    if (noError())
        Elements->optimizeOrdering();
    if (noError())
        FaceElements->optimizeOrdering();
    if (noError())
        Points->optimizeOrdering();
    if (noError())
        ContactElements->optimizeOrdering();
}

/// regenerates list of tags in use for node file and element files
void Mesh::updateTagList()
{
    if (noError()) Nodes->updateTagList();
    if (noError()) Elements->updateTagList();
    if (noError()) FaceElements->updateTagList();
    if (noError()) Points->updateTagList();
    if (noError()) ContactElements->updateTagList();
}

/// assigns new node reference numbers to all element files
void Mesh::relabelElementNodes(const std::vector<index_t>& newNode, index_t offset)
{
    Elements->relabelNodes(newNode, offset);
    FaceElements->relabelNodes(newNode, offset);
    ContactElements->relabelNodes(newNode, offset);
    Points->relabelNodes(newNode, offset);
}

void Mesh::resolveNodeIds()
{
    // Initially the element nodes refer to the numbering defined by the global
    // id assigned to the nodes in the NodeFile. It is also not ensured that
    // all nodes referred by an element are actually available on the process.
    // At the output, a local node labeling is used and all nodes are
    // available. In particular the numbering of the element nodes is between
    // 0 and NodeFile->numNodes.
    // The function does not create a distribution of the degrees of freedom.

    // find the minimum and maximum id used by elements
    index_t min_id=std::numeric_limits<index_t>::max();
    index_t max_id=std::numeric_limits<index_t>::min();
    std::pair<index_t,index_t> range(Elements->getNodeRange());
    max_id=std::max(max_id,range.second);
    min_id=std::min(min_id,range.first);
    range=FaceElements->getNodeRange();
    max_id=std::max(max_id,range.second);
    min_id=std::min(min_id,range.first);
    range=ContactElements->getNodeRange();
    max_id=std::max(max_id,range.second);
    min_id=std::min(min_id,range.first);
    range=Points->getNodeRange();
    max_id=std::max(max_id,range.second);
    min_id=std::min(min_id,range.first);
#ifdef Finley_TRACE
    index_t global_min_id, global_max_id;
#ifdef ESYS_MPI
    index_t id_range[2], global_id_range[2];
    id_range[0]=-min_id;
    id_range[1]=max_id;
    MPI_Allreduce(id_range, global_id_range, 2, MPI_DIM_T, MPI_MAX, MPIInfo->comm);
    global_min_id=-global_id_range[0];
    global_max_id=global_id_range[1];
#else
    global_min_id=min_id;
    global_max_id=max_id;
#endif
    printf("Node id range used by elements is %d:%d\n",global_min_id,global_max_id);
#endif
    if (min_id>max_id) {
        max_id=-1;
        min_id=0;
    }
  
    // allocate mappings for new local node labeling to global node labeling
    // (newLocalToGlobalNodeLabels) and global node labeling to the new local
    // node labeling (globalToNewLocalNodeLabels[i-min_id] is the new local id
    // of global node i)
    index_t len=(max_id>=min_id) ? max_id-min_id+1 : 0;

    // mark the nodes referred by elements in usedMask
    std::vector<short> usedMask(len, -1);
    markNodes(usedMask, min_id, false);

    // create a local labeling newLocalToGlobalNodeLabels of the local nodes
    // by packing the mask usedMask
    std::vector<index_t> newLocalToGlobalNodeLabels=util::packMask(usedMask);
    const dim_t newNumNodes=newLocalToGlobalNodeLabels.size();
    usedMask.clear();

    // invert the new labeling and shift the index newLocalToGlobalNodeLabels
    // to global node ids
    std::vector<index_t> globalToNewLocalNodeLabels(len, -1);

#pragma omp parallel for
    for (index_t n=0; n<newNumNodes; n++) {
#ifdef BOUNDS_CHECK
        if (newLocalToGlobalNodeLabels[n] >= len || newLocalToGlobalNodeLabels[n] < 0) {
            printf("BOUNDS_CHECK %s %d n=%d\n", __FILE__, __LINE__, n);
            exit(1);
        }
#endif
        globalToNewLocalNodeLabels[newLocalToGlobalNodeLabels[n]]=n;
        newLocalToGlobalNodeLabels[n]+=min_id;
    }

    // create a new table
    NodeFile *newNodeFile=new NodeFile(getDim(), MPIInfo);
    if (noError()) {
        newNodeFile->allocTable(newNumNodes);
    }
    if (noError()) {
        if (len)
            newNodeFile->gather_global(&newLocalToGlobalNodeLabels[0], Nodes);
        else
            newNodeFile->gather_global(NULL, Nodes);
    }
    if (noError()) {
        delete Nodes;
        Nodes=newNodeFile;
        // relabel nodes of the elements
        relabelElementNodes(globalToNewLocalNodeLabels, min_id);
    } else
        throw FinleyAdapterException("Errors occurred during node resolution");
}

/// sets new coordinates for the nodes
void Mesh::setCoordinates(const escript::Data& newX)
{
    Nodes->setCoordinates(newX);
}

void Mesh::addTagMap(const char* name, int tag_key) 
{
   tagMap[std::string(name)]=tag_key;
}

int Mesh::getTag(const char* name) const
{
    TagMap::const_iterator it = tagMap.find(name);
    if (it == tagMap.end()) {
        std::stringstream ss;
        ss << "getTag: unknown tag name " << name << ".";
        const std::string errorMsg(ss.str());
        setError(VALUE_ERROR, errorMsg.c_str());
        return -1;
    }
    return it->second;
}

bool Mesh::isValidTagName(const char* name) const
{
   return (tagMap.count(std::string(name)) > 0);
}

} // namespace finley

