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


/****************************************************************************

  Finley: Mesh

*****************************************************************************/

#include "Mesh.h"
#include "IndexList.h"

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
    Points(NULL),
    FullFullPattern(NULL),
    FullReducedPattern(NULL),
    ReducedFullPattern(NULL),
    ReducedReducedPattern(NULL)
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
    Paso_SystemMatrixPattern_free(FullFullPattern);
    Paso_SystemMatrixPattern_free(FullReducedPattern);
    Paso_SystemMatrixPattern_free(ReducedFullPattern);
    Paso_SystemMatrixPattern_free(ReducedReducedPattern);
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
void Mesh::createMappings(const std::vector<int>& dofDistribution,
                          const std::vector<int>& nodeDistribution)
{
    std::vector<short> maskReducedNodes(Nodes->numNodes, -1);
    markNodes(maskReducedNodes, 0, true);
    std::vector<int> indexReducedNodes = util::packMask(maskReducedNodes);
    if (noError())
        Nodes->createNodeMappings(indexReducedNodes, dofDistribution,
                                  nodeDistribution);
}

/// redistributes the Nodes and Elements including overlap
/// according to the DOF distribution. It will create an element colouring
/// but will not create any mappings.
void Mesh::distributeByRankOfDOF(const std::vector<int>& dof_distribution)
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
    const std::pair<int,int> dof_range(Nodes->getDOFRange());
    const int len=dof_range.second-dof_range.first+1;
    // local mask for used nodes
    std::vector<int> localDOF_mask(len, -1);
    std::vector<int> localDOF_map(Nodes->numNodes, -1);

#pragma omp parallel for
    for (int n=0; n<Nodes->numNodes; n++) {
#ifdef BOUNDS_CHECK
        if ((Nodes->globalDegreesOfFreedom[n]-dof_range.first) >= len ||
                (Nodes->globalDegreesOfFreedom[n]-dof_range.first) < 0) {
            printf("BOUNDS_CHECK %s %d\n", __FILE__, __LINE__);
            exit(1);
        }
#endif
        localDOF_mask[Nodes->globalDegreesOfFreedom[n]-dof_range.first]=n;
    }

    int numDOFs=0;
    for (int n=0; n<len; n++) {
        const int k=localDOF_mask[n];
        if (k>=0) {
             localDOF_mask[n]=numDOFs;
             numDOFs++;
          }
    }
#pragma omp parallel for
    for (int n=0; n<Nodes->numNodes; n++) {
        const int k=localDOF_mask[Nodes->globalDegreesOfFreedom[n]-dof_range.first];
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
        printf("=== %s:\nnumber of elements=%d\ncolor range=[%d,%d]\n",
               Elements->referenceElementSet->referenceElement->Type->Name,
               Elements->numElements, Elements->minColor, Elements->maxColor);
        if (Elements->numElements > 0) {
            const int NN=Elements->referenceElementSet->referenceElement->Type->numNodes;
            const int NN2=Elements->numNodes;
            printf("Id,Tag,Owner,Color,Nodes\n");
            for (int i=0; i<Elements->numElements; i++) {
                printf("%d,%d,%d,%d,", Elements->Id[i], Elements->Tag[i],
                        Elements->Owner[i], Elements->Color[i]);
                for (int j=0; j<NN; j++)
                    printf(" %d", Nodes->Id[Elements->Nodes[INDEX2(j,i,NN2)]]);
                printf("\n");
            }
        }
    }

    // write face elements
    if (FaceElements) {
        printf("=== %s:\nnumber of elements=%d\ncolor range=[%d,%d]\n",
               FaceElements->referenceElementSet->referenceElement->Type->Name,
               FaceElements->numElements, FaceElements->minColor,
               FaceElements->maxColor);
        if (FaceElements->numElements > 0) {
            const int NN=FaceElements->referenceElementSet->referenceElement->Type->numNodes;
            const int NN2=FaceElements->numNodes;
            printf("Id,Tag,Owner,Color,Nodes\n");
            for (int i=0; i<FaceElements->numElements; i++) {
                printf("%d,%d,%d,%d,", FaceElements->Id[i],
                        FaceElements->Tag[i], FaceElements->Owner[i],
                        FaceElements->Color[i]);
                for (int j=0; j<NN; j++)
                    printf(" %d", Nodes->Id[FaceElements->Nodes[INDEX2(j,i,NN2)]]);
                printf("\n");
            }
        }
    }

    // write Contact elements
    if (ContactElements) {
        printf("=== %s:\nnumber of elements=%d\ncolor range=[%d,%d]\n",
               ContactElements->referenceElementSet->referenceElement->Type->Name,
               ContactElements->numElements, ContactElements->minColor,
               ContactElements->maxColor);
        if (ContactElements->numElements > 0) {
            const int NN=ContactElements->referenceElementSet->referenceElement->Type->numNodes;
            const int NN2=ContactElements->numNodes;
            printf("Id,Tag,Owner,Color,Nodes\n");
            for (int i=0; i<ContactElements->numElements; i++) {
                printf("%d,%d,%d,%d,", ContactElements->Id[i],
                        ContactElements->Tag[i], ContactElements->Owner[i],
                        ContactElements->Color[i]);
                for (int j=0; j<NN; j++)
                    printf(" %d", Nodes->Id[ContactElements->Nodes[INDEX2(j,i,NN2)]]);
                printf("\n");
            }
        }
    }
  
    // write points
    if (Points) {
        printf("=== %s:\nnumber of elements=%d\ncolor range=[%d,%d]\n",
               Points->referenceElementSet->referenceElement->Type->Name,
               Points->numElements, Points->minColor, Points->maxColor);
        if (Points->numElements > 0) {
            const int NN=Points->referenceElementSet->referenceElement->Type->numNodes;
            const int NN2=Points->numNodes;
            printf("Id,Tag,Owner,Color,Nodes\n");
            for (int i=0; i<Points->numElements; i++) {
                printf("%d,%d,%d,%d,", Points->Id[i], Points->Tag[i],
                        Points->Owner[i], Points->Color[i]);
                for (int j=0; j<NN; j++)
                    printf(" %d", Nodes->Id[Points->Nodes[INDEX2(j,i,NN2)]]);
                printf("\n");
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
                                    int firstDOF, int lastDOF, bool useLinear)
{
    const int *dofIndex = (useLinear ? Nodes->globalReducedDOFIndex
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
void Mesh::optimizeDOFLabeling(const std::vector<int>& distribution)
{
    const int myRank=MPIInfo->rank;
    const int mpiSize=MPIInfo->size;
    const int myFirstVertex=distribution[myRank];
    const int myLastVertex=distribution[myRank+1];
    const int myNumVertices=myLastVertex-myFirstVertex;
    int len=0;
    for (int p=0; p<mpiSize; ++p)
        len=std::max(len, distribution[p+1]-distribution[p]);

    IndexList* index_list=new IndexList[myNumVertices];
    std::vector<int> newGlobalDOFID(len);
    // create the adjacency structure xadj and adjncy
#pragma omp parallel
    {
        // insert contributions from element matrices into columns index
        IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list,
                myFirstVertex, myLastVertex, Elements,
                Nodes->globalDegreesOfFreedom,
                Nodes->globalDegreesOfFreedom);
        IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list,
                myFirstVertex, myLastVertex, FaceElements,
                Nodes->globalDegreesOfFreedom,
                Nodes->globalDegreesOfFreedom);
        IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list,
                myFirstVertex, myLastVertex, ContactElements,
                Nodes->globalDegreesOfFreedom,
                Nodes->globalDegreesOfFreedom);
        IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list,
                myFirstVertex, myLastVertex, Points,
                Nodes->globalDegreesOfFreedom,
                Nodes->globalDegreesOfFreedom);
    }
    // create the local matrix pattern
    Paso_Pattern *pattern=IndexList_createPattern(0, myNumVertices,
            index_list, myFirstVertex, myLastVertex, -myFirstVertex);

    if (noError())
        Paso_Pattern_reduceBandwidth(pattern, &newGlobalDOFID[0]); 

    Paso_Pattern_free(pattern);
    delete[] index_list;
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
            const int firstVertex=distribution[current_rank];
            const int lastVertex=distribution[current_rank+1];
#pragma omp parallel for
            for (int i=0; i<Nodes->numNodes; ++i) {
                const int k=Nodes->globalDegreesOfFreedom[i];
                if (firstVertex<=k && k<lastVertex) {
                    Nodes->globalDegreesOfFreedom[i]=newGlobalDOFID[k-firstVertex];
                }
            }
   
            if (p<mpiSize-1) { // the final send can be skipped
#ifdef ESYS_MPI
                MPI_Status status;
                MPI_Sendrecv_replace(&newGlobalDOFID[0], len, MPI_INT,
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
    std::vector<int> distribution(MPIInfo->size+1);

    // first we create dense labeling for the DOFs
    int newGlobalNumDOFs=Nodes->createDenseDOFLabeling();

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
        std::vector<int> nodeDistribution(MPIInfo->size+1);
        markNodes(maskReducedNodes, 0, true);
        std::vector<int> indexReducedNodes = util::packMask(maskReducedNodes);

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
void Mesh::createColoring(const std::vector<int>& dofMap)
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
void Mesh::relabelElementNodes(const std::vector<int>& newNode, int offset)
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
    int min_id=std::numeric_limits<int>::max();
    int max_id=std::numeric_limits<int>::min();
    std::pair<int,int> range(Elements->getNodeRange());
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
    int global_min_id, global_max_id;
#ifdef ESYS_MPI
    int id_range[2], global_id_range[2];
    id_range[0]=-min_id;
    id_range[1]=max_id;
    MPI_Allreduce(id_range, global_id_range, 2, MPI_INT, MPI_MAX, MPIInfo->comm);
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
    int len=(max_id>=min_id) ? max_id-min_id+1 : 0;

    // mark the nodes referred by elements in usedMask
    std::vector<short> usedMask(len, -1);
    markNodes(usedMask, min_id, false);

    // create a local labeling newLocalToGlobalNodeLabels of the local nodes
    // by packing the mask usedMask
    std::vector<int> newLocalToGlobalNodeLabels=util::packMask(usedMask);
    const int newNumNodes=newLocalToGlobalNodeLabels.size();
    usedMask.clear();

    // invert the new labeling and shift the index newLocalToGlobalNodeLabels
    // to global node ids
    std::vector<int> globalToNewLocalNodeLabels(len, -1);

#pragma omp parallel for
    for (int n=0; n<newNumNodes; n++) {
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
        newNodeFile->gather_global(newLocalToGlobalNodeLabels, Nodes);
    }
    if (noError()) {
        delete Nodes;
        Nodes=newNodeFile;
        // relabel nodes of the elements
        relabelElementNodes(globalToNewLocalNodeLabels, min_id);
    }
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

