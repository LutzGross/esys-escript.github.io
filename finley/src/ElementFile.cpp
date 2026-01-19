
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

#include "ElementFile.h"

#include <escript/Data.h>
#include <escript/index.h>

#include <algorithm> // std::swap

namespace finley {

/// constructor
/// use ElementFile::allocTable to allocate the element table
ElementFile::ElementFile(const_ReferenceElementSet_ptr refSet,
                         escript::JMPI mpiInfo) :
    MPIInfo(mpiInfo),
    referenceElementSet(refSet),
    numElements(0),
    Id(NULL),
    Tag(NULL),
    Owner(NULL),
    Nodes(NULL),
    Color(NULL),
    minColor(0),
    maxColor(-1)
{
    jacobians = new ElementFile_Jacobians(
            referenceElementSet->referenceElement->BasisFunctions);
    jacobians_reducedQ = new ElementFile_Jacobians(
            referenceElementSet->referenceElementReducedQuadrature->BasisFunctions);
    jacobians_reducedS = new ElementFile_Jacobians(
            referenceElementSet->referenceElement->LinearBasisFunctions);
    jacobians_reducedS_reducedQ = new ElementFile_Jacobians(
            referenceElementSet->referenceElementReducedQuadrature->LinearBasisFunctions);

    numNodes = referenceElementSet->getNumNodes();
}

/// destructor
ElementFile::~ElementFile()
{
    freeTable();   
    delete jacobians;
    delete jacobians_reducedS;
    delete jacobians_reducedQ;
    delete jacobians_reducedS_reducedQ;
}

/// allocates the element table within this element file to hold NE elements.
void ElementFile::allocTable(dim_t NE) 
{
    if (numElements > 0)
        freeTable();

    numElements = NE;
    Owner = new int[numElements];
    Id = new index_t[numElements];
    Nodes = new index_t[numElements * numNodes];
    Tag = new int[numElements];
    Color = new index_t[numElements];

    // this initialization makes sure that data are located on the right
    // processor
#pragma omp parallel for
    for (index_t e = 0; e < numElements; e++) {
        for (int i = 0; i < numNodes; i++)
            Nodes[INDEX2(i, e, numNodes)] = -1;
        Owner[e] = -1;
        Id[e] = -1;
        Tag[e] = -1;
        Color[e] = -1;
    }
}

/// deallocates the element table within this element file
void ElementFile::freeTable()
{
    delete[] Owner;
    delete[] Id;
    delete[] Nodes;
    delete[] Tag;
    delete[] Color;
    tagsInUse.clear();
    numElements = 0;
    maxColor = -1;
    minColor = 0;
}

/// copies element file 'in' into this element file starting from 'offset'.
/// The elements offset to in->numElements+offset-1 will be overwritten
void ElementFile::copyTable(index_t offset, index_t nodeOffset,
                            index_t idOffset, const ElementFile* in)
{
    const int NN_in = in->numNodes;
    if (NN_in > numNodes) {
        throw escript::ValueError("ElementFile::copyTable: dimensions of element files don't match.");
    }

#pragma omp parallel for
    for (index_t n = 0; n < in->numElements; n++) {
        Owner[offset + n] = in->Owner[n];
        Id[offset + n] = in->Id[n] + idOffset;
        Tag[offset + n] = in->Tag[n];
        for (int i = 0; i < numNodes; i++)
            Nodes[INDEX2(i, offset + n, numNodes)] =
                            in->Nodes[INDEX2(i, n, NN_in)] + nodeOffset;
    }
}

void ElementFile::gather(const index_t* index, const ElementFile* in)
{
    const int NN_in = in->numNodes;
#pragma omp parallel for
    for (index_t e = 0; e < numElements; e++) {
        const index_t k = index[e];
        Id[e] = in->Id[k];
        Tag[e] = in->Tag[k];
        Owner[e] = in->Owner[k];
        Color[e] = in->Color[k] + maxColor + 1;
        for (int j = 0; j < std::min(numNodes, NN_in); j++)
            Nodes[INDEX2(j, e, numNodes)] = in->Nodes[INDEX2(j, k, NN_in)];
    }
    minColor = std::min(minColor, in->minColor+maxColor+1);
    maxColor = std::max(maxColor, in->maxColor+maxColor+1);
}

/// scatters the ElementFile in into this ElementFile.
/// A conservative assumption on the coloring is made.
void ElementFile::scatter(index_t* index, const ElementFile* in)
{
    const int NN_in = in->numNodes;
#pragma omp parallel for
    for (index_t e = 0; e < in->numElements; e++) {
        const index_t k = index[e];
        Owner[k] = in->Owner[e];
        Id[k] = in->Id[e];
        Tag[k] = in->Tag[e];
        Color[k] = in->Color[e]+maxColor+1;
        for (int j = 0; j < std::min(numNodes,NN_in); j++)
            Nodes[INDEX2(j,k,numNodes)] = in->Nodes[INDEX2(j,e,NN_in)];
    }
    minColor = std::min(minColor, in->minColor+maxColor+1);
    maxColor = std::max(maxColor, in->maxColor+maxColor+1);
}

void ElementFile::swapTable(ElementFile* other)
{
    std::swap(numElements, other->numElements);
    std::swap(Owner, other->Owner);
    std::swap(Id, other->Id);
    std::swap(Nodes, other->Nodes);
    std::swap(Tag, other->Tag);
    std::swap(Color, other->Color);
    std::swap(minColor, other->minColor);
    std::swap(maxColor, other->maxColor);
    std::swap(tagsInUse, other->tagsInUse);
}

void ElementFile::optimizeOrdering()
{
    if (numElements < 1)
        return;

    const int NN = referenceElementSet->getNumNodes();
    util::ValueAndIndexList item_list(numElements);
    index_t* index = new index_t[numElements];
    ElementFile* out = new ElementFile(referenceElementSet, MPIInfo);
    out->allocTable(numElements);
#pragma omp parallel for
    for (index_t e = 0; e < numElements; e++) {
        std::pair<index_t,index_t> entry(Nodes[INDEX2(0, e, NN)], e);
        for (int i = 1; i < NN; i++)
            entry.first = std::min(entry.first, Nodes[INDEX2(i, e, NN)]);
        item_list[e] = entry;
    }
    util::sortValueAndIndex(item_list);

#pragma omp parallel for
    for (index_t e = 0; e < numElements; e++)
        index[e] = item_list[e].second;

    out->gather(index, this);
    swapTable(out);
    delete out;
    delete[] index;
}

/// assigns new node reference numbers to the elements.
/// If k is the old node, the new node is newNode[k-offset].
void ElementFile::relabelNodes(const std::vector<index_t>& newNode, index_t offset)
{
#pragma omp parallel for
    for (index_t j=0; j<numElements; j++) {
        for (int i=0; i<numNodes; i++) {
            Nodes[INDEX2(i,j,numNodes)]=
                        newNode[Nodes[INDEX2(i,j,numNodes)]-offset];
        }
    }
}

void ElementFile::setTags(int newTag, const escript::Data& mask)
{
    const int numQuad = referenceElementSet->borrowReferenceElement(
            util::hasReducedIntegrationOrder(mask))
            ->Parametrization->numQuadNodes; 
    if (1 != mask.getDataPointSize()) {
        throw escript::ValueError("ElementFile::setTags: number of components of mask must be 1.");
    } else if (mask.getNumDataPointsPerSample() != numQuad ||
            mask.getNumSamples() != numElements) {
        throw escript::ValueError("ElementFile::setTags: illegal number of samples of mask Data object");
    }

    if (mask.actsExpanded()) {
#pragma omp parallel for
        for (index_t n = 0; n < numElements; n++) {
            if (mask.getSampleDataRO(n)[0] > 0)
                Tag[n] = newTag;
        }
    } else {
#pragma omp parallel for
        for (index_t n = 0; n < numElements; n++) {
            const double* mask_array = mask.getSampleDataRO(n);
            bool check = false;
            for (int q = 0; q < numQuad; q++)
                check = check || mask_array[q];
            if (check)
                Tag[n] = newTag;
        }
    }
    updateTagList();
}

/// Tries to reduce the number of colours used to colour the elements
void ElementFile::createColoring(const IndexVector& dofMap)
{
    if (numElements < 1)
        return;

    const int NN = numNodes;
    const std::pair<index_t,index_t> idRange(util::getMinMaxInt(
                                            1, dofMap.size(), &dofMap[0]));
    const index_t len = idRange.second-idRange.first+1;

    // reset color vector
#pragma omp parallel for
    for (index_t e = 0; e < numElements; e++)
        Color[e] = -1;

    index_t numUncoloredElements = numElements;
    minColor = 0;
    maxColor = -1;
    while (numUncoloredElements>0) {
        // initialize the mask marking nodes used by a color
        std::vector<index_t> maskDOF(len, -1);
        numUncoloredElements=0;

        // TODO: OMP
        for (index_t e=0; e<numElements; e++) {
            if (Color[e] < 0) {
                // find out if element e is independent from the elements
                // already coloured:
                bool independent = true; 
                for (int i=0; i<NN; i++) {
#ifdef BOUNDS_CHECK
                    ESYS_ASSERT(Nodes[INDEX2(i, e, NN)] >= 0, "BOUNDS_CHECK");
                    ESYS_ASSERT(Nodes[INDEX2(i, e, NN)] < dofMap.size(), "BOUNDS_CHECK");
                    ESYS_ASSERT(dofMap[Nodes[INDEX2(i, e, NN)]] - idRange.first < len, "BOUNDS_CHECK");
                    ESYS_ASSERT(dofMap[Nodes[INDEX2(i, e, NN)]] - idRange.first >= 0, "BOUNDS_CHECK");
#endif
                    if (maskDOF[dofMap[Nodes[INDEX2(i,e,NN)]]-idRange.first]>0) {
                        independent=false;
                        break;
                    }
                }
                // if e is independent a new color is assigned and the nodes
                // are marked as being used
                if (independent) {
                    for (int i=0; i<NN; i++)
                        maskDOF[dofMap[Nodes[INDEX2(i,e,NN)]]-idRange.first] = 1;
                    Color[e]=maxColor+1;
                } else {
                    numUncoloredElements++;
                }
            } // if no colour yet
        } // for all elements
        maxColor++;
    } // end of while loop
}

void ElementFile::markNodes(std::vector<short>& mask, int offset, bool useLinear)
{
    const_ReferenceElement_ptr refElement(referenceElementSet->
                                            borrowReferenceElement(false));
    if (useLinear) {
        const int NN=refElement->numLinearNodes;
        const int *lin_nodes=refElement->Type->linearNodes;
#pragma omp parallel for
        for (index_t e=0; e<numElements; e++) {
            for (int i=0; i<NN; i++) {
                mask[Nodes[INDEX2(lin_nodes[i],e,numNodes)]-offset]=1;
            }
        }
    } else {
        const int NN=refElement->Type->numNodes;
#pragma omp parallel for
        for (index_t e=0; e<numElements; e++) {
            for (int i=0; i<NN; i++) {
                mask[Nodes[INDEX2(i,e,numNodes)]-offset]=1;
            }
        }
    }
}

/// redistributes the elements including overlap by rank
void ElementFile::distributeByRankOfDOF(const std::vector<int>& mpiRankOfDOF, index_t* index)
{
    const int size=MPIInfo->size;

    if (size > 1) {
#ifdef ESYS_MPI
        int myRank=MPIInfo->rank;
        int numRequests=0;
        std::vector<MPI_Request> mpi_requests(8*size);
        std::vector<MPI_Status> mpi_stati(8*size);

        // count the number of elements that have to be sent to each processor
        // (send_count) and define a new element owner as the processor with
        // the largest number of DOFs and the smallest id
        std::vector<dim_t> send_count(size);
        std::vector<dim_t> recv_count(size);
        std::vector<int> newOwner(numElements);
#pragma omp parallel
        {
            std::vector<dim_t> loc_send_count(size);
#pragma omp for
            for (index_t e=0; e<numElements; e++) {
                if (Owner[e] == myRank) {
                    newOwner[e]=myRank;
                    std::vector<dim_t> loc_proc_mask(size);
                    for (int j=0; j<numNodes; j++) {
                        const int p=mpiRankOfDOF[Nodes[INDEX2(j,e,numNodes)]];
                        loc_proc_mask[p]++;
                    }
                    dim_t loc_proc_mask_max=0;
                    for (int p=0; p<size; ++p) {
                        if (loc_proc_mask[p] > 0)
                            loc_send_count[p]++;
                        if (loc_proc_mask[p] > loc_proc_mask_max) {
                            newOwner[e]=p;
                            loc_proc_mask_max=loc_proc_mask[p];
                        }
                    }
                } else {
                    newOwner[e]=-1;
                }
            }
#pragma omp critical
            {
                for (int p=0; p<size; ++p)
                    send_count[p]+=loc_send_count[p];
            }
        }
        MPI_Alltoall(&send_count[0], 1, MPI_DIM_T, &recv_count[0], 1, MPI_DIM_T,
                     MPIInfo->comm);
        // get the new number of elements for this processor
        dim_t newNumElements=0;
        dim_t numElementsInBuffer=0;
        for (int p=0; p<size; ++p) {
            newNumElements+=recv_count[p];
            numElementsInBuffer+=send_count[p];
        }

        std::vector<index_t> Id_buffer(numElementsInBuffer);
        std::vector<int> Tag_buffer(numElementsInBuffer);
        std::vector<int> Owner_buffer(numElementsInBuffer);
        std::vector<index_t> Nodes_buffer(numElementsInBuffer*numNodes);
        std::vector<index_t> send_offset(size);
        std::vector<index_t> recv_offset(size);
        std::vector<unsigned char> proc_mask(size);

        // calculate the offsets for the processor buffers
        for (int p=0; p<size-1; ++p) {
            recv_offset[p+1]=recv_offset[p]+recv_count[p];
            send_offset[p+1]=send_offset[p]+send_count[p];
        }

        send_count.assign(size, 0);
        // copy element into buffers. proc_mask makes sure that an element is
        // copied once only for each processor
        for (index_t e=0; e<numElements; e++) {
            if (Owner[e] == myRank) {
                proc_mask.assign(size, true);
                for (int j=0; j<numNodes; j++) {
                    const int p=mpiRankOfDOF[Nodes[INDEX2(j,e,numNodes)]];
                    if (proc_mask[p]) {
                        index_t k=send_offset[p]+send_count[p];
                        Id_buffer[k]=Id[e];
                        Tag_buffer[k]=Tag[e];
                        Owner_buffer[k]=newOwner[e];
                        for (int i=0; i<numNodes; i++)
                            Nodes_buffer[INDEX2(i,k,numNodes)]=
                                    index[Nodes[INDEX2(i,e,numNodes)]];
                        send_count[p]++;
                        proc_mask[p]=false;
                    }
                }
            }
        }
        // allocate new tables
        allocTable(newNumElements);

        // start to receive new elements
        for (int p=0; p<size; ++p) {
            if (recv_count[p] > 0) {
                MPI_Irecv(&Id[recv_offset[p]], recv_count[p], MPI_DIM_T, p,
                        MPIInfo->counter()+myRank, MPIInfo->comm,
                        &mpi_requests[numRequests]);
                numRequests++;
                MPI_Irecv(&Tag[recv_offset[p]], recv_count[p], MPI_INT, p,
                        MPIInfo->counter()+size+myRank, MPIInfo->comm,
                        &mpi_requests[numRequests]);
                numRequests++;
                MPI_Irecv(&Owner[recv_offset[p]], recv_count[p], MPI_INT, p,
                        MPIInfo->counter()+2*size+myRank, MPIInfo->comm,
                        &mpi_requests[numRequests]);
                numRequests++;
                MPI_Irecv(&Nodes[recv_offset[p]*numNodes],
                        recv_count[p]*numNodes, MPI_DIM_T, p,
                        MPIInfo->counter()+3*size+myRank, MPIInfo->comm,
                        &mpi_requests[numRequests]);
                numRequests++;
            }
        }
        // now the buffers can be sent away
        for (int p=0; p<size; ++p) {
            if (send_count[p] > 0) {
                MPI_Issend(&Id_buffer[send_offset[p]], send_count[p], MPI_DIM_T,
                        p, MPIInfo->counter()+p, MPIInfo->comm,
                        &mpi_requests[numRequests]);
                numRequests++;
                MPI_Issend(&Tag_buffer[send_offset[p]], send_count[p], MPI_INT,
                        p, MPIInfo->counter()+size+p, MPIInfo->comm,
                        &mpi_requests[numRequests]);
                numRequests++;
                MPI_Issend(&Owner_buffer[send_offset[p]], send_count[p],
                        MPI_INT, p, MPIInfo->counter()+2*size+p,
                        MPIInfo->comm, &mpi_requests[numRequests]);
                numRequests++;
                MPI_Issend(&Nodes_buffer[send_offset[p]*numNodes],
                        send_count[p]*numNodes, MPI_DIM_T, p,
                        MPIInfo->counter()+3*size+p, MPIInfo->comm,
                        &mpi_requests[numRequests]);
                numRequests++;
            }
        }
        MPIInfo->incCounter(4*size);
        // wait for the requests to be finalized
        MPI_Waitall(numRequests, &mpi_requests[0], &mpi_stati[0]);
#endif
    } else { // single rank
#pragma omp parallel for
        for (index_t e=0; e<numElements; e++) {
            Owner[e]=0;
            for (int i=0; i<numNodes; i++)
                Nodes[INDEX2(i,e,numNodes)]=index[Nodes[INDEX2(i,e,numNodes)]];
        }
    }
}

#ifdef ESYS_HAVE_HDF5
void ElementFile::dump(H5::Group h5_grp) const
{
    #ifdef ESYS_INDEXTYPE_LONG
        H5::DataType h5_type_index = H5::PredType::NATIVE_LONG;
    #else
        H5::DataType h5_type_index = H5::PredType::NATIVE_INT;
    #endif
    hsize_t h5_dims_numElements[1] = { 1 };
    uint h5_values_numElements[1] = { static_cast<uint>(numElements) };
    uint h5_values_numNodes[1] = { static_cast<uint>(numNodes) };
    int h5_values_typeId[1] = { static_cast<int>(referenceElementSet->referenceElement->Type->TypeId) };
    H5::Attribute h5_attr_numElements = h5_grp.createAttribute("numElements", H5::PredType::NATIVE_UINT, H5::DataSpace(1, h5_dims_numElements ) );
    h5_attr_numElements.write( H5::PredType::NATIVE_UINT, h5_values_numElements );
    H5::Attribute h5_attr_numNodes = h5_grp.createAttribute("numNodes", H5::PredType::NATIVE_UINT, H5::DataSpace(1, h5_dims_numElements ) );
    h5_attr_numNodes.write( H5::PredType::NATIVE_UINT, h5_values_numNodes );
    H5::Attribute h5_attr_typeId = h5_grp.createAttribute("TypeId", H5::PredType::NATIVE_INT, H5::DataSpace(1, h5_dims_numElements ) );
    h5_attr_typeId.write( H5::PredType::NATIVE_INT, h5_values_typeId );

    hsize_t h5_dims_ids[1] = { static_cast<hsize_t>(numElements)  };
    H5::DataSet h5_ds_ids = h5_grp.createDataSet("Ids", H5::DataType(h5_type_index), H5::DataSpace(1, h5_dims_ids ) );
    h5_ds_ids.write( &Id[0], H5::DataType(h5_type_index));
    H5::DataSet h5_ds_tags = h5_grp.createDataSet("Tags", H5::DataType(h5_type_index), H5::DataSpace(1, h5_dims_ids ) );
    h5_ds_tags.write( &Tag[0], H5::DataType(h5_type_index));
    H5::DataSet h5_ds_owners = h5_grp.createDataSet("Owners", H5::DataType(h5_type_index), H5::DataSpace(1, h5_dims_ids ) );
    h5_ds_owners.write( &Owner[0], H5::DataType(h5_type_index));
    H5::DataSet h5_ds_colors = h5_grp.createDataSet("Colors", H5::DataType(h5_type_index), H5::DataSpace(1, h5_dims_ids ) );
    h5_ds_colors.write( &Color[0], H5::DataType(h5_type_index));
    hsize_t h5_dims_nodes[1] = { static_cast<hsize_t>(numElements * numNodes) };
    H5::DataSet h5_ds_nodes = h5_grp.createDataSet("Nodes", H5::DataType(h5_type_index), H5::DataSpace(1, h5_dims_nodes ) );
    h5_ds_nodes.write( &Nodes[0], H5::DataType(h5_type_index));
}
#endif


#ifdef ESYS_HAVE_HDF5
finley::ElementFile* loadElements_hdf5(const H5::Group h5_grp, const int integration_order, const int reduced_integration_order, const escript::JMPI mpiInfo)
{
        #ifdef ESYS_INDEXTYPE_LONG
            H5::DataType h5_type_index = H5::PredType::NATIVE_LONG;
        #else
            H5::DataType h5_type_index = H5::PredType::NATIVE_INT;
        #endif

        uint h5_numElements=0, h5_numNodes=0;
        int h5_typeId=0;


        H5::Attribute h5_attr_ne(h5_grp.openAttribute("numElements"));
        H5::DataType h5_type_ne(h5_attr_ne.getDataType());
        if (  h5_type_ne != H5::PredType::NATIVE_UINT  ) {
             throw FinleyException("Error - finley.load: illegal data type for number of elements in HDF5 file.");
        }
        if ( h5_attr_ne.getStorageSize() !=  1 * h5_type_ne.getSize() ) {
             throw FinleyException("Error - finley.load: number of elements in HDF5 file needs to be a single value.");
        }
        h5_attr_ne.read(h5_type_ne, &h5_numElements);


        H5::Attribute h5_attr_nn(h5_grp.openAttribute("numNodes"));
        H5::DataType h5_type_nn(h5_attr_nn.getDataType());
        if (  h5_type_nn != H5::PredType::NATIVE_UINT  ) {
             throw FinleyException("Error - finley.load: illegal data type for number of nodes per elements in HDF5 file.");
        }
        if ( h5_attr_nn.getStorageSize() !=  1 * h5_type_nn.getSize() ) {
             throw FinleyException("Error - finley.load: number of nodes per elements in HDF5 file needs to be a single value.");
        }
        h5_attr_nn.read(h5_type_nn, &h5_numNodes);

        H5::Attribute h5_attr_etype(h5_grp.openAttribute("TypeId"));
        H5::DataType h5_type_etype(h5_attr_etype.getDataType());
        if (  h5_type_etype != H5::PredType::NATIVE_INT  ) {
             throw FinleyException("Error - finley.load: illegal data type for element type in HDF5 file.");
        }
        if ( h5_attr_etype.getStorageSize() !=  1 * h5_type_etype.getSize() ) {
             throw FinleyException("Error - finley.load: element type in HDF5 file needs to be a single value.");
        }
        h5_attr_etype.read(h5_type_etype, &h5_typeId);

        // ... create reference element and element file:
        const_ReferenceElementSet_ptr refElements(new ReferenceElementSet( static_cast<ElementTypeId>(h5_typeId), integration_order, reduced_integration_order));
        ElementFile* elements = new ElementFile(refElements, mpiInfo);


        const dim_t num_Elements = static_cast<dim_t>(h5_numElements);
        elements->allocTable(num_Elements);

        H5::DataSet h5_ds_ids = h5_grp.openDataSet("Ids");
        H5::DataType h5_type_ids(h5_ds_ids.getDataType());
        if (  h5_type_ids != H5::DataType(h5_type_index)  ) {
             throw FinleyException("Error - finley.load: illegal data type for element ids in HDF5 file.");
        }
        if ( h5_ds_ids.getStorageSize() !=  num_Elements * h5_type_ids.getSize() ) {
             throw FinleyException("Error - finley.load: number of element ids in HDF5 file is incorrect.");
        }
        h5_ds_ids.read(&elements->Id[0], h5_type_ids);

        H5::DataSet h5_ds_tags = h5_grp.openDataSet("Tags");
        H5::DataType h5_type_tags(h5_ds_tags.getDataType());
        if (  h5_type_tags != H5::DataType(h5_type_index)  ) {
             throw FinleyException("Error - finley.load: illegal data type for element tags in HDF5 file.");
        }
        if ( h5_ds_tags.getStorageSize() !=  num_Elements * h5_type_tags.getSize() ) {
             throw FinleyException("Error - finley.load: number of element tags in HDF5 file is incorrect.");
        }
        h5_ds_tags.read(&elements->Tag[0], h5_type_tags);

        H5::DataSet h5_ds_ows = h5_grp.openDataSet("Owners");
        H5::DataType h5_type_ows(h5_ds_ows.getDataType());
        if (  h5_type_ows != H5::DataType(h5_type_index)  ) {
             throw FinleyException("Error - finley.load: illegal data type for element owners in HDF5 file.");
        }
        if ( h5_ds_ows.getStorageSize() !=  num_Elements * h5_type_ows.getSize() ) {
             throw FinleyException("Error - finley.load: number of element owners in HDF5 file is incorrect.");
        }
        h5_ds_ows.read(&elements->Owner[0], h5_type_ows);

        H5::DataSet h5_ds_cls = h5_grp.openDataSet("Colors");
        H5::DataType h5_type_cls(h5_ds_cls.getDataType());
        if (  h5_type_cls != H5::DataType(h5_type_index)  ) {
             throw FinleyException("Error - finley.load: illegal data type for element colors in HDF5 file.");
        }
        if ( h5_ds_cls.getStorageSize() !=  num_Elements * h5_type_cls.getSize() ) {
             throw FinleyException("Error - finley.load: number of element colors in HDF5 file is incorrect.");
        }
        h5_ds_cls.read(&elements->Owner[0], h5_type_cls);

        H5::DataSet h5_ds_nds = h5_grp.openDataSet("Nodes");
        H5::DataType h5_type_nds(h5_ds_nds.getDataType());
        if (  h5_type_nds != H5::DataType(h5_type_index)  ) {
             throw FinleyException("Error - finley.load: illegal data type for element nodes in HDF5 file.");
        }
std::cout <<  h5_ds_nds.getStorageSize() << " " <<  num_Elements << " " <<  (elements->numNodes) << " " <<  h5_type_nds.getSize()  << " " <<  h5_numNodes << std::endl;
        if ( h5_ds_nds.getStorageSize() !=  num_Elements * (elements->numNodes) * h5_type_nds.getSize() ) {
             throw FinleyException("Error - finley.load: number of element nodes in HDF5 file is incorrect.");
        }
        h5_ds_nds.read(&elements->Nodes[0], h5_type_nds);

       // .... Now we need to adjust maxColor
       elements->minColor = 0;
       elements->maxColor = num_Elements-1;
       if  ( num_Elements > 0 ) {
            index_t mc = elements->Color[0];
            #pragma omp parallel
            for (index_t i = 1; i < num_Elements; ++i) {
                if (mc < elements->Color[i]) {
                    #pragma omp critical
                    {
                        mc = elements->Color[i];
                    }
                }
            }
            elements->maxColor = mc;
       }
       elements->updateTagList();
       return elements;
}
#endif // ESYS_HAVE_HDF5
} // namespace finley

