
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

  Finley: ElementFile

*****************************************************************************/

#include "ElementFile.h"
#include <escript/Data.h>

#include <algorithm> // std::swap

namespace finley {

/// constructor
/// use ElementFile::allocTable to allocate the element table
ElementFile::ElementFile(const_ReferenceElementSet_ptr refSet,
                         esysUtils::JMPI& mpiInfo) :
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
    MPIInfo=mpiInfo;
 
    jacobians=new ElementFile_Jacobians(
            referenceElementSet->referenceElement->BasisFunctions);
    jacobians_reducedQ=new ElementFile_Jacobians(
            referenceElementSet->referenceElementReducedQuadrature->BasisFunctions);
    jacobians_reducedS=new ElementFile_Jacobians(
            referenceElementSet->referenceElement->LinearBasisFunctions);
    jacobians_reducedS_reducedQ=new ElementFile_Jacobians(
            referenceElementSet->referenceElementReducedQuadrature->LinearBasisFunctions);

    numNodes=referenceElementSet->getNumNodes();
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
void ElementFile::allocTable(int NE) 
{
    if (numElements>0)
        freeTable();

    numElements=NE;
    Owner=new int[numElements];
    Id=new int[numElements];
    Nodes=new int[numElements*numNodes];
    Tag=new int[numElements];
    Color=new int[numElements];
  
    // this initialization makes sure that data are located on the right
    // processor
#pragma omp parallel for
    for (int e=0; e<numElements; e++) {
        for (int i=0; i<numNodes; i++)
            Nodes[INDEX2(i,e,numNodes)]=-1;
        Owner[e]=-1;
        Id[e]=-1;
        Tag[e]=-1;
        Color[e]=-1;
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
    numElements=0;
    maxColor=-1;
    minColor=0;
}

/// copies element file 'in' into this element file starting from 'offset'.
/// The elements offset to in->numElements+offset-1 will be overwritten
void ElementFile::copyTable(int offset, int nodeOffset, int idOffset,
                            const ElementFile* in)
{
    const int NN_in=in->numNodes;
    if (NN_in > numNodes) {
        setError(TYPE_ERROR, "ElementFile::copyTable: dimensions of element files don't match.");
        return;
    }

#pragma omp parallel for
    for (int n=0; n<in->numElements; n++) {
          Owner[offset+n]=in->Owner[n];
          Id[offset+n]=in->Id[n]+idOffset;
          Tag[offset+n]=in->Tag[n];
          for (int i=0; i<numNodes; i++)
              Nodes[INDEX2(i,offset+n,numNodes)] =
                            in->Nodes[INDEX2(i,n,NN_in)]+nodeOffset;
    }
}

void ElementFile::gather(int* index, const ElementFile* in)
{
    const int NN_in=in->numNodes;
#pragma omp parallel for
    for (int e=0; e<numElements; e++) {
        const int k=index[e];
        Id[e]=in->Id[k];
        Tag[e]=in->Tag[k];
        Owner[e]=in->Owner[k];
        Color[e]=in->Color[k]+maxColor+1;
        for (int j=0; j<std::min(numNodes,NN_in); j++)
            Nodes[INDEX2(j,e,numNodes)]=in->Nodes[INDEX2(j,k,NN_in)];
    }
    minColor=std::min(minColor, in->minColor+maxColor+1);
    maxColor=std::max(maxColor, in->maxColor+maxColor+1);
}

/// scatters the ElementFile in into this ElementFile.
/// A conservative assumption on the coloring is made.
void ElementFile::scatter(int* index, const ElementFile* in)
{
    const int NN_in=in->numNodes;
#pragma omp parallel for
    for (int e=0; e<in->numElements; e++) {
        const int k=index[e];
        Owner[k]=in->Owner[e];
        Id[k]=in->Id[e];
        Tag[k]=in->Tag[e];
        Color[k]=in->Color[e]+maxColor+1;
        for (int j=0; j<std::min(numNodes,NN_in); j++)
            Nodes[INDEX2(j,k,numNodes)]=in->Nodes[INDEX2(j,e,NN_in)];
    }
    minColor=std::min(minColor, in->minColor+maxColor+1);
    maxColor=std::max(maxColor, in->maxColor+maxColor+1);
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
    if (numElements<1)
        return;

    const int NN=referenceElementSet->getNumNodes();
    util::ValueAndIndexList item_list(numElements);
    int *index=new int[numElements];
    ElementFile* out=new ElementFile(referenceElementSet, MPIInfo);
    out->allocTable(numElements);
    if (noError()) {
#pragma omp parallel for
        for (int e=0; e<numElements; e++) {
            std::pair<int,int> entry(Nodes[INDEX2(0,e,NN)], e);
            for (int i=1; i<NN; i++)
                entry.first=std::min(entry.first, Nodes[INDEX2(i,e,NN)]);
            item_list[e] = entry;
        }
        util::sortValueAndIndex(item_list);
#pragma omp parallel for
        for (int e=0; e<numElements; e++)
            index[e]=item_list[e].second;
        out->gather(index, this);
        swapTable(out);
    }
    delete out;
    delete[] index;
}

/// assigns new node reference numbers to the elements.
/// If k is the old node, the new node is newNode[k-offset].
void ElementFile::relabelNodes(const std::vector<int>& newNode, int offset)
{
#pragma omp parallel for
    for (int j=0; j<numElements; j++) {
        for (int i=0; i<numNodes; i++) {
            Nodes[INDEX2(i,j,numNodes)]=
                        newNode[Nodes[INDEX2(i,j,numNodes)]-offset];
        }
    }
}

void ElementFile::setTags(const int newTag, const escript::Data& mask)
{
    resetError();

    const int numQuad=referenceElementSet->borrowReferenceElement(
            util::hasReducedIntegrationOrder(mask))
            ->Parametrization->numQuadNodes; 
    if (1 != mask.getDataPointSize()) {
        setError(TYPE_ERROR, "ElementFile::setTags: number of components of mask must be 1.");
        return;
    } else if (mask.getNumDataPointsPerSample() != numQuad ||
            mask.getNumSamples() != numElements) {
        setError(TYPE_ERROR, "ElementFile::setTags: illegal number of samples of mask Data object");
        return;
    }

    if (mask.actsExpanded()) {
#pragma omp parallel for
        for (int n=0; n<numElements; n++) {
            if (mask.getSampleDataRO(n)[0] > 0)
                Tag[n]=newTag;
        }
    } else {
#pragma omp parallel for
        for (int n=0; n<numElements; n++) {
            const double *mask_array=mask.getSampleDataRO(n);
            bool check=false;
            for (int q=0; q<numQuad; q++)
                check = (check || mask_array[q]);
            if (check)
                Tag[n]=newTag;
        }
    }
    updateTagList();
}

/// Tries to reduce the number of colours used to colour the elements
void ElementFile::createColoring(const std::vector<int>& dofMap)
{
    if (numElements < 1)
        return;

    const int NN = numNodes;
    const std::pair<int,int> idRange(util::getMinMaxInt(
                                            1, dofMap.size(), &dofMap[0]));
    const int len=idRange.second-idRange.first+1;

    // reset color vector
#pragma omp parallel for
    for (int e=0; e<numElements; e++)
        Color[e]=-1;

    int numUncoloredElements=numElements;
    minColor=0;
    maxColor=-1;
    while (numUncoloredElements>0) {
        // initialize the mask marking nodes used by a color
        std::vector<int> maskDOF(len, -1);
        numUncoloredElements=0;

        // TODO: OMP
        for (int e=0; e<numElements; e++) {
            if (Color[e] < 0) {
                // find out if element e is independent from the elements
                // already coloured:
                bool independent = true; 
                for (int i=0; i<NN; i++) {
#ifdef BOUNDS_CHECK
if (Nodes[INDEX2(i,e,NN)] < 0 || Nodes[INDEX2(i,e,NN)] >= dofMap.size()) {
    printf("BOUNDS_CHECK %s %d i=%d e=%d NN=%d min_id=%d Nodes[INDEX2...]=%d\n",
            __FILE__, __LINE__, i, e, NN, idRange.first, Nodes[INDEX2(i,e,NN)]);
    exit(1);
}
if ((dofMap[Nodes[INDEX2(i,e,NN)]]-idRange.first) >= len ||
        (dofMap[Nodes[INDEX2(i,e,NN)]]-idRange.first) < 0) {
    printf("BOUNDS_CHECK %s %d i=%d e=%d NN=%d min_id=%d dof=%d\n",
            __FILE__, __LINE__, i, e, NN, idRange.first, dofMap[Nodes[INDEX2(i,e,NN)]]-idRange.first);
    exit(1);
}
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
        for (int e=0; e<numElements; e++) {
            for (int i=0; i<NN; i++) {
                mask[Nodes[INDEX2(lin_nodes[i],e,numNodes)]-offset]=1;
            }
        }
    } else {
        const int NN=refElement->Type->numNodes;
#pragma omp parallel for
        for (int e=0; e<numElements; e++) {
            for (int i=0; i<NN; i++) {
                mask[Nodes[INDEX2(i,e,numNodes)]-offset]=1;
            }
        }
    }
}

void ElementFile::markDOFsConnectedToRange(int* mask, int offset, int marker,
        int firstDOF, int lastDOF, const int *dofIndex, bool useLinear) 
{
    const_ReferenceElement_ptr refElement(referenceElementSet->
                                            borrowReferenceElement(false));
    if (useLinear) {
        const int NN=refElement->numLinearNodes;
        const int *lin_nodes=refElement->Type->linearNodes;
        for (int color=minColor; color<=maxColor; color++) {
#pragma omp parallel for
            for (int e=0; e<numElements; e++) {
                if (Color[e]==color) {
                    for (int i=0; i<NN; i++) {
                        const int k=dofIndex[Nodes[INDEX2(lin_nodes[i],e,numNodes)]];
                        if (firstDOF<=k && k<lastDOF) {
                            for (int j=0; j<NN; j++)
                                mask[dofIndex[Nodes[INDEX2(lin_nodes[j],e,numNodes)]]-offset]=marker;
                            break;
                        }
                    }
                }
            }
        }
    } else {
        const int NN=refElement->Type->numNodes;
        for (int color=minColor; color<=maxColor; color++) {
#pragma omp parallel for
            for (int e=0; e<numElements; e++) {
                if (Color[e]==color) {
                    for (int i=0; i<NN; i++) {
                        const int k=dofIndex[Nodes[INDEX2(i,e,numNodes)]];
                        if (firstDOF<=k && k<lastDOF) {
                            for (int j=0; j<NN; j++)
                                mask[dofIndex[Nodes[INDEX2(j,e,numNodes)]]-offset]=marker;
                            break;
                        }
                    }
                }
            }
        }
    }
}

/// redistributes the elements including overlap by rank
void ElementFile::distributeByRankOfDOF(const std::vector<int>& mpiRankOfDOF, int* index)
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
        std::vector<int> send_count(size);
        std::vector<int> recv_count(size);
        std::vector<int> newOwner(numElements);
#pragma omp parallel
        {
            std::vector<int> loc_send_count(size);
#pragma omp for
            for (int e=0; e<numElements; e++) {
                if (Owner[e] == myRank) {
                    newOwner[e]=myRank;
                    std::vector<int> loc_proc_mask(size);
                    for (int j=0; j<numNodes; j++) {
                        const int p=mpiRankOfDOF[Nodes[INDEX2(j,e,numNodes)]];
                        loc_proc_mask[p]++;
                    }
                    int loc_proc_mask_max=0;
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
        MPI_Alltoall(&send_count[0], 1, MPI_INT, &recv_count[0], 1, MPI_INT,
                     MPIInfo->comm);
        // get the new number of elements for this processor
        int newNumElements=0;
        int numElementsInBuffer=0;
        for (int p=0; p<size; ++p) {
            newNumElements+=recv_count[p];
            numElementsInBuffer+=send_count[p];
        }

        std::vector<int> Id_buffer(numElementsInBuffer);
        std::vector<int> Tag_buffer(numElementsInBuffer);
        std::vector<int> Owner_buffer(numElementsInBuffer);
        std::vector<int> Nodes_buffer(numElementsInBuffer*numNodes);
        std::vector<int> send_offset(size);
        std::vector<int> recv_offset(size);
        std::vector<unsigned char> proc_mask(size);

        // calculate the offsets for the processor buffers
        for (int p=0; p<size-1; ++p) {
            recv_offset[p+1]=recv_offset[p]+recv_count[p];
            send_offset[p+1]=send_offset[p]+send_count[p];
        }

        send_count.assign(size, 0);
        // copy element into buffers. proc_mask makes sure that an element is
        // copied once only for each processor
        for (int e=0; e<numElements; e++) {
            if (Owner[e] == myRank) {
                proc_mask.assign(size, TRUE);
                for (int j=0; j<numNodes; j++) {
                    const int p=mpiRankOfDOF[Nodes[INDEX2(j,e,numNodes)]];
                    if (proc_mask[p]) {
                        int k=send_offset[p]+send_count[p];
                        Id_buffer[k]=Id[e];
                        Tag_buffer[k]=Tag[e];
                        Owner_buffer[k]=newOwner[e];
                        for (int i=0; i<numNodes; i++)
                            Nodes_buffer[INDEX2(i,k,numNodes)]=
                                    index[Nodes[INDEX2(i,e,numNodes)]];
                        send_count[p]++;
                        proc_mask[p]=FALSE;
                    }
                }
            }
        }
        // allocate new tables
        allocTable(newNumElements);

        // start to receive new elements
        for (int p=0; p<size; ++p) {
            if (recv_count[p] > 0) {
                MPI_Irecv(&Id[recv_offset[p]], recv_count[p], MPI_INT, p,
                        MPIInfo->msg_tag_counter+myRank, MPIInfo->comm,
                        &mpi_requests[numRequests]);
                numRequests++;
                MPI_Irecv(&Tag[recv_offset[p]], recv_count[p], MPI_INT, p,
                        MPIInfo->msg_tag_counter+size+myRank, MPIInfo->comm,
                        &mpi_requests[numRequests]);
                numRequests++;
                MPI_Irecv(&Owner[recv_offset[p]], recv_count[p], MPI_INT, p,
                        MPIInfo->msg_tag_counter+2*size+myRank, MPIInfo->comm,
                        &mpi_requests[numRequests]);
                numRequests++;
                MPI_Irecv(&Nodes[recv_offset[p]*numNodes],
                        recv_count[p]*numNodes, MPI_INT, p,
                        MPIInfo->msg_tag_counter+3*size+myRank, MPIInfo->comm,
                        &mpi_requests[numRequests]);
                numRequests++;
            }
        }
        // now the buffers can be sent away
        for (int p=0; p<size; ++p) {
            if (send_count[p] > 0) {
                MPI_Issend(&Id_buffer[send_offset[p]], send_count[p], MPI_INT,
                        p, MPIInfo->msg_tag_counter+p, MPIInfo->comm,
                        &mpi_requests[numRequests]);
                numRequests++;
                MPI_Issend(&Tag_buffer[send_offset[p]], send_count[p], MPI_INT,
                        p, MPIInfo->msg_tag_counter+size+p, MPIInfo->comm,
                        &mpi_requests[numRequests]);
                numRequests++;
                MPI_Issend(&Owner_buffer[send_offset[p]], send_count[p],
                        MPI_INT, p, MPIInfo->msg_tag_counter+2*size+p,
                        MPIInfo->comm, &mpi_requests[numRequests]);
                numRequests++;
                MPI_Issend(&Nodes_buffer[send_offset[p]*numNodes],
                        send_count[p]*numNodes, MPI_INT, p,
                        MPIInfo->msg_tag_counter+3*size+p, MPIInfo->comm,
                        &mpi_requests[numRequests]);
                numRequests++;
            }
        }
        ESYS_MPI_INC_COUNTER(*MPIInfo, 4*size);
        // wait for the requests to be finalized
        MPI_Waitall(numRequests, &mpi_requests[0], &mpi_stati[0]);
#endif
    } else { // single rank
#pragma omp parallel for
        for (int e=0; e<numElements; e++) {
            Owner[e]=0;
            for (int i=0; i<numNodes; i++)
                Nodes[INDEX2(i,e,numNodes)]=index[Nodes[INDEX2(i,e,numNodes)]];
        }
    }
}

} // namespace finley

