
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

#include "ElementFile.h"
#include "ShapeTable.h"

#include <escript/index.h>

namespace dudley {

ElementFile::ElementFile(ElementTypeId type, escript::JMPI mpiInfo) :
    MPIInfo(mpiInfo),
    numElements(0),
    Id(NULL),
    Tag(NULL),
    Owner(NULL),
    Nodes(NULL),
    Color(NULL),
    minColor(0),
    maxColor(-1),
    etype(type)
{
    jacobians = new ElementFile_Jacobians();
    jacobians_reducedQ = new ElementFile_Jacobians();

    numDim = Dims[type];
    numNodes = numDim + 1;
    numLocalDim = localDims[type];
    numShapes = numLocalDim + 1;
    ename = getElementName(type);
}

ElementFile::~ElementFile()
{
    freeTable();
    delete jacobians;
    delete jacobians_reducedQ;
}

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
    maxColor = -1;
    minColor = 0;
}

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

void ElementFile::copyTable(index_t offset, index_t nodeOffset,
                            index_t idOffset, const ElementFile* in)
{
    const int NN_in = in->numNodes;
    if (NN_in > numNodes) {
        throw DudleyException("ElementFile::copyTable: dimensions of element files don't match.");
    }

    if (MPIInfo->comm != in->MPIInfo->comm) {
        throw DudleyException("ElementFile::copyTable: MPI communicators of element files don't match.");
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

void ElementFile::print(const index_t* nodesId) const
{
    std::cout << "=== " << ename << ":\nnumber of elements=" << numElements
              << "\ncolor range=[" << minColor << "," << maxColor << "]\n";

    if (numElements > 0) {
        std::cout << "Id,Tag,Owner,Color,Nodes" << std::endl;
        for (index_t i = 0; i < numElements; i++) {
            std::cout << Id[i] << "," << Tag[i] << ","
                << Owner[i] << "," << Color[i] << ",";
            for (int j = 0; j < numNodes; j++)
                std::cout << " " << nodesId[Nodes[INDEX2(j, i, numNodes)]];
            std::cout << std::endl;
        }
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
    minColor = std::min(minColor, in->minColor + maxColor + 1);
    maxColor = std::max(maxColor, in->maxColor + maxColor + 1);
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

    util::ValueAndIndexList item_list(numElements);
    index_t* index = new index_t[numElements];
    ElementFile* out = new ElementFile(etype, MPIInfo);
    out->allocTable(numElements);

#pragma omp parallel for
    for (index_t e = 0; e < numElements; e++) {
        std::pair<index_t,index_t> entry(Nodes[INDEX2(0, e, numNodes)], e);
        for (int i = 1; i < numNodes; i++)
            entry.first = std::min(entry.first, Nodes[INDEX2(i, e, numNodes)]);
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

void ElementFile::setTags(int newTag, const escript::Data& mask)
{
    if (mask.isComplex())
    {
      throw DudleyException("ElementFile::setTags: mask argument must not be complex.");
    }
    const int numQuad = hasReducedIntegrationOrder(mask) ? 1 : numNodes;

    if (1 != mask.getDataPointSize()) {
        throw DudleyException("ElementFile::setTags: number of components of mask must be 1.");
    } else if (!mask.numSamplesEqual(numQuad, numElements)) {
        throw DudleyException("ElementFile::setTags: illegal number of samples of mask Data object");
    }

    escript::DataTypes::real_t wantreal=0;
    if (mask.actsExpanded()) {
#pragma omp parallel for
        for (index_t n = 0; n < numElements; n++) {
            if (mask.getSampleDataRO(n, wantreal)[0] > 0)
                Tag[n] = newTag;
        }
    } else {
#pragma omp parallel for
        for (index_t n = 0; n < numElements; n++) {
            const double* mask_array = mask.getSampleDataRO(n, wantreal);
            bool check = false;
            for (int q = 0; q < numQuad; q++)
                check = check || mask_array[q];
            if (check)
                Tag[n] = newTag;
        }
    }
    updateTagList();
}

void ElementFile::markNodes(std::vector<short>& mask, index_t offset) const
{
#pragma omp parallel for
    for (index_t e = 0; e < numElements; e++) {
        for (int i = 0; i < numNodes; i++) {
            mask[Nodes[INDEX2(i, e, numNodes)] - offset] = 1;
        }
    }
}

void ElementFile::relabelNodes(const index_t* newNode, index_t offset)
{
#pragma omp parallel for
    for (index_t j = 0; j < numElements; j++) {
        for (int i = 0; i < numNodes; i++) {
            Nodes[INDEX2(i, j, numNodes)] =
                              newNode[Nodes[INDEX2(i, j, numNodes)] - offset];
        }
    }
}

} // namespace dudley

