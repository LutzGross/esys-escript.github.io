
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

#ifndef __DUDLEY_ELEMENTFILE_H__
#define __DUDLEY_ELEMENTFILE_H__

#include "Dudley.h"
#include "NodeFile.h"
#include "ElementType.h"
#include "Util.h"

namespace dudley {

struct ElementFile_Jacobians
{
    ElementFile_Jacobians();
    ~ElementFile_Jacobians();

    /// status of mesh when jacobians were updated last time
    int status;
    /// number of spatial dimensions
    int numDim;
    /// number of quadrature nodes used to calculate jacobians
    int numQuad;
    /// number of shape functions
    int numShapes;
    /// number of elements
    dim_t numElements;
    /// used to compute volume
    double *absD;
    /// used to compute volume
    double quadweight;
    /// derivatives of shape functions in global coordinates at quadrature
    /// points
    double* DSDX;
};

class ElementFile
{
public:
    ElementFile(ElementTypeId etype, escript::JMPI mpiInfo);
    ~ElementFile();

    /// allocates the element table within an element file to hold NE elements
    void allocTable(dim_t NE);

    /// deallocates the element table within an element file
    void freeTable();

    /// copies element file `in` into this element file starting from `offset`.
    /// The elements `offset` to in->numElements+offset-1 will be overwritten.
    void copyTable(index_t offset, index_t nodeOffset, index_t idOffset,
                   const ElementFile* in);

    /// prints information about this element file to stdout
    void print(const index_t* nodesId) const;

    /// redistributes the elements including overlap by rank
    void distributeByRankOfDOF(const int* mpiRankOfDOF,
                               const index_t* nodesId);

    /// Tries to reduce the number of colors used to color elements in this
    /// ElementFile
    void createColoring(dim_t numNodes, const index_t* degreeOfFreedom);

    /// reorders the elements so that they are stored close to the nodes
    void optimizeOrdering();

    /// assigns new node reference numbers to the elements.
    /// If k is the old node, the new node is newNode[k-offset].
    void relabelNodes(const index_t* newNode, index_t offset);

    void markNodes(std::vector<short>& mask, index_t offset) const;

    /// gathers the elements from the element file `in` using
    /// index[0:out->elements-1]. `index` has to be between 0 and
    /// in->numElements-1. A conservative assumption on the colouring is made.
    void gather(const index_t* index, const ElementFile* in);

    /// sets element tags to newTag where mask > 0
    void setTags(int newTag, const escript::Data& mask);

    ElementFile_Jacobians* borrowJacobians(const NodeFile* nodes,
                                           bool reducedOrder) const;

    /// returns the minimum and maximum reference number of nodes
    /// describing the elements
    inline std::pair<index_t,index_t> getNodeRange() const;

    inline void updateTagList();

private:
    void swapTable(ElementFile* other);

public:
    escript::JMPI MPIInfo;

    /// number of elements
    dim_t numElements;

    /// Id[i] is the id number of node i. This number is used when elements
    /// are resorted. In the entire code the term 'element id' refers to i and
    /// not to Id[i] unless explicitly stated otherwise.
    index_t* Id;

    /// Tag[i] is the tag of element i
    int* Tag;

    /// Owner[i] contains the rank that owns element i
    int* Owner;

    /// array of tags which are actually used
    std::vector<int> tagsInUse;

    /// number of nodes per element
    int numNodes;

    /// Nodes[INDEX(k, i, numNodes)] is the k-th node in the i-th element.
    index_t* Nodes;

    /// assigns each element a color. Elements with the same color don't share
    /// a node so they can be processed simultaneously. At any time Color must
    /// provide a valid value. In any case one can set Color[e]=e for all e.
    index_t* Color;

    /// minimum color value
    index_t minColor;

    /// maximum color value
    index_t maxColor;

    /// number of spatial dimensions of the domain
    int numDim;

    /// dimension of the element e.g. 2 for a line in 2D or 3D
    int numLocalDim;

    /// element type ID
    ElementTypeId etype;

    /// name of element type
    const char *ename;

    /// number of shape functions
    int numShapes;

private:
    /// jacobians of the shape function used for solution approximation
    ElementFile_Jacobians* jacobians;

    /// jacobians of the shape function used for solution approximation for
    /// reduced integration order
    ElementFile_Jacobians* jacobians_reducedQ;
};

inline std::pair<index_t,index_t> ElementFile::getNodeRange() const
{
    return util::getMinMaxInt(numNodes, numElements, Nodes);
}

inline void ElementFile::updateTagList()
{
    util::setValuesInUse(Tag, numElements, tagsInUse, MPIInfo);
}


} // namespace dudley

#endif // __DUDLEY_ELEMENTFILE_H__

