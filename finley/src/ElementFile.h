
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __FINLEY_ELEMENTFILE_H__
#define __FINLEY_ELEMENTFILE_H__

#include "Finley.h"
#include "NodeFile.h"
#include "ReferenceElementSets.h"
#include "Util.h"

#ifdef ESYS_HAVE_HDF5
  #include <H5Cpp.h>
#endif

namespace finley {

struct ElementFile_Jacobians
{
    ElementFile_Jacobians(const_ShapeFunction_ptr basis);
    ~ElementFile_Jacobians();

    /// status of mesh when jacobians were updated last time
    int status;
    /// number of spatial dimensions
    int numDim;
    /// basis function used
    const_ShapeFunction_ptr BasisFunctions;
    /// total number of quadrature nodes used to calculate jacobians
    /// = numSub * BasisFunctions->numQuadNodes
    int numQuadTotal;
    /// number of sides (=1 normal, =2 contact)
    int numSides;
    /// offset to sides (borrowed reference)
    const int* offsets;
    /// number of subelements
    int numSub;
    /// total number of shape functions = BasisFunctions->numShapes * numSides
    int numShapesTotal;
    /// local node selection list of length numSub*numShapesTotal
    /// (borrowed reference)
    const int* node_selection;
    /// number of elements
    dim_t numElements;
    /// local volume
    double* volume;
    /// derivatives of shape functions in global coordinates at quadrature
    /// points
    double* DSDX;
};

class ElementFile
{
public:
    ElementFile(const_ReferenceElementSet_ptr refElementSet,
                escript::JMPI mpiInfo);
    ~ElementFile();

    /// allocates the element table within an element file to hold NE elements
    void allocTable(dim_t NE);

    /// deallocates the element table within an element file
    void freeTable();

    /// copies element file `in` into this element file starting from `offset`.
    /// The elements `offset` to in->numElements+offset-1 will be overwritten.
    void copyTable(index_t offset, index_t nodeOffset, index_t idOffset,
                   const ElementFile* in);

    /// redistributes the elements including overlap by rank
    void distributeByRankOfDOF(const std::vector<int>& mpiRankOfDOF,
                               index_t* nodesId);

    /// Tries to reduce the number of colors used to color elements in this
    /// ElementFile
    void createColoring(const IndexVector& dofMap);

    /// reorders the elements so that they are stored close to the nodes
    void optimizeOrdering();

    /// assigns new node reference numbers to the elements.
    /// If k is the old node, the new node is newNode[k-offset].
    void relabelNodes(const IndexVector& newNode, index_t offset);

    void markNodes(std::vector<short>& mask, int offset, bool useLinear);

    void gather(const index_t* index, const ElementFile* in);

    void scatter(index_t* index, const ElementFile* in);

    void setTags(const int newTag, const escript::Data& mask);

    ElementFile_Jacobians* borrowJacobians(const NodeFile*, bool, bool) const;

    /// returns the minimum and maximum reference number of nodes describing
    /// the elements
    inline std::pair<index_t,index_t> getNodeRange() const;

    /// updates the list of tags in use. This method must be called by all
    /// ranks.
    inline void updateTagList();

    /// dump element data to HDF5 file
    #ifdef ESYS_HAVE_HDF5
    void dump(H5::Group h5_grp) const;
    #endif
private:
    void swapTable(ElementFile* other);

public:
    escript::JMPI MPIInfo;

    /// the reference element to be used
    const_ReferenceElementSet_ptr referenceElementSet;
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
    /// Note that in the way the nodes are ordered Nodes[INDEX(k, i, numNodes)
    /// is the k-th node of element i when referring to the linear version of
    /// the mesh.
    index_t* Nodes;

    /// assigns each element a color. Elements with the same color don't share
    /// a node so they can be processed simultaneously.
    /// At any time Color must provide a valid value. In any case one can set
    /// Color[e]=e for all e
    index_t* Color;

    /// minimum color value
    index_t minColor;

    /// maximum color value
    index_t maxColor;

    /// jacobians of the shape function used for solution approximation
    ElementFile_Jacobians* jacobians;

    /// jacobians of the shape function used for solution approximation for
    /// reduced order of shape function
    ElementFile_Jacobians* jacobians_reducedS;

    /// jacobians of the shape function used for solution approximation for
    /// reduced integration order
    ElementFile_Jacobians* jacobians_reducedQ;

    /// jacobians of the shape function used for solution approximation for
    /// reduced integration order and reduced order of shape function
    ElementFile_Jacobians* jacobians_reducedS_reducedQ;
};

inline std::pair<index_t,index_t> ElementFile::getNodeRange() const
{
    return util::getMinMaxInt(numNodes, numElements, Nodes);
}

inline void ElementFile::updateTagList()
{
    util::setValuesInUse(Tag, numElements, tagsInUse, MPIInfo);
}

#ifdef ESYS_HAVE_HDF5
FINLEY_DLL_API
ElementFile* loadElements_hdf5(const H5::Group h5_grp,
                                const int integration_order,
                                const int reduced_integration_order,
                                const escript::JMPI mpiInfo);
#endif // ESYS_HAVE_HDF5

} // namespace finley

#endif // __FINLEY_ELEMENTFILE_H__

