
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

#ifndef __FINLEY_ELEMENTFILE_H__
#define __FINLEY_ELEMENTFILE_H__

#include "Finley.h"
#include "NodeFile.h"
#include "ReferenceElementSets.h"
#include "Util.h"
#include "esysUtils/Esys_MPI.h"

namespace finley {

struct ElementFile_Jacobians {
    ElementFile_Jacobians(Finley_ShapeFunction*);
    ~ElementFile_Jacobians();

    /// status of mesh when jacobians were updated last time
    Finley_Status_t status;
    /// number of spatial dimensions
    int numDim;
    /// basis function used
    Finley_ShapeFunction* BasisFunctions;
    /// total number of quadrature nodes used to calculate jacobians
    /// = numSub * BasisFunctions->numQuadNodes
    int numQuadTotal;
    /// number of sides (=1 normal, =2 contact)
    int numSides;
    /// offset to sides (borrowed reference)
    int* offsets;
    /// number of subelements
    int numSub;
    /// total number of shape functions = BasisFunctions->numShapes * numSides
    int numShapesTotal;
    /// local node selection list of length numSub*numShapesTotal
    /// (borrowed reference)
    int* node_selection;
    /// number of elements
    int numElements;
    /// local volume
    double* volume;
    /// derivatives of shape functions in global coordinates at quadrature
    /// points
    double* DSDX;
};

class ElementFile
{
public:
    ElementFile(Finley_ReferenceElementSet* referenceElementSet,
                Esys_MPIInfo *mpiInfo);
    ~ElementFile();

    void allocTable(int numElements);
    void freeTable();

    void distributeByRankOfDOF(int* mpiRankOfDOF, int *Id);
    void createColoring(int nodeCount, int* degreeOfFreedom);
    /// reorders the elements so that they are stored close to the nodes
    void optimizeOrdering();
    /// assigns new node reference numbers to the elements
    void relabelNodes(int* newNode, int offset);
    void markNodes(int* mask, int offset, bool useLinear);
    void scatter(int* index, const ElementFile* in);
    void gather(int* index, const ElementFile* in);
    void copyTable(int offset, int nodeOffset, int idOffset,
                   const ElementFile* in);

    void markDOFsConnectedToRange(int* mask, int offset, int marker,
                                  int firstDOF, int lastDOF, int *dofIndex,
                                  bool useLinear);

    void setTags(const int newTag, const escript::Data& mask);
    ElementFile_Jacobians* borrowJacobians(NodeFile*, bool_t, bool_t);
    /// returns the minimum and maximum reference number of nodes describing
    /// the elements
    inline std::pair<int,int> getNodeRange() const;
    inline void updateTagList();

private:
    void swapTable(ElementFile* other);

public:
    Esys_MPIInfo *MPIInfo;

    /// the reference element to be used
    Finley_ReferenceElementSet *referenceElementSet;
    /// number of elements
    int numElements;
    /// Id[i] is the id number of node i. This number is used when elements
    /// are resorted. In the entire code the term 'element id' refers to i and
    /// not to Id[i] unless explicitly stated otherwise.
    int *Id;
    /// Tag[i] is the tag of element i
    int *Tag;
    /// Owner[i] contains the rank that owns element i
    int *Owner;
    /// array of tags which are actually used
    std::vector<int> tagsInUse;
    /// number of nodes per element
    int numNodes;
    /// Nodes[INDEX(k, i, numNodes)] is the k-th node in the i-th element.
    /// Note that in the way the nodes are ordered Nodes[INDEX(k, i, numNodes)
    /// is the k-th node of element i when referring to the linear version of
    /// the mesh.
    int *Nodes;
    /// minimum color
    int minColor;
    /// maximum color
    int maxColor;
    /// assigns each element a color. Elements with the same color don't share
    /// a node so they can be processed simultaneously.
    /// At any time Color must provide a valid value. In any case one can set
    /// Color[e]=e for all e
    int *Color;
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

inline std::pair<int,int> ElementFile::getNodeRange() const
{
    return util::getMinMaxInt(numNodes, numElements, Nodes);
}


inline void ElementFile::updateTagList()
{
    util::setValuesInUse(Tag, numElements, tagsInUse, MPIInfo);
}

} // namespace finley

#endif // __FINLEY_ELEMENTFILE_H__

