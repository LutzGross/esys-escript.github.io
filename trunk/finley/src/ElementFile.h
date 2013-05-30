
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


#ifndef INC_FINLEY_ELEMENTFILE
#define INC_FINLEY_ELEMENTFILE

#include "Finley.h"
#include "NodeFile.h"
#include "ReferenceElementSets.h"
#include "escript/DataC.h"

#ifdef ESYS_MPI
#include "esysUtils/Esys_MPI.h"
#endif

#include <vector>

struct Finley_ElementFile_Jacobians {
  Finley_Status_t status;               /* status of mesh when jacobians were updated last time */
  dim_t numDim;                         /* spatial dimension */
  Finley_ShapeFunction* BasisFunctions; /* basis function used */
  dim_t numQuadTotal;           /* total number of quadrature nodes used to calculate jacobians = numSub * BasisFunctions->numQuadNodes*/
  dim_t numSides;                   /* number of sides (=1 normal, =2 contact) */
  index_t* offsets;         /* offset to sides (borrowed reference) */
  dim_t numSub;         /* number of subelements        */
  dim_t numShapesTotal;         /* total number of shape functions =  BasisFunctions->numShapes * numSides */
  index_t* node_selection;      /* local node selection list of length numSub * numShapesTotal  (borrowed reference)  */
  dim_t numElements;            /* number of elements */
  double* volume;                       /* local volume */
  double* DSDX;                         /* derivatives of shape functions in global coordinates at quadrature points*/
};

typedef struct Finley_ElementFile_Jacobians Finley_ElementFile_Jacobians;

struct Finley_ElementFile {
  Esys_MPIInfo *MPIInfo;
  Esys_MPI_rank *Owner;

  Finley_ReferenceElementSet *referenceElementSet; /* the reference element to be used */

  dim_t numElements;                             /* number of elements. */
  
  index_t *Id;                                 /* Id[i] is the id number of
                                                node i. this number is not
                                                used but useful when
                                                elements are resorted. In
                                                the entire code the term
                                                'element id' refers to i
                                                and not to Id[i] unless
                                                explicitly stated
                                                otherwise. */

  index_t *Tag;                                /* Tag[i] is the tag of element i. */

  std::vector<int> tagsInUse;                  /* array of tags which are actually used */

  dim_t numNodes;                              /* number of nodes per element */
  index_t *Nodes;                              /* Nodes[INDEX(k, i, numNodes)]
                                                is the k-the node in the
                                                i-the element. Note that
                                                in the way the nodes are
                                                ordered Nodes[INDEX(k, i, numNodes)
                                                is the k-th node of element i
                                                when referring to the
                                                linear version of the mesh. */
  index_t minColor;                           /* minimum color */
  index_t maxColor;                           /* maximum color */
  index_t *Color;                             /* assigns each element a color. Elements with the same color     
                                                 don't share a node so they can be processed simultaneously.
                                                 At anytime Color must provide a valid value. In any case one can set  
                                                 Color[e]=e for all e */

  Finley_ElementFile_Jacobians* jacobians;           /* jacobians of the shape function used for solution approximation */
  Finley_ElementFile_Jacobians* jacobians_reducedS;  /* jacobians of the shape function used for solution approximation for reduced order of shape function*/
  Finley_ElementFile_Jacobians* jacobians_reducedQ;  /* jacobians of the shape function used for solution approximation for reduced integration order*/
  Finley_ElementFile_Jacobians* jacobians_reducedS_reducedQ;  /* jacobians of the shape function used for solution approximation for reduced integration order and reduced order of shape function*/

};

typedef struct Finley_ElementFile Finley_ElementFile;
Finley_ElementFile* Finley_ElementFile_alloc(Finley_ReferenceElementSet* referenceElementSet, Esys_MPIInfo *MPIInfo);
void Finley_ElementFile_free(Finley_ElementFile*);
void Finley_ElementFile_allocTable(Finley_ElementFile*,dim_t);
void Finley_ElementFile_freeTable(Finley_ElementFile*);
void Finley_ElementFile_setElementDistribution(Finley_ElementFile* in, dim_t* distribution);
dim_t Finley_ElementFile_getGlobalNumElements(Finley_ElementFile* in);
dim_t Finley_ElementFile_getMyNumElements(Finley_ElementFile* in);
index_t Finley_ElementFile_getFirstElement(Finley_ElementFile* in); 
void Finley_ElementFile_distributeByRankOfDOF(Finley_ElementFile* self, Esys_MPI_rank* mpiRankOfDOF, index_t *Id);

void Finley_ElementFile_createColoring(Finley_ElementFile* in,dim_t numNodes,dim_t* degreeOfFreedom);
void Finley_ElementFile_optimizeOrdering(Finley_ElementFile** in);
void Finley_ElementFile_setNodeRange(dim_t*,dim_t*,Finley_ElementFile*);
void Finley_ElementFile_relableNodes(dim_t*,dim_t,Finley_ElementFile*);
void Finley_ElementFile_markNodes(dim_t*,dim_t,dim_t,Finley_ElementFile*,dim_t);
void Finley_ElementFile_scatter(dim_t*,Finley_ElementFile*,Finley_ElementFile*);
void Finley_ElementFile_gather(dim_t*,Finley_ElementFile*,Finley_ElementFile*);
void Finley_ElementFile_copyTable(dim_t,Finley_ElementFile*,dim_t,dim_t,Finley_ElementFile*);
void Finley_ElementFile_markDOFsConnectedToRange(index_t* mask,index_t offset,index_t marker,index_t firstDOF,index_t lastDOF,index_t *dofIndex,Finley_ElementFile*in ,bool_t useLinear);

void Finley_ElementFile_setTags(Finley_ElementFile* ,const int, escriptDataC*);
Finley_ElementFile_Jacobians* Finley_ElementFile_Jacobians_alloc(Finley_ShapeFunction* );
void Finley_ElementFile_Jacobians_dealloc(Finley_ElementFile_Jacobians*);
Finley_ElementFile_Jacobians* Finley_ElementFile_borrowJacobians(Finley_ElementFile*, finley::NodeFile*, bool_t, bool_t);
void Finley_ElementFile_setTagsInUse(Finley_ElementFile* in);


#endif /* #ifndef INC_FINLEY_ELEMENTFILE */

