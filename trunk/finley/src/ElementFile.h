
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#ifndef INC_FINLEY_ELEMENTFILE
#define INC_FINLEY_ELEMENTFILE

#include "Finley.h"
#include "NodeFile.h"
#include "ReferenceElementSets.h"
#include "escript/DataC.h"

#ifdef PASO_MPI
#include "paso/Paso_MPI.h"
#endif


struct Finley_ElementFile_Jacobeans {
  Finley_Status_t status;               /* status of mesh when jacobeans where updated last time */
  dim_t numDim;                         /* spatial dimension */
  Finley_ShapeFunction* BasisFunctions; /* basis function used */
  dim_t numQuadTotal;           /* total number of quadrature nodes used to calculate jacobeans = numSub * BasisFunctions->numQuadNodes*/
  dim_t numSides;                   /* number of sides (=1 normal, =2 contact) */
  index_t* offsets;         /* offset to sides (borrowed reference) */
  dim_t numSub;         /* number of subelements        */
  dim_t numShapesTotal;         /* total number of shape functions =  BasisFunctions->numShapes * numSides */
  index_t* node_selection;      /* local node selection list of length numSub * numShapesTotal  (borrowed reference)  */
  dim_t numElements;            /* number of elements */
  double* volume;                       /* local volume */
  double* DSDX;                         /* derivatives of shape functions in global coordinates at quadrature points*/
};

typedef struct Finley_ElementFile_Jacobeans Finley_ElementFile_Jacobeans;

struct Finley_ElementFile {
  Paso_MPIInfo *MPIInfo;
  Paso_MPI_rank *Owner;

  Finley_ReferenceElementSet *referenceElementSet; /* the reference element to be used */

  dim_t numElements;                             /* number of elements. */
  
  index_t *Id;                                 /* Id[i] is the id nmber of
                                                node i. this number is not
                                                used but useful when
                                                elements are resorted. in
                                                the entire code the term
                                                'element id' refers to i
                                                but nor to Id[i] if not
                                                explicitly stated
                                                otherwise. */

  index_t *Tag;                                /* Tag[i] is the tag of element i. */

  index_t *tagsInUse;                  /* array of tags which are actually used */
  dim_t     numTagsInUse;               /* number of tags used */


  dim_t numNodes;                              /* number of nodes per element */
  index_t *Nodes;                              /* Nodes[INDEX(k, i, numNodes)]
                                                is the k-the node in the
                                                i-the element. note that
                                                in the way the nodes are
                                                ordered Nodes[INDEX(k, i, numNodes)
                                                is k-the node of element i
                                                when refering to the
                                                linear version of the
                                                mesh. */
  index_t minColor;                           /* minimum color */
  index_t maxColor;                           /* maximum color */
  index_t *Color;                             /* assigns each element a color. elements with the same color     
                                                 are don't share a node so they can be processed simultaneously 
                                                 at anytime Color must provide a valid value. In any case one can set  
                                                 Color[e]=e  for all e */
  index_t order;                                /* order of the element integration scheme*/
  index_t reduced_order;                        /* order of the reduced element integration scheme*/

  Finley_ElementFile_Jacobeans* jacobeans;           /* jacobeans of the shape function used for solution approximation */
  Finley_ElementFile_Jacobeans* jacobeans_reducedS;  /* jacobeans of the shape function used for solution approximation for reduced order of shape function*/
  Finley_ElementFile_Jacobeans* jacobeans_reducedQ;  /* jacobeans of the shape function used for solution approximation for reduced integration order*/
  Finley_ElementFile_Jacobeans* jacobeans_reducedS_reducedQ;  /* jacobeans of the shape function used for solution approximation for reduced integration order and  reduced order of shape function*/

};

typedef struct Finley_ElementFile Finley_ElementFile;
Finley_ElementFile* Finley_ElementFile_alloc(Finley_ReferenceElementSet* referenceElementSet, Paso_MPIInfo *MPIInfo);
void Finley_ElementFile_free(Finley_ElementFile*);
void Finley_ElementFile_allocTable(Finley_ElementFile*,dim_t);
void Finley_ElementFile_freeTable(Finley_ElementFile*);
void Finley_ElementFile_setElementDistribution(Finley_ElementFile* in, dim_t* distribution);
dim_t Finley_ElementFile_getGlobalNumElements(Finley_ElementFile* in);
dim_t Finley_ElementFile_getMyNumElements(Finley_ElementFile* in);
index_t Finley_ElementFile_getFirstElement(Finley_ElementFile* in); 
void Finley_ElementFile_distributeByRankOfDOF(Finley_ElementFile* self, Paso_MPI_rank* mpiRankOfDOF, index_t *Id);

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
Finley_ElementFile_Jacobeans* Finley_ElementFile_Jacobeans_alloc(Finley_ShapeFunction* );
void Finley_ElementFile_Jacobeans_dealloc(Finley_ElementFile_Jacobeans*);
Finley_ElementFile_Jacobeans* Finley_ElementFile_borrowJacobeans(Finley_ElementFile*, Finley_NodeFile*, bool_t, bool_t);
void Finley_ElementFile_setTagsInUse(Finley_ElementFile* in);


#endif /* #ifndef INC_FINLEY_ELEMENTFILE */

