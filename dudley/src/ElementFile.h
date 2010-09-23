
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#ifndef INC_DUDLEY_ELEMENTFILE
#define INC_DUDLEY_ELEMENTFILE

#include "Dudley.h"
#include "NodeFile.h"
#include "ReferenceElementSets.h"
#include "escript/DataC.h"

#ifdef PASO_MPI
#include "paso/Paso_MPI.h"
#endif


struct Dudley_ElementFile_Jacobeans {
  Dudley_Status_t status;               /* status of mesh when jacobeans where updated last time */
  dim_t numDim;                         /* spatial dimension */
  Dudley_ShapeFunction* BasisFunctions; /* basis function used */
  dim_t numQuadTotal;           /* total number of quadrature nodes used to calculate jacobeans = BasisFunctions->numQuadNodes*/

  dim_t numShapes;
//  dim_t numShapesTotal;         /* total number of shape functions =  BasisFunctions->numShapes */
  dim_t numElements;            /* number of elements */
  double* absD;			/* used to compute volume */
  double quadweight;		/* used to compute volume */
  double* DSDX;                         /* derivatives of shape functions in global coordinates at quadrature points*/
};

typedef struct Dudley_ElementFile_Jacobeans Dudley_ElementFile_Jacobeans;

struct Dudley_ElementFile {
  Paso_MPIInfo *MPIInfo;
  Paso_MPI_rank *Owner;

  Dudley_ReferenceElementSet *referenceElementSet; /* the reference element to be used */

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

  Dudley_ElementFile_Jacobeans* jacobeans;           /* jacobeans of the shape function used for solution approximation */
  Dudley_ElementFile_Jacobeans* jacobeans_reducedQ;  /* jacobeans of the shape function used for solution approximation for reduced integration order*/
  dim_t numDim;
};

typedef struct Dudley_ElementFile Dudley_ElementFile;
Dudley_ElementFile* Dudley_ElementFile_alloc(Dudley_ReferenceElementSet* referenceElementSet, Paso_MPIInfo *MPIInfo);
void Dudley_ElementFile_free(Dudley_ElementFile*);
void Dudley_ElementFile_allocTable(Dudley_ElementFile*,dim_t);
void Dudley_ElementFile_freeTable(Dudley_ElementFile*);
void Dudley_ElementFile_setElementDistribution(Dudley_ElementFile* in, dim_t* distribution);
dim_t Dudley_ElementFile_getGlobalNumElements(Dudley_ElementFile* in);
dim_t Dudley_ElementFile_getMyNumElements(Dudley_ElementFile* in);
index_t Dudley_ElementFile_getFirstElement(Dudley_ElementFile* in); 
void Dudley_ElementFile_distributeByRankOfDOF(Dudley_ElementFile* self, Paso_MPI_rank* mpiRankOfDOF, index_t *Id);

void Dudley_ElementFile_createColoring(Dudley_ElementFile* in,dim_t numNodes,dim_t* degreeOfFreedom);
void Dudley_ElementFile_optimizeOrdering(Dudley_ElementFile** in);
void Dudley_ElementFile_setNodeRange(dim_t*,dim_t*,Dudley_ElementFile*);
void Dudley_ElementFile_relableNodes(dim_t*,dim_t,Dudley_ElementFile*);
void Dudley_ElementFile_markNodes(dim_t*,dim_t,dim_t,Dudley_ElementFile*,dim_t);
void Dudley_ElementFile_scatter(dim_t*,Dudley_ElementFile*,Dudley_ElementFile*);
void Dudley_ElementFile_gather(dim_t*,Dudley_ElementFile*,Dudley_ElementFile*);
void Dudley_ElementFile_copyTable(dim_t,Dudley_ElementFile*,dim_t,dim_t,Dudley_ElementFile*);
void Dudley_ElementFile_markDOFsConnectedToRange(index_t* mask,index_t offset,index_t marker,index_t firstDOF,index_t lastDOF,index_t *dofIndex,Dudley_ElementFile*in ,bool_t useLinear);

void Dudley_ElementFile_setTags(Dudley_ElementFile* ,const int, escriptDataC*);
Dudley_ElementFile_Jacobeans* Dudley_ElementFile_Jacobeans_alloc(Dudley_ShapeFunction* );
void Dudley_ElementFile_Jacobeans_dealloc(Dudley_ElementFile_Jacobeans*);
Dudley_ElementFile_Jacobeans* Dudley_ElementFile_borrowJacobeans(Dudley_ElementFile*, Dudley_NodeFile*, bool_t, bool_t);
void Dudley_ElementFile_setTagsInUse(Dudley_ElementFile* in);


#endif /* #ifndef INC_DUDLEY_ELEMENTFILE */

