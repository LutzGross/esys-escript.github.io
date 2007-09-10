/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/
/* Version: $Id$ */


#ifndef INC_FINLEY_MESH
#define INC_FINLEY_MESH

/**************************************************************/

/* Finley: Mesh */

/* A mesh is built from nodes and elements which are describing the
   domain, the surface and point sources. (the latter are needed to
   establish links with other codes, in particular to particle
   codes). The nodes are stored a Finley_NodeFile and elements in a
   Finley_ElementFile. A Finley_NodeFile and three Finley_ElementFile
   containing the elements describing the domain, surface and point
   sources respectively. Notice that the surface elements do not
   necessaryly cover the entire surface of the domain. */

/* The element type is fixed by the reference element, see
   ReferenceElement.h. The numbering of the nodes starts with 0. */

/* Important: it is assumed that every node is appearing in at least
   one element or surface element and that any node used in an
   element, surface element or as a point is specified in the
   Finley_Node, see also Finley_resolveNodeIds. */

/* In some cases it is useful to refer to a mesh entirly built from
   order 1 (=linear) elements. The linear version of the mesh can be
   accessed by referning to the first few nodes of each element
   (thanks to the way the nodes are ordered). As the numbering of
   these nodes is not continuous a relabeling vectors are introduced
   in the Finley_NodeFile. This feature is not fully implemented
   yet. */

/* allnodes and elements are tagged. the tag allows to group nodes and
   elements. A typical application is to mark surface elements on a
   certain portion of the domain with the same tag. All these surface
   elements can then assigned the same value eg. for the pressure. */

/* Thespacial dimension is determined by the type of elements
   used. The spacial dimension should be accessed by the function
   Finley_Mesh_getDim. Notice that the element type also determines
   the type of surface elements to be used. */

/**************************************************************/

#include "Finley.h"
#include "NodeFile.h"
#include "ElementFile.h"
#include "TagMap.h"
#include "Util.h"
#include "paso/SystemMatrixPattern.h"
#include "escript/DataC.h"

#ifdef PASO_MPI
#include "paso/Paso_MPI.h"
#endif

/**************************************************************/

/*  this struct holds a mesh: */

struct Finley_Mesh {
  char* Name;                           /* the name of the mesh */
  index_t order;                        /* integration order */
  index_t reduced_order;                /* reduced integration order */
  dim_t reference_counter;              /* counts the number of references to the mesh; */
  Finley_NodeFile* Nodes;               /* the table of the nodes */
  Finley_ElementFile* Elements;         /* the table of the elements */
  Finley_ElementFile* FaceElements;     /* the table of the face elements */
  Finley_ElementFile* ContactElements;  /* the table of the contact elements */
  Finley_ElementFile* Points;           /* the table of points (treated as elements of dimension 0) */
  Finley_TagMap* TagMap;                /* the tag map mapping names to tag keys */

  /* pointer to the sparse matrix pattern */

  Paso_SystemMatrixPattern *FullFullPattern;
  Paso_SystemMatrixPattern *FullReducedPattern;
  Paso_SystemMatrixPattern *ReducedFullPattern;
  Paso_SystemMatrixPattern *ReducedReducedPattern;
  Paso_MPIInfo *MPIInfo;
};

typedef struct Finley_Mesh Finley_Mesh;

/* these structures are used for matching surfaces elements: */

struct Finley_Mesh_findMatchingFaces_center{
   index_t refId;
   double x[MAX_numDim];
};
typedef struct Finley_Mesh_findMatchingFaces_center Finley_Mesh_findMatchingFaces_center;

/**************************************************************/

/*  interfaces: */
Finley_Mesh* Finley_Mesh_alloc(char* name,dim_t numDim, index_t order, index_t reduced_order, Paso_MPIInfo *mpi_info);
Finley_Mesh* Finley_Mesh_reference(Finley_Mesh*);
dim_t Finley_Mesh_getDim(Finley_Mesh*);
void Finley_Mesh_free(Finley_Mesh*);

void Finley_Mesh_addTagMap(Finley_Mesh *mesh_p,const char* name, index_t tag_key);
index_t Finley_Mesh_getTag(Finley_Mesh *mesh_p,const char* name);
bool_t Finley_Mesh_isValidTagName(Finley_Mesh *mesh_p,const char* name);
void Finley_Mesh_distributeByRankOfDOF(Finley_Mesh* in, dim_t *distribution);
Paso_SystemMatrixPattern* Finley_getPattern(Finley_Mesh *mesh,bool_t reduce_row_order, bool_t reduce_col_order);
Paso_SystemMatrixPattern* Finley_makePattern(Finley_Mesh *mesh,bool_t reduce_row_order, bool_t reduce_col_order);
void Finley_Mesh_write(Finley_Mesh*,char*);
void Finley_Mesh_dump(Finley_Mesh *in,char* fname);
Finley_Mesh* Finley_Mesh_load(char* fname);
Finley_Mesh* Finley_Mesh_read(char*,index_t, index_t, bool_t);
Finley_Mesh* Finley_Mesh_readGmsh(char*,index_t, index_t, index_t, bool_t);
void Finley_Mesh_setCoordinates(Finley_Mesh*,escriptDataC*);
void Finley_Mesh_setElements(Finley_Mesh* self,Finley_ElementFile *elements);
void Finley_Mesh_setFaceElements(Finley_Mesh* self,Finley_ElementFile *elements);
void Finley_Mesh_setContactElements(Finley_Mesh* self,Finley_ElementFile *elements);
void Finley_Mesh_setPoints(Finley_Mesh* self,Finley_ElementFile *elements);
void Finley_Mesh_optimizeDOFDistribution(Finley_Mesh* in,dim_t *distribution);
void Finley_Mesh_prepare(Finley_Mesh* in, bool_t optimize);
void Finley_Mesh_createColoring(Finley_Mesh* in, index_t *node_localDOF_map);
void Finley_Mesh_optimizeElementOrdering(Finley_Mesh* in);
void Finley_Mesh_resolveNodeIds(Finley_Mesh*);
void Finley_Mesh_createMappings(Finley_Mesh* in, index_t *distribution);
void Finley_Mesh_createNodeFileMappings(Finley_Mesh* in, dim_t numReducedNodes, index_t* indexReducedNodes, index_t* dof_first_component);
void Finley_Mesh_markDOFsConnectedToRange(index_t* mask, index_t offset, index_t marker,index_t firstDOF,index_t lastDOF,Finley_Mesh* in, bool_t useLinear);

void Finley_Mesh_optimizeDOFLabeling(Finley_Mesh*,dim_t *);


Finley_Mesh* Finley_Mesh_merge(dim_t, Finley_Mesh**);

void Finley_Mesh_relableElementNodes(int*,int,Finley_Mesh*);
void Finley_Mesh_markNodes(int*,int,Finley_Mesh*,int);

void Finley_Mesh_glueFaces(Finley_Mesh* self,double safety_factor,double tolerance, bool_t);
void Finley_Mesh_joinFaces(Finley_Mesh* self,double safety_factor,double tolerance, bool_t);

int Finley_Mesh_findMatchingFaces_compar(const void*,const void*);
void Finley_Mesh_findMatchingFaces(Finley_NodeFile*,Finley_ElementFile *,double,double, int*, int*,int*,int*);
void Finley_Mesh_print(Finley_Mesh *in);
void Finley_Mesh_saveDX(const char * filename_p, Finley_Mesh *mesh_p, const dim_t num_data,char* *names_p,escriptDataC* *data_pp);
void Finley_Mesh_optimizeNodeLabeling(Finley_Mesh* mesh_p);
dim_t Finley_Mesh_FindMinDegreeNode(Paso_SystemMatrixPattern* pattern_p,index_t* available,index_t indicator);
index_t Finley_Mesh_getDegree(Paso_SystemMatrixPattern* pattern_p, index_t *label);


void Finley_Mesh_saveVTK(const char * filename_p, Finley_Mesh *mesh_p, const dim_t num_data,char* *names_p,escriptDataC* *data_pp);

#endif /* #ifndef INC_FINLEY_MESH */

