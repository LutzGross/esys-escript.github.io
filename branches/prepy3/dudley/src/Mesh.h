
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

#ifndef INC_DUDLEY_MESH
#define INC_DUDLEY_MESH

/**************************************************************/

/* Dudley: Mesh */

/* A mesh is built from nodes and elements which are describing the
   domain, the surface and point sources. (the latter are needed to
   establish links with other codes, in particular to particle
   codes). The nodes are stored a Dudley_NodeFile and elements in a
   Dudley_ElementFile. A Dudley_NodeFile and three Dudley_ElementFile
   containing the elements describing the domain, surface and point
   sources respectively. Notice that the surface elements do not
   necessaryly cover the entire surface of the domain. */

/* The element type is fixed by the reference element, see
   ReferenceElement.h. The numbering of the nodes starts with 0. */

/* Important: it is assumed that every node is appearing in at least
   one element or surface element and that any node used in an
   element, surface element or as a point is specified in the
   Dudley_Node, see also Dudley_resolveNodeIds. */

/* In some cases it is useful to refer to a mesh entirly built from
   order 1 (=linear) elements. The linear version of the mesh can be
   accessed by referning to the first few nodes of each element
   (thanks to the way the nodes are ordered). As the numbering of
   these nodes is not continuous a relabeling vectors are introduced
   in the Dudley_NodeFile. This feature is not fully implemented
   yet. */

/* allnodes and elements are tagged. the tag allows to group nodes and
   elements. A typical application is to mark surface elements on a
   certain portion of the domain with the same tag. All these surface
   elements can then assigned the same value eg. for the pressure. */

/* Thespacial dimension is determined by the type of elements
   used. The spacial dimension should be accessed by the function
   Dudley_Mesh_getDim. Notice that the element type also determines
   the type of surface elements to be used. */

/**************************************************************/

#include "Dudley.h"
#include "NodeFile.h"
#include "ElementFile.h"
#include "TagMap.h"
#include "Util.h"
#include "paso/SystemMatrixPattern.h"
#include "escript/DataC.h"

#ifdef ESYS_MPI
#include "esysUtils/Esys_MPI.h"
#endif

/**************************************************************/

/*  this struct holds a mesh: */

struct Dudley_Mesh {
    char *Name;			/* the name of the mesh */
    dim_t reference_counter;	/* counts the number of references to the mesh; */
    dim_t approximationOrder;
    dim_t reducedApproximationOrder;
    dim_t integrationOrder;
    dim_t reducedIntegrationOrder;
    Dudley_NodeFile *Nodes;	/* the table of the nodes */
    Dudley_ElementFile *Elements;	/* the table of the elements */
    Dudley_ElementFile *FaceElements;	/* the table of the face elements */
    Dudley_ElementFile *Points;	/* the table of points (treated as elements of dimension 0) */
    Dudley_TagMap *TagMap;	/* the tag map mapping names to tag keys */

    /* pointer to the sparse matrix pattern */

    Paso_SystemMatrixPattern *FullFullPattern;
    Paso_SystemMatrixPattern *FullReducedPattern;
    Paso_SystemMatrixPattern *ReducedFullPattern;
    Paso_SystemMatrixPattern *ReducedReducedPattern;
    Esys_MPIInfo *MPIInfo;
};

typedef struct Dudley_Mesh Dudley_Mesh;

/* these structures are used for matching surfaces elements: */

struct Dudley_Mesh_findMatchingFaces_center {
    index_t refId;
    double x[MAX_numDim];
};
typedef struct Dudley_Mesh_findMatchingFaces_center Dudley_Mesh_findMatchingFaces_center;

/**************************************************************/

/*  interfaces: */
Dudley_Mesh *Dudley_Mesh_alloc(char *name, dim_t numDim, Esys_MPIInfo * mpi_info);
Dudley_Mesh *Dudley_Mesh_reference(Dudley_Mesh *);
dim_t Dudley_Mesh_getDim(Dudley_Mesh *);
void Dudley_Mesh_free(Dudley_Mesh *);

void Dudley_Mesh_addTagMap(Dudley_Mesh * mesh_p, const char *name, index_t tag_key);
index_t Dudley_Mesh_getTag(Dudley_Mesh * mesh_p, const char *name);
bool_t Dudley_Mesh_isValidTagName(Dudley_Mesh * mesh_p, const char *name);
void Dudley_Mesh_distributeByRankOfDOF(Dudley_Mesh * in, dim_t * distribution);
Paso_SystemMatrixPattern *Dudley_getPattern(Dudley_Mesh * mesh, bool_t reduce_row_order, bool_t reduce_col_order);
Paso_SystemMatrixPattern *Dudley_makePattern(Dudley_Mesh * mesh, bool_t reduce_row_order, bool_t reduce_col_order);
void Dudley_Mesh_write(Dudley_Mesh *, char *);
void Dudley_Mesh_dump(Dudley_Mesh * in, char *fname);
void Dudley_PrintMesh_Info(Dudley_Mesh *, bool_t);
Dudley_Mesh *Dudley_Mesh_load(char *fname);
Dudley_Mesh *Dudley_Mesh_read(char *, index_t, index_t, bool_t);
Dudley_Mesh *Dudley_Mesh_readGmsh(char *, index_t, index_t, index_t, bool_t, bool_t);
void Dudley_Mesh_setOrders(Dudley_Mesh * in);

void Dudley_Mesh_setCoordinates(Dudley_Mesh *, escriptDataC *);
void Dudley_Mesh_setElements(Dudley_Mesh * self, Dudley_ElementFile * elements);
void Dudley_Mesh_setFaceElements(Dudley_Mesh * self, Dudley_ElementFile * elements);
void Dudley_Mesh_setPoints(Dudley_Mesh * self, Dudley_ElementFile * elements);

void Dudley_Mesh_optimizeDOFDistribution(Dudley_Mesh * in, dim_t * distribution);
void Dudley_Mesh_prepare(Dudley_Mesh * in, bool_t optimize);
void Dudley_Mesh_createColoring(Dudley_Mesh * in, index_t * node_localDOF_map);
void Dudley_Mesh_optimizeElementOrdering(Dudley_Mesh * in);
void Dudley_Mesh_resolveNodeIds(Dudley_Mesh *);
void Dudley_Mesh_createMappings(Dudley_Mesh * in, index_t * dof_distribution, index_t * node_distribution);
void Dudley_Mesh_createNodeFileMappings(Dudley_Mesh * in, dim_t numReducedNodes, index_t * indexReducedNodes,
					index_t * dof_first_component, index_t * nodes_first_component);
void Dudley_Mesh_markDOFsConnectedToRange(index_t * mask, index_t offset, index_t marker, index_t firstDOF,
					  index_t lastDOF, Dudley_Mesh * in, bool_t useLinear);

void Dudley_Mesh_optimizeDOFLabeling(Dudley_Mesh *, dim_t *);

Dudley_Mesh *Dudley_Mesh_merge(dim_t, Dudley_Mesh **);

void Dudley_Mesh_relableElementNodes(int *, int, Dudley_Mesh *);
void Dudley_Mesh_markNodes(int *, int, Dudley_Mesh *, int);

void Dudley_Mesh_glueFaces(Dudley_Mesh * self, double safety_factor, double tolerance, bool_t);
void Dudley_Mesh_joinFaces(Dudley_Mesh * self, double safety_factor, double tolerance, bool_t);

int Dudley_Mesh_findMatchingFaces_compar(const void *, const void *);
void Dudley_Mesh_findMatchingFaces(Dudley_NodeFile *, Dudley_ElementFile *, double, double, int *, int *, int *, int *);
void Dudley_Mesh_print(Dudley_Mesh * in);
void Dudley_Mesh_saveDX(const char *filename_p, Dudley_Mesh * mesh_p, const dim_t num_data, char **names_p,
			escriptDataC * *data_pp);
void Dudley_Mesh_optimizeNodeLabeling(Dudley_Mesh * mesh_p);
dim_t Dudley_Mesh_FindMinDegreeNode(Paso_SystemMatrixPattern * pattern_p, index_t * available, index_t indicator);
index_t Dudley_Mesh_getDegree(Paso_SystemMatrixPattern * pattern_p, index_t * label);

void Dudley_Mesh_saveVTK(const char *filename_p, Dudley_Mesh * mesh_p, const dim_t num_data, char **names_p,
			 escriptDataC * *data_pp, const char *metadata, const char *metadata_schema);
void Dudley_Mesh_setTagsInUse(Dudley_Mesh * in);

int Dudley_Mesh_getStatus(Dudley_Mesh * in);

#endif				/* #ifndef INC_DUDLEY_MESH */