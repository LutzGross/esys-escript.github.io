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

#ifndef INC_FINLEY_NODEFILE
#define INC_FINLEY_NODEFILE

#define MAX_numDim 3

#include "Finley.h"
#include "escript/DataC.h"

#ifdef PASO_MPI
#include "Distribution.h"
#include "./paso/CommBuffer.h"
#endif

struct Finley_NodeFile {
#ifdef PASO_MPI
  Paso_MPIInfo *MPIInfo;              /* MPI information */
  Finley_NodeDistribution *degreeOfFreedomDistribution;  /* information about the distribution of degrees of freedom
                                              on this subdomain and over other subdomains */
  Finley_NodeDistribution *reducedDegreeOfFreedomDistribution;
  Paso_CommBuffer *CommBuffer;
  Paso_CommBuffer *reducedCommBuffer;

	index_t *Dom;													/* flags whether a node references an internal/boundary/external DOF. Is one of
																					 either NODE_INTERNAL, NODE_BOUNDARY or NODE_EXTERNAL, as defined in
																					 Mesh.h	*/ 
#endif
  index_t isPrepared;                   /* UNKNOWN,  UNPREPARED, PREPAED  to indicate that the Nodes are ready for calculation */
  dim_t numNodes;                      /* number of nodes */
  dim_t numDim;                        /* spatial dimension */
  index_t *Id;                         /* Id[i] is the id number of node i. this number is not really
					      used but useful when the nodes are
					      relabled. see also Finley_resolveNodeIds. */
                                       /* in the entire code the term 'node id' refers to i but nor to Id[i] */
                                       /* if not explicitly stated otherwise. */
  index_t *Tag;                        /* Tag[i] is the tag of node i. */
  double *Coordinates;                 /* Coordinates[INDEX2(k,i,numDim)] is the k-th coordinate of the */
                                       /* node i. */
  index_t* degreeOfFreedom;            /* degreeOfFreedom[i] is the degree of freedom assigned to node i */
                                       /* this index is used to consider periodic boundary conditions by assigning */
                                       /* the same degreeOfFreedom to the same node */
  index_t* degreeOfFreedomId;          /* the ids of the degreeOfFreedom */
  /* the following data are set by Finley_NodeFile_preparNodes */
  dim_t numDegreesOfFreedom;           /* number of degrees of freedom in the mesh (<=numNodes)*/
                                       /* notice that numDegreesOfFreedom=max(degreeOfFreedom[:]) */
  index_t* reducedDegreeOfFreedom;     /* reducedDegreeOfFreedom[i] is the degree of freedom in the reduced version of the mesh */
                                       /* assigned to node i. reducedDegreeOfFreedom[i]<0 indicates that the node is not appearing */
                                       /* as a degree of freedom in the reduced version of the mesh. */
  index_t* reducedDegreeOfFreedomId;  /* the ids of the reducedDegreeOfFreedom */
  dim_t reducedNumDegreesOfFreedom;    /* number of degrees of freedom in the reduced version of the mesh */
                                       /* notice that reducedNumDegreesOfFreedom=max(reducedDegreeOfFreedom[:]) */
  /* this information is used by interfaces to postprocessing tools that can deal with linear elements only: */
  dim_t reducedNumNodes;               /* returnes the number of nodes in the mesh when looking at the corners */
                                       /* of elements only                                                     */
  index_t *toReduced;                  /* toReduced[i] is the node id in the reduced mesh. if toReduced[i]<0 it means that the node does not appear in the reduced mesh */

  int status; /* the status counts the updates done on the node coordinates */
              /* the value of status is increased by when the node coordinates are updated.*/
                                                                                                                                                                                                 
                                                                                                                                                                                                 
};



typedef struct Finley_NodeFile Finley_NodeFile;

#ifdef PASO_MPI
Finley_NodeFile* Finley_NodeFile_alloc(dim_t, Paso_MPIInfo *MPIInfo);
void Finley_NodeFile_allocTable(Finley_NodeFile*,dim_t);
#else
Finley_NodeFile* Finley_NodeFile_alloc(dim_t numDim);
void Finley_NodeFile_allocTable(Finley_NodeFile*,dim_t);
#endif

/* Finley_NodeFile * Finley_NodeFile_alloc(dim_t); */
void Finley_NodeFile_dealloc(Finley_NodeFile*);
void Finley_NodeFile_setIdRange(index_t*,index_t*,Finley_NodeFile*);
void Finley_NodeFile_scatter(index_t*,Finley_NodeFile*,Finley_NodeFile*);
void Finley_NodeFile_gather(index_t*,Finley_NodeFile*,Finley_NodeFile*);
void Finley_NodeFile_setCoordinates(Finley_NodeFile*,escriptDataC*);
void Finley_NodeFile_copyTable(dim_t,Finley_NodeFile*,dim_t,dim_t,Finley_NodeFile*);
void Finley_NodeFile_setTags(Finley_NodeFile*,const int,escriptDataC*);
void Finley_NodeFile_deallocTable(Finley_NodeFile*);

#endif

/*
 * $Log$
 * Revision 1.3  2005/09/15 03:44:23  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.2.2.1  2005/09/07 06:26:20  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.2  2005/07/08 04:07:55  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.1  2005/06/29 02:34:54  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.2  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 * Revision 1.1.1.1  2004/06/24 04:00:40  johng
 * Initial version of eys using boost-python.
 *
 *
 */
