/* $Id$ */

#ifndef INC_FINLEY_NODEFILE
#define INC_FINLEY_NODEFILE

#define MAX_numDim 3

#include "Common.h"
#include "escript/Data/DataC.h"

struct Finley_NodeFile {
  maybelong numNodes;                      /* number of nodes */
  maybelong numDim;                        /* spatial dimension */
  maybelong *Id;                           /* Id[i] is the id number of node i. this number is not really
					      used but useful when the nodes are
					      relabled. see also Finley_resolveNodeIds. */
                                           /* in the entire code the term 'node id' refers to i but nor to Id[i] */
                                           /* if not explicitly stated otherwise. */
  maybelong *Tag;                          /* Tag[i] is the tag of node i. */
  double *Coordinates;                     /* Coordinates[INDEX2(k,i,numDim)] is the k-th coordinate of the */
                                           /* node i. */
  maybelong* degreeOfFreedom;              /* degreeOfFreedom[i] is the degree of freedom assigned to node i */
                                           /* this index is used to consider periodic boundary conditions by assigning */
                                           /* the same degreeOfFreedom to the same node */
  /* the following data are set by Finley_NodeFile_? */
  maybelong numDegreesOfFreedom;           /* number of degrees of freedom in the mesh (<=numNodes)*/
                                           /* notice that numDegreesOfFreedom=max(degreeOfFreedom[:]) */
  maybelong* reducedDegreeOfFreedom;       /* reducedDegreeOfFreedom[i] is the degree of freedom in the reduced version of the mesh */
                                           /* assigned to node i. reducedDegreeOfFreedom[i]<0 indicates that the node is not appearing */
                                           /* as a degree of freedom in the reduced version of the mesh. */
  maybelong reducedNumDegreesOfFreedom;    /* number of degrees of freedom in the reduced version of the mesh */
                                           /* notice that reducedNumDegreesOfFreedom=max(reducedDegreeOfFreedom[:]) */
  /* this information is used by interfaces to postprocessing tools that can deal with linear elements only: */
  maybelong reducedNumNodes;               /* returnes the number of nodes in the mesh when looking at the corners */
                                           /* of elements only                                                     */
  maybelong *toReduced;                    /* toReduced[i] is the node id in the reduced mesh. if toReduced[i]<0 it means that the node does not appear in the reduced mesh */
};

typedef struct Finley_NodeFile Finley_NodeFile;

Finley_NodeFile * Finley_NodeFile_alloc(int);
void Finley_NodeFile_dealloc(Finley_NodeFile*);
void Finley_NodeFile_setIdRange(int*,int*,Finley_NodeFile*);
void Finley_NodeFile_scatter(int*,Finley_NodeFile*,Finley_NodeFile*);
void Finley_NodeFile_gather(int*,Finley_NodeFile*,Finley_NodeFile*);
void Finley_NodeFile_setCoordinates(Finley_NodeFile*,escriptDataC*);
void Finley_NodeFile_copyTable(int,Finley_NodeFile*,int,int,Finley_NodeFile*);
void Finley_NodeFile_allocTable(Finley_NodeFile*,int);
void Finley_NodeFile_deallocTable(Finley_NodeFile*);

#endif

/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:57  jgs
 * Initial revision
 *
 * Revision 1.2  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 * Revision 1.1.1.1  2004/06/24 04:00:40  johng
 * Initial version of eys using boost-python.
 *
 *
 */
