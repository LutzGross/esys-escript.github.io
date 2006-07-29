/* created by Ben Cumming on 26/04/2006 */

#ifndef INC_DISTRIBUTION
#define INC_DISTRIBUTION


#include "Finley.h"

#ifdef PASO_MPI

#include "./paso/CommBuffer.h"

#define NODE_INTERNAL 1
#define NODE_BOUNDARY 2
#define NODE_EXTERNAL 3

#define ELEMENT_INTERNAL 1
#define ELEMENT_BOUNDARY 2

/* DATA TYPES */

/**************************************************** 
  describes the set of nodes shared between two
  neighbouring domains. Each of the domains has
  a local copy that show where to map information
  received from/sent to the other neighbour
****************************************************/
struct Finley_NodeGhostEdge
{
  dim_t reference_counter;
  index_t numForward;  /* number of local nodes referenced by neighbour */
  index_t numBackward; /* number of neighbours nodes referenced by this process */
  /* local indices of the forward and backward referenced nodes */
  index_t *indexForward;   
  index_t *indexBackward;
};

typedef struct Finley_NodeGhostEdge Finley_NodeGhostEdge;

/**************************************************** 
  describes the distribution of the nodes/DOF stored
  on the local process, along with their connections
  with nodes/DOF in neighbouring domains.
****************************************************/
struct Finley_NodeDistribution
{
  dim_t reference_counter;
  Paso_MPIInfo *MPIInfo;
  index_t numLocal;       /* total number of nodes on local domain
                             numLocal = numBoundary + numInternal */
  index_t numInternal;    /* number of local nodes that are internal, 
                             that is, not dependent on nodes in other domains*/
  index_t numBoundary;    /* number of local nodes that are dependant on nodes
                             in other domains */
  index_t numExternal;    /* number of nodes belonging to other subdomains that
                             share elements with local boundary nodes */
  index_t *indexExternal;  /* global indices of the external nodes stored on this Pid */
  index_t *vtxdist;        /* process i has nodes with global indices
                             vtxdist[i] to vtxdist[i]-1. */
  index_t numGlobal;      /* total number of nodes in the global domain */
  index_t numNeighbours;  /* number of neighbour domains */
  index_t *neighbours;    /* list of ranks of neighbours */
  Finley_NodeGhostEdge **edges; /* ghost edges shared with each neighbour */
};

typedef struct Finley_NodeDistribution Finley_NodeDistribution;

/* not used at the moment, but could be used in the future for more efficient
   calculation of integrals etc on boundary elements... */
struct Finley_ElementGhostEdge
{
  dim_t reference_counter;
  index_t numShared;  /* number of elements shared with neighbour */
  index_t *shared;    /* local indices of the elements shared with the neighbour */
  index_t *pointers;  /* Yale-style description of dependancies on neighbour elements */
  index_t *index;
};

typedef struct Finley_ElementGhostEdge Finley_ElementGhostEdge;

struct Finley_ElementDistribution
{
  dim_t reference_counter;
  Paso_MPIInfo *MPIInfo;
  index_t numLocal;       /* total number of elements on local domain
                             numLocal = numBoundary + numInternal */
  index_t numInternal;    /* number of local elements that are internal, 
                             that is, not dependent on elements in other domains*/
  index_t numBoundary;    /* number of local elements that are dependant on elements
                             in other domains */
  dim_t *vtxdist;         /* process i has elements with Id's between vtxdist[i] up to vtxdist[i+1]-1 inclusive */
  /* there will be further stuff here, as the need for it in domain decomposition arises */
  /* ... */
};

typedef struct Finley_ElementDistribution Finley_ElementDistribution;

/***************************************
 Function prototypes 
***************************************/

/* Finley_NodeDistribution */
Finley_NodeDistribution*  Finley_NodeDistribution_alloc( Paso_MPIInfo *MPIInfo );
void                      Finley_NodeDistribution_dealloc( Finley_NodeDistribution *in );
Finley_NodeDistribution*  Finley_NodeDistribution_getReference( Finley_NodeDistribution *in );
void                      Finley_NodeDistribution_allocTable( Finley_NodeDistribution *in, dim_t numLocal, dim_t numExternal, dim_t numNeighbours );
void                      Finley_NodeDistribution_deallocTable( Finley_NodeDistribution *in );
void                      Finley_NodeDistribution_addForward( Finley_NodeDistribution *in, index_t domain, dim_t numForward, index_t* indexLocal  );
void                      Finley_NodeDistribution_addBackward( Finley_NodeDistribution *in, index_t domain, dim_t numBackward, index_t* indexLocal  );
void                      Finley_NodeDistribution_calculateIndexExternal( Finley_NodeDistribution *Distribution, Paso_CommBuffer *CommBuffer );
void                      Finley_NodeDistribution_formCommBuffer( Finley_NodeDistribution *in, Paso_CommBuffer *CommBuffer );
void                      Finley_NodeDistribution_print( Finley_NodeDistribution *in, FILE *fid );

/* Finley_NodeGhostEdge */
Finley_NodeGhostEdge* Finley_NodeGhostEdge_alloc( void );
void                  Finley_NodeGhostEdge_dealloc( Finley_NodeGhostEdge *in );
Finley_NodeGhostEdge* Finley_NodeGhostEdge_getReference( Finley_NodeGhostEdge *in );
void                  Finley_NodeGhostEdge_allocTable( Finley_NodeGhostEdge *in, dim_t numForward, dim_t numBackward );
void                  Finley_NodeGhostEdge_deallocTable( Finley_NodeGhostEdge *in );

/* Finley_ElementDistribution */
Finley_ElementDistribution* Finley_ElementDistribution_alloc( Paso_MPIInfo *MPIInfo );
void                        Finley_ElementDistribution_dealloc( Finley_ElementDistribution* in );
Finley_ElementDistribution* Finley_ElementDistribution_getReference( Finley_ElementDistribution* in );
void 												Finley_ElementDistribution_allocTable( Finley_ElementDistribution *in, dim_t numElements, dim_t numElementsThis );

#endif
#endif
