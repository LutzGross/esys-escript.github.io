/* $Id$ */

#ifndef INC_FINLEY_ELEMENTFILE
#define INC_FINLEY_ELEMENTFILE

#include "Common.h"
#include "ReferenceElements.h"

struct Finley_ElementFile {

  Finley_RefElement* ReferenceElement;           /* the reference element, see
						    Reference element.c */

  Finley_RefElement* LinearReferenceElement;     /* the reference element for
						    the linear mesh. it is
						    vital that both are using
						    the same quadrature
						    scheme */

  maybelong numElements;                         /* number of elements. */

  maybelong *Id;                                 /* Id[i] is the id nmber of
						    node i. this number is not
						    used but useful when
						    elements are resorted. in
						    the entire code the term
						    'element id' refers to i
						    but nor to Id[i] if not
						    explicitly stated
						    otherwise. */

  maybelong *Tag;                                /* Tag[i] is the tag of
						    element i. */

  maybelong *Nodes;                              /* Nodes[INDEX(k, i,
						    ReferenceElement.Type.numNodes)
						    is the k-the node in the
						    i-the element. note that
						    in the way the nodes are
						    ordered Nodes[INDEX(k, i,
						    LinearReferenceElement.Type.numNodes)
						    is k-the node of element i
						    when refering to the
						    linear version of the
						    mesh. */
  maybelong numColors;                           /* number of colours (sould be as small as possible) */
  maybelong *Color;                              /* assigns each element a color. elements with the same color     */
						 /* are don't share a node so they can be processed simultaneously */
                                                 /* at anytime Color must provide a valid value. In any case one can set  */
                                                 /* Color[e]=e  for all e */
  int order;					 /* order of the element */
};

typedef struct Finley_ElementFile Finley_ElementFile;

Finley_ElementFile* Finley_ElementFile_alloc(ElementTypeId,int);
void Finley_ElementFile_dealloc(Finley_ElementFile*);
void Finley_ElementFile_improveColoring(Finley_ElementFile* in,maybelong numNodes,maybelong* degreeOfFreedom);
void Finley_ElementFile_optimizeDistribution(Finley_ElementFile** in);
void Finley_ElementFile_setNodeRange(int*,int*,Finley_ElementFile*);
void Finley_ElementFile_relableNodes(int*,int,Finley_ElementFile*);
void Finley_ElementFile_markNodes(int*,int,Finley_ElementFile*,int);
void Finley_ElementFile_scatter(int*,Finley_ElementFile*,Finley_ElementFile*);
void Finley_ElementFile_gather(int*,Finley_ElementFile*,Finley_ElementFile*);
void Finley_ElementFile_copyTable(int,Finley_ElementFile*,int,int,Finley_ElementFile*);
void Finley_ElementFile_allocTable(Finley_ElementFile*,int);
void Finley_ElementFile_deallocTable(Finley_ElementFile*);
void Finley_ElementFile_prepare(Finley_ElementFile** in,maybelong numNodes,maybelong* degreeOfFreedom);

#endif /* #ifndef INC_FINLEY_ELEMENTFILE */

/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:57  jgs
 * Initial revision
 *
 * Revision 1.1.1.1  2004/06/24 04:00:40  johng
 * Initial version of eys using boost-python.
 *
 *
*/
