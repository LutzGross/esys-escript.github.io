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

#ifndef INC_FINLEY_ELEMENTFILE
#define INC_FINLEY_ELEMENTFILE

#include "Finley.h"
#include "NodeFile.h"
#include "ReferenceElements.h"
#include "escript/DataC.h"

#ifdef PASO_MPI
#include "paso/Paso_MPI.h"
#include "Distribution.h"
#endif

struct Finley_ElementFile_Jacobeans {
  Finley_Status_t status;               /* status of mesh when jacobeans where updated last time */
  dim_t numDim;                         /* spatial dimension */
  Finley_RefElement* ReferenceElement;  /* reference elemnt used to calculate jacobeans (this is a borrowd reference) */
  double* volume;                       /* local volume */
  double* DSDX;                         /* derivatives of shape functions in global coordinates at quadrature points*/
};

typedef struct Finley_ElementFile_Jacobeans Finley_ElementFile_Jacobeans;

struct Finley_ElementFile {
#ifdef PASO_MPI
  Paso_MPIInfo *MPIInfo;
  Finley_ElementDistribution *elementDistribution;
	index_t *Dom;
#endif

  Finley_RefElement* ReferenceElement;           /* the reference element, see
						    Reference element.c */

  Finley_RefElement* LinearReferenceElement;     /* the reference element for
						    the linear mesh. it is
						    vital that both are using
						    the same quadrature
						    scheme */

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

  index_t *Tag;                                /* Tag[i] is the tag of
						    element i. */

  index_t *Nodes;                              /* Nodes[INDEX(k, i,
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
  index_t minColor;                           /* minimum color */
  index_t maxColor;                           /* maximum color */
  index_t *Color;                             /* assigns each element a color. elements with the same color     */
				              /* are don't share a node so they can be processed simultaneously */
                                              /* at anytime Color must provide a valid value. In any case one can set  */
                                              /* Color[e]=e  for all e */
  index_t order;			       /* order of the element */

  Finley_ElementFile_Jacobeans* jacobeans;           /* element jacobeans */
  Finley_ElementFile_Jacobeans* jacobeans_reducedS;  /* element jacobeans for reduced order of shape function*/
  Finley_ElementFile_Jacobeans* jacobeans_reducedQ;  /* element jacobeans for reduced integration order*/
  Finley_ElementFile_Jacobeans* jacobeans_reducedS_reducedQ;  /* element jacobeans for reduced integration order and  reduced order of shape function*/

};

typedef struct Finley_ElementFile Finley_ElementFile;

#ifndef PASO_MPI
Finley_ElementFile* Finley_ElementFile_alloc(ElementTypeId,dim_t);
#else
Finley_ElementFile* Finley_ElementFile_alloc( ElementTypeId, dim_t, Paso_MPIInfo* );
void Finley_ElementFile_setDomainFlags( Finley_ElementFile *in  );
#endif

void Finley_ElementFile_dealloc(Finley_ElementFile*);
void Finley_ElementFile_improveColoring(Finley_ElementFile* in,dim_t numNodes,dim_t* degreeOfFreedom);
void Finley_ElementFile_optimizeDistribution(Finley_ElementFile** in);
void Finley_ElementFile_setNodeRange(dim_t*,dim_t*,Finley_ElementFile*);
void Finley_ElementFile_relableNodes(dim_t*,dim_t,Finley_ElementFile*);
void Finley_ElementFile_markNodes(dim_t*,dim_t,Finley_ElementFile*,dim_t);
void Finley_ElementFile_scatter(dim_t*,Finley_ElementFile*,Finley_ElementFile*);
void Finley_ElementFile_gather(dim_t*,Finley_ElementFile*,Finley_ElementFile*);
void Finley_ElementFile_copyTable(dim_t,Finley_ElementFile*,dim_t,dim_t,Finley_ElementFile*);
void Finley_ElementFile_allocTable(Finley_ElementFile*,dim_t);
void Finley_ElementFile_deallocTable(Finley_ElementFile*);
void Finley_ElementFile_prepare(Finley_ElementFile** in,dim_t numNodes,dim_t* degreeOfFreedom);
void Finley_ElementFile_setTags(Finley_ElementFile*,const int,escriptDataC*);
Finley_ElementFile_Jacobeans* Finley_ElementFile_Jacobeans_alloc(Finley_RefElement*);
void Finley_ElementFile_Jacobeans_dealloc(Finley_ElementFile_Jacobeans*);
Finley_ElementFile_Jacobeans* Finley_ElementFile_borrowJacobeans(Finley_ElementFile*, Finley_NodeFile*, bool_t, bool_t);


#endif /* #ifndef INC_FINLEY_ELEMENTFILE */

/*
 * $Log$
 * Revision 1.3  2005/09/15 03:44:21  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.2.2.1  2005/09/07 06:26:18  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.2  2005/07/08 04:07:49  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.1  2005/06/29 02:34:49  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.1.1.1  2004/06/24 04:00:40  johng
 * Initial version of eys using boost-python.
 *
 *
*/
