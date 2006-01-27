/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2003,2004,2005 -  All Rights Reserved              *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

/**************************************************************/

/*   Finley: integration schemes for element shapes Tri, Quad, Hex, Tet, Line, Point */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#ifndef INC_FINLEY_QUADRATURE
#define INC_FINLEY_QUADRATURE

/**************************************************************/

#include "Finley.h"

/**************************************************************/

#define MAX_numQuadNodesLine 10

/**************************************************************/

/*     Interfaces: */

typedef void (Finley_Quad_getNodes) (dim_t,double*,double*);
Finley_Quad_getNodes Finley_Quad_getNodesTri;
Finley_Quad_getNodes Finley_Quad_getNodesTet;
Finley_Quad_getNodes Finley_Quad_getNodesRec;
Finley_Quad_getNodes Finley_Quad_getNodesHex;
Finley_Quad_getNodes Finley_Quad_getNodesLine;
Finley_Quad_getNodes Finley_Quad_getNodesPoint;
Finley_Quad_getNodes Finley_Quad_getNodesTriOnFace;
Finley_Quad_getNodes Finley_Quad_getNodesRecOnFace;
Finley_Quad_getNodes Finley_Quad_getNodesLineOnFace;
Finley_Quad_getNodes Finley_Quad_getNodesPointOnFace;

typedef dim_t (Finley_Quad_getNumNodes) (dim_t);
Finley_Quad_getNumNodes Finley_Quad_getNumNodesPoint;
Finley_Quad_getNumNodes Finley_Quad_getNumNodesLine;
Finley_Quad_getNumNodes Finley_Quad_getNumNodesTri;
Finley_Quad_getNumNodes Finley_Quad_getNumNodesRec;
Finley_Quad_getNumNodes Finley_Quad_getNumNodesTet;
Finley_Quad_getNumNodes Finley_Quad_getNumNodesHex;

void Finley_Quad_makeNodesOnFace(dim_t, dim_t,double*,double*, Finley_Quad_getNodes);

#endif /* #ifndef INC_FINLEY_QUADRATURE */

/*
 * $Log$
 * Revision 1.3  2005/09/15 03:44:23  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.2.2.1  2005/09/07 06:26:20  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.2  2005/07/08 04:07:56  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.1  2005/06/29 02:34:55  gross
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
