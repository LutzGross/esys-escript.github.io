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

/*   Finley: Mesh: prepares the mesh for further calculations  */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Finley_Mesh_prepare(Finley_Mesh* in) {

       /* set the labeling vectors in node files: */
       Finley_Mesh_prepareNodes(in);

       /* rearrange elements: */
       Finley_Mesh_optimizeElementDistribution(in);

       /* improve coloring */
       Finley_Mesh_improveColoring(in);
}
/*                                                      */
/*  tries to reduce the coloring for all element files: */
/*                                                      */
void Finley_Mesh_improveColoring(Finley_Mesh* in) {
  Finley_ElementFile_improveColoring(in->Elements,in->Nodes->numNodes,in->Nodes->degreeOfFreedom);
  Finley_ElementFile_improveColoring(in->FaceElements,in->Nodes->numNodes,in->Nodes->degreeOfFreedom);
  Finley_ElementFile_improveColoring(in->Points,in->Nodes->numNodes,in->Nodes->degreeOfFreedom);
  Finley_ElementFile_improveColoring(in->ContactElements,in->Nodes->numNodes,in->Nodes->degreeOfFreedom);
}
/*                                                                    */
/*  redistribute elements to minimize communication during assemblage */
/*                                                                    */
void Finley_Mesh_optimizeElementDistribution(Finley_Mesh* in) {
  Finley_ElementFile_optimizeDistribution(&(in->Elements));
  Finley_ElementFile_optimizeDistribution(&(in->FaceElements));
  Finley_ElementFile_optimizeDistribution(&(in->Points));
  Finley_ElementFile_optimizeDistribution(&(in->ContactElements));
}

/*
* $Log$
* Revision 1.2  2005/09/15 03:44:22  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.1.1.1.6.1  2005/09/07 06:26:19  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/

