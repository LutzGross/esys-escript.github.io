/**************************************************************/

/*   Finley: Mesh: prepares the mesh for further calculations  */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003 */
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
* Revision 1.1  2004/10/26 06:53:57  jgs
* Initial revision
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/

