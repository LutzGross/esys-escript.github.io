/* $Id$ */
/**************************************************************/
                                                                                                                                               
/*   Finley: Mesh: prepares the mesh for further calculations  */
                                                                                                                                               
/**************************************************************/
                                                                                                                                               
/*   Copyrights by ACcESS Australia 2003 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */
                                                                                                                                               
/**************************************************************/
                                                                                                                                               
                                                                                                                                               
#include "Common.h"
#include "ElementFile.h"

/**************************************************************/
                                                                                                                                               
void Finley_ElementFile_prepare(Finley_ElementFile** in,maybelong numNodes, maybelong *degreeOfFreedom) {
                                                                                                                                               
       /* rearrange elements: */
       Finley_ElementFile_optimizeDistribution(in);
                                                                                                                                               
       /* improve coloring */
       Finley_ElementFile_improveColoring(*in,numNodes,degreeOfFreedom);
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
