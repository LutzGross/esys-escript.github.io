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
                                                                                                                                               
void Finley_ElementFile_prepare(Finley_ElementFile** in,dim_t numNodes, index_t *degreeOfFreedom) {
                                                                                                                                               
       /* rearrange elements: */
       Finley_ElementFile_optimizeDistribution(in);
                                                                                                                                               
       /* improve coloring */
       Finley_ElementFile_improveColoring(*in,numNodes,degreeOfFreedom);
}
/* 
* $Log$
* Revision 1.2  2005/07/08 04:07:50  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.1.1.1.2.1  2005/06/29 02:34:50  gross
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
