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

/**************************************************************/
                                                                                                                                               
/*   Finley: Mesh: prepares the mesh for further calculations  */
                                                                                                                                               
/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Vesrion: $Id$ */

/**************************************************************/
                                                                                                                                               
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
* Revision 1.3  2005/09/15 03:44:22  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.2.2.1  2005/09/07 06:26:18  gross
* the solver from finley are put into the standalone package paso now
*
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
