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
