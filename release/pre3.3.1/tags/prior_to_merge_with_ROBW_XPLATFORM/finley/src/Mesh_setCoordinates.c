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

/*   Finley: Mesh: sets new coordinates for nodes */


/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Mesh.h"
#include "Util.h"

/**************************************************************/


void Finley_Mesh_setCoordinates(Finley_Mesh* self,escriptDataC* newX) {
  Finley_NodeFile_setCoordinates(self->Nodes,newX);
  Finley_ElementFile_setCoordinates(self->Elements,newX);
  Finley_ElementFile_setCoordinates(self->FaceElements,newX);
  Finley_ElementFile_setCoordinates(self->ContactElements,newX);
  Finley_ElementFile_setCoordinates(self->Points,newX);
}
/*
* $Log$
*/
