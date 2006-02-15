/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2003,2004,2005,2006 -  All Rights Reserved              *
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
