
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

/**************************************************************/

/*   Dudley: Mesh: sets new coordinates for nodes */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Dudley_Mesh_setCoordinates(Dudley_Mesh * self, escriptDataC * newX)
{
    Dudley_NodeFile_setCoordinates(self->Nodes, newX);
}

/*
* $Log$
*/
