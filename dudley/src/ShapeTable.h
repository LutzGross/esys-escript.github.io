
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/* Shape Function info
These tables are a much simplified version of content from finley's ShapeFunctions files

This file is not to be included in .h files - only .c files should have any use for it
*/

#ifndef SHAPETABLE_DUDLEY
#define SHAPETABLE_DUDLEY

/* This optimisation assumes that loop unrolling is enabled for this file */
static const double DTDV_2D[3][2]={{-1,-1}, {1,0}, {0,1}};
static const double DTDV_3D[4][3]={{-1, -1, -1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};



#endif 



