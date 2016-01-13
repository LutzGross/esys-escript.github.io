
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


/**************************************************************/

/*    Paso finite element solver */

/**************************************************************/

/*  Copyrights by ACcESS Australia, 2003,2004,2005 */
/*  Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#ifndef INC_PASO
#define INC_PASO

#include "Common.h"
#include "esysUtils/error.h"

#define MATRIX_FORMAT_DEFAULT 1
#define MATRIX_FORMAT_CSC 2
#define MATRIX_FORMAT_BLK1 4
#define MATRIX_FORMAT_OFFSET1 8
#define MATRIX_FORMAT_TRILINOS_CRS 16
#define MATRIX_FORMAT_DIAGONAL_BLOCK 32

#define PASO_ONE (double)(1.0)
#define PASO_ZERO (double)(0.0)

#endif /* #ifndef INC_PASO */
