
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/****************************************************************************/

/* Paso: SystemMatrix :                                       */
/*  sets the values of the system matrix to a value           */

/****************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

/****************************************************************************/

void  Paso_SystemMatrix_setValues(Paso_SystemMatrix* in,double value)
{
    if (in!=NULL) {
        in->mainBlock->setValues(value);
        in->col_coupleBlock->setValues(value);
        in->row_coupleBlock->setValues(value);
        in->is_balanced = false;
    }
}
