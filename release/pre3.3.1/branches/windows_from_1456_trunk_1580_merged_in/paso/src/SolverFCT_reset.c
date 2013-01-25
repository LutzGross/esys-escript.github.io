
/* $Id$ */

/*******************************************************
 *
 *       Copyright 2008 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/* Paso: Transort problem: reset operators */

/**************************************************************/

/* Author: l.gross@uq.edu.au */

/**************************************************************/

#include "SolverFCT.h"

/**************************************************************/

void  Paso_FCTransportProblem_reset(Paso_FCTransportProblem* in) 
{
    Paso_SystemMatrix_setValues(in->transport_matrix, 0.);
    Paso_SystemMatrix_setValues(in->mass_matrix, 0.);
    Paso_solve_free(in->iteration_matrix);
    in->valid_matrices=FALSE;
}
