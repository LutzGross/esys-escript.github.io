/* $Id: $ */

/*******************************************************
 *
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#ifndef INC_SOLVERFCT
#define INC_SOLVERFCT

#include "SystemMatrix.h"

struct Paso_Solver_FluxControl {
    Paso_SystemMatrix * matrix;
    dim_t num_colors;
    index_t *colorOf;
    index_t *main_iptr;
};
typedef struct Paso_Solver_FluxControl Paso_Solver_FluxControl;


void Paso_Solver_FluxControl_free(Paso_Solver_FluxControl* in);
Paso_Solver_FluxControl* Paso_SolverFCT_getFluxControl(Paso_SystemMatrix * A);
void Paso_Solver_FluxControl_setAntiDiffusiveFlux(Paso_Solver_FluxControl * fc, double * u, double* fa);
void Paso_Solver_FluxControl_addDiffusion(Paso_Solver_FluxControl * fc, double alpha, Paso_SystemMatrix * B);



#endif /* #ifndef INC_SOLVERFCT */
