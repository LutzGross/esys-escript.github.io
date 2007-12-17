/* $Id:$ */

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

/**************************************************************/

/* Paso: Flux correction transport solver
 *
 * solves Mu_t=Du+Ku+q
 *
 *  where is D is diffusive (not checked)
 *  		- D is symmetric
 *  		- row sums are equal to zero.
 *  and  K is the advective part.
 *
 *        u(0) >= 0
 *
 * intially fctp->transport_matrix defines the diffusive part 
 * but the matrix is updated by the adevctive part + artificial diffusion 
 *
*/
/**************************************************************/

/* Author: l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "Solver.h"
#include "SolverFCT.h"
#include "escript/blocktimer.h"

/***********************************************************************************/


void Paso_SolverFCT_solve(Paso_FCTransportProblem* fctp, double* u, double dt, Paso_Options* options,Paso_Performance* pp) {

   if (dt<=0.) {
       Paso_setError(TYPE_ERROR,"Paso_SolverFCT_solve: dt must be positive.");
   }
   if (! fctp->valid_matrices) {

        /* extract the row sum of the advective part */
        Paso_SystemMatrix_rowSum(fctp->flux_matrix,fctp->row_sum_flux_matrix);

        /* add the advective part + artificial diffusion to the diffusive part */
        Paso_FCTransportProblem_addAdvectivePart(fctp,1.);

        if (Paso_noError()) fctp->valid_matrices=TRUE;
   }

}
