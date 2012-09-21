
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


#ifndef INC_PASO_TESTS
#define INC_PASO_TESTS

/************************************************************************************/

/*   Some testing routines: */

/************************************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004,2005 */
/*   author: artak@uq.edu.au */

/************************************************************************************/

#include "paso/Common.h"

/************************************************************************************/

void Paso_test_run(Paso_SystemMatrix* A,double* b,dim_t level) ;
void Paso_test_matrix(Paso_SystemMatrix* A, double* b, Paso_Options* options );
void Paso_test_data(char *fileName_p, double* b, Paso_Options* options );
#endif /* #ifndef INC_PASO_TESTS */
