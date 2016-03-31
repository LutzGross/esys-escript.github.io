
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
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
