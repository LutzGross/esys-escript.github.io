
/*******************************************************
*
* Copyright (c) 2010-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

/* Author: Joel Fenwick */

/* An MPI aware Maximal independent set algorithm for SystemMatrix */

#ifndef INC_SYSTEMMATRIX_MIS

#include "SystemMatrix.h"

/* Returns the number of local nodes in the MIS.
   The second param will store a pointer to a list of the nodes in the MIS.
   Deallocating the value sent back in the second param is caller's responsibility.
   
   Sets paso error on failure.
*/
index_t Paso_SystemMatrix_getMIS(Paso_SystemMatrix*, index_t**);

#endif
