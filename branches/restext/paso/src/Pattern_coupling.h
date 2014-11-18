
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   Paso: Pattern_coupling                                   */

/**************************************************************/

/*   Author: Artak Amirbekyan */

/**************************************************************/

#ifndef INC_PASO_PATTERN_COUPLING
#define INC_PASO_PATTERN_COUPLING

#include "SparseMatrix.h"

/**************************************************************/

void Paso_Pattern_coup(Paso_SparseMatrix* A, index_t* mis_marker, double thershold);
void Paso_Pattern_RS(Paso_SparseMatrix* A, index_t* mis_marker, double theta);
void Paso_Pattern_Aggregiation(Paso_SparseMatrix* A, index_t* mis_marker, double theta);
void Paso_Pattern_greedy(Paso_Pattern* pattern, index_t* mis_marker);
void Paso_Pattern_greedy_color(Paso_Pattern* pattern, index_t* mis_marker);


#endif 