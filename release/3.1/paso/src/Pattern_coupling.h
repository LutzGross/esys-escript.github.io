
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

/*   Paso: Pattern_coupling                                   */

/**************************************************************/

/*   Author: Artak Amirbekyan */

/**************************************************************/

#ifndef INC_PASO_PATTERN_COUPLING
#define INC_PASO_PATTERN_COUPLING

#include "SparseMatrix.h"

/**************************************************************/

void Paso_Pattern_YS(Paso_SparseMatrix* A, index_t* mis_marker, double thershold);
void Paso_Pattern_RS(Paso_SparseMatrix* A, index_t* mis_marker, double theta);
void Paso_Pattern_Aggregiation(Paso_SparseMatrix* A, index_t* mis_marker, double theta);
void Paso_Pattern_greedy(Paso_Pattern* pattern, index_t* mis_marker);
void Paso_Pattern_greedy_color(Paso_Pattern* pattern, index_t* mis_marker);
void Paso_Pattern_greedy_diag(Paso_SparseMatrix* A, index_t* mis_marker, double thershold);

void Paso_Pattern_YS_plus(Paso_SparseMatrix* A, index_t* mis_marker, double alpha, double taw, double delta);
void Paso_Pattern_Standard(Paso_SparseMatrix* A, index_t* mis_marker, double theta);
void Paso_Pattern_greedy_RS(Paso_SparseMatrix* A, index_t* mis_marker, double theta);
void Paso_Pattern_greedy_Agg(Paso_SparseMatrix* A, index_t* mis_marker, double theta);
/*dim_t how_many(dim_t n,dim_t* S_i, int value1, dim_t* addedSet, int value2);*/

dim_t how_many(dim_t i,Paso_Pattern * S, bool_t transpose);
dim_t arg_max(dim_t n, dim_t* lambda, dim_t mask);
Paso_Pattern* Paso_Pattern_getTranspose(Paso_Pattern* P);
void Paso_Pattern_getReport(dim_t n,index_t* mis_marker);
void Paso_Pattern_Read(char *fileName,dim_t n,index_t* mis_marker);

#endif 
