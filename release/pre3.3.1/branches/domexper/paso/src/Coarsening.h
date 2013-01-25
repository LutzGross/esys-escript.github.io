
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

#ifndef INC_COARSE
#define INC_COARSE

#include "SystemMatrix.h"

#define PASO_COARSENING_IN_F 1
#define PASO_COARSENING_IN_C 2

void Paso_Coarsening_Local(index_t* mis_marker, Paso_SparseMatrix* A, double threshold, const index_t coarsening_method);

void Paso_Coarsening_Local_YS(Paso_SparseMatrix* A, index_t* mis_marker, double thershold);
void Paso_Coarsening_Local_RS(Paso_SparseMatrix* A, index_t* mis_marker, double theta);
void Paso_Coarsening_Local_Aggregiation(Paso_SparseMatrix* A, index_t* mis_marker, double theta);


/*===============*/
void Paso_Coarsening_Local_greedy(Paso_Pattern* pattern, index_t* mis_marker);
void Paso_Coarsening_Local_greedy_color(Paso_Pattern* pattern, index_t* mis_marker);
void Paso_Coarsening_Local_greedy_diag(Paso_SparseMatrix* A, index_t* mis_marker, double thershold);

void Paso_Coarsening_Local_YS_plus(Paso_SparseMatrix* A, index_t* mis_marker, double alpha, double taw, double delta);
void Paso_Coarsening_Local_Standard(Paso_SparseMatrix* A, index_t* mis_marker, double theta);
void Paso_Coarsening_Local_greedy_RS(Paso_SparseMatrix* A, index_t* mis_marker, double theta);
void Paso_Coarsening_Local_greedy_Agg(Paso_SparseMatrix* A, index_t* mis_marker, double theta);
/*dim_t how_many(dim_t n,dim_t* S_i, int value1, dim_t* addedSet, int value2);*/
void Paso_Coarsening_Local_Standard_Block(Paso_SparseMatrix* A, index_t* mis_marker, double theta);

dim_t how_many(dim_t i,Paso_Pattern * S, bool_t transpose);
dim_t arg_max(dim_t n, dim_t* lambda, dim_t mask);
Paso_Pattern* Paso_Coarsening_Local_getTranspose(Paso_Pattern* P);
void Paso_Coarsening_Local_getReport(dim_t n,index_t* mis_marker);
void Paso_Coarsening_Local_Read(char *fileName,dim_t n,index_t* mis_marker);
void Paso_Coarsening_Local_Write(char *fileName,dim_t n,index_t* mis_marker);


#endif
