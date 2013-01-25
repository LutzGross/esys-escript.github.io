
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#if !defined  escript_UtilC_20040611_H
#define escript_UtilC_20040611_H
#include "system_dep.h"

#define ESCRIPT_MAX_DATA_RANK 4
/**
   \brief
   solver options, they have to be consistent with LinearPDE.py settings 

   Description:
   solver options, they have to be consistent with LinearPDE.py settings 

*/
#define ESCRIPT_DEFAULT  0
#define ESCRIPT_DIRECT  1
#define ESCRIPT_CHOLEVSKY  2
#define ESCRIPT_PCG  3
#define ESCRIPT_CR  4
#define ESCRIPT_CGS  5
#define ESCRIPT_BICGSTAB  6
#define ESCRIPT_SSOR  7
#define ESCRIPT_ILU0  8
#define ESCRIPT_ILUT  9
#define ESCRIPT_JACOBI  10
#define ESCRIPT_GMRES  11
#define ESCRIPT_PRES20  12
#define ESCRIPT_LUMPING  13
#define ESCRIPT_NO_REORDERING  17
#define ESCRIPT_MINIMUM_FILL_IN  18
#define ESCRIPT_NESTED_DISSECTION  19
#define ESCRIPT_MKL  15
#define ESCRIPT_UMFPACK  16
#define ESCRIPT_ITERATIVE  20
#define ESCRIPT_PASO  21
#define ESCRIPT_AMG  22
#define ESCRIPT_REC_ILU   23
#define ESCRIPT_TRILINOS   24
#define ESCRIPT_NONLINEAR_GMRES   25
#define ESCRIPT_TFQMR   26
#define ESCRIPT_MINRES   27
#define ESCRIPT_GAUSS_SEIDEL 28
#define ESCRIPT_RILU 29
#define ESCRIPT_DEFAULT_REORDERING 30
#define ESCRIPT_SUPER_LU 31
#define ESCRIPT_PASTIX 32
#define ESCRIPT_YAIR_SHAPIRA_COARSENING 33
#define ESCRIPT_RUGE_STUEBEN_COARSENING 34
#define ESCRIPT_AGGREGATION_COARSENING 35
#define ESCRIPT_NO_PRECONDITIONER 36
#define ESCRIPT_MIN_COARSE_MATRIX_SIZE 37
#define ESCRIPT_AMLI  38
#define ESCRIPT_STANDARD_COARSENING 39
#define ESCRIPT_CLASSIC_INTERPOLATION_WITH_FF_COUPLING 50
#define ESCRIPT_CLASSIC_INTERPOLATION 51
#define ESCRIPT_DIRECT_INTERPOLATION 52
#define ESCRIPT_BOOMERAMG 60
#define ESCRIPT_CIJP_FIXED_RANDOM_COARSENING 61
#define ESCRIPT_CIJP_COARSENING 62
#define ESCRIPT_FALGOUT_COARSENING 63
#define ESCRIPT_PMIS_COARSENING 64
#define ESCRIPT_HMIS_COARSENING 65
#define ESCRIPT_LINEAR_CRANK_NICOLSON 66
#define ESCRIPT_CRANK_NICOLSON 67
#define ESCRIPT_BACKWARD_EULER 68


#endif
