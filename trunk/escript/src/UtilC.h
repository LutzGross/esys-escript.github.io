/* $Id$ */
/* 
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/
                                                                           
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
#define ESCRIPT_DEFAULT 0
#define ESCRIPT_DIRECT 1
#define ESCRIPT_CHOLEVSKY 2
#define ESCRIPT_PCG 3
#define ESCRIPT_CR 4
#define ESCRIPT_CGS 5
#define ESCRIPT_BICGSTAB 6
#define ESCRIPT_SSOR 7
#define ESCRIPT_ILU0 8
#define ESCRIPT_ILUT 9
#define ESCRIPT_JACOBI 10
#define ESCRIPT_GMRES 11
#define ESCRIPT_PRES20 12
#define ESCRIPT_LUMPING 13
#define ESCRIPT_NO_REORDERING 17
#define ESCRIPT_MINIMUM_FILL_IN 18
#define ESCRIPT_NESTED_DISSECTION 19
#define ESCRIPT_SCSL 14
#define ESCRIPT_MKL 15
#define ESCRIPT_UMFPACK 16
#define ESCRIPT_ITERATIVE 20
#define ESCRIPT_PASO 21
#define ESCRIPT_AMG 22
#define ESCRIPT_RILU 23

/* keyword in solver options relevant for lineraPDE class */
#define ESCRIPT_TOLERANCE_KEY "tolerance"
#define ESCRIPT_METHOD_KEY "method"
#define ESCRIPT_SYMMETRY_KEY "symmetric"
#define ESCRIPT_PACKAGE_KEY "package"

#endif
