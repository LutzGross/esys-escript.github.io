
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   Paso: SystemMatrix and SystemVector */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003,2004,2005 */
/*   Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_OPTIONS
#define INC_PASO_OPTIONS

#include "Paso_MPI.h"

/* solver options */

#define PASO_DEFAULT 0
#define PASO_DIRECT 1
#define PASO_CHOLEVSKY 2
#define PASO_PCG 3
#define PASO_CR 4
#define PASO_CGS 5
#define PASO_BICGSTAB 6
#define PASO_SSOR 7
#define PASO_ILU0 8
#define PASO_ILUT 9
#define PASO_JACOBI 10
#define PASO_GMRES 11
#define PASO_PRES20 12
#define PASO_LUMPING 13
#define PASO_NO_REORDERING 17
#define PASO_MINIMUM_FILL_IN 18
#define PASO_NESTED_DISSECTION 19
#define PASO_MKL 15
#define PASO_UMFPACK 16
#define PASO_ITERATIVE 20
#define PASO_PASO 21
#define PASO_AMG 22
#define PASO_REC_ILU  23
#define PASO_TRILINOS  24
#define PASO_NONLINEAR_GMRES  25
#define PASO_TFQMR  26
#define PASO_MINRES  27
#define PASO_GAUSS_SEIDEL 28
#define PASO_GS PASO_GAUSS_SEIDEL
#define PASO_RILU 29
#define PASO_DEFAULT_REORDERING 30
#define PASO_SUPER_LU 31
#define PASO_PASTIX 32
#define PASO_YAIR_SHAPIRA_COARSENING 33
#define PASO_RUGE_STUEBEN_COARSENING 34
#define PASO_AGGREGATION_COARSENING 35
#define PASO_NO_PRECONDITIONER 36

typedef struct {
    index_t method;
    index_t package;
    bool_t symmetric;
    double tolerance;
    double absolute_tolerance;
    double inner_tolerance;
    bool_t adapt_inner_tolerance;
    bool_t verbose;
    bool_t reordering;
    index_t preconditioner;
    dim_t iter_max;
    dim_t inner_iter_max;
    double drop_tolerance;
    double drop_storage;
    index_t truncation;
    index_t restart;
    dim_t sweeps;
    dim_t pre_sweeps;
    dim_t post_sweeps;
    dim_t level_max;
    double coarsening_threshold;
    bool_t accept_failed_convergence;
    index_t coarsening_method;
    double relaxation_factor;

    /* diagnostic values */
    dim_t num_iter;
    dim_t num_level;
    dim_t num_inner_iter;
    double time;
    double set_up_time;
    double residual_norm;
    bool_t converged;

} Paso_Options;

/*  interfaces: */


PASO_DLL_API
void Paso_Options_setDefaults(Paso_Options* in);

PASO_DLL_API
void Paso_Options_show(const Paso_Options* options);

PASO_DLL_API
void Paso_Options_showDiagnostics(const Paso_Options* options);
const char* Paso_Options_name(const index_t key);
index_t Paso_Options_getPackage(index_t solver,index_t package, bool_t symmetry, Paso_MPIInfo *mpi_info);

index_t Paso_Options_getSolver(index_t solver,index_t package, bool_t symmetry, Paso_MPIInfo *mpi_info);

#define Paso_Options_copy(in,out) memcpy((Paso_Options*)out,(Paso_Options*)in,sizeof(Paso_Options))

#endif /* #ifndef INC_PASO_OPTIONS */
