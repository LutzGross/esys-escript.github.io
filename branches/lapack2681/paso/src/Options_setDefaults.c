
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

/*   Paso: solver options */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004 */
/*   author: l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "Options.h"

/**************************************************************/

/* set the default values for solver options */

void Paso_Options_setDefaults(Paso_Options* options) {
  options->verbose=FALSE;
  options->method=PASO_DEFAULT;
  options->package=PASO_DEFAULT;
  options->symmetric=FALSE;
  options->reordering=PASO_NO_REORDERING;
  options->tolerance=1.E-8;
  options->absolute_tolerance=0.;
  options->inner_tolerance=0.9;
  options->adapt_inner_tolerance=TRUE;
  options->preconditioner=PASO_JACOBI;
  options->iter_max=10000;
  options->inner_iter_max=10;
  options->drop_tolerance=0.01;
  options->drop_storage=2.;
  options->restart=-1;
  options->truncation=20;
  options->sweeps=2;
  options->pre_sweeps=2;
  options->post_sweeps=2;
  options->coarsening_threshold=0.05;
  options->min_coarse_matrix_size=500;
  options->level_max=5;
  options->accept_failed_convergence=FALSE;
  options->coarsening_method=PASO_DEFAULT;
  options->relaxation_factor=0.95;

  /* diagnostic values */
  options->num_iter=-1;
  options->num_level=-1;
  options->num_inner_iter=-1;
  options->time=-1.;
  options->set_up_time=-1.;
  options->net_time=-1.;
  options->residual_norm=-1.;
  options->converged=FALSE;
}
void Paso_Options_showDiagnostics(const Paso_Options* options) {
	printf("Paso diagonsitics:\n");
	printf("\tnum_iter = %d\n",options->num_iter);
	printf("\tnum_level = %d\n",options->num_level);
	printf("\tnum_inner_iter = %d\n",options->num_inner_iter);
	printf("\ttime = %e\n",options->time);
	printf("\tset_up_time = %e\n",options->set_up_time);
        printf("\tnet_time = %e\n",options->net_time);
	printf("\tresidual_norm = %e\n",options->residual_norm);
	printf("\tconverged = %d\n",options->converged);
}
const char* Paso_Options_name(const index_t key){
    switch (key) {
       case  PASO_DEFAULT:
          return "DEFAULT";
       case  PASO_DIRECT:
          return "DIRECT";
       case  PASO_CHOLEVSKY:
          return "CHOLEVSKY";
       case  PASO_PCG:
          return "PCG";
       case  PASO_CR:
          return "CR";
       case  PASO_CGS:
          return "CGS";
       case  PASO_BICGSTAB:
          return "BICGSTAB";
       case  PASO_SSOR:
          return "SSOR";
       case  PASO_ILU0:
          return "ILU0";
       case  PASO_ILUT:
          return "ILUT";
       case  PASO_JACOBI:
          return "JACOBI";
       case  PASO_GMRES:
          return "GMRES";
       case  PASO_PRES20:
          return "PRES20";
       case  PASO_LUMPING:
          return "LUMPING";
       case  PASO_NO_REORDERING:
          return "NO_REORDERING";
       case  PASO_MINIMUM_FILL_IN:
          return "MINIMUM_FILL_IN";
       case  PASO_NESTED_DISSECTION:
          return "NESTED_DISSECTION";
       case  PASO_MKL:
          return "MKL";
       case  PASO_UMFPACK:
          return "UMFPACK";
       case  PASO_ITERATIVE:
          return "ITERATIVE";
       case  PASO_PASO:
          return "PASO";
       case  PASO_AMG:
          return "AMG";
       case  PASO_REC_ILU:
          return "REC_ILU";
       case  PASO_TRILINOS:
          return "TRILINOS";
       case  PASO_NONLINEAR_GMRES:
          return "NONLINEAR_GMRES";
       case  PASO_TFQMR :
          return "TFQMR";
       case  PASO_MINRES:
          return "MINRES";
       case  PASO_GAUSS_SEIDEL:
          return "GAUSS_SEIDEL";
       case  PASO_RILU:
          return "RILU";
       case  PASO_DEFAULT_REORDERING:
          return "DEFAULT_REORDERING";
       case  PASO_SUPER_LU:
          return "SUPER_LU";
       case  PASO_PASTIX:
          return "PASTIX";
       case  PASO_YAIR_SHAPIRA_COARSENING:
          return "YAIR_SHAPIRA_COARSENING";
       case  PASO_RUGE_STUEBEN_COARSENING:
          return "RUGE_STUEBEN_COARSENING";
       case  PASO_AGGREGATION_COARSENING:
          return "AGGREGATION_COARSENING";
       case  PASO_NO_PRECONDITIONER:
          return "NO_PRECONDITIONER";
       case  PASO_MIN_COARSE_MATRIX_SIZE:
          return "MIN_COARSE_MATRIX_SIZE";
       default:
	  return "<unknown>";
    }
}
void Paso_Options_show(const Paso_Options* options ) {
	printf("Paso options settings:\n");
	printf("\tverbose = %d\n",options->verbose);
	printf("\tmethod = %s (%d)\n",Paso_Options_name(options->method),options->method);
	printf("\tpackage = %s (%d)\n",Paso_Options_name(options->package),options->package);
	printf("\tsymmetric = %d\n",options->symmetric);
	printf("\treordering = %s (%d)\n",Paso_Options_name(options->reordering),options->reordering);
	printf("\ttolerance = %e\n",options->tolerance);
	printf("\tabsolute_tolerance = %e\n",options->absolute_tolerance);
	printf("\tinner_tolerance = %e\n",options->inner_tolerance);
	printf("\tadapt_inner_tolerance = %d\n",options->adapt_inner_tolerance);
	printf("\tpreconditioner =  %s (%d)\n",Paso_Options_name(options->preconditioner),options->preconditioner);
	printf("\titer_max = %d\n",options->iter_max);
	printf("\tinner_iter_max = %d\n",options->inner_iter_max);
	printf("\tdrop_tolerance = %e\n",options->drop_tolerance);
	printf("\tdrop_storage = %e\n",options->drop_storage);
	printf("\trestart = %d\n",options->restart);
	printf("\ttruncation = %d\n",options->truncation);
	printf("\tsweeps = %d\n",options->sweeps);
	printf("\tpre_sweeps = %d\n",options->pre_sweeps);
	printf("\tpost_sweeps = %d\n",options->post_sweeps);
	printf("\tcoarsening_threshold = %e\n",options->coarsening_threshold);
	printf("\tlevel_max = %d\n",options->level_max);
	printf("\taccept_failed_convergence = %d\n",options->accept_failed_convergence);
	printf("\tcoarsening_method = %s (%d)\n",Paso_Options_name(options->coarsening_method), options->coarsening_method);
	printf("\trelaxation_factor = %e\n",options->relaxation_factor);
}
