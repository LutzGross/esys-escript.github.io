
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
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


/************************************************************************************/

/* Paso: SystemMatrix: interface to HYPRE, a software library of
   high performance preconditioners and solvers. We use the 
   BoomerAMG provided in this library */

/************************************************************************************/

/* Copyrights by ACcESS Australia 2011 */
/* Author: Lin Gao, l.gao@uq.edu.au */

/************************************************************************************/

#include "Paso.h"
#include "BOOMERAMG.h"
#include "Preconditioner.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/************************************************************************************/

/*  free any extra stuff possibly used by the BoomerAMG library */
void Paso_Preconditioner_BoomerAMG_free(Paso_Preconditioner_BoomerAMG * in) {
#ifdef BOOMERAMG
     if (in != NULL) {
	HYPRE_IJMatrixDestroy(in->pt->A);
	HYPRE_IJVectorDestroy(in->pt->b);
	HYPRE_IJVectorDestroy(in->pt->x);

        /* Destroy solver */
        HYPRE_BoomerAMGDestroy(in->pt->solver);
	MEMFREE(in->pt);
     }
#endif
}


/* allocate necessary space for the BoomerAMG library */
Paso_Preconditioner_BoomerAMG* Paso_Preconditioner_BoomerAMG_alloc(
			Paso_SystemMatrix * A, Paso_Options* options)
{
#ifdef BOOMERAMG
     index_t ilower; /* first row in current processer, number is given by
                           the global indices. Can be 0- or 1-based indexing */
     index_t iupper; /* last row in current processer, number is given by
                           the global indices. Row partitioning must be 
                           contiguous, i.e., iupper for proc i must equal 
                           ilower-1 for proc i+1 */
     index_t jlower; /* jlower and jupper define a column partitioning */
     index_t jupper;
     index_t nrows;
     dim_t offset;
     index_t i;
     index_t *ncols, *rows, *cols;
     double *val;
     Paso_Preconditioner_BoomerAMG* out=NULL;
     Paso_BOOMERAMG_Handler* pt=NULL;

     if (A->type & MATRIX_FORMAT_OFFSET1 ||
         A->mainBlock->col_block_size != 1 ||
         A->mainBlock->row_block_size != 1 ||
         ((A->type != MATRIX_FORMAT_DEFAULT) &&
          (A->type != (MATRIX_FORMAT_DEFAULT+MATRIX_FORMAT_BLK1))
         )
        ) {
        printf("Paso_Preconditioner_BoomerAMG: BoomerAMG requires CSR format with index offset 0 and block size 1.");
        return NULL;
     }

     out=MEMALLOC(1,Paso_Preconditioner_BoomerAMG);
     if (Esys_checkPtr(out)) return NULL;

     pt=(MEMALLOC(1,Paso_BOOMERAMG_Handler));
     if (Esys_checkPtr(pt)) return NULL;

     out->pt = pt;

     /* set up inputs for BoomerAMG */
     nrows = A->mainBlock->numRows;

     ilower = A->row_distribution->first_component[A->mpi_info->rank];
     iupper = ilower + nrows - 1;

     jlower = A->col_distribution->first_component[A->mpi_info->rank];
     jupper = ilower + A->mainBlock->numCols - 1;

     rows = MEMALLOC(nrows,index_t);
     #pragma omp parallel for schedule(static) private(i)
     for (i=0; i<nrows; i++) rows[i]=i+ilower;

     /* merge MainBlock and CoupleBlock to get a complete column 
	entries for each row allocated to current processor */
     Paso_SystemMatrix_mergeMainAndCouple(A, &ncols, &cols, &val);

     /* ncols returned from mergeMainAndCouple() is a vector of 
	locations that start a row. we transfer it into a vector
	of the number of columns in each row */
     for (i=0; i<nrows; i++) ncols[i]=ncols[i+1]-ncols[i];

     /* create and initialize an empty matrix object. _IJMatrixCreate()  
	is a collective call with each process passing its own row extents, 
	ilower and iupper */
     HYPRE_IJMatrixCreate(A->mpi_info->comm, ilower, iupper, jlower, jupper, &(out->pt->A));
     HYPRE_IJMatrixSetObjectType(out->pt->A, HYPRE_PARCSR);
     HYPRE_IJMatrixInitialize(out->pt->A);

     /* set matrix coefficients. set matrix values for some number of 
	rows (nrows) and some number of columns in each row (ncols).
	The actual row and column numbers of the matrix values to be
	set are given by "rows" and "cols" */
     HYPRE_IJMatrixSetValues(out->pt->A, nrows, ncols, rows, cols, val);

     /* finalize the matrix assembly and make the matrix "ready to
	use". _IJMatrixAssemble() is a collective call */
     HYPRE_IJMatrixAssemble(out->pt->A);
     HYPRE_IJMatrixGetObject(out->pt->A, (void **) &(out->pt->parcsr_A));

     /* Create the rhs and solution */
     HYPRE_IJVectorCreate(A->mpi_info->comm, ilower, iupper,&(out->pt->b));
     HYPRE_IJVectorSetObjectType(out->pt->b, HYPRE_PARCSR);
     HYPRE_IJVectorInitialize(out->pt->b);
     HYPRE_IJVectorCreate(A->mpi_info->comm, ilower, iupper,&(out->pt->x));
     HYPRE_IJVectorSetObjectType(out->pt->x, HYPRE_PARCSR);
     HYPRE_IJVectorInitialize(out->pt->x);

     /* create solver */
     HYPRE_BoomerAMGCreate(&(out->pt->solver));

     /* set some parameters */
     /* define the smoother to be used (default is 3): 
	0 -- Jacobi
	1 -- Gauss-Seidel, sequential (very slow!)
	2 -- Gauss-Seidel, points in parallel, boundary sequential (slow!)
	3 -- hybrid Gauss-Seidel or SOR, forward solve
	4 -- hybrid Gauss-Seidel or SOR, backward solve
	5 -- hybrid chaotic Gauss-Seidel (works only with OpenMP)
	6 -- hybrid symmetric Gauss-Seidel or SSOR
	9 -- Gaussian elimination (only on coarsest level) */
     switch (options->smoother) {
	case PASO_GS:
	default:
	  HYPRE_BoomerAMGSetRelaxType(out->pt->solver, 3);
	  break;
        case PASO_JACOBI:
          HYPRE_BoomerAMGSetRelaxType(out->pt->solver, 0);
          break;
     }
     /* parallel coarsening algorithm used (default is 6):
	0 -- CLJP-coarsening (a parallel coarsening algorithm using 
	     independent sets.
	1 -- classical Ruge-Stueben coarsening on each processor, no 
	     boundary treatment (not recommended!)
	3 -- classical Ruge-Stueben coarsening on each processor, followed 
	     by a third pass, which adds coarse points on the boundaries
	6 -- Falgout coarsening (uses 1 first, followed by CLJP using the 
	     interior coarse points generated by 1 as its first 
	     independent set)
	7 -- CLJP-coarsening (using a fixed random vector, for debugging 
	     purposes only)
	8 -- PMIS-coarsening (a parallel coarsening algorithm using 
	     independent sets, generating lower complexities than CLJP, 
	     might also lead to slower convergence)
	9 -- PMIS-coarsening (using a fixed random vector, for debugging 
	     purposes only)
	10 -- HMIS-coarsening (uses one pass Ruge-Stueben on each processor 
	     independently, followed by PMIS using the interior C-points 
	     generated as its first independent set)
	11 -- one-pass Ruge-Stueben coarsening on each processor, no 
	     boundary treatment (not recommended!) */
     switch (options->coarsening_method) {
	case PASO_FALGOUT_COARSENING:
        default:
          HYPRE_BoomerAMGSetCoarsenType(out->pt->solver, 6);
          break;
	case PASO_RUGE_STUEBEN_COARSENING:
	  HYPRE_BoomerAMGSetCoarsenType(out->pt->solver, 3);
	  break;
	case PASO_CIJP_COARSENING:
	  HYPRE_BoomerAMGSetCoarsenType(out->pt->solver, 0);
	  break;
	case PASO_CIJP_FIXED_RANDOM_COARSENING:
	  HYPRE_BoomerAMGSetCoarsenType(out->pt->solver, 7);
	  break;
	case PASO_PMIS_COARSENING:
	  HYPRE_BoomerAMGSetCoarsenType(out->pt->solver, 8);
	  break;
	case PASO_HMIS_COARSENING:
	  HYPRE_BoomerAMGSetCoarsenType(out->pt->solver, 10);
	  break;
     }

     if (options->verbose)
        HYPRE_BoomerAMGSetPrintLevel(out->pt->solver, 3);  
     /* maximum number of levels. Default is 25 */
     if (options->level_max > 0) 
	HYPRE_BoomerAMGSetMaxLevels(out->pt->solver, options->level_max);
     /* number of is sweeps */
     if (options->sweeps > 0)
	HYPRE_BoomerAMGSetNumSweeps(out->pt->solver, options->sweeps);
     /* 1 for V-cycle and 2 for W-cycle */
     if (options->cycle_type == 1 || options->cycle_type == 2)
	HYPRE_BoomerAMGSetCycleType(out->pt->solver, options->cycle_type);
     /* AMG coarsening strength threshold. Default is 0.25. For 2D Laplace 
	operations, 0.25 is a good value. For 3D Laplace operations, 0.5 or 
	0.6 is a better value. For elasticity problems, a large strength 
	threshold, such as 0.9, is oftern better */
     if (options->coarsening_threshold > 0)
	HYPRE_BoomerAMGSetStrongThreshold(out->pt->solver, options->coarsening_threshold); 
     /* threshold for diagonal dominant rows. Default is 0.9 */
     if (options->diagonal_dominance_threshold > 0)
	HYPRE_BoomerAMGSetMaxRowSum(out->pt->solver, options->diagonal_dominance_threshold); 
     /* conv. tolerance */
     /* HYPRE_BoomerAMGSetTol(out->pt->solver, 1e-7);*/ 

     /* free memory allocation*/
     MEMFREE(rows);
     MEMFREE(ncols);
     MEMFREE(cols);
     MEMFREE(val);
     return out;
#else
     return NULL;
#endif
}


/*  call the solver: */
void Paso_Preconditioner_BoomerAMG_solve(Paso_SystemMatrix* A,
			  Paso_Preconditioner_BoomerAMG * amg,
                          double* out,
                          double* in)
{
#ifdef BOOMERAMG
     index_t ilower; /* first row in current processer, number is given by
                           the global indices. Can be 0- or 1-based indexing */
     index_t iupper; /* last row in current processer, number is given by
			   the global indices. Row partitioning must be 
			   contiguous, i.e., iupper for proc i must equal 
			   ilower-1 for proc i+1 */
     index_t i, nrows;
     dim_t offset;
     index_t *rows;

     if (A->type & MATRIX_FORMAT_OFFSET1 ||
	 A->mainBlock->col_block_size != 1 ||
	 A->mainBlock->row_block_size != 1 ||
	 ((A->type != MATRIX_FORMAT_DEFAULT) && 
	  (A->type != (MATRIX_FORMAT_DEFAULT+MATRIX_FORMAT_BLK1))
	 )
	) {
        printf("Paso_Preconditioner_BoomerAMG: BoomerAMG requires CSR format with index offset 0 and block size 1.");
        return;
     }

     {
	/* set up inputs for BoomerAMG */
	nrows = A->mainBlock->numRows;
	ilower = A->row_distribution->first_component[A->mpi_info->rank];
	iupper = ilower + nrows - 1;
	rows = MEMALLOC(nrows,index_t);
	#pragma omp parallel for schedule(static) private(i)
	for (i=0; i<nrows; i++) rows[i]=i+ilower;

	/* set values for the rhs and solution */
	HYPRE_IJVectorSetValues(amg->pt->b, nrows, rows, in);
	HYPRE_IJVectorAssemble(amg->pt->b);
	HYPRE_IJVectorGetObject(amg->pt->b, (void **) &(amg->pt->par_b));
	HYPRE_IJVectorSetValues(amg->pt->x, nrows, rows, out);
	HYPRE_IJVectorAssemble(amg->pt->x);
        HYPRE_IJVectorGetObject(amg->pt->x, (void **) &(amg->pt->par_x));

	/* setup and solve */
	HYPRE_BoomerAMGSetup(amg->pt->solver, amg->pt->parcsr_A, amg->pt->par_b, amg->pt->par_x);
	HYPRE_BoomerAMGSolve(amg->pt->solver, amg->pt->parcsr_A, amg->pt->par_b, amg->pt->par_x);

        HYPRE_IJVectorGetValues(amg->pt->x, nrows, rows, out);

	/* free memory allocation*/
	MEMFREE(rows); 
     }
#endif
}
