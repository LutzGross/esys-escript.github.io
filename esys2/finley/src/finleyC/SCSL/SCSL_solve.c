/* $Id$ */

/**************************************************************/

/* Finley: interface to the direct solvers                    */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "System.h"
#include "SCSL.h"
#if DIRECT_SOLVER == SGI_SCSL
#include <scsl_sparse.h>
static int TokenList[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
static int Token[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
#endif

/**************************************************************/

/*  free any extra stuff possibly used by the SCSL library */

void Finley_SCSL_solve_free(Finley_SystemMatrix* A) {
#if DIRECT_SOLVER == SGI_SCSL
      if (A->solve!=NULL) {
	int token=*(int*)(A->solve);
          TokenList[token]=0;
          if (A->symmetric) {
             DPSLDLT_Destroy(token);
          } else {
             DPSLDU_Destroy(token);
          }
          A->solve=NULL;
      }
#endif
}
/*  call the iterative solver: */

void Finley_SCSL_solve(Finley_SystemMatrix* A,
                               double* out,
                               double* in,
                               Finley_SolverOptions* options) {
#if DIRECT_SOLVER == SGI_SCSL
  double time0;
  int k,l=sizeof(TokenList)/sizeof(int),method;
  long long non_zeros;
  double ops;
  if (A->col_block_size!=1) {
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"SCSL linear solver can only be applied to block size 1.");
      return;
  }
  if (A->type!=CSC) {
      Taipan_ErrorCode=TYPE_ERROR;
      sprintf(Taipan_ErrorMsg,"SGI SCSL direct solver can only be applied to CSC format.");
      return;
  }
  #ifdef Finley_TRACE
  printf("solver is SCSL\n");
  #endif
  /* if no token has been assigned a free token must be found*/
  if (A->solve==NULL) {
     /* find the next available token */
    for(k=0;k<l && TokenList[k]!=0;k++);
     if (k==l) {
         Finley_ErrorCode=TYPE_ERROR;
         sprintf(Finley_ErrorMsg,"SCSL linear solver can only deal with %d matrices simultaneously.",l);
         return;
     } else {
       TokenList[k] = 1;
       A->solve=&Token[k];
       #ifdef Finley_TRACE
       printf("solving system %d\n",*(int*)(A->solve));
       #endif
       /* map the reordering method onto SCSL */
       switch (options->reordering) {
             case NO_REORDERING:
                  method=0;
                  break;
             case MINIMUM_FILL_IN:
                  method=1;
                  break;
             default:
                  method=2;
                  break;
       }
       time0=Finley_timer();
       if (A->symmetric) {
            DPSLDLT_Ordering(*(int*)(A->solve),method);
            DPSLDLT_Preprocess(*(int*)(A->solve),A->num_rows,A->ptr,A->index,&non_zeros,&ops);
            DPSLDLT_Factor(*(int*)(A->solve),A->num_rows,A->ptr,A->index,A->val);
       } else {
            DPSLDU_Ordering(*(int*)(A->solve),method);
            DPSLDU_Preprocess(*(int*)(A->solve),A->num_rows,A->ptr,A->index,&non_zeros,&ops);
            DPSLDU_Factor(*(int*)(A->solve),A->num_rows,A->ptr,A->index,A->val);
       }
       printf("timing SCSL: factorization: %.4e sec\n",Finley_timer()-time0);
    }
  }
  time0=Finley_timer();
  if (A->symmetric) {
     DPSLDLT_Solve(*(int*)(A->solve),out,in);
  } else {
     DPSLDU_Solve(*(int*)(A->solve),out,in);
  }
  printf("timing SCSL: solve: %.4e sec\n",Finley_timer()-time0);
#endif
}
/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:57  jgs
 * Initial revision
 *
 * Revision 1.1.2.2  2004/10/26 06:36:44  jgs
 * *** empty log message ***
 *
 * Revision 1.2  2004/10/13 01:53:43  gross
 * bug in CSC assembling fixed
 *
 * Revision 1.1  2004/07/02 04:21:14  gross
 * Finley C code has been included
 *
 *
 */
