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
static int TokenSym[]={FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE};
#endif

/**************************************************************/

/*  free any extra stuff possibly used by the SCSL library */

void Finley_SCSL_direct_free(Finley_SystemMatrix* A) {
#if DIRECT_SOLVER == SGI_SCSL
      if (A->direct!=NULL) {
	int token=*(int*)(A->direct);
        TokenList[token]=0;
        if (TokenSym[token]) {
            DPSLDLT_Destroy(token);
        } else {
            DPSLDU_Destroy(token);
        }
        A->direct=NULL;
     }
#endif
}
/*  call the iterative solver: */

void Finley_SCSL_direct(Finley_SystemMatrix* A,
                               double* out,
                               double* in,
                               Finley_SolverOptions* options) {
#if DIRECT_SOLVER == SGI_SCSL
  double time0;
  int token,l=sizeof(TokenList)/sizeof(int),reordering_method,method;
  long long non_zeros;
  double ops;
  if (options->symmetric) {
      method=ESCRIPT_CHOLEVSKY;
  } else {
      method=ESCRIPT_DIRECT;
  }
  if (A->col_block_size!=1) {
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"SCSL linear solver can only be applied to block size 1.");
      return;
  }
  if (method==ESCRIPT_CHOLEVSKY) {
      if (A->type!=CSC_SYM) {
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"SGI SCSL direct solver for symmetric matrices can only be applied to symmetric CSC format.");
          return;
      }
  } else {
      if (A->type!=CSC) {
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"SGI SCSL direct solver can only be applied to CSC format.");
          return;
      }
  }
  /* if no token has been assigned a free token must be found*/
  if (A->direct==NULL) {
     /* find the next available token */
    for(token=0;token<l && TokenList[token]!=0;token++);
     if (token==l) {
         Finley_ErrorCode=TYPE_ERROR;
         sprintf(Finley_ErrorMsg,"SCSL linear solver can only deal with %d matrices simultaneously.",l);
         return;
     } else {
       TokenList[token] = 1;
       TokenSym[token]=(method==ESCRIPT_CHOLEVSKY);
       A->direct=&Token[token];
       /* map the reordering method onto SCSL */
       switch (options->reordering) {
             case ESCRIPT_NO_REORDERING:
                  reordering_method=0;
                  break;
             case ESCRIPT_MINIMUM_FILL_IN:
                  reordering_method=1;
                  break;
             default:
                  reordering_method=2;
                  break;
       }
       time0=Finley_timer();
       if (TokenSym[token]) {
            DPSLDLT_Ordering(token,reordering_method);
            DPSLDLT_Preprocess(token,A->num_rows,A->pattern->ptr,A->pattern->index,&non_zeros,&ops);
            DPSLDLT_Factor(token,A->num_rows,A->pattern->ptr,A->pattern->index,A->val);
            if (options->verbose) printf("timing SCSL: Cholevsky factorization: %.4e sec (token = %d)\n",Finley_timer()-time0,token);
       } else {
            DPSLDU_Ordering(token,reordering_method);
            DPSLDU_Preprocess(token,A->num_rows,A->pattern->ptr,A->pattern->index,&non_zeros,&ops);
            DPSLDU_Factor(token,A->num_rows,A->pattern->ptr,A->pattern->index,A->val);
            if (options->verbose) printf("timing SCSL: LDU factorization: %.4e sec (token = %d)\n",Finley_timer()-time0,token);
       }
    }
  }
  time0=Finley_timer();
  token=*(int*)(A->direct);
  if (TokenSym[token]) {
     DPSLDLT_Solve(token,out,in);
  } else {
     DPSLDU_Solve(token,out,in);
  }
  if (options->verbose) printf("timing SCSL: solve: %.4e sec (token = %d)\n",Finley_timer()-time0,token);
#else
    Finley_ErrorCode=SYSTEM_ERROR;
    sprintf(Finley_ErrorMsg,"No SCSL direct solver available.");
#endif
}
/*
 * $Log$
 * Revision 1.2  2004/12/15 07:08:34  jgs
 * *** empty log message ***
 *
 * Revision 1.1.2.1  2004/11/12 06:58:20  gross
 * a lot of changes to get the linearPDE class running: most important change is that there is no matrix format exposed to the user anymore. the format is chosen by the Domain according to the solver and symmetry
 *
 *
 */
