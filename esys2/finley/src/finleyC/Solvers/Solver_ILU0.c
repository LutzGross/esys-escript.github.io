/* $Id$ */

/**************************************************************/

/* Finley: ILU(0) preconditioner: */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "System.h"
#include "Solver.h"

void Finley_Solver_setILU0_1(int,maybelong*,maybelong*,double*,maybelong,maybelong*,maybelong*);
void Finley_Solver_solveILU0_1(int,double*,maybelong*,maybelong*,double*,maybelong,maybelong*,maybelong*);


/**************************************************************/

/* inverse preconditioner setup */

void Finley_Solver_setILU0(Finley_SystemMatrix * A_p) {
#if ITERATIVE_SOLVER == NO_LIB
  Finley_Solver_Preconditioner* prec=(Finley_Solver_Preconditioner*) (A_p->iterative);
  maybelong n=A_p->num_rows;
  /* allocate arrays */
  prec->mainDiag =MEMALLOC(n,maybelong);
  prec->color = MEMALLOC(n,maybelong);
  prec->values = MEMALLOC((A_p->len)*(size_t) n,double);
  if (! (Finley_checkPtr(prec->mainDiag) || Finley_checkPtr(prec->color) || Finley_checkPtr(prec->values))) {
     /* find the main diagonal: */
     Finley_Solver_getMainDiagonal(A_p,prec->mainDiag);
     /* color the rows : */
     Finley_Solver_coloring(A_p->pattern,&(prec->numColors),prec->color);
     if (Finley_ErrorCode==NO_ERROR) {
         /* allocate vector to hold main diagonal entries: */
         Finley_SystemMatrix_copy(A_p,prec->values);
         /* allocate vector to hold main diagonal entries: */
         if (A_p->row_block_size > 1) {
           Finley_ErrorCode = TYPE_ERROR;
           sprintf(Finley_ErrorMsg, "ILU preconditioner requires block size =1 .");
           return;
         } else {
           Finley_Solver_setILU0_1(n,A_p->pattern->ptr,A_p->pattern->index,prec->values,prec->numColors,prec->color,prec->mainDiag);
         }
     }
  }
  if (Finley_ErrorCode!=NO_ERROR) Finley_Preconditioner_free(A_p->iterative);
  return;
#endif
} 

/**************************************************************/

/* inverse preconditioner solution using raw pointers */

/* should be called within a parallel region                                              */
/* barrier synconization should be performed to make sure that the input vector available */

void Finley_Solver_solveILU0(Finley_SystemMatrix * A_p, double * x, double * b) {
#if ITERATIVE_SOLVER == NO_LIB
     Finley_Solver_Preconditioner* prec=(Finley_Solver_Preconditioner*) (A_p->iterative);
     maybelong n=A_p->num_rows;
     maybelong i;
     #pragma omp for private(i)
     for (i=0;i<n*A_p->row_block_size;i++) x[i]=b[i];
     if (A_p->row_block_size > 1) {
        /* Finley_Solver_solveILU0_blk(n,x,A_p->pattern->ptr,A_p->pattern->index,prec->values,prec->numColors,prec->color,prec->mainDiag); */
     } else {
         Finley_Solver_solveILU0_1(n,x,A_p->pattern->ptr,A_p->pattern->index,prec->values,prec->numColors,prec->color,prec->mainDiag);
     }
     return;
#endif
}
/**************************************************************/

/* applies ILU preconditioner on x for block size 1 (in place) :  */

void Finley_Solver_solveILU0_1(int n,double* x, maybelong* ptr,maybelong* index,double* values,
                                     maybelong numColors,maybelong* color, maybelong* mainDiag) {
#if ITERATIVE_SOLVER == NO_LIB
      int iColor,i,k,ik;
      /* forward substitution */
      iColor=0;
      while (iColor<numColors) {
         #pragma omp for private(i,k,ik)
         for (i=0;i<n;i++) {
             if (color[i]==iColor) {
                 for (ik=ptr[i]; ik < ptr[i+1]; ik++) {
                    k=index[ik]; 
                    if (color[k]<color[i]) x[i]-=values[ik]*x[k];
                 } /* end of ik-loop */
             } 
         } /* end of i-loop */
         iColor++;
      }
      /* backwards substitution */
      iColor=numColors;
      while (0<iColor) {
         iColor--;
         #pragma omp for private(i,k,ik)
         for (i=0;i<n;i++) {
             if (color[i]==iColor) {
                 for (ik=ptr[i]; ik < ptr[i+1]; ik++) {
                    k=index[ik]; 
                    if (color[k]>color[i]) x[i]-=values[ik]*x[k];
                 } /* end of ik-loop */
                 x[i]*=values[mainDiag[i]];
             }
         }
      }
#endif
}
/**************************************************************/

/* ILU preconditioner for block size 1 :  */

void Finley_Solver_setILU0_1(int n,maybelong* ptr,maybelong* index,double* values,
                                     maybelong numColors,maybelong* color, maybelong* mainDiag) {
#if ITERATIVE_SOLVER == NO_LIB
   double value_ik;
   int iColor,i,k,j,ij,kj,ik,kColor;
   #pragma omp parallel private(iColor)
   {
      for (iColor=0;iColor<numColors;iColor++) {
         #pragma omp for private(i,k,j,ij,kj,ik,kColor,value_ik)
         for (i=0;i<n;i++) {
             if (color[i]==iColor) {
/* printf("processing row = %d\n",i); */
                 for (kColor=0;kColor<color[i];kColor++) {
                    /* run through the row elements a_i,k considering rows have been eleminated only */
                    for (ik=ptr[i]; ik < ptr[i+1]; ik++) {
                       k=index[ik]; 
                       if (color[k]==kColor) {
                          /* up-date the row entry (i,k) */
                          value_ik=values[ik]*values[mainDiag[k]];
                          values[ik]=value_ik;
/* printf("  working on element (%d,%d):%e\n",i,k,values[ik]); */
                          /* run through the entries a_k,j in row k considering rows have not been eleminated only */
                          for (kj=ptr[k];kj<ptr[k+1]; kj++) {
                              j=index[kj]; 
                              if (color[j]>kColor) {
                                 /* is (i,j) in the shape of the matrix? */
                                 /* this means that index[ij]=j */
                                 for (ij=ptr[i];ij<ptr[i+1];ij++) {
                                    if (index[ij]==j) {
                                       values[ij]-=value_ik*values[kj];
/* printf("    updating (%d,%d) by (%d,%d): %e\n",i,j,k,j,values[ij]); */
                                       break;
                                    }
                                 }
                             }
                          } /* end of kj-loop */
                       }
                    } /* end of ik-loop */
                }
                /* invert the main diagonal element */
                if (ABS(values[mainDiag[i]])>0.) {
                   values[mainDiag[i]]=1./values[mainDiag[i]];
                } else {
                   Finley_ErrorCode = ZERO_DIVISION_ERROR;
                   sprintf(Finley_ErrorMsg, "ILU factorization failed in row %d.",i);
                }
/* printf("finalizing row %d : entry %e\n",i,values[mainDiag[i]]);  */
             } /* end of color if */
         } /* end of row loop */
      } /* end of color loop */
   }
#endif
}

/*
 * $Log$
 * Revision 1.2  2004/12/14 05:39:32  jgs
 * *** empty log message ***
 *
 * Revision 1.1.1.1.2.2  2004/11/24 01:37:17  gross
 * some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
 *
 * Revision 1.1.1.1.2.1  2004/11/12 06:58:21  gross
 * a lot of changes to get the linearPDE class running: most important change is that there is no matrix format exposed to the user anymore. the format is chosen by the Domain according to the solver and symmetry
 *
 * Revision 1.1.1.1  2004/10/26 06:53:58  jgs
 * initial import of project esys2
 *
 * Revision 1.1  2004/07/02 04:21:14  gross
 * Finley C code has been included
 *
 *
 */
