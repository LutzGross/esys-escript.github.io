/* $Id$ */

/**********************************************************************/

/* Finley: Solver: set the pointer of the main diagonal in A_p->value */

/**********************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "System.h"
#include "Util.h"
#include "Solver.h"


/* used to generate pseudo random numbers: */
static double Finley_Solver_coloring_seed=.4142135623730951;

#define CONNECTED_TO_CURRENT_COLOR -1
#define NOT_COLORED_NOW -2
#define MAYBE_COLORED_NOW -3

/**************************************************************/

/* inverse preconditioner setup */

void Finley_Solver_coloring(Finley_SystemMatrix* A_p,maybelong* numColors, maybelong* color) {
#if ITERATIVE_SOLVER == NO_LIB
  int irow,iptr,currentColor,naib;
  maybelong n=A_p->num_rows;
  double *value=(double*) TMPMEMALLOC(n*sizeof(double));
  if (!Finley_checkPtr(value)) {
       /* no colors used yet */
       currentColor=0;
       /* initialize color */
       #pragma omp parallel for private(irow) schedule(static)
       for (irow=0;irow<n;irow++) color[irow]=CONNECTED_TO_CURRENT_COLOR;

       while (Finley_Util_isAny(n,color,CONNECTED_TO_CURRENT_COLOR)) {
          /* all rows which are not colored yet (color[irow]<0) are indicated 
               to be free for the next color assigment*/
          #pragma omp parallel for private(irow) schedule(static)
          for (irow=0;irow<n;irow++) if (color[irow]==CONNECTED_TO_CURRENT_COLOR) color[irow]=MAYBE_COLORED_NOW;
          /* find rows marked with MAYBE_COLORED_NOW which are not connected to the colored rows but could be connected */
          /* to CONNECTED_TO_CURRENT_COLOR rows which are the rows connected to colored rows */
          /* the following while loop finds MAYBE_COLORED_NOW rows which could be linked to CONNECTED_TO_CURRENT_COLOR rows but */
          /* not connected rows already colors.                                                                      */
          while (Finley_Util_isAny(n,color,MAYBE_COLORED_NOW)) {
	     #pragma omp parallel
	     {
                /* assign a random number to all MAYBE_COLORED_NOW or currentColor colored row: */
                /* non-MAYBE_COLORED_NOW rows get value 2. which is bigger than all other values */
                /* which means that they have no effect when checking for minimum values in the naigbourhood */
                #pragma omp for private(irow) schedule(static)
                for (irow=0;irow<n;irow++) {
                   if (color[irow]==MAYBE_COLORED_NOW) {
                       value[irow]=fmod(Finley_Solver_coloring_seed*(irow+1),1.);
                   } else {
                       value[irow]=2.;
                   }
                }
                /* update the seed */
                #pragma omp master
                Finley_Solver_coloring_seed=fmod(sqrt(Finley_Solver_coloring_seed*(n+1)),1.);
                /* if a row is MAYBE_COLORED_NOW and its value is smaller than all the values of its naigbours */
                /* the row is independent from the un-colored rows */
                /* otherwise the row is connected to another row and the row is marked as NOT_COLORED_NOW */
                /* i.e. it may use in the next choice for value if it is not marked as CONNECTED_TO_CURRENT_COLOR */
                #pragma omp for private(naib,irow,iptr) schedule(static)
                for (irow=0;irow<n;irow++) {
                   if (color[irow]==MAYBE_COLORED_NOW) {
                     for (iptr=A_p->ptr[irow];iptr<A_p->ptr[irow+1]; iptr++) {
                          naib=A_p->index[iptr];
                          if (naib!=irow && value[naib]<=value[irow]) {
                               color[irow]=NOT_COLORED_NOW;
                               break;
                          }
                     }
                   }
                }
                /* all naigbours of MAYBE_COLORED_NOW row which are not colored yet are marked as CONNECTED_TO_CURRENT_COLOR */
                /* They cannot be used for this coloring round */
                #pragma omp for private(naib,iptr,irow) schedule(static)
                for (irow=0;irow<n;irow++) {
                   if (color[irow]==MAYBE_COLORED_NOW) {
                     for (iptr=A_p->ptr[irow];iptr<A_p->ptr[irow+1]; iptr++) {
                          naib=A_p->index[iptr];
                          if (naib!=irow && color[naib]<0) color[naib]=CONNECTED_TO_CURRENT_COLOR;
                     }
                     color[irow]=currentColor;
                   } else if (color[irow]==NOT_COLORED_NOW) {
                     color[irow]=MAYBE_COLORED_NOW;
                   }
                }
            }
	  } 
          /* get to the next color */
          currentColor++;
       }
  }
  *numColors=currentColor;

  #ifdef FINLEY_SOLVER_TRACE
  printf("row coloring finalized. %d colors used.\n",*numColors);
  #endif
  TMPMEMFREE(value);
#endif
}
#undef CONNECTED_TO_CURRENT_COLOR
#undef MAYBE_COLORED_NOW
#undef NOT_COLORED_NOW

/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:58  jgs
 * Initial revision
 *
 * Revision 1.1  2004/07/02 04:21:14  gross
 * Finley C code has been included
 *
 *
 */
