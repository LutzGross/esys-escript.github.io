/* $Id$ */

/**********************************************************************/

/* Finley: SystemMatrixPattern: Finley_SystemMatrixPattern_mis 

   searches for a maximal independent set MIS in the matrix pattern 
   vertices in the maximal independent set are marked in mis_marker

*/
/**********************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005              */
/* Author: gross@access.edu.au                                */

/**************************************************************/

#include "Finley.h"
#include "System.h"
#include "Util.h"


/* used to generate pseudo random numbers: */

static double Finley_SystemMatrixPattern_mis_seed=.4142135623730951;


/***************************************************************/
 
#define IS_AVAILABLE -1
#define IS_IN_MIS_NOW -2
#define IS_IN_MIS -3
#define IS_CONNECTED_TO_MIS -4

void Finley_SystemMatrixPattern_mis(Finley_SystemMatrixPattern* pattern_p, maybelong* mis_marker) {

  maybelong i,naib,iptr;
  int flag;
  maybelong n=pattern_p->n_ptr;
  double *value=TMPMEMALLOC(n,double);
  if (!Finley_checkPtr(value)) {
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<n;++i) mis_marker[i]=IS_AVAILABLE;
   
     /* is there any vertex available ?*/
     while (Finley_Util_isAny(n,mis_marker,IS_AVAILABLE)) {
        /* step 1: assign a random number in [0,1] to each vertex */
        /* step 2: is the vertex is available, check if its value is the smaller than all values of its naigbours 
                   if the answer is yes, the vertex is put into the independend set and all 
                   its naighbours are removed from the graph by setting it mis_marker to FALSE */
      
        /* step 2: is the vertex is available, check if its value is the smaller than all values of its naigbours */

           /* assign random number in [0,1] to each vertex */
           #pragma omp parallel for private(i) schedule(static)
           for (i=0;i<n;++i) {
                 if (mis_marker[i]==IS_AVAILABLE) {
                    value[i]=fmod(Finley_SystemMatrixPattern_mis_seed*(i+1),1.);
                 } else {
                    value[i]=2.;
                 }
           }
           /* update the seed */
           Finley_SystemMatrixPattern_mis_seed=fmod(sqrt(Finley_SystemMatrixPattern_mis_seed*(n+1)),1.);
           /* detect independent vertices as those vertices that have a value less than all values of its naigbours */
           #pragma omp parallel for private(naib,i,iptr,flag) schedule(static) 
           for (i=0;i<n;++i) {
              if (mis_marker[i]==IS_AVAILABLE) {
                 flag=IS_IN_MIS_NOW;
                 for (iptr=pattern_p->ptr[i];iptr<pattern_p->ptr[i+1]; ++iptr) {
                     naib=pattern_p->index[iptr];
                     if (naib!=i && value[naib]<=value[i]) {
                        flag=IS_AVAILABLE;
                        break;
                     }
                 }
                 mis_marker[i]=flag;
              }
           }
           /* detect independent vertices as those vertices that have a value less than all values of its naigbours */
           #pragma omp parallel for private(naib,i,iptr) schedule(static)
           for (i=0;i<n;i++) {
              if (mis_marker[i]==IS_IN_MIS_NOW) {
                 for (iptr=pattern_p->ptr[i];iptr<pattern_p->ptr[i+1]; ++iptr) {
                     naib=pattern_p->index[iptr];
                     if (naib!=i) mis_marker[naib]=IS_CONNECTED_TO_MIS;
                 }
                 mis_marker[i]=IS_IN_MIS;
              }
           }
     }
     /* swap to TRUE/FALSE in mis_marker */
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<n;i++) mis_marker[i]=(mis_marker[i]==IS_IN_MIS);
  }
  TMPMEMFREE(value);
}
#undef IS_AVAILABLE 
#undef IS_IN_MIS_NOW 
#undef IS_IN_MIS 
#undef IS_CONNECTED_TO_MIS 

/*
 * $Log$
 * Revision 1.2  2005/03/04 07:12:46  jgs
 * *** empty log message ***
 *
 * Revision 1.1.2.1  2005/03/02 23:35:06  gross
 * reimplementation of the ILU in Finley. block size>1 still needs some testing
 *
 *
 */
