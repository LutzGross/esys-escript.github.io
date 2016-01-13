/* $Id$ */

/*
********************************************************************************
*               Copyright   2006 by ACcESS MNRF                                *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

/**********************************************************************/

/* Paso: SystemMatrixPattern: Paso_SystemMatrixPattern_mis 

   searches for a maximal independent set MIS in the matrix pattern 
   vertices in the maximal independent set are marked in mis_marker
   nodes to be considered are marked by -1 on the input in mis_marker

*/
/**********************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005              */
/* Author: gross@access.edu.au                                */

/**************************************************************/

#include "Paso.h"
#include "PasoUtil.h"
#include "SystemMatrixPattern.h"


/* used to generate pseudo random numbers: */

static double Paso_SystemMatrixPattern_mis_seed=.4142135623730951;


/***************************************************************/
 
#define IS_AVAILABLE -1
#define IS_IN_MIS_NOW -2
#define IS_IN_MIS -3
#define IS_CONNECTED_TO_MIS -4

void Paso_SystemMatrixPattern_mis(Paso_SystemMatrixPattern* pattern_p, index_t* mis_marker) {

  index_t index_offset=(pattern_p->type & PATTERN_FORMAT_OFFSET1 ? 1:0);
  dim_t i;
  index_t naib,iptr;
  bool_t flag;
  dim_t n=pattern_p->n_ptr;
  if (pattern_p->type & PATTERN_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"Paso_SystemMatrixPattern_mis: symmetric matrix pattern is not supported yet");
    return;
  }
  double *value=TMPMEMALLOC(n,double);
  if (!Paso_checkPtr(value)) {

   
     /* is there any vertex available ?*/
     while (Paso_Util_isAny(n,mis_marker,IS_AVAILABLE)) {
        /* step 1: assign a random number in [0,1] to each vertex */
        /* step 2: is the vertex is available, check if its value is the smaller than all values of its naigbours 
                   if the answer is yes, the vertex is put into the independend set and all 
                   its naighbours are removed from the graph by setting it mis_marker to FALSE */
      
        /* step 2: is the vertex is available, check if its value is the smaller than all values of its naigbours */

           /* assign random number in [0,1] to each vertex */
           #pragma omp parallel for private(i) schedule(static)
           for (i=0;i<n;++i) {
                 if (mis_marker[i]==IS_AVAILABLE) {
                    value[i]=fmod(Paso_SystemMatrixPattern_mis_seed*(i+1),1.);
                 } else {
                    value[i]=2.;
                 }
           }
           /* update the seed */
           /* Paso_SystemMatrixPattern_mis_seed=fmod(sqrt(Paso_SystemMatrixPattern_mis_seed*(n+1)),1.); */
           /* detect independent vertices as those vertices that have a value less than all values of its naigbours */
           #pragma omp parallel for private(naib,i,iptr,flag) schedule(static) 
           for (i=0;i<n;++i) {
              if (mis_marker[i]==IS_AVAILABLE) {
                 flag=IS_IN_MIS_NOW;
                 for (iptr=pattern_p->ptr[i]-index_offset;iptr<pattern_p->ptr[i+1]-index_offset; ++iptr) {
                     naib=pattern_p->index[iptr]-index_offset;
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
                 for (iptr=pattern_p->ptr[i]-index_offset;iptr<pattern_p->ptr[i+1]-index_offset; ++iptr) {
                     naib=pattern_p->index[iptr]-index_offset;
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
