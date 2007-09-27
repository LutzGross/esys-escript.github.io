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

/**************************************************************/

/*   Some utility routines: */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004,2005 */

/**************************************************************/

#include "Common.h"
#include "PasoUtil.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/* returns true if array contains value */
bool_t Paso_Util_isAny(dim_t N,index_t* array,index_t value) {
   bool_t out=FALSE;
   dim_t i;
   #pragma omp parallel for private(i) schedule(static) reduction(||:out)
   for (i=0;i<N;i++) out = out || (array[i]==value);
   return out;
}


/* copies source to taget: */
void Paso_copyDouble(dim_t n,double* source, double* target) {
  dim_t i;
  for (i=0;i<n;i++) target[i]=source[i];
}

/**************************************************************/

/* calculates the cummultative sum in array and returns the total sum */

/**************************************************************/
index_t Paso_Util_cumsum(dim_t N,index_t* array) {
   index_t out=0,tmp;
   dim_t i;
   #ifdef _OPENMP
      index_t partial_sums[omp_get_max_threads()],sum;
      #pragma omp parallel private(sum,i,tmp)
      {
        sum=0;
        #pragma omp for schedule(static)
        for (i=0;i<N;++i) sum+=array[i];
        partial_sums[omp_get_thread_num()]=sum;
        #pragma omp barrier
        #pragma omp master
        {
          out=0;
          for (i=0;i<omp_get_max_threads();++i) {
             tmp=out;
             out+=partial_sums[i];
             partial_sums[i]=tmp;
           } 
        }
        #pragma omp barrier
        sum=partial_sums[omp_get_thread_num()];
        #pragma omp for schedule(static)
        for (i=0;i<N;++i) {
          tmp=sum;
          sum+=array[i];
          array[i]=tmp;
        } 
      }
   #else 
      for (i=0;i<N;++i) {
         tmp=out;
         out+=array[i];
         array[i]=tmp;
      }
   #endif
   return out;
}

/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:39  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.1  2005/09/05 06:29:48  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
