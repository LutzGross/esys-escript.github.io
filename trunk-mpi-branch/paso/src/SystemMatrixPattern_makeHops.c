/* $Id:*/

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

/* Paso: SystemMatrixPatternPattern_makeHops : sets the hop vector to hand over to Pattern object` */

/**************************************************************/
 
/* Copyrights by ACcESS Australia 2003, 2004,2005, 2006, 2007 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrixPattern.h"

/**************************************************************/

/* allocates a SystemMatrixPattern  */

void Paso_SystemMatrixPattern_makeHops(int type,
                                       Paso_Distribution* distribution,
                                       index_t* ptr,
                                       index_t* index,
                                       dim_t* numHops,
                                       index_t **hop) {
  dim_t i, p, j, d, *mask=NULL, *mask_g=NULL, r;
  register dim_t itmp;
  Paso_MPIInfo* mpi_info=distribution->mpi_info;
  dim_t s=mpi_info->size;
  dim_t myProc=mpi_info->rank;
  dim_t myN=distribution->myNumComponents;
  *numHops=0;
  *hop=MEMALLOC(s,index_t);
  if (Paso_checkPtr(*hop)) return;
  if (s < 2) {
     *numHops=1;
     *hop[0]=1;
  } else {
     #ifdef PASO_MPI
        mask=TMPMEMALLOC(s,index_t);
        mask_g=TMPMEMALLOC(s,index_t);
        if (Paso_checkPtr(mask) || Paso_checkPtr(mask_g)) {
           MEMFREE(*hop);
        } else {
           for (i=0;i<s;++i) {
               mask[i]=0;
               (*hop)[i]=0;
           }
           mask[0]=1;
           if (type & PATTERN_FORMAT_OFFSET1) {
              #pragma omp parallel private(p,i,j,itmp,r)
              {
                 for (p=0; p<s; ++p) {
                   r=0;
                   #pragma omp for schedule(static)
                   for (i=0;i<myN;++i) {
                       for (j=ptr[i]; j<ptr[i+1]; ++j) {
                           itmp=index[j];
                           if ( (distribution->first_component[p]<itmp) && 
                                (itmp<=distribution->first_component[p+1]) ) {
                               r=1;
                               break;
                           }
                       }
                    }
                    #pragma omp critical
                    {
                       mask[PASO_MPI_mod(p-myProc,s)]=MAX(mask[PASO_MPI_mod(p-myProc,s)],r);
                    }
                 }
              }
          } else {
              #pragma omp parallel private(p,i,j,itmp,r)
              {
                 for (p=0; p<s; ++p) {
                   r=0;
                   #pragma omp for schedule(static)
                   for (i=0;i<myN;++i) {
                       for (j=ptr[i]; j<ptr[i+1]; ++j) {
                           itmp=index[j];
                           if ( (distribution->first_component[p]<=itmp) && 
                                (itmp<distribution->first_component[p+1]) ) {
                               r=1;
                               break;
                           }
                       }
                    }
                    #pragma omp critical
                    {
                       mask[PASO_MPI_mod(p-myProc,s)]=MAX(mask[PASO_MPI_mod(p-myProc,s)],r);
                    }
                 }
              }
          }
          /* make global*/
          MPI_Allreduce(mask, mask_g, s, MPI_INT, MPI_MAX, mpi_info->comm);
          /* get the hops now */
          *numHops=1;
          (*hop)[0]=1;
          for (d=1; d<s; ++d) {
             if (mask_g[d]==1) {
                (*hop)[*numHops]=1;
                (*numHops)++;
             } else {
                (*hop)[*numHops-1]++;
             }
          }
{ int q;
printf("numHops: ");
for (q=0;q<*numHops;++q) printf(" %d",(*hop)[q]);
printf("\n");
}
        }
        TMPMEMFREE(mask);
        TMPMEMFREE(mask_g);
     #else
        *numHops=1;
        (*hop)[0]=1;
     #endif
   }
}
