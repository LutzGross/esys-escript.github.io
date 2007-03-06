/* $Id$ */

/*
********************************************************************************
*               Copyright  2006 by ACcESS MNRF                                 *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

/**************************************************************/

/*   Paso: system matrix pattern                            */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2004,2005 */
/*   Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_SYSTEMMATRIXPATTERN
#define INC_PASO_SYSTEMMATRIXPATTERN

#include "Distribution.h"
#include "Common.h"
#include "Paso_MPI.h"

/**************************************************************/

#define PATTERN_FORMAT_DEFAULT 0
#define PATTERN_FORMAT_SYM 1
#define PATTERN_FORMAT_OFFSET1 2

typedef struct Paso_SystemMatrixPattern {
  int type;
  dim_t myNumOutput;
  index_t numOutput;
  dim_t myNumInput;
  index_t numInput;
  index_t* ptr;
  index_t* index;
  dim_t myLen;
  dim_t reference_counter;
  Paso_MPIInfo *mpi_info;
  Paso_Distribution *output_distribution; 
  Paso_Distribution *input_distribution; 
} Paso_SystemMatrixPattern;


/*  interfaces: */

Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_alloc(int type, Paso_Distribution* ptr_distribution, Paso_Distribution* index_range_distribution, index_t* ptr, index_t* index);
Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_reference(Paso_SystemMatrixPattern*);
void Paso_SystemMatrixPattern_dealloc(Paso_SystemMatrixPattern*);
int Paso_comparIndex(const void *,const void *);
Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_getSubpattern(Paso_SystemMatrixPattern*,dim_t,dim_t,index_t*,index_t*);
void Paso_SystemMatrixPattern_mis(Paso_SystemMatrixPattern* pattern_p, index_t* mis_marker);
Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_unrollBlocks(Paso_SystemMatrixPattern*,int, dim_t,dim_t);

#endif /* #ifndef INC_PASO_SYSTEMPATTERN */
