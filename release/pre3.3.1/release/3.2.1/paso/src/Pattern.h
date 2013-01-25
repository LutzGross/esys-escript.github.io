
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   Paso: pattern                                            */

/**************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#ifndef INC_PASO_PATTERN
#define INC_PASO_PATTERN

#include "Paso.h"
#include "IndexList.h"

/**************************************************************/

typedef struct Paso_Pattern {
  int type;
  dim_t numOutput;	/* Number of rows the ptr array [CSR] for CSC it's the number of cols*/
  dim_t numInput;	/* Number of cols [CSR] */
  dim_t len;		/* number of non-zeros */
  index_t* ptr;		/* ptr[n] to ptr[n+1] lists indicies (in index) of non-zeros in row n*/
  index_t* index;	/* Non-major indicies of non-zeros (in CSR this will be col numbers) */ 
  index_t *main_iptr;  /* pointer to main diagonal entry */
  dim_t numColors;    /* number of colors */
  index_t* coloring;     /* coloring index: input with the same color are not connected */
  dim_t reference_counter;
} Paso_Pattern;

PASO_DLL_API
Paso_Pattern* Paso_Pattern_alloc(int type, dim_t numOutput, dim_t numInput, index_t* ptr, index_t* index);

PASO_DLL_API

PASO_DLL_API
Paso_Pattern* Paso_Pattern_getReference(Paso_Pattern*);

PASO_DLL_API
void Paso_Pattern_free(Paso_Pattern*);

PASO_DLL_API
int Paso_comparIndex(const void *,const void *);

PASO_DLL_API
Paso_Pattern* Paso_Pattern_unrollBlocks(Paso_Pattern*,int, dim_t,dim_t);

PASO_DLL_API
Paso_Pattern* Paso_Pattern_getSubpattern(Paso_Pattern*,dim_t,dim_t,index_t*,index_t*);

PASO_DLL_API
bool_t Paso_Pattern_isEmpty(Paso_Pattern* in);

PASO_DLL_API
void Paso_Pattern_mis(Paso_Pattern* pattern_p, index_t* mis_marker);

PASO_DLL_API
void Paso_Pattern_reduceBandwidth(Paso_Pattern* self,index_t* oldToNew);

PASO_DLL_API
void Paso_Pattern_color(Paso_Pattern* patter, index_t* num_colors, index_t* colorOf);
Paso_Pattern* Paso_Pattern_multiply(int type, Paso_Pattern* A, Paso_Pattern* B);

PASO_DLL_API
Paso_Pattern* Paso_Pattern_binop(int type, Paso_Pattern* A, Paso_Pattern* B);

PASO_DLL_API
index_t* Paso_Pattern_borrowMainDiagonalPointer(Paso_Pattern* A);

PASO_DLL_API
Paso_Pattern* Paso_Pattern_fromIndexListArray(dim_t n0, Paso_IndexListArray* index_list_array,index_t range_min,index_t range_max, index_t index_offset);

PASO_DLL_API
dim_t Paso_Pattern_getNumColors(Paso_Pattern* A);

PASO_DLL_API
index_t* Paso_Pattern_borrowColoringPointer(Paso_Pattern* A);

PASO_DLL_API
dim_t Paso_Pattern_maxDeg(Paso_Pattern* A);

#endif /* #ifndef INC_PASO_PATTERN */
