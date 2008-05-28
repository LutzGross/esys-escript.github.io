
/* $Id: Pattern.h 1306 2007-09-18 05:51:09Z ksteube $ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/*   Paso: pattern                                            */

/**************************************************************/

/*   Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_PATTERN
#define INC_PASO_PATTERN

#include "Common.h"

/**************************************************************/

#define PATTERN_FORMAT_DEFAULT 0
#define PATTERN_FORMAT_SYM 1
#define PATTERN_FORMAT_OFFSET1 2

typedef struct Paso_Pattern {
  int type;
  dim_t numOutput;
  dim_t numInput;
  dim_t input_block_size; /* logical block size in the input */
  dim_t output_block_size; /* logical block size in the output */
  dim_t block_size; /* = input_block_size * output_block_size */

  dim_t len;
  index_t* ptr;
  index_t* index;
  dim_t reference_counter;
} Paso_Pattern;

/*  interfaces: */

Paso_Pattern* Paso_Pattern_alloc(int type, dim_t input_block_size, dim_t output_block_size, dim_t numOutput, index_t* ptr, index_t* index);
Paso_Pattern* Paso_Pattern_getReference(Paso_Pattern*);
void Paso_Pattern_free(Paso_Pattern*);
int Paso_comparIndex(const void *,const void *);
Paso_Pattern* Paso_Pattern_unrollBlocks(Paso_Pattern*,int, dim_t,dim_t);
Paso_Pattern* Paso_Pattern_getSubpattern(Paso_Pattern*,dim_t,dim_t,index_t*,index_t*);
bool_t Paso_Pattern_isEmpty(Paso_Pattern* in);
void Paso_Pattern_mis(Paso_Pattern* pattern_p, index_t* mis_marker);
void Paso_Pattern_reduceBandwidth(Paso_Pattern* self,index_t* oldToNew);
void Paso_Pattern_color(Paso_Pattern* patter, index_t* num_colors, index_t* colorOf);

#endif /* #ifndef INC_PASO_SYSTEMPATTERN */
