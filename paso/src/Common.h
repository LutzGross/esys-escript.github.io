
/* $Id$ */

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

#ifndef INC_PASO_COMMON
#define INC_PASO_COMMON

/**************************************************************/

/*    Finley finite element solver: common include file       */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003 */
/*   Version: $Id$ */

/**************************************************************/

/* some system values */
/* FIXME: This is not satisfactory. _ECC, __INTEL_COMPILER, and other               */
/*        intel compiler pre-defines need to be handled (__ICL, __ICC come to mind) */
#if ( defined __INTEL_COMPILER )
#include <mathimf.h>
#else
#include <math.h>
#endif


#include <float.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#define LenString_MAX FILENAME_MAX*2
#define LenErrorMsg_MAX LenString_MAX
 
/* on some arcitectures it could be a good idea to use long rather than int */
/* this has not really been tested */

typedef int dim_t;
typedef int index_t;
typedef int bool_t;
typedef int type_t;
typedef int err_t;

#define INDEX_T_MAX INT_MAX
#define EPSILON DBL_EPSILON
#define LARGE_POSITIVE_FLOAT DBL_MAX

/**************************************************************/

/*   some useful functions: */

#define FALSE 0
#define TRUE 1
#define UNKNOWN -1
#define DBLE(_x_) (double)(_x_)
#define INDEX1(_X1_) (_X1_)
#define INDEX2(_X1_,_X2_,_N1_) ((_X1_)+(_N1_)*(_X2_))
#define INDEX3(_X1_,_X2_,_X3_,_N1_,_N2_) ((_X1_)+(_N1_)*INDEX2(_X2_,_X3_,_N2_))
#define INDEX4(_X1_,_X2_,_X3_,_X4_,_N1_,_N2_,_N3_) ((_X1_)+(_N1_)*INDEX3(_X2_,_X3_,_X4_,_N2_,_N3_))
#define INDEX5(_X1_,_X2_,_X3_,_X4_,_X5_,_N1_,_N2_,_N3_,_N4_) ((_X1_)+(_N1_)*INDEX4(_X2_,_X3_,_X4_,_X5_,_N2_,_N3_,_N4_))
#define INDEX6(_X1_,_X2_,_X3_,_X4_,_X5_,_X6_,_N1_,_N2_,_N3_,_N4_,_N5_) ((_X1_)+(_N1_)*INDEX5(_X2_,_X3_,_X4_,_X5_,_X6_,_N2_,_N3_,_N4_,_N5_))

#define MAX(_arg1_,_arg2_) ((_arg1_)>(_arg2_) ?  (_arg1_) : (_arg2_))
#define MAX3(_arg1_,_arg2_,_arg3_) MAX(_arg1_,MAX(_arg2_,_arg3_))
#define MIN(_arg1_,_arg2_) ((_arg1_)>(_arg2_) ?  (_arg2_) : (_arg1_)) 
#define MIN3(_arg1_,_arg2_,_arg3_) MIN(_arg1_,MIN(_arg2_,_arg3_))
#define ABS(_arg_) MAX((_arg_),-(_arg_))
#define SIGN(_arg_) ((_arg_)>0 ?  1  : ((_arg_)<0 ? -1 : 0 ))
#define SWAP(_a0_,_a1_,_type_) { \
                                _type_ s; \
                                s=(_a0_); \
                                _a0_= (_a1_); \
                                _a1_=s; \
                               }
/**************************************************************/
/*    memory allocation:                                      */
/*    Wise to not use PASO_MALLOC/FREE/REALLOC and            */
/*    PASO_THREAD... directly. These are only for tailoring   */
/*    the main macros that follow                             */
/**************************************************************/

#if defined(_WIN32) // Use python for memory management on windows.

  #include <python.h>

  #define PASO_MALLOC PyMem_Malloc
  #define PASO_FREE PyMem_Free
  #define PASO_REALLOC PyMem_Realloc

#else

  #include <malloc.h>

  #define PASO_MALLOC malloc
  #define PASO_FREE free
  #define PASO_REALLOC realloc

#endif


/* FIXME: This is not satisfactory. _ECC, __INTEL_COMPILER, and other               */
/*        intel compiler pre-defines need to be handled (__ICL, __ICC come to mind) */
/*        Also, _WIN32 may take this branch one day...                              */
/*        SO KEEP ALL THREAD_MEMALLOC/FREEs CONFINED TO THE PASO LIBRARY.           */

#if defined(__ECC) && defined(_OPENMP) // ECC version of intel compiler with openmp.
  #include <omp.h>
  #define PASO_THREAD_MALLOC kmp_malloc
  #define PASO_THREAD_FREE kmp_free
#else
  #define PASO_THREAD_MALLOC PASO_MALLOC
  #define PASO_THREAD_FREE PASO_FREE
#endif

/******************The main macros ************************************/ 

#define MEMALLOC(_LENGTH_,_TYPE_)                                     \
  (_TYPE_*) PASO_MALLOC(((size_t)(_LENGTH_))*sizeof(_TYPE_))

/* do {} while(0) -  an old trick for bracketing a macro that */
/* makes sure a semi-colon does no harm.                      */

#define MEMFREE(_PTR_)                                                  \
do                                                                      \
{                                                                       \
  if ((void *)(_PTR_) != NULL ) { PASO_FREE(_PTR_); (_PTR_) = NULL; }   \
} while(0)

#define MEMREALLOC(_POINTER_,_LENGTH_,_TYPE_)                             \
do                                                                        \
{                                                                         \
   if( (_POINTER_)!=NULL )                                                \
   {                                                                      \
      _POINTER_ = (_TYPE_*)PASO_REALLOC((void*)(_POINTER_),               \
                                   ((size_t)(_LENGTH_))*sizeof(_TYPE_) ); \
   }                                                                      \
   else                                                                   \
   {                                                                      \
      _POINTER_ = (_TYPE_*)PASO_MALLOC( ((size_t)(_LENGTH_))*sizeof(_TYPE_) ); \
   }                                                                      \
} while(0)

#define TMPMEMALLOC MEMALLOC
#define TMPMEMFREE MEMFREE

#define THREAD_MEMALLOC(_LENGTH_,_TYPE_)                          \
   (_TYPE_*) PASO_THREAD_MALLOC/**/(((size_t)(_LENGTH_))*sizeof(_TYPE_))

#define THREAD_MEMFREE(_PTR_)                                                \
do                                                                           \
{                                                                            \
  if ((void *)(_PTR_) != NULL ) { PASO_THREAD_FREE(_PTR_); (_PTR_) = NULL; } \
} while(0)


/*
	This was useful for seeing what memory is being allocated if you get an "Out of memory" error
        if used again, bracket with do {} while(0)
#define TMPMEMALLOC(_LENGTH_,_TYPE_) (printf("TMPMEMALLOC at %s %d #bytes=%d*%d=%d\n", __FILE__, __LINE__, _LENGTH_, sizeof(_TYPE_), ((size_t)(_LENGTH_) * sizeof(_TYPE_))) , (_TYPE_*) malloc(((size_t)(_LENGTH_))*sizeof(_TYPE_)))

	This was useful for seeing where memory is being freed if you get an "glibc detected...free" error
        if used again, bracket with do {} while(0)
#define TMPMEMFREE(_PTR_) if ((void *)(_PTR_) != NULL ) { printf("TMPMEMFREE AAA %s %d\n", __FILE__, __LINE__); free(_PTR_); (_PTR_) = NULL; printf("TMPMEMFREE BBB %s %d\n", __FILE__, __LINE__); }
*/

#endif /* #ifndef INC_PASO_COMMON */
