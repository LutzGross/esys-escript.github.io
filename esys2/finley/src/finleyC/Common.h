/* $Id$ */

#ifndef INC_FINLEY_COMMON
#define INC_FINLEY_COMMON

/**************************************************************/

/*    Finley finite element solver: common include file       */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003 */
/*   Version: $Id$ */

/**************************************************************/

/* some system values */

#include <float.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>

#define LenString_MAX FILENAME_MAX*2
#define LenErrorMsg_MAX LenString_MAX
 
/* on some arcitectures it could be a good idea to use long rather than int */
/* this has not really been tested */

typedef int maybelong;
#define MAYBELONG_MAX INT_MAX
#define HUGE_VAL DBL_MAX
#define EPSILON DBL_EPSILON

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
#define MIN(_arg1_,_arg2_) ((_arg1_)>(_arg2_) ?  (_arg2_) : (_arg1_)) 
#define ABS(_arg_) MAX((_arg_),-(_arg_))
/**************************************************************/

/*    memory allocation:                                      */

#define TMPMEMALLOC(_LENGTH_,_TYPE_) (_TYPE_*) malloc(((size_t)(_LENGTH_))*sizeof(_TYPE_))
#define TMPMEMFREE(_PTR_) if ((void *)(_PTR_) != NULL ) { free(_PTR_); (_PTR_) = NULL; }
#ifdef __ECC
  #define MEMALLOC(_LENGTH_,_TYPE_) (_TYPE_*) malloc(((size_t)(_LENGTH_))*sizeof(_TYPE_))
  #define MEMFREE(_PTR_) if ((void *)(_PTR_) != NULL ) { free(_PTR_); (_PTR_) = NULL; }
  #ifdef _OPENMP
      #define THREAD_MEMALLOC(_LENGTH_,_TYPE_) (_TYPE_*) kmp_malloc(((size_t)(_LENGTH_))*sizeof(_TYPE_))
      #define THREAD_MEMFREE(_PTR_) if ((void *)(_PTR_) != NULL ) { kmp_free(_PTR_); (_PTR_) = NULL; }
  #else
     #define THREAD_MEMALLOC(_LENGTH_,_TYPE_) TMPMEMALLOC(_LENGTH_,_TYPE_)
     #define THREAD_MEMFREE(_PTR_) TMPMEMFREE(_PTR_)
  #endif
#else
  #define MEMALLOC(_LENGTH_,_TYPE_) (_TYPE_*) malloc(((size_t)(_LENGTH_))*sizeof(_TYPE_))
  #define MEMFREE(_PTR_) if ((void *)(_PTR_) != NULL ) { free(_PTR_); (_PTR_) = NULL; }
  #define THREAD_MEMALLOC(_LENGTH_,_TYPE_) TMPMEMALLOC(_LENGTH_,_TYPE_)
  #define THREAD_MEMFREE(_PTR_) TMPMEMFREE(_PTR_)
#endif

#endif /* #ifndef INC_FINLEY_COMMON */

/*
 * $Log$
 * Revision 1.2  2004/12/14 05:39:30  jgs
 * *** empty log message ***
 *
 * Revision 1.1.1.1.2.1  2004/11/24 01:37:13  gross
 * some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.3  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 * Revision 1.2  2004/07/01 23:49:17  gross
 * macro for copying of a double array added
 *
 * Revision 1.1.1.1  2004/06/24 04:00:40  johng
 * Initial version of eys using boost-python.
 *
 *
 */
