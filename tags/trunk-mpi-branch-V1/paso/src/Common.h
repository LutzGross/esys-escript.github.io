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

#ifndef INC_PASO_COMMON
#define INC_PASO_COMMON

/**************************************************************/

/*    Finley finite element solver: common include file       */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003 */
/*   Version: $Id$ */

/**************************************************************/

/* some system values */
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
  //#define MEMALLOC(_LENGTH_,_TYPE_) (_TYPE_*) malloc(((size_t)(_LENGTH_))*sizeof(_TYPE_)); if( !(_LENGTH_) ) printf( "Trying malloc(0) at %s %d\n", __FILE__, __LINE__ ); //else printf( "malloc(%d) at %s %d\n", _LENGTH_, __FILE__, __LINE__ );
  #define MEMALLOC(_LENGTH_,_TYPE_) (_TYPE_*) malloc(((size_t)(_LENGTH_))*sizeof(_TYPE_))
  #define MEMFREE(_PTR_) if ((void *)(_PTR_) != NULL ) { free(_PTR_); (_PTR_) = NULL; }
  #ifdef _OPENMP
      #define THREAD_MEMALLOC(_LENGTH_,_TYPE_) (_TYPE_*) kmp_malloc(((size_t)(_LENGTH_))*sizeof(_TYPE_))
      #define THREAD_MEMFREE(_PTR_) if ((void *)(_PTR_) != NULL ) { kmp_free(_PTR_); (_PTR_) = NULL; }
  #else
     #define THREAD_MEMALLOC(_LENGTH_,_TYPE_) TMPMEMALLOC(_LENGTH_,_TYPE_)
     #define THREAD_MEMFREE(_PTR_) TMPMEMFREE(_PTR_)
  #endif
  #define MEMREALLOC(_POINTER_,_LENGTH_,_TYPE_) if( (_POINTER_)!=NULL ){ _POINTER_ = (_TYPE_*)realloc((void*)(_POINTER_),((size_t)(_LENGTH_))*sizeof(_TYPE_) );  }else{ _POINTER_ = (_TYPE_*)malloc( ((size_t)(_LENGTH_))*sizeof(_TYPE_) ); }

#else
  #define MEMALLOC(_LENGTH_,_TYPE_) (_TYPE_*) malloc(((size_t)(_LENGTH_))*sizeof(_TYPE_))
  #define MEMFREE(_PTR_) if ((void *)(_PTR_) != NULL ) { free(_PTR_); (_PTR_) = NULL; }
  #define THREAD_MEMALLOC(_LENGTH_,_TYPE_) TMPMEMALLOC(_LENGTH_,_TYPE_)
  #define THREAD_MEMFREE(_PTR_) TMPMEMFREE(_PTR_)
  #define MEMREALLOC(_POINTER_,_LENGTH_,_TYPE_) if( (_POINTER_)!=NULL ){ _POINTER_ = (_TYPE_*)realloc((void*)(_POINTER_),((size_t)(_LENGTH_))*sizeof(_TYPE_) );  }else{ _POINTER_ = (_TYPE_*)malloc( ((size_t)(_LENGTH_))*sizeof(_TYPE_) ); }
#endif

#endif /* #ifndef INC_PASO_COMMON */
