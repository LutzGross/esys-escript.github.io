
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


#ifndef INC_ESYS_MEM
#define INC_ESYS_MEM

/************************************************************************************/
/*   Macros to deal with memory management */
/********************************************/


/************************************************************************************/
/*    memory allocation:                                      */
/*    Wise to not use PASO_MALLOC/FREE/REALLOC and            */
/*    PASO_THREAD... directly. These are only for tailoring   */
/*    the main macros that follow                             */
/************************************************************************************/


/*#if defined(_WIN32) */ /* Use python for memory management on windows. */
/*
  #include <python.h>

  #define PASO_MALLOC PyMem_Malloc
  #define PASO_FREE PyMem_Free
  #define PASO_REALLOC PyMem_Realloc

#else
*/
  #include <stdlib.h>

  #define PASO_MALLOC malloc
  #define PASO_FREE free
  #define PASO_REALLOC realloc

/*#endif */

/* FIXME: This is not satisfactory.                                */
/* _ECC, __INTEL_COMPILER, and other                               */
/* intel compiler pre-defines need to be handled                   */
/* (__ICL, __ICC come to mind)                                     */
/* Also, _WIN32 may take this branch one day...                    */
/* SO KEEP ALL THREAD_MEMALLOC/FREEs CONFINED TO THE PASO LIBRARY. */

#if defined(__ECC) && defined(_OPENMP) /* ECC version of intel compiler with openmp. */
  #include <omp.h>
  #define PASO_THREAD_MALLOC kmp_malloc
  #define PASO_THREAD_FREE kmp_free
#else
  #define PASO_THREAD_MALLOC PASO_MALLOC
  #define PASO_THREAD_FREE PASO_FREE
#endif


/* Prepare for the day that this becomes sharable. */
/* and we wish to do multi-file optimisations on windows */

#define PASO_DLL_API

#ifdef _WIN32
#   ifndef PASO_STATIC_LIB
#      undef PASO_DLL_API
#      ifdef PASO_EXPORTS
#         define PASO_DLL_API __declspec(dllexport)
#      else
#         define PASO_DLL_API __declspec(dllimport)
#      endif
#   endif
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

#define MEMREALLOC(_RETP_,_POINTER_,_LENGTH_,_TYPE_)                    \
do                                                                        \
{                                                                         \
   if( (_POINTER_)!=NULL )                                                \
   {                                                                      \
      _RETP_ = (_TYPE_*)PASO_REALLOC((void*)(_POINTER_),               \
                                   ((size_t)(_LENGTH_))*sizeof(_TYPE_) ); \
   }                                                                      \
   else                                                                   \
   {                                                                      \
      _RETP_ = (_TYPE_*)PASO_MALLOC( ((size_t)(_LENGTH_))*sizeof(_TYPE_) ); \
   }                                                                      \
} while(0)

#define TMPMEMALLOC MEMALLOC
#define TMPMEMFREE MEMFREE
#define TMPMEMREALLOC MEMREALLOC

#define THREAD_MEMALLOC(_LENGTH_,_TYPE_)                          \
   (_TYPE_*) PASO_THREAD_MALLOC/**/(((size_t)(_LENGTH_))*sizeof(_TYPE_))

#define THREAD_MEMFREE(_PTR_)                                                \
do                                                                           \
{                                                                            \
  if ((void *)(_PTR_) != NULL ) { PASO_THREAD_FREE(_PTR_); (_PTR_) = NULL; } \
} while(0)


#endif 
