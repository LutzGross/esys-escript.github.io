
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/**
\file esys_malloc.h
\ingroup Other
 */
//
// @(#) esys_malloc.h
//

#ifndef esys_malloc_h
#define esys_malloc_h

#ifdef _WIN32

#   include <python.h>

#   define ESYS_MALLOC PyMem_Malloc
#   define ESYS_FREE PyMem_Free
#   define ESYS_REALLOC PyMem_Realloc

#else

#   include <stdlib.h>

#   define ESYS_MALLOC ::malloc
#   define ESYS_FREE ::free
#   define ESYS_REALLOC ::realloc

#endif

namespace esysUtils
{

   inline
   void *malloc(size_t len)
   {
      return ESYS_MALLOC(len);
   }

   inline
   void free(void *ptr)
   {
      ESYS_FREE(ptr);
      return;
   }

   inline
   void *realloc(void *ptr, size_t len)
   {
      return ESYS_REALLOC(ptr,len);
   }
}

#undef ESYS_MALLOC
#undef ESYS_FREE
#undef ESYS_REALLOC

#endif
