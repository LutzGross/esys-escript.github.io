/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

#include "Utils.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace escript {

void setNumberOfThreads(const int num_threads) 
{

   #ifdef _OPENMP
   omp_set_num_threads(num_threads);
   #endif

}

int getNumberOfThreads() 
{
   #ifdef _OPENMP
   return omp_get_max_threads();
   #else
   return 1;
   #endif

}

}  // end of namespace
