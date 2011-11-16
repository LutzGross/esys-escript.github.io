
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

#ifndef __RIPLEY_UTIL_H__
#define __RIPLEY_UTIL_H__

#include <ripley/Ripley.h>

struct Esys_MPIInfo;

namespace ripley {

/**
   \brief
   Returns the minimum and maximum value in an IndexVector as a pair
   disregarding the value 'ignore'.
*/
IndexPair getFlaggedMinMax(const IndexVector &values, index_t ignore);

/**
   \brief
   Returns the minimum and maximum value in an IndexVector over all processes
   as a pair.
*/
IndexPair getGlobalMinMax(const IndexVector &values, Esys_MPIInfo *mpiInfo);

/**
   \brief
   Returns the minimum and maximum value in an IndexVector as a pair.
*/
IndexPair getMinMax(const IndexVector &values);

/**
   \brief
   Returns all unique indices contained in an IndexVector over all processes.
*/
IndexVector getUniqueValues(const IndexVector &values, Esys_MPIInfo *mpiinfo);

/**
   \brief
   Returns an IndexVector that contains all values >=0 contained within mask.
*/
IndexVector packMask(const IndexVector &mask);

} // end of namespace ripley

#endif // __RIPLEY_UTIL_H__

