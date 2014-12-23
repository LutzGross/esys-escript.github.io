
/*****************************************************************************
*
* Copyright (c) 2010-2014 by University of Queensland
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


#ifndef __ESYS_TYPES_H__
#define __ESYS_TYPES_H__
 
#ifdef ESYS_INDEXTYPE_LONG
typedef long index_t;
#else
typedef int index_t;
#endif

typedef index_t dim_t;
typedef int type_t;
typedef int err_t;

#endif // __ESYS_TYPES_H__

