
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#ifndef POINTERS_H_2008
#define POINTERS_H_2008

/** \file Pointers.h 
  \brief Typedefs and macros for reference counted storage.
*/

// The idea is that we should be able to easily switch between shared_ptr
// and intrusive_ptr if required

// Where to find the base class which supplies refcounting
#define REFCOUNT_BASE_FILE <boost/enable_shared_from_this.hpp>
// The name of the class to extend
#define REFCOUNT_BASE_CLASS(x) boost::enable_shared_from_this<x>

#define POINTER_WRAPPER_CLASS(x) boost::shared_ptr<x>

#define REFCOUNTNS	boost

#include REFCOUNT_BASE_FILE

#endif

