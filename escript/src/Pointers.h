
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#ifndef POINTERS_H_2008
#define POINTERS_H_2008

/** \file Pointers.h 
\brief Typedefs and macros for reference counted storage.
*/

// The idea is that we should be able to easily switch between shared_ptr and intrusive_ptr if required

// Where to find the base class which supplies refcounting
#define REFCOUNT_BASE_FILE <boost/enable_shared_from_this.hpp>
// The name of the class to extend
#define REFCOUNT_BASE_CLASS(x) boost::enable_shared_from_this<x>


#define POINTER_WRAPPER_CLASS(x) boost::shared_ptr<x>


#include REFCOUNT_BASE_FILE


#endif
