
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

#if !defined bruce_BruceFactory_20050901_H
#define bruce_BruceFactory_20050901_H
#include "system_dep.h"
#include "escript/AbstractContinuousDomain.h"

namespace bruce {

  /**
     \brief
     A suite of factory methods for creating various Bruces.

     Description:
     A suite of factory methods for creating various Bruces.
  */
 
  /**
     \brief
     Creates a rectangular mesh with n0 x n1 x n2 elements over the brick 
     [l0,0,0] x [0,l1,0] x [0,0,l2].

     \param n0,n1,n2 Input - number of elements in each dimension
     \param l0,l1,l2 Input - length of each side of brick
  */
  BRUCE_DLL_API escript::AbstractContinuousDomain* brick(int n0=2, int n1=2, int n2=2,
		                           double l0=1.0, double l1=1.0, double l2=1.0);

  /**
     \brief
     Creates a rectangular mesh with n0 x n1 elements over the rectangle
     [l0,0] x [0,l1].

     \param n0,n1 Input - number of elements in each dimension
     \param l0,l1 Input - length of each side of rectangle
  */
  BRUCE_DLL_API escript::AbstractContinuousDomain* rectangle(int n0=2, int n1=2,
				               double l0=1.0, double l1=1.0);
 
} // end of namespace

#endif
