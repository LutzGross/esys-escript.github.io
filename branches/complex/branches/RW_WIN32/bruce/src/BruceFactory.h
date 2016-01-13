/* $Id$ */
/* 
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2005 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

#if !defined bruce_BruceFactory_20050901_H
#define bruce_BruceFactory_20050901_H
#ifdef MSVC
#ifdef BRUCE_EXPORTS
#define BRUCE_DLL __declspec(dllexport)
#else
#define BRUCE_DLL __declspec(dllimport)
#endif
#else
#define BRUCE_DLL
#endif
#include "escript/Data/AbstractContinuousDomain.h"
#include "Bruce/Bruce.h"

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
  BRUCE_DLL escript::AbstractContinuousDomain* brick(int n0=2, int n1=2, int n2=2,
		                           double l0=1.0, double l1=1.0, double l2=1.0);

  /**
     \brief
     Creates a rectangular mesh with n0 x n1 elements over the rectangle
     [l0,0] x [0,l1].

     \param n0,n1 Input - number of elements in each dimension
     \param l0,l1 Input - length of each side of rectangle
  */
  BRUCE_DLL escript::AbstractContinuousDomain* rectangle(int n0=2, int n1=2,
				               double l0=1.0, double l1=1.0);
 
} // end of namespace

#endif
