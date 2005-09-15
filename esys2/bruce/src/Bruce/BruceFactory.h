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

#include "escript/Data/AbstractContinuousDomain.h"
#include "bruce/Bruce/Bruce.h"

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
     [0,l0] x [0,l1] x [0,l2].

     \param n0,n1,n2 Input - number of elements in each dimension
     \param l0,l1,l2 Input - length of each side of brick
  */
  escript::AbstractContinuousDomain* brick(int n0=1,int n1=1,int n2=1,
		                           double l0=1.0,double l1=1.0,double l2=1.0);

  /**
     \brief
     Creates a rectangular mesh with n0 x n1 elements over the rectangle
     [0,l0] x [0,l1].

     \param n0,n1 Input - number of elements in each dimension
     \param l0,l1 Input - length of each side of brick
  */
  escript::AbstractContinuousDomain* rectangle(int n0=1,int n1=1,
				      double l0=1.0, double l1=1.0);
 
} // end of namespace

#endif
