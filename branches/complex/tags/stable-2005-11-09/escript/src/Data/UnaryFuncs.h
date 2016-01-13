// $Id$
/* 
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/
                                                                           
#if !defined escript_UnaryFuncs_20041124_H
#define escript_UnaryFuncs_20041124_H

namespace escript {

inline
double
fsign(double x)
{
  if (x == 0) {
    return 0;
  } else {
    return x/fabs(x);
  }
}

} // end of namespace
#endif
