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

#if !defined  escript_Utils_H
#define escript_Utils_H

namespace escript {

  /**
     \brief
     some functions

  */

  /**
     \brief
     set the number of threads 
  */
  void setNumberOfThreads(const int num_threads);

  /**
     \brief
     returns  the number of threads 
  */
  int getNumberOfThreads();

} // end of namespace
#endif
