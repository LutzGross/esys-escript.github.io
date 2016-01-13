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
                                                                           
#if !defined  esysUtils_esysExceptionTranslator_20040419_H
#define esysUtils_esysExceptionTranslator_20040419_H

#include "EsysException.h"
#include "boost/python/errors.hpp"

namespace esysUtils {
  /**
     \brief
     Function which translates an EsysException into a python exception
  */
	void esysExceptionTranslator(EsysException const& e);
} // end of namespace
#endif
