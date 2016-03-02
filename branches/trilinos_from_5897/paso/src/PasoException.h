
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

#ifndef __PASO_EXCEPTION_H__
#define __PASO_EXCEPTION_H__

#include <escript/EsysException.h>

namespace paso {

/**
  \brief
  PasoException exception class.

  Description:
  PasoException exception class.
  The class provides a public function returning the exception name
*/
class PasoException : public escript::EsysException
{
public:
    PasoException(const std::string& str) : EsysException(str) {}
    virtual ~PasoException() throw() {}
};

void checkPasoError(); 


} // end of namespace

#endif // __PASO_EXCEPTION_H__

