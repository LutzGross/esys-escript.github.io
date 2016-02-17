
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

#include "PasoException.h"
#include <esysUtils/error.h>

namespace paso
{

const std::string PasoException::exceptionNameValue("PasoException");

const std::string& PasoException::exceptionName() const
{
  return exceptionNameValue;
}

void checkPasoError() 
{
  if (Esys_noError()) {
    return;
  } else {
    //
    // reset the error code to no error otherwise the next call to
    // this function may resurrect a previous error
    Esys_resetError();
    throw PasoException(Esys_getErrorMessage());
  }
}

} // namespace paso

