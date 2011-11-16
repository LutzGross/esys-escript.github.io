
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include <ripley/RipleyError.h>
#include <ripley/RipleyException.h>
extern "C" {
#include <esysUtils/error.h>
}

namespace ripley {

void checkPasoError() 
{
    if (!Esys_noError()) {
        // reset the error code to no error otherwise the next call to
        // this function may resurrect a previous error
        Esys_resetError();
        throw RipleyException(Esys_getErrorMessage());
    }
}

}  // end of namespace ripley

