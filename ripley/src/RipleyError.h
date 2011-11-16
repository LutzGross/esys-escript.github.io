
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

#ifndef __RIPLEY_ERROR_H__
#define __RIPLEY_ERROR_H__

#include <ripley/system_dep.h>

namespace ripley {

/**
   \brief
   Converts a C paso error into a C++ exception.
*/
RIPLEY_DLL_API
void checkPasoError();

} // namespace ripley

#endif // __RIPLEY_ERROR_H__

