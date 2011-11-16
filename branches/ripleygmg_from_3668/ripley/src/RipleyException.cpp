
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

#include <ripley/RipleyException.h>

namespace ripley {

const std::string RipleyException::exceptionNameValue("RipleyException");

const std::string& RipleyException::exceptionName() const
{
    return exceptionNameValue;
}

} // namespace ripley

