
/*****************************************************************************
*
* Copyright (c) 2014-2016 by The University of Queensland
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

// Adapted from FunctionSpaceException.h

#ifndef  __ESCRIPT_SPLITWORLDEXCEPTION_H__
#define  __ESCRIPT_SPLITWORLDEXCEPTION_H__

#include "EsysException.h"

namespace escript
{

class SplitWorldException : public EsysException
{
public:
    SplitWorldException(const std::string& str) : EsysException(str) {}
    virtual ~SplitWorldException() throw() {}
};

} // end of namespace

#endif // __ESCRIPT_SPLITWORLDEXCEPTION_H__

