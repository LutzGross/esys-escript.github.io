
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


#if !defined  escript_FunctionSpaceException_20040602_H
#define escript_FunctionSpaceException_20040602_H
#include "system_dep.h"

#include "EsysException.h"

namespace escript
{

class FunctionSpaceException : public EsysException
{
public:
    FunctionSpaceException(const std::string& str) : EsysException(str) {}
    virtual ~FunctionSpaceException() throw() {}
};

} // end of namespace

#endif

