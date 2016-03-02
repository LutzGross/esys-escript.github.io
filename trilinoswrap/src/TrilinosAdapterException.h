
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

#ifndef __TRILINOSADAPTEREXCEPTION_H__
#define __TRILINOSADAPTEREXCEPTION_H__

#include <escript/EsysException.h>

namespace esys_trilinos {

class TrilinosAdapterException : public escript::EsysException
{
public:
    TrilinosAdapterException(const std::string& str)
        : escript::EsysException(str) {}
};


} // namespace esys_trilinos

#endif // __TRILINOSADAPTEREXCEPTION_H__

