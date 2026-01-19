
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#if !defined  escript_TransportProblemException_20040608_H
#define escript_TransportProblemException_20040608_H

#include "system_dep.h"
#include "EsysException.h"

namespace escript
{

class TransportProblemException : public EsysException
{
public:
    TransportProblemException(const std::string& str) : EsysException(str) {}
};

} // end of namespace

#endif
