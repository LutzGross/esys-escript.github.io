
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


#if !defined  escript_SystemMatrixException_20040608_H
#define escript_SystemMatrixException_20040608_H

#include "system_dep.h"
#include "EsysException.h"

namespace escript
{
  /**
  \brief
  SystemMatrixException exception class.

  Description:
  SystemMatrixException exception class.
  The class provides a public function returning the exception name
  */
class SystemMatrixException : public EsysException
{
public:
    SystemMatrixException(const std::string& str) : EsysException(str) {}

    /// Destructor
    virtual ~SystemMatrixException() throw() {}
};

} // end of namespace

#endif

