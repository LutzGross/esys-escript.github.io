/*****************************************************************************
*
* Copyright (c) 2003-2019 by The University of Queensland
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

#ifndef __OXLEY_EXCEPTION_H__
#define __OXLEY_EXCEPTION_H__

#include <escript/EsysException.h>

#include <exception>

namespace oxley {

/**
   \brief
   OxleyException exception class.
*/
class OxleyException : public escript::EsysException
{
public:
    OxleyException(const std::string& str) : escript::EsysException(str) {}
};

} // end of namespace oxley


/**
  \brief
  An exception class that signals an invalid argument value
*/
class ValueError : public escript::EsysException
{
public:
    ValueError(const std::string& str) : EsysException(str) {}
};

#endif
