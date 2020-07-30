
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
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

#ifndef __ESCRIPT_ESYSEXCEPTION_H__
#define __ESCRIPT_ESYSEXCEPTION_H__

#include <exception>
#include <string>

namespace escript {

/**
   \brief The base class for escript exceptions.
*/
class EsysException : public std::exception
{
public:
    /**
      \brief
      Constructor which creates an Exception with the given message
     
      @param message - Exception message.
    */
    EsysException(const std::string& message) : msg(message) {}

    /// Destructor
    virtual ~EsysException() throw() {}

    /**
      \brief
      Returns a description of the exception.
    */
    inline virtual const char* what() const throw() { return msg.c_str(); }

private:
    //
    // the exception message
    std::string msg;
};

/**
  \brief
  An exception class for assertions within escript
*/
class AssertException : public EsysException
{
public:
    AssertException(const std::string& str) : EsysException(str) {}
};

/**
  \brief
  An exception class for Input/Output errors
*/
class IOError : public EsysException
{
public:
    IOError(const std::string& str) : EsysException(str) {}
};

/**
  \brief
  An exception class for features which are not (yet) implemented
*/
class NotImplementedError : public EsysException
{
public:
    NotImplementedError(const std::string& str) : EsysException(str) {}
};

/**
  \brief
  An exception class that signals an invalid argument value
*/
class ValueError : public EsysException
{
public:
    ValueError(const std::string& str) : EsysException(str) {}
};

} // namespace escript

#endif // __ESCRIPT_ESYSEXCEPTION_H__

