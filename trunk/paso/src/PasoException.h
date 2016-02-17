
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

#ifndef __PASO_EXCEPTION_H__
#define __PASO_EXCEPTION_H__

#include <esysUtils/EsysException.h>

namespace paso {

/**
  \brief
  PasoException exception class.

  Description:
  PasoException exception class.
  The class provides a public function returning the exception name
*/
class PasoException : public esysUtils::EsysException
{

protected:
    typedef EsysException Parent;

public:
    /**
    \brief
    Default constructor for the exception.
    */
    PasoException() : Parent() { updateMessage(); }

    /**
    \brief
    Constructor for the exception.
    */
    PasoException(const char *cstr) : Parent(cstr) { updateMessage(); }

    /**
    \brief
    Constructor for the exception.
    */
    PasoException(const std::string &str) : Parent(str) { updateMessage(); }

    /**
    \brief
    Copy Constructor for the exception.
    */
    PasoException(const PasoException &other) : Parent(other)
    {
        updateMessage();
    }

    /// Destructor
    virtual ~PasoException() THROW(NO_ARG) {}

    /**
    \brief
    Assignment operator.
    */
    inline PasoException& operator=(const PasoException &other ) THROW(NO_ARG)
    {
        Parent::operator=(other);
        updateMessage();
        return *this;
    }

    /**
    \brief
    Returns the name of the exception.
    */
    virtual const std::string & exceptionName() const;

private:
    //
    // the exception name is immutable and class-wide.
    // Inheritor note; you need one of these too.
    // and an overloaded exceptionName() in your .cpp implementation file. 
    static const std::string exceptionNameValue;
};

void checkPasoError(); 


} // end of namespace

#endif // __PASO_EXCEPTION_H__

