
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/

#ifndef __RIPLEY_EXCEPTION_H__
#define __RIPLEY_EXCEPTION_H__

#include <ripley/system_dep.h>
#include <esysUtils/EsysException.h>

namespace ripley {

/**
   \brief
   RipleyException exception class.
*/
class RIPLEY_DLL_API RipleyException : public esysUtils::EsysException
{
protected:
    typedef EsysException Parent;

public:
    /**
       \brief
       Default constructor for the exception.
    */
    RipleyException() : Parent() { updateMessage(); }

    /**
       \brief
       Constructor with message.
    */
    RipleyException(const char *cstr) : Parent(cstr) { updateMessage(); }

    /**
       \brief
       Constructor with message.
    */
    RipleyException(const std::string &str) : Parent(str) { updateMessage(); }

    /**
       \brief
       Copy Constructor.
    */
    RipleyException(const RipleyException &other) : Parent(other)
    {
        updateMessage();
    }

    /// Destructor
    virtual ~RipleyException() THROW(NO_ARG) {}

    /**
       \brief
       Assignment operator.
    */
    inline RipleyException& operator=(const RipleyException &other ) THROW(NO_ARG)
    {
        Parent::operator=(other);
        updateMessage();
        return *this;
    }

    /**
       \brief
       Returns the name of the exception.
    */
    virtual const std::string& exceptionName() const;

private:
    //
    // the exception name is immutable and class-wide.
    static const std::string exceptionNameValue;
};

} // end of namespace ripley

#endif // __RIPLEY_EXCEPTION_H__

