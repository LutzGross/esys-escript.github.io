
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

#ifndef __RIPLEY_ADAPTEREXCEPTION_H__
#define __RIPLEY_ADAPTEREXCEPTION_H__

#include <ripley/system_dep.h>
#include <esysUtils/EsysException.h>

namespace ripley {

/**
   \brief
   RipleyException exception class.
*/
class RipleyException : public esysUtils::EsysException
{
protected:
    typedef EsysException Parent;

public:
    /**
       \brief
       Default constructor for the exception.
    */
    RIPLEY_DLL_API
    RipleyException() : Parent() { updateMessage(); }

    /**
       \brief
       Constructor with message.
    */
    RIPLEY_DLL_API
    RipleyException(const char *cstr) : Parent(cstr) { updateMessage(); }

    /**
       \brief
       Constructor with message.
    */
    RIPLEY_DLL_API
    RipleyException(const std::string &str) : Parent(str) { updateMessage(); }

    /**
       \brief
       Copy Constructor.
    */
    RIPLEY_DLL_API
    RipleyException(const RipleyException &other) : Parent(other)
    {
        updateMessage();
    }

    /// Destructor
    RIPLEY_DLL_API
    virtual ~RipleyException() THROW(NO_ARG) {}

    /**
       \brief
       Assignment operator.
    */
    RIPLEY_DLL_API
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
    RIPLEY_DLL_API
    virtual const std::string& exceptionName() const;

private:
    //
    // the exception name is immutable and class-wide.
    static const std::string exceptionNameValue;
};

} // end of namespace ripley

#endif // __RIPLEY_ADAPTEREXCEPTION_H__

