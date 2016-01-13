
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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

#ifndef __SOLVEROPTIONS_EXCEPTION_H__
#define __SOLVEROPTIONS_EXCEPTION_H__

#include "system_dep.h"
#include <esysUtils/EsysException.h>

namespace escript {

/**
   \brief
   SolverOptionsException exception class.
*/
class ESCRIPT_DLL_API SolverOptionsException : public esysUtils::EsysException
{
protected:
    typedef EsysException Parent;

public:
    /**
       \brief
       Default constructor for the exception.
    */
    SolverOptionsException() : Parent() { updateMessage(); }

    /**
       \brief
       Constructor with message.
    */
    SolverOptionsException(const char *cstr) : Parent(cstr) { updateMessage(); }

    /**
       \brief
       Constructor with message.
    */
    SolverOptionsException(const std::string &str) : Parent(str) { updateMessage(); }

    /**
       \brief
       Copy Constructor.
    */
    SolverOptionsException(const SolverOptionsException &other) : Parent(other)
    {
        updateMessage();
    }

    /// Destructor
    virtual ~SolverOptionsException() THROW(NO_ARG) {}

    /**
       \brief
       Assignment operator.
    */
    inline SolverOptionsException& operator=(const SolverOptionsException &other ) THROW(NO_ARG)
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

}

#endif

