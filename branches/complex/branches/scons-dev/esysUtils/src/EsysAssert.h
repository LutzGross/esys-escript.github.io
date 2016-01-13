
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#if !defined escript_EsysAssert_20040330_H
#define escript_EsysAssert_20040330_H
#include "system_dep.h"
/**
   \brief
   EsysAssert is a MACRO that will throw an exception if the boolean
   condition specified is false.

   Description:
   EsysAssert is conditionally compiled into code only when DOASSERT is
   defined.  When DOASSERT is not defined, the EsysAssert statement is
   entirely removed from code.
*/

//
// Note that the ANSI C Standard requires all headers to be idempotent except
// <assert.h> which is explicitly required not to be idempotent (section 4.1.2).
// This version of EsysAssert follows this requirement, consequently this
// part of the header is intentionally outside the single pass guard.
//

#undef EsysAssert

#if defined DOASSERT

//
// DOASSERT is defined, replace EsysAssert with Exception throw
//

#include "EsysAssertException.h"
#include <sstream>

namespace esysUtils {

  class ErrStream
  {
    public:
    template <typename Tmpl>
    ErrStream& operator<<(Tmpl t)
    {
      std::stringstream str;
      str << t;
      m_msg += str.str();
      
      return *this;
    }
    
    inline
    const std::string &toString() const
    {
      return m_msg;
    }

    private:
      std::string m_msg;
  };

  inline
  std::ostream& operator<<(std::ostream& oStream,
                                  const ErrStream& errStream)
  {
    oStream << errStream.toString();
    return oStream;
  }

}

#define EsysAssert(AssertTest,AssertMessage) \
   (void)((AssertTest) || \
           ((esysUtils::EsysAssertException::assertFailure(#AssertTest, __DATE__, __FILE__, __LINE__, \
             (esysUtils::ErrStream()<<AssertMessage).toString())),0),0)

#else

//
// DOASSERT os not defined, replace EsysAssert with "NO-OP"
//

#define EsysAssert(AssertTest,AssertMessage) ((void)0)

#endif

#endif
