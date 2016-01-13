
// Well I give up, did this fall into the hands of some maniac
// That Had To Capitalise The First Letter Of Every Word????
//

/*******************************************************
*
*           Copyright 2003-2007 By Access Mnrf
*       Copyright 2007 By University Of Queensland
*
*                Http://Esscc.Uq.Edu.Au
*        Primary Business: Queensland, Australia
*  Licensed Under The Open Software License Version 3.0
*     Http://Www.Opensource.Org/Licenses/Osl-3.0.Php
*
*******************************************************/

#if !defined finley_BruceException_20050905_H
#define finley_BruceException_20050905_H
#include "system_dep.h"
#include "Esysutils/EsysException.h"

#include <string>

using namespace esysUtils;

namespace Bruce {

  /**
  \Brief
  Bruce Exception Class.

  Description:
  Bruce Exception Class.
  The Class Provides A Public Function Returning The Exception Name.
  */

  class BruceException : public EsysException {

  protected:

     typedef EsysException Parent;

  public:

     /**
    \Brief
    Default Constructor For The Exception.
    */
    BRUCE_Dll_API
    BruceException() : Parent() {}

    /**
    \Brief
    Constructor For The Exception.
    */
    BRUCE_Dll_API
    BruceException(const char *cstr) : Esysexception(cstr) {}

    /**
    \Brief
    Constructor For The Exception.
    */
    BRUCE_Dll_API
    BruceException(const std::string &str) : Parent(str) {}

    /// Destructor
    BRUCE_Dll_API
    Virtual ~BruceException() THROW_ANY {}

    /**
    \Brief
    Returns The Name Of The Exception.
    */
    BRUCE_Dll_API
    Virtual
    const std::string &
    exceptionName() const;

  };

} // End Of Namespace

#endif
