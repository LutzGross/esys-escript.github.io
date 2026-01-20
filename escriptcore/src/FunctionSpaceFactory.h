
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __ESCRIPT_FUNCTIONSPACEFACTORY_H__
#define __ESCRIPT_FUNCTIONSPACEFACTORY_H__

#include "system_dep.h"
#include "AbstractDomain.h"
#include "FunctionSpace.h"

namespace escript {

  /**
     \brief
     Create function space objects.

     Description:
     Create function space objects.

  */

  /**
     \brief
     Return a continuous FunctionSpace (overlapped node values)
  */
  ESCRIPT_DLL_API FunctionSpace continuousFunction(const AbstractDomain& domain);

  /**
     \brief
     Return a continuous with reduced order FunctionSpace (overlapped node values on reduced element order)
  */
  ESCRIPT_DLL_API FunctionSpace reducedContinuousFunction(const AbstractDomain& domain);

  /**
     \brief
     Return a function FunctionSpace
  */
  ESCRIPT_DLL_API FunctionSpace function(const AbstractDomain& domain);
  /**
     \brief
     Return a function FunctionSpace with reduced integration order
  */
  ESCRIPT_DLL_API FunctionSpace reducedFunction(const AbstractDomain& domain);
  /**
     \brief
     Return a function on boundary FunctionSpace
  */
  ESCRIPT_DLL_API FunctionSpace functionOnBoundary(const AbstractDomain& domain);
  /**
     \brief
     Return a function on boundary FunctionSpace with reduced integration order
  */
  ESCRIPT_DLL_API FunctionSpace reducedFunctionOnBoundary(const AbstractDomain& domain);
  /**
     \brief
     Return a FunctionSpace on left side of contact
  */
  ESCRIPT_DLL_API FunctionSpace functionOnContactZero(const AbstractDomain& domain);
  /**
     \brief
     Return a FunctionSpace  on left side of contact with reduced integration order
  */
  ESCRIPT_DLL_API FunctionSpace reducedFunctionOnContactZero(const AbstractDomain& domain);
  /**
     \brief
     Return a FunctionSpace on right side of contact
  */
  ESCRIPT_DLL_API FunctionSpace functionOnContactOne(const AbstractDomain& domain);
  /**
     \brief
     Return a FunctionSpace on right side of contact with reduced integration order
  */
  ESCRIPT_DLL_API FunctionSpace reducedFunctionOnContactOne(const AbstractDomain& domain);
  /**
     \brief
     Return a FunctionSpace
  */
  ESCRIPT_DLL_API FunctionSpace solution(const AbstractDomain& domain);
  /**
     \brief
     Return a FunctionSpace with reduced integration order
  */
  ESCRIPT_DLL_API FunctionSpace reducedSolution(const AbstractDomain& domain);
  /**
     \brief
     Return a FunctionSpace
  */
  ESCRIPT_DLL_API FunctionSpace diracDeltaFunctions(const AbstractDomain& domain);

} // end of namespace

#endif // __ESCRIPT_FUNCTIONSPACEFACTORY_H__

