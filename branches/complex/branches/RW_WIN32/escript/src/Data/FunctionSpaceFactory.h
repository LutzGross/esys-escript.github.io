/* 
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/
                                                                           
#if !defined  escript_FunctionSpaceFactory_20040604_H
#define escript_FunctionSpaceFactory_20040604_H
#ifdef MSVC
#ifdef ESCRIPT_EXPORTS
#define ESCRIPT_DLL __declspec(dllexport)
#else
#define ESCRIPT_DLL __declspec(dllimport)
#endif
#else
#define ESCRIPT_DLL
#endif
#include "escript/Data/AbstractDomain.h"
#include "escript/Data/FunctionSpace.h"

namespace escript {

  /**
     \brief
     Create function space objects.

     Description:
     Create function space objects.

	 TODO: Should be a class using the Factory Pattern in order to allow extension
  */

  /**
     \brief
     Return a continuous FunctionSpace
  */
  ESCRIPT_DLL FunctionSpace continuousFunction(const AbstractDomain& domain);

  /**
     \brief
     Return a functon FunctionSpace
  */
  ESCRIPT_DLL FunctionSpace function(const AbstractDomain& domain);
  /**
     \brief
     Return a function on boundary FunctionSpace
  */
  ESCRIPT_DLL FunctionSpace functionOnBoundary(const AbstractDomain& domain);
  /**
     \brief
     Return a FunctionSpace
  */
  ESCRIPT_DLL FunctionSpace functionOnContactZero(const AbstractDomain& domain);
  /**
     \brief
     Return a FunctionSpace
  */
  ESCRIPT_DLL FunctionSpace functionOnContactOne(const AbstractDomain& domain);
  /**
     \brief
     Return a FunctionSpace
  */
  ESCRIPT_DLL FunctionSpace solution(const AbstractDomain& domain);
  /**
     \brief
     Return a FunctionSpace
  */
  ESCRIPT_DLL FunctionSpace reducedSolution(const AbstractDomain& domain);
  /**
     \brief
     Return a FunctionSpace
  */
  ESCRIPT_DLL FunctionSpace diracDeltaFunction(const AbstractDomain& domain);

} // end of namespace
#endif
