//$Id$
/*=============================================================================
 
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT ACcESS 2004 -  All Rights Reserved                         *
 *                                                                            *
 * This software is the property of ACcESS.  No part of this code             *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that                          *
 * person has a software license agreement with ACcESS.                       *
 *                                                                            *
 ******************************************************************************
 
******************************************************************************/

#include "escript/Data/Data.h"
#include "escript/Data/FunctionSpace.h"
#include "escript/Data/FunctionSpaceFactory.h"
#include "escript/Data/DataFactory.h"
#include "escript/Data/AbstractContinuousDomain.h"
#include "escript/Data/AbstractDomain.h"
#include "esysUtils/esysExceptionTranslator.h"

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/object.hpp>
#include <boost/python/tuple.hpp>
                                                                                                                                          
using namespace boost::python;

/**
   @memo
   escript is the python module that contains the interfaces
   to the C++ side of escript.

   @version 1.0.0 

   @doc

   Class Description:
   Data

   Class Limitations:
   None

   Class Conditions of Use:
   None

   Throws:
   None

*/

BOOST_PYTHON_MODULE(escriptcpp)
{

  //
  // Interface for AbstractDomain
  //
  class_<escript::AbstractDomain>("Domain",no_init)
     .def("getX",&escript::AbstractDomain::getX)
     .def("getNormal",&escript::AbstractDomain::getNormal)
     .def("getSize",&escript::AbstractDomain::getSize)
     .def(self == self)
     .def(self != self);

  //
  // Interface for AbstractContinuousDomain
  //
  class_<escript::AbstractContinuousDomain, bases<escript::AbstractDomain> >("ContinuousDomain",no_init)
       .def("getSystemMatrixTypeId",&escript::AbstractContinuousDomain::getSystemMatrixTypeId);

  //
  // Interface for FunctionSpace
  //
  class_<escript::FunctionSpace>("FunctionSpace",init<>())
     .def("getDim",&escript::FunctionSpace::getDim)
     .def("getX",&escript::FunctionSpace::getX)
     .def("getNormal",&escript::FunctionSpace::getNormal)
     .def("getSize",&escript::FunctionSpace::getSize)
     .def("toString",&escript::FunctionSpace::toString)
     .def(self == self)
     .def(self != self);

  //
  // Interface for Data
  //
  class_<escript::Data>("Data",init<>())
    // various constructors for Data objects
    .def(init<const numeric::array&, optional<const escript::FunctionSpace&, bool> >(args("value","what","expand")))
    .def(init<const object&, optional<const escript::FunctionSpace&, bool> >(args("value","what","expand")))
    .def(init<const double, const tuple&, optional<const escript::FunctionSpace&, bool> >(args("value","shape","what","expand")))
    .def(init<const escript::Data&, const escript::FunctionSpace&>(args("value","what")))
    .def(init<const escript::Data&>())
    // Note for Lutz, Need to specify the call policy in order to return a
    // reference. In this case return_internal_reference.
    .def("toString",&escript::Data::toString)
    .def("getDomain",&escript::Data::getDomain,return_internal_reference<>())
    .def("getFunctionSpace",&escript::Data::getFunctionSpace,return_internal_reference<>())
    .def("isEmpty",&escript::Data::isEmpty)
    .def("getShape",&escript::Data::getShapeTuple)
    .def("getRank",&escript::Data::getDataPointRank)
    .def("copyWithMask",&escript::Data::copyWithMask)
    .def("setTaggedValue",&escript::Data::setTaggedValue)
    .def("setRefValue",&escript::Data::setRefValue)
    .def("getRefValue",&escript::Data::getRefValue)
    //.def("expand",&escript::Data::expand)
    //.def("tag",&escript::Data::tag)
    .def("saveDX",&escript::Data::saveDX)
    .def("saveVTK",&escript::Data::saveVTK)
    .def("wherePositive",&escript::Data::wherePositive)
    .def("whereNegative",&escript::Data::whereNegative)
    .def("whereNonNegative",&escript::Data::whereNonNegative)
    .def("whereNonPositive",&escript::Data::whereNonPositive)
    .def("whereZero",&escript::Data::whereZero)
    .def("whereNonZero",&escript::Data::whereNonZero)
    // Unary functions for Data
    .def("interpolate",&escript::Data::interpolate)
    .def("grad",&escript::Data::gradOn)
    .def("grad",&escript::Data::grad)
    .def("integrate",&escript::Data::integrate)
    .def("transpose",&escript::Data::transpose)
    .def("trace",&escript::Data::trace)
    .def("sin",&escript::Data::sin)
    .def("cos",&escript::Data::cos)
    .def("tan",&escript::Data::tan)
    .def("log",&escript::Data::log)
    .def("ln",&escript::Data::ln)
    .def("Lsup",&escript::Data::Lsup)
    .def("sup",&escript::Data::sup)
    .def("inf",&escript::Data::inf)
    .def("__abs__",&escript::Data::abs)
    .def("exp",&escript::Data::exp)
    .def("sqrt",&escript::Data::sqrt)
    .def("maxval",&escript::Data::maxval)
    .def("minval",&escript::Data::minval)
    .def("length",&escript::Data::length)
    .def("sign",&escript::Data::sign)
    // following implements the python "-" negation operator
    .def("__neg__",&escript::Data::neg)
    // following implements the python "+" identity operator
    .def("__pos__",&escript::Data::pos)
    // following two functions implement the python [] operator
    .def("__getitem__",&escript::Data::getItem)
    .def("__setitem__",&escript::Data::setItemO)
    .def("__setitem__",&escript::Data::setItemD)
    // following two functions implement the python ** operator
    .def("__pow__",&escript::Data::powO)
    .def("__pow__",&escript::Data::powD)
    // NOTE:: The order of these declarations is important. Anything
    // declared before the generic declaration isn't found so the generic
    // version will be called. 
    .def(self += other<object>())
    .def(self += self)
    .def(self + other<object>())
    .def(other<object>() + self)
    .def(self + self)
    .def(self -= other<object>())
    .def(self -= self)
    .def(self - other<object>())
    .def(other<object>() - self)
    .def(self - self)
    .def(self *= other<object>())
    .def(self *= self)
    .def(self * other<object>())
    .def(other<object>() * self)
    .def(self * self)
    .def(self /= other<object>())
    .def(self /= self)
    .def(self / other<object>())
    .def(other<object>() / self)
    .def(self / self)
    // Need scope resolution due to a bug either in the compiler or
    // the boost code. This calls operator << for Data.
    .def(self_ns::str(self));

  //
  // Factory methods for function space
  //
  def("ContinuousFunction",escript::continuousFunction);
  def("Function",escript::function);
  def("FunctionOnBoundary",escript::functionOnBoundary);
  def("FunctionOnContactZero",escript::functionOnContactZero);
  def("FunctionOnContactOne",escript::functionOnContactOne);
  def("Solution",escript::solution);
  def("ReducedSolution",escript::reducedSolution);
  def("DiracDeltaFunction",escript::diracDeltaFunction);

  //
  // Factory methods for Data
  //
  def("Scalar",escript::Scalar,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false));
  def("Vector",escript::Vector,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false));
  def("Tensor",escript::Tensor,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false));
  def("Tensor3",escript::Tensor3,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false));
  def("Tensor4",escript::Tensor4,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false));

  //
  // Interface for AbstractSystemMatrix
  //
  class_<escript::AbstractSystemMatrix>("Operator",init<>())
     .def("isEmpty",&escript::AbstractSystemMatrix::isEmpty)
     .def("solve",&escript::AbstractSystemMatrix::solve)
     .def("of",&escript::AbstractSystemMatrix::vectorMultiply)
     .def("saveMM",&escript::AbstractSystemMatrix::saveMM)
     .def("setValue",&escript::AbstractSystemMatrix::setValue)
     .def("resetSolver",&escript::AbstractSystemMatrix::resetSolver)
     .def(self*other<escript::Data>());

  //
  // Register esysExceptionTranslator
  //
  register_exception_translator<esysUtils::EsysException>(&esysUtils::esysExceptionTranslator);

}
