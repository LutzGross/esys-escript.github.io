//$Id$
/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

#include "Data.h"
#include "FunctionSpace.h"
#include "FunctionSpaceFactory.h"
#include "DataFactory.h"
#include "AbstractContinuousDomain.h"
#include "AbstractDomain.h"
#include "Utils.h"
#include "AbstractSystemMatrix.h"
#include "DataVector.h"

#include "esysUtils/esysExceptionTranslator.h"

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/object.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/numeric.hpp>

using namespace boost::python;

/*! \mainpage Esys Documentation
 *
 * \version 1.0.0
 *
 * - \ref escript
 *
 * - \ref esys_exception "Esys Exception"
 *
 * - \ref finley
 *
 * - <a href=http://iservo.edu.au/esys/epydoc/index.html>Python module documentation (epydoc generated)</a>
 *
 */

/*! \page escript Escript
 * Escript is the python module that contains the interfaces
 * to the C++ side of escript.
 *
 * \version 1.0.0 
 *
 * \section class_desc Class Description:
 * Data
 *
 * \section class_limits Class Limitations:
 * None
 *
 * \section class_conds Class Conditions of Use:
 * None
 *
 * \section class_throws Throws:
 * None
 *
 */

BOOST_PYTHON_MODULE(escriptcpp)
{
  def("setNumberOfThreads",escript::setNumberOfThreads);
  def("getNumberOfThreads",escript::getNumberOfThreads);
  def("releaseUnusedMemory",escript::releaseUnusedMemory);


  //
  // Interface for AbstractDomain
  //
  class_<escript::AbstractDomain>("Domain",no_init)
     .def("getX",&escript::AbstractDomain::getX)
     .def("getNormal",&escript::AbstractDomain::getNormal)
     .def("getSize",&escript::AbstractDomain::getSize)
     .def("saveVTK",&escript::AbstractDomain::saveVTK)
     .def("saveDX",&escript::AbstractDomain::saveDX)
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
     .def("getDomain",&escript::FunctionSpace::getDomain,return_internal_reference<>())
     .def("getX",&escript::FunctionSpace::getX)
     .def("getNormal",&escript::FunctionSpace::getNormal)
     .def("getSize",&escript::FunctionSpace::getSize)
     .def("setTags",&escript::FunctionSpace::setTags)
     .def("getTagFromDataPointNo",&escript::FunctionSpace::getTagFromDataPointNo)
     .def("__str__",&escript::FunctionSpace::toString)
     .def(self == self)
     .def(self != self);
  //
  // Interface for Data
  //
  class_<escript::Data>("Data","TEST DOCUMENTATION",init<>())
    // various constructors for Data objects
    .def(init<const numeric::array&, optional<const escript::FunctionSpace&, bool> >(args("value","what","expand")))
    .def(init<const object&, optional<const escript::FunctionSpace&, bool> >(args("value","what","expand")))
    .def(init<const double, const tuple&, optional<const escript::FunctionSpace&, bool> >(args("value","shape","what","expand")))
    .def(init<const escript::Data&, const escript::FunctionSpace&>(args("value","what")))
    .def(init<const escript::Data&>())
    // Note for Lutz, Need to specify the call policy in order to return a
    // reference. In this case return_internal_reference.
    .def("__str__",&escript::Data::toString)
    .def("getDomain",&escript::Data::getDomain,return_internal_reference<>())
    .def("getFunctionSpace",&escript::Data::getFunctionSpace,return_internal_reference<>())
    .def("isEmpty",&escript::Data::isEmpty)
    .def("isProtected",&escript::Data::isProtected)
    .def("setProtection",&escript::Data::setProtection)
    .def("getShape",&escript::Data::getShapeTuple)
    .def("getRank",&escript::Data::getDataPointRank)
    .def("copyWithMask",&escript::Data::copyWithMask)
    .def("setTaggedValue",&escript::Data::setTaggedValue)
    .def("setRefValue",&escript::Data::setRefValue)
    .def("getRefValue",&escript::Data::getRefValue)
    .def("expand",&escript::Data::expand)
    .def("tag",&escript::Data::tag)
    .def("copy",&escript::Data::copy)
    .def("convertToNumArray",&escript::Data::convertToNumArray)
    .def("getNumberOfSamples",&escript::Data::getNumSamples)
    .def("convertToNumArrayFromSampleNo",&escript::Data::convertToNumArrayFromSampleNo)
    .def("convertToNumArrayFromDPNo",&escript::Data::convertToNumArrayFromDPNo)
    .def("fillFromNumArray",&escript::Data::fillFromNumArray)
    .def("interpolate",&escript::Data::interpolate)
    .def("mindp",&escript::Data::mindp)
    .def("saveDX",&escript::Data::saveDX)
    .def("saveVTK",&escript::Data::saveVTK)
    .def("getTagNumber",&escript::Data::getTagNumber)
    .def("archiveData",&escript::Data::archiveData)
    .def("extractData",&escript::Data::extractData)
    // Unary functions for Data
    .def("_interpolate",&escript::Data::interpolate)
    .def("_grad",&escript::Data::gradOn)
    .def("_grad",&escript::Data::grad)
    .def("_transpose",&escript::Data::transpose)
    .def("_trace",&escript::Data::trace)
    .def("_maxval",&escript::Data::maxval)
    .def("_minval",&escript::Data::minval)
    .def("_wherePositive",&escript::Data::wherePositive)
    .def("_whereNegative",&escript::Data::whereNegative)
    .def("_whereNonNegative",&escript::Data::whereNonNegative)
    .def("_whereNonPositive",&escript::Data::whereNonPositive)
    .def("_whereZero",&escript::Data::whereZero,(arg("tol")=0.0))
    .def("_whereNonZero",&escript::Data::whereNonZero,(arg("tol")=0.0))
    .def("_erf",&escript::Data::erf)
    .def("_sin",&escript::Data::sin)
    .def("_cos",&escript::Data::cos)
    .def("_tan",&escript::Data::tan)
    .def("_asin",&escript::Data::asin)
    .def("_acos",&escript::Data::acos)
    .def("_atan",&escript::Data::atan)
    .def("_sinh",&escript::Data::sinh)
    .def("_cosh",&escript::Data::cosh)
    .def("_tanh",&escript::Data::tanh)
    .def("_asinh",&escript::Data::asinh)
    .def("_acosh",&escript::Data::acosh)
    .def("_atanh",&escript::Data::atanh)
    .def("_exp",&escript::Data::exp)
    .def("_sqrt",&escript::Data::sqrt)
    .def("_log10",&escript::Data::log10)
    .def("_log",&escript::Data::log)
    .def("_sign",&escript::Data::sign)
    .def("_symmetric",&escript::Data::symmetric)
    .def("_nonsymmetric",&escript::Data::nonsymmetric)
    .def("_trace",&escript::Data::trace)
    .def("_swap_axes",&escript::Data::swapaxes)
    .def("_eigenvalues",&escript::Data::eigenvalues)
    .def("_eigenvalues_and_eigenvectors",&escript::Data::eigenvalues_and_eigenvectors,(arg("tol")=1.e-13))
    // functions returning a single real number:
    .def("_Lsup",&escript::Data::Lsup)
    .def("_sup",&escript::Data::sup)
    .def("_inf",&escript::Data::inf)
    .def("_integrate",&escript::Data::integrate)

    // following implements the python abs operator
    .def("__abs__",&escript::Data::abs)
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
    .def("__rpow__",&escript::Data::rpowO)
    // NOTE:: The order of these declarations is important. Anything
    // declared before the generic declaration isn't found so the generic
    // version will be called. 
    .def(self + other<object>())
    .def(other<object>() + self)
    .def(self + self)
    .def(self += other<object>())
    .def(self += self)

    .def(self - other<object>())
    .def(other<object>() - self)
    .def(self - self)
    .def(self -= other<object>())
    .def(self -= self)

    .def(self * other<object>())
    .def(other<object>() * self)
    .def(self * self)
    .def(self *= other<object>())
    .def(self *= self)

    .def(self / other<object>())
    .def(other<object>() / self)
    .def(self / self)
    .def(self /= other<object>())
    .def(self /= self)
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
  // Binary operators
  //
  def("C_GeneralTensorProduct",escript::C_GeneralTensorProduct,
      (arg("arg0")=escript::Data(),
       arg("arg1")=escript::Data(),
       arg("axis_offset")=0,
       arg("transpose")=0));

  //
  // Interface for AbstractSystemMatrix
  //
  class_<escript::AbstractSystemMatrix>("Operator",init<>())
     .def("isEmpty",&escript::AbstractSystemMatrix::isEmpty)
     .def("solve",&escript::AbstractSystemMatrix::solve)
     .def("of",&escript::AbstractSystemMatrix::vectorMultiply)
     .def("saveMM",&escript::AbstractSystemMatrix::saveMM)
     .def("saveHB",&escript::AbstractSystemMatrix::saveHB)
     .def("resetValues",&escript::AbstractSystemMatrix::resetValues)
     .def(self*other<escript::Data>());

  //
  // Register esysExceptionTranslator
  //
  register_exception_translator<esysUtils::EsysException>(&esysUtils::esysExceptionTranslator);

}
