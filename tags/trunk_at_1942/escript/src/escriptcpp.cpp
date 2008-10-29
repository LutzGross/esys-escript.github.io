
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "Data.h"
#include "FunctionSpace.h"
#include "FunctionSpaceFactory.h"
#include "DataFactory.h"
#include "AbstractContinuousDomain.h"
#include "AbstractDomain.h"
#include "Utils.h"
#include "AbstractSystemMatrix.h"
#include "AbstractTransportProblem.h"
#include "DataVector.h"
#include "paso/Paso_MPI.h"
#include "EscriptParams.h"



extern "C" {
#include "escript/blocktimer.h"
}

#include "esysUtils/esysExceptionTranslator.h"

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/object.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/numeric.hpp>
#include <boost/smart_ptr.hpp>

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
  def("blocktimer_initialize",blocktimer_initialize);
  def("blocktimer_reportSortByName",blocktimer_reportSortByName);
  def("blocktimer_reportSortByTime",blocktimer_reportSortByTime);
  def("blocktimer_increment",blocktimer_increment);
  def("blocktimer_time",blocktimer_time);
  def("getVersion",escript::getSvnVersion);
  def("printParallelThreadCounts",escript::printParallelThreadCnt);
  def("getMPISizeWorld",escript::getMPISizeWorld);
  def("getMPIRankWorld",escript::getMPIRankWorld);


  //
  // Interface for AbstractDomain
  //
  class_<escript::AbstractDomain, escript::Domain_ptr>("Domain",no_init)
     .def("setTagMap",&escript::AbstractDomain::setTagMap)
     .def("getTag",&escript::AbstractDomain::getTag)
     .def("isValidTagName",&escript::AbstractDomain::isValidTagName)
     .def("showTagNames",&escript::AbstractDomain::showTagNames)
     .def("getX",&escript::AbstractDomain::getX)
     .def("getDim",&escript::AbstractDomain::getDim)
     .def("getNormal",&escript::AbstractDomain::getNormal)
     .def("getSize",&escript::AbstractDomain::getSize)
     .def("saveVTK",&escript::AbstractDomain::saveVTK)
     .def("saveDX",&escript::AbstractDomain::saveDX)
     .def("getMPISize",&escript::AbstractDomain::getMPISize)
     .def("getMPIRank",&escript::AbstractDomain::getMPIRank)
     .def("MPIBarrier",&escript::AbstractDomain::MPIBarrier)
     .def("onMasterProcessor",&escript::AbstractDomain::onMasterProcessor)

     .def(self == self)
     .def(self != self);

  //
  // Interface for AbstractContinuousDomain
  //
  class_<escript::AbstractContinuousDomain, bases<escript::AbstractDomain> >("ContinuousDomain",no_init)
       .def("getSystemMatrixTypeId",&escript::AbstractContinuousDomain::getSystemMatrixTypeId)
       .def("getTransportTypeId",&escript::AbstractContinuousDomain::getTransportTypeId);

  //
  // Interface for FunctionSpace
  //
  class_<escript::FunctionSpace> fs_definer("FunctionSpace",init<>());
  fs_definer.def("getDim",&escript::FunctionSpace::getDim);
//   fs_definer.def("getDomain",&escript::FunctionSpace::getDomain,
//                  return_internal_reference<>());
  fs_definer.def("getDomain",&escript::FunctionSpace::getDomainPython);
  fs_definer.def("getX",&escript::FunctionSpace::getX);
  fs_definer.def("getNormal",&escript::FunctionSpace::getNormal);
  fs_definer.def("getSize",&escript::FunctionSpace::getSize);
  fs_definer.def("setTags",&escript::FunctionSpace::setTags);
  fs_definer.def("getTagFromDataPointNo",
                 &escript::FunctionSpace::getTagFromDataPointNo);
  fs_definer.def("getReferenceIDFromDataPointNo", &escript::FunctionSpace::getReferenceIDFromDataPointNo);
  fs_definer.def("getListOfTags",&escript::FunctionSpace::getListOfTags);
  fs_definer.def("__str__", &escript::FunctionSpace::toString);
  fs_definer.def(self == self);
  fs_definer.def(self != self);
  //
  // Interface for Data
  //
  class_<escript::Data>("Data","TEST DOCUMENTATION",init<>() )
    // various constructors for Data objects
    .def(init<const numeric::array&, optional<const escript::FunctionSpace&, bool> >(args("value","what","expand")))
    .def(init<const object&, optional<const escript::FunctionSpace&, bool> >(args("value","what","expand")))
    .def(init<const double, const tuple&, optional<const escript::FunctionSpace&, bool> >(args("value","shape","what","expand")))
    .def(init<const escript::Data&, const escript::FunctionSpace&>(args("value","what")))
    .def(init<const escript::Data&>())
    // Note for Lutz, Need to specify the call policy in order to return a
    // reference. In this case return_internal_reference.
    .def("__str__",&escript::Data::toString)
//     .def("getDomain",&escript::Data::getDomain,return_internal_reference<>())
    .def("getDomain",&escript::Data::getDomainPython)
    .def("getFunctionSpace",&escript::Data::getFunctionSpace,return_value_policy<copy_const_reference>())
    .def("isEmpty",&escript::Data::isEmpty)
    .def("isProtected",&escript::Data::isProtected)
    .def("setProtection",&escript::Data::setProtection)
    .def("getShape",&escript::Data::getShapeTuple)
    .def("getRank",&escript::Data::getDataPointRank)
    .def("dump",&escript::Data::dump)
    .def("copyWithMask",&escript::Data::copyWithMask)
    .def("setTaggedValue",&escript::Data::setTaggedValue)
    .def("setTaggedValue",&escript::Data::setTaggedValueByName)
    .def("getNumberOfDataPoints",&escript::Data::getNumDataPoints)
    .def("isExpanded",&escript::Data::isExpanded)
    .def("isTagged",&escript::Data::isTagged)
    .def("expand",&escript::Data::expand)
    .def("tag",&escript::Data::tag)
    .def("copy",&escript::Data::copy)
    .def("copy",&escript::Data::copySelf,return_value_policy<manage_new_object>())
    .def("setValueOfDataPoint",&escript::Data::setValueOfDataPointToPyObject)
    .def("setValueOfDataPoint",&escript::Data::setValueOfDataPointToArray)
    .def("setValueOfDataPoint",&escript::Data::setValueOfDataPoint)
    .def("getValueOfDataPoint",&escript::Data::getValueOfDataPoint)
    .def("getValueOfGlobalDataPoint",&escript::Data::getValueOfGlobalDataPoint)
    .def("setToZero",&escript::Data::setToZero)
    .def("interpolate",&escript::Data::interpolate)
    .def("minGlobalDataPoint",&escript::Data::minGlobalDataPoint)
    .def("saveDX",&escript::Data::saveDX)
    .def("saveVTK",&escript::Data::saveVTK)
    .def("getTagNumber",&escript::Data::getTagNumber)
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
  def("ReducedContinuousFunction",escript::reducedContinuousFunction);
  def("Function",escript::function);
  def("ReducedFunction",escript::reducedFunction);
  def("FunctionOnBoundary",escript::functionOnBoundary);
  def("ReducedFunctionOnBoundary",escript::reducedFunctionOnBoundary);
  def("FunctionOnContactZero",escript::functionOnContactZero);
  def("ReducedFunctionOnContactZero",escript::reducedFunctionOnContactZero);
  def("FunctionOnContactOne",escript::functionOnContactOne);
  def("ReducedFunctionOnContactOne",escript::reducedFunctionOnContactOne);
  def("Solution",escript::solution);
  def("ReducedSolution",escript::reducedSolution);
  def("DiracDeltaFunction",escript::diracDeltaFunction);

  //
  // Factory methods for Data
  //
  def("load",escript::load);
  def("loadIsConfigured",escript::loadConfigured);
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
  // Interface for AbstractTransportProblem
  //
  class_<escript::AbstractTransportProblem>("TransportProblem",init<>())
     .def("isEmpty",&escript::AbstractTransportProblem::isEmpty)
     .def("solve",&escript::AbstractTransportProblem::solve)
     .def("setInitialValue",&escript::AbstractTransportProblem::setInitialValue)
     .def("insertConstraint",&escript::AbstractTransportProblem::insertConstraint)
     .def("reset",&escript::AbstractTransportProblem::resetTransport)
     .def("resetValues",&escript::AbstractTransportProblem::resetTransport)
     .def("getSafeTimeStepSize",&escript::AbstractTransportProblem::getSafeTimeStepSize)
     .def("getUnlimitedTimeStepSize",&escript::AbstractTransportProblem::getUnlimitedTimeStepSize);

  // Functions to modify global parameters
  def("setEscriptParamInt",escript::setEscriptParamInt,
      (arg("value")=0));
  def("getEscriptParamInt",escript::getEscriptParamInt,
      (arg("sentinel")=0));


  //
  // Register esysExceptionTranslator
  //
  register_exception_translator<esysUtils::EsysException>(&esysUtils::esysExceptionTranslator);

}
