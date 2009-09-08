
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
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
#include "TestDomain.h"


extern "C" {
#include "esysUtils/blocktimer.h"
}

#include "esysUtils/esysExceptionTranslator.h"

#include <boost/version.hpp>
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/object.hpp>
#include <boost/python/tuple.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/version.hpp>

using namespace boost::python;

/*! \mainpage Esys Documentation
 *
 * \version 3.0.0
 *
 * - \ref escript
 *
 * - \ref esys_exception "Esys Exception"
 *
 * - \ref finley
 *
 * - <a href="../../epydoc/index.html">Python module documentation (epydoc generated)</a>
 *
 */

/*! \page escript Escript
 * Escript is the python module that contains the interfaces
 * to the C++ side of escript.
 *
 * 
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

  #if BOOST_VERSION >= 103500
// params are: bool show_user_defined, bool show_py_signatures, bool show_cpp_signatures
  docstring_options docopt(true,true,false);
  #endif

  def("setNumberOfThreads",escript::setNumberOfThreads,"Use of this method is strongly discouraged.");
  def("getNumberOfThreads",escript::getNumberOfThreads,"Return the maximum number of threads"
" available to OpenMP.");
  def("releaseUnusedMemory",escript::releaseUnusedMemory);
  def("blocktimer_initialize",blocktimer_initialize);
  def("blocktimer_reportSortByName",blocktimer_reportSortByName);
  def("blocktimer_reportSortByTime",blocktimer_reportSortByTime);
  def("blocktimer_increment",blocktimer_increment);
  def("blocktimer_time",blocktimer_time);
  def("getVersion",escript::getSvnVersion,"This method will only report accurate version numbers for clean checkouts.");
  def("printParallelThreadCounts",escript::printParallelThreadCnt);
  def("getMPISizeWorld",escript::getMPISizeWorld,"Return number of MPI processes in the job.");
  def("getMPIRankWorld",escript::getMPIRankWorld,"Return the rank of this process in the MPI World.");
  def("MPIBarrierWorld",escript::MPIBarrierWorld,"Wait until all MPI processes have reached this point.");
  def("getMPIWorldMax",escript::getMPIWorldMax,"\nEach MPI process calls this function with a"
" value for arg1. The maximum value is computed and returned.\n\n:rtype: int");
  def("getMPIWorldSum",escript::getMPIWorldSum,"\nEach MPI process calls this function with a"
" value for arg1. The values are added up and the total value is returned.\n\n:rtype: int");
  def("getMachinePrecision",escript::getMachinePrecision);
  def("getMaxFloat",escript::getMaxFloat);
  def("_saveDataCSV",escript::saveDataCSV, (args("filename","arg","sep","csep"), arg("append")=false),
"Saves data objects passed in a python dictionary to a file.\n"
"The data objects must be over the same domain and be able to be interpolated to the same FunctionSpace.\n"
"If one of the dictionary keys is named ``mask``, then only samples where ``mask`` has a positive\n"
"value will be written to the file.\n\n"
"A header line giving the names of each column will be output first.\n"
"The keys given in the dictionary will be used to name columns.\n"
"Then the data will be output, one line per sample (for all data).\n"
"\n:param filename:\n:type filename: ``string``\n"
":param arg: dictionary of named `Data` objects. If one is called ``mask`` it must be scalar data."
"\n:type arg: ``dict``\n"
":param sep:separator for columns (defaults to \", \"\n:type sep:``string``\n"
":param csep:separator for fields within data object (defaults to \"_\"\n:type csep:``string``\n"
":param append: If True, write to the end of ``filename``\n:type append: ``string``\n");


  //
  // Interface for AbstractDomain
  //
  class_<escript::AbstractDomain, escript::Domain_ptr>("Domain","Base class for all domains.",no_init)
     .def("getStatus",&escript::AbstractDomain::getStatus,"The status of a domain changes whenever the domain is modified\n\n:rtype: int")
     .def("setTagMap",&escript::AbstractDomain::setTagMap,args("name","tag"),
"Give a tag number a name.\n\n:param name: Name for the tag\n:type name: ``string``\n"
":param tag: numeric id\n:type tag: ``int``\n:note: Tag names must be unique within a domain")
     .def("getTag",&escript::AbstractDomain::getTag,args("name"),":return: tag id for "
"``name``\n:rtype: ``string``")
     .def("isValidTagName",&escript::AbstractDomain::isValidTagName,args("name"),
":return: True is ``name`` corresponds to a tag\n:rtype: ``bool``")
     .def("showTagNames",&escript::AbstractDomain::showTagNames,":return: A space separated list of tag names\n:rtype: ``string``")
     .def("getX",&escript::AbstractDomain::getX,":rtype: `Data`\n:return: Locations in the"
"`Domain`. FunctionSpace is chosen appropriately")
     .def("getDim",&escript::AbstractDomain::getDim,":rtype: `int`\n:return: Spatial dimension of the `Domain`")
     .def("getNormal",&escript::AbstractDomain::getNormal,":rtype: `escript`\n:return: Boundary normals")
     .def("getSize",&escript::AbstractDomain::getSize,":return: the local size of samples. The function space is chosen appropriately\n:rtype: `Data`")
     .def("saveVTK",&escript::AbstractDomain::saveVTK,args("filename","arg", "metadata"
, "metadata_schema"),
":param filename: \n:type filename: ``string``\n:param arg: items to be added"
"\n:type arg: ``dict``\n:param metadata: string representing some meta data to be added\n:type metadata: ``string``"
"\n:param metadata_schema: schema type for metadata\n:type metadata_schema: ``string``")
     .def("dump",&escript::AbstractDomain::dump,args("filename"),"Dumps the domain to a file"
":param filename:\n:type filename: string")
     .def("saveDX",&escript::AbstractDomain::saveDX,args("filename","arg"),"Saves a dictonary of Data objects to an OpenDX input file.\n\n:param filename:\n:type filename: ``string``"
"\n:param arg:\n:type arg: `dict`")
     .def("getMPISize",&escript::AbstractDomain::getMPISize,":return: the number of processes used for this `Domain`\n:rtype: ``int``")
     .def("getMPIRank",&escript::AbstractDomain::getMPIRank,":return: the rank of this process\n:rtype: ``int``")
     .def("MPIBarrier",&escript::AbstractDomain::MPIBarrier,"Wait until all processes have reached this point")
     .def("onMasterProcessor",&escript::AbstractDomain::onMasterProcessor,":return: True if this code is executing on the master process\n:rtype: `bool`")

     .def(self == self)
     .def(self != self);

  //
  // Interface for AbstractContinuousDomain
  //
  class_<escript::AbstractContinuousDomain, bases<escript::AbstractDomain> >("ContinuousDomain","Class representing continuous domains",no_init)
       .def("getSystemMatrixTypeId",&escript::AbstractContinuousDomain::getSystemMatrixTypeId,
args("solver", "preconditioner", "package", "symmetry"),
":return: the identifier of the matrix type to be used for the global stiffness matrix "
"when a particular solver package and symmetric matrix is used.\n"
":rtype: int")
       .def("getTransportTypeId",&escript::AbstractContinuousDomain::getTransportTypeId,
args("solver", "preconditioner", "package", "symmetry"));


  //
  // Interface for TestDomain
  //
  class_ <escript::TestDomain, bases<escript::AbstractDomain> >("TestDomain", "Test Class for domains with no structure. May be removed from future releases without notice.", init<int,int>());

  // This is the only python visible way to get a TestDomain
  def("getTestDomainFunctionSpace",&escript::getTestDomainFunctionSpace, args("dpps", 
 "samples"),
"For testing only. May be removed without notice.");

  //
  // Interface for FunctionSpace
  //
  class_<escript::FunctionSpace> fs_definer("FunctionSpace","A FunctionSpace describes which points from the `Domain` to use to represent functions.",init<>());	// Doco goes in the empty string param
  fs_definer.def("getDim",&escript::FunctionSpace::getDim,":return: the spatial dimension of the underlying domain.\n:rtype: int");
//   fs_definer.def("getDomain",&escript::FunctionSpace::getDomain,
//                  return_internal_reference<>());
  fs_definer.def("getDomain",&escript::FunctionSpace::getDomainPython,":return: the underlying `Domain` for this FunctionSpace.\n:rtype: `Domain`");
  fs_definer.def("getX",&escript::FunctionSpace::getX,"\n:return: a function whose values are its input coordinates. ie an identity function.\n:rtype: `Data`");
  fs_definer.def("getNormal",&escript::FunctionSpace::getNormal,":return: the surface normal field.\n\n:rtype: `Data`");
  fs_definer.def("getSize",&escript::FunctionSpace::getSize,":return: sample size\n:rtype: `Data`");
  fs_definer.def("setTags",&escript::FunctionSpace::setTags,args("newtag","mask"),
"Set tags according to a mask\n\n:param newtag: tag number to set\n:type newtag: non-zero ``int``\n:param mask: Samples which correspond to positive values in the mask will be set to ``newtag``.\n:type mask: scalar `Data`");
  fs_definer.def("getTagFromDataPointNo",
                 &escript::FunctionSpace::getTagFromDataPointNo,":return: the tag associated with the given sample number.\n:rtype: int");
  fs_definer.def("getReferenceIDFromDataPointNo", &escript::FunctionSpace::getReferenceIDFromDataPointNo,args("dataPointNo"),":return: the reference number associated with ``dataPointNo``\n:rtype: int ");
  fs_definer.def("getListOfTags",&escript::FunctionSpace::getListOfTags,":return: a list of the tags used in this function space\n:rtype: ``list``");
  fs_definer.def("__str__", &escript::FunctionSpace::toString);
  fs_definer.def(self == self);
  fs_definer.def(self != self);
  //
  // Interface for Data
  //
  class_<escript::Data>("Data","Represents a collection of datapoints. It is used to store the values of a function. For more details please consult the c++ class documentation.",init<>() )
    // various constructors for Data objects
    .def(init<const object&, optional<const escript::FunctionSpace&, bool> >(args("value","what","expand")))
    .def(init<const double, const tuple&, optional<const escript::FunctionSpace&, bool> >(args("value","shape","what","expand")))
    .def(init<const escript::Data&, const escript::FunctionSpace&>(args("value","what")))
    .def(init<const escript::Data&>())
    // Note for Lutz, Need to specify the call policy in order to return a
    // reference. In this case return_internal_reference.
    .def("__str__",&escript::Data::toString)
    .def("getDomain",&escript::Data::getDomainPython,":rtype: `Domain`")
    .def("getFunctionSpace",&escript::Data::getFunctionSpace,return_value_policy<copy_const_reference>(),":rtype: `FunctionSpace`")
    .def("isEmpty",&escript::Data::isEmpty,"Is this object an instance of ``DataEmpty``\n\n:rtype: ``bool``\n:note: This is not the same thing as asking if the object contains datapoints.")
    .def("isProtected",&escript::Data::isProtected,"Can this instance be modified.\n:rtype: ``bool``")
    .def("setProtection",&escript::Data::setProtection,"Disallow modifications to this data object\n\n:note: This method does not allow you to undo protection.")
    .def("getShape",&escript::Data::getShapeTuple,"\nReturns the shape of the datapoints in this object as a python tuple. Scalar data has the shape ``()``\n\n:rtype: ``tuple``")
    .def("getRank",&escript::Data::getDataPointRank,":return: the number of indicies required to address a component of a datapoints\n:rtype: positive ``int``")
    .def("dump",&escript::Data::dump,args("fileName"),"Save the data as a netCDF file\n\n:param fileName: \n:type fileName: ``string``")
    .def("toListOfTuples",&escript::Data::toListOfTuples, (arg("scalarastuple")=false),
"Return the datapoints of this object in a list. Each datapoint is stored as a tuple.\n\n"
":param scalarastuple: if True, scalar data will be wrapped as a tuple."
" True => [(0), (1), (2)]; False => [0, 1, 2]")
    .def("copyWithMask",&escript::Data::copyWithMask,args("other","mask"),
"Selectively copy values from ``other`` `Data`."
"Datapoints which correspond to postive values in ``mask`` will be copied from ``other``\n"
"\n:param other: source of values\n"
":type other: `Data`\n:param mask:\n:type mask: Scalar `Data`")
    .def("setTaggedValue",&escript::Data::setTaggedValue,args("tagKey","value"),
"Set the value of tagged Data.\n\n:param tagKey: tag to update\n:type tagKey: ``int``\n")
    .def("setTaggedValue",&escript::Data::setTaggedValueByName,args("name","value"),":param name: tag to update\n:type name: ``string``\n"
":param value: value to set tagged data to\n:type value: ``object`` which acts like an array, ``tuple`` or ``list``\n")
    .def("getNumberOfDataPoints",&escript::Data::getNumDataPoints,":rtype: ``int``\n:return: Number of datapoints in the object")
    .def("isExpanded",&escript::Data::isExpanded,":rtype: ``bool``\n:return: True if this ``Data`` is expanded.")
    .def("isTagged",&escript::Data::isTagged,":rtype: ``bool``\n:return: True if this ``Data`` is expanded.")
    .def("isConstant",&escript::Data::isConstant,":rtype: ``bool``\n:return: True if this ``Data`` is an instance of ``DataConstant``\n:note: This does not mean the data is immutable.")
    .def("isLazy",&escript::Data::isLazy,":rtype: ``bool``\n:return: True if this ``Data`` is lazy.")
    .def("isReady",&escript::Data::isReady,":rtype: ``bool``\n:return: True if this ``Data`` is not lazy.")
    .def("expand",&escript::Data::expand,"Convert the data to expanded representation if it is not expanded already.")
    .def("tag",&escript::Data::tag,"Convert data to tagged representation if it is not already tagged or expanded")
    .def("resolve",&escript::Data::resolve,"Convert the data to non-lazy representation.")
    .def("copy",&escript::Data::copy,args("other"),"Make this object a copy of ``other``\n"
"\n:note: The two objects will act independently from now on. That is, changing ``other`` "
"after this call will not change this object and vice versa.")
    .def("copy",&escript::Data::copySelf,":note: In the no argument form, a new object will be returned which is an independent copy of this object.")
    .def("delay",&escript::Data::delay,"Convert this object into lazy representation")
    .def("setValueOfDataPoint",&escript::Data::setValueOfDataPointToPyObject,args("dataPointNo","value"))
    .def("setValueOfDataPoint",&escript::Data::setValueOfDataPointToArray)
    .def("setValueOfDataPoint",&escript::Data::setValueOfDataPoint,"\nModify the value of a single datapoint.\n\n:param dataPointNo:\n"
":type dataPointNo: int\n:param value: \n:type value: ``float`` or an object which acts like an array, ``tuple`` or ``list``\n:warning: Use of this operation is discouraged. It prevents some optimisations from operating.")
    .def("getTupleForDataPoint",&escript::Data::getValueOfDataPointAsTuple,args("dataPointNo"),
":return: Value of the specified datapoint\n:rtype: ``tuple``\n:param dataPointNo: datapoint to access\n:type dataPointNo: ``int``")
    .def("getTupleForGlobalDataPoint",&escript::Data::getValueOfGlobalDataPointAsTuple,args("procNo","dataPointNo"),"Get a specific datapoint from a specific process\n\n"
":rtype: ``tuple``\n:param procNo: MPI rank of the process\n:type procNo: positive ``int``"
"\n:param dataPointNo: datapoint to access\n:type dataPointNo: int")
    .def("setToZero",&escript::Data::setToZero,"After this call the object will store values of the same shape as before but all components will be zero.")
    .def("interpolate",&escript::Data::interpolate,args("functionspace"),"Interpolate this object's values into a new functionspace.")
    .def("interpolateTable", &escript::Data::interpolateFromTable2DP, 
(arg("table"),arg("Amin"),arg("Astep"), arg("B"), arg("Bmin"), arg("Bstep"), arg("undef")=1.e50),
"Creates a new Data object by interpolating using the source data (which are\n"
"looked up in ``table``)\n``A`` must be the outer dimension on the table\n\n"
":param table: two dimensional collection of values\n"
":param Amin: The base of locations in table\n:type Amin: float\n"
":param Astep: size of gap between each item in the table\n:type Astep: float\n"
":param undef: upper bound on interpolated values\n:type undef: float\n"
":param B: Scalar representing the second coordinate to be mapped into the table\n"
":type B: `Data`\n"
":param Bmin: The base of locations in table for 2nd dimension\n:type Bmin: float\n"
":param Bstep: size of gap between each item in the table for 2nd dimension\n:type Bstep: float\n"
":raise RuntimeError(DataException): if the cordinates do not map into the table or if the interpolated value is above ``undef``"
"\n:rtype: `Data`"
)
    .def("interpolateTable", &escript::Data::interpolateFromTable1DP, 
(arg("table"),arg("Amin"),arg("Astep"), arg("undef")=1.e50),
"Creates a new Data object by interpolating using the source data (which are\n"
"looked up in ``table``)\n\n"
":param table: one dimensional collection of values\n"
":param Amin: The base of locations in table\n:type Amin: float\n"
":param Astep: size of gap between each item in the table\n:type Astep: float\n"
":param undef: upper bound on interpolated values\n:type undef: float\n"
":raise RuntimeError(DataException): if the cordinates do not map into the table or if the interpolated value is above ``undef``"
"\n:rtype: `Data`"
)


    .def("minGlobalDataPoint",&escript::Data::minGlobalDataPoint)
    .def("maxGlobalDataPoint",&escript::Data::maxGlobalDataPoint)
    .def("saveDX",&escript::Data::saveDX,args("fileName"),"Save the object in DX format.\n\n"
":param fileName: filename\n:type fileName: string")
    .def("saveVTK",&escript::Data::saveVTK, args("fileName"),"Save the object in VTK format.\n\n"
":param fileName: \n:type fileName: string")
    .def("getTagNumber",&escript::Data::getTagNumber,args("dpno"),"Return tag number for the specified datapoint\n\n:rtype: int\n:param dpno: datapoint number\n:type dpno: int")
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
    .def("_Lsup",&escript::Data::Lsup,":return: the Lsup-norm of the object\n:rtype: float\n:note: If the ``Data`` contains no values, zero will be returned instead.")
    .def("_sup",&escript::Data::sup,":return: the maximum value over all data points.\n:rtype: float\n:note: If the ``Data`` contains no values a large negative value will be returned instead.")
    .def("_inf",&escript::Data::inf,":return: minimum value over all components and all data points\n:rtype: float\n:note: If the ``Data`` contains no values a large positive value will be returned instead.")
    .def("_integrateToTuple",&escript::Data::integrateToTuple,":return: Calculate the integral over the function space domain as a python tuple\n:rtype: tuple")
    // following implements the python abs operator
    .def("__abs__",&escript::Data::abs,":return: absolute value\n\n:rtype: `Data`")
    // following implements the python "-" negation operator
    .def("__neg__",&escript::Data::neg, ":return: negation of the values in this object\n:rtype: `Data`")
    // following implements the python "+" identity operator
    .def("__pos__",&escript::Data::pos, "\nThe unary + operator\n\n:rtype: `Data`")
    // following two functions implement the python [] operator
    .def("__getitem__",&escript::Data::getItem,"Used by the python [] operator\n\n:rtype: `Data`")
    .def("__setitem__",&escript::Data::setItemO,"Used by the python [] operator")
    .def("__setitem__",&escript::Data::setItemD,"Used by the python [] operator")
    // following two functions implement the python ** operator
    .def("__pow__",&escript::Data::powO,"Used by the python ** operator\n\n:rtype: `Data`")
    .def("__pow__",&escript::Data::powD)
    .def("__rpow__",&escript::Data::rpowO,"\nUsed by the python ** operator\n\n:rtype: `Data`")
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
  def("ContinuousFunction",escript::continuousFunction,args("domain"),
":return: a continuous FunctionSpace (overlapped node values)\n"
":rtype: `FunctionSpace`");
  def("ReducedContinuousFunction",escript::reducedContinuousFunction,args("domain"),
":return: a continuous with reduced order FunctionSpace (overlapped node values on reduced element order)\n"
":rtype: `FunctionSpace`");
  def("Function",escript::function,args("domain"),":return: a function `FunctionSpace`\n"
":rtype: `FunctionSpace`");
  def("ReducedFunction",escript::reducedFunction, args("domain"),":return: a function FunctionSpace with reduced integration order\n:rtype: `FunctionSpace`");
  def("FunctionOnBoundary",escript::functionOnBoundary, args("domain"), ":return: a function on boundary FunctionSpace\n:rtype: `FunctionSpace`");
  def("ReducedFunctionOnBoundary",escript::reducedFunctionOnBoundary, args("domain"), 
":return: a function on boundary FunctionSpace with reduced integration order\n"
":rtype: `FunctionSpace`");
  def("FunctionOnContactZero",escript::functionOnContactZero, args("domain"), ":return: Return a FunctionSpace on left side of contact\n:rtype: `FunctionSpace`");
  def("ReducedFunctionOnContactZero",escript::reducedFunctionOnContactZero, args("domain"),
 ":return: a FunctionSpace  on left side of contact with reduced integration order\n:rtype: `FunctionSpace`");
  def("FunctionOnContactOne",escript::functionOnContactOne, args("domain"), ":return: Return a FunctionSpace on right side of contact\n:rtype: `FunctionSpace`");
  def("ReducedFunctionOnContactOne",escript::reducedFunctionOnContactOne, args("domain"),
 ":return: Return a FunctionSpace on right side of contact with reduced integration order\n"
":rtype: `FunctionSpace`");
  def("Solution",escript::solution, args("domain"), ":rtype: `FunctionSpace`");
  def("ReducedSolution",escript::reducedSolution, args("domain"), ":rtype: `FunctionSpace`");
  def("DiracDeltaFunction",escript::diracDeltaFunction, args("domain"), ":rtype: `FunctionSpace`");





  //
  // Factory methods for Data
  //
  def("load",escript::load, args("fileName","domain"), "reads Data on domain from file in netCDF format\n\n:param fileName:\n:type fileName: ``string``\n:param domain:\n:type domain: `Domain`");
  def("loadIsConfigured",escript::loadConfigured,":return: True if the load function is configured.");
  def("Scalar",escript::Scalar,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false),
"Construct a Data object containing scalar data-points.\n\n:param value: scalar value for all points\n"
"\n:rtype: `Data`\n:type value: float\n:param what: FunctionSpace for Data\n:type what: `FunctionSpace`\n"
":param expanded: If True, a value is stored for each point. If False, more efficient representations may be used\n"
":type expanded: ``bool``");
  def("Vector",escript::Vector,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false),
"Construct a Data object containing rank1 data-points.\n\n:param value: scalar value for all points\n"
"\n:rtype: `Data`\n:type value: float\n:param what: FunctionSpace for Data\n:type what: `FunctionSpace`\n"
":param expanded: If True, a value is stored for each point. If False, more efficient representations may be used\n"
":type expanded: ``bool``");
  def("Tensor",escript::Tensor,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false),
"Construct a Data object containing rank2 data-points.\n\n:param value: scalar value for all points\n"
"\n:rtype: `Data`\n:type value: float\n:param what: FunctionSpace for Data\n:type what: `FunctionSpace`\n"
":param expanded: If True, a value is stored for each point. If False, more efficient representations may be used\n"
":type expanded: ``bool``");
  def("Tensor3",escript::Tensor3,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false),
"Construct a Data object containing rank3 data-points.\n\n:param value: scalar value for all points\n"
"\n:rtype: `Data`\n:type value: float\n:param what: FunctionSpace for Data\n:type what: `FunctionSpace`\n"
":param expanded: If True, a value is stored for each point. If False, more efficient representations may be used\n"
":type expanded: ``bool``"
);
  def("Tensor4",escript::Tensor4,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false),
"Construct a Data object containing rank4 data-points.\n\n:param value: scalar value for all points\n"
"\n:rtype: `Data`\n:type value: float\n:param what: FunctionSpace for Data\n:type what: `FunctionSpace`\n"
":param expanded: If True, a value is stored for each point. If False, more efficient representations may be used\n"
":type expanded: ``bool``"
);

  //
  // Binary operators
  //
  def("C_GeneralTensorProduct",escript::C_GeneralTensorProduct,
      (arg("arg0"),
       arg("arg1"),
       arg("axis_offset")=0,
       arg("transpose")=0),
"Compute a tensor product of two Data objects.\n\n:rtype: `Data`\n:param arg0:\n"
":param arg1:\n:param axis_offset:\n:type axis_offset: ``int``\n"
":param transpose: 0: transpose neither, 1: transpose arg0, 2: transpose arg1\n"
":type transpose: int");

  //
  // Interface for AbstractSystemMatrix
  //
  class_<escript::AbstractSystemMatrix>("Operator","",init<>())    // Doco goes in the empty string param
     .def("isEmpty",&escript::AbstractSystemMatrix::isEmpty,":rtype: ``bool``\n"
":return: True if matrix is empty")
     .def("solve",&escript::AbstractSystemMatrix::solve, args("in","options"),
":return: the solution *u* of the linear system *this*u=in*\n\n:param in:\n:type in: `Data`")
     .def("of",&escript::AbstractSystemMatrix::vectorMultiply,args("right"),
"matrix*vector multiplication")
     .def("saveMM",&escript::AbstractSystemMatrix::saveMM, args("fileName"), 
"writes the matrix to a file using the Matrix Market file format")
     .def("saveHB",&escript::AbstractSystemMatrix::saveHB, args("filename"),
"writes the matrix to a file using the Harwell-Boeing file format")
     .def("resetValues",&escript::AbstractSystemMatrix::resetValues, "resets the matrix entries")
     .def(self*other<escript::Data>());
  //
  // Interface for AbstractTransportProblem
  //
  class_<escript::AbstractTransportProblem>("TransportProblem","",init<>())    // Doco goes in the empty string param
     .def("isEmpty",&escript::AbstractTransportProblem::isEmpty,":rtype: ``int``")
     .def("solve",&escript::AbstractTransportProblem::solve, args("source","dt", "options"),
"returns the solution *u* for a time step *dt>0*\n\n:rtype: `Data`\n"
":param source:\n:type source: `Data`")
     .def("setInitialValue",&escript::AbstractTransportProblem::setInitialValue, args("u"),
"sets the value for *u* at time *t=0*.")
     .def("insertConstraint",&escript::AbstractTransportProblem::insertConstraint,
args("source", "q", "r","factor"),
"inserts constraint *u_{,t}=r* where *q>0*  into the problem using a weighting factor")
     .def("reset",&escript::AbstractTransportProblem::resetTransport,
"resets the transport operator typically as they have been updated.")
     .def("resetValues",&escript::AbstractTransportProblem::resetTransport)
     .def("getSafeTimeStepSize",&escript::AbstractTransportProblem::getSafeTimeStepSize)
     .def("getUnlimitedTimeStepSize",&escript::AbstractTransportProblem::getUnlimitedTimeStepSize);

  // Functions to modify global parameters
  def("setEscriptParamInt",escript::setEscriptParamInt,
      (arg("name"), arg("value")=0), "Modify the value of an escript tuning parameter\n\n"
":param name:\n:type name: ``string``\n:param value:\n:type value: ``int``");
  def("getEscriptParamInt",escript::getEscriptParamInt,
      (arg("name"),arg("sentinel")=0), "Read the value of an escript tuning paramter\n\n"
":param name: parameter to lookup\n:type name: ``string``\n:param sentinel: Value to be returned if ``name`` is not a known parameter\n"
":type sentinel: ``int``");
  def("listEscriptParams",escript::listEscriptParams,":return: A list of pairs (p,d) where p is the name of a parameter for escript and d is a description.");

  //
  // Register esysExceptionTranslator
  //
  register_exception_translator<esysUtils::EsysException>(&esysUtils::esysExceptionTranslator);
}
