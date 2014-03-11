
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
#include "esysUtils/Esys_MPI.h"
#include "EscriptParams.h"
#include "TestDomain.h"
#include "SubWorld.h"
#include "SplitWorld.h"
#include "Crates.h"

#include "esysUtils/blocktimer.h"

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
 * \version 3.3.1
 *
 * Main modules/namespaces:
 *
 * - \ref escript
 *
 * - \ref paso
 *
 * - \ref finley
 *
 * - \ref dudley
 *
 * - \ref ripley
 *
 * - \ref weipa
 *
 * Depending on your system, only one of the following will work:
 *
 * - <a href="../../sphinx_api/index.html">Python module documentation (sphinx generated)</a>
 *
 * - <a href="../../epydoc/index.html">Python module documentation (epydoc generated)</a>
 *
 */

/*
namespace escript
{
  
  // Note: not virtual because it calls the virtual probeInterpolationOnDomain
  ESCRIPT_DLL_API
  bool canInterpolate(FunctionSpace src, FunctionSpace dest)
  {
      return src.getDomain()->probeInterpolationOnDomain(src.getTypeCode(), dest.getTypeCode());
  }  
  
  
}
*/

#include <boost/python/raw_function.hpp>


bool test1(double d)
{
    std::cout << "Line " << __LINE__ << std::endl;
    return true;
}


bool test2(boost::python::object o, double x, double y, double z)
{
    std::cout << "Line " << __LINE__ << std::endl;
    return true;
}

bool test3(double x, double y, double z)
{
    std::cout << "Line " << __LINE__ << std::endl;
    return true;
}

bool test4(double x, double y, double z, bool bozo1=true, bool bozo2=false)
{
    std::cout << "Line " << __LINE__ << std::endl;
    return true;
}


namespace
{

object raw1(boost::python::tuple t, boost::python::dict kwargs)
{
    std::cout << "In raw1\n";
    if (len(t)<2)
    {
        return object(false);
    }
    t[0](*t, **kwargs);
    return object();
}

object raw2(boost::python::tuple t, boost::python::dict kwargs)
{
    std::cout << "In raw2\n";
//    raw1(t, kwargs);
    if (len(t)<1)
    {
        return object(false);
    }
    object target=t[0];
    object zz=t.slice(1,len(t));
    tuple t2=tuple(zz);
    target(*t2, **kwargs);
    //target(t2[0], t2[1], t2[2], t2[3], **kwargs);   // This doesn't work
    return object();
}


}



tuple raw(tuple args, dict kw)
{
    return make_tuple(args, kw);
}

BOOST_PYTHON_MODULE(escriptcpp)
{

  #if BOOST_VERSION >= 103500
// params are: bool show_user_defined, bool show_py_signatures, bool show_cpp_signatures
  docstring_options docopt(true,true,false);
  #endif

  scope().attr("__doc__") = "To use this module, please import esys.escript";      
  
  
    def("jf1", &test1, "Here is a docstring");
  def("jf2", &test2);
  def("jf3", &test3);
  def("jf4", &test4, args("x","y","z","bozo1","bozo2"));

  def("passthrough", raw_function(raw1));
  def("pass2", raw_function(raw2,2));

  
  
  class_<escript::AbstractCrate, escript::crate_ptr, boost::noncopyable>("AbstractCrate", "", no_init);
  class_<escript::DataCrate, bases<escript::AbstractCrate> >("DataCrate","Stores and merges Data objects for transport between worlds.",
    init<std::string, std::string>(args("label", "operation"))
  );

/* begin SubWorld things */
  // Why doesn't this have a doc-string?   Because it doesn't compile if you try to add one  
  def("buildDomains", raw_function(escript::raw_buildDomains,2));
      
  class_<escript::SplitWorld, boost::noncopyable>("SplitWorld", "Manages a group of sub worlds", init<unsigned int>(args("num_worlds")))
    .def("registerCrate", &escript::SplitWorld::registerCrate, arg("crate")) 
    .def("runJobs", &escript::SplitWorld::runJobs, arg("tuplelist"), "Create a set of jobs and execute them on the subworlds.");

  // This class has no methods. This is deliberate - at this stage, I would like this to be an opaque type  
  class_ <escript::SubWorld, escript::SubWorld_ptr, boost::noncopyable>("SubWorld", "Information about a group of workers.", no_init);
/* end SubWorld things */
  
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
  def("runMPIProgram",escript::runMPIProgram,"Spawns an external MPI program using a separate communicator.");
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
"That is, items in each column will be printed in the same order.\n"
"So you can be sure that values in the same row correspond to the same input value.\n\n"
"\n:param filename:\n:type filename: ``string``\n"
":param arg: dictionary of named `Data` objects. If one is called ``mask`` it must be scalar data.\n"
":type arg: ``dict``\n"
":param sep: separator for columns (defaults to ',')\n"
":type sep: ``string``\n"
":param csep: separator for fields within data object (defaults to \"_\")\n:type csep: ``string``\n"
":param append: If True, write to the end of ``filename``\n:type append: ``string``\n");
   def("canInterpolate", &escript::canInterpolate, args("src", "dest"),":param src: Source FunctionSpace\n:param dest: Destination FunctionSpace\n:return: True if src can be interpolated to dest\n:rtype: `bool`");

  //
  // Interface for AbstractDomain
  //
  class_<escript::AbstractDomain, escript::Domain_ptr, boost::noncopyable >("Domain","Base class for all domains.",no_init)
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
     .def("dump",&escript::AbstractDomain::dump,args("filename"),"Dumps the domain to a file"
":param filename:\n:type filename: string")
     .def("getMPISize",&escript::AbstractDomain::getMPISize,":return: the number of processes used for this `Domain`\n:rtype: ``int``")
     .def("getMPIRank",&escript::AbstractDomain::getMPIRank,":return: the rank of this process\n:rtype: ``int``")
     .def("MPIBarrier",&escript::AbstractDomain::MPIBarrier,"Wait until all processes have reached this point")
     .def("onMasterProcessor",&escript::AbstractDomain::onMasterProcessor,":return: True if this code is executing on the master process\n:rtype: `bool`")
     .def("supportsContactElements", &escript::AbstractDomain::supportsContactElements,"Does this domain support contact elements.")
     .def(self == self)
     .def(self != self);

  //
  // Interface for AbstractContinuousDomain
  //
  class_<escript::AbstractContinuousDomain, bases<escript::AbstractDomain>, boost::noncopyable >("ContinuousDomain","Class representing continuous domains",no_init)
       .def("getSystemMatrixTypeId",&escript::AbstractContinuousDomain::getSystemMatrixTypeId,
args("solver", "preconditioner", "package", "symmetry"),
":return: the identifier of the matrix type to be used for the global stiffness matrix "
"when a particular solver package and symmetric matrix is used.\n"
":rtype: int")
       .def("getTransportTypeId",&escript::AbstractContinuousDomain::getTransportTypeId,
args("solver", "preconditioner", "package", "symmetry"))

      .def("addPDEToSystem",&escript::AbstractContinuousDomain::addPDEToSystem,
args("mat", "rhs","A", "B", "C", "D", "X", "Y", "d", "y", "d_contact", "y_contact", "d_dirac", "y_dirac"),
"adds a PDE onto the stiffness matrix mat and a rhs\n\n"
":param mat:\n:type mat: `OperatorAdapter`\n:param rhs:\n:type rhs: `Data`\n"
":param A:\n:type A: `Data`\n"
":param B:\n:type B: `Data`\n"
":param C:\n:type C: `Data`\n"
":param D:\n:type D: `Data`\n"
":param X:\n:type X: `Data`\n"
":param Y:\n:type Y: `Data`\n"
":param d:\n:type d: `Data`\n"
":param d_contact:\n:type d_contact: `Data`\n"
":param y_contact:\n:type y_contact: `Data`\n"
":param d_dirac:\n:type d_dirac: `Data`\n"
":param y_dirac:\n:type y_dirac: `Data`\n"
)
      .def("addPDEToRHS",&escript::AbstractContinuousDomain::addPDEToRHS, 
args("rhs", "X", "Y", "y", "y_contact", "y_dirac"),
"adds a PDE onto the stiffness matrix mat and a rhs\n\n"
":param rhs:\n:type rhs: `Data`\n"
":param X:\n:type X: `Data`\n"
":param Y:\n:type Y: `Data`\n"
":param y:\n:type y: `Data`\n"
":param y_contact:\n:type y_contact: `Data`\n"
":param y_dirac:\n:type y_dirac: `Data`"
)
      .def("addPDEToTransportProblem",&escript::AbstractContinuousDomain::addPDEToTransportProblem,
args( "tp", "source", "M", "A", "B", "C", "D", "X", "Y", "d", "y", "d_contact", "y_contact", "d_dirac", "y_dirac"),
":param tp:\n:type tp: `TransportProblemAdapter`\n"
":param source:\n:type source: `Data`\n"
":param M:\n:type M: `Data`\n"
":param A:\n:type A: `Data`\n"
":param B:\n:type B: `Data`\n"
":param C:\n:type C: `Data`\n"
":param D:\n:type D: `Data`\n"
":param X:\n:type X: `Data`\n"
":param Y:\n:type Y: `Data`\n"
":param d:\n:type d: `Data`\n"
":param y:\n:type y: `Data`\n"
":param d_contact:\n:type d_contact: `Data`\n"
":param y_contact:\n:type y_contact: `Data`\n"
":param d_dirac:\n:type d_dirac: `Data`\n"
":param y_dirac:\n:type y_dirac: `Data`\n"
)
      .def("newOperator",&escript::AbstractContinuousDomain::newSystemMatrix,
args("row_blocksize", "row_functionspace", "column_blocksize", "column_functionspace", "type"),
"creates a SystemMatrixAdapter stiffness matrix and initializes it with zeros\n\n"
":param row_blocksize:\n:type row_blocksize: ``int``\n"
":param row_functionspace:\n:type row_functionspace: `FunctionSpace`\n"
":param column_blocksize:\n:type column_blocksize: ``int``\n"
":param column_functionspace:\n:type column_functionspace: `FunctionSpace`\n"
":param type:\n:type type: ``int``\n"
)
      .def("newTransportProblem",&escript::AbstractContinuousDomain::newTransportProblem,
args("theta", "blocksize", "functionspace", "type"),
"creates a TransportProblemAdapter\n\n"
":param theta:\n:type theta: ``float``\n"
":param blocksize:\n:type blocksize: ``int``\n"
":param functionspace:\n:type functionspace: `FunctionSpace`\n"
":param type:\n:type type: ``int``\n"
)
      .def("getDataShape",&escript::AbstractContinuousDomain::getDataShape, args("functionSpaceCode"),
":return: a pair (dps, ns) where dps=the number of data points per sample, and ns=the number of samples\n:rtype: ``tuple``")
      .def("print_mesh_info",&escript::AbstractContinuousDomain::Print_Mesh_Info,(arg("full")=false),
":param full:\n:type full: ``bool``")
      .def("getDescription",&escript::AbstractContinuousDomain::getDescription,
":return: a description for this domain\n:rtype: ``string``")
      .def("setX",&escript::AbstractContinuousDomain::setNewX,
args("arg"), "assigns new location to the domain\n\n:param arg:\n:type arg: `Data`")
      .def("getNumDataPointsGlobal",&escript::AbstractContinuousDomain::getNumDataPointsGlobal,
":return: the number of data points summed across all MPI processes\n"
":rtype: ``int``");




  //
  // Interface for TestDomain
  //
  class_ <escript::TestDomain, bases<escript::AbstractDomain> >("TestDomain", 
	"Test Class for domains with no structure. May be removed from future releases without notice.", no_init);

  // This is the only python visible way to get a TestDomain
  def("getTestDomainFunctionSpace",&escript::getTestDomainFunctionSpace, (arg("dpps"), 
 arg("samples"), arg("size")=1),
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
"Set tags according to a mask\n\n:param newtag: tag number to set\n:type newtag: string, non-zero ``int``\n:param mask: Samples which correspond to positive values in the mask will be set to ``newtag``.\n:type mask: scalar `Data`");
  fs_definer.def("setTags",&escript::FunctionSpace::setTagsByString,args("newtag","mask"));
  fs_definer.def("getTagFromDataPointNo",
                 &escript::FunctionSpace::getTagFromDataPointNo,":return: the tag associated with the given sample number.\n:rtype: int");
  fs_definer.def("getReferenceIDFromDataPointNo", &escript::FunctionSpace::getReferenceIDFromDataPointNo,args("dataPointNo"),":return: the reference number associated with ``dataPointNo``\n:rtype: int ");
  fs_definer.def("getListOfTags",&escript::FunctionSpace::getListOfTags,":return: a list of the tags used in this function space\n:rtype: ``list``");
  fs_definer.def("getApproximationOrder", &escript::FunctionSpace::getApproximationOrder,":return: the approximation order referring to the maximum degree of a polynomial which can be represented exactly in interpolation and/or integration.\n:rtype: ``int``");
  fs_definer.def("__str__", &escript::FunctionSpace::toString);
  fs_definer.def(self == self);
  fs_definer.def(self != self);
  //
  // Interface for Data
  //
  class_<escript::Data>("Data"/*,shared_ptr<Data>*/, "Represents a collection of datapoints. It is used to store the values of a function. For more details please consult the c++ class documentation.",init<>() )
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
    .def("getRank",&escript::Data::getDataPointRank,":return: the number of indices required to address a component of a datapoint\n:rtype: positive ``int``")
    .def("dump",&escript::Data::dump,args("fileName"),"Save the data as a netCDF file\n\n:param fileName: \n:type fileName: ``string``")
    .def("toListOfTuples",&escript::Data::toListOfTuples, (arg("scalarastuple")=false),
"Return the datapoints of this object in a list. Each datapoint is stored as a tuple.\n\n"
":param scalarastuple: if True, scalar data will be wrapped as a tuple."
" True => [(0), (1), (2)]; False => [0, 1, 2]")
    .def("copyWithMask",&escript::Data::copyWithMask,args("other","mask"),
"Selectively copy values from ``other`` `Data`."
"Datapoints which correspond to positive values in ``mask`` will be copied from ``other``\n"
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
    .def("_setTupleForGlobalDataPoint", &escript::Data::setTupleForGlobalDataPoint)
    .def("setValueOfDataPoint",&escript::Data::setValueOfDataPoint,"\nModify the value of a single datapoint.\n\n:param dataPointNo:\n"
":type dataPointNo: int\n:param value: \n:type value: ``float`` or an object which acts like an array, ``tuple`` or ``list``\n:warning: Use of this operation is discouraged. It prevents some optimisations from operating.")
    .def("getTupleForDataPoint",&escript::Data::getValueOfDataPointAsTuple,args("dataPointNo"),
":return: Value of the specified datapoint\n:rtype: ``tuple``\n:param dataPointNo: datapoint to access\n:type dataPointNo: ``int``")
    .def("getTupleForGlobalDataPoint",&escript::Data::getValueOfGlobalDataPointAsTuple,args("procNo","dataPointNo"),"Get a specific datapoint from a specific process\n\n"
":rtype: ``tuple``\n:param procNo: MPI rank of the process\n:type procNo: positive ``int``"
"\n:param dataPointNo: datapoint to access\n:type dataPointNo: int")
    .def("setToZero",&escript::Data::setToZero,"After this call the object will store values of the same shape as before but all components will be zero.")
    .def("interpolate",&escript::Data::interpolate,args("functionspace"),"Interpolate this object's values into a new functionspace.")
    .def("_interpolateTable3d", &escript::Data::interpolateFromTable3DP, 
(arg("table"),arg("Amin"),arg("Astep"), arg("B"), arg("Bmin"), arg("Bstep"), arg("C"), arg("Cmin"), arg("Cstep"), arg("undef")=1.e50, arg("check_boundaries")=false, "For internal use only. Please use the interpolateTable function.")
)

    .def("interpolateTable", &escript::Data::interpolateFromTable2DP, 
(arg("table"),arg("Amin"),arg("Astep"), arg("B"), arg("Bmin"), arg("Bstep"), arg("undef")=1.e50, arg("check_boundaries")=false),
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
":param check_boundaries: if true, then values outside the boundaries will be rejected. If false, then boundary values will be used.\n"
":raise RuntimeError(DataException): if the coordinates do not map into the table or if the interpolated value is above ``undef``"
"\n:rtype: `Data`"
)
    .def("interpolateTable", &escript::Data::interpolateFromTable1DP, 
(arg("table"),arg("Amin"),arg("Astep"), arg("undef")=1.e50, arg("check_boundaries")=false)/*,
"Creates a new Data object by interpolating using the source data (which are\n"
"looked up in ``table``)\n\n"
":param table: one dimensional collection of values\n"
":param Amin: The base of locations in table\n:type Amin: float\n"
":param Astep: size of gap between each item in the table\n:type Astep: float\n"
":param undef: upper bound on interpolated values\n:type undef: float\n"
":param check_boundaries: if true, then values outside the boundaries will be rejected. If false, then boundary values will be used.\n"
":raise RuntimeError(DataException): if the coordinates do not map into the table or if the interpolated value is above ``undef``"
"\n:rtype: `Data`"
*/
)
    .def("nonuniformInterpolate", &escript::Data::nonuniforminterp, "1D interpolation with non equally spaced points",
      (arg("in"), arg("out"), arg("check_boundaries")),
"Creates a Data object by linear interpolation of the function F(in)->out\n\n"
":param in: input values of interpolation function\n"
":param out: corresponding output values of interpolation function\n"
":param check_boundaries: If True, an exception will the thrown if the data object contains values"
"outside the range given by ``in``.\n"
    )
    .def("nonuniformSlope", &escript::Data::nonuniformslope, "1D interpolation of slope with non equally spaced points",
      (arg("in"), arg("out"), arg("check_boundaries")),
"Creates a Data object by computing the slope of the function F(in)->out\n\n"
":param in: input values of interpolation function\n"
":param out: corresponding output values of interpolation function\n"
":param check_boundaries: If True, an exception will the thrown if the data object contains values"
"outside the range given by ``in``.\n"
    )	 
    .def("minGlobalDataPoint",&escript::Data::minGlobalDataPoint,"Please consider using getInfLocator() from pdetools instead.")
    .def("maxGlobalDataPoint",&escript::Data::maxGlobalDataPoint, "Please consider using getSupLocator() from pdetools instead.")
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
    // following three functions implement the python [] operator
    .def("__getitem__",&escript::Data::getItem,"Used by the python [] operator\n\n:rtype: `Data`")
    .def("__setitem__",&escript::Data::setItemO,"Used by the python [] operator")
    .def("__setitem__",&escript::Data::setItemD,"Used by the python [] operator")
    // following three functions implement the python ** operator
    .def("__pow__",&escript::Data::powO,"Used by the python ** operator\n\n:rtype: `Data`")
    .def("__pow__",&escript::Data::powD)
    .def("__rpow__",&escript::Data::rpowO,"\nUsed by the python ** operator\n\n:rtype: `Data`")
    // following two functions implement the newer python / operator
    .def("__truediv__",&escript::Data::truedivO)
    .def("__truediv__",&escript::Data::truedivD)
    .def("__rtruediv__",&escript::Data::rtruedivO)
    // NOTE:: The order of these declarations is important. Anything
    // declared before the generic declaration isn't found so the generic
    // version will be called. 
//    .def(self + other<object>())
//    .def(other<object>() + self)
//    .def(self + self)
    .def(self += other<object>())
    .def(self += self)

//     .def(self - other<object>())
//     .def(other<object>() - self)
//     .def(self - self)
    .def(self -= other<object>())
    .def(self -= self)

//     .def(self * other<object>())
//     .def(other<object>() * self)
//     .def(self * self)
    .def(self *= other<object>())
    .def(self *= self)

//     .def(self / other<object>())
//     .def(other<object>() / self)
//     .def(self / self)
    .def(self /= other<object>())
    .def(self /= self)
    // Need scope resolution due to a bug either in the compiler or
    // the boost code. This calls operator << for Data.
    .def(self_ns::str(self))
    .def("_inverse", &escript::Data::matrixInverse, ":return: inverse of square matrices\n")
//    .def("__add__", &escript::Data::addOperatorD)
    .def("__add__", &escript::Data::__add__)
    .def("__radd__", &escript::Data::__add__)  // its the same coz + is commutative
    .def("__sub__", &escript::Data::__sub__)
    .def("__rsub__", &escript::Data::__rsub__)
    .def("__mul__", &escript::Data::__mul__)   
    .def("__rmul__", &escript::Data::__mul__)   // commutative
    .def("__div__", &escript::Data::__div__)   
    .def("__rdiv__", &escript::Data::__rdiv__)   // commutative
    
    ;

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
  def("DiracDeltaFunctions",escript::diracDeltaFunctions, args("domain"), ":rtype: `FunctionSpace`");





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
 def("Vector", escript::VectorFromObj,
      (arg("value"),
	arg("what")=escript::FunctionSpace(),
	arg("expanded")=false));
  def("Tensor",escript::Tensor,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false),
"Construct a Data object containing rank2 data-points.\n\n:param value: scalar value for all points\n"
"\n:rtype: `Data`\n:type value: float\n:param what: FunctionSpace for Data\n:type what: `FunctionSpace`\n"
":param expanded: If True, a value is stored for each point. If False, more efficient representations may be used\n"
":type expanded: ``bool``");
 def("Tensor", escript::TensorFromObj,
      (arg("value"),
	arg("what")=escript::FunctionSpace(),
	arg("expanded")=false));
  def("Tensor3",escript::Tensor3,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false),
"Construct a Data object containing rank3 data-points.\n\n:param value: scalar value for all points\n"
"\n:rtype: `Data`\n:type value: float\n:param what: FunctionSpace for Data\n:type what: `FunctionSpace`\n"
":param expanded: If True, a value is stored for each point. If False, more efficient representations may be used\n"
":type expanded: ``bool``"
);
 def("Tensor3", escript::Tensor3FromObj,
      (arg("value"),
	arg("what")=escript::FunctionSpace(),
	arg("expanded")=false));
  def("Tensor4",escript::Tensor4,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false),
"Construct a Data object containing rank4 data-points.\n\n:param value: scalar value for all points\n"
"\n:rtype: `Data`\n:type value: float\n:param what: FunctionSpace for Data\n:type what: `FunctionSpace`\n"
":param expanded: If True, a value is stored for each point. If False, more efficient representations may be used\n"
":type expanded: ``bool``"
);
 def("Tensor4", escript::Tensor4FromObj,
      (arg("value"),
	arg("what")=escript::FunctionSpace(),
	arg("expanded")=false));


 def("RandomData", escript::randomData, (arg("shape"), arg("fs"), arg("seed")=0, arg("filter")=boost::python::tuple()),
"Creates a new expanded Data object containing pseudo-random values.\n\n"
":param shape: datapoint shape\n:type shape: tuple\n"
":param fs: function space for data object.\n:type fs: `FunctionSpace`\n"
":param seed: seed for random number generator.\n:type seed: long\n");

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
  class_<escript::AbstractSystemMatrix,escript::ASM_ptr, boost::noncopyable>("Operator","",init<>())    // Doco goes in the empty string param
     .def("isEmpty",&escript::AbstractSystemMatrix::isEmpty,":rtype: ``bool``\n"
":return: True if matrix is empty")
     .def("solve",&escript::AbstractSystemMatrix::solve, args("in","options"),
":return: the solution *u* of the linear system *this*u=in*\n\n:param in:\n:type in: `Data`")
     .def("of",&escript::AbstractSystemMatrix::vectorMultiply,args("right"),
"matrix*vector multiplication")
     .def("nullifyRowsAndCols",&escript::AbstractSystemMatrix::nullifyRowsAndCols)
     .def("saveMM",&escript::AbstractSystemMatrix::saveMM, args("fileName"), 
"writes the matrix to a file using the Matrix Market file format")
     .def("saveHB",&escript::AbstractSystemMatrix::saveHB, args("filename"),
"writes the matrix to a file using the Harwell-Boeing file format")
     .def("resetValues",&escript::AbstractSystemMatrix::resetValues, "resets the matrix entries")   
     .def(self*other<escript::Data>());
  //
  // Interface for AbstractTransportProblem
  //
  class_<escript::AbstractTransportProblem, escript::ATP_ptr, boost::noncopyable >("TransportProblem","",init<>())    // Doco goes in the empty string param
     .def("isEmpty",&escript::AbstractTransportProblem::isEmpty,":rtype: ``int``")
     .def("solve",&escript::AbstractTransportProblem::solve, args("u0","source","dt", "options"),
"returns the solution *u* for a time step *dt>0* with initial value u0\n\n:rtype: `Data`\n"
":param source:\n:type source: `Data`")
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
      (arg("name"),arg("sentinel")=0), "Read the value of an escript tuning parameter\n\n"
":param name: parameter to lookup\n:type name: ``string``\n:param sentinel: Value to be returned if ``name`` is not a known parameter\n"
":type sentinel: ``int``");
  def("listEscriptParams",escript::listEscriptParams,":return: A list of pairs (p,d) where p is the name of a parameter for escript and d is a description.");


  def("resolveGroup", escript::resolveGroup);

#ifdef IKNOWWHATIMDOING

  def("applyBinaryCFunction", escript::applyBinaryCFunction, (arg("function"), arg("outshape"),
arg("in1"), 
arg("in2"))
);
#endif

  def("_condEval", escript::condEval, (arg("mask"), arg("trueval"), arg("falseval")));

  //
  // Register esysExceptionTranslator
  //
  register_exception_translator<esysUtils::EsysException>(&esysUtils::esysExceptionTranslator);
}
