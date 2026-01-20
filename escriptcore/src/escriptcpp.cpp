
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

#include "AbstractContinuousDomain.h"
#include "AbstractSystemMatrix.h"
#include "AbstractTransportProblem.h"
#include "Data.h"
#include "DataFactory.h"
#include "DataVector.h"
#include "EscriptParams.h"
#include "ExceptionTranslators.h"
#include "FunctionSpace.h"
#include "FunctionSpaceFactory.h"
#include "SolverOptions.h"
#include "TestDomain.h"
#include "Utils.h"

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/module.hpp>
#include <boost/python/object.hpp>
#include <boost/python/tuple.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/version.hpp>

#ifdef ESYS_HAVE_BOOST_NUMPY
#include <boost/python/numpy.hpp>
#endif

using namespace boost::python;

/*! \mainpage Esys Documentation
 *
 * \version 4.2
 *
 * Main modules/namespaces:
 *
 * - \ref escript
 *
 * - \ref paso
 *
 * - \ref finley
 *
 * - \ref ripley
 *
 * - \ref speckley
 *
 * - \ref weipa
 *
 * This documentation describes the C++ layer of escript and related libraries.
 * For documentation of the python API, please see:
 * <a href="../../sphinx_api/index.html">Here</a>
 *
 */

#include <boost/python/raw_function.hpp>

namespace
{

bool block_cmp_data(const escript::Data&, boost::python::object o)
{
    PyErr_SetString(PyExc_TypeError,"Python relational operators are not defined for Data objects.");
    boost::python::throw_error_already_set();
    return false;
}


bool block_eq_data(const escript::Data&, boost::python::object o)
{
    PyErr_SetString(PyExc_TypeError,"The Python == and != operators are not defined for Data objects. "
      "To check for object identity use 'is'.  To check for numerical similarity of x and y, use Lsup(x-y)<TOL "
      "for a suitable tolerance.");
    boost::python::throw_error_already_set();
    return false;
}

bool block_cmp_functionspace(const escript::FunctionSpace&, boost::python::object o)
{
    PyErr_SetString(PyExc_TypeError,"Python relational operators are not defined for FunctionSpaces.");
    boost::python::throw_error_already_set();
    return false;
}

bool block_cmp_domains(const escript::AbstractDomain&, boost::python::object o)
{
    PyErr_SetString(PyExc_TypeError,"Python relational operators are not defined for Domains.");
    boost::python::throw_error_already_set();
    return false;
}
}

BOOST_PYTHON_MODULE(escriptcpp)
{
#if BOOST_VERSION >= 103500
// params are: bool show_user_defined, bool show_py_signatures, bool show_cpp_signatures
    docstring_options docopt(true,true,false);
#endif
#ifdef ESYS_HAVE_BOOST_NUMPY
    numpy::initialize();
#endif

    scope().attr("__doc__") = "To use this module, please import esys.escript";

    // register escript's default translators
    REGISTER_ESCRIPT_EXCEPTION_TRANSLATORS;
  def("setNumberOfThreads",escript::setNumberOfThreads,"Use of this method is strongly discouraged.");
  def("getNumberOfThreads",escript::getNumberOfThreads,"Return the maximum number of threads"
        " available to OpenMP.");
  def("releaseUnusedMemory",escript::DataTypes::releaseUnusedMemory);
#ifdef GIT_BUILD
  def("getVersion",escript::getGitVersion,"This method will only report accurate version numbers for clean checkouts.");
#else
  def("getVersion",escript::getSvnVersion,"This method will only report accurate version numbers for clean checkouts.");
#endif
  def("printParallelThreadCounts",escript::printParallelThreadCnt);
  def("getMPISizeWorld",escript::getMPISizeWorld,"Return number of MPI processes in the job.");
  def("getMPIRankWorld",escript::getMPIRankWorld,"Return the rank of this process in the MPI World.");
  def("MPIBarrierWorld",escript::MPIBarrierWorld,"Wait until all MPI processes have reached this point.");
  def("getMPIWorldMax",escript::getMPIWorldMax,"\n"
        "Each MPI process calls this function with a"
        " value for arg1. The maximum value is computed and returned.\n\n"
        ":rtype: int");
  def("getMPIWorldSum",escript::getMPIWorldSum,"\n"
        "Each MPI process calls this function with a"
        " value for arg1. The values are added up and the total value is returned.\n\n"
        ":rtype: int");
  def("runMPIProgram",escript::runMPIProgram,"Spawns an external MPI program using a separate communicator.");
  def("getMachinePrecision",escript::getMachinePrecision);
  def("getMaxFloat",escript::getMaxFloat);
  def("_saveDataCSV",escript::saveDataCSV, (args("filename","arg","sep","csep"), arg("refid")=false, arg("append")=false),
        "Saves data objects passed in a python dictionary to a file.\n"
        "The data objects must be over the same domain and be able to be interpolated to the same FunctionSpace.\n"
        "If one of the dictionary keys is named ``mask``, then only samples where ``mask`` has a positive\n"
        "value will be written to the file.\n\n"
        "A header line giving the names of each column will be output first.\n"
        "The keys given in the dictionary will be used to name columns.\n"
        "Then the data will be output, one line per sample (for all data).\n"
        "That is, items in each column will be printed in the same order.\n"
        "So you can be sure that values in the same row correspond to the same input value.\n\n"
        "\n"
        ":param filename:\n"
        ":type filename: ``string``\n"
        ":param arg: dictionary of named `Data` objects. If one is called ``mask`` it must be scalar data.\n"
        ":type arg: ``dict``\n"
        ":param sep: separator for columns (defaults to ',')\n"
        ":type sep: ``string``\n"
        ":param csep: separator for fields within data object (defaults to \"_\")\n"
        ":type csep: ``string``\n"
        ":param refid: If True, includes a column containing the element id numbers \n"
        ":type refid: ``string``\n"
        ":param append: If True, write to the end of ``filename``\n"
        ":type append: ``string``\n"
        "");
  def("_getNumpy",escript::getNumpy, arg("arg"),
        "Takes in a data object (or objects) and returns a numpy array\n"
        ":param arg: Data object\n"
        ":rtype: numpy ndarray\n"
        "");
  def("_convertToNumpy",escript::convertToNumpy, arg("arg"),
        "Takes in a data object (or objects) and returns a numpy array\n"
        ":param arg: Data object\n"
        ":rtype: numpy ndarray\n"
        "");
// #ifdef ESYS_HAVE_BOOST_NUMPY
//   def("_numpyToData", escript::numpyToData,(arg("array"), arg("isComplex"), arg("functionspace")),
//         "Takes in a numpy ndarray and function space and returns a Data object\n"
//         ":param array: A numpy ndarray\n"
//         ":param isComplex: boolean. True for complex data \n"
//         ":param functionspace: A FunctionSpace\n"
//         ":rtype: Data object\n"
//         "");
// #endif
  def("canInterpolate", &escript::canInterpolate, args("src", "dest"),":param src: Source FunctionSpace\n"
        ":param dest: Destination FunctionSpace\n"
        ":return: True if src can be interpolated to dest\n"
        ":rtype: `bool`");

  //
  // Interface for AbstractDomain
  //
  class_<escript::AbstractDomain, escript::Domain_ptr, boost::noncopyable >("Domain","Base class for all domains.",no_init)
     .def("getStatus",&escript::AbstractDomain::getStatus,"The status of a domain changes whenever the domain is modified\n\n"
        ":rtype: int")
     .def("setTagMap",&escript::AbstractDomain::setTagMap,args("name","tag"),
        "Give a tag number a name.\n\n"
        ":param name: Name for the tag\n"
        ":type name: ``string``\n"
        ":param tag: numeric id\n"
        ":type tag: ``int``\n"
        ":note: Tag names must be unique within a domain")
     .def("getTag",&escript::AbstractDomain::getTag,args("name"),":return: tag id for "
        "``name``\n"
        ":rtype: ``string``")
     .def("isValidTagName",&escript::AbstractDomain::isValidTagName,args("name"),
        ":return: True is ``name`` corresponds to a tag\n"
        ":rtype: ``bool``")
     .def("showTagNames",&escript::AbstractDomain::showTagNames,":return: A space separated list of tag names\n"
        ":rtype: ``string``")
     .def("getX",&escript::AbstractDomain::getX,":rtype: `Data`\n"
        ":return: Locations in the"
        "`Domain`. FunctionSpace is chosen appropriately")
#ifdef ESYS_HAVE_BOOST_NUMPY
     .def("getNumpyX",&escript::AbstractDomain::getNumpyX,":rtype: `numpy ndarray`\n"
        ":return: Locations in the"
        "`Domain`. FunctionSpace is chosen appropriately")
     .def("isCellOriented", &escript::AbstractDomain::isCellOriented, arg("functionSpaceCode"),
        ":return: true is the data is cell centered.\n"
        ":rtype: ``int``")
#endif
#if defined(ESYS_MPI) && defined(ESYS_HAVE_MPI4PY)
      // .def("setMPIComm", &escript::AbstractDomain::setMPIComm, 
      //       (arg("COMM")),
      //       "Sets the escript MPI comm to the COMM"
      //       ":param COMM: An MPI communicator"
      //       )
#endif
     .def("getDim",&escript::AbstractDomain::getDim,":rtype: `int`\n"
        ":return: Spatial dimension of the `Domain`")
     .def("getNormal",&escript::AbstractDomain::getNormal,":rtype: `escript`\n"
        ":return: Boundary normals")
     .def("getSize",&escript::AbstractDomain::getSize,":return: the local size of samples. The function space is chosen appropriately\n"
        ":rtype: `Data`")
     .def("dump",&escript::AbstractDomain::dump,args("filename"),"Dumps the domain to a file\n\n"
        ":param filename: \n"
        ":type filename: string")

     .def("getMPISize",&escript::AbstractDomain::getMPISize,":return: the number of processes used for this `Domain`\n"
        ":rtype: ``int``")
     .def("getMPIRank",&escript::AbstractDomain::getMPIRank,":return: the rank of this process\n"
        ":rtype: ``int``")
     .def("getMPIComm",&escript::AbstractDomain::getMPICommPython,":return: the MPI communicator for this domain as an mpi4py.MPI.Comm object (or None if MPI/mpi4py not enabled)\n"
        ":rtype: mpi4py.MPI.Comm or None")
     .def("MPIBarrier",&escript::AbstractDomain::MPIBarrier,"Wait until all processes have reached this point")
     .def("isRootRank",&escript::AbstractDomain::isRootRank,":return: True if this code is executing on the root rank (rank 0)\n"
        ":rtype: `bool`")
     .def("supportsContactElements", &escript::AbstractDomain::supportsContactElements,"Does this domain support contact elements.")
     .def(self == self)
     .def(self != self)
     .def("__lt__", block_cmp_domains)
     .def("__le__", block_cmp_domains)
     .def("__gt__", block_cmp_domains)
     .def("__ge__", block_cmp_domains);

  //
  // Interface for AbstractContinuousDomain
  //
  class_<escript::AbstractContinuousDomain, bases<escript::AbstractDomain>, boost::noncopyable >("ContinuousDomain","Class representing continuous domains",no_init)
       .def("getSystemMatrixTypeId", &escript::AbstractContinuousDomain::getSystemMatrixTypeId,
args("options"),
        ":return: the identifier of the matrix type to be used for the global stiffness matrix "
        "when particular solver options are used.\n"
        ":rtype: int")
       .def("getTransportTypeId",&escript::AbstractContinuousDomain::getTransportTypeId,
args("solver", "preconditioner", "package", "symmetry"))

      .def("addPDEToSystem",&escript::AbstractContinuousDomain::addPDEToSystem,
args("mat", "rhs","A", "B", "C", "D", "X", "Y", "d", "y", "d_contact", "y_contact", "d_dirac", "y_dirac"),
        "adds a PDE onto the stiffness matrix mat and a rhs\n\n"
        ":param mat:\n"
        ":type mat: `OperatorAdapter`\n"
        ":param rhs:\n"
        ":type rhs: `Data`\n"
        ":param A:\n"
        ":type A: `Data`\n"
        ":param B:\n"
        ":type B: `Data`\n"
        ":param C:\n"
        ":type C: `Data`\n"
        ":param D:\n"
        ":type D: `Data`\n"
        ":param X:\n"
        ":type X: `Data`\n"
        ":param Y:\n"
        ":type Y: `Data`\n"
        ":param d:\n"
        ":type d: `Data`\n"
        ":param d_contact:\n"
        ":type d_contact: `Data`\n"
        ":param y_contact:\n"
        ":type y_contact: `Data`\n"
        ":param d_dirac:\n"
        ":type d_dirac: `Data`\n"
        ":param y_dirac:\n"
        ":type y_dirac: `Data`\n"
)
      .def("addPDEToRHS",&escript::AbstractContinuousDomain::addPDEToRHS,
args("rhs", "X", "Y", "y", "y_contact", "y_dirac"),
        "adds a PDE onto the stiffness matrix mat and a rhs\n\n"
        ":param rhs:\n"
        ":type rhs: `Data`\n"
        ":param X:\n"
        ":type X: `Data`\n"
        ":param Y:\n"
        ":type Y: `Data`\n"
        ":param y:\n"
        ":type y: `Data`\n"
        ":param y_contact:\n"
        ":type y_contact: `Data`\n"
        ":param y_dirac:\n"
        ":type y_dirac: `Data`"
)
      .def("addPDEToTransportProblem",&escript::AbstractContinuousDomain::addPDEToTransportProblem,
args( "tp", "source", "M", "A", "B", "C", "D", "X", "Y", "d", "y", "d_contact", "y_contact", "d_dirac", "y_dirac"),
        ":param tp:\n"
        ":type tp: `TransportProblemAdapter`\n"
        ":param source:\n"
        ":type source: `Data`\n"
        ":param M:\n"
        ":type M: `Data`\n"
        ":param A:\n"
        ":type A: `Data`\n"
        ":param B:\n"
        ":type B: `Data`\n"
        ":param C:\n"
        ":type C: `Data`\n"
        ":param D:\n"
        ":type D: `Data`\n"
        ":param X:\n"
        ":type X: `Data`\n"
        ":param Y:\n"
        ":type Y: `Data`\n"
        ":param d:\n"
        ":type d: `Data`\n"
        ":param y:\n"
        ":type y: `Data`\n"
        ":param d_contact:\n"
        ":type d_contact: `Data`\n"
        ":param y_contact:\n"
        ":type y_contact: `Data`\n"
        ":param d_dirac:\n"
        ":type d_dirac: `Data`\n"
        ":param y_dirac:\n"
        ":type y_dirac: `Data`\n"
)
      .def("newOperator",&escript::AbstractContinuousDomain::newSystemMatrix,
args("row_blocksize", "row_functionspace", "column_blocksize", "column_functionspace", "type"),
        "creates a SystemMatrixAdapter stiffness matrix and initializes it with zeros\n\n"
        ":param row_blocksize:\n"
        ":type row_blocksize: ``int``\n"
        ":param row_functionspace:\n"
        ":type row_functionspace: `FunctionSpace`\n"
        ":param column_blocksize:\n"
        ":type column_blocksize: ``int``\n"
        ":param column_functionspace:\n"
        ":type column_functionspace: `FunctionSpace`\n"
        ":param type:\n"
        ":type type: ``int``\n"
)
      .def("newTransportProblem",&escript::AbstractContinuousDomain::newTransportProblem,
args("theta", "blocksize", "functionspace", "type"),
        "creates a TransportProblemAdapter\n\n"
        ":param theta:\n"
        ":type theta: ``float``\n"
        ":param blocksize:\n"
        ":type blocksize: ``int``\n"
        ":param functionspace:\n"
        ":type functionspace: `FunctionSpace`\n"
        ":param type:\n"
        ":type type: ``int``\n"
)
      .def("getDataShape",&escript::AbstractContinuousDomain::getDataShape, args("functionSpaceCode"),
        ":return: a pair (dps, ns) where dps=the number of data points per sample, and ns=the number of samples\n"
        ":rtype: ``tuple``")
      .def("print_mesh_info",&escript::AbstractContinuousDomain::Print_Mesh_Info,(arg("full")=false),
        ":param full:\n"
        ":type full: ``bool``")
      .def("getDescription",&escript::AbstractContinuousDomain::getDescription,
        ":return: a description for this domain\n"
        ":rtype: ``string``")
      .def("setX",&escript::AbstractContinuousDomain::setNewX,
args("arg"), "assigns new location to the domain\n\n"
        ":param arg:\n"
        ":type arg: `Data`")
      .def("getNumDataPointsGlobal",&escript::AbstractContinuousDomain::getNumDataPointsGlobal,
        ":return: the number of data points summed across all MPI processes\n"
        ":rtype: ``int``")
#ifdef ESYS_HAVE_TRILINOS
      .def("finaliseA",&escript::AbstractContinuousDomain::finaliseA, args("mat", "iz"),
        "finalises the matrix system before it is passed to the solver\n\n"
        ":param mat:\n"
        ":type mat: ``AbstractSystemMatrix``\n"
        ":param iz:\n"
        ":type iz: ``CRSMatrix``\n"
        ":rtype ``AbstractSystemMatrix``")
      .def("finaliseRhs",&escript::AbstractContinuousDomain::finaliseRhs, args("rhs", "z"),
        "finalises the matrix system before it is passed to the solver\n\n"
        ":param rhs:\n"
        ":type rhs: ``AbstractSystemMatrix``\n"
        ":param z:\n"
        ":type z: ``CRSMatrix``\n"
        ":rtype ``AbstractSystemMatrix``");
#else
        ;
#endif //ESYS_HAVE_TRILINOS



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
  fs_definer.def("getDim",&escript::FunctionSpace::getDim,":return: the spatial dimension of the underlying domain.\n"
        ":rtype: int");
//   fs_definer.def("getDomain",&escript::FunctionSpace::getDomain,
//                  return_internal_reference<>());
  fs_definer.def("getDomain",&escript::FunctionSpace::getDomainPython,":return: the underlying `Domain` for this FunctionSpace.\n"
        ":rtype: `Domain`");
  fs_definer.def("getMPIComm",&escript::FunctionSpace::getMPIComm,":return: the MPI communicator for this function space's domain as an mpi4py.MPI.Comm object (or None if MPI/mpi4py not enabled)\n"
        ":rtype: mpi4py.MPI.Comm or None");
  fs_definer.def("getX",&escript::FunctionSpace::getX,"\n"
        ":return: a function whose values are its input coordinates. ie an identity function.\n"
        ":rtype: `Data`");
  fs_definer.def("getNormal",&escript::FunctionSpace::getNormal,":return: the surface normal field.\n\n"
        ":rtype: `Data`");
  fs_definer.def("getSize",&escript::FunctionSpace::getSize,":return: sample size\n"
        ":rtype: `Data`");
  fs_definer.def("setTags",&escript::FunctionSpace::setTags,args("newtag","mask"),
        "Set tags according to a mask\n\n"
        ":param newtag: tag number to set\n"
        ":type newtag: string, non-zero ``int``\n"
        ":param mask: Samples which correspond to positive values in the mask will be set to ``newtag``.\n"
        ":type mask: scalar `Data`");
  fs_definer.def("setTags",&escript::FunctionSpace::setTagsByString,args("newtag","mask"));
  fs_definer.def("getTagFromDataPointNo",
                 &escript::FunctionSpace::getTagFromDataPointNo,":return: the tag associated with the given sample number.\n"
        ":rtype: int");
  fs_definer.def("getReferenceIDFromDataPointNo", &escript::FunctionSpace::getReferenceIDFromDataPointNo,args("dataPointNo"),":return: the reference number associated with ``dataPointNo``\n"
        ":rtype: int ");
  fs_definer.def("getListOfTags",&escript::FunctionSpace::getListOfTags,":return: a list of the tags used in this function space\n"
        ":rtype: ``list``");
  fs_definer.def("getApproximationOrder", &escript::FunctionSpace::getApproximationOrder,":return: the approximation order referring to the maximum degree of a polynomial which can be represented exactly in interpolation and/or integration.\n"
        ":rtype: ``int``");
  fs_definer.def("getTypeCode",&escript::FunctionSpace::getTypeCode,":rtype: `int`");
  fs_definer.def("__str__", &escript::FunctionSpace::toString);
  fs_definer.def("__lt__",block_cmp_functionspace);
  fs_definer.def("__le__",block_cmp_functionspace);
  fs_definer.def("__gt__",block_cmp_functionspace);
  fs_definer.def("__ge__",block_cmp_functionspace);
  fs_definer.def(self == self);
  fs_definer.def(self != self);



  //
  // Interface for Data
  //
  class_<escript::Data>("Data"/*,shared_ptr<Data>*/, "Represents a collection of datapoints. It is used to store the values of a function. For more details please consult the c++ class documentation.",init<>())
    // various constructors for Data objects
    .def(init<object, optional<object, object, object>>(args("value", "p2", "p3", "p4")))
    // Note for Lutz, Need to specify the call policy in order to return a
    // reference. In this case return_internal_reference.
    .def("__str__",&escript::Data::toString)
    .def("getDomain",&escript::Data::getDomainPython,":rtype: `Domain`")
    .def("getFunctionSpace",&escript::Data::getFunctionSpace,return_value_policy<copy_const_reference>(),":rtype: `FunctionSpace`")
    .def("getMPIComm",&escript::Data::getMPIComm,":return: the MPI communicator for this data's domain as an mpi4py.MPI.Comm object (or None if MPI/mpi4py not enabled)\n"
        ":rtype: mpi4py.MPI.Comm or None")
    .def("getX",&escript::Data::getXFromFunctionSpace,
        "Returns the spatial coordinates of the spatial nodes.\n"
        ":rtype: `Data`")
    .def("isEmpty",&escript::Data::isEmpty,"Is this object an instance of ``DataEmpty``\n\n"
        ":rtype: ``bool``\n"
        ":note: This is not the same thing as asking if the object contains datapoints.")
    .def("isProtected",&escript::Data::isProtected,"Can this instance be modified.\n"
        ":rtype: ``bool``")
    .def("setProtection",&escript::Data::setProtection,"Disallow modifications to this data object\n\n"
        ":note: This method does not allow you to undo protection.")
    .def("getShape",&escript::Data::getShapeTuple,"\n"
        "Returns the shape of the datapoints in this object as a python tuple. Scalar data has the shape ``()``\n\n"
        ":rtype: ``tuple``")
    .def("getRank",&escript::Data::getDataPointRank,":return: the number of indices required to address a component of a datapoint\n"
        ":rtype: positive ``int``")
    .def("dump",&escript::Data::dump,args("fileName"),"Save the data as a HDF5 file\n\n"
        ":param fileName: \n"
        ":type fileName: ``string``")
    .def("toListOfTuples",&escript::Data::toListOfTuples, (arg("scalarastuple")=false),
        "Return the datapoints of this object in a list. Each datapoint is stored as a tuple.\n\n"
        ":param scalarastuple: if True, scalar data will be wrapped as a tuple."
        " True => [(0), (1), (2)]; False => [0, 1, 2]")
    .def("copyWithMask",&escript::Data::copyWithMask,args("other","mask"),
        "Selectively copy values from ``other`` `Data`."
        "Datapoints which correspond to positive values in ``mask`` will be copied from ``other``\n"
        "\n"
        ":param other: source of values\n"
        ":type other: `Data`\n"
        ":param mask:\n"
        ":type mask: Scalar `Data`")
    .def("setTaggedValue",&escript::Data::setTaggedValue,args("tagKey","value"),
        "Set the value of tagged Data.\n\n"
        ":param tagKey: tag to update\n"
        ":type tagKey: ``int``\n"
        "")
    .def("setTaggedValue",&escript::Data::setTaggedValueByName,args("name","value"),":param name: tag to update\n"
        ":type name: ``string``\n"
        ":param value: value to set tagged data to\n"
        ":type value: ``object`` which acts like an array, ``tuple`` or ``list``\n"
        "")
    .def("getNumberOfDataPoints",&escript::Data::getNumDataPoints,":rtype: ``int``\n"
        ":return: Number of datapoints in the object")
    .def("isExpanded",&escript::Data::isExpanded,":rtype: ``bool``\n"
        ":return: True if this ``Data`` is expanded.")
    .def("isTagged",&escript::Data::isTagged,":rtype: ``bool``\n"
        ":return: True if this ``Data`` is expanded.")
    .def("isConstant",&escript::Data::isConstant,":rtype: ``bool``\n"
        ":return: True if this ``Data`` is an instance of ``DataConstant``\n"
        ":note: This does not mean the data is immutable.")
    .def("isLazy",&escript::Data::isLazy,":rtype: ``bool``\n"
        ":return: True if this ``Data`` is lazy.")
    .def("isReady",&escript::Data::isReady,":rtype: ``bool``\n"
        ":return: True if this ``Data`` is not lazy.")
    .def("isComplex", &escript::Data::isComplex,":rtype: ``bool``\n"
	":return: True if this ``Data`` stores complex values.")
    .def("expand",&escript::Data::expand,"Convert the data to expanded representation if it is not expanded already.")
    .def("hasNaN",&escript::Data::hasNaN,"Returns return true if data contains NaN. [Note that for complex values, hasNaN and hasInf are not mutually exclusive.]")
    .def("replaceNaN",&escript::Data::replaceNaNPython,args("value"),"Replaces NaN values with value. [Note, for complex Data, both real and imaginary components are replaced even if only one part is NaN].")
    .def("hasInf",&escript::Data::hasInf,"Returns return true if data contains +-Inf.  [Note that for complex values, hasNaN and hasInf are not mutually exclusive.]")
    .def("replaceInf",&escript::Data::replaceInfPython,args("value"),"Replaces +-Inf values with value. [Note, for complex Data, both real and imaginary components are replaced even if only one part is Inf].")
    .def("tag",&escript::Data::tag,"Convert data to tagged representation if it is not already tagged or expanded")
    .def("resolve",&escript::Data::resolve,"Convert the data to non-lazy representation.")
    .def("copy",&escript::Data::copy,args("other"),"Make this object a copy of ``other``\n"
        "\n"
        ":note: The two objects will act independently from now on. That is, changing ``other`` "
        "after this call will not change this object and vice versa.")
    .def("copy",&escript::Data::copySelf,":note: In the no argument form, a new object will be returned which is an independent copy of this object.")
    .def("delay",&escript::Data::delay,"Convert this object into lazy representation")
    .def("setValueOfDataPoint",&escript::Data::setValueOfDataPointToPyObject,args("dataPointNo","value"))
    .def("setValueOfDataPoint",&escript::Data::setValueOfDataPointToArray)
    .def("_setTupleForGlobalDataPoint", &escript::Data::setTupleForGlobalDataPoint)
    .def("setValueOfDataPoint",&escript::Data::setValueOfDataPoint,"\n"
        "Modify the value of a single datapoint.\n\n"
        ":param dataPointNo:\n"
        ":type dataPointNo: int\n"
        ":param value: \n"
        ":type value: ``float`` or an object which acts like an array, ``tuple`` or ``list``\n"
        ":warning: Use of this operation is discouraged. It prevents some optimisations from operating.")
    .def("getTupleForDataPoint",&escript::Data::getValueOfDataPointAsTuple,args("dataPointNo"),
        ":return: Value of the specified datapoint\n"
        ":rtype: ``tuple``\n"
        ":param dataPointNo: datapoint to access\n"
        ":type dataPointNo: ``int``")
    .def("getTupleForGlobalDataPoint",&escript::Data::getValueOfGlobalDataPointAsTuple,args("procNo","dataPointNo"),"Get a specific datapoint from a specific process\n\n"
        ":rtype: ``tuple``\n"
        ":param procNo: MPI rank of the process\n"
        ":type procNo: positive ``int``"
        "\n"
        ":param dataPointNo: datapoint to access\n"
        ":type dataPointNo: int")
    .def("setToZero",&escript::Data::setToZero,"After this call the object will store values of the same shape as before but all components will be zero.")
    .def("interpolate",&escript::Data::interpolate,args("functionspace"),"Interpolate this object's values into a new functionspace.")
    .def("_interpolateTable3d", &escript::Data::interpolateFromTable3DP,
(arg("table"),arg("Amin"),arg("Astep"), arg("B"), arg("Bmin"), arg("Bstep"), arg("C"), arg("Cmin"), arg("Cstep"), arg("undef")=1.e50, arg("check_boundaries")=false, "For internal use only. Please use the interpolateTable function.")
)

    .def("interpolateTable", &escript::Data::interpolateFromTable2DP,
(arg("table"),arg("Amin"),arg("Astep"), arg("B"), arg("Bmin"), arg("Bstep"), arg("undef")=1.e50, arg("check_boundaries")=false),
        "Creates a new Data object by interpolating using the source data (which are\n"
        "looked up in ``table``)\n"
        "``A`` must be the outer dimension on the table\n\n"
        ":param table: two dimensional collection of values\n"
        ":param Amin: The base of locations in table\n"
        ":type Amin: float\n"
        ":param Astep: size of gap between each item in the table\n"
        ":type Astep: float\n"
        ":param undef: upper bound on interpolated values\n"
        ":type undef: float\n"
        ":param B: Scalar representing the second coordinate to be mapped into the table\n"
        ":type B: `Data`\n"
        ":param Bmin: The base of locations in table for 2nd dimension\n"
        ":type Bmin: float\n"
        ":param Bstep: size of gap between each item in the table for 2nd dimension\n"
        ":type Bstep: float\n"
        ":param check_boundaries: if true, then values outside the boundaries will be rejected. If false, then boundary values will be used.\n"
        ":raise RuntimeError(DataException): if the coordinates do not map into the table or if the interpolated value is above ``undef``"
        "\n"
        ":rtype: `Data`"
)
    .def("interpolateTable", &escript::Data::interpolateFromTable1DP,
(arg("table"),arg("Amin"),arg("Astep"), arg("undef")=1.e50, arg("check_boundaries")=false)/*,
        "Creates a new Data object by interpolating using the source data (which are\n"
        "looked up in ``table``)\n\n"
        ":param table: one dimensional collection of values\n"
        ":param Amin: The base of locations in table\n"
        ":type Amin: float\n"
        ":param Astep: size of gap between each item in the table\n"
        ":type Astep: float\n"
        ":param undef: upper bound on interpolated values\n"
        ":type undef: float\n"
        ":param check_boundaries: if true, then values outside the boundaries will be rejected. If false, then boundary values will be used.\n"
        ":raise RuntimeError(DataException): if the coordinates do not map into the table or if the interpolated value is above ``undef``"
        "\n"
        ":rtype: `Data`"
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
    .def("internal_minGlobalDataPoint",&escript::Data::minGlobalDataPoint,"Please consider using getInfLocator() from pdetools instead.")
    .def("internal_maxGlobalDataPoint",&escript::Data::maxGlobalDataPoint, "Please consider using getSupLocator() from pdetools instead.")
    .def("getTagNumber",&escript::Data::getTagNumber,args("dpno"),"Return tag number for the specified datapoint\n\n"
        ":rtype: int\n"
        ":param dpno: datapoint number\n"
        ":type dpno: int")
    // Unary functions for Data
    .def("conjugate", &escript::Data::conjugate)
    .def("real", &escript::Data::real)
    .def("imag", &escript::Data::imag)
    .def("promote", &escript::Data::complicate)
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
    .def("_besselFirstKind",&escript::Data::besselFirstKind,arg("order"))
    .def("_besselSecondKind",&escript::Data::besselSecondKind,arg("order"))
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
    .def("_antisymmetric",&escript::Data::antisymmetric)
    .def("_hermitian",&escript::Data::hermitian)
    .def("_antihermitian",&escript::Data::antihermitian)
    .def("_trace",&escript::Data::trace)
    .def("_swap_axes",&escript::Data::swapaxes)
    .def("_eigenvalues",&escript::Data::eigenvalues)
    .def("_eigenvalues_and_eigenvectors",&escript::Data::eigenvalues_and_eigenvectors,(arg("tol")=1.e-13))
    // functions returning a single real number:
    .def("_Lsup",&escript::Data::Lsup,":return: the Lsup-norm of the object\n"
        ":rtype: float\n"
        ":note: If the ``Data`` contains no values, zero will be returned instead.")
    .def("_sup",&escript::Data::sup,":return: the maximum value over all data points.\n"
        ":rtype: float\n"
        ":note: If the ``Data`` contains no values a large negative value will be returned instead.")
    .def("_inf",&escript::Data::inf,":return: minimum value over all components and all data points\n"
        ":rtype: float\n"
        ":note: If the ``Data`` contains no values a large positive value will be returned instead.")
    .def("_integrateToTuple",&escript::Data::integrateToTuple,":return: Calculate the integral over the function space domain as a python tuple\n"
        ":rtype: tuple")
    // following implements the python abs operator
    .def("__abs__",&escript::Data::abs,":return: absolute value\n\n"
        ":rtype: `Data`")
    // following implements the python "-" negation operator
    .def("__neg__",&escript::Data::neg, ":return: negation of the values in this object\n"
        ":rtype: `Data`")
    // following implements the python "+" identity operator
    .def("__pos__",&escript::Data::pos, "\n"
        "The unary + operator\n\n"
        ":rtype: `Data`")
    // following three functions implement the python [] operator
    .def("__getitem__",&escript::Data::getItem,"Used by the python [] operator\n\n"
        ":rtype: `Data`")
    .def("__setitem__",&escript::Data::setItemO,"Used by the python [] operator")
    .def("__setitem__",&escript::Data::setItemD,"Used by the python [] operator")
    // following three functions implement the python ** operator
    .def("__pow__",&escript::Data::powO,"Used by the python ** operator\n\n"
        ":rtype: `Data`")
    .def("__pow__",&escript::Data::powD)
    .def("__rpow__",&escript::Data::rpowO,"\n"
        "Used by the python ** operator\n\n"
        ":rtype: `Data`")
    // following two functions implement the newer python / operator
    .def("__truediv__",&escript::Data::truedivO)
    .def("__truediv__",&escript::Data::truedivD)
    .def("__rtruediv__",&escript::Data::rtruedivO)
    .def("__gt__",block_cmp_data)
    .def("__lt__",block_cmp_data)
    .def("__le__",block_cmp_data)
    .def("__ge__",block_cmp_data)
    .def("phase",&escript::Data::phase)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wself-assign-overloaded"
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
#pragma clang diagnostic pop
    .def("_inverse", &escript::Data::matrixInverse, ":return: inverse of square matrices\n"
        "")
//    .def("__add__", &escript::Data::addOperatorD)
    .def("__add__", &escript::Data::__add__)
    .def("__radd__", &escript::Data::__add__)  // its the same coz + is commutative
    .def("__sub__", &escript::Data::__sub__)
    .def("__rsub__", &escript::Data::__rsub__)
    .def("__mul__", &escript::Data::__mul__)
    .def("__rmul__", &escript::Data::__mul__)   // commutative
    .def("__div__", &escript::Data::__div__)
    .def("__rdiv__", &escript::Data::__rdiv__)   // commutative
    .def("__eq__", block_eq_data)		// stop people from using ==
    .def("__ne__", block_eq_data)		// stop people from using !=
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
  def("ReducedFunction",escript::reducedFunction, args("domain"),":return: a function FunctionSpace with reduced integration order\n"
        ":rtype: `FunctionSpace`");
  def("FunctionOnBoundary",escript::functionOnBoundary, args("domain"), ":return: a function on boundary FunctionSpace\n"
        ":rtype: `FunctionSpace`");
  def("ReducedFunctionOnBoundary",escript::reducedFunctionOnBoundary, args("domain"),
        ":return: a function on boundary FunctionSpace with reduced integration order\n"
        ":rtype: `FunctionSpace`");
  def("FunctionOnContactZero",escript::functionOnContactZero, args("domain"), ":return: Return a FunctionSpace on left side of contact\n"
        ":rtype: `FunctionSpace`");
  def("ReducedFunctionOnContactZero",escript::reducedFunctionOnContactZero, args("domain"),
        ":return: a FunctionSpace  on left side of contact with reduced integration order\n"
        ":rtype: `FunctionSpace`");
  def("FunctionOnContactOne",escript::functionOnContactOne, args("domain"), ":return: Return a FunctionSpace on right side of contact\n"
        ":rtype: `FunctionSpace`");
  def("ReducedFunctionOnContactOne",escript::reducedFunctionOnContactOne, args("domain"),
        ":return: Return a FunctionSpace on right side of contact with reduced integration order\n"
        ":rtype: `FunctionSpace`");
  def("Solution",escript::solution, args("domain"), ":rtype: `FunctionSpace`");
  def("ReducedSolution",escript::reducedSolution, args("domain"), ":rtype: `FunctionSpace`");
  def("DiracDeltaFunctions",escript::diracDeltaFunctions, args("domain"), ":rtype: `FunctionSpace`");





  //
  // Factory methods for Data
  //
 def("load",escript::load_hdf5, args("fileName","domain"), "reads real Data on domain from file in HDF5 format\n\n"
        ":param fileName:\n"
        ":type fileName: ``string``\n"
        ":param domain:\n"
        ":type domain: `Domain`");
  def("loadIsConfigured",escript::loadConfigured,":return: True if the load function is configured.");
  def("Scalar",escript::ScalarFromObj,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false),
        "Construct a Data object containing scalar data-points.\n\n"
        ":param value: scalar value for all points\n"
        "\n"
        ":rtype: `Data`\n"
        ":type value: float\n"
        ":param what: FunctionSpace for Data\n"
        ":type what: `FunctionSpace`\n"
        ":param expanded: If True, a value is stored for each point. If False, more efficient representations may be used\n"
        ":type expanded: ``bool``");
  def("ComplexScalar",escript::ComplexScalarFromObj,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false),
        "Construct a Data object containing scalar data-points.\n\n"
        ":param value: scalar value for all points\n"
        "\n"
        ":rtype: `Data`\n"
        ":type value: float\n"
        ":param what: FunctionSpace for Data\n"
        ":type what: `FunctionSpace`\n"
        ":param expanded: If True, a value is stored for each point. If False, more efficient representations may be used\n"
        ":type expanded: ``bool``");
  def("Vector",escript::Vector,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false),
        "Construct a Data object containing rank1 data-points.\n\n"
        ":param value: scalar value for all points\n"
        "\n"
        ":rtype: `Data`\n"
        ":type value: float\n"
        ":param what: FunctionSpace for Data\n"
        ":type what: `FunctionSpace`\n"
        ":param expanded: If True, a value is stored for each point. If False, more efficient representations may be used\n"
        ":type expanded: ``bool``");
  def("ComplexVector",escript::ComplexVector,
    (arg("value")=0.0,
     arg("what")=escript::FunctionSpace(),
     arg("expanded")=false),
      "Construct a Data object containing rank1 data-points.\n\n"
      ":param value: scalar value for all points\n"
      "\n"
      ":rtype: `Data`\n"
      ":type value: float\n"
      ":param what: FunctionSpace for Data\n"
      ":type what: `FunctionSpace`\n"
      ":param expanded: If True, a value is stored for each point. If False, more efficient representations may be used\n"
      ":type expanded: ``bool``");
 def("Vector", escript::VectorFromObj,
      (arg("value"),
	arg("what")=escript::FunctionSpace(),
	arg("expanded")=false));
  def("ComplexVector", escript::ComplexVectorFromObj,
         (arg("value"),
   	arg("what")=escript::FunctionSpace(),
   	arg("expanded")=false));
  def("Tensor",escript::Tensor,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false),
        "Construct a Data object containing rank2 data-points.\n\n"
        ":param value: scalar value for all points\n"
        "\n"
        ":rtype: `Data`\n"
        ":type value: float\n"
        ":param what: FunctionSpace for Data\n"
        ":type what: `FunctionSpace`\n"
        ":param expanded: If True, a value is stored for each point. If False, more efficient representations may be used\n"
        ":type expanded: ``bool``");
    def("ComplexTensor",escript::ComplexTensor,
        (arg("value")=0.0,
         arg("what")=escript::FunctionSpace(),
         arg("expanded")=false),
          "Construct a Data object containing rank2 data-points.\n\n"
          ":param value: scalar value for all points\n"
          "\n"
          ":rtype: `Data`\n"
          ":type value: float\n"
          ":param what: FunctionSpace for Data\n"
          ":type what: `FunctionSpace`\n"
          ":param expanded: If True, a value is stored for each point. If False, more efficient representations may be used\n"
          ":type expanded: ``bool``");
 def("Tensor", escript::TensorFromObj,
      (arg("value"),
	arg("what")=escript::FunctionSpace(),
	arg("expanded")=false));
def("ComplexTensor", escript::ComplexTensorFromObj,
         (arg("value"),
   	arg("what")=escript::FunctionSpace(),
   	arg("expanded")=false));
  def("Tensor3",escript::Tensor3,
      (arg("value")=0.0,
       arg("what")=escript::FunctionSpace(),
       arg("expanded")=false),
        "Construct a Data object containing rank3 data-points.\n\n"
        ":param value: scalar value for all points\n"
        "\n"
        ":rtype: `Data`\n"
        ":type value: float\n"
        ":param what: FunctionSpace for Data\n"
        ":type what: `FunctionSpace`\n"
        ":param expanded: If True, a value is stored for each point. If False, more efficient representations may be used\n"
        ":type expanded: ``bool``"
);
def("ComplexTensor3",escript::ComplexTensor3,
    (arg("value")=0.0,
     arg("what")=escript::FunctionSpace(),
     arg("expanded")=false),
      "Construct a Data object containing rank3 data-points.\n\n"
      ":param value: scalar value for all points\n"
      "\n"
      ":rtype: `Data`\n"
      ":type value: float\n"
      ":param what: FunctionSpace for Data\n"
      ":type what: `FunctionSpace`\n"
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
        "Construct a Data object containing rank4 data-points.\n\n"
        ":param value: scalar value for all points\n"
        "\n"
        ":rtype: `Data`\n"
        ":type value: float\n"
        ":param what: FunctionSpace for Data\n"
        ":type what: `FunctionSpace`\n"
        ":param expanded: If True, a value is stored for each point. If False, more efficient representations may be used\n"
        ":type expanded: ``bool``"
);
def("ComplexTensor3", escript::ComplexTensor3FromObj,
     (arg("value"),
   arg("what")=escript::FunctionSpace(),
   arg("expanded")=false));
 def("ComplexTensor4",escript::ComplexTensor4,
     (arg("value")=0.0,
      arg("what")=escript::FunctionSpace(),
      arg("expanded")=false),
       "Construct a Data object containing rank4 data-points.\n\n"
       ":param value: scalar value for all points\n"
       "\n"
       ":rtype: `Data`\n"
       ":type value: float\n"
       ":param what: FunctionSpace for Data\n"
       ":type what: `FunctionSpace`\n"
       ":param expanded: If True, a value is stored for each point. If False, more efficient representations may be used\n"
       ":type expanded: ``bool``"
);
 def("Tensor4", escript::Tensor4FromObj,
      (arg("value"),
	arg("what")=escript::FunctionSpace(),
	arg("expanded")=false));
def("ComplexTensor4", escript::ComplexTensor4FromObj,
         (arg("value"),
   	arg("what")=escript::FunctionSpace(),
   	arg("expanded")=false));
def("ComplexData", escript::ComplexData,
         (arg("value"),
   	arg("what")=escript::FunctionSpace(),
   	arg("expanded")=false));

 def("RandomData", escript::randomData, (arg("shape"), arg("fs"), arg("seed")=0, arg("filter")=boost::python::tuple()),
        "Creates a new expanded Data object containing pseudo-random values. With no filter, values are drawn uniformly at random from [0,1].\n\n"
        ":param shape: datapoint shape\n"
        ":type shape: tuple\n"
        ":param fs: function space for data object.\n"
        ":type fs: `FunctionSpace`\n"
        ":param seed: seed for random number generator.\n"
        ":type seed: long\n"
        "");

  //
  // Binary operators
  //
  def("C_GeneralTensorProduct",escript::C_GeneralTensorProduct,
      (arg("arg0"),
       arg("arg1"),
       arg("axis_offset")=0,
       arg("transpose")=0),
        "Compute a tensor product of two Data objects.\n\n"
        ":rtype: `Data`\n"
        ":param arg0:\n"
        ":param arg1:\n"
        ":param axis_offset:\n"
        ":type axis_offset: ``int``\n"
        ":param transpose: 0: transpose neither, 1: transpose arg0, 2: transpose arg1\n"
        ":type transpose: int");

  //
  // Interface for AbstractSystemMatrix
  //
  class_<escript::AbstractSystemMatrix,escript::ASM_ptr, boost::noncopyable>("Operator","",init<>())    // Doco goes in the empty string param
     .def("isEmpty",&escript::AbstractSystemMatrix::isEmpty,":rtype: ``bool``\n"
        ":return: True if matrix is empty")
     .def("solve",&escript::AbstractSystemMatrix::solve, args("in","options"),
        ":return: the solution *u* of the linear system *this*u=in*\n\n"
        ":param in:\n"
        ":type in: `Data`")
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
        "returns the solution *u* for a time step *dt>0* with initial value u0\n\n"
        ":rtype: `Data`\n"
        ":param source:\n"
        ":type source: `Data`")
     .def("insertConstraint",&escript::AbstractTransportProblem::insertConstraint,
args("source", "q", "r","factor"),
        "inserts constraint *u_{,t}=r* where *q>0*  into the problem using a weighting factor")
     .def("reset",&escript::AbstractTransportProblem::resetTransport,
        "resets the transport operator typically as they have been updated.")
     .def("resetValues",&escript::AbstractTransportProblem::resetTransport)
     .def("getSafeTimeStepSize",&escript::AbstractTransportProblem::getSafeTimeStepSize)
     .def("getUnlimitedTimeStepSize",&escript::AbstractTransportProblem::getUnlimitedTimeStepSize);

  enum_<escript::SolverOptions>("SolverOptions")
    .value("DEFAULT", escript::SO_DEFAULT)

    .value("MKL", escript::SO_PACKAGE_MKL)
    .value("PASO", escript::SO_PACKAGE_PASO)
    .value("TRILINOS", escript::SO_PACKAGE_TRILINOS)
    .value("UMFPACK", escript::SO_PACKAGE_UMFPACK)
    .value("MUMPS", escript::SO_PACKAGE_MUMPS)

    .value("BICGSTAB", escript::SO_METHOD_BICGSTAB)
    .value("CGLS", escript::SO_METHOD_CGLS)
    .value("CGS", escript::SO_METHOD_CGS)
    .value("CHOLEVSKY", escript::SO_METHOD_CHOLEVSKY)
    .value("CR", escript::SO_METHOD_CR)
    .value("DIRECT", escript::SO_METHOD_DIRECT)
    .value("DIRECT_MUMPS", escript::SO_METHOD_DIRECT_MUMPS)
    .value("DIRECT_PARDISO", escript::SO_METHOD_DIRECT_PARDISO)
    .value("DIRECT_SUPERLU", escript::SO_METHOD_DIRECT_SUPERLU)
    .value("DIRECT_TRILINOS", escript::SO_METHOD_DIRECT_TRILINOS)
    .value("GMRES", escript::SO_METHOD_GMRES)
    .value("HRZ_LUMPING", escript::SO_METHOD_HRZ_LUMPING)
    .value("ITERATIVE", escript::SO_METHOD_ITERATIVE)
    .value("LSQR", escript::SO_METHOD_LSQR)
    .value("LUMPING", escript::SO_METHOD_ROWSUM_LUMPING)
    .value("MINRES", escript::SO_METHOD_MINRES)
    .value("NONLINEAR_GMRES", escript::SO_METHOD_NONLINEAR_GMRES)
    .value("PCG", escript::SO_METHOD_PCG)
    .value("PRES20", escript::SO_METHOD_PRES20)
    .value("ROWSUM_LUMPING", escript::SO_METHOD_ROWSUM_LUMPING)
    .value("TFQMR", escript::SO_METHOD_TFQMR)

    .value("AMG", escript::SO_PRECONDITIONER_AMG)
    .value("GAUSS_SEIDEL", escript::SO_PRECONDITIONER_GAUSS_SEIDEL)
    .value("ILU0", escript::SO_PRECONDITIONER_ILU0)
    .value("ILUT", escript::SO_PRECONDITIONER_ILUT)
    .value("JACOBI", escript::SO_PRECONDITIONER_JACOBI)
    .value("NO_PRECONDITIONER", escript::SO_PRECONDITIONER_NONE)
    .value("REC_ILU", escript::SO_PRECONDITIONER_REC_ILU)
    .value("RILU", escript::SO_PRECONDITIONER_RILU)

    .value("BACKWARD_EULER", escript::SO_ODESOLVER_BACKWARD_EULER)
    .value("CRANK_NICOLSON", escript::SO_ODESOLVER_CRANK_NICOLSON)
    .value("LINEAR_CRANK_NICOLSON", escript::SO_ODESOLVER_LINEAR_CRANK_NICOLSON)

    .value("CLASSIC_INTERPOLATION", escript::SO_INTERPOLATION_CLASSIC)
    .value("CLASSIC_INTERPOLATION_WITH_FF_COUPLING", escript::SO_INTERPOLATION_CLASSIC_WITH_FF_COUPLING)
    .value("DIRECT_INTERPOLATION", escript::SO_INTERPOLATION_DIRECT)

    .value("DEFAULT_REORDERING", escript::SO_REORDERING_DEFAULT)
    .value("MINIMUM_FILL_IN", escript::SO_REORDERING_MINIMUM_FILL_IN)
    .value("NESTED_DISSECTION", escript::SO_REORDERING_NESTED_DISSECTION)
    .value("NO_REORDERING", escript::SO_REORDERING_NONE);


  class_<escript::SolverBuddy, escript::SB_ptr >("SolverBuddy","",init<>())
    .def("setTrilinosParameter", &escript::SolverBuddy::setTrilinosParameter,
            "Sets a Trilinos preconditioner/solver parameter.\n\n"
        ":note Escript does not check for validity of the parameter name\n"
        "(e.g. spelling mistakes). Parameters are passed 1:1 to escript's\n"
        "Trilinos wrapper and from there to the relevant Trilinos package.\n"
        "See the relevant Trilinos documentation for valid parameter strings\n"
        "and values."
        ":note This method does nothing in a non-Trilinos build.")
    .def("getTrilinosParameters", &escript::SolverBuddy::getTrilinosParameters,
            "Returns a dictionary of set Trilinos parameters.\n\n"
        ":note This method returns an empty dictionary in a non-Trilinos build.")
    .def("getSummary", &escript::SolverBuddy::getSummary,"Returns a string reporting the current settings")
    .def("__str__", &escript::SolverBuddy::getSummary)
    .def("getName", &escript::SolverBuddy::getName, args("key"),"Returns the name of a given key\n\n"
        ":param key: a valid key")
    .def("resetDiagnostics", &escript::SolverBuddy::resetDiagnostics, args("all")=false,"Resets the diagnostics\n\n"
        ":param all: if ``all`` is ``True`` all diagnostics including accumulative counters are reset.\n"
        ":type all: ``bool``")
    .def("_updateDiagnostics", &escript::SolverBuddy::updateDiagnosticsPy, args("key", "value"),"Updates diagnostic information\n\n"
        ":param name: name of  diagnostic information\n"
        ":type name: ``str`` in the list 'num_iter', 'num_level',\n"
        "'num_inner_iter', 'time', 'set_up_time', 'net_time',\n"
        "'residual_norm', 'converged'.\n"
        ":param value: new value of the diagnostic information\n"
        ":note: this function is used by a solver to report diagnostics\n"
        "information.")
    .def("getDiagnostics", &escript::SolverBuddy::getDiagnostics, args("name"),"Returns the diagnostic information ``name``. Possible values are:\n\n"
        "- 'num_iter': the number of iteration steps\n"
        "- 'cum_num_iter': the cumulative number of iteration steps\n"
        "- 'num_level': the number of level in multi level solver\n"
        "- 'num_inner_iter': the number of inner iteration steps\n"
        "- 'cum_num_inner_iter': the cumulative number of inner iteration steps\n"
        "- 'time': execution time\n"
        "- 'cum_time': cumulative execution time\n"
        "- 'set_up_time': time to set up of the solver, typically this includes factorization and reordering\n"
        "- 'cum_set_up_time': cumulative time to set up of the solver\n"
        "- 'net_time': net execution time, excluding setup time for the solver and execution time for preconditioner\n"
        "- 'cum_net_time': cumulative net execution time\n"
        "- 'preconditioner_size': size of preconditioner [Bytes]\n"
        "- 'converged': return True if solution has converged.\n"
        "- 'time_step_backtracking_used': returns True if time step back tracking has been used.\n"
        "- 'coarse_level_sparsity': returns the sparsity of the matrix on the coarsest level\n"
        "- 'num_coarse_unknowns': returns the number of unknowns on the coarsest level\n\n\n"
        ":param name: name of diagnostic information to return\n"
        ":type name: ``str`` in the list above.\n"
        ":return: requested value. 0 is returned if the value is yet to be defined.\n"
        ":note: If the solver has thrown an exception diagnostic values have an undefined status.")
    .def("hasConverged", &escript::SolverBuddy::hasConverged,"Returns ``True`` if the last solver call has been finalized successfully.\n\n"
        ":note: if an exception has been thrown by the solver the status of this"
        "flag is undefined.\n")
    .def("setPreconditioner", &escript::SolverBuddy::setPreconditioner, args("preconditioner"),"Sets the preconditioner to be used.\n\n"
        ":param preconditioner: key of the preconditioner to be used.\n"
        ":type preconditioner: in `ILU0`, `ILUT`, `JACOBI`, `AMG`, , `REC_ILU`, `GAUSS_SEIDEL`, `RILU`, `NO_PRECONDITIONER`\n"
        ":note: Not all packages support all preconditioner. It can be assumed that a package makes a reasonable choice if it encounters an unknown"
        "preconditioner.\n")
    .def("getPreconditioner", &escript::SolverBuddy::getPreconditioner,"Returns the key of the preconditioner to be used.\n\n"
        ":rtype: in the list `ILU0`, `ILUT`, `JACOBI`, `AMG`, `REC_ILU`, `GAUSS_SEIDEL`, `RILU`,  `NO_PRECONDITIONER`")
    .def("setSolverMethod", &escript::SolverBuddy::setSolverMethod, args("method"),"Sets the solver method to be used. Use ``method``=``DIRECT`` to indicate that a direct rather than an iterative solver should be used and use ``method``=``ITERATIVE`` to indicate that an iterative rather than a direct solver should be used.\n\n"
        ":param method: key of the solver method to be used.\n"
        ":type method: in `DEFAULT`, `DIRECT`, `CHOLEVSKY`, `PCG`, `CR`, `CGS`, `BICGSTAB`, `GMRES`, `PRES20`, `ROWSUM_LUMPING`, `HRZ_LUMPING`, `ITERATIVE`, `NONLINEAR_GMRES`, `TFQMR`, `MINRES`\n"
        ":note: Not all packages support all solvers. It can be assumed that a package makes a reasonable choice if it encounters an unknown solver method.")
    .def("getSolverMethod", &escript::SolverBuddy::getSolverMethod,"Returns key of the solver method to be used.\n\n"
        ":rtype: in the list `DEFAULT`, `DIRECT`, `CHOLEVSKY`, `PCG`, `CR`, `CGS`, `BICGSTAB`, `GMRES`, `PRES20`, `ROWSUM_LUMPING`, `HRZ_LUMPING`, `MINRES`, `ITERATIVE`, `NONLINEAR_GMRES`, `TFQMR`")
    .def("setPackage", &escript::SolverBuddy::setPackage, args("package"),"Sets the solver package to be used as a solver.\n\n"
        ":param package: key of the solver package to be used.\n"
        ":type package: in `DEFAULT`, `PASO`, `CUSP`, `MKL`, `UMFPACK`, `MUMPS`, `TRILINOS`\n"
        ":note: Not all packages are support on all implementation. An exception may be thrown on some platforms if a particular package is requested.")
    .def("getPackage", &escript::SolverBuddy::getPackage,"Returns the solver package key\n\n"
        ":rtype: in the list `DEFAULT`, `PASO`, `CUSP`, `MKL`, `UMFPACK`, `MUMPS`, `TRILINOS`")
    .def("setReordering", &escript::SolverBuddy::setReordering, args("ordering"),"Sets the key of the reordering method to be applied if supported by the solver. Some direct solvers support reordering to optimize compute time and storage use during elimination.\n\n"
        ":param ordering: selects the reordering strategy.\n"
        ":type ordering: in 'NO_REORDERING', 'MINIMUM_FILL_IN', 'NESTED_DISSECTION', 'DEFAULT_REORDERING'")
    .def("getReordering", &escript::SolverBuddy::getReordering,"Returns the key of the reordering method to be applied if supported by the solver.\n\n"
        ":rtype: in `NO_REORDERING`, `MINIMUM_FILL_IN`, `NESTED_DISSECTION`, `DEFAULT_REORDERING`")
    .def("setRestart", &escript::SolverBuddy::setRestart, args("restart"),"Sets the number of iterations steps after which GMRES performs a restart.\n\n"
        ":param restart: number of iteration steps after which to perform a restart. If 0 no restart is performed.\n"
        ":type restart: ``int``")
    .def("getRestart", &escript::SolverBuddy::getRestart,"Returns the number of iterations steps after which GMRES performs a restart. If 0 is returned no restart is performed.\n\n"
        ":rtype: ``int``")
    .def("setTruncation", &escript::SolverBuddy::setTruncation, args("truncation"),"Sets the number of residuals in GMRES to be stored for orthogonalization. The more residuals are stored the faster GMRES converged\n\n"
        ":param truncation: truncation\n"
        ":type truncation: ``int``")
    .def("getTruncation", &escript::SolverBuddy::getTruncation,"Returns the number of residuals in GMRES to be stored for orthogonalization\n\n"
        ":rtype: ``int``")
    .def("setInnerIterMax", &escript::SolverBuddy::setInnerIterMax, args("iter_max"),"Sets the maximum number of iteration steps for the inner iteration.\n\n"
        ":param iter_max: maximum number of inner iterations\n"
        ":type iter_max: ``int``")
    .def("getInnerIterMax", &escript::SolverBuddy::getInnerIterMax,"Returns maximum number of inner iteration steps\n\n"
        ":rtype: ``int``")
    .def("setIterMax", &escript::SolverBuddy::setIterMax, args("iter_max"),"Sets the maximum number of iteration steps\n\n"
        ":param iter_max: maximum number of iteration steps\n"
        ":type iter_max: ``int``")
    .def("getIterMax", &escript::SolverBuddy::getIterMax,"Returns maximum number of iteration steps\n\n"
        ":rtype: ``int``")
    .def("setNumSweeps", &escript::SolverBuddy::setNumSweeps, args("sweeps"),"Sets the number of sweeps in a Jacobi or Gauss-Seidel/SOR preconditioner.\n\n"
        ":param sweeps: number of sweeps\n"
        ":type sweeps: positive ``int``")
    .def("getNumSweeps", &escript::SolverBuddy::getNumSweeps,"Returns the number of sweeps in a Jacobi or Gauss-Seidel/SOR preconditioner.\n\n"
        ":rtype: ``int``")
    .def("setTolerance", &escript::SolverBuddy::setTolerance, args("rtol"),"Sets the relative tolerance for the solver\n\n"
        ":param rtol: relative tolerance\n"
        ":type rtol: non-negative ``float``")
    .def("getTolerance", &escript::SolverBuddy::getTolerance,"Returns the relative tolerance for the solver\n\n"
        ":rtype: ``float``")
    .def("setAbsoluteTolerance", &escript::SolverBuddy::setAbsoluteTolerance, args("atol"),"Sets the absolute tolerance for the solver\n\n"
        ":param atol:  absolute tolerance\n"
        ":type atol: non-negative ``float``")
    .def("getAbsoluteTolerance", &escript::SolverBuddy::getAbsoluteTolerance,"Returns the absolute tolerance for the solver\n\n"
        ":rtype: ``float``")
    .def("setInnerTolerance", &escript::SolverBuddy::setInnerTolerance, args("rtol"),"Sets the relative tolerance for an inner iteration scheme, for instance on the coarsest level in a multi-level scheme.\n\n"
        ":param rtol: inner relative tolerance\n"
        ":type rtol: positive ``float``")
    .def("getInnerTolerance", &escript::SolverBuddy::getInnerTolerance,"Returns the relative tolerance for an inner iteration scheme\n\n"
        ":rtype: ``float``")
    .def("setDropTolerance", &escript::SolverBuddy::setDropTolerance, args("drop_tol"),"Sets the relative drop tolerance in ILUT\n\n"
        ":param drop_tol: drop tolerance\n"
        ":type drop_tol: positive ``float``")
    .def("getDropTolerance", &escript::SolverBuddy::getDropTolerance,"Returns the relative drop tolerance in ILUT\n\n"
        ":rtype: ``float``")
    .def("setDropStorage", &escript::SolverBuddy::setDropStorage, args("drop"),"Sets the maximum allowed increase in storage for ILUT. ``storage`` =2 would mean that a doubling of the storage needed for the coefficient matrix is\n"
        "allowed in the ILUT factorization.\n\n"
        ":param storage: allowed storage increase\n"
        ":type storage: ``float``")
    .def("getDropStorage", &escript::SolverBuddy::getDropStorage,"Returns the maximum allowed increase in storage for ILUT\n\n"
        ":rtype: ``float``")
    .def("setRelaxationFactor", &escript::SolverBuddy::setRelaxationFactor, args("relaxation"),"Sets the relaxation factor used to add dropped elements in RILU to the main diagonal.\n\n"
        ":param factor: relaxation factor\n"
        ":type factor: ``float``\n"
        ":note: RILU with a relaxation factor 0 is identical to ILU0")
    .def("getRelaxationFactor", &escript::SolverBuddy::getRelaxationFactor,"Returns the relaxation factor used to add dropped elements in RILU to the main diagonal.\n\n"
        ":rtype: ``float``")
    .def("isComplex", &escript::SolverBuddy::isComplex,"Checks if the coefficient matrix is set to be complex-valued.\n\n"
        ":return: True if a complex-valued PDE is indicated, False otherwise\n"
        ":rtype: ``bool``")
    .def("setComplex", &escript::SolverBuddy::setComplex, args("complex"),"Sets the complex flag for the coefficient matrix to ``flag``.\n\n"
        ":param flag: If True, the complex flag is set otherwise reset.\n"
        ":type flag: ``bool``")
    .def("setDim", &escript::SolverBuddy::setDim, args("dim"),"Sets the dimension of the problem.\n\n"
        ":param dim: Either 2 or 3.\n"
        ":rtype: ``int``")
    .def("getDim", &escript::SolverBuddy::getDim, "Returns the dimension of the problem.\n\n"
        ":rtype: ``int``")
    .def("isSymmetric", &escript::SolverBuddy::isSymmetric,"Checks if symmetry of the coefficient matrix is indicated.\n\n"
        ":return: True if a symmetric PDE is indicated, False otherwise\n"
        ":rtype: ``bool``")
    .def("setSymmetryOn", &escript::SolverBuddy::setSymmetryOn,"Sets the symmetry flag to indicate that the coefficient matrix is symmetric.")
    .def("setSymmetryOff", &escript::SolverBuddy::setSymmetryOff,"Clears the symmetry flag for the coefficient matrix.")
    .def("setSymmetry", &escript::SolverBuddy::setSymmetry, args("symmetry"),"Sets the symmetry flag for the coefficient matrix to ``flag``.\n\n"
        ":param flag: If True, the symmetry flag is set otherwise reset.\n"
        ":type flag: ``bool``")
    .def("isHermitian", &escript::SolverBuddy::isHermitian,"Checks if the coefficient matrix is indicated to be Hermitian.\n\n"
        ":return: True if a hermitian PDE is indicated, False otherwise\n"
        ":rtype: ``bool``")
    .def("setHermitianOn", &escript::SolverBuddy::setHermitianOn,"Sets the hermitian flag to indicate that the coefficient matrix is hermitian.")
    .def("setHermitianOff", &escript::SolverBuddy::setHermitianOff,"Clears the hermitian flag for the coefficient matrix.")
    .def("setHermitian", &escript::SolverBuddy::setHermitian, args("hermitian"),"Sets the hermitian flag for the coefficient matrix to ``flag``.\n\n"
        ":param flag: If True, the hermitian flag is set otherwise reset.\n"
        ":type flag: ``bool``")
    // .def("useDirectSolver",&escript::SolverBuddy::useDirect,":rtype: `int`\n"
    //     ":sets a Direct Solver")
    .def("isVerbose", &escript::SolverBuddy::isVerbose,"Returns ``True`` if the solver is expected to be verbose.\n\n"
        ":return: True if verbosity of switched on.\n"
        ":rtype: ``bool``")
    .def("setVerbosityOn", &escript::SolverBuddy::setVerbosityOn,"Switches the verbosity of the solver on.")
    .def("setVerbosityOff", &escript::SolverBuddy::setVerbosityOff,"Switches the verbosity of the solver off.")
    .def("setVerbosity", &escript::SolverBuddy::setVerbosity, args("verbosity"),"Sets the verbosity flag for the solver to ``flag``.\n\n"
        ":param verbose: If ``True``, the verbosity of the solver is switched on.\n"
        ":type verbose: ``bool``")
    .def("adaptInnerTolerance", &escript::SolverBuddy::adaptInnerTolerance,"Returns ``True`` if the tolerance of the inner solver is selected automatically. Otherwise the inner tolerance set by `setInnerTolerance` is used.\n\n"
        ":return: ``True`` if inner tolerance adaption is chosen.\n"
        ":rtype: ``bool``")
    .def("setInnerToleranceAdaptionOn", &escript::SolverBuddy::setInnerToleranceAdaptionOn,"Switches the automatic selection of inner tolerance on")
    .def("setInnerToleranceAdaptionOff", &escript::SolverBuddy::setInnerToleranceAdaptionOff,"Switches the automatic selection of inner tolerance off.")
    .def("setInnerToleranceAdaption", &escript::SolverBuddy::setInnerToleranceAdaption, args("adapt"),"Sets the flag to indicate automatic selection of the inner tolerance.\n\n"
        ":param adapt: If ``True``, the inner tolerance is selected automatically.\n"
        ":type adapt: ``bool``")
    .def("acceptConvergenceFailure", &escript::SolverBuddy::acceptConvergenceFailure,"Returns ``True`` if a failure to meet the stopping criteria within the given number of iteration steps is not raising in exception. This is useful\n"
        "if a solver is used in a non-linear context where the non-linear solver can continue even if the returned the solution does not necessarily meet the stopping criteria. One can use the `hasConverged` method to check if the\n"
        "last call to the solver was successful.\n\n"
        ":return: ``True`` if a failure to achieve convergence is accepted.\n"
        ":rtype: ``bool``")
    .def("setAcceptanceConvergenceFailureOn", &escript::SolverBuddy::setAcceptanceConvergenceFailureOn,"Switches the acceptance of a failure of convergence on")
    .def("setAcceptanceConvergenceFailureOff", &escript::SolverBuddy::setAcceptanceConvergenceFailureOff,"Switches the acceptance of a failure of convergence off.")
    .def("setAcceptanceConvergenceFailure", &escript::SolverBuddy::setAcceptanceConvergenceFailure, args("accept"),"Sets the flag to indicate the acceptance of a failure of convergence.\n\n"
        ":param accept: If ``True``, any failure to achieve convergence is accepted.\n"
        ":type accept: ``bool``")
    .def("useLocalPreconditioner", &escript::SolverBuddy::useLocalPreconditioner,"Returns ``True`` if the preconditoner is applied locally on each MPI. This reduces communication costs and speeds up the application of the preconditioner but at the costs of more iteration steps. This can be an advantage on clusters with slower interconnects.\n\n"
        ":return: ``True`` if local preconditioning is applied\n"
        ":rtype: ``bool``")
    .def("setLocalPreconditionerOn", &escript::SolverBuddy::setLocalPreconditionerOn,"Sets the flag to use  local preconditioning to on")
    .def("setLocalPreconditionerOff", &escript::SolverBuddy::setLocalPreconditionerOff,"Sets the flag to use  local preconditioning to off")
    .def("setLocalPreconditioner", &escript::SolverBuddy::setLocalPreconditioner, args("local"),"Sets the flag to use  local preconditioning\n\n"
        ":param use: If ``True``, local preconditioning on each MPI rank is applied\n"
        ":type use: ``bool``")
    .def("setNumRefinements", &escript::SolverBuddy::setNumRefinements, args("refinements"),"Sets the number of refinement steps to refine the solution when a direct solver is applied.\n\n"
        ":param refinements: number of refinements\n"
        ":type refinements: non-negative ``int``")
    .def("getNumRefinements", &escript::SolverBuddy::getNumRefinements,"Returns the number of refinement steps to refine the solution when a direct solver is applied.\n\n"
        ":rtype: non-negative ``int``")
    .def("setOxleyDomain", &escript::SolverBuddy::setOxleyDomain, args("using_oxley"), "Sets the parameter Oxley_Domain.\n\n"
        ":rtype: non-negative ``int``")
    .def("getOxleyDomain", &escript::SolverBuddy::getOxleyDomain,"True if we are using an Oxley domain, False otherwise.\n\n"
        ":rtype: non-negative ``int``")
    .def("setODESolver", &escript::SolverBuddy::setODESolver, args("solver"),"Set the solver method for ODEs.\n\n"
        ":param method: key of the ODE solver method to be used.\n"
        ":type method: in `CRANK_NICOLSON`, `BACKWARD_EULER`, `LINEAR_CRANK_NICOLSON`")
    .def("getODESolver", &escript::SolverBuddy::getODESolver,"Returns key of the solver method for ODEs.\n\n"
        ":param method: key of the ODE solver method to be used.\n"
        ":type method: in `CRANK_NICOLSON`, `BACKWARD_EULER`, `LINEAR_CRANK_NICOLSON`");


  // Functions to get/modify global parameters/features
  def("setEscriptParamInt", escript::setEscriptParamInt,
      (arg("name"), arg("value")=0), "Modify the value of an escript tuning parameter\n\n"
        ":param name:\n"
        ":type name: ``string``\n"
        ":param value:\n"
        ":type value: ``int``");

  def("getEscriptParamInt", escript::getEscriptParamInt,
      (arg("name"),arg("sentinel")=0), "Read the value of an escript tuning parameter\n\n"
        ":param name: parameter to lookup\n"
        ":type name: ``string``\n"
        ":param sentinel: Value to be returned if ``name`` is not a known parameter\n"
        ":type sentinel: ``int``");

  def("listEscriptParams", escript::listEscriptParams,
        ":return: A list of tuples (p,v,d) where p is the name of a parameter "
        "for escript, v is its current value, and d is a description.");

  def("hasFeature", escript::hasFeature,
      (arg("name")), "Check if escript was compiled with a certain feature\n\n"
        ":param name: feature to lookup\n"
        ":type name: ``string``");
  def("listFeatures", escript::listFeatures,
        ":return: A list of strings representing the features escript supports.");

  def("resolveGroup", escript::resolveGroup);

#ifdef IKNOWWHATIMDOING
  def("applyBinaryCFunction", escript::applyBinaryCFunction,
          (arg("function"), arg("outshape"), arg("in1"), arg("in2"))
  );
#endif

  def("_condEval", escript::condEval, (arg("mask"), arg("trueval"), arg("falseval")));
}
