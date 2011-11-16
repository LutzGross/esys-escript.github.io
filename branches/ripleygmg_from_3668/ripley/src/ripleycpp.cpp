
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include <ripley/RipleyDomainFactory.h>
#include <esysUtils/esysExceptionTranslator.h>

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/detail/defaults_gen.hpp>
#include <boost/version.hpp>

using namespace boost::python;

/**
   \page ripley Ripley
   Ripley is the python module name that contains the interfaces
   to the C++ wrapper to ripley.

   \version 1.0.0 

   \section class_desc Class Description:
   None

   \section class_limits Class Limitations:
   None

   \section class_conds Class Conditions of Use:
   None

   \section throws Throws:
   None
*/

BOOST_PYTHON_MODULE(ripleycpp)
{
// This feature was added in boost v1.34
#if ((BOOST_VERSION/100)%1000 > 34) || (BOOST_VERSION/100000 >1)
  // params are: bool show_user_defined, bool show_py_signatures, bool show_cpp_signatures
  docstring_options docopt(true, true, false);
#endif

  register_exception_translator<ripley::RipleyException>(&(esysUtils::esysExceptionTranslator));

  def("LoadMesh",ripley::loadMesh, (arg("fileName")), ":rtype: `Domain`");

  def("ReadMesh",ripley::readMesh, (arg("fileName"), arg("optimize")=true)
        ,"Reads a mesh from a file. For MPI parallel runs the mesh is fanned out to multiple processes.\n\n"
        ":rtype: `Domain`\n:param fileName:\n:type fileName: ``string``\n"
        ":param optimize: Whether to optimize node labels\n:type optimize: ``bool``");

  def("ReadGmsh",ripley::readGmsh,
      (arg("fileName"), arg("numDim"), arg("optimize")=true)
        ,"Reads a gmsh mesh file\n\n"
        ":rtype: `Domain`\n:param fileName:\n:type fileName: ``string``\n"
        ":param optimize: Whether to optimize node labels\n:type optimize: ``bool``\n");

  def ("Brick", ripley::brick,
      (arg("n0")=1, arg("n1")=1, arg("n2")=1,
      arg("l0")=1.0, arg("l1")=1.0, arg("l2")=1.0,
      arg("optimize")=false)
        ,"Creates a rectangular mesh with n0 x n1 x n2 elements over the brick [0,l0] x [0,l1] x [0,l2]."
        "\n\n:param n0:\n:type n0:\n:param n1:\n:type n1:\n:param n2:\n:type n2:\n"
        ":param l0: length of side 0\n:param l1:\n:param l2:\n"
        ":param optimize: Whether to optimize node labels\n:type optimize: ``bool``");

  def ("Rectangle", ripley::rectangle,
      (arg("n0")=1, arg("n1")=1, arg("l0")=1.0, arg("l1")=1.0,
      arg("optimize")=false)
        ,"Creates a rectangular mesh with n0 x n1 elements over the rectangle [0,l0] x [0,l1]."
        "\n\n:param n0:\n:type n0:\n:param n1:\n:type n1:\n"
        ":param l0: length of side 0\n:param l1:\n"
        ":param optimize: Whether to optimize node labels\n:type optimize: ``bool``");

  class_<ripley::RipleyDomain, bases<escript::AbstractContinuousDomain> >
      ("RipleyDomain","A concrete class representing a domain. For more details, please consult the c++ documentation.", no_init)
      .def(init<const ripley::RipleyDomain&>())
      .def("write",&ripley::RipleyDomain::write,args("filename"),
"Write the current mesh to a file with the given name.")
      .def("print_mesh_info",&ripley::RipleyDomain::Print_Mesh_Info,(arg("full")=false),
":param full:\n:type full: ``bool``")
      .def("dump",&ripley::RipleyDomain::dump,args("fileName")
,"dumps the mesh to a file with the given name.")
      .def("getDescription",&ripley::RipleyDomain::getDescription,
":return: a description for this domain\n:rtype: ``string``")
      .def("getDim",&ripley::RipleyDomain::getDim,":rtype: ``int``")
      .def("getDataShape",&ripley::RipleyDomain::getDataShape, args("functionSpaceCode"),
":return: a pair (dps, ns) where dps=the number of data points per sample, and ns=the number of samples\n:rtype: ``tuple``")
      .def("getNumDataPointsGlobal",&ripley::RipleyDomain::getNumDataPointsGlobal,
":return: the number of data points summed across all MPI processes\n"
":rtype: ``int``")
      .def("addPDEToSystem",&ripley::RipleyDomain::addPDEToSystem,
args("mat", "rhs","A", "B", "C", "D", "X", "Y", "d", "y", "d_contact", "y_contact"),
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
)
      .def("addPDEToLumpedSystem",&ripley::RipleyDomain::addPDEToLumpedSystem,
args("mat", "D", "d"),
"adds a PDE onto the lumped stiffness matrix\n\n"
":param mat:\n:type mat: `Data`\n"
":param D:\n:type D: `Data`\n"
":param d:\n:type d: `Data`\n"
":param useHRZ:\n:type useHRZ: bool\n"
)
      .def("addPDEToRHS",&ripley::RipleyDomain::addPDEToRHS, 
args("rhs", "X", "Y", "y", "y_contact"),
"adds a PDE onto the stiffness matrix mat and a rhs\n\n"
":param rhs:\n:type rhs: `Data`\n"
":param X:\n:type X: `Data`\n"
":param Y:\n:type Y: `Data`\n"
":param y:\n:type y: `Data`\n"
":param y_contact:\n:type y_contact: `Data`"
)
      .def("addPDEToTransportProblem",&ripley::RipleyDomain::addPDEToTransportProblem,
args( "tp", "source", "M", "A", "B", "C", "D", "X", "Y", "d", "y", "d_contact", "y_contact"),
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
)
      .def("newOperator",&ripley::RipleyDomain::newSystemMatrix,
args("row_blocksize", "row_functionspace", "column_blocksize", "column_functionspace", "type"),
"creates a SystemMatrixAdapter stiffness matrix and initializes it with zeros\n\n"
":param row_blocksize:\n:type row_blocksize: ``int``\n"
":param row_functionspace:\n:type row_functionspace: `FunctionSpace`\n"
":param column_blocksize:\n:type column_blocksize: ``int``\n"
":param column_functionspace:\n:type column_functionspace: `FunctionSpace`\n"
":param type:\n:type type: ``int``\n"
)
      .def("newTransportProblem",&ripley::RipleyDomain::newTransportProblem,
args("theta", "blocksize", "functionspace", "type"),
"creates a TransportProblemAdapter\n\n"
":param theta:\n:type theta: ``float``\n"
":param blocksize:\n:type blocksize: ``int``\n"
":param functionspace:\n:type functionspace: `FunctionSpace`\n"
":param type:\n:type type: ``int``\n"
)
      .def("getSystemMatrixTypeId",&ripley::RipleyDomain::getSystemMatrixTypeId,
args("solver", "preconditioner", "package", "symmetry"),
":return: the identifier of the matrix type to be used for the global stiffness matrix when a particular solver, package, perconditioner, and symmetric matrix is used.\n"
":rtype: ``int``\n"
":param solver:\n:type solver: ``int``\n"
":param preconditioner:\n:type preconditioner: ``int``\n"
":param package:\n:type package: ``int``\n"
":param symmetry:\n:type symmetry: ``int``\n"
)
      .def("getTransportTypeId",&ripley::RipleyDomain::getTransportTypeId,
args("solver", "preconditioner", "package", "symmetry"),
":return: the identifier of the transport problem type to be used when a particular solver, perconditioner, package and symmetric matrix is used.\n"
":rtype: ``int``\n"
":param solver:\n:type solver: ``int``\n"
":param preconditioner:\n:type preconditioner: ``int``\n"
":param package:\n:type package: ``int``\n"
":param symmetry:\n:type symmetry: ``int``\n"
)
      .def("setX",&ripley::RipleyDomain::setNewX,
args("arg"), "assigns new location to the domain\n\n:param arg:\n:type arg: `Data`")
      .def("getX",&ripley::RipleyDomain::getX, ":return: locations in the FEM nodes\n\n"
":rtype: `Data`")
      .def("getNormal",&ripley::RipleyDomain::getNormal,
":return: boundary normals at the quadrature point on the face elements\n"
":rtype: `Data`")
      .def("getSize",&ripley::RipleyDomain::getSize,":return: the element size\n"
":rtype: `Data`")
      .def("setTagMap",&ripley::RipleyDomain::setTagMap,args("name","tag"),
"Give a tag number a name.\n\n:param name: Name for the tag\n:type name: ``string``\n"
":param tag: numeric id\n:type tag: ``int``\n:note: Tag names must be unique within a domain")
      .def("getTag",&ripley::RipleyDomain::getTag,args("name"),":return: tag id for "
"``name``\n:rtype: ``string``")
      .def("isValidTagName",&ripley::RipleyDomain::isValidTagName,args("name"),
":return: True is ``name`` corresponds to a tag\n:rtype: ``bool``")
      .def("showTagNames",&ripley::RipleyDomain::showTagNames,":return: A space separated list of tag names\n:rtype: ``string``")
      .def("getMPISize",&ripley::RipleyDomain::getMPISize,":return: the number of processes used for this `Domain`\n:rtype: ``int``")
      .def("getMPIRank",&ripley::RipleyDomain::getMPIRank,":return: the rank of this process\n:rtype: ``int``")
      .def("MPIBarrier",&ripley::RipleyDomain::MPIBarrier,"Wait until all processes have reached this point")
      .def("onMasterProcessor",&ripley::RipleyDomain::onMasterProcessor,":return: True if this code is executing on the master process\n:rtype: `bool`");

  class_<ripley::SystemMatrixAdapter, bases<escript::AbstractSystemMatrix> >
      ("OperatorAdapter","A concrete class representing an operator. For more details, please see the c++ documentation.", no_init)
      .def("print_matrix_info",&ripley::SystemMatrixAdapter::Print_Matrix_Info,(arg("full")=false),"prints information about a system matrix")
      .def("nullifyRowsAndCols",&ripley::SystemMatrixAdapter::nullifyRowsAndCols)
      .def("resetValues",&ripley::SystemMatrixAdapter::resetValues, "resets the matrix entries")
      .def("saveMM",&ripley::SystemMatrixAdapter::saveMM,args("fileName"), 
"writes the matrix to a file using the Matrix Market file format")
      .def("saveHB",&ripley::SystemMatrixAdapter::saveHB, args("filename"),
"writes the matrix to a file using the Harwell-Boeing file format");

  class_<ripley::TransportProblemAdapter, bases<escript::AbstractTransportProblem> >
      ("TransportProblemAdapter","",no_init)
      .def("getSafeTimeStepSize",&ripley::TransportProblemAdapter::getSafeTimeStepSize)
      .def("getUnlimitedTimeStepSize",&ripley::TransportProblemAdapter::getUnlimitedTimeStepSize)
      .def("resetTransport",&ripley::TransportProblemAdapter::resetTransport,
"resets the transport operator typically as they have been updated");
}

