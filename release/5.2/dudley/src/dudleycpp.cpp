
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include <dudley/Dudley.h>
#include <dudley/DomainFactory.h>
#include <dudley/DudleyDomain.h>

#include <escript/ExceptionTranslators.h>

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/module.hpp>
#include <boost/python/detail/defaults_gen.hpp>
#include <boost/version.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(dudleycpp)
{
// This feature was added in boost v1.34
#if ((BOOST_VERSION/100)%1000 > 34) || (BOOST_VERSION/100000 >1)
    // params are: bool show_user_defined, bool show_py_signatures, bool show_cpp_signatures
    docstring_options docopt(true, true, false);
#endif

    scope().attr("__doc__") = "To use this module, please import esys.dudley";

    // register escript's default translators
    REGISTER_ESCRIPT_EXCEPTION_TRANSLATORS;
    register_exception_translator<dudley::DudleyException>(&escript::RuntimeErrorTranslator);

  def("LoadMesh", dudley::DudleyDomain::load,
      (arg("fileName") = "file.nc"), ":rtype: `DudleyDomain`");

  def("ReadMesh", dudley::readMesh,
      (arg("fileName")="file.fly", arg("integrationOrder")=-1, arg("reducedIntegrationOrder")=-1, arg("optimize")=true)
	,"Read a mesh from a fly file. For MPI parallel runs fan out the mesh to multiple processes.\n\n"
":rtype: `Domain`\n:param fileName:\n:type fileName: ``string``\n"
":param integrationOrder: order of the quadrature scheme. Ignored.\n"
":type integrationOrder: ``int``\n"
":param reducedIntegrationOrder: order of reduced quadrature scheme. Ignored.\n"
":param optimize: Enable optimisation of node labels\n:type optimize: ``bool``");

  def("ReadGmsh", dudley::readGmsh,
      (arg("fileName") = "file.msh",
       arg("numDim"), 
       arg("integrationOrder") = -1, 
       arg("reducedIntegrationOrder") = -1, 
       arg("optimize") = true)
,"Read a gmsh mesh file\n\n"
":rtype: `Domain`\n:param fileName:\n:type fileName: ``string``\n"
":param integrationOrder: order of the quadrature scheme. Always 2.\n"
":type integrationOrder: ``int``\n"
":param reducedIntegrationOrder: order of reduced quadrature scheme. Always 0.\n"
":param optimize: Enable optimisation of node labels\n:type optimize: ``bool``\n");

  def ("__Brick_driver", dudley::brick_driver, arg("args"));
  def ("__Rectangle_driver", dudley::rectangle_driver, arg("args"));

  class_<dudley::DudleyDomain, bases<escript::AbstractContinuousDomain> >
      ("DudleyDomain", "A concrete class representing a dudley domain. For more details, please consult the C++ documentation.", no_init)
      .def(init<const dudley::DudleyDomain&>())
      .def("write", &dudley::DudleyDomain::write, args("filename"),
"Write the current mesh to a file with the given name.")
      .def("print_mesh_info", &dudley::DudleyDomain::Print_Mesh_Info, (arg("full")=false),
":param full:\n:type full: ``bool``")
      .def("dump", &dudley::DudleyDomain::dump, args("fileName")
,"dumps the mesh to a file with the given name.")
      .def("getDescription", &dudley::DudleyDomain::getDescription,
":return: a description for this domain\n:rtype: ``string``")
      .def("getDim", &dudley::DudleyDomain::getDim,":rtype: ``int``")
      .def("getDataShape", &dudley::DudleyDomain::getDataShape, args("functionSpaceCode"),
":return: a pair (dps, ns) where dps=the number of data points per sample, and ns=the number of samples\n:rtype: ``tuple``")
      .def("getNumDataPointsGlobal", &dudley::DudleyDomain::getNumDataPointsGlobal,
":return: the number of data points summed across all MPI processes\n"
":rtype: ``int``")
      .def("addPDEToSystem", &dudley::DudleyDomain::addPDEToSystem,
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
      .def("addPDEToLumpedSystem", &dudley::DudleyDomain::addPDEToLumpedSystem,
args("mat", "D", "d"),
"adds a PDE onto the lumped stiffness matrix\n\n"
":param mat:\n:type mat: `Data`\n"
":param D:\n:type D: `Data`\n"
":param d:\n:type d: `Data`\n"
":param useHRZ:\n:type useHRZ: bool\n"
)
      .def("addPDEToRHS", &dudley::DudleyDomain::addPDEToRHS, 
args("rhs", "X", "Y", "y", "y_contact"),
"adds a PDE onto the stiffness matrix mat and a rhs\n\n"
":param rhs:\n:type rhs: `Data`\n"
":param X:\n:type X: `Data`\n"
":param Y:\n:type Y: `Data`\n"
":param y:\n:type y: `Data`\n"
":param y_contact:\n:type y_contact: `Data`"
)
      .def("addPDEToTransportProblem", &dudley::DudleyDomain::addPDEToTransportProblem,
args( "tp", "source", "M", "A", "B", "C", "D", "X", "Y", "d", "y", "d_contact", "y_contact"),
":param tp:\n:type tp: `AbstractTransportProblem`\n"
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
      .def("newOperator", &dudley::DudleyDomain::newSystemMatrix,
args("row_blocksize", "row_functionspace", "column_blocksize", "column_functionspace", "type"),
"creates a stiffness matrix and initializes it with zeros\n\n"
":param row_blocksize:\n:type row_blocksize: ``int``\n"
":param row_functionspace:\n:type row_functionspace: `FunctionSpace`\n"
":param column_blocksize:\n:type column_blocksize: ``int``\n"
":param column_functionspace:\n:type column_functionspace: `FunctionSpace`\n"
":param type:\n:type type: ``int``\n"
)
      .def("newTransportProblem", &dudley::DudleyDomain::newTransportProblem,
args("theta", "blocksize", "functionspace", "type"),
"creates a TransportProblem\n\n"
":param theta:\n:type theta: ``float``\n"
":param blocksize:\n:type blocksize: ``int``\n"
":param functionspace:\n:type functionspace: `FunctionSpace`\n"
":param type:\n:type type: ``int``\n"
)
      .def("getSystemMatrixTypeId", &dudley::DudleyDomain::getSystemMatrixTypeId,
args("options"),
":return: the identifier of the matrix type to be used for the global stiffness matrix when particular solver options are used.\n"
":rtype: ``int``\n"
":param options:\n:type options: `SolverBuddy`\n"
)
      .def("getTransportTypeId", &dudley::DudleyDomain::getTransportTypeId,
args("solver", "preconditioner", "package", "symmetry"),
":return: the identifier of the transport problem type to be used when a particular solver, preconditioner, package and symmetric matrix is used.\n"
":rtype: ``int``\n"
":param solver:\n:type solver: ``int``\n"
":param preconditioner:\n:type preconditioner: ``int``\n"
":param package:\n:type package: ``int``\n"
":param symmetry:\n:type symmetry: ``int``\n"
)
      .def("setX", &dudley::DudleyDomain::setNewX,
args("arg"), "assigns new location to the domain\n\n:param arg:\n:type arg: `Data`")
      .def("getX", &dudley::DudleyDomain::getX, ":return: locations in the FEM nodes\n\n"
":rtype: `Data`")
      .def("getNormal", &dudley::DudleyDomain::getNormal,
":return: boundary normals at the quadrature point on the face elements\n"
":rtype: `Data`")
      .def("getSize", &dudley::DudleyDomain::getSize,":return: the element size\n"
":rtype: `Data`")
      .def("setTagMap", &dudley::DudleyDomain::setTagMap,args("name","tag"),
"Give a tag number a name.\n\n:param name: Name for the tag\n:type name: ``string``\n"
":param tag: numeric id\n:type tag: ``int``\n:note: Tag names must be unique within a domain")
      .def("getTag", &dudley::DudleyDomain::getTag,args("name"),":return: tag id for "
"``name``\n:rtype: ``string``")
      .def("isValidTagName", &dudley::DudleyDomain::isValidTagName,args("name"),
":return: True is ``name`` corresponds to a tag\n:rtype: ``bool``")
      .def("showTagNames", &dudley::DudleyDomain::showTagNames,":return: A space separated list of tag names\n:rtype: ``string``")
      .def("getMPISize", &dudley::DudleyDomain::getMPISize,":return: the number of processes used for this `Domain`\n:rtype: ``int``")
      .def("getMPIRank", &dudley::DudleyDomain::getMPIRank,":return: the rank of this process\n:rtype: ``int``")
      .def("MPIBarrier", &dudley::DudleyDomain::MPIBarrier,"Wait until all processes have reached this point")
      .def("onMasterProcessor", &dudley::DudleyDomain::onMasterProcessor,":return: True if this code is executing on the master process\n:rtype: `bool`");
}

