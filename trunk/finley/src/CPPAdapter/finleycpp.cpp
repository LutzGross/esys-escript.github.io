
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


#ifdef PASO_MPI
#include <mpi.h>
#include "paso/Paso_MPI.h"
#endif
extern "C" {
#include "../Finley.h"
}

#include "MeshAdapter.h"
#include "MeshAdapterFactory.h"
#include "SystemMatrixAdapter.h"
#include "TransportProblemAdapter.h"

#include "FinleyAdapterException.h"
// #include "esysUtils/EsysException.h"
#include "esysUtils/esysExceptionTranslator.h"

#include "escript/AbstractContinuousDomain.h"

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/detail/defaults_gen.hpp>
#include <boost/version.hpp>

using namespace boost::python;

/**
   \page finley Finley
   Finley is the python module name that contains the interfaces
   to the C++ wrapper to finley.

   \version 1.0.0 

   \section class_desc Class Description:
   Data

   \section class_limits Class Limitations:
   None

   \section class_conds Class Conditions of Use:
   None

   \section throws Throws:
   None

*/

//
// The BOOST_PYTHON_FUNCTION_OVERLOADS macro generates function overloads for optional
// arguments to the respective finley functions.
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// NOTE: If the number of arguments to the finley functions change
// the magic numbers in the BOOST_PYTHON_FUNCTION_OVERLOADS call 
// must change.
//
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// BOOST_PYTHON_FUNCTION_OVERLOADS(readMesh_overloads,finley::readMesh,1,2)
// BOOST_PYTHON_FUNCTION_OVERLOADS(brick_overloads,finley::brick,0,12)
// BOOST_PYTHON_FUNCTION_OVERLOADS(rectangle_overloads,finley::rectangle,0,9)
// BOOST_PYTHON_FUNCTION_OVERLOADS(interval_overloads,finley::interval,0,6)
// BOOST_PYTHON_FUNCTION_OVERLOADS(glueFaces_overloads,finley::glueFaces,1,3)
// BOOST_PYTHON_FUNCTION_OVERLOADS(joinFaces_overloads,finley::joinFaces,1,3)

BOOST_PYTHON_MODULE(finleycpp)
{
// This feature was added in boost v1.34
#if ((BOOST_VERSION/100)%1000 > 34) || (BOOST_VERSION/100000 >1)
  // params are: bool show_user_defined, bool show_py_signatures, bool show_cpp_signatures
  docstring_options docopt(true, true, false);
#endif

  //
  // NOTE: The return_value_policy is necessary for functions that
  // return pointers.
  //
  register_exception_translator<finley::FinleyAdapterException>(&(esysUtils::esysExceptionTranslator));

  def("LoadMesh",finley::loadMesh,
      (arg("fileName")="file.nc"),":rtype: `Domain`"
/*      ,return_value_policy<manage_new_object>());*/
      );

  def("ReadMesh",finley::readMesh,
      (arg("fileName")="file.fly",arg("integrationOrder")=-1,  arg("reducedIntegrationOrder")=-1,  arg("optimize")=true)
/*      ,return_value_policy<manage_new_object>());*/
	,"Read a mesh from a file. For MPI parallel runs fan out the mesh to multiple processes.\n\n"
":rtype: `Domain`\n:param fileName:\n:type fileName: ``string``\n"
":param integrationOrder: order of the quadrature scheme. If *integrationOrder<0* the integration order is selected independently.\n"
":type integrationOrder: ``int``\n"
":param reducedIntegrationOrder: order of the quadrature scheme. If *reducedIntegrationOrder<0* the integration order is selected independently.\n"
":param optimize: Enable optimisation of node labels\n:type optimize: ``bool``");

  def("ReadGmsh",finley::readGmsh,
      (arg("fileName")="file.msh",arg("numDim"), arg("integrationOrder")=-1, arg("reducedIntegrationOrder")=-1, arg("optimize")=true)
//       ,return_value_policy<manage_new_object>());
,"Read a gmsh mesh file\n\n"
":rtype: `Domain`\n:param fileName:\n:type fileName: ``string``\n"
":param integrationOrder: order of the quadrature scheme. If *integrationOrder<0* the integration order is selected independently.\n"
":type integrationOrder: ``int``\n"
":param reducedIntegrationOrder: order of the quadrature scheme. If *reducedIntegrationOrder<0* the integration order is selected independently.\n"
":param optimize: Enable optimisation of node labels\n:type optimize: ``bool``");

  def ("Brick",finley::brick,
      (arg("n0")=1,arg("n1")=1,arg("n2")=1,
      arg("order")=1,
      arg("l0")=1.0,arg("l1")=1.0,arg("l2")=1.0,
      arg("periodic0")=false,arg("periodic1")=false,arg("periodic2")=false,
      arg("integrationOrder")=-1,  arg("reducedIntegrationOrder")=-1,
      arg("useElementsOnFace")=false,
      arg("useFullElementOrder")=false,
      arg("optimize")=false)
//       ,return_value_policy<manage_new_object>());
,"Creates a rectangular mesh with n0 x n1 x n2 elements over the brick [0,l0] x [0,l1] x [0,l2]."
"\n\n:param n0:\n:type n0:\n:param n1:\n:type n1:\n:param n2:\n:type n2:\n"
":param order: =1 or =2 gives the order of shape function\n"
":param l0: length of side 0\n:param l1:\n:param l2:\n"
":param integrationOrder: order of the quadrature scheme. If integrationOrder<0 the integration order is selected independently.\n"
":param reducedIntegrationOrder: order of the quadrature scheme. If reducedIntegrationOrder<0 the integration order is selected independently.\n"
":param useElementsOnFace:  whether or not to use elements on face\n"
":type useElementsOnFace: ``int``"
":param periodic0:  whether or not boundary conditions are periodic\n"
":param periodic1:\n:param periodic2:\n"
":param useFullElementOrder:\n:param optimize:\n"
);

  def ("Rectangle",finley::rectangle,
      (arg("n0")=1,arg("n1")=1,arg("order")=1,
      arg("l0")=1.0,arg("l1")=1.0,
      arg("periodic0")=false,arg("periodic1")=false,
      arg("integrationOrder")=-1,  arg("reducedIntegrationOrder")=-1,
      arg("useElementsOnFace")=false,
      arg("useFullElementOrder")=false,
      arg("optimize")=false)
//       ,return_value_policy<manage_new_object>());
,"Creates a rectangular mesh with n0 x n1 elements over the brick [0,l0] x [0,l1]."
"\n\n:param n0:\n:type n0:\n:param n1:\n:type n1:\n"
":param order: =1 or =2 gives the order of shape function\n"
":param l0: length of side 0\n:param l1:\n"
":param integrationOrder: order of the quadrature scheme. If integrationOrder<0 the integration order is selected independently.\n"
":param reducedIntegrationOrder: order of the quadrature scheme. If reducedIntegrationOrder<0 the integration order is selected independently.\n"
":param useElementsOnFace:  whether or not to use elements on face\n"
":type useElementsOnFace: ``int``"
":param periodic0:  whether or not boundary conditions are periodic\n"
":param periodic1:\n"
":param useFullElementOrder:\n:param optimize:\n");

  def("Merge",finley::meshMerge,args("meshList")
//       ,return_value_policy<manage_new_object>());
,"Merges a list of meshes into one mesh.\n\n:rtype: `Domain`"
  );

  def("GlueFaces",finley::glueFaces,
      (arg("meshList"),arg("safetyFactor")=0.2,
      arg("tolerance")=1.e-8,
      arg("optimize")=true)
//       ,return_value_policy<manage_new_object>());
,"Detects matching faces in the mesh, removes them from the mesh and joins the elements touched by the face elements."
	);

  def("JoinFaces",finley::joinFaces,
      (arg("meshList"), arg("safetyFactor")=0.2,
      arg("tolerance")=1.e-8,
      arg("optimize")=true)
//       ,return_value_policy<manage_new_object>());
,"Detects matching faces in the mesh and replaces them by joint elements."
	);


  class_<finley::MeshAdapter, bases<escript::AbstractContinuousDomain> >
      ("MeshAdapter","A concrete class representing a domain. For more details, please consult the c++ documentation.",init<optional <Finley_Mesh*> >())
      .def(init<const finley::MeshAdapter&>())
      .def("write",&finley::MeshAdapter::write,args("filename"),
"Write the current mesh to a file with the given name.")
      .def("print_mesh_info",&finley::MeshAdapter::Print_Mesh_Info,(arg("full")=false),
":param full:\n:type full: ``bool``")
      .def("dump",&finley::MeshAdapter::dump,args("fileName")
,"dumps the mesh to a file with the given name.")
      .def("getDescription",&finley::MeshAdapter::getDescription,
":return: a description for this domain\n:rtype: ``string``")
      .def("getDim",&finley::MeshAdapter::getDim,":rtype: ``int``")
      .def("getDataShape",&finley::MeshAdapter::getDataShape, args("functionSpaceCode"),
":return: a pair (dps, ns) where dps=the number of data points per sample, and ns=the number of samples\n:rtype: ``tuple``")
      .def("getNumDataPointsGlobal",&finley::MeshAdapter::getNumDataPointsGlobal,
":return: the number of data points summed across all MPI processes\n"
":rtype: ``int``")
      .def("addPDEToSystem",&finley::MeshAdapter::addPDEToSystem,
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
      .def("addPDEToLumpedSystem",&finley::MeshAdapter::addPDEToLumpedSystem,
args("mat", "D", "d"),
"adds a PDE onto the lumped stiffness matrix\n\n"
":param mat:\n:type mat: `Data`\n"
":param D:\n:type D: `Data`\n"
":param d:\n:type d: `Data`\n")
      .def("addPDEToRHS",&finley::MeshAdapter::addPDEToRHS, 
args("rhs", "X", "Y", "y", "y_contact"),
"adds a PDE onto the stiffness matrix mat and a rhs\n\n"
":param rhs:\n:type rhs: `Data`\n"
":param X:\n:type X: `Data`\n"
":param Y:\n:type Y: `Data`\n"
":param y:\n:type y: `Data`\n"
":param y_contact:\n:type y_contact: `Data`"
)
      .def("addPDEToTransportProblem",&finley::MeshAdapter::addPDEToTransportProblem,
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
      .def("newOperator",&finley::MeshAdapter::newSystemMatrix,
args("row_blocksize", "row_functionspace", "column_blocksize", "column_functionspace", "type"),
"creates a SystemMatrixAdapter stiffness matrix and initializes it with zeros\n\n"
":param row_blocksize:\n:type row_blocksize: ``int``\n"
":param row_functionspace:\n:type row_functionspace: `FunctionSpace`\n"
":param column_blocksize:\n:type column_blocksize: ``int``\n"
":param column_functionspace:\n:type column_functionspace: `FunctionSpace`\n"
":param type:\n:type type: ``int``\n"
)
      .def("newTransportProblem",&finley::MeshAdapter::newTransportProblem,
args("theta", "blocksize", "functionspace", "type"),
"creates a TransportProblemAdapter\n\n"
":param theta:\n:type theta: ``float``\n"
":param blocksize:\n:type blocksize: ``int``\n"
":param functionspace:\n:type functionspace: `FunctionSpace`\n"
":param type:\n:type type: ``int``\n"
)
      .def("getSystemMatrixTypeId",&finley::MeshAdapter::getSystemMatrixTypeId,
args("solver", "preconditioner", "package", "symmetry"),
":return: the identifier of the matrix type to be used for the global stiffness matrix when a particular solver, package, perconditioner, and symmetric matrix is used.\n"
":rtype: ``int``\n"
":param solver:\n:type solver: ``int``\n"
":param preconditioner:\n:type preconditioner: ``int``\n"
":param package:\n:type package: ``int``\n"
":param symmetry:\n:type symmetry: ``int``\n"
)
      .def("getTransportTypeId",&finley::MeshAdapter::getTransportTypeId,
args("solver", "preconditioner", "package", "symmetry"),
":return: the identifier of the transport problem type to be used when a particular solver, perconditioner, package and symmetric matrix is used.\n"
":rtype: ``int``\n"
":param solver:\n:type solver: ``int``\n"
":param preconditioner:\n:type preconditioner: ``int``\n"
":param package:\n:type package: ``int``\n"
":param symmetry:\n:type symmetry: ``int``\n"
)
      .def("setX",&finley::MeshAdapter::setNewX,
args("arg"), "assigns new location to the domain\n\n:param arg:\n:type arg: `Data`")
      .def("getX",&finley::MeshAdapter::getX, ":return: locations in the FEM nodes\n\n"
":rtype: `Data`")
      .def("getNormal",&finley::MeshAdapter::getNormal,
":return: boundary normals at the quadrature point on the face elements\n"
":rtype: `Data`")
      .def("getSize",&finley::MeshAdapter::getSize,":return: the element size\n"
":rtype: `Data`")
      .def("saveDX",&finley::MeshAdapter::saveDX,args("filename" ,"arg"),
"Saves a dictonary of Data objects to an OpenDX input file. The keywords are used as identifier"
"\n\n:param filename: \n:type filename: ``string``\n"
"\n:param arg: \n:type arg: ``dict``\n")
      .def("saveVTK",&finley::MeshAdapter::saveVTK,
args("filename" ,"arg",  "metadata", "metadata_schema"),
"Saves a dictonary of Data objects to an VTK XML input file. The keywords are used as identifier"
"\n\n:param filename:\n:type filename: ``string``\n"
":param arg:\n:type arg: ``dict``\n"
":param metadata:\n:type metadata: ``string``\n"
":param metadata_schema:\n:type metadata_schema: ``string``\n"
)
      .def("setTagMap",&finley::MeshAdapter::setTagMap,args("name","tag"),
"Give a tag number a name.\n\n:param name: Name for the tag\n:type name: ``string``\n"
":param tag: numeric id\n:type tag: ``int``\n:note: Tag names must be unique within a domain")
      .def("getTag",&finley::MeshAdapter::getTag,args("name"),":return: tag id for "
"``name``\n:rtype: ``string``")
      .def("isValidTagName",&finley::MeshAdapter::isValidTagName,args("name"),
":return: True is ``name`` corresponds to a tag\n:rtype: ``bool``")
      .def("showTagNames",&finley::MeshAdapter::showTagNames,":return: A space separated list of tag names\n:rtype: ``string``")
      .def("getMPISize",&finley::MeshAdapter::getMPISize,":return: the number of processes used for this `Domain`\n:rtype: ``int``")
      .def("getMPIRank",&finley::MeshAdapter::getMPIRank,":return: the rank of this process\n:rtype: ``int``")
      .def("MPIBarrier",&finley::MeshAdapter::MPIBarrier,"Wait until all processes have reached this point")
      .def("onMasterProcessor",&finley::MeshAdapter::onMasterProcessor,":return: True if this code is executing on the master process\n:rtype: `bool`");

  class_<finley::SystemMatrixAdapter, bases<escript::AbstractSystemMatrix> >
      ("OperatorAdapter","A concrete class representing an operator. For more details, please see the c++ documentation.", no_init)
      .def("print_matrix_info",&finley::SystemMatrixAdapter::Print_Matrix_Info,(arg("full")=false),"prints information about a system matrix")
      .def("nullifyRowsAndCols",&finley::SystemMatrixAdapter::nullifyRowsAndCols)
      .def("resetValues",&finley::SystemMatrixAdapter::resetValues, "resets the matrix entries")
      .def("saveMM",&finley::SystemMatrixAdapter::saveMM,args("fileName"), 
"writes the matrix to a file using the Matrix Market file format")
      .def("saveHB",&finley::SystemMatrixAdapter::saveHB, args("filename"),
"writes the matrix to a file using the Harwell-Boeing file format");

  class_<finley::TransportProblemAdapter, bases<escript::AbstractTransportProblem> >
      ("TransportProblemAdapter","",no_init)
      .def("getSafeTimeStepSize",&finley::TransportProblemAdapter::getSafeTimeStepSize)
      .def("getUnlimitedTimeStepSize",&finley::TransportProblemAdapter::getUnlimitedTimeStepSize)
      .def("resetTransport",&finley::TransportProblemAdapter::resetTransport,
"resets the transport operator typically as they have been updated");
}
