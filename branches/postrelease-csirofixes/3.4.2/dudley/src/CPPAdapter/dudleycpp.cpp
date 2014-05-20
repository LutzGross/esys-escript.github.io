
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


#ifdef ESYS_MPI
#include "esysUtils/Esys_MPI.h"
#endif
#include "../Dudley.h"

#include <pasowrap/SystemMatrixAdapter.h>
#include <pasowrap/TransportProblemAdapter.h>

#include "MeshAdapter.h"
#include "MeshAdapterFactory.h"

#include "DudleyAdapterException.h"
// #include "esysUtils/EsysException.h"
#include "esysUtils/esysExceptionTranslator.h"

#include "escript/AbstractContinuousDomain.h"

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
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
  
  
  //
  // NOTE: The return_value_policy is necessary for functions that
  // return pointers.
  //
  register_exception_translator<dudley::DudleyAdapterException>(&(esysUtils::RuntimeErrorTranslator));

  def("LoadMesh",dudley::loadMesh,
      (arg("fileName")="file.nc"),":rtype: `Domain`"
/*      ,return_value_policy<manage_new_object>());*/
      );

  def("ReadMesh",dudley::readMesh,
      (arg("fileName")="file.fly",arg("integrationOrder")=-1,  arg("reducedIntegrationOrder")=-1,  arg("optimize")=true)
/*      ,return_value_policy<manage_new_object>());*/
	,"Read a mesh from a file. For MPI parallel runs fan out the mesh to multiple processes.\n\n"
":rtype: `Domain`\n:param fileName:\n:type fileName: ``string``\n"
":param integrationOrder: order of the quadrature scheme. If *integrationOrder<0* the integration order is selected independently.\n"
":type integrationOrder: ``int``\n"
":param reducedIntegrationOrder: order of the quadrature scheme. If *reducedIntegrationOrder<0* the integration order is selected independently.\n"
":param optimize: Enable optimisation of node labels\n:type optimize: ``bool``");

  def("ReadGmsh",dudley::readGmsh,
      (arg("fileName")="file.msh",
       arg("numDim"), 
       arg("integrationOrder")=-1, 
       arg("reducedIntegrationOrder")=-1, 
       arg("optimize")=true,  
       arg("useMacroElements")=false)
//       ,return_value_policy<manage_new_object>());
,"Read a gmsh mesh file\n\n"
":rtype: `Domain`\n:param fileName:\n:type fileName: ``string``\n"
":param integrationOrder: order of the quadrature scheme. If *integrationOrder<0* the integration order is selected independently.\n"
":type integrationOrder: ``int``\n"
":param reducedIntegrationOrder: order of the quadrature scheme. If *reducedIntegrationOrder<0* the integration order is selected independently.\n"
":param optimize: Enable optimisation of node labels\n:type optimize: ``bool``\n"
":param useMacroElements: Enable the usage of macro elements instead of second order elements.\n:type useMacroElements: ``bool``"
);

  def ("Brick",dudley::brick,
      (arg("n0")=1,arg("n1")=1,arg("n2")=1,
      arg("order")=1,
      arg("l0")=1.0,arg("l1")=1.0,arg("l2")=1.0,
      arg("periodic0")=false,arg("periodic1")=false,arg("periodic2")=false,
      arg("integrationOrder")=-1,  arg("reducedIntegrationOrder")=-1,
      arg("useElementsOnFace")=false,
      arg("useFullElementOrder")=false,
      arg("optimize")=false)
,"Creates a tetrahedral mesh by subdividing n0 x n1 x n2 rectangular elements over the brick [0,l0] x [0,l1] x [0,l2]."
"We accept floating point values for n0, n1 only to ease transition of scripts to python3 when the time comes."
"\n\n:param n0:\n:type n0:\n:param n1:\n:type n1:\n:param n2:\n:type n2:\n"
":param order: =1, =-1 or =2 gives the order of shape function. If -1 macro elements of order 1 are used.\n"
":param l0: length of side 0\n:param l1:\n:param l2:\n"
":param integrationOrder: order of the quadrature scheme. If integrationOrder<0 the integration order is selected independently.\n"
":param reducedIntegrationOrder: order of the quadrature scheme. If reducedIntegrationOrder<0 the integration order is selected independently.\n"
":param useElementsOnFace:  whether or not to use elements on face\n"
":type useElementsOnFace: ``int``"
":param periodic0:  whether or not boundary conditions are periodic\n"
":param periodic1:\n:param periodic2:\n"
":param useFullElementOrder:\n:param optimize:\n"
":param optimize: Enable optimisation of node labels\n:type optimize: ``bool``"
);

  def ("Rectangle",dudley::rectangle,
      (arg("n0")=1,arg("n1")=1,arg("order")=1,
      arg("l0")=1.0,arg("l1")=1.0,
      arg("periodic0")=false,arg("periodic1")=false,
      arg("integrationOrder")=-1,  arg("reducedIntegrationOrder")=-1,
      arg("useElementsOnFace")=false,
      arg("useFullElementOrder")=false,
      arg("optimize")=false)
,"Creates a triangular mesh by subdividing n0 x n1 rectangular elements over the brick [0,l0] x [0,l1]."
"We accept floating point values for n0, n1 only to ease transition of scripts to python3 when the time comes."
"\n\n:param n0:\n:type n0:\n:param n1:\n:type n1:\n"
":param order: =1, =-1 or =2 gives the order of shape function. If -1 macro elements of order 1 are used.\n"
":param l0: length of side 0\n:param l1:\n"
":param integrationOrder: order of the quadrature scheme. If integrationOrder<0 the integration order is selected independently.\n"
":param reducedIntegrationOrder: order of the quadrature scheme. If reducedIntegrationOrder<0 the integration order is selected independently.\n"
":param useElementsOnFace:  whether or not to use elements on face\n"
":type useElementsOnFace: ``int``"
":param periodic0:  whether or not boundary conditions are periodic\n"
":param periodic1:\n"
":param useFullElementOrder:\n:param optimize:\n"
":param useMacroElements: Enable the usage of first order macro elements.\n:type useMacroElements: ``bool``\n"
":param optimize: Enable optimisation of node labels\n:type optimize: ``bool``"
);

  class_<dudley::MeshAdapter, bases<escript::AbstractContinuousDomain> >
      ("MeshAdapter","A concrete class representing a domain. For more details, please consult the c++ documentation.",init<optional <Dudley_Mesh*> >())
      .def(init<const dudley::MeshAdapter&>())
      .def("write",&dudley::MeshAdapter::write,args("filename"),
"Write the current mesh to a file with the given name.")
      .def("print_mesh_info",&dudley::MeshAdapter::Print_Mesh_Info,(arg("full")=false),
":param full:\n:type full: ``bool``")
      .def("dump",&dudley::MeshAdapter::dump,args("fileName")
,"dumps the mesh to a file with the given name.")
      .def("getDescription",&dudley::MeshAdapter::getDescription,
":return: a description for this domain\n:rtype: ``string``")
      .def("getDim",&dudley::MeshAdapter::getDim,":rtype: ``int``")
      .def("getDataShape",&dudley::MeshAdapter::getDataShape, args("functionSpaceCode"),
":return: a pair (dps, ns) where dps=the number of data points per sample, and ns=the number of samples\n:rtype: ``tuple``")
      .def("getNumDataPointsGlobal",&dudley::MeshAdapter::getNumDataPointsGlobal,
":return: the number of data points summed across all MPI processes\n"
":rtype: ``int``")
      .def("addPDEToSystem",&dudley::MeshAdapter::addPDEToSystem,
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
      .def("addPDEToLumpedSystem",&dudley::MeshAdapter::addPDEToLumpedSystem,
args("mat", "D", "d"),
"adds a PDE onto the lumped stiffness matrix\n\n"
":param mat:\n:type mat: `Data`\n"
":param D:\n:type D: `Data`\n"
":param d:\n:type d: `Data`\n"
":param useHRZ:\n:type useHRZ: bool\n"
)
      .def("addPDEToRHS",&dudley::MeshAdapter::addPDEToRHS, 
args("rhs", "X", "Y", "y", "y_contact"),
"adds a PDE onto the stiffness matrix mat and a rhs\n\n"
":param rhs:\n:type rhs: `Data`\n"
":param X:\n:type X: `Data`\n"
":param Y:\n:type Y: `Data`\n"
":param y:\n:type y: `Data`\n"
":param y_contact:\n:type y_contact: `Data`"
)
      .def("addPDEToTransportProblem",&dudley::MeshAdapter::addPDEToTransportProblem,
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
      .def("newOperator",&dudley::MeshAdapter::newSystemMatrix,
args("row_blocksize", "row_functionspace", "column_blocksize", "column_functionspace", "type"),
"creates a SystemMatrixAdapter stiffness matrix and initializes it with zeros\n\n"
":param row_blocksize:\n:type row_blocksize: ``int``\n"
":param row_functionspace:\n:type row_functionspace: `FunctionSpace`\n"
":param column_blocksize:\n:type column_blocksize: ``int``\n"
":param column_functionspace:\n:type column_functionspace: `FunctionSpace`\n"
":param type:\n:type type: ``int``\n"
)
      .def("newTransportProblem",&dudley::MeshAdapter::newTransportProblem,
args("theta", "blocksize", "functionspace", "type"),
"creates a TransportProblemAdapter\n\n"
":param theta:\n:type theta: ``float``\n"
":param blocksize:\n:type blocksize: ``int``\n"
":param functionspace:\n:type functionspace: `FunctionSpace`\n"
":param type:\n:type type: ``int``\n"
)
      .def("getSystemMatrixTypeId",&dudley::MeshAdapter::getSystemMatrixTypeId,
args("solver", "preconditioner", "package", "symmetry"),
":return: the identifier of the matrix type to be used for the global stiffness matrix when a particular solver, package, preconditioner, and symmetric matrix is used.\n"
":rtype: ``int``\n"
":param solver:\n:type solver: ``int``\n"
":param preconditioner:\n:type preconditioner: ``int``\n"
":param package:\n:type package: ``int``\n"
":param symmetry:\n:type symmetry: ``int``\n"
)
      .def("getTransportTypeId",&dudley::MeshAdapter::getTransportTypeId,
args("solver", "preconditioner", "package", "symmetry"),
":return: the identifier of the transport problem type to be used when a particular solver, preconditioner, package and symmetric matrix is used.\n"
":rtype: ``int``\n"
":param solver:\n:type solver: ``int``\n"
":param preconditioner:\n:type preconditioner: ``int``\n"
":param package:\n:type package: ``int``\n"
":param symmetry:\n:type symmetry: ``int``\n"
)
      .def("setX",&dudley::MeshAdapter::setNewX,
args("arg"), "assigns new location to the domain\n\n:param arg:\n:type arg: `Data`")
      .def("getX",&dudley::MeshAdapter::getX, ":return: locations in the FEM nodes\n\n"
":rtype: `Data`")
      .def("getNormal",&dudley::MeshAdapter::getNormal,
":return: boundary normals at the quadrature point on the face elements\n"
":rtype: `Data`")
      .def("getSize",&dudley::MeshAdapter::getSize,":return: the element size\n"
":rtype: `Data`")
      .def("setTagMap",&dudley::MeshAdapter::setTagMap,args("name","tag"),
"Give a tag number a name.\n\n:param name: Name for the tag\n:type name: ``string``\n"
":param tag: numeric id\n:type tag: ``int``\n:note: Tag names must be unique within a domain")
      .def("getTag",&dudley::MeshAdapter::getTag,args("name"),":return: tag id for "
"``name``\n:rtype: ``string``")
      .def("isValidTagName",&dudley::MeshAdapter::isValidTagName,args("name"),
":return: True is ``name`` corresponds to a tag\n:rtype: ``bool``")
      .def("showTagNames",&dudley::MeshAdapter::showTagNames,":return: A space separated list of tag names\n:rtype: ``string``")
      .def("getMPISize",&dudley::MeshAdapter::getMPISize,":return: the number of processes used for this `Domain`\n:rtype: ``int``")
      .def("getMPIRank",&dudley::MeshAdapter::getMPIRank,":return: the rank of this process\n:rtype: ``int``")
      .def("MPIBarrier",&dudley::MeshAdapter::MPIBarrier,"Wait until all processes have reached this point")
      .def("onMasterProcessor",&dudley::MeshAdapter::onMasterProcessor,":return: True if this code is executing on the master process\n:rtype: `bool`");

}
