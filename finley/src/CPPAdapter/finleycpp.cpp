
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "finley/Finley.h"
#include "finley/FinleyException.h"
#include "MeshAdapter.h"
#include "MeshAdapterFactory.h"

#include <escript/AbstractContinuousDomain.h>
#include <escript/ExceptionTranslators.h>

#include "esysUtils/Esys_MPI.h"

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/detail/defaults_gen.hpp>
#include <boost/version.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(finleycpp)
{
// This feature was added in boost v1.34
#if ((BOOST_VERSION/100)%1000 > 34) || (BOOST_VERSION/100000 >1)
  // params are: bool show_user_defined, bool show_py_signatures, bool show_cpp_signatures
  docstring_options docopt(true, true, false);
#endif

  scope().attr("__doc__") = "To use this module, please import esys.finley";    

    // register escript's default translators
    REGISTER_ESCRIPT_EXCEPTION_TRANSLATORS;
    register_exception_translator<finley::FinleyException>(&escript::RuntimeErrorTranslator);

  //
  // NOTE: The return_value_policy is necessary for functions that
  // return pointers.
  //
  def("LoadMesh",finley::loadMesh,
      (arg("fileName")="file.nc"),":rtype: `Domain`"
/*      ,return_value_policy<manage_new_object>());*/
      );

  
  def("__ReadMesh_driver", finley::readMesh_driver,
      (arg("params"))
//  def("ReadMesh",finley::readMesh,
//      (arg("fileName")="file.fly",arg("integrationOrder")=-1,  arg("reducedIntegrationOrder")=-1,  arg("optimize")=true)
	,"Read a mesh from a file. For MPI parallel runs fan out the mesh to multiple processes.\n\n"
":rtype: `Domain`\n:param fileName:\n:type fileName: ``string``\n"
":param integrationOrder: order of the quadrature scheme. If *integrationOrder<0* the integration order is selected independently.\n"
":type integrationOrder: ``int``\n"
":param reducedIntegrationOrder: order of the quadrature scheme. If *reducedIntegrationOrder<0* the integration order is selected independently.\n"
":param optimize: Enable optimisation of node labels\n:type optimize: ``bool``");

  def("__ReadGmsh_driver", finley::readGmsh_driver,
      (arg("params"))  
//  def("ReadGmsh",finley::readGmsh,
//      (arg("fileName")="file.msh",
//       arg("numDim"), 
//       arg("integrationOrder")=-1, 
//       arg("reducedIntegrationOrder")=-1, 
//       arg("optimize")=true,  
//       arg("useMacroElements")=false)
,"Read a gmsh mesh file\n\n"
":rtype: `Domain`\n:param fileName:\n:type fileName: ``string``\n"
":param integrationOrder: order of the quadrature scheme. If *integrationOrder<0* the integration order is selected independently.\n"
":type integrationOrder: ``int``\n"
":param reducedIntegrationOrder: order of the quadrature scheme. If *reducedIntegrationOrder<0* the integration order is selected independently.\n"
":param optimize: Enable optimisation of node labels\n:type optimize: ``bool``\n"
":param useMacroElements: Enable the usage of macro elements instead of second order elements.\n:type useMacroElements: ``bool``"
);

  def ("__Brick_driver",finley::brick_driver,
      (arg("params"))
,"Creates a rectangular mesh with n0 x n1 x n2 elements over the brick [0,l0] x [0,l1] x [0,l2]."
"\n\n:param n0: number of elements in direction 0\n:type n0: ``int``\n:param n1: number of elements in direction 1\n:type n1: ``int``\n"
":param n2:number of elements in direction 2\n:type n2: ``int``\n"
":param order: =1, =-1 or =2 gives the order of shape function. If -1 macro elements of order 1 are used.\n"
":param l0: length of side 0\n"
":type  l0: ``float``\n"
":param l1: length of side 1\n"
":type  l1: ``float``\n"
":param l2: length of side 2\n"
":type  l2: ``float``\n"
":param periodic0: whether or not boundary conditions are periodic in direction 0\n:type periodic0: ``bool``\n"
":param periodic1: whether or not boundary conditions are periodic in direction 1\n:type periodic1: ``bool``\n"
":param periodic2: whether or not boundary conditions are periodic in direction 2\n:type periodic2: ``bool``\n"
":param integrationOrder: order of the quadrature scheme. If integrationOrder<0 the integration order is selected independently.\n"
":param reducedIntegrationOrder: order of the quadrature scheme. If reducedIntegrationOrder<0 the integration order is selected independently.\n"
":param useElementsOnFace:  whether or not to use elements on face\n"
":type useElementsOnFace: ``int``\n"
":param useFullElementOrder: Whether or not to use Hex27 elements\n"":type useFullElementOrder: ``bool``\n"
":param optimize: Enable optimisation of node labels\n:type optimize: ``bool``"
);

  def ("__Rectangle_driver",finley::rectangle_driver,
      (arg("args")) 
,"Creates a rectangular mesh with n0 x n1 elements over the brick [0,l0] x [0,l1]."
"\n\n:param n0:\n:type n0:\n:param n1:\n:type n1:\n"
":param order: =1, =-1 or =2 gives the order of shape function. If -1 macro elements of order 1 are used.\n"
":param l0: length of side 0\n:param l1:\n"
":param integrationOrder: order of the quadrature scheme. If integrationOrder<0 the integration order is selected independently.\n"
":param reducedIntegrationOrder: order of the quadrature scheme. If reducedIntegrationOrder<0 the integration order is selected independently.\n"
":param useElementsOnFace:  whether or not to use elements on face\n"
":type useElementsOnFace: ``int``"
":param periodic0:  whether or not boundary conditions are periodic\n"
":param periodic1:\n"
":param useFullElementOrder: Whether or not to use Rec9 elements\n"":type useFullElementOrder: ``bool``\n"
":param optimize: Enable optimisation of node labels\n:type optimize: ``bool``"
);

  def("Merge",finley::meshMerge,args("meshList")
,"Merges a list of meshes into one mesh.\n\n:rtype: `Domain`"
  );

  def("GlueFaces",finley::glueFaces,
      (arg("meshList"),arg("safetyFactor")=0.2,
      arg("tolerance")=1.e-8,
      arg("optimize")=true)
,"Detects matching faces in the mesh, removes them from the mesh and joins the elements touched by the face elements."
	);

  def("JoinFaces",finley::joinFaces,
      (arg("meshList"), arg("safetyFactor")=0.2,
      arg("tolerance")=1.e-8,
      arg("optimize")=true)
,"Detects matching faces in the mesh and replaces them by joint elements."
	);


  class_<finley::MeshAdapter, bases<escript::AbstractContinuousDomain> >
      ("MeshAdapter","A concrete class representing a domain. For more details, please consult the c++ documentation.",init<optional <finley::Mesh*> >())
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
":param d:\n:type d: `Data`\n"
":param useHRZ:\n:type useHRZ: bool\n")
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
":param tp:\n:type tp: `TransportProblem`\n"
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
"creates a stiffness matrix and initializes it with zeros\n\n"
":param row_blocksize:\n:type row_blocksize: ``int``\n"
":param row_functionspace:\n:type row_functionspace: `FunctionSpace`\n"
":param column_blocksize:\n:type column_blocksize: ``int``\n"
":param column_functionspace:\n:type column_functionspace: `FunctionSpace`\n"
":param type:\n:type type: ``int``\n"
)
      .def("newTransportProblem",&finley::MeshAdapter::newTransportProblem,
args("theta", "blocksize", "functionspace", "type"),
"creates a TransportProblem\n\n"
":param theta:\n:type theta: ``float``\n"
":param blocksize:\n:type blocksize: ``int``\n"
":param functionspace:\n:type functionspace: `FunctionSpace`\n"
":param type:\n:type type: ``int``\n"
)
      .def("getSystemMatrixTypeId",&finley::MeshAdapter::getSystemMatrixTypeId,
args("options"),
":return: the identifier of the matrix type to be used for the global stiffness matrix when particular solver options are used.\n"
":rtype: ``int``\n"
":param options:\n:type options: `SolverBuddy`\n"
)
      .def("getTransportTypeId",&finley::MeshAdapter::getTransportTypeId,
args("solver", "preconditioner", "package", "symmetry"),
":return: the identifier of the transport problem type to be used when a particular solver, preconditioner, package and symmetric matrix is used.\n"
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
      .def("onMasterProcessor",&finley::MeshAdapter::onMasterProcessor,":return: True if this code is executing on the master process\n:rtype: `bool`")
//       .def("addDiracPoints", &finley::MeshAdapter::addDiracPoints,(arg("points"), arg("tags")=finley::EmptyPythonList),
//       "adds points to support more Dirac delta function.\n\n"
//        	":param points: list of points where Dirac delta function are to be defined. The location will be adjusted to\n"
// 	"               to match the nearest node of the finite element mesh.\n:type points: ``list`` of ``list``s of ``floats``\n"
// 	":param tags: list of tags to be assigned to the points.\n:type tags: ``list`` of ``int``s or ``strings``\n"     
//       )
//       .def("addDiracPoint", &finley::MeshAdapter::addDiracPoint,(arg("point"), arg("tag")=-1),
//         "adds a point to support more Dirac delta function.\n\n"
//        	":param point: point where Dirac delta function are to be defined. The location will be adjusted to\n"
// 	"               to match the nearest node of the finite element mesh.\n:type point: ``list`` of ``floats``\n"
// 	":param tags: tag to be assigned to the point.\n:type tag: ``int``s or ``strings``\n"     
//       ) 
//       .def("addDiracPoint", &finley::MeshAdapter::addDiracPointWithTagName,(arg("point"), arg("tag")),
//         "adds a point to support more Dirac delta function.\n\n"
//        	":param point: point where Dirac delta function are to be defined. The location will be adjusted to\n"
// 	"               to match the nearest node of the finite element mesh.\n:type point: ``list`` of ``floats``\n"
// 	":param tags: tag to be assigned to the point.\n:type tag: ``int``s or ``strings``\n"     
//       ) 
 ;
}

