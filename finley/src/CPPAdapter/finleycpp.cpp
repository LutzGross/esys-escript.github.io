
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
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

  //
  // NOTE: The return_value_policy is necessary for functions that
  // return pointers.
  //
  register_exception_translator<finley::FinleyAdapterException>(&(esysUtils::esysExceptionTranslator));

  def("LoadMesh",finley::loadMesh,
      (arg("fileName")="file.nc"),
      return_value_policy<manage_new_object>());
  def("ReadMesh",finley::readMesh,
      (arg("fileName")="file.fly",arg("integrationOrder")=-1,  arg("reducedIntegrationOrder")=-1,  arg("optimize")=true),
      return_value_policy<manage_new_object>());
  def("ReadMeshMPI",finley::readMeshMPI,
      (arg("fileName")="file.fly",arg("integrationOrder")=-1,  arg("reducedIntegrationOrder")=-1,  arg("optimize")=true),
      return_value_policy<manage_new_object>());

  def("ReadGmsh",finley::readGmsh,
      (arg("fileName")="file.msh",arg("numDim"), arg("integrationOrder")=-1, arg("reducedIntegrationOrder")=-1, arg("optimize")=true),
      return_value_policy<manage_new_object>());

  def ("Brick",finley::brick,
      (arg("n0")=1,arg("n1")=1,arg("n2")=1,
      arg("order")=1,
      arg("l0")=1.0,arg("l1")=1.0,arg("l2")=1.0,
      arg("periodic0")=false,arg("periodic1")=false,arg("periodic2")=false,
      arg("integrationOrder")=-1,  arg("reducedIntegrationOrder")=-1,
      arg("useElementsOnFace")=false,
      arg("useFullElementOrder")=false,
      arg("optimize")=false),
      return_value_policy<manage_new_object>());

  def ("Rectangle",finley::rectangle,
      (arg("n0")=1,arg("n1")=1,arg("order")=1,
      arg("l0")=1.0,arg("l1")=1.0,
      arg("periodic0")=false,arg("periodic1")=false,
      arg("integrationOrder")=-1,  arg("reducedIntegrationOrder")=-1,
      arg("useElementsOnFace")=false,
      arg("useFullElementOrder")=false,
      arg("optimize")=false),
      return_value_policy<manage_new_object>());

  def("Merge",finley::meshMerge,
      return_value_policy<manage_new_object>());

  def("GlueFaces",finley::glueFaces,
      (arg("safetyFactor")=0.2,
      arg("tolerance")=1.e-8,
      arg("optimize")=true),
      return_value_policy<manage_new_object>());

  def("JoinFaces",finley::joinFaces,
      (arg("safetyFactor")=0.2,
      arg("tolerance")=1.e-8,
      arg("optimize")=true),
      return_value_policy<manage_new_object>());



  class_<finley::MeshAdapter, bases<escript::AbstractContinuousDomain> >
      ("MeshAdapter",init<optional <Finley_Mesh*> >())
      .def(init<const finley::MeshAdapter&>())
      .def("write",&finley::MeshAdapter::write)
      .def("print_mesh_info",&finley::MeshAdapter::Print_Mesh_Info,(arg("full")=false))
      .def("dump",&finley::MeshAdapter::dump)
      .def("getDescription",&finley::MeshAdapter::getDescription)
      .def("getDim",&finley::MeshAdapter::getDim)
      .def("getDataShape",&finley::MeshAdapter::getDataShape)
      .def("getNumDataPointsGlobal",&finley::MeshAdapter::getNumDataPointsGlobal)
      .def("addPDEToSystem",&finley::MeshAdapter::addPDEToSystem)
      .def("addPDEToLumpedSystem",&finley::MeshAdapter::addPDEToLumpedSystem)
      .def("addPDEToRHS",&finley::MeshAdapter::addPDEToRHS)
      .def("addPDEToTransportProblem",&finley::MeshAdapter::addPDEToTransportProblem)
      .def("newOperator",&finley::MeshAdapter::newSystemMatrix)
      .def("newTransportProblem",&finley::MeshAdapter::newTransportProblem)
      .def("getSystemMatrixTypeId",&finley::MeshAdapter::getSystemMatrixTypeId)
      .def("setX",&finley::MeshAdapter::setNewX)
      .def("getX",&finley::MeshAdapter::getX)
      .def("getNormal",&finley::MeshAdapter::getNormal)
      .def("getSize",&finley::MeshAdapter::getSize)
      .def("saveDX",&finley::MeshAdapter::saveDX)
      .def("saveVTK",&finley::MeshAdapter::saveVTK)
      .def("setTagMap",&finley::MeshAdapter::setTagMap)
      .def("getTag",&finley::MeshAdapter::getTag)
      .def("isValidTagName",&finley::MeshAdapter::isValidTagName)
      .def("showTagNames",&finley::MeshAdapter::showTagNames)
      .def("getMPISize",&finley::MeshAdapter::getMPISize)
      .def("getMPIRank",&finley::MeshAdapter::getMPIRank)
      .def("getMPIComm",&finley::MeshAdapter::getMPIComm);

  class_<finley::SystemMatrixAdapter, bases<escript::AbstractSystemMatrix> >
      ("OperatorAdapter",no_init)
      .def("print_matrix_info",&finley::SystemMatrixAdapter::Print_Matrix_Info,(arg("full")=false))
      .def("nullifyRowsAndCols",&finley::SystemMatrixAdapter::nullifyRowsAndCols)
      .def("resetValues",&finley::SystemMatrixAdapter::resetValues)
      .def("saveMM",&finley::SystemMatrixAdapter::saveMM)
      .def("saveHB",&finley::SystemMatrixAdapter::saveHB);

  class_<finley::TransportProblemAdapter, bases<escript::AbstractTransportProblem> >
      ("TransportProblemAdapter",no_init)
      .def("getSafeTimeStepSize",&finley::TransportProblemAdapter::getSafeTimeStepSize)
      .def("resetTransport",&finley::TransportProblemAdapter::resetTransport);
}
