
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

#include <ripley/Brick.h>
#include <ripley/Rectangle.h>
#include <esysUtils/esysExceptionTranslator.h>

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/detail/defaults_gen.hpp>
#include <boost/version.hpp>

using namespace boost::python;

namespace ripley {

// These wrappers are required to make the shared pointers work through the
// Python wrapper

escript::Domain_ptr _brick(int n0, int n1, int n2, double l0, double l1, double l2, int d0, int d1, int d2)
{
    return escript::Domain_ptr(new Brick(n0, n1, n2, l0, l1, l2, d0, d1, d2));
}

escript::Domain_ptr _rectangle(int n0, int n1, double l0, double l1, int d0, int d1)
{
    return escript::Domain_ptr(new Rectangle(n0, n1, l0, l1, d0, d1));
}

}

/**
    \page ripley Ripley
    Ripley is the python module name that contains the interfaces
    to the C++ wrapper to ripley.
*/

BOOST_PYTHON_MODULE(ripleycpp)
{
// This feature was added in boost v1.34
#if ((BOOST_VERSION/100)%1000 > 34) || (BOOST_VERSION/100000 >1)
    // params are: bool show_user_defined, bool show_py_signatures, bool show_cpp_signatures
    docstring_options docopt(true, true, false);
#endif

    register_exception_translator<ripley::RipleyException>(&(esysUtils::esysExceptionTranslator));

    def("Brick", ripley::_brick, (arg("n0"),arg("n1"),arg("n2"),arg("l0")=1.0,arg("l1")=1.0,arg("l2")=1.0,arg("d0")=1,arg("d1")=1,arg("d2")=1),
"Creates a hexagonal mesh with n0 x n1 x n2 elements over the brick [0,l0] x [0,l1] x [0,l2].\n\n"
":param n0: number of elements in direction 0\n:type n0: ``int``\n"
":param n1: number of elements in direction 1\n:type n1: ``int``\n"
":param n2: number of elements in direction 2\n:type n2: ``int``\n"
":param l0: length of side 0\n:type l0: ``float``\n"
":param l1: length of side 1\n:type l1: ``float``\n"
":param l2: length of side 2\n:type l2: ``float``\n"
":param d0: number of subdivisions in direction 0\n:type d0: ``int``\n"
":param d1: number of subdivisions in direction 1\n:type d1: ``int``\n"
":param d2: number of subdivisions in direction 2\n:type d2: ``int``");

    def("Rectangle", ripley::_rectangle, (arg("n0"),arg("n1"),arg("l0")=1.0,arg("l1")=1.0,arg("d0")=1,arg("d1")=1),
"Creates a rectangular mesh with n0 x n1 elements over the rectangle [0,l0] x [0,l1].\n\n"
":param n0: number of elements in direction 0\n:type n0: ``int``\n"
":param n1: number of elements in direction 1\n:type n1: ``int``\n"
":param l0: length of side 0\n:type l0: ``float``\n"
":param l1: length of side 1\n:type l1: ``float``\n"
":param d0: number of subdivisions in direction 0\n:type d0: ``int``\n"
":param d1: number of subdivisions in direction 1\n:type d1: ``int``");

    def("LoadMesh", ripley::RipleyDomain::loadMesh, (arg("filename")),
           "Loads a ripley domain from a dump file" ":rtype: `Domain`");

    def("ReadMesh", ripley::RipleyDomain::readMesh, (arg("filename")),
            "Reads a ripley domain from a file created by write().\n\n"
            ":rtype: `RipleyDomain`\n:param filename:\n:type filename: ``string``\n");

    class_<ripley::RipleyDomain, bases<escript::AbstractContinuousDomain> >
        ("RipleyDomain", "", no_init)
        .def("write", &ripley::RipleyDomain::write, args("filename"),
                "Writes the current mesh to a file with the given name. It can subsequently be recovered using ReadMesh().")
        .def("print_mesh_info", &ripley::RipleyDomain::Print_Mesh_Info, (arg("full")=false),
                "Prints out a summary about the mesh.\n"
                ":param full: whether to output additional data\n:type full: ``bool``")
        .def("dump", &ripley::RipleyDomain::dump, args("filename"),
                "Dumps the mesh to a file with the given name.")
        .def("getDescription", &ripley::RipleyDomain::getDescription,
":return: a description for this domain\n:rtype: ``string``")
        .def("getDim", &ripley::RipleyDomain::getDim, ":rtype: ``int``")
        .def("getDataShape", &ripley::RipleyDomain::getDataShape, args("functionSpaceCode"),
":return: a pair (dps, ns) where dps=the number of data points per sample, and ns=the number of samples\n:rtype: ``tuple``")
        .def("getNumDataPointsGlobal", &ripley::RipleyDomain::getNumDataPointsGlobal,
":return: the number of data points summed across all MPI processes\n"
":rtype: ``int``")
        .def("addPDEToSystem",&ripley::RipleyDomain::addPDEToSystem,
args("mat", "rhs", "A", "B", "C", "D", "X", "Y", "d", "y", "d_contact", "y_contact"),
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
//      .def("setX",&ripley::RipleyDomain::setNewX,
//args("arg"), "assigns new location to the domain\n\n:param arg:\n:type arg: `Data`")
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
":return: True if ``name`` corresponds to a tag, otherwise False\n:rtype: ``bool``")
        .def("showTagNames",&ripley::RipleyDomain::showTagNames,":return: A space separated list of tag names\n:rtype: ``string``")
        .def("getMPISize",&ripley::RipleyDomain::getMPISize,":return: the number of processes used for this `Domain`\n:rtype: ``int``")
        .def("getMPIRank",&ripley::RipleyDomain::getMPIRank,":return: the rank of this process\n:rtype: ``int``")
        .def("MPIBarrier",&ripley::RipleyDomain::MPIBarrier,"Wait until all processes have reached this point")
        .def("onMasterProcessor",&ripley::RipleyDomain::onMasterProcessor,":return: True if this code is executing on the master process\n:rtype: `bool`");

    class_<ripley::Brick, bases<ripley::RipleyDomain> >("RipleyBrick", "", no_init);

    class_<ripley::Rectangle, bases<ripley::RipleyDomain> >("RipleyRectangle", "", no_init);

}

