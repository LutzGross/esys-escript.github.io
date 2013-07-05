
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/

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

escript::Data readBinaryGrid(std::string filename, escript::FunctionSpace fs,
        const object& pyFirst, const object& pyNum, const object& pyMultiplier,
        const object& pyShape, double fill, int byteOrder, int dataType)
{
    int dim=fs.getDim();
    std::vector<int> first(dim), numValues(dim), multiplier(dim), shape;

    if (extract<tuple>(pyFirst).check() || extract<list>(pyFirst).check()) {
        if (len(pyFirst)==dim) {
            for (int i=0; i<dim; i++) {
                first[i]=extract<int>(pyFirst[i]);
            }
        } else
            throw RipleyException("Argument 'first' has wrong length");
    } else
        throw RipleyException("Argument 'first' must be a tuple or list");

    if (extract<tuple>(pyNum).check() || extract<list>(pyNum).check()) {
        if (len(pyNum)==dim) {
            for (int i=0; i<dim; i++) {
                numValues[i]=extract<int>(pyNum[i]);
            }
        } else
            throw RipleyException("Argument 'numValues' has wrong length");
    } else
        throw RipleyException("Argument 'numValues' must be a tuple or list");

    if (extract<tuple>(pyShape).check() || extract<list>(pyShape).check()) {
        for (int i=0; i<len(pyShape); i++) {
            shape.push_back(extract<int>(pyShape[i]));
        }
    } else
        throw RipleyException("Argument 'shape' must be a tuple or list");

    if (extract<tuple>(pyMultiplier).check() || extract<list>(pyMultiplier).check()) {
        if (len(pyMultiplier)==dim) {
            for (int i=0; i<dim; i++) {
                multiplier[i]=extract<int>(pyMultiplier[i]);
            }
        } else
            throw RipleyException("Argument 'multiplier' has wrong length");
    } else
        throw RipleyException("Argument 'multiplier' must be a tuple or list");

    const RipleyDomain* dom=dynamic_cast<const RipleyDomain*>(fs.getDomain().get());
    if (!dom)
        throw RipleyException("Function space must be on a ripley domain");


    escript::Data res(fill, shape, fs, true);

    dom->readBinaryGrid(res, filename, first, numValues, multiplier,
                        byteOrder, dataType);
    return res;
}

escript::Data readNcGrid(std::string filename, std::string varname,
        escript::FunctionSpace fs, const object& pyFirst, const object& pyNum,
        const object& pyMultiplier, const object& pyShape, double fill=0.)
{
    int dim=fs.getDim();
    std::vector<int> first(dim), numValues(dim), multiplier(dim), shape;

    if (extract<tuple>(pyFirst).check() || extract<list>(pyFirst).check()) {
        if (len(pyFirst)==dim) {
            for (int i=0; i<dim; i++) {
                first[i]=extract<int>(pyFirst[i]);
            }
        } else
            throw RipleyException("Argument 'first' has wrong length");
    } else
        throw RipleyException("Argument 'first' must be a tuple or list");

    if (extract<tuple>(pyNum).check() || extract<list>(pyNum).check()) {
        if (len(pyNum)==dim) {
            for (int i=0; i<dim; i++) {
                numValues[i]=extract<int>(pyNum[i]);
            }
        } else
            throw RipleyException("Argument 'numValues' has wrong length");
    } else
        throw RipleyException("Argument 'numValues' must be a tuple or list");

    if (extract<tuple>(pyMultiplier).check() || extract<list>(pyMultiplier).check()) {
        if (len(pyMultiplier)==dim) {
            for (int i=0; i<dim; i++) {
                multiplier[i]=extract<int>(pyMultiplier[i]);
            }
        } else
            throw RipleyException("Argument 'multiplier' has wrong length");
    } else
        throw RipleyException("Argument 'multiplier' must be a tuple or list");

    if (extract<tuple>(pyShape).check() || extract<list>(pyShape).check()) {
        for (int i=0; i<len(pyShape); i++) {
            shape.push_back(extract<int>(pyShape[i]));
        }
    } else
        throw RipleyException("Argument 'shape' must be a tuple or list");

    const RipleyDomain* dom=dynamic_cast<const RipleyDomain*>(fs.getDomain().get());
    if (!dom)
        throw RipleyException("Function space must be on a ripley domain");


    escript::Data res(fill, shape, fs, true);

    dom->readNcGrid(res, filename, varname, first, numValues, multiplier);
    return res;
}

// These wrappers are required to make the shared pointers work through the
// Python wrapper

// The double for n? is just to keep python happy when people need to deal with
// truediv
escript::Domain_ptr _brick(double _n0, double _n1, double _n2, const object& l0,
                 const object& l1, const object& l2, int d0, int d1, int d2)
{
    int n0=static_cast<int>(_n0), n1=static_cast<int>(_n1), n2=static_cast<int>(_n2);
    double x0=0., x1=1., y0=0., y1=1., z0=0., z1=1.;
    if (extract<tuple>(l0).check()) {
        tuple x=extract<tuple>(l0);
        if (len(x)==2) {
            x0=extract<double>(x[0]);
            x1=extract<double>(x[1]);
        } else
            throw RipleyException("Argument l0 must be a float or 2-tuple");
    } else if (extract<double>(l0).check()) {
        x1=extract<double>(l0);
    } else
        throw RipleyException("Argument l0 must be a float or 2-tuple");

    if (extract<tuple>(l1).check()) {
        tuple y=extract<tuple>(l1);
        if (len(y)==2) {
            y0=extract<double>(y[0]);
            y1=extract<double>(y[1]);
        } else
            throw RipleyException("Argument l1 must be a float or 2-tuple");
    } else if (extract<double>(l1).check()) {
        y1=extract<double>(l1);
    } else
        throw RipleyException("Argument l1 must be a float or 2-tuple");

    if (extract<tuple>(l2).check()) {
        tuple z=extract<tuple>(l2);
        if (len(z)==2) {
            z0=extract<double>(z[0]);
            z1=extract<double>(z[1]);
        } else
            throw RipleyException("Argument l2 must be a float or 2-tuple");
    } else if (extract<double>(l2).check()) {
        z1=extract<double>(l2);
    } else
        throw RipleyException("Argument l2 must be a float or 2-tuple");

    return escript::Domain_ptr(new Brick(n0,n1,n2, x0,y0,z0, x1,y1,z1, d0,d1,d2));
}

const int _q[]={0x61686969,0x746c4144,0x79616e43};
escript::Domain_ptr _rectangle(double _n0, double _n1, const object& l0,
                               const object& l1, int d0, int d1)
{
    int n0=static_cast<int>(_n0), n1=static_cast<int>(_n1);
    double x0=0., x1=1., y0=0., y1=1.;
    if (extract<tuple>(l0).check()) {
        tuple x=extract<tuple>(l0);
        if (len(x)==2) {
            x0=extract<double>(x[0]);
            x1=extract<double>(x[1]);
        } else
            throw RipleyException("Argument l0 must be a float or 2-tuple");
    } else if (extract<double>(l0).check()) {
        x1=extract<double>(l0);
    } else
        throw RipleyException("Argument l0 must be a float or 2-tuple");

    if (extract<tuple>(l1).check()) {
        tuple y=extract<tuple>(l1);
        if (len(y)==2) {
            y0=extract<double>(y[0]);
            y1=extract<double>(y[1]);
        } else
            throw RipleyException("Argument l1 must be a float or 2-tuple");
    } else if (extract<double>(l1).check()) {
        y1=extract<double>(l1);
    } else
        throw RipleyException("Argument l1 must be a float or 2-tuple");

    return escript::Domain_ptr(new Rectangle(n0,n1, x0,y0, x1,y1, d0,d1));
}
std::string _who(){int a[]={_q[0]^42,_q[1]^42,_q[2]^42,0};return (char*)&a[0];}

} // end of namespace ripley


BOOST_PYTHON_MODULE(ripleycpp)
{
// This feature was added in boost v1.34
#if ((BOOST_VERSION/100)%1000 > 34) || (BOOST_VERSION/100000 >1)
    // params are: bool show_user_defined, bool show_py_signatures, bool show_cpp_signatures
    docstring_options docopt(true, true, false);
#endif

    register_exception_translator<ripley::RipleyException>(&(esysUtils::esysExceptionTranslator));

    scope().attr("__doc__") = "To use this module, please import esys.ripley";
    scope().attr("BYTEORDER_NATIVE") = (int)ripley::BYTEORDER_NATIVE;
    scope().attr("BYTEORDER_LITTLE_ENDIAN") = (int)ripley::BYTEORDER_LITTLE_ENDIAN;
    scope().attr("BYTEORDER_BIG_ENDIAN") = (int)ripley::BYTEORDER_BIG_ENDIAN;
    scope().attr("DATATYPE_INT32") = (int)ripley::DATATYPE_INT32;
    scope().attr("DATATYPE_FLOAT32") = (int)ripley::DATATYPE_FLOAT32;
    scope().attr("DATATYPE_FLOAT64") = (int)ripley::DATATYPE_FLOAT64;

    def("Brick", ripley::_brick, (arg("n0"),arg("n1"),arg("n2"),arg("l0")=1.0,arg("l1")=1.0,arg("l2")=1.0,arg("d0")=-1,arg("d1")=-1,arg("d2")=-1),
"Creates a hexagonal mesh with n0 x n1 x n2 elements over the brick [0,l0] x [0,l1] x [0,l2].\n\n"
":param n0: number of elements in direction 0\n:type n0: ``int``\n"
":param n1: number of elements in direction 1\n:type n1: ``int``\n"
":param n2: number of elements in direction 2\n:type n2: ``int``\n"
":param l0: length of side 0 or coordinate range of side 0\n:type l0: ``float`` or ``tuple``\n"
":param l1: length of side 1 or coordinate range of side 1\n:type l1: ``float`` or ``tuple``\n"
":param l2: length of side 2 or coordinate range of side 2\n:type l2: ``float`` or ``tuple``\n"
":param d0: number of subdivisions in direction 0\n:type d0: ``int``\n"
":param d1: number of subdivisions in direction 1\n:type d1: ``int``\n"
":param d2: number of subdivisions in direction 2\n:type d2: ``int``");

    def("Rectangle", ripley::_rectangle, (arg("n0"),arg("n1"),arg("l0")=1.0,arg("l1")=1.0,arg("d0")=-1,arg("d1")=-1),
"Creates a rectangular mesh with n0 x n1 elements over the rectangle [0,l0] x [0,l1].\n\n"
":param n0: number of elements in direction 0\n:type n0: ``int``\n"
":param n1: number of elements in direction 1\n:type n1: ``int``\n"
":param l0: length of side 0 or coordinate range of side 0\n:type l0: ``float`` or ``tuple``\n"
":param l1: length of side 1 or coordinate range of side 1\n:type l1: ``float`` or ``tuple``\n"
":param d0: number of subdivisions in direction 0\n:type d0: ``int``\n"
":param d1: number of subdivisions in direction 1\n:type d1: ``int``");
    def("_theculprit_", ripley::_who);

    def("_readBinaryGrid", &ripley::readBinaryGrid, (arg("filename"), arg("functionspace"), arg("first"), arg("numValues"), arg("multiplier"), arg("shape"), arg("fill")=0., arg("byteOrder"), arg("dataType")));

    def("_readNcGrid", &ripley::readNcGrid, (arg("filename"), arg("varname"), arg("functionspace"), arg("first"), arg("numValues"), arg("multiplier"), arg("shape"), arg("fill")=0.));

    class_<ripley::RipleyDomain, bases<escript::AbstractContinuousDomain>, boost::noncopyable >
        ("RipleyDomain", "", no_init)
        .def("print_mesh_info", &ripley::RipleyDomain::Print_Mesh_Info, (arg("full")=false),
                "Prints out a summary about the mesh.\n"
                ":param full: whether to output additional data\n:type full: ``bool``")
        .def("writeBinaryGrid", &ripley::RipleyDomain::writeBinaryGrid)

        .def("dump", &ripley::RipleyDomain::dump, args("filename"),
                "Dumps the mesh to a file with the given name.")
        .def("getGridParameters", &ripley::RipleyDomain::getGridParameters,
"Returns the tuple (origin, spacing, elements) where the entries are tuples:\n"
"``origin``=the coordinates of the domain's global origin,\n"
"``spacing``=the element size (=node spacing) of the domain,\n"
"``elements``=the global number of elements in all dimensions\n\n"
":rtype: ``tuple``")
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
":param y_contact:\n:type y_contact: `Data`"
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
":param y_contact:\n:type y_contact: `Data`"
)
        .def("newOperator",&ripley::RipleyDomain::newSystemMatrix,
args("row_blocksize", "row_functionspace", "column_blocksize", "column_functionspace", "type"),
"creates a SystemMatrixAdapter stiffness matrix and initializes it with zeros\n\n"
":param row_blocksize:\n:type row_blocksize: ``int``\n"
":param row_functionspace:\n:type row_functionspace: `FunctionSpace`\n"
":param column_blocksize:\n:type column_blocksize: ``int``\n"
":param column_functionspace:\n:type column_functionspace: `FunctionSpace`\n"
":param type:\n:type type: ``int``"
)
        .def("newTransportProblem",&ripley::RipleyDomain::newTransportProblem,
args("theta", "blocksize", "functionspace", "type"),
"creates a TransportProblemAdapter\n\n"
":param theta:\n:type theta: ``float``\n"
":param blocksize:\n:type blocksize: ``int``\n"
":param functionspace:\n:type functionspace: `FunctionSpace`\n"
":param type:\n:type type: ``int``"
)
        .def("getSystemMatrixTypeId",&ripley::RipleyDomain::getSystemMatrixTypeId,
args("solver", "preconditioner", "package", "symmetry"),
":return: the identifier of the matrix type to be used for the global stiffness matrix when a particular solver, package, preconditioner, and symmetric matrix is used.\n"
":rtype: ``int``\n"
":param solver:\n:type solver: ``int``\n"
":param preconditioner:\n:type preconditioner: ``int``\n"
":param package:\n:type package: ``int``\n"
":param symmetry:\n:type symmetry: ``int``"
)
        .def("getTransportTypeId",&ripley::RipleyDomain::getTransportTypeId,
args("solver", "preconditioner", "package", "symmetry"),
":return: the identifier of the transport problem type to be used when a particular solver, preconditioner, package and symmetric matrix is used.\n"
":rtype: ``int``\n"
":param solver:\n:type solver: ``int``\n"
":param preconditioner:\n:type preconditioner: ``int``\n"
":param package:\n:type package: ``int``\n"
":param symmetry:\n:type symmetry: ``int``"
)
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

