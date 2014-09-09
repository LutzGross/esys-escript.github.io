
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

#include <ripley/AbstractAssembler.h>
#include <ripley/Brick.h>
#include <ripley/Rectangle.h>
#include <esysUtils/esysExceptionTranslator.h>

#include <boost/python.hpp> 
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/detail/defaults_gen.hpp>
#include <boost/version.hpp>

#include "escript/SubWorld.h"

using namespace boost::python;

namespace ripley {

template<typename T>
std::vector<T> extractPyArray(const object& obj, const std::string& name,
                              int expectedLength=0)
{
    std::vector<T> result;
    if (extract<tuple>(obj).check() || extract<list>(obj).check()) {
        if (expectedLength==0 || len(obj)==expectedLength) {
            for (int i=0; i<len(obj); i++) {
                result.push_back(extract<T>(obj[i]));
            }
        } else {
            std::stringstream ssError;
            ssError << "argument '" << name << "' has wrong length";
            std::string error(ssError.str());
            throw RipleyException(error.c_str());
        }
    } else {
        std::stringstream ssError;
        ssError << "argument '" << name << "' must be a tuple or list";
        std::string error(ssError.str());
        throw RipleyException(error.c_str());
    }

    return result;
}

escript::Data readBinaryGrid(std::string filename, escript::FunctionSpace fs,
        const object& pyShape, double fill, int byteOrder, int dataType,
        const object& pyFirst, const object& pyNum, const object& pyMultiplier,
        const object& pyReverse)
{
    int dim=fs.getDim();
    ReaderParameters params;

    params.first = extractPyArray<int>(pyFirst, "first", dim);
    params.numValues = extractPyArray<int>(pyNum, "numValues", dim);
    params.multiplier = extractPyArray<int>(pyMultiplier, "multiplier", dim);
    params.reverse = extractPyArray<int>(pyReverse, "reverse", dim);
    params.byteOrder = byteOrder;
    params.dataType = dataType;
    std::vector<int> shape(extractPyArray<int>(pyShape, "shape"));

    const RipleyDomain* dom=dynamic_cast<const RipleyDomain*>(fs.getDomain().get());
    if (!dom)
        throw RipleyException("Function space must be on a ripley domain");

    escript::Data res(fill, shape, fs, true);
    dom->readBinaryGrid(res, filename, params);
    return res;
}

#ifdef USE_BOOSTIO
escript::Data readBinaryGridFromZipped(std::string filename, escript::FunctionSpace fs,
        const object& pyShape, double fill, int byteOrder, int dataType,
        const object& pyFirst, const object& pyNum, const object& pyMultiplier,
        const object& pyReverse)
{
    int dim=fs.getDim();
    ReaderParameters params;

    params.first = extractPyArray<int>(pyFirst, "first", dim);
    params.numValues = extractPyArray<int>(pyNum, "numValues", dim);
    params.multiplier = extractPyArray<int>(pyMultiplier, "multiplier", dim);
    params.reverse = extractPyArray<int>(pyReverse, "reverse", dim);
    params.byteOrder = byteOrder;
    params.dataType = dataType;
    std::vector<int> shape(extractPyArray<int>(pyShape, "shape"));

    const RipleyDomain* dom=dynamic_cast<const RipleyDomain*>(fs.getDomain().get());
    if (!dom)
        throw RipleyException("Function space must be on a ripley domain");

    escript::Data res(fill, shape, fs, true);
    dom->readBinaryGridFromZipped(res, filename, params);
    return res;
}
#endif

escript::Data readNcGrid(std::string filename, std::string varname,
        escript::FunctionSpace fs, const object& pyShape, double fill,
        const object& pyFirst, const object& pyNum, const object& pyMultiplier,
        const object& pyReverse)
{
    int dim=fs.getDim();
    ReaderParameters params;

    params.first = extractPyArray<int>(pyFirst, "first", dim);
    params.numValues = extractPyArray<int>(pyNum, "numValues", dim);
    params.multiplier = extractPyArray<int>(pyMultiplier, "multiplier", dim);
    params.reverse = extractPyArray<int>(pyReverse, "reverse", dim);
    std::vector<int> shape(extractPyArray<int>(pyShape, "shape"));

    const RipleyDomain* dom=dynamic_cast<const RipleyDomain*>(fs.getDomain().get());
    if (!dom)
        throw RipleyException("Function space must be on a ripley domain");

    escript::Data res(fill, shape, fs, true);
    dom->readNcGrid(res, filename, varname, params);
    return res;
}

// These wrappers are required to make the shared pointers work through the
// Python wrapper

// The double for n? is just to keep python happy when people need to deal with
// truediv
escript::Domain_ptr _brick(double _n0, double _n1, double _n2, const object& l0,
                 const object& l1, const object& l2, int d0, int d1, int d2,
                 const object& objpoints, const object& objtags, escript::SubWorld_ptr world)
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
    boost::python::list pypoints=extract<boost::python::list>(objpoints);
    boost::python::list pytags=extract<boost::python::list>(objtags);
    int numpts=extract<int>(pypoints.attr("__len__")());
    int numtags=extract<int>(pytags.attr("__len__")());
    std::vector<double> points;
    std::vector<int> tags;
    tags.resize(numtags, -1);
    for (int i=0;i<numpts;++i) {
        tuple temp = extract<tuple>(pypoints[i]);
        int l=extract<int>(temp.attr("__len__")());
        if (l != 3)
            throw RipleyException("Number of coordinates for each dirac point must match dimensions.");
        for (int k=0;k<l;++k) {
            points.push_back(extract<double>(temp[k]));
        }
    }
    std::map<std::string, int> tagstonames;
    int curmax=40;
    // but which order to assign tags to names?????
    for (int i=0;i<numtags;++i) {
        extract<int> ex_int(pytags[i]);
        extract<std::string> ex_str(pytags[i]);
        if (ex_int.check()) {
            tags[i]=ex_int();
            if (tags[i]>= curmax) {
                curmax=tags[i]+1;
            }
        } else if (ex_str.check()) {
            std::string s=ex_str();
            std::map<std::string, int>::iterator it=tagstonames.find(s);
            if (it!=tagstonames.end()) {
                // we have the tag already so look it up
                tags[i]=it->second;
            } else {
                tagstonames[s]=curmax;
                tags[i]=curmax;
                curmax++;
            }
        } else {
            throw RipleyException("Error - Unable to extract tag value.");
        }
    }
    if (numtags != numpts)
        throw RipleyException("Number of tags does not match number of points.");
    return escript::Domain_ptr(new Brick(n0,n1,n2, x0,y0,z0, x1,y1,z1, d0,d1,d2,
                                            points, tags, tagstonames, world));
}

//const int _q[]={0x61686969,0x746c4144,0x79616e43};
const int _q[]={0x62207363, 0x6574735F, 0x2020214e};
escript::Domain_ptr _rectangle(double _n0, double _n1, const object& l0,
                               const object& l1, int d0, int d1, 
                               const object& objpoints, const object& objtags,
			      escript::SubWorld_ptr world
			      )
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
    boost::python::list pypoints=extract<boost::python::list>(objpoints);
    boost::python::list pytags=extract<boost::python::list>(objtags);
    int numpts=extract<int>(pypoints.attr("__len__")());
    int numtags=extract<int>(pytags.attr("__len__")());
    std::vector<double> points;
    std::vector<int> tags;
    tags.resize(numtags, -1);
    for (int i=0;i<numpts;++i) {
        tuple temp = extract<tuple>(pypoints[i]);
        int l=extract<int>(temp.attr("__len__")());
        if (l != 2)
            throw RipleyException("Number of coordinates for each dirac point must match dimensions.");
        for (int k=0;k<l;++k) {
            points.push_back(extract<double>(temp[k]));
        }
    }
    std::map<std::string, int> tagstonames;
    int curmax=40;
    // but which order to assign tags to names?????
    for (int i=0;i<numtags;++i) {
        extract<int> ex_int(pytags[i]);
        extract<std::string> ex_str(pytags[i]);
        if (ex_int.check()) {
            tags[i]=ex_int();
            if (tags[i]>= curmax) {
                curmax=tags[i]+1;
            }
        } else if (ex_str.check()) {
            std::string s=ex_str();
            std::map<std::string, int>::iterator it=tagstonames.find(s);
            if (it!=tagstonames.end()) {
                // we have the tag already so look it up
                tags[i]=it->second;
            } else {
                tagstonames[s]=curmax;
                tags[i]=curmax;
                curmax++;
            }
        } else {
            throw RipleyException("Error - Unable to extract tag value.");
        }
    }
    if (numtags != numpts)
        throw RipleyException("Number of tags does not match number of points.");
    return escript::Domain_ptr(new Rectangle(n0,n1, x0,y0, x1,y1, d0,d1,
                                             points, tags, tagstonames, world));
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

    register_exception_translator<ripley::RipleyException>(&(esysUtils::RuntimeErrorTranslator));

    scope().attr("__doc__") = "To use this module, please import esys.ripley";
    scope().attr("BYTEORDER_NATIVE") = (int)ripley::BYTEORDER_NATIVE;
    scope().attr("BYTEORDER_LITTLE_ENDIAN") = (int)ripley::BYTEORDER_LITTLE_ENDIAN;
    scope().attr("BYTEORDER_BIG_ENDIAN") = (int)ripley::BYTEORDER_BIG_ENDIAN;
    scope().attr("DATATYPE_INT32") = (int)ripley::DATATYPE_INT32;
    scope().attr("DATATYPE_FLOAT32") = (int)ripley::DATATYPE_FLOAT32;
    scope().attr("DATATYPE_FLOAT64") = (int)ripley::DATATYPE_FLOAT64;

    def("Brick", ripley::_brick, (arg("n0"),arg("n1"),arg("n2"),arg("l0")=1.0,arg("l1")=1.0,arg("l2")=1.0,
        arg("d0")=-1,arg("d1")=-1,arg("d2")=-1,arg("diracPoints")=list(),arg("diracTags")=list(), arg("escriptworld")=escript::SubWorld_ptr()),
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

    def("Rectangle", ripley::_rectangle, (arg("n0"),arg("n1"),arg("l0")=1.0,arg("l1")=1.0,arg("d0")=-1,arg("d1")=-1,arg("diracPoints")=list(),arg("diracTags")=list(), arg("escriptworld")=escript::SubWorld_ptr()),
"Creates a rectangular mesh with n0 x n1 elements over the rectangle [0,l0] x [0,l1].\n\n"
":param n0: number of elements in direction 0\n:type n0: ``int``\n"
":param n1: number of elements in direction 1\n:type n1: ``int``\n"
":param l0: length of side 0 or coordinate range of side 0\n:type l0: ``float`` or ``tuple``\n"
":param l1: length of side 1 or coordinate range of side 1\n:type l1: ``float`` or ``tuple``\n"
":param d0: number of subdivisions in direction 0\n:type d0: ``int``\n"
":param d1: number of subdivisions in direction 1\n:type d1: ``int``");
    def("_theculprit_", ripley::_who);

    def("readBinaryGrid", &ripley::readBinaryGrid, (arg("filename"),
                arg("functionspace"), arg("shape"), arg("fill")=0.,
                arg("byteOrder"), arg("dataType"), arg("first"),
                arg("numValues"), arg("multiplier"), arg("reverse")),
"Reads a binary Grid");
#ifdef USE_BOOSTIO
    def("_readBinaryGridFromZipped", &ripley::readBinaryGridFromZipped, (arg("filename"),
                arg("functionspace"), arg("shape"), arg("fill")=0.,
                arg("byteOrder"), arg("dataType"), arg("first"),
                arg("numValues"), arg("multiplier"), arg("reverse")),
"Reads a binary Grid");
#endif
    def("_readNcGrid", &ripley::readNcGrid, (arg("filename"), arg("varname"),
                arg("functionspace"), arg("shape"), arg("fill"), arg("first"),
                arg("numValues"), arg("multiplier"), arg("reverse")),
"Reads a grid from a netCDF file");

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
        .def("addToSystem",&ripley::RipleyDomain::addToSystemFromPython,
            args("mat", "rhs", "data"),
            "adds a PDE to the system, results depend on domain\n\n"
            ":param mat:\n:type mat: `OperatorAdapter`\n"
            ":param rhs:\n:type rhs: `Data`\n"
            ":param data:\ntype data: `list`")
        .def("addToRHS",&ripley::RipleyDomain::addToRHSFromPython,
            args("rhs", "data"),
            "adds a PDE onto the stiffness matrix mat and a rhs, "
            "results depends on domain\n\n"
            ":param rhs:\n:type rhs: `Data`\n"
            ":param data:\ntype data: `list`")
        .def("createAssembler", &ripley::RipleyDomain::createAssemblerFromPython,
            args("typename", "options"),
            "request from the domain an assembler of the specified type, if "
            "supported, using the supplied options (if provided)"
            ":param typename:\n:type typename: `string`\n"
            ":param options:\n:type options: `list`\n")
        .def("addPDEToTransportProblem",&ripley::RipleyDomain::addPDEToTransportProblemFromPython,
            args("tp", "source", "data"),
            ":param tp:\n:type tp: `TransportProblemAdapter`\n"
            ":param source:\n:type source: `Data`\n"
            ":param data:\ntype data: `list`")
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
args("options"),
":return: the identifier of the matrix type to be used for the global stiffness matrix when particular solver options are used.\n"
":rtype: ``int``\n"
":param options:\n:type options: `SolverBuddy`\n"
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
    /* These two class exports are necessary to ensure that the extra methods added by ripley make it to python.
     * This change became necessary when the Brick and Rectangle constructors turned into factories instead of classes */
    class_<ripley::Brick, bases<ripley::RipleyDomain> >("RipleyBrick", "", no_init);
    class_<ripley::Rectangle, bases<ripley::RipleyDomain> >("RipleyRectangle", "", no_init);
    class_<ripley::AbstractAssembler, ripley::Assembler_ptr, boost::noncopyable >
        ("AbstractAssembler", "", no_init);
}

