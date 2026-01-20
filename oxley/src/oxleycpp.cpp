/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/

#include <oxley/Brick.h>
#include <oxley/Rectangle.h>
#include <oxley/OtherAlgorithms.h>
#include <oxley/OxleyDomain.h>
#include <oxley/RefinementZone.h>

#include <boost/python.hpp>
#ifdef ESYS_HAVE_BOOST_NUMPY
#include <boost/python/numpy.hpp>
#include <boost/python/numpy/dtype.hpp>
#endif

// using namespace boost::python::numpy;

namespace oxley {

escript::Domain_ptr _rectangle(double _n0, double _n1,
                        const object& l0, const object& l1, int d0, int d1,
                        const object& objpoints, const object& objtags,
                        int periodic0, int periodic1, const object& py_comm)
{
	// // Integration Order
 //    if (order < 2 || order > 10)
 //        throw OxleyException("Order must be in the range 2 to 10");

    // Number of nodes in each direction
    dim_t n0=static_cast<dim_t>(_n0), n1=static_cast<dim_t>(_n1);
    double x0=0., x1=1., y0=0., y1=1.;

    // Length of the domain in each direction
    if (extract<tuple>(l0).check()) {
        tuple x=extract<tuple>(l0);
        if (len(x)==2) {
            x0=extract<double>(x[0]);
            x1=extract<double>(x[1]);
        } else
            throw OxleyException("Argument l0 must be a float or 2-tuple");
    } else if (extract<double>(l0).check()) {
        x1=extract<double>(l0);
    } else
        throw OxleyException("Argument l0 must be a float or 2-tuple");

    if (extract<tuple>(l1).check()) {
        tuple y=extract<tuple>(l1);
        if (len(y)==2) {
            y0=extract<double>(y[0]);
            y1=extract<double>(y[1]);
        } else
            throw OxleyException("Argument l1 must be a float or 2-tuple");
    } else if (extract<double>(l1).check()) {
        y1=extract<double>(l1);
    } else
        throw OxleyException("Argument l1 must be a float or 2-tuple");

    int order = 1;

    // process tags and points
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
            throw OxleyException("Number of coordinates for each dirac point must match dimensions.");
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
            if (tags[i] >= curmax) {
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
            throw OxleyException("Error - Unable to extract tag value.");
        }
    }
    if (numtags != numpts)
        throw OxleyException("Number of tags does not match number of points.");

    // Handle optional MPI communicator
    escript::JMPI jmpi = escript::makeInfoFromPyComm(py_comm);

    return escript::Domain_ptr(new Rectangle(jmpi, order, n0,n1, x0,y0, x1,y1,
                                d0,d1, points, tags, tagstonames, periodic0, periodic1));
}


escript::Domain_ptr _brick(double _n0, double _n1, double _n2,
                        const object& l0, const object& l1, const object& l2,
                        int d0, int d1, int d2,
                        const object& objpoints, const object& objtags,
                        int periodic0, int periodic1, int periodic2, const object& py_comm)
{
    // Integration Order
    int order=2;
    if (order < 2 || order > 10)
        throw OxleyException("Order must be in the range 2 to 10");

    // Number of nodes in each direction
    dim_t n0=static_cast<dim_t>(_n0), n1=static_cast<dim_t>(_n1), n2=static_cast<dim_t>(_n2);;
    double x0=0., x1=1., y0=0., y1=1., z0=0., z1=1.;

    // Length of the domain in each direction
    if (extract<tuple>(l0).check()) {
        tuple x=extract<tuple>(l0);
        if (len(x)==2) {
            x0=extract<double>(x[0]);
            x1=extract<double>(x[1]);
        } else
            throw OxleyException("Argument l0 must be a float or 2-tuple");
    } else if (extract<double>(l0).check()) {
        x1=extract<double>(l0);
    } else
        throw OxleyException("Argument l0 must be a float or 2-tuple");

    if (extract<tuple>(l1).check()) {
        tuple y=extract<tuple>(l1);
        if (len(y)==2) {
            y0=extract<double>(y[0]);
            y1=extract<double>(y[1]);
        } else
            throw OxleyException("Argument l1 must be a float or 2-tuple");
    } else if (extract<double>(l1).check()) {
        y1=extract<double>(l1);
    } else
        throw OxleyException("Argument l1 must be a float or 2-tuple");

    if (extract<tuple>(l2).check()) {
        tuple z=extract<tuple>(l2);
        if (len(z)==2) {
            z0=extract<double>(z[0]);
            z1=extract<double>(z[1]);
        } else
            throw OxleyException("Argument l2 must be a float or 2-tuple");
    } else if (extract<double>(l2).check()) {
        z1=extract<double>(l2);
    } else
        throw OxleyException("Argument l2 must be a float or 2-tuple");

    // process tags and points
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
            throw OxleyException("Number of coordinates for each dirac point must match dimensions.");
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
            throw OxleyException("Error - Unable to extract tag value.");
        }
    }
    if (numtags != numpts)
        throw OxleyException("Number of tags does not match number of points.");

    // Handle optional MPI communicator
    escript::JMPI jmpi = escript::makeInfoFromPyComm(py_comm);

    return escript::Domain_ptr(new Brick(jmpi, order, n0,n1,n2, x0,y0,z0, x1,y1,z1,
                                d0,d1,d2, points, tags, tagstonames, periodic0,periodic1,periodic2));
}

// //tmp
// oxley::RefinementZone_Ptr _refinementZone()
// {
//     return oxley::RefinementZone_Ptr(new RefinementZone());
//     // return oxley::RefinementZone2D_Ptr(new RefinementZone2D());
// }


oxley::RefinementZone2D_Ptr _refinementZone2D()
{
    return oxley::RefinementZone2D_Ptr(new RefinementZone2D());
}

oxley::RefinementZone3D_Ptr _refinementZone3D()
{
    return oxley::RefinementZone3D_Ptr(new RefinementZone3D());
}

BOOST_PYTHON_MODULE(oxleycpp)
{

    def("Rectangle", oxley::_rectangle, (
    arg("n0"),arg("n1"),
    arg("l0")=1.0,arg("l1")=1.0,
    arg("d0")=-1,arg("d1")=-1,
    arg("diracPoints")=list(), arg("diracTags")=list(),
    arg("periodic0")=0,arg("periodic1")=0,
    arg("comm")=object()),
    "Creates a rectangular p4est mesh with n0 x n1 elements over the rectangle [0,l0] x [0,l1].\n\n"
    ":param n0: number of elements in direction 0\n:type n0: ``int``\n"
    ":param n1: number of elements in direction 1\n:type n1: ``int``\n"
    ":param l0: length of side 0 or coordinate range of side 0\n:type l0: ``float`` or ``tuple``\n"
    ":param l1: length of side 1 or coordinate range of side 1\n:type l1: ``float`` or ``tuple``\n"
    ":param d0: number of subdivisions in direction 0\n:type d0: ``int``\n"
    ":param d1: number of subdivisions in direction 1\n:type d1: ``int``\n"
    ":param comm: MPI communicator (optional, from mpi4py)\n:type comm: ``mpi4py.MPI.Comm``");

    def("Brick", oxley::_brick, (
    arg("n0"),arg("n1"),arg("n2"),
    arg("l0")=1.0,arg("l1")=1.0,arg("l2")=1.0,
    arg("d0")=-1,arg("d1")=-1,arg("d2")=-1,
    arg("diracPoints")=list(), arg("diracTags")=list(),
    arg("periodic0")=0,arg("periodic1")=0,arg("periodic2")=0,
    arg("comm")=object()),
    "Creates a brick p4est mesh with n0 x n1 x n2 elements over the rectangle [0,l0] x [0,l1] x [0,l2].\n\n"
    ":param order: order of the elements: ``int``\n"
    ":param n0: number of elements in direction 0\n:type n0: ``int``\n"
    ":param n1: number of elements in direction 1\n:type n1: ``int``\n"
    ":param n2: number of elements in direction 2\n:type n2: ``int``\n"
    ":param l0: length of side 0 or coordinate range of side 0\n:type l0: ``float`` or ``tuple``\n"
    ":param l1: length of side 1 or coordinate range of side 1\n:type l1: ``float`` or ``tuple``\n"
    ":param l2: length of side 2 or coordinate range of side 1\n:type l2: ``float`` or ``tuple``\n"
    ":param d0: number of subdivisions in direction 0\n:type d0: ``int``\n"
    ":param d1: number of subdivisions in direction 1\n:type d1: ``int``\n"
    ":param d2: number of subdivisions in direction 2\n:type d2: ``int``\n"
    ":param comm: MPI communicator (optional, from mpi4py)\n:type comm: ``mpi4py.MPI.Comm``");

    // def("RefinementZone", oxley::_refinementZone, 
    //     "Creates a refinement zone parent class.\n\n"
    //     );

    def("RefinementZone2D", oxley::_refinementZone2D, 
        "Creates a refinement zone of dimension 2.\n\n"
        );

    def("RefinementZone3D", oxley::_refinementZone3D, 
        "Creates a refinement zone of dimension 3.\n\n"
        );

// #ifdef ESYS_HAVE_BOOST_NUMPY
//     // def("addSurface", oxley::_addSurface, (arg("domain"),arg("x"),arg("y"),arg("z")));
//     def("addSurface", oxley::_addCurve, (arg("domain"),arg("x"),arg("y")));
//     def("addSurface", oxley::_addSurface, (arg("domain"),arg("x"),arg("y"),arg("z")));
// #endif

    class_<oxley::OxleyDomain, bases<escript::AbstractContinuousDomain>, boost::noncopyable >
        ("OxleyDomain", "", no_init)
        .def("addToRHS",&oxley::OxleyDomain::addToRHSFromPython,
            args("rhs", "data"),
            "adds a PDE onto the stiffness matrix mat and a rhs, "
            "results depends on domain\n\n"
            ":param rhs:\n:type rhs: `Data`\n"
            ":param data:\n:type data: `list`\n")
        .def("addToSystem",&oxley::OxleyDomain::addToSystemFromPython,
            args("mat", "rhs", "data"),
            "adds a PDE to the system, results depend on domain\n\n"
            ":param mat:\n:type mat: `OperatorAdapter`\n"
            ":param rhs:\n:type rhs: `Data`\n"
            ":param data:\n:type data: `list`\n")
        .def("createAssembler", &oxley::OxleyDomain::createAssemblerFromPython,
            args("typename", "options"),
            "request from the domain an assembler of the specified type, if "
            "supported, using the supplied options (if provided)"
            ":param typename:\n:type typename: `string`\n"
            ":param options:\n:type options: `list`\n")
        .def("dump", &oxley::OxleyDomain::dump, args("filename"),
            "Dumps the mesh to a silo file with the name `filename`."
            ":param filename:\n:type typename: `string`\n")
        .def("getDataShape", &oxley::OxleyDomain::getDataShape, args("functionSpaceCode"),
            ":return: a pair (dps, ns) where dps is the number of data points per sample, and ns is the number of samples\n"
            ":rtype: ``tuple``")        
        .def("getDescription", &oxley::OxleyDomain::getDescription,
                "Prints out a description of the mesh.")
        .def("getDim", &oxley::OxleyDomain::getDim, ":rtype: ``int``")
        .def("getNormal",&oxley::OxleyDomain::getNormal,
            ":return: boundary normals at the quadrature point on the face elements\n"
            ":rtype: `Data`")
        .def("getNumVertices", &oxley::OxleyDomain::getNumVertices,
            "Returns the number of corners in the mesh.\n"
            ":rtype: ``int``")
        .def("getSystemMatrixTypeId",&oxley::OxleyDomain::getSystemMatrixTypeId,
            args("options"),
            ":return: the identifier of the matrix type to be used for the global stiffness matrix when particular solver options are used.\n"
            ":rtype: ``int``\n"
            ":param options:\n:type options: `SolverBuddy`\n")
        .def("setTagMap",&oxley::OxleyDomain::setTagMap,args("name","tag"),
            "Give a tag number a name.\n\n:param name: Name for the tag\n:type name: ``string``\n"
            ":param tag: numeric id\n:type tag: ``int``\n:note: Tag names must be unique within a domain")
        .def("getTransportTypeId",&oxley::OxleyDomain::getTransportTypeId,
            args("solver", "preconditioner", "package", "symmetry"),
            ":return: the identifier of the transport problem type to be used when a particular solver, preconditioner, package and symmetric matrix is used.\n"
            ":rtype: ``int``\n"
            ":param solver:\n:type solver: ``int``\n"
            ":param preconditioner:\n:type preconditioner: ``int``\n"
            ":param package:\n:type package: ``int``\n"
            ":param symmetry:\n:type symmetry: ``int``")
        .def("getX",&oxley::OxleyDomain::getX, ":return: locations in the FEM nodes\n\n"
            ":rtype: `Data`")
        #ifdef ESYS_HAVE_TRILINOS
        .def("makeZ",&oxley::OxleyDomain::makeZ, arg("complex"), "creates the matrix Z")
        .def("makeIZ",&oxley::OxleyDomain::makeIZ, arg("complex"), "creates the matrix IZ")
        .def("finaliseA",&oxley::OxleyDomain::finaliseA, args("mat","isComplex"), "finalisesLHS")
        .def("finaliseRhs",&oxley::OxleyDomain::finaliseRhs, arg("rhs"), "finalisesRHS")
        #endif
        .def("saveFsType",&oxley::OxleyDomain::saveFsType, arg("rhs"), "saves the fs type")
        .def("getOrigFsType",&oxley::OxleyDomain::getOrigFsType, "returns the fs type")
        .def("loadMesh", &oxley::OxleyDomain::loadMesh, (arg("filename")),
                "Loads a mesh (in p4est format)\n"
                ":param filename: The name of the file to load\n")
        .def("newOperator",&oxley::OxleyDomain::newSystemMatrix,
            args("row_blocksize", "row_functionspace", "column_blocksize", "column_functionspace", "type"),
            "creates a SystemMatrixAdapter stiffness matrix and initializes it with zeros\n\n"
            ":param row_blocksize:\n:type row_blocksize: ``int``\n"
            ":param row_functionspace:\n:type row_functionspace: `FunctionSpace`\n"
            ":param column_blocksize:\n:type column_blocksize: ``int``\n"
            ":param column_functionspace:\n:type column_functionspace: `FunctionSpace`\n"
            ":param type:\n:type type: ``int``")
        .def("refineMesh", &oxley::OxleyDomain::refineMesh, (args("RefinementAlgorithm")),
                "Refines the mesh.\n\n"
                ":param RefinementAlgorithm: The refinement algorithm. "
                "Accepted values are ``uniform``, ``MARE2DEM``.\n"
                ":type RefinementAlgorithm: ``str``")
        // .def("resetRhs",&oxley::OxleyDomain::resetRhs, arg("rhs"), "resets the RHS")
        .def("saveMesh", &oxley::OxleyDomain::saveMesh, (arg("filename")),
                "Saves the mesh to file using p4est format\n"
                ":param filename: The name of the output file\n")
        .def("setRefinementLevel", &oxley::OxleyDomain::setRefinementLevels, (arg("refinementlevels")),
                "Sets the number of levels of refinement\n"
                ":param refinementLevels:\ntype int: `Maximum number of levels of refinement,`\n")
        .def("setAdaptiveRefinement", &oxley::OxleyDomain::setAdaptiveRefinement, (arg("on/off")), 
                "Sets adaptive refinement on or off\n"
                ":param on/off:\n:type bool: true or false")
        .def("showTagNames",&oxley::OxleyDomain::showTagNames,
                ":return: A space separated list of tag names\n:rtype: ``string``")
        .def("updateSolutionInformation", &oxley::OxleyDomain::updateSolutionInformation, (arg("solution")),
                "Internal use. Updates the mesh with the latest solution\n")
        .def("updateMeshInformation", &oxley::OxleyDomain::updateMeshInformation, 
                "Internal use. Refines the mesh based on the solution information\n")
        .def("writeToVTK", &oxley::OxleyDomain::writeToVTK, (arg("filename"), arg("writeMesh")=false),
                "Writes the mesh to a VTK file.\n"
                ":param filename: The name of the output file\n"
                ":param writeMesh: Boolean: Only writes the mesh to file")

        ;

    // These two class exports are necessary to ensure that the extra methods
    // added by oxley make it to python. 
    class_<oxley::Brick, bases<oxley::OxleyDomain> >("OxleyBrick", "", no_init)
        .def("applyRefinementZone",&oxley::Brick::apply_refinementzone, (args("RefinementZone")),
                "Applies a RefinementZone to the Brick.\n"
                ":param RefinementZone:\n:type RefinementZone: A RefinementZone. \n")
        .def("refineBoundary", &oxley::Brick::refineBoundary, (args("boundary","dx")),
                "Refines the mesh near a boundary.\n\n"
                ":param boundary: The boundary name (top, bottom, right, left, north, south, east, west).\n"
                ":type boundary: ``str``\n"
                ":param dx: All quadrants closer to the boundary than dx will be refined.\n"
                ":type dx: ``float``")
        .def("refineRegion", &oxley::Brick::refineRegion, (arg("x0")=-1,arg("x1")=-1,arg("y0")=-1,arg("y1")=-1,arg("z0")=-1,arg("z1")=-1),
                "Refines the mesh within the interior of a region.\n"
                ":param x0:\n:type double: boundary of the region.\n"
                ":param x1:\n:type double: boundary of the region.\n"
                ":param y0:\n:type double: boundary of the region.\n"
                ":param y1:\n:type double: boundary of the region.\n")
        .def("refinePoint", &oxley::Brick::refinePoint, (arg("x0")=-1,arg("y0")=-1,arg("z0")=-1),
                "Refines the mesh around the point (x0,y0) to the level of refinement.\n"
                "set by setRefinementLevel\n"
                ":param x0:\n:type double: x coordinate of the point to be refined.\n"
                ":param y0:\n:type double: y coordinate of the point to be refined.\n")
        .def("refineSphere", &oxley::Brick::refineSphere, (arg("x0")=-1,arg("y0")=-1,arg("z0")=-1,arg("r")=-1),
                "Refines the mesh around the point (x0,y0) to the level of refinement.\n"
                "set by setRefinementLevel\n"
                ":param x0:\n:type double: x coordinate of the point to be refined.\n"
                ":param y0:\n:type double: y coordinate of the point to be refined.\n"
                ":param r: \n:type double: radius of the circle.\n")
        ;

    class_<oxley::Rectangle, bases<oxley::OxleyDomain> >("OxleyRectangle", "", no_init)
        .def("applyRefinementZone",&oxley::Rectangle::apply_refinementzone, (args("RefinementZone")),
                "Applies a RefinementZone to the Rectangle.\n"
                ":param RefinementZone:\n:type RefinementZone: A RefinementZone. \n")
        .def("refineBoundary", &oxley::Rectangle::refineBoundary, (args("boundary","dx")),
                "Refines the mesh near a boundary.\n\n"
                ":param boundary: The boundary name (top, bottom, right, left).\n"
                ":type boundary: ``str``\n"
                ":param dx: All quadrants closer to the boundary than dx will be refined.\n"
                ":type dx: ``float``")
        .def("refineRegion", &oxley::Rectangle::refineRegion, (arg("x0")=-1,arg("x1")=-1,arg("y0")=-1,arg("y1")=-1),
                "Refines the mesh within the interior of a region.\n"
                ":param x0:\n:type double: boundary of the region.\n"
                ":param x1:\n:type double: boundary of the region.\n"
                ":param y0:\n:type double: boundary of the region.\n"
                ":param y1:\n:type double: boundary of the region.\n")
        .def("refinePoint", &oxley::Rectangle::refinePoint, (arg("x0")=-1,arg("y0")=-1),
                "Refines the mesh around the point (x0,y0) to the level of refinement.\n"
                "set by setRefinementLevel\n"
                ":param x0:\n:type double: x coordinate of the point to be refined.\n"
                ":param y0:\n:type double: y coordinate of the point to be refined.\n")
        .def("refineCircle", &oxley::Rectangle::refineCircle, (arg("x0")=-1,arg("y0")=-1,arg("r")=-1),
                "Refines the mesh around the point (x0,y0) to the level of refinement.\n"
                "set by setRefinementLevel\n"
                ":param x0:\n:type double: x coordinate of the point to be refined.\n"
                ":param y0:\n:type double: y coordinate of the point to be refined.\n"
                ":param r: \n:type double: radius of the circle.\n")
        .def("interpolate", &oxley::Rectangle::interpolateAcross, (arg("target"),arg("source")),
                "Interpolates source to target\n"
                ":param source:\n:type Data: The source Data object. \n"
                ":param target:\n:type Data: The target Data object. \n")
        ;

    class_<oxley::RefinementZone>("RefinementZone", "")

    ;

    class_<oxley::RefinementZone2D, bases<oxley::RefinementZone>>("RefinementZone2D")
        .def("setRefinementLevel", &oxley::RefinementZone2D::setRefinementLevel, (arg("level")),
                "Sets the level of refinement\n"
                ":param level:\n:type int: the level of the refinement.\n")
        .def("print", &oxley::RefinementZone2D::print,
                "Prints the current queue to console\n")
        .def("remove", &oxley::RefinementZone2D::deleteFromQueue, (arg("n")),
                "Removes the n^th item from the queue\n"
                ":param n:\n:type int: the refinement to remove.\n")
        .def("refinePoint", &oxley::RefinementZone2D::refinePoint, (arg("x0"),arg("y0"),arg("level")=-1),
                "Refines the mesh around the point (x0,y0) to the level of refinement"
                "set by setRefinementLevel \n"
                ":param x0:\n:type float: x coordinate of the point to be refined.\n"
                ":param y0:\n:type float: y coordinate of the point to be refined.\n"
                ":param level:\n:type float: the level of refinement.\n")
        .def("refineRegion", &oxley::RefinementZone2D::refineRegion, (arg("x0"),arg("y0"),arg("x1"),arg("y1"),arg("level")=-1),
                "Refines the mesh around a rectangular region bound by the points (x0,y0)"
                "and (x1,y1) to the level of refinement set by setRefinementLevel\n"
                ":param x0:\n:type float: x coordinate of the upper left coordinate.\n"
                ":param y0:\n:type float: y coordinate of the upper left coordinate.\n"
                ":param x1:\n:type float: x coordinate of the lower right coordinate.\n"
                ":param y1:\n:type float: y coordinate of the lower right coordinate.\n"
                ":param level:\n:type float: the level of refinement.\n")
        .def("refineCircle", &oxley::RefinementZone2D::refineCircle, (arg("x0"),arg("y0"),arg("r"),arg("level")=-1),
                "Refines the mesh around a circular region with radius r and center"
                "and (x0,y0) to the level of refinement set by setRefinementLevel\n"
                ":param x0:\n:type float: x coordinate of the center of the circle.\n"
                ":param y0:\n:type float: y coordinate of the center of the circle.\n"
                ":param r :\n:type float: the radius of the circle.\n"
                ":param level:\n:type float: the level of refinement.\n")
        .def("refineBorder", &oxley::RefinementZone2D::refineBorder, (arg("Border"),arg("dx"),arg("level")=-1),
                "Refines the border of the mesh to depth dx to the level of refinement"
                "set by setRefinementLevel\n"
                ":param Border:\n:type string: The border to refine (top,bottom,right,left).\n"
                ":param dx:\n:type float: the depth of the refinement.\n"
                ":param level:\n:type float: the level of refinement.\n")
        .def("refineMask", &oxley::RefinementZone2D::refineMask, (args("mask")),
                "Refines the mesh in regions defined by a mask\n"
                ":param mask:\n:type Data: a mask.\n")
        ;

    class_<oxley::RefinementZone3D, bases<oxley::RefinementZone>>("RefinementZone3D")
        .def("setRefinementLevel", &oxley::RefinementZone3D::setRefinementLevel, (args("level")),
                "Sets the level of refinement\n"
                ":param level:\n:type int: the level of the refinement.\n")
        .def("print", &oxley::RefinementZone3D::print, (arg("level")),
                "Prints the current queue to console\n")
        .def("remove", &oxley::RefinementZone3D::deleteFromQueue, (arg("n")),
                "Removes the n^th item from the queue\n"
                ":param n:\n:type int: the refinement to remove.\n")
        .def("refinePoint", &oxley::RefinementZone3D::refinePoint, (arg("x0"),arg("y0"),arg("z0"),arg("level")=-1),
                "Refines the mesh around the point (x0,y0,z0) to the level of refinement"
                "set by setRefinementLevel \n"
                ":param x0:\n:type float: x coordinate of the point to be refined.\n"
                ":param y0:\n:type float: y coordinate of the point to be refined.\n"
                ":param z0:\n:type float: z coordinate of the point to be refined.\n"
                ":param level:\n:type float: the level of refinement.\n")
        .def("refineRegion", &oxley::RefinementZone3D::refineRegion, (arg("x0"),arg("y0"),arg("z0"),arg("x1"),arg("y1"),arg("z1"),arg("level")=-1),
                "Refines the mesh around a rectangular region bound by the points (x0,y0,z0)"
                "and (x1,y1,z1) to the level of refinement set by setRefinementLevel\n"
                ":param x0:\n:type float: x coordinate of the upper left coordinate.\n"
                ":param y0:\n:type float: y coordinate of the upper left coordinate.\n"
                ":param z0:\n:type float: z coordinate of the upper left coordinate.\n"
                ":param x1:\n:type float: x coordinate of the lower right coordinate.\n"
                ":param y1:\n:type float: y coordinate of the lower right coordinate.\n"
                ":param z1:\n:type float: z coordinate of the lower right coordinate.\n"
                ":param level:\n:type float: the level of refinement.\n")
        .def("refineSphere", &oxley::RefinementZone3D::refineSphere, (arg("x0"),arg("y0"),arg("z0"),arg("r"),arg("level")=-1),
                "Refines the mesh around a spherical region with radius r and center"
                "and (x0,y0,z0) to the level of refinement set by setRefinementLevel\n"
                ":param x0:\n:type float: x coordinate of the center of the circle.\n"
                ":param y0:\n:type float: y coordinate of the center of the circle.\n"
                ":param z0:\n:type float: z coordinate of the center of the circle.\n"
                ":param r :\n:type float: the radius of the circle.\n"
                ":param level:\n:type float: the level of refinement.\n")
        .def("refineBorder", &oxley::RefinementZone3D::refineBorder, (arg("Border"),arg("dx"),arg("level")=-1),
                "Refines the border of the mesh to depth dx to the level of refinement"
                "set by setRefinementLevel\n"
                ":param Border:\n:type string: The border to refine (top,bottom,right,left).\n"
                ":param dx:\n:type float: the depth of the refinement.\n"
                ":param level:\n:type float: the level of refinement.\n")
        .def("refineMask", &oxley::RefinementZone3D::refineMask, (arg("mask"),arg("level")=-1),
                "Refines the mesh in regions defined by a mask\n"
                ":param mask:\n:type Data: a mask.\n")
        ;

    class_<oxley::AbstractAssembler, oxley::Assembler_ptr, boost::noncopyable >  ("AbstractAssembler", "", no_init);

    // register_ptr_to_python<boost::shared_ptr<RefinementZone>>();
    // register_ptr_to_python<boost::shared_ptr<RefinementZone2D>>();
    // register_ptr_to_python<boost::shared_ptr<RefinementZone3D>>();

}

} //namespace oxley
