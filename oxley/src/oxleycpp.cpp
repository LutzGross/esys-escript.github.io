/*****************************************************************************
*
* Copyright (c) 2003-2019 by The University of Queensland
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

#include <oxley/Brick.h>
#include <oxley/Rectangle.h>
#include <oxley/OtherAlgorithms.h>
#include <oxley/OxleyDomain.h>

#include <boost/python.hpp>
#ifdef ESYS_HAVE_BOOST_NUMPY
#include <boost/python/numpy.hpp>
#include <boost/python/numpy/dtype.hpp>
#endif

// #include <p4est_algorithms.h> //aeae

// using namespace boost::python::numpy;

namespace oxley {

escript::Domain_ptr _rectangle(double _n0, double _n1,
                        const object& l0, const object& l1, int d0, int d1,
                        int periodic0, int periodic1)
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

    return escript::Domain_ptr(new Rectangle(order, n0,n1, x0,y0, x1,y1, d0,d1, periodic0, periodic1));
}


escript::Domain_ptr _brick(int order, double _n0, double _n1, double _n2,
                        const object& l0, const object& l1, const object& l2,
                        int d0, int d1, int d2,
                        int periodic0, int periodic1, int periodic2)
{
    // Integration Order
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

    return escript::Domain_ptr(new Brick(order, n0,n1,n2, x0,y0,z0, x1,y1,z1, d0,d1,d2, periodic0,periodic1,periodic2));
}

#ifdef ESYS_HAVE_BOOST_NUMPY
void _addCurve(OxleyDomainRect_ptr domain,
    boost::python::numpy::ndarray &x,
    boost::python::numpy::ndarray &y)
{
    // Check that the x and z ndarrays are of reasonable dimensions
    int ndx = x.get_nd();
    int ndy = y.get_nd();
    const long int* lx = x.get_shape();
    const long int* ly = y.get_shape();

    if(x.get_dtype() != boost::python::numpy::dtype::get_builtin<double>())
        throw OxleyException("x must be an array of doubles.");
    if(y.get_dtype() != boost::python::numpy::dtype::get_builtin<double>())
        throw OxleyException("y must be an array of doubles.");
    if(ndx != 1)
        throw OxleyException("x has invalid dimensions.");
    if(ndy != 1)
        throw OxleyException("y has invalid dimensions");
    if(lx[0] != ly[0])
        throw OxleyException("x and y are different lengths.");

    // This structure is used to store the information used during the algorithm
    // a pointer is attached to p4est and then passed around between functions
    addSurfaceData * surfacedata = new addSurfaceData;
    p4estData *forestData = (p4estData *) domain->borrow_forestData();
    forestData->assign_info(surfacedata);
    domain->set_temp_data(surfacedata);
    surfacedata->x.clear();
    surfacedata->y.clear();
    surfacedata->x.resize(lx[0],-1.0);
    surfacedata->y.resize(ly[0],-1.0);
    double* px = reinterpret_cast<double*>(x.get_data());
    double* py = reinterpret_cast<double*>(y.get_data());
    // Note: boost::python::numpy::ndarray is not threadsafe
    for(long i = 0; i < *lx; i++)
        surfacedata->x[i] = *(px + i);
    for(long i = 0; i < *ly; i++)
        surfacedata->y[i] = *(py + i);

    if(!domain->getDescription().compare("oxley::rectangle"))
    {
        addSurface(domain);
        delete surfacedata;

    }
    else
    {
        delete surfacedata;
        std::string message = "Invalid domain. Was expecting an oxley::rectangle, not a " + domain->getDescription();
        throw OxleyException(message);
    }
}
#endif

#ifdef ESYS_HAVE_BOOST_NUMPY
void _addSurface(OxleyDomainBrick_ptr domain,
    boost::python::numpy::ndarray &x,
    boost::python::numpy::ndarray &y,
    boost::python::numpy::ndarray &z)
{
    // Check that the x and z ndarrays are of reasonable dimensions
    int ndx = x.get_nd();
    int ndy = y.get_nd();
    int ndz = z.get_nd();
    const long int* lx = x.get_shape();
    const long int* ly = y.get_shape();
    const long int* lz = z.get_shape();
    if(x.get_dtype() != boost::python::numpy::dtype::get_builtin<double>())
        throw OxleyException("x must be an array of doubles.");
    if(y.get_dtype() != boost::python::numpy::dtype::get_builtin<double>())
        throw OxleyException("y must be an array of doubles.");
    if(z.get_dtype() != boost::python::numpy::dtype::get_builtin<double>())
        throw OxleyException("z must be an array of doubles.");
    if(ndx != 1)
        throw OxleyException("x has invalid dimensions.");
    if(ndy != 1)
        throw OxleyException("y has invalid dimensions");
    if(ndz != 1)
        throw OxleyException("z has invalid dimensions");
    if(lz[0] != lx[0]*ly[0])
        throw OxleyException("z is not of length x*y.");

    // This structure is used to store the information used during the algorithm
    // a pointer is attached to p4est and then passed around between functions
    addSurfaceData * surfacedata = new addSurfaceData;
    p4estData *forestData = (p4estData *) domain->borrow_forestData();
    // forestData->info = surfacedata;
    forestData->assign_info(surfacedata);
    domain->set_temp_data(surfacedata);
    surfacedata->x.clear();
    surfacedata->y.clear();
    surfacedata->z.clear();
    surfacedata->x.resize(lx[0],-1.0);
    surfacedata->y.resize(ly[0],-1.0);
    surfacedata->z.resize(lz[0],-1.0);
    double* px = reinterpret_cast<double*>(x.get_data());
    double* py = reinterpret_cast<double*>(y.get_data());
    double* pz = reinterpret_cast<double*>(z.get_data());
    std::vector<double> vx(*lx);
    std::vector<double> vy(*ly);
    std::vector<double> vz(*lz);
    // Note: boost::python::numpy::ndarray is not threadsafe
    for(long i = 0; i < *lx; i++)
        surfacedata->x[i] = *(px + i);
    for(long i = 0; i < *ly; i++)
        surfacedata->y[i] = *(py + i);
    for(long i = 0; i < *lz; i++)
        surfacedata->z[i] = *(pz + i);

    if(!domain->getDescription().compare("oxley::brick"))
    {
        addSurface(domain);
        delete surfacedata;
    }
    else
    {
        delete surfacedata;
        std::string message = "Invalid domain. Was expecting an oxley::brick, not a " + domain->getDescription();
        throw OxleyException(message);
    }
}
#endif

BOOST_PYTHON_MODULE(oxleycpp)
{

    def("Rectangle", oxley::_rectangle, (
    arg("n0"),arg("n1"),
    arg("l0")=1.0,arg("l1")=1.0,
    arg("d0")=-1,arg("d1")=-1,
    arg("periodic0")=0,arg("periodic1")=0),
    "Creates a rectangular p4est mesh with n0 x n1 elements over the rectangle [0,l0] x [0,l1].\n\n"
    ":param n0: number of elements in direction 0\n:type n0: ``int``\n"
    ":param n1: number of elements in direction 1\n:type n1: ``int``\n"
    ":param l0: length of side 0 or coordinate range of side 0\n:type l0: ``float`` or ``tuple``\n"
    ":param l1: length of side 1 or coordinate range of side 1\n:type l1: ``float`` or ``tuple``\n"
    ":param d0: number of subdivisions in direction 0\n:type d0: ``int``\n"
    ":param d1: number of subdivisions in direction 1\n:type d1: ``int``");

    def("Brick", oxley::_brick, (arg("order"),
    arg("n0"),arg("n1"),arg("n2"),
    arg("l0")=1.0,arg("l1")=1.0,arg("l2")=1.0,
    arg("d0")=-1,arg("d1")=-1,arg("d2")=-1,
    arg("periodic0")=0,arg("periodic1")=0,arg("periodic2")=0),
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
    ":param d2: number of subdivisions in direction 2\n:type d2: ``int``");

#ifdef ESYS_HAVE_BOOST_NUMPY
    // def("addSurface", oxley::_addSurface, (arg("domain"),arg("x"),arg("y"),arg("z")));
    def("addSurface", oxley::_addCurve, (arg("domain"),arg("x"),arg("y")));
    def("addSurface", oxley::_addSurface, (arg("domain"),arg("x"),arg("y"),arg("z")));
#endif

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
                "Returns the number of corners in the mesh.")
        .def("getSystemMatrixTypeId",&oxley::OxleyDomain::getSystemMatrixTypeId,
            args("options"),
            ":return: the identifier of the matrix type to be used for the global stiffness matrix when particular solver options are used.\n"
            ":rtype: ``int``\n"
            ":param options:\n:type options: `SolverBuddy`\n")
        .def("getTransportTypeId",&oxley::OxleyDomain::getTransportTypeId,
            args("solver", "preconditioner", "package", "symmetry"),
            ":return: the identifier of the transport problem type to be used when a particular solver, preconditioner, package and symmetric matrix is used.\n"
            ":rtype: ``int``\n"
            ":param solver:\n:type solver: ``int``\n"
            ":param preconditioner:\n:type preconditioner: ``int``\n"
            ":param package:\n:type package: ``int``\n"
            ":param symmetry:\n:type symmetry: ``int``"
            )
        .def("getX",&oxley::OxleyDomain::getX, ":return: locations in the FEM nodes\n\n"
            ":rtype: `Data`")
        .def("newOperator",&oxley::OxleyDomain::newSystemMatrix,
            args("row_blocksize", "row_functionspace", "column_blocksize", "column_functionspace", "type"),
            "creates a SystemMatrixAdapter stiffness matrix and initializes it with zeros\n\n"
            ":param row_blocksize:\n:type row_blocksize: ``int``\n"
            ":param row_functionspace:\n:type row_functionspace: `FunctionSpace`\n"
            ":param column_blocksize:\n:type column_blocksize: ``int``\n"
            ":param column_functionspace:\n:type column_functionspace: `FunctionSpace`\n"
            ":param type:\n:type type: ``int``")
        .def("refineMesh", &oxley::OxleyDomain::refineMesh, (args("RefineLevel","RefinementAlgorithm")),
                "Refines the mesh.\n"
                ":param RefineLevel:\n:type int: `Maximum levels of refinement,`\n"
                ":param RefinementAlgorithm:\n:type string: `The refinement algorithm \n"
                "       accepted values are \"uniform\"")
        .def("refineBoundary", &oxley::OxleyDomain::refineBoundary, (args("boundary","dx")),
                "Refines the mesh near a boundary.\n"
                ":param boundary:\n:type string: `The boundary (n,s,e,w) \n"
                ":param dx:\n:type double: all quadrants closer to the boundary than dx will be refined. ")
        .def("setRefinementLevels", &oxley::OxleyDomain::setRefinementLevels, (arg("refinementlevels")),
                "Sets the number of levels of refinement\n"
                ":param refinementLevels:\ntype int: `Maximum number of levels of refinement,`\n")
        .def("writeToVTK", &oxley::OxleyDomain::writeToVTK, (arg("filename"), arg("writeMesh")=false),
                "Writes the mesh to a VTK file.\n"
                ":param filename: The name of the output file\n"
                ":param writeMesh: Boolean: Only writes the mesh to file")
        ;

    // These two class exports are necessary to ensure that the extra methods
    // added by oxley make it to python. 
    class_<oxley::Brick, bases<oxley::OxleyDomain> >("OxleyBrick", "", no_init);
    class_<oxley::Rectangle, bases<oxley::OxleyDomain> >("OxleyRectangle", "", no_init);
    class_<oxley::AbstractAssembler, oxley::Assembler_ptr, boost::noncopyable >  ("AbstractAssembler", "", no_init);
}

} //namespace oxley
