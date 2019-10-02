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

#include <iostream> //ae: temporary

#include <oxley/Brick.h>
#include <oxley/Rectangle.h>
#include <oxley/OtherAlgorithms.h>
#include <oxley/OxleyDomain.h>

#include <boost/python.hpp>
#ifdef ESYS_HAVE_BOOST_NUMPY
#include <boost/python/numpy.hpp>
#include <boost/python/numpy/dtype.hpp>
#endif

// using namespace boost::python::numpy;

namespace oxley {

escript::Domain_ptr _rectangle(int order, double _n0, double _n1,
                        const object& l0, const object& l1, int d0, int d1,
                        int periodic0, int periodic1)
{
	// Integration Order
    if (order < 2 || order > 10)
        throw OxleyException("Order must be in the range 2 to 10");

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
    if(ndx != 1)
        throw OxleyException("x has invalid dimensions.");
    if(ndy != 1)
        throw OxleyException("y has invalid dimensions");
    if(lx[0] != ly[0])
        throw OxleyException("x and y are different lengths.");

    // This structure is used to store the information used during the algorithm
    // a pointer is attached to p4est and then passed around between functions
    p4estData *forestData = (p4estData *) domain->p4est->user_pointer;
    addSurfaceData * surfacedata = new addSurfaceData;
    forestData->info = surfacedata;
    surfacedata->x.clear();
    surfacedata->y.clear();
    surfacedata->x.resize(lx[0],-1.0);
    surfacedata->y.resize(ly[0],-1.0);
    double* px = reinterpret_cast<double*>(x.get_data());
    double* py = reinterpret_cast<double*>(y.get_data());
    // Note: boost::python::numpy::ndarray is not threadsafe
    for(long xx = 0; xx < *lx; xx++)
    {
        surfacedata->x[xx] = *(px + xx);
        surfacedata->y[xx] = *(py + xx);
    }

    if(!domain->getDescription().compare("oxley::rectangle"))
    {
        addSurface(domain);
        delete surfacedata;
    }
    else
    {
        delete surfacedata;
        throw OxleyException("Invalid domain. Was expecting an oxley::rectangle.");
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
    if(ndx != 1)
        throw OxleyException("x has invalid dimensions.");
    if(ndy != 1)
        throw OxleyException("y has invalid dimensions");
    if(ndz != 1)
        throw OxleyException("z has invalid dimensions");
    if(lz[0] != lx[0]*ly[0])
        throw OxleyException("z i snot of length x*y.");

    // This structure is used to store the information used during the algorithm
    // a pointer is attached to p4est and then passed around between functions
    p8estData *forestData = (p8estData *) domain->p8est->user_pointer;
    addSurfaceData * surfacedata = new addSurfaceData;
    forestData->info = surfacedata;
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
    for(long xx = 0; xx < *lx; xx++)
    {
        surfacedata->x[xx] = *(px + xx);
        surfacedata->y[xx] = *(py + xx);
    }
    for(long xx = 0; xx < *lz; xx++){
        surfacedata->z[xx] = *(pz + xx);
    }

    if(!domain->getDescription().compare("oxley::rectangle"))
    {
        addSurface(domain);
        delete surfacedata;
    }
    else
    {
        delete surfacedata;
        throw OxleyException("Invalid domain. Was expecting an oxley::rectangle.");
    }
}
#endif

BOOST_PYTHON_MODULE(oxleycpp)
{

    def("Rectangle", oxley::_rectangle, (arg("order"),
    arg("n0"),arg("n1"),
    arg("l0")=1.0,arg("l1")=1.0,
    arg("d0")=-1,arg("d1")=-1,
    arg("periodic0")=0,arg("periodic1")=0),
    "Creates a rectangular p4est mesh with n0 x n1 elements over the rectangle [0,l0] x [0,l1].\n\n"
    ":param order: order of the elements: ``int``\n"
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
        .def("getDescription", &oxley::OxleyDomain::getDescription,
                "Prints out a description of the mesh.")
        .def("writeToVTK", &oxley::OxleyDomain::writeToVTK, (arg("filename"), arg("writeTagInfo")=false),
                "Writes the mesh to a VTK file.\n"
                ":param filename: The name of the output file\n"
                ":param writeTagInfo: Boolean: Whether to output tag info")
        .def("setRefinementLevels", &oxley::OxleyDomain::setRefinementLevels, (arg("refinementlevels")),
                "Sets the number of levels of refinement\n"
                ":param refinementLevels:\ntype int: `Maximum number of levels of refinement,`\n")
        .def("refine", &oxley::OxleyDomain::refineMesh, (args("maxRecursion","RefinementAlgorithm")),
                "Refines the mesh.\n"
                ":param maxRecursion:\n:type int: `Maximum number of levels of refinement,`\n"
                ":param RefinementAlgorithm:\n:type string: `The refinement algorithm \n"
                "       accepted values are \"uniform\"")
        .def("getNumVertices", &oxley::OxleyDomain::getNumVertices,
                "Returns the number of corners in the mesh.");

    class_<oxley::Rectangle, bases<oxley::OxleyDomain> > ("OxleyRectangle", "", no_init);
    class_<oxley::Brick, bases<oxley::OxleyDomain> > ("OxleyBrick", "", no_init);
}

} //namespace oxley
