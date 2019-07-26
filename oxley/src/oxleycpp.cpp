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

#include <boost/python.hpp>
#ifdef ESYS_HAVE_BOOST_NUMPY
#include <boost/python/numpy.hpp>
#include <boost/python/numpy/dtype.hpp>
#endif

using namespace boost::python;

namespace oxley {

escript::Domain_ptr _rectangle(int order, double _n0, double _n1,
                        const object& l0, const object& l1, int d0, int d1)
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

    return escript::Domain_ptr(new Rectangle(order, n0,n1, x0,y0, x1,y1, d0,d1));
}

#ifdef ESYS_HAVE_BOOST_NUMPY
escript::Domain_ptr _rectangle_numpy(int order, 
    boost::python::numpy::ndarray ndx, boost::python::numpy::ndarray ndy, int d0, int d1)
{
    // Integration Order
    if(order < 2 || order > 10)
        throw OxleyException("Order must be in the range 2 to 10.");

    // boost::python::list xtemp = ndx;
    // boost::python::list ytemp = ndy;

    // int n0 = ndx.get_nd();
    // int n1 = ndy.get_nd();
    // boost::python::numpy::dtype xtype = ndx.get_dtype();
    // boost::python::numpy::dtype ytype = ndy.get_dtype();

    // // ndx.squeeze();
    // char* xdatapointer = ndx.get_data();
    // int* xpoint = reinterpret_cast<int*>(ndx.get_data());
    // for(int i = 0; i < n0; i++){
    //     std::cout << xpoint[i] << std::endl;
    // }

    int n0=10,n1=10,x0=0, y0=0, x1=1, y1=1; // ae: This is temporary
    return escript::Domain_ptr(new Rectangle(order, n0,n1, x0,y0, x1,y1, d0,d1));
}
#else
escript::Domain_ptr _rectangle_numpy(int order, 
    boost::python::dict ndx, boost::python::dict ndy, int d0, int d1)
{
    throw OxleyException("This requires the boost::numpy libraries.");
}
#endif

escript::Domain_ptr _brick(int order, double _n0, double _n1, double _n2,
                        const object& l0, const object& l1, const object& l2, 
                        int d0, int d1, int d2)
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

    return escript::Domain_ptr(new Brick(order, n0,n1,n2, x0,y0,z0, x1,y1,z1, d0,d1,d2));
}

#ifdef ESYS_HAVE_BOOST_NUMPY
escript::Domain_ptr _brick_numpy(int order, 
    boost::python::numpy::ndarray ndx, boost::python::numpy::ndarray ndy, boost::python::numpy::ndarray ndz,  int d0, int d1, int d2)
{
    // Integration Order
    if(order < 2 || order > 10)
        throw OxleyException("Order must be in the range 2 to 10.");

    int n0=10,n1=10,n2=10, x0=0,y0=0,z0=0, x1=1,y1=1,z1=1; // ae: This is temporary
    return escript::Domain_ptr(new Brick(order, n0,n1,n2, x0,y0,z0, x1,y1,z1, d0,d1,d2));
}
#else
escript::Domain_ptr _rectangle_numpy(int order, 
    boost::python::numpy::ndarray ndx, boost::python::numpy::ndarray ndy, boost::python::numpy::ndarray ndz,  int d0, int d1, int d2)
{
    throw OxleyException("This requires the boost::numpy libraries.");
}
#endif

BOOST_PYTHON_MODULE(oxleycpp)
{
    // Initialise numpy, if it wasn't initialised somewhere else
    boost::python::numpy::initialize();

    def("Rectangle", oxley::_rectangle, (arg("order"),arg("n0"),arg("n1"),arg("l0")=1.0,arg("l1")=1.0,arg("d0")=-1,arg("d1")=-1),
"Creates a rectangular p4est mesh with n0 x n1 elements over the rectangle [0,l0] x [0,l1].\n\n"
":param order: order of the elements: ``int``\n"
":param n0: number of elements in direction 0\n:type n0: ``int``\n"
":param n1: number of elements in direction 1\n:type n1: ``int``\n"
":param l0: length of side 0 or coordinate range of side 0\n:type l0: ``float`` or ``tuple``\n"
":param l1: length of side 1 or coordinate range of side 1\n:type l1: ``float`` or ``tuple``\n"
":param d0: number of subdivisions in direction 0\n:type d0: ``int``\n"
":param d1: number of subdivisions in direction 1\n:type d1: ``int``");

    def("Rectangle", oxley::_rectangle_numpy, (arg("order"),arg("ndx"),arg("ndy"),arg("d0")=-1,arg("d1")=-1),
"Creates a rectangular p4est mesh with n0 x n1 elements over the rectangle [0,l0] x [0,l1].\n\n"
":param order: order of the elements: ``int``\n"
":param ndx: a numpy array with the coordinates of the nodes in direction 0\n:type ndx: ``numpy::ndarray``\n"
":param ndy: a numpy array with the coordinates of the nodes in direction 1\n:type ndy: ``numpy::ndarray``\n"
":param d0: number of subdivisions in direction 0\n:type d0: ``int``\n"
":param d1: number of subdivisions in direction 1\n:type d1: ``int``");

    def("Brick", oxley::_brick, (arg("order"),arg("n0"),arg("n1"),arg("n2"),arg("l0")=1.0,arg("l1")=1.0,arg("l2")=1.0,arg("d0")=-1,arg("d1")=-1,arg("d2")=-1),
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

    def("Rectangle", oxley::_brick_numpy, (arg("order"),arg("ndx"),arg("ndy"),arg("ndz"),arg("d0")=-1,arg("d1")=-1,arg("d2")=-1),
"Creates a rectangular p4est mesh with n0 x n1 x n2 elements over the rectangle [0,l0] x [0,l1] x [0,l2].\n\n"
":param order: order of the elements: ``int``\n"
":param ndx: a numpy array with the coordinates of the nodes in direction 0\n:type ndx: ``numpy::ndarray``\n"
":param ndy: a numpy array with the coordinates of the nodes in direction 1\n:type ndy: ``numpy::ndarray``\n"
":param ndz: a numpy array with the coordinates of the nodes in direction 1\n:type ndz: ``numpy::ndarray``\n"
":param d0: number of subdivisions in direction 0\n:type d0: ``int``\n"
":param d1: number of subdivisions in direction 1\n:type d1: ``int``\n"
":param d1: number of subdivisions in direction 1\n:type d2: ``int``");

    class_<oxley::OxleyDomain, bases<escript::AbstractContinuousDomain>, boost::noncopyable >
        ("OxleyDomain", "", no_init)
        .def("getDescription", &oxley::OxleyDomain::getDescription,
                "Prints out a description of the mesh.")
        .def("writeToVTK", &oxley::OxleyDomain::writeToVTK, (arg("filename")),
                "Writes the mesh to a VTK file.")
        .def("refine", &oxley::OxleyDomain::refineMesh, (args("maxRecursion","RefinementAlgorithm")),
                "Refines the mesh.\n"
                ":param maxRecursion:\n:type int: `Maximum number of levels of refinement,`\n"
                ":param RefinementAlgorithm:\n:type string: `The refinement algorithm ('GridCorners')`\n");

    class_<oxley::Rectangle, bases<oxley::OxleyDomain> >("OxleyRectangle", "", no_init);
    class_<oxley::Brick, bases<oxley::OxleyDomain> >("OxleyBrick", "", no_init);
}

} //namespace oxley

