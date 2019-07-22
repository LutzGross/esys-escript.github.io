
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
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

// #include <oxley/OxleyDomain.h>
#include <oxley/Rectangle.h>

// #include <escript/ExceptionTranslators.h>

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/module.hpp>
#include <boost/python/detail/defaults_gen.hpp>
#include <boost/version.hpp>

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


BOOST_PYTHON_MODULE(oxleycpp)
{
    // class_<oxley::OxleyDomain, bases<escript::AbstractContinuousDomain>, boost::noncopyable >

    def("Rectangle", oxley::_rectangle, (arg("order"),arg("n0"),arg("n1"),arg("l0")=1.0,arg("l1")=1.0,arg("d0")=-1,arg("d1")=-1),
"Creates a rectangular p4est mesh with n0 x n1 elements over the rectangle [0,l0] x [0,l1].\n\n"
":param order: order of the elements: ``int``\n"
":param n0: number of elements in direction 0\n:type n0: ``int``\n"
":param n1: number of elements in direction 1\n:type n1: ``int``\n"
":param l0: length of side 0 or coordinate range of side 0\n:type l0: ``float`` or ``tuple``\n"
":param l1: length of side 1 or coordinate range of side 1\n:type l1: ``float`` or ``tuple``\n"
":param d0: number of subdivisions in direction 0\n:type d0: ``int``\n"
":param d1: number of subdivisions in direction 1\n:type d1: ``int``");

	def("getDescription", &oxley::OxleyDomain::getDescription,
            ":return: a description for this domain\n:rtype: ``string``");
        
}

} //namespace oxley