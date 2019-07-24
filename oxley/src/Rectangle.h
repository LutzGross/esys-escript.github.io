
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

#ifndef __OXLEY_RECTANGLE_H__
#define __OXLEY_RECTANGLE_H__

#include <escript/EsysMPI.h>
#include <escript/SubWorld.h>

#include <oxley/OxleyDomain.h>

#include <p4est.h>
#include <p4est_connectivity.h>

#include <boost/python.hpp>
#ifdef ESYS_HAVE_BOOST_NUMPY
#include <boost/python/numpy.hpp>
#endif

using namespace boost::python;

namespace oxley {

/**
   \brief
   Rectangle is the 2-dimensional implementation of a Oxleydomain.
*/
class Rectangle: public OxleyDomain
{
public:

    /**
       \brief creates a rectangular mesh with n0 x n1 elements over the
              rectangle [x0,x1] x [y0,y1].
       \param 
    */
    Rectangle(int order, dim_t n0, dim_t n1, double x0, double y0, double x1, double y1, int d0, int d1);

    /**
       \brief creates a rectangular mesh from numpy arrays [x,y].
            Requires boost numpy 
       \param 
    */
#ifdef ESYS_HAVE_BOOST_NUMPY
    Rectangle(int order, dim_t n0, dim_t n1, boost::python::numpy::ndarray x, boost::python::numpy::ndarray y);
#endif
    
    /**
       \brief
       Destructor.
    */
    ~Rectangle();

    /**
       \brief
       returns a description for this domain
    */
    virtual std::string getDescription() const;

    /**
       \brief
       dumps the mesh to a file with the given name
       \param filename The name of the output file
    */
    virtual void dump(const std::string& filename) const;

    /**
       \brief
       writes the current mesh to a file with the given name
       \param filename The name of the file to write to
    */
    virtual void write(const std::string& filename) const;
    
    /**
       \brief
       interpolates data given on source onto target where source and target
       are given on different domains
    */
    virtual void interpolateAcross(escript::Data& target,
                                   const escript::Data& source) const;

    /**
       \brief
       determines whether interpolation from source to target is possible
    */
    virtual bool probeInterpolationAcross(int, const escript::AbstractDomain&,
            int) const;

    /**
       \brief
       returns true if this rank owns the sample id.
    */
    virtual bool ownSample(int fs_code, index_t id) const;

    /**
       \brief
       copies the surface normals at data points into out. The actual function
       space to be considered is defined by out. out has to be defined on this
       domain.
    */
    virtual void setToNormal(escript::Data& out) const;

    /**
       \brief
       copies the size of samples into out. The actual function space to be
       considered is defined by out. out has to be defined on this domain.
    */
    virtual void setToSize(escript::Data& out) const;

    /**
     * \brief 
       Returns a Data object filled with random data passed through filter.
    */ 
    virtual escript::Data randomFill(const escript::DataTypes::ShapeType& shape,
       const escript::FunctionSpace& what, long seed, const boost::python::tuple& filter) const;

    /**
       \brief
       returns the array of reference numbers for a function space type
       \param fsType The function space type
    */
    const dim_t* borrowSampleReferenceIDs(int fsType) const;


private:

    // This object records the connectivity of the p4est quadrants
    p4est_connectivity_t *connectivity;

    // A p4est
    p4est_t *p4est;

    // origin of domain
    double m_origin[2];

    // side lengths of domain
    double m_length[2];

    // number of spatial subdivisions
    int m_NX[2];

    // total number of elements in each dimension
    dim_t m_gNE[2];

    // number of elements for this rank in each dimension including shared
    dim_t m_NE[2];

    // periodic boundary conditions
    int periodic[2] {false, false};

};

////////////////////////////// inline methods ////////////////////////////////

} //end namespace


#endif //__OXLEY_RECTANGLE_H__
