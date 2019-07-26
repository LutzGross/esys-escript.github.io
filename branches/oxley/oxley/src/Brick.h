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

#ifndef __OXLEY_BRICK_H__
#define __OXLEY_BRICK_H__

#include <escript/EsysMPI.h>
#include <escript/SubWorld.h>

#include <oxley/OxleyDomain.h>

#include <p8est.h>
#include <p8est_connectivity.h>

#include <boost/python.hpp>

using namespace boost::python;

namespace oxley {

/**
   \brief
   Brick is the 2-dimensional implementation of an Oxleydomain.
*/
class Brick: public OxleyDomain
{
public:

    /**
       \brief creates a rectangular mesh with n0 x n1 x n2 elements over the
              rectangle [x0,x1] x [y0,y1] x [z0,z1].
       \param 
    */
    // Brick();

    Brick(int order, dim_t n0, dim_t n1, dim_t n2, double x0, double y0, double z0,
      double x1, double y1, double z1, int d0, int d1, int d2);

    /**
       \brief
       Destructor.
    */
    ~Brick();

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

    /**
       \brief
       writes the mesh to a VTK file
       \param filename The file name
    */
    virtual void writeToVTK(std::string filename) const;

    /**
       \brief
       refines the mesh using enum RefinementAlgorithm
    */
    virtual void refineMesh(int maxRecursion, std::string RefinementAlgorithm);

private:

    // This object records the connectivity of the p8est quadrants
    p8est_connectivity_t *connectivity;

    // A p8est
    p8est_t *p8est;

    // origin of domain
    double m_origin[3];

    // side lengths of domain
    double m_length[3];

    // number of spatial subdivisions
    int m_NX[3];

    // total number of elements in each dimension
    dim_t m_gNE[3];

    // number of elements for this rank in each dimension including shared
    dim_t m_NE[3];

    // periodic boundary conditions
    int periodic[3] {false, false, false};


};

////////////////////////////// inline methods ////////////////////////////////

} //end namespace


#endif //__OXLEY_BRICK_H__
