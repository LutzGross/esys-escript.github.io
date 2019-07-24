
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

namespace bp = boost::python;

namespace oxley {

    /**
       \brief
       Constructor
    */
Brick::Brick(int order, dim_t n0, dim_t n1, dim_t n2, double x0, double y0, double z0,
      double x1, double y1, double z1, int d0, int d1, int d2): OxleyDomain(2, order){


    // Possible error: User passes invalid values for the dimensions
    if(n0 <= 0 || n1 <= 0 || n2 <= 0)
        throw OxleyException("Number of elements in each spatial dimension must be positive");

    // Ignore d0 and d1 if we are running in serial
    m_mpiInfo = escript::makeInfo(MPI_COMM_WORLD);
    if(m_mpiInfo->size == 1) {
        d0=1;
        d1=1;
        d2=1;
    }

    // ensure number of subdivisions is valid and nodes can be distributed
    // among number of ranks
    if(d0*d1*d2 != m_mpiInfo->size)
        throw OxleyException("Invalid number of spatial subdivisions");

    //Create a connectivity
    connectivity = p8est_connectivity_new_brick((int) n0, (int) n1, (int) n2, 
        periodic[0], periodic[1], periodic[2]);

    // Create a forest that is not refined; it consists of the root octant.
    p8est = p8est_new(m_mpiInfo->comm, connectivity, 0, NULL, NULL);

    // Record the physical dimensions of the domain and the location of the origin
    m_origin[0] = x0;
    m_origin[1] = y0;
    m_origin[2] = z0;
    m_length[0] = x1-x0;
    m_length[1] = y1-y0;
    m_length[2] = z1-z0;

    // number of elements in each dimension 
    m_gNE[0] = n0;
    m_gNE[1] = n1;
    m_gNE[2] = n2;

    //number of elements for this rank in each dimension including shared
    m_NX[0] = d0;
    m_NX[1] = d1;
    m_NX[2] = d2;

    // local number of elements
    m_NE[0] = m_gNE[0] / d0;
    m_NE[1] = m_gNE[1] / d1;
    m_NE[2] = m_gNE[2] / d2;

    }

    /**
       \brief
       Destructor.
    */
Brick::~Brick(){

    p8est_destroy(p8est);
    p8est_connectivity_destroy(connectivity);

    }

    /**
       \brief
       returns a description for this domain
    */
std::string Brick::getDescription() const{
        return "oxley::Brick";
    }


/**
   \brief
   dumps the mesh to a file with the given name
   \param filename The name of the output file
*/
void Brick::dump(const std::string& filename) const
{
    throw OxleyException("dump: not supported");
}

/**
   \brief
   writes the current mesh to a file with the given name
   \param filename The name of the file to write to
*/
void Brick::write(const std::string& filename) const
{
    throw OxleyException("write: not supported");
}

bool Brick::probeInterpolationAcross(int fsType_source,
        const escript::AbstractDomain& domain, int fsType_target) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

void Brick::interpolateAcross(escript::Data& target, const escript::Data& source) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

void Brick::setToNormal(escript::Data& out) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

void Brick::setToSize(escript::Data& out) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

bool Brick::ownSample(int fsType, index_t id) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}


/* This is a wrapper for filtered (and non-filtered) randoms
 * For detailed doco see randomFillWorker
*/
escript::Data Brick::randomFill(const escript::DataTypes::ShapeType& shape,
                                    const escript::FunctionSpace& fs,
                                    long seed, const bp::tuple& filter) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}


const dim_t* Brick::borrowSampleReferenceIDs(int fsType) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}



} // end of namespace oxley

