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

#include <oxley/Rectangle.h>

#include <p4est_vtk.h>

// ae - this is temporary
inline void xx(std::string message){
    std::cout << message << std::endl;
}


namespace bp = boost::python;



namespace oxley {

    /**
       \brief
       Constructor
    */
Rectangle::Rectangle(int order, dim_t n0, dim_t n1, 
    double x0, double y0, double x1, double y1, int d0, int d1): 
    OxleyDomain(2, order){

    // Possible error: User passes invalid values for the dimensions
    if(n0 <= 0 || n1 <= 0)
        throw OxleyException("Number of elements in each spatial dimension must be positive");

    // Ignore d0 and d1 if we are running in serial
    m_mpiInfo = escript::makeInfo(MPI_COMM_WORLD);
    if(m_mpiInfo->size == 1) {
        d0=1;
        d1=1;
    }

    // ensure number of subdivisions is valid and nodes can be distributed
    // among number of ranks
    if(d0*d1 != m_mpiInfo->size)
        throw OxleyException("Invalid number of spatial subdivisions");

    //Create a connectivity
    connectivity = p4est_connectivity_new_brick((int) n0, (int) n1, periodic[0], periodic[1]);

    // Create a forest that is not refined; it consists of the root octant. 
    p4est = p4est_new(m_mpiInfo->comm, connectivity, 0, NULL, NULL);

    // Record the physical dimensions of the domain and the location of the origin
    m_origin[0] = x0;
    m_origin[1] = y0;
    m_length[0] = x1-x0;
    m_length[1] = y1-y0;

    // number of elements in each dimension 
    m_gNE[0] = n0;
    m_gNE[1] = n1;

    //number of elements for this rank in each dimension including shared
    m_NX[0] = d0;
    m_NX[1] = d1;

    // local number of elements
    m_NE[0] = m_gNE[0] / d0;
    m_NE[1] = m_gNE[1] / d1;

    }

    /**
       \brief
       Destructor.
    */
Rectangle::~Rectangle(){

    p4est_destroy(p4est);
    p4est_connectivity_destroy(connectivity);

    }

    /**
       \brief
       returns a description for this domain
    */
std::string Rectangle::getDescription() const{
        return "oxley::Rectangle";
    }


/**
   \brief
   dumps the mesh to a file with the given name
   \param filename The name of the output file
*/
void Rectangle::dump(const std::string& filename) const
{
    throw OxleyException("dump: not supported");
}

/**
   \brief
   writes the current mesh to a file with the given name
   \param filename The name of the file to write to
*/
void Rectangle::write(const std::string& filename) const
{
    throw OxleyException("write: not supported");
}

bool Rectangle::probeInterpolationAcross(int fsType_source,
        const escript::AbstractDomain& domain, int fsType_target) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

void Rectangle::interpolateAcross(escript::Data& target, const escript::Data& source) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

void Rectangle::setToNormal(escript::Data& out) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

void Rectangle::setToSize(escript::Data& out) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

bool Rectangle::ownSample(int fsType, index_t id) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}


/* This is a wrapper for filtered (and non-filtered) randoms
 * For detailed doco see randomFillWorker
*/
escript::Data Rectangle::randomFill(const escript::DataTypes::ShapeType& shape,
                                    const escript::FunctionSpace& fs,
                                    long seed, const bp::tuple& filter) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}


const dim_t* Rectangle::borrowSampleReferenceIDs(int fsType) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

void Rectangle::writeToVTK(std::string filename) const
{
    // Check that the filename is reasonable
    int stringlength=filename.length();
    if(stringlength>4){
        if(filename.compare(stringlength-4,4,".vtk") != 0)
            filename.append(".vtk");
    } else {
        filename.append(".vtk");
    }

    // Write to file
    const char * name = filename.c_str();
    p4est_vtk_write_file(p4est, NULL, name);
}

} // end of namespace oxley
