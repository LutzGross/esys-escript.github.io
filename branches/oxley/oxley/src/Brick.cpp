
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

    }

    /**
       \brief
       Destructor.
    */
Brick::~Brick(){

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

