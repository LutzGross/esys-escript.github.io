
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

namespace bp = boost::python;

namespace oxley {

    /**
       \brief
       Constructor
    */
Rectangle::Rectangle(int order, dim_t n0, dim_t n1, double x0, double y0, double x1, double y1, 
    int d0, int d1): OxleyDomain(2, order){

    }

    /**
       \brief
       Destructor.
    */
Rectangle::~Rectangle(){

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

} // end of namespace oxley

