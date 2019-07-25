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

#include <string>

#include <oxley/OxleyDomain.h>

#include <escript/Data.h>

namespace bp = boost::python;

using namespace std;

namespace oxley {

    /**
       \brief
       Constructor 
    */
    OxleyDomain::OxleyDomain(dim_t dim, int order)
    {
        
    }

    /**
       \brief
       Destructor
    */
    OxleyDomain::~OxleyDomain(){

    }


    std::string OxleyDomain::functionSpaceTypeAsString(int fsType) const 
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    bool OxleyDomain::operator==(const AbstractDomain& other) const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    int OxleyDomain::getTagFromSampleNo(int fsType, index_t sampleNo) const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    string OxleyDomain::showTagNames() const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    void OxleyDomain::interpolateOnDomain(escript::Data& target, const escript::Data& in) const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    bool OxleyDomain::probeInterpolationOnDomain(int fsType_source, int fsType_target) const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    signed char OxleyDomain::preferredInterpolationOnDomain(int fsType_source, int fsType_target) const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    bool OxleyDomain::commonFunctionSpace(const vector<int>& fs, int& resultcode) const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    escript::Data OxleyDomain::getX() const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    std::string OxleyDomain::getDescription() const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    escript::Data OxleyDomain::getNormal() const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    escript::Data OxleyDomain::getSize() const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    void OxleyDomain::setToX(escript::Data& arg) const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    void OxleyDomain::setToGradient(escript::Data& grad, const escript::Data& arg) const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    void OxleyDomain::setTags(const int fsType, const int newTag, const escript::Data& mask) const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    bool OxleyDomain::isCellOriented(int fsType) const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    int OxleyDomain::getNumberOfTagsInUse(int fsType) const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    const int* OxleyDomain::borrowListOfTagsInUse(int fsType) const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    bool OxleyDomain::canTag(int fsType) const
    {
        throw OxleyException("currently not implemented"); // ae: This is temporary
    }

    void OxleyDomain::writeToVTK(std::string filename) const
    {
        throw OxleyException("unknown error");
    }

} // end of namespace oxley

