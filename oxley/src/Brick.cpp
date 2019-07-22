
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

#include <oxley/Brick.h>

using namespace boost::python;

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
       writes information about the mesh to standard output
       \param full whether to print additional data
    */
void Brick::Print_Mesh_Info(const bool full) const {


    }



} // end of namespace oxley

