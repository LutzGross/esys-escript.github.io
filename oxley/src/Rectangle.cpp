
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

#include <oxley/Rectangle.h>

namespace oxley {


Rectangle::Rectangle(){



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
       writes information about the mesh to standard output
       \param full whether to print additional data
    */
void Rectangle::Print_Mesh_Info(const bool full) const {


    }



} // end of namespace oxley

