
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

#include <oxley/OxleyDomain.h>
#include <p4est/p4est.h>

namespace oxley {

/**
   \brief
   Rectangle is the 2-dimensional implementation of a SpeckleyDomain.
*/
class Rectangle: public OxleyDomain
{
public:

    /**
       \brief creates a rectangular mesh with n0 x n1 elements over the
              rectangle [x0,x1] x [y0,y1].
       \param 
    */
    Rectangle();

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
       writes information about the mesh to standard output
       \param full whether to print additional data
    */
    void Print_Mesh_Info(const bool full=false) const;


};

////////////////////////////// inline methods ////////////////////////////////

} //end namespace


