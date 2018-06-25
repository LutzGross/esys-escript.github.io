
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

#ifndef __DUDLEY_DOMAINFACTORY_H__
#define __DUDLEY_DOMAINFACTORY_H__

#include <dudley/DudleyDomain.h>

#include <boost/python/list.hpp>

#include <sstream>

/**
    \brief
    A suite of factory methods for creating 2D and 3D dudley domains.
*/

namespace dudley {

/**
    \brief
    reads a mesh from a fly file. For MPI parallel runs fans out the mesh to
    multiple processes.
    \param fileName the name of the file
    \param integrationOrder ignored
    \param reducedIntegrationOrder ignored
    \param optimize whether to optimize the node labels
*/
escript::Domain_ptr readMesh(const std::string& fileName,
                             int integrationOrder = -1,
                             int reducedIntegrationOrder = -1,
                             bool optimize = false);

/**
    \brief
    reads a gmsh mesh file
    \param fileName the name of the file
    \param numDim spatial dimensionality
    \param integrationOrder ignored
    \param reducedIntegrationOrder ignored
    \param optimize whether to optimize the node labels 
*/
escript::Domain_ptr readGmsh(const std::string& fileName, int numDim,
                             int integrationOrder = -1,
                             int reducedIntegrationOrder = -1,
                             bool optimize = false);

/**
    \brief
    Creates a rectangular mesh with n0 x n1 x n2 elements over the brick 
    [0,l0] x [0,l1] x [0,l2].

    \param jmpi pointer to MPI world information structure
    \param n0,n1,n2 number of elements in each dimension
    \param order ignored
    \param l0,l1,l2 length of each side of brick
    \param integrationOrder ignored
    \param reducedIntegrationOrder ignored
    \param optimize
*/
escript::Domain_ptr brick(escript::JMPI jmpi,
                    dim_t n0=1, dim_t n1=1, dim_t n2=1, int order=1,
                    double l0=1.0, double l1=1.0, double l2=1.0,
                    bool periodic0=false, bool periodic1=false, bool periodic2=false,
                    int integrationOrder=-1, int reducedIntegrationOrder=-1,
                    bool useElementsOnFace=false, bool useFullElementOrder=false,
                    bool optimize=false);

/**
    \brief Python driver for brick()
    \param args see brick() definition for order of params
*/
escript::Domain_ptr brick_driver(const boost::python::list& args);

/**
    \brief
    Creates a 2-dimensional rectangular mesh with n0 x n1 x 2 Tri3 elements
    over the rectangle [0,l0] x [0,l1]. The doubling of elements is due to
    splitting of rectangular elements.

    \param jmpi pointer to MPI world information structure
    \param n0,n1 number of elements in each dimension
    \param order ignored
    \param l0,l1 length of each side of rectangle
    \param periodic0,periodic1 ignored
    \param integrationOrder ignored
    \param reducedIntegrationOrder ignored
    \param useElementsOnFace ignored
    \param useFullElementOrder ignored
    \param optimize whether to optimize labelling
*/
escript::Domain_ptr rectangle(escript::JMPI jmpi,
                              dim_t n0 = 1, dim_t n1 = 1, int order = 1,
                              double l0 = 1.0, double l1 = 1.0,
                              bool periodic0 = false, bool periodic1 = false,
                              int integrationOrder = -1,
                              int reducedIntegrationOrder = -1,
                              bool useElementsOnFace = false,
                              bool useFullElementOrder = false,
                              bool optimize = false);

/**
    \brief Python driver for rectangle()
    \param args see rectangle() definition for order of params
*/
escript::Domain_ptr rectangle_driver(const boost::python::list& args);


} // end of namespace

#endif // __DUDLEY_DOMAINFACTORY_H__

