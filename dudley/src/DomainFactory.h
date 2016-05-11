
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

namespace dudley {
  /**
     \brief
     A suite of factory methods for creating 2D and 3D dudley domains.
  */
 
  /**
     \brief
     recovers domain from a dump file
     \param fileName Input - The name of the file
  */
  escript::Domain_ptr loadMesh(const std::string& fileName);

  /**
     \brief
     Read a mesh from a file. For MPI parallel runs fan out the mesh to
     multiple processes.
     \param fileName Input -  The name of the file
     \param integrationOrder Input - ignored
     \param reducedIntegrationOrder Input - ignored
     \param optimize Input - switches on the optimization of node labels 
  */
  escript::Domain_ptr readMesh(const std::string& fileName,
                               int integrationOrder=-1,
                               int reducedIntegrationOrder=-1,
                               bool optimize=false);

  /**
     \brief
     Read a gmsh mesh file
     \param fileName Input -  The name of the file
     \param numDim Input -  spatial dimensionality
     \param integrationOrder Input - ignored
     \param reducedIntegrationOrder Input - ignored
     \param optimize Input - switches on the optimization of node labels 
  */
  escript::Domain_ptr readGmsh(const std::string& fileName,
                               int numDim, 
                               int integrationOrder = -1,
                               int reducedIntegrationOrder = -1, 
                               bool optimize = false);
                                     
                                     
   /**
   \brief Python driver for brick()
   \param args see brick() definition for order of params
   */
   escript::Domain_ptr brick_driver(const boost::python::list& args);

   /**
   \brief Python driver for rectangle()
   \param args see rectangle() definition for order of params
   */
   escript::Domain_ptr rectangle_driver(const boost::python::list& args);   
   
  /**
     \brief
     Creates a rectangular mesh with n0 x n1 x n2 elements over the brick 
     [0,l0] x [0,l1] x [0,l2].

     \param n0,n1,n2 Input - number of elements in each dimension
     \param order Input - ignored
     \param l0,l1,l2 Input - length of each side of brick
     \param integrationOrder Input - ignored
     \param reducedIntegrationOrder Input - ignored
     \param optimize
  */
  escript::Domain_ptr brick(escript::JMPI mpi_info, dim_t n0=1, dim_t n1=1,
                    dim_t n2=1, int order=1,
                    double l0=1.0, double l1=1.0, double l2=1.0,
                    int periodic0=0, int periodic1=0, int periodic2=0,
                    int integrationOrder=-1, int reducedIntegrationOrder=-1,
                    int useElementsOnFace=0, int useFullElementOrder=0,
                    bool optimize=false);

  /**
     \brief
     Creates a 2-dimensional rectangular mesh with n0 x n1 x 2 Tri3 elements
     over the rectangle [0,l0] x [0,l1]. The doubling of elements is due to
     splitting of rectangular elements.

     \param n0,n1 Input - number of elements in each dimension
     \param order Input - ignored
     \param l0,l1 Input - length of each side of rectangle
     \param integrationOrder Input - ignored
     \param reducedIntegrationOrder Input - ignored
     \param optimize - whether to optimize labelling
  */
  escript::Domain_ptr rectangle(escript::JMPI mpi_info, dim_t n0 = 1,
                                dim_t n1 = 1, int order = 1,
                                double l0 = 1.0, double l1 = 1.0,
                                int periodic0 = false, int periodic1 = false,
                                int integrationOrder = -1,
                                int reducedIntegrationOrder = -1,
                                int useElementsOnFace = 0,
                                int useFullElementOrder = 0,
                                bool optimize = false);

} // end of namespace

#endif // __DUDLEY_DOMAINFACTORY_H__

