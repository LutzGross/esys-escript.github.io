
/*****************************************************************************
*
* Copyright (c) 2003-2015 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


#if !defined  dudley_MeshAdapterFactory_20040526_H
#define dudley_MeshAdapterFactory_20040526_H
#include "system_dep.h"

#include "MeshAdapter.h"

#include "escript/AbstractContinuousDomain.h"

#include <boost/python/list.hpp>

#include <sstream>

namespace dudley {
  /**
     \brief
     A suite of factory methods for creating various MeshAdapters.

     Description:
     A suite of factory methods for creating various MeshAdapters.
  */
 
  /**
     \brief
     recovers mesg from a dump file
     \param fileName Input -  The name of the file.
  */
  DUDLEY_DLL_API
/*  escript::AbstractContinuousDomain* loadMesh(const std::string& fileName);*/
  escript::Domain_ptr loadMesh(const std::string& fileName);
  /**
     \brief
     Read a mesh from a file. For MPI parallel runs fan out the mesh to multiple processes.
     \param fileName Input -  The name of the file.
     \param integrationOrder Input - order of the quadrature scheme.  
     If integrationOrder<0 the integration order is selected independently.
     \param reducedIntegrationOrder Input - order of the reduced quadrature scheme.  
     If reducedIntegrationOrder<0 the integration order is selected independently.
     \param optimize Input - switches on the optimization of node labels 
  */
  DUDLEY_DLL_API
//   escript::AbstractContinuousDomain* readMesh(const std::string& fileName,
   escript::Domain_ptr readMesh(const std::string& fileName,
				     int integrationOrder=-1,
                                     int reducedIntegrationOrder=-1,
                                     int optimize=0);
  /**
     \brief
     Read a gmsh mesh file
     \param fileName Input -  The name of the file.
     \param numDim Input -  spatial dimension
     \param integrationOrder Input - order of the quadrature scheme.  
     If integrationOrder<0 the integration order is selected independently.
     \param reducedIntegrationOrder Input - order of the reduced quadrature scheme.  
     If reducedIntegrationOrder<0 the integration order is selected independently.
     \param optimize Input - switches on the optimization of node labels 
     \param useMacroElements
  */
  DUDLEY_DLL_API
//   escript::AbstractContinuousDomain* readGmsh(const std::string& fileName,
  escript::Domain_ptr readGmsh(const std::string& fileName,
				     int numDim, 
				     int integrationOrder=-1,
				     int reducedIntegrationOrder=-1, 
				     int optimize=0,
				     int useMacroElements=0);
				     
				     
   /**
   \brief Python driver for brick()
   \param args see brick() definition for order of params
   */
   DUDLEY_DLL_API
   escript::Domain_ptr brick_driver(const boost::python::list& args);

   /**
   \brief Python driver for rectangle()
   \param args see rectangle() definition for order of params
   */
   DUDLEY_DLL_API
   escript::Domain_ptr rectangle_driver(const boost::python::list& args);   
   
  /**
     \brief
     Creates a rectangular mesh with n0 x n1 x n2 elements over the brick 
     [0,l0] x [0,l1] x [0,l2].

     \param n0,n1,n2 Input - number of elements in each dimension
     \param order Input - =1, =-1 or =2 gives the order of shape function (-1= macro elements of order 1)
     \param l0,l1,l2 Input - length of each side of brick
     \param integrationOrder Input - order of the quadrature scheme.  
     If integrationOrder<0 the integration order is selected independently.
     \param reducedIntegrationOrder Input - order of the reduced quadrature scheme.  
     If reducedIntegrationOrder<0 the integration order is selected independently.
     \param useElementsOnFace Input - whether or not to use elements on face
     \param periodic0, periodic1, periodic2 Input - whether or not boundary 
     conditions of the dimension are periodic
     \param useFullElementOrder
     \param optimize
  */
  DUDLEY_DLL_API
  escript::Domain_ptr brick(esysUtils::JMPI& mpi_info, double n0=1,double n1=1,double n2=1,int order=1,
                    double l0=1.0,double l1=1.0,double l2=1.0,
                    int periodic0=0,int periodic1=0,
                    int periodic2=0,
                    int integrationOrder=-1,
                    int reducedIntegrationOrder=-1,
                    int useElementsOnFace=0,
                    int useFullElementOrder=0,
                    int optimize=0);

  /**
     \brief
     Creates a rectangular mesh with n0 x n1 elements over the brick 
     [0,l0] x [0,l1].

     \param n0,n1 Input - number of elements in each dimension [We only except floats for py transition]
     \param order Input - =1, =-1 or =2 gives the order of shape function (-1= macro elements of order 1)
     \param l0,l1 Input - length of each side of brick
     \param integrationOrder Input - order of the quadrature scheme. 
     If integrationOrder<0 the integration order is selected 
     independently.
     \param reducedIntegrationOrder Input - order of the reduced quadrature scheme.  
     If reducedIntegrationOrder<0 the integration order is selected independently.
     \param periodic0, periodic1 Input - whether or not the boundary
     conditions of the dimension are periodic
     \param useElementsOnFace Input - whether or not to use elements on face
     \param useFullElementOrder
     \param optimize
  */
  DUDLEY_DLL_API
  escript::Domain_ptr rectangle(esysUtils::JMPI& mpi_info, double n0=1,double n1=1,int order=1,
                                      double l0=1.0, double l1=1.0,
                                      int periodic0=false,int periodic1=false,
                                      int integrationOrder=-1,
                                      int reducedIntegrationOrder=-1,
                                      int useElementsOnFace=0,
                                      int useFullElementOrder=0,
                                      int optimize=0);

//  /**
//     \brief
//     Merges a list of meshes into one list.
//     \param meshList Input - The list of meshes.
//  */
//  DUDLEY_DLL_API
// //   escript::AbstractContinuousDomain* meshMerge(const boost::python::list& meshList);
//  escript::Domain_ptr meshMerge(const boost::python::list& meshList);

 
} // end of namespace
#endif
