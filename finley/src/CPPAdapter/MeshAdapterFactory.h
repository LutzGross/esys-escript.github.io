
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


#if !defined  finley_MeshAdapterFactory_20040526_H
#define finley_MeshAdapterFactory_20040526_H
#include "system_dep.h"

#include "finley/Finley.h"
#include "finley/Mesh.h"
#include "finley/RectangularMesh.h"

#include "MeshAdapter.h"

#include "escript/AbstractContinuousDomain.h"

#include <boost/python/list.hpp>

#include <sstream>

namespace finley {
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
  FINLEY_DLL_API
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
  FINLEY_DLL_API
//   escript::AbstractContinuousDomain* readMesh(const std::string& fileName,
   escript::Domain_ptr readMesh(const std::string& fileName,
                                     int integrationOrder=-1,
                                     int reducedIntegrationOrder=-1,
                                     bool optimize=false);
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
  FINLEY_DLL_API
//   escript::AbstractContinuousDomain* readGmsh(const std::string& fileName,
  escript::Domain_ptr readGmsh(const std::string& fileName,
                               int numDim, 
                               int integrationOrder=-1,
                               int reducedIntegrationOrder=-1, 
                               bool optimize=false,
                               bool useMacroElements=false);

  /**
     \brief
     Creates a rectangular mesh with n0 x n1 x n2 elements over the brick 
     [0,l0] x [0,l1] x [0,l2].

     \param n0,n1,n2 number of elements in each dimension
     \param order =1, =-1 or =2 gives the order of shape function
                  (-1= macro elements of order 1)
     \param l0,l1,l2 length of each side of brick
     \param periodic0, periodic1, periodic2 whether or not boundary 
            conditions of the dimension are periodic
     \param integrationOrder order of the quadrature scheme.  
          If integrationOrder<0 the integration order is selected independently.
     \param reducedIntegrationOrder order of the reduced quadrature scheme.  
          If reducedIntegrationOrder<0 the integration order is selected independently.
     \param useElementsOnFace whether or not to use elements on face
     \param useFullElementOrder whether to use second order elements
     \param optimize whether to apply optimization
     \param points
     \param tags
     \param tagnamestonums
  */
FINLEY_DLL_API
escript::Domain_ptr brick(int n0=1, int n1=1, int n2=1, int order=1,
                    double l0=1.0, double l1=1.0, double l2=1.0,
                    bool periodic0=false, bool periodic1=false,
                    bool periodic2=false,
                    int integrationOrder=-1,
                    int reducedIntegrationOrder=-1,
                    bool useElementsOnFace=false,
                    bool useFullElementOrder=false,
                    bool optimize=false,
                    const std::vector<double>& points=std::vector<double>(),
                    const std::vector<int>& tags=std::vector<int>(),
                    const std::map<std::string, int>& tagnamestonums=std::map<std::string, int>()
                    );
                    
   /**
   \brief Python driver for brick()
   \param args see brick() definition for order of params
   */
   FINLEY_DLL_API
   escript::Domain_ptr brick_driver(const boost::python::list& args);

   /**
   \brief Python driver for rectangle()
   \param args see rectangle() definition for order of params
   */
   FINLEY_DLL_API
   escript::Domain_ptr rectangle_driver(const boost::python::list& args);   
   
   
  /**
     \brief
     Creates a rectangular mesh with n0 x n1 elements over the brick 
     [0,l0] x [0,l1].

     \param n0,n1 number of elements in each dimension
     \param l0,l1 length of each side of brick
     \param order =1, =-1 or =2 gives the order of shape function
                  (-1= macro elements of order 1)
     \param periodic0, periodic1 whether or not the boundary conditions of the
            dimension are periodic
     \param integrationOrder order of the quadrature scheme. 
            If integrationOrder<0 the integration order is selected 
            independently.
     \param reducedIntegrationOrder order of the reduced quadrature scheme.  
            If reducedIntegrationOrder<0 the integration order is selected
            independently.
     \param useElementsOnFace whether or not to use elements on face
     \param useFullElementOrder
     \param optimize
     \param points
     \param tags
     \param tagnamestonums
  */
FINLEY_DLL_API
  escript::Domain_ptr rectangle(int n0=1, int n1=1, int order=1,
                                double l0=1.0, double l1=1.0,
                                bool periodic0=false, bool periodic1=false,
                                int integrationOrder=-1,
                                int reducedIntegrationOrder=-1, 
                                bool useElementsOnFace=false,
                                bool useFullElementOrder=false,
                                bool optimize=false,
                                const std::vector<double>& points=std::vector<double>(),
                                const std::vector<int>& tags=std::vector<int>(),
                                const std::map<std::string, int>& tagnamestonums=std::map<std::string, int>()
                        );
  /**
     \brief
     Merges a list of meshes into one list.
     \param meshList Input - The list of meshes.
  */
  FINLEY_DLL_API
//   escript::AbstractContinuousDomain* meshMerge(const boost::python::list& meshList);
  escript::Domain_ptr meshMerge(const boost::python::list& meshList);
  /**
     \brief
     Detects matching faces in the mesh, removes them from the mesh 
     and joins the elements touched by the face elements.
     \param meshList Input - The list of meshes.
     \param safetyFactor Input - ??
     \param tolerance Input - ??
     \param optimize Input - switches on the optimization of node labels 
  */
  FINLEY_DLL_API
//   escript::AbstractContinuousDomain* glueFaces(const boost::python::list& meshList,
  escript::Domain_ptr glueFaces(const boost::python::list& meshList,
                           double safetyFactor=0.2, 
                           double tolerance=1.e-8,
                           bool optimize=false);
  /**
     \brief
     Detects matching faces in the mesh and replaces them by joint elements.
     \param meshList Input - The list of meshes.
     \param safetyFactor Input - ??
     \param tolerance Input - ??
     \param optimize Input - switches on the optimization of node labels 
  */
  FINLEY_DLL_API
//   escript::AbstractContinuousDomain* joinFaces(const boost::python::list& meshList,
  escript::Domain_ptr joinFaces(const boost::python::list& meshList,
                        double safetyFactor=0.2, 
                        double tolerance=1.e-8,
                        bool optimize=false);
 
} // end of namespace
#endif
