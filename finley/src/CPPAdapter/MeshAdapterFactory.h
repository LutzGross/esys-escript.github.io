/* $Id$ */
/* 
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

#if !defined  finley_MeshAdapterFactory_20040526_H
#define finley_MeshAdapterFactory_20040526_H
#include "system_dep.h"

extern "C" {
#include "../Finley.h"
#include "../Mesh.h"
#include "../RectangularMesh.h"
}

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
     Read a mesh from a file
     \param fileName Input -  The name of the file.
     \param integrationOrder Input - order of the quadrature scheme.  
     If integrationOrder<0 the integration order is selected independently.
  */
  FINLEY_DLL_API
  escript::AbstractContinuousDomain* readMesh(const std::string& fileName,
				     int integrationOrder=-1);
  /**
     \brief
     Creates a rectangular mesh with n0 x n1 x n2 elements over the brick 
     [0,l0] x [0,l1] x [0,l2].

     \param n0,n1,n2 Input - number of elements in each dimension
     \param order Input - =1 or =2 gives the order of shape function
     \param l0,l1,l2 Input - length of each side of brick
     \param integrationOrder Input - order of the quadrature scheme.
     \param useElementsOnFace Input - whether or not to use elements on face
     \param periodic0, periodic1, periodic2 Input - whether or not boundary 
     conditions of the dimension are periodic
  */
  FINLEY_DLL_API
  escript::AbstractContinuousDomain* brick(int n0=1,int n1=1,int n2=1,int order=1,
		    double l0=1.0,double l1=1.0,double l2=1.0,
		    int periodic0=0,int periodic1=0,
		    int periodic2=0,
		    int integrationOrder=-1,
		    int useElementsOnFace=0);
  /**
     \brief
     Creates a rectangular mesh with n0 x n1 elements over the brick 
     [0,l0] x [0,l1].

     \param n0,n1 Input - number of elements in each dimension
     \param order Input - =1 or =2 gives the order of shape function
     \param l0,l1 Input - length of each side of brick
     \param integrationOrder Input - order of the quadrature scheme. 
     If integrationOrder<0 the integration order is selected 
     independently.
     \param periodic0, periodic1 Input - whether or not the boundary
     conditions of the dimension are periodic
     \param useElementsOnFace Input - whether or not to use elements on face
  */
  FINLEY_DLL_API
  escript::AbstractContinuousDomain* rectangle(int n0=1,int n1=1,int order=1,
				      double l0=1.0, double l1=1.0,
				      int periodic0=false,int periodic1=false,
				      int integrationOrder=-1,
				      int useElementsOnFace=false);
  /**
     \brief
     Creates an equidistant mesh with n elements over the interval [0,l].
     \param n0 Input - number of elements
     \param order Input - =1 or =2 gives the order of shape function.
     \param l0 Input - length of the brick
     \param integrationOrder Input - order of the quadrature scheme. 
     If integrationOrder<0 the integration order is selected 
     independently.
     \param periodic0 Input - whether or not the boundary conditions are
     periodic
     \param useElementsOnFace Input - whether or not to use the elements
     on the face
  */
  FINLEY_DLL_API
  escript::AbstractContinuousDomain* interval(int n0=1,int order=1,double l0=1.0,
				     int periodic0=false,
				     int integrationOrder=-1,
				     int useElementsOnFace=false);
  /**
     \brief
     Merges a list of meshes into one list.
     \param meshList Input - The list of meshes.
  */
  FINLEY_DLL_API
  escript::AbstractContinuousDomain* meshMerge(const boost::python::list& meshList);
  /**
     \brief
     Detects matching faces in the mesh, removes them from the mesh 
     and joins the elements touched by the face elements.
     \param meshList Input - The list of meshes.
     \param safetyFactor Input - ??
     \param tolerance Input - ??
  */
  FINLEY_DLL_API
  escript::AbstractContinuousDomain* glueFaces(const boost::python::list& meshList,
			   double safetyFactor=0.2, 
			   double tolerance=1.e-8);
  /**
     \brief
     Detects matching faces in the mesh and replaces them by joint elements.
     \param meshList Input - The list of meshes.
     \param safetyFactor Input - ??
     \param tolerance Input - ??
  */
  FINLEY_DLL_API
  escript::AbstractContinuousDomain* joinFaces(const boost::python::list& meshList,
			double safetyFactor=0.2, 
			double tolerance=1.e-8);
 
} // end of namespace
#endif
