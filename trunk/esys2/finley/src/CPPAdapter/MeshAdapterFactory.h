/* $Id$
/* 
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/
                                                                           
#if !defined  finley_MeshAdapterFactory_20040526_H
#define finley_MeshAdapterFactory_20040526_H

#include "finley/CPPAdapter/MeshAdapter.h"

#include "escript/Data/AbstractContinuousDomain.h"

#include <boost/python/list.hpp>

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
  */
  escript::AbstractContinuousDomain* readMesh(const std::string& fileName,
				     int integrationOrder=-1);
  /**
     \brief
     Creates a rectangular mesh with n1 x n2 x n3 elements over the brick 
     [0,l1] x [0,l2] x [0,l3].

     \param n1,n2,n3 Input - 
     \param order Input - =1 or =2 gives the order of shape function
     \param l1,l2,l3 Input - 
     \param integrationOrder Input - order of the quadrature scheme.
  */
  escript::AbstractContinuousDomain* brick(int n1=1,int n2=1,int n3=1,int order=1,
		    double l1=1.0,double l2=1.0,double l3=1.0,
		    int periodic0=0,int periodic1=0,
		    int periodic2=0,
		    int integrationOrder=-1,
		    int useElementsOnFace=0);
  /**
     \brdirief
     Creates a rectangular mesh with n1 x n2 elements over the brick 
     [0,l1] x [0,l2].

     \param n1,n2 Input - 
     \param order Input - =1 or =2 gives the order of shape function
     \param l1,l2 Input - 
     \param integrationOrder Input - order of the quadrature scheme. 
     If integrationOrder<0 the integration order is selected 
     independently.
  */
  escript::AbstractContinuousDomain* rectangle(int n1=1,int n2=1,int order=1,
				      double l1=1.0, double l2=1.0,
				      int periodic0=false,int periodic1=false,
				      int integrationOrder=-1,
				      int useElementsOnFace=false);
  /**
     \brief
     Creates an equidistant mesh with n elements over the interval [0,l].
     \param n1 Input -
     \param order Input - =1 or =2 gives the order of shape function.
     \param integrationOrder Input - order of the quadrature scheme. 
     If integrationOrder<0 the integration order is selected 
     independently.
  */
  escript::AbstractContinuousDomain* interval(int n1=1,int order=1,double l1=1.0,
				     int periodic0=false,
				     int integrationOrder=-1,
				     int useElementsOnFace=false);
  /**
     \brief
     Merges a list of meshes into one list.
     \param meshList Input - The list of meshes.
  */
  escript::AbstractContinuousDomain* meshMerge(const boost::python::list& meshList);
  /**
     \brief
     Detects matching faces in the mesh, removes them from the mesh 
     and joins the elements touched by the face elements.
     \param meshList Input - The list of meshes.
     \param safetyFacto Input - ??
     \param tolerance Input - ??
  */
  escript::AbstractContinuousDomain* glueFaces(const boost::python::list& meshList,
			   double safetyFactor=0.2, 
			   double tolerance=100.*std::numeric_limits<double>::epsilon());
  /**
     \brief
     Detects matching faces in the mesh and replaces them by joint elements.
     \param meshList Input - The list of meshes.
     \param safetyFacto Input - ??
     \param tolerance Input - ??
  */
  escript::AbstractContinuousDomain* joinFaces(const boost::python::list& meshList,
			double safety_factor=0.2, 
			double tolerance=100.*std::numeric_limits<double>::epsilon());
 
} // end of namespace
#endif
