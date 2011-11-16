
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#ifndef __RIPLEY_DOMAINFACTORY_H__
#define __RIPLEY_DOMAINFACTORY_H__

#include <ripley/RipleyDomain.h>
#include <escript/AbstractContinuousDomain.h>

namespace ripley {
  /**
     \brief
     A suite of factory methods for creating ripley meshes.

     Description:
     A suite of factory methods for creating ripley meshes.
  */
 
  /**
     \brief recovers mesh from a file created with the dump() method.
     \param fileName the name of the file
  */
  RIPLEY_DLL_API
  escript::Domain_ptr loadMesh(const std::string& fileName);

  /**
     \brief reads a mesh from a file. For MPI parallel runs fans out the
            mesh to multiple processes.
     \param fileName the name of the file
     \param optimize whether to optimize node labels 
  */
  RIPLEY_DLL_API
  escript::Domain_ptr readMesh(const std::string& fileName, bool optimize=false);

  /**
     \brief reads a mesh from a gmsh file
     \param fileName the name of the file
     \param numDim spatial dimension
     \param optimize whether to optimize node labels 
  */
  RIPLEY_DLL_API
  escript::Domain_ptr readGmsh(const std::string& fileName, int numDim, bool optimize=false);

  /**
     \brief creates a rectangular mesh with n0 x n1 x n2 elements over the
            brick [0,l0] x [0,l1] x [0,l2].
     \param n0,n1,n2 number of elements in each dimension
     \param l0,l1,l2 length of each side of brick
     \param optimize whether to optimize node labels
  */
  RIPLEY_DLL_API
  escript::Domain_ptr brick(int n0=1, int n1=1, int n2=1, double l0=1.0,
                            double l1=1.0, double l2=1.0, bool optimize=false);

  /**
     \brief creates a rectangular mesh with n0 x n1 elements over the rectangle
            [0,l0] x [0,l1].
     \param n0,n1 number of elements in each dimension
     \param l0,l1 length of each side of rectangle
     \param optimize whether to optimize node labels
  */
  RIPLEY_DLL_API
  escript::Domain_ptr rectangle(int n0=1, int n1=1, double l0=1.0,
                                double l1=1.0, bool optimize=false);

 
} // end of namespace ripley

#endif // __RIPLEY_DOMAINFACTORY_H__

