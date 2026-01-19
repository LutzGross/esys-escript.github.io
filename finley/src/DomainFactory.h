
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __FINLEY_DOMAINFACTORY_H__
#define __FINLEY_DOMAINFACTORY_H__

#include <finley/FinleyDomain.h>

#include <boost/python/list.hpp>


#include <sstream>
namespace bp = boost::python;
/**
    \brief
    A suite of factory methods for creating various finley domains.
*/

namespace finley {

/**
    \brief Python driver for readMesh()
    \param fileName Path to mesh file
    \param integrationOrder Order of quadrature scheme (-1 for automatic)
    \param reducedIntegrationOrder Order of reduced quadrature scheme (-1 for automatic)
    \param optimize Enable optimization of node labels
    \param diracPoints Dirac point coordinates
    \param diracTags Dirac point tags
    \param comm MPI communicator (optional, defaults to MPI_COMM_WORLD if None)
*/
FINLEY_DLL_API
escript::Domain_ptr readMesh_driver(const std::string& fileName,
                                     int integrationOrder,
                                     int reducedIntegrationOrder,
                                     bool optimize,
                                     const boost::python::list& diracPoints,
                                     const boost::python::list& diracTags,
                                     const boost::python::object& comm);

/**
    \brief Python driver for readGMesh()
    \param fileName Path to gmsh file
    \param numDim Number of spatial dimensions
    \param integrationOrder Order of quadrature scheme (-1 for automatic)
    \param reducedIntegrationOrder Order of reduced quadrature scheme (-1 for automatic)
    \param optimize Enable optimization of node labels
    \param useMacroElements Enable usage of macro elements instead of second order elements
    \param diracPoints Dirac point coordinates
    \param diracTags Dirac point tags
    \param comm MPI communicator (optional, defaults to MPI_COMM_WORLD if None)
*/
FINLEY_DLL_API
escript::Domain_ptr readGmsh_driver(const std::string& fileName,
                                     int numDim,
                                     int integrationOrder,
                                     int reducedIntegrationOrder,
                                     bool optimize,
                                     bool useMacroElements,
                                     const boost::python::list& diracPoints,
                                     const boost::python::list& diracTags,
                                     const boost::python::object& comm);

/**
    \brief
    Creates a rectangular mesh with n0 x n1 x n2 elements over the brick 
    [0,l0] x [0,l1] x [0,l2].

    \param jmpi pointer to MPI world information structure
    \param n0,n1,n2 number of elements in each dimension
    \param order order of shape functions (1, 2, or -1 for macro
                 elements of order 1)
    \param l0,l1,l2 length of each side of brick
    \param periodic0,periodic1,periodic2 whether or not boundary 
           conditions of the dimension are periodic
    \param integrationOrder order of the quadrature scheme.
                            If <0 the order is selected automatically.
    \param reducedIntegrationOrder order of the reduced quadrature scheme.
                                   If <0 the order is selected automatically.
    \param useElementsOnFace whether or not to use elements on face
    \param useFullElementOrder whether to use second order elements
    \param optimize whether to apply optimization of node labels
    \param points dirac points to add
    \param tags
    \param tagNamesToNums
*/
FINLEY_DLL_API
escript::Domain_ptr brick(escript::JMPI jmpi,
                    dim_t n0=1, dim_t n1=1, dim_t n2=1, int order=1,
                    double l0=1.0, double l1=1.0, double l2=1.0,
                    bool periodic0=false, bool periodic1=false, bool periodic2=false,
                    int integrationOrder=-1, int reducedIntegrationOrder=-1,
                    bool useElementsOnFace=false,
                    bool useFullElementOrder=false, bool optimize=false,
                    const std::vector<double>& points=std::vector<double>(),
                    const std::vector<int>& tags=std::vector<int>(),
                    const std::map<std::string, int>& tagNamesToNums=std::map<std::string, int>());

/**
    \brief Python driver for brick()
    \param n_tuple tuple of (n0, n1, n2) - number of elements in each direction
    \param order order of shape functions (1, 2, or -1 for macro elements)
    \param l_tuple tuple of (l0, l1, l2) - length of each side
    \param periodic_tuple tuple of (periodic0, periodic1, periodic2) - periodic boundary conditions
    \param integrationOrder order of quadrature scheme (-1 for automatic)
    \param reducedIntegrationOrder order of reduced quadrature scheme (-1 for automatic)
    \param useElementsOnFace whether to use elements on face
    \param useFullElementOrder whether to use Hex27 elements
    \param optimize enable optimization of node labels
    \param diracPoints Dirac point coordinates
    \param diracTags Dirac point tags
    \param comm MPI communicator (optional, defaults to MPI_COMM_WORLD if None)
*/
FINLEY_DLL_API
escript::Domain_ptr brick_driver(const boost::python::tuple& n_tuple,
                                 int order,
                                 const boost::python::tuple& l_tuple,
                                 const boost::python::tuple& periodic_tuple,
                                 int integrationOrder, int reducedIntegrationOrder,
                                 bool useElementsOnFace, bool useFullElementOrder,
                                 bool optimize,
                                 const boost::python::list& diracPoints,
                                 const boost::python::list& diracTags,
                                 const boost::python::object& comm);

/**
    \brief
    Creates a 2-dimensional rectangular mesh with n0 x n1 elements over the
    rectangle [0,l0] x [0,l1].

    \param jmpi pointer to MPI world information structure
    \param n0,n1 number of elements in each dimension
    \param order order of shape functions (1, 2, or -1 for macro
                 elements of order 1)
    \param l0,l1 length of each side of rectangle
    \param periodic0,periodic1 whether or not the boundary conditions of the
                               dimension are periodic
    \param integrationOrder order of the quadrature scheme.
                            If <0 the order is selected automatically.
    \param reducedIntegrationOrder order of the reduced quadrature scheme.
                                   If <0 the order is selected automatically.
    \param useElementsOnFace whether or not to use elements on face
    \param useFullElementOrder
    \param optimize whether to optimize labelling
    \param points
    \param tags
    \param tagNamesToNums
*/
FINLEY_DLL_API
escript::Domain_ptr rectangle(escript::JMPI jmpi,
                              dim_t n0 = 1, dim_t n1 = 1, int order = 1,
                              double l0 = 1.0, double l1 = 1.0,
                              bool periodic0 = false, bool periodic1 = false,
                              int integrationOrder = -1,
                              int reducedIntegrationOrder = -1,
                              bool useElementsOnFace = false,
                              bool useFullElementOrder = false,
                              bool optimize = false,
                              const std::vector<double>& points = std::vector<double>(),
                              const std::vector<int>& tags = std::vector<int>(),
                              const std::map<std::string, int>& tagNamesToNums = std::map<std::string, int>());

/**
    \brief Python driver for rectangle()
    \param n_tuple tuple of (n0, n1) - number of elements in each direction
    \param order order of shape functions (1, 2, or -1 for macro elements)
    \param l_tuple tuple of (l0, l1) - length of each side
    \param periodic_tuple tuple of (periodic0, periodic1) - periodic boundary conditions
    \param integrationOrder order of quadrature scheme (-1 for automatic)
    \param reducedIntegrationOrder order of reduced quadrature scheme (-1 for automatic)
    \param useElementsOnFace whether to use elements on face
    \param useFullElementOrder whether to use Rec9 elements
    \param optimize enable optimization of node labels
    \param diracPoints Dirac point coordinates
    \param diracTags Dirac point tags
    \param comm MPI communicator (optional, defaults to MPI_COMM_WORLD if None)
*/
FINLEY_DLL_API
escript::Domain_ptr rectangle_driver(const boost::python::tuple& n_tuple,
                                     int order,
                                     const boost::python::tuple& l_tuple,
                                     const boost::python::tuple& periodic_tuple,
                                     int integrationOrder, int reducedIntegrationOrder,
                                     bool useElementsOnFace, bool useFullElementOrder,
                                     bool optimize,
                                     const boost::python::list& diracPoints,
                                     const boost::python::list& diracTags,
                                     const boost::python::object& comm);

/**
    \brief
    Merges a list of meshes into one list.
    \param meshList Input - The list of meshes.
*/
FINLEY_DLL_API
escript::Domain_ptr meshMerge(const boost::python::list& meshList);

/**
    \brief
    Detects matching faces in the mesh, removes them from the mesh 
    and joins the elements touched by the face elements.
    \param meshList The list of meshes.
    \param safetyFactor
    \param tolerance
    \param optimize switches on the optimization of node labels 
*/
FINLEY_DLL_API
escript::Domain_ptr glueFaces(const boost::python::list& meshList,
                              double safetyFactor = 0.2, double tolerance = 1.e-8,
                              bool optimize = false);

/**
    \brief
    Detects matching faces in the mesh and replaces them by joint elements.
    \param meshList The list of meshes
    \param safetyFactor
    \param tolerance
    \param optimize switches on the optimization of node labels 
*/
FINLEY_DLL_API
escript::Domain_ptr joinFaces(const boost::python::list& meshList,
                              double safetyFactor = 0.2, double tolerance = 1.e-8,
                              bool optimize = false);

 
} // end of namespace

#endif // __FINLEY_DOMAINFACTORY_H__

