
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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

#include "MeshAdapterFactory.h"
#include "esysUtils/blocktimer.h"
#ifdef ESYS_MPI
#include "esysUtils/Esys_MPI.h"
#endif
#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

#include <boost/python/extract.hpp>
#include <boost/scoped_array.hpp>

#include <sstream>


using namespace std;
using namespace escript;

namespace finley {

#ifdef USE_NETCDF
  // A convenience method to retrieve an integer attribute from a NetCDF file
  int NetCDF_Get_Int_Attribute(NcFile *dataFile, const std::string &fName, char *attr_name)
  {
    NcAtt *attr;
    char error_msg[LenErrorMsg_MAX];
    if (! (attr=dataFile->get_att(attr_name)) ) {
      sprintf(error_msg,"loadMesh: Error retrieving integer attribute '%s' from NetCDF file '%s'", attr_name, fName.c_str());
      throw DataException(error_msg);
    }
    int temp = attr->as_int(0);
    delete attr;
    return(temp);
  }
#endif

  inline void cleanupAndThrow(Mesh* mesh, string msg)
  {
      delete mesh;
      string msgPrefix("loadMesh: NetCDF operation failed - ");
      throw DataException(msgPrefix+msg);
  }

//   AbstractContinuousDomain* loadMesh(const std::string& fileName)
  Domain_ptr loadMesh(const std::string& fileName)
  {
#ifdef USE_NETCDF
    esysUtils::JMPI mpi_info = esysUtils::makeInfo( MPI_COMM_WORLD );
    Mesh *mesh_p=NULL;
    char error_msg[LenErrorMsg_MAX];

    const std::string fName(esysUtils::appendRankToFileName(fileName,
                                              mpi_info->size, mpi_info->rank));

    double blocktimer_start = blocktimer_time();
    resetError();

    // Open NetCDF file for reading
    NcAtt *attr;
    NcVar *nc_var_temp;
    // netCDF error handler
    NcError err(NcError::silent_nonfatal);
    // Create the NetCDF file.
    NcFile dataFile(fName.c_str(), NcFile::ReadOnly);
    if (!dataFile.is_valid()) {
      sprintf(error_msg,"loadMesh: Opening NetCDF file '%s' for reading failed.", fName.c_str());
      setError(IO_ERROR,error_msg);
      throw DataException(error_msg);
    }

    // Read NetCDF integer attributes
    int mpi_size                        = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"mpi_size");
    int mpi_rank                        = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"mpi_rank");
    int numDim                          = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"numDim");
    int order                           = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"order");
    int reduced_order                   = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"reduced_order");
    int numNodes                        = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"numNodes");
    int num_Elements                    = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_Elements");
    int num_FaceElements                = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_FaceElements");
    int num_ContactElements             = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_ContactElements");
    int num_Points                      = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_Points");
    int num_Elements_numNodes           = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_Elements_numNodes");
    int Elements_TypeId                 = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"Elements_TypeId");
    int num_FaceElements_numNodes       = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_FaceElements_numNodes");
    int FaceElements_TypeId             = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"FaceElements_TypeId");
    int num_ContactElements_numNodes    = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_ContactElements_numNodes");
    int ContactElements_TypeId          = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"ContactElements_TypeId");
    int Points_TypeId                   = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"Points_TypeId");
    int num_Tags                        = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_Tags");

    // Verify size and rank
    if (mpi_info->size != mpi_size) {
      sprintf(error_msg, "loadMesh: The NetCDF file '%s' can only be read on %d CPUs instead of %d", fName.c_str(), mpi_size, mpi_info->size);
      throw DataException(error_msg);
    }
    if (mpi_info->rank != mpi_rank) {
      sprintf(error_msg, "loadMesh: The NetCDF file '%s' should be read on CPU #%d instead of %d", fName.c_str(), mpi_rank, mpi_info->rank);
      throw DataException(error_msg);
    }

    // Read mesh name
    if (! (attr=dataFile.get_att("Name")) ) {
      sprintf(error_msg,"loadMesh: Error retrieving mesh name from NetCDF file '%s'", fName.c_str());
      throw DataException(error_msg);
    }
    boost::scoped_array<char> name(attr->as_string(0));
    delete attr;

    /* allocate mesh */
    mesh_p = new Mesh(name.get(), numDim, mpi_info);
    if (noError()) {

        /* read nodes */
        mesh_p->Nodes->allocTable(numNodes);
        // Nodes_Id
        if (! ( nc_var_temp = dataFile.get_var("Nodes_Id")) )
            cleanupAndThrow(mesh_p, "get_var(Nodes_Id)");
        if (! nc_var_temp->get(&mesh_p->Nodes->Id[0], numNodes) )
            cleanupAndThrow(mesh_p, "get(Nodes_Id)");
        // Nodes_Tag
        if (! ( nc_var_temp = dataFile.get_var("Nodes_Tag")) )
            cleanupAndThrow(mesh_p, "get_var(Nodes_Tag)");
        if (! nc_var_temp->get(&mesh_p->Nodes->Tag[0], numNodes) )
            cleanupAndThrow(mesh_p, "get(Nodes_Tag)");
        // Nodes_gDOF
        if (! ( nc_var_temp = dataFile.get_var("Nodes_gDOF")) )
            cleanupAndThrow(mesh_p, "get_var(Nodes_gDOF)");
        if (! nc_var_temp->get(&mesh_p->Nodes->globalDegreesOfFreedom[0], numNodes) )
            cleanupAndThrow(mesh_p, "get(Nodes_gDOF)");
        // Nodes_gNI
        if (! ( nc_var_temp = dataFile.get_var("Nodes_gNI")) )
            cleanupAndThrow(mesh_p, "get_var(Nodes_gNI)");
        if (! nc_var_temp->get(&mesh_p->Nodes->globalNodesIndex[0], numNodes) )
            cleanupAndThrow(mesh_p, "get(Nodes_gNI)");
        // Nodes_grDfI
        if (! ( nc_var_temp = dataFile.get_var("Nodes_grDfI")) )
            cleanupAndThrow(mesh_p, "get_var(Nodes_grDfI)");
        if (! nc_var_temp->get(&mesh_p->Nodes->globalReducedDOFIndex[0], numNodes) )
            cleanupAndThrow(mesh_p, "get(Nodes_grDfI)");
        // Nodes_grNI
        if (! ( nc_var_temp = dataFile.get_var("Nodes_grNI")) )
            cleanupAndThrow(mesh_p, "get_var(Nodes_grNI)");
        if (! nc_var_temp->get(&mesh_p->Nodes->globalReducedNodesIndex[0], numNodes) )
            cleanupAndThrow(mesh_p, "get(Nodes_grNI)");
        // Nodes_Coordinates
        if (!(nc_var_temp = dataFile.get_var("Nodes_Coordinates")))
            cleanupAndThrow(mesh_p, "get_var(Nodes_Coordinates)");
        if (! nc_var_temp->get(&(mesh_p->Nodes->Coordinates[0]), numNodes, numDim) )
            cleanupAndThrow(mesh_p, "get(Nodes_Coordinates)");
        mesh_p->Nodes->updateTagList();

        /* read elements */
        if (noError()) {
            const_ReferenceElementSet_ptr refElements(new ReferenceElementSet(
                        (ElementTypeId)Elements_TypeId, order, reduced_order));
            if (noError())  {
                mesh_p->Elements=new ElementFile(refElements, mpi_info);
            }
            if (noError())
                mesh_p->Elements->allocTable(num_Elements);
            if (noError()) {
                mesh_p->Elements->minColor=0;
                mesh_p->Elements->maxColor=num_Elements-1;
                if (num_Elements>0) {
                   // Elements_Id
                   if (! ( nc_var_temp = dataFile.get_var("Elements_Id")) )
                       cleanupAndThrow(mesh_p, "get_var(Elements_Id)");
                   if (! nc_var_temp->get(&mesh_p->Elements->Id[0], num_Elements) )
                       cleanupAndThrow(mesh_p, "get(Elements_Id)");
                   // Elements_Tag
                   if (! ( nc_var_temp = dataFile.get_var("Elements_Tag")) )
                       cleanupAndThrow(mesh_p, "get_var(Elements_Tag)");
                   if (! nc_var_temp->get(&mesh_p->Elements->Tag[0], num_Elements) )
                       cleanupAndThrow(mesh_p, "get(Elements_Tag)");
                   // Elements_Owner
                   if (! ( nc_var_temp = dataFile.get_var("Elements_Owner")) )
                       cleanupAndThrow(mesh_p, "get_var(Elements_Owner)");
                   if (! nc_var_temp->get(&mesh_p->Elements->Owner[0], num_Elements) )
                       cleanupAndThrow(mesh_p, "get(Elements_Owner)");
                   // Elements_Color
                   if (! ( nc_var_temp = dataFile.get_var("Elements_Color")) )
                       cleanupAndThrow(mesh_p, "get_var(Elements_Color)");
                   if (! nc_var_temp->get(&mesh_p->Elements->Color[0], num_Elements) )
                       cleanupAndThrow(mesh_p, "get(Elements_Color)");
                   // Now we need to adjust maxColor
                   index_t mc=mesh_p->Elements->Color[0];
                   for (index_t i=1;i<num_Elements;++i) {
                       if (mc<mesh_p->Elements->Color[i]) {
                           mc = mesh_p->Elements->Color[i];
                       }
                   }
                   mesh_p->Elements->maxColor=mc;
                   // Elements_Nodes
                   int *Elements_Nodes = TMPMEMALLOC(num_Elements*num_Elements_numNodes,int);
                   if (!(nc_var_temp = dataFile.get_var("Elements_Nodes"))) {
                       TMPMEMFREE(Elements_Nodes);
                       cleanupAndThrow(mesh_p, "get_var(Elements_Nodes)");
                   }
                   if (! nc_var_temp->get(&(Elements_Nodes[0]), num_Elements, num_Elements_numNodes) ) {
                       TMPMEMFREE(Elements_Nodes);
                       cleanupAndThrow(mesh_p, "get(Elements_Nodes)");
                   }

                   // Copy temp array into mesh_p->Elements->Nodes
                   for (int i=0; i<num_Elements; i++) {
                       for (int j=0; j<num_Elements_numNodes; j++) {
                           mesh_p->Elements->Nodes[INDEX2(j,i,num_Elements_numNodes)]
                                = Elements_Nodes[INDEX2(j,i,num_Elements_numNodes)];
                       }
                   }
                   TMPMEMFREE(Elements_Nodes);
                } /* num_Elements>0 */
                mesh_p->Elements->updateTagList();
            }
        }

        /* get the face elements */
        if (noError()) {
            const_ReferenceElementSet_ptr refFaceElements(
                    new ReferenceElementSet((ElementTypeId)FaceElements_TypeId,
                        order, reduced_order));
            if (noError())  {
                mesh_p->FaceElements=new ElementFile(refFaceElements, mpi_info);
            }
            if (noError())
                mesh_p->FaceElements->allocTable(num_FaceElements);
            if (noError()) {
                mesh_p->FaceElements->minColor=0;
                mesh_p->FaceElements->maxColor=num_FaceElements-1;
                if (num_FaceElements>0) {
                   // FaceElements_Id
                   if (! ( nc_var_temp = dataFile.get_var("FaceElements_Id")) )
                       cleanupAndThrow(mesh_p, "get_var(FaceElements_Id)");
                   if (! nc_var_temp->get(&mesh_p->FaceElements->Id[0], num_FaceElements) )
                       cleanupAndThrow(mesh_p, "get(FaceElements_Id)");
                   // FaceElements_Tag
                   if (! ( nc_var_temp = dataFile.get_var("FaceElements_Tag")) )
                       cleanupAndThrow(mesh_p, "get_var(FaceElements_Tag)");
                   if (! nc_var_temp->get(&mesh_p->FaceElements->Tag[0], num_FaceElements) )
                       cleanupAndThrow(mesh_p, "get(FaceElements_Tag)");
                   // FaceElements_Owner
                   if (! ( nc_var_temp = dataFile.get_var("FaceElements_Owner")) )
                       cleanupAndThrow(mesh_p, "get_var(FaceElements_Owner)");
                   if (! nc_var_temp->get(&mesh_p->FaceElements->Owner[0], num_FaceElements) )
                       cleanupAndThrow(mesh_p, "get(FaceElements_Owner)");
                   // FaceElements_Color
                   if (! ( nc_var_temp = dataFile.get_var("FaceElements_Color")) )
                       cleanupAndThrow(mesh_p, "get_var(FaceElements_Color)");
                   if (! nc_var_temp->get(&mesh_p->FaceElements->Color[0], num_FaceElements) )
                       cleanupAndThrow(mesh_p, "get(FaceElements_Color)");
                   // Now we need to adjust maxColor
                   index_t mc=mesh_p->FaceElements->Color[0];
                   for (index_t i=1;i<num_FaceElements;++i) {
                       if (mc<mesh_p->FaceElements->Color[i]) {
                           mc = mesh_p->FaceElements->Color[i];
                       }
                   }
                   mesh_p->FaceElements->maxColor=mc;
                   // FaceElements_Nodes
                   int *FaceElements_Nodes = TMPMEMALLOC(num_FaceElements*num_FaceElements_numNodes,int);
                   if (!(nc_var_temp = dataFile.get_var("FaceElements_Nodes"))) {
                       TMPMEMFREE(FaceElements_Nodes);
                       cleanupAndThrow(mesh_p, "get_var(FaceElements_Nodes)");
                   }
                   if (! nc_var_temp->get(&(FaceElements_Nodes[0]), num_FaceElements, num_FaceElements_numNodes) ) {
                       TMPMEMFREE(FaceElements_Nodes);
                       cleanupAndThrow(mesh_p, "get(FaceElements_Nodes)");
                   }
                   // Copy temp array into mesh_p->FaceElements->Nodes
                   for (int i=0; i<num_FaceElements; i++) {
                       for (int j=0; j<num_FaceElements_numNodes; j++) {
                           mesh_p->FaceElements->Nodes[INDEX2(j,i,num_FaceElements_numNodes)] = FaceElements_Nodes[INDEX2(j,i,num_FaceElements_numNodes)];
                       }
                   }
                   TMPMEMFREE(FaceElements_Nodes);
                } /* num_FaceElements>0 */
                mesh_p->FaceElements->updateTagList();
            }
        }

        /* get the Contact elements */
        if (noError()) {
            const_ReferenceElementSet_ptr refContactElements(
                 new ReferenceElementSet((ElementTypeId)ContactElements_TypeId,
                     order, reduced_order));
            if (noError()) {
                mesh_p->ContactElements=new ElementFile(refContactElements, mpi_info);
            }
            if (noError())
                mesh_p->ContactElements->allocTable(num_ContactElements);
            if (noError()) {
                mesh_p->ContactElements->minColor=0;
                mesh_p->ContactElements->maxColor=num_ContactElements-1;
                if (num_ContactElements>0) {
                   // ContactElements_Id
                   if (! ( nc_var_temp = dataFile.get_var("ContactElements_Id")) )
                       cleanupAndThrow(mesh_p, "get_var(ContactElements_Id)");
                   if (! nc_var_temp->get(&mesh_p->ContactElements->Id[0], num_ContactElements) )
                       cleanupAndThrow(mesh_p, "get(ContactElements_Id)");
                   // ContactElements_Tag
                   if (! ( nc_var_temp = dataFile.get_var("ContactElements_Tag")) )
                       cleanupAndThrow(mesh_p, "get_var(ContactElements_Tag)");
                   if (! nc_var_temp->get(&mesh_p->ContactElements->Tag[0], num_ContactElements) )
                       cleanupAndThrow(mesh_p, "get(ContactElements_Tag)");
                   // ContactElements_Owner
                   if (! ( nc_var_temp = dataFile.get_var("ContactElements_Owner")) )
                       cleanupAndThrow(mesh_p, "get_var(ContactElements_Owner)");
                   if (! nc_var_temp->get(&mesh_p->ContactElements->Owner[0], num_ContactElements) )
                       cleanupAndThrow(mesh_p, "get(ContactElements_Owner)");
                   // ContactElements_Color
                   if (! ( nc_var_temp = dataFile.get_var("ContactElements_Color")) )
                       cleanupAndThrow(mesh_p, "get_var(ContactElements_Color)");
                   if (! nc_var_temp->get(&mesh_p->ContactElements->Color[0], num_ContactElements) )
                       cleanupAndThrow(mesh_p, "get(ContactElements_Color)");
                   // Now we need to adjust maxColor
                   index_t mc=mesh_p->ContactElements->Color[0];
                   for (index_t i=1;i<num_ContactElements;++i) {
                       if (mc<mesh_p->ContactElements->Color[i]) {
                           mc = mesh_p->ContactElements->Color[i];
                       }
                   }
                   mesh_p->ContactElements->maxColor=mc;
                   // ContactElements_Nodes
                   int *ContactElements_Nodes = TMPMEMALLOC(num_ContactElements*num_ContactElements_numNodes,int);
                   if (!(nc_var_temp = dataFile.get_var("ContactElements_Nodes"))) {
                       TMPMEMFREE(ContactElements_Nodes);
                       cleanupAndThrow(mesh_p, "get_var(ContactElements_Nodes)");
                   }
                   if (! nc_var_temp->get(&(ContactElements_Nodes[0]), num_ContactElements, num_ContactElements_numNodes) ) {
                       TMPMEMFREE(ContactElements_Nodes);
                       cleanupAndThrow(mesh_p, "get(ContactElements_Nodes)");
                   }
                   // Copy temp array into mesh_p->ContactElements->Nodes
                   for (int i=0; i<num_ContactElements; i++) {
                       for (int j=0; j<num_ContactElements_numNodes; j++) {
                           mesh_p->ContactElements->Nodes[INDEX2(j,i,num_ContactElements_numNodes)]= ContactElements_Nodes[INDEX2(j,i,num_ContactElements_numNodes)];
                       }
                   }
                   TMPMEMFREE(ContactElements_Nodes);
               } /* num_ContactElements>0 */
               mesh_p->ContactElements->updateTagList();
           }
        }

        /* get the Points (nodal elements) */
        if (noError()) {
            const_ReferenceElementSet_ptr refPoints(new ReferenceElementSet(
                        (ElementTypeId)Points_TypeId, order, reduced_order));
            if (noError())  {
                mesh_p->Points=new ElementFile(refPoints, mpi_info);
            }
            if (noError())
                mesh_p->Points->allocTable(num_Points);
            if (noError()) {
                mesh_p->Points->minColor=0;
                mesh_p->Points->maxColor=num_Points-1;
                if (num_Points>0) {
                   // Points_Id
                   if (! ( nc_var_temp = dataFile.get_var("Points_Id")))
                       cleanupAndThrow(mesh_p, "get_var(Points_Id)");
                   if (! nc_var_temp->get(&mesh_p->Points->Id[0], num_Points))
                       cleanupAndThrow(mesh_p, "get(Points_Id)");
                   // Points_Tag
                   if (! ( nc_var_temp = dataFile.get_var("Points_Tag")))
                       cleanupAndThrow(mesh_p, "get_var(Points_Tag)");
                   if (! nc_var_temp->get(&mesh_p->Points->Tag[0], num_Points))
                       cleanupAndThrow(mesh_p, "get(Points_Tag)");
                   // Points_Owner
                   if (! ( nc_var_temp = dataFile.get_var("Points_Owner")))
                       cleanupAndThrow(mesh_p, "get_var(Points_Owner)");
                   if (!nc_var_temp->get(&mesh_p->Points->Owner[0], num_Points))
                       cleanupAndThrow(mesh_p, "get(Points_Owner)");
                   // Points_Color
                   if (! ( nc_var_temp = dataFile.get_var("Points_Color")))
                       cleanupAndThrow(mesh_p, "get_var(Points_Color)");
                   if (!nc_var_temp->get(&mesh_p->Points->Color[0], num_Points))
                       cleanupAndThrow(mesh_p, "get(Points_Color)");
                   // Now we need to adjust maxColor
                   index_t mc=mesh_p->Points->Color[0];
                   for (index_t i=1;i<num_Points;++i) {
                       if (mc<mesh_p->Points->Color[i]) {
                           mc = mesh_p->Points->Color[i];
                       }
                   }
                   mesh_p->Points->maxColor=mc;
                   // Points_Nodes
                   int *Points_Nodes = TMPMEMALLOC(num_Points,int);
                   if (!(nc_var_temp = dataFile.get_var("Points_Nodes"))) {
                       TMPMEMFREE(Points_Nodes);
                       cleanupAndThrow(mesh_p, "get_var(Points_Nodes)");
                   }
                   if (! nc_var_temp->get(&(Points_Nodes[0]), num_Points) ) {
                       TMPMEMFREE(Points_Nodes);
                       cleanupAndThrow(mesh_p, "get(Points_Nodes)");
                   }
                   // Copy temp array into mesh_p->Points->Nodes
                   for (int i=0; i<num_Points; i++) {
                       mesh_p->Points->Id[mesh_p->Points->Nodes[INDEX2(0,i,1)]] = Points_Nodes[i];
                   }
                   TMPMEMFREE(Points_Nodes);
                } /* num_Points>0 */
                mesh_p->Points->updateTagList();
            }
        }

        /* get the tags */
        if (noError()) {
          if (num_Tags>0) {
            // Temp storage to gather node IDs
            int *Tags_keys = TMPMEMALLOC(num_Tags, int);
            char name_temp[4096];
            int i;

            // Tags_keys
            if (! ( nc_var_temp = dataFile.get_var("Tags_keys")) ) {
                TMPMEMFREE(Tags_keys);
                cleanupAndThrow(mesh_p, "get_var(Tags_keys)");
            }
            if (! nc_var_temp->get(&Tags_keys[0], num_Tags) ) {
                TMPMEMFREE(Tags_keys);
                cleanupAndThrow(mesh_p, "get(Tags_keys)");
            }
            for (i=0; i<num_Tags; i++) {
              // Retrieve tag name
              sprintf(name_temp, "Tags_name_%d", i);
              if (! (attr=dataFile.get_att(name_temp)) ) {
                  TMPMEMFREE(Tags_keys);
                  sprintf(error_msg,"get_att(%s)", name_temp);
                  cleanupAndThrow(mesh_p, error_msg);
              }
              boost::scoped_array<char> name(attr->as_string(0));
              delete attr;
              mesh_p->addTagMap(name.get(), Tags_keys[i]);
            }
            TMPMEMFREE(Tags_keys);
          }
        }
   
        if (noError()) {
            // Nodes_DofDistribution
            std::vector<int> first_DofComponent(mpi_size+1);
            if (! (nc_var_temp = dataFile.get_var("Nodes_DofDistribution")) ) {
                cleanupAndThrow(mesh_p, "get_var(Nodes_DofDistribution)");
            }
            if (!nc_var_temp->get(&first_DofComponent[0], mpi_size+1)) {
                cleanupAndThrow(mesh_p, "get(Nodes_DofDistribution)");
            }

            // Nodes_NodeDistribution
            std::vector<int> first_NodeComponent(mpi_size+1);
            if (! (nc_var_temp = dataFile.get_var("Nodes_NodeDistribution")) ) {
                cleanupAndThrow(mesh_p, "get_var(Nodes_NodeDistribution)");
            }
            if (!nc_var_temp->get(&first_NodeComponent[0], mpi_size+1)) {
                cleanupAndThrow(mesh_p, "get(Nodes_NodeDistribution)");
            }
            mesh_p->createMappings(first_DofComponent, first_NodeComponent);
        }

    } /* noError() after new Mesh() */

    checkFinleyError();
    Domain_ptr dom(new MeshAdapter(mesh_p));

    blocktimer_increment("LoadMesh()", blocktimer_start);
    return dom;
#else
    throw DataException("loadMesh: not compiled with NetCDF. Please contact your installation manager.");
#endif /* USE_NETCDF */
  }

  Domain_ptr readMesh(const std::string& fileName,
                      int integrationOrder,
                      int reducedIntegrationOrder,
                      bool optimize)
  {
    if (fileName.size() == 0 )
        throw DataException("Null file name!");

    double blocktimer_start = blocktimer_time();
    Mesh* fMesh=Mesh::read(fileName, integrationOrder, reducedIntegrationOrder, optimize);
    checkFinleyError();
    blocktimer_increment("ReadMesh()", blocktimer_start);
    return Domain_ptr(new MeshAdapter(fMesh));
  }

  Domain_ptr readGmsh(const std::string& fileName,
                                     int numDim,
                                     int integrationOrder,
                                     int reducedIntegrationOrder,
                                     bool optimize,
                                     bool useMacroElements)
  {
    if (fileName.size() == 0 )
        throw DataException("Null file name!");

    double blocktimer_start = blocktimer_time();
    Mesh* fMesh=Mesh::readGmsh(fileName, numDim, integrationOrder, reducedIntegrationOrder, optimize, useMacroElements);
    checkFinleyError();
    blocktimer_increment("ReadGmsh()", blocktimer_start);
    return Domain_ptr(new MeshAdapter(fMesh));
  }

/*  AbstractContinuousDomain* brick(int n0,int n1,int n2,int order,*/
  Domain_ptr brick(esysUtils::JMPI& info, int n0, int n1, int n2, int order,
                   double l0, double l1, double l2,
                   bool periodic0, bool periodic1, bool periodic2,
                   int integrationOrder, int reducedIntegrationOrder,
                   bool useElementsOnFace, bool useFullElementOrder,
                   bool optimize, const std::vector<double>& points,
                   const std::vector<int>& tags,
                   const std::map<std::string, int>& tagnamestonums
		  )
  {
    const int numElements[] = {n0, n1, n2};
    const double length[] = {l0, l1, l2};
    const bool periodic[] = {periodic0, periodic1, periodic2};

    Mesh* fMesh = NULL;
    if (order==1) {
        fMesh=RectangularMesh_Hex8(numElements, length, periodic,
                integrationOrder, reducedIntegrationOrder,
                useElementsOnFace, useFullElementOrder, optimize,
		info);
    } else if (order==2) {
        fMesh=RectangularMesh_Hex20(numElements, length, periodic,
                integrationOrder, reducedIntegrationOrder,
                useElementsOnFace, useFullElementOrder, false, optimize,
		info);
    } else if (order==-1) {
        fMesh=RectangularMesh_Hex20(numElements, length, periodic,
                integrationOrder, reducedIntegrationOrder,
                useElementsOnFace, useFullElementOrder, true, optimize,
		info);
    } else {
        stringstream message;
        message << "Illegal interpolation order " << order;
        throw FinleyAdapterException(message.str());
    }

    // Convert any finley errors into a C++ exception
    checkFinleyError();
    MeshAdapter* dom = new MeshAdapter(fMesh);
    dom->addDiracPoints(points, tags);
    Mesh* out=dom->getMesh().get();     
    for (map<string, int>::const_iterator it=tagnamestonums.begin();it!=tagnamestonums.end();++it)
    {
        out->addTagMap(it->first.c_str(), it->second);
    }
    out->Points->updateTagList();
    return Domain_ptr(dom);
  }

  Domain_ptr brick_driver(const boost::python::list& args)
  {
      using boost::python::extract;

      // we need to convert lists to stl vectors
      boost::python::list pypoints=extract<boost::python::list>(args[15]);
      boost::python::list pytags=extract<boost::python::list>(args[16]);
      int numpts=extract<int>(pypoints.attr("__len__")());
      int numtags=extract<int>(pytags.attr("__len__")());
      vector<double> points;
      vector<int> tags;
      tags.resize(numtags, -1);
      for (int i=0;i<numpts;++i) {
          boost::python::object temp=pypoints[i];
          int l=extract<int>(temp.attr("__len__")());
          for (int k=0;k<l;++k) {
              points.push_back(extract<double>(temp[k]));           
          }
      }
      map<string, int> namestonums;
      int curmax=40; // bricks use up to 30
      for (int i=0;i<numtags;++i) {
          extract<int> ex_int(pytags[i]);
          extract<string> ex_str(pytags[i]);
          if (ex_int.check()) {
              tags[i]=ex_int();
              if (tags[i]>= curmax) {
                  curmax=tags[i]+1;
              }
          } else if (ex_str.check()) {
              string s=ex_str();
              map<string, int>::iterator it=namestonums.find(s);
              if (it!=namestonums.end()) {
                  // we have the tag already so look it up
                  tags[i]=it->second;
              } else {
                  namestonums[s]=curmax;
                  tags[i]=curmax;
                  curmax++;
              }
          } else {
              throw FinleyAdapterException("Error - Unable to extract tag value.");
          }
        
      }
      boost::python::object pworld=args[17];
      esysUtils::JMPI info;
      if (!pworld.is_none())
      {
	  extract<SubWorld_ptr> ex(pworld);
	  if (!ex.check())
	  {
	      throw FinleyAdapterException("Invalid escriptWorld parameter.");
	  }
	  info=ex()->getMPI();
      }
      else
      {
	  info=esysUtils::makeInfo(MPI_COMM_WORLD);

      }
      return brick(info, static_cast<int>(extract<float>(args[0])),
                   static_cast<int>(extract<float>(args[1])),
                   static_cast<int>(extract<float>(args[2])),
                   extract<int>(args[3]), extract<double>(args[4]),
                   extract<double>(args[5]), extract<double>(args[6]),
                   extract<int>(args[7]), extract<int>(args[8]),
                   extract<int>(args[9]), extract<int>(args[10]),
                   extract<int>(args[11]), extract<int>(args[12]),
                   extract<int>(args[13]), extract<int>(args[14]),
                   points, tags, namestonums
		  );
  }

  Domain_ptr rectangle(esysUtils::JMPI& info, int n0, int n1, int order,
                       double l0, double l1,
                       bool periodic0, bool periodic1,
                       int integrationOrder,
                       int reducedIntegrationOrder,
                       bool useElementsOnFace,
                       bool useFullElementOrder,
                       bool optimize,
                       const vector<double>& points,
                       const vector<int>& tags,
                       const std::map<std::string, int>& tagnamestonums
		      )
  {
    const int numElements[] = {n0, n1};
    const double length[] = {l0, l1};
    const bool periodic[] = {periodic0, periodic1};

    Mesh* fMesh = NULL;
    if (order==1) {
        fMesh=RectangularMesh_Rec4(numElements, length, periodic,
                integrationOrder, reducedIntegrationOrder,
                useElementsOnFace, useFullElementOrder, optimize,
		info);
    } else if (order==2) {
        fMesh=RectangularMesh_Rec8(numElements, length, periodic,
                integrationOrder, reducedIntegrationOrder,
                useElementsOnFace,useFullElementOrder, false, optimize,
		info);
    } else if (order==-1) {
        fMesh=RectangularMesh_Rec8(numElements, length, periodic,
                integrationOrder, reducedIntegrationOrder,
                useElementsOnFace, useFullElementOrder, true, optimize,
		info);
    } else {
        stringstream message;
        message << "Illegal interpolation order " << order;
        throw FinleyAdapterException(message.str());
    }

    // Convert any finley errors into a C++ exception
    checkFinleyError();
    MeshAdapter* dom = new MeshAdapter(fMesh);
    dom->addDiracPoints(points, tags);
    Mesh* out=dom->getMesh().get();     
    for (map<string, int>::const_iterator it=tagnamestonums.begin();it!=tagnamestonums.end();++it)
    {
        out->addTagMap(it->first.c_str(), it->second);
    }
    out->Points->updateTagList();
    return Domain_ptr(dom);
  }

  Domain_ptr meshMerge(const boost::python::list& meshList)
  {
    // extract the meshes from meshList
    int num=boost::python::extract<int>(meshList.attr("__len__")());
    vector<Mesh*> meshes(num);
    for (int i=0; i<num; ++i) {
        AbstractContinuousDomain& meshListMember=boost::python::extract<AbstractContinuousDomain&>(meshList[i]);
        const MeshAdapter* finley_meshListMember=static_cast<const MeshAdapter*>(&meshListMember);
        meshes[i]=finley_meshListMember->getFinley_Mesh();
    }

    // merge the meshes
    Mesh* fMesh=Mesh_merge(meshes);

    // Convert any finley errors into a C++ exception
    checkFinleyError();
    return Domain_ptr(new MeshAdapter(fMesh));
  }

  Domain_ptr rectangle_driver(const boost::python::list& args)
  {
      using boost::python::extract;

      // we need to convert lists to stl vectors
      boost::python::list pypoints=extract<boost::python::list>(args[12]);
      boost::python::list pytags=extract<boost::python::list>(args[13]);
      int numpts=extract<int>(pypoints.attr("__len__")());
      int numtags=extract<int>(pytags.attr("__len__")());
      vector<double> points;
      vector<int> tags;
      tags.resize(numtags, -1);
      for (int i=0;i<numpts;++i)
      {
          boost::python::object temp=pypoints[i];
          int l=extract<int>(temp.attr("__len__")());
          for (int k=0;k<l;++k)
          {
              points.push_back(extract<double>(temp[k]));           
          }
      }
      map<string, int> tagstonames;
      int curmax=40;
      // but which order to assign tags to names?????
      for (int i=0;i<numtags;++i)
      {
          extract<int> ex_int(pytags[i]);
          extract<string> ex_str(pytags[i]);
          if (ex_int.check())
          {
              tags[i]=ex_int();
              if (tags[i]>= curmax)
              {
                  curmax=tags[i]+1;
              }
          } 
          else if (ex_str.check())
          {
              string s=ex_str();
              map<string, int>::iterator it=tagstonames.find(s);
              if (it!=tagstonames.end())
              {
                  // we have the tag already so look it up
                  tags[i]=it->second;
              }
              else
              {
                  tagstonames[s]=curmax;
                  tags[i]=curmax;
                  curmax++;
              }
          }
          else
          {
              throw FinleyAdapterException("Error - Unable to extract tag value.");
          }
      }
      boost::python::object pworld=args[14];
      esysUtils::JMPI info;
      if (!pworld.is_none())
      {
          extract<SubWorld_ptr> ex(pworld);
	  if (!ex.check())
	  {
	      throw FinleyAdapterException("Invalid escriptWorld parameter.");
          }
          info=ex()->getMPI();
      }
      else
      {
          info=esysUtils::makeInfo(MPI_COMM_WORLD);
      }

      return rectangle(info, static_cast<int>(extract<float>(args[0])),
                       static_cast<int>(extract<float>(args[1])),
                       extract<int>(args[2]), extract<double>(args[3]),
                       extract<double>(args[4]), extract<int>(args[5]),
                       extract<int>(args[6]), extract<int>(args[7]),
                       extract<int>(args[8]), extract<int>(args[9]),
                       extract<int>(args[10]), extract<int>(args[11]), 
                       points, tags, tagstonames
		       );
  }  


  Domain_ptr glueFaces(const boost::python::list& meshList,
                       double safety_factor, 
                       double tolerance,
                       bool optimize)
  {
    Mesh* fMesh=0;
    //
    // merge the meshes:
    Domain_ptr merged_meshes=meshMerge(meshList);

    //
    // glue the faces:
    const MeshAdapter* merged_finley_meshes=dynamic_cast<const MeshAdapter*>(merged_meshes.get());
    fMesh=merged_finley_meshes->getFinley_Mesh();
    fMesh->glueFaces(safety_factor, tolerance, optimize);

    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    return merged_meshes;
  }

  Domain_ptr joinFaces(const boost::python::list& meshList,
                       double safety_factor, 
                       double tolerance,
                       bool optimize)
  {
    Mesh* fMesh=0;
    //
    // merge the meshes:
    Domain_ptr merged_meshes=meshMerge(meshList);
    //
    // join the faces:
    const MeshAdapter* merged_finley_meshes=static_cast<const MeshAdapter*>(merged_meshes.get());
    fMesh=merged_finley_meshes->getFinley_Mesh();
    fMesh->joinFaces(safety_factor, tolerance, optimize);
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    return merged_meshes;
  }

  // end of namespace

}

