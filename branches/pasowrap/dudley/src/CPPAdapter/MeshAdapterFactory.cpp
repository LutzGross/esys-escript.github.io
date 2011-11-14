
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include "MeshAdapterFactory.h"
#include "DudleyError.h"
extern "C" {
#include "esysUtils/blocktimer.h"
#include "dudley/Dudley.h"
#include "dudley/Mesh.h"
#include "dudley/TriangularMesh.h"
#ifdef ESYS_MPI
#include "esysUtils/Esys_MPI.h"
#endif
}

#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

#include <boost/python/extract.hpp>

#include <sstream>

using namespace std;
using namespace escript;

namespace dudley {

#ifdef USE_NETCDF
  // A convenience method to retrieve an integer attribute from a NetCDF file
  int NetCDF_Get_Int_Attribute(NcFile *dataFile, char *fName, char *attr_name) {
    NcAtt *attr;
    char error_msg[LenErrorMsg_MAX];
    if (! (attr=dataFile->get_att(attr_name)) ) {
      sprintf(error_msg,"loadMesh: Error retrieving integer attribute '%s' from NetCDF file '%s'", attr_name, fName);
      throw DataException(error_msg);
    }
    int temp = attr->as_int(0);
    delete attr;
    return(temp);
  }
#endif

  inline void cleanupAndThrow(Dudley_Mesh* mesh, Esys_MPIInfo* info, string msg)
  {
      Dudley_Mesh_free(mesh);
      Esys_MPIInfo_free(info);
      string msgPrefix("loadMesh: NetCDF operation failed - ");
      throw DataException(msgPrefix+msg);
  }

//   AbstractContinuousDomain* loadMesh(const std::string& fileName)
  Domain_ptr loadMesh(const std::string& fileName)
  {
#ifdef USE_NETCDF
    Esys_MPIInfo *mpi_info = Esys_MPIInfo_alloc( MPI_COMM_WORLD );
    Dudley_Mesh *mesh_p=NULL;
    char error_msg[LenErrorMsg_MAX];

    char *fName = Esys_MPI_appendRankToFileName(fileName.c_str(),
                                                mpi_info->size,
                                                mpi_info->rank);

    double blocktimer_start = blocktimer_time();
    Dudley_resetError();
    int *first_DofComponent, *first_NodeComponent;

    // Open NetCDF file for reading
    NcAtt *attr;
    NcVar *nc_var_temp;
    // netCDF error handler
    NcError err(NcError::silent_nonfatal);
    // Create the NetCDF file.
    NcFile dataFile(fName, NcFile::ReadOnly);
    if (!dataFile.is_valid()) {
      sprintf(error_msg,"loadMesh: Opening NetCDF file '%s' for reading failed.", fName);
      Dudley_setError(IO_ERROR,error_msg);
      Esys_MPIInfo_free( mpi_info );
      throw DataException(error_msg);
    }

    // Read NetCDF integer attributes
    int mpi_size                        = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"mpi_size");
    int mpi_rank                        = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"mpi_rank");
    int numDim                          = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"numDim");
    int numNodes                        = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"numNodes");
    int num_Elements                    = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_Elements");
    int num_FaceElements                = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_FaceElements");
    int num_Points                      = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_Points");
    int num_Elements_numNodes           = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_Elements_numNodes");
    int Elements_TypeId                 = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"Elements_TypeId");
    int num_FaceElements_numNodes       = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_FaceElements_numNodes");
    int FaceElements_TypeId             = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"FaceElements_TypeId");
    int Points_TypeId                   = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"Points_TypeId");
    int num_Tags                        = NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_Tags");

    // Verify size and rank
    if (mpi_info->size != mpi_size) {
      sprintf(error_msg, "loadMesh: The NetCDF file '%s' can only be read on %d CPUs instead of %d", fName, mpi_size, mpi_info->size);
      throw DataException(error_msg);
    }
    if (mpi_info->rank != mpi_rank) {
      sprintf(error_msg, "loadMesh: The NetCDF file '%s' should be read on CPU #%d instead of %d", fName, mpi_rank, mpi_info->rank);
      throw DataException(error_msg);
    }

    // Read mesh name
    if (! (attr=dataFile.get_att("Name")) ) {
      sprintf(error_msg,"loadMesh: Error retrieving mesh name from NetCDF file '%s'", fName);
      throw DataException(error_msg);
    }
    char *name = attr->as_string(0);
    delete attr;

    TMPMEMFREE(fName);

    /* allocate mesh */
    mesh_p = Dudley_Mesh_alloc(name,numDim,mpi_info);
    if (Dudley_noError()) {

        /* read nodes */
        Dudley_NodeFile_allocTable(mesh_p->Nodes, numNodes);
        // Nodes_Id
        if (! ( nc_var_temp = dataFile.get_var("Nodes_Id")) )
            cleanupAndThrow(mesh_p, mpi_info, "get_var(Nodes_Id)");
        if (! nc_var_temp->get(&mesh_p->Nodes->Id[0], numNodes) )
            cleanupAndThrow(mesh_p, mpi_info, "get(Nodes_Id)");
        // Nodes_Tag
        if (! ( nc_var_temp = dataFile.get_var("Nodes_Tag")) )
            cleanupAndThrow(mesh_p, mpi_info, "get_var(Nodes_Tag)");
        if (! nc_var_temp->get(&mesh_p->Nodes->Tag[0], numNodes) )
            cleanupAndThrow(mesh_p, mpi_info, "get(Nodes_Tag)");
        // Nodes_gDOF
        if (! ( nc_var_temp = dataFile.get_var("Nodes_gDOF")) )
            cleanupAndThrow(mesh_p, mpi_info, "get_var(Nodes_gDOF)");
        if (! nc_var_temp->get(&mesh_p->Nodes->globalDegreesOfFreedom[0], numNodes) )
            cleanupAndThrow(mesh_p, mpi_info, "get(Nodes_gDOF)");
        // Nodes_gNI
        if (! ( nc_var_temp = dataFile.get_var("Nodes_gNI")) )
            cleanupAndThrow(mesh_p, mpi_info, "get_var(Nodes_gNI)");
        if (! nc_var_temp->get(&mesh_p->Nodes->globalNodesIndex[0], numNodes) )
            cleanupAndThrow(mesh_p, mpi_info, "get(Nodes_gNI)");
        // Nodes_grDfI
        if (! ( nc_var_temp = dataFile.get_var("Nodes_grDfI")) )
            cleanupAndThrow(mesh_p, mpi_info, "get_var(Nodes_grDfI)");
        if (! nc_var_temp->get(&mesh_p->Nodes->globalReducedDOFIndex[0], numNodes) )
            cleanupAndThrow(mesh_p, mpi_info, "get(Nodes_grDfI)");
        // Nodes_grNI
        if (! ( nc_var_temp = dataFile.get_var("Nodes_grNI")) )
            cleanupAndThrow(mesh_p, mpi_info, "get_var(Nodes_grNI)");
        if (! nc_var_temp->get(&mesh_p->Nodes->globalReducedNodesIndex[0], numNodes) )
            cleanupAndThrow(mesh_p, mpi_info, "get(Nodes_grNI)");
        // Nodes_Coordinates
        if (!(nc_var_temp = dataFile.get_var("Nodes_Coordinates")))
            cleanupAndThrow(mesh_p, mpi_info, "get_var(Nodes_Coordinates)");
        if (! nc_var_temp->get(&(mesh_p->Nodes->Coordinates[0]), numNodes, numDim) )
            cleanupAndThrow(mesh_p, mpi_info, "get(Nodes_Coordinates)");

        /* read elements */
        if (Dudley_noError()) {
            mesh_p->Elements=Dudley_ElementFile_alloc((Dudley_ElementTypeId)Elements_TypeId, mpi_info);
            if (Dudley_noError())
                Dudley_ElementFile_allocTable(mesh_p->Elements, num_Elements);
            if (Dudley_noError()) {
                mesh_p->Elements->minColor=0;
                mesh_p->Elements->maxColor=num_Elements-1;
                if (num_Elements>0) {
                   // Elements_Id
                   if (! ( nc_var_temp = dataFile.get_var("Elements_Id")) )
                       cleanupAndThrow(mesh_p, mpi_info, "get_var(Elements_Id)");
                   if (! nc_var_temp->get(&mesh_p->Elements->Id[0], num_Elements) )
                       cleanupAndThrow(mesh_p, mpi_info, "get(Elements_Id)");
                   // Elements_Tag
                   if (! ( nc_var_temp = dataFile.get_var("Elements_Tag")) )
                       cleanupAndThrow(mesh_p, mpi_info, "get_var(Elements_Tag)");
                   if (! nc_var_temp->get(&mesh_p->Elements->Tag[0], num_Elements) )
                       cleanupAndThrow(mesh_p, mpi_info, "get(Elements_Tag)");
                   // Elements_Owner
                   if (! ( nc_var_temp = dataFile.get_var("Elements_Owner")) )
                       cleanupAndThrow(mesh_p, mpi_info, "get_var(Elements_Owner)");
                   if (! nc_var_temp->get(&mesh_p->Elements->Owner[0], num_Elements) )
                       cleanupAndThrow(mesh_p, mpi_info, "get(Elements_Owner)");
                   // Elements_Color
                   if (! ( nc_var_temp = dataFile.get_var("Elements_Color")) )
                       cleanupAndThrow(mesh_p, mpi_info, "get_var(Elements_Color)");
                   if (! nc_var_temp->get(&mesh_p->Elements->Color[0], num_Elements) )
                       cleanupAndThrow(mesh_p, mpi_info, "get(Elements_Color)");
                   // Elements_Nodes
                   int *Elements_Nodes = TMPMEMALLOC(num_Elements*num_Elements_numNodes,int);
                   if (!(nc_var_temp = dataFile.get_var("Elements_Nodes"))) {
                       TMPMEMFREE(Elements_Nodes);
                       cleanupAndThrow(mesh_p, mpi_info, "get_var(Elements_Nodes)");
                   }
                   if (! nc_var_temp->get(&(Elements_Nodes[0]), num_Elements, num_Elements_numNodes) ) {
                       TMPMEMFREE(Elements_Nodes);
                       cleanupAndThrow(mesh_p, mpi_info, "get(Elements_Nodes)");
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
            }
        }

        /* get the face elements */
        if (Dudley_noError()) {
            mesh_p->FaceElements=Dudley_ElementFile_alloc((Dudley_ElementTypeId)FaceElements_TypeId, mpi_info);
            if (Dudley_noError())
                Dudley_ElementFile_allocTable(mesh_p->FaceElements, num_FaceElements);
            if (Dudley_noError()) {
                mesh_p->FaceElements->minColor=0;
                mesh_p->FaceElements->maxColor=num_FaceElements-1;
                if (num_FaceElements>0) {
                   // FaceElements_Id
                   if (! ( nc_var_temp = dataFile.get_var("FaceElements_Id")) )
                       cleanupAndThrow(mesh_p, mpi_info, "get_var(FaceElements_Id)");
                   if (! nc_var_temp->get(&mesh_p->FaceElements->Id[0], num_FaceElements) )
                       cleanupAndThrow(mesh_p, mpi_info, "get(FaceElements_Id)");
                   // FaceElements_Tag
                   if (! ( nc_var_temp = dataFile.get_var("FaceElements_Tag")) )
                       cleanupAndThrow(mesh_p, mpi_info, "get_var(FaceElements_Tag)");
                   if (! nc_var_temp->get(&mesh_p->FaceElements->Tag[0], num_FaceElements) )
                       cleanupAndThrow(mesh_p, mpi_info, "get(FaceElements_Tag)");
                   // FaceElements_Owner
                   if (! ( nc_var_temp = dataFile.get_var("FaceElements_Owner")) )
                       cleanupAndThrow(mesh_p, mpi_info, "get_var(FaceElements_Owner)");
                   if (! nc_var_temp->get(&mesh_p->FaceElements->Owner[0], num_FaceElements) )
                       cleanupAndThrow(mesh_p, mpi_info, "get(FaceElements_Owner)");
                   // FaceElements_Color
                   if (! ( nc_var_temp = dataFile.get_var("FaceElements_Color")) )
                       cleanupAndThrow(mesh_p, mpi_info, "get_var(FaceElements_Color)");
                   if (! nc_var_temp->get(&mesh_p->FaceElements->Color[0], num_FaceElements) )
                       cleanupAndThrow(mesh_p, mpi_info, "get(FaceElements_Color)");
                   // FaceElements_Nodes
                   int *FaceElements_Nodes = TMPMEMALLOC(num_FaceElements*num_FaceElements_numNodes,int);
                   if (!(nc_var_temp = dataFile.get_var("FaceElements_Nodes"))) {
                       TMPMEMFREE(FaceElements_Nodes);
                       cleanupAndThrow(mesh_p, mpi_info, "get_var(FaceElements_Nodes)");
                   }
                   if (! nc_var_temp->get(&(FaceElements_Nodes[0]), num_FaceElements, num_FaceElements_numNodes) ) {
                       TMPMEMFREE(FaceElements_Nodes);
                       cleanupAndThrow(mesh_p, mpi_info, "get(FaceElements_Nodes)");
                   }
                   // Copy temp array into mesh_p->FaceElements->Nodes
                   for (int i=0; i<num_FaceElements; i++) {
                       for (int j=0; j<num_FaceElements_numNodes; j++) {
                           mesh_p->FaceElements->Nodes[INDEX2(j,i,num_FaceElements_numNodes)] = FaceElements_Nodes[INDEX2(j,i,num_FaceElements_numNodes)];
                       }
                   }
                   TMPMEMFREE(FaceElements_Nodes);
                } /* num_FaceElements>0 */
            }
        }

        /* get the Points (nodal elements) */
        if (Dudley_noError()) {
            mesh_p->Points=Dudley_ElementFile_alloc((Dudley_ElementTypeId)Points_TypeId, mpi_info);
            if (Dudley_noError())
                Dudley_ElementFile_allocTable(mesh_p->Points, num_Points);
            if (Dudley_noError()) {
                mesh_p->Points->minColor=0;
                mesh_p->Points->maxColor=num_Points-1;
                if (num_Points>0) {
                   // Points_Id
                   if (! ( nc_var_temp = dataFile.get_var("Points_Id")))
                       cleanupAndThrow(mesh_p, mpi_info, "get_var(Points_Id)");
                   if (! nc_var_temp->get(&mesh_p->Points->Id[0], num_Points))
                       cleanupAndThrow(mesh_p, mpi_info, "get(Points_Id)");
                   // Points_Tag
                   if (! ( nc_var_temp = dataFile.get_var("Points_Tag")))
                       cleanupAndThrow(mesh_p, mpi_info, "get_var(Points_Tag)");
                   if (! nc_var_temp->get(&mesh_p->Points->Tag[0], num_Points))
                       cleanupAndThrow(mesh_p, mpi_info, "get(Points_Tag)");
                   // Points_Owner
                   if (! ( nc_var_temp = dataFile.get_var("Points_Owner")))
                       cleanupAndThrow(mesh_p, mpi_info, "get_var(Points_Owner)");
                   if (!nc_var_temp->get(&mesh_p->Points->Owner[0], num_Points))
                       cleanupAndThrow(mesh_p, mpi_info, "get(Points_Owner)");
                   // Points_Color
                   if (! ( nc_var_temp = dataFile.get_var("Points_Color")))
                       cleanupAndThrow(mesh_p, mpi_info, "get_var(Points_Color)");
                   if (!nc_var_temp->get(&mesh_p->Points->Color[0], num_Points))
                       cleanupAndThrow(mesh_p, mpi_info, "get(Points_Color)");
                   // Points_Nodes
                   int *Points_Nodes = TMPMEMALLOC(num_Points,int);
                   if (!(nc_var_temp = dataFile.get_var("Points_Nodes"))) {
                       TMPMEMFREE(Points_Nodes);
                       cleanupAndThrow(mesh_p, mpi_info, "get_var(Points_Nodes)");
                   }
                   if (! nc_var_temp->get(&(Points_Nodes[0]), num_Points) ) {
                       TMPMEMFREE(Points_Nodes);
                       cleanupAndThrow(mesh_p, mpi_info, "get(Points_Nodes)");
                   }
                   // Copy temp array into mesh_p->Points->Nodes
                   for (int i=0; i<num_Points; i++) {
                       mesh_p->Nodes->Id[mesh_p->Points->Nodes[INDEX2(0,i,1)]] = Points_Nodes[i];
                   }
                   TMPMEMFREE(Points_Nodes);
                } /* num_Points>0 */
            }
        }

        /* get the tags */
        if (Dudley_noError()) {
          if (num_Tags>0) {
            // Temp storage to gather node IDs
            int *Tags_keys = TMPMEMALLOC(num_Tags, int);
            char name_temp[4096];
            int i;

            // Tags_keys
            if (! ( nc_var_temp = dataFile.get_var("Tags_keys")) ) {
                TMPMEMFREE(Tags_keys);
                cleanupAndThrow(mesh_p, mpi_info, "get_var(Tags_keys)");
            }
            if (! nc_var_temp->get(&Tags_keys[0], num_Tags) ) {
                TMPMEMFREE(Tags_keys);
                cleanupAndThrow(mesh_p, mpi_info, "get(Tags_keys)");
            }
            for (i=0; i<num_Tags; i++) {
              // Retrieve tag name
              sprintf(name_temp, "Tags_name_%d", i);
              if (! (attr=dataFile.get_att(name_temp)) ) {
                  TMPMEMFREE(Tags_keys);
                  sprintf(error_msg,"get_att(%s)", name_temp);
                  cleanupAndThrow(mesh_p, mpi_info, error_msg);
              }
              char *name = attr->as_string(0);
              delete attr;
              Dudley_Mesh_addTagMap(mesh_p, name, Tags_keys[i]);
            }
            TMPMEMFREE(Tags_keys);
          }
        }
   
        if (Dudley_noError()) {
            // Nodes_DofDistribution
            first_DofComponent = TMPMEMALLOC(mpi_size+1,index_t);
            if (! ( nc_var_temp = dataFile.get_var("Nodes_DofDistribution")) ) {
                TMPMEMFREE(first_DofComponent);
                cleanupAndThrow(mesh_p, mpi_info, "get_var(Nodes_DofDistribution)");
            }
            if (! nc_var_temp->get(&first_DofComponent[0], mpi_size+1) ) {
                TMPMEMFREE(first_DofComponent);
                cleanupAndThrow(mesh_p, mpi_info, "get(Nodes_DofDistribution)");
            }

            // Nodes_NodeDistribution
            first_NodeComponent = TMPMEMALLOC(mpi_size+1,index_t);
            if (! ( nc_var_temp = dataFile.get_var("Nodes_NodeDistribution")) ) {
                TMPMEMFREE(first_DofComponent);
                TMPMEMFREE(first_NodeComponent);
                cleanupAndThrow(mesh_p, mpi_info, "get_var(Nodes_NodeDistribution)");
            }
            if (! nc_var_temp->get(&first_NodeComponent[0], mpi_size+1) ) {
                TMPMEMFREE(first_DofComponent);
                TMPMEMFREE(first_NodeComponent);
                cleanupAndThrow(mesh_p, mpi_info, "get(Nodes_NodeDistribution)");
            }
            Dudley_Mesh_createMappings(mesh_p, first_DofComponent, first_NodeComponent);
            TMPMEMFREE(first_DofComponent);
            TMPMEMFREE(first_NodeComponent);
        }

    } /* Dudley_noError() after Dudley_Mesh_alloc() */

    checkDudleyError();
    AbstractContinuousDomain* dom=new MeshAdapter(mesh_p);

    if (! Dudley_noError()) {
        Dudley_Mesh_free(mesh_p);
    }

    blocktimer_increment("LoadMesh()", blocktimer_start);
    return dom->getPtr();
#else
    throw DataException("loadMesh: not compiled with NetCDF. Please contact your installation manager.");
#endif /* USE_NETCDF */
  }

  Domain_ptr readMesh(const std::string& fileName,
                      int integrationOrder,
                      int reducedIntegrationOrder,
                      int optimize)
  {
    //
    // create a copy of the filename to overcome the non-constness of call
    // to Dudley_Mesh_read
    Dudley_Mesh* fMesh=0;
    // Win32 refactor
    if( fileName.size() == 0 )
    {
       throw DataException("Null file name!");
    }

    char *fName = TMPMEMALLOC(fileName.size()+1,char);
        
    strcpy(fName,fileName.c_str());
    double blocktimer_start = blocktimer_time();

    fMesh=Dudley_Mesh_read(fName,integrationOrder, reducedIntegrationOrder, (optimize ? TRUE : FALSE));
    checkDudleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    
    /* win32 refactor */
    TMPMEMFREE(fName);
    
    blocktimer_increment("ReadMesh()", blocktimer_start);
    return temp->getPtr();
  }

  Domain_ptr readGmsh(const std::string& fileName,
                                     int numDim,
                                     int integrationOrder,
                                     int reducedIntegrationOrder,
                                     int optimize,
                                     int useMacroElements)
  {
    //
    // create a copy of the filename to overcome the non-constness of call
    // to Dudley_Mesh_read
    Dudley_Mesh* fMesh=0;
    // Win32 refactor
    if( fileName.size() == 0 )
    {
       throw DataException("Null file name!");
    }

    char *fName = TMPMEMALLOC(fileName.size()+1,char);
        
    strcpy(fName,fileName.c_str());
    double blocktimer_start = blocktimer_time();

    fMesh=Dudley_Mesh_readGmsh(fName, numDim, integrationOrder, reducedIntegrationOrder, (optimize ? TRUE : FALSE), (useMacroElements ? TRUE : FALSE));
    checkDudleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    
    /* win32 refactor */
    TMPMEMFREE(fName);
    
    blocktimer_increment("ReadGmsh()", blocktimer_start);
    return temp->getPtr();
  }

  Domain_ptr brick(int n0,int n1,int n2,int order,
                   double l0,double l1,double l2,
                   int periodic0,int periodic1,
                   int periodic2,
                   int integrationOrder,
                   int reducedIntegrationOrder,
                   int useElementsOnFace,
                   int useFullElementOrder,
                   int optimize)
  {
    int numElements[]={n0,n1,n2};
    double length[]={l0,l1,l2};

    if (periodic0 || periodic1) // we don't support periodic boundary conditions
    {
        throw DudleyAdapterException("Dudley does not support periodic boundary conditions.");
    }
    else if (integrationOrder>3 || reducedIntegrationOrder>1)
    {
        throw DudleyAdapterException("Dudley does not support the requested integrationOrders.");
    }
    else if (useElementsOnFace || useFullElementOrder)
    {
        throw DudleyAdapterException("Dudley does not support useElementsOnFace or useFullElementOrder.");
    }
    if (order>1)
    {
        throw DudleyAdapterException("Dudley does not support element order greater than 1.");
    }

    //
    // linearInterpolation
    Dudley_Mesh* fMesh=NULL;

    fMesh=Dudley_TriangularMesh_Tet4(numElements, length, integrationOrder,
                        reducedIntegrationOrder, (optimize ? TRUE : FALSE));

    //
    // Convert any dudley errors into a C++ exception
    checkDudleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    return temp->getPtr();
  }

  Domain_ptr rectangle(int n0,int n1,int order,
                       double l0, double l1,
                       int periodic0,int periodic1,
                       int integrationOrder,
                       int reducedIntegrationOrder,
                       int useElementsOnFace,
                       int useFullElementOrder,
                       int optimize)
  {
    int numElements[]={n0,n1};
    double length[]={l0,l1};

    if (periodic0 || periodic1) // we don't support periodic boundary conditions
    {
        throw DudleyAdapterException("Dudley does not support periodic boundary conditions.");
    }
    else if (integrationOrder>3 || reducedIntegrationOrder>1)
    {
        throw DudleyAdapterException("Dudley does not support the requested integrationOrders.");
    }
    else if (useElementsOnFace || useFullElementOrder)
    {
        throw DudleyAdapterException("Dudley does not support useElementsOnFace or useFullElementOrder.");
    }

    if (order>1)
    {
        throw DudleyAdapterException("Dudley does not support element order greater than 1.");
    }
    Dudley_Mesh* fMesh=Dudley_TriangularMesh_Tri3(numElements, length,
          integrationOrder, reducedIntegrationOrder, (optimize ? TRUE : FALSE));
    //
    // Convert any dudley errors into a C++ exception
    checkDudleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);

    return temp->getPtr();
  }

  // end of namespace

}

