
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
#include "FinleyError.h"
extern "C" {
#include "esysUtils/blocktimer.h"
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

namespace finley {

#ifdef USE_NETCDF
  // A convenience method to retrieve an integer attribute from a NetCDF file
  int NetCDF_Get_Int_Attribute(NcFile *dataFile, char *fName, char *attr_name) {
    NcAtt *attr;
    char error_msg[LenErrorMsg_MAX];
    if (! (attr=dataFile->get_att(attr_name)) ) {
      sprintf(error_msg,"Error retrieving integer attribute '%s' from NetCDF file '%s'", attr_name, fName);
      throw DataException(error_msg);
    }
    int temp = attr->as_int(0);
    delete attr;
    return(temp);
  }
#endif

//   AbstractContinuousDomain* loadMesh(const std::string& fileName)
  Domain_ptr loadMesh(const std::string& fileName)
  {
#ifdef USE_NETCDF
    Esys_MPIInfo *mpi_info = Esys_MPIInfo_alloc( MPI_COMM_WORLD );
    AbstractContinuousDomain* temp;
    Finley_Mesh *mesh_p=NULL;
    char error_msg[LenErrorMsg_MAX];

    char *fName = Esys_MPI_appendRankToFileName(fileName.c_str(),
                                                mpi_info->size,
                                                mpi_info->rank);

    double blocktimer_start = blocktimer_time();
    Finley_resetError();
    int *first_DofComponent, *first_NodeComponent;

    // Open NetCDF file for reading
    NcAtt *attr;
    NcVar *nc_var_temp;
    // netCDF error handler
    NcError err(NcError::silent_nonfatal);
    // Create the NetCDF file.
    NcFile dataFile(fName, NcFile::ReadOnly);
    if (!dataFile.is_valid()) {
      sprintf(error_msg,"loadMesh: Opening file NetCDF %s for reading failed.", fName);
      Finley_setError(IO_ERROR,error_msg);
      Esys_MPIInfo_free( mpi_info );
      throw DataException(error_msg);
    }

    // Read NetCDF integer attributes
    int mpi_size			= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"mpi_size");
    int mpi_rank			= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"mpi_rank");
    int numDim				= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"numDim");
    int order				= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"order");
    int reduced_order			= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"reduced_order");
    int numNodes			= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"numNodes");
    int num_Elements			= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_Elements");
    int num_FaceElements		= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_FaceElements");
    int num_ContactElements		= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_ContactElements");
    int num_Points			= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_Points");
    int num_Elements_numNodes		= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_Elements_numNodes");
    int Elements_TypeId			= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"Elements_TypeId");
    int num_FaceElements_numNodes	= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_FaceElements_numNodes");
    int FaceElements_TypeId		= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"FaceElements_TypeId");
    int num_ContactElements_numNodes	= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_ContactElements_numNodes");
    int ContactElements_TypeId		= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"ContactElements_TypeId");
    int Points_TypeId			= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"Points_TypeId");
    int num_Tags			= NetCDF_Get_Int_Attribute(&dataFile, fName, (char *)"num_Tags");

    // Verify size and rank
    if (mpi_info->size != mpi_size) {
      sprintf(error_msg, "Error loadMesh - The NetCDF file '%s' can only be read on %d CPUs instead of %d", fName, mpi_size, mpi_info->size);
      throw DataException(error_msg);
    }
    if (mpi_info->rank != mpi_rank) {
      sprintf(error_msg, "Error loadMesh - The NetCDF file '%s' should be read on CPU #%d instead of %d", fName, mpi_rank, mpi_info->rank);
      throw DataException(error_msg);
    }

    // Read mesh name
    if (! (attr=dataFile.get_att("Name")) ) {
      sprintf(error_msg,"Error retrieving mesh name from NetCDF file '%s'", fName);
      throw DataException(error_msg);
    }
    char *name = attr->as_string(0);
    delete attr;

    string msgPrefix("Error in loadMesh: NetCDF operation failed - ");

    /* allocate mesh */
    mesh_p = Finley_Mesh_alloc(name,numDim,mpi_info);
    if (Finley_noError()) {

        /* read nodes */
        Finley_NodeFile_allocTable(mesh_p->Nodes, numNodes);
        // Nodes_Id
        if (! ( nc_var_temp = dataFile.get_var("Nodes_Id")) )
          throw DataException(msgPrefix+"get_var(Nodes_Id)");
        if (! nc_var_temp->get(&mesh_p->Nodes->Id[0], numNodes) ) {
          TMPMEMFREE(mesh_p->Nodes->Id);
          throw DataException("get(Nodes_Id)");
        }
        // Nodes_Tag
        if (! ( nc_var_temp = dataFile.get_var("Nodes_Tag")) )
          throw DataException("get_var(Nodes_Tag)");
        if (! nc_var_temp->get(&mesh_p->Nodes->Tag[0], numNodes) ) {
          TMPMEMFREE(mesh_p->Nodes->Tag);
          throw DataException("get(Nodes_Tag)");
        }
        // Nodes_gDOF
        if (! ( nc_var_temp = dataFile.get_var("Nodes_gDOF")) )
          throw DataException("get_var(Nodes_gDOF)");
        if (! nc_var_temp->get(&mesh_p->Nodes->globalDegreesOfFreedom[0], numNodes) ) {
          TMPMEMFREE(mesh_p->Nodes->globalDegreesOfFreedom);
          throw DataException("get(Nodes_gDOF)");
        }
        // Nodes_gNI
        if (! ( nc_var_temp = dataFile.get_var("Nodes_gNI")) )
          throw DataException("get_var(Nodes_gNI)");
        if (! nc_var_temp->get(&mesh_p->Nodes->globalNodesIndex[0], numNodes) ) {
          TMPMEMFREE(mesh_p->Nodes->globalNodesIndex);
          throw DataException("get(Nodes_gNI)");
        }
        // Nodes_grDfI
        if (! ( nc_var_temp = dataFile.get_var("Nodes_grDfI")) )
          throw DataException("get_var(Nodes_grDfI)");
        if (! nc_var_temp->get(&mesh_p->Nodes->globalReducedDOFIndex[0], numNodes) ) {
          TMPMEMFREE(mesh_p->Nodes->globalReducedDOFIndex);
          throw DataException("get(Nodes_grDfI)");
        }
        // Nodes_grNI
        if (! ( nc_var_temp = dataFile.get_var("Nodes_grNI")) )
          throw DataException("get_var(Nodes_grNI)");
        if (! nc_var_temp->get(&mesh_p->Nodes->globalReducedNodesIndex[0], numNodes) ) {
          TMPMEMFREE(mesh_p->Nodes->globalReducedNodesIndex);
          throw DataException("get(Nodes_grNI)");
        }
        // Nodes_Coordinates
        if (!(nc_var_temp = dataFile.get_var("Nodes_Coordinates"))) {
          TMPMEMFREE(mesh_p->Nodes->Coordinates);
          throw DataException("get_var(Nodes_Coordinates)");
        }
        if (! nc_var_temp->get(&(mesh_p->Nodes->Coordinates[0]), numNodes, numDim) ) {
          TMPMEMFREE(mesh_p->Nodes->Coordinates);
          throw DataException("get(Nodes_Coordinates)");
        }
        // Nodes_DofDistribution
        first_DofComponent = TMPMEMALLOC(mpi_size+1,index_t);
        if (! ( nc_var_temp = dataFile.get_var("Nodes_DofDistribution")) )
          throw DataException("get_var(Nodes_DofDistribution)");
        if (! nc_var_temp->get(&first_DofComponent[0], mpi_size+1) ) {
          throw DataException("get(Nodes_DofDistribution)");
        }

        // Nodes_NodeDistribution
        first_NodeComponent = TMPMEMALLOC(mpi_size+1,index_t);
        if (! ( nc_var_temp = dataFile.get_var("Nodes_NodeDistribution")) )
          throw DataException("get_var(Nodes_NodeDistribution)");
        if (! nc_var_temp->get(&first_NodeComponent[0], mpi_size+1) ) {
          throw DataException("get(Nodes_NodeDistribution)");
        }

        /* read elements */
        if (Finley_noError()) {
		  Finley_ReferenceElementSet  *refElements=	Finley_ReferenceElementSet_alloc((Finley_ElementTypeId)Elements_TypeId,order, reduced_order);
		  if (Finley_noError())  {
			  mesh_p->Elements=Finley_ElementFile_alloc(refElements, mpi_info);
		  }
		  Finley_ReferenceElementSet_dealloc(refElements);
          if (Finley_noError()) Finley_ElementFile_allocTable(mesh_p->Elements, num_Elements);
		  if (Finley_noError()) {
			  mesh_p->Elements->minColor=0;
			  mesh_p->Elements->maxColor=num_Elements-1;
			  if (num_Elements>0) {
                   // Elements_Id
                   if (! ( nc_var_temp = dataFile.get_var("Elements_Id")) )
                     throw DataException("get_var(Elements_Id)");
                   if (! nc_var_temp->get(&mesh_p->Elements->Id[0], num_Elements) ) {
                     TMPMEMFREE(mesh_p->Elements->Id);
                     throw DataException("get(Elements_Id)");
                   }
                   // Elements_Tag
                   if (! ( nc_var_temp = dataFile.get_var("Elements_Tag")) )
                     throw DataException("get_var(Elements_Tag)");
                   if (! nc_var_temp->get(&mesh_p->Elements->Tag[0], num_Elements) ) {
                     TMPMEMFREE(mesh_p->Elements->Tag);
                     throw DataException("get(Elements_Tag)");
                   }
                   // Elements_Owner
                   if (! ( nc_var_temp = dataFile.get_var("Elements_Owner")) )
                     throw DataException("get_var(Elements_Owner)");
                   if (! nc_var_temp->get(&mesh_p->Elements->Owner[0], num_Elements) ) {
                     TMPMEMFREE(mesh_p->Elements->Owner);
                     throw DataException("get(Elements_Owner)");
                   }
                   // Elements_Color
                   if (! ( nc_var_temp = dataFile.get_var("Elements_Color")) )
                     throw DataException("get_var(Elements_Color)");
                   if (! nc_var_temp->get(&mesh_p->Elements->Color[0], num_Elements) ) {
                     TMPMEMFREE(mesh_p->Elements->Color);
                     throw DataException("get(Elements_Color)");
                   }
		   // Now we need to adjust the maxColor
		   index_t mc=mesh_p->Elements->Color[0];
		   for (index_t i=1;i<num_Elements;++i)
		   {
			if (mc<mesh_p->Elements->Color[i])
			{
			   mc= mesh_p->Elements->Color[i];
			}
		   }
		   mesh_p->Elements->maxColor=mc;
                   // Elements_Nodes
                   int *Elements_Nodes = TMPMEMALLOC(num_Elements*num_Elements_numNodes,int);
                   if (!(nc_var_temp = dataFile.get_var("Elements_Nodes"))) {
                     TMPMEMFREE(mesh_p->Elements->Nodes);
                     throw DataException("get_var(Elements_Nodes)");
                   }
                   if (! nc_var_temp->get(&(Elements_Nodes[0]), num_Elements, num_Elements_numNodes) ) {
                     TMPMEMFREE(Elements_Nodes);
                     throw DataException("get(Elements_Nodes)");
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
        if (Finley_noError()) {
		  Finley_ReferenceElementSet *refFaceElements=	Finley_ReferenceElementSet_alloc((Finley_ElementTypeId)FaceElements_TypeId	,order, reduced_order);
		  if (Finley_noError())  {
			  mesh_p->FaceElements=Finley_ElementFile_alloc(refFaceElements, mpi_info);
		  }
		  Finley_ReferenceElementSet_dealloc(refFaceElements);	
          if (Finley_noError()) Finley_ElementFile_allocTable(mesh_p->FaceElements, num_FaceElements);
		  if (Finley_noError()) {
			  mesh_p->FaceElements->minColor=0;
			  mesh_p->FaceElements->maxColor=num_FaceElements-1;
			  if (num_FaceElements>0) {
                   // FaceElements_Id
                   if (! ( nc_var_temp = dataFile.get_var("FaceElements_Id")) )
                     throw DataException("get_var(FaceElements_Id)");
                   if (! nc_var_temp->get(&mesh_p->FaceElements->Id[0], num_FaceElements) ) {
                     TMPMEMFREE(mesh_p->FaceElements->Id);
                     throw DataException("get(FaceElements_Id)");
                   }
                   // FaceElements_Tag
                   if (! ( nc_var_temp = dataFile.get_var("FaceElements_Tag")) )
                     throw DataException("get_var(FaceElements_Tag)");
                   if (! nc_var_temp->get(&mesh_p->FaceElements->Tag[0], num_FaceElements) ) {
                     TMPMEMFREE(mesh_p->FaceElements->Tag);
                     throw DataException("get(FaceElements_Tag)");
                   }
                   // FaceElements_Owner
                   if (! ( nc_var_temp = dataFile.get_var("FaceElements_Owner")) )
                     throw DataException("get_var(FaceElements_Owner)");
                   if (! nc_var_temp->get(&mesh_p->FaceElements->Owner[0], num_FaceElements) ) {
                     TMPMEMFREE(mesh_p->FaceElements->Owner);
                     throw DataException("get(FaceElements_Owner)");
                   }
                   // FaceElements_Color
                   if (! ( nc_var_temp = dataFile.get_var("FaceElements_Color")) )
                     throw DataException("get_var(FaceElements_Color)");
                   if (! nc_var_temp->get(&mesh_p->FaceElements->Color[0], num_FaceElements) ) {
                     TMPMEMFREE(mesh_p->FaceElements->Color);
                     throw DataException("get(FaceElements_Color)");
                   }
		   // Now we need to adjust the maxColor
		   index_t mc=mesh_p->FaceElements->Color[0];
		   for (index_t i=1;i<num_FaceElements;++i)
		   {
			if (mc<mesh_p->FaceElements->Color[i])
			{
			   mc= mesh_p->FaceElements->Color[i];
			}
		   }
		   mesh_p->FaceElements->maxColor=mc;
                   // FaceElements_Nodes
                   int *FaceElements_Nodes = TMPMEMALLOC(num_FaceElements*num_FaceElements_numNodes,int);
                   if (!(nc_var_temp = dataFile.get_var("FaceElements_Nodes"))) {
                     TMPMEMFREE(mesh_p->FaceElements->Nodes);
                     throw DataException("get_var(FaceElements_Nodes)");
                   }
                   if (! nc_var_temp->get(&(FaceElements_Nodes[0]), num_FaceElements, num_FaceElements_numNodes) ) {
                     TMPMEMFREE(FaceElements_Nodes);
                     throw DataException("get(FaceElements_Nodes)");
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

        /* get the Contact elements */
        if (Finley_noError()) {
		  Finley_ReferenceElementSet *refContactElements=	Finley_ReferenceElementSet_alloc((Finley_ElementTypeId)ContactElements_TypeId,order, reduced_order);
		  if (Finley_noError())  {
			  mesh_p->ContactElements=Finley_ElementFile_alloc(refContactElements, mpi_info);
		  }
		  Finley_ReferenceElementSet_dealloc(refContactElements);	
          if (Finley_noError()) Finley_ElementFile_allocTable(mesh_p->ContactElements, num_ContactElements);
		  if (Finley_noError()) {
			  mesh_p->ContactElements->minColor=0;
			  mesh_p->ContactElements->maxColor=num_ContactElements-1;
			  if (num_ContactElements>0) {
	           // ContactElements_Id
                   if (! ( nc_var_temp = dataFile.get_var("ContactElements_Id")) )
                     throw DataException("get_var(ContactElements_Id)");
                   if (! nc_var_temp->get(&mesh_p->ContactElements->Id[0], num_ContactElements) ) {
                     TMPMEMFREE(mesh_p->ContactElements->Id);
                     throw DataException("get(ContactElements_Id)");
                   }
	           // ContactElements_Tag
                   if (! ( nc_var_temp = dataFile.get_var("ContactElements_Tag")) )
                     throw DataException("get_var(ContactElements_Tag)");
                   if (! nc_var_temp->get(&mesh_p->ContactElements->Tag[0], num_ContactElements) ) {
                     TMPMEMFREE(mesh_p->ContactElements->Tag);
                     throw DataException("get(ContactElements_Tag)");
                   }
	           // ContactElements_Owner
                   if (! ( nc_var_temp = dataFile.get_var("ContactElements_Owner")) )
                     throw DataException("get_var(ContactElements_Owner)");
                   if (! nc_var_temp->get(&mesh_p->ContactElements->Owner[0], num_ContactElements) ) {
                     TMPMEMFREE(mesh_p->ContactElements->Owner);
                     throw DataException("get(ContactElements_Owner)");
                   }
	           // ContactElements_Color
                   if (! ( nc_var_temp = dataFile.get_var("ContactElements_Color")) )
                     throw DataException("get_var(ContactElements_Color)");
                   if (! nc_var_temp->get(&mesh_p->ContactElements->Color[0], num_ContactElements) ) {
                     TMPMEMFREE(mesh_p->ContactElements->Color);
                     throw DataException("get(ContactElements_Color)");
                   }
		   // Now we need to adjust the maxColor
		   index_t mc=mesh_p->ContactElements->Color[0];
		   for (index_t i=1;i<num_ContactElements;++i)
		   {
			if (mc<mesh_p->ContactElements->Color[i])
			{
			   mc= mesh_p->ContactElements->Color[i];
			}
		   }
		   mesh_p->ContactElements->maxColor=mc;
	           // ContactElements_Nodes
			   		int *ContactElements_Nodes = TMPMEMALLOC(num_ContactElements*num_ContactElements_numNodes,int);
                   if (!(nc_var_temp = dataFile.get_var("ContactElements_Nodes"))) {
                     TMPMEMFREE(mesh_p->ContactElements->Nodes);
                     throw DataException("get_var(ContactElements_Nodes)");
                   }
                   if (! nc_var_temp->get(&(ContactElements_Nodes[0]), num_ContactElements, num_ContactElements_numNodes) ) {
                     TMPMEMFREE(ContactElements_Nodes);
                     throw DataException("get(ContactElements_Nodes)");
                   }
				   // Copy temp array into mesh_p->ContactElements->Nodes
				   	for (int i=0; i<num_ContactElements; i++) {
						for (int j=0; j<num_ContactElements_numNodes; j++) {
							mesh_p->ContactElements->Nodes[INDEX2(j,i,num_ContactElements_numNodes)]= ContactElements_Nodes[INDEX2(j,i,num_ContactElements_numNodes)];
						}
					}
					TMPMEMFREE(ContactElements_Nodes);
				 } /* num_ContactElements>0 */
			  } 
		  
	}

        /* get the Points (nodal elements) */
        if (Finley_noError()) {
		  Finley_ReferenceElementSet *refPoints=	Finley_ReferenceElementSet_alloc((Finley_ElementTypeId)Points_TypeId,order, reduced_order);
		  if (Finley_noError())  {
			  mesh_p->Points=Finley_ElementFile_alloc(refPoints, mpi_info);
		  }
		  Finley_ReferenceElementSet_dealloc(refPoints);
          if (Finley_noError()) Finley_ElementFile_allocTable(mesh_p->Points, num_Points);
		  if (Finley_noError()) {
			  mesh_p->Points->minColor=0;
			  mesh_p->Points->maxColor=num_Points-1;
			  if (num_Points>0) {
	           // Points_Id
                   if (! ( nc_var_temp = dataFile.get_var("Points_Id")) )
                     throw DataException("get_var(Points_Id)");
                   if (! nc_var_temp->get(&mesh_p->Points->Id[0], num_Points) ) {
                     TMPMEMFREE(mesh_p->Points->Id);
                     throw DataException("get(Points_Id)");
                   }
	           // Points_Tag
                   if (! ( nc_var_temp = dataFile.get_var("Points_Tag")) )
                     throw DataException("get_var(Points_Tag)");
                   if (! nc_var_temp->get(&mesh_p->Points->Tag[0], num_Points) ) {
                     TMPMEMFREE(mesh_p->Points->Tag);
                     throw DataException("get(Points_Tag)");
                   }
	           // Points_Owner
                   if (! ( nc_var_temp = dataFile.get_var("Points_Owner")) )
                     throw DataException("get_var(Points_Owner)");
                   if (! nc_var_temp->get(&mesh_p->Points->Owner[0], num_Points) ) {
                     TMPMEMFREE(mesh_p->Points->Owner);
                     throw DataException("get(Points_Owner)");
                   }
	           // Points_Color
                   if (! ( nc_var_temp = dataFile.get_var("Points_Color")) )
                     throw DataException("get_var(Points_Color)");
                   if (! nc_var_temp->get(&mesh_p->Points->Color[0], num_Points) ) {
                     TMPMEMFREE(mesh_p->Points->Color);
                     throw DataException("get(Points_Color)");
                   }
		   // Now we need to adjust the maxColor
		   index_t mc=mesh_p->Points->Color[0];
		   for (index_t i=1;i<num_Points;++i)
		   {
			if (mc<mesh_p->Points->Color[i])
			{
			   mc= mesh_p->Points->Color[i];
			}
		   }
		   mesh_p->Points->maxColor=mc;
	           // Points_Nodes
		   int *Points_Nodes = TMPMEMALLOC(num_Points,int);
                   if (!(nc_var_temp = dataFile.get_var("Points_Nodes"))) {
                     TMPMEMFREE(mesh_p->Points->Nodes);
                     throw DataException("get_var(Points_Nodes)");
                   }
                   if (! nc_var_temp->get(&(Points_Nodes[0]), num_Points) ) {
                     TMPMEMFREE(Points_Nodes);
                     throw DataException("get(Points_Nodes)");
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
        if (Finley_noError()) {
          if (num_Tags>0) {
            // Temp storage to gather node IDs
            int *Tags_keys = TMPMEMALLOC(num_Tags, int);
            char name_temp[4096];
	    int i;

	    // Tags_keys
            if (! ( nc_var_temp = dataFile.get_var("Tags_keys")) )
              throw DataException("get_var(Tags_keys)");
            if (! nc_var_temp->get(&Tags_keys[0], num_Tags) ) {
              TMPMEMFREE(Tags_keys);
              throw DataException("get(Tags_keys)");
            }
	    for (i=0; i<num_Tags; i++) {
              // Retrieve tag name
              sprintf(name_temp, "Tags_name_%d", i);
              if (! (attr=dataFile.get_att(name_temp)) ) {
                sprintf(error_msg,"Error retrieving tag name from NetCDF file '%s'", fName);
                throw DataException(error_msg);
              }
              char *name = attr->as_string(0);
              delete attr;
              Finley_Mesh_addTagMap(mesh_p, name, Tags_keys[i]);
	    }
	  }
	}
   
        if (Finley_noError()) Finley_Mesh_createMappings(mesh_p, first_DofComponent, first_NodeComponent);
        TMPMEMFREE(first_DofComponent);
        TMPMEMFREE(first_NodeComponent);

    } /* Finley_noError() after Finley_Mesh_alloc() */

    checkFinleyError();
    temp=new MeshAdapter(mesh_p);

    if (! Finley_noError()) {
      Finley_Mesh_free(mesh_p);
    }

    /* win32 refactor */
    TMPMEMFREE(fName);

    blocktimer_increment("LoadMesh()", blocktimer_start);
    return temp->getPtr();
#else
    throw DataException("Error - loadMesh: is not compiled with NetCDF. Please contact your installation manager.");
#endif /* USE_NETCDF */
  }

  Domain_ptr readMesh(const std::string& fileName,
  				     int integrationOrder,
                                     int reducedIntegrationOrder,
                                     int optimize)
  {
    //
    // create a copy of the filename to overcome the non-constness of call
    // to Finley_Mesh_read
    Finley_Mesh* fMesh=0;
    // Win32 refactor
    if( fileName.size() == 0 )
    {
       throw DataException("Null file name!");
    }

    char *fName = TMPMEMALLOC(fileName.size()+1,char);
	
    strcpy(fName,fileName.c_str());
    double blocktimer_start = blocktimer_time();

    fMesh=Finley_Mesh_read(fName,integrationOrder, reducedIntegrationOrder, (optimize ? TRUE : FALSE));
    checkFinleyError();
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
    // to Finley_Mesh_read
    Finley_Mesh* fMesh=0;
    // Win32 refactor
    if( fileName.size() == 0 )
    {
       throw DataException("Null file name!");
    }

    char *fName = TMPMEMALLOC(fileName.size()+1,char);
	
    strcpy(fName,fileName.c_str());
    double blocktimer_start = blocktimer_time();

    fMesh=Finley_Mesh_readGmsh(fName, numDim, integrationOrder, reducedIntegrationOrder, (optimize ? TRUE : FALSE), (useMacroElements ? TRUE : FALSE));
    checkFinleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    
    /* win32 refactor */
    TMPMEMFREE(fName);
    
    blocktimer_increment("ReadGmsh()", blocktimer_start);
    return temp->getPtr();
  }

/*  AbstractContinuousDomain* brick(int n0,int n1,int n2,int order,*/
  Domain_ptr brick(int n0,int n1,int n2,int order,
		    double l0,double l1,double l2,
		    int periodic0,int periodic1,
		    int periodic2,
		    int integrationOrder,
                    int reducedIntegrationOrder,
		    int useElementsOnFace,
		    int useFullElementOrder,
                     int optimize,
		    const std::vector<double>& points,
		    const std::vector<int>& tags, const std::map<std::string, int>& tagnamestonums
		    )
  {
    int numElements[]={n0,n1,n2};
    double length[]={l0,l1,l2};
    int periodic[]={periodic0, periodic1, periodic2};

    //
    // linearInterpolation
    Finley_Mesh* fMesh=NULL;

    if (order==1) {
           fMesh=Finley_RectangularMesh_Hex8(numElements,length,periodic,integrationOrder,reducedIntegrationOrder,
					useElementsOnFace,useFullElementOrder,(optimize ? TRUE : FALSE)) ;
    } else if (order==2) {
           fMesh=Finley_RectangularMesh_Hex20(numElements,length,periodic,integrationOrder,reducedIntegrationOrder,
					      useElementsOnFace,useFullElementOrder,FALSE, (optimize ? TRUE : FALSE)) ;
    } else if (order==-1) {
           fMesh=Finley_RectangularMesh_Hex20(numElements,length,periodic,integrationOrder,reducedIntegrationOrder,
					      useElementsOnFace,useFullElementOrder,TRUE,(optimize ? TRUE : FALSE)) ;
    } else {
      stringstream temp;
      temp << "Illegal interpolation order: " << order;
      setFinleyError(VALUE_ERROR,temp.str().c_str());
    }
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    MeshAdapter* temp=new MeshAdapter(fMesh);
    temp->addDiracPoints(points, tags);
    Finley_Mesh* out=temp->getMesh().get();	
    for (map<string, int>::const_iterator it=tagnamestonums.begin();it!=tagnamestonums.end();++it)
    {
	  Finley_Mesh_addTagMap(out, it->first.c_str(), it->second);
    }        
    return temp->getPtr();
  }



  Domain_ptr brick_driver(const boost::python::list& args)
  {
      using boost::python::extract;

      // we need to convert lists to stl vectors
      boost::python::list pypoints=boost::python::extract<boost::python::list>(args[15]);
      boost::python::list pytags=boost::python::extract<boost::python::list>(args[16]);
      int numpts=boost::python::extract<int>(pypoints.attr("__len__")());
      int numtags=boost::python::extract<int>(pytags.attr("__len__")());
      vector<double> points;
      vector<int> tags;
      tags.resize(numtags, -1);
      for (int i=0;i<numpts;++i)
      {
	  boost::python::object temp=pypoints[i];
	  int l=boost::python::extract<int>(temp.attr("__len__")());
	  for (int k=0;k<l;++k)
	  {
	      points.push_back(extract<double>(temp[k]));	    
	  }
      }
      map<string, int> namestonums;
      int curmax=40;				// bricks use up to 30
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
	      map<string, int>::iterator it=namestonums.find(s);
	      if (it!=namestonums.end())
	      {
		  // we have the tag already so look it up
		  tags[i]=it->second;
	      }
	      else
	      {
		  namestonums[s]=curmax;
		  tags[i]=curmax;
		  curmax++;
	      }
	  }
	  else
	  {
	      throw FinleyAdapterException("Error - Unable to extract tag value.");
	  }
	
      }
      
      Domain_ptr res=brick(boost::python::extract<int>(args[0]), boost::python::extract<int>(args[1]),
			   boost::python::extract<int>(args[2]), 
		   boost::python::extract<int>(args[3]), boost::python::extract<double>(args[4]),
		   boost::python::extract<double>(args[5]),
		   boost::python::extract<double>(args[6]), boost::python::extract<int>(args[7]),
		   boost::python::extract<int>(args[8]),
		   boost::python::extract<int>(args[9]),
		   boost::python::extract<int>(args[10]), boost::python::extract<int>(args[11]),
		   boost::python::extract<int>(args[12]), 
		   boost::python::extract<int>(args[13]), boost::python::extract<int>(args[14]),
		   points,
		   tags, namestonums);

      return res;
  }  
    
/*  AbstractContinuousDomain*  rectangle(int n0,int n1,int order,*/
  Domain_ptr  rectangle(int n0,int n1,int order,
			double l0, double l1,
			int periodic0,int periodic1,
			int integrationOrder,
                        int reducedIntegrationOrder,
			int useElementsOnFace,
		        int useFullElementOrder,
                        int optimize,
			const vector<double>& points,
			const vector<int>& tags,
			const std::map<std::string, int>& tagnamestonums)
  {
    int numElements[]={n0,n1};
    double length[]={l0,l1};
    int periodic[]={periodic0, periodic1};

    Finley_Mesh* fMesh=0;
    if (order==1) {
            fMesh=Finley_RectangularMesh_Rec4(numElements, length,periodic,integrationOrder,reducedIntegrationOrder,
					useElementsOnFace,useFullElementOrder,(optimize ? TRUE : FALSE));
    } else if (order==2) {
            fMesh=Finley_RectangularMesh_Rec8(numElements,length,periodic,integrationOrder,reducedIntegrationOrder,
					      useElementsOnFace,useFullElementOrder,FALSE,(optimize ? TRUE : FALSE));
    } else if (order==-1) {
            fMesh=Finley_RectangularMesh_Rec8(numElements,length,periodic,integrationOrder,reducedIntegrationOrder,
					      useElementsOnFace,useFullElementOrder,TRUE,(optimize ? TRUE : FALSE));
    } else {
      stringstream temp;
      temp << "Illegal interpolation order: " << order;
      setFinleyError(VALUE_ERROR,temp.str().c_str());
    }
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    MeshAdapter* temp=new MeshAdapter(fMesh);
    temp->addDiracPoints(points, tags);
    Finley_Mesh* out=temp->getMesh().get();	
    for (map<string, int>::const_iterator it=tagnamestonums.begin();it!=tagnamestonums.end();++it)
    {
	  Finley_Mesh_addTagMap(out, it->first.c_str(), it->second);
    }
    Finley_ElementFile_setTagsInUse(out->Points);
    return temp->getPtr();
  }

  Domain_ptr meshMerge(const boost::python::list& meshList)
  {
    Finley_Mesh* fMesh=0;
    //
    // extract the meshes from meshList
    int numMsh=boost::python::extract<int>(meshList.attr("__len__")());
    Finley_Mesh **mshes = (numMsh) ? TMPMEMALLOC(numMsh,Finley_Mesh*) : (Finley_Mesh**)NULL;
    for (int i=0;i<numMsh;++i) {
         AbstractContinuousDomain& meshListMember=boost::python::extract<AbstractContinuousDomain&>(meshList[i]);
         const MeshAdapter* finley_meshListMember=static_cast<const MeshAdapter*>(&meshListMember);
         mshes[i]=finley_meshListMember->getFinley_Mesh();
    }
    //
    // merge the meshes:
    fMesh=Finley_Mesh_merge(numMsh,mshes);
	  TMPMEMFREE(mshes);
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);

    return temp->getPtr();
  }

  Domain_ptr rectangle_driver(const boost::python::list& args)
  {
      using boost::python::extract;

      // we need to convert lists to stl vectors
      boost::python::list pypoints=boost::python::extract<boost::python::list>(args[12]);
      boost::python::list pytags=boost::python::extract<boost::python::list>(args[13]);
      int numpts=boost::python::extract<int>(pypoints.attr("__len__")());
      int numtags=boost::python::extract<int>(pytags.attr("__len__")());
      vector<double> points;
      vector<int> tags;
      tags.resize(numtags, -1);
      for (int i=0;i<numpts;++i)
      {
	  boost::python::object temp=pypoints[i];
	  int l=boost::python::extract<int>(temp.attr("__len__")());
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
      
      return rectangle(boost::python::extract<int>(args[0]), boost::python::extract<int>(args[1]), boost::python::extract<int>(args[2]), 
		   boost::python::extract<double>(args[3]),
		   boost::python::extract<double>(args[4]),
		   boost::python::extract<int>(args[5]), boost::python::extract<int>(args[6]),
		   boost::python::extract<int>(args[7]),
		   boost::python::extract<int>(args[8]),
		   boost::python::extract<int>(args[9]), boost::python::extract<int>(args[10]),
		   boost::python::extract<int>(args[11]), 
		   points,
		   tags, tagstonames);
  }  


  Domain_ptr  glueFaces(const boost::python::list& meshList,
                 	               double safety_factor, 
			               double tolerance,
                                       int optimize)
  {
    Finley_Mesh* fMesh=0;
    //
    // merge the meshes:
    Domain_ptr merged_meshes=meshMerge(meshList);

    //
    // glue the faces:
    const MeshAdapter* merged_finley_meshes=dynamic_cast<const MeshAdapter*>(merged_meshes.get());
    fMesh=merged_finley_meshes->getFinley_Mesh();
    Finley_Mesh_glueFaces(fMesh,safety_factor,tolerance,(optimize ? TRUE : FALSE));

    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    return merged_meshes->getPtr();
  }
  Domain_ptr  joinFaces(const boost::python::list& meshList,
			double safety_factor, 
			double tolerance,
                        int optimize)
  {
    Finley_Mesh* fMesh=0;
    //
    // merge the meshes:
    Domain_ptr merged_meshes=meshMerge(meshList);
    //
    // join the faces:
    const MeshAdapter* merged_finley_meshes=static_cast<const MeshAdapter*>(merged_meshes.get());
    fMesh=merged_finley_meshes->getFinley_Mesh();
    Finley_Mesh_joinFaces(fMesh,safety_factor,tolerance, (optimize ? TRUE : FALSE));
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    return merged_meshes->getPtr();
  }

  // end of namespace

}
