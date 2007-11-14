
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#ifdef PASO_MPI
#include <mpi.h>
#endif
#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif
#include "MeshAdapterFactory.h"
#include "FinleyError.h"
extern "C" {
#include "escript/blocktimer.h"
}

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

  AbstractContinuousDomain* loadMesh(const std::string& fileName)
  {
#ifdef USE_NETCDF
    Paso_MPIInfo *mpi_info = Paso_MPIInfo_alloc( MPI_COMM_WORLD );
    AbstractContinuousDomain* temp;
    Finley_Mesh *mesh_p=NULL;
    char error_msg[LenErrorMsg_MAX];
    // create a copy of the filename to overcome the non-constness of call
    // to Finley_Mesh_read
    // Win32 refactor
    char *fName = ((fileName.size()+1)>0) ? TMPMEMALLOC((fileName.size()+1),char) : (char*)NULL;
    strcpy(fName,fileName.c_str());

    printf("ksteube finley::loadMesh %s\n", fName);

    double blocktimer_start = blocktimer_time();
    Finley_resetError();

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
      Paso_MPIInfo_free( mpi_info );
      throw DataException("Error - loadMesh:: Could not read NetCDF file.");
    }

    // Read NetCDF integer attributes
    int mpi_size		= NetCDF_Get_Int_Attribute(&dataFile, fName, "mpi_size");
    int mpi_rank		= NetCDF_Get_Int_Attribute(&dataFile, fName, "mpi_rank");
    int numDim			= NetCDF_Get_Int_Attribute(&dataFile, fName, "numDim");
    int order			= NetCDF_Get_Int_Attribute(&dataFile, fName, "order");
    int reduced_order		= NetCDF_Get_Int_Attribute(&dataFile, fName, "reduced_order");
    int numNodes		= NetCDF_Get_Int_Attribute(&dataFile, fName, "numNodes");
    int num_Elements		= NetCDF_Get_Int_Attribute(&dataFile, fName, "num_Elements");
    int num_FaceElements	= NetCDF_Get_Int_Attribute(&dataFile, fName, "num_FaceElements");
    int num_ContactElements	= NetCDF_Get_Int_Attribute(&dataFile, fName, "num_ContactElements");
    int num_Points		= NetCDF_Get_Int_Attribute(&dataFile, fName, "num_Points");

    // Read mesh name
    if (! (attr=dataFile.get_att("Name")) ) {
      sprintf(error_msg,"Error retrieving mesh name from NetCDF file '%s'", fName);
      throw DataException(error_msg);
    }
    char *name = attr->as_string(0);
    delete attr;

    /* allocate mesh */
    mesh_p = Finley_Mesh_alloc(name,numDim,order,reduced_order,mpi_info);
    if (Finley_noError()) {

        /* read nodes */
        Finley_NodeFile_allocTable(mesh_p->Nodes, numNodes);
	// Nodes_Id
        if (! ( nc_var_temp = dataFile.get_var("Nodes_Id")) )
          throw DataException("Error - loadMesh:: unable to read Nodes_Id from netCDF file: " + *fName);
        if (! nc_var_temp->get(&mesh_p->Nodes->Id[0], numNodes) ) {
          free(&mesh_p->Nodes->Id);
          throw DataException("Error - loadMesh:: unable to recover Nodes_Id from NetCDF file: " + *fName);
        }
// printf("ksteube Nodes_Id: "); for (int i=0; i<numNodes; i++) { printf(" %d", mesh_p->Nodes->Id[i]); } printf("\n");
	// Nodes_Tag
        if (! ( nc_var_temp = dataFile.get_var("Nodes_Tag")) )
          throw DataException("Error - loadMesh:: unable to read Nodes_Tag from netCDF file: " + *fName);
        if (! nc_var_temp->get(&mesh_p->Nodes->Tag[0], numNodes) ) {
          free(&mesh_p->Nodes->Tag);
          throw DataException("Error - loadMesh:: unable to recover Nodes_Tag from NetCDF file: " + *fName);
        }
	// Nodes_gDOF
        if (! ( nc_var_temp = dataFile.get_var("Nodes_gDOF")) )
          throw DataException("Error - loadMesh:: unable to read Nodes_gDOF from netCDF file: " + *fName);
        if (! nc_var_temp->get(&mesh_p->Nodes->globalDegreesOfFreedom[0], numNodes) ) {
          free(&mesh_p->Nodes->globalDegreesOfFreedom);
          throw DataException("Error - loadMesh:: unable to recover Nodes_gDOF from NetCDF file: " + *fName);
        }
	// Nodes_gNI
        if (! ( nc_var_temp = dataFile.get_var("Nodes_gNI")) )
          throw DataException("Error - loadMesh:: unable to read Nodes_gNI from netCDF file: " + *fName);
        if (! nc_var_temp->get(&mesh_p->Nodes->globalNodesIndex[0], numNodes) ) {
          free(&mesh_p->Nodes->globalNodesIndex);
          throw DataException("Error - loadMesh:: unable to recover Nodes_gNI from NetCDF file: " + *fName);
        }
	// Nodes_grDfI
        if (! ( nc_var_temp = dataFile.get_var("Nodes_grDfI")) )
          throw DataException("Error - loadMesh:: unable to read Nodes_grDfI from netCDF file: " + *fName);
        if (! nc_var_temp->get(&mesh_p->Nodes->globalReducedDOFIndex[0], numNodes) ) {
          free(&mesh_p->Nodes->globalReducedDOFIndex);
          throw DataException("Error - loadMesh:: unable to recover Nodes_grDfI from NetCDF file: " + *fName);
        }
	// Nodes_grNI
        if (! ( nc_var_temp = dataFile.get_var("Nodes_grNI")) )
          throw DataException("Error - loadMesh:: unable to read Nodes_grNI from netCDF file: " + *fName);
        if (! nc_var_temp->get(&mesh_p->Nodes->globalReducedNodesIndex[0], numNodes) ) {
          free(&mesh_p->Nodes->globalReducedNodesIndex);
          throw DataException("Error - loadMesh:: unable to recover Nodes_grNI from NetCDF file: " + *fName);
        }
	// Nodes_Coordinates
        if (!(nc_var_temp = dataFile.get_var("Nodes_Coordinates"))) {
          free(&mesh_p->Nodes->Coordinates);
          throw DataException("Error - loadMesh:: unable to read Nodes_Coordinates from netCDF file: " + *fName);
        }
        if (! nc_var_temp->get(&(mesh_p->Nodes->Coordinates[0]), numNodes, numDim) ) {
          free(&mesh_p->Nodes->Coordinates);
          throw DataException("Error - load:: unable to recover Nodes_Coordinates from netCDF file: " + *fName);
        }

#if 0 /* Not yet...finish the above first */
        /* read elements */
        if (Finley_noError()) {
             mesh_p->Elements=Finley_ElementFile_alloc(typeID,mesh_p->order, mesh_p->reduced_order, mpi_info);
             if (Finley_noError()) {
                 Finley_ElementFile_allocTable(mesh_p->Elements, numEle);
                 mesh_p->Elements->minColor=0;
                 mesh_p->Elements->maxColor=numEle-1;
                 if (Finley_noError()) {
		 }
	     }
	}
#endif

        /* get the face elements */

        /* get the Contact face element */

        /* get the nodal element */

        /* get the name tags */

    } /* Finley_noError() after Finley_Mesh_alloc() */
   
#if 0 /* Not yet...finish the above first */
    if (Finley_noError()) Finley_Mesh_resolveNodeIds(mesh_p);
    if (Finley_noError()) Finley_Mesh_prepare(mesh_p, optimize);
#endif

    checkFinleyError();
    temp=new MeshAdapter(mesh_p);

    if (! Finley_noError()) {
      Finley_Mesh_free(mesh_p);
    }

    /* win32 refactor */
    TMPMEMFREE(fName);

    blocktimer_increment("LoadMesh()", blocktimer_start);
    return temp;
#else
    throw DataException("Error - loadMesh: is not compiled with NetCFD. Please contact your installation manager.");
#endif /* USE_NETCDF */
  }

  AbstractContinuousDomain* readMesh(const std::string& fileName,
  				     int integrationOrder,
                                     int reducedIntegrationOrder,
                                     int optimize)
  {
    //
    // create a copy of the filename to overcome the non-constness of call
    // to Finley_Mesh_read
    Finley_Mesh* fMesh=0;
    // Win32 refactor
    char *fName = ((fileName.size()+1)>0) ? TMPMEMALLOC((fileName.size()+1),char) : (char*)NULL;
    strcpy(fName,fileName.c_str());
    double blocktimer_start = blocktimer_time();

    fMesh=Finley_Mesh_read(fName,integrationOrder, reducedIntegrationOrder, (optimize ? TRUE : FALSE));
    checkFinleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    
    /* win32 refactor */
    TMPMEMFREE(fName);
    
    blocktimer_increment("ReadMesh()", blocktimer_start);
    return temp;
  }

  AbstractContinuousDomain* readGmsh(const std::string& fileName,
                                     int numDim,
                                     int integrationOrder,
                                     int reducedIntegrationOrder,
                                     int optimize)
  {
    //
    // create a copy of the filename to overcome the non-constness of call
    // to Finley_Mesh_read
    Finley_Mesh* fMesh=0;
    // Win32 refactor
    char *fName = ((fileName.size()+1)>0) ? TMPMEMALLOC((fileName.size()+1),char) : (char*)NULL;
    strcpy(fName,fileName.c_str());
    double blocktimer_start = blocktimer_time();

    fMesh=Finley_Mesh_readGmsh(fName, numDim, integrationOrder, reducedIntegrationOrder, (optimize ? TRUE : FALSE));
    checkFinleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    
    /* win32 refactor */
    TMPMEMFREE(fName);
    
    blocktimer_increment("ReadGmsh()", blocktimer_start);
    return temp;
  }

  AbstractContinuousDomain* brick(int n0,int n1,int n2,int order,
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
    int periodic[]={periodic0, periodic1, periodic2};

    //
    // linearInterpolation
    Finley_Mesh* fMesh=NULL;

    if (order==1) {
      fMesh=Finley_RectangularMesh_Hex8(numElements,length,periodic,integrationOrder,reducedIntegrationOrder,
					useElementsOnFace,useFullElementOrder,(optimize ? TRUE : FALSE)) ;
    } 
		else if (order==2) {
      fMesh=Finley_RectangularMesh_Hex20(numElements,length,periodic,integrationOrder,reducedIntegrationOrder,
					 useElementsOnFace,useFullElementOrder,(optimize ? TRUE : FALSE)) ;
    } else {
      stringstream temp;
      temp << "Illegal interpolation order: " << order;
      setFinleyError(VALUE_ERROR,temp.str().c_str());
    }
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    return temp;
  }
  AbstractContinuousDomain*  rectangle(int n0,int n1,int order,
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
    int periodic[]={periodic0, periodic1};

    Finley_Mesh* fMesh=0;
    if (order==1) {
      fMesh=Finley_RectangularMesh_Rec4(numElements, length,periodic,integrationOrder,reducedIntegrationOrder,
					useElementsOnFace,useFullElementOrder,(optimize ? TRUE : FALSE));
    }
    else if (order==2) {
      fMesh=Finley_RectangularMesh_Rec8(numElements,length,periodic,integrationOrder,reducedIntegrationOrder,
					useElementsOnFace,useFullElementOrder,(optimize ? TRUE : FALSE));
    }
    else {
      stringstream temp;
      temp << "Illegal interpolation order: " << order;
      setFinleyError(VALUE_ERROR,temp.str().c_str());
    }
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    return temp;
  }

  AbstractContinuousDomain* meshMerge(const boost::python::list& meshList)
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

    return temp;
  }
  AbstractContinuousDomain*  glueFaces(const boost::python::list& meshList,
                 	               double safety_factor, 
			               double tolerance,
                                       int optimize)
  {
    Finley_Mesh* fMesh=0;
    //
    // merge the meshes:
    AbstractContinuousDomain* merged_meshes=meshMerge(meshList);
    //
    // glue the faces:
    const MeshAdapter* merged_finley_meshes=static_cast<const MeshAdapter*>(merged_meshes);
    fMesh=merged_finley_meshes->getFinley_Mesh();
    Finley_Mesh_glueFaces(fMesh,safety_factor,tolerance,(optimize ? TRUE : FALSE));

    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    return merged_meshes;
  }
  AbstractContinuousDomain*  joinFaces(const boost::python::list& meshList,
			double safety_factor, 
			double tolerance,
                        int optimize)
  {
    Finley_Mesh* fMesh=0;
    //
    // merge the meshes:
    AbstractContinuousDomain* merged_meshes=meshMerge(meshList);
    //
    // join the faces:
    const MeshAdapter* merged_finley_meshes=static_cast<const MeshAdapter*>(merged_meshes);
    fMesh=merged_finley_meshes->getFinley_Mesh();
    Finley_Mesh_joinFaces(fMesh,safety_factor,tolerance, (optimize ? TRUE : FALSE));
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    return merged_meshes;
  }

  // end of namespace

}
