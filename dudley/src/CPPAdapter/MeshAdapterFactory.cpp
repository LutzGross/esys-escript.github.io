
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

#define ESNEEDPYTHON
#include "esysUtils/first.h"


#include "MeshAdapterFactory.h"
#include "DudleyError.h"
#include "esysUtils/blocktimer.h"
#include "dudley/Dudley.h"
#include "dudley/Mesh.h"
#include "dudley/TriangularMesh.h"
#ifdef ESYS_MPI
#include "esysUtils/Esys_MPI.h"
#endif

#include "escript/SubWorld.h"

#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

#include <boost/python/extract.hpp>
#include <boost/scoped_array.hpp>

#include <sstream>

using namespace std;
using namespace escript;

namespace dudley {

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

  inline void cleanupAndThrow(Dudley_Mesh* mesh, string msg)
  {
      Dudley_Mesh_free(mesh);
      string msgPrefix("loadMesh: NetCDF operation failed - ");
      throw DataException(msgPrefix+msg);
  }

//   AbstractContinuousDomain* loadMesh(const std::string& fileName)
  Domain_ptr loadMesh(const std::string& fileName)
  {
#ifdef USE_NETCDF
    esysUtils::JMPI mpi_info = esysUtils::makeInfo( MPI_COMM_WORLD );
    Dudley_Mesh *mesh_p=NULL;
    char error_msg[LenErrorMsg_MAX];

    std::string fName(esysUtils::appendRankToFileName(fileName, mpi_info->size,
                                                      mpi_info->rank));

    double blocktimer_start = blocktimer_time();
    Dudley_resetError();
    int *first_DofComponent, *first_NodeComponent;

    // Open NetCDF file for reading
    NcAtt *attr;
    NcVar *nc_var_temp;
    // netCDF error handler
    NcError err(NcError::silent_nonfatal);
    // Create the NetCDF file.
    NcFile dataFile(fName.c_str(), NcFile::ReadOnly);
    if (!dataFile.is_valid()) {
      sprintf(error_msg,"loadMesh: Opening NetCDF file '%s' for reading failed.", fName.c_str());
      Dudley_setError(IO_ERROR,error_msg);
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
    mesh_p = Dudley_Mesh_alloc(name.get(), numDim, mpi_info);
    if (Dudley_noError()) {

        /* read nodes */
        Dudley_NodeFile_allocTable(mesh_p->Nodes, numNodes);
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

        Dudley_NodeFile_setTagsInUse(mesh_p->Nodes);

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
                   // Elements_Nodes
                   int *Elements_Nodes = new int[num_Elements*num_Elements_numNodes];
                   if (!(nc_var_temp = dataFile.get_var("Elements_Nodes"))) {
                       delete[] Elements_Nodes;
                       cleanupAndThrow(mesh_p, "get_var(Elements_Nodes)");
                   }
                   if (! nc_var_temp->get(&(Elements_Nodes[0]), num_Elements, num_Elements_numNodes) ) {
                       delete[] Elements_Nodes;
                       cleanupAndThrow(mesh_p, "get(Elements_Nodes)");
                   }

                   // Copy temp array into mesh_p->Elements->Nodes
                   for (int i=0; i<num_Elements; i++) {
                       for (int j=0; j<num_Elements_numNodes; j++) {
                           mesh_p->Elements->Nodes[INDEX2(j,i,num_Elements_numNodes)]
                                = Elements_Nodes[INDEX2(j,i,num_Elements_numNodes)];
                       }
                   }
                   delete[] Elements_Nodes;
                   Dudley_ElementFile_setTagsInUse(mesh_p->Elements);
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
                   // FaceElements_Nodes
                   int *FaceElements_Nodes = new int[num_FaceElements*num_FaceElements_numNodes];
                   if (!(nc_var_temp = dataFile.get_var("FaceElements_Nodes"))) {
                       delete[] FaceElements_Nodes;
                       cleanupAndThrow(mesh_p, "get_var(FaceElements_Nodes)");
                   }
                   if (! nc_var_temp->get(&(FaceElements_Nodes[0]), num_FaceElements, num_FaceElements_numNodes) ) {
                       delete[] FaceElements_Nodes;
                       cleanupAndThrow(mesh_p, "get(FaceElements_Nodes)");
                   }
                   // Copy temp array into mesh_p->FaceElements->Nodes
                   for (int i=0; i<num_FaceElements; i++) {
                       for (int j=0; j<num_FaceElements_numNodes; j++) {
                           mesh_p->FaceElements->Nodes[INDEX2(j,i,num_FaceElements_numNodes)] = FaceElements_Nodes[INDEX2(j,i,num_FaceElements_numNodes)];
                       }
                   }
                   delete[] FaceElements_Nodes;
                   Dudley_ElementFile_setTagsInUse(mesh_p->FaceElements);
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
                   // Points_Nodes
                   int *Points_Nodes = new int[num_Points];
                   if (!(nc_var_temp = dataFile.get_var("Points_Nodes"))) {
                       delete[] Points_Nodes;
                       cleanupAndThrow(mesh_p, "get_var(Points_Nodes)");
                   }
                   if (! nc_var_temp->get(&(Points_Nodes[0]), num_Points) ) {
                       delete[] Points_Nodes;
                       cleanupAndThrow(mesh_p, "get(Points_Nodes)");
                   }
                   // Copy temp array into mesh_p->Points->Nodes
                   for (int i=0; i<num_Points; i++) {
                       mesh_p->Points->Id[mesh_p->Points->Nodes[INDEX2(0,i,1)]] = Points_Nodes[i];
                   }
                   delete[] Points_Nodes;
                   Dudley_ElementFile_setTagsInUse(mesh_p->Points);
                } /* num_Points>0 */
            }
        }

        /* get the tags */
        if (Dudley_noError()) {
          if (num_Tags>0) {
            // Temp storage to gather node IDs
            int *Tags_keys = new int[num_Tags];
            char name_temp[4096];
            int i;

            // Tags_keys
            if (! ( nc_var_temp = dataFile.get_var("Tags_keys")) ) {
                delete[] Tags_keys;
                cleanupAndThrow(mesh_p, "get_var(Tags_keys)");
            }
            if (! nc_var_temp->get(&Tags_keys[0], num_Tags) ) {
                delete[] Tags_keys;
                cleanupAndThrow(mesh_p, "get(Tags_keys)");
            }
            for (i=0; i<num_Tags; i++) {
              // Retrieve tag name
              sprintf(name_temp, "Tags_name_%d", i);
              if (! (attr=dataFile.get_att(name_temp)) ) {
                  delete[] Tags_keys;
                  sprintf(error_msg,"get_att(%s)", name_temp);
                  cleanupAndThrow(mesh_p, error_msg);
              }
              boost::scoped_array<char> name(attr->as_string(0));
              delete attr;
              Dudley_Mesh_addTagMap(mesh_p, name.get(), Tags_keys[i]);
            }
            delete[] Tags_keys;
          }
        }
   
        if (Dudley_noError()) {
            // Nodes_DofDistribution
            first_DofComponent = new index_t[mpi_size+1];
            if (! ( nc_var_temp = dataFile.get_var("Nodes_DofDistribution")) ) {
                delete[] first_DofComponent;
                cleanupAndThrow(mesh_p, "get_var(Nodes_DofDistribution)");
            }
            if (! nc_var_temp->get(&first_DofComponent[0], mpi_size+1) ) {
                delete[] first_DofComponent;
                cleanupAndThrow(mesh_p, "get(Nodes_DofDistribution)");
            }

            // Nodes_NodeDistribution
            first_NodeComponent = new index_t[mpi_size+1];
            if (! ( nc_var_temp = dataFile.get_var("Nodes_NodeDistribution")) ) {
                delete[] first_DofComponent;
                delete[] first_NodeComponent;
                cleanupAndThrow(mesh_p, "get_var(Nodes_NodeDistribution)");
            }
            if (! nc_var_temp->get(&first_NodeComponent[0], mpi_size+1) ) {
                delete[] first_DofComponent;
                delete[] first_NodeComponent;
                cleanupAndThrow(mesh_p, "get(Nodes_NodeDistribution)");
            }
            Dudley_Mesh_createMappings(mesh_p, first_DofComponent, first_NodeComponent);
            delete[] first_DofComponent;
            delete[] first_NodeComponent;
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

    char *fName = new char[fileName.size()+1];
        
    strcpy(fName,fileName.c_str());
    double blocktimer_start = blocktimer_time();

    fMesh=Dudley_Mesh_read(fName,integrationOrder, reducedIntegrationOrder, (optimize ? TRUE : FALSE));
    checkDudleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    
    delete[] fName;
    
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

    char *fName = new char[fileName.size()+1];
        
    strcpy(fName,fileName.c_str());
    double blocktimer_start = blocktimer_time();

    fMesh=Dudley_Mesh_readGmsh(fName, numDim, integrationOrder, reducedIntegrationOrder, (optimize ? TRUE : FALSE), (useMacroElements ? TRUE : FALSE));
    checkDudleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    
    delete[] fName;
    
    blocktimer_increment("ReadGmsh()", blocktimer_start);
    return temp->getPtr();
  }

  Domain_ptr brick(esysUtils::JMPI& mpi_info, double n0, double n1,double n2,int order,
                   double l0,double l1,double l2,
                   int periodic0,int periodic1,
                   int periodic2,
                   int integrationOrder,
                   int reducedIntegrationOrder,
                   int useElementsOnFace,
                   int useFullElementOrder,
                   int optimize)
  {
    int numElements[]={static_cast<int>(n0),static_cast<int>(n1),static_cast<int>(n2)};
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
                        reducedIntegrationOrder, (optimize ? TRUE : FALSE),
                        mpi_info);

    //
    // Convert any dudley errors into a C++ exception
    checkDudleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    return temp->getPtr();
  }

  Domain_ptr brick_driver(const boost::python::list& args)
  {
      using boost::python::extract;

//       // we need to convert lists to stl vectors
//       boost::python::list pypoints=extract<boost::python::list>(args[15]);
//       boost::python::list pytags=extract<boost::python::list>(args[16]);
//       int numpts=extract<int>(pypoints.attr("__len__")());
//       int numtags=extract<int>(pytags.attr("__len__")());
//       vector<double> points;
//       vector<int> tags;
//       tags.resize(numtags, -1);
//       for (int i=0;i<numpts;++i) {
//           boost::python::object temp=pypoints[i];
//           int l=extract<int>(temp.attr("__len__")());
//           for (int k=0;k<l;++k) {
//               points.push_back(extract<double>(temp[k]));           
//           }
//       }
//       map<string, int> namestonums;
//       int curmax=40; // bricks use up to 30
//       for (int i=0;i<numtags;++i) {
//           extract<int> ex_int(pytags[i]);
//           extract<string> ex_str(pytags[i]);
//           if (ex_int.check()) {
//               tags[i]=ex_int();
//               if (tags[i]>= curmax) {
//                   curmax=tags[i]+1;
//               }
//           } else if (ex_str.check()) {
//               string s=ex_str();
//               map<string, int>::iterator it=namestonums.find(s);
//               if (it!=namestonums.end()) {
//                   // we have the tag already so look it up
//                   tags[i]=it->second;
//               } else {
//                   namestonums[s]=curmax;
//                   tags[i]=curmax;
//                   curmax++;
//               }
//           } else {
//               throw DudleyAdapterException("Error - Unable to extract tag value.");
//           }
//         
//       }
      boost::python::object pworld=args[15];
      esysUtils::JMPI info;
      if (!pworld.is_none())
      {
	  extract<SubWorld_ptr> ex(pworld);
	  if (!ex.check())
	  {	  
	      throw DudleyAdapterException("Invalid escriptworld parameter.");
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
                   extract<int>(args[13]), extract<int>(args[14])
                   );
  }  
  
  
  Domain_ptr rectangle_driver(const boost::python::list& args)
  {
      using boost::python::extract;
/*
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
              throw DudleyAdapterException("Error - Unable to extract tag value.");
          }
      }*/
      boost::python::object pworld=args[12];
      esysUtils::JMPI info;
      if (!pworld.is_none())
      {
          extract<SubWorld_ptr> ex(pworld);
	  if (!ex.check())
	  {
	      throw DudleyAdapterException("Invalid escriptworld parameter.");
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
                       extract<int>(args[10]), extract<int>(args[11]) 
		       );
  }  
  
  
  
  Domain_ptr rectangle(esysUtils::JMPI& mpi_info, double n0, double n1, int order,
                       double l0, double l1,
                       int periodic0,int periodic1,
                       int integrationOrder,
                       int reducedIntegrationOrder,
                       int useElementsOnFace,
                       int useFullElementOrder,
                       int optimize)
  {
    int numElements[]={static_cast<int>(n0), static_cast<int>(n1)};
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
          integrationOrder, reducedIntegrationOrder, (optimize ? TRUE : FALSE),
          mpi_info);
    //
    // Convert any dudley errors into a C++ exception
    checkDudleyError();
    MeshAdapter* ma=new MeshAdapter(fMesh);
    return Domain_ptr(ma);
  }

  // end of namespace

}

