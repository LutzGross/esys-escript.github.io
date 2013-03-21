
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

#include <pasowrap/PasoException.h>
#include <pasowrap/TransportProblemAdapter.h>
#include "MeshAdapter.h"
#include "escript/Data.h"
#include "escript/DataFactory.h"
#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

#include "esysUtils/blocktimer.h"


#include <boost/python/import.hpp>

using namespace std;
using namespace escript;
using namespace paso;

namespace finley {

//
// define the static constants
MeshAdapter::FunctionSpaceNamesMapType MeshAdapter::m_functionSpaceTypeNames;
const int MeshAdapter::DegreesOfFreedom=FINLEY_DEGREES_OF_FREEDOM;
const int MeshAdapter::ReducedDegreesOfFreedom=FINLEY_REDUCED_DEGREES_OF_FREEDOM;
const int MeshAdapter::Nodes=FINLEY_NODES;
const int MeshAdapter::ReducedNodes=FINLEY_REDUCED_NODES;
const int MeshAdapter::Elements=FINLEY_ELEMENTS;
const int MeshAdapter::ReducedElements=FINLEY_REDUCED_ELEMENTS;
const int MeshAdapter::FaceElements=FINLEY_FACE_ELEMENTS;
const int MeshAdapter::ReducedFaceElements=FINLEY_REDUCED_FACE_ELEMENTS;
const int MeshAdapter::Points=FINLEY_POINTS;
const int MeshAdapter::ContactElementsZero=FINLEY_CONTACT_ELEMENTS_1;
const int MeshAdapter::ReducedContactElementsZero=FINLEY_REDUCED_CONTACT_ELEMENTS_1;
const int MeshAdapter::ContactElementsOne=FINLEY_CONTACT_ELEMENTS_2;
const int MeshAdapter::ReducedContactElementsOne=FINLEY_REDUCED_CONTACT_ELEMENTS_2;

MeshAdapter::MeshAdapter(Finley_Mesh* finleyMesh)
{
   setFunctionSpaceTypeNames();
   //
   // need to use a null_deleter as Finley_Mesh_free deletes the pointer
   // for us.
   m_finleyMesh.reset(finleyMesh,null_deleter());
}

//
// The copy constructor should just increment the use count
MeshAdapter::MeshAdapter(const MeshAdapter& in):
m_finleyMesh(in.m_finleyMesh)
{
   setFunctionSpaceTypeNames();
}

MeshAdapter::~MeshAdapter()
{
   //
   // I hope the case for the pointer being zero has been taken care of.
   //  cout << "In MeshAdapter destructor." << endl;
   if (m_finleyMesh.unique()) {
      Finley_Mesh_free(m_finleyMesh.get());
   }
}

int MeshAdapter::getMPISize() const
{
   return m_finleyMesh.get()->MPIInfo->size;
}
int MeshAdapter::getMPIRank() const
{
   return m_finleyMesh.get()->MPIInfo->rank;
}
void MeshAdapter::MPIBarrier() const
{
#ifdef ESYS_MPI
   MPI_Barrier(m_finleyMesh.get()->MPIInfo->comm);
#endif
   return;
}
bool MeshAdapter::onMasterProcessor() const
{
   return m_finleyMesh.get()->MPIInfo->rank == 0;
}


#ifdef ESYS_MPI
  MPI_Comm
#else
  unsigned int
#endif
MeshAdapter::getMPIComm() const
{
#ifdef ESYS_MPI
	return m_finleyMesh->MPIInfo->comm;
#else
	return 0;
#endif
}


Finley_Mesh* MeshAdapter::getFinley_Mesh() const {
   return m_finleyMesh.get();
}

void MeshAdapter::write(const string& fileName) const
{
   char *fName = (fileName.size()+1>0) ? new char[fileName.size()+1] : (char*)NULL;
   strcpy(fName,fileName.c_str());
   Finley_Mesh_write(m_finleyMesh.get(),fName);
   checkFinleyError();
   delete[] fName;
}

void MeshAdapter::Print_Mesh_Info(const bool full) const
{
   Finley_PrintMesh_Info(m_finleyMesh.get(), full);
}

void MeshAdapter::dump(const string& fileName) const
{
#ifdef USE_NETCDF
   const NcDim* ncdims[12] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
   NcVar *ids;
   int *int_ptr;
   Finley_Mesh *mesh = m_finleyMesh.get();
   Finley_TagMap* tag_map;
   int num_Tags = 0;
   int mpi_size				= mesh->MPIInfo->size;
   int mpi_rank				= mesh->MPIInfo->rank;
   int numDim				= mesh->Nodes->numDim;
   int numNodes				= mesh->Nodes->numNodes;
   int num_Elements			= mesh->Elements->numElements;
   int num_FaceElements			= mesh->FaceElements->numElements;
   int num_ContactElements		= mesh->ContactElements->numElements;
   int num_Points			= mesh->Points->numElements;
   int num_Elements_numNodes		= mesh->Elements->numNodes;
   int num_FaceElements_numNodes	= mesh->FaceElements->numNodes;
   int num_ContactElements_numNodes	= mesh->ContactElements->numNodes;
#ifdef ESYS_MPI
   MPI_Status status;
#endif

/* Incoming token indicates it's my turn to write */
#ifdef ESYS_MPI
   if (mpi_rank>0) MPI_Recv(&num_Tags, 0, MPI_INT, mpi_rank-1, 81800, mesh->MPIInfo->comm, &status);
#endif

   char *newFileName = Esys_MPI_appendRankToFileName(fileName.c_str(),
                                                     mpi_size, mpi_rank);

   /* Figure out how much storage is required for tags */
   tag_map = mesh->TagMap;
   num_Tags = 0;
   while (tag_map) {
      num_Tags++;
      tag_map=tag_map->next;
   }

   // NetCDF error handler
   NcError err(NcError::verbose_nonfatal);
   // Create the file.
   NcFile dataFile(newFileName, NcFile::Replace);
   string msgPrefix("Error in MeshAdapter::dump: NetCDF operation failed - ");
   // check if writing was successful
   if (!dataFile.is_valid())
      throw DataException(msgPrefix+"Open file for output");

   // Define dimensions (num_Elements and dim_Elements are identical,
   // dim_Elements only appears if > 0)
   if (! (ncdims[0] = dataFile.add_dim("numNodes", numNodes)) )
      throw DataException(msgPrefix+"add_dim(numNodes)");
   if (! (ncdims[1] = dataFile.add_dim("numDim", numDim)) )
      throw DataException(msgPrefix+"add_dim(numDim)");
   if (! (ncdims[2] = dataFile.add_dim("mpi_size_plus_1", mpi_size+1)) )
      throw DataException(msgPrefix+"add_dim(mpi_size)");
   if (num_Elements>0)
      if (! (ncdims[3] = dataFile.add_dim("dim_Elements", num_Elements)) )
         throw DataException(msgPrefix+"add_dim(dim_Elements)");
   if (num_FaceElements>0)
      if (! (ncdims[4] = dataFile.add_dim("dim_FaceElements", num_FaceElements)) )
         throw DataException(msgPrefix+"add_dim(dim_FaceElements)");
   if (num_ContactElements>0)
      if (! (ncdims[5] = dataFile.add_dim("dim_ContactElements", num_ContactElements)) )
         throw DataException(msgPrefix+"add_dim(dim_ContactElements)");
   if (num_Points>0)
      if (! (ncdims[6] = dataFile.add_dim("dim_Points", num_Points)) )
         throw DataException(msgPrefix+"add_dim(dim_Points)");
   if (num_Elements>0)
      if (! (ncdims[7] = dataFile.add_dim("dim_Elements_Nodes", num_Elements_numNodes)) )
         throw DataException(msgPrefix+"add_dim(dim_Elements_Nodes)");
   if (num_FaceElements>0)
      if (! (ncdims[8] = dataFile.add_dim("dim_FaceElements_numNodes", num_FaceElements_numNodes)) )
         throw DataException(msgPrefix+"add_dim(dim_FaceElements_numNodes)");
   if (num_ContactElements>0)
      if (! (ncdims[9] = dataFile.add_dim("dim_ContactElements_numNodes", num_ContactElements_numNodes)) )
         throw DataException(msgPrefix+"add_dim(dim_ContactElements_numNodes)");
   if (num_Tags>0)
      if (! (ncdims[10] = dataFile.add_dim("dim_Tags", num_Tags)) )
         throw DataException(msgPrefix+"add_dim(dim_Tags)");

   // Attributes: MPI size, MPI rank, Name, order, reduced_order
   if (!dataFile.add_att("mpi_size", mpi_size) )
      throw DataException(msgPrefix+"add_att(mpi_size)");
   if (!dataFile.add_att("mpi_rank", mpi_rank) )
      throw DataException(msgPrefix+"add_att(mpi_rank)");
   if (!dataFile.add_att("Name",mesh->Name) )
      throw DataException(msgPrefix+"add_att(Name)");
   if (!dataFile.add_att("numDim",numDim) )
      throw DataException(msgPrefix+"add_att(order)");
   if (!dataFile.add_att("order",mesh->integrationOrder) )
      throw DataException(msgPrefix+"add_att(order)");
   if (!dataFile.add_att("reduced_order",mesh->reducedIntegrationOrder) )
      throw DataException(msgPrefix+"add_att(reduced_order)");
   if (!dataFile.add_att("numNodes",numNodes) )
      throw DataException(msgPrefix+"add_att(numNodes)");
   if (!dataFile.add_att("num_Elements",num_Elements) )
      throw DataException(msgPrefix+"add_att(num_Elements)");
   if (!dataFile.add_att("num_FaceElements",num_FaceElements) )
      throw DataException(msgPrefix+"add_att(num_FaceElements)");
   if (!dataFile.add_att("num_ContactElements",num_ContactElements) )
      throw DataException(msgPrefix+"add_att(num_ContactElements)");
   if (!dataFile.add_att("num_Points",num_Points) )
      throw DataException(msgPrefix+"add_att(num_Points)");
   if (!dataFile.add_att("num_Elements_numNodes",num_Elements_numNodes) )
      throw DataException(msgPrefix+"add_att(num_Elements_numNodes)");
   if (!dataFile.add_att("num_FaceElements_numNodes",num_FaceElements_numNodes) )
      throw DataException(msgPrefix+"add_att(num_FaceElements_numNodes)");
   if (!dataFile.add_att("num_ContactElements_numNodes",num_ContactElements_numNodes) )
      throw DataException(msgPrefix+"add_att(num_ContactElements_numNodes)");
   if (!dataFile.add_att("Elements_TypeId", mesh->Elements->referenceElementSet->referenceElement->Type->TypeId) )
      throw DataException(msgPrefix+"add_att(Elements_TypeId)");
   if (!dataFile.add_att("FaceElements_TypeId", mesh->FaceElements->referenceElementSet->referenceElement->Type->TypeId) )
      throw DataException(msgPrefix+"add_att(FaceElements_TypeId)");
   if (!dataFile.add_att("ContactElements_TypeId", mesh->ContactElements->referenceElementSet->referenceElement->Type->TypeId) )
      throw DataException(msgPrefix+"add_att(ContactElements_TypeId)");
   if (!dataFile.add_att("Points_TypeId", mesh->Points->referenceElementSet->referenceElement->Type->TypeId) )
      throw DataException(msgPrefix+"add_att(Points_TypeId)");
   if (!dataFile.add_att("num_Tags", num_Tags) )
      throw DataException(msgPrefix+"add_att(num_Tags)");

   // // // // // Nodes // // // // //

   // Nodes nodeDistribution
   if (! ( ids = dataFile.add_var("Nodes_NodeDistribution", ncInt, ncdims[2])) )
      throw DataException(msgPrefix+"add_var(Nodes_NodeDistribution)");
   int_ptr = &mesh->Nodes->nodesDistribution->first_component[0];
   if (! (ids->put(int_ptr, mpi_size+1)) )
      throw DataException(msgPrefix+"put(Nodes_NodeDistribution)");

   // Nodes degreesOfFreedomDistribution
   if (! ( ids = dataFile.add_var("Nodes_DofDistribution", ncInt, ncdims[2])) )
      throw DataException(msgPrefix+"add_var(Nodes_DofDistribution)");
   int_ptr = &mesh->Nodes->degreesOfFreedomDistribution->first_component[0];
   if (! (ids->put(int_ptr, mpi_size+1)) )
      throw DataException(msgPrefix+"put(Nodes_DofDistribution)");

   // Only write nodes if non-empty because NetCDF doesn't like empty arrays
   // (it treats them as NC_UNLIMITED)
   if (numNodes>0) {

      // Nodes Id
      if (! ( ids = dataFile.add_var("Nodes_Id", ncInt, ncdims[0])) )
         throw DataException(msgPrefix+"add_var(Nodes_Id)");
      int_ptr = &mesh->Nodes->Id[0];
      if (! (ids->put(int_ptr, numNodes)) )
         throw DataException(msgPrefix+"put(Nodes_Id)");

      // Nodes Tag
      if (! ( ids = dataFile.add_var("Nodes_Tag", ncInt, ncdims[0])) )
         throw DataException(msgPrefix+"add_var(Nodes_Tag)");
      int_ptr = &mesh->Nodes->Tag[0];
      if (! (ids->put(int_ptr, numNodes)) )
         throw DataException(msgPrefix+"put(Nodes_Tag)");

      // Nodes gDOF
      if (! ( ids = dataFile.add_var("Nodes_gDOF", ncInt, ncdims[0])) )
         throw DataException(msgPrefix+"add_var(Nodes_gDOF)");
      int_ptr = &mesh->Nodes->globalDegreesOfFreedom[0];
      if (! (ids->put(int_ptr, numNodes)) )
         throw DataException(msgPrefix+"put(Nodes_gDOF)");

      // Nodes global node index
      if (! ( ids = dataFile.add_var("Nodes_gNI", ncInt, ncdims[0])) )
         throw DataException(msgPrefix+"add_var(Nodes_gNI)");
      int_ptr = &mesh->Nodes->globalNodesIndex[0];
      if (! (ids->put(int_ptr, numNodes)) )
         throw DataException(msgPrefix+"put(Nodes_gNI)");

      // Nodes grDof
      if (! ( ids = dataFile.add_var("Nodes_grDfI", ncInt, ncdims[0])) )
         throw DataException(msgPrefix+"add_var(Nodes_grDfI)");
      int_ptr = &mesh->Nodes->globalReducedDOFIndex[0];
      if (! (ids->put(int_ptr, numNodes)) )
         throw DataException(msgPrefix+"put(Nodes_grDfI)");

      // Nodes grNI
      if (! ( ids = dataFile.add_var("Nodes_grNI", ncInt, ncdims[0])) )
         throw DataException(msgPrefix+"add_var(Nodes_grNI)");
      int_ptr = &mesh->Nodes->globalReducedNodesIndex[0];
      if (! (ids->put(int_ptr, numNodes)) )
         throw DataException(msgPrefix+"put(Nodes_grNI)");

      // Nodes Coordinates
      if (! ( ids = dataFile.add_var("Nodes_Coordinates", ncDouble, ncdims[0], ncdims[1]) ) )
         throw DataException(msgPrefix+"add_var(Nodes_Coordinates)");
      if (! (ids->put(&(mesh->Nodes->Coordinates[INDEX2(0,0,numDim)]), numNodes, numDim)) )
         throw DataException(msgPrefix+"put(Nodes_Coordinates)");

   }

   // // // // // Elements // // // // //

   if (num_Elements>0) {

      // Elements_Id
      if (! ( ids = dataFile.add_var("Elements_Id", ncInt, ncdims[3])) )
         throw DataException(msgPrefix+"add_var(Elements_Id)");
      int_ptr = &mesh->Elements->Id[0];
      if (! (ids->put(int_ptr, num_Elements)) )
         throw DataException(msgPrefix+"put(Elements_Id)");

      // Elements_Tag
      if (! ( ids = dataFile.add_var("Elements_Tag", ncInt, ncdims[3])) )
         throw DataException(msgPrefix+"add_var(Elements_Tag)");
      int_ptr = &mesh->Elements->Tag[0];
      if (! (ids->put(int_ptr, num_Elements)) )
         throw DataException(msgPrefix+"put(Elements_Tag)");

      // Elements_Owner
      if (! ( ids = dataFile.add_var("Elements_Owner", ncInt, ncdims[3])) )
         throw DataException(msgPrefix+"add_var(Elements_Owner)");
      int_ptr = &mesh->Elements->Owner[0];
      if (! (ids->put(int_ptr, num_Elements)) )
         throw DataException(msgPrefix+"put(Elements_Owner)");

      // Elements_Color
      if (! ( ids = dataFile.add_var("Elements_Color", ncInt, ncdims[3])) )
         throw DataException(msgPrefix+"add_var(Elements_Color)");
      int_ptr = &mesh->Elements->Color[0];
      if (! (ids->put(int_ptr, num_Elements)) )
         throw DataException(msgPrefix+"put(Elements_Color)");

      // Elements_Nodes
      if (! ( ids = dataFile.add_var("Elements_Nodes", ncInt, ncdims[3], ncdims[7]) ) )
         throw DataException(msgPrefix+"add_var(Elements_Nodes)");
      if (! (ids->put(&(mesh->Elements->Nodes[0]), num_Elements, num_Elements_numNodes)) )
         throw DataException(msgPrefix+"put(Elements_Nodes)");

   }

   // // // // // Face_Elements // // // // //

   if (num_FaceElements>0) {

      // FaceElements_Id
      if (! ( ids = dataFile.add_var("FaceElements_Id", ncInt, ncdims[4])) )
         throw DataException(msgPrefix+"add_var(FaceElements_Id)");
      int_ptr = &mesh->FaceElements->Id[0];
      if (! (ids->put(int_ptr, num_FaceElements)) )
         throw DataException(msgPrefix+"put(FaceElements_Id)");

      // FaceElements_Tag
      if (! ( ids = dataFile.add_var("FaceElements_Tag", ncInt, ncdims[4])) )
         throw DataException(msgPrefix+"add_var(FaceElements_Tag)");
      int_ptr = &mesh->FaceElements->Tag[0];
      if (! (ids->put(int_ptr, num_FaceElements)) )
         throw DataException(msgPrefix+"put(FaceElements_Tag)");

      // FaceElements_Owner
      if (! ( ids = dataFile.add_var("FaceElements_Owner", ncInt, ncdims[4])) )
         throw DataException(msgPrefix+"add_var(FaceElements_Owner)");
      int_ptr = &mesh->FaceElements->Owner[0];
      if (! (ids->put(int_ptr, num_FaceElements)) )
         throw DataException(msgPrefix+"put(FaceElements_Owner)");

      // FaceElements_Color
      if (! ( ids = dataFile.add_var("FaceElements_Color", ncInt, ncdims[4])) )
         throw DataException(msgPrefix+"add_var(FaceElements_Color)");
      int_ptr = &mesh->FaceElements->Color[0];
      if (! (ids->put(int_ptr, num_FaceElements)) )
         throw DataException(msgPrefix+"put(FaceElements_Color)");

      // FaceElements_Nodes
      if (! ( ids = dataFile.add_var("FaceElements_Nodes", ncInt, ncdims[4], ncdims[8]) ) )
         throw DataException(msgPrefix+"add_var(FaceElements_Nodes)");
      if (! (ids->put(&(mesh->FaceElements->Nodes[0]), num_FaceElements, num_FaceElements_numNodes)) )
         throw DataException(msgPrefix+"put(FaceElements_Nodes)");

   }

   // // // // // Contact_Elements // // // // //

   if (num_ContactElements>0) {

      // ContactElements_Id
      if (! ( ids = dataFile.add_var("ContactElements_Id", ncInt, ncdims[5])) )
         throw DataException(msgPrefix+"add_var(ContactElements_Id)");
      int_ptr = &mesh->ContactElements->Id[0];
      if (! (ids->put(int_ptr, num_ContactElements)) )
         throw DataException(msgPrefix+"put(ContactElements_Id)");

      // ContactElements_Tag
      if (! ( ids = dataFile.add_var("ContactElements_Tag", ncInt, ncdims[5])) )
         throw DataException(msgPrefix+"add_var(ContactElements_Tag)");
      int_ptr = &mesh->ContactElements->Tag[0];
      if (! (ids->put(int_ptr, num_ContactElements)) )
         throw DataException(msgPrefix+"put(ContactElements_Tag)");

      // ContactElements_Owner
      if (! ( ids = dataFile.add_var("ContactElements_Owner", ncInt, ncdims[5])) )
         throw DataException(msgPrefix+"add_var(ContactElements_Owner)");
      int_ptr = &mesh->ContactElements->Owner[0];
      if (! (ids->put(int_ptr, num_ContactElements)) )
         throw DataException(msgPrefix+"put(ContactElements_Owner)");

      // ContactElements_Color
      if (! ( ids = dataFile.add_var("ContactElements_Color", ncInt, ncdims[5])) )
         throw DataException(msgPrefix+"add_var(ContactElements_Color)");
      int_ptr = &mesh->ContactElements->Color[0];
      if (! (ids->put(int_ptr, num_ContactElements)) )
         throw DataException(msgPrefix+"put(ContactElements_Color)");

      // ContactElements_Nodes
      if (! ( ids = dataFile.add_var("ContactElements_Nodes", ncInt, ncdims[5], ncdims[9]) ) )
         throw DataException(msgPrefix+"add_var(ContactElements_Nodes)");
      if (! (ids->put(&(mesh->ContactElements->Nodes[0]), num_ContactElements, num_ContactElements_numNodes)) )
         throw DataException(msgPrefix+"put(ContactElements_Nodes)");

   }

   // // // // // Points // // // // //

   if (num_Points>0) {

      fprintf(stderr, "\n\n\nWARNING: MeshAdapter::dump has not been tested with Point elements\n\n\n");

      // Points_Id
      if (! ( ids = dataFile.add_var("Points_Id", ncInt, ncdims[6])) )
         throw DataException(msgPrefix+"add_var(Points_Id)");
      int_ptr = &mesh->Points->Id[0];
      if (! (ids->put(int_ptr, num_Points)) )
         throw DataException(msgPrefix+"put(Points_Id)");

      // Points_Tag
      if (! ( ids = dataFile.add_var("Points_Tag", ncInt, ncdims[6])) )
         throw DataException(msgPrefix+"add_var(Points_Tag)");
      int_ptr = &mesh->Points->Tag[0];
      if (! (ids->put(int_ptr, num_Points)) )
         throw DataException(msgPrefix+"put(Points_Tag)");

      // Points_Owner
      if (! ( ids = dataFile.add_var("Points_Owner", ncInt, ncdims[6])) )
         throw DataException(msgPrefix+"add_var(Points_Owner)");
      int_ptr = &mesh->Points->Owner[0];
      if (! (ids->put(int_ptr, num_Points)) )
         throw DataException(msgPrefix+"put(Points_Owner)");

      // Points_Color
      if (! ( ids = dataFile.add_var("Points_Color", ncInt, ncdims[6])) )
         throw DataException(msgPrefix+"add_var(Points_Color)");
      int_ptr = &mesh->Points->Color[0];
      if (! (ids->put(int_ptr, num_Points)) )
         throw DataException(msgPrefix+"put(Points_Color)");

      // Points_Nodes
      // mesh->Nodes->Id[mesh->Points->Nodes[INDEX2(0,i,1)]]
      if (! ( ids = dataFile.add_var("Points_Nodes", ncInt, ncdims[6]) ) )
         throw DataException(msgPrefix+"add_var(Points_Nodes)");
      if (! (ids->put(&(mesh->Points->Nodes[0]), num_Points)) )
         throw DataException(msgPrefix+"put(Points_Nodes)");

   }

   // // // // // TagMap // // // // //

   if (num_Tags>0) {

      // Temp storage to gather node IDs
      int *Tags_keys = new  int[num_Tags];
      char name_temp[4096];

      /* Copy tag data into temp arrays */
      tag_map = mesh->TagMap;
      if (tag_map) {
         int i = 0;
         while (tag_map) {
            Tags_keys[i++] = tag_map->tag_key;
            tag_map=tag_map->next;
         }
      }

      // Tags_keys
      if (! ( ids = dataFile.add_var("Tags_keys", ncInt, ncdims[10])) )
         throw DataException(msgPrefix+"add_var(Tags_keys)");
      int_ptr = &Tags_keys[0];
      if (! (ids->put(int_ptr, num_Tags)) )
         throw DataException(msgPrefix+"put(Tags_keys)");

      // Tags_names_*
      // This is an array of strings, it should be stored as an array but
      // instead I have hacked in one attribute per string because the NetCDF
      // manual doesn't tell how to do an array of strings
      tag_map = mesh->TagMap;
      if (tag_map) {
         int i = 0;
         while (tag_map) {
            sprintf(name_temp, "Tags_name_%d", i);
            if (!dataFile.add_att(name_temp, tag_map->name) )
               throw DataException(msgPrefix+"add_att(Tags_names_XX)");
            tag_map=tag_map->next;
            i++;
         }
      }

      delete[] Tags_keys;
   }

/* Send token to next MPI process so he can take his turn */
#ifdef ESYS_MPI
   if (mpi_rank<mpi_size-1) MPI_Send(&num_Tags, 0, MPI_INT, mpi_rank+1, 81800, mesh->MPIInfo->comm);
#endif

   // NetCDF file is closed by destructor of NcFile object

#else
   Finley_setError(IO_ERROR, "MeshAdapter::dump: not configured with NetCDF. Please contact your installation manager.");
#endif	/* USE_NETCDF */
   checkFinleyError();
}

string MeshAdapter::getDescription() const
{
   return "FinleyMesh";
}

string MeshAdapter::functionSpaceTypeAsString(int functionSpaceType) const
{
   FunctionSpaceNamesMapType::iterator loc;
   loc=m_functionSpaceTypeNames.find(functionSpaceType);
   if (loc==m_functionSpaceTypeNames.end()) {
      return "Invalid function space type code.";
   } else {
      return loc->second;
   }
}

bool MeshAdapter::isValidFunctionSpaceType(int functionSpaceType) const
{
   FunctionSpaceNamesMapType::iterator loc;
   loc=m_functionSpaceTypeNames.find(functionSpaceType);
   return (loc!=m_functionSpaceTypeNames.end());
}

void MeshAdapter::setFunctionSpaceTypeNames()
{
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(DegreesOfFreedom,"Finley_DegreesOfFreedom"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(ReducedDegreesOfFreedom,"Finley_ReducedDegreesOfFreedom"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(Nodes,"Finley_Nodes"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(ReducedNodes,"Finley_Reduced_Nodes"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(Elements,"Finley_Elements"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(ReducedElements,"Finley_Reduced_Elements"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(FaceElements,"Finley_Face_Elements"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(ReducedFaceElements,"Finley_Reduced_Face_Elements"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(Points,"Finley_Points"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(ContactElementsZero,"Finley_Contact_Elements_0"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(ReducedContactElementsZero,"Finley_Reduced_Contact_Elements_0"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(ContactElementsOne,"Finley_Contact_Elements_1"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(ReducedContactElementsOne,"Finley_Reduced_Contact_Elements_1"));
}

int MeshAdapter::getContinuousFunctionCode() const
{
   return Nodes;
}
int MeshAdapter::getReducedContinuousFunctionCode() const
{
   return ReducedNodes;
}

int MeshAdapter::getFunctionCode() const
{
   return Elements;
}
int MeshAdapter::getReducedFunctionCode() const
{
   return ReducedElements;
}

int MeshAdapter::getFunctionOnBoundaryCode() const
{
   return FaceElements;
}
int MeshAdapter::getReducedFunctionOnBoundaryCode() const
{
   return ReducedFaceElements;
}

int MeshAdapter::getFunctionOnContactZeroCode() const
{
   return ContactElementsZero;
}
int MeshAdapter::getReducedFunctionOnContactZeroCode() const
{
   return ReducedContactElementsZero;
}

int MeshAdapter::getFunctionOnContactOneCode() const
{
   return ContactElementsOne;
}
int MeshAdapter::getReducedFunctionOnContactOneCode() const
{
   return ReducedContactElementsOne;
}

int MeshAdapter::getSolutionCode() const
{
   return DegreesOfFreedom;
}

int MeshAdapter::getReducedSolutionCode() const
{
   return ReducedDegreesOfFreedom;
}

int MeshAdapter::getDiracDeltaFunctionsCode() const
{
   return Points;
}

//
// return the spatial dimension of the Mesh:
//
int MeshAdapter::getDim() const
{
   Finley_Mesh* mesh=m_finleyMesh.get();
   int numDim=Finley_Mesh_getDim(mesh);
   checkFinleyError();
   return numDim;
}

//
// Return the number of data points summed across all MPI processes
//
int MeshAdapter::getNumDataPointsGlobal() const
{
   return Finley_NodeFile_getGlobalNumNodes(m_finleyMesh.get()->Nodes);
}

//
// return the number of data points per sample and the number of samples
// needed to represent data on a parts of the mesh.
//
pair<int,int> MeshAdapter::getDataShape(int functionSpaceCode) const
{
   int numDataPointsPerSample=0;
   int numSamples=0;
   Finley_Mesh* mesh=m_finleyMesh.get();
   switch (functionSpaceCode) {
   case(Nodes):
   numDataPointsPerSample=1;
   numSamples=Finley_NodeFile_getNumNodes(mesh->Nodes);
   break;
   case(ReducedNodes):
   numDataPointsPerSample=1;
   numSamples=Finley_NodeFile_getNumReducedNodes(mesh->Nodes);
   break;
   case(Elements):
   if (mesh->Elements!=NULL) {
      numSamples=mesh->Elements->numElements;
      numDataPointsPerSample=mesh->Elements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
   }
   break;
   case(ReducedElements):
   if (mesh->Elements!=NULL) {
      numSamples=mesh->Elements->numElements;
      numDataPointsPerSample=mesh->Elements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
   }
   break;
   case(FaceElements):
   if (mesh->FaceElements!=NULL) {
      numDataPointsPerSample=mesh->FaceElements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
      numSamples=mesh->FaceElements->numElements;
   }
   break;
   case(ReducedFaceElements):
   if (mesh->FaceElements!=NULL) {
      numDataPointsPerSample=mesh->FaceElements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
      numSamples=mesh->FaceElements->numElements;
   }
   break;
   case(Points):
   if (mesh->Points!=NULL) {
      numDataPointsPerSample=1;
      numSamples=mesh->Points->numElements;
   }
   break;
   case(ContactElementsZero):
   if (mesh->ContactElements!=NULL) {
      numDataPointsPerSample=mesh->ContactElements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
      numSamples=mesh->ContactElements->numElements;
   }
   break;
   case(ReducedContactElementsZero):
   if (mesh->ContactElements!=NULL) {
      numDataPointsPerSample=mesh->ContactElements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
      numSamples=mesh->ContactElements->numElements;
   }
   break;
   case(ContactElementsOne):
   if (mesh->ContactElements!=NULL) {
      numDataPointsPerSample=mesh->ContactElements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
      numSamples=mesh->ContactElements->numElements;
   }
   break;
   case(ReducedContactElementsOne):
   if (mesh->ContactElements!=NULL) {
      numDataPointsPerSample=mesh->ContactElements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
      numSamples=mesh->ContactElements->numElements;
   }
   break;
   case(DegreesOfFreedom):
   if (mesh->Nodes!=NULL) {
      numDataPointsPerSample=1;
      numSamples=Finley_NodeFile_getNumDegreesOfFreedom(mesh->Nodes);
   }
   break;
   case(ReducedDegreesOfFreedom):
   if (mesh->Nodes!=NULL) {
      numDataPointsPerSample=1;
      numSamples=Finley_NodeFile_getNumReducedDegreesOfFreedom(mesh->Nodes);
   }
   break;
   default:
      stringstream temp;
      temp << "Error - Invalid function space type: " << functionSpaceCode << " for domain: " << getDescription();
      throw FinleyAdapterException(temp.str());
      break;
   }
   return pair<int,int>(numDataPointsPerSample,numSamples);
}

//
// adds linear PDE of second order into a given stiffness matrix and right hand side:
//
void MeshAdapter::addPDEToSystem(
                                 AbstractSystemMatrix& mat, escript::Data& rhs,
                                 const escript::Data& A, const escript::Data& B, const escript::Data& C,const  escript::Data& D,const  escript::Data& X,const  escript::Data& Y,
                                 const escript::Data& d, const escript::Data& y, 
                                 const escript::Data& d_contact,const escript::Data& y_contact, 
                                 const escript::Data& d_dirac,const escript::Data& y_dirac) const
{
   SystemMatrixAdapter* smat=dynamic_cast<SystemMatrixAdapter*>(&mat);
   if (smat==0)
   {
	throw FinleyAdapterException("finley only supports Paso system matrices.");
   }
   escriptDataC _rhs=rhs.getDataC();
   escriptDataC _A  =A.getDataC();
   escriptDataC _B=B.getDataC();
   escriptDataC _C=C.getDataC();
   escriptDataC _D=D.getDataC();
   escriptDataC _X=X.getDataC();
   escriptDataC _Y=Y.getDataC();
   escriptDataC _d=d.getDataC();
   escriptDataC _y=y.getDataC();
   escriptDataC _d_contact=d_contact.getDataC();
   escriptDataC _y_contact=y_contact.getDataC();
   escriptDataC _d_dirac=d_dirac.getDataC();
   escriptDataC _y_dirac=y_dirac.getDataC();

   Finley_Mesh* mesh=m_finleyMesh.get();

   Finley_Assemble_PDE(mesh->Nodes,mesh->Elements,smat->getPaso_SystemMatrix(), &_rhs, &_A, &_B, &_C, &_D, &_X, &_Y );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->FaceElements, smat->getPaso_SystemMatrix(), &_rhs, 0, 0, 0, &_d, 0, &_y );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->ContactElements, smat->getPaso_SystemMatrix(), &_rhs , 0, 0, 0, &_d_contact, 0, &_y_contact );
   checkFinleyError();

    Finley_Assemble_PDE(mesh->Nodes,mesh->Points, smat->getPaso_SystemMatrix(), &_rhs , 0, 0, 0, &_d_dirac, 0, &_y_dirac );
   checkFinleyError();
}

void  MeshAdapter::addPDEToLumpedSystem(
                                        escript::Data& mat,
                                        const escript::Data& D,
                                        const escript::Data& d,
                                        const escript::Data& d_dirac,
                                        const bool useHRZ) const
{
   escriptDataC _mat=mat.getDataC();
   escriptDataC _D=D.getDataC();
   escriptDataC _d=d.getDataC();
   escriptDataC _d_dirac=d_dirac.getDataC();

   Finley_Mesh* mesh=m_finleyMesh.get();

   Finley_Assemble_LumpedSystem(mesh->Nodes,mesh->Elements,&_mat, &_D, useHRZ);
   checkFinleyError();
   
   Finley_Assemble_LumpedSystem(mesh->Nodes,mesh->FaceElements,&_mat, &_d, useHRZ);
   checkFinleyError();

   Finley_Assemble_LumpedSystem(mesh->Nodes,mesh->Points,&_mat, &_d_dirac, useHRZ);
   checkFinleyError();

}


//
// adds linear PDE of second order into the right hand side only
//
void MeshAdapter::addPDEToRHS( escript::Data& rhs, const  escript::Data& X,const  escript::Data& Y, const escript::Data& y, const escript::Data& y_contact,  const escript::Data& y_dirac) const
{
   Finley_Mesh* mesh=m_finleyMesh.get();

   escriptDataC _rhs=rhs.getDataC();
   escriptDataC _X=X.getDataC();
   escriptDataC _Y=Y.getDataC();
   escriptDataC _y=y.getDataC();
   escriptDataC _y_contact=y_contact.getDataC();
   escriptDataC _y_dirac=y_dirac.getDataC();

   Finley_Assemble_PDE(mesh->Nodes,mesh->Elements, 0, &_rhs, 0, 0, 0, 0, &_X, &_Y );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->FaceElements, 0, &_rhs, 0, 0, 0, 0, 0, &_y );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->ContactElements, 0, &_rhs , 0, 0, 0, 0, 0, &_y_contact );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->Points, 0, &_rhs , 0, 0, 0, 0, 0, &_y_dirac );
   checkFinleyError();

}
//
// adds PDE of second order into a transport problem
//
void MeshAdapter::addPDEToTransportProblem(
                                           AbstractTransportProblem& tp, escript::Data& source, const escript::Data& M,
                                           const escript::Data& A, const escript::Data& B, const escript::Data& C,
                                           const  escript::Data& D,const  escript::Data& X,const  escript::Data& Y,
                                           const escript::Data& d, const escript::Data& y, 
                                           const escript::Data& d_contact,const escript::Data& y_contact,
                                           const escript::Data& d_dirac, const escript::Data& y_dirac) const
{
   TransportProblemAdapter* tpa=dynamic_cast<TransportProblemAdapter*>(&tp);
   if (tpa==0)
   {
	throw FinleyAdapterException("finley only supports Paso transport problems.");
   }


   DataTypes::ShapeType shape;
   source.expand();
   escriptDataC _source=source.getDataC();
   escriptDataC _M=M.getDataC();
   escriptDataC _A=A.getDataC();
   escriptDataC _B=B.getDataC();
   escriptDataC _C=C.getDataC();
   escriptDataC _D=D.getDataC();
   escriptDataC _X=X.getDataC();
   escriptDataC _Y=Y.getDataC();
   escriptDataC _d=d.getDataC();
   escriptDataC _y=y.getDataC();
   escriptDataC _d_contact=d_contact.getDataC();
   escriptDataC _y_contact=y_contact.getDataC();
   escriptDataC _d_dirac=d_dirac.getDataC();
   escriptDataC _y_dirac=y_dirac.getDataC();


   Finley_Mesh* mesh=m_finleyMesh.get();
   Paso_TransportProblem* _tp = tpa->getPaso_TransportProblem();

   Finley_Assemble_PDE(mesh->Nodes,mesh->Elements,_tp->mass_matrix, &_source, 0, 0, 0, &_M, 0, 0 );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->Elements,_tp->transport_matrix, &_source, &_A, &_B, &_C, &_D, &_X, &_Y );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->FaceElements, _tp->transport_matrix, &_source, 0, 0, 0, &_d, 0, &_y );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->ContactElements, _tp->transport_matrix, &_source , 0, 0, 0, &_d_contact, 0, &_y_contact );
   checkFinleyError();
 
   Finley_Assemble_PDE(mesh->Nodes,mesh->Points, _tp->transport_matrix, &_source , 0, 0, 0, &_d_dirac, 0, &_y_dirac );
   checkFinleyError();

}

//
// interpolates data between different function spaces:
//
void MeshAdapter::interpolateOnDomain(escript::Data& target,const escript::Data& in) const
{
   const MeshAdapter& inDomain=dynamic_cast<const MeshAdapter&>(*(in.getFunctionSpace().getDomain()));
   const MeshAdapter& targetDomain=dynamic_cast<const MeshAdapter&>(*(target.getFunctionSpace().getDomain()));
   if (inDomain!=*this)  
      throw FinleyAdapterException("Error - Illegal domain of interpolant.");
   if (targetDomain!=*this) 
      throw FinleyAdapterException("Error - Illegal domain of interpolation target.");

   Finley_Mesh* mesh=m_finleyMesh.get();
   escriptDataC _target=target.getDataC();
   escriptDataC _in=in.getDataC();
   switch(in.getFunctionSpace().getTypeCode()) {
   case(Nodes):
      switch(target.getFunctionSpace().getTypeCode()) {
      case(Nodes):
      case(ReducedNodes):
      case(DegreesOfFreedom):
      case(ReducedDegreesOfFreedom):
      Finley_Assemble_CopyNodalData(mesh->Nodes,&_target,&_in);
      break;
      case(Elements):
      case(ReducedElements):
      Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&_in,&_target);
      break;
      case(FaceElements):
      case(ReducedFaceElements):
      Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&_in,&_target);
      break;
      case(Points):
      Finley_Assemble_interpolate(mesh->Nodes,mesh->Points,&_in,&_target);
      break;
      case(ContactElementsZero):
      case(ReducedContactElementsZero):
      Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&_in,&_target);
      break;
      case(ContactElementsOne):
      case(ReducedContactElementsOne):
      Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&_in,&_target);
      break;
      default:
         stringstream temp;
         temp << "Error - Interpolation on Domain: Finley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
         throw FinleyAdapterException(temp.str());
         break;
      }
      break;
   case(ReducedNodes):
      switch(target.getFunctionSpace().getTypeCode()) {
      case(Nodes):
      case(ReducedNodes):
      case(DegreesOfFreedom):
      case(ReducedDegreesOfFreedom):
      Finley_Assemble_CopyNodalData(mesh->Nodes,&_target,&_in);
      break;
      case(Elements):
      case(ReducedElements):
      Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&_in,&_target);
      break;
      case(FaceElements):
      case(ReducedFaceElements):
      Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&_in,&_target);
      break;
      case(Points):
      Finley_Assemble_interpolate(mesh->Nodes,mesh->Points,&_in,&_target);
      break;
      case(ContactElementsZero):
      case(ReducedContactElementsZero):
      Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&_in,&_target);
      break;
      case(ContactElementsOne):
      case(ReducedContactElementsOne):
      Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&_in,&_target);
      break;
      default:
         stringstream temp;
         temp << "Error - Interpolation on Domain: Finley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
         throw FinleyAdapterException(temp.str());
         break;
      }
      break;
   case(Elements):
      if (target.getFunctionSpace().getTypeCode()==Elements) {
         Finley_Assemble_CopyElementData(mesh->Elements,&_target,&_in);
      } else if (target.getFunctionSpace().getTypeCode()==ReducedElements) {
         Finley_Assemble_AverageElementData(mesh->Elements,&_target,&_in);
      } else {
         throw FinleyAdapterException("Error - No interpolation with data on elements possible.");
      }
      break;
   case(ReducedElements):
      if (target.getFunctionSpace().getTypeCode()==ReducedElements) {
         Finley_Assemble_CopyElementData(mesh->Elements,&_target,&_in);
      } else {
         throw FinleyAdapterException("Error - No interpolation with data on elements with reduced integration order possible.");
      }
      break;
   case(FaceElements):
      if (target.getFunctionSpace().getTypeCode()==FaceElements) {
         Finley_Assemble_CopyElementData(mesh->FaceElements,&_target,&_in);
      } else if (target.getFunctionSpace().getTypeCode()==ReducedFaceElements) {
         Finley_Assemble_AverageElementData(mesh->FaceElements,&_target,&_in);
      } else {
         throw FinleyAdapterException("Error - No interpolation with data on face elements possible.");
      }
      break;
   case(ReducedFaceElements):
      if (target.getFunctionSpace().getTypeCode()==ReducedFaceElements) {
         Finley_Assemble_CopyElementData(mesh->FaceElements,&_target,&_in);
      } else {
         throw FinleyAdapterException("Error - No interpolation with data on face elements with reduced integration order possible.");
      }
      break;
   case(Points):
      if (target.getFunctionSpace().getTypeCode()==Points) {
         Finley_Assemble_CopyElementData(mesh->Points,&_target,&_in);
      } else {
         throw FinleyAdapterException("Error - No interpolation with data on points possible.");
      }
      break;
   case(ContactElementsZero):
   case(ContactElementsOne):
      if (target.getFunctionSpace().getTypeCode()==ContactElementsZero || target.getFunctionSpace().getTypeCode()==ContactElementsOne) {
         Finley_Assemble_CopyElementData(mesh->ContactElements,&_target,&_in);
      } else if (target.getFunctionSpace().getTypeCode()==ReducedContactElementsZero || target.getFunctionSpace().getTypeCode()==ReducedContactElementsOne) {
         Finley_Assemble_AverageElementData(mesh->ContactElements,&_target,&_in);
      } else {
         throw FinleyAdapterException("Error - No interpolation with data on contact elements possible.");
      }
      break;
   case(ReducedContactElementsZero):
   case(ReducedContactElementsOne):
      if (target.getFunctionSpace().getTypeCode()==ReducedContactElementsZero || target.getFunctionSpace().getTypeCode()==ReducedContactElementsOne) {
         Finley_Assemble_CopyElementData(mesh->ContactElements,&_target,&_in);
      } else {
         throw FinleyAdapterException("Error - No interpolation with data on contact elements with reduced integration order possible.");
      }
      break;
   case(DegreesOfFreedom):      
      switch(target.getFunctionSpace().getTypeCode()) {
      case(ReducedDegreesOfFreedom):
      case(DegreesOfFreedom):
      Finley_Assemble_CopyNodalData(mesh->Nodes,&_target,&_in);
      break;
   
      case(Nodes):
      case(ReducedNodes):
      if (getMPISize()>1) {
         escript::Data temp=escript::Data(in);
         temp.expand();
         escriptDataC _in2 = temp.getDataC();
         Finley_Assemble_CopyNodalData(mesh->Nodes,&_target,&_in2);
      } else {
         Finley_Assemble_CopyNodalData(mesh->Nodes,&_target,&_in);
      }
      break;
      case(Elements):
      case(ReducedElements):
      if (getMPISize()>1) {
         escript::Data temp=escript::Data( in,  continuousFunction(*this) );
         escriptDataC _in2 = temp.getDataC();
         Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&_in2,&_target);
      } else {
         Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&_in,&_target);
      }
      break;
      case(FaceElements):
      case(ReducedFaceElements):
      if (getMPISize()>1) {
         escript::Data temp=escript::Data( in,  continuousFunction(*this) );
         escriptDataC _in2 = temp.getDataC();
         Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&_in2,&_target);
   
      } else {
         Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&_in,&_target);
      }
      break;
      case(Points):
      if (getMPISize()>1) {
         //escript::Data temp=escript::Data( in,  continuousFunction(*this) );
         //escriptDataC _in2 = temp.getDataC();
      } else {
         Finley_Assemble_interpolate(mesh->Nodes,mesh->Points,&_in,&_target);
      }
      break;
      case(ContactElementsZero):
      case(ContactElementsOne):
      case(ReducedContactElementsZero):
      case(ReducedContactElementsOne):
      if (getMPISize()>1) {
         escript::Data temp=escript::Data( in,  continuousFunction(*this) );
         escriptDataC _in2 = temp.getDataC();
         Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&_in2,&_target);
      } else {
         Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&_in,&_target);
      }
      break;
      default:
         stringstream temp;
         temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
         throw FinleyAdapterException(temp.str());
         break;
      }
      break;
   case(ReducedDegreesOfFreedom):
      switch(target.getFunctionSpace().getTypeCode()) {
      case(Nodes):
      throw FinleyAdapterException("Error - Finley does not support interpolation from reduced degrees of freedom to mesh nodes.");
      break;
      case(ReducedNodes):
      if (getMPISize()>1) {
         escript::Data temp=escript::Data(in);
         temp.expand();
         escriptDataC _in2 = temp.getDataC();
         Finley_Assemble_CopyNodalData(mesh->Nodes,&_target,&_in2);
      } else {
         Finley_Assemble_CopyNodalData(mesh->Nodes,&_target,&_in);
      }
      break;
      case(DegreesOfFreedom):
      throw FinleyAdapterException("Error - Finley does not support interpolation from reduced degrees of freedom to degrees of freedom");
      break;
      case(ReducedDegreesOfFreedom):
      Finley_Assemble_CopyNodalData(mesh->Nodes,&_target,&_in);
      break;
      case(Elements):
      case(ReducedElements):
      if (getMPISize()>1) {
         escript::Data temp=escript::Data( in,  reducedContinuousFunction(*this) );
         escriptDataC _in2 = temp.getDataC();
         Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&_in2,&_target);
      } else {
         Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&_in,&_target);
      }
      break;
      case(FaceElements):
      case(ReducedFaceElements):
      if (getMPISize()>1) {
         escript::Data temp=escript::Data( in,  reducedContinuousFunction(*this) );
         escriptDataC _in2 = temp.getDataC();
         Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&_in2,&_target);
      } else {
         Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&_in,&_target);
      }
      break;
      case(Points):
      if (getMPISize()>1) {
         escript::Data temp=escript::Data( in,  reducedContinuousFunction(*this) );
         escriptDataC _in2 = temp.getDataC();
         Finley_Assemble_interpolate(mesh->Nodes,mesh->Points,&_in2,&_target);
      } else {
         Finley_Assemble_interpolate(mesh->Nodes,mesh->Points,&_in,&_target);
      }
      break;
      case(ContactElementsZero):
      case(ContactElementsOne):
      case(ReducedContactElementsZero):
      case(ReducedContactElementsOne):
      if (getMPISize()>1) {
         escript::Data temp=escript::Data( in,  reducedContinuousFunction(*this) );
         escriptDataC _in2 = temp.getDataC();
         Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&_in2,&_target);
      } else {
         Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&_in,&_target);
      }
      break;
      default:
         stringstream temp;
         temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
         throw FinleyAdapterException(temp.str());
         break;
      }
      break;
   default:
      stringstream temp;
      temp << "Error - Interpolation On Domain: Finley does not know anything about function space type %d" << in.getFunctionSpace().getTypeCode();
      throw FinleyAdapterException(temp.str());
      break;
   }
   checkFinleyError();
}

//
// copies the locations of sample points into x:
//
void MeshAdapter::setToX(escript::Data& arg) const
{
   const MeshAdapter& argDomain=dynamic_cast<const MeshAdapter&>(*(arg.getFunctionSpace().getDomain()));
   if (argDomain!=*this) 
      throw FinleyAdapterException("Error - Illegal domain of data point locations");
   Finley_Mesh* mesh=m_finleyMesh.get();
   // in case of values node coordinates we can do the job directly:
   if (arg.getFunctionSpace().getTypeCode()==Nodes) {
      escriptDataC _arg=arg.getDataC();
      Finley_Assemble_NodeCoordinates(mesh->Nodes,&_arg);
   } else {
      escript::Data tmp_data=Vector(0.0,continuousFunction(*this),true);
      escriptDataC _tmp_data=tmp_data.getDataC();
      Finley_Assemble_NodeCoordinates(mesh->Nodes,&_tmp_data);
      // this is then interpolated onto arg:
      interpolateOnDomain(arg,tmp_data);
   }
   checkFinleyError();
}

//
// return the normal vectors at the location of data points as a Data object:
//
void MeshAdapter::setToNormal(escript::Data& normal) const
{
/*   const MeshAdapter& normalDomain=dynamic_cast<const MeshAdapter&>(normal.getFunctionSpace().getDomain());*/
   const MeshAdapter& normalDomain=dynamic_cast<const MeshAdapter&>(*(normal.getFunctionSpace().getDomain()));
   if (normalDomain!=*this) 
      throw FinleyAdapterException("Error - Illegal domain of normal locations");
   Finley_Mesh* mesh=m_finleyMesh.get();
   escriptDataC _normal=normal.getDataC();
   switch(normal.getFunctionSpace().getTypeCode()) {
   case(Nodes):
   throw FinleyAdapterException("Error - Finley does not support surface normal vectors for nodes");
   break;
   case(ReducedNodes):
   throw FinleyAdapterException("Error - Finley does not support surface normal vectors for reduced nodes");
   break;
   case(Elements):
   throw FinleyAdapterException("Error - Finley does not support surface normal vectors for elements");
   break;
   case(ReducedElements):
   throw FinleyAdapterException("Error - Finley does not support surface normal vectors for elements with reduced integration order");
   break;
   case (FaceElements):
   Finley_Assemble_setNormal(mesh->Nodes,mesh->FaceElements,&_normal);
   break;
   case (ReducedFaceElements):
   Finley_Assemble_setNormal(mesh->Nodes,mesh->FaceElements,&_normal);
   break;
   case(Points):
   throw FinleyAdapterException("Error - Finley does not support surface normal vectors for point elements");
   break;
   case (ContactElementsOne):
   case (ContactElementsZero):
   Finley_Assemble_setNormal(mesh->Nodes,mesh->ContactElements,&_normal);
   break;
   case (ReducedContactElementsOne):
   case (ReducedContactElementsZero):
   Finley_Assemble_setNormal(mesh->Nodes,mesh->ContactElements,&_normal);
   break;
   case(DegreesOfFreedom):
   throw FinleyAdapterException("Error - Finley does not support surface normal vectors for degrees of freedom.");
   break;
   case(ReducedDegreesOfFreedom):
   throw FinleyAdapterException("Error - Finley does not support surface normal vectors for reduced degrees of freedom.");
   break;
   default:
      stringstream temp;
      temp << "Error - Normal Vectors: Finley does not know anything about function space type " << normal.getFunctionSpace().getTypeCode();
      throw FinleyAdapterException(temp.str());
      break;
   }
   checkFinleyError();
}

//
// interpolates data to other domain:
//
void MeshAdapter::interpolateACross(escript::Data& target,const escript::Data& source) const
{
   const_Domain_ptr targetDomain_p=target.getFunctionSpace().getDomain();
   const MeshAdapter* targetDomain=dynamic_cast<const MeshAdapter*>(targetDomain_p.get());
   if (targetDomain!=this) 
      throw FinleyAdapterException("Error - Illegal domain of interpolation target");

   throw FinleyAdapterException("Error - Finley does not allow interpolation across domains yet.");
}

//
// calculates the integral of a function defined of arg:
//
void MeshAdapter::setToIntegrals(vector<double>& integrals,const escript::Data& arg) const
{
   const MeshAdapter& argDomain=dynamic_cast<const MeshAdapter&>(*(arg.getFunctionSpace().getDomain()));
   if (argDomain!=*this) 
      throw FinleyAdapterException("Error - Illegal domain of integration kernel");

   double blocktimer_start = blocktimer_time();
   Finley_Mesh* mesh=m_finleyMesh.get();
   escriptDataC _temp;
   escript::Data temp;
   escriptDataC _arg=arg.getDataC();
   switch(arg.getFunctionSpace().getTypeCode()) {
   case(Nodes):
   temp=escript::Data( arg, escript::function(*this) );
   _temp=temp.getDataC();
   Finley_Assemble_integrate(mesh->Nodes,mesh->Elements,&_temp,&integrals[0]);
   break;
   case(ReducedNodes):
   temp=escript::Data( arg, escript::function(*this) );
   _temp=temp.getDataC();
   Finley_Assemble_integrate(mesh->Nodes,mesh->Elements,&_temp,&integrals[0]);
   break;
   case(Elements):
   Finley_Assemble_integrate(mesh->Nodes,mesh->Elements,&_arg,&integrals[0]);
   break;
   case(ReducedElements):
   Finley_Assemble_integrate(mesh->Nodes,mesh->Elements,&_arg,&integrals[0]);
   break;
   case(FaceElements):
   Finley_Assemble_integrate(mesh->Nodes,mesh->FaceElements,&_arg,&integrals[0]);
   break;
   case(ReducedFaceElements):
   Finley_Assemble_integrate(mesh->Nodes,mesh->FaceElements,&_arg,&integrals[0]);
   break;
   case(Points):
   throw FinleyAdapterException("Error - Integral of data on points is not supported.");
   break;
   case(ContactElementsZero):
   Finley_Assemble_integrate(mesh->Nodes,mesh->ContactElements,&_arg,&integrals[0]);
   break;
   case(ReducedContactElementsZero):
   Finley_Assemble_integrate(mesh->Nodes,mesh->ContactElements,&_arg,&integrals[0]);
   break;
   case(ContactElementsOne):
   Finley_Assemble_integrate(mesh->Nodes,mesh->ContactElements,&_arg,&integrals[0]);
   break;
   case(ReducedContactElementsOne):
   Finley_Assemble_integrate(mesh->Nodes,mesh->ContactElements,&_arg,&integrals[0]);
   break;
   case(DegreesOfFreedom):
   temp=escript::Data( arg, escript::function(*this) );
   _temp=temp.getDataC();
   Finley_Assemble_integrate(mesh->Nodes,mesh->Elements,&_temp,&integrals[0]);
   break;
   case(ReducedDegreesOfFreedom):
   temp=escript::Data( arg, escript::function(*this) );
   _temp=temp.getDataC();
   Finley_Assemble_integrate(mesh->Nodes,mesh->Elements,&_temp,&integrals[0]);
   break;
   default:
      stringstream temp;
      temp << "Error - Integrals: Finley does not know anything about function space type " << arg.getFunctionSpace().getTypeCode();
      throw FinleyAdapterException(temp.str());
      break;
   }
   checkFinleyError();
   blocktimer_increment("integrate()", blocktimer_start);
}

//
// calculates the gradient of arg:
//
void MeshAdapter::setToGradient(escript::Data& grad,const escript::Data& arg) const
{
   const MeshAdapter& argDomain=dynamic_cast<const MeshAdapter&>(*(arg.getFunctionSpace().getDomain()));
   if (argDomain!=*this)
      throw FinleyAdapterException("Error - Illegal domain of gradient argument");
   const MeshAdapter& gradDomain=dynamic_cast<const MeshAdapter&>(*(grad.getFunctionSpace().getDomain()));
   if (gradDomain!=*this)
      throw FinleyAdapterException("Error - Illegal domain of gradient");

   Finley_Mesh* mesh=m_finleyMesh.get();
   escriptDataC _grad=grad.getDataC();
   escriptDataC nodeDataC;
   escript::Data temp;
   if (getMPISize()>1) {
      if( arg.getFunctionSpace().getTypeCode() == DegreesOfFreedom ) {
         temp=escript::Data( arg,  continuousFunction(*this) );
         nodeDataC = temp.getDataC();
      } else if( arg.getFunctionSpace().getTypeCode() == ReducedDegreesOfFreedom ) {
         temp=escript::Data( arg,  reducedContinuousFunction(*this) );
         nodeDataC = temp.getDataC();
      } else {
         nodeDataC = arg.getDataC();
      }
   } else {
      nodeDataC = arg.getDataC();
   }
   switch(grad.getFunctionSpace().getTypeCode()) {
   case(Nodes):
   throw FinleyAdapterException("Error - Gradient at nodes is not supported.");
   break;
   case(ReducedNodes):
   throw FinleyAdapterException("Error - Gradient at reduced nodes is not supported.");
   break;
   case(Elements):
   Finley_Assemble_gradient(mesh->Nodes,mesh->Elements,&_grad,&nodeDataC);
   break;
   case(ReducedElements):
   Finley_Assemble_gradient(mesh->Nodes,mesh->Elements,&_grad,&nodeDataC);
   break;
   case(FaceElements):
   Finley_Assemble_gradient(mesh->Nodes,mesh->FaceElements,&_grad,&nodeDataC);
   break;
   case(ReducedFaceElements):
   Finley_Assemble_gradient(mesh->Nodes,mesh->FaceElements,&_grad,&nodeDataC);
   break;
   case(Points):
   throw FinleyAdapterException("Error - Gradient at points is not supported.");
   break;
   case(ContactElementsZero):
   Finley_Assemble_gradient(mesh->Nodes,mesh->ContactElements,&_grad,&nodeDataC);
   break;
   case(ReducedContactElementsZero):
   Finley_Assemble_gradient(mesh->Nodes,mesh->ContactElements,&_grad,&nodeDataC);
   break;
   case(ContactElementsOne):
   Finley_Assemble_gradient(mesh->Nodes,mesh->ContactElements,&_grad,&nodeDataC);
   break;
   case(ReducedContactElementsOne):
   Finley_Assemble_gradient(mesh->Nodes,mesh->ContactElements,&_grad,&nodeDataC);
   break;
   case(DegreesOfFreedom):
   throw FinleyAdapterException("Error - Gradient at degrees of freedom is not supported.");
   break;
   case(ReducedDegreesOfFreedom):
   throw FinleyAdapterException("Error - Gradient at reduced degrees of freedom is not supported.");
   break;
   default:
      stringstream temp;
      temp << "Error - Gradient: Finley does not know anything about function space type " << arg.getFunctionSpace().getTypeCode();
      throw FinleyAdapterException(temp.str());
      break;
   }
   checkFinleyError();
}

//
// returns the size of elements:
//
void MeshAdapter::setToSize(escript::Data& size) const
{
   Finley_Mesh* mesh=m_finleyMesh.get();
   escriptDataC tmp=size.getDataC();
   switch(size.getFunctionSpace().getTypeCode()) {
   case(Nodes):
   throw FinleyAdapterException("Error - Size of nodes is not supported.");
   break;
   case(ReducedNodes):
   throw FinleyAdapterException("Error - Size of reduced nodes is not supported.");
   break;
   case(Elements):
   Finley_Assemble_getSize(mesh->Nodes,mesh->Elements,&tmp);
   break;
   case(ReducedElements):
   Finley_Assemble_getSize(mesh->Nodes,mesh->Elements,&tmp);
   break;
   case(FaceElements):
   Finley_Assemble_getSize(mesh->Nodes,mesh->FaceElements,&tmp);
   break;
   case(ReducedFaceElements):
   Finley_Assemble_getSize(mesh->Nodes,mesh->FaceElements,&tmp);
   break;
   case(Points):
   throw FinleyAdapterException("Error - Size of point elements is not supported.");
   break;
   case(ContactElementsZero):
   case(ContactElementsOne):
   Finley_Assemble_getSize(mesh->Nodes,mesh->ContactElements,&tmp);
   break;
   case(ReducedContactElementsZero):
   case(ReducedContactElementsOne):
   Finley_Assemble_getSize(mesh->Nodes,mesh->ContactElements,&tmp);
   break;
   case(DegreesOfFreedom):
   throw FinleyAdapterException("Error - Size of degrees of freedom is not supported.");
   break;
   case(ReducedDegreesOfFreedom):
   throw FinleyAdapterException("Error - Size of reduced degrees of freedom is not supported.");
   break;
   default:
      stringstream temp;
      temp << "Error - Element size: Finley does not know anything about function space type " << size.getFunctionSpace().getTypeCode();
      throw FinleyAdapterException(temp.str());
      break;
   }
   checkFinleyError();
}

//
// sets the location of nodes
//
void MeshAdapter::setNewX(const escript::Data& new_x)
{
   Finley_Mesh* mesh=m_finleyMesh.get();
   escriptDataC tmp;
   const MeshAdapter& newDomain=dynamic_cast<const MeshAdapter&>(*(new_x.getFunctionSpace().getDomain()));
   if (newDomain!=*this) 
      throw FinleyAdapterException("Error - Illegal domain of new point locations");
   if ( new_x.getFunctionSpace() == continuousFunction(*this) ) {
       tmp = new_x.getDataC();
       Finley_Mesh_setCoordinates(mesh,&tmp);
   } else {
       throw FinleyAdapterException("As of version escript3.3 - SetNewX only accepts ContinuousFunction arguments please interpolate.");     
/*       escript::Data new_x_inter=escript::Data( new_x,  continuousFunction(*this) );
       tmp = new_x_inter.getDataC();
       Finley_Mesh_setCoordinates(mesh,&tmp);*/
   }
   checkFinleyError();
}

bool MeshAdapter::ownSample(int fs_code, index_t id) const
{
    if (getMPISize() > 1) {
#ifdef ESYS_MPI
        index_t myFirstNode=0, myLastNode=0, k=0;
        index_t* globalNodeIndex=0;
        Finley_Mesh* mesh_p=m_finleyMesh.get();
        /*
         * this method is only used by saveDataCSV which would use the returned
         * values for reduced nodes wrongly so this case is disabled for now
        if (fs_code == FINLEY_REDUCED_NODES) 
        {
            myFirstNode = Finley_NodeFile_getFirstReducedNode(mesh_p->Nodes);
            myLastNode = Finley_NodeFile_getLastReducedNode(mesh_p->Nodes);
            globalNodeIndex = Finley_NodeFile_borrowGlobalReducedNodesIndex(mesh_p->Nodes);
        }
        else
        */
        if (fs_code == FINLEY_NODES)
        {
            myFirstNode = Finley_NodeFile_getFirstNode(mesh_p->Nodes);
            myLastNode = Finley_NodeFile_getLastNode(mesh_p->Nodes);
            globalNodeIndex = Finley_NodeFile_borrowGlobalNodesIndex(mesh_p->Nodes);
        }
        else
        {
            throw FinleyAdapterException("Unsupported function space type for ownSample()");
        }

        k=globalNodeIndex[id];
        return static_cast<bool>( (myFirstNode <= k) && (k < myLastNode) );
#endif
    }
    return true;
}


//
// creates a SystemMatrixAdapter stiffness matrix an initializes it with zeros
//
ASM_ptr MeshAdapter::newSystemMatrix(
                                                 const int row_blocksize,
                                                 const escript::FunctionSpace& row_functionspace,
                                                 const int column_blocksize,
                                                 const escript::FunctionSpace& column_functionspace,
                                                 const int type) const
{
   int reduceRowOrder=0;
   int reduceColOrder=0;
   // is the domain right?
   const MeshAdapter& row_domain=dynamic_cast<const MeshAdapter&>(*(row_functionspace.getDomain()));
   if (row_domain!=*this) 
      throw FinleyAdapterException("Error - domain of row function space does not match the domain of matrix generator.");
   const MeshAdapter& col_domain=dynamic_cast<const MeshAdapter&>(*(column_functionspace.getDomain()));
   if (col_domain!=*this) 
      throw FinleyAdapterException("Error - domain of columnn function space does not match the domain of matrix generator.");
   // is the function space type right 
   if (row_functionspace.getTypeCode()==DegreesOfFreedom) {
      reduceRowOrder=0;
   } else if (row_functionspace.getTypeCode()==ReducedDegreesOfFreedom) {
      reduceRowOrder=1;
   } else {
      throw FinleyAdapterException("Error - illegal function space type for system matrix rows.");
   }
   if (column_functionspace.getTypeCode()==DegreesOfFreedom) {
      reduceColOrder=0;
   } else if (column_functionspace.getTypeCode()==ReducedDegreesOfFreedom) {
      reduceColOrder=1;
   } else {
      throw FinleyAdapterException("Error - illegal function space type for system matrix columns.");
   }
   // generate matrix:
 
   Paso_SystemMatrixPattern* fsystemMatrixPattern=Finley_getPattern(getFinley_Mesh(),reduceRowOrder,reduceColOrder);
   checkFinleyError();
   Paso_SystemMatrix* fsystemMatrix;
   int trilinos = 0;
   if (trilinos) {
#ifdef TRILINOS
      /* Allocation Epetra_VrbMatrix here */
#endif
   }
   else {
      fsystemMatrix=Paso_SystemMatrix_alloc(type,fsystemMatrixPattern,row_blocksize,column_blocksize,FALSE);
   }
   checkPasoError();
   Paso_SystemMatrixPattern_free(fsystemMatrixPattern);
   SystemMatrixAdapter* sma=new SystemMatrixAdapter(fsystemMatrix, row_blocksize, row_functionspace, column_blocksize, column_functionspace);
   return ASM_ptr(sma);
//   return SystemMatrixAdapter(fsystemMatrix,row_blocksize,row_functionspace,column_blocksize,column_functionspace);
}

//
// creates a TransportProblemAdapter
//
ATP_ptr MeshAdapter::newTransportProblem(
                                                         const int blocksize,
                                                         const escript::FunctionSpace& functionspace,
                                                         const int type) const
{
   int reduceOrder=0;
   // is the domain right?
   const MeshAdapter& domain=dynamic_cast<const MeshAdapter&>(*(functionspace.getDomain()));
   if (domain!=*this) 
      throw FinleyAdapterException("Error - domain of function space does not match the domain of transport problem generator.");
   // is the function space type right 
   if (functionspace.getTypeCode()==DegreesOfFreedom) {
      reduceOrder=0;
   } else if (functionspace.getTypeCode()==ReducedDegreesOfFreedom) {
      reduceOrder=1;
   } else {
      throw FinleyAdapterException("Error - illegal function space type for system matrix rows.");
   }
   // generate matrix:
 
   Paso_SystemMatrixPattern* fsystemMatrixPattern=Finley_getPattern(getFinley_Mesh(),reduceOrder,reduceOrder);
   checkFinleyError();
   Paso_TransportProblem* transportProblem;
   transportProblem=Paso_TransportProblem_alloc(fsystemMatrixPattern,blocksize);
   checkPasoError();
   Paso_SystemMatrixPattern_free(fsystemMatrixPattern);
   TransportProblemAdapter* tpa=new TransportProblemAdapter(transportProblem,blocksize,functionspace);
   return ATP_ptr(tpa);
//   return TransportProblemAdapter(transportProblem, blocksize,functionspace);
}

//
// vtkObject MeshAdapter::createVtkObject() const
// TODO:
//

//
// returns true if data at the atom_type is considered as being cell centered:
bool MeshAdapter::isCellOriented(int functionSpaceCode) const
{
   switch(functionSpaceCode) {
   case(Nodes):
   case(DegreesOfFreedom):
   case(ReducedDegreesOfFreedom):
   return false;
   break;
   case(Elements):
   case(FaceElements):
   case(Points):
   case(ContactElementsZero):
   case(ContactElementsOne):
   case(ReducedElements):
   case(ReducedFaceElements):
   case(ReducedContactElementsZero):
   case(ReducedContactElementsOne):
   return true;
   break;
   default:
      stringstream temp;
      temp << "Error - Cell: Finley does not know anything about function space type " << functionSpaceCode;
      throw FinleyAdapterException(temp.str());
      break;
   }
   checkFinleyError();
   return false;
}

bool
MeshAdapter::commonFunctionSpace(const vector<int>& fs, int& resultcode) const
{
   /* The idea is to use equivalence classes. [Types which can be interpolated back and forth]
	class 1: DOF <-> Nodes
	class 2: ReducedDOF <-> ReducedNodes
	class 3: Points
	class 4: Elements
	class 5: ReducedElements
	class 6: FaceElements
	class 7: ReducedFaceElements
	class 8: ContactElementZero <-> ContactElementOne
	class 9: ReducedContactElementZero <-> ReducedContactElementOne

   There is also a set of lines. Interpolation is possible down a line but not between lines.
   class 1 and 2 belong to all lines so aren't considered.
	line 0: class 3
	line 1: class 4,5
	line 2: class 6,7
	line 3: class 8,9

   For classes with multiple members (eg class 2) we have vars to record if there is at least one instance.
   eg hasnodes is true if we have at least one instance of Nodes.
   */
    if (fs.empty())
    {
        return false;
    }
    vector<int> hasclass(10);
    vector<int> hasline(4);	
    bool hasnodes=false;
    bool hasrednodes=false;
    bool hascez=false;
    bool hasrcez=false;
    for (int i=0;i<fs.size();++i)
    {
	switch(fs[i])
	{
	case(Nodes):   hasnodes=true;	// no break is deliberate
	case(DegreesOfFreedom):
		hasclass[1]=1;
		break;
	case(ReducedNodes):    hasrednodes=true;	// no break is deliberate
	case(ReducedDegreesOfFreedom):
		hasclass[2]=1;
		break;
	case(Points):
		hasline[0]=1;
		hasclass[3]=1;
		break;
	case(Elements):
		hasclass[4]=1;
		hasline[1]=1;
		break;
	case(ReducedElements):
		hasclass[5]=1;
		hasline[1]=1;
		break;
	case(FaceElements):
		hasclass[6]=1;
		hasline[2]=1;
		break;
	case(ReducedFaceElements):
		hasclass[7]=1;
		hasline[2]=1;
		break;
	case(ContactElementsZero):  hascez=true;	// no break is deliberate
	case(ContactElementsOne):
		hasclass[8]=1;
		hasline[3]=1;
		break;
	case(ReducedContactElementsZero):   hasrcez=true;   // no break is deliberate
	case(ReducedContactElementsOne):
		hasclass[9]=1;
		hasline[3]=1;
		break;
	default:
		return false;
	}
    }
    int totlines=hasline[0]+hasline[1]+hasline[2]+hasline[3];
    // fail if we have more than one leaf group

    if (totlines>1)
    {
	return false;	// there are at least two branches we can't interpolate between
    }
    else if (totlines==1)
    {
	if (hasline[0]==1)		// we have points
	{
	    resultcode=Points;
	}
	else if (hasline[1]==1)
	{
	    if (hasclass[5]==1)
	    {
		resultcode=ReducedElements;
	    }
	    else
	    {
		resultcode=Elements;
	    }
	}
	else if (hasline[2]==1)
	{
	    if (hasclass[7]==1)
	    {
		resultcode=ReducedFaceElements;
	    }
	    else
	    {
		resultcode=FaceElements;
	    }
	}
	else	// so we must be in line3
	{
	    if (hasclass[9]==1)
	    {
		// need something from class 9
		resultcode=(hasrcez?ReducedContactElementsZero:ReducedContactElementsOne);
	    }
	    else
	    {
		// something from class 8
		resultcode=(hascez?ContactElementsZero:ContactElementsOne);
	    }
	}
    }
    else	// totlines==0
    {
	if (hasclass[2]==1)
	{
		// something from class 2
		resultcode=(hasrednodes?ReducedNodes:ReducedDegreesOfFreedom);
	}
	else
	{	// something from class 1
		resultcode=(hasnodes?Nodes:DegreesOfFreedom);
	}
    }
    return true;
}

bool MeshAdapter::probeInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const
{
   switch(functionSpaceType_source) {
   case(Nodes):
   	switch(functionSpaceType_target) {
	case(Nodes):
	case(ReducedNodes):
	case(ReducedDegreesOfFreedom):
	case(DegreesOfFreedom):
	case(Elements):
	case(ReducedElements):
	case(FaceElements):
	case(ReducedFaceElements):
	case(Points):
	case(ContactElementsZero):
	case(ReducedContactElementsZero):
	case(ContactElementsOne):
	case(ReducedContactElementsOne):
	return true;
	default:
	      stringstream temp;
	      temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << functionSpaceType_target;
	      throw FinleyAdapterException(temp.str());
   }
   break;
   case(ReducedNodes):
	switch(functionSpaceType_target) {
	case(ReducedNodes):
	case(ReducedDegreesOfFreedom):
	case(Elements):
	case(ReducedElements):
	case(FaceElements):
	case(ReducedFaceElements):
	case(Points):
	case(ContactElementsZero):
	case(ReducedContactElementsZero):
	case(ContactElementsOne):
	case(ReducedContactElementsOne):
	return true;
	case(Nodes):
	case(DegreesOfFreedom):
	return false;
	default:
		stringstream temp;
		temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << functionSpaceType_target;
		throw FinleyAdapterException(temp.str());
   }
   break;
   case(Elements):
	if (functionSpaceType_target==Elements) {
	  return true;
	} else if (functionSpaceType_target==ReducedElements) {
	  return true;
        } else {
          return false;
        }
   case(ReducedElements):
	if (functionSpaceType_target==ReducedElements) {
	  return true;
	} else {
          return false;
	}
   case(FaceElements):
	if (functionSpaceType_target==FaceElements) {
      		return true;
	} else if (functionSpaceType_target==ReducedFaceElements) {
      		return true;
	} else {
      		return false;
	}
   case(ReducedFaceElements):
	if (functionSpaceType_target==ReducedFaceElements) {
      		return true;
	} else {
		return false;
	}
   case(Points):
	if (functionSpaceType_target==Points) {
      		return true;
	} else {
      		return false;
	}
   case(ContactElementsZero):
   case(ContactElementsOne):
	if (functionSpaceType_target==ContactElementsZero || functionSpaceType_target==ContactElementsOne) {
      		return true;
	} else if (functionSpaceType_target==ReducedContactElementsZero || functionSpaceType_target==ReducedContactElementsOne) {
      		return true;
	} else {
      		return false;
	}
   case(ReducedContactElementsZero):
   case(ReducedContactElementsOne):
	if (functionSpaceType_target==ReducedContactElementsZero || functionSpaceType_target==ReducedContactElementsOne) {
      		return true;
	} else {
      		return false;
	}
   case(DegreesOfFreedom):
	switch(functionSpaceType_target) {
	case(ReducedDegreesOfFreedom):
	case(DegreesOfFreedom):
	case(Nodes):
	case(ReducedNodes):
	case(Elements):
	case(ReducedElements):
	case(Points):
	case(FaceElements):
	case(ReducedFaceElements):
	case(ContactElementsZero):
	case(ReducedContactElementsZero):
	case(ContactElementsOne):
	case(ReducedContactElementsOne):
	return true;
	default:
		stringstream temp;
		temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << functionSpaceType_target;
		throw FinleyAdapterException(temp.str());
	}
   	break;
   case(ReducedDegreesOfFreedom):
   switch(functionSpaceType_target) {
	case(ReducedDegreesOfFreedom):
	case(ReducedNodes):
	case(Elements):
	case(ReducedElements):
	case(FaceElements):
	case(ReducedFaceElements):
	case(Points):
	case(ContactElementsZero):
	case(ReducedContactElementsZero):
	case(ContactElementsOne):
	case(ReducedContactElementsOne):
	return true;
	case(Nodes):
	case(DegreesOfFreedom):
	return false;
	default:
		stringstream temp;
		temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << functionSpaceType_target;
		throw FinleyAdapterException(temp.str());
	}
	break;
   default:
      stringstream temp;
      temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << functionSpaceType_source;
      throw FinleyAdapterException(temp.str());
      break;
   }
   checkFinleyError();
   return false;
}

signed char MeshAdapter::preferredInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const
{
    if (probeInterpolationOnDomain(functionSpaceType_source, functionSpaceType_target))
    {
        return 1;
    }
    else if (probeInterpolationOnDomain(functionSpaceType_target, functionSpaceType_source))
    {
        return -1;
    }
    return 0;
}  
  
bool MeshAdapter::probeInterpolationACross(int functionSpaceType_source,const AbstractDomain& targetDomain, int functionSpaceType_target) const
{
   return false;
}

bool MeshAdapter::operator==(const AbstractDomain& other) const
{
   const MeshAdapter* temp=dynamic_cast<const MeshAdapter*>(&other);
   if (temp!=0) {
      return (m_finleyMesh==temp->m_finleyMesh);
   } else {
      return false;
   }
}

bool MeshAdapter::operator!=(const AbstractDomain& other) const
{
   return !(operator==(other));
}

int MeshAdapter::getSystemMatrixTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const
{
   Finley_Mesh* mesh=m_finleyMesh.get();
   return SystemMatrixAdapter::getSystemMatrixTypeId(solver, preconditioner,
           package, symmetry, mesh->MPIInfo);
}

int MeshAdapter::getTransportTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const
{
   Finley_Mesh* mesh=m_finleyMesh.get();
   return TransportProblemAdapter::getTransportTypeId(solver, preconditioner,
           package, symmetry, mesh->MPIInfo);
}

escript::Data MeshAdapter::getX() const
{
   return continuousFunction(*this).getX();
}

escript::Data MeshAdapter::getNormal() const
{
   return functionOnBoundary(*this).getNormal();
}

escript::Data MeshAdapter::getSize() const
{
   return escript::function(*this).getSize();
}

const int* MeshAdapter::borrowSampleReferenceIDs(int functionSpaceType) const
{
   int *out = NULL;
   Finley_Mesh* mesh=m_finleyMesh.get();
   switch (functionSpaceType) {
   case(Nodes):
   out=mesh->Nodes->Id;
   break;
   case(ReducedNodes):
   out=mesh->Nodes->reducedNodesId;
   break;
   case(Elements):
   out=mesh->Elements->Id;
   break;
   case(ReducedElements):
   out=mesh->Elements->Id;
   break;
   case(FaceElements):
   out=mesh->FaceElements->Id;
   break;
   case(ReducedFaceElements):
   out=mesh->FaceElements->Id;
   break;
   case(Points):
   out=mesh->Points->Id;
   break;
   case(ContactElementsZero):
   out=mesh->ContactElements->Id;
   break;
   case(ReducedContactElementsZero):
   out=mesh->ContactElements->Id;
   break;
   case(ContactElementsOne):
   out=mesh->ContactElements->Id;
   break;
   case(ReducedContactElementsOne):
   out=mesh->ContactElements->Id;
   break;
   case(DegreesOfFreedom):
   out=mesh->Nodes->degreesOfFreedomId;
   break;
   case(ReducedDegreesOfFreedom):
   out=mesh->Nodes->reducedDegreesOfFreedomId;
   break;
   default:
      stringstream temp;
      temp << "Error - Invalid function space type: " << functionSpaceType << " for domain: " << getDescription();
      throw FinleyAdapterException(temp.str());
      break;
   }
   return out;
}
int MeshAdapter::getTagFromSampleNo(int functionSpaceType, int sampleNo) const
{
   int out=0;
   Finley_Mesh* mesh=m_finleyMesh.get();
   switch (functionSpaceType) {
   case(Nodes):
   out=mesh->Nodes->Tag[sampleNo];
   break;
   case(ReducedNodes):
   throw FinleyAdapterException(" Error - ReducedNodes does not support tags.");
   break;
   case(Elements):
   out=mesh->Elements->Tag[sampleNo];
   break;
   case(ReducedElements):
   out=mesh->Elements->Tag[sampleNo];
   break;
   case(FaceElements):
   out=mesh->FaceElements->Tag[sampleNo];
   break;
   case(ReducedFaceElements):
   out=mesh->FaceElements->Tag[sampleNo];
   break;
   case(Points):
   out=mesh->Points->Tag[sampleNo];
   break;
   case(ContactElementsZero):
   out=mesh->ContactElements->Tag[sampleNo];
   break;
   case(ReducedContactElementsZero):
   out=mesh->ContactElements->Tag[sampleNo];
   break;
   case(ContactElementsOne):
   out=mesh->ContactElements->Tag[sampleNo];
   break;
   case(ReducedContactElementsOne):
   out=mesh->ContactElements->Tag[sampleNo];
   break;
   case(DegreesOfFreedom):
   throw FinleyAdapterException(" Error - DegreesOfFreedom does not support tags.");
   break;
   case(ReducedDegreesOfFreedom):
   throw FinleyAdapterException(" Error - ReducedDegreesOfFreedom does not support tags.");
   break;
   default:
      stringstream temp;
      temp << "Error - Invalid function space type: " << functionSpaceType << " for domain: " << getDescription();
      throw FinleyAdapterException(temp.str());
      break;
   }
   return out;
}


void MeshAdapter::setTags(const int functionSpaceType, const int newTag, const escript::Data& mask) const
{
   Finley_Mesh* mesh=m_finleyMesh.get();
   escriptDataC tmp=mask.getDataC();
   switch(functionSpaceType) {
   case(Nodes):
   Finley_NodeFile_setTags(mesh->Nodes,newTag,&tmp);
   break;
   case(ReducedNodes):
   throw FinleyAdapterException("Error - ReducedNodes does not support tags");
   break;
   case(DegreesOfFreedom):
   throw FinleyAdapterException("Error - DegreesOfFreedom does not support tags");
   break;
   case(ReducedDegreesOfFreedom):
   throw FinleyAdapterException("Error - ReducedDegreesOfFreedom does not support tags");
   break;
   case(Elements):
   Finley_ElementFile_setTags(mesh->Elements,newTag,&tmp);
   break;
   case(ReducedElements):
   Finley_ElementFile_setTags(mesh->Elements,newTag,&tmp);
   break;
   case(FaceElements):
   Finley_ElementFile_setTags(mesh->FaceElements,newTag,&tmp);
   break;
   case(ReducedFaceElements):
   Finley_ElementFile_setTags(mesh->FaceElements,newTag,&tmp);
   break;
   case(Points):
   Finley_ElementFile_setTags(mesh->Points,newTag,&tmp);
   break;
   case(ContactElementsZero):
   Finley_ElementFile_setTags(mesh->ContactElements,newTag,&tmp);
   break;
   case(ReducedContactElementsZero):
   Finley_ElementFile_setTags(mesh->ContactElements,newTag,&tmp);
   break;
   case(ContactElementsOne):
   Finley_ElementFile_setTags(mesh->ContactElements,newTag,&tmp);
   break;
   case(ReducedContactElementsOne):
   Finley_ElementFile_setTags(mesh->ContactElements,newTag,&tmp);
   break;
   default:
      stringstream temp;
      temp << "Error - Finley does not know anything about function space type " << functionSpaceType;
      throw FinleyAdapterException(temp.str());
   }
   checkFinleyError();
   return;
}

void MeshAdapter::setTagMap(const string& name,  int tag)
{
   Finley_Mesh* mesh=m_finleyMesh.get();
   Finley_Mesh_addTagMap(mesh, name.c_str(),tag);
   checkFinleyError();
   // throwStandardException("MeshAdapter::set TagMap is not implemented.");
}

int MeshAdapter::getTag(const string& name) const
{
   Finley_Mesh* mesh=m_finleyMesh.get();
   int tag=0;
   tag=Finley_Mesh_getTag(mesh, name.c_str());
   checkFinleyError();
   // throwStandardException("MeshAdapter::getTag is not implemented.");
   return tag;
}

bool MeshAdapter::isValidTagName(const string& name) const
{
   Finley_Mesh* mesh=m_finleyMesh.get();
   return Finley_Mesh_isValidTagName(mesh,name.c_str());
}

string MeshAdapter::showTagNames() const
{
   stringstream temp;
   Finley_Mesh* mesh=m_finleyMesh.get();
   Finley_TagMap* tag_map=mesh->TagMap;
   while (tag_map) {
      temp << tag_map->name;
      tag_map=tag_map->next;
      if (tag_map) temp << ", ";
   }
   return temp.str();
}

int MeshAdapter::getNumberOfTagsInUse(int functionSpaceCode) const
{
  Finley_Mesh* mesh=m_finleyMesh.get();
  dim_t numTags=0;
  switch(functionSpaceCode) {
   case(Nodes):
          numTags=mesh->Nodes->numTagsInUse;
          break;
   case(ReducedNodes):
          throw FinleyAdapterException("Error - ReducedNodes does not support tags");
          break;
   case(DegreesOfFreedom):
          throw FinleyAdapterException("Error - DegreesOfFreedom does not support tags");
          break;
   case(ReducedDegreesOfFreedom):
          throw FinleyAdapterException("Error - ReducedDegreesOfFreedom does not support tags");
          break;
   case(Elements):
   case(ReducedElements):
          numTags=mesh->Elements->numTagsInUse;
          break;
   case(FaceElements):
   case(ReducedFaceElements):
          numTags=mesh->FaceElements->numTagsInUse;
          break;
   case(Points):
          numTags=mesh->Points->numTagsInUse;
          break;
   case(ContactElementsZero):
   case(ReducedContactElementsZero):
   case(ContactElementsOne):
   case(ReducedContactElementsOne):
          numTags=mesh->ContactElements->numTagsInUse;
          break;
   default:
      stringstream temp;
      temp << "Error - Finley does not know anything about function space type " << functionSpaceCode;
      throw FinleyAdapterException(temp.str());
  }
  return numTags;
}

const int* MeshAdapter::borrowListOfTagsInUse(int functionSpaceCode) const
{
  Finley_Mesh* mesh=m_finleyMesh.get();
  index_t* tags=NULL;
  switch(functionSpaceCode) {
   case(Nodes):
          tags=mesh->Nodes->tagsInUse;
          break;
   case(ReducedNodes):
          throw FinleyAdapterException("Error - ReducedNodes does not support tags");
          break;
   case(DegreesOfFreedom):
          throw FinleyAdapterException("Error - DegreesOfFreedom does not support tags");
          break;
   case(ReducedDegreesOfFreedom):
          throw FinleyAdapterException("Error - ReducedDegreesOfFreedom does not support tags");
          break;
   case(Elements):
   case(ReducedElements):
          tags=mesh->Elements->tagsInUse;
          break;
   case(FaceElements):
   case(ReducedFaceElements):
          tags=mesh->FaceElements->tagsInUse;
          break;
   case(Points):
          tags=mesh->Points->tagsInUse;
          break;
   case(ContactElementsZero):
   case(ReducedContactElementsZero):
   case(ContactElementsOne):
   case(ReducedContactElementsOne):
          tags=mesh->ContactElements->tagsInUse;
          break;
   default:
      stringstream temp;
      temp << "Error - Finley does not know anything about function space type " << functionSpaceCode;
      throw FinleyAdapterException(temp.str());
  }
  return tags;
}


bool MeshAdapter::canTag(int functionSpaceCode) const
{
  switch(functionSpaceCode) {
   case(Nodes):
   case(Elements):
   case(ReducedElements):
   case(FaceElements):
   case(ReducedFaceElements):
   case(Points):
   case(ContactElementsZero):
   case(ReducedContactElementsZero):
   case(ContactElementsOne):
   case(ReducedContactElementsOne):
          return true;
   case(ReducedNodes):
   case(DegreesOfFreedom):
   case(ReducedDegreesOfFreedom):
	  return false;
   default:
	return false;
  }
}

AbstractDomain::StatusType MeshAdapter::getStatus() const
{
  Finley_Mesh* mesh=m_finleyMesh.get();
  return Finley_Mesh_getStatus(mesh);
}

int MeshAdapter::getApproximationOrder(const int functionSpaceCode) const
{
   
  Finley_Mesh* mesh=m_finleyMesh.get();
  int order =-1;
  switch(functionSpaceCode) {
   case(Nodes):
   case(DegreesOfFreedom):
          order=mesh->approximationOrder;
          break;
   case(ReducedNodes):
   case(ReducedDegreesOfFreedom):
          order=mesh->reducedApproximationOrder;
          break;
   case(Elements):
   case(FaceElements):
   case(Points):
   case(ContactElementsZero):
   case(ContactElementsOne):
          order=mesh->integrationOrder;
          break;
   case(ReducedElements):
   case(ReducedFaceElements):
   case(ReducedContactElementsZero):
   case(ReducedContactElementsOne):
          order=mesh->reducedIntegrationOrder;
          break;
   default:
      stringstream temp;
      temp << "Error - Finley does not know anything about function space type " << functionSpaceCode;
      throw FinleyAdapterException(temp.str());
  }
  return order;
}

bool MeshAdapter::supportsContactElements() const
{
  return true;
}

ReferenceElementSetWrapper::ReferenceElementSetWrapper(Finley_ElementTypeId id, index_t order, index_t reducedOrder)
{
  m_refSet = Finley_ReferenceElementSet_alloc(id, order, reducedOrder);
}

ReferenceElementSetWrapper::~ReferenceElementSetWrapper()
{
  Finley_ReferenceElementSet_dealloc(m_refSet);
}

// points will be flattened
void MeshAdapter:: addDiracPoints(const std::vector<double>& points, const std::vector<int>& tags) const
{
      const int dim = getDim();
      int numPoints=points.size()/dim;
      int numTags=tags.size();
      Finley_Mesh* mesh=m_finleyMesh.get();
      
      if ( points.size() % dim != 0 )
      {
	throw FinleyAdapterException("Error - number of coords does not appear to be a multiple of dimension.");
      }
      
      if  ( (numTags > 0) && ( numPoints !=  numTags ) )
	 throw FinleyAdapterException("Error - if tags are given number of tags and points must match.");
      
      double* points_ptr=new  double[numPoints * dim];
      int*    tags_ptr= new  int[numPoints];
      
      for (int i=0;i<numPoints;++i) {
	   points_ptr[ i * dim     ] = points[i * dim ];
	   if ( dim > 1 ) points_ptr[ i * dim + 1 ] = points[i * dim + 1];
	   if ( dim > 2 ) points_ptr[ i * dim + 2 ] = points[i * dim + 2];	  
           tags_ptr[i]=tags[i];
      }
      
      Finley_Mesh_addPoints(mesh, numPoints, points_ptr, tags_ptr);
      checkFinleyError();
      
      delete[] points_ptr;
      delete[] tags_ptr;
}


// void MeshAdapter:: addDiracPoints(const boost::python::list& points, const boost::python::list& tags) const
// {
//       const int dim = getDim();
//       int numPoints=boost::python::extract<int>(points.attr("__len__")());
//       int numTags=boost::python::extract<int>(tags.attr("__len__")());
//       Finley_Mesh* mesh=m_finleyMesh.get();
//       
//       if  ( (numTags > 0) && ( numPoints !=  numTags ) )
// 	 throw FinleyAdapterException("Error - if tags are given number of tags and points must match.");
//       
//       double* points_ptr=TMPMEMALLOC(numPoints * dim, double);
//       int*    tags_ptr= TMPMEMALLOC(numPoints, int);
//       
//       for (int i=0;i<numPoints;++i) {
// 	   int tag_id=-1;
// 	   int numComps=boost::python::extract<int>(points[i].attr("__len__")());
// 	   if  ( numComps !=   dim ) {
//                stringstream temp;	        
//                temp << "Error - illegal number of components " << numComps << " for point " << i;               
//                throw FinleyAdapterException(temp.str());
// 	   }
// 	   points_ptr[ i * dim     ] = boost::python::extract<double>(points[i][0]);
// 	   if ( dim > 1 ) points_ptr[ i * dim + 1 ] = boost::python::extract<double>(points[i][1]);
// 	   if ( dim > 2 ) points_ptr[ i * dim + 2 ] = boost::python::extract<double>(points[i][2]);
// 	   
// 	   if ( numTags > 0) {
// 	          boost::python::extract<string> ex_str(tags[i]);
// 		  if  ( ex_str.check() ) {
// 		      tag_id=getTag( ex_str());
// 		  } else {
// 		       boost::python::extract<int> ex_int(tags[i]);
// 		       if ( ex_int.check() ) {
// 		           tag_id=ex_int();
// 		       } else { 
// 			    stringstream temp;	         
// 		            temp << "Error - unable to extract tag for point " << i; 
// 			    throw FinleyAdapterException(temp.str());
// 		      }
// 		  }
// 	   }	  
//            tags_ptr[i]=tag_id;
//       }
//       
//       Finley_Mesh_addPoints(mesh, numPoints, points_ptr, tags_ptr);
//       checkPasoError();
//       
//       TMPMEMFREE(points_ptr);
//       TMPMEMFREE(tags_ptr);
// }

/*
void MeshAdapter:: addDiracPoint( const boost::python::list& point, const int tag) const
{   
    boost::python::list points =  boost::python::list();
    boost::python::list tags =  boost::python::list();
    points.append(point);
    tags.append(tag);
    addDiracPoints(points, tags);
}
*/

/*
void MeshAdapter:: addDiracPointWithTagName( const boost::python::list& point, const std::string& tag) const
{
        boost::python::list points =   boost::python::list();
        boost::python::list tags =   boost::python::list();
        points.append(point);
        tags.append(tag);
        addDiracPoints(points, tags);
} 
*/
}  // end of namespace
