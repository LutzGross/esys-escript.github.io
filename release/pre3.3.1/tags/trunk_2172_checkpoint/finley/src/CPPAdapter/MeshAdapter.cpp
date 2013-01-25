
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "MeshAdapter.h"
#include "escript/Data.h"
#include "escript/DataFactory.h"
#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif
#ifdef PASO_MPI
#include <mpi.h>
#include "paso/Paso_MPI.h"
#endif
extern "C" {
#include "esysUtils/blocktimer.h"
}

using namespace std;
using namespace escript;

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
#ifdef PASO_MPI
   MPI_Barrier(m_finleyMesh.get()->MPIInfo->comm);
#endif
   return;
}
bool MeshAdapter::onMasterProcessor() const
{
   return m_finleyMesh.get()->MPIInfo->rank == 0;
}


Finley_Mesh* MeshAdapter::getFinley_Mesh() const {
   return m_finleyMesh.get();
}

void MeshAdapter::write(const std::string& fileName) const
{
   char *fName = (fileName.size()+1>0) ? TMPMEMALLOC(fileName.size()+1,char) : (char*)NULL;
   strcpy(fName,fileName.c_str());
   Finley_Mesh_write(m_finleyMesh.get(),fName);
   checkFinleyError();
   TMPMEMFREE(fName);
}

void MeshAdapter::Print_Mesh_Info(const bool full=false) const
{
   Finley_PrintMesh_Info(m_finleyMesh.get(), full);
}

void MeshAdapter::dump(const std::string& fileName) const
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
#ifdef PASO_MPI
   MPI_Status status;
#endif

/* Incoming token indicates it's my turn to write */
#ifdef PASO_MPI
   if (mpi_rank>0) MPI_Recv(&num_Tags, 0, MPI_INT, mpi_rank-1, 81800, mesh->MPIInfo->comm, &status);
#endif

   char *newFileName = Paso_MPI_appendRankToFileName(fileName.c_str(),
                                                     mpi_size, mpi_rank);

   /* Figure out how much storage is required for tags */
   tag_map = mesh->TagMap;
   num_Tags = 0;
   if (tag_map) {
      while (tag_map) {
         num_Tags++;
         tag_map=tag_map->next;
      }
   }

   // NetCDF error handler
   NcError err(NcError::verbose_nonfatal);
   // Create the file.
   NcFile dataFile(newFileName, NcFile::Replace);
   // check if writing was successful
   if (!dataFile.is_valid())
      throw DataException("Error - MeshAdapter::dump: opening of NetCDF file for output failed: " + *newFileName);

   // Define dimensions (num_Elements and dim_Elements are identical, dim_Elements only appears if > 0)
   if (! (ncdims[0] = dataFile.add_dim("numNodes", numNodes)) )
      throw DataException("Error - MeshAdapter::dump: appending dimension numNodes to netCDF file failed: " + *newFileName);
   if (! (ncdims[1] = dataFile.add_dim("numDim", numDim)) )
      throw DataException("Error - MeshAdapter::dump: appending dimension numDim to netCDF file failed: " + *newFileName);
   if (! (ncdims[2] = dataFile.add_dim("mpi_size_plus_1", mpi_size+1)) )
      throw DataException("Error - MeshAdapter::dump: appending dimension mpi_size to netCDF file failed: " + *newFileName);
   if (num_Elements>0)
      if (! (ncdims[3] = dataFile.add_dim("dim_Elements", num_Elements)) )
         throw DataException("Error - MeshAdapter::dump: appending dimension dim_Elements to netCDF file failed: " + *newFileName);
   if (num_FaceElements>0)
      if (! (ncdims[4] = dataFile.add_dim("dim_FaceElements", num_FaceElements)) )
         throw DataException("Error - MeshAdapter::dump: appending dimension dim_FaceElements to netCDF file failed: " + *newFileName);
   if (num_ContactElements>0)
      if (! (ncdims[5] = dataFile.add_dim("dim_ContactElements", num_ContactElements)) )
         throw DataException("Error - MeshAdapter::dump: appending dimension dim_ContactElements to netCDF file failed: " + *newFileName);
   if (num_Points>0)
      if (! (ncdims[6] = dataFile.add_dim("dim_Points", num_Points)) )
         throw DataException("Error - MeshAdapter::dump: appending dimension dim_Points to netCDF file failed: " + *newFileName);
   if (num_Elements>0)
      if (! (ncdims[7] = dataFile.add_dim("dim_Elements_Nodes", num_Elements_numNodes)) )
         throw DataException("Error - MeshAdapter::dump: appending dimension dim_Elements_Nodes to netCDF file failed: " + *newFileName);
   if (num_FaceElements>0)
      if (! (ncdims[8] = dataFile.add_dim("dim_FaceElements_numNodes", num_FaceElements_numNodes)) )
         throw DataException("Error - MeshAdapter::dump: appending dimension dim_FaceElements_numNodes to netCDF file failed: " + *newFileName);
   if (num_ContactElements>0)
      if (! (ncdims[9] = dataFile.add_dim("dim_ContactElements_numNodes", num_ContactElements_numNodes)) )
         throw DataException("Error - MeshAdapter::dump: appending dimension dim_ContactElements_numNodes to netCDF file failed: " + *newFileName);
   if (num_Tags>0)
      if (! (ncdims[10] = dataFile.add_dim("dim_Tags", num_Tags)) )
         throw DataException("Error - MeshAdapter::dump: appending dimension dim_Tags to netCDF file failed: " + *newFileName);

   // Attributes: MPI size, MPI rank, Name, order, reduced_order
   if (!dataFile.add_att("mpi_size", mpi_size) )
      throw DataException("Error - MeshAdapter::dump: appending mpi_size to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("mpi_rank", mpi_rank) )
      throw DataException("Error - MeshAdapter::dump: appending mpi_rank to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("Name",mesh->Name) )
      throw DataException("Error - MeshAdapter::dump: appending Name to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("numDim",numDim) )
      throw DataException("Error - MeshAdapter::dump: appending order to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("order",mesh->order) )
      throw DataException("Error - MeshAdapter::dump: appending order to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("reduced_order",mesh->reduced_order) )
      throw DataException("Error - MeshAdapter::dump: appending reduced_order to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("numNodes",numNodes) )
      throw DataException("Error - MeshAdapter::dump: appending numNodes to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("num_Elements",num_Elements) )
      throw DataException("Error - MeshAdapter::dump: appending num_Elements to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("num_FaceElements",num_FaceElements) )
      throw DataException("Error - MeshAdapter::dump: appending num_FaceElements to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("num_ContactElements",num_ContactElements) )
      throw DataException("Error - MeshAdapter::dump: appending num_ContactElements to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("num_Points",num_Points) )
      throw DataException("Error - MeshAdapter::dump: appending num_Points to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("num_Elements_numNodes",num_Elements_numNodes) )
      throw DataException("Error - MeshAdapter::dump: appending num_Elements_numNodes to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("num_FaceElements_numNodes",num_FaceElements_numNodes) )
      throw DataException("Error - MeshAdapter::dump: appending num_FaceElements_numNodes to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("num_ContactElements_numNodes",num_ContactElements_numNodes) )
      throw DataException("Error - MeshAdapter::dump: appending num_ContactElements_numNodes to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("Elements_TypeId", mesh->Elements->ReferenceElement->Type->TypeId) )
      throw DataException("Error - MeshAdapter::dump: appending Elements_TypeId to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("FaceElements_TypeId", mesh->FaceElements->ReferenceElement->Type->TypeId) )
      throw DataException("Error - MeshAdapter::dump: appending FaceElements_TypeId to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("ContactElements_TypeId", mesh->ContactElements->ReferenceElement->Type->TypeId) )
      throw DataException("Error - MeshAdapter::dump: appending ContactElements_TypeId to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("Points_TypeId", mesh->Points->ReferenceElement->Type->TypeId) )
      throw DataException("Error - MeshAdapter::dump: appending Points_TypeId to NetCDF file failed: " + *newFileName);
   if (!dataFile.add_att("num_Tags", num_Tags) )
      throw DataException("Error - MeshAdapter::dump: appending num_Tags to NetCDF file failed: " + *newFileName);

   // // // // // Nodes // // // // //

   // Only write nodes if non-empty because NetCDF doesn't like empty arrays (it treats them as NC_UNLIMITED)
   if (numNodes>0) {

      // Nodes Id
      if (! ( ids = dataFile.add_var("Nodes_Id", ncInt, ncdims[0])) )
         throw DataException("Error - MeshAdapter::dump: appending Nodes_Id to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->Nodes->Id[0];
      if (! (ids->put(int_ptr, numNodes)) )
         throw DataException("Error - MeshAdapter::dump: copy Nodes_Id to netCDF buffer failed: " + *newFileName);

      // Nodes Tag
      if (! ( ids = dataFile.add_var("Nodes_Tag", ncInt, ncdims[0])) )
         throw DataException("Error - MeshAdapter::dump: appending Nodes_Tag to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->Nodes->Tag[0];
      if (! (ids->put(int_ptr, numNodes)) )
         throw DataException("Error - MeshAdapter::dump: copy Nodes_Tag to netCDF buffer failed: " + *newFileName);

      // Nodes gDOF
      if (! ( ids = dataFile.add_var("Nodes_gDOF", ncInt, ncdims[0])) )
         throw DataException("Error - MeshAdapter::dump: appending Nodes_gDOF to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->Nodes->globalDegreesOfFreedom[0];
      if (! (ids->put(int_ptr, numNodes)) )
         throw DataException("Error - MeshAdapter::dump: copy Nodes_gDOF to netCDF buffer failed: " + *newFileName);

      // Nodes global node index
      if (! ( ids = dataFile.add_var("Nodes_gNI", ncInt, ncdims[0])) )
         throw DataException("Error - MeshAdapter::dump: appending Nodes_gNI to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->Nodes->globalNodesIndex[0];
      if (! (ids->put(int_ptr, numNodes)) )
         throw DataException("Error - MeshAdapter::dump: copy Nodes_gNI to netCDF buffer failed: " + *newFileName);

      // Nodes grDof
      if (! ( ids = dataFile.add_var("Nodes_grDfI", ncInt, ncdims[0])) )
         throw DataException("Error - MeshAdapter::dump: appending Nodes_grDfI to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->Nodes->globalReducedDOFIndex[0];
      if (! (ids->put(int_ptr, numNodes)) )
         throw DataException("Error - MeshAdapter::dump: copy Nodes_grDfI to netCDF buffer failed: " + *newFileName);

      // Nodes grNI
      if (! ( ids = dataFile.add_var("Nodes_grNI", ncInt, ncdims[0])) )
         throw DataException("Error - MeshAdapter::dump: appending Nodes_grNI to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->Nodes->globalReducedNodesIndex[0];
      if (! (ids->put(int_ptr, numNodes)) )
         throw DataException("Error - MeshAdapter::dump: copy Nodes_grNI to netCDF buffer failed: " + *newFileName);

      // Nodes Coordinates
      if (! ( ids = dataFile.add_var("Nodes_Coordinates", ncDouble, ncdims[0], ncdims[1]) ) )
         throw DataException("Error - MeshAdapter::dump: appending Nodes_Coordinates to netCDF file failed: " + *newFileName);
      if (! (ids->put(&(mesh->Nodes->Coordinates[INDEX2(0,0,numDim)]), numNodes, numDim)) )
         throw DataException("Error - MeshAdapter::dump: copy Nodes_Coordinates to netCDF buffer failed: " + *newFileName);

      // Nodes degreesOfFreedomDistribution
      if (! ( ids = dataFile.add_var("Nodes_DofDistribution", ncInt, ncdims[2])) )
         throw DataException("Error - MeshAdapter::dump: appending Nodes_DofDistribution to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->Nodes->degreesOfFreedomDistribution->first_component[0];
      if (! (ids->put(int_ptr, mpi_size+1)) )
         throw DataException("Error - MeshAdapter::dump: copy Nodes_DofDistribution to netCDF buffer failed: " + *newFileName);

      // Nodes nodeDistribution
      if (! ( ids = dataFile.add_var("Nodes_NodeDistribution", ncInt, ncdims[2])) )
         throw DataException("Error - MeshAdapter::dump: appending Nodes_NodeDistribution to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->Nodes->nodesDistribution->first_component[0];
      if (! (ids->put(int_ptr, mpi_size+1)) )
         throw DataException("Error - MeshAdapter::dump: copy Nodes_NodeDistribution to netCDF buffer failed: " + *newFileName);

   }

   // // // // // Elements // // // // //

   if (num_Elements>0) {

      // Elements_Id
      if (! ( ids = dataFile.add_var("Elements_Id", ncInt, ncdims[3])) )
         throw DataException("Error - MeshAdapter::dump: appending Elements_Id to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->Elements->Id[0];
      if (! (ids->put(int_ptr, num_Elements)) )
         throw DataException("Error - MeshAdapter::dump: copy Elements_Id to netCDF buffer failed: " + *newFileName);

      // Elements_Tag
      if (! ( ids = dataFile.add_var("Elements_Tag", ncInt, ncdims[3])) )
         throw DataException("Error - MeshAdapter::dump: appending Elements_Tag to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->Elements->Tag[0];
      if (! (ids->put(int_ptr, num_Elements)) )
         throw DataException("Error - MeshAdapter::dump: copy Elements_Tag to netCDF buffer failed: " + *newFileName);

      // Elements_Owner
      if (! ( ids = dataFile.add_var("Elements_Owner", ncInt, ncdims[3])) )
         throw DataException("Error - MeshAdapter::dump: appending Elements_Owner to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->Elements->Owner[0];
      if (! (ids->put(int_ptr, num_Elements)) )
         throw DataException("Error - MeshAdapter::dump: copy Elements_Owner to netCDF buffer failed: " + *newFileName);

      // Elements_Color
      if (! ( ids = dataFile.add_var("Elements_Color", ncInt, ncdims[3])) )
         throw DataException("Error - MeshAdapter::dump: appending Elements_Color to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->Elements->Color[0];
      if (! (ids->put(int_ptr, num_Elements)) )
         throw DataException("Error - MeshAdapter::dump: copy Elements_Color to netCDF buffer failed: " + *newFileName);

      // Elements_Nodes
      if (! ( ids = dataFile.add_var("Elements_Nodes", ncInt, ncdims[3], ncdims[7]) ) )
         throw DataException("Error - MeshAdapter::dump: appending Elements_Nodes to netCDF file failed: " + *newFileName);
      if (! (ids->put(&(mesh->Elements->Nodes[0]), num_Elements, num_Elements_numNodes)) )
         throw DataException("Error - MeshAdapter::dump: copy Elements_Nodes to netCDF buffer failed: " + *newFileName);

   }

   // // // // // Face_Elements // // // // //

   if (num_FaceElements>0) {

      // FaceElements_Id
      if (! ( ids = dataFile.add_var("FaceElements_Id", ncInt, ncdims[4])) )
         throw DataException("Error - MeshAdapter::dump: appending FaceElements_Id to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->FaceElements->Id[0];
      if (! (ids->put(int_ptr, num_FaceElements)) )
         throw DataException("Error - MeshAdapter::dump: copy FaceElements_Id to netCDF buffer failed: " + *newFileName);

      // FaceElements_Tag
      if (! ( ids = dataFile.add_var("FaceElements_Tag", ncInt, ncdims[4])) )
         throw DataException("Error - MeshAdapter::dump: appending FaceElements_Tag to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->FaceElements->Tag[0];
      if (! (ids->put(int_ptr, num_FaceElements)) )
         throw DataException("Error - MeshAdapter::dump: copy FaceElements_Tag to netCDF buffer failed: " + *newFileName);

      // FaceElements_Owner
      if (! ( ids = dataFile.add_var("FaceElements_Owner", ncInt, ncdims[4])) )
         throw DataException("Error - MeshAdapter::dump: appending FaceElements_Owner to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->FaceElements->Owner[0];
      if (! (ids->put(int_ptr, num_FaceElements)) )
         throw DataException("Error - MeshAdapter::dump: copy FaceElements_Owner to netCDF buffer failed: " + *newFileName);

      // FaceElements_Color
      if (! ( ids = dataFile.add_var("FaceElements_Color", ncInt, ncdims[4])) )
         throw DataException("Error - MeshAdapter::dump: appending FaceElements_Color to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->FaceElements->Color[0];
      if (! (ids->put(int_ptr, num_FaceElements)) )
         throw DataException("Error - MeshAdapter::dump: copy FaceElements_Color to netCDF buffer failed: " + *newFileName);

      // FaceElements_Nodes
      if (! ( ids = dataFile.add_var("FaceElements_Nodes", ncInt, ncdims[4], ncdims[8]) ) )
         throw DataException("Error - MeshAdapter::dump: appending FaceElements_Nodes to netCDF file failed: " + *newFileName);
      if (! (ids->put(&(mesh->FaceElements->Nodes[0]), num_FaceElements, num_FaceElements_numNodes)) )
         throw DataException("Error - MeshAdapter::dump: copy FaceElements_Nodes to netCDF buffer failed: " + *newFileName);

   }

   // // // // // Contact_Elements // // // // //

   if (num_ContactElements>0) {

      // ContactElements_Id
      if (! ( ids = dataFile.add_var("ContactElements_Id", ncInt, ncdims[5])) )
         throw DataException("Error - MeshAdapter::dump: appending ContactElements_Id to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->ContactElements->Id[0];
      if (! (ids->put(int_ptr, num_ContactElements)) )
         throw DataException("Error - MeshAdapter::dump: copy ContactElements_Id to netCDF buffer failed: " + *newFileName);

      // ContactElements_Tag
      if (! ( ids = dataFile.add_var("ContactElements_Tag", ncInt, ncdims[5])) )
         throw DataException("Error - MeshAdapter::dump: appending ContactElements_Tag to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->ContactElements->Tag[0];
      if (! (ids->put(int_ptr, num_ContactElements)) )
         throw DataException("Error - MeshAdapter::dump: copy ContactElements_Tag to netCDF buffer failed: " + *newFileName);

      // ContactElements_Owner
      if (! ( ids = dataFile.add_var("ContactElements_Owner", ncInt, ncdims[5])) )
         throw DataException("Error - MeshAdapter::dump: appending ContactElements_Owner to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->ContactElements->Owner[0];
      if (! (ids->put(int_ptr, num_ContactElements)) )
         throw DataException("Error - MeshAdapter::dump: copy ContactElements_Owner to netCDF buffer failed: " + *newFileName);

      // ContactElements_Color
      if (! ( ids = dataFile.add_var("ContactElements_Color", ncInt, ncdims[5])) )
         throw DataException("Error - MeshAdapter::dump: appending ContactElements_Color to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->ContactElements->Color[0];
      if (! (ids->put(int_ptr, num_ContactElements)) )
         throw DataException("Error - MeshAdapter::dump: copy ContactElements_Color to netCDF buffer failed: " + *newFileName);

      // ContactElements_Nodes
      if (! ( ids = dataFile.add_var("ContactElements_Nodes", ncInt, ncdims[5], ncdims[9]) ) )
         throw DataException("Error - MeshAdapter::dump: appending ContactElements_Nodes to netCDF file failed: " + *newFileName);
      if (! (ids->put(&(mesh->ContactElements->Nodes[0]), num_ContactElements, num_ContactElements_numNodes)) )
         throw DataException("Error - MeshAdapter::dump: copy ContactElements_Nodes to netCDF buffer failed: " + *newFileName);

   }

   // // // // // Points // // // // //

   if (num_Points>0) {

      fprintf(stderr, "\n\n\nWARNING: MeshAdapter::dump has not been tested with Point elements\n\n\n");

      // Points_Id
      if (! ( ids = dataFile.add_var("Points_Id", ncInt, ncdims[6])) )
         throw DataException("Error - MeshAdapter::dump: appending Points_Id to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->Points->Id[0];
      if (! (ids->put(int_ptr, num_Points)) )
         throw DataException("Error - MeshAdapter::dump: copy Points_Id to netCDF buffer failed: " + *newFileName);

      // Points_Tag
      if (! ( ids = dataFile.add_var("Points_Tag", ncInt, ncdims[6])) )
         throw DataException("Error - MeshAdapter::dump: appending Points_Tag to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->Points->Tag[0];
      if (! (ids->put(int_ptr, num_Points)) )
         throw DataException("Error - MeshAdapter::dump: copy Points_Tag to netCDF buffer failed: " + *newFileName);

      // Points_Owner
      if (! ( ids = dataFile.add_var("Points_Owner", ncInt, ncdims[6])) )
         throw DataException("Error - MeshAdapter::dump: appending Points_Owner to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->Points->Owner[0];
      if (! (ids->put(int_ptr, num_Points)) )
         throw DataException("Error - MeshAdapter::dump: copy Points_Owner to netCDF buffer failed: " + *newFileName);

      // Points_Color
      if (! ( ids = dataFile.add_var("Points_Color", ncInt, ncdims[6])) )
         throw DataException("Error - MeshAdapter::dump: appending Points_Color to netCDF file failed: " + *newFileName);
      int_ptr = &mesh->Points->Color[0];
      if (! (ids->put(int_ptr, num_Points)) )
         throw DataException("Error - MeshAdapter::dump: copy Points_Color to netCDF buffer failed: " + *newFileName);

      // Points_Nodes
      // mesh->Nodes->Id[mesh->Points->Nodes[INDEX2(0,i,1)]]
      if (! ( ids = dataFile.add_var("Points_Nodes", ncInt, ncdims[6]) ) )
         throw DataException("Error - MeshAdapter::dump: appending Points_Nodes to netCDF file failed: " + *newFileName);
      if (! (ids->put(&(mesh->Points->Nodes[0]), num_Points)) )
         throw DataException("Error - MeshAdapter::dump: copy Points_Nodes to netCDF buffer failed: " + *newFileName);

   }

   // // // // // TagMap // // // // //

   if (num_Tags>0) {

      // Temp storage to gather node IDs
      int *Tags_keys = TMPMEMALLOC(num_Tags, int);
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
         throw DataException("Error - MeshAdapter::dump: appending Tags_keys to netCDF file failed: " + *newFileName);
      int_ptr = &Tags_keys[0];
      if (! (ids->put(int_ptr, num_Tags)) )
         throw DataException("Error - MeshAdapter::dump: copy Tags_keys to netCDF buffer failed: " + *newFileName);

      // Tags_names_*
      // This is an array of strings, it should be stored as an array but instead I have hacked in one attribute per string
      // because the NetCDF manual doesn't tell how to do an array of strings
      tag_map = mesh->TagMap;
      if (tag_map) {
         int i = 0;
         while (tag_map) {
            sprintf(name_temp, "Tags_name_%d", i);
            if (!dataFile.add_att(name_temp, tag_map->name) )
               throw DataException("Error - MeshAdapter::dump: appending Tags_names_ to NetCDF file failed: " + *newFileName);
            tag_map=tag_map->next;
            i++;
         }
      }

      TMPMEMFREE(Tags_keys);

   }

/* Send token to next MPI process so he can take his turn */
#ifdef PASO_MPI
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

int MeshAdapter::getDiracDeltaFunctionCode() const
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
      numDataPointsPerSample=mesh->Elements->ReferenceElement->numQuadNodes;
   }
   break;
   case(ReducedElements):
   if (mesh->Elements!=NULL) {
      numSamples=mesh->Elements->numElements;
      numDataPointsPerSample=mesh->Elements->ReferenceElementReducedOrder->numQuadNodes;
   }
   break;
   case(FaceElements):
   if (mesh->FaceElements!=NULL) {
      numDataPointsPerSample=mesh->FaceElements->ReferenceElement->numQuadNodes;
      numSamples=mesh->FaceElements->numElements;
   }
   break;
   case(ReducedFaceElements):
   if (mesh->FaceElements!=NULL) {
      numDataPointsPerSample=mesh->FaceElements->ReferenceElementReducedOrder->numQuadNodes;
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
      numDataPointsPerSample=mesh->ContactElements->ReferenceElement->numQuadNodes;
      numSamples=mesh->ContactElements->numElements;
   }
   break;
   case(ReducedContactElementsZero):
   if (mesh->ContactElements!=NULL) {
      numDataPointsPerSample=mesh->ContactElements->ReferenceElementReducedOrder->numQuadNodes;
      numSamples=mesh->ContactElements->numElements;
   }
   break;
   case(ContactElementsOne):
   if (mesh->ContactElements!=NULL) {
      numDataPointsPerSample=mesh->ContactElements->ReferenceElement->numQuadNodes;
      numSamples=mesh->ContactElements->numElements;
   }
   break;
   case(ReducedContactElementsOne):
   if (mesh->ContactElements!=NULL) {
      numDataPointsPerSample=mesh->ContactElements->ReferenceElementReducedOrder->numQuadNodes;
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
                                 SystemMatrixAdapter& mat, escript::Data& rhs,
                                 const escript::Data& A, const escript::Data& B, const escript::Data& C,const  escript::Data& D,const  escript::Data& X,const  escript::Data& Y,
                                 const escript::Data& d, const escript::Data& y, 
                                 const escript::Data& d_contact,const escript::Data& y_contact) const
{
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

   Finley_Mesh* mesh=m_finleyMesh.get();

   Finley_Assemble_PDE(mesh->Nodes,mesh->Elements,mat.getPaso_SystemMatrix(), &_rhs, &_A, &_B, &_C, &_D, &_X, &_Y );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->FaceElements, mat.getPaso_SystemMatrix(), &_rhs, 0, 0, 0, &_d, 0, &_y );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->ContactElements, mat.getPaso_SystemMatrix(), &_rhs , 0, 0, 0, &_d_contact, 0, &_y_contact );
   checkFinleyError();
}

void  MeshAdapter::addPDEToLumpedSystem(
                                        escript::Data& mat,
                                        const escript::Data& D,
                                        const escript::Data& d) const
{
   escriptDataC _mat=mat.getDataC();
   escriptDataC _D=D.getDataC();
   escriptDataC _d=d.getDataC();

   Finley_Mesh* mesh=m_finleyMesh.get();

   Finley_Assemble_LumpedSystem(mesh->Nodes,mesh->Elements,&_mat, &_D);
   Finley_Assemble_LumpedSystem(mesh->Nodes,mesh->FaceElements,&_mat, &_d);

   checkFinleyError();
}


//
// adds linear PDE of second order into the right hand side only
//
void MeshAdapter::addPDEToRHS( escript::Data& rhs, const  escript::Data& X,const  escript::Data& Y, const escript::Data& y, const escript::Data& y_contact) const
{
   Finley_Mesh* mesh=m_finleyMesh.get();

   escriptDataC _rhs=rhs.getDataC();
   escriptDataC _X=X.getDataC();
   escriptDataC _Y=Y.getDataC();
   escriptDataC _y=y.getDataC();
   escriptDataC _y_contact=y_contact.getDataC();

   Finley_Assemble_PDE(mesh->Nodes,mesh->Elements, 0, &_rhs, 0, 0, 0, 0, &_X, &_Y );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->FaceElements, 0, &_rhs, 0, 0, 0, 0, 0, &_y );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->ContactElements, 0, &_rhs , 0, 0, 0, 0, 0, &_y_contact );
   checkFinleyError();
}
//
// adds PDE of second order into a transport problem
//
void MeshAdapter::addPDEToTransportProblem(
                                           TransportProblemAdapter& tp, escript::Data& source, const escript::Data& M,
                                           const escript::Data& A, const escript::Data& B, const escript::Data& C,
                                           const  escript::Data& D,const  escript::Data& X,const  escript::Data& Y,
                                           const escript::Data& d, const escript::Data& y, 
                                           const escript::Data& d_contact,const escript::Data& y_contact) const
{
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

   Finley_Mesh* mesh=m_finleyMesh.get();
   Paso_FCTransportProblem* _tp = tp.getPaso_FCTransportProblem();

   Finley_Assemble_PDE(mesh->Nodes,mesh->Elements,_tp->mass_matrix, &_source, 0, 0, 0, &_M, 0, 0 );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->Elements,_tp->transport_matrix, &_source, &_A, &_B, &_C, &_D, &_X, &_Y );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->FaceElements, _tp->transport_matrix, &_source, 0, 0, 0, &_d, 0, &_y );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->ContactElements, _tp->transport_matrix, &_source , 0, 0, 0, &_d_contact, 0, &_y_contact );
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
      escript::Data temp=escript::Data( in,  continuousFunction(asAbstractContinuousDomain()) );
      escriptDataC _in2 = temp.getDataC();
      Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&_in2,&_target);
   } else {
      Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&_in,&_target);
   }
   break;
   case(FaceElements):
   case(ReducedFaceElements):
   if (getMPISize()>1) {
      escript::Data temp=escript::Data( in,  continuousFunction(asAbstractContinuousDomain()) );
      escriptDataC _in2 = temp.getDataC();
      Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&_in2,&_target);

   } else {
      Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&_in,&_target);
   }
   break;
   case(Points):
   if (getMPISize()>1) {
      escript::Data temp=escript::Data( in,  continuousFunction(asAbstractContinuousDomain()) );
      escriptDataC _in2 = temp.getDataC();
   } else {
      Finley_Assemble_interpolate(mesh->Nodes,mesh->Points,&_in,&_target);
   }
   break;
   case(ContactElementsZero):
   case(ContactElementsOne):
   case(ReducedContactElementsZero):
   case(ReducedContactElementsOne):
   if (getMPISize()>1) {
      escript::Data temp=escript::Data( in,  continuousFunction(asAbstractContinuousDomain()) );
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
      escript::Data temp=escript::Data( in,  reducedContinuousFunction(asAbstractContinuousDomain()) );
      escriptDataC _in2 = temp.getDataC();
      Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&_in2,&_target);
   } else {
      Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&_in,&_target);
   }
   break;
   case(FaceElements):
   case(ReducedFaceElements):
   if (getMPISize()>1) {
      escript::Data temp=escript::Data( in,  reducedContinuousFunction(asAbstractContinuousDomain()) );
      escriptDataC _in2 = temp.getDataC();
      Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&_in2,&_target);
   } else {
      Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&_in,&_target);
   }
   break;
   case(Points):
   if (getMPISize()>1) {
      escript::Data temp=escript::Data( in,  reducedContinuousFunction(asAbstractContinuousDomain()) );
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
      escript::Data temp=escript::Data( in,  reducedContinuousFunction(asAbstractContinuousDomain()) );
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
      escript::Data tmp_data=Vector(0.0,continuousFunction(asAbstractContinuousDomain()),true);
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
void MeshAdapter::setToIntegrals(std::vector<double>& integrals,const escript::Data& arg) const
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
   temp=escript::Data( arg, escript::function(asAbstractContinuousDomain()) );
   _temp=temp.getDataC();
   Finley_Assemble_integrate(mesh->Nodes,mesh->Elements,&_temp,&integrals[0]);
   break;
   case(ReducedNodes):
   temp=escript::Data( arg, escript::function(asAbstractContinuousDomain()) );
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
   temp=escript::Data( arg, escript::function(asAbstractContinuousDomain()) );
   _temp=temp.getDataC();
   Finley_Assemble_integrate(mesh->Nodes,mesh->Elements,&_temp,&integrals[0]);
   break;
   case(ReducedDegreesOfFreedom):
   temp=escript::Data( arg, escript::function(asAbstractContinuousDomain()) );
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
         temp=escript::Data( arg,  continuousFunction(asAbstractContinuousDomain()) );
         nodeDataC = temp.getDataC();
      } else if( arg.getFunctionSpace().getTypeCode() == ReducedDegreesOfFreedom ) {
         temp=escript::Data( arg,  reducedContinuousFunction(asAbstractContinuousDomain()) );
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
   tmp = new_x.getDataC();
   Finley_Mesh_setCoordinates(mesh,&tmp);
   checkFinleyError();
}

//
// Helper for the save* methods. Extracts optional data variable names and the
// corresponding pointers from python dictionary. Caller must free arrays.
//
void MeshAdapter::extractArgsFromDict(const boost::python::dict& arg, int& numData, char**& names, escriptDataC*& data, escriptDataC**& dataPtr) const
{
   numData = boost::python::extract<int>(arg.attr("__len__")());
   /* win32 refactor */
   names = (numData>0) ? TMPMEMALLOC(numData, char*) : (char**)NULL;
   data = (numData>0) ? TMPMEMALLOC(numData,escriptDataC) : (escriptDataC*)NULL;
   dataPtr = (numData>0) ? TMPMEMALLOC(numData,escriptDataC*) : (escriptDataC**)NULL;

   boost::python::list keys=arg.keys();
   for (int i=0; i<numData; ++i) {
      std::string n=boost::python::extract<std::string>(keys[i]);
      escript::Data& d=boost::python::extract<escript::Data&>(arg[keys[i]]);
      if (dynamic_cast<const MeshAdapter&>(*(d.getFunctionSpace().getDomain())) !=*this) 
         throw FinleyAdapterException("Error: Data must be defined on same Domain");
      data[i] = d.getDataC();
      dataPtr[i] = &(data[i]);
      names[i] = TMPMEMALLOC(n.length()+1, char);
      strcpy(names[i], n.c_str());
   }
}

//
// saves mesh and optionally data arrays in openDX format
//
void MeshAdapter::saveDX(const std::string& filename,const boost::python::dict& arg) const
{
   int num_data;
   char **names;
   escriptDataC *data;
   escriptDataC **ptr_data;

   extractArgsFromDict(arg, num_data, names, data, ptr_data);
   Finley_Mesh_saveDX(filename.c_str(), m_finleyMesh.get(), num_data, names, ptr_data);
   checkFinleyError();
 
   /* win32 refactor */
   TMPMEMFREE(data);
   TMPMEMFREE(ptr_data);
   for(int i=0; i<num_data; i++)
   {
      TMPMEMFREE(names[i]);
   }
   TMPMEMFREE(names);

   return;
}

//
// saves mesh and optionally data arrays in VTK format
//
void MeshAdapter::saveVTK(const std::string& filename,const boost::python::dict& arg) const
{
   int num_data;
   char **names;
   escriptDataC *data;
   escriptDataC **ptr_data;

   extractArgsFromDict(arg, num_data, names, data, ptr_data);
   Finley_Mesh_saveVTK(filename.c_str(), m_finleyMesh.get(), num_data, names, ptr_data);
   checkFinleyError();

   /* win32 refactor */
   TMPMEMFREE(data);
   TMPMEMFREE(ptr_data);
   for(int i=0; i<num_data; i++)
   {
      TMPMEMFREE(names[i]);
   }
   TMPMEMFREE(names);
}

//
// creates a SystemMatrixAdapter stiffness matrix an initializes it with zeros
//
SystemMatrixAdapter MeshAdapter::newSystemMatrix(
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
      fsystemMatrix=Paso_SystemMatrix_alloc(type,fsystemMatrixPattern,row_blocksize,column_blocksize);
   }
   checkPasoError();
   Paso_SystemMatrixPattern_free(fsystemMatrixPattern);
   return SystemMatrixAdapter(fsystemMatrix,row_blocksize,row_functionspace,column_blocksize,column_functionspace);
}

//
// creates a TransportProblemAdapter
//
TransportProblemAdapter MeshAdapter::newTransportProblem(
                                                         const double theta,
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
   Paso_FCTransportProblem* transportProblem;
   transportProblem=Paso_FCTransportProblem_alloc(theta,fsystemMatrixPattern,blocksize);
   checkPasoError();
   Paso_SystemMatrixPattern_free(fsystemMatrixPattern);
   return TransportProblemAdapter(transportProblem,theta,blocksize,functionspace);
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
   int out=Paso_SystemMatrix_getSystemMatrixTypeId(SystemMatrixAdapter::mapOptionToPaso(solver),SystemMatrixAdapter::mapOptionToPaso(preconditioner), SystemMatrixAdapter::mapOptionToPaso(package),symmetry?1:0);
   checkPasoError();
   return out;
}
int MeshAdapter::getTransportTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const
{
   int out=Paso_FCTransportProblem_getTypeId(SystemMatrixAdapter::mapOptionToPaso(solver),SystemMatrixAdapter::mapOptionToPaso(preconditioner), SystemMatrixAdapter::mapOptionToPaso(package),symmetry?1:0);
   checkPasoError();
   return out;
}

escript::Data MeshAdapter::getX() const
{
   return continuousFunction(asAbstractContinuousDomain()).getX();
}

escript::Data MeshAdapter::getNormal() const
{
   return functionOnBoundary(asAbstractContinuousDomain()).getNormal();
}

escript::Data MeshAdapter::getSize() const
{
   return escript::function(asAbstractContinuousDomain()).getSize();
}

int* MeshAdapter::borrowSampleReferenceIDs(int functionSpaceType) const
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

void MeshAdapter::setTagMap(const std::string& name,  int tag)
{
   Finley_Mesh* mesh=m_finleyMesh.get();
   Finley_Mesh_addTagMap(mesh, name.c_str(),tag);
   checkPasoError();
   // throwStandardException("MeshAdapter::set TagMap is not implemented.");
}

int MeshAdapter::getTag(const std::string& name) const
{
   Finley_Mesh* mesh=m_finleyMesh.get();
   int tag=0;
   tag=Finley_Mesh_getTag(mesh, name.c_str());
   checkPasoError();
   // throwStandardException("MeshAdapter::getTag is not implemented.");
   return tag;
}

bool MeshAdapter::isValidTagName(const std::string& name) const
{
   Finley_Mesh* mesh=m_finleyMesh.get();
   return Finley_Mesh_isValidTagName(mesh,name.c_str());
}

std::string MeshAdapter::showTagNames() const
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
int* MeshAdapter::borrowListOfTagsInUse(int functionSpaceCode) const
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


}  // end of namespace
