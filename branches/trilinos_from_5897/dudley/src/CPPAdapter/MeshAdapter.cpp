
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

#include "MeshAdapter.h"
#include "esysUtils/blocktimer.h"
#include "esysUtils/EsysRandom.h"
#include <escript/Data.h>
#include <escript/DataFactory.h>
#include <escript/SolverOptions.h>

#include <paso/SystemMatrix.h>
#include <paso/Transport.h>

#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif
#include <boost/python/tuple.hpp>

using namespace std;
using namespace escript;
namespace bp = boost::python;

namespace dudley {

//
// define the static constants
MeshAdapter::FunctionSpaceNamesMapType MeshAdapter::m_functionSpaceTypeNames;
const int MeshAdapter::DegreesOfFreedom=DUDLEY_DEGREES_OF_FREEDOM;
const int MeshAdapter::ReducedDegreesOfFreedom=DUDLEY_REDUCED_DEGREES_OF_FREEDOM;
const int MeshAdapter::Nodes=DUDLEY_NODES;
const int MeshAdapter::ReducedNodes=DUDLEY_REDUCED_NODES;
const int MeshAdapter::Elements=DUDLEY_ELEMENTS;
const int MeshAdapter::ReducedElements=DUDLEY_REDUCED_ELEMENTS;
const int MeshAdapter::FaceElements=DUDLEY_FACE_ELEMENTS;
const int MeshAdapter::ReducedFaceElements=DUDLEY_REDUCED_FACE_ELEMENTS;
const int MeshAdapter::Points=DUDLEY_POINTS;

MeshAdapter::MeshAdapter(Dudley_Mesh* dudleyMesh)
{
   setFunctionSpaceTypeNames();
   //
   // need to use a null_deleter as Dudley_Mesh_free deletes the pointer
   // for us.
   m_dudleyMesh.reset(dudleyMesh,null_deleter());
}

//
// The copy constructor should just increment the use count
MeshAdapter::MeshAdapter(const MeshAdapter& in):
m_dudleyMesh(in.m_dudleyMesh)
{
   setFunctionSpaceTypeNames();
}

MeshAdapter::~MeshAdapter()
{
   //
   // I hope the case for the pointer being zero has been taken care of.
   //  cout << "In MeshAdapter destructor." << endl;
   if (m_dudleyMesh.unique()) {
      Dudley_Mesh_free(m_dudleyMesh.get());
   }
}

int MeshAdapter::getMPISize() const
{
   return m_dudleyMesh.get()->MPIInfo->size;
}
int MeshAdapter::getMPIRank() const
{
   return m_dudleyMesh.get()->MPIInfo->rank;
}
void MeshAdapter::MPIBarrier() const
{
#ifdef ESYS_MPI
   MPI_Barrier(m_dudleyMesh.get()->MPIInfo->comm);
#endif
   return;
}
bool MeshAdapter::onMasterProcessor() const
{
   return m_dudleyMesh.get()->MPIInfo->rank == 0;
}

MPI_Comm MeshAdapter::getMPIComm() const
{
    return m_dudleyMesh->MPIInfo->comm;
}


Dudley_Mesh* MeshAdapter::getDudley_Mesh() const
{
   return m_dudleyMesh.get();
}

void MeshAdapter::write(const string& fileName) const
{
   char *fName = (fileName.size()+1>0) ? new char[fileName.size()+1] : (char*)NULL;
   strcpy(fName,fileName.c_str());
   Dudley_Mesh_write(m_dudleyMesh.get(),fName);
   checkDudleyError();
   delete[] fName;
}

void MeshAdapter::Print_Mesh_Info(const bool full) const
{
   Dudley_PrintMesh_Info(m_dudleyMesh.get(), full);
}

void MeshAdapter::dump(const string& fileName) const
{
#ifdef USE_NETCDF
   const NcDim* ncdims[12] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
   NcVar *ids;
   int *int_ptr;
   Dudley_Mesh *mesh = m_dudleyMesh.get();
   Dudley_TagMap* tag_map;
   int num_Tags = 0;
   int mpi_size                         = mesh->MPIInfo->size;
   int mpi_rank                         = mesh->MPIInfo->rank;
   int numDim                           = mesh->Nodes->numDim;
   int numNodes                         = mesh->Nodes->numNodes;
   int num_Elements                     = mesh->Elements->numElements;
   int num_FaceElements                 = mesh->FaceElements->numElements;
   int num_Points                       = mesh->Points->numElements;
   int num_Elements_numNodes            = mesh->Elements->numNodes;
   int num_FaceElements_numNodes        = mesh->FaceElements->numNodes;
#ifdef ESYS_MPI
   MPI_Status status;
#endif

/* Incoming token indicates it's my turn to write */
#ifdef ESYS_MPI
   if (mpi_rank>0) MPI_Recv(&num_Tags, 0, MPI_INT, mpi_rank-1, 81800, mesh->MPIInfo->comm, &status);
#endif

   string newFileName(esysUtils::appendRankToFileName(
                                            fileName, mpi_size, mpi_rank));

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
   NcFile dataFile(newFileName.c_str(), NcFile::Replace);
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
   if (num_Points>0)
      if (! (ncdims[6] = dataFile.add_dim("dim_Points", num_Points)) )
         throw DataException(msgPrefix+"add_dim(dim_Points)");
   if (num_Elements>0)
      if (! (ncdims[7] = dataFile.add_dim("dim_Elements_Nodes", num_Elements_numNodes)) )
         throw DataException(msgPrefix+"add_dim(dim_Elements_Nodes)");
   if (num_FaceElements>0)
      if (! (ncdims[8] = dataFile.add_dim("dim_FaceElements_numNodes", num_FaceElements_numNodes)) )
         throw DataException(msgPrefix+"add_dim(dim_FaceElements_numNodes)");
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
   if (!dataFile.add_att("num_Points",num_Points) )
      throw DataException(msgPrefix+"add_att(num_Points)");
   if (!dataFile.add_att("num_Elements_numNodes",num_Elements_numNodes) )
      throw DataException(msgPrefix+"add_att(num_Elements_numNodes)");
   if (!dataFile.add_att("num_FaceElements_numNodes",num_FaceElements_numNodes) )
      throw DataException(msgPrefix+"add_att(num_FaceElements_numNodes)");
   if (!dataFile.add_att("Elements_TypeId", mesh->Elements->etype) )
      throw DataException(msgPrefix+"add_att(Elements_TypeId)");
   if (!dataFile.add_att("FaceElements_TypeId", mesh->FaceElements->etype) )
      throw DataException(msgPrefix+"add_att(FaceElements_TypeId)");
   if (!dataFile.add_att("Points_TypeId", mesh->Points->etype) )
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
      int *Tags_keys = new int[num_Tags];
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
   Dudley_setError(IO_ERROR, "MeshAdapter::dump: not configured with NetCDF. Please contact your installation manager.");
#endif  /* USE_NETCDF */
   checkDudleyError();
}

string MeshAdapter::getDescription() const
{
   return "DudleyMesh";
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
   (FunctionSpaceNamesMapType::value_type(DegreesOfFreedom,"Dudley_DegreesOfFreedom [Solution(domain)]"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(ReducedDegreesOfFreedom,"Dudley_ReducedDegreesOfFreedom [ReducedSolution(domain)]"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(Nodes,"Dudley_Nodes [ContinuousFunction(domain)]"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(ReducedNodes,"Dudley_Reduced_Nodes [ReducedContinuousFunction(domain)]"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(Elements,"Dudley_Elements [Function(domain)]"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(ReducedElements,"Dudley_Reduced_Elements [ReducedFunction(domain)]"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(FaceElements,"Dudley_Face_Elements [FunctionOnBoundary(domain)]"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(ReducedFaceElements,"Dudley_Reduced_Face_Elements [ReducedFunctionOnBoundary(domain)]"));
   m_functionSpaceTypeNames.insert
   (FunctionSpaceNamesMapType::value_type(Points,"Dudley_Points [DiracDeltaFunctions(domain)]"));
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
   throw DudleyAdapterException("Dudley does not support contact elements.");
}

int MeshAdapter::getReducedFunctionOnContactZeroCode() const
{
   throw DudleyAdapterException("Dudley does not support contact elements.");
}

int MeshAdapter::getFunctionOnContactOneCode() const
{
   throw DudleyAdapterException("Dudley does not support contact elements.");
}

int MeshAdapter::getReducedFunctionOnContactOneCode() const
{
   throw DudleyAdapterException("Dudley does not support contact elements.");
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
   Dudley_Mesh* mesh=m_dudleyMesh.get();
   int numDim=Dudley_Mesh_getDim(mesh);
   checkDudleyError();
   return numDim;
}

//
// Return the number of data points summed across all MPI processes
//
int MeshAdapter::getNumDataPointsGlobal() const
{
   return Dudley_NodeFile_getGlobalNumNodes(m_dudleyMesh.get()->Nodes);
}

//
// return the number of data points per sample and the number of samples
// needed to represent data on a parts of the mesh.
//
pair<int,int> MeshAdapter::getDataShape(int functionSpaceCode) const
{
   int numDataPointsPerSample=0;
   int numSamples=0;
   Dudley_Mesh* mesh=m_dudleyMesh.get();
   switch (functionSpaceCode) {
   case(Nodes):
   numDataPointsPerSample=1;
   numSamples=Dudley_NodeFile_getNumNodes(mesh->Nodes);
   break;
   case(ReducedNodes):
   numDataPointsPerSample=1;
   numSamples=Dudley_NodeFile_getNumReducedNodes(mesh->Nodes);
   break;
   case(Elements):
   if (mesh->Elements!=NULL) {
      numSamples=mesh->Elements->numElements;
      numDataPointsPerSample=mesh->Elements->numLocalDim+1/*referenceElementSet->referenceElement->BasisFunctions->numQuadNodes*/;
   }
   break;
   case(ReducedElements):
   if (mesh->Elements!=NULL) {
      numSamples=mesh->Elements->numElements;
      numDataPointsPerSample=(mesh->Elements->numLocalDim==0)?0:1;
   }
   break;
   case(FaceElements):
   if (mesh->FaceElements!=NULL) {
      numDataPointsPerSample=mesh->FaceElements->numLocalDim+1/*referenceElementSet->referenceElement->BasisFunctions->numQuadNodes*/;
      numSamples=mesh->FaceElements->numElements;
   }
   break;
   case(ReducedFaceElements):
   if (mesh->FaceElements!=NULL) {
      numDataPointsPerSample=(mesh->FaceElements->numLocalDim==0)?0:1/*referenceElementSet->referenceElementReducedQuadrature->BasisFunctions->numQuadNodes*/;
      numSamples=mesh->FaceElements->numElements;
   }
   break;
   case(Points):
   if (mesh->Points!=NULL) {
      numDataPointsPerSample=1;
      numSamples=mesh->Points->numElements;
   }
   break;
   case(DegreesOfFreedom):
   if (mesh->Nodes!=NULL) {
      numDataPointsPerSample=1;
      numSamples=Dudley_NodeFile_getNumDegreesOfFreedom(mesh->Nodes);
   }
   break;
   case(ReducedDegreesOfFreedom):
   if (mesh->Nodes!=NULL) {
      numDataPointsPerSample=1;
      numSamples=Dudley_NodeFile_getNumReducedDegreesOfFreedom(mesh->Nodes);
   }
   break;
   default:
      stringstream temp;
      temp << "Error - Invalid function space type: " << functionSpaceCode << " for domain: " << getDescription();
      throw DudleyAdapterException(temp.str());
      break;
   }
   return pair<int,int>(numDataPointsPerSample,numSamples);
}

//
// adds linear PDE of second order into a given stiffness matrix and right hand side:
//
void MeshAdapter::addPDEToSystem(AbstractSystemMatrix& mat, escript::Data& rhs,
                                 const escript::Data& A, const escript::Data& B,
                                 const escript::Data& C, const escript::Data& D,
                                 const escript::Data& X, const escript::Data& Y,
                                 const escript::Data& d, const escript::Data& y,
                                 const escript::Data& d_contact,
                                 const escript::Data& y_contact,
                                 const escript::Data& d_dirac,
                                 const escript::Data& y_dirac) const
{
    if (!d_contact.isEmpty() || !y_contact.isEmpty())
        throw DudleyAdapterException("Dudley does not support contact elements");

    paso::SystemMatrix* smat = dynamic_cast<paso::SystemMatrix*>(&mat);
    if (smat) {
        paso::SystemMatrix_ptr S(smat->shared_from_this());
        Dudley_Mesh* mesh = m_dudleyMesh.get();

        Dudley_Assemble_PDE(mesh->Nodes, mesh->Elements, S, &rhs, &A, &B, &C, &D, &X, &Y);
        checkDudleyError();

        Dudley_Assemble_PDE(mesh->Nodes, mesh->FaceElements, S, &rhs, 0, 0, 0, &d, 0, &y);
        checkDudleyError();

        Dudley_Assemble_PDE(mesh->Nodes,mesh->Points, S, &rhs, 0, 0, 0, &d_dirac, 0, &y_dirac);
        checkDudleyError();
        return;
    }
    throw DudleyAdapterException("Dudley only accepts Paso system matrices");
}

void  MeshAdapter::addPDEToLumpedSystem(escript::Data& mat,
                                        const escript::Data& D,
                                        const escript::Data& d,
                                        const escript::Data& d_dirac,
                                        bool useHRZ) const
{
   Dudley_Mesh* mesh=m_dudleyMesh.get();

   Dudley_Assemble_LumpedSystem(mesh->Nodes,mesh->Elements, &mat, &D, useHRZ);
   checkDudleyError();
   
   Dudley_Assemble_LumpedSystem(mesh->Nodes,mesh->FaceElements, &mat, &d, useHRZ);
   checkDudleyError();

   Dudley_Assemble_LumpedSystem(mesh->Nodes,mesh->FaceElements, &mat, &d_dirac, useHRZ);
   checkDudleyError();

}


//
// adds linear PDE of second order into the right hand side only
//
void MeshAdapter::addPDEToRHS( escript::Data& rhs, const escript::Data& X,const  escript::Data& Y, const escript::Data& y, const escript::Data& y_contact, const escript::Data& y_dirac) const
{
   if (!y_contact.isEmpty())
   {
        throw DudleyAdapterException("Dudley does not support y_contact");
   }
   Dudley_Mesh* mesh=m_dudleyMesh.get();

   Dudley_Assemble_PDE(mesh->Nodes,mesh->Elements, escript::ASM_ptr(), &rhs,
                       NULL, NULL, NULL, NULL, &X, &Y);
   checkDudleyError();

   Dudley_Assemble_PDE(mesh->Nodes,mesh->FaceElements, escript::ASM_ptr(),
                       &rhs, NULL, NULL, NULL, NULL, NULL, &y);
   checkDudleyError();

   Dudley_Assemble_PDE(mesh->Nodes,mesh->Points, escript::ASM_ptr(), &rhs,
                       NULL, NULL, NULL, NULL, NULL, &y_dirac);
   checkDudleyError();
}
//
// adds PDE of second order into a transport problem
//
void MeshAdapter::addPDEToTransportProblem(
        AbstractTransportProblem& tp, escript::Data& source,
        const escript::Data& M, const escript::Data& A, const escript::Data& B,
        const escript::Data& C, const  escript::Data& D,const escript::Data& X,
        const escript::Data& Y, const escript::Data& d, const escript::Data& y,
        const escript::Data& d_contact,const escript::Data& y_contact,
        const escript::Data& d_dirac, const escript::Data& y_dirac) const
{
    if (!d_contact.isEmpty())
    {
        throw DudleyAdapterException("Dudley does not support d_contact");
    }
    if (!y_contact.isEmpty())
    {
        throw DudleyAdapterException("Dudley does not support y_contact");
    }   
    paso::TransportProblem* ptp = dynamic_cast<paso::TransportProblem*>(&tp);
    if (!ptp)
    {
        throw DudleyAdapterException("Dudley only accepts Paso transport problems");
    }
    DataTypes::ShapeType shape;
    source.expand();

   Dudley_Mesh* mesh=m_dudleyMesh.get();

   Dudley_Assemble_PDE(mesh->Nodes, mesh->Elements, ptp->borrowMassMatrix(),
                       &source, 0, 0, 0, &M, 0, 0);
   checkDudleyError();

   Dudley_Assemble_PDE(mesh->Nodes, mesh->Elements, ptp->borrowTransportMatrix(),
                       &source, &A, &B, &C, &D, &X, &Y);
   checkDudleyError();

   Dudley_Assemble_PDE(mesh->Nodes, mesh->FaceElements,
                       ptp->borrowTransportMatrix(), &source, NULL, NULL, NULL,
                       &d, NULL, &y);
   checkDudleyError();

   Dudley_Assemble_PDE(mesh->Nodes, mesh->Points, ptp->borrowTransportMatrix(),
                       &source, NULL, NULL, NULL, &d_dirac, NULL, &y_dirac);
   checkDudleyError();
}

//
// interpolates data between different function spaces:
//
void MeshAdapter::interpolateOnDomain(escript::Data& target,const escript::Data& in) const
{
   const MeshAdapter& inDomain=dynamic_cast<const MeshAdapter&>(*(in.getFunctionSpace().getDomain()));
   const MeshAdapter& targetDomain=dynamic_cast<const MeshAdapter&>(*(target.getFunctionSpace().getDomain()));
   if (inDomain!=*this)  
      throw DudleyAdapterException("Error - Illegal domain of interpolant.");
   if (targetDomain!=*this) 
      throw DudleyAdapterException("Error - Illegal domain of interpolation target.");

   Dudley_Mesh* mesh=m_dudleyMesh.get();
   switch(in.getFunctionSpace().getTypeCode()) {
   case(Nodes):
      switch(target.getFunctionSpace().getTypeCode()) {
      case(Nodes):
      case(ReducedNodes):
      case(DegreesOfFreedom):
      case(ReducedDegreesOfFreedom):
      Dudley_Assemble_CopyNodalData(mesh->Nodes,&target,&in);
      break;
      case(Elements):
      case(ReducedElements):
      Dudley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&in,&target);
      break;
      case(FaceElements):
      case(ReducedFaceElements):
      Dudley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&in,&target);
      break;
      case(Points):
      Dudley_Assemble_interpolate(mesh->Nodes,mesh->Points,&in,&target);
      break;
      default:
         stringstream temp;
         temp << "Error - Interpolation on Domain: Dudley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
         throw DudleyAdapterException(temp.str());
         break;
      }
      break;
   case(ReducedNodes):
      switch(target.getFunctionSpace().getTypeCode()) {
      case(Nodes):
      case(ReducedNodes):
      case(DegreesOfFreedom):
      case(ReducedDegreesOfFreedom):
      Dudley_Assemble_CopyNodalData(mesh->Nodes,&target,&in);
      break;
      case(Elements):
      case(ReducedElements):
      Dudley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&in,&target);
      break;
      case(FaceElements):
      case(ReducedFaceElements):
      Dudley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&in,&target);
      break;
      case(Points):
      Dudley_Assemble_interpolate(mesh->Nodes,mesh->Points,&in,&target);
      break;
      default:
         stringstream temp;
         temp << "Error - Interpolation on Domain: Dudley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
         throw DudleyAdapterException(temp.str());
         break;
      }
      break;
   case(Elements):
      if (target.getFunctionSpace().getTypeCode()==Elements) {
         Dudley_Assemble_CopyElementData(mesh->Elements,&target,&in);
      } else if (target.getFunctionSpace().getTypeCode()==ReducedElements) {
         Dudley_Assemble_AverageElementData(mesh->Elements,&target,&in);
      } else {
         throw DudleyAdapterException("Error - No interpolation with data on elements possible.");
      }
      break;
   case(ReducedElements):
      if (target.getFunctionSpace().getTypeCode()==ReducedElements) {
         Dudley_Assemble_CopyElementData(mesh->Elements,&target,&in);
      } else {
         throw DudleyAdapterException("Error - No interpolation with data on elements with reduced integration order possible.");
      }
      break;
   case(FaceElements):
      if (target.getFunctionSpace().getTypeCode()==FaceElements) {
         Dudley_Assemble_CopyElementData(mesh->FaceElements,&target,&in);
      } else if (target.getFunctionSpace().getTypeCode()==ReducedFaceElements) {
         Dudley_Assemble_AverageElementData(mesh->FaceElements,&target,&in);
      } else {
         throw DudleyAdapterException("Error - No interpolation with data on face elements possible.");
      }
      break;
   case(ReducedFaceElements):
      if (target.getFunctionSpace().getTypeCode()==ReducedFaceElements) {
         Dudley_Assemble_CopyElementData(mesh->FaceElements,&target,&in);
      } else {
         throw DudleyAdapterException("Error - No interpolation with data on face elements with reduced integration order possible.");
      }
      break;
   case(Points):
      if (target.getFunctionSpace().getTypeCode()==Points) {
         Dudley_Assemble_CopyElementData(mesh->Points,&target,&in);
      } else {
         throw DudleyAdapterException("Error - No interpolation with data on points possible.");
      }
      break;
   case(DegreesOfFreedom):      
      switch(target.getFunctionSpace().getTypeCode()) {
      case(ReducedDegreesOfFreedom):
      case(DegreesOfFreedom):
      Dudley_Assemble_CopyNodalData(mesh->Nodes,&target,&in);
      break;
   
      case(Nodes):
      case(ReducedNodes):
      if (getMPISize()>1) {
         escript::Data temp=escript::Data(in);
         temp.expand();
         Dudley_Assemble_CopyNodalData(mesh->Nodes,&target,&temp);
      } else {
         Dudley_Assemble_CopyNodalData(mesh->Nodes,&target,&in);
      }
      break;
      case(Elements):
      case(ReducedElements):
      if (getMPISize()>1) {
         escript::Data temp=escript::Data( in,  continuousFunction(*this) );
         Dudley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&temp,&target);
      } else {
         Dudley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&in,&target);
      }
      break;
      case(FaceElements):
      case(ReducedFaceElements):
      if (getMPISize()>1) {
         escript::Data temp=escript::Data( in,  continuousFunction(*this) );
         Dudley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&temp,&target);
   
      } else {
         Dudley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&in,&target);
      }
      break;
      case(Points):
      if (getMPISize()>1) {
         //escript::Data temp=escript::Data( in,  continuousFunction(*this) );
         //escriptDataC _in2 = temp.getDataC();
      } else {
         Dudley_Assemble_interpolate(mesh->Nodes,mesh->Points,&in,&target);
      }
      break;
      default:
         stringstream temp;
         temp << "Error - Interpolation On Domain: Dudley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
         throw DudleyAdapterException(temp.str());
         break;
      }
      break;
   case(ReducedDegreesOfFreedom):
      switch(target.getFunctionSpace().getTypeCode()) {
      case(Nodes):
      throw DudleyAdapterException("Error - Dudley does not support interpolation from reduced degrees of freedom to mesh nodes.");
      break;
      case(ReducedNodes):
      if (getMPISize()>1) {
         escript::Data temp=escript::Data(in);
         temp.expand();
         Dudley_Assemble_CopyNodalData(mesh->Nodes,&target,&temp);
      } else {
         Dudley_Assemble_CopyNodalData(mesh->Nodes,&target,&in);
      }
      break;
      case(DegreesOfFreedom):
      throw DudleyAdapterException("Error - Dudley does not support interpolation from reduced degrees of freedom to degrees of freedom");
      break;
      case(ReducedDegreesOfFreedom):
      Dudley_Assemble_CopyNodalData(mesh->Nodes,&target,&in);
      break;
      case(Elements):
      case(ReducedElements):
      if (getMPISize()>1) {
         escript::Data temp=escript::Data( in,  reducedContinuousFunction(*this) );
         Dudley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&temp,&target);
      } else {
         Dudley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&in,&target);
      }
      break;
      case(FaceElements):
      case(ReducedFaceElements):
      if (getMPISize()>1) {
         escript::Data temp=escript::Data( in,  reducedContinuousFunction(*this) );
         Dudley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&temp,&target);
      } else {
         Dudley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&in,&target);
      }
      break;
      case(Points):
      if (getMPISize()>1) {
         escript::Data temp=escript::Data( in,  reducedContinuousFunction(*this) );
         Dudley_Assemble_interpolate(mesh->Nodes,mesh->Points,&temp,&target);
      } else {
         Dudley_Assemble_interpolate(mesh->Nodes,mesh->Points,&in,&target);
      }
      break;
      default:
         stringstream temp;
         temp << "Error - Interpolation On Domain: Dudley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
         throw DudleyAdapterException(temp.str());
         break;
      }
      break;
   default:
      stringstream temp;
      temp << "Error - Interpolation On Domain: Dudley does not know anything about function space type %d" << in.getFunctionSpace().getTypeCode();
      throw DudleyAdapterException(temp.str());
      break;
   }
   checkDudleyError();
}

//
// copies the locations of sample points into x:
//
void MeshAdapter::setToX(escript::Data& arg) const
{
   const MeshAdapter& argDomain=dynamic_cast<const MeshAdapter&>(*(arg.getFunctionSpace().getDomain()));
   if (argDomain!=*this) 
      throw DudleyAdapterException("Error - Illegal domain of data point locations");
   Dudley_Mesh* mesh=m_dudleyMesh.get();
   // in case of values node coordinates we can do the job directly:
   if (arg.getFunctionSpace().getTypeCode()==Nodes) {
      Dudley_Assemble_NodeCoordinates(mesh->Nodes,&arg);
   } else {
      escript::Data tmp_data=Vector(0.0,continuousFunction(*this),true);
      Dudley_Assemble_NodeCoordinates(mesh->Nodes,&tmp_data);
      // this is then interpolated onto arg:
      interpolateOnDomain(arg,tmp_data);
   }
   checkDudleyError();
}

//
// return the normal vectors at the location of data points as a Data object:
//
void MeshAdapter::setToNormal(escript::Data& normal) const
{
/*   const MeshAdapter& normalDomain=dynamic_cast<const MeshAdapter&>(normal.getFunctionSpace().getDomain());*/
   const MeshAdapter& normalDomain=dynamic_cast<const MeshAdapter&>(*(normal.getFunctionSpace().getDomain()));
   if (normalDomain!=*this) 
      throw DudleyAdapterException("Error - Illegal domain of normal locations");
   Dudley_Mesh* mesh=m_dudleyMesh.get();
   switch(normal.getFunctionSpace().getTypeCode()) {
   case(Nodes):
   throw DudleyAdapterException("Error - Dudley does not support surface normal vectors for nodes");
   break;
   case(ReducedNodes):
   throw DudleyAdapterException("Error - Dudley does not support surface normal vectors for reduced nodes");
   break;
   case(Elements):
   throw DudleyAdapterException("Error - Dudley does not support surface normal vectors for elements");
   break;
   case(ReducedElements):
   throw DudleyAdapterException("Error - Dudley does not support surface normal vectors for elements with reduced integration order");
   break;
   case (FaceElements):
   Dudley_Assemble_setNormal(mesh->Nodes,mesh->FaceElements,&normal);
   break;
   case (ReducedFaceElements):
   Dudley_Assemble_setNormal(mesh->Nodes,mesh->FaceElements,&normal);
   break;
   case(Points):
   throw DudleyAdapterException("Error - Dudley does not support surface normal vectors for point elements");
   break;
   case(DegreesOfFreedom):
   throw DudleyAdapterException("Error - Dudley does not support surface normal vectors for degrees of freedom.");
   break;
   case(ReducedDegreesOfFreedom):
   throw DudleyAdapterException("Error - Dudley does not support surface normal vectors for reduced degrees of freedom.");
   break;
   default:
      stringstream temp;
      temp << "Error - Normal Vectors: Dudley does not know anything about function space type " << normal.getFunctionSpace().getTypeCode();
      throw DudleyAdapterException(temp.str());
      break;
   }
   checkDudleyError();
}

//
// interpolates data to other domain:
//
void MeshAdapter::interpolateAcross(escript::Data& target,const escript::Data& source) const
{
   const_Domain_ptr targetDomain_p=target.getFunctionSpace().getDomain();
   const MeshAdapter* targetDomain=dynamic_cast<const MeshAdapter*>(targetDomain_p.get());
   if (targetDomain!=this) 
      throw DudleyAdapterException("Error - Illegal domain of interpolation target");

   throw DudleyAdapterException("Error - Dudley does not allow interpolation across domains yet.");
}

//
// calculates the integral of a function defined of arg:
//
void MeshAdapter::setToIntegrals(vector<double>& integrals,const escript::Data& arg) const
{
   const MeshAdapter& argDomain=dynamic_cast<const MeshAdapter&>(*(arg.getFunctionSpace().getDomain()));
   if (argDomain!=*this) 
      throw DudleyAdapterException("Error - Illegal domain of integration kernel");

   double blocktimer_start = blocktimer_time();
   Dudley_Mesh* mesh=m_dudleyMesh.get();
   escript::Data temp;
   switch(arg.getFunctionSpace().getTypeCode()) {
   case(Nodes):
   temp=escript::Data( arg, escript::function(*this) );
   Dudley_Assemble_integrate(mesh->Nodes,mesh->Elements,&temp,&integrals[0]);
   break;
   case(ReducedNodes):
   temp=escript::Data( arg, escript::function(*this) );
   Dudley_Assemble_integrate(mesh->Nodes,mesh->Elements,&temp,&integrals[0]);
   break;
   case(Elements):
   Dudley_Assemble_integrate(mesh->Nodes,mesh->Elements,&arg,&integrals[0]);
   break;
   case(ReducedElements):
   Dudley_Assemble_integrate(mesh->Nodes,mesh->Elements,&arg,&integrals[0]);
   break;
   case(FaceElements):
   Dudley_Assemble_integrate(mesh->Nodes,mesh->FaceElements,&arg,&integrals[0]);
   break;
   case(ReducedFaceElements):
   Dudley_Assemble_integrate(mesh->Nodes,mesh->FaceElements,&arg,&integrals[0]);
   break;
   case(Points):
   throw DudleyAdapterException("Error - Integral of data on points is not supported.");
   break;
   case(DegreesOfFreedom):
   temp=escript::Data( arg, escript::function(*this) );
   Dudley_Assemble_integrate(mesh->Nodes,mesh->Elements,&temp,&integrals[0]);
   break;
   case(ReducedDegreesOfFreedom):
   temp=escript::Data( arg, escript::function(*this) );
   Dudley_Assemble_integrate(mesh->Nodes,mesh->Elements,&temp,&integrals[0]);
   break;
   default:
      stringstream temp;
      temp << "Error - Integrals: Dudley does not know anything about function space type " << arg.getFunctionSpace().getTypeCode();
      throw DudleyAdapterException(temp.str());
      break;
   }
   checkDudleyError();
   blocktimer_increment("integrate()", blocktimer_start);
}

//
// calculates the gradient of arg:
//
void MeshAdapter::setToGradient(escript::Data& grad,const escript::Data& arg) const
{
   const MeshAdapter& argDomain=dynamic_cast<const MeshAdapter&>(*(arg.getFunctionSpace().getDomain()));
   if (argDomain!=*this)
      throw DudleyAdapterException("Error - Illegal domain of gradient argument");
   const MeshAdapter& gradDomain=dynamic_cast<const MeshAdapter&>(*(grad.getFunctionSpace().getDomain()));
   if (gradDomain!=*this)
      throw DudleyAdapterException("Error - Illegal domain of gradient");

   Dudley_Mesh* mesh=m_dudleyMesh.get();
   const escript::Data* nodeData=0;
   escript::Data temp;
   if (getMPISize()>1) {
      if( arg.getFunctionSpace().getTypeCode() == DegreesOfFreedom ) {
         temp=escript::Data( arg,  continuousFunction(*this) );
         nodeData = &temp;
      } else if( arg.getFunctionSpace().getTypeCode() == ReducedDegreesOfFreedom ) {
         temp=escript::Data( arg,  reducedContinuousFunction(*this) );
         nodeData = &temp;
      } else {
         nodeData = &arg;
      }
   } else {
      nodeData = &arg;
   }
   switch(grad.getFunctionSpace().getTypeCode()) {
   case(Nodes):
   throw DudleyAdapterException("Error - Gradient at nodes is not supported.");
   break;
   case(ReducedNodes):
   throw DudleyAdapterException("Error - Gradient at reduced nodes is not supported.");
   break;
   case(Elements):
   Dudley_Assemble_gradient(mesh->Nodes,mesh->Elements,&grad, nodeData);
   break;
   case(ReducedElements):
   Dudley_Assemble_gradient(mesh->Nodes,mesh->Elements,&grad, nodeData);
   break;
   case(FaceElements):
   Dudley_Assemble_gradient(mesh->Nodes,mesh->FaceElements,&grad, nodeData);
   break;
   case(ReducedFaceElements):
   Dudley_Assemble_gradient(mesh->Nodes,mesh->FaceElements,&grad, nodeData);
   break;
   case(Points):
   throw DudleyAdapterException("Error - Gradient at points is not supported.");
   break;
   case(DegreesOfFreedom):
   throw DudleyAdapterException("Error - Gradient at degrees of freedom is not supported.");
   break;
   case(ReducedDegreesOfFreedom):
   throw DudleyAdapterException("Error - Gradient at reduced degrees of freedom is not supported.");
   break;
   default:
      stringstream temp;
      temp << "Error - Gradient: Dudley does not know anything about function space type " << arg.getFunctionSpace().getTypeCode();
      throw DudleyAdapterException(temp.str());
      break;
   }
   checkDudleyError();
}

//
// returns the size of elements:
//
void MeshAdapter::setToSize(escript::Data& size) const
{
   Dudley_Mesh* mesh=m_dudleyMesh.get();
   switch(size.getFunctionSpace().getTypeCode()) {
   case(Nodes):
   throw DudleyAdapterException("Error - Size of nodes is not supported.");
   break;
   case(ReducedNodes):
   throw DudleyAdapterException("Error - Size of reduced nodes is not supported.");
   break;
   case(Elements):
   Dudley_Assemble_getSize(mesh->Nodes,mesh->Elements,&size);
   break;
   case(ReducedElements):
   Dudley_Assemble_getSize(mesh->Nodes,mesh->Elements,&size);
   break;
   case(FaceElements):
   Dudley_Assemble_getSize(mesh->Nodes,mesh->FaceElements,&size);
   break;
   case(ReducedFaceElements):
   Dudley_Assemble_getSize(mesh->Nodes,mesh->FaceElements,&size);
   break;
   case(Points):
   throw DudleyAdapterException("Error - Size of point elements is not supported.");
   break;
   case(DegreesOfFreedom):
   throw DudleyAdapterException("Error - Size of degrees of freedom is not supported.");
   break;
   case(ReducedDegreesOfFreedom):
   throw DudleyAdapterException("Error - Size of reduced degrees of freedom is not supported.");
   break;
   default:
      stringstream temp;
      temp << "Error - Element size: Dudley does not know anything about function space type " << size.getFunctionSpace().getTypeCode();
      throw DudleyAdapterException(temp.str());
      break;
   }
   checkDudleyError();
}

//
// sets the location of nodes
//
void MeshAdapter::setNewX(const escript::Data& new_x)
{
   Dudley_Mesh* mesh=m_dudleyMesh.get();
   const MeshAdapter& newDomain=dynamic_cast<const MeshAdapter&>(*(new_x.getFunctionSpace().getDomain()));
   if (newDomain!=*this) 
      throw DudleyAdapterException("Error - Illegal domain of new point locations");
   if ( new_x.getFunctionSpace() == continuousFunction(*this) ) {
       Dudley_Mesh_setCoordinates(mesh,&new_x);
   } else {
       throw DudleyAdapterException("As of version escript3.3 - SetNewX only accepts ContinuousFunction arguments please interpolate.");      

   }
   checkDudleyError();
}

bool MeshAdapter::ownSample(int fs_code, index_t id) const
{
    if (getMPISize()>1) {
#ifdef ESYS_MPI
        index_t myFirstNode=0, myLastNode=0, k=0;
        index_t* globalNodeIndex=0;
        Dudley_Mesh* mesh_p=m_dudleyMesh.get();
        if (fs_code == DUDLEY_REDUCED_NODES) 
        {
            myFirstNode = Dudley_NodeFile_getFirstReducedNode(mesh_p->Nodes);
            myLastNode = Dudley_NodeFile_getLastReducedNode(mesh_p->Nodes);
            globalNodeIndex = Dudley_NodeFile_borrowGlobalReducedNodesIndex(mesh_p->Nodes);
        }
        else if (fs_code == DUDLEY_NODES)
        {
            myFirstNode = Dudley_NodeFile_getFirstNode(mesh_p->Nodes);
            myLastNode = Dudley_NodeFile_getLastNode(mesh_p->Nodes);
            globalNodeIndex = Dudley_NodeFile_borrowGlobalNodesIndex(mesh_p->Nodes);
        }
        else
        {
            throw DudleyAdapterException("unsupported function space type for ownSample()");
        }
        k=globalNodeIndex[id];
        return static_cast<bool>( (myFirstNode <= k) && (k < myLastNode) );
#endif
    }
    return true;
}



//
// creates a stiffness matrix an initializes it with zeros
//
ASM_ptr MeshAdapter::newSystemMatrix(const int row_blocksize,
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
      throw DudleyAdapterException("Error - domain of row function space does not match the domain of matrix generator.");
   const MeshAdapter& col_domain=dynamic_cast<const MeshAdapter&>(*(column_functionspace.getDomain()));
   if (col_domain!=*this) 
      throw DudleyAdapterException("Error - domain of column function space does not match the domain of matrix generator.");
   // is the function space type right 
   if (row_functionspace.getTypeCode()==DegreesOfFreedom) {
      reduceRowOrder=0;
   } else if (row_functionspace.getTypeCode()==ReducedDegreesOfFreedom) {
      reduceRowOrder=1;
   } else {
      throw DudleyAdapterException("Error - illegal function space type for system matrix rows.");
   }
   if (column_functionspace.getTypeCode()==DegreesOfFreedom) {
      reduceColOrder=0;
   } else if (column_functionspace.getTypeCode()==ReducedDegreesOfFreedom) {
      reduceColOrder=1;
   } else {
      throw DudleyAdapterException("Error - illegal function space type for system matrix columns.");
   }
   // generate matrix:
 
   paso::SystemMatrixPattern_ptr fsystemMatrixPattern(Dudley_getPattern(getDudley_Mesh(),reduceRowOrder,reduceColOrder));
   checkDudleyError();
   paso::SystemMatrix_ptr sm;
   int trilinos = 0;
   if (trilinos) {
#ifdef TRILINOS
      /* Allocation Epetra_VrbMatrix here */
#endif
   }
   else {
      sm.reset(new paso::SystemMatrix(type, fsystemMatrixPattern,
                  row_blocksize, column_blocksize, false, row_functionspace,
                  column_functionspace));
   }
   checkPasoError();
   return sm;
}

//
// creates a TransportProblem
//
ATP_ptr MeshAdapter::newTransportProblem(int blocksize,
                                         const escript::FunctionSpace& fs,
                                         int type) const
{
   int reduceOrder=0;
   // is the domain right?
   const MeshAdapter& domain=dynamic_cast<const MeshAdapter&>(*(fs.getDomain()));
   if (domain!=*this) 
      throw DudleyAdapterException("Error - domain of function space does not match the domain of transport problem generator.");
   // is the function space type right 
   if (fs.getTypeCode()==DegreesOfFreedom) {
      reduceOrder=0;
   } else if (fs.getTypeCode()==ReducedDegreesOfFreedom) {
      reduceOrder=1;
   } else {
      throw DudleyAdapterException("Error - illegal function space type for system matrix rows.");
   }
   // generate matrix:
 
   paso::SystemMatrixPattern_ptr fsystemMatrixPattern(Dudley_getPattern(
               getDudley_Mesh(),reduceOrder,reduceOrder));
   checkDudleyError();
   paso::TransportProblem_ptr transportProblem(new paso::TransportProblem(
                                            fsystemMatrixPattern, blocksize,
                                            fs));
   checkPasoError();
   return transportProblem;
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
   case(ReducedElements):
   case(ReducedFaceElements):
   return true;
   break;
   default:
      stringstream temp;
      temp << "Error - Cell: Dudley does not know anything about function space type " << functionSpaceCode;
      throw DudleyAdapterException(temp.str());
      break;
   }
   checkDudleyError();
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
    for (int i=0;i<fs.size();++i)
    {
        switch(fs[i])
        {
        case(Nodes):   hasnodes=true;   // no break is deliberate
        case(DegreesOfFreedom):
                hasclass[1]=1;
                break;
        case(ReducedNodes):    hasrednodes=true;        // no break is deliberate
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
        default:
                return false;
        }
    }
    int totlines=hasline[0]+hasline[1]+hasline[2]+hasline[3];
    // fail if we have more than one leaf group

    if (totlines>1)
    {
        return false;   // there are at least two branches we can't interpolate between
    }
    else if (totlines==1)
    {
        if (hasline[0]==1)              // we have points
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
        else    // so we must be in line3
        {

            throw DudleyAdapterException("Programmer Error - choosing between contact elements - we should never get here.");

        }
    }
    else        // totlines==0
    {
        if (hasclass[2]==1)
        {
                // something from class 2
                resultcode=(hasrednodes?ReducedNodes:ReducedDegreesOfFreedom);
        }
        else
        {       // something from class 1
                resultcode=(hasnodes?Nodes:DegreesOfFreedom);
        }
    }
    return true;
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
        return true;
        default:
              stringstream temp;
              temp << "Error - Interpolation On Domain: Dudley does not know anything about function space type " << functionSpaceType_target;
              throw DudleyAdapterException(temp.str());
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
        return true;
        case(Nodes):
        case(DegreesOfFreedom):
        return false;
        default:
                stringstream temp;
                temp << "Error - Interpolation On Domain: Dudley does not know anything about function space type " << functionSpaceType_target;
                throw DudleyAdapterException(temp.str());
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
        return true;
        default:
                stringstream temp;
                temp << "Error - Interpolation On Domain: Dudley does not know anything about function space type " << functionSpaceType_target;
                throw DudleyAdapterException(temp.str());
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
        return true;
        case(Nodes):
        case(DegreesOfFreedom):
        return false;
        default:
                stringstream temp;
                temp << "Error - Interpolation On Domain: Dudley does not know anything about function space type " << functionSpaceType_target;
                throw DudleyAdapterException(temp.str());
        }
        break;
   default:
      stringstream temp;
      temp << "Error - Interpolation On Domain: Dudley does not know anything about function space type " << functionSpaceType_source;
      throw DudleyAdapterException(temp.str());
      break;
   }
   checkDudleyError();
   return false;
}

bool MeshAdapter::probeInterpolationAcross(int functionSpaceType_source,const AbstractDomain& targetDomain, int functionSpaceType_target) const
{
   return false;
}

bool MeshAdapter::operator==(const AbstractDomain& other) const
{
   const MeshAdapter* temp=dynamic_cast<const MeshAdapter*>(&other);
   if (temp!=0) {
      return (m_dudleyMesh==temp->m_dudleyMesh);
   } else {
      return false;
   }
}

bool MeshAdapter::operator!=(const AbstractDomain& other) const
{
   return !(operator==(other));
}

int MeshAdapter::getSystemMatrixTypeId(const boost::python::object& options) const
{
    const escript::SolverBuddy& sb = bp::extract<escript::SolverBuddy>(options);

    return paso::SystemMatrix::getSystemMatrixTypeId(sb.getSolverMethod(),
                sb.getPreconditioner(), sb.getPackage(), sb.isSymmetric(),
                m_dudleyMesh->MPIInfo);
}

int MeshAdapter::getTransportTypeId(int solver, int preconditioner,
                                    int package, bool symmetry) const
{
   Dudley_Mesh* mesh=m_dudleyMesh.get();
   return paso::TransportProblem::getTypeId(solver, preconditioner, package,
                                            symmetry, mesh->MPIInfo);
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
   Dudley_Mesh* mesh=m_dudleyMesh.get();
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
   case(DegreesOfFreedom):
   out=mesh->Nodes->degreesOfFreedomId;
   break;
   case(ReducedDegreesOfFreedom):
   out=mesh->Nodes->reducedDegreesOfFreedomId;
   break;
   default:
      stringstream temp;
      temp << "Error - Invalid function space type: " << functionSpaceType << " for domain: " << getDescription();
      throw DudleyAdapterException(temp.str());
      break;
   }
   return out;
}
int MeshAdapter::getTagFromSampleNo(int functionSpaceType, int sampleNo) const
{
   int out=0;
   Dudley_Mesh* mesh=m_dudleyMesh.get();
   switch (functionSpaceType) {
   case(Nodes):
   out=mesh->Nodes->Tag[sampleNo];
   break;
   case(ReducedNodes):
   throw DudleyAdapterException(" Error - ReducedNodes does not support tags.");
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
   case(DegreesOfFreedom):
   throw DudleyAdapterException(" Error - DegreesOfFreedom does not support tags.");
   break;
   case(ReducedDegreesOfFreedom):
   throw DudleyAdapterException(" Error - ReducedDegreesOfFreedom does not support tags.");
   break;
   default:
      stringstream temp;
      temp << "Error - Invalid function space type: " << functionSpaceType << " for domain: " << getDescription();
      throw DudleyAdapterException(temp.str());
      break;
   }
   return out;
}


void MeshAdapter::setTags(const int functionSpaceType, const int newTag, const escript::Data& mask) const
{
   Dudley_Mesh* mesh=m_dudleyMesh.get();
   switch(functionSpaceType) {
   case(Nodes):
   Dudley_NodeFile_setTags(mesh->Nodes,newTag,&mask);
   break;
   case(ReducedNodes):
   throw DudleyAdapterException("Error - ReducedNodes does not support tags");
   break;
   case(DegreesOfFreedom):
   throw DudleyAdapterException("Error - DegreesOfFreedom does not support tags");
   break;
   case(ReducedDegreesOfFreedom):
   throw DudleyAdapterException("Error - ReducedDegreesOfFreedom does not support tags");
   break;
   case(Elements):
   Dudley_ElementFile_setTags(mesh->Elements,newTag,&mask);
   break;
   case(ReducedElements):
   Dudley_ElementFile_setTags(mesh->Elements,newTag,&mask);
   break;
   case(FaceElements):
   Dudley_ElementFile_setTags(mesh->FaceElements,newTag,&mask);
   break;
   case(ReducedFaceElements):
   Dudley_ElementFile_setTags(mesh->FaceElements,newTag,&mask);
   break;
   case(Points):
   Dudley_ElementFile_setTags(mesh->Points,newTag,&mask);
   break;
   default:
      stringstream temp;
      temp << "Error - Dudley does not know anything about function space type " << functionSpaceType;
      throw DudleyAdapterException(temp.str());
   }
   checkDudleyError();
   return;
}

void MeshAdapter::setTagMap(const string& name,  int tag)
{
   Dudley_Mesh* mesh=m_dudleyMesh.get();
   Dudley_Mesh_addTagMap(mesh, name.c_str(),tag);
   checkDudleyError();
   // throwStandardException("MeshAdapter::set TagMap is not implemented.");
}

int MeshAdapter::getTag(const string& name) const
{
   Dudley_Mesh* mesh=m_dudleyMesh.get();
   int tag=0;
   tag=Dudley_Mesh_getTag(mesh, name.c_str());
   checkDudleyError();
   // throwStandardException("MeshAdapter::getTag is not implemented.");
   return tag;
}

bool MeshAdapter::isValidTagName(const string& name) const
{
   Dudley_Mesh* mesh=m_dudleyMesh.get();
   return Dudley_Mesh_isValidTagName(mesh,name.c_str());
}

string MeshAdapter::showTagNames() const
{
   stringstream temp;
   Dudley_Mesh* mesh=m_dudleyMesh.get();
   Dudley_TagMap* tag_map=mesh->TagMap;
   while (tag_map) {
      temp << tag_map->name;
      tag_map=tag_map->next;
      if (tag_map) temp << ", ";
   }
   return temp.str();
}

int MeshAdapter::getNumberOfTagsInUse(int functionSpaceCode) const
{
  Dudley_Mesh* mesh=m_dudleyMesh.get();
  dim_t numTags=0;
  switch(functionSpaceCode) {
   case(Nodes):
          numTags=mesh->Nodes->numTagsInUse;
          break;
   case(ReducedNodes):
          throw DudleyAdapterException("Error - ReducedNodes does not support tags");
          break;
   case(DegreesOfFreedom):
          throw DudleyAdapterException("Error - DegreesOfFreedom does not support tags");
          break;
   case(ReducedDegreesOfFreedom):
          throw DudleyAdapterException("Error - ReducedDegreesOfFreedom does not support tags");
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
   default:
      stringstream temp;
      temp << "Error - Dudley does not know anything about function space type " << functionSpaceCode;
      throw DudleyAdapterException(temp.str());
  }
  return numTags;
}

const int* MeshAdapter::borrowListOfTagsInUse(int functionSpaceCode) const
{
  Dudley_Mesh* mesh=m_dudleyMesh.get();
  index_t* tags=NULL;
  switch(functionSpaceCode) {
   case(Nodes):
          tags=mesh->Nodes->tagsInUse;
          break;
   case(ReducedNodes):
          throw DudleyAdapterException("Error - ReducedNodes does not support tags");
          break;
   case(DegreesOfFreedom):
          throw DudleyAdapterException("Error - DegreesOfFreedom does not support tags");
          break;
   case(ReducedDegreesOfFreedom):
          throw DudleyAdapterException("Error - ReducedDegreesOfFreedom does not support tags");
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
   default:
      stringstream temp;
      temp << "Error - Dudley does not know anything about function space type " << functionSpaceCode;
      throw DudleyAdapterException(temp.str());
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
  Dudley_Mesh* mesh=m_dudleyMesh.get();
  return Dudley_Mesh_getStatus(mesh);
}

int MeshAdapter::getApproximationOrder(const int functionSpaceCode) const
{
   
  Dudley_Mesh* mesh=m_dudleyMesh.get();
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
          order=mesh->integrationOrder;
          break;
   case(ReducedElements):
   case(ReducedFaceElements):
          order=mesh->reducedIntegrationOrder;
          break;
   default:
      stringstream temp;
      temp << "Error - Dudley does not know anything about function space type " << functionSpaceCode;
      throw DudleyAdapterException(temp.str());
  }
  return order;
}


bool MeshAdapter::supportsContactElements() const
{
    return false;
}

escript::Data MeshAdapter::randomFill(const escript::DataTypes::ShapeType& shape,
       const escript::FunctionSpace& what, long seed,
       const boost::python::tuple& filter) const
{
    Data towipe(0, shape, what, true);
    // since we just made this object, no sharing is possible and we don't need to check for
    // exlusive write
    escript::DataTypes::ValueType& dv=towipe.getExpandedVectorReference();
    const size_t dvsize=dv.size();
    esysUtils::randomFillArray(seed, &(dv[0]), dvsize);
    return towipe;       
}


}  // end of namespace
