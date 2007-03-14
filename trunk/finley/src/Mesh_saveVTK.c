/*
************************************************************
*          Copyright 2006 by ACcESS MNRF                   *
*                                                          *
*              http://www.access.edu.au                    *
*       Primary Business: Queensland, Australia            *
*  Licensed under the Open Software License version 3.0    *
*     http://www.opensource.org/licenses/osl-3.0.php       *
*                                                          *
************************************************************
*/


/**************************************************************/

/*   writes data and mesh in a vtk file */

/**************************************************************/

/*   Author: Paul Cochrane, cochrane@esscc.uq.edu.au */
/*   MPI-IO version: Derick Hawcroft, d.hawcroft@uq.edu.au         */

/*   Version: $Id$ */

/**************************************************************/


#include "Mesh.h"
#include "vtkCellType.h"  /* copied from vtk source directory !!! */

/*
 MPI version notes:
 
 ******************************************************************************
 ***                                                                       ****
 *** WARNING: Won't work for meshes with peridodic boundary conditions yet **** 
 ***                                                                       ****  
 ******************************************************************************
 
 In this version, the rank==0 process writes *all* opening and closing 
 XML tags.
 Individual process data is copied to a buffer before being written
 out. The  routines are collectively called and will be called in the natural
 ordering i.e 0 to maxProcs-1.
 
 Notable Notables:
 the struct localIndexCache: stores local domain indices for faster  reference
*/

#ifdef PASO_MPI


//#define MPIO_HINTS



#define MPIO_DEBUG(str) \
{ \
	if(myRank == 0) \
	printf("==== MPI-IO => %s \n", str); \
}

void Finley_Mesh_saveVTK_MPIO(const char * filename_p, Finley_Mesh *mesh_p, const dim_t num_data,char* *names_p, escriptDataC* *data_pp)
{
  int    numPoints,
  numCells = -1,
             myRank,comm,gsize,
             numLocal,
             nDim,
             shape;
  size_t __n;
  int i,j,k,m,n,count;
  int numGlobalCells = 0;
  index_t  *nodesGlobal=NULL;   // used to get the connectivity  right for VTK

  /* variables associatted with write_celldata/pointdata */
  int numPointsPerSample,
  nComp,
  nCompReqd;
  double* values, rtmp;

  // Local element info (for debugging)
  size_t numLocalCells,
  numInternalCells,
  numBoundaryCells;

  int rank;

  int amode = MPI_MODE_CREATE | MPI_MODE_WRONLY |  MPI_MODE_SEQUENTIAL;

  comm   = mesh_p->Nodes->MPIInfo->comm;
  myRank = mesh_p->Nodes->MPIInfo->rank;
  gsize  = mesh_p->Nodes->MPIInfo->size;

  MPI_File fh;
  MPI_Status status;
  MPI_Request req;
  MPI_Info infoHints;

  char error_msg[LenErrorMsg_MAX];

  int i_data;

  int nodetype=FINLEY_DEGREES_OF_FREEDOM;
  int elementtype=FINLEY_UNKNOWN;
  bool_t isCellCentered[num_data],write_celldata=FALSE,write_pointdata=FALSE;

  ElementTypeId TypeId;

  int numVTKNodesPerElement;
  int cellType;
  char elemTypeStr[32];

  Finley_ElementFile* elements=NULL;


  // Local node info
  int numInternalNodes,
  numLocalNodes,
  numBoundaryNodes,
  localDOF;  // equals to  (DOF of Internal Nodes) +  (DOF of Boundary Nodes) of local domain


  nDim  = mesh_p->Nodes->numDim;

#ifdef MPIO_HINTS
  MPI_Info_create(&infoHints);
  //  MPI_Info_set(infoHints, "striping_unit",        "424288");
  //  MPI_Info_set(infoHints, "striping_factor",      "16");
  //  MPI_Info_set(infoHints, "collective_buffering", "true");
  //  MPI_Info_set(infoHints, "cb_block_size",        "131072");
  //  MPI_Info_set(infoHints, "cb_buffer_size",       "1048567");
  //  MPI_Info_set(infoHints, "cb_nodes",             "8");
  //    MPI_Info_set(infoHints, "access_style", "write_once, sequential");

  //XFS only
  //   MPI_Info_set(infoHints, "direct_write",          "true");
#else
  infoHints = MPI_INFO_NULL;
#endif

  // Holds a local node/element index into the global array
  struct localIndexCache
  {
    index_t *values;
    int size;
  };

  struct localIndexCache nodeCache,
  		  elementCache;

  // Collective Call
  MPI_File_open(mesh_p->Nodes->MPIInfo->comm, (char*)filename_p, amode,infoHints, &fh);
  MPI_File_set_view(fh,MPI_DISPLACEMENT_CURRENT,MPI_CHAR, MPI_CHAR, "native" , infoHints);

  MPIO_DEBUG(" ***** Enter saveVTK ******")

  for (i_data=0;i_data<num_data;++i_data)
  {
    if (! isEmpty(data_pp[i_data]) )
    {
      switch(getFunctionSpaceType(data_pp[i_data]))
      {
      case FINLEY_DEGREES_OF_FREEDOM:
        nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
        if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS)
        {
          elementtype=FINLEY_ELEMENTS;
        }
        else
        {
          Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
          return ;
        }
        isCellCentered[i_data]=FALSE;
        break;
      case FINLEY_REDUCED_DEGREES_OF_FREEDOM:
        nodetype = FINLEY_REDUCED_DEGREES_OF_FREEDOM;
        if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS)
        {
          elementtype=FINLEY_ELEMENTS;
        }
        else
        {
          Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
          return;
        }
        isCellCentered[i_data]=FALSE;
        break;
      case FINLEY_NODES:
        nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
        if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS)
        {
          elementtype=FINLEY_ELEMENTS;
        }
        else
        {
          Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
          return;
        }
        isCellCentered[i_data]=FALSE;
        break;
      case FINLEY_ELEMENTS:
        nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
        if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS)
        {
          elementtype=FINLEY_ELEMENTS;
        }
        else
        {
          Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
          return;
        }
        isCellCentered[i_data]=TRUE;
        break;
      case FINLEY_FACE_ELEMENTS:
        nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
        if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_FACE_ELEMENTS)
        {
          elementtype=FINLEY_FACE_ELEMENTS;
        }
        else
        {
          Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
          return;

        }
        isCellCentered[i_data]=TRUE;
        break;
      case FINLEY_POINTS:
        nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
        if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_POINTS)
        {
          elementtype=FINLEY_POINTS;
        }
        else
        {
          Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
          return;

        }
        isCellCentered[i_data]=TRUE;
        break;
      case FINLEY_CONTACT_ELEMENTS_1:
        nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
        if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_CONTACT_ELEMENTS_1)
        {
          elementtype=FINLEY_CONTACT_ELEMENTS_1;
        }
        else
        {
          Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
          return;

        }
        isCellCentered[i_data]=TRUE;
        break;
      case FINLEY_CONTACT_ELEMENTS_2:
        nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
        if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_CONTACT_ELEMENTS_1)
        {
          elementtype=FINLEY_CONTACT_ELEMENTS_1;
        }
        else
        {
          Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
          return;

        }
        isCellCentered[i_data]=TRUE;
        break;
      default:
        sprintf(error_msg,"saveVTK: Finley does not know anything about function space type %d",getFunctionSpaceType(data_pp[i_data]));
        Finley_setError(TYPE_ERROR,error_msg);
        return;

      }

      if (isCellCentered[i_data])
      {
        write_celldata=TRUE;
      }
      else
      {
        write_pointdata=TRUE;
      }
    }
  }

  Finley_NodeDistribution *dist;
  if (nodetype==FINLEY_REDUCED_DEGREES_OF_FREEDOM)
  {
    dist = mesh_p->Nodes->reducedDegreeOfFreedomDistribution;
  }
  else
  {
    dist = mesh_p->Nodes->degreeOfFreedomDistribution;
  }

  numInternalNodes = dist->numInternal;
  numBoundaryNodes = dist->numBoundary;

  localDOF =  dist->numLocal;

  numPoints        = dist->numGlobal;

  if (elementtype==FINLEY_UNKNOWN) elementtype=FINLEY_ELEMENTS;
  switch(elementtype)
  {
  case FINLEY_ELEMENTS:
    elements=mesh_p->Elements;
    break;
  case FINLEY_FACE_ELEMENTS:
    elements=mesh_p->FaceElements;
    break;
  case FINLEY_POINTS:
    elements=mesh_p->Points;
    break;
  case FINLEY_CONTACT_ELEMENTS_1:
    elements=mesh_p->ContactElements;
    break;
  }
  if (elements==NULL)
  {
    Finley_setError(SYSTEM_ERROR,"saveVTK: undefined element file");
    return ;
  }

  numCells =  elements->numElements;
  numGlobalCells = elements->elementDistribution->vtxdist[gsize];
  numLocalCells    = elements->elementDistribution->numLocal;
  numInternalCells = elements->elementDistribution->numInternal;
  numBoundaryCells = elements->elementDistribution->numBoundary;

  if (nodetype==FINLEY_REDUCED_DEGREES_OF_FREEDOM)
  {
    TypeId = elements->LinearReferenceElement->Type->TypeId;
  }
  else
  {
    TypeId = elements->ReferenceElement->Type->TypeId;
  }

  switch(TypeId)
  {
  case Point1:
  case Line2Face:
  case Line3Face:
  case Point1_Contact:
  case Line2Face_Contact:
  case Line3Face_Contact:
    cellType = VTK_VERTEX;
    numVTKNodesPerElement = 1;
    strcpy(elemTypeStr, "VTK_VERTEX");
    break;

  case Line2:
  case Tri3Face:
  case Rec4Face:
  case Line2_Contact:
  case Tri3_Contact:
  case Tri3Face_Contact:
  case Rec4Face_Contact:
    cellType = VTK_LINE;
    numVTKNodesPerElement = 2;
    strcpy(elemTypeStr, "VTK_LINE");
    break;

  case Tri3:
  case Tet4Face:
  case Tet4Face_Contact:
    cellType = VTK_TRIANGLE;
    numVTKNodesPerElement = 3;
    strcpy(elemTypeStr, "VTK_TRIANGLE");
    break;

  case Rec4:
  case Hex8Face:
  case Rec4_Contact:
  case Hex8Face_Contact:
    cellType = VTK_QUAD;
    numVTKNodesPerElement = 4;
    strcpy(elemTypeStr, "VTK_QUAD");
    break;

  case Tet4:
    cellType = VTK_TETRA;
    numVTKNodesPerElement = 4;
    strcpy(elemTypeStr, "VTK_TETRA");
    break;

  case Hex8:
    cellType = VTK_HEXAHEDRON;
    numVTKNodesPerElement = 8;
    strcpy(elemTypeStr, "VTK_HEXAHEDRON");
    break;

  case Line3:
  case Tri6Face:
  case Rec8Face:
  case Line3_Contact:
  case Tri6Face_Contact:
  case Rec8Face_Contact:
    cellType = VTK_QUADRATIC_EDGE;
    numVTKNodesPerElement = 3;
    strcpy(elemTypeStr, "VTK_QUADRATIC_EDGE");
    break;

  case Tri6:
  case Tet10Face:
  case Tri6_Contact:
  case Tet10Face_Contact:
    cellType = VTK_QUADRATIC_TRIANGLE;
    numVTKNodesPerElement = 6;
    strcpy(elemTypeStr, "VTK_QUADRATIC_TRIANGLE");
    break;

  case Rec8:
  case Hex20Face:
  case Rec8_Contact:
  case Hex20Face_Contact:
    cellType = VTK_QUADRATIC_QUAD;
    numVTKNodesPerElement = 8;
    strcpy(elemTypeStr, "VTK_QUADRATIC_QUAD");
    break;

  case Tet10:
    cellType = VTK_QUADRATIC_TETRA;
    numVTKNodesPerElement = 10;
    strcpy(elemTypeStr, "VTK_QUADRATIC_TETRA");
    break;

  case Hex20:
    cellType = VTK_QUADRATIC_HEXAHEDRON;
    numVTKNodesPerElement = 20;
    strcpy(elemTypeStr, "VTK_QUADRATIC_HEXAHEDRON");
    break;

  default:
    sprintf(error_msg, "saveVTK: Element type %s is not supported by VTK",elements->ReferenceElement->Type->Name);
    Finley_setError(VALUE_ERROR,error_msg);
    return;
  }

  /* Write XML Header */
  if(myRank == 0)
  {
    char header[400];

    sprintf(header,"<?xml version=\"1.0\"?>\n" \
            "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n" \
            "<UnstructuredGrid>\n" \
            "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n" \
            "<Points>\n" \
            "<DataArray NumberOfComponents=\"%d\" type=\"Float64\" format=\"ascii\">\n"
            ,numPoints,numGlobalCells,MAX(3,nDim));


    MPI_File_iwrite_shared(fh,header,strlen(header),MPI_CHAR,&req);
    MPI_Wait(&req,&status);
  }

  MPIO_DEBUG(" Writing Coordinate Points... ")

  numLocalNodes=localDOF;

  //  we will be writing values which vary from 13-15 chars hence the strlen()
  char* largebuf = MEMALLOC( numLocalNodes*15*nDim + numLocalNodes*2 + 1 ,char);
  largebuf[0] = '\0';
  char tmpbuf[15];
  int tsz=0;
  int numNodesOutput=0;
  index_t pos=0;

  index_t *vtxdist = NULL, *DOFNodes=NULL,*forwardBuffer=NULL,*backwardBuffer=NULL;

  DOFNodes   = MEMALLOC(mesh_p->Nodes->numNodes,index_t);
  nodeCache.values = MEMALLOC( numLocalNodes, index_t);
  index_t bc_pos = 0;

  // Custom string concat:  avoids expensive strlen(3) call by strcat(3)
  // Note the implicit assumption on the variable "tsz"
  int __len,__j;
  char  *zero = "0.000000e+00 ";
  char  *newline = "\n";
  
#define __STRCAT(dest,chunk,tsz)  \
{                  \
   __len = strlen(chunk); \
   __j = -1;      \
   while(__j++ < __len)  \
    *(dest+tsz+__j)=*(chunk+__j); \
   tsz+=__len;              \
}
  
  // Loop over all nodes    
  for (i = 0; i < mesh_p->Nodes->numNodes; i++)
  {
    // This is the bit that will break for periodic BCs because it assumes that there is a one to one
    // correspondance between nodes and Degrees of freedom
    //TODO: handle periodic BC's 
    DOFNodes[mesh_p->Nodes->degreeOfFreedom[i]] = i;

    // Is this node local to the domain ?
    if( mesh_p->Nodes->degreeOfFreedom[i] < localDOF )
    {
      for (j = 0; j < nDim; j++)
      {
        sprintf(tmpbuf,"%e ", mesh_p->Nodes->Coordinates[INDEX2(j, i, nDim)] );
        __STRCAT(largebuf,tmpbuf,tsz)
      }
      for (k=0; k<3-nDim; k++)
      {
        __STRCAT(largebuf,zero,tsz)
      }
      __STRCAT(largebuf,newline,tsz)
      nodeCache.values[numNodesOutput++]=i;
    }
  }

  nodeCache.size=numNodesOutput;

  largebuf[tsz] = '\0';
  MPI_File_write_ordered(fh, largebuf,tsz, MPI_CHAR, &status);
  MEMFREE(largebuf);

  nodesGlobal = MEMALLOC(mesh_p->Nodes->numNodes ,index_t);

  // form distribution info on who output which nodes
  vtxdist = MEMALLOC( gsize+1, index_t );
  vtxdist[0]=0;
  MPI_Allgather(&numNodesOutput,1,MPI_INT,vtxdist+1,1,MPI_INT,comm);
  for( i=0; i<gsize; i++ )
    vtxdist[i+1]+=vtxdist[i];

  // will not work for periodic boundary conditions
  // calculate the local nodes file positions
  pos = 0;
  for( i=0; i<mesh_p->Nodes->numNodes; i++ )
  {
    if( mesh_p->Nodes->degreeOfFreedom[i]< localDOF )
    {
      nodesGlobal[i] = vtxdist[myRank] + pos++;
    }
    else
      nodesGlobal[i] = -1;
  }

  // communicate the local Nodes file position to the interested parties
  // send local info
  forwardBuffer = MEMALLOC( mesh_p->Nodes->numNodes, index_t );
  for( n=0; n < dist->numNeighbours; n++ )
  {
    if(  dist->edges[n]->numForward)
    {
      for( i=0; i < dist->edges[n]->numForward; i++ )
        forwardBuffer[i] = nodesGlobal[DOFNodes[dist->edges[n]->indexForward[i] ]];
      Paso_CommBuffer_pack( mesh_p->Nodes->CommBuffer, dist->neighbours[n], NULL, forwardBuffer, sizeof(index_t), 0 );
      Paso_CommBuffer_send( mesh_p->Nodes->CommBuffer, dist->neighbours[n], sizeof(index_t) );
    }
  }
  // receive external info
  backwardBuffer = MEMALLOC( mesh_p->Nodes->numNodes, index_t );
  for( n=0; n < dist->numNeighbours; n++ )
  {
    if( dist->edges[n]->numBackward )
    {
      Paso_CommBuffer_recv(mesh_p->Nodes->CommBuffer, dist->neighbours[n], sizeof(index_t));
      Paso_CommBuffer_unpack(mesh_p->Nodes->CommBuffer, dist->neighbours[n], NULL, backwardBuffer, sizeof(index_t), 0 );
      /* TODO: voodoo to handle perdiodic  BC's */
      for( i=0; i<dist->edges[n]->numBackward; i++ )
        nodesGlobal[DOFNodes[dist->edges[n]->indexBackward[i] ]] = backwardBuffer[i];
    }
  }
  

  
  MEMFREE(vtxdist);
  MEMFREE(DOFNodes);
  MEMFREE(backwardBuffer);
  MEMFREE(forwardBuffer);

  if( myRank == 0)
  {
    char* tags = "</DataArray>\n</Points>\n<Cells>\n<DataArray Name=\"connectivity\" type=\"Int32\" " \
                 "format=\"ascii\">\n" ;
    MPI_File_iwrite_shared(fh,tags,strlen(tags),MPI_CHAR,&req);
    MPI_Wait(&req,&status);
  }
  MPIO_DEBUG(" Done Writing Coordinate Points ")

  /* BEGIN CONNECTIVITY */

  int NN = elements->ReferenceElement->Type->numNodes; /* num Nodes holding ref-element */

  // Collective
  MPIO_DEBUG(" Writing Connectivity... ")

  size_t sz = numLocalCells*6*numVTKNodesPerElement + numLocalCells;
  largebuf = MEMALLOC(sz,char);
  largebuf[0] = '\0';
  tsz=0;
  pos = 0;
  // numCells?
  elementCache.values = MEMALLOC(numLocalCells,index_t);
  if (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM)
  {
    for (i = 0; i < numCells; i++)
    {

      if (elements->Id[i] >= elements->elementDistribution->vtxdist[i] &&  elements->Id[i] <= elements->elementDistribution->vtxdist[i+1] - 1 )
      {
        for (j = 0; j < numVTKNodesPerElement; j++)
        {
          sprintf(tmpbuf,"%d ",nodesGlobal[mesh_p->Nodes->toReduced[elements->Nodes[INDEX2(elements->ReferenceElement->Type->linearNodes[j], i, NN)]]]);
          __STRCAT(largebuf,tmpbuf,tsz)
        }
        __STRCAT(largebuf,newline,tsz)
        elementCache.values[pos++]=i;
      }
    }
  }
  else if (VTK_QUADRATIC_HEXAHEDRON==cellType)
  {
    char tmpbuf2[20*20*2];

    for (i = 0; i < numCells; i++)
    {
      if( elements->Id[i] >= elements->elementDistribution->vtxdist[myRank] && elements->Id[i] <= elements->elementDistribution->vtxdist[myRank+1]-1)
      {
        sprintf(tmpbuf2,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
                nodesGlobal[elements->Nodes[INDEX2(0, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(1, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(2, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(3, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(4, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(5, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(6, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(7, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(8, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(9, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(10, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(11, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(16, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(17, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(18, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(19, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(12, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(13, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(14, i, NN)]],
                nodesGlobal[elements->Nodes[INDEX2(15, i, NN)]]);
        __STRCAT(largebuf,tmpbuf2,tsz)
        elementCache.values[pos++]=i;
      }
    }
  }
  else if (numVTKNodesPerElement!=NN)
  {

    for (i = 0; i < numCells; i++)
    {
      if (elements->Id[i] >= elements->elementDistribution->vtxdist[i] &&  elements->Id[i] <= elements->elementDistribution->vtxdist[i+1] - 1 )
      {
        for (j = 0; j < numVTKNodesPerElement; j++)
        {
          sprintf(tmpbuf,"%d ", nodesGlobal[elements->Nodes[INDEX2(elements->ReferenceElement->Type->geoNodes[j], i, NN)]]);
          __STRCAT(largebuf,tmpbuf,tsz)
        }
        __STRCAT(largebuf,newline,tsz)
        elementCache.values[pos++]=i;
      }
    }
  }
  else
  {
    for(i = 0;i  < numCells ; i++)
    {
      if( elements->Id[i] >= elements->elementDistribution->vtxdist[myRank] && elements->Id[i] <= elements->elementDistribution->vtxdist[myRank+1]-1)
      {
        for (j = 0; j < numVTKNodesPerElement; j++)
        {
          sprintf(tmpbuf,"%d ", nodesGlobal[ elements->Nodes[INDEX2(j, i, NN) ] ] );
          __STRCAT(largebuf,tmpbuf,tsz)
        }
        __STRCAT(largebuf,newline,tsz)
        elementCache.values[pos++]=i;
      }
    }
  }

  elementCache.size = pos;

  largebuf[tsz] = '\0';
  MPI_File_write_ordered(fh,largebuf,tsz, MPI_CHAR, &status);
  MEMFREE(largebuf);
  MPIO_DEBUG(" Done Writing Connectivity ")
  MPIO_DEBUG(" Writing Offsets & Types... ")

  // Non-Collective
  if( myRank == 0)
  {
    // write out the DataArray element for the offsets
    char* tag1 = "</DataArray>\n<DataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\">\n";
    char* tag2 = "</DataArray>\n";
    char *tag3 =  "<DataArray Name=\"types\" type=\"UInt8\" format=\"ascii\">\n";
    char *tag4 = "</DataArray>\n</Cells>\n";

    int n = numVTKNodesPerElement;

    // allocate an upper bound on number of bytes needed  
    int sz=0;
    int lg = log10(numGlobalCells * n) + 1;
    sz += numGlobalCells*lg; 
    sz += numGlobalCells;   
    tsz = 0;

    char* largebuf = MEMALLOC(sz  + strlen(tag1) + strlen(tag2) + strlen(tag3) + strlen(tag4),char);
    largebuf[0] ='\0';
    char tmp[10];
    __STRCAT(largebuf,tag1,tsz)
    for (i=numVTKNodesPerElement; i<=numGlobalCells*numVTKNodesPerElement; i+=numVTKNodesPerElement)
    {
      sprintf(tmp,"%d\n", i);
      __STRCAT(largebuf,tmp,tsz)
    }
    __STRCAT(largebuf,tag2,tsz)
    largebuf[tsz] = '\0';
    MPI_File_iwrite_shared(fh,largebuf, tsz,MPI_CHAR,&req);
    MPI_Wait(&req,&status);

    // re-using buffer!
    largebuf[0] = '\0';
    tsz = 0;
    __STRCAT(largebuf,tag3,tsz)
    for (i=0; i<numGlobalCells; i++)
    {
      sprintf(tmp, "%d\n", cellType);
      __STRCAT(largebuf,tmp,tsz)
    }
    __STRCAT(largebuf,tag4,tsz)
    largebuf[tsz] = '\0';
    MPI_File_iwrite_shared(fh,largebuf,tsz,MPI_CHAR,&req);
    MPI_Wait(&req,&status);
    MEMFREE(largebuf);
  }

  MPIO_DEBUG(" Done Writing Offsets & Types ")

  // Write Cell data header Tags
  if(myRank == 0)
  {
    MPIO_DEBUG(" Writing Cell Data ...")
    if( write_celldata)
    {
      char tmpBuf[80];
      char header[600];
      // mark the active data arrays
      bool_t set_scalar=FALSE,set_vector=FALSE, set_tensor=FALSE;
      sprintf(tmpBuf, "<CellData");
      strcat(header,tmpBuf);
      for (i_data =0 ;i_data<num_data;++i_data)
      {
        if (! isEmpty(data_pp[i_data]) && isCellCentered[i_data])
        {
          // if the rank == 0:   --> scalar data
          // if the rank == 1:   --> vector data
          // if the rank == 2:   --> tensor data

          switch(getDataPointRank(data_pp[i_data]))
          {
          case 0:
            if (! set_scalar)
            {
              sprintf(tmpBuf," Scalars=\"%s\"",names_p[i_data]);
              strcat(header,tmpBuf);
              set_scalar=TRUE;
            }
            break;
          case 1:
            if (! set_vector)
            {
              sprintf(tmpBuf," Vectors=\"%s\"",names_p[i_data]);
	      strcat(header,tmpBuf);
              set_vector=TRUE;
            }
            break;
          case 2:
            if (! set_tensor)
            {
              sprintf(tmpBuf," Tensors=\"%s\"",names_p[i_data]);
	      strcat(header,tmpBuf);
              set_tensor=TRUE;
            }
            break;
          default:
            sprintf(error_msg, "saveVTK: data %s: Vtk can't handle objects with rank greater than 2.",names_p[i_data]);
            Finley_setError(VALUE_ERROR,error_msg);
            return;
          }
        }
      }
      strcat(header, ">\n");
      MPI_File_iwrite_shared(fh,header,strlen(header),MPI_CHAR,&req);
      MPI_Wait(&req,&status);
    }
  }

  // write actual data (collective)
  if(write_celldata)
  {
    for (i_data =0 ;i_data<num_data;++i_data)
    {
      if (! isEmpty(data_pp[i_data]) && isCellCentered[i_data])
      {
        numPointsPerSample = elements->ReferenceElement->numQuadNodes;
        rank = getDataPointRank(data_pp[i_data]);
        nComp = getDataPointSize(data_pp[i_data]);
        nCompReqd=1;   // the number of components required by vtk
        shape=0;
        if (rank == 0)
        {
          nCompReqd = 1;
        }
        else if (rank == 1)
        {
          shape=getDataPointShape(data_pp[i_data], 0);
          if  (shape>3)
          {
            Finley_setError(VALUE_ERROR, "saveVTK: rank 1 object must have less then 4 components");
            return;
          }
          nCompReqd = 3;
        }
        else
        {
          shape=getDataPointShape(data_pp[i_data], 0);
          if  (shape>3 || shape != getDataPointShape(data_pp[i_data], 1))
          {
            Finley_setError(VALUE_ERROR, "saveVTK: rank 2 object must have less then 4x4 components and must have a square shape");
            return;
          }
          nCompReqd = 9;
        }

        if( myRank == 0)
        {
          char header[250];
          sprintf(header,"<DataArray Name=\"%s\" type=\"Float64\" NumberOfComponents=\"%d\" format=\"ascii\">\n",names_p[i_data], nCompReqd);
          MPI_File_iwrite_shared(fh,header,strlen(header),MPI_CHAR,&req);
          MPI_Wait(&req,&status);
        }

        // Write the actual data 
        char tmpbuf[15];
        char* largebuf = MEMALLOC(nCompReqd*numLocalCells*15 + numLocalCells*nCompReqd + nCompReqd + 15,char);
        largebuf[0] = '\0';
        size_t tsz = 0;

        double sampleAvg[nComp];

        for (k=0; k<elementCache.size; k++)
        {
          i = elementCache.values[k];

          values = getSampleData(data_pp[i_data], i);
          // averaging over the number of points in the sample
          for (n=0; n<nComp; n++)
          {
            if (isExpanded(data_pp[i_data])) {
               rtmp = 0.;
               for (j=0; j<numPointsPerSample; j++) rtmp += values[INDEX2(n,j,nComp)];
               sampleAvg[n] = rtmp/numPointsPerSample;
            } else {
               sampleAvg[n] = values[n];
            }
          }
          // if the number of required components is more than the number
          // of actual components, pad with zeros

          // probably only need to get shape of first element
          // write the data different ways for scalar, vector and tensor
          if (nCompReqd == 1)
          {
            sprintf(tmpbuf, " %e", sampleAvg[0]);
            __STRCAT(largebuf,tmpbuf,tsz)
          }
          else if (nCompReqd == 3)
          {
            // write out the data
            for (m=0; m<shape; m++)
            {
              sprintf(tmpbuf, " %e", sampleAvg[m]);
              __STRCAT(largebuf,tmpbuf,tsz)
            }
            for (m=0; m<nCompReqd-shape; m++)
            {
              __STRCAT(largebuf,zero,tsz)
            }
          }
          else if (nCompReqd == 9)
          {
            // tensor data, so have a 3x3 matrix to output as a row
            // of 9 data points
            count = 0;
            for (m=0; m<shape; m++)
            {
              for (n=0; n<shape; n++)
              {
                sprintf(tmpbuf, " %e", sampleAvg[count]);
                __STRCAT(largebuf,tmpbuf,tsz)
                count++;
              }
              for (n=0; n<3-shape; n++)
              {
                __STRCAT(largebuf,zero,tsz)
              }
            }
            for (m=0; m<3-shape; m++)
              for (n=0; n<3; n++)
              {
                __STRCAT(largebuf,zero,tsz)
              }
          }
          __STRCAT(largebuf,newline,tsz)
        }
        largebuf[tsz] = '\0';
        MPI_File_write_ordered(fh,largebuf,tsz,MPI_CHAR,&status);
        MEMFREE(largebuf);
        if( myRank == 0)
        {
          char *tag = "</DataArray>\n";
          MPI_File_iwrite_shared(fh,tag,strlen(tag),MPI_CHAR,&req);
          MPI_Wait(&req,&status);
        }

      }
    }
    // closing celldata tag
    if(myRank == 0)
    {
      char* tag =  "</CellData>\n";
      MPI_File_iwrite_shared(fh,tag,strlen(tag),MPI_CHAR,&req);
      MPI_Wait(&req,&status);
    }

    MPIO_DEBUG(" Done Writing Cell Data ")
  }


  // Write Point Data Header Tags
  if( myRank == 0)
  {
    char header[600];
    char tmpBuf[50];

    if (write_pointdata)
    {
      MPIO_DEBUG(" Writing Pointdata... ")
      // mark the active data arrays
      bool_t set_scalar=FALSE,set_vector=FALSE, set_tensor=FALSE;
      sprintf(header, "<PointData");
      for (i_data =0 ;i_data<num_data;++i_data)
      {
        if (! isEmpty(data_pp[i_data]) && !isCellCentered[i_data])
        {
          // if the rank == 0:   --> scalar data
          // if the rank == 1:   --> vector data
          // if the rank == 2:   --> tensor data

          switch(getDataPointRank(data_pp[i_data]))
          {
          case 0:
            if (! set_scalar)
            {
              sprintf(tmpBuf," Scalars=\"%s\"",names_p[i_data]);
              strcat(header,tmpBuf);
              set_scalar=TRUE;
            }
            break;
          case 1:
            if (! set_vector)
            {
              sprintf(tmpBuf," Vectors=\"%s\"",names_p[i_data]);
              strcat(header,tmpBuf);
              set_vector=TRUE;
            }
            break;
          case 2:
            if (! set_tensor)
            {
              sprintf(tmpBuf," Tensors=\"%s\"",names_p[i_data]);
              strcat(header,tmpBuf);
              set_tensor=TRUE;
            }
            break;
          default:
            sprintf(error_msg, "saveVTK: data %s: Vtk can't handle objects with rank greater than 2.",names_p[i_data]);
            Finley_setError(VALUE_ERROR,error_msg);
            return;
          }
        }
      }
      strcat(header, ">\n");
      MPI_File_iwrite_shared(fh,header,strlen(header),MPI_CHAR,&req);
      MPI_Wait(&req,&status);
    }
  }

  // write actual data
  if(write_pointdata)
  {
    for (i_data =0 ;i_data<num_data;++i_data)
    {
      if (! isEmpty(data_pp[i_data]) && !isCellCentered[i_data])
      {
        numPointsPerSample = elements->ReferenceElement->numQuadNodes;
        rank = getDataPointRank(data_pp[i_data]);
        nComp = getDataPointSize(data_pp[i_data]);
        nCompReqd=1;   // the number of components required by vtk
        shape=0;
        if (rank == 0)
        {
          nCompReqd = 1;
        }
        else if (rank == 1)
        {
          shape=getDataPointShape(data_pp[i_data], 0);
          if  (shape>3)
          {
            Finley_setError(VALUE_ERROR, "saveVTK: rank 1 object must have less then 4 components");
            return;
          }
          nCompReqd = 3;
        }
        else
        {
          shape=getDataPointShape(data_pp[i_data], 0);
          if  (shape>3 || shape != getDataPointShape(data_pp[i_data], 1))
          {
            Finley_setError(VALUE_ERROR, "saveVTK: rank 2 object must have less then 4x4 components and must have a square shape");
            return;
          }
          nCompReqd = 9;
        }

        if( myRank == 0)
        {
          char header[250];
          sprintf(header, "<DataArray Name=\"%s\" type=\"Float64\" NumberOfComponents=\"%d\" format=\"ascii\">\n",names_p[i_data], nCompReqd);
          MPI_File_iwrite_shared(fh,header,strlen(header),MPI_CHAR,&req);
          MPI_Wait(&req,&status);
        }
        // write out the data
        // if the number of required components is more than the number
        // of actual components, pad with zeros

        char tmpbuf[15];
        char* largebuf = MEMALLOC(nCompReqd*numLocalNodes*15 + numLocal*nCompReqd + nCompReqd + 15,char);
        largebuf[0] = '\0';
        bool_t do_write=TRUE;
        size_t tsz = 0;

        for(k=0;k < nodeCache.size;k++)
        {
          i = nodeCache.values[k];

          if (nodetype==FINLEY_REDUCED_DEGREES_OF_FREEDOM)
          {
            if (mesh_p->Nodes->toReduced[i]>=0)
            {
              switch(getFunctionSpaceType(data_pp[i_data]))
              {
              case FINLEY_DEGREES_OF_FREEDOM:
                values = getSampleData(data_pp[i_data],mesh_p->Nodes->degreeOfFreedom[i]);
                break;
              case FINLEY_REDUCED_DEGREES_OF_FREEDOM:
                values = getSampleData(data_pp[i_data],mesh_p->Nodes->reducedDegreeOfFreedom[i]);
                break;
              case FINLEY_NODES:
                values = getSampleData(data_pp[i_data],i);
                break;
              }
              do_write=TRUE;
            }
            else
            {
              do_write=FALSE;
            }
          }
          else
          {
            do_write=TRUE;
            switch(getFunctionSpaceType(data_pp[i_data]))
            {
            case FINLEY_DEGREES_OF_FREEDOM:
              values = getSampleData(data_pp[i_data],mesh_p->Nodes->degreeOfFreedom[i]);
              break;
            case FINLEY_NODES:
              values = getSampleData(data_pp[i_data],i);
              break;
            }
          }
          // write the data different ways for scalar, vector and tensor
          if (do_write)
          {
            if (nCompReqd == 1)
            {
              sprintf(tmpbuf," %e", values[0]);
              __STRCAT(largebuf,tmpbuf,tsz)
            }
            else if (nCompReqd == 3)
            {
              for (m=0; m<shape; m++)
              {

                sprintf(tmpbuf," %e",values[m]);
                __STRCAT(largebuf,tmpbuf,tsz)
              }
              for (m=0; m<nCompReqd-shape; m++)
              {
                __STRCAT(largebuf,zero,tsz)
              }
            }
            else if (nCompReqd == 9)
            {
              // tensor data, so have a 3x3 matrix to output as a row
              //  of 9 data points
              count = 0;
              for (m=0; m<shape; m++)
              {
                for (n=0; n<shape; n++)
                {
                  sprintf(tmpbuf," %e", values[count]);
                  __STRCAT(largebuf,tmpbuf,tsz)
                  count++;
                }
                for (n=0; n<3-shape; n++)
                {
                  __STRCAT(largebuf,zero,tsz)
                }
              }
              for (m=0; m<3-shape; m++)
              {
                for (n=0; n<3; n++)
                {
                  __STRCAT(largebuf,zero,tsz)
                }
              }
            }
            __STRCAT(largebuf,newline,tsz)
          }

        }
        // Write out local data

        largebuf[tsz] = '\0';
        MPI_File_write_ordered(fh,largebuf,tsz,MPI_CHAR,&status);
        MEMFREE(largebuf);
        if( myRank == 0)
        {
          char *tag = "</DataArray>\n";
          MPI_File_iwrite_shared(fh,tag,strlen(tag),MPI_CHAR,&req);
          MPI_Wait(&req,&status);
        }
      }
    }
    // Finish off with closing tag
    if(myRank == 0)
    {
      char* tag =  "</PointData>\n";
      MPI_File_iwrite_shared(fh,tag,strlen(tag),MPI_CHAR,&req);
      MPI_Wait(&req,&status);
    }
    MPIO_DEBUG(" Done Writing Pointdata ")
  }
  // end write_pointdata

  // tag and bag...  
  if (myRank == 0)
  {
    char *footer = "</Piece>\n</UnstructuredGrid>\n</VTKFile>";
    MPI_File_iwrite_shared(fh,footer,strlen(footer),MPI_CHAR,&req);
    MPI_Wait(&req,&status);
  }

  MEMFREE(nodesGlobal);
  MEMFREE(nodeCache.values);
  MEMFREE(elementCache.values);
#ifdef MPIO_HINTS
  MPI_Info_free(&infoHints);
#undef MPIO_HINTS
#endif
  MPI_File_close(&fh);
  MPIO_DEBUG(" ***** Exit saveVTK ***** ")
#undef __STRCAT
}

#undef MPIO_DEBUG
#else




void Finley_Mesh_saveVTK(const char * filename_p, Finley_Mesh *mesh_p, const dim_t num_data,char* *names_p, escriptDataC* *data_pp)
{
  #define NCOMP_MAX 9
  char error_msg[LenErrorMsg_MAX];
  double sampleAvg[NCOMP_MAX];
  /* if there is no mesh we just return */

  int i, j, k, numVTKNodesPerElement,i_data,m, count, n, rank,shape, numPoints, cellType, numCells,
  nDim, numPointsPerSample, nComp, nCompReqd, NN;
  index_t j2;
  int nodetype=FINLEY_DEGREES_OF_FREEDOM;
  int elementtype=FINLEY_UNKNOWN;
  double* values, rtmp;
  char elemTypeStr[32];
  FILE * fileHandle_p = NULL;
  bool_t do_write, *isCellCentered=NULL,write_celldata=FALSE,write_pointdata=FALSE;
  Finley_ElementFile* elements=NULL;
  ElementTypeId TypeId;

  /* open the file and check handle */
  if (mesh_p==NULL) return;

  fileHandle_p = fopen(filename_p, "w");
  if (fileHandle_p==NULL)
  {
    sprintf(error_msg, "saveVTK: File %s could not be opened for writing.", filename_p);
    Finley_setError(IO_ERROR,error_msg);
    return;
  }
  /* find the mesh type to be written */
  isCellCentered=TMPMEMALLOC(num_data,bool_t);


  if (Finley_checkPtr(isCellCentered)) {
     fclose(fileHandle_p);
     return;
  }
  for (i_data=0;i_data<num_data;++i_data)
  {
    if (! isEmpty(data_pp[i_data]))
    {
      switch(getFunctionSpaceType(data_pp[i_data]))
      {
      case FINLEY_DEGREES_OF_FREEDOM:
        nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
        if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS)
        {
          elementtype=FINLEY_ELEMENTS;
        }
        else
        {
          Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
          fclose(fileHandle_p);
          TMPMEMFREE(isCellCentered);
          return;
        }
        isCellCentered[i_data]=FALSE;
        break;
      case FINLEY_REDUCED_DEGREES_OF_FREEDOM:
        nodetype = FINLEY_REDUCED_DEGREES_OF_FREEDOM;
        if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS)
        {
          elementtype=FINLEY_ELEMENTS;
        }
        else
        {
          Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
          fclose(fileHandle_p);
          TMPMEMFREE(isCellCentered);
          return;
        }
        isCellCentered[i_data]=FALSE;
        break;
      case FINLEY_NODES:
        nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
        if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS)
        {
          elementtype=FINLEY_ELEMENTS;
        }
        else
        {
          Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
          fclose(fileHandle_p);
          TMPMEMFREE(isCellCentered);
          return;
        }
        isCellCentered[i_data]=FALSE;
        break;
      case FINLEY_ELEMENTS:
        nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
        if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS)
        {
          elementtype=FINLEY_ELEMENTS;
        }
        else
        {
          Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
          fclose(fileHandle_p);
          return;
          TMPMEMFREE(isCellCentered);
        }
        isCellCentered[i_data]=TRUE;
        break;
      case FINLEY_FACE_ELEMENTS:
        nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
        if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_FACE_ELEMENTS)
        {
          elementtype=FINLEY_FACE_ELEMENTS;
        }
        else
        {
          Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
          fclose(fileHandle_p);
          TMPMEMFREE(isCellCentered);
          return;
        }
        isCellCentered[i_data]=TRUE;
        break;
      case FINLEY_POINTS:
        nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
        if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_POINTS)
        {
          elementtype=FINLEY_POINTS;
        }
        else
        {
          Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
          fclose(fileHandle_p);
          TMPMEMFREE(isCellCentered);
          return;
        }
        isCellCentered[i_data]=TRUE;
        break;
      case FINLEY_CONTACT_ELEMENTS_1:
        nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
        if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_CONTACT_ELEMENTS_1)
        {
          elementtype=FINLEY_CONTACT_ELEMENTS_1;
        }
        else
        {
          Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
          fclose(fileHandle_p);
          TMPMEMFREE(isCellCentered);
          return;
        }
        isCellCentered[i_data]=TRUE;
        break;
      case FINLEY_CONTACT_ELEMENTS_2:
        nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
        if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_CONTACT_ELEMENTS_1)
        {
          elementtype=FINLEY_CONTACT_ELEMENTS_1;
        }
        else
        {
          Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
          fclose(fileHandle_p);
          TMPMEMFREE(isCellCentered);
          return;
        }
        isCellCentered[i_data]=TRUE;
        break;
      default:
        sprintf(error_msg,"saveVTK: Finley does not know anything about function space type %d",getFunctionSpaceType(data_pp[i_data]));
        Finley_setError(TYPE_ERROR,error_msg);
        fclose(fileHandle_p);
        TMPMEMFREE(isCellCentered);
        return;
      }
      if (isCellCentered[i_data])
      {
        write_celldata=TRUE;
      }
      else
      {
        write_pointdata=TRUE;
      }
    }
  }
  /* select nomber of points and the mesh component */
  numPoints = mesh_p->Nodes->numNodes;
  if (nodetype==FINLEY_REDUCED_DEGREES_OF_FREEDOM)
  {
    numPoints = mesh_p->Nodes->reducedNumNodes;
  }
  else
  {
    numPoints = mesh_p->Nodes->numNodes;
  }
  if (elementtype==FINLEY_UNKNOWN) elementtype=FINLEY_ELEMENTS;
  switch(elementtype)
  {
  case FINLEY_ELEMENTS:
    elements=mesh_p->Elements;
    break;
  case FINLEY_FACE_ELEMENTS:
    elements=mesh_p->FaceElements;
    break;
  case FINLEY_POINTS:
    elements=mesh_p->Points;
    break;
  case FINLEY_CONTACT_ELEMENTS_1:
    elements=mesh_p->ContactElements;
    break;
  }
  if (elements==NULL)
  {
    Finley_setError(SYSTEM_ERROR,"saveVTK: undefined element file");
    fclose(fileHandle_p);
    TMPMEMFREE(isCellCentered);
    return;
  }
  /* map finley element type to VTK element type */
  numCells = elements->numElements;
  if (nodetype==FINLEY_REDUCED_DEGREES_OF_FREEDOM)
  {
    TypeId = elements->LinearReferenceElement->Type->TypeId;
  }
  else
  {
    TypeId = elements->ReferenceElement->Type->TypeId;
  }

  switch(TypeId)
  {
  case Point1:
  case Line2Face:
  case Line3Face:
  case Point1_Contact:
  case Line2Face_Contact:
  case Line3Face_Contact:
    cellType = VTK_VERTEX;
    numVTKNodesPerElement = 1;
    strcpy(elemTypeStr, "VTK_VERTEX");
    break;

  case Line2:
  case Tri3Face:
  case Rec4Face:
  case Line2_Contact:
  case Tri3_Contact:
  case Tri3Face_Contact:
  case Rec4Face_Contact:
    cellType = VTK_LINE;
    numVTKNodesPerElement = 2;
    strcpy(elemTypeStr, "VTK_LINE");
    break;

  case Tri3:
  case Tet4Face:
  case Tet4Face_Contact:
    cellType = VTK_TRIANGLE;
    numVTKNodesPerElement = 3;
    strcpy(elemTypeStr, "VTK_TRIANGLE");
    break;

  case Rec4:
  case Hex8Face:
  case Rec4_Contact:
  case Hex8Face_Contact:
    cellType = VTK_QUAD;
    numVTKNodesPerElement = 4;
    strcpy(elemTypeStr, "VTK_QUAD");
    break;

  case Tet4:
    cellType = VTK_TETRA;
    numVTKNodesPerElement = 4;
    strcpy(elemTypeStr, "VTK_TETRA");
    break;

  case Hex8:
    cellType = VTK_HEXAHEDRON;
    numVTKNodesPerElement = 8;
    strcpy(elemTypeStr, "VTK_HEXAHEDRON");
    break;

  case Line3:
  case Tri6Face:
  case Rec8Face:
  case Line3_Contact:
  case Tri6Face_Contact:
  case Rec8Face_Contact:
    cellType = VTK_QUADRATIC_EDGE;
    numVTKNodesPerElement = 3;
    strcpy(elemTypeStr, "VTK_QUADRATIC_EDGE");
    break;

  case Tri6:
  case Tet10Face:
  case Tri6_Contact:
  case Tet10Face_Contact:
    cellType = VTK_QUADRATIC_TRIANGLE;
    numVTKNodesPerElement = 6;
    strcpy(elemTypeStr, "VTK_QUADRATIC_TRIANGLE");
    break;

  case Rec8:
  case Hex20Face:
  case Rec8_Contact:
  case Hex20Face_Contact:
    cellType = VTK_QUADRATIC_QUAD;
    numVTKNodesPerElement = 8;
    strcpy(elemTypeStr, "VTK_QUADRATIC_QUAD");
    break;

  case Tet10:
    cellType = VTK_QUADRATIC_TETRA;
    numVTKNodesPerElement = 10;
    strcpy(elemTypeStr, "VTK_QUADRATIC_TETRA");
    break;

  case Hex20:
    cellType = VTK_QUADRATIC_HEXAHEDRON;
    numVTKNodesPerElement = 20;
    strcpy(elemTypeStr, "VTK_QUADRATIC_HEXAHEDRON");
    break;

  default:
    sprintf(error_msg, "saveVTK: Element type %s is not supported by VTK",elements->ReferenceElement->Type->Name);
    Finley_setError(VALUE_ERROR,error_msg);
    fclose(fileHandle_p);
    TMPMEMFREE(isCellCentered);
    return;
  }
  /* xml header */
  fprintf(fileHandle_p, "<?xml version=\"1.0\"?>\n");
  fprintf(fileHandle_p, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n");

  /* finley uses an unstructured mesh, so UnstructuredGrid *should* work */
  fprintf(fileHandle_p, "<UnstructuredGrid>\n");

  /* is there only one "piece" to the data?? */
  fprintf(fileHandle_p, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",numPoints, numCells);
  /* now for the points; equivalent to positions section in saveDX() */
  /* "The points element explicitly defines coordinates for each point
  * individually.  It contains one DataArray element describing an array
  * with three components per value, each specifying the coordinates of one
  * point" - from Vtk User's Guide
  */
  fprintf(fileHandle_p, "<Points>\n");
  /*
  * the reason for this if statement is explained in the long comment below
  */
  nDim = mesh_p->Nodes->numDim;
  fprintf(fileHandle_p, "<DataArray NumberOfComponents=\"%d\" type=\"Float64\" format=\"ascii\">\n",MAX(3,nDim));
  /* vtk/mayavi doesn't like 2D data, it likes 3D data with a degenerate
  * third dimension to handle 2D data (like a sheet of paper).  So, if
  * nDim is 2, we have to append zeros to the array to get this third
  * dimension, and keep the visualisers happy.
  * Indeed, if nDim is less than 3, must pad all empty dimensions, so
  * that the total number of dims is 3.
  */
  if (nodetype==FINLEY_REDUCED_DEGREES_OF_FREEDOM)
  {
    for (i = 0; i < mesh_p->Nodes->numNodes; i++)
    {
      if (mesh_p->Nodes->toReduced[i]>=0)
      {
        for (j = 0; j < nDim; j++) fprintf(fileHandle_p, " %e",mesh_p->Nodes->Coordinates[INDEX2(j, i, nDim)]);
        for (k=0; k<3-nDim; k++) fprintf(fileHandle_p, " %e",0.);
        fprintf(fileHandle_p, "\n");
      }
    }
  }
  else
  {
    for (i = 0; i < mesh_p->Nodes->numNodes; i++)
    {

      for (j = 0; j < nDim; j++) fprintf(fileHandle_p, " %e",mesh_p->Nodes->Coordinates[INDEX2(j, i, nDim)]);
      for (k=0; k<3-nDim; k++) fprintf(fileHandle_p, " %e",0.);
      fprintf(fileHandle_p, "\n");
    }

  }
  fprintf(fileHandle_p, "</DataArray>\n");
  fprintf(fileHandle_p, "</Points>\n");

  /* write out the DataArray element for the connectivity */

  NN = elements->ReferenceElement->Type->numNodes;
  fprintf(fileHandle_p, "<Cells>\n");
  fprintf(fileHandle_p, "<DataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\">\n");

  if (nodetype==FINLEY_REDUCED_DEGREES_OF_FREEDOM)
  {
    for (i = 0; i < numCells; i++)
    {
      for (j = 0; j < numVTKNodesPerElement; j++)
        fprintf(fileHandle_p,"%d ",mesh_p->Nodes->toReduced[elements->Nodes[INDEX2(elements->ReferenceElement->Type->linearNodes[j], i, NN)]]);
      fprintf(fileHandle_p, "\n");
    }
  }
  else if (VTK_QUADRATIC_HEXAHEDRON==cellType)
  {
    for (i = 0; i < numCells; i++)
    {
      fprintf(fileHandle_p,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
              elements->Nodes[INDEX2(0, i, NN)],
              elements->Nodes[INDEX2(1, i, NN)],
              elements->Nodes[INDEX2(2, i, NN)],
              elements->Nodes[INDEX2(3, i, NN)],
              elements->Nodes[INDEX2(4, i, NN)],
              elements->Nodes[INDEX2(5, i, NN)],
              elements->Nodes[INDEX2(6, i, NN)],
              elements->Nodes[INDEX2(7, i, NN)],
              elements->Nodes[INDEX2(8, i, NN)],
              elements->Nodes[INDEX2(9, i, NN)],
              elements->Nodes[INDEX2(10, i, NN)],
              elements->Nodes[INDEX2(11, i, NN)],
              elements->Nodes[INDEX2(16, i, NN)],
              elements->Nodes[INDEX2(17, i, NN)],
              elements->Nodes[INDEX2(18, i, NN)],
              elements->Nodes[INDEX2(19, i, NN)],
              elements->Nodes[INDEX2(12, i, NN)],
              elements->Nodes[INDEX2(13, i, NN)],
              elements->Nodes[INDEX2(14, i, NN)],
              elements->Nodes[INDEX2(15, i, NN)]);
    }
  }
  else if (numVTKNodesPerElement!=NN)
  {
    for (i = 0; i < numCells; i++)
    {
      for (j = 0; j < numVTKNodesPerElement; j++) fprintf(fileHandle_p,"%d ", elements->Nodes[INDEX2(elements->ReferenceElement->Type->geoNodes[j], i, NN)]);
      fprintf(fileHandle_p, "\n");
    }
  }
  else
  {
    for (i = 0; i < numCells; i++)
    {
      for (j = 0; j < numVTKNodesPerElement; j++) fprintf(fileHandle_p,"%d ", elements->Nodes[INDEX2(j, i, NN)]);
      fprintf(fileHandle_p, "\n");
    }
  }
  fprintf(fileHandle_p, "</DataArray>\n");

  /* write out the DataArray element for the offsets */
  fprintf(fileHandle_p, "<DataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\">\n");
  for (i=numVTKNodesPerElement; i<=numCells*numVTKNodesPerElement; i+=numVTKNodesPerElement) fprintf(fileHandle_p, "%d\n", i);
  fprintf(fileHandle_p, "</DataArray>\n");

  /* write out the DataArray element for the types */
  fprintf(fileHandle_p, "<DataArray Name=\"types\" type=\"UInt8\" format=\"ascii\">\n");
  for (i=0; i<numCells; i++) fprintf(fileHandle_p, "%d\n", cellType);
  fprintf(fileHandle_p, "</DataArray>\n");

  /* finish off the <Cells> element */
  fprintf(fileHandle_p, "</Cells>\n");

  /* cell data */
  if (write_celldata)
  {
    /* mark the active data arrays */
    bool_t set_scalar=FALSE,set_vector=FALSE, set_tensor=FALSE;
    fprintf(fileHandle_p, "<CellData");
    for (i_data =0 ;i_data<num_data;++i_data)
    {
      if (! isEmpty(data_pp[i_data]) && isCellCentered[i_data])
      {
        /* if the rank == 0:   --> scalar data
        * if the rank == 1:   --> vector data
        * if the rank == 2:   --> tensor data
        */
        switch(getDataPointRank(data_pp[i_data]))
        {
        case 0:
          if (! set_scalar)
          {
            fprintf(fileHandle_p," Scalars=\"%s\"",names_p[i_data]);
            set_scalar=TRUE;
          }
          break;
        case 1:
          if (! set_vector)
          {
            fprintf(fileHandle_p," Vectors=\"%s\"",names_p[i_data]);
            set_vector=TRUE;
          }
          break;
        case 2:
          if (! set_tensor)
          {
            fprintf(fileHandle_p," Tensors=\"%s\"",names_p[i_data]);
            set_tensor=TRUE;
          }
          break;
        default:
          sprintf(error_msg, "saveVTK: data %s: Vtk can't handle objects with rank greater than 2.",names_p[i_data]);
          Finley_setError(VALUE_ERROR,error_msg);
          fclose(fileHandle_p);
          TMPMEMFREE(isCellCentered);
          return;
        }
      }
    }
    fprintf(fileHandle_p, ">\n");
    /* write the arrays */
    for (i_data =0 ;i_data<num_data;++i_data)
    {
      if (! isEmpty(data_pp[i_data]) && isCellCentered[i_data])
      {
        numPointsPerSample = elements->ReferenceElement->numQuadNodes;
        rank = getDataPointRank(data_pp[i_data]);
        nComp = getDataPointSize(data_pp[i_data]);
        nCompReqd=1;   /* the number of components required by vtk */
        shape=0;
        if (rank == 0)
        {
          nCompReqd = 1;
        }
        else if (rank == 1)
        {
          shape=getDataPointShape(data_pp[i_data], 0);
          if  (shape>3)
          {
            Finley_setError(VALUE_ERROR, "saveVTK: rank 1 object must have less then 4 components");
            fclose(fileHandle_p);
            TMPMEMFREE(isCellCentered);
            return;
          }
          nCompReqd = 3;
        }
        else
        {
          shape=getDataPointShape(data_pp[i_data], 0);
          if  (shape>3 || shape != getDataPointShape(data_pp[i_data], 1))
          {
            Finley_setError(VALUE_ERROR, "saveVTK: rank 2 object must have less then 4x4 components and must have a square shape");
            fclose(fileHandle_p);
            TMPMEMFREE(isCellCentered);
            return;
          }
          nCompReqd = 9;
        }
        fprintf(fileHandle_p, "<DataArray Name=\"%s\" type=\"Float64\" NumberOfComponents=\"%d\" format=\"ascii\">\n",names_p[i_data], nCompReqd);

        for (i=0; i<numCells; i++)
        {
          values = getSampleData(data_pp[i_data], i);
          /* averaging over the number of points in the sample */
          for (k=0; k<MIN(nComp,NCOMP_MAX); k++)
          {
            if (isExpanded(data_pp[i_data])) {
               rtmp = 0.;
               for (j=0; j<numPointsPerSample; j++) rtmp += values[INDEX2(k,j,nComp)];
               sampleAvg[k] = rtmp/numPointsPerSample;
            } else {
               sampleAvg[k] = values[k];
            }
 
          }
          /* if the number of required components is more than the number
          * of actual components, pad with zeros
          */
          /* probably only need to get shape of first element */
          /* write the data different ways for scalar, vector and tensor */
          if (nCompReqd == 1)
          {
            fprintf(fileHandle_p, " %e", sampleAvg[0]);
          }
          else if (nCompReqd == 3)
          {
            /* write out the data */
            for (m=0; m<shape; m++) fprintf(fileHandle_p, " %e", sampleAvg[m]);
            for (m=0; m<nCompReqd-shape; m++) fprintf(fileHandle_p, " %e", 0.);
          }
          else if (nCompReqd == 9)
          {
            /* tensor data, so have a 3x3 matrix to output as a row
            * of 9 data points */
            count = 0;
            for (m=0; m<shape; m++)
            {
              for (n=0; n<shape; n++)
              {
                fprintf(fileHandle_p, " %e", sampleAvg[count]);
                count++;
              }
              for (n=0; n<3-shape; n++) fprintf(fileHandle_p, " %e", 0.);
            }
            for (m=0; m<3-shape; m++)
              for (n=0; n<3; n++) fprintf(fileHandle_p, " %e", 0.);
          }
          fprintf(fileHandle_p, "\n");
        }
        fprintf(fileHandle_p, "</DataArray>\n");
      }
    }
    fprintf(fileHandle_p, "</CellData>\n");
  }
  /* point data */
  if (write_pointdata)
  {
    /* mark the active data arrays */
    bool_t set_scalar=FALSE,set_vector=FALSE, set_tensor=FALSE;
    fprintf(fileHandle_p, "<PointData");
    for (i_data =0 ;i_data<num_data;++i_data)
    {
      if (! isEmpty(data_pp[i_data]) && !isCellCentered[i_data])
      {
        /* if the rank == 0:   --> scalar data
        * if the rank == 1:   --> vector data
        * if the rank == 2:   --> tensor data
        */
        switch(getDataPointRank(data_pp[i_data]))
        {
        case 0:
          if (! set_scalar)
          {
            fprintf(fileHandle_p," Scalars=\"%s\"",names_p[i_data]);
            set_scalar=TRUE;
          }
          break;
        case 1:
          if (! set_vector)
          {
            fprintf(fileHandle_p," Vectors=\"%s\"",names_p[i_data]);
            set_vector=TRUE;
          }
          break;
        case 2:
          if (! set_tensor)
          {
            fprintf(fileHandle_p," Tensors=\"%s\"",names_p[i_data]);
            set_tensor=TRUE;
          }
          break;
        default:
          sprintf(error_msg, "saveVTK: data %s: Vtk can't handle objects with rank greater than 2.",names_p[i_data]);
          Finley_setError(VALUE_ERROR,error_msg);
          fclose(fileHandle_p);
          TMPMEMFREE(isCellCentered);
          return;
        }
      }
    }
    fprintf(fileHandle_p, ">\n");
    /* write the arrays */
    for (i_data =0 ;i_data<num_data;++i_data)
    {
      if (! isEmpty(data_pp[i_data]) && !isCellCentered[i_data])
      {
        numPointsPerSample = elements->ReferenceElement->numQuadNodes;
        rank = getDataPointRank(data_pp[i_data]);
        nComp = getDataPointSize(data_pp[i_data]);
        nCompReqd=1;   /* the number of components required by vtk */
        shape=0;
        if (rank == 0)
        {
          nCompReqd = 1;
        }
        else if (rank == 1)
        {
          shape=getDataPointShape(data_pp[i_data], 0);
          if  (shape>3)
          {
            Finley_setError(VALUE_ERROR, "saveVTK: rank 1 object must have less then 4 components");
            fclose(fileHandle_p);
            TMPMEMFREE(isCellCentered);
            return;
          }
          nCompReqd = 3;
        }
        else
        {
          shape=getDataPointShape(data_pp[i_data], 0);
          if  (shape>3 || shape != getDataPointShape(data_pp[i_data], 1))
          {
            Finley_setError(VALUE_ERROR, "saveVTK: rank 2 object must have less then 4x4 components and must have a square shape");
            fclose(fileHandle_p);
            TMPMEMFREE(isCellCentered);
            return;
          }
          nCompReqd = 9;
        }
        fprintf(fileHandle_p, "<DataArray Name=\"%s\" type=\"Float64\" NumberOfComponents=\"%d\" format=\"ascii\">\n",names_p[i_data], nCompReqd);
        /* write out the data */
        /* if the number of required components is more than the number
        * of actual components, pad with zeros
        */
        do_write=TRUE;
        for (i=0; i<mesh_p->Nodes->numNodes; i++)
        {
          if (nodetype==FINLEY_REDUCED_DEGREES_OF_FREEDOM)
          {
            if (mesh_p->Nodes->toReduced[i]>=0)
            {
              switch(getFunctionSpaceType(data_pp[i_data]))
              {
              case FINLEY_DEGREES_OF_FREEDOM:
                values = getSampleData(data_pp[i_data],mesh_p->Nodes->degreeOfFreedom[i]);
                break;
              case FINLEY_REDUCED_DEGREES_OF_FREEDOM:
                values = getSampleData(data_pp[i_data],mesh_p->Nodes->reducedDegreeOfFreedom[i]);
                break;
              case FINLEY_NODES:
                values = getSampleData(data_pp[i_data],i);
                break;
              }
              do_write=TRUE;
            }
            else
            {
              do_write=FALSE;
            }
          }
          else
          {
            do_write=TRUE;
            switch(getFunctionSpaceType(data_pp[i_data]))
            {
            case FINLEY_DEGREES_OF_FREEDOM:
              values = getSampleData(data_pp[i_data],mesh_p->Nodes->degreeOfFreedom[i]);
              break;
            case FINLEY_NODES:
              values = getSampleData(data_pp[i_data],i);
              break;
            }
          }
          /* write the data different ways for scalar, vector and tensor */
          if (do_write)
          {
            if (nCompReqd == 1)
            {
              fprintf(fileHandle_p, " %e", values[0]);
            }
            else if (nCompReqd == 3)
            {
              for (m=0; m<shape; m++) fprintf(fileHandle_p, " %e", values[m]);
              for (m=0; m<nCompReqd-shape; m++) fprintf(fileHandle_p, " %e", 0.);
            }
            else if (nCompReqd == 9)
            {
              /* tensor data, so have a 3x3 matrix to output as a row
              * of 9 data points */
              count = 0;
              for (m=0; m<shape; m++)
              {
                for (n=0; n<shape; n++)
                {
                  fprintf(fileHandle_p, " %e", values[count]);
                  count++;
                }
                for (n=0; n<3-shape; n++) fprintf(fileHandle_p, " %e", 0.);
              }
              for (m=0; m<3-shape; m++)
                for (n=0; n<3; n++) fprintf(fileHandle_p, " %e", 0.);
            }
            fprintf(fileHandle_p, "\n");
          }
        }
        fprintf(fileHandle_p, "</DataArray>\n");
      }
    }
    fprintf(fileHandle_p, "</PointData>\n");
  }
  /* finish off the piece */
  fprintf(fileHandle_p, "</Piece>\n");

  fprintf(fileHandle_p, "</UnstructuredGrid>\n");
  /* write the xml footer */
  fprintf(fileHandle_p, "</VTKFile>\n");
  /* close the file */
  fclose(fileHandle_p);
  TMPMEMFREE(isCellCentered);
  return;
}
#endif

