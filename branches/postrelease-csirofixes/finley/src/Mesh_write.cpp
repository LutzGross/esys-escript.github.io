
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


/************************************************************************************/

/*   Finley: write Mesh in finley file format */

/************************************************************************************/

#include "Mesh.h"

namespace finley {

/// writes the mesh to the external file fname using the Finley file format
void Mesh::write(const std::string fname) const
{
  char error_msg[LenErrorMsg_MAX];
  FILE *f;
  int NN,i,j,numDim;

  if (MPIInfo->size >1 ) {
    setError(IO_ERROR,"Mesh_write: only single processor runs are supported.");
    return;

  }
  /* open file */
  f=fopen(fname.c_str(), "w");
  if (f==NULL) {
    sprintf(error_msg,"Mesh_write: Opening file %s for writing failed.",fname.c_str());
    setError(IO_ERROR,error_msg);
    return;
  }

  /* write header */

  fprintf(f, "%s\n", m_name.c_str());
  
  /*  write nodes: */
  
  if (Nodes!=NULL) {
    numDim=getDim();
    fprintf(f,"%1dD-Nodes %d\n", numDim, Nodes->numNodes);
    for (i=0;i<Nodes->numNodes;i++) {
      fprintf(f,"%d %d %d",Nodes->Id[i],Nodes->globalDegreesOfFreedom[i],Nodes->Tag[i]);
      for (j=0;j<numDim;j++) fprintf(f," %20.15e",Nodes->Coordinates[INDEX2(j,i,numDim)]);
      fprintf(f,"\n");
    }
  } else {
    fprintf(f,"0D-Nodes 0\n");
  }
  
  /*  write elements: */

  if (Elements!=NULL) {
    fprintf(f, "%s %d\n",Elements->referenceElementSet->referenceElement->Type->Name,Elements->numElements);
    NN=Elements->numNodes; 
    for (i=0;i<Elements->numElements;i++) {
      fprintf(f,"%d %d",Elements->Id[i],Elements->Tag[i]);
      for (j=0;j<NN;j++) fprintf(f," %d",Nodes->Id[Elements->Nodes[INDEX2(j,i,NN)]]);
      fprintf(f,"\n");
    }
  } else {
    fprintf(f,"Tet4 0\n");
  }

  /*  write face elements: */
  if (FaceElements!=NULL) {
    fprintf(f, "%s %d\n", FaceElements->referenceElementSet->referenceElement->Type->Name,FaceElements->numElements);
    NN=FaceElements->numNodes;
    for (i=0;i<FaceElements->numElements;i++) {
      fprintf(f,"%d %d",FaceElements->Id[i],FaceElements->Tag[i]);
      for (j=0;j<NN;j++) fprintf(f," %d",Nodes->Id[FaceElements->Nodes[INDEX2(j,i,NN)]]);
      fprintf(f,"\n");
    }
  } else {
    fprintf(f,"Tri3 0\n");
  }

  /*  write Contact elements : */
  if (ContactElements!=NULL) {
    fprintf(f, "%s %d\n",ContactElements->referenceElementSet->referenceElement->Type->Name,ContactElements->numElements);
    NN=ContactElements->numNodes;
    for (i=0;i<ContactElements->numElements;i++) {
      fprintf(f,"%d %d",ContactElements->Id[i],ContactElements->Tag[i]);
      for (j=0;j<NN;j++) fprintf(f," %d",Nodes->Id[ContactElements->Nodes[INDEX2(j,i,NN)]]);
      fprintf(f,"\n");
    }
  } else {
    fprintf(f,"Tri3_Contact 0\n");
  }
  
  /*  write points: */
  if (Points!=NULL) {
    fprintf(f, "%s %d\n",Points->referenceElementSet->referenceElement->Type->Name,Points->numElements);
    for (i=0;i<Points->numElements;i++) {
      fprintf(f,"%d %d %d\n",Points->Id[i],Points->Tag[i],Nodes->Id[Points->Nodes[INDEX2(0,i,1)]]);
    }
  } else {
    fprintf(f,"Point1 0\n");
  }

    /*  write tags:*/
    if (tagMap.size()>0) {
        fprintf(f, "Tags\n");
        TagMap::const_iterator it;
        for (it=tagMap.begin(); it!=tagMap.end(); it++) {
            fprintf(f, "%s %d\n", it->first.c_str(), it->second);
        }
    }
    fclose(f);
#ifdef Finley_TRACE
    printf("mesh %s has been written to file %s\n", m_name, fname.c_str());
#endif
}

void Mesh::printInfo(bool full)
{
  int NN,i,j,numDim;

  fprintf(stdout, "PrintMesh_Info running on CPU %d of %d\n",MPIInfo->rank, MPIInfo->size);
  fprintf(stdout, "\tMesh name '%s'\n", m_name.c_str());
  fprintf(stdout, "\tApproximation order %d\n",approximationOrder);
  fprintf(stdout, "\tReduced Approximation order %d\n",reducedApproximationOrder);
  fprintf(stdout, "\tIntegration order %d\n",integrationOrder);
  fprintf(stdout, "\tReduced Integration order %d\n",reducedIntegrationOrder);

  /* write nodes: */
  if (Nodes!=NULL) {
    numDim=getDim();
    fprintf(stdout, "\tNodes: %1dD-Nodes %d\n", numDim, Nodes->numNodes);
    if (full) {
      fprintf(stdout, "\t     Id   Tag  gDOF   gNI grDfI  grNI:  Coordinates\n");
      for (i=0;i<Nodes->numNodes;i++) {
        fprintf(stdout, "\t  %5d %5d %5d %5d %5d %5d: ", Nodes->Id[i], Nodes->Tag[i], Nodes->globalDegreesOfFreedom[i], Nodes->globalNodesIndex[i], Nodes->globalReducedDOFIndex[i], Nodes->globalReducedNodesIndex[i]);
        for (j=0;j<numDim;j++) fprintf(stdout," %20.15e",Nodes->Coordinates[INDEX2(j,i,numDim)]);
        fprintf(stdout,"\n");
      }
    }
  } else {
    fprintf(stdout, "\tNodes: 0D-Nodes 0\n");
  }

  /* write elements: */
  if (Elements!=NULL) {
    int mine=0, overlap=0;
    for (i=0;i<Elements->numElements;i++) {
      if (Elements->Owner[i] == MPIInfo->rank) mine++;
      else overlap++;
    }
    fprintf(stdout, "\tElements: %s %d (TypeId=%d) owner=%d overlap=%d\n",Elements->referenceElementSet->referenceElement->Type->Name,Elements->numElements,Elements->referenceElementSet->referenceElement->Type->TypeId, mine, overlap);
    NN=Elements->numNodes;
    if (full) {
      fprintf(stdout, "\t     Id   Tag Owner Color:  Nodes\n");
      for (i=0;i<Elements->numElements;i++) {
        fprintf(stdout, "\t  %5d %5d %5d %5d: ",Elements->Id[i],Elements->Tag[i],Elements->Owner[i],Elements->Color[i]);
        for (j=0;j<NN;j++) fprintf(stdout," %5d",Nodes->Id[Elements->Nodes[INDEX2(j,i,NN)]]);
        fprintf(stdout,"\n");
      }
    }
  } else {
    fprintf(stdout, "\tElements: Tet4 0\n");
  }

  /* write face elements: */
  if (FaceElements!=NULL) {
    int mine=0, overlap=0;
    for (i=0;i<FaceElements->numElements;i++) {
      if (FaceElements->Owner[i] == MPIInfo->rank) mine++;
      else overlap++;
    }
    fprintf(stdout, "\tFace elements: %s %d (TypeId=%d) owner=%d overlap=%d\n", FaceElements->referenceElementSet->referenceElement->Type->Name,FaceElements->numElements,FaceElements->referenceElementSet->referenceElement->Type->TypeId, mine, overlap);
    NN=FaceElements->numNodes;
    if (full) {
      fprintf(stdout, "\t     Id   Tag Owner Color:  Nodes\n");
      for (i=0;i<FaceElements->numElements;i++) {
        fprintf(stdout, "\t  %5d %5d %5d %5d: ",FaceElements->Id[i],FaceElements->Tag[i],FaceElements->Owner[i],FaceElements->Color[i]);
        for (j=0;j<NN;j++) fprintf(stdout," %5d",Nodes->Id[FaceElements->Nodes[INDEX2(j,i,NN)]]);
        fprintf(stdout,"\n");
      }
    }
  } else {
    fprintf(stdout, "\tFace elements: Tri3 0\n");
  }

  /* write Contact elements : */
  if (ContactElements!=NULL) {
    int mine=0, overlap=0;
    for (i=0;i<ContactElements->numElements;i++) {
      if (ContactElements->Owner[i] == MPIInfo->rank) mine++;
      else overlap++;
    }
    fprintf(stdout, "\tContact elements: %s %d (TypeId=%d) owner=%d overlap=%d\n",ContactElements->referenceElementSet->referenceElement->Type->Name,ContactElements->numElements,ContactElements->referenceElementSet->referenceElement->Type->TypeId, mine, overlap);
    NN=ContactElements->numNodes;
    if (full) {
      fprintf(stdout, "\t     Id   Tag Owner Color:  Nodes\n");
      for (i=0;i<ContactElements->numElements;i++) {
        fprintf(stdout, "\t  %5d %5d %5d %5d: ",ContactElements->Id[i],ContactElements->Tag[i],ContactElements->Owner[i],ContactElements->Color[i]);
        for (j=0;j<NN;j++) fprintf(stdout," %5d",Nodes->Id[ContactElements->Nodes[INDEX2(j,i,NN)]]);
        fprintf(stdout,"\n");
      }
    }
  } else {
    fprintf(stdout, "\tContact elements: Tri3_Contact 0\n");
  }

  /* write points: */
  if (Points!=NULL) {
    int mine=0, overlap=0;
    for (i=0;i<Points->numElements;i++) {
      if (Points->Owner[i] == MPIInfo->rank) mine++;
      else overlap++;
    }
    fprintf(stdout, "\tPoints: %s %d (TypeId=%d) owner=%d overlap=%d\n",Points->referenceElementSet->referenceElement->Type->Name,Points->numElements,Points->referenceElementSet->referenceElement->Type->TypeId, mine, overlap);
    if (full) {
      fprintf(stdout, "\t     Id   Tag Owner Color:  Nodes\n");
      for (i=0;i<Points->numElements;i++) {
        fprintf(stdout, "\t  %5d %5d %5d %5d %5d\n",Points->Id[i],Points->Tag[i],Points->Owner[i],Points->Color[i],Nodes->Id[Points->Nodes[INDEX2(0,i,1)]]);
      }
    }
  } else {
    fprintf(stdout, "\tPoints: Point1 0\n");
  }

    // write tags
    if (tagMap.size()>0) {
        fprintf(stdout, "\tTags:\n");
        TagMap::const_iterator it;
        for (it=tagMap.begin(); it!=tagMap.end(); it++) {
            fprintf(stdout, "\t  %5d %s\n", it->second, it->first.c_str());
        }
    }
}

} // namespace finley

