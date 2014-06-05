
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


/****************************************************************************

  Finley: read mesh from gmsh file

*****************************************************************************/

#include "Mesh.h"
#include <cstdio>

namespace finley {

#define FSCANF_CHECK(scan_ret, reason) { if (scan_ret == EOF) { perror(reason); setError(IO_ERROR,"scan error while reading finley file"); return NULL;} }

#define MAX_numNodes_gmsh 20

Mesh* Mesh::readGmsh(esysUtils::JMPI& mpi_info, const std::string fname, int numDim, int order,
                     int reduced_order, bool optimize, bool useMacroElements)
{
    double version = 1.0;
    int format = 0, size = sizeof(double), scan_ret;
    int numNodes, totalNumElements=0, numTags=0, numNodesPerElement=0, numNodesPerElement2, element_dim=0;
    int e, i0, j, gmsh_type, partition_id, itmp, elementary_id, tag_key;
    int numElements=0, numFaceElements=0, *id=NULL, *tag=NULL, *vertices=NULL;
    char line[LenString_MAX+1], name[LenString_MAX+1];
    char error_msg[LenErrorMsg_MAX];
    double rtmp0, rtmp1;
    const_ReferenceElementSet_ptr refPoints, refContactElements, refFaceElements, refElements;
#ifdef Finley_TRACE
    double time0=timer();
#endif
    FILE * fileHandle_p = NULL;
    ElementTypeId* element_type=NULL;

    resetError();
    if (mpi_info->size>1) {
        setError(IO_ERROR, "reading gmsh files with MPI is not supported yet.");
        return NULL;
    }

    // allocate mesh
    Mesh* mesh_p = new Mesh(fname, numDim, mpi_info);

    // get file handle
    fileHandle_p = fopen(fname.c_str(), "r");
    if (fileHandle_p==NULL) {
        sprintf(error_msg, "Opening Gmsh file %s for reading failed.", fname.c_str());
        setError(IO_ERROR, error_msg);
        return NULL;
    }

    /* start reading */
    while(1) {
       if (!noError()) break;
       /* find line staring with $ */
       do {
         if( ! fgets(line, sizeof(line), fileHandle_p) ) break;
         if(feof(fileHandle_p)) break;
       } while(line[0] != '$');

       if (feof(fileHandle_p)) break;

       /* format */
       if (!strncmp(&line[1], "MeshFormat", 10)) {
         scan_ret = fscanf(fileHandle_p, "%lf %d %d\n", &version, &format, &size);
     FSCANF_CHECK(scan_ret, "fscanf: Mesh_readGmsh");
       }
       /* nodes are read */
       if ( !strncmp(&line[1], "NOD", 3)   ||
            !strncmp(&line[1], "NOE", 3)   ||
            !strncmp(&line[1], "Nodes", 5)    ) {

         scan_ret = fscanf(fileHandle_p, "%d", &numNodes);
     FSCANF_CHECK(scan_ret, "fscanf: Mesh_readGmsh");
         if (!noError()) break;
         mesh_p->Nodes->allocTable(numNodes);
         if (!noError()) break;
         for (i0 = 0; i0 < numNodes; i0++) {
            if (1 == numDim) {
           scan_ret = fscanf(fileHandle_p, "%d %le %le %le\n", &mesh_p->Nodes->Id[i0],
                  &mesh_p->Nodes->Coordinates[INDEX2(0,i0,numDim)],
                         &rtmp0,
                         &rtmp1);
           FSCANF_CHECK(scan_ret, "fscanf: Mesh_readGmsh");
            } else if (2 == numDim) {
           scan_ret = fscanf(fileHandle_p, "%d %le %le %le\n", &mesh_p->Nodes->Id[i0],
                  &mesh_p->Nodes->Coordinates[INDEX2(0,i0,numDim)],
                  &mesh_p->Nodes->Coordinates[INDEX2(1,i0,numDim)],
                         &rtmp0);
           FSCANF_CHECK(scan_ret, "fscanf: Mesh_readGmsh");
            } else if (3 == numDim) {
           scan_ret = fscanf(fileHandle_p, "%d %le %le %le\n", &mesh_p->Nodes->Id[i0],
                  &mesh_p->Nodes->Coordinates[INDEX2(0,i0,numDim)],
                  &mesh_p->Nodes->Coordinates[INDEX2(1,i0,numDim)],
                  &mesh_p->Nodes->Coordinates[INDEX2(2,i0,numDim)]);
           FSCANF_CHECK(scan_ret, "fscanf: Mesh_readGmsh");
            }
            mesh_p->Nodes->globalDegreesOfFreedom[i0]=mesh_p->Nodes->Id[i0];
            mesh_p->Nodes->Tag[i0]=0;
         }
       }
       /* elements */
       else if(!strncmp(&line[1], "ELM", 3) ||
        !strncmp(&line[1], "Elements", 8)) {

         ElementTypeId final_element_type = NoRef;
         ElementTypeId final_face_element_type = NoRef;
         ElementTypeId contact_element_type = NoRef;
         numElements=0;
         numFaceElements=0;
         scan_ret = fscanf(fileHandle_p, "%d", &totalNumElements);
     FSCANF_CHECK(scan_ret, "fscanf: Mesh_readGmsh");

         id=new int[totalNumElements];
         tag=new int[totalNumElements];

         element_type=new ElementTypeId[totalNumElements];
         vertices=new int[totalNumElements*MAX_numNodes_gmsh];
            /* read all in */
            for(e = 0; e < totalNumElements; e++) {
              scan_ret = fscanf(fileHandle_p, "%d %d", &id[e], &gmsh_type);
          FSCANF_CHECK(scan_ret, "fscanf: Mesh_readGmsh");
              switch (gmsh_type) {
                  case 1:  /* line order 1 */
                      element_type[e]=Line2;
                      element_dim=1;
                      numNodesPerElement=2;
                      break;
                  case 2:  /* triangle order 1 */
                      element_type[e]=Tri3;
                      numNodesPerElement= 3;
                      element_dim=2;
                      break;
                  case 3:  /* quadrilateral order 1 */
                      element_type[e]=Rec4;
                      numNodesPerElement= 4;
                      element_dim=2;
                      break;
                  case 4:  /* tetrahedron order 1 */
                      element_type[e]=Tet4;
                      numNodesPerElement= 4;
                      element_dim=3;
                      break;
                  case 5:  /* hexahedron order 1 */
                      element_type[e]=Hex8;
                      numNodesPerElement= 8;
                      element_dim=3;
                      break;
                  case 8:  /* line order 2 */
                      if (useMacroElements) {
                          element_type[e]=Line3Macro;
                      } else {
                          element_type[e]=Line3;
                      }
                      numNodesPerElement= 3;
                      element_dim=1;
                      break;
                  case 9:  /* triangle order 2 */
                      if (useMacroElements) {
                           element_type[e]=Tri6Macro;
                      } else {
                           element_type[e]=Tri6;
                      }
                      numNodesPerElement= 6;
                      element_dim=2;
                      break;
                  case 10:  /* quadrilateral order 2 */
                      if (useMacroElements) {
                          element_type[e]=Rec9Macro;
                      } else {
                          element_type[e]=Rec9;
                      }
                      numNodesPerElement= 9;
                      element_dim=2;
                      break;
                  case 11:  /* tetrahedron order 2 */
                      if (useMacroElements) {
                          element_type[e]=Tet10Macro;
                      } else {
                          element_type[e]=Tet10;
                      }
                      numNodesPerElement= 10;
                      element_dim=3;
                      break;
                  case 16:  /* rectangular order 2 */
                      element_type[e]=Rec8;
                      numNodesPerElement= 8;
                      element_dim=2;
                      break;
                  case 17:  /* hexahedron order 2 */
                      element_type[e]=Hex20;
                      numNodesPerElement= 20;
                      element_dim=3;
                      break;
                  case 15 :  /* point */
                      element_type[e]=Point1;
                      numNodesPerElement= 1;
                      element_dim=0;
                      break;
                  default:
                     element_type[e]=NoRef;
                     sprintf(error_msg,"Unexpected gmsh element type %d in mesh file %s.", gmsh_type, fname.c_str());
                     setError(IO_ERROR,error_msg);
              }
              if (element_dim == numDim) {
                 if (final_element_type == NoRef) {
                    final_element_type = element_type[e];
                 } else if (final_element_type != element_type[e]) {
                     sprintf(error_msg,"Finley can handle a single type of internal elements only.");
                     setError(IO_ERROR,error_msg);
                     break;
                 }
                 numElements++;
              } else if (element_dim == numDim-1) {
                 if (final_face_element_type == NoRef) {
                    final_face_element_type = element_type[e];
                 } else if (final_face_element_type != element_type[e]) {
                     sprintf(error_msg,"Finley can handle a single type of face elements only.");
                     setError(IO_ERROR,error_msg);
                     break;
                 }
                 numFaceElements++;
              }

       if(version <= 1.0){
         scan_ret = fscanf(fileHandle_p, "%d %d %d", &tag[e], &elementary_id, &numNodesPerElement2);
         FSCANF_CHECK(scan_ret, "fscanf: Mesh_readGmsh");
         partition_id = 1;
                if (numNodesPerElement2 != numNodesPerElement) {
                     sprintf(error_msg,"Illegal number of nodes for element %d in mesh file %s.", id[e], fname.c_str());
                     setError(IO_ERROR,error_msg);
                }
       } else {
         scan_ret = fscanf(fileHandle_p, "%d", &numTags);
         FSCANF_CHECK(scan_ret, "fscanf: Mesh_readGmsh");
         elementary_id = tag[e] = partition_id = 1;
                numNodesPerElement2=-1;
         for(j = 0; j < numTags; j++){
           scan_ret = fscanf(fileHandle_p, "%d", &itmp);
           FSCANF_CHECK(scan_ret, "fscanf: Mesh_readGmsh");
           if (j == 0) {
             tag[e] = itmp;
           } else if (j == 1) {
             elementary_id = itmp;
           } else if (j == 2) {
             partition_id = itmp;
                  }
           /* ignore any other tags */
         }
       }
              if (!noError()) break;
              for(j = 0; j < numNodesPerElement; j++) {
        scan_ret = fscanf(fileHandle_p, "%d", &vertices[INDEX2(j,e,MAX_numNodes_gmsh)]);
            FSCANF_CHECK(scan_ret, "fscanf: Mesh_readGmsh");
          }
              /* for tet10 the last two nodes need to be swapped */
              if ((element_type[e]==Tet10) || (element_type[e]==Tet10Macro)) {
                   itmp=vertices[INDEX2(9,e,MAX_numNodes_gmsh)];
                   vertices[INDEX2(9,e,MAX_numNodes_gmsh)]=vertices[INDEX2(8,e,MAX_numNodes_gmsh)];
                   vertices[INDEX2(8,e,MAX_numNodes_gmsh)]=itmp;
              }
            }

            // all elements have been read, now we have to identify the
            // elements for finley
            if (noError()) {
              /* first we have to identify the elements to define Elements and FaceElements */
              if (final_element_type == NoRef) {
                 if (numDim==1) {
                    final_element_type=Line2;
                 } else if (numDim==2) {
                    final_element_type=Tri3;
                 } else if (numDim==3) {
                    final_element_type=Tet4;
                 }
              }
              if (final_face_element_type == NoRef) {
                 if (numDim==1) {
                    final_face_element_type=Point1;
                 } else if (numDim==2) {
                    final_face_element_type=Line2;
                 } else if (numDim==3) {
                    final_face_element_type=Tri3;
                 }
              }
              if (final_face_element_type == Line2) {
                  contact_element_type=Line2_Contact;
              } else  if ( (final_face_element_type == Line3) || (final_face_element_type == Line3Macro) ) {
                  contact_element_type=Line3_Contact;
              } else  if (final_face_element_type == Tri3) {
                  contact_element_type=Tri3_Contact;
              } else  if ( (final_face_element_type == Tri6) || (final_face_element_type == Tri6Macro)) {
                  contact_element_type=Tri6_Contact;
              } else {
                  contact_element_type=Point1_Contact;
              }
              refElements.reset(new ReferenceElementSet(final_element_type, order, reduced_order));
              refFaceElements.reset(new ReferenceElementSet(final_face_element_type, order, reduced_order));
              refContactElements.reset(new ReferenceElementSet(contact_element_type, order, reduced_order));
              refPoints.reset(new ReferenceElementSet(Point1, order, reduced_order));
              mesh_p->Elements=new ElementFile(refElements, mpi_info);
              mesh_p->FaceElements=new ElementFile(refFaceElements, mpi_info);
              mesh_p->ContactElements=new ElementFile(refContactElements, mpi_info);
              mesh_p->Points=new ElementFile(refPoints, mpi_info);
              if (noError()) {
                  mesh_p->Elements->allocTable(numElements);
                  mesh_p->FaceElements->allocTable(numFaceElements);
                  mesh_p->ContactElements->allocTable(0);
                  mesh_p->Points->allocTable(0);
                  if (noError()) {
                      mesh_p->Elements->minColor=0;
                      mesh_p->Elements->maxColor=numElements-1;
                      mesh_p->FaceElements->minColor=0;
                      mesh_p->FaceElements->maxColor=numFaceElements-1;
                      mesh_p->ContactElements->minColor=0;
                      mesh_p->ContactElements->maxColor=0;
                      mesh_p->Points->minColor=0;
                      mesh_p->Points->maxColor=0;
                      numElements=0;
                      numFaceElements=0;
                      for(e = 0; e < totalNumElements; e++) {
                         if (element_type[e] == final_element_type) {
                            mesh_p->Elements->Id[numElements]=id[e];
                            mesh_p->Elements->Tag[numElements]=tag[e];
                            mesh_p->Elements->Color[numElements]=numElements;
                            mesh_p->Elements->Owner[numElements]=0;
                            for (j = 0; j<  mesh_p->Elements->numNodes; ++j)  {
                                  mesh_p->Elements->Nodes[INDEX2(j, numElements, mesh_p->Elements->numNodes)]=vertices[INDEX2(j,e,MAX_numNodes_gmsh)];
                            }
                            numElements++;
                         } else if (element_type[e] == final_face_element_type) {
                            mesh_p->FaceElements->Id[numFaceElements]=id[e];
                            mesh_p->FaceElements->Tag[numFaceElements]=tag[e];
                            mesh_p->FaceElements->Color[numFaceElements]=numFaceElements;
                            mesh_p->FaceElements->Owner[numFaceElements]=0;
                            for (j=0; j<mesh_p->FaceElements->numNodes; ++j) {
                                     mesh_p->FaceElements->Nodes[INDEX2(j, numFaceElements, mesh_p->FaceElements->numNodes)]=vertices[INDEX2(j,e,MAX_numNodes_gmsh)];
                            }
                            numFaceElements++;
                         }
                      }
                 }
              }
         }
         /* and clean up */
         delete[] id;
         delete[] tag;
         delete[] element_type;
         delete[] vertices;
      }
     /* name tags (thanks to Antoine Lefebvre, antoine.lefebvre2@mail.mcgill.ca ) */
     else if (!strncmp(&line[1], "PhysicalNames", 13)) {
         scan_ret = fscanf(fileHandle_p, "%d", &numTags);
         FSCANF_CHECK(scan_ret, "fscanf: Mesh_readGmsh");
         if (! noError()) break;
         for (i0 = 0; i0 < numTags; i0++) {
            scan_ret = fscanf(fileHandle_p, "%d %d %s\n", &itmp, &tag_key, name);
            FSCANF_CHECK(scan_ret, "fscanf: Mesh_readGmsh");
            //if (! (itmp == 2)) setError(IO_ERROR,"Mesh_readGmsh: expecting two entries per physical name.");
            if ( strlen(name) < 3 ) setError(IO_ERROR,"Mesh_readGmsh: illegal tagname (\" missing?)");
            if (! noError()) break;
            name[strlen(name)-1]='\0';
            mesh_p->addTagMap(&name[1], tag_key);
         }
      }

      /* search for end of data block */
      do {
         if (!fgets(line, sizeof(line), fileHandle_p)) {
            sprintf(error_msg,"Unexpected end of file in %s", fname.c_str());
            setError(IO_ERROR,error_msg);
         }
         if (feof(fileHandle_p)) {
            sprintf(error_msg,"Unexpected end of file in %s", fname.c_str());
            setError(IO_ERROR,error_msg);
         }
         if (! noError()) break;
       } while(line[0] != '$');

    } // end of read

    // close file
    fclose(fileHandle_p);
    // clean up
    if (!noError()) {
        delete mesh_p;
        return NULL;
    }
    // resolve id's
    if (noError()) mesh_p->resolveNodeIds();
    // rearrange elements
    if (noError()) mesh_p->prepare(optimize);
    if (!noError()) {
        delete mesh_p;
        return NULL;
    }
    return mesh_p;
}

} // namespace finley

