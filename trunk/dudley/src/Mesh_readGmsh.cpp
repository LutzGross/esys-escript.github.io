
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include "Mesh.h"

#include <cstdio>

namespace dudley {

#define MAX_numNodes_gmsh 20

/*  reads a mesh from a Dudley file of name fname */
Dudley_Mesh *Dudley_Mesh_readGmsh(char *fname, index_t numDim, index_t order, index_t reduced_order, bool optimize)
{
    double version = 1.0;
    int format = 0, size = sizeof(double), scan_ret;
    dim_t numNodes, totalNumElements = 0, numTags = 0, numNodesPerElement = 0, numNodesPerElement2, element_dim = 0;
    index_t e, i0, j, gmsh_type, partition_id, itmp, elementary_id;
    index_t numElements = 0, numFaceElements = 0, *id = NULL, *tag = NULL, *vertices = NULL;
    Dudley_Mesh *mesh_p = NULL;
    char line[1024];
    char error_msg[1024];
    double rtmp0, rtmp1;
    FILE *fileHandle_p = NULL;
    Dudley_ElementTypeId *element_type = NULL;

    /* No! Bad! take a parameter for this */
    escript::JMPI mpi_info = escript::makeInfo(MPI_COMM_WORLD);
    if (mpi_info->size > 1)
        throw DudleyException("reading GMSH with MPI is not supported yet.");

    /* allocate mesh */
    mesh_p = Dudley_Mesh_alloc(fname, numDim, mpi_info);

    /* get file handle */
    fileHandle_p = fopen(fname, "r");
    if (fileHandle_p == NULL) {
        sprintf(error_msg, "Opening Gmsh file %s for reading failed.", fname);
        throw DudleyException(error_msg);
    }

    /* start reading */
    while (1) {
        /* find line staring with $ */
        do {
            if (!fgets(line, sizeof(line), fileHandle_p))
                break;
            if (feof(fileHandle_p))
                break;
        }
        while (line[0] != '$');

        if (feof(fileHandle_p))
            break;

        /* format */
        if (!strncmp(&line[1], "MeshFormat", 10)) {
            scan_ret = fscanf(fileHandle_p, "%lf %d %d\n", &version, &format, &size);
            if (scan_ret == EOF)
                throw DudleyException("readGmsh: early EOF while reading file");
        }
        /* nodes are read */
        if (!strncmp(&line[1], "NOD", 3) || !strncmp(&line[1], "NOE", 3) || !strncmp(&line[1], "Nodes", 5)) {

            scan_ret = fscanf(fileHandle_p, "%d", &numNodes);
            if (scan_ret == EOF)
                throw DudleyException("readGmsh: early EOF while reading file");
            Dudley_NodeFile_allocTable(mesh_p->Nodes, numNodes);
            for (i0 = 0; i0 < numNodes; i0++)
            {
                if (1 == numDim)
                {
                    scan_ret = fscanf(fileHandle_p, "%d %le %le %le\n", &mesh_p->Nodes->Id[i0],
                                      &mesh_p->Nodes->Coordinates[INDEX2(0, i0, numDim)], &rtmp0, &rtmp1);
                    if (scan_ret == EOF)
                        throw DudleyException("readGmsh: early EOF while reading file");
                }
                else if (2 == numDim)
                {
                    scan_ret = fscanf(fileHandle_p, "%d %le %le %le\n", &mesh_p->Nodes->Id[i0],
                                      &mesh_p->Nodes->Coordinates[INDEX2(0, i0, numDim)],
                                      &mesh_p->Nodes->Coordinates[INDEX2(1, i0, numDim)], &rtmp0);
                    if (scan_ret == EOF)
                        throw DudleyException("readGmsh: early EOF while reading file");
                }
                else if (3 == numDim)
                {
                    scan_ret = fscanf(fileHandle_p, "%d %le %le %le\n", &mesh_p->Nodes->Id[i0],
                                      &mesh_p->Nodes->Coordinates[INDEX2(0, i0, numDim)],
                                      &mesh_p->Nodes->Coordinates[INDEX2(1, i0, numDim)],
                                      &mesh_p->Nodes->Coordinates[INDEX2(2, i0, numDim)]);
                    if (scan_ret == EOF)
                        throw DudleyException("readGmsh: early EOF while reading file");
                }
                mesh_p->Nodes->globalDegreesOfFreedom[i0] = mesh_p->Nodes->Id[i0];
                mesh_p->Nodes->Tag[i0] = 0;
            }
        }
        /* elements */
        else if (!strncmp(&line[1], "ELM", 3) || !strncmp(&line[1], "Elements", 8))
        {
            Dudley_ElementTypeId final_element_type = Dudley_NoRef;
            Dudley_ElementTypeId final_face_element_type = Dudley_NoRef;
            numElements = 0;
            numFaceElements = 0;
            scan_ret = fscanf(fileHandle_p, "%d", &totalNumElements);
            if (scan_ret == EOF)
                throw DudleyException("readGmsh: early EOF while reading file");

            id = new index_t[totalNumElements];
            tag = new index_t[totalNumElements];

            element_type = new Dudley_ElementTypeId[totalNumElements];
            vertices = new index_t[totalNumElements * MAX_numNodes_gmsh];
            /* read all in */
            for (e = 0; e < totalNumElements; e++)
            {
                scan_ret = fscanf(fileHandle_p, "%d %d", &id[e], &gmsh_type);
                if (scan_ret == EOF)
                    throw DudleyException("readGmsh: early EOF while reading file");
                switch (gmsh_type)
                {
                case 1: /* line order 1 */
                    element_type[e] = Dudley_Line2;
                    element_dim = 1;
                    numNodesPerElement = 2;
                    break;
                case 2: /* triangle order 1 */
                    element_type[e] = Dudley_Tri3;
                    numNodesPerElement = 3;
                    element_dim = 2;
                    break;
                case 4: /* tetrahedron order 1 */
                    element_type[e] = Dudley_Tet4;
                    numNodesPerElement = 4;
                    element_dim = 3;
                    break;
                case 15:        /* point */
                    element_type[e] = Dudley_Point1;
                    numNodesPerElement = 1;
                    element_dim = 0;
                    break;
                default:
                    element_type[e] = Dudley_NoRef;
                    sprintf(error_msg, "Unexpected gmsh element type %d in mesh file %s.", gmsh_type, fname);
                    throw DudleyException(error_msg);
                }
                if (element_dim == numDim)
                {
                    if (final_element_type == Dudley_NoRef)
                    {
                        final_element_type = element_type[e];
                    }
                    else if (final_element_type != element_type[e])
                    {
                        sprintf(error_msg, "Dudley can handle a single type of internal elements only.");
                        throw DudleyException(error_msg);
                    }
                    numElements++;
                }
                else if (element_dim == numDim - 1)
                {
                    if (final_face_element_type == Dudley_NoRef) {
                        final_face_element_type = element_type[e];
                    } else if (final_face_element_type != element_type[e]) {
                        sprintf(error_msg, "Dudley can handle a single type of face elements only.");
                        throw DudleyException(error_msg);
                    }
                    numFaceElements++;
                }

                if (version <= 1.0) {
                    scan_ret = fscanf(fileHandle_p, "%d %d %d", &tag[e], &elementary_id, &numNodesPerElement2);
                    if (scan_ret == EOF)
                        throw DudleyException("readGmsh: early EOF while reading file");
                    partition_id = 1;
                    if (numNodesPerElement2 != numNodesPerElement) {
                        sprintf(error_msg, "Illegal number of nodes for element %d in mesh file %s.", id[e],
                                fname);
                        throw DudleyException(error_msg);
                    }
                } else {
                    scan_ret = fscanf(fileHandle_p, "%d", &numTags);
                    if (scan_ret == EOF)
                        throw DudleyException("readGmsh: early EOF while reading file");
                    elementary_id = tag[e] = partition_id = 1;
                    numNodesPerElement2 = -1;
                    for (j = 0; j < numTags; j++) {
                        scan_ret = fscanf(fileHandle_p, "%d", &itmp);
                        if (scan_ret == EOF)
                            throw DudleyException("readGmsh: early EOF while reading file");
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
                for (j = 0; j < numNodesPerElement; j++)
                {
                    scan_ret = fscanf(fileHandle_p, "%d", &vertices[INDEX2(j, e, MAX_numNodes_gmsh)]);
                    if (scan_ret == EOF)
                        throw DudleyException("readGmsh: early EOF while reading file");
                }
            }
            /* all elements have been read, now we have to identify the elements for dudley */

            /* first we have to identify the elements to define Elements and FaceElements */
            if (final_element_type == Dudley_NoRef) {
                if (numDim == 1) {
                    final_element_type = Dudley_Line2;
                } else if (numDim == 2) {
                    final_element_type = Dudley_Tri3;
                } else if (numDim == 3) {
                    final_element_type = Dudley_Tet4;
                }
            }
            if (final_face_element_type == Dudley_NoRef) {
                if (numDim == 1) {
                    final_face_element_type = Dudley_Point1;
                } else if (numDim == 2) {
                    final_face_element_type = Dudley_Line2;
                } else if (numDim == 3) {
                    final_face_element_type = Dudley_Tri3;
                }
            }
            mesh_p->Elements = Dudley_ElementFile_alloc(final_element_type, mpi_info);
            mesh_p->FaceElements = Dudley_ElementFile_alloc(final_face_element_type, mpi_info);
            mesh_p->Points = Dudley_ElementFile_alloc(Dudley_Point1, mpi_info);
            Dudley_ElementFile_allocTable(mesh_p->Elements, numElements);
            Dudley_ElementFile_allocTable(mesh_p->FaceElements, numFaceElements);
            Dudley_ElementFile_allocTable(mesh_p->Points, 0);
            mesh_p->Elements->minColor = 0;
            mesh_p->Elements->maxColor = numElements - 1;
            mesh_p->FaceElements->minColor = 0;
            mesh_p->FaceElements->maxColor = numFaceElements - 1;
            mesh_p->Points->minColor = 0;
            mesh_p->Points->maxColor = 0;
            numElements = 0;
            numFaceElements = 0;
            for (e = 0; e < totalNumElements; e++)
            {
                if (element_type[e] == final_element_type)
                {
                    mesh_p->Elements->Id[numElements] = id[e];
                    mesh_p->Elements->Tag[numElements] = tag[e];
                    mesh_p->Elements->Color[numElements] = numElements;
                    mesh_p->Elements->Owner[numElements] = 0;
                    for (j = 0; j < mesh_p->Elements-> /*referenceElementSet-> */ numNodes; ++j)
                    {
                        mesh_p->Elements->Nodes[INDEX2
                                                (j, numElements,
                                                 mesh_p->
                                                 Elements-> /*referenceElementSet-> */ numNodes)] =
                            vertices[INDEX2(j, e, MAX_numNodes_gmsh)];
                    }
                    numElements++;
                }
                else if (element_type[e] == final_face_element_type)
                {
                    mesh_p->FaceElements->Id[numFaceElements] = id[e];
                    mesh_p->FaceElements->Tag[numFaceElements] = tag[e];
                    mesh_p->FaceElements->Color[numFaceElements] = numFaceElements;
                    mesh_p->FaceElements->Owner[numFaceElements] = 0;
                    for (j = 0; j < mesh_p->FaceElements-> /*referenceElementSet-> */ numNodes; ++j)
                    {
                        mesh_p->FaceElements->Nodes[INDEX2
                                                    (j, numFaceElements,
                                                     mesh_p->
                                                     FaceElements-> /*referenceElementSet-> */
                                                     numNodes)] =
                            vertices[INDEX2(j, e, MAX_numNodes_gmsh)];
                    }
                    numFaceElements++;
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
            char name[1024];
            index_t tag_key;
            scan_ret = fscanf(fileHandle_p, "%d", &numTags);
            if (scan_ret == EOF)
                throw DudleyException("readGmsh: early EOF while reading file");
            for (i0 = 0; i0 < numTags; i0++) {
                scan_ret = fscanf(fileHandle_p, "%d %d %s\n", &itmp, &tag_key, name);
                if (scan_ret == EOF)
                    throw DudleyException("readGmsh: early EOF while reading file");
                if (itmp != 2)
                    throw DudleyException("Dudley_Mesh_readGmsh: expecting two entries per physical name.");
                if (strlen(name) < 3)
                    throw DudleyException("Dudley_Mesh_readGmsh: illegal tagname (\" missing?)");
                name[strlen(name)-1]='\0';
                Dudley_Mesh_addTagMap(mesh_p,&name[1],tag_key);
            }
          }
        /* search for end of data block */
        do
        {
            if (!fgets(line, sizeof(line), fileHandle_p)) {
                sprintf(error_msg, "Unexpected end of file in %s", fname);
                throw DudleyException(error_msg);
            }
            if (feof(fileHandle_p)) {
                sprintf(error_msg, "Unexpected end of file in %s", fname);
                throw DudleyException(error_msg);
            }
        }
        while (line[0] != '$');
    }

    fclose(fileHandle_p);
    Dudley_Mesh_resolveNodeIds(mesh_p);
    Dudley_Mesh_prepare(mesh_p, optimize);
    return mesh_p;
}

} // namespace dudley

