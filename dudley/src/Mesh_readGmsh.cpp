
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

using escript::IOError;

namespace dudley {

#define MAX_numNodes_gmsh 20

/// reads a mesh from a gmsh file of name filename
Mesh* Mesh::readGmsh(escript::JMPI mpiInfo, const std::string& filename,
                     int numDim, bool optimize)
{
    double version = 1.0;
    int format = 0, size = sizeof(double), scan_ret;
    dim_t numNodes, totalNumElements = 0, numTags = 0;
    int numNodesPerElement = 0, numNodesPerElement2, element_dim = 0;
    index_t gmsh_type, partition_id, itmp, elementary_id;
    index_t numElements = 0, numFaceElements = 0, *id = NULL, *tag = NULL, *vertices = NULL;
    char line[1024];
    double rtmp0, rtmp1;

    if (mpiInfo->size > 1)
        throw DudleyException("reading GMSH with MPI is not supported yet.");

    // allocate mesh
    Mesh* mesh = new Mesh(filename, numDim, mpiInfo);

    // open file
    FILE* fileHandle = fopen(filename.c_str(), "r");
    if (fileHandle == NULL) {
        std::stringstream ss;
        ss << "Opening gmsh file " << filename << " for reading failed.";
        throw IOError(ss.str());
    }

    // start reading
    while (1) {
        // find line staring with $
        do {
            if (!fgets(line, sizeof(line), fileHandle))
                break;
            if (feof(fileHandle))
                break;
        } while (line[0] != '$');

        if (feof(fileHandle))
            break;

        // format
        if (!strncmp(&line[1], "MeshFormat", 10)) {
            scan_ret = fscanf(fileHandle, "%lf %d %d\n", &version, &format, &size);
            if (scan_ret == EOF)
                throw IOError("readGmsh: early EOF while reading file");
        }
        // nodes
        if (!strncmp(&line[1], "NOD", 3) || !strncmp(&line[1], "NOE", 3) || !strncmp(&line[1], "Nodes", 5)) {
            scan_ret = fscanf(fileHandle, "%d", &numNodes);
            if (scan_ret == EOF)
                throw IOError("readGmsh: early EOF while reading file");
            mesh->Nodes->allocTable(numNodes);
            for (index_t i0 = 0; i0 < numNodes; i0++) {
                if (1 == numDim) {
                    scan_ret = fscanf(fileHandle, "%d %le %le %le\n",
                            &mesh->Nodes->Id[i0],
                            &mesh->Nodes->Coordinates[INDEX2(0, i0, numDim)],
                            &rtmp0, &rtmp1);
                    if (scan_ret == EOF)
                        throw IOError("readGmsh: early EOF while reading file");
                } else if (2 == numDim) {
                    scan_ret = fscanf(fileHandle, "%d %le %le %le\n",
                            &mesh->Nodes->Id[i0],
                            &mesh->Nodes->Coordinates[INDEX2(0, i0, numDim)],
                            &mesh->Nodes->Coordinates[INDEX2(1, i0, numDim)],
                            &rtmp0);
                    if (scan_ret == EOF)
                        throw IOError("readGmsh: early EOF while reading file");
                } else if (3 == numDim) {
                    scan_ret = fscanf(fileHandle, "%d %le %le %le\n",
                            &mesh->Nodes->Id[i0],
                            &mesh->Nodes->Coordinates[INDEX2(0, i0, numDim)],
                            &mesh->Nodes->Coordinates[INDEX2(1, i0, numDim)],
                            &mesh->Nodes->Coordinates[INDEX2(2, i0, numDim)]);
                    if (scan_ret == EOF)
                        throw IOError("readGmsh: early EOF while reading file");
                }
                mesh->Nodes->globalDegreesOfFreedom[i0] = mesh->Nodes->Id[i0];
                mesh->Nodes->Tag[i0] = 0;
            }
        }
        // elements
        else if (!strncmp(&line[1], "ELM", 3) || !strncmp(&line[1], "Elements", 8)) {
            ElementTypeId final_element_type = Dudley_NoRef;
            ElementTypeId final_face_element_type = Dudley_NoRef;
            numElements = 0;
            numFaceElements = 0;
            scan_ret = fscanf(fileHandle, "%d", &totalNumElements);
            if (scan_ret == EOF)
                throw IOError("readGmsh: early EOF while reading file");

            id = new index_t[totalNumElements];
            tag = new index_t[totalNumElements];

            ElementTypeId* element_type = new ElementTypeId[totalNumElements];
            vertices = new index_t[totalNumElements * MAX_numNodes_gmsh];
            // read all in
            for (index_t e = 0; e < totalNumElements; e++) {
                scan_ret = fscanf(fileHandle, "%d %d", &id[e], &gmsh_type);
                if (scan_ret == EOF)
                    throw IOError("readGmsh: early EOF while reading file");
                switch (gmsh_type) {
                    case 1: // line order 1
                        element_type[e] = Dudley_Line2;
                        element_dim = 1;
                        numNodesPerElement = 2;
                        break;
                    case 2: // triangle order 1
                        element_type[e] = Dudley_Tri3;
                        numNodesPerElement = 3;
                        element_dim = 2;
                        break;
                    case 4: // tetrahedron order 1
                        element_type[e] = Dudley_Tet4;
                        numNodesPerElement = 4;
                        element_dim = 3;
                        break;
                    case 15: // point
                        element_type[e] = Dudley_Point1;
                        numNodesPerElement = 1;
                        element_dim = 0;
                        break;
                    default: {
                        std::stringstream ss;
                        ss << "Unexpected gmsh element type " << gmsh_type
                           << " in mesh file " << filename;
                        throw IOError(ss.str());
                    }
                }
                if (element_dim == numDim) {
                    if (final_element_type == Dudley_NoRef) {
                        final_element_type = element_type[e];
                    } else if (final_element_type != element_type[e]) {
                        throw IOError("Dudley can handle a single type of "
                                      "internal elements only.");
                    }
                    numElements++;
                } else if (element_dim == numDim - 1) {
                    if (final_face_element_type == Dudley_NoRef) {
                        final_face_element_type = element_type[e];
                    } else if (final_face_element_type != element_type[e]) {
                        throw IOError("Dudley can handle a single type of "
                                      "face elements only.");
                    }
                    numFaceElements++;
                }

                if (version <= 1.0) {
                    scan_ret = fscanf(fileHandle, "%d %d %d", &tag[e],
                                      &elementary_id, &numNodesPerElement2);
                    if (scan_ret == EOF)
                        throw IOError("readGmsh: early EOF while reading file");
                    partition_id = 1;
                    if (numNodesPerElement2 != numNodesPerElement) {
                        std::stringstream ss;
                        ss << "Illegal number of nodes for element " << id[e]
                            << " in mesh file " << filename;
                        throw IOError(ss.str());
                    }
                } else {
                    scan_ret = fscanf(fileHandle, "%d", &numTags);
                    if (scan_ret == EOF)
                        throw IOError("readGmsh: early EOF while reading file");
                    elementary_id = tag[e] = partition_id = 1;
                    numNodesPerElement2 = -1;
                    for (int j = 0; j < numTags; j++) {
                        scan_ret = fscanf(fileHandle, "%d", &itmp);
                        if (scan_ret == EOF)
                            throw IOError("readGmsh: early EOF while reading file");
                        if (j == 0) {
                            tag[e] = itmp;
                        } else if (j == 1) {
                            elementary_id = itmp;
                        } else if (j == 2) {
                            partition_id = itmp;
                        }
                        // ignore any other tags
                    }
                }
                for (int j = 0; j < numNodesPerElement; j++) {
                    scan_ret = fscanf(fileHandle, "%d",
                            &vertices[INDEX2(j, e, MAX_numNodes_gmsh)]);
                    if (scan_ret == EOF)
                        throw IOError("readGmsh: early EOF while reading file");
                }
            }
            // all elements have been read, now we have to identify the
            // dudley elements to define Elements and FaceElements
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
            mesh->Elements = new ElementFile(final_element_type, mpiInfo);
            mesh->FaceElements = new ElementFile(final_face_element_type, mpiInfo);
            mesh->Points = new ElementFile(Dudley_Point1, mpiInfo);
            mesh->Elements->allocTable(numElements);
            mesh->FaceElements->allocTable(numFaceElements);
            mesh->Points->allocTable(0);
            mesh->Elements->minColor = 0;
            mesh->Elements->maxColor = numElements - 1;
            mesh->FaceElements->minColor = 0;
            mesh->FaceElements->maxColor = numFaceElements - 1;
            mesh->Points->minColor = 0;
            mesh->Points->maxColor = 0;
            numElements = 0;
            numFaceElements = 0;
            for (index_t e = 0; e < totalNumElements; e++) {
                if (element_type[e] == final_element_type) {
                    mesh->Elements->Id[numElements] = id[e];
                    mesh->Elements->Tag[numElements] = tag[e];
                    mesh->Elements->Color[numElements] = numElements;
                    mesh->Elements->Owner[numElements] = 0;
                    for (int j = 0; j < mesh->Elements->numNodes; ++j) {
                        mesh->Elements->Nodes[INDEX2(j, numElements,
                                                 mesh->Elements->numNodes)] =
                            vertices[INDEX2(j, e, MAX_numNodes_gmsh)];
                    }
                    numElements++;
                } else if (element_type[e] == final_face_element_type) {
                    mesh->FaceElements->Id[numFaceElements] = id[e];
                    mesh->FaceElements->Tag[numFaceElements] = tag[e];
                    mesh->FaceElements->Color[numFaceElements] = numFaceElements;
                    mesh->FaceElements->Owner[numFaceElements] = 0;
                    for (int j = 0; j < mesh->FaceElements->numNodes; ++j)
                    {
                        mesh->FaceElements->Nodes[INDEX2(j, numFaceElements,
                                              mesh->FaceElements->numNodes)] =
                            vertices[INDEX2(j, e, MAX_numNodes_gmsh)];
                    }
                    numFaceElements++;
                }
            }
            // and clean up
            delete[] id;
            delete[] tag;
            delete[] element_type;
            delete[] vertices;
        } // end reading elements  
        // name tags (thanks to Antoine Lefebvre, antoine.lefebvre2@mail.mcgill.ca )
        else if (!strncmp(&line[1], "PhysicalNames", 13)) {
            char name[1024];
            index_t tag_key;
            scan_ret = fscanf(fileHandle, "%d", &numTags);
            if (scan_ret == EOF)
                throw IOError("readGmsh: early EOF while reading file");
            for (int i0 = 0; i0 < numTags; i0++) {
                scan_ret = fscanf(fileHandle, "%d %d %s\n", &itmp, &tag_key, name);
                if (scan_ret == EOF)
                    throw IOError("readGmsh: early EOF while reading file");
                if (itmp != 2)
                    throw IOError("readGmsh: expecting two entries per physical name.");
                if (strlen(name) < 3)
                    throw IOError("readGmsh: illegal tagname (\" missing?)");
                name[strlen(name)-1]='\0';
                mesh->addTagMap(&name[1], tag_key);
            }
        }
        // search for end of data block
        do {
            if (!fgets(line, sizeof(line), fileHandle)) {
                std::stringstream ss;
                ss << "Unexpected end of file in " << filename;
                throw IOError(ss.str());
            }
            if (feof(fileHandle)) {
                std::stringstream ss;
                ss << "Unexpected end of file in " << filename;
                throw IOError(ss.str());
            }
        } while (line[0] != '$');
    }

    fclose(fileHandle);
    mesh->resolveNodeIds();
    mesh->prepare(optimize);
    return mesh;
}

} // namespace dudley

