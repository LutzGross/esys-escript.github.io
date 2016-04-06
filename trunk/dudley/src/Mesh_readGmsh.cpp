
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
    int format = 0, size = sizeof(double);
    dim_t numNodes, totalNumElements = 0;
    int numTags = 0;
    int numNodesPerElement = 0, numNodesPerElement2, element_dim = 0;
    int gmsh_type, partition_id, itmp, elementary_id;
    index_t numElements = 0, numFaceElements = 0, *id = NULL, *vertices = NULL;
    int* tag = NULL;
    std::string line;

    if (mpiInfo->size > 1)
        throw DudleyException("reading gmsh with MPI is not supported yet.");

    // allocate mesh
    Mesh* mesh = new Mesh(filename, numDim, mpiInfo);

    // open file
    std::ifstream fileHandle(filename);
    if (!fileHandle.good()) {
        std::stringstream ss;
        ss << "Opening gmsh file " << filename << " for reading failed.";
        throw IOError(ss.str());
    }

    // start reading
    while (1) {
        // find line staring with $
        do {
            std::getline(fileHandle, line);
            if (!fileHandle.good())
                break;
        } while (line[0] != '$');

        if (fileHandle.eof())
            break;

        // format
        if (line.substr(1,10) == "MeshFormat") {
            std::getline(fileHandle, line);
            if (fileHandle.eof())
                throw IOError("readGmsh: early EOF while reading file");
            std::stringstream ss(line);
            ss >> version >> format >> size;
        } else if (line.substr(1,3) == "NOD" || line.substr(1,3) == "NOE" ||
                   line.substr(1,5) == "Nodes") {
            // nodes
            std::getline(fileHandle, line);
            if (fileHandle.eof())
                throw IOError("readGmsh: early EOF while reading file");
            numNodes = std::stol(line);
            mesh->Nodes->allocTable(numNodes);
            for (index_t i0 = 0; i0 < numNodes; i0++) {
                std::getline(fileHandle, line);
                if (!fileHandle.good())
                    throw IOError("readGmsh: early EOF while reading file");
                std::stringstream ss(line);
                ss >> mesh->Nodes->Id[i0]
                   >> mesh->Nodes->Coordinates[INDEX2(0, i0, numDim)];
                if (numDim > 1)
                    ss >> mesh->Nodes->Coordinates[INDEX2(1, i0, numDim)];
                if (numDim > 2)
                    ss >> mesh->Nodes->Coordinates[INDEX2(2, i0, numDim)];
                mesh->Nodes->globalDegreesOfFreedom[i0] = mesh->Nodes->Id[i0];
                mesh->Nodes->Tag[i0] = 0;
            }
        } else if (line.substr(1,3) == "ELM" || line.substr(1,8) == "Elements") {
            // elements
            ElementTypeId final_element_type = Dudley_NoRef;
            ElementTypeId final_face_element_type = Dudley_NoRef;
            numElements = 0;
            numFaceElements = 0;
            std::getline(fileHandle, line);
            if (fileHandle.eof())
                throw IOError("readGmsh: early EOF while reading file");
            totalNumElements = std::stol(line);

            id = new index_t[totalNumElements];
            tag = new int[totalNumElements];

            ElementTypeId* element_type = new ElementTypeId[totalNumElements];
            vertices = new index_t[totalNumElements * MAX_numNodes_gmsh];
            // read all in
            for (index_t e = 0; e < totalNumElements; e++) {
                std::getline(fileHandle, line);
                if (fileHandle.eof())
                    throw IOError("readGmsh: early EOF while reading file");
                std::stringstream ss(line);
                ss >> id[e];
                ss >> gmsh_type;
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
                        std::stringstream ss2;
                        ss2 << "Unexpected gmsh element type " << gmsh_type
                           << " in mesh file " << filename;
                        throw IOError(ss2.str());
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
                    ss >> tag[e] >> elementary_id >> numNodesPerElement2;
                    partition_id = 1;
                    if (numNodesPerElement2 != numNodesPerElement) {
                        std::stringstream ss;
                        ss << "Illegal number of nodes for element " << id[e]
                            << " in mesh file " << filename;
                        throw IOError(ss.str());
                    }
                } else {
                    ss >> numTags;
                    elementary_id = tag[e] = partition_id = 1;
                    numNodesPerElement2 = -1;
                    for (int j = 0; j < numTags; j++) {
                        ss >> itmp;
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
                    ss >> vertices[INDEX2(j, e, MAX_numNodes_gmsh)];
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
        } else if (line.substr(1,13) == "PhysicalNames") {
            // name tags (thanks to Antoine Lefebvre,
            // antoine.lefebvre2@mail.mcgill.ca)
            std::getline(fileHandle, line);
            if (fileHandle.eof())
                throw IOError("readGmsh: early EOF while reading file");
            numTags = std::stoi(line);
            for (int i0 = 0; i0 < numTags; i0++) {
                std::getline(fileHandle, line);
                if (fileHandle.eof())
                    throw IOError("readGmsh: early EOF while reading file");
                std::stringstream ss(line);
                int tag_key;
                ss >> itmp >> tag_key;
                std::string name = line.substr((int)ss.tellg()+1);
                if (itmp != 2)
                    throw IOError("readGmsh: expecting two entries per physical name.");
                if (name.length() < 3)
                    throw IOError("readGmsh: illegal tagname (\" missing?)");
                name = name.substr(1, name.length()-2);
                mesh->addTagMap(name, tag_key);
            }
        }
        // search for end of data block
        do {
            std::getline(fileHandle, line);
            if (fileHandle.eof()) {
                std::stringstream ss;
                ss << "Unexpected end of file in " << filename;
                throw IOError(ss.str());
            }
        } while (line[0] != '$');
    }

    fileHandle.close();
    mesh->resolveNodeIds();
    mesh->prepare(optimize);
    return mesh;
}

} // namespace dudley

