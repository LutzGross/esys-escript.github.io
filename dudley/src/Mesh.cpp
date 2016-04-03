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

namespace dudley {

/// constructor
Mesh::Mesh(const std::string name, int numDim, escript::JMPI mpi_info) :
    m_name(name),
    // order of shapeFunctions is always 1 in Dudley
    approximationOrder(1),
    integrationOrder(2),
    reducedIntegrationOrder(0),
    Elements(NULL),
    FaceElements(NULL),
    Points(NULL),
    MPIInfo(mpi_info)
{
    // allocate node table
    Nodes = new NodeFile(numDim, mpi_info);
}

/// destructor
Mesh::~Mesh()
{
    delete Nodes;
    delete FaceElements;
    delete Elements;
    delete Points;
}

void Mesh::setElements(ElementFile* elements)
{
    delete Elements;
    Elements = elements;
}

void Mesh::setFaceElements(ElementFile* elements)
{
    delete FaceElements;
    FaceElements = elements;
}

void Mesh::setPoints(ElementFile* elements)
{
    delete Points;
    Points = elements;
}

void Mesh::createMappings(const std::vector<index_t>& dofDist,
                          const std::vector<index_t>& nodeDist)
{
    Nodes->createNodeMappings(dofDist, nodeDist);
}

void Mesh::markNodes(std::vector<short>& mask, index_t offset) const
{
    Elements->markNodes(mask, offset);
    FaceElements->markNodes(mask, offset);
    Points->markNodes(mask, offset);
}

void Mesh::relabelElementNodes(const index_t* newNode, index_t offset)
{
    Elements->relabelNodes(newNode, offset);
    FaceElements->relabelNodes(newNode, offset);
    Points->relabelNodes(newNode, offset);
}

void Mesh::setCoordinates(const escript::Data& newX)
{
    Nodes->setCoordinates(newX);
}

void Mesh::addTagMap(const std::string& name, int tag_key)
{
    tagMap[name] = tag_key;
}

int Mesh::getTag(const std::string& name) const
{
    TagMap::const_iterator it = tagMap.find(name);
    if (it == tagMap.end()) {
        std::stringstream ss;
        ss << "getTag: unknown tag name " << name << ".";
        const std::string errorMsg(ss.str());
        throw escript::ValueError(errorMsg);
    }
    return it->second;
}

bool Mesh::isValidTagName(const std::string& name) const
{
    return (tagMap.count(name) > 0);
}

} // namespace dudley

