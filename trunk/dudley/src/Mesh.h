
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

#ifndef __DUDLEY_MESH_H__
#define __DUDLEY_MESH_H__

/****************************************************************************

   Dudley: Mesh

   A mesh is built from nodes and elements which are describing the
   domain, the surface and point sources (the latter are needed to
   establish links with other codes, in particular to particle
   codes). The nodes are stored in a NodeFile and elements in an
   ElementFile. A NodeFile and three ElementFiles containing the
   elements describe the domain, surface and point sources,
   respectively. Notice that the surface elements do not necessarily
   cover the entire surface of the domain.

   The element type is fixed by the reference element, see
   ReferenceElement.h. The numbering of the nodes starts with 0.

   Important: it is assumed that every node appears in at least
   one element or surface element and that any node used in an
   element, surface element or as a point is specified in the
   NodeFile, see also Dudley_resolveNodeIds.

   In some cases it is useful to refer to a mesh entirly built from
   order 1 (=linear) elements. The linear version of the mesh can be
   accessed by referning to the first few nodes of each element
   (thanks to the way the nodes are ordered). As the numbering of
   these nodes is not continuous a relabeling vector is introduced
   in the NodeFile. This feature is not fully implemented yet.

   All nodes and elements are tagged. The tag allows to group nodes and
   elements. A typical application is to mark surface elements on a
   certain portion of the domain with the same tag. All these surface
   elements can then assigned the same value e.g. for the pressure.

   The spatial dimensionality is determined by the type of elements
   used and can be queried using the function Dudley_Mesh_getDim.
   Notice that the element type also determines the type of surface
   elements to be used.

*****************************************************************************/

#include "Dudley.h"
#include "ElementFile.h"
#include "NodeFile.h"
#include "Util.h"

#include <paso/SystemMatrixPattern.h>

#ifdef USE_TRILINOS
#include <trilinoswrap/types.h>
#endif

namespace dudley {

typedef std::map<std::string, int> TagMap;

class Mesh
{
public:
    Mesh(const std::string name, int numDim, escript::JMPI mpi_info);
    ~Mesh();

    static Mesh* read(escript::JMPI mpiInfo, const std::string& filename,
                      bool optimize);
    static Mesh* readGmsh(escript::JMPI mpiInfo, const std::string& filename,
                          int numDim, bool optimize);

    /// writes the mesh to the external file filename using the Dudley file
    /// format
    void write(const std::string& filename) const;

    int getDim() const { return Nodes->numDim; }
    int getStatus() const { return Nodes->status; }

    void addPoints(int numPoints, const double *points_ptr, const int *tags_ptr);
    void addTagMap(const std::string& name, int tag_key);
    int getTag(const std::string& name) const;
    bool isValidTagName(const std::string& name) const;

    /// returns a reference to the paso matrix pattern
    paso::SystemMatrixPattern_ptr getPasoPattern();

    void printInfo(bool);

    /// prints the mesh details to standard output
    void print();

    /// assigns new coordinates to the nodes
    void setCoordinates(const escript::Data& newX);
    void setElements(ElementFile* elements);
    void setFaceElements(ElementFile* elements);
    void setPoints(ElementFile* elements);

    void prepare(bool optimize);

    /// Initially the element nodes refer to the numbering defined by the
    /// global id assigned to the nodes in the NodeFile. It is also not ensured
    /// that all nodes referred by an element are actually available on the
    /// process. At the output, a local node labeling is used and all nodes are
    /// available. In particular the numbering of the element nodes is between
    /// 0 and Nodes->numNodes.
    /// The function does not create a distribution of the degrees of freedom.
    void resolveNodeIds();

    void createMappings(const std::vector<index_t>& dofDistribution,
                        const std::vector<index_t>& nodeDistribution);

    /// assigns new node reference numbers to all element files.
    /// If k is the old node, the new node is newNode[k-offset].
    void relabelElementNodes(const index_t* newNode, index_t offset);

#ifdef USE_TRILINOS
    /// creates and returns a Trilinos CRS graph suitable to build a sparse
    /// matrix
    esys_trilinos::const_TrilinosGraph_ptr createTrilinosGraph() const;
#endif

private:
    paso::SystemMatrixPattern_ptr makePasoPattern() const;
    void createColoring(const index_t* dofMap);
    void distributeByRankOfDOF(const std::vector<index_t>& distribution);
    void markNodes(std::vector<short>& mask, index_t offset) const;
    void optimizeDOFDistribution(std::vector<index_t>& distribution);
    void optimizeDOFLabeling(const std::vector<index_t>& distribution);
    void optimizeElementOrdering();
    void updateTagList();
    void printElementInfo(const ElementFile* e, const std::string& title,
                          const std::string& defaultType, bool full) const;

    void writeElementInfo(std::ostream& stream, const ElementFile* e,
                          const std::string& defaultType) const;

public:
    /// the name of the mesh
    std::string m_name;
    int approximationOrder;
    int integrationOrder;
    int reducedIntegrationOrder;
    // the table of the nodes
    NodeFile* Nodes;
    // the table of the elements
    ElementFile* Elements;
    // the table of the face elements
    ElementFile* FaceElements;
    // the table of points (treated as elements of dimension 0)
    ElementFile* Points;
    // the tag map mapping names to tag keys
    TagMap tagMap;

    // pointer to the sparse matrix pattern
    paso::SystemMatrixPattern_ptr pasoPattern;
    escript::JMPI MPIInfo;
};

} // namespace dudley

#endif // __DUDLEY_MESH_H__

