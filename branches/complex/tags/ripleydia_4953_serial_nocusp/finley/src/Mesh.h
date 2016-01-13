
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

#ifndef __FINLEY_MESH_H__
#define __FINLEY_MESH_H__

/****************************************************************************

   Finley: Mesh

   A mesh is built from nodes and elements which are describing the
   domain, the surface and point sources (the latter are needed to
   establish links with other codes, in particular to particle
   codes). The nodes are stored in a NodeFile and elements in an
   ElementFile. Four ElementFiles containing the elements
   describe the domain, surface, contact and point sources, respectively.
   Notice that the surface elements do not necessarily cover the entire
   surface of the domain.

   The element type is fixed by the reference element, see
   ReferenceElement.h. The numbering of the nodes starts with 0.

   Important: it is assumed that every node appears in at least
   one element or surface element and that any node used in an
   element, surface element or as a point is specified in the
   NodeFile, see also resolveNodeIds.

   In some cases it is useful to refer to a mesh entirely built from
   order 1 (=linear) elements. The linear version of the mesh can be
   accessed by referring to the first few nodes of each element
   (thanks to the way the nodes are ordered). As the numbering of
   these nodes is not continuous a relabeling vector is introduced
   in the NodeFile. This feature is not fully implemented yet.

   All nodes and elements are tagged. The tag allows to group nodes and
   elements. A typical application is to mark surface elements on a
   certain portion of the domain with the same tag. All these surface
   elements can then be assigned the same value e.g. for the pressure.

   The spatial dimensionality is determined by the type of elements
   used and can be queried using getDim(). Notice that the element type
   also determines the type of surface elements to be used.

*****************************************************************************/

#include "Finley.h"
#include "NodeFile.h"
#include "ElementFile.h"
#include "Util.h"
#include "paso/SystemMatrixPattern.h"

#include <map>
#include <string>

namespace escript {
    class Data;
}

namespace finley {

typedef std::map<std::string, int> TagMap;

/****************************************************************************/

class Mesh
{
public:
    Mesh(const std::string name, int numDim, esysUtils::JMPI& mpi_info);
    ~Mesh();

    static Mesh* load(esysUtils::JMPI& mpi_info, const std::string fname);
    static Mesh* read(esysUtils::JMPI& mpi_info, const std::string fname, int order, int reducedOrder,
                      bool optimize);
    static Mesh* readGmsh(esysUtils::JMPI& mpi_info, const std::string fname, int numDim, int order,
                          int reducedOrder, bool optimize,
                          bool useMacroElements);

    void write(const std::string fname) const;

    int getDim() const { return Nodes->numDim; }
    int getStatus() const { return Nodes->status; }

    void addPoints(int numPoints, const double *points_ptr, const int *tags_ptr);
    void addTagMap(const char* name, int tag_key);
    int getTag(const char* name) const;
    bool isValidTagName(const char* name) const;
    paso::SystemMatrixPattern_ptr getPattern(bool reduce_row_order, bool reduce_col_order);
    paso::SystemMatrixPattern_ptr makePattern(bool reduce_row_order, bool reduce_col_order);
    void printInfo(bool);

    void setCoordinates(const escript::Data& newX);
    void setElements(ElementFile *elements);
    void setFaceElements(ElementFile *elements);
    void setContactElements(ElementFile *elements);
    void setPoints(ElementFile *elements);

    void prepare(bool optimize);
    void resolveNodeIds();
    void createMappings(const std::vector<int>& dofDistribution,
                        const std::vector<int>& nodeDistribution);
    void markDOFsConnectedToRange(int* mask, int offset, int marker,
                                  int firstDOF, int lastDOF, bool useLinear);
    
    void relabelElementNodes(const std::vector<int>&, int offset);

    void glueFaces(double safetyFactor, double tolerance, bool);
    void joinFaces(double safetyFactor, double tolerance, bool);

    void findMatchingFaces(double, double, int*, int*, int*, int*);
    void print();

private:
    void createColoring(const std::vector<int>& dofMap);
    void distributeByRankOfDOF(const std::vector<int>& distribution);
    void markNodes(std::vector<short>& mask, int offset, bool useLinear);
    void optimizeDOFDistribution(std::vector<int>& distribution);
    void optimizeDOFLabeling(const std::vector<int>& distribution);
    void optimizeElementOrdering();
    void setOrders();
    void updateTagList();

public:
    // the name of the mesh
    std::string m_name;
    int approximationOrder;
    int reducedApproximationOrder;
    int integrationOrder;
    int reducedIntegrationOrder;
    // the table of the nodes
    NodeFile* Nodes;
    // the table of the elements
    ElementFile* Elements;
    // the table of the face elements
    ElementFile* FaceElements;
    // the table of the contact elements
    ElementFile* ContactElements;
    // the table of points (treated as elements of dimension 0)
    ElementFile* Points;
    // the tag map mapping names to tag keys
    TagMap tagMap;

    // pointers to the sparse matrix patterns
    paso::SystemMatrixPattern_ptr FullFullPattern;
    paso::SystemMatrixPattern_ptr FullReducedPattern;
    paso::SystemMatrixPattern_ptr ReducedFullPattern;
    paso::SystemMatrixPattern_ptr ReducedReducedPattern;
    esysUtils::JMPI MPIInfo;
};

// this structure is used for matching surface elements
struct FaceCenter
{
   int refId;
   std::vector<double> x;
};


Mesh* Mesh_merge(const std::vector<Mesh*>& meshes);


} // namespace finley

#endif // __FINLEY_MESH_H__

