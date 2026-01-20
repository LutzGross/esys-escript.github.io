
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __WEIPA_SPECKLEYELEMENTS_H__
#define __WEIPA_SPECKLEYELEMENTS_H__

#include <weipa/ElementData.h>
#include <weipa/SpeckleyNodes.h>

class DBfile;

namespace speckley {
class SpeckleyDomain;
}


namespace weipa {
 
class SpeckleyElements;
typedef boost::shared_ptr<SpeckleyElements> SpeckleyElements_ptr;

/// \brief Stores and manipulates one type of speckley mesh elements (cells,
///        faces).
///
/// This class provides functionality to manipulate a speckley domain.
/// It is able to load element data from NetCDF files or retrieve it from
/// a SpeckleyDomain instance.
///
/// \note The corresponding mesh nodes are not part of this class but are
///       stored in a NodeData instance.
class SpeckleyElements : public ElementData
{
public:

    /// \brief Constructor with name and accompanying NodeData object.
    SpeckleyElements(const std::string& elementName, SpeckleyNodes_ptr nodes);

    /// \brief Copy constructor
    SpeckleyElements(const SpeckleyElements& e);

    /// \brief Destructor
    virtual ~SpeckleyElements() {}

    /// \brief Initialises with data from a SpeckleyDomain instance.
    bool initFromSpeckley(const speckley::SpeckleyDomain* speckleyDomain, int fsType);

    /// \brief Moves "ghost" elements (whose owner does not match ownIndex)
    ///        and the corresponding data to the end of the arrays.
    void reorderGhostZones(int ownIndex);

    /// \brief Removes "ghost" elements.
    void removeGhostZones(int ownIndex);

    /// \brief Writes connectivity data to a stream in VTK text format.
    virtual void writeConnectivityVTK(std::ostream& os);

    /// \brief Writes element data into given directory in given Silo file.
    ///
    /// Since the mesh depends on element information this method also writes
    /// the node mesh itself. If Silo was not available at compile time or if
    /// a Silo function fails this method returns false.
    bool writeToSilo(DBfile* dbfile, const std::string& siloPath,
                     const StringVec& labels, const StringVec& units,
                     bool writeMeshData);

    /// \brief Returns the names of the meshes associated with the elements.
    virtual StringVec getMeshNames() const;

    /// \brief Returns a vector with the mesh variable names.
    virtual StringVec getVarNames() const;

    /// \brief Returns the number of elements.
    virtual int getNumElements() const { return numElements; }

    /// \brief Returns the number of nodes per element.
    virtual int getNodesPerElement() const { return nodesPerElement; }

    /// \brief Returns the number of "ghost" elements.
    virtual int getGhostCount() const { return numGhostElements; }

    /// \brief Returns the type of the elements.
    virtual ZoneType getType() const { return type; }

    /// \brief Returns a vector of the node IDs used by the elements.
    virtual const IntVec& getNodeList() const { return nodes; }

    /// \brief Returns a vector of element IDs.
    virtual const IntVec& getIDs() const { return ID; }

    /// \brief Returns an array of data values for the name provided.
    ///
    /// The name must be one of the names returned from getVarNames().
    virtual const IntVec& getVarDataByName(const std::string varName) const;

    /// \brief Returns the node mesh instance used by the elements.
    virtual NodeData_ptr getNodes() const { return nodeMesh; }

    /// \brief Returns the reduced elements.
    virtual ElementData_ptr getReducedElements() const { return ElementData_ptr(); }
 
    /// \brief Returns a QuadMaskInfo structure for given functionspace code.
    virtual const QuadMaskInfo& getQuadMask(int functionSpace) const { return quadMask; }
 
    /// \brief If the original element type is not supported they are
    ///        subdivided into N smaller elements (e.g. one Rec9 -> four Rec4)
    ///        and this method returns the multiplication factor N.
    virtual int getElementFactor() const { return 1; }

private:
    SpeckleyElements() {}
    void buildMeshes();
    IntVec prepareGhostIndices(int ownIndex);
    void reorderArray(IntVec& v, const IntVec& idx, int elementsPerIndex);

    SpeckleyNodes_ptr nodeMesh;
    SpeckleyNodes_ptr originalMesh;
    std::string name;
    int numElements;
    int numGhostElements;
    int nodesPerElement;
    ZoneType type;
    IntVec nodes;
    IntVec ID, tag, owner;
    QuadMaskInfo quadMask; // dummy
};

} // namespace weipa

#endif // __WEIPA_SPECKLEYELEMENTS_H__

