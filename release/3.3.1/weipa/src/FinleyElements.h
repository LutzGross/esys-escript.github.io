
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/

#ifndef __WEIPA_FINLEYELEMENTS_H__
#define __WEIPA_FINLEYELEMENTS_H__

#include <weipa/ElementData.h>
#include <weipa/FinleyNodes.h>

extern "C" {
#include <dudley/ElementType.h> // for Dudley_ElementTypeId
#include <finley/ReferenceElements.h> // for Finley_ElementTypeId
}

class DBfile;
class NcFile;

struct Dudley_ElementFile;
struct Finley_ElementFile;

namespace weipa {
 
/// \brief Holds information that is used to convert from finley element types
///        to elements supported by Silo and VTK.
struct FinleyElementInfo
{
    ZoneType elementType, reducedElementType;
    int elementFactor;
    int elementSize, reducedElementSize;
    const size_t* multiCellIndices;
    bool useQuadNodes;
    int quadDim;
};

class FinleyElements;
typedef boost::shared_ptr<FinleyElements> FinleyElements_ptr;

/// \brief Stores and manipulates one type of finley mesh elements (cells,
///        faces or contacts).
///
/// This class provides functionality to manipulate a finley element file.
/// It is able to load element data from NetCDF files or retrieve it from
/// a Finley_ElementFile instance.
///
/// \note The corresponding mesh nodes are not part of this class but are
///       stored in a NodeData instance.
class FinleyElements : public ElementData
{
public:

    /// \brief Constructor with name and accompanying NodeData object.
    FinleyElements(const std::string& elementName, FinleyNodes_ptr nodes);

    /// \brief Copy constructor
    FinleyElements(const FinleyElements& e);

    /// \brief Destructor
    virtual ~FinleyElements() {}

    /// \brief Initialises with data from a Dudley_ElementFile instance.
    bool initFromDudley(const Dudley_ElementFile* dudleyFile);

    /// \brief Initialises with data from a Finley_ElementFile instance.
    bool initFromFinley(const Finley_ElementFile* finleyFile);

    /// \brief Reads element data from escript/finley NetCDF file.
    bool readFromNc(NcFile* ncfile);

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

    /// \brief Returns the original type id of the Finley reference elements.
    Finley_ElementTypeId getFinleyTypeId() const { return finleyTypeId; }

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
    virtual ElementData_ptr getReducedElements() const { return reducedElements; }
 
    /// \brief Returns a QuadMaskInfo structure for given functionspace code.
    virtual const QuadMaskInfo& getQuadMask(int functionSpace) const;
 
    /// \brief If the original element type is not supported they are
    ///        subdivided into N smaller elements (e.g. one Rec9 -> four Rec4)
    ///        and this method returns the multiplication factor N.
    virtual int getElementFactor() const { return elementFactor; }

private:
    FinleyElements() {}
    FinleyElementInfo getDudleyTypeInfo(Dudley_ElementTypeId typeId);
    FinleyElementInfo getFinleyTypeInfo(Finley_ElementTypeId typeId);
    void buildMeshes();
    void buildReducedElements(const FinleyElementInfo& f);
    IntVec prepareGhostIndices(int ownIndex);
    void reorderArray(IntVec& v, const IntVec& idx, int elementsPerIndex);
    QuadMaskInfo buildQuadMask(const CoordArray& quadNodes, int numQNodes);

    FinleyElements_ptr reducedElements;
    FinleyNodes_ptr nodeMesh;
    FinleyNodes_ptr originalMesh;
    std::string name;
    int numElements;
    int numGhostElements;
    int nodesPerElement;
    ZoneType type;
    Finley_ElementTypeId finleyTypeId;
    IntVec nodes;
    IntVec color, ID, tag;
    IntVec owner;
    QuadMaskInfo quadMask, reducedQuadMask;
    int elementFactor;
};

} // namespace weipa

#endif // __WEIPA_FINLEYELEMENTS_H__

