
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

#ifndef __WEIPA_DATAVAR_H__
#define __WEIPA_DATAVAR_H__

#include <weipa/DomainChunk.h>
#include <ostream>

class DBfile;
//class NcFile;

namespace escript {
    class Data;
}

namespace weipa {

/// \brief A class that provides functionality to read an escript data object
/// from a dump file or an escript::Data instance and write that data in Silo
/// or VTK XML formats.
class DataVar
{
public:
    /// \brief Constructor with variable name
    DataVar(const std::string& name);

    /// \brief Copy constructor. Performs a deep copy of the data values.
    DataVar(const DataVar& d);

    /// \brief Destructor
    ~DataVar();

    /// \brief Initialises values and IDs from an escript::Data instance.
    ///
    /// \return true if a valid instance of expanded data was supplied,
    ///         false if the initialisation was unsuccessful.
    bool initFromEscript(escript::Data& escriptData, const_DomainChunk_ptr dom);

    /// \brief Initialises with integral mesh data like IDs or tags.
    bool initFromMeshData(const_DomainChunk_ptr dom, const IntVec& data,
            int fsCode, Centering c, NodeData_ptr nodes, const IntVec& id);

    /// \brief Reads values and IDs for this variable from an escript dump
    ///        file.
    ///
    /// \return true if the file was found and contains valid escript data
    ///         with at least one sample, false otherwise.
    /// \note Only expanded data with rank <=2 is supported at the moment.
    bool initFromFile(const std::string& filename, const_DomainChunk_ptr dom);

    /// \brief Writes the data into given directory within a Silo file.
    ///
    /// If Silo was not available at compile time this method returns false.
    bool writeToSilo(DBfile* dbfile, const std::string& siloPath,
                     const std::string& units);

    /// \brief Writes the data values to ostream in VTK text format.
    void writeToVTK(std::ostream& os, int ownIndex);

    /// \brief Returns the rank of the data.
    int getRank() const { return rank; }

    /// \brief Returns true if the variable data is node centered, false if
    ///        zone centered.
    bool isNodeCentered() const;

    /// \brief Returns the name of the associated mesh.
    ///
    /// The returned name is one of the sub-meshes of the mesh set with
    /// setMesh() determined on the basis of the function space type and
    /// whether reduced elements are used or not.
    std::string getMeshName() const { return meshName; }

    /// \brief Returns the shape vector of the data.
    ///
    /// The shape vector has as many elements as the rank of this variable.
    const IntVec& getShape() const { return shape; }

    /// \brief Returns the variable name.
    std::string getName() const { return varName; }

    /// \brief Returns the Silo tensor definition for this tensor.
    ///
    /// If the data is tensor data then the components of the tensor are stored
    /// separately in the Silo file. This method then returns a string that
    /// contains the proper Silo expression to put the tensor together again.
    /// For non-tensor data this method returns an empty string.
    std::string getTensorDef() const;

    /// \brief Returns the number of data values.
    int getNumberOfSamples() const { return numSamples; }

    /// \brief Returns the array of data values where array[i] is the i-th
    ///        component of the data.
    const CoordArray& getData() const { return dataArray; }

    /// \brief Returns a flattened array of data values, i.e. the ordering is
    ///        s0c0 s0c1 s0c2 s1c0 s1c1 s1c2 s2c0 ...
    ///        where s denotes the sample number and c the component.
    float* getDataFlat() const;

    /// \brief Returns the total number of components (sum of shape elements).
    int getNumberOfComponents() const;

private:
    void cleanup();

    /// \brief Averages and filters data.
    ///
    /// If a sample consists of more than one data point then this method
    /// averages over the data points and returns the resulting array.
    /// In any case, the data is filtered according to the stride value.
    float* averageData(const float* src, size_t stride);

    /// \brief Prepares a sample ID -> index mapping which is used to reorder
    ///        data.
    IndexMap buildIndexMap();

    /// \brief Reorders the samples according to the corresponding node or
    ///        element IDs.
    ///
    /// \return true if the function space is supported and the number of
    ///         elements or nodes matches the number of data samples.
    bool reorderSamples();

    /// \brief Outputs sample at index to output stream in VTK XML format
    void sampleToStream(std::ostream& os, int index);

    bool initialized;
    const_DomainChunk_ptr domain;
    std::string varName;
    int numSamples, rank, ptsPerSample, funcSpace;
    Centering centering;
    IntVec shape;
    IntVec sampleID;
    CoordArray dataArray;
    std::string meshName, siloMeshName;
};

inline IndexMap DataVar::buildIndexMap()
{
    IndexMap sampleID2idx;
    int idx = sampleID.size()-1;
    // see this thread for why this is done the way it's done:
    // http://www.tech-archive.net/Archive/VC/microsoft.public.vc.stl/2005-01/0075.html
    IntVec::const_reverse_iterator idIt = sampleID.rbegin();
    IntVec::const_reverse_iterator endIt = sampleID.rend();
    for (; idIt != endIt; idIt++, idx--)
        sampleID2idx[*idIt] = idx;

    return sampleID2idx;
}

} // namespace weipa

#endif // __WEIPA_DATAVAR_H__

