
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

//
// DataVar.h
//
#ifndef __DATAVAR_H__
#define __DATAVAR_H__

#include <escriptreader/common.h>

class DBfile;
class NcFile;

namespace EscriptReader {

class MeshWithElements;

/// \brief A class that provides functionality to read an escript data object in
/// the native NetCDF format and to write that data in Silo format (if available).
class DataVar
{
public:
    /// \brief Constructor with variable name
    DataVar(const std::string& name);

    /// \brief Copy constructor. Performs a deep copy of the data values.
    DataVar(const DataVar& d);

    /// \brief Special constructor for integral mesh variables like node IDs
    DataVar(const std::string& name, const IntVec& data,
            MeshWithElements* mesh);

    /// \brief Destructor
    ~DataVar();

    /// \brief Appends raw data and IDs from rhs.
    ///
    /// Reordered data becomes invalid, i.e. you have to call setMesh() again
    /// before getData().
    /// \return true if rhs is compatible and appending succeeded, false
    /// otherwise.
    bool append(const DataVar& rhs);

    /// \brief Reads values and IDs for this variable from escript NetCDF file.
    ///
    /// \return true if the file was found and contains valid escript data with
    /// at least one sample, false otherwise.
    /// \note Only expanded data of up to rank 2 is supported at the
    /// moment.
    bool readFromNc(const std::string& filename);

    /// \brief Associates the data with the given mesh and reorders the samples
    /// according to the corresponding node IDs.
    ///
    /// \return true if the function space is supported and the number of
    /// elements or nodes corresponds to the number of data samples.
    /// \note This method must be called before writeToSilo() or getData().
    bool setMesh(MeshWithElements* mesh);

    /// \brief Writes the data into given directory in given Silo file.
    ///
    /// If Silo was not available at compile time or the mesh was not set
    /// beforehand using setMesh() or if a Silo function fails this method
    /// returns false.
    bool writeToSilo(DBfile* dbfile, const std::string& siloPath);

    /// \brief Returns the rank of the data.
    int getRank() const { return rank; }

    /// \brief Returns true if the variable is node centered, false if zone
    /// centered.
    bool isNodeCentered() const;

    /// \brief Returns the name of the associated mesh.
    ///
    /// The returned name is one of the meshes in mainMesh depending on function
    /// space type and whether reduced elements are used or not.
    std::string getMeshName(MeshWithElements* mainMesh) const;

    /// \brief Returns the shape vector of the data.
    ///
    /// The shape vector has as many elements as the rank of this variable.
    const IntVec& getShape() const { return shape; }

    /// \brief Returns the variable name.
    const std::string& getName() const { return varName; }

    /// \brief Returns the Silo tensor definition for this tensor.
    ///
    /// If the data is tensor data then the components of the tensor are stored
    /// separately in the Silo file. This method then returns a string that
    /// contains the proper Silo expression to put the tensor together again.
    /// For non-tensor data this method returns an empty string.
    std::string getTensorDef() const;

    /// \brief Returns the number of data values.
    int getNumberOfSamples() const { return reorderedNumSamples; }

    /// \brief Returns the reordered (not raw!) data values.
    const CoordArray& getData() const { return reorderedData; }

private:
    /// \brief Averages and filters data.
    ///
    /// If a sample consists of more than one data point then this method
    /// averages over the data points and returns the resulting array.
    /// In any case, the data is filtered according to the stride value.
    float* averageData(const float* src, size_t stride);

    /// \brief Prepares a sample ID -> index mapping which is used to reorder
    /// data.
    void buildIndexMap();

    void reorderSamples(const IndexMap& id2idxMap, const IntVec& requiredIDs);

    void handleGhostZones(const IntVec& reorderArray);

    /// \brief Reads scalar data from NetCDF file.
    void readRank0Data(NcFile* ncfile);

    /// \brief Reads vector data from NetCDF file.
    void readRank1Data(NcFile* ncfile);

    /// \brief Reads tensor data from NetCDF file.
    void readRank2Data(NcFile* ncfile);

    std::string varName;
    int numSamples, rank, ptsPerSample, centering, funcSpace;
    IntVec shape;
    IntVec sampleID;
    IndexMap sampleID2idx;
    CoordArray rawData;       /// data as read from NetCDF file
    CoordArray reorderedData; /// reordered and filtered data
    int reorderedNumSamples;
    std::string siloMeshName;
    MeshWithElements* fullMesh;
};

inline void DataVar::buildIndexMap()
{
    sampleID2idx.clear();
    int idx = 0;
    IntVec::const_iterator idIt;
    for (idIt = sampleID.begin(); idIt != sampleID.end(); idIt++, idx++)
        sampleID2idx[*idIt] = idx;
}

} // namespace EscriptReader

#endif // __DATAVAR_H__

