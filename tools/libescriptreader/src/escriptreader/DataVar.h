
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
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

//
// A class that provides functionality to read an escript data object in the
// native NetCDF format and to write that data in Silo format (if available).
//
class DataVar
{
public:
    /// Constructor with variable name
    DataVar(const std::string& name);

    /// Copy constructor. Performs a deep copy of the data values.
    DataVar(const DataVar& d);

    /// Special constructor for integral mesh variables like node IDs
    DataVar(const std::string& name, const IntVec& data,
            MeshWithElements* mesh);

    /// Destructor
    ~DataVar();

    /// Appends raw data and IDs from rhs. Reordered data becomes invalid,
    /// i.e. you have to call setMesh() again before getData().
    /// Returns true if rhs is compatible and appending succeeded, false
    /// otherwise.
    bool append(const DataVar& rhs);

    /// Reads values and IDs for this variable from escript NetCDF file.
    /// Returns true if the file was found and contains valid escript data
    /// with at least one sample, false otherwise.
    /// Note that only expanded data of up to rank 2 is supported at the
    /// moment.
    bool readFromNc(const std::string& filename);

    /// Associates the data with the given mesh and reorders the samples
    /// according to the corresponding node IDs.
    /// Returns true if the function space is supported and the number of
    /// elements or nodes corresponds to the number of data samples.
    /// Note: This method must be called before writeToSilo() or getData().
    bool setMesh(MeshWithElements* mesh);

    /// Writes the data into given directory in given Silo file.
    /// If Silo was not available at compile time or the mesh was not set
    /// beforehand using setMesh() or if a Silo function fails this method
    /// returns false.
    bool writeToSilo(DBfile* dbfile, const std::string& siloPath);

    /// Returns the rank of the data.
    int getRank() const { return rank; }

    /// Returns true if the variable is node centered, false if zone centered.
    bool isNodeCentered() const;

    /// Returns the name of the associated mesh which is one of the meshes
    /// in mainMesh depending on function space type and whether reduced
    /// elements are used or not.
    std::string getMeshName(MeshWithElements* mainMesh) const;

    /// Returns the shape vector of the data. The shape vector has as many
    /// elements as the rank of this variable.
    const IntVec& getShape() const { return shape; }

    /// Returns the variable name.
    const std::string& getName() const { return varName; }

    /// If the data is tensor data then the components of the tensor are stored
    /// separately in the Silo file. This method then returns a string that
    /// contains the proper Silo expression to put the tensor together again.
    /// For non-tensor data this method returns an empty string.
    std::string getTensorDef() const;

    /// Returns the number of data values.
    int getNumberOfSamples() const { return reorderedNumSamples; }

    /// Returns the reordered (not raw!) data values.
    const CoordArray& getData() const { return reorderedData; }

private:
    /// If a sample consists of more than one data point then this method
    /// averages over the data points and returns the resulting array.
    /// In any case, the data is filtered according to the stride value.
    float* averageData(const float* src, size_t stride);

    /// Prepares a sample ID -> index mapping which is used to reorder data.
    void buildIndexMap();

    void reorderSamples(const IndexMap& id2idxMap, const IntVec& requiredIDs);

    void handleGhostZones(const IntVec& reorderArray);

    /// Reads scalar data from NetCDF file.
    void readRank0Data(NcFile* ncfile);

    /// Reads vector data from NetCDF file.
    void readRank1Data(NcFile* ncfile);

    /// Reads tensor data from NetCDF file.
    void readRank2Data(NcFile* ncfile);

    std::string varName;
    int numSamples, rank, ptsPerSample, centering, funcSpace;
    IntVec shape;
    IntVec sampleID;
    IndexMap sampleID2idx;
    CoordArray rawData;       // data as read from NetCDF file
    CoordArray reorderedData; // reordered and filtered data
    int reorderedNumSamples;
    std::string siloMeshName;
    MeshWithElements* fullMesh;
};

//
//
//
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

