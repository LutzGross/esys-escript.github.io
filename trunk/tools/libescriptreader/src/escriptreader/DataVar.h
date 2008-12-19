
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
class MeshWithElements;

//
//
//
class DataVar
{
public:
    /// Constructor with variable name
    DataVar(const std::string& name);

    /// Copy constructor
    DataVar(const DataVar& d);

    /// Special constructor for integral mesh variables
    DataVar(const std::string& name, const IntVec& data,
            MeshWithElements* mesh);

    /// Destructor
    ~DataVar();

    /// Appends data from rhs
    bool append(const DataVar& rhs);

    /// Reads values and IDs for this variable from NetCDF file
    bool readFromNc(const std::string& ncFile);

    /// Associates the data with the given mesh and reorders the samples
    /// according to the corresponding node IDs.
    /// This method must be called before writeToSilo() or getData()
    bool setMesh(MeshWithElements* mesh);

    /// Writes the data into given directory in given Silo file using
    /// provided underlying mesh
    bool writeToSilo(DBfile* dbfile, const std::string& siloPath);

    /// Returns the rank of the data
    int getRank() const { return rank; }

    /// Returns true if the variable is node centered, false if zone centered
    bool isNodeCentered() const;

    /// Returns the name of the associated mesh which is one of the meshes
    /// in mainMesh
    std::string getMeshName(MeshWithElements* mainMesh) const;

    /// Returns the shape vector of the data
    const IntVec& getShape() const { return shape; }

    /// Returns the variable name
    const std::string& getName() const { return varName; }
    std::string getTensorDef() const;

    /// Returns the number of data values
    int getNumberOfSamples() const { return reorderedNumSamples; }

    const CoordArray& getData() const { return reorderedData; }

private:
    /// If a sample consists of more than one data point then this method
    /// averages over the data points and returns the resulting array.
    /// In any case, the data is filtered according to the stride value.
    float* averageData(const float* src, size_t stride);

    void buildIndexMap();

    void reorderSamples(const IndexMap& id2idxMap, const IntVec& requiredIDs);
    void handleGhostZones(const IntVec& reorderArray);

    /// Reads scalar data from NetCDF file
    void readRank0Data(NcFile* ncfile);

    /// Reads vector data from NetCDF file
    void readRank1Data(NcFile* ncfile);

    /// Reads tensor data from NetCDF file
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

#endif // __DATAVAR_H__

