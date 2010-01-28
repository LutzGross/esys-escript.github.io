
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

#ifndef __DATAVAR_H__
#define __DATAVAR_H__

#include <escriptexport/escriptexport.h>

class DBfile;
class NcFile;

namespace escript {
    class Data;
}

namespace escriptexport {

class FinleyMesh;

/// \brief A class that provides functionality to read an escript data object
/// from NetCDF file or an escript::Data instance and write that data in Silo
/// or VTK XML formats.
class DataVar
{
public:
    /// \brief Constructor with variable name
    ESCRIPTEXPORT_DLL_API
    DataVar(const std::string& name);

    /// \brief Copy constructor. Performs a deep copy of the data values.
    ESCRIPTEXPORT_DLL_API
    DataVar(const DataVar& d);

    /// \brief Destructor
    ESCRIPTEXPORT_DLL_API
    ~DataVar();

    /// \brief Initialises values and IDs from an escript::Data instance.
    ///
    /// \return true if a valid instance of expanded data was supplied,
    ///         false if the initialisation was unsuccessful.
    ESCRIPTEXPORT_DLL_API
    bool initFromEscript(escript::Data& escriptData, FinleyMesh_ptr mesh);

    /// \brief Initialises with integral mesh data like IDs or tags.
    ///
    /// The data is retrieved from the mesh using the variable name.
    ESCRIPTEXPORT_DLL_API
    bool initFromMesh(FinleyMesh_ptr mesh);

    /// \brief Reads values and IDs for this variable from escript NetCDF file.
    ///
    /// \return true if the file was found and contains valid escript data
    ///         with at least one sample, false otherwise.
    /// \note Only expanded data with rank <=2 is supported at the moment.
    ESCRIPTEXPORT_DLL_API
    bool initFromNetCDF(const std::string& filename, FinleyMesh_ptr mesh);

    /// \brief Writes the data into given directory in given Silo file.
    ///
    /// If Silo was not available at compile time or the mesh was not set
    /// beforehand using setMesh() or if a Silo function fails this method
    /// returns false.
    ESCRIPTEXPORT_DLL_API
    bool writeToSilo(DBfile* dbfile, const std::string& siloPath);

    /// \brief Returns the rank of the data.
    ESCRIPTEXPORT_DLL_API
    int getRank() const { return rank; }

    /// \brief Returns true if the variable data is node centered, false if
    ///        zone centered.
    ESCRIPTEXPORT_DLL_API
    bool isNodeCentered() const;

    /// \brief Returns the name of the associated mesh.
    ///
    /// The returned name is one of the sub-meshes of the mesh set with
    /// setMesh() determined on the basis of the function space type and
    /// whether reduced elements are used or not.
    ESCRIPTEXPORT_DLL_API
    std::string getMeshName() const { return meshName; }

    /// \brief Returns the shape vector of the data.
    ///
    /// The shape vector has as many elements as the rank of this variable.
    ESCRIPTEXPORT_DLL_API
    const IntVec& getShape() const { return shape; }

    /// \brief Returns the variable name.
    ESCRIPTEXPORT_DLL_API
    std::string getName() const { return varName; }

    /// \brief Returns the Silo tensor definition for this tensor.
    ///
    /// If the data is tensor data then the components of the tensor are stored
    /// separately in the Silo file. This method then returns a string that
    /// contains the proper Silo expression to put the tensor together again.
    /// For non-tensor data this method returns an empty string.
    ESCRIPTEXPORT_DLL_API
    std::string getTensorDef() const;

    /// \brief Returns the number of data values.
    ESCRIPTEXPORT_DLL_API
    int getNumberOfSamples() const { return numSamples; }

    /// \brief Returns the array of data values where array[i] is the i-th
    ///        component of the data.
    ESCRIPTEXPORT_DLL_API
    const CoordArray& getData() const { return dataArray; }

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

    bool initialized;
    FinleyMesh_ptr finleyMesh;
    std::string varName;
    int numSamples, rank, ptsPerSample, centering, funcSpace;
    IntVec shape;
    IntVec sampleID;
    CoordArray dataArray;
    std::string meshName, siloMeshName;
};

inline IndexMap DataVar::buildIndexMap()
{
    IndexMap sampleID2idx;
    int idx = sampleID.size()-1;
    IntVec::const_reverse_iterator idIt;
    for (idIt = sampleID.rbegin(); idIt != sampleID.rend(); idIt++, idx--)
        sampleID2idx[*idIt] = idx;

    return sampleID2idx;
}

} // namespace escriptexport

#endif // __DATAVAR_H__

