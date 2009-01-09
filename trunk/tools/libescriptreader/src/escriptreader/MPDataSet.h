
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
// MPDataSet.h
//
#ifndef __MPDATASET_H__
#define __MPDATASET_H__

#include <escriptreader/common.h>

class DBfile;

namespace EscriptReader {

class DataVar;
class MeshWithElements;

typedef std::vector<MeshWithElements*> MeshBlocks;

struct VarInfo {
    std::string fileName;
    std::string varName;
    DataVar* dataVar;
    bool valid;
};
typedef std::vector<VarInfo> VarVector;


/// \brief Represents an escript/finley dataset including mesh and variables for
/// one timestep. The domain may consist of multiple blocks.
class MPDataSet
{
public:
    /// \brief Constructor.
    MPDataSet();

    /// \brief Destructor.
    ~MPDataSet();

    /// \brief Loads mesh and variables.
    bool load(const std::string meshFile, const StringVec& varFiles,
              const StringVec& varNames, int nBlocks);

    /// \brief Loads only variables using given external mesh.
    bool load(const MeshBlocks& m, const std::string siloFile,
              const StringVec& varFiles, const StringVec& varNames);

    /// \brief Saves the dataset to a Silo file.
    bool saveAsSilo(const std::string siloFile, bool useMultiMesh=true);

    /// \brief Returns the dataset's mesh.
    ///
    /// \note The caller is responsible for freeing the memory of the returned
    /// mesh.
    MeshBlocks extractMesh() { keepMesh = true; return meshBlocks; }

    /// \brief Returns the dataset's mesh.
    const MeshBlocks& getMesh() const { return meshBlocks; }

    /// \brief Returns a vector with the dataset's variables.
    const VarVector& getVariables() const { return variables; }

    /// \brief Returns a vector with the mesh variables.
    const VarVector& getMeshVariables() const { return meshVariables; }

private:
    typedef std::vector<DataVar*> DataParts;

    bool readMeshes();
    bool readVariables();
    void convertMeshVariables();
    void assembleVariable(VarInfo& vi, const DataParts& parts);
    void putSiloMultiMesh(DBfile* dbfile, std::string meshName);
    void putSiloMultiTensor(DBfile* dbfile, const DataVar* var);
    void putSiloMultiVar(DBfile* dbfile, std::string varName,
                         bool useMeshFile = false);

    int numParts;
    bool keepMesh;
    std::string meshFmt;
    std::string siloMeshFile;
    MeshBlocks meshBlocks;
    VarVector variables, meshVariables;
};

} // namespace EscriptReader

#endif // __MPDATASET_H__

