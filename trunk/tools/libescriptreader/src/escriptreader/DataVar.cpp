
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
// DataVar.cpp
//
#include <escriptreader/DataVar.h>
#include <escriptreader/ElementData.h>
#include <escriptreader/MeshWithElements.h>
#include <netcdf.hh>
#if HAVE_SILO
#include <silo.h>
#endif

using namespace std;

namespace EscriptReader {
    
enum {
    NODE_CENTERED = 1,
    ZONE_CENTERED = 2
};

//
// Constructor
//
DataVar::DataVar(const string& name) :
    varName(name), numSamples(0), rank(0), ptsPerSample(0), centering(0),
    reorderedNumSamples(0), fullMesh(NULL)
{
}

//
// Destructor
//
DataVar::~DataVar()
{
    CoordArray::iterator it;
    for (it = reorderedData.begin(); it != reorderedData.end(); it++)
        delete[] *it;
    for (it = rawData.begin(); it != rawData.end(); it++)
        delete[] *it;
}

//
// Copy constructor
//
DataVar::DataVar(const DataVar& d) :
    varName(d.varName), numSamples(d.numSamples),
    rank(d.rank), ptsPerSample(d.ptsPerSample), centering(d.centering),
    funcSpace(d.funcSpace), shape(d.shape), sampleID(d.sampleID),
    reorderedNumSamples(d.reorderedNumSamples), fullMesh(d.fullMesh)
{
    CoordArray::const_iterator it;
    for (it = d.rawData.begin(); it != d.rawData.end(); it++) {
        float* c = new float[numSamples];
        copy(*it, (*it)+numSamples, c);
        rawData.push_back(c);
    }
    for (it = d.reorderedData.begin(); it != d.reorderedData.end(); it++) {
        float* c = new float[reorderedNumSamples];
        copy(*it, (*it)+reorderedNumSamples, c);
        reorderedData.push_back(c);
    }
}

//
// Special constructor for mesh data
//
DataVar::DataVar(const string& name, const IntVec& data,
                 MeshWithElements* mesh) :
    varName(name)
{
    numSamples = data.size();

    float* c = new float[numSamples];
    rawData.push_back(c);
    IntVec::const_iterator it;
    for (it=data.begin(); it != data.end(); it++)
        *c++ = static_cast<float>(*it);

    rank = 0;
    ptsPerSample = 1;
    if (name.compare(0, 6, "Nodes_") == 0) {
        funcSpace = FINLEY_NODES;
        centering = NODE_CENTERED;
        sampleID.insert(sampleID.end(), mesh->getNodeIDs().begin(),
                mesh->getNodeIDs().end());
    } else if (name.compare(0, 9, "Elements_") == 0) {
        funcSpace = FINLEY_ELEMENTS;
        centering = ZONE_CENTERED;
        sampleID.insert(sampleID.end(), mesh->getElements()->getIDs().begin(),
                mesh->getElements()->getIDs().end());
    } else if (name.compare(0, 13, "FaceElements_") == 0) {
        funcSpace = FINLEY_FACE_ELEMENTS;
        centering = ZONE_CENTERED;
        sampleID.insert(sampleID.end(),
                mesh->getFaceElements()->getIDs().begin(),
                mesh->getFaceElements()->getIDs().end());
    } else if (name.compare(0, 16, "ContactElements_") == 0) {
        funcSpace = FINLEY_CONTACT_ELEMENTS_1;
        centering = ZONE_CENTERED;
        sampleID.insert(sampleID.end(),
                mesh->getContactElements()->getIDs().begin(),
                mesh->getContactElements()->getIDs().end());
    } else if (name.compare(0, 7, "Points_") == 0) {
        funcSpace = FINLEY_POINTS;
        centering = NODE_CENTERED;
        sampleID.insert(sampleID.end(), mesh->getPoints()->getIDs().begin(),
                mesh->getPoints()->getIDs().end());
    }

    shape.clear();
    reorderedNumSamples = 0;
}

//
// Appends raw data including IDs from rhs.
//
bool DataVar::append(const DataVar& rhs)
{
    // check if variables are compatible
    if (varName != rhs.varName || ptsPerSample != rhs.ptsPerSample || 
            rank != rhs.rank || shape.size() != rhs.shape.size() ||
            centering != rhs.centering)
        return false;

    for (size_t i=0; i<shape.size(); i++)
        if (shape[i] != rhs.shape[i])
            return false;

    sampleID.insert(sampleID.end(), rhs.sampleID.begin(), rhs.sampleID.end());
    for (size_t i=0; i<rawData.size(); i++) {
        float* c = new float[numSamples+rhs.numSamples];
        copy(rawData[i], rawData[i]+numSamples, c);
        copy(rhs.rawData[i], rhs.rawData[i]+rhs.numSamples, c+numSamples);
        delete[] rawData[i];
        rawData[i] = c;
    }
    numSamples += rhs.numSamples;

    // invalidate previously reordered data
    CoordArray::iterator it;
    for (it = reorderedData.begin(); it != reorderedData.end(); it++)
        delete[] *it;
    reorderedData.clear();
    reorderedNumSamples = 0;
    
    return true;
}

//
// Returns a subset of the src array according to stride parameter.
// If samples consist of multiple values they are averaged beforehand.
// Used to separate (x0,y0,z0,x1,y1,z1,...) into (x0,x1,...), (y0,y1,...) and
// (z0,z1,...)
//
float* DataVar::averageData(const float* src, size_t stride)
{
    float* res = new float[numSamples];

    if (ptsPerSample == 1) {
        float* dest = res;
        for (int i=0; i<numSamples; i++, src+=stride)
            *dest++ = *src;
    } else {
        float* dest = res;
        for (int i=0; i<numSamples; i++) {
            double tmpVal = 0.0;
            for (int j=0; j<ptsPerSample; j++, src+=stride)
                tmpVal += *src;
            *dest++ = (float)(tmpVal / ptsPerSample);
        }
    }
    return res;
}

//
// Reads scalar data (rank=0) from NetCDF file and stores the values
// after averaging.
//
void DataVar::readRank0Data(NcFile* ncfile)
{
    shape.clear();
    float* tempData = new float[ptsPerSample*numSamples];
    NcVar* var = ncfile->get_var("data");
    var->get(tempData, ptsPerSample, numSamples);

    float* c = averageData(tempData, 1);
    rawData.push_back(c);

    delete[] tempData;
}

//
// Reads vector data (rank=1) from NetCDF file and stores the components
// separately after averaging.
//
void DataVar::readRank1Data(NcFile* ncfile)
{
    shape.clear();
    NcDim* dim = ncfile->get_dim("d0");
    shape.push_back(dim->size());

    float* tempData = new float[shape[0]*ptsPerSample*numSamples];
    NcVar* var = ncfile->get_var("data");
    var->get(tempData, shape[0], ptsPerSample, numSamples);

    for (int i=0; i<shape[0]; i++) {
        const float* src = tempData;
        src+=i;
        float* c = averageData(src, shape[0]);
        rawData.push_back(c);
    }
    delete[] tempData;
}

//
// Like readRank1Data() but for tensor data (rank=2).
//
void DataVar::readRank2Data(NcFile* ncfile)
{
    shape.clear();
    NcDim* dim = ncfile->get_dim("d0");
    shape.push_back(dim->size());
    dim = ncfile->get_dim("d1");
    shape.push_back(dim->size());

    float* tempData = new float[shape[0]*shape[1]*ptsPerSample*numSamples];
    NcVar* var = ncfile->get_var("data");
    var->get(tempData, shape[0], shape[1], ptsPerSample, numSamples);

    for (int i=0; i<shape[1]; i++) {
        for (int j=0; j<shape[0]; j++) {
            const float* src = tempData;
            src+=i*shape[0]+j;
            float* c = averageData(src, shape[0]*shape[1]);
            rawData.push_back(c);
        }
    }
    delete[] tempData;
}

//
// Reads a NetCDF file in escript/finley format.
//
bool DataVar::readFromNc(const string& filename)
{
    NcError ncerr(NcError::silent_nonfatal);    
    NcFile* input = new NcFile(filename.c_str());
    if (!input->is_valid()) {
        cerr << "Could not open input file " << filename << "." << endl;
        delete input;
        return false;
    }

    NcDim* dim;
    NcAtt* att;

    dim = input->get_dim("num_samples");
    numSamples = dim->size();

    // if there are no data samples bail out
    if (numSamples == 0) {
        delete input;
        return false;
    }

    att = input->get_att("type_id");
    int typeID = att->as_int(0);
    if (typeID != 2) {
        cerr << "WARNING: Only expanded data supported at the moment!" << endl;
        delete input;
        return false;
    }

    att = input->get_att("rank");
    rank = att->as_int(0);

    dim = input->get_dim("num_data_points_per_sample");
    ptsPerSample = dim->size();

    att = input->get_att("function_space_type");
    funcSpace = att->as_int(0);

#ifdef _DEBUG
    cout << varName << ":\t" << numSamples << " samples,  "
        << ptsPerSample << " pts/s,  rank: " << rank << endl;
#endif

    sampleID.clear();
    sampleID.insert(sampleID.end(), numSamples, 0);
    NcVar* var = input->get_var("id");
    var->get(&sampleID[0], numSamples);

    switch (rank) {
        case 0:
            readRank0Data(input);
            break;
        case 1:
            readRank1Data(input);
            break;
        case 2:
            readRank2Data(input);
            break;
        default:
            cerr << "WARNING: Rank " << rank << " data is not supported!\n";
            delete input;
            return false;
    }

    delete input;
    return true;
}

//
// Returns one of the mesh names provided by mainMesh that matches the
// data variable's function space type and reduced/unreduced state.
//
string DataVar::getMeshName(MeshWithElements* mainMesh) const
{
    string name;

    switch (funcSpace) {
        case FINLEY_REDUCED_NODES:
        case FINLEY_REDUCED_DEGREES_OF_FREEDOM:
        case FINLEY_REDUCED_ELEMENTS:
        case FINLEY_ELEMENTS:
            if (mainMesh->getElements()->reducedCount > 0) {
                name = mainMesh->getElements()->reducedMesh->getName();
            } else {
                name = mainMesh->getElements()->fullMesh->getName();
            }
            break;

        case FINLEY_REDUCED_FACE_ELEMENTS:
        case FINLEY_FACE_ELEMENTS:
            if (mainMesh->getFaceElements()->reducedCount > 0) {
                name = mainMesh->getFaceElements()->reducedMesh->getName();
            } else {
                name = mainMesh->getFaceElements()->fullMesh->getName();
            }
            break;

        case FINLEY_REDUCED_CONTACT_ELEMENTS_1:
        case FINLEY_REDUCED_CONTACT_ELEMENTS_2:
        case FINLEY_CONTACT_ELEMENTS_1:
        case FINLEY_CONTACT_ELEMENTS_2:
            if (mainMesh->getContactElements()->reducedCount > 0) {
                name = mainMesh->getContactElements()->reducedMesh->getName();
            } else {
                name = mainMesh->getContactElements()->fullMesh->getName();
            }
            break;

        case FINLEY_NODES:
        case FINLEY_DEGREES_OF_FREEDOM:
            name = mainMesh->getElements()->fullMesh->getName();
            break;

        case FINLEY_POINTS:
            name = mainMesh->getPoints()->fullMesh->getName();
            break;
    }
    return name;
}

//
// Returns true if the data values are nodal, false if they are zonal.
//
bool DataVar::isNodeCentered() const
{
    return (funcSpace == FINLEY_REDUCED_NODES ||
            funcSpace == FINLEY_REDUCED_DEGREES_OF_FREEDOM ||
            funcSpace == FINLEY_NODES ||
            funcSpace == FINLEY_DEGREES_OF_FREEDOM ||
            funcSpace == FINLEY_POINTS);
}

//
// Filters and reorders the raw sample values according to the IDs provided
// in 'requiredIDs'. This is used to have data arrays ordered according to
// the underlying mesh (i.e. DataID[i]==MeshNodeID[i])
//
void DataVar::reorderSamples(const IndexMap& id2idxMap,
                             const IntVec& requiredIDs)
{
    CoordArray::iterator it;
    for (it = reorderedData.begin(); it != reorderedData.end(); it++)
        delete[] *it;
    reorderedData.clear();

    buildIndexMap();
    for (size_t i=0; i < rawData.size(); i++) {
        float* c = new float[reorderedNumSamples];
        reorderedData.push_back(c);
        const float* src = rawData[i];
        IntVec::const_iterator idIt = requiredIDs.begin();
        for (; idIt != requiredIDs.end(); idIt++) {
            size_t srcIdx = sampleID2idx.find(*idIt)->second;
            size_t destIdx = id2idxMap.find(*idIt)->second;
            c[destIdx] = src[srcIdx];
        }
    }
}

//
// For zonal data this method reorders the values according to the indices
// given in reorderArray. This is used to move ghost zones to the end of
// the arrays which conforms to Silo's expected format.
// Nodal data is not changed by this method.
//
void DataVar::handleGhostZones(const IntVec& reorderArray)
{
    if (centering == NODE_CENTERED)
        return;

    vector<float*>::iterator it;
    for (it = reorderedData.begin(); it!=reorderedData.end(); it++) {
        float* original = *it;
        float* reordered = new float[reorderedNumSamples];
        float* arrIt = reordered;
        IntVec::const_iterator idxIt;
        for (idxIt=reorderArray.begin(); idxIt!=reorderArray.end(); idxIt++)
            *arrIt++ = original[*idxIt];

        delete[] *it;
        *it = reordered;
    }
}

//
// Makes the top-level mesh known to this data variable. The mesh is used
// to reorder and filter the samples and inquire the number of ghost zones.
//
bool DataVar::setMesh(MeshWithElements* mesh)
{
    if (fullMesh == mesh)
        return true;

    const IndexMap* id2idxMap;
    const IntVec* reqIDs;
    const IntVec* reorderArray = NULL;

    switch (funcSpace) {
        case FINLEY_REDUCED_NODES:
        case FINLEY_REDUCED_DEGREES_OF_FREEDOM:
            {
                centering = NODE_CENTERED;
                ElementData* cells = mesh->getElements();
                if (cells->reducedCount > 0) {
                    if (cells->getReducedGhostCount())
                        reorderArray = &cells->reducedIndexArray;
                    siloMeshName = cells->reducedMesh->getFullSiloName();
                    id2idxMap = &cells->reducedMesh->getIndexMap();
                    reqIDs = &cells->reducedMesh->getNodeIDs();
                    reorderedNumSamples = cells->reducedMesh->getNumNodes();
                } else {
                    if (cells->getGhostCount())
                        reorderArray = &cells->indexArray;
                    siloMeshName = cells->fullMesh->getFullSiloName();
                    id2idxMap = &cells->fullMesh->getIndexMap();
                    reqIDs = &cells->fullMesh->getNodeIDs();
                    reorderedNumSamples = cells->fullMesh->getNumNodes();
                }
            }
            break;

        case FINLEY_NODES:
        case FINLEY_DEGREES_OF_FREEDOM:
            {
                centering = NODE_CENTERED;
                ElementData* cells = mesh->getElements();
                if (cells->getGhostCount())
                    reorderArray = &cells->indexArray;
                siloMeshName = cells->fullMesh->getFullSiloName();
                id2idxMap = &cells->fullMesh->getIndexMap();
                reqIDs = &cells->fullMesh->getNodeIDs();
                reorderedNumSamples = cells->fullMesh->getNumNodes();
            }
            break;

        case FINLEY_REDUCED_ELEMENTS:
        case FINLEY_ELEMENTS:
            {
                centering = ZONE_CENTERED;
                ElementData* cells = mesh->getElements();
                id2idxMap = &cells->ID2idx;
                reqIDs = &cells->getIDs();
                if (cells->reducedCount > 0) {
                    if (cells->getReducedGhostCount())
                        reorderArray = &cells->reducedIndexArray;
                    reorderedNumSamples = cells->reducedCount;
                    siloMeshName = cells->reducedMesh->getFullSiloName();
                } else {
                    if (cells->getGhostCount())
                        reorderArray = &cells->indexArray;
                    reorderedNumSamples = cells->count;
                    siloMeshName = cells->fullMesh->getFullSiloName();
                }
            }
            break;

        case FINLEY_REDUCED_FACE_ELEMENTS:
        case FINLEY_FACE_ELEMENTS:
            {
                centering = ZONE_CENTERED;
                ElementData* cells = mesh->getFaceElements();
                id2idxMap = &cells->ID2idx;
                reqIDs = &cells->getIDs();
                if (cells->reducedCount > 0) {
                    if (cells->getReducedGhostCount())
                        reorderArray = &cells->reducedIndexArray;
                    reorderedNumSamples = cells->reducedCount;
                    siloMeshName = cells->reducedMesh->getFullSiloName();
                } else {
                    if (cells->getGhostCount())
                        reorderArray = &cells->indexArray;
                    reorderedNumSamples = cells->count;
                    siloMeshName = cells->fullMesh->getFullSiloName();
                }
            }
            break;

        case FINLEY_REDUCED_CONTACT_ELEMENTS_1:
        case FINLEY_REDUCED_CONTACT_ELEMENTS_2:
        case FINLEY_CONTACT_ELEMENTS_1:
        case FINLEY_CONTACT_ELEMENTS_2:
            {
                centering = ZONE_CENTERED;
                ElementData* cells = mesh->getContactElements();
                id2idxMap = &cells->ID2idx;
                reqIDs = &cells->getIDs();
                if (cells->reducedCount > 0) {
                    if (cells->getReducedGhostCount())
                        reorderArray = &cells->reducedIndexArray;
                    reorderedNumSamples = cells->reducedCount;
                    siloMeshName = cells->reducedMesh->getFullSiloName();
                } else {
                    if (cells->getGhostCount())
                        reorderArray = &cells->indexArray;
                    reorderedNumSamples = cells->count;
                    siloMeshName = cells->fullMesh->getFullSiloName();
                }
            }
            break;

        case FINLEY_POINTS:
            {
                centering = NODE_CENTERED;
                ElementData* cells = mesh->getPoints();
                if (cells->getGhostCount())
                    reorderArray = &cells->indexArray;
                siloMeshName = cells->fullMesh->getFullSiloName();
                id2idxMap = &cells->ID2idx;
                reqIDs = &cells->getIDs();
                reorderedNumSamples = cells->count;
            }
            break;

        default:
            cerr << "Unknown function space type " << funcSpace << "!\n";
            return false;
    }

    if (reorderedNumSamples > numSamples) {
        cerr << "WARNING: " << varName << " has " << numSamples
            << " instead of " << reorderedNumSamples << " samples!" << endl;
        return false;
    }

    fullMesh = mesh;

    reorderSamples(*id2idxMap, *reqIDs);
    if (reorderArray)
        handleGhostZones(*reorderArray);
    return true;
}

///////////////////////////////
// SILO related methods follow
///////////////////////////////

//
// If the data is tensor data then the components of the tensor are stored
// separately in the Silo file. This method then returns a string that
// contains the proper Silo expression to put the tensor together again.
// For non-tensor data this method returns an empty string.
//
string DataVar::getTensorDef() const
{
    if (rank < 2)
        return string();
    
    /// Format string for Silo 2x2 tensor
    const string tensor2DefFmt =
        "{{ <%sa_00>, <%sa_01> },"
        " { <%sa_10>, <%sa_11> }}";

    /// Format string for Silo 3x3 tensor
    const string tensor3DefFmt =
        "{{ <%sa_00>, <%sa_01>, <%sa_02> },"
        " { <%sa_10>, <%sa_11>, <%sa_12> },"
        " { <%sa_20>, <%sa_21>, <%sa_22> }}";

    string tensorDef;
    string tensorDir = varName+string("_comps/");
    if (shape[1] == 3) {
        char* tDef = new char[tensor3DefFmt.length()+9*tensorDir.length()];
        sprintf(tDef, tensor3DefFmt.c_str(),
                tensorDir.c_str(), tensorDir.c_str(), tensorDir.c_str(),
                tensorDir.c_str(), tensorDir.c_str(), tensorDir.c_str(),
                tensorDir.c_str(), tensorDir.c_str(), tensorDir.c_str());
        tensorDef = tDef;
        delete[] tDef;
    } else {
        char* tDef = new char[tensor2DefFmt.length()+4*tensorDir.length()];
        sprintf(tDef, tensor2DefFmt.c_str(),
                tensorDir.c_str(), tensorDir.c_str(),
                tensorDir.c_str(), tensorDir.c_str(),
                tensorDir.c_str(), tensorDir.c_str());
        tensorDef = tDef;
        delete[] tDef;
    }
    return tensorDef;
}

//
// Writes the data to given Silo file under the virtual path provided.
// The corresponding mesh must have been written already and made known
// to this variable by a call to setMesh().
//
bool DataVar::writeToSilo(DBfile* dbfile, const string& siloPath)
{
#if HAVE_SILO
    if (numSamples == 0)
        return true;

    // have to set mesh first
    if (!fullMesh)
        return false;

    int ret;

    if (siloPath != "") {
        ret = DBSetDir(dbfile, siloPath.c_str());
        if (ret != 0)
            return false;
    }
    
    char* meshName = (char*)siloMeshName.c_str();
    int dcenter = (centering == NODE_CENTERED ? DB_NODECENT : DB_ZONECENT);

    if (rank == 0) {
        ret = DBPutUcdvar1(dbfile, varName.c_str(), meshName, reorderedData[0],
                reorderedNumSamples, NULL, 0, DB_FLOAT, dcenter, NULL);
    }
    else if (rank == 1) {
        const string comps[3] = {
            varName+string("_x"), varName+string("_y"), varName+string("_z")
        };
        const char* varnames[3] = {
            comps[0].c_str(), comps[1].c_str(), comps[2].c_str()
        };

        ret = DBPutUcdvar(dbfile, varName.c_str(), meshName, shape[0],
                (char**)varnames, &reorderedData[0], reorderedNumSamples, NULL,
                0, DB_FLOAT, dcenter, NULL);
    }
    else {
        string tensorDir = varName+string("_comps/");
        ret = DBMkdir(dbfile, tensorDir.c_str());
        if (ret == 0) {
            int one = 1;
            DBoptlist* optList = DBMakeOptlist(1);
            DBAddOption(optList, DBOPT_HIDE_FROM_GUI, &one);

            for (int i=0; i<shape[1]; i++) {
                for (int j=0; j<shape[0]; j++) {
                    ostringstream varname;
                    varname << tensorDir << "a_" << i << j;
                    ret = DBPutUcdvar1(dbfile, varname.str().c_str(), meshName,
                            reorderedData[i*shape[0]+j], reorderedNumSamples,
                            NULL, 0, DB_FLOAT, dcenter, optList);
                    if (ret != 0) break;
                }
                if (ret != 0) break;
            }
            DBFreeOptlist(optList);
        } // ret==0
    } // rank

    DBSetDir(dbfile, "/");
    return (ret == 0);

#else // !HAVE_SILO
    return false;
#endif
}

} // namespace EscriptReader

