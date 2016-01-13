
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

#include <escriptexport/DataVar.h>
#include <escriptexport/ElementData.h>
#include <escriptexport/FinleyMesh.h>
#include <escriptexport/NodeData.h>
#include <escript/Data.h>

#if USE_NETCDF
#include <netcdf.hh>
#endif

#if USE_SILO
#include <silo.h>
#endif

using namespace std;

namespace escriptexport {
    
enum {
    NODE_CENTERED = 1,
    ZONE_CENTERED = 2
};

//
// Constructor
//
DataVar::DataVar(const string& name) :
    initialized(false), varName(name),
    numSamples(0), rank(0), ptsPerSample(0), centering(0)
{
}

//
// Copy constructor
//
DataVar::DataVar(const DataVar& d) :
    varName(d.varName), numSamples(d.numSamples),
    rank(d.rank), ptsPerSample(d.ptsPerSample), centering(d.centering),
    funcSpace(d.funcSpace), shape(d.shape), sampleID(d.sampleID)
{
    if (numSamples > 0) {
        CoordArray::const_iterator it;
        for (it = d.dataArray.begin(); it != d.dataArray.end(); it++) {
            float* c = new float[numSamples];
            copy(*it, (*it)+numSamples, c);
            dataArray.push_back(c);
        }
    }
    initialized = d.initialized;
}

//
// Destructor
//
DataVar::~DataVar()
{
    cleanup();
}

//
//
//
void DataVar::cleanup()
{
    CoordArray::iterator it;
    for (it = dataArray.begin(); it != dataArray.end(); it++)
        delete[] *it;
    dataArray.clear();
    shape.clear();
    sampleID.clear();
    numSamples = 0;
    initialized = false;
}

//
//
//
bool DataVar::initFromEscript(escript::Data& escriptData, FinleyMesh_ptr mesh)
{
    cleanup();

    if (!escriptData.actsExpanded()) {
        cerr << "WARNING: Only expanded data supported!" << endl;
        return false;
    }

    rank = escriptData.getDataPointRank();
    ptsPerSample = escriptData.getNumDataPointsPerSample();
    shape = escriptData.getDataPointShape();
    funcSpace = escriptData.getFunctionSpace().getTypeCode();
    numSamples = escriptData.getNumSamples();

    if (funcSpace == FINLEY_REDUCED_NODES || funcSpace == FINLEY_NODES) {
        centering = NODE_CENTERED;
    } else {
        centering = ZONE_CENTERED;
    }

#ifdef _DEBUG
    cout << varName << ":\t" << numSamples << " samples,  "
        << ptsPerSample << " pts/s,  rank: " << rank << endl;
#endif

    initialized = true;

    if (numSamples == 0)
        return true;

    const int* iPtr = escriptData.getFunctionSpace().borrowSampleReferenceIDs();
    sampleID.insert(sampleID.end(), numSamples, 0);
    copy(iPtr, iPtr+numSamples, sampleID.begin());

    size_t dimSize = 1;
    if (rank > 0)
        dimSize *= shape[0];
    if (rank > 1)
        dimSize *= shape[1];
    if (rank > 2) {
        cerr << "WARNING: Rank " << rank << " data is not supported!\n";
        initialized = false;
    }

    if (initialized) {
        size_t dataSize = dimSize * ptsPerSample;
        float* tempData = new float[dataSize*numSamples];
        float* destPtr = tempData;
        for (int sampleNo=0; sampleNo<numSamples; sampleNo++) {
            const escript::DataAbstract::ValueType::value_type* values =
                escriptData.getSampleDataRO(sampleNo);
            copy(values, values+dataSize, destPtr);
            destPtr += dataSize;
        }

        const float* srcPtr = tempData;
        for (int i=0; i < dimSize; i++, srcPtr++) {
            float* c = averageData(srcPtr, dimSize);
            dataArray.push_back(c);
        }
        delete[] tempData;

        initialized = filterSamples(mesh);
    }

    return initialized;
}

//
// Initialise with mesh data
//
bool DataVar::initFromMesh(FinleyMesh_ptr mesh)
{
    cleanup();
    
    const IntVec& data = mesh->getVarDataByName(varName);
    rank = 0;
    ptsPerSample = 1;
    numSamples = data.size();

    if (numSamples > 0) {
        float* c = new float[numSamples];
        dataArray.push_back(c);
        IntVec::const_iterator it;
        for (it=data.begin(); it != data.end(); it++)
            *c++ = static_cast<float>(*it);

        if (varName.compare(0, 6, "Nodes_") == 0) {
            funcSpace = FINLEY_NODES;
            centering = NODE_CENTERED;
            sampleID.insert(sampleID.end(),
                    mesh->getNodes()->getNodeIDs().begin(),
                    mesh->getNodes()->getNodeIDs().end());
        } else if (varName.compare(0, 9, "Elements_") == 0) {
            funcSpace = FINLEY_ELEMENTS;
            centering = ZONE_CENTERED;
            sampleID.insert(sampleID.end(),
                    mesh->getElements()->getIDs().begin(),
                    mesh->getElements()->getIDs().end());
        } else if (varName.compare(0, 13, "FaceElements_") == 0) {
            funcSpace = FINLEY_FACE_ELEMENTS;
            centering = ZONE_CENTERED;
            sampleID.insert(sampleID.end(),
                    mesh->getFaceElements()->getIDs().begin(),
                    mesh->getFaceElements()->getIDs().end());
        } else if (varName.compare(0, 16, "ContactElements_") == 0) {
            funcSpace = FINLEY_CONTACT_ELEMENTS_1;
            centering = ZONE_CENTERED;
            sampleID.insert(sampleID.end(),
                    mesh->getContactElements()->getIDs().begin(),
                    mesh->getContactElements()->getIDs().end());
        } else {
            return false;
        }

        NodeData_ptr nodes = mesh->getMeshForFinleyFS(funcSpace);
        meshName = nodes->getName();
        siloMeshName = nodes->getFullSiloName();
    }
    initialized = true;

    return initialized;
}

//
// Reads variable data from NetCDF file
//
bool DataVar::initFromNetCDF(const string& filename, FinleyMesh_ptr mesh)
{
    cleanup();
    
#if USE_NETCDF
    NcError ncerr(NcError::silent_nonfatal);    
    NcFile* input = new NcFile(filename.c_str());
    if (!input->is_valid()) {
        cerr << "Could not open input file " << filename << "." << endl;
        delete input;
        return false;
    }

    NcDim* dim;
    NcAtt* att;

    att = input->get_att("type_id");
    int typeID = att->as_int(0);
    if (typeID != 2) {
        cerr << "WARNING: Only expanded data supported!" << endl;
        delete input;
        return false;
    }

    att = input->get_att("rank");
    rank = att->as_int(0);

    dim = input->get_dim("num_data_points_per_sample");
    ptsPerSample = dim->size();

    att = input->get_att("function_space_type");
    funcSpace = att->as_int(0);

    if (funcSpace == FINLEY_REDUCED_NODES || funcSpace == FINLEY_NODES) {
        centering = NODE_CENTERED;
    } else {
        centering = ZONE_CENTERED;
    }

    dim = input->get_dim("num_samples");
    numSamples = dim->size();

#ifdef _DEBUG
    cout << varName << ":\t" << numSamples << " samples,  "
        << ptsPerSample << " pts/s,  rank: " << rank << endl;
#endif

    initialized = true;

    // if there are no data samples we're done
    if (numSamples == 0) {
        delete input;
        return true;
    }

    sampleID.insert(sampleID.end(), numSamples, 0);
    NcVar* var = input->get_var("id");
    var->get(&sampleID[0], numSamples);

    size_t dimSize = 1;
    vector<long> counts;

    if (rank > 0) {
        dim = input->get_dim("d0");
        int d = dim->size();
        shape.push_back(d);
        counts.push_back(d);
        dimSize *= d;
    }
    if (rank > 1) {
        dim = input->get_dim("d1");
        int d = dim->size();
        shape.push_back(d);
        counts.push_back(d);
        dimSize *= d;
    }
    if (rank > 2) {
        cerr << "WARNING: Rank " << rank << " data is not supported!\n";
        initialized = false;
    }
 
    if (initialized) {
        size_t dataSize = dimSize*numSamples*ptsPerSample;
        counts.push_back(ptsPerSample);
        counts.push_back(numSamples);
        float* tempData = new float[dataSize];
        NcVar* var = input->get_var("data");
        var->get(tempData, &counts[0]);

        const float* srcPtr = tempData;
        for (int i=0; i < dimSize; i++, srcPtr++) {
            float* c = averageData(srcPtr, dimSize);
            dataArray.push_back(c);
        }
        delete[] tempData;

        initialized = filterSamples(mesh);
    }

    delete input;
#endif // USE_NETCDF

    return initialized;
}

//
// Returns true if the data values are nodal, false if they are zonal.
//
bool DataVar::isNodeCentered() const
{
    return (centering == NODE_CENTERED);
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
// Filters and reorders the raw sample values according to the IDs provided
// in 'requiredIDs'. This is used to have data arrays ordered according to
// the underlying mesh (i.e. DataID[i]==MeshNodeID[i])
//
bool DataVar::filterSamples(FinleyMesh_ptr finleyMesh)
{
    if (numSamples == 0)
        return true;

    IndexMap id2idxMap;
    const IntVec* requiredIDs = NULL;

    NodeData_ptr nodes = finleyMesh->getMeshForFinleyFS(funcSpace);
    if (nodes == NULL)
        return false;

    int requiredNumSamples = 0;

    if (centering == NODE_CENTERED) {
        id2idxMap = nodes->getIndexMap();
        requiredIDs = &nodes->getNodeIDs();
        requiredNumSamples = nodes->getNumNodes();
    } else {
        ElementData_ptr cells = finleyMesh->getElementsForFinleyFS(funcSpace);
        if (cells == NULL)
            return false;

        id2idxMap = cells->getIndexMap();
        requiredIDs = &cells->getIDs();
        if (cells->getReducedNumElements() > 0) {
            requiredNumSamples = cells->getReducedNumElements();
        } else {
            requiredNumSamples = cells->getNumElements();
        }
    }

    if (requiredNumSamples > numSamples) {
        cerr << "ERROR: " << varName << " has " << numSamples
            << " instead of " << requiredNumSamples << " samples!" << endl;
        return false;
    }

    meshName = nodes->getName();
    siloMeshName = nodes->getFullSiloName();

    IndexMap sampleID2idx = buildIndexMap();
    numSamples = requiredNumSamples;

    // now filter the data
    for (size_t i=0; i < dataArray.size(); i++) {
        float* c = new float[numSamples];
        const float* src = dataArray[i];
        IntVec::const_iterator idIt = requiredIDs->begin();
        for (; idIt != requiredIDs->end(); idIt++) {
            size_t srcIdx = sampleID2idx.find(*idIt)->second;
            size_t destIdx = id2idxMap.find(*idIt)->second;
            c[destIdx] = src[srcIdx];
        }
        delete[] dataArray[i];
        dataArray[i] = c;
    }
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
    if (rank < 2 || !initialized)
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
#if USE_SILO
    if (!initialized)
        return false;

    if (numSamples == 0)
        return true;

    int ret;

    if (siloPath != "") {
        ret = DBSetDir(dbfile, siloPath.c_str());
        if (ret != 0)
            return false;
    }
 
    char* siloMesh = const_cast<char*>(siloMeshName.c_str());
    int dcenter = (centering == NODE_CENTERED ? DB_NODECENT : DB_ZONECENT);

    if (rank == 0) {
        ret = DBPutUcdvar1(dbfile, varName.c_str(), siloMesh, dataArray[0],
                numSamples, NULL, 0, DB_FLOAT, dcenter, NULL);
    }
    else if (rank == 1) {
        const string comps[3] = {
            varName+string("_x"), varName+string("_y"), varName+string("_z")
        };
        const char* varnames[3] = {
            comps[0].c_str(), comps[1].c_str(), comps[2].c_str()
        };

        ret = DBPutUcdvar(dbfile, varName.c_str(), siloMesh, shape[0],
                (char**)varnames, &dataArray[0], numSamples, NULL,
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
                    ret = DBPutUcdvar1(dbfile, varname.str().c_str(), siloMesh,
                            dataArray[i*shape[0]+j], numSamples,
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

#else // !USE_SILO
    return false;
#endif
}

} // namespace escriptexport

