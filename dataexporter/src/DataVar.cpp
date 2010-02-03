
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
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
#ifndef VISIT_PLUGIN
#include <escript/Data.h>
#endif

#if USE_NETCDF
#include <netcdfcpp.h>
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
    initialized(d.initialized), finleyMesh(d.finleyMesh),
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
#ifndef VISIT_PLUGIN
    cleanup();

    if (!escriptData.actsExpanded()) {
        cerr << "WARNING: Only expanded data supported!" << endl;
        return false;
    }

    finleyMesh = mesh;
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

    NodeData_ptr nodes = finleyMesh->getMeshForFinleyFS(funcSpace);
    if (nodes == NULL)
        return false;

    meshName = nodes->getName();
    siloMeshName = nodes->getFullSiloName();
    initialized = true;

    // no samples? Nothing more to do.
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

        initialized = reorderSamples();
    }

    return initialized;

#else // VISIT_PLUGIN
    return false;
#endif
}

//
// Initialise with mesh data
//
bool DataVar::initFromMesh(FinleyMesh_ptr mesh)
{
    cleanup();
    
    finleyMesh = mesh;
    rank = 0;
    ptsPerSample = 1;
    NodeData_ptr nodes;

    if (varName.find("ContactElements_") != varName.npos) {
        funcSpace = FINLEY_CONTACT_ELEMENTS_1;
        centering = ZONE_CENTERED;
        string elementName = varName.substr(0, varName.find('_'));
        ElementData_ptr elements = mesh->getElementsByName(elementName);
        nodes = elements->getNodeMesh();
        sampleID = elements->getIDs();
    } else if (varName.find("FaceElements_") != varName.npos) {
        funcSpace = FINLEY_FACE_ELEMENTS;
        centering = ZONE_CENTERED;
        string elementName = varName.substr(0, varName.find('_'));
        ElementData_ptr elements = mesh->getElementsByName(elementName);
        nodes = elements->getNodeMesh();
        sampleID = elements->getIDs();
    } else if (varName.find("Elements_") != varName.npos) {
        funcSpace = FINLEY_ELEMENTS;
        centering = ZONE_CENTERED;
        string elementName = varName.substr(0, varName.find('_'));
        ElementData_ptr elements = mesh->getElementsByName(elementName);
        nodes = elements->getNodeMesh();
        sampleID = elements->getIDs();
    } else if (varName.find("Nodes_") != varName.npos) {
        funcSpace = FINLEY_NODES;
        centering = NODE_CENTERED;
        nodes = mesh->getNodes();
        sampleID = nodes->getNodeIDs();
    } else {
        cerr << "WARNING: Unrecognized mesh variable '" << varName << "'\n";
        return false;
    }

    meshName = nodes->getName();
    siloMeshName = nodes->getFullSiloName();

    const IntVec& data = mesh->getVarDataByName(varName);
    numSamples = data.size();

    if (numSamples > 0) {
        float* c = new float[numSamples];
        dataArray.push_back(c);
        IntVec::const_iterator it;
        for (it=data.begin(); it != data.end(); it++)
            *c++ = static_cast<float>(*it);
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

    finleyMesh = mesh;
    NodeData_ptr nodes = finleyMesh->getMeshForFinleyFS(funcSpace);
    if (nodes == NULL) {
        delete input;
        return false;
    }

    meshName = nodes->getName();
    siloMeshName = nodes->getFullSiloName();
    initialized = true;

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
 
    if (initialized && numSamples > 0) {
        sampleID.insert(sampleID.end(), numSamples, 0);
        NcVar* var = input->get_var("id");
        var->get(&sampleID[0], numSamples);

        size_t dataSize = dimSize*numSamples*ptsPerSample;
        counts.push_back(ptsPerSample);
        counts.push_back(numSamples);
        float* tempData = new float[dataSize];
        var = input->get_var("data");
        var->get(tempData, &counts[0]);

        const float* srcPtr = tempData;
        for (int i=0; i < dimSize; i++, srcPtr++) {
            float* c = averageData(srcPtr, dimSize);
            dataArray.push_back(c);
        }
        delete[] tempData;

        initialized = reorderSamples();
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
    float* res;

    if (ptsPerSample == 1) {
        res = new float[numSamples];
        float* dest = res;
        for (int i=0; i<numSamples; i++, src+=stride)
            *dest++ = *src;
    } else {
        ElementData_ptr cells = finleyMesh->getElementsForFinleyFS(funcSpace);
        int cellFactor = cells->getElementFactor();
        res = new float[cellFactor * numSamples];
        float* dest = res;
        QuadMaskInfo qmi = cells->getQuadMask(funcSpace);
        if (qmi.mask.size() > 0) {
            const float* tmpSrc = src;
            for (int i=0; i<numSamples; i++, tmpSrc+=stride*ptsPerSample) {
                for (int l=0; l<cellFactor; l++) {
                    double tmpVal = 0.0;
                    for (int j=0; j<ptsPerSample; j++) {
                        if (qmi.mask[l][j] != 0) {
                            tmpVal += *(tmpSrc+stride*j);
                        }
                    }
                    *dest++ = (float)(tmpVal / qmi.factor[l]);
                }
            }
        } else {
            for (int i=0; i<numSamples; i++) {
                double tmpVal = 0.0;
                for (int j=0; j<ptsPerSample; j++, src+=stride) {
                    tmpVal += *src;
                }
                tmpVal /= ptsPerSample;
                for (int l=0; l<cellFactor; l++) {
                    *dest++ = static_cast<float>(tmpVal);
                }
            }
        }
    }
    return res;
}

//
// Filters and reorders the raw sample values according to the node/element
// IDs. This is used to have data arrays ordered according to the underlying
// mesh (i.e. DataID[i]==MeshNodeID[i])
//
bool DataVar::reorderSamples()
{
    if (numSamples == 0)
        return true;

    const IntVec* requiredIDs = NULL;
    int requiredNumSamples = 0;
    int cellFactor = 1;

    if (centering == NODE_CENTERED) {
        NodeData_ptr nodes = finleyMesh->getMeshForFinleyFS(funcSpace);
        requiredIDs = &nodes->getNodeIDs();
        requiredNumSamples = nodes->getNumNodes();
    } else {
        ElementData_ptr cells = finleyMesh->getElementsForFinleyFS(funcSpace);
        if (cells == NULL)
            return false;

        requiredIDs = &cells->getIDs();
        requiredNumSamples = cells->getNumElements();
        cellFactor = cells->getElementFactor();
        if (cellFactor > 1) {
            numSamples *= cellFactor;
            // update sample IDs
            IntVec newSampleID(numSamples);
            IntVec::const_iterator idIt = sampleID.begin();
            IntVec::iterator newIDit = newSampleID.begin();
            for (; idIt != sampleID.end(); idIt++, newIDit+=cellFactor) {
                fill(newIDit, newIDit+cellFactor, *idIt);
            }
            sampleID.swap(newSampleID);
        }
    }

    if (requiredNumSamples > numSamples) {
        cerr << "ERROR: " << varName << " has " << numSamples
            << " instead of " << requiredNumSamples << " samples!" << endl;
        return false;
    }

    IndexMap sampleID2idx = buildIndexMap();
    numSamples = requiredNumSamples;

    // now filter the data
    for (size_t i=0; i < dataArray.size(); i++) {
        float* c = new float[numSamples];
        const float* src = dataArray[i];
        IntVec::const_iterator idIt = requiredIDs->begin();
        size_t destIdx = 0;
        for (; idIt != requiredIDs->end(); idIt+=cellFactor, destIdx+=cellFactor) {
            size_t srcIdx = sampleID2idx.find(*idIt)->second;
            copy(&src[srcIdx], &src[srcIdx+cellFactor], &c[destIdx]);
        }
        delete[] dataArray[i];
        dataArray[i] = c;
    }

    // sample IDs now = mesh node/element IDs
    sampleID = *requiredIDs;

    return true;
}

//
//
//
void DataVar::sampleToStream(ostream& os, int index)
{
    if (rank == 0) {
        os << dataArray[0][index];
    } else if (rank == 1) {
        if (shape[0] < 3)
            os << dataArray[0][index] << " " << dataArray[1][index]
                << " " << 0.;
        else
            os << dataArray[0][index] << " " << dataArray[1][index]
                << " " << dataArray[2][index];
    } else if (rank == 2) {
        if (shape[1] < 3) {
            os << dataArray[0][index] << " " << dataArray[1][index]
                << " " << 0. << " ";
            os << dataArray[2][index] << " " << dataArray[3][index]
                << " " << 0. << " ";
            os << 0. << " " << 0. << " " << 0.;
        } else {
            os << dataArray[0][index] << " " << dataArray[1][index]
                << " " << dataArray[2][index] << " ";
            os << dataArray[3][index] << " " << dataArray[4][index]
                << " " << dataArray[5][index] << " ";
            os << dataArray[6][index] << " " << dataArray[7][index]
                << " " << dataArray[8][index];
        }
    }
    os << endl;
}

//
//
//
void DataVar::writeToVTK(ostream& os, int ownIndex)
{
    if (numSamples == 0)
        return;

    if (isNodeCentered()) {
        // data was reordered in reorderSamples() but for VTK we write the
        // original node mesh and thus need the original ordering...
        const IntVec& requiredIDs = finleyMesh->getNodes()->getNodeIDs();
        const IntVec& nodeGNI = finleyMesh->getNodes()->getGlobalNodeIndices();
        const IntVec& nodeDist = finleyMesh->getNodes()->getNodeDistribution();
        int firstId = nodeDist[ownIndex];
        int lastId = nodeDist[ownIndex+1];
        IndexMap sampleID2idx = buildIndexMap();
        for (int i=0; i<nodeGNI.size(); i++) {
            if (firstId <= nodeGNI[i] && nodeGNI[i] < lastId) {
                int idx = sampleID2idx[requiredIDs[i]];
                sampleToStream(os, idx);
            }
        }
    } else {
        // cell data: ghost cells have been removed so do not write ghost
        // samples (which are the last elements in the arrays)
        int toWrite =
            finleyMesh->getElementsByName(meshName)->getNumElements();
        for (int i=0; i<toWrite; i++) {
            sampleToStream(os, i);
        }
    }
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

