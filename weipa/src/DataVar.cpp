
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <weipa/DataVar.h>
#include <weipa/DomainChunk.h>
#include <weipa/ElementData.h>
#include <weipa/NodeData.h>

#include <escript/Utils.h>

#ifndef VISIT_PLUGIN
#include <escript/Data.h>
#endif


#ifdef ESYS_HAVE_SILO
#include <silo.h>
#endif

#ifdef ESYS_HAVE_NETCDF4
#include <ncFile.h>
#include <ncVar.h>
#include <ncDim.h>
#include <ncGroupAtt.h>
#endif

#include <numeric> // for accumulate
#include <iostream> // for cerr
#include <sstream>
#include <stdio.h>

using namespace std;
#ifdef ESYS_HAVE_NETCDF4
using namespace netCDF;
#endif

namespace weipa {
    
//
// Constructor
//
DataVar::DataVar(const string& name) :
    initialized(false), varName(name),
    numSamples(0), rank(0), ptsPerSample(0)
{
}

//
// Copy constructor
//
DataVar::DataVar(const DataVar& d) :
    initialized(d.initialized), domain(d.domain),
    varName(d.varName), numSamples(d.numSamples),
    rank(d.rank), ptsPerSample(d.ptsPerSample), funcSpace(d.funcSpace),
    centering(d.centering), shape(d.shape), sampleID(d.sampleID)
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
bool DataVar::initFromEscript(escript::Data& escriptData, const_DomainChunk_ptr dom)
{
#ifndef VISIT_PLUGIN
    cleanup();

    if (!escriptData.isConstant() && !escriptData.actsExpanded()) {
        cerr << "WARNING: Weipa only supports constant & expanded data, "
            << "not initializing " << varName << endl;
        return false;
    }

    domain = dom;
    rank = escriptData.getDataPointRank();
    ptsPerSample = escriptData.getNumDataPointsPerSample();
    shape = escriptData.getDataPointShape();
    funcSpace = escriptData.getFunctionSpace().getTypeCode();
    numSamples = escriptData.getNumSamples();
    centering = domain->getCenteringForFunctionSpace(funcSpace);

#ifdef _DEBUG
    cout << varName << ":\t" << numSamples << " samples,  "
        << ptsPerSample << " pts/s,  rank: " << rank << endl;
#endif

    NodeData_ptr nodes = domain->getMeshForFunctionSpace(funcSpace);
    if (nodes == NULL)
        return false;

    meshName = nodes->getName();
    siloMeshName = nodes->getFullSiloName();
    initialized = true;

    // no samples? Nothing more to do.
    if (numSamples == 0)
        return true;

    const escript::DataTypes::dim_t* iPtr = escriptData.getFunctionSpace().borrowSampleReferenceIDs();
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

    // special case: shape=(1,) or shape=(1,1) -> convert to scalar
    if (dimSize==1 && rank>0) {
        rank=0;
        shape.clear();
    }

    if (initialized) {
        size_t dataSize = dimSize * ptsPerSample;
        float* tempData = new float[dataSize*numSamples];
        float* destPtr = tempData;
        if (escriptData.isConstant()) {
            const escript::DataTypes::real_t* values =
                escriptData.getDataRO();
            for (int pointNo=0; pointNo<numSamples*ptsPerSample; pointNo++) {
                copy(values, values+dimSize, destPtr);
                destPtr += dimSize;
            }
        } else {
            for (int sampleNo=0; sampleNo<numSamples; sampleNo++) {
                const escript::DataTypes::real_t* values =
                    escriptData.getSampleDataRO(sampleNo);
                copy(values, values+dataSize, destPtr);
                destPtr += dataSize;
            }
        }

        const float* srcPtr = tempData;
        for (size_t i=0; i < dimSize; i++, srcPtr++) {
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
// Initialise with domain data
//
bool DataVar::initFromMeshData(const_DomainChunk_ptr dom, const IntVec& data,
        int fsCode, Centering c, NodeData_ptr nodes, const IntVec& id)
{
    cleanup();
    
    domain = dom;
    rank = 0;
    ptsPerSample = 1;
    centering = c;
    sampleID = id;
    meshName = nodes->getName();
    siloMeshName = nodes->getFullSiloName();
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
// Reads variable data from dump file
//
bool DataVar::initFromFile(const string& filename, const_DomainChunk_ptr dom)
{
    cleanup();

#ifdef ESYS_HAVE_NETCDF4

    NcFile input;
    if (!escript::openNcFile(input, filename))
    {
        cerr << "Could not open input file " << filename << "." << endl;
        return false;
    }      

    NcDim dim;
    NcGroupAtt att;

    att = input.getAtt("type_id");
    int typeID;
    att.getValues(&typeID);
    if (typeID != 2) {
        cerr << "WARNING: Only expanded data supported!" << endl;
        return false;
    }

    att = input.getAtt("rank");
    att.getValues(&rank);

    dim = input.getDim("num_data_points_per_sample");
    ptsPerSample = dim.getSize();

    att = input.getAtt("function_space_type");
    att.getValues(&funcSpace);

    centering = dom->getCenteringForFunctionSpace(funcSpace);

    dim = input.getDim("num_samples");
    numSamples = dim.getSize();

#ifdef _DEBUG
    cout << varName << ":\t" << numSamples << " samples,  "
        << ptsPerSample << " pts/s,  rank: " << rank << endl;
#endif

    domain = dom;
    NodeData_ptr nodes = domain->getMeshForFunctionSpace(funcSpace);
    if (nodes == NULL) {
        return false;
    }

    meshName = nodes->getName();
    siloMeshName = nodes->getFullSiloName();
    initialized = true;

    size_t dimSize = 1;
    vector<long> counts;

    if (rank > 0) {
        dim = input.getDim("d0");
        int d = dim.getSize();
        shape.push_back(d);
        counts.push_back(d);
        dimSize *= d;
    }
    if (rank > 1) {
        dim = input.getDim("d1");
        int d = dim.getSize();
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
        NcVar var = input.getVar("id");
        var.getVar(&sampleID[0]);   // numSamples

        size_t dataSize = dimSize*numSamples*ptsPerSample;
        counts.push_back(ptsPerSample);
        counts.push_back(numSamples);
        float* tempData = new float[dataSize];
        var = input.getVar("data");
        var.getVar(tempData);   // &counts[0]

        const float* srcPtr = tempData;
        for (size_t i=0; i < dimSize; i++, srcPtr++) {
            float* c = averageData(srcPtr, dimSize);
            dataArray.push_back(c);
        }
        delete[] tempData;

        initialized = reorderSamples();
    }
#endif // ESYS_HAVE_NETCDF4

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
        ElementData_ptr cells = domain->getElementsForFunctionSpace(funcSpace);
        int cellFactor = cells->getElementFactor();
        res = new float[cellFactor * numSamples];
        float* dest = res;
        QuadMaskInfo qmi = cells->getQuadMask(funcSpace);
        if (!qmi.mask.empty()) {
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
        NodeData_ptr nodes = domain->getMeshForFunctionSpace(funcSpace);
        requiredIDs = &nodes->getNodeIDs();
        requiredNumSamples = nodes->getNumNodes();
    } else {
        ElementData_ptr cells = domain->getElementsForFunctionSpace(funcSpace);
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
int DataVar::getNumberOfComponents() const
{
    return (rank == 0 ? 1 : accumulate(shape.begin(), shape.end(), 0));
}

//
//
//
float* DataVar::getDataFlat() const
{
    int totalSize = numSamples * getNumberOfComponents();
    float* res = new float[totalSize];
    if (rank == 0) {
        copy(dataArray[0], dataArray[0]+numSamples, res);
    } else if (rank == 1) {
        float *dest = res;
        for (size_t c=0; c<numSamples; c++) {
            for (size_t i=0; i<shape[0]; i++) {
                *dest++ = dataArray[i][c];
            }
        }
    } else if (rank == 2) {
        float *dest = res;
        for (size_t c=0; c<numSamples; c++) {
            for (int i=0; i<shape[1]; i++) {
                for (int j=0; j<shape[0]; j++) {
                    *dest++ = dataArray[i*shape[0]+j][c];
                }
            }
        }
    }

    return res;
}

//
//
//
void DataVar::sampleToStream(ostream& os, int index)
{
    // index is -1 for dummy samples, i.e. if writing the full mesh but
    // only a reduced number of samples is required
    if (rank == 0) {
        if (index < 0) {
            os << 0.;
        } else {
            os << dataArray[0][index];
        }
    } else if (rank == 1) {
        if (index < 0) {
            os << 0. << " " << 0.  << " " << 0.;
        } else if (shape[0] < 3) {
            os << dataArray[0][index] << " " << dataArray[1][index]
                << " " << 0.;
        } else {
            os << dataArray[0][index] << " " << dataArray[1][index]
                << " " << dataArray[2][index];
        }
    } else if (rank == 2) {
        if (index < 0) {
            os << 0. << " " << 0. << " " << 0. << " ";
            os << 0. << " " << 0. << " " << 0. << " ";
            os << 0. << " " << 0. << " " << 0.;
        } else if (shape[1] < 3) {
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
        // original node mesh and thus need the original ordering.
        // Note, that this also means we may not have samples for all nodes
        // in which case we set idx to -1 and write a dummy sample
        const IntVec& requiredIDs = domain->getNodes()->getNodeIDs();
        const IntVec& nodeGNI = domain->getNodes()->getGlobalNodeIndices();
        const IntVec& nodeDist = domain->getNodes()->getNodeDistribution();
        int firstId = nodeDist[ownIndex];
        int lastId = nodeDist[ownIndex+1];
        IndexMap sampleID2idx = buildIndexMap();
        for (size_t i=0; i<nodeGNI.size(); i++) {
            if (firstId <= nodeGNI[i] && nodeGNI[i] < lastId) {
                IndexMap::const_iterator it = sampleID2idx.find(requiredIDs[i]);
                int idx = (it==sampleID2idx.end() ? -1 : (int)it->second);
                sampleToStream(os, idx);
            }
        }
    } else {
        // cell data: ghost cells have been removed so do not write ghost
        // samples (which are the last elements in the arrays)
        int toWrite = domain->getElementsByName(meshName)->getNumElements();
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
    const int ltDef = tensor3DefFmt.length()+9*tensorDir.length();
    char* tDef = new char[ltDef];

    if (shape[1] == 3) {
        snprintf(tDef, ltDef, tensor3DefFmt.c_str(),
                tensorDir.c_str(), tensorDir.c_str(), tensorDir.c_str(),
                tensorDir.c_str(), tensorDir.c_str(), tensorDir.c_str(),
                tensorDir.c_str(), tensorDir.c_str(), tensorDir.c_str());
        tensorDef = tDef;
        delete[] tDef;
    } else {
        snprintf(tDef, ltDef, tensor2DefFmt.c_str(),
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
// The corresponding mesh must have been written already.
//
bool DataVar::writeToSilo(DBfile* dbfile, const string& siloPath,
                          const string& units)
{
#ifdef ESYS_HAVE_SILO
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
    DBoptlist* optList = DBMakeOptlist(2);
    if (units.length()>0) {
        DBAddOption(optList, DBOPT_UNITS, (void*)units.c_str());
    }

    if (rank == 0) {
        ret = DBPutUcdvar1(dbfile, varName.c_str(), siloMesh, dataArray[0],
                numSamples, NULL, 0, DB_FLOAT, dcenter, optList);
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
                0, DB_FLOAT, dcenter, optList);
    }
    else {
        string tensorDir = varName+string("_comps/");
        ret = DBMkdir(dbfile, tensorDir.c_str());
        if (ret == 0) {
            int one = 1;
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
        } // ret==0
    } // rank

    DBFreeOptlist(optList);
    DBSetDir(dbfile, "/");
    return (ret == 0);

#else // !ESYS_HAVE_SILO
    return false;
#endif
}

} // namespace weipa

