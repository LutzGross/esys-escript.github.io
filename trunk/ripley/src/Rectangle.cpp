
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/

#include <ripley/Rectangle.h>
#include <paso/SystemMatrix.h>
#include <esysUtils/esysFileWriter.h>

#include <boost/scoped_array.hpp>
#include "esysUtils/EsysRandom.h"

#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

#if USE_SILO
#include <silo.h>
#ifdef ESYS_MPI
#include <pmpio.h>
#endif
#endif

#include <iomanip>

using namespace std;
using esysUtils::FileWriter;

namespace ripley {

Rectangle::Rectangle(int n0, int n1, double x0, double y0, double x1,
                     double y1, int d0, int d1) :
    RipleyDomain(2)
{
    // ignore subdivision parameters for serial run
    if (m_mpiInfo->size == 1) {
        d0=1;
        d1=1;
    }

    bool warn=false;
    // if number of subdivisions is non-positive, try to subdivide by the same
    // ratio as the number of elements
    if (d0<=0 && d1<=0) {
        warn=true;
        d0=max(1, (int)(sqrt(m_mpiInfo->size*(n0+1)/(float)(n1+1))));
        d1=m_mpiInfo->size/d0;
        if (d0*d1 != m_mpiInfo->size) {
            // ratios not the same so subdivide side with more elements only
            if (n0>n1) {
                d0=0;
                d1=1;
            } else {
                d0=1;
                d1=0;
            }
        }
    }
    if (d0<=0) {
        // d1 is preset, determine d0 - throw further down if result is no good
        d0=m_mpiInfo->size/d1;
    } else if (d1<=0) {
        // d0 is preset, determine d1 - throw further down if result is no good
        d1=m_mpiInfo->size/d0;
    }

    // ensure number of subdivisions is valid and nodes can be distributed
    // among number of ranks
    if (d0*d1 != m_mpiInfo->size)
        throw RipleyException("Invalid number of spatial subdivisions");

    if (warn) {
        cout << "Warning: Automatic domain subdivision (d0=" << d0 << ", d1="
            << d1 << "). This may not be optimal!" << endl;
    }

    double l0 = x1-x0;
    double l1 = y1-y0;
    m_dx[0] = l0/n0;
    m_dx[1] = l1/n1;

    if ((n0+1)%d0 > 0) {
        n0=(int)round((float)(n0+1)/d0+0.5)*d0-1;
        l0=m_dx[0]*n0;
        cout << "Warning: Adjusted number of elements and length. N0="
            << n0 << ", l0=" << l0 << endl;
    }
    if ((n1+1)%d1 > 0) {
        n1=(int)round((float)(n1+1)/d1+0.5)*d1-1;
        l1=m_dx[1]*n1;
        cout << "Warning: Adjusted number of elements and length. N1="
            << n1 << ", l1=" << l1 << endl;
    }

    if ((d0 > 1 && (n0+1)/d0<2) || (d1 > 1 && (n1+1)/d1<2))
        throw RipleyException("Too few elements for the number of ranks");

    m_gNE[0] = n0;
    m_gNE[1] = n1;
    m_origin[0] = x0;
    m_origin[1] = y0;
    m_length[0] = l0;
    m_length[1] = l1;
    m_NX[0] = d0;
    m_NX[1] = d1;

    // local number of elements (with and without overlap)
    m_NE[0] = m_ownNE[0] = (d0>1 ? (n0+1)/d0 : n0);
    if (m_mpiInfo->rank%d0>0 && m_mpiInfo->rank%d0<d0-1)
        m_NE[0]++;
    else if (d0>1 && m_mpiInfo->rank%d0==d0-1)
        m_ownNE[0]--;

    m_NE[1] = m_ownNE[1] = (d1>1 ? (n1+1)/d1 : n1);
    if (m_mpiInfo->rank/d0>0 && m_mpiInfo->rank/d0<d1-1)
        m_NE[1]++;
    else if (d1>1 && m_mpiInfo->rank/d0==d1-1)
        m_ownNE[1]--;

    // local number of nodes
    m_NN[0] = m_NE[0]+1;
    m_NN[1] = m_NE[1]+1;

    // bottom-left node is at (offset0,offset1) in global mesh
    m_offset[0] = (n0+1)/d0*(m_mpiInfo->rank%d0);
    if (m_offset[0] > 0)
        m_offset[0]--;
    m_offset[1] = (n1+1)/d1*(m_mpiInfo->rank/d0);
    if (m_offset[1] > 0)
        m_offset[1]--;

    populateSampleIds();
    createPattern();
}

Rectangle::~Rectangle()
{
    Paso_SystemMatrixPattern_free(m_pattern);
    Paso_Connector_free(m_connector);
}

string Rectangle::getDescription() const
{
    return "ripley::Rectangle";
}

bool Rectangle::operator==(const AbstractDomain& other) const
{
    const Rectangle* o=dynamic_cast<const Rectangle*>(&other);
    if (o) {
        return (RipleyDomain::operator==(other) &&
                m_gNE[0]==o->m_gNE[0] && m_gNE[1]==o->m_gNE[1]
                && m_origin[0]==o->m_origin[0] && m_origin[1]==o->m_origin[1]
                && m_length[0]==o->m_length[0] && m_length[1]==o->m_length[1]
                && m_NX[0]==o->m_NX[0] && m_NX[1]==o->m_NX[1]);
    }

    return false;
}

void Rectangle::readNcGrid(escript::Data& out, string filename, string varname,
            const vector<int>& first, const vector<int>& numValues,
            const vector<int>& multiplier) const
{
#ifdef USE_NETCDF
    // check destination function space
    int myN0, myN1;
    if (out.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NN[0];
        myN1 = m_NN[1];
    } else if (out.getFunctionSpace().getTypeCode() == Elements ||
                out.getFunctionSpace().getTypeCode() == ReducedElements) {
        myN0 = m_NE[0];
        myN1 = m_NE[1];
    } else
        throw RipleyException("readNcGrid(): invalid function space for output data object");

    if (first.size() != 2)
        throw RipleyException("readNcGrid(): argument 'first' must have 2 entries");

    if (numValues.size() != 2)
        throw RipleyException("readNcGrid(): argument 'numValues' must have 2 entries");

    if (multiplier.size() != 2)
        throw RipleyException("readNcGrid(): argument 'multiplier' must have 2 entries");
    for (size_t i=0; i<multiplier.size(); i++)
        if (multiplier[i]<1)
            throw RipleyException("readNcGrid(): all multipliers must be positive");

    // check file existence and size
    NcFile f(filename.c_str(), NcFile::ReadOnly);
    if (!f.is_valid())
        throw RipleyException("readNcGrid(): cannot open file");

    NcVar* var = f.get_var(varname.c_str());
    if (!var)
        throw RipleyException("readNcGrid(): invalid variable");

    // TODO: rank>0 data support
    const int numComp = out.getDataPointSize();
    if (numComp > 1)
        throw RipleyException("readNcGrid(): only scalar data supported");

    const int dims = var->num_dims();
    boost::scoped_array<long> edges(var->edges());

    // is this a slice of the data object (dims!=2)?
    // note the expected ordering of edges (as in numpy: y,x)
    if ( (dims==2 && (numValues[1] > edges[0] || numValues[0] > edges[1]))
            || (dims==1 && numValues[1]>1) ) {
        throw RipleyException("readNcGrid(): not enough data in file");
    }

    // check if this rank contributes anything
    if (first[0] >= m_offset[0]+myN0 || first[0]+numValues[0]*multiplier[0] <= m_offset[0] ||
            first[1] >= m_offset[1]+myN1 || first[1]+numValues[1]*multiplier[1] <= m_offset[1])
        return;

    // now determine how much this rank has to write

    // first coordinates in data object to write to
    const int first0 = max(0, first[0]-m_offset[0]);
    const int first1 = max(0, first[1]-m_offset[1]);
    // indices to first value in file
    const int idx0 = max(0, m_offset[0]-first[0]);
    const int idx1 = max(0, m_offset[1]-first[1]);
    // number of values to read
    const int num0 = min(numValues[0]-idx0, myN0-first0);
    const int num1 = min(numValues[1]-idx1, myN1-first1);

    vector<double> values(num0*num1);
    if (dims==2) {
        var->set_cur(idx1, idx0);
        var->get(&values[0], num1, num0);
    } else {
        var->set_cur(idx0);
        var->get(&values[0], num0);
    }

    const int dpp = out.getNumDataPointsPerSample();
    out.requireWrite();

    for (index_t y=0; y<num1; y++) {
#pragma omp parallel for
        for (index_t x=0; x<num0; x++) {
            const int baseIndex = first0+x*multiplier[0]
                                  +(first1+y*multiplier[1])*myN0;
            const int srcIndex = y*num0+x;
            if (!isnan(values[srcIndex])) {
                for (index_t m1=0; m1<multiplier[1]; m1++) {
                    for (index_t m0=0; m0<multiplier[0]; m0++) {
                        const int dataIndex = baseIndex+m0+m1*myN0;
                        double* dest = out.getSampleDataRW(dataIndex);
                        for (index_t q=0; q<dpp; q++) {
                            *dest++ = values[srcIndex];
                        }
                    }
                }
            }
        }
    }
#else
    throw RipleyException("readNcGrid(): not compiled with netCDF support");
#endif
}

void Rectangle::readBinaryGrid(escript::Data& out, string filename,
                               const vector<int>& first,
                               const vector<int>& numValues,
                               const std::vector<int>& multiplier,
                               int byteOrder, int dataType) const
{
    // the mapping is not universally correct but should work on our
    // supported platforms
    switch (dataType) {
        case DATATYPE_INT32:
            readBinaryGridImpl<int>(out, filename, first, numValues,
                                    multiplier, byteOrder);
            break;
        case DATATYPE_FLOAT32:
            readBinaryGridImpl<float>(out, filename, first, numValues,
                                      multiplier, byteOrder);
            break;
        case DATATYPE_FLOAT64:
            readBinaryGridImpl<double>(out, filename, first, numValues,
                                       multiplier, byteOrder);
            break;
        default:
            throw RipleyException("readBinaryGrid(): invalid or unsupported datatype");
    }
}

template<typename ValueType>
void Rectangle::readBinaryGridImpl(escript::Data& out, const string& filename,
                                   const vector<int>& first,
                                   const vector<int>& numValues,
                                   const std::vector<int>& multiplier,
                                   int byteOrder) const
{
    // check destination function space
    int myN0, myN1;
    if (out.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NN[0];
        myN1 = m_NN[1];
    } else if (out.getFunctionSpace().getTypeCode() == Elements ||
                out.getFunctionSpace().getTypeCode() == ReducedElements) {
        myN0 = m_NE[0];
        myN1 = m_NE[1];
    } else
        throw RipleyException("readBinaryGrid(): invalid function space for output data object");

    // check file existence and size
    ifstream f(filename.c_str(), ifstream::binary);
    if (f.fail()) {
        throw RipleyException("readBinaryGrid(): cannot open file");
    }
    f.seekg(0, ios::end);
    const int numComp = out.getDataPointSize();
    const int filesize = f.tellg();
    const int reqsize = numValues[0]*numValues[1]*numComp*sizeof(ValueType);
    if (filesize < reqsize) {
        f.close();
        throw RipleyException("readBinaryGrid(): not enough data in file");
    }

    // check if this rank contributes anything
    if (first[0] >= m_offset[0]+myN0 || first[0]+numValues[0] <= m_offset[0] ||
            first[1] >= m_offset[1]+myN1 || first[1]+numValues[1] <= m_offset[1]) {
        f.close();
        return;
    }

    // now determine how much this rank has to write

    // first coordinates in data object to write to
    const int first0 = max(0, first[0]-m_offset[0]);
    const int first1 = max(0, first[1]-m_offset[1]);
    // indices to first value in file
    const int idx0 = max(0, m_offset[0]-first[0]);
    const int idx1 = max(0, m_offset[1]-first[1]);
    // number of values to read
    const int num0 = min(numValues[0]-idx0, myN0-first0);
    const int num1 = min(numValues[1]-idx1, myN1-first1);

    out.requireWrite();
    vector<ValueType> values(num0*numComp);
    const int dpp = out.getNumDataPointsPerSample();

    for (int y=0; y<num1; y++) {
        const int fileofs = numComp*(idx0+(idx1+y)*numValues[0]);
        f.seekg(fileofs*sizeof(ValueType));
        f.read((char*)&values[0], num0*numComp*sizeof(ValueType));
        for (int x=0; x<num0; x++) {
            const int baseIndex = first0+x*multiplier[0]
                                    +(first1+y*multiplier[1])*myN0;
            for (int m1=0; m1<multiplier[1]; m1++) {
                for (int m0=0; m0<multiplier[0]; m0++) {
                    const int dataIndex = baseIndex+m0+m1*myN0;
                    double* dest = out.getSampleDataRW(dataIndex);
                    for (int c=0; c<numComp; c++) {
                        ValueType val = values[x*numComp+c];

                        if (byteOrder != BYTEORDER_NATIVE) {
                            char* cval = reinterpret_cast<char*>(&val);
                            // this will alter val!!
                            byte_swap32(cval);
                        }
                        if (!std::isnan(val)) {
                            for (int q=0; q<dpp; q++) {
                                *dest++ = static_cast<double>(val);
                            }
                        }
                    }
                }
            }
        }
    }

    f.close();
}

void Rectangle::writeBinaryGrid(const escript::Data& in, string filename,
                                int byteOrder, int dataType) const
{
    // the mapping is not universally correct but should work on our
    // supported platforms
    switch (dataType) {
        case DATATYPE_INT32:
            writeBinaryGridImpl<int>(in, filename, byteOrder);
            break;
        case DATATYPE_FLOAT32:
            writeBinaryGridImpl<float>(in, filename, byteOrder);
            break;
        case DATATYPE_FLOAT64:
            writeBinaryGridImpl<double>(in, filename, byteOrder);
            break;
        default:
            throw RipleyException("writeBinaryGrid(): invalid or unsupported datatype");
    }
}

template<typename ValueType>
void Rectangle::writeBinaryGridImpl(const escript::Data& in,
                                    const string& filename, int byteOrder) const
{
    // check function space and determine number of points
    int myN0, myN1;
    int totalN0, totalN1;
    if (in.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NN[0];
        myN1 = m_NN[1];
        totalN0 = m_gNE[0]+1;
        totalN1 = m_gNE[1]+1;
    } else if (in.getFunctionSpace().getTypeCode() == Elements ||
                in.getFunctionSpace().getTypeCode() == ReducedElements) {
        myN0 = m_NE[0];
        myN1 = m_NE[1];
        totalN0 = m_gNE[0];
        totalN1 = m_gNE[1];
    } else
        throw RipleyException("writeBinaryGrid(): invalid function space of data object");

    const int numComp = in.getDataPointSize();
    const int dpp = in.getNumDataPointsPerSample();

    if (numComp > 1 || dpp > 1)
        throw RipleyException("writeBinaryGrid(): only scalar, single-value data supported");

    escript::Data* _in = const_cast<escript::Data*>(&in);
    const int fileSize = sizeof(ValueType)*numComp*dpp*totalN0*totalN1;

    // from here on we know that each sample consists of one value
    FileWriter fw;
    fw.openFile(filename, fileSize);
    MPIBarrier();

    for (index_t y=0; y<myN1; y++) {
        const int fileofs = (m_offset[0]+(m_offset[1]+y)*totalN0)*sizeof(ValueType);
        ostringstream oss;

        for (index_t x=0; x<myN0; x++) {
            const double* sample = _in->getSampleDataRO(y*myN0+x);
            ValueType fvalue = static_cast<ValueType>(*sample);
            if (byteOrder == BYTEORDER_NATIVE) {
                oss.write((char*)&fvalue, sizeof(fvalue));
            } else {
                char* value = reinterpret_cast<char*>(&fvalue);
                oss.write(byte_swap32(value), sizeof(fvalue));
            }
        }
        fw.writeAt(oss, fileofs);
    }
    fw.close();
}

void Rectangle::dump(const string& fileName) const
{
#if USE_SILO
    string fn(fileName);
    if (fileName.length() < 6 || fileName.compare(fileName.length()-5, 5, ".silo") != 0) {
        fn+=".silo";
    }

    int driver=DB_HDF5;    
    DBfile* dbfile = NULL;
    const char* blockDirFmt = "/block%04d";

#ifdef ESYS_MPI
    PMPIO_baton_t* baton = NULL;
    const int NUM_SILO_FILES = 1;
#endif

    if (m_mpiInfo->size > 1) {
#ifdef ESYS_MPI
        baton = PMPIO_Init(NUM_SILO_FILES, PMPIO_WRITE, m_mpiInfo->comm,
                    0x1337, PMPIO_DefaultCreate, PMPIO_DefaultOpen,
                    PMPIO_DefaultClose, (void*)&driver);
        // try the fallback driver in case of error
        if (!baton && driver != DB_PDB) {
            driver = DB_PDB;
            baton = PMPIO_Init(NUM_SILO_FILES, PMPIO_WRITE, m_mpiInfo->comm,
                        0x1338, PMPIO_DefaultCreate, PMPIO_DefaultOpen,
                        PMPIO_DefaultClose, (void*)&driver);
        }
        if (baton) {
            char siloPath[64];
            snprintf(siloPath, 64, blockDirFmt, PMPIO_RankInGroup(baton, m_mpiInfo->rank));
            dbfile = (DBfile*) PMPIO_WaitForBaton(baton, fn.c_str(), siloPath);
        }
#endif
    } else {
        dbfile = DBCreate(fn.c_str(), DB_CLOBBER, DB_LOCAL,
                getDescription().c_str(), driver);
        // try the fallback driver in case of error
        if (!dbfile && driver != DB_PDB) {
            driver = DB_PDB;
            dbfile = DBCreate(fn.c_str(), DB_CLOBBER, DB_LOCAL,
                    getDescription().c_str(), driver);
        }
        char siloPath[64];
        snprintf(siloPath, 64, blockDirFmt, 0);
        DBMkDir(dbfile, siloPath);
        DBSetDir(dbfile, siloPath);
    }

    if (!dbfile)
        throw RipleyException("dump: Could not create Silo file");

    /*
    if (driver==DB_HDF5) {
        // gzip level 1 already provides good compression with minimal
        // performance penalty. Some tests showed that gzip levels >3 performed
        // rather badly on escript data both in terms of time and space
        DBSetCompression("ERRMODE=FALLBACK METHOD=GZIP LEVEL=1");
    }
    */

    boost::scoped_ptr<double> x(new double[m_NN[0]]);
    boost::scoped_ptr<double> y(new double[m_NN[1]]);
    double* coords[2] = { x.get(), y.get() };
#pragma omp parallel
    {
#pragma omp for nowait
        for (dim_t i0 = 0; i0 < m_NN[0]; i0++) {
            coords[0][i0]=getLocalCoordinate(i0, 0);
        }
#pragma omp for nowait
        for (dim_t i1 = 0; i1 < m_NN[1]; i1++) {
            coords[1][i1]=getLocalCoordinate(i1, 1);
        }
    }
    int* dims = const_cast<int*>(getNumNodesPerDim());

    // write mesh
    DBPutQuadmesh(dbfile, "mesh", NULL, coords, dims, 2, DB_DOUBLE,
            DB_COLLINEAR, NULL);

    // write node ids
    DBPutQuadvar1(dbfile, "nodeId", "mesh", (void*)&m_nodeId[0], dims, 2,
            NULL, 0, DB_INT, DB_NODECENT, NULL);

    // write element ids
    dims = const_cast<int*>(getNumElementsPerDim());
    DBPutQuadvar1(dbfile, "elementId", "mesh", (void*)&m_elementId[0],
            dims, 2, NULL, 0, DB_INT, DB_ZONECENT, NULL);

    // rank 0 writes multimesh and multivar
    if (m_mpiInfo->rank == 0) {
        vector<string> tempstrings;
        vector<char*> names;
        for (dim_t i=0; i<m_mpiInfo->size; i++) {
            stringstream path;
            path << "/block" << setw(4) << setfill('0') << right << i << "/mesh";
            tempstrings.push_back(path.str());
            names.push_back((char*)tempstrings.back().c_str());
        }
        vector<int> types(m_mpiInfo->size, DB_QUAD_RECT);
        DBSetDir(dbfile, "/");
        DBPutMultimesh(dbfile, "multimesh", m_mpiInfo->size, &names[0],
               &types[0], NULL);
        tempstrings.clear();
        names.clear();
        for (dim_t i=0; i<m_mpiInfo->size; i++) {
            stringstream path;
            path << "/block" << setw(4) << setfill('0') << right << i << "/nodeId";
            tempstrings.push_back(path.str());
            names.push_back((char*)tempstrings.back().c_str());
        }
        types.assign(m_mpiInfo->size, DB_QUADVAR);
        DBPutMultivar(dbfile, "nodeId", m_mpiInfo->size, &names[0],
               &types[0], NULL);
        tempstrings.clear();
        names.clear();
        for (dim_t i=0; i<m_mpiInfo->size; i++) {
            stringstream path;
            path << "/block" << setw(4) << setfill('0') << right << i << "/elementId";
            tempstrings.push_back(path.str());
            names.push_back((char*)tempstrings.back().c_str());
        }
        DBPutMultivar(dbfile, "elementId", m_mpiInfo->size, &names[0],
               &types[0], NULL);
    }

    if (m_mpiInfo->size > 1) {
#ifdef ESYS_MPI
        PMPIO_HandOffBaton(baton, dbfile);
        PMPIO_Finish(baton);
#endif
    } else {
        DBClose(dbfile);
    }

#else // USE_SILO
    throw RipleyException("dump: no Silo support");
#endif
}

const int* Rectangle::borrowSampleReferenceIDs(int fsType) const
{
    switch (fsType) {
        case Nodes:
        case ReducedNodes: // FIXME: reduced
            return &m_nodeId[0];
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom: // FIXME: reduced
            return &m_dofId[0];
        case Elements:
        case ReducedElements:
            return &m_elementId[0];
        case FaceElements:
        case ReducedFaceElements:
            return &m_faceId[0];
        default:
            break;
    }

    stringstream msg;
    msg << "borrowSampleReferenceIDs: invalid function space type " << fsType;
    throw RipleyException(msg.str());
}

bool Rectangle::ownSample(int fsType, index_t id) const
{
    if (getMPISize()==1)
        return true;

    switch (fsType) {
        case Nodes:
        case ReducedNodes: // FIXME: reduced
            return (m_dofMap[id] < getNumDOF());
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
            return true;
        case Elements:
        case ReducedElements:
            // check ownership of element's bottom left node
            return (m_dofMap[id%m_NE[0]+m_NN[0]*(id/m_NE[0])] < getNumDOF());
        case FaceElements:
        case ReducedFaceElements:
            {
                // determine which face the sample belongs to before
                // checking ownership of corresponding element's first node
                dim_t n=0;
                for (size_t i=0; i<4; i++) {
                    n+=m_faceCount[i];
                    if (id<n) {
                        index_t k;
                        if (i==1)
                            k=m_NN[0]-2;
                        else if (i==3)
                            k=m_NN[0]*(m_NN[1]-2);
                        else
                            k=0;
                        // determine whether to move right or up
                        const index_t delta=(i/2==0 ? m_NN[0] : 1);
                        return (m_dofMap[k+(id-n+m_faceCount[i])*delta] < getNumDOF());
                    }
                }
                return false;
            }
        default:
            break;
    }

    stringstream msg;
    msg << "ownSample: invalid function space type " << fsType;
    throw RipleyException(msg.str());
}

void Rectangle::setToNormal(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == FaceElements) {
        out.requireWrite();
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k1);
                    // set vector at two quadrature points
                    *o++ = -1.;
                    *o++ = 0.;
                    *o++ = -1.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k1);
                    // set vector at two quadrature points
                    *o++ = 1.;
                    *o++ = 0.;
                    *o++ = 1.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k0);
                    // set vector at two quadrature points
                    *o++ = 0.;
                    *o++ = -1.;
                    *o++ = 0.;
                    *o = -1.;
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k0);
                    // set vector at two quadrature points
                    *o++ = 0.;
                    *o++ = 1.;
                    *o++ = 0.;
                    *o = 1.;
                }
            }
        } // end of parallel section
    } else if (out.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        out.requireWrite();
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k1);
                    *o++ = -1.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k1);
                    *o++ = 1.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k0);
                    *o++ = 0.;
                    *o = -1.;
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k0);
                    *o++ = 0.;
                    *o = 1.;
                }
            }
        } // end of parallel section

    } else {
        stringstream msg;
        msg << "setToNormal: invalid function space type "
            << out.getFunctionSpace().getTypeCode();
        throw RipleyException(msg.str());
    }
}

void Rectangle::setToSize(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == Elements
            || out.getFunctionSpace().getTypeCode() == ReducedElements) {
        out.requireWrite();
        const dim_t numQuad=out.getNumDataPointsPerSample();
        const double size=sqrt(m_dx[0]*m_dx[0]+m_dx[1]*m_dx[1]);
#pragma omp parallel for
        for (index_t k = 0; k < getNumElements(); ++k) {
            double* o = out.getSampleDataRW(k);
            fill(o, o+numQuad, size);
        }
    } else if (out.getFunctionSpace().getTypeCode() == FaceElements
            || out.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        out.requireWrite();
        const dim_t numQuad=out.getNumDataPointsPerSample();
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k1);
                    fill(o, o+numQuad, m_dx[1]);
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k1);
                    fill(o, o+numQuad, m_dx[1]);
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k0);
                    fill(o, o+numQuad, m_dx[0]);
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k0);
                    fill(o, o+numQuad, m_dx[0]);
                }
            }
        } // end of parallel section

    } else {
        stringstream msg;
        msg << "setToSize: invalid function space type "
            << out.getFunctionSpace().getTypeCode();
        throw RipleyException(msg.str());
    }
}

void Rectangle::Print_Mesh_Info(const bool full) const
{
    RipleyDomain::Print_Mesh_Info(full);
    if (full) {
        cout << "     Id  Coordinates" << endl;
        cout.precision(15);
        cout.setf(ios::scientific, ios::floatfield);
        for (index_t i=0; i < getNumNodes(); i++) {
            cout << "  " << setw(5) << m_nodeId[i]
                << "  " << getLocalCoordinate(i%m_NN[0], 0)
                << "  " << getLocalCoordinate(i/m_NN[0], 1) << endl;
        }
    }
}


//protected
void Rectangle::assembleCoordinates(escript::Data& arg) const
{
    escriptDataC x = arg.getDataC();
    int numDim = m_numDim;
    if (!isDataPointShapeEqual(&x, 1, &numDim))
        throw RipleyException("setToX: Invalid Data object shape");
    if (!numSamplesEqual(&x, 1, getNumNodes()))
        throw RipleyException("setToX: Illegal number of samples in Data object");

    arg.requireWrite();
#pragma omp parallel for
    for (dim_t i1 = 0; i1 < m_NN[1]; i1++) {
        for (dim_t i0 = 0; i0 < m_NN[0]; i0++) {
            double* point = arg.getSampleDataRW(i0+m_NN[0]*i1);
            point[0] = getLocalCoordinate(i0, 0);
            point[1] = getLocalCoordinate(i1, 1);
        }
    }
}

//protected
void Rectangle::assembleGradient(escript::Data& out, escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    const double cx0 = .21132486540518711775/m_dx[0];
    const double cx1 = .78867513459481288225/m_dx[0];
    const double cx2 = 1./m_dx[0];
    const double cy0 = .21132486540518711775/m_dx[1];
    const double cy1 = .78867513459481288225/m_dx[1];
    const double cy2 = 1./m_dx[1];

    if (out.getFunctionSpace().getTypeCode() == Elements) {
        out.requireWrite();
#pragma omp parallel
        {
            vector<double> f_00(numComp);
            vector<double> f_01(numComp);
            vector<double> f_10(numComp);
            vector<double> f_11(numComp);
#pragma omp for
            for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,k1+1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,k1+1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(INDEX2(k0,k1,m_NE[0]));
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx1 + (f_11[i]-f_01[i])*cx0;
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy1 + (f_11[i]-f_10[i])*cy0;
                        o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx1 + (f_11[i]-f_01[i])*cx0;
                        o[INDEX3(i,1,1,numComp,2)] = (f_01[i]-f_00[i])*cy0 + (f_11[i]-f_10[i])*cy1;
                        o[INDEX3(i,0,2,numComp,2)] = (f_10[i]-f_00[i])*cx0 + (f_11[i]-f_01[i])*cx1;
                        o[INDEX3(i,1,2,numComp,2)] = (f_01[i]-f_00[i])*cy1 + (f_11[i]-f_10[i])*cy0;
                        o[INDEX3(i,0,3,numComp,2)] = (f_10[i]-f_00[i])*cx0 + (f_11[i]-f_01[i])*cx1;
                        o[INDEX3(i,1,3,numComp,2)] = (f_01[i]-f_00[i])*cy0 + (f_11[i]-f_10[i])*cy1;
                    } // end of component loop i
                } // end of k0 loop
            } // end of k1 loop
        } // end of parallel section
    } else if (out.getFunctionSpace().getTypeCode() == ReducedElements) {
        out.requireWrite();
#pragma omp parallel
        {
            vector<double> f_00(numComp);
            vector<double> f_01(numComp);
            vector<double> f_10(numComp);
            vector<double> f_11(numComp);
#pragma omp for
            for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,k1+1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,k1+1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(INDEX2(k0,k1,m_NE[0]));
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i] + f_11[i] - f_00[i] - f_01[i])*cx2/2;
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i] + f_11[i] - f_00[i] - f_10[i])*cy2/2;
                    } // end of component loop i
                } // end of k0 loop
            } // end of k1 loop
        } // end of parallel section
    } else if (out.getFunctionSpace().getTypeCode() == FaceElements) {
        out.requireWrite();
#pragma omp parallel
        {
            vector<double> f_00(numComp);
            vector<double> f_01(numComp);
            vector<double> f_10(numComp);
            vector<double> f_11(numComp);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(0,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(0,k1+1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(1,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(1,k1+1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx1 + (f_11[i]-f_01[i])*cx0;
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy2;
                        o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx0 + (f_11[i]-f_01[i])*cx1;
                        o[INDEX3(i,1,1,numComp,2)] = (f_01[i]-f_00[i])*cy2;
                    } // end of component loop i
                } // end of k1 loop
            } // end of face 0
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(m_NN[0]-2,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(m_NN[0]-2,k1+1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1+1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx1 + (f_11[i]-f_01[i])*cx0;
                        o[INDEX3(i,1,0,numComp,2)] = (f_11[i]-f_10[i])*cy2;
                        o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx0 + (f_11[i]-f_01[i])*cx1;
                        o[INDEX3(i,1,1,numComp,2)] = (f_11[i]-f_10[i])*cy2;
                    } // end of component loop i
                } // end of k1 loop
            } // end of face 1
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,0, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,0, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx2;
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy1 + (f_11[i]-f_10[i])*cy0;
                        o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx2;
                        o[INDEX3(i,1,1,numComp,2)] = (f_01[i]-f_00[i])*cy0 + (f_11[i]-f_10[i])*cy1;
                    } // end of component loop i
                } // end of k0 loop
            } // end of face 2
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,m_NN[1]-2, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,m_NN[1]-1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,m_NN[1]-2, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,m_NN[1]-1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_11[i]-f_01[i])*cx2;
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy1 + (f_11[i]-f_10[i])*cy0;
                        o[INDEX3(i,0,1,numComp,2)] = (f_11[i]-f_01[i])*cx2;
                        o[INDEX3(i,1,1,numComp,2)] = (f_01[i]-f_00[i])*cy0 + (f_11[i]-f_10[i])*cy1;
                    } // end of component loop i
                } // end of k0 loop
            } // end of face 3
        } // end of parallel section

    } else if (out.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        out.requireWrite();
#pragma omp parallel
        {
            vector<double> f_00(numComp);
            vector<double> f_01(numComp);
            vector<double> f_10(numComp);
            vector<double> f_11(numComp);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(0,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(0,k1+1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(1,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(1,k1+1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i] + f_11[i] - f_00[i] - f_01[i])*cx2/2;
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy2;
                    } // end of component loop i
                } // end of k1 loop
            } // end of face 0
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(m_NN[0]-2,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(m_NN[0]-2,k1+1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1+1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i] + f_11[i] - f_00[i] - f_01[i])*cx2/2;
                        o[INDEX3(i,1,0,numComp,2)] = (f_11[i]-f_10[i])*cy2;
                    } // end of component loop i
                } // end of k1 loop
            } // end of face 1
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,0, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,0, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx2;
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i] + f_11[i] - f_00[i] - f_10[i])*cy2/2;
                    } // end of component loop i
                } // end of k0 loop
            } // end of face 2
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,m_NN[1]-2, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,m_NN[1]-1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,m_NN[1]-2, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,m_NN[1]-1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_11[i]-f_01[i])*cx2;
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i] + f_11[i] - f_00[i] - f_10[i])*cy2/2;
                    } // end of component loop i
                } // end of k0 loop
            } // end of face 3
        } // end of parallel section
    }
}

//protected
void Rectangle::assembleIntegrate(vector<double>& integrals, escript::Data& arg) const
{
    const dim_t numComp = arg.getDataPointSize();
    const index_t left = (m_offset[0]==0 ? 0 : 1);
    const index_t bottom = (m_offset[1]==0 ? 0 : 1);
    const int fs=arg.getFunctionSpace().getTypeCode();
    if (fs == Elements && arg.actsExpanded()) {
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
            const double w = m_dx[0]*m_dx[1]/4.;
#pragma omp for nowait
            for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                    const double* f = arg.getSampleDataRO(INDEX2(k0, k1, m_NE[0]));
                    for (index_t i=0; i < numComp; ++i) {
                        const double f0 = f[INDEX2(i,0,numComp)];
                        const double f1 = f[INDEX2(i,1,numComp)];
                        const double f2 = f[INDEX2(i,2,numComp)];
                        const double f3 = f[INDEX2(i,3,numComp)];
                        int_local[i]+=(f0+f1+f2+f3)*w;
                    }  // end of component loop i
                } // end of k0 loop
            } // end of k1 loop
#pragma omp critical
            for (index_t i=0; i<numComp; i++)
                integrals[i]+=int_local[i];
        } // end of parallel section

    } else if (fs==ReducedElements || (fs==Elements && !arg.actsExpanded())) {
        const double w = m_dx[0]*m_dx[1];
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
#pragma omp for nowait
            for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                    const double* f = arg.getSampleDataRO(INDEX2(k0, k1, m_NE[0]));
                    for (index_t i=0; i < numComp; ++i) {
                        int_local[i]+=f[i]*w;
                    }
                }
            }
#pragma omp critical
            for (index_t i=0; i<numComp; i++)
                integrals[i]+=int_local[i];
        } // end of parallel section

    } else if (fs == FaceElements && arg.actsExpanded()) {
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
            const double w0 = m_dx[0]/2.;
            const double w1 = m_dx[1]/2.;
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    const double* f = arg.getSampleDataRO(m_faceOffset[0]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        const double f0 = f[INDEX2(i,0,numComp)];
                        const double f1 = f[INDEX2(i,1,numComp)];
                        int_local[i]+=(f0+f1)*w1;
                    }  // end of component loop i
                } // end of k1 loop
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    const double* f = arg.getSampleDataRO(m_faceOffset[1]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        const double f0 = f[INDEX2(i,0,numComp)];
                        const double f1 = f[INDEX2(i,1,numComp)];
                        int_local[i]+=(f0+f1)*w1;
                    }  // end of component loop i
                } // end of k1 loop
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                    const double* f = arg.getSampleDataRO(m_faceOffset[2]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        const double f0 = f[INDEX2(i,0,numComp)];
                        const double f1 = f[INDEX2(i,1,numComp)];
                        int_local[i]+=(f0+f1)*w0;
                    }  // end of component loop i
                } // end of k0 loop
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                    const double* f = arg.getSampleDataRO(m_faceOffset[3]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        const double f0 = f[INDEX2(i,0,numComp)];
                        const double f1 = f[INDEX2(i,1,numComp)];
                        int_local[i]+=(f0+f1)*w0;
                    }  // end of component loop i
                } // end of k0 loop
            }
#pragma omp critical
            for (index_t i=0; i<numComp; i++)
                integrals[i]+=int_local[i];
        } // end of parallel section

    } else if (fs==ReducedFaceElements || (fs==FaceElements && !arg.actsExpanded())) {
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    const double* f = arg.getSampleDataRO(m_faceOffset[0]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        int_local[i]+=f[i]*m_dx[1];
                    }
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    const double* f = arg.getSampleDataRO(m_faceOffset[1]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        int_local[i]+=f[i]*m_dx[1];
                    }
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                    const double* f = arg.getSampleDataRO(m_faceOffset[2]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        int_local[i]+=f[i]*m_dx[0];
                    }
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                    const double* f = arg.getSampleDataRO(m_faceOffset[3]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        int_local[i]+=f[i]*m_dx[0];
                    }
                }
            }

#pragma omp critical
            for (index_t i=0; i<numComp; i++)
                integrals[i]+=int_local[i];
        } // end of parallel section
    } // function space selector
}

//protected
dim_t Rectangle::insertNeighbourNodes(IndexVector& index, index_t node) const
{
    const dim_t nDOF0 = (m_gNE[0]+1)/m_NX[0];
    const dim_t nDOF1 = (m_gNE[1]+1)/m_NX[1];
    const int x=node%nDOF0;
    const int y=node/nDOF0;
    dim_t num=0;
    // loop through potential neighbours and add to index if positions are
    // within bounds
    for (int i1=-1; i1<2; i1++) {
        for (int i0=-1; i0<2; i0++) {
            // skip node itself
            if (i0==0 && i1==0)
                continue;
            // location of neighbour node
            const int nx=x+i0;
            const int ny=y+i1;
            if (nx>=0 && ny>=0 && nx<nDOF0 && ny<nDOF1) {
                index.push_back(ny*nDOF0+nx);
                num++;
            }
        }
    }

    return num;
}

//protected
void Rectangle::nodesToDOF(escript::Data& out, escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    out.requireWrite();

    const index_t left = (m_offset[0]==0 ? 0 : 1);
    const index_t bottom = (m_offset[1]==0 ? 0 : 1);
    const dim_t nDOF0 = (m_gNE[0]+1)/m_NX[0];
    const dim_t nDOF1 = (m_gNE[1]+1)/m_NX[1];
#pragma omp parallel for
    for (index_t i=0; i<nDOF1; i++) {
        for (index_t j=0; j<nDOF0; j++) {
            const index_t n=j+left+(i+bottom)*m_NN[0];
            const double* src=in.getSampleDataRO(n);
            copy(src, src+numComp, out.getSampleDataRW(j+i*nDOF0));
        }
    }
}

//protected
void Rectangle::dofToNodes(escript::Data& out, escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    Paso_Coupler* coupler = Paso_Coupler_alloc(m_connector, numComp);
    in.requireWrite();
    Paso_Coupler_startCollect(coupler, in.getSampleDataRW(0));

    const dim_t numDOF = getNumDOF();
    out.requireWrite();
    const double* buffer = Paso_Coupler_finishCollect(coupler);

#pragma omp parallel for
    for (index_t i=0; i<getNumNodes(); i++) {
        const double* src=(m_dofMap[i]<numDOF ?
                in.getSampleDataRO(m_dofMap[i])
                : &buffer[(m_dofMap[i]-numDOF)*numComp]);
        copy(src, src+numComp, out.getSampleDataRW(i));
    }
    Paso_Coupler_free(coupler);
}

//private
void Rectangle::populateSampleIds()
{
    // degrees of freedom are numbered from left to right, bottom to top in
    // each rank, continuing on the next rank (ranks also go left-right,
    // bottom-top).
    // This means rank 0 has id 0...n0-1, rank 1 has id n0...n1-1 etc. which
    // helps when writing out data rank after rank.

    // build node distribution vector first.
    // rank i owns m_nodeDistribution[i+1]-nodeDistribution[i] nodes which is
    // constant for all ranks in this implementation
    m_nodeDistribution.assign(m_mpiInfo->size+1, 0);
    const dim_t numDOF=getNumDOF();
    for (dim_t k=1; k<m_mpiInfo->size; k++) {
        m_nodeDistribution[k]=k*numDOF;
    }
    m_nodeDistribution[m_mpiInfo->size]=getNumDataPointsGlobal();
    m_nodeId.resize(getNumNodes());
    m_dofId.resize(numDOF);
    m_elementId.resize(getNumElements());

    // populate face element counts
    //left
    if (m_offset[0]==0)
        m_faceCount[0]=m_NE[1];
    else
        m_faceCount[0]=0;
    //right
    if (m_mpiInfo->rank%m_NX[0]==m_NX[0]-1)
        m_faceCount[1]=m_NE[1];
    else
        m_faceCount[1]=0;
    //bottom
    if (m_offset[1]==0)
        m_faceCount[2]=m_NE[0];
    else
        m_faceCount[2]=0;
    //top
    if (m_mpiInfo->rank/m_NX[0]==m_NX[1]-1)
        m_faceCount[3]=m_NE[0];
    else
        m_faceCount[3]=0;

    m_faceId.resize(getNumFaceElements());

    const index_t left = (m_offset[0]==0 ? 0 : 1);
    const index_t bottom = (m_offset[1]==0 ? 0 : 1);
    const dim_t nDOF0 = (m_gNE[0]+1)/m_NX[0];
    const dim_t nDOF1 = (m_gNE[1]+1)/m_NX[1];

#define globalNodeId(x,y) \
    ((m_offset[0]+x)/nDOF0)*nDOF0*nDOF1+(m_offset[0]+x)%nDOF0 \
    + ((m_offset[1]+y)/nDOF1)*nDOF0*nDOF1*m_NX[0]+((m_offset[1]+y)%nDOF1)*nDOF0

    // set corner id's outside the parallel region
    m_nodeId[0] = globalNodeId(0, 0);
    m_nodeId[m_NN[0]-1] = globalNodeId(m_NN[0]-1, 0);
    m_nodeId[m_NN[0]*(m_NN[1]-1)] = globalNodeId(0, m_NN[1]-1);
    m_nodeId[m_NN[0]*m_NN[1]-1] = globalNodeId(m_NN[0]-1,m_NN[1]-1);
#undef globalNodeId

#pragma omp parallel
    {
        // populate degrees of freedom and own nodes (identical id)
#pragma omp for nowait
        for (dim_t i=0; i<nDOF1; i++) {
            for (dim_t j=0; j<nDOF0; j++) {
                const index_t nodeIdx=j+left+(i+bottom)*m_NN[0];
                const index_t dofIdx=j+i*nDOF0;
                m_dofId[dofIdx] = m_nodeId[nodeIdx]
                    = m_nodeDistribution[m_mpiInfo->rank]+dofIdx;
            }
        }

        // populate the rest of the nodes (shared with other ranks)
        if (m_faceCount[0]==0) { // left column
#pragma omp for nowait
            for (dim_t i=0; i<nDOF1; i++) {
                const index_t nodeIdx=(i+bottom)*m_NN[0];
                const index_t dofId=(i+1)*nDOF0-1;
                m_nodeId[nodeIdx]
                    = m_nodeDistribution[m_mpiInfo->rank-1]+dofId;
            }
        }
        if (m_faceCount[1]==0) { // right column
#pragma omp for nowait
            for (dim_t i=0; i<nDOF1; i++) {
                const index_t nodeIdx=(i+bottom+1)*m_NN[0]-1;
                const index_t dofId=i*nDOF0;
                m_nodeId[nodeIdx]
                    = m_nodeDistribution[m_mpiInfo->rank+1]+dofId;
            }
        }
        if (m_faceCount[2]==0) { // bottom row
#pragma omp for nowait
            for (dim_t i=0; i<nDOF0; i++) {
                const index_t nodeIdx=i+left;
                const index_t dofId=nDOF0*(nDOF1-1)+i;
                m_nodeId[nodeIdx]
                    = m_nodeDistribution[m_mpiInfo->rank-m_NX[0]]+dofId;
            }
        }
        if (m_faceCount[3]==0) { // top row
#pragma omp for nowait
            for (dim_t i=0; i<nDOF0; i++) {
                const index_t nodeIdx=m_NN[0]*(m_NN[1]-1)+i+left;
                const index_t dofId=i;
                m_nodeId[nodeIdx]
                    = m_nodeDistribution[m_mpiInfo->rank+m_NX[0]]+dofId;
            }
        }

        // populate element id's
#pragma omp for nowait
        for (dim_t i1=0; i1<m_NE[1]; i1++) {
            for (dim_t i0=0; i0<m_NE[0]; i0++) {
                m_elementId[i0+i1*m_NE[0]]=(m_offset[1]+i1)*m_gNE[0]+m_offset[0]+i0;
            }
        }

        // face elements
#pragma omp for
        for (dim_t k=0; k<getNumFaceElements(); k++)
            m_faceId[k]=k;
    } // end parallel section

    m_nodeTags.assign(getNumNodes(), 0);
    updateTagsInUse(Nodes);

    m_elementTags.assign(getNumElements(), 0);
    updateTagsInUse(Elements);

    // generate face offset vector and set face tags
    const index_t LEFT=1, RIGHT=2, BOTTOM=10, TOP=20;
    const index_t faceTag[] = { LEFT, RIGHT, BOTTOM, TOP };
    m_faceOffset.assign(4, -1);
    m_faceTags.clear();
    index_t offset=0;
    for (size_t i=0; i<4; i++) {
        if (m_faceCount[i]>0) {
            m_faceOffset[i]=offset;
            offset+=m_faceCount[i];
            m_faceTags.insert(m_faceTags.end(), m_faceCount[i], faceTag[i]);
        }
    }
    setTagMap("left", LEFT);
    setTagMap("right", RIGHT);
    setTagMap("bottom", BOTTOM);
    setTagMap("top", TOP);
    updateTagsInUse(FaceElements);
}

//private
void Rectangle::createPattern()
{
    const dim_t nDOF0 = (m_gNE[0]+1)/m_NX[0];
    const dim_t nDOF1 = (m_gNE[1]+1)/m_NX[1];
    const index_t left = (m_offset[0]==0 ? 0 : 1);
    const index_t bottom = (m_offset[1]==0 ? 0 : 1);

    // populate node->DOF mapping with own degrees of freedom.
    // The rest is assigned in the loop further down
    m_dofMap.assign(getNumNodes(), 0);
#pragma omp parallel for
    for (index_t i=bottom; i<bottom+nDOF1; i++) {
        for (index_t j=left; j<left+nDOF0; j++) {
            m_dofMap[i*m_NN[0]+j]=(i-bottom)*nDOF0+j-left;
        }
    }

    // build list of shared components and neighbours by looping through
    // all potential neighbouring ranks and checking if positions are
    // within bounds
    const dim_t numDOF=nDOF0*nDOF1;
    vector<IndexVector> colIndices(numDOF); // for the couple blocks
    RankVector neighbour;
    IndexVector offsetInShared(1,0);
    IndexVector sendShared, recvShared;
    int numShared=0;
    const int x=m_mpiInfo->rank%m_NX[0];
    const int y=m_mpiInfo->rank/m_NX[0];
    for (int i1=-1; i1<2; i1++) {
        for (int i0=-1; i0<2; i0++) {
            // skip this rank
            if (i0==0 && i1==0)
                continue;
            // location of neighbour rank
            const int nx=x+i0;
            const int ny=y+i1;
            if (nx>=0 && ny>=0 && nx<m_NX[0] && ny<m_NX[1]) {
                neighbour.push_back(ny*m_NX[0]+nx);
                if (i0==0) {
                    // sharing top or bottom edge
                    const int firstDOF=(i1==-1 ? 0 : numDOF-nDOF0);
                    const int firstNode=(i1==-1 ? left : m_NN[0]*(m_NN[1]-1)+left);
                    offsetInShared.push_back(offsetInShared.back()+nDOF0);
                    for (dim_t i=0; i<nDOF0; i++, numShared++) {
                        sendShared.push_back(firstDOF+i);
                        recvShared.push_back(numDOF+numShared);
                        if (i>0)
                            colIndices[firstDOF+i-1].push_back(numShared);
                        colIndices[firstDOF+i].push_back(numShared);
                        if (i<nDOF0-1)
                            colIndices[firstDOF+i+1].push_back(numShared);
                        m_dofMap[firstNode+i]=numDOF+numShared;
                    }
                } else if (i1==0) {
                    // sharing left or right edge
                    const int firstDOF=(i0==-1 ? 0 : nDOF0-1);
                    const int firstNode=(i0==-1 ? bottom*m_NN[0] : (bottom+1)*m_NN[0]-1);
                    offsetInShared.push_back(offsetInShared.back()+nDOF1);
                    for (dim_t i=0; i<nDOF1; i++, numShared++) {
                        sendShared.push_back(firstDOF+i*nDOF0);
                        recvShared.push_back(numDOF+numShared);
                        if (i>0)
                            colIndices[firstDOF+(i-1)*nDOF0].push_back(numShared);
                        colIndices[firstDOF+i*nDOF0].push_back(numShared);
                        if (i<nDOF1-1)
                            colIndices[firstDOF+(i+1)*nDOF0].push_back(numShared);
                        m_dofMap[firstNode+i*m_NN[0]]=numDOF+numShared;
                    }
                } else {
                    // sharing a node
                    const int dof=(i0+1)/2*(nDOF0-1)+(i1+1)/2*(numDOF-nDOF0);
                    const int node=(i0+1)/2*(m_NN[0]-1)+(i1+1)/2*m_NN[0]*(m_NN[1]-1);
                    offsetInShared.push_back(offsetInShared.back()+1);
                    sendShared.push_back(dof);
                    recvShared.push_back(numDOF+numShared);
                    colIndices[dof].push_back(numShared);
                    m_dofMap[node]=numDOF+numShared;
                    ++numShared;
                }
            }
        }
    }

    // create connector
    Paso_SharedComponents *snd_shcomp = Paso_SharedComponents_alloc(
            numDOF, neighbour.size(), &neighbour[0], &sendShared[0],
            &offsetInShared[0], 1, 0, m_mpiInfo);
    Paso_SharedComponents *rcv_shcomp = Paso_SharedComponents_alloc(
            numDOF, neighbour.size(), &neighbour[0], &recvShared[0],
            &offsetInShared[0], 1, 0, m_mpiInfo);
    m_connector = Paso_Connector_alloc(snd_shcomp, rcv_shcomp);
    Paso_SharedComponents_free(snd_shcomp);
    Paso_SharedComponents_free(rcv_shcomp);

    // create main and couple blocks
    Paso_Pattern *mainPattern = createMainPattern();
    Paso_Pattern *colPattern, *rowPattern;
    createCouplePatterns(colIndices, numShared, &colPattern, &rowPattern);

    // allocate paso distribution
    Paso_Distribution* distribution = Paso_Distribution_alloc(m_mpiInfo,
            const_cast<index_t*>(&m_nodeDistribution[0]), 1, 0);

    // finally create the system matrix
    m_pattern = Paso_SystemMatrixPattern_alloc(MATRIX_FORMAT_DEFAULT,
            distribution, distribution, mainPattern, colPattern, rowPattern,
            m_connector, m_connector);

    Paso_Distribution_free(distribution);

    // useful debug output
    /*
    cout << "--- rcv_shcomp ---" << endl;
    cout << "numDOF=" << numDOF << ", numNeighbors=" << neighbour.size() << endl;
    for (size_t i=0; i<neighbour.size(); i++) {
        cout << "neighbor[" << i << "]=" << neighbour[i]
            << " offsetInShared[" << i+1 << "]=" << offsetInShared[i+1] << endl;
    }
    for (size_t i=0; i<recvShared.size(); i++) {
        cout << "shared[" << i << "]=" << recvShared[i] << endl;
    }
    cout << "--- snd_shcomp ---" << endl;
    for (size_t i=0; i<sendShared.size(); i++) {
        cout << "shared[" << i << "]=" << sendShared[i] << endl;
    }
    cout << "--- dofMap ---" << endl;
    for (size_t i=0; i<m_dofMap.size(); i++) {
        cout << "m_dofMap[" << i << "]=" << m_dofMap[i] << endl;
    }
    cout << "--- colIndices ---" << endl;
    for (size_t i=0; i<colIndices.size(); i++) {
        cout << "colIndices[" << i << "].size()=" << colIndices[i].size() << endl;
    }
    */

    /*
    cout << "--- main_pattern ---" << endl;
    cout << "M=" << mainPattern->numOutput << ", N=" << mainPattern->numInput << endl;
    for (size_t i=0; i<mainPattern->numOutput+1; i++) {
        cout << "ptr[" << i << "]=" << mainPattern->ptr[i] << endl;
    }
    for (size_t i=0; i<mainPattern->ptr[mainPattern->numOutput]; i++) {
        cout << "index[" << i << "]=" << mainPattern->index[i] << endl;
    }
    */

    /*
    cout << "--- colCouple_pattern ---" << endl;
    cout << "M=" << colPattern->numOutput << ", N=" << colPattern->numInput << endl;
    for (size_t i=0; i<colPattern->numOutput+1; i++) {
        cout << "ptr[" << i << "]=" << colPattern->ptr[i] << endl;
    }
    for (size_t i=0; i<colPattern->ptr[colPattern->numOutput]; i++) {
        cout << "index[" << i << "]=" << colPattern->index[i] << endl;
    }
    */

    /*
    cout << "--- rowCouple_pattern ---" << endl;
    cout << "M=" << rowPattern->numOutput << ", N=" << rowPattern->numInput << endl;
    for (size_t i=0; i<rowPattern->numOutput+1; i++) {
        cout << "ptr[" << i << "]=" << rowPattern->ptr[i] << endl;
    }
    for (size_t i=0; i<rowPattern->ptr[rowPattern->numOutput]; i++) {
        cout << "index[" << i << "]=" << rowPattern->index[i] << endl;
    }
    */

    Paso_Pattern_free(mainPattern);
    Paso_Pattern_free(colPattern);
    Paso_Pattern_free(rowPattern);
}

//private
void Rectangle::addToMatrixAndRHS(Paso_SystemMatrix* S, escript::Data& F,
         const vector<double>& EM_S, const vector<double>& EM_F, bool addS,
         bool addF, index_t firstNode, dim_t nEq, dim_t nComp) const
{
    IndexVector rowIndex;
    rowIndex.push_back(m_dofMap[firstNode]);
    rowIndex.push_back(m_dofMap[firstNode+1]);
    rowIndex.push_back(m_dofMap[firstNode+m_NN[0]]);
    rowIndex.push_back(m_dofMap[firstNode+m_NN[0]+1]);
    if (addF) {
        double *F_p=F.getSampleDataRW(0);
        for (index_t i=0; i<rowIndex.size(); i++) {
            if (rowIndex[i]<getNumDOF()) {
                for (index_t eq=0; eq<nEq; eq++) {
                    F_p[INDEX2(eq, rowIndex[i], nEq)]+=EM_F[INDEX2(eq,i,nEq)];
                }
            }
        }
    }
    if (addS) {
        addToSystemMatrix(S, rowIndex, nEq, rowIndex, nComp, EM_S);
    }
}

//protected
void Rectangle::interpolateNodesOnElements(escript::Data& out,
                                        escript::Data& in, bool reduced) const
{
    const dim_t numComp = in.getDataPointSize();
    if (reduced) {
        out.requireWrite();
        const double c0 = 0.25;
#pragma omp parallel
        {
            vector<double> f_00(numComp);
            vector<double> f_01(numComp);
            vector<double> f_10(numComp);
            vector<double> f_11(numComp);
#pragma omp for
            for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,k1+1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,k1+1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(INDEX2(k0,k1,m_NE[0]));
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = c0*(f_00[i] + f_01[i] + f_10[i] + f_11[i]);
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of k1 loop */
        } /* end of parallel section */
    } else {
        out.requireWrite();
        const double c0 = 0.16666666666666666667;
        const double c1 = 0.044658198738520451079;
        const double c2 = 0.62200846792814621559;
#pragma omp parallel
        {
            vector<double> f_00(numComp);
            vector<double> f_01(numComp);
            vector<double> f_10(numComp);
            vector<double> f_11(numComp);
#pragma omp for
            for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,k1+1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,k1+1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(INDEX2(k0,k1,m_NE[0]));
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = c0*(f_01[i] + f_10[i]) + c1*f_11[i] + c2*f_00[i];
                        o[INDEX2(i,numComp,1)] = c0*(f_00[i] + f_11[i]) + c1*f_01[i] + c2*f_10[i];
                        o[INDEX2(i,numComp,2)] = c0*(f_00[i] + f_11[i]) + c1*f_10[i] + c2*f_01[i];
                        o[INDEX2(i,numComp,3)] = c0*(f_01[i] + f_10[i]) + c1*f_00[i] + c2*f_11[i];
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of k1 loop */
        } /* end of parallel section */
    }
}

//protected
void Rectangle::interpolateNodesOnFaces(escript::Data& out, escript::Data& in,
                                        bool reduced) const
{
    const dim_t numComp = in.getDataPointSize();
    if (reduced) {
        out.requireWrite();
#pragma omp parallel
        {
            vector<double> f_00(numComp);
            vector<double> f_01(numComp);
            vector<double> f_10(numComp);
            vector<double> f_11(numComp);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(0,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(0,k1+1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = (f_00[i] + f_01[i])/2;
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of face 0 */
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1+1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = (f_10[i] + f_11[i])/2;
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of face 1 */
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,0, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,0, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = (f_00[i] + f_10[i])/2;
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of face 2 */
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,m_NN[1]-1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,m_NN[1]-1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = (f_01[i] + f_11[i])/2;
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of face 3 */
        } /* end of parallel section */
    } else {
        out.requireWrite();
        const double c0 = 0.21132486540518711775;
        const double c1 = 0.78867513459481288225;
#pragma omp parallel
        {
            vector<double> f_00(numComp);
            vector<double> f_01(numComp);
            vector<double> f_10(numComp);
            vector<double> f_11(numComp);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(0,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(0,k1+1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = c0*f_01[i] + c1*f_00[i];
                        o[INDEX2(i,numComp,1)] = c0*f_00[i] + c1*f_01[i];
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of face 0 */
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1+1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = c1*f_10[i] + c0*f_11[i];
                        o[INDEX2(i,numComp,1)] = c1*f_11[i] + c0*f_10[i];
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of face 1 */
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,0, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,0, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = c0*f_10[i] + c1*f_00[i];
                        o[INDEX2(i,numComp,1)] = c0*f_00[i] + c1*f_10[i];
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of face 2 */
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,m_NN[1]-1, m_NN[0])), numComp*sizeof(double));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,m_NN[1]-1, m_NN[0])), numComp*sizeof(double));
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = c0*f_11[i] + c1*f_01[i];
                        o[INDEX2(i,numComp,1)] = c0*f_01[i] + c1*f_11[i];
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of face 3 */
        } /* end of parallel section */
    }
}

//protected
void Rectangle::assemblePDESingle(Paso_SystemMatrix* mat,
        escript::Data& rhs, const escript::Data& A, const escript::Data& B,
        const escript::Data& C, const escript::Data& D,
        const escript::Data& X, const escript::Data& Y) const
{
    const double SQRT3 = 1.73205080756887719318;
    const double w1 = 1.0/24.0;
    const double w5 = -SQRT3/24 + 1.0/12;
    const double w2 = -SQRT3/24 - 1.0/12;
    const double w19 = -m_dx[0]/12;
    const double w11 = w19*(SQRT3 + 3)/12;
    const double w14 = w19*(-SQRT3 + 3)/12;
    const double w16 = w19*(5*SQRT3 + 9)/12;
    const double w17 = w19*(-5*SQRT3 + 9)/12;
    const double w27 = w19*(-SQRT3 - 3)/2;
    const double w28 = w19*(SQRT3 - 3)/2;
    const double w18 = -m_dx[1]/12;
    const double w12 = w18*(5*SQRT3 + 9)/12;
    const double w13 = w18*(-5*SQRT3 + 9)/12;
    const double w10 = w18*(SQRT3 + 3)/12;
    const double w15 = w18*(-SQRT3 + 3)/12;
    const double w25 = w18*(-SQRT3 - 3)/2;
    const double w26 = w18*(SQRT3 - 3)/2;
    const double w22 = m_dx[0]*m_dx[1]/144;
    const double w20 = w22*(SQRT3 + 2);
    const double w21 = w22*(-SQRT3 + 2);
    const double w23 = w22*(4*SQRT3 + 7);
    const double w24 = w22*(-4*SQRT3 + 7);
    const double w3 = m_dx[0]/(24*m_dx[1]);
    const double w7 = w3*(SQRT3 + 2);
    const double w8 = w3*(-SQRT3 + 2);
    const double w6 = -m_dx[1]/(24*m_dx[0]);
    const double w0 = w6*(SQRT3 + 2);
    const double w4 = w6*(-SQRT3 + 2);

    rhs.requireWrite();
#pragma omp parallel
    {
        for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
            for (index_t k1=k1_0; k1<m_NE[1]; k1+=2) {
                for (index_t k0=0; k0<m_NE[0]; ++k0)  {
                    bool add_EM_S=false;
                    bool add_EM_F=false;
                    vector<double> EM_S(4*4, 0);
                    vector<double> EM_F(4, 0);
                    const index_t e = k0 + m_NE[0]*k1;
                    ///////////////
                    // process A //
                    ///////////////
                    if (!A.isEmpty()) {
                        add_EM_S = true;
                        const double* A_p = const_cast<escript::Data*>(&A)->getSampleDataRO(e);
                        if (A.actsExpanded()) {
                            const double A_00_0 = A_p[INDEX3(0,0,0,2,2)];
                            const double A_01_0 = A_p[INDEX3(0,1,0,2,2)];
                            const double A_10_0 = A_p[INDEX3(1,0,0,2,2)];
                            const double A_11_0 = A_p[INDEX3(1,1,0,2,2)];
                            const double A_00_1 = A_p[INDEX3(0,0,1,2,2)];
                            const double A_01_1 = A_p[INDEX3(0,1,1,2,2)];
                            const double A_10_1 = A_p[INDEX3(1,0,1,2,2)];
                            const double A_11_1 = A_p[INDEX3(1,1,1,2,2)];
                            const double A_00_2 = A_p[INDEX3(0,0,2,2,2)];
                            const double A_01_2 = A_p[INDEX3(0,1,2,2,2)];
                            const double A_10_2 = A_p[INDEX3(1,0,2,2,2)];
                            const double A_11_2 = A_p[INDEX3(1,1,2,2,2)];
                            const double A_00_3 = A_p[INDEX3(0,0,3,2,2)];
                            const double A_01_3 = A_p[INDEX3(0,1,3,2,2)];
                            const double A_10_3 = A_p[INDEX3(1,0,3,2,2)];
                            const double A_11_3 = A_p[INDEX3(1,1,3,2,2)];
                            const double tmp0 = w3*(A_11_0 + A_11_1 + A_11_2 + A_11_3);
                            const double tmp1 = w1*(A_01_0 + A_01_3 - A_10_1 - A_10_2);
                            const double tmp2 = w4*(A_00_2 + A_00_3);
                            const double tmp3 = w0*(A_00_0 + A_00_1);
                            const double tmp4 = w5*(A_01_2 - A_10_3);
                            const double tmp5 = w2*(-A_01_1 + A_10_0);
                            const double tmp6 = w5*(A_01_3 + A_10_0);
                            const double tmp7 = w3*(-A_11_0 - A_11_1 - A_11_2 - A_11_3);
                            const double tmp8 = w6*(A_00_0 + A_00_1 + A_00_2 + A_00_3);
                            const double tmp9 = w1*(A_01_1 + A_01_2 + A_10_1 + A_10_2);
                            const double tmp10 = w2*(-A_01_0 - A_10_3);
                            const double tmp11 = w4*(A_00_0 + A_00_1);
                            const double tmp12 = w0*(A_00_2 + A_00_3);
                            const double tmp13 = w5*(A_01_1 - A_10_0);
                            const double tmp14 = w2*(-A_01_2 + A_10_3);
                            const double tmp15 = w7*(A_11_0 + A_11_2);
                            const double tmp16 = w4*(-A_00_2 - A_00_3);
                            const double tmp17 = w0*(-A_00_0 - A_00_1);
                            const double tmp18 = w5*(A_01_3 + A_10_3);
                            const double tmp19 = w8*(A_11_1 + A_11_3);
                            const double tmp20 = w2*(-A_01_0 - A_10_0);
                            const double tmp21 = w7*(A_11_1 + A_11_3);
                            const double tmp22 = w4*(-A_00_0 - A_00_1);
                            const double tmp23 = w0*(-A_00_2 - A_00_3);
                            const double tmp24 = w5*(A_01_0 + A_10_0);
                            const double tmp25 = w8*(A_11_0 + A_11_2);
                            const double tmp26 = w2*(-A_01_3 - A_10_3);
                            const double tmp27 = w5*(-A_01_1 - A_10_2);
                            const double tmp28 = w1*(-A_01_0 - A_01_3 - A_10_0 - A_10_3);
                            const double tmp29 = w2*(A_01_2 + A_10_1);
                            const double tmp30 = w7*(-A_11_1 - A_11_3);
                            const double tmp31 = w1*(-A_01_1 - A_01_2 + A_10_0 + A_10_3);
                            const double tmp32 = w5*(-A_01_0 + A_10_2);
                            const double tmp33 = w8*(-A_11_0 - A_11_2);
                            const double tmp34 = w6*(-A_00_0 - A_00_1 - A_00_2 - A_00_3);
                            const double tmp35 = w2*(A_01_3 - A_10_1);
                            const double tmp36 = w5*(A_01_0 + A_10_3);
                            const double tmp37 = w2*(-A_01_3 - A_10_0);
                            const double tmp38 = w7*(-A_11_0 - A_11_2);
                            const double tmp39 = w5*(-A_01_3 + A_10_1);
                            const double tmp40 = w8*(-A_11_1 - A_11_3);
                            const double tmp41 = w2*(A_01_0 - A_10_2);
                            const double tmp42 = w5*(A_01_1 - A_10_3);
                            const double tmp43 = w2*(-A_01_2 + A_10_0);
                            const double tmp44 = w5*(A_01_2 - A_10_0);
                            const double tmp45 = w2*(-A_01_1 + A_10_3);
                            const double tmp46 = w5*(-A_01_0 + A_10_1);
                            const double tmp47 = w2*(A_01_3 - A_10_2);
                            const double tmp48 = w5*(-A_01_1 - A_10_1);
                            const double tmp49 = w2*(A_01_2 + A_10_2);
                            const double tmp50 = w5*(-A_01_3 + A_10_2);
                            const double tmp51 = w2*(A_01_0 - A_10_1);
                            const double tmp52 = w5*(-A_01_2 - A_10_1);
                            const double tmp53 = w2*(A_01_1 + A_10_2);
                            const double tmp54 = w5*(-A_01_2 - A_10_2);
                            const double tmp55 = w2*(A_01_1 + A_10_1);
                            EM_S[INDEX2(0,0,4)]+=tmp15 + tmp16 + tmp17 + tmp18 + tmp19 + tmp20 + tmp9;
                            EM_S[INDEX2(0,1,4)]+=tmp0 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5;
                            EM_S[INDEX2(0,2,4)]+=tmp31 + tmp34 + tmp38 + tmp39 + tmp40 + tmp41;
                            EM_S[INDEX2(0,3,4)]+=tmp28 + tmp52 + tmp53 + tmp7 + tmp8;
                            EM_S[INDEX2(1,0,4)]+=tmp0 + tmp2 + tmp3 + tmp31 + tmp50 + tmp51;
                            EM_S[INDEX2(1,1,4)]+=tmp16 + tmp17 + tmp21 + tmp25 + tmp28 + tmp54 + tmp55;
                            EM_S[INDEX2(1,2,4)]+=tmp10 + tmp6 + tmp7 + tmp8 + tmp9;
                            EM_S[INDEX2(1,3,4)]+=tmp1 + tmp30 + tmp33 + tmp34 + tmp44 + tmp45;
                            EM_S[INDEX2(2,0,4)]+=tmp1 + tmp34 + tmp38 + tmp40 + tmp42 + tmp43;
                            EM_S[INDEX2(2,1,4)]+=tmp36 + tmp37 + tmp7 + tmp8 + tmp9;
                            EM_S[INDEX2(2,2,4)]+=tmp15 + tmp19 + tmp22 + tmp23 + tmp28 + tmp48 + tmp49;
                            EM_S[INDEX2(2,3,4)]+=tmp0 + tmp11 + tmp12 + tmp31 + tmp46 + tmp47;
                            EM_S[INDEX2(3,0,4)]+=tmp27 + tmp28 + tmp29 + tmp7 + tmp8;
                            EM_S[INDEX2(3,1,4)]+=tmp30 + tmp31 + tmp32 + tmp33 + tmp34 + tmp35;
                            EM_S[INDEX2(3,2,4)]+=tmp0 + tmp1 + tmp11 + tmp12 + tmp13 + tmp14;
                            EM_S[INDEX2(3,3,4)]+=tmp21 + tmp22 + tmp23 + tmp24 + tmp25 + tmp26 + tmp9;
                        } else { // constant data
                            const double A_00 = A_p[INDEX2(0,0,2)];
                            const double A_01 = A_p[INDEX2(0,1,2)];
                            const double A_10 = A_p[INDEX2(1,0,2)];
                            const double A_11 = A_p[INDEX2(1,1,2)];
                            const double tmp0 = 6*w1*(A_01 - A_10);
                            const double tmp1 = 6*w1*(A_01 + A_10);
                            const double tmp2 = 6*w1*(-A_01 - A_10);
                            const double tmp3 = 6*w1*(-A_01 + A_10);
                            EM_S[INDEX2(0,0,4)]+=-8*A_00*w6 + 8*A_11*w3 + tmp1;
                            EM_S[INDEX2(0,1,4)]+=8*A_00*w6 + 4*A_11*w3 + tmp0;
                            EM_S[INDEX2(0,2,4)]+=-4*A_00*w6 - 8*A_11*w3 + tmp3;
                            EM_S[INDEX2(0,3,4)]+=4*A_00*w6 - 4*A_11*w3 + tmp2;
                            EM_S[INDEX2(1,0,4)]+=8*A_00*w6 + 4*A_11*w3 + tmp3;
                            EM_S[INDEX2(1,1,4)]+=-8*A_00*w6 + 8*A_11*w3 + tmp2;
                            EM_S[INDEX2(1,2,4)]+=4*A_00*w6 - 4*A_11*w3 + tmp1;
                            EM_S[INDEX2(1,3,4)]+=-4*A_00*w6 - 8*A_11*w3 + tmp0;
                            EM_S[INDEX2(2,0,4)]+=-4*A_00*w6 - 8*A_11*w3 + tmp0;
                            EM_S[INDEX2(2,1,4)]+=4*A_00*w6 - 4*A_11*w3 + tmp1;
                            EM_S[INDEX2(2,2,4)]+=-8*A_00*w6 + 8*A_11*w3 + tmp2;
                            EM_S[INDEX2(2,3,4)]+=8*A_00*w6 + 4*A_11*w3 + tmp3;
                            EM_S[INDEX2(3,0,4)]+=4*A_00*w6 - 4*A_11*w3 + tmp2;
                            EM_S[INDEX2(3,1,4)]+=-4*A_00*w6 - 8*A_11*w3 + tmp3;
                            EM_S[INDEX2(3,2,4)]+=8*A_00*w6 + 4*A_11*w3 + tmp0;
                            EM_S[INDEX2(3,3,4)]+=-8*A_00*w6 + 8*A_11*w3 + tmp1;
                        }
                    }
                    ///////////////
                    // process B //
                    ///////////////
                    if (!B.isEmpty()) {
                        add_EM_S=true;
                        const double* B_p=const_cast<escript::Data*>(&B)->getSampleDataRO(e);
                        if (B.actsExpanded()) {
                            const double B_0_0 = B_p[INDEX2(0,0,2)];
                            const double B_1_0 = B_p[INDEX2(1,0,2)];
                            const double B_0_1 = B_p[INDEX2(0,1,2)];
                            const double B_1_1 = B_p[INDEX2(1,1,2)];
                            const double B_0_2 = B_p[INDEX2(0,2,2)];
                            const double B_1_2 = B_p[INDEX2(1,2,2)];
                            const double B_0_3 = B_p[INDEX2(0,3,2)];
                            const double B_1_3 = B_p[INDEX2(1,3,2)];
                            const double tmp0 = w11*(B_1_0 + B_1_1);
                            const double tmp1 = w14*(B_1_2 + B_1_3);
                            const double tmp2 = w15*(-B_0_1 - B_0_3);
                            const double tmp3 = w10*(-B_0_0 - B_0_2);
                            const double tmp4 = w11*(B_1_2 + B_1_3);
                            const double tmp5 = w14*(B_1_0 + B_1_1);
                            const double tmp6 = w11*(-B_1_2 - B_1_3);
                            const double tmp7 = w14*(-B_1_0 - B_1_1);
                            const double tmp8 = w11*(-B_1_0 - B_1_1);
                            const double tmp9 = w14*(-B_1_2 - B_1_3);
                            const double tmp10 = w10*(-B_0_1 - B_0_3);
                            const double tmp11 = w15*(-B_0_0 - B_0_2);
                            const double tmp12 = w15*(B_0_0 + B_0_2);
                            const double tmp13 = w10*(B_0_1 + B_0_3);
                            const double tmp14 = w10*(B_0_0 + B_0_2);
                            const double tmp15 = w15*(B_0_1 + B_0_3);
                            EM_S[INDEX2(0,0,4)]+=B_0_0*w12 + B_0_1*w10 + B_0_2*w15 + B_0_3*w13 + B_1_0*w16 + B_1_1*w14 + B_1_2*w11 + B_1_3*w17;
                            EM_S[INDEX2(0,1,4)]+=B_0_0*w10 + B_0_1*w12 + B_0_2*w13 + B_0_3*w15 + tmp0 + tmp1;
                            EM_S[INDEX2(0,2,4)]+=B_1_0*w11 + B_1_1*w17 + B_1_2*w16 + B_1_3*w14 + tmp14 + tmp15;
                            EM_S[INDEX2(0,3,4)]+=tmp12 + tmp13 + tmp4 + tmp5;
                            EM_S[INDEX2(1,0,4)]+=-B_0_0*w12 - B_0_1*w10 - B_0_2*w15 - B_0_3*w13 + tmp0 + tmp1;
                            EM_S[INDEX2(1,1,4)]+=-B_0_0*w10 - B_0_1*w12 - B_0_2*w13 - B_0_3*w15 + B_1_0*w14 + B_1_1*w16 + B_1_2*w17 + B_1_3*w11;
                            EM_S[INDEX2(1,2,4)]+=tmp2 + tmp3 + tmp4 + tmp5;
                            EM_S[INDEX2(1,3,4)]+=B_1_0*w17 + B_1_1*w11 + B_1_2*w14 + B_1_3*w16 + tmp10 + tmp11;
                            EM_S[INDEX2(2,0,4)]+=-B_1_0*w16 - B_1_1*w14 - B_1_2*w11 - B_1_3*w17 + tmp14 + tmp15;
                            EM_S[INDEX2(2,1,4)]+=tmp12 + tmp13 + tmp8 + tmp9;
                            EM_S[INDEX2(2,2,4)]+=B_0_0*w15 + B_0_1*w13 + B_0_2*w12 + B_0_3*w10 - B_1_0*w11 - B_1_1*w17 - B_1_2*w16 - B_1_3*w14;
                            EM_S[INDEX2(2,3,4)]+=B_0_0*w13 + B_0_1*w15 + B_0_2*w10 + B_0_3*w12 + tmp6 + tmp7;
                            EM_S[INDEX2(3,0,4)]+=tmp2 + tmp3 + tmp8 + tmp9;
                            EM_S[INDEX2(3,1,4)]+=-B_1_0*w14 - B_1_1*w16 - B_1_2*w17 - B_1_3*w11 + tmp10 + tmp11;
                            EM_S[INDEX2(3,2,4)]+=-B_0_0*w15 - B_0_1*w13 - B_0_2*w12 - B_0_3*w10 + tmp6 + tmp7;
                            EM_S[INDEX2(3,3,4)]+=-B_0_0*w13 - B_0_1*w15 - B_0_2*w10 - B_0_3*w12 - B_1_0*w17 - B_1_1*w11 - B_1_2*w14 - B_1_3*w16;
                        } else { // constant data
                            const double B_0 = B_p[0];
                            const double B_1 = B_p[1];
                            EM_S[INDEX2(0,0,4)]+= 2*B_0*w18 + 2*B_1*w19;
                            EM_S[INDEX2(0,1,4)]+= 2*B_0*w18 +   B_1*w19;
                            EM_S[INDEX2(0,2,4)]+=   B_0*w18 + 2*B_1*w19;
                            EM_S[INDEX2(0,3,4)]+=   B_0*w18 +   B_1*w19;
                            EM_S[INDEX2(1,0,4)]+=-2*B_0*w18 +   B_1*w19;
                            EM_S[INDEX2(1,1,4)]+=-2*B_0*w18 + 2*B_1*w19;
                            EM_S[INDEX2(1,2,4)]+=  -B_0*w18 +   B_1*w19;
                            EM_S[INDEX2(1,3,4)]+=  -B_0*w18 + 2*B_1*w19;
                            EM_S[INDEX2(2,0,4)]+=   B_0*w18 - 2*B_1*w19;
                            EM_S[INDEX2(2,1,4)]+=   B_0*w18 -   B_1*w19;
                            EM_S[INDEX2(2,2,4)]+= 2*B_0*w18 - 2*B_1*w19;
                            EM_S[INDEX2(2,3,4)]+= 2*B_0*w18 -   B_1*w19;
                            EM_S[INDEX2(3,0,4)]+=  -B_0*w18 -   B_1*w19;
                            EM_S[INDEX2(3,1,4)]+=  -B_0*w18 - 2*B_1*w19;
                            EM_S[INDEX2(3,2,4)]+=-2*B_0*w18 -   B_1*w19;
                            EM_S[INDEX2(3,3,4)]+=-2*B_0*w18 - 2*B_1*w19;
                        }
                    }
                    ///////////////
                    // process C //
                    ///////////////
                    if (!C.isEmpty()) {
                        add_EM_S=true;
                        const double* C_p=const_cast<escript::Data*>(&C)->getSampleDataRO(e);
                        if (C.actsExpanded()) {
                            const double C_0_0 = C_p[INDEX2(0,0,2)];
                            const double C_1_0 = C_p[INDEX2(1,0,2)];
                            const double C_0_1 = C_p[INDEX2(0,1,2)];
                            const double C_1_1 = C_p[INDEX2(1,1,2)];
                            const double C_0_2 = C_p[INDEX2(0,2,2)];
                            const double C_1_2 = C_p[INDEX2(1,2,2)];
                            const double C_0_3 = C_p[INDEX2(0,3,2)];
                            const double C_1_3 = C_p[INDEX2(1,3,2)];
                            const double tmp0 = w11*(C_1_0 + C_1_1);
                            const double tmp1 = w14*(C_1_2 + C_1_3);
                            const double tmp2 = w15*(C_0_0 + C_0_2);
                            const double tmp3 = w10*(C_0_1 + C_0_3);
                            const double tmp4 = w11*(-C_1_0 - C_1_1);
                            const double tmp5 = w14*(-C_1_2 - C_1_3);
                            const double tmp6 = w11*(-C_1_2 - C_1_3);
                            const double tmp7 = w14*(-C_1_0 - C_1_1);
                            const double tmp8 = w11*(C_1_2 + C_1_3);
                            const double tmp9 = w14*(C_1_0 + C_1_1);
                            const double tmp10 = w10*(-C_0_1 - C_0_3);
                            const double tmp11 = w15*(-C_0_0 - C_0_2);
                            const double tmp12 = w15*(-C_0_1 - C_0_3);
                            const double tmp13 = w10*(-C_0_0 - C_0_2);
                            const double tmp14 = w10*(C_0_0 + C_0_2);
                            const double tmp15 = w15*(C_0_1 + C_0_3);
                            EM_S[INDEX2(0,0,4)]+=C_0_0*w12 + C_0_1*w10 + C_0_2*w15 + C_0_3*w13 + C_1_0*w16 + C_1_1*w14 + C_1_2*w11 + C_1_3*w17;
                            EM_S[INDEX2(0,1,4)]+=-C_0_0*w12 - C_0_1*w10 - C_0_2*w15 - C_0_3*w13 + tmp0 + tmp1;
                            EM_S[INDEX2(0,2,4)]+=-C_1_0*w16 - C_1_1*w14 - C_1_2*w11 - C_1_3*w17 + tmp14 + tmp15;
                            EM_S[INDEX2(0,3,4)]+=tmp12 + tmp13 + tmp4 + tmp5;
                            EM_S[INDEX2(1,0,4)]+=C_0_0*w10 + C_0_1*w12 + C_0_2*w13 + C_0_3*w15 + tmp0 + tmp1;
                            EM_S[INDEX2(1,1,4)]+=-C_0_0*w10 - C_0_1*w12 - C_0_2*w13 - C_0_3*w15 + C_1_0*w14 + C_1_1*w16 + C_1_2*w17 + C_1_3*w11;
                            EM_S[INDEX2(1,2,4)]+=tmp2 + tmp3 + tmp4 + tmp5;
                            EM_S[INDEX2(1,3,4)]+=-C_1_0*w14 - C_1_1*w16 - C_1_2*w17 - C_1_3*w11 + tmp10 + tmp11;
                            EM_S[INDEX2(2,0,4)]+=C_1_0*w11 + C_1_1*w17 + C_1_2*w16 + C_1_3*w14 + tmp14 + tmp15;
                            EM_S[INDEX2(2,1,4)]+=tmp12 + tmp13 + tmp8 + tmp9;
                            EM_S[INDEX2(2,2,4)]+=C_0_0*w15 + C_0_1*w13 + C_0_2*w12 + C_0_3*w10 - C_1_0*w11 - C_1_1*w17 - C_1_2*w16 - C_1_3*w14;
                            EM_S[INDEX2(2,3,4)]+=-C_0_0*w15 - C_0_1*w13 - C_0_2*w12 - C_0_3*w10 + tmp6 + tmp7;
                            EM_S[INDEX2(3,0,4)]+=tmp2 + tmp3 + tmp8 + tmp9;
                            EM_S[INDEX2(3,1,4)]+=C_1_0*w17 + C_1_1*w11 + C_1_2*w14 + C_1_3*w16 + tmp10 + tmp11;
                            EM_S[INDEX2(3,2,4)]+=C_0_0*w13 + C_0_1*w15 + C_0_2*w10 + C_0_3*w12 + tmp6 + tmp7;
                            EM_S[INDEX2(3,3,4)]+=-C_0_0*w13 - C_0_1*w15 - C_0_2*w10 - C_0_3*w12 - C_1_0*w17 - C_1_1*w11 - C_1_2*w14 - C_1_3*w16;
                        } else { // constant data
                            const double C_0 = C_p[0];
                            const double C_1 = C_p[1];
                            EM_S[INDEX2(0,0,4)]+= 2*C_0*w18 + 2*C_1*w19;
                            EM_S[INDEX2(0,1,4)]+=-2*C_0*w18 +   C_1*w19;
                            EM_S[INDEX2(0,2,4)]+=   C_0*w18 - 2*C_1*w19;
                            EM_S[INDEX2(0,3,4)]+=  -C_0*w18 -   C_1*w19;
                            EM_S[INDEX2(1,0,4)]+= 2*C_0*w18 +   C_1*w19;
                            EM_S[INDEX2(1,1,4)]+=-2*C_0*w18 + 2*C_1*w19;
                            EM_S[INDEX2(1,2,4)]+=   C_0*w18 -   C_1*w19;
                            EM_S[INDEX2(1,3,4)]+=  -C_0*w18 - 2*C_1*w19;
                            EM_S[INDEX2(2,0,4)]+=   C_0*w18 + 2*C_1*w19;
                            EM_S[INDEX2(2,1,4)]+=  -C_0*w18 +   C_1*w19;
                            EM_S[INDEX2(2,2,4)]+= 2*C_0*w18 - 2*C_1*w19;
                            EM_S[INDEX2(2,3,4)]+=-2*C_0*w18 -   C_1*w19;
                            EM_S[INDEX2(3,0,4)]+=   C_0*w18 +   C_1*w19;
                            EM_S[INDEX2(3,1,4)]+=  -C_0*w18 + 2*C_1*w19;
                            EM_S[INDEX2(3,2,4)]+= 2*C_0*w18 -   C_1*w19;
                            EM_S[INDEX2(3,3,4)]+=-2*C_0*w18 - 2*C_1*w19;
                        }
                    }
                    ///////////////
                    // process D //
                    ///////////////
                    if (!D.isEmpty()) {
                        add_EM_S=true;
                        const double* D_p=const_cast<escript::Data*>(&D)->getSampleDataRO(e);
                        if (D.actsExpanded()) {
                            const double D_0 = D_p[0];
                            const double D_1 = D_p[1];
                            const double D_2 = D_p[2];
                            const double D_3 = D_p[3];
                            const double tmp0 = w21*(D_2 + D_3);
                            const double tmp1 = w20*(D_0 + D_1);
                            const double tmp2 = w22*(D_0 + D_1 + D_2 + D_3);
                            const double tmp3 = w21*(D_0 + D_1);
                            const double tmp4 = w20*(D_2 + D_3);
                            const double tmp5 = w22*(D_1 + D_2);
                            const double tmp6 = w21*(D_0 + D_2);
                            const double tmp7 = w20*(D_1 + D_3);
                            const double tmp8 = w21*(D_1 + D_3);
                            const double tmp9 = w20*(D_0 + D_2);
                            const double tmp10 = w22*(D_0 + D_3);
                            EM_S[INDEX2(0,0,4)]+=D_0*w23 + D_3*w24 + tmp5;
                            EM_S[INDEX2(0,1,4)]+=tmp0 + tmp1;
                            EM_S[INDEX2(0,2,4)]+=tmp8 + tmp9;
                            EM_S[INDEX2(0,3,4)]+=tmp2;
                            EM_S[INDEX2(1,0,4)]+=tmp0 + tmp1;
                            EM_S[INDEX2(1,1,4)]+=D_1*w23 + D_2*w24 + tmp10;
                            EM_S[INDEX2(1,2,4)]+=tmp2;
                            EM_S[INDEX2(1,3,4)]+=tmp6 + tmp7;
                            EM_S[INDEX2(2,0,4)]+=tmp8 + tmp9;
                            EM_S[INDEX2(2,1,4)]+=tmp2;
                            EM_S[INDEX2(2,2,4)]+=D_1*w24 + D_2*w23 + tmp10;
                            EM_S[INDEX2(2,3,4)]+=tmp3 + tmp4;
                            EM_S[INDEX2(3,0,4)]+=tmp2;
                            EM_S[INDEX2(3,1,4)]+=tmp6 + tmp7;
                            EM_S[INDEX2(3,2,4)]+=tmp3 + tmp4;
                            EM_S[INDEX2(3,3,4)]+=D_0*w24 + D_3*w23 + tmp5;
                        } else { // constant data
                            const double D_0 = D_p[0];
                            EM_S[INDEX2(0,0,4)]+=16*D_0*w22;
                            EM_S[INDEX2(0,1,4)]+=8*D_0*w22;
                            EM_S[INDEX2(0,2,4)]+=8*D_0*w22;
                            EM_S[INDEX2(0,3,4)]+=4*D_0*w22;
                            EM_S[INDEX2(1,0,4)]+=8*D_0*w22;
                            EM_S[INDEX2(1,1,4)]+=16*D_0*w22;
                            EM_S[INDEX2(1,2,4)]+=4*D_0*w22;
                            EM_S[INDEX2(1,3,4)]+=8*D_0*w22;
                            EM_S[INDEX2(2,0,4)]+=8*D_0*w22;
                            EM_S[INDEX2(2,1,4)]+=4*D_0*w22;
                            EM_S[INDEX2(2,2,4)]+=16*D_0*w22;
                            EM_S[INDEX2(2,3,4)]+=8*D_0*w22;
                            EM_S[INDEX2(3,0,4)]+=4*D_0*w22;
                            EM_S[INDEX2(3,1,4)]+=8*D_0*w22;
                            EM_S[INDEX2(3,2,4)]+=8*D_0*w22;
                            EM_S[INDEX2(3,3,4)]+=16*D_0*w22;
                        }
                    }
                    ///////////////
                    // process X //
                    ///////////////
                    if (!X.isEmpty()) {
                        add_EM_F=true;
                        const double* X_p=const_cast<escript::Data*>(&X)->getSampleDataRO(e);
                        if (X.actsExpanded()) {
                            const double X_0_0 = X_p[INDEX2(0,0,2)];
                            const double X_1_0 = X_p[INDEX2(1,0,2)];
                            const double X_0_1 = X_p[INDEX2(0,1,2)];
                            const double X_1_1 = X_p[INDEX2(1,1,2)];
                            const double X_0_2 = X_p[INDEX2(0,2,2)];
                            const double X_1_2 = X_p[INDEX2(1,2,2)];
                            const double X_0_3 = X_p[INDEX2(0,3,2)];
                            const double X_1_3 = X_p[INDEX2(1,3,2)];
                            const double tmp0 = 6*w15*(X_0_2 + X_0_3);
                            const double tmp1 = 6*w10*(X_0_0 + X_0_1);
                            const double tmp2 = 6*w11*(X_1_0 + X_1_2);
                            const double tmp3 = 6*w14*(X_1_1 + X_1_3);
                            const double tmp4 = 6*w11*(X_1_1 + X_1_3);
                            const double tmp5 = w25*(X_0_0 + X_0_1);
                            const double tmp6 = w26*(X_0_2 + X_0_3);
                            const double tmp7 = 6*w14*(X_1_0 + X_1_2);
                            const double tmp8 = w27*(X_1_0 + X_1_2);
                            const double tmp9 = w28*(X_1_1 + X_1_3);
                            const double tmp10 = w25*(-X_0_2 - X_0_3);
                            const double tmp11 = w26*(-X_0_0 - X_0_1);
                            const double tmp12 = w27*(X_1_1 + X_1_3);
                            const double tmp13 = w28*(X_1_0 + X_1_2);
                            const double tmp14 = w25*(X_0_2 + X_0_3);
                            const double tmp15 = w26*(X_0_0 + X_0_1);
                            EM_F[0]+=tmp0 + tmp1 + tmp2 + tmp3;
                            EM_F[1]+=tmp4 + tmp5 + tmp6 + tmp7;
                            EM_F[2]+=tmp10 + tmp11 + tmp8 + tmp9;
                            EM_F[3]+=tmp12 + tmp13 + tmp14 + tmp15;
                        } else { // constant data
                            const double X_0 = X_p[0];
                            const double X_1 = X_p[1];
                            EM_F[0]+= 6*X_0*w18 + 6*X_1*w19;
                            EM_F[1]+=-6*X_0*w18 + 6*X_1*w19;
                            EM_F[2]+= 6*X_0*w18 - 6*X_1*w19;
                            EM_F[3]+=-6*X_0*w18 - 6*X_1*w19;
                        }
                    }
                    ///////////////
                    // process Y //
                    ///////////////
                    if (!Y.isEmpty()) {
                        add_EM_F=true;
                        const double* Y_p=const_cast<escript::Data*>(&Y)->getSampleDataRO(e);
                        if (Y.actsExpanded()) {
                            const double Y_0 = Y_p[0];
                            const double Y_1 = Y_p[1];
                            const double Y_2 = Y_p[2];
                            const double Y_3 = Y_p[3];
                            const double tmp0 = 6*w22*(Y_1 + Y_2);
                            const double tmp1 = 6*w22*(Y_0 + Y_3);
                            EM_F[0]+=6*Y_0*w20 + 6*Y_3*w21 + tmp0;
                            EM_F[1]+=6*Y_1*w20 + 6*Y_2*w21 + tmp1;
                            EM_F[2]+=6*Y_1*w21 + 6*Y_2*w20 + tmp1;
                            EM_F[3]+=6*Y_0*w21 + 6*Y_3*w20 + tmp0;
                        } else { // constant data
                            EM_F[0]+=36*Y_p[0]*w22;
                            EM_F[1]+=36*Y_p[0]*w22;
                            EM_F[2]+=36*Y_p[0]*w22;
                            EM_F[3]+=36*Y_p[0]*w22;
                        }
                    }

                    // add to matrix (if add_EM_S) and RHS (if add_EM_F)
                    const index_t firstNode=m_NN[0]*k1+k0;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F, firstNode);
                } // end k0 loop
            } // end k1 loop
        } // end of colouring
    } // end of parallel region
}

//protected
void Rectangle::assemblePDESingleReduced(Paso_SystemMatrix* mat,
        escript::Data& rhs, const escript::Data& A, const escript::Data& B,
        const escript::Data& C, const escript::Data& D,
        const escript::Data& X, const escript::Data& Y) const
{
    const double w0 = 1./4;
    const double w1 = m_dx[0]/8;
    const double w2 = m_dx[1]/8;
    const double w3 = m_dx[0]*m_dx[1]/16;
    const double w4 = m_dx[0]/(4*m_dx[1]);
    const double w5 = m_dx[1]/(4*m_dx[0]);

    rhs.requireWrite();
#pragma omp parallel
    {
        for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
            for (index_t k1=k1_0; k1<m_NE[1]; k1+=2) {
                for (index_t k0=0; k0<m_NE[0]; ++k0)  {
                    bool add_EM_S=false;
                    bool add_EM_F=false;
                    vector<double> EM_S(4*4, 0);
                    vector<double> EM_F(4, 0);
                    const index_t e = k0 + m_NE[0]*k1;
                    ///////////////
                    // process A //
                    ///////////////
                    if (!A.isEmpty()) {
                        add_EM_S=true;
                        const double* A_p=const_cast<escript::Data*>(&A)->getSampleDataRO(e);
                        const double A_00 = A_p[INDEX2(0,0,2)];
                        const double A_10 = A_p[INDEX2(1,0,2)];
                        const double A_01 = A_p[INDEX2(0,1,2)];
                        const double A_11 = A_p[INDEX2(1,1,2)];
                        const double tmp0 = (A_01 + A_10)*w0;
                        const double tmp1 = A_00*w5;
                        const double tmp2 = A_01*w0;
                        const double tmp3 = A_10*w0;
                        const double tmp4 = A_11*w4;
                        EM_S[INDEX2(0,0,4)]+=tmp4 + tmp0 + tmp1;
                        EM_S[INDEX2(1,0,4)]+=tmp4 - tmp1 + tmp3 - tmp2;
                        EM_S[INDEX2(2,0,4)]+=tmp2 - tmp3 - tmp4 + tmp1;
                        EM_S[INDEX2(3,0,4)]+=-tmp1 - tmp4 - tmp0;
                        EM_S[INDEX2(0,1,4)]+=tmp4 - tmp1 + tmp2 - tmp3;
                        EM_S[INDEX2(1,1,4)]+=tmp4 + tmp1 - tmp0;
                        EM_S[INDEX2(2,1,4)]+=-tmp1 + tmp0 - tmp4;
                        EM_S[INDEX2(3,1,4)]+=-tmp4 + tmp1 + tmp3 - tmp2;
                        EM_S[INDEX2(0,2,4)]+=-tmp4 + tmp1 + tmp3 - tmp2;
                        EM_S[INDEX2(1,2,4)]+=-tmp1 + tmp0 - tmp4;
                        EM_S[INDEX2(2,2,4)]+=tmp4 + tmp1 - tmp0;
                        EM_S[INDEX2(3,2,4)]+=tmp4 - tmp1 + tmp2 - tmp3;
                        EM_S[INDEX2(0,3,4)]+=-tmp1 - tmp4 - tmp0;
                        EM_S[INDEX2(1,3,4)]+=tmp2 - tmp3 - tmp4 + tmp1;
                        EM_S[INDEX2(2,3,4)]+=tmp4 - tmp1 + tmp3 - tmp2;
                        EM_S[INDEX2(3,3,4)]+=tmp4 + tmp0 + tmp1;
                    }
                    ///////////////
                    // process B //
                    ///////////////
                    if (!B.isEmpty()) {
                        add_EM_S=true;
                        const double* B_p=const_cast<escript::Data*>(&B)->getSampleDataRO(e);
                        const double tmp0 = B_p[0]*w2;
                        const double tmp1 = B_p[1]*w1;
                        EM_S[INDEX2(0,0,4)]+=-tmp0 - tmp1;
                        EM_S[INDEX2(1,0,4)]+= tmp0 - tmp1;
                        EM_S[INDEX2(2,0,4)]+= tmp1 - tmp0;
                        EM_S[INDEX2(3,0,4)]+= tmp0 + tmp1;
                        EM_S[INDEX2(0,1,4)]+=-tmp0 - tmp1;
                        EM_S[INDEX2(1,1,4)]+= tmp0 - tmp1;
                        EM_S[INDEX2(2,1,4)]+= tmp1 - tmp0;
                        EM_S[INDEX2(3,1,4)]+= tmp0 + tmp1;
                        EM_S[INDEX2(0,2,4)]+=-tmp0 - tmp1;
                        EM_S[INDEX2(1,2,4)]+= tmp0 - tmp1;
                        EM_S[INDEX2(2,2,4)]+= tmp1 - tmp0;
                        EM_S[INDEX2(3,2,4)]+= tmp0 + tmp1;
                        EM_S[INDEX2(0,3,4)]+=-tmp0 - tmp1;
                        EM_S[INDEX2(1,3,4)]+= tmp0 - tmp1;
                        EM_S[INDEX2(2,3,4)]+= tmp1 - tmp0;
                        EM_S[INDEX2(3,3,4)]+= tmp0 + tmp1;
                    }
                    ///////////////
                    // process C //
                    ///////////////
                    if (!C.isEmpty()) {
                        add_EM_S=true;
                        const double* C_p=const_cast<escript::Data*>(&C)->getSampleDataRO(e);
                        const double tmp0 = C_p[0]*w2;
                        const double tmp1 = C_p[1]*w1;
                        EM_S[INDEX2(0,0,4)]+=-tmp1 - tmp0;
                        EM_S[INDEX2(1,0,4)]+=-tmp1 - tmp0;
                        EM_S[INDEX2(2,0,4)]+=-tmp1 - tmp0;
                        EM_S[INDEX2(3,0,4)]+=-tmp1 - tmp0;
                        EM_S[INDEX2(0,1,4)]+= tmp0 - tmp1;
                        EM_S[INDEX2(1,1,4)]+= tmp0 - tmp1;
                        EM_S[INDEX2(2,1,4)]+= tmp0 - tmp1;
                        EM_S[INDEX2(3,1,4)]+= tmp0 - tmp1;
                        EM_S[INDEX2(0,2,4)]+= tmp1 - tmp0;
                        EM_S[INDEX2(1,2,4)]+= tmp1 - tmp0;
                        EM_S[INDEX2(2,2,4)]+= tmp1 - tmp0;
                        EM_S[INDEX2(3,2,4)]+= tmp1 - tmp0;
                        EM_S[INDEX2(0,3,4)]+= tmp0 + tmp1;
                        EM_S[INDEX2(1,3,4)]+= tmp0 + tmp1;
                        EM_S[INDEX2(2,3,4)]+= tmp0 + tmp1;
                        EM_S[INDEX2(3,3,4)]+= tmp0 + tmp1;
                    }
                    ///////////////
                    // process D //
                    ///////////////
                    if (!D.isEmpty()) {
                        add_EM_S=true;
                        const double* D_p=const_cast<escript::Data*>(&D)->getSampleDataRO(e);
                        EM_S[INDEX2(0,0,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(1,0,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(2,0,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(3,0,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(0,1,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(1,1,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(2,1,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(3,1,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(0,2,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(1,2,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(2,2,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(3,2,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(0,3,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(1,3,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(2,3,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(3,3,4)]+=D_p[0]*w3;
                    }
                    ///////////////
                    // process X //
                    ///////////////
                    if (!X.isEmpty()) {
                        add_EM_F=true;
                        const double* X_p=const_cast<escript::Data*>(&X)->getSampleDataRO(e);
                        const double wX0 = 4*X_p[0]*w2;
                        const double wX1 = 4*X_p[1]*w1;
                        EM_F[0]+=-wX0 - wX1;
                        EM_F[1]+=-wX1 + wX0;
                        EM_F[2]+=-wX0 + wX1;
                        EM_F[3]+= wX0 + wX1;
                    }
                    ///////////////
                    // process Y //
                    ///////////////
                    if (!Y.isEmpty()) {
                        add_EM_F=true;
                        const double* Y_p=const_cast<escript::Data*>(&Y)->getSampleDataRO(e);
                        EM_F[0]+=4*Y_p[0]*w3;
                        EM_F[1]+=4*Y_p[0]*w3;
                        EM_F[2]+=4*Y_p[0]*w3;
                        EM_F[3]+=4*Y_p[0]*w3;
                    }

                    // add to matrix (if add_EM_S) and RHS (if add_EM_F)
                    const index_t firstNode=m_NN[0]*k1+k0;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F, firstNode);
                } // end k0 loop
            } // end k1 loop
        } // end of colouring
    } // end of parallel region
}

//protected
void Rectangle::assemblePDESystem(Paso_SystemMatrix* mat,
        escript::Data& rhs, const escript::Data& A, const escript::Data& B,
        const escript::Data& C, const escript::Data& D,
        const escript::Data& X, const escript::Data& Y) const
{
    dim_t numEq, numComp;
    if (!mat)
        numEq=numComp=(rhs.isEmpty() ? 1 : rhs.getDataPointSize());
    else {
        numEq=mat->logical_row_block_size;
        numComp=mat->logical_col_block_size;
    }
    const double SQRT3 = 1.73205080756887719318;
    const double w1 = 1.0/24;
    const double w5 = -SQRT3/24 + 1.0/12;
    const double w2 = -SQRT3/24 - 1.0/12;
    const double w19 = -m_dx[0]/12;
    const double w11 = w19*(SQRT3 + 3)/12;
    const double w14 = w19*(-SQRT3 + 3)/12;
    const double w16 = w19*(5*SQRT3 + 9)/12;
    const double w17 = w19*(-5*SQRT3 + 9)/12;
    const double w27 = w19*(-SQRT3 - 3)/2;
    const double w28 = w19*(SQRT3 - 3)/2;
    const double w18 = -m_dx[1]/12;
    const double w10 = w18*(SQRT3 + 3)/12;
    const double w15 = w18*(-SQRT3 + 3)/12;
    const double w12 = w18*(5*SQRT3 + 9)/12;
    const double w13 = w18*(-5*SQRT3 + 9)/12;
    const double w25 = w18*(-SQRT3 - 3)/2;
    const double w26 = w18*(SQRT3 - 3)/2;
    const double w22 = m_dx[0]*m_dx[1]/144;
    const double w20 = w22*(SQRT3 + 2);
    const double w21 = w22*(-SQRT3 + 2);
    const double w23 = w22*(4*SQRT3 + 7);
    const double w24 = w22*(-4*SQRT3 + 7);
    const double w3 = m_dx[0]/(24*m_dx[1]);
    const double w7 = w3*(SQRT3 + 2);
    const double w8 = w3*(-SQRT3 + 2);
    const double w6 = -m_dx[1]/(24*m_dx[0]);
    const double w0 = w6*(SQRT3 + 2);
    const double w4 = w6*(-SQRT3 + 2);

    rhs.requireWrite();
#pragma omp parallel
    {
        for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
            for (index_t k1=k1_0; k1<m_NE[1]; k1+=2) {
                for (index_t k0=0; k0<m_NE[0]; ++k0)  {
                    bool add_EM_S=false;
                    bool add_EM_F=false;
                    vector<double> EM_S(4*4*numEq*numComp, 0);
                    vector<double> EM_F(4*numEq, 0);
                    const index_t e = k0 + m_NE[0]*k1;
                    ///////////////
                    // process A //
                    ///////////////
                    if (!A.isEmpty()) {
                        add_EM_S = true;
                        const double* A_p = const_cast<escript::Data*>(&A)->getSampleDataRO(e);
                        if (A.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double A_00_0 = A_p[INDEX5(k,0,m,0,0,numEq,2,numComp,2)];
                                    const double A_01_0 = A_p[INDEX5(k,0,m,1,0,numEq,2,numComp,2)];
                                    const double A_10_0 = A_p[INDEX5(k,1,m,0,0,numEq,2,numComp,2)];
                                    const double A_11_0 = A_p[INDEX5(k,1,m,1,0,numEq,2,numComp,2)];
                                    const double A_00_1 = A_p[INDEX5(k,0,m,0,1,numEq,2,numComp,2)];
                                    const double A_01_1 = A_p[INDEX5(k,0,m,1,1,numEq,2,numComp,2)];
                                    const double A_10_1 = A_p[INDEX5(k,1,m,0,1,numEq,2,numComp,2)];
                                    const double A_11_1 = A_p[INDEX5(k,1,m,1,1,numEq,2,numComp,2)];
                                    const double A_00_2 = A_p[INDEX5(k,0,m,0,2,numEq,2,numComp,2)];
                                    const double A_01_2 = A_p[INDEX5(k,0,m,1,2,numEq,2,numComp,2)];
                                    const double A_10_2 = A_p[INDEX5(k,1,m,0,2,numEq,2,numComp,2)];
                                    const double A_11_2 = A_p[INDEX5(k,1,m,1,2,numEq,2,numComp,2)];
                                    const double A_00_3 = A_p[INDEX5(k,0,m,0,3,numEq,2,numComp,2)];
                                    const double A_01_3 = A_p[INDEX5(k,0,m,1,3,numEq,2,numComp,2)];
                                    const double A_10_3 = A_p[INDEX5(k,1,m,0,3,numEq,2,numComp,2)];
                                    const double A_11_3 = A_p[INDEX5(k,1,m,1,3,numEq,2,numComp,2)];
                                    const double tmp0 = w3*(A_11_0 + A_11_1 + A_11_2 + A_11_3);
                                    const double tmp1 = w1*(A_01_0 + A_01_3 - A_10_1 - A_10_2);
                                    const double tmp2 = w4*(A_00_2 + A_00_3);
                                    const double tmp3 = w0*(A_00_0 + A_00_1);
                                    const double tmp4 = w5*(A_01_2 - A_10_3);
                                    const double tmp5 = w2*(-A_01_1 + A_10_0);
                                    const double tmp6 = w5*(A_01_3 + A_10_0);
                                    const double tmp7 = w3*(-A_11_0 - A_11_1 - A_11_2 - A_11_3);
                                    const double tmp8 = w6*(A_00_0 + A_00_1 + A_00_2 + A_00_3);
                                    const double tmp9 = w1*(A_01_1 + A_01_2 + A_10_1 + A_10_2);
                                    const double tmp10 = w2*(-A_01_0 - A_10_3);
                                    const double tmp11 = w4*(A_00_0 + A_00_1);
                                    const double tmp12 = w0*(A_00_2 + A_00_3);
                                    const double tmp13 = w5*(A_01_1 - A_10_0);
                                    const double tmp14 = w2*(-A_01_2 + A_10_3);
                                    const double tmp15 = w7*(A_11_0 + A_11_2);
                                    const double tmp16 = w4*(-A_00_2 - A_00_3);
                                    const double tmp17 = w0*(-A_00_0 - A_00_1);
                                    const double tmp18 = w5*(A_01_3 + A_10_3);
                                    const double tmp19 = w8*(A_11_1 + A_11_3);
                                    const double tmp20 = w2*(-A_01_0 - A_10_0);
                                    const double tmp21 = w7*(A_11_1 + A_11_3);
                                    const double tmp22 = w4*(-A_00_0 - A_00_1);
                                    const double tmp23 = w0*(-A_00_2 - A_00_3);
                                    const double tmp24 = w5*(A_01_0 + A_10_0);
                                    const double tmp25 = w8*(A_11_0 + A_11_2);
                                    const double tmp26 = w2*(-A_01_3 - A_10_3);
                                    const double tmp27 = w5*(-A_01_1 - A_10_2);
                                    const double tmp28 = w1*(-A_01_0 - A_01_3 - A_10_0 - A_10_3);
                                    const double tmp29 = w2*(A_01_2 + A_10_1);
                                    const double tmp30 = w7*(-A_11_1 - A_11_3);
                                    const double tmp31 = w1*(-A_01_1 - A_01_2 + A_10_0 + A_10_3);
                                    const double tmp32 = w5*(-A_01_0 + A_10_2);
                                    const double tmp33 = w8*(-A_11_0 - A_11_2);
                                    const double tmp34 = w6*(-A_00_0 - A_00_1 - A_00_2 - A_00_3);
                                    const double tmp35 = w2*(A_01_3 - A_10_1);
                                    const double tmp36 = w5*(A_01_0 + A_10_3);
                                    const double tmp37 = w2*(-A_01_3 - A_10_0);
                                    const double tmp38 = w7*(-A_11_0 - A_11_2);
                                    const double tmp39 = w5*(-A_01_3 + A_10_1);
                                    const double tmp40 = w8*(-A_11_1 - A_11_3);
                                    const double tmp41 = w2*(A_01_0 - A_10_2);
                                    const double tmp42 = w5*(A_01_1 - A_10_3);
                                    const double tmp43 = w2*(-A_01_2 + A_10_0);
                                    const double tmp44 = w5*(A_01_2 - A_10_0);
                                    const double tmp45 = w2*(-A_01_1 + A_10_3);
                                    const double tmp46 = w5*(-A_01_0 + A_10_1);
                                    const double tmp47 = w2*(A_01_3 - A_10_2);
                                    const double tmp48 = w5*(-A_01_1 - A_10_1);
                                    const double tmp49 = w2*(A_01_2 + A_10_2);
                                    const double tmp50 = w5*(-A_01_3 + A_10_2);
                                    const double tmp51 = w2*(A_01_0 - A_10_1);
                                    const double tmp52 = w5*(-A_01_2 - A_10_1);
                                    const double tmp53 = w2*(A_01_1 + A_10_2);
                                    const double tmp54 = w5*(-A_01_2 - A_10_2);
                                    const double tmp55 = w2*(A_01_1 + A_10_1);
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp15 + tmp16 + tmp17 + tmp18 + tmp19 + tmp20 + tmp9;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp31 + tmp34 + tmp38 + tmp39 + tmp40 + tmp41;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp28 + tmp52 + tmp53 + tmp7 + tmp8;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp0 + tmp2 + tmp3 + tmp31 + tmp50 + tmp51;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp16 + tmp17 + tmp21 + tmp25 + tmp28 + tmp54 + tmp55;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp10 + tmp6 + tmp7 + tmp8 + tmp9;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp1 + tmp30 + tmp33 + tmp34 + tmp44 + tmp45;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp1 + tmp34 + tmp38 + tmp40 + tmp42 + tmp43;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp36 + tmp37 + tmp7 + tmp8 + tmp9;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp15 + tmp19 + tmp22 + tmp23 + tmp28 + tmp48 + tmp49;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp0 + tmp11 + tmp12 + tmp31 + tmp46 + tmp47;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp27 + tmp28 + tmp29 + tmp7 + tmp8;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp30 + tmp31 + tmp32 + tmp33 + tmp34 + tmp35;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp0 + tmp1 + tmp11 + tmp12 + tmp13 + tmp14;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp21 + tmp22 + tmp23 + tmp24 + tmp25 + tmp26 + tmp9;
                                }
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double A_00 = A_p[INDEX4(k,0,m,0, numEq,2, numComp)];
                                    const double A_01 = A_p[INDEX4(k,0,m,1, numEq,2, numComp)];
                                    const double A_10 = A_p[INDEX4(k,1,m,0, numEq,2, numComp)];
                                    const double A_11 = A_p[INDEX4(k,1,m,1, numEq,2, numComp)];
                                    const double tmp0 = 6*w1*(A_01 - A_10);
                                    const double tmp1 = 6*w1*(A_01 + A_10);
                                    const double tmp2 = 6*w1*(-A_01 - A_10);
                                    const double tmp3 = 6*w1*(-A_01 + A_10);
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=-8*A_00*w6 + 8*A_11*w3 + tmp1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=8*A_00*w6 + 4*A_11*w3 + tmp0;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=-4*A_00*w6 - 8*A_11*w3 + tmp3;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=4*A_00*w6 - 4*A_11*w3 + tmp2;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=8*A_00*w6 + 4*A_11*w3 + tmp3;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=-8*A_00*w6 + 8*A_11*w3 + tmp2;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=4*A_00*w6 - 4*A_11*w3 + tmp1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=-4*A_00*w6 - 8*A_11*w3 + tmp0;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=-4*A_00*w6 - 8*A_11*w3 + tmp0;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=4*A_00*w6 - 4*A_11*w3 + tmp1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=-8*A_00*w6 + 8*A_11*w3 + tmp2;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=8*A_00*w6 + 4*A_11*w3 + tmp3;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=4*A_00*w6 - 4*A_11*w3 + tmp2;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=-4*A_00*w6 - 8*A_11*w3 + tmp3;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=8*A_00*w6 + 4*A_11*w3 + tmp0;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=-8*A_00*w6 + 8*A_11*w3 + tmp1;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process B //
                    ///////////////
                    if (!B.isEmpty()) {
                        add_EM_S=true;
                        const double* B_p=const_cast<escript::Data*>(&B)->getSampleDataRO(e);
                        if (B.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double B_0_0 = B_p[INDEX4(k,0,m,0, numEq,2,numComp)];
                                    const double B_1_0 = B_p[INDEX4(k,1,m,0, numEq,2,numComp)];
                                    const double B_0_1 = B_p[INDEX4(k,0,m,1, numEq,2,numComp)];
                                    const double B_1_1 = B_p[INDEX4(k,1,m,1, numEq,2,numComp)];
                                    const double B_0_2 = B_p[INDEX4(k,0,m,2, numEq,2,numComp)];
                                    const double B_1_2 = B_p[INDEX4(k,1,m,2, numEq,2,numComp)];
                                    const double B_0_3 = B_p[INDEX4(k,0,m,3, numEq,2,numComp)];
                                    const double B_1_3 = B_p[INDEX4(k,1,m,3, numEq,2,numComp)];
                                    const double tmp0 = w11*(B_1_0 + B_1_1);
                                    const double tmp1 = w14*(B_1_2 + B_1_3);
                                    const double tmp2 = w15*(-B_0_1 - B_0_3);
                                    const double tmp3 = w10*(-B_0_0 - B_0_2);
                                    const double tmp4 = w11*(B_1_2 + B_1_3);
                                    const double tmp5 = w14*(B_1_0 + B_1_1);
                                    const double tmp6 = w11*(-B_1_2 - B_1_3);
                                    const double tmp7 = w14*(-B_1_0 - B_1_1);
                                    const double tmp8 = w11*(-B_1_0 - B_1_1);
                                    const double tmp9 = w14*(-B_1_2 - B_1_3);
                                    const double tmp10 = w10*(-B_0_1 - B_0_3);
                                    const double tmp11 = w15*(-B_0_0 - B_0_2);
                                    const double tmp12 = w15*(B_0_0 + B_0_2);
                                    const double tmp13 = w10*(B_0_1 + B_0_3);
                                    const double tmp14 = w10*(B_0_0 + B_0_2);
                                    const double tmp15 = w15*(B_0_1 + B_0_3);
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=B_0_0*w12 + B_0_1*w10 + B_0_2*w15 + B_0_3*w13 + B_1_0*w16 + B_1_1*w14 + B_1_2*w11 + B_1_3*w17;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=B_0_0*w10 + B_0_1*w12 + B_0_2*w13 + B_0_3*w15 + tmp0 + tmp1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=B_1_0*w11 + B_1_1*w17 + B_1_2*w16 + B_1_3*w14 + tmp14 + tmp15;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp12 + tmp13 + tmp4 + tmp5;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=-B_0_0*w12 - B_0_1*w10 - B_0_2*w15 - B_0_3*w13 + tmp0 + tmp1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=-B_0_0*w10 - B_0_1*w12 - B_0_2*w13 - B_0_3*w15 + B_1_0*w14 + B_1_1*w16 + B_1_2*w17 + B_1_3*w11;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp2 + tmp3 + tmp4 + tmp5;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=B_1_0*w17 + B_1_1*w11 + B_1_2*w14 + B_1_3*w16 + tmp10 + tmp11;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=-B_1_0*w16 - B_1_1*w14 - B_1_2*w11 - B_1_3*w17 + tmp14 + tmp15;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp12 + tmp13 + tmp8 + tmp9;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=B_0_0*w15 + B_0_1*w13 + B_0_2*w12 + B_0_3*w10 - B_1_0*w11 - B_1_1*w17 - B_1_2*w16 - B_1_3*w14;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=B_0_0*w13 + B_0_1*w15 + B_0_2*w10 + B_0_3*w12 + tmp6 + tmp7;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp2 + tmp3 + tmp8 + tmp9;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=-B_1_0*w14 - B_1_1*w16 - B_1_2*w17 - B_1_3*w11 + tmp10 + tmp11;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=-B_0_0*w15 - B_0_1*w13 - B_0_2*w12 - B_0_3*w10 + tmp6 + tmp7;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=-B_0_0*w13 - B_0_1*w15 - B_0_2*w10 - B_0_3*w12 - B_1_0*w17 - B_1_1*w11 - B_1_2*w14 - B_1_3*w16;
                                }
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double wB0 = B_p[INDEX3(k,0,m,numEq,2)]*w18;
                                    const double wB1 = B_p[INDEX3(k,1,m,numEq,2)]*w19;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+= 2*wB0 + 2*wB1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+= 2*wB0 +   wB1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=   wB0 + 2*wB1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=   wB0 +   wB1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=-2*wB0 +   wB1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=-2*wB0 + 2*wB1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=  -wB0 +   wB1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=  -wB0 + 2*wB1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=   wB0 - 2*wB1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=   wB0 -   wB1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+= 2*wB0 - 2*wB1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+= 2*wB0 -   wB1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=  -wB0 -   wB1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=  -wB0 - 2*wB1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=-2*wB0 -   wB1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=-2*wB0 - 2*wB1;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process C //
                    ///////////////
                    if (!C.isEmpty()) {
                        add_EM_S=true;
                        const double* C_p=const_cast<escript::Data*>(&C)->getSampleDataRO(e);
                        if (C.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double C_0_0 = C_p[INDEX4(k,m,0, 0, numEq,numComp,2)];
                                    const double C_1_0 = C_p[INDEX4(k,m,1, 0, numEq,numComp,2)];
                                    const double C_0_1 = C_p[INDEX4(k,m,0, 1, numEq,numComp,2)];
                                    const double C_1_1 = C_p[INDEX4(k,m,1, 1, numEq,numComp,2)];
                                    const double C_0_2 = C_p[INDEX4(k,m,0, 2, numEq,numComp,2)];
                                    const double C_1_2 = C_p[INDEX4(k,m,1, 2, numEq,numComp,2)];
                                    const double C_0_3 = C_p[INDEX4(k,m,0, 3, numEq,numComp,2)];
                                    const double C_1_3 = C_p[INDEX4(k,m,1, 3, numEq,numComp,2)];
                                    const double tmp0 = w11*(C_1_0 + C_1_1);
                                    const double tmp1 = w14*(C_1_2 + C_1_3);
                                    const double tmp2 = w15*(C_0_0 + C_0_2);
                                    const double tmp3 = w10*(C_0_1 + C_0_3);
                                    const double tmp4 = w11*(-C_1_0 - C_1_1);
                                    const double tmp5 = w14*(-C_1_2 - C_1_3);
                                    const double tmp6 = w11*(-C_1_2 - C_1_3);
                                    const double tmp7 = w14*(-C_1_0 - C_1_1);
                                    const double tmp8 = w11*(C_1_2 + C_1_3);
                                    const double tmp9 = w14*(C_1_0 + C_1_1);
                                    const double tmp10 = w10*(-C_0_1 - C_0_3);
                                    const double tmp11 = w15*(-C_0_0 - C_0_2);
                                    const double tmp12 = w15*(-C_0_1 - C_0_3);
                                    const double tmp13 = w10*(-C_0_0 - C_0_2);
                                    const double tmp14 = w10*(C_0_0 + C_0_2);
                                    const double tmp15 = w15*(C_0_1 + C_0_3);
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=C_0_0*w12 + C_0_1*w10 + C_0_2*w15 + C_0_3*w13 + C_1_0*w16 + C_1_1*w14 + C_1_2*w11 + C_1_3*w17;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=-C_0_0*w12 - C_0_1*w10 - C_0_2*w15 - C_0_3*w13 + tmp0 + tmp1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=-C_1_0*w16 - C_1_1*w14 - C_1_2*w11 - C_1_3*w17 + tmp14 + tmp15;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp12 + tmp13 + tmp4 + tmp5;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=C_0_0*w10 + C_0_1*w12 + C_0_2*w13 + C_0_3*w15 + tmp0 + tmp1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=-C_0_0*w10 - C_0_1*w12 - C_0_2*w13 - C_0_3*w15 + C_1_0*w14 + C_1_1*w16 + C_1_2*w17 + C_1_3*w11;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp2 + tmp3 + tmp4 + tmp5;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=-C_1_0*w14 - C_1_1*w16 - C_1_2*w17 - C_1_3*w11 + tmp10 + tmp11;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=C_1_0*w11 + C_1_1*w17 + C_1_2*w16 + C_1_3*w14 + tmp14 + tmp15;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp12 + tmp13 + tmp8 + tmp9;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=C_0_0*w15 + C_0_1*w13 + C_0_2*w12 + C_0_3*w10 - C_1_0*w11 - C_1_1*w17 - C_1_2*w16 - C_1_3*w14;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=-C_0_0*w15 - C_0_1*w13 - C_0_2*w12 - C_0_3*w10 + tmp6 + tmp7;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp2 + tmp3 + tmp8 + tmp9;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=C_1_0*w17 + C_1_1*w11 + C_1_2*w14 + C_1_3*w16 + tmp10 + tmp11;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=C_0_0*w13 + C_0_1*w15 + C_0_2*w10 + C_0_3*w12 + tmp6 + tmp7;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=-C_0_0*w13 - C_0_1*w15 - C_0_2*w10 - C_0_3*w12 - C_1_0*w17 - C_1_1*w11 - C_1_2*w14 - C_1_3*w16;
                                }
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double wC0 = C_p[INDEX3(k,m,0,numEq,numComp)]*w18;
                                    const double wC1 = C_p[INDEX3(k,m,1,numEq,numComp)]*w19;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+= 2*wC0 + 2*wC1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=-2*wC0 +   wC1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=   wC0 - 2*wC1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=  -wC0 -   wC1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+= 2*wC0 +   wC1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=-2*wC0 + 2*wC1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=   wC0 -   wC1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=  -wC0 - 2*wC1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=   wC0 + 2*wC1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=  -wC0 +   wC1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+= 2*wC0 - 2*wC1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=-2*wC0 -   wC1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=   wC0 +   wC1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=  -wC0 + 2*wC1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+= 2*wC0 -   wC1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=-2*wC0 - 2*wC1;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process D //
                    ///////////////
                    if (!D.isEmpty()) {
                        add_EM_S=true;
                        const double* D_p=const_cast<escript::Data*>(&D)->getSampleDataRO(e);
                        if (D.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double D_0 = D_p[INDEX3(k,m,0,numEq,numComp)];
                                    const double D_1 = D_p[INDEX3(k,m,1,numEq,numComp)];
                                    const double D_2 = D_p[INDEX3(k,m,2,numEq,numComp)];
                                    const double D_3 = D_p[INDEX3(k,m,3,numEq,numComp)];
                                    const double tmp0 = w21*(D_2 + D_3);
                                    const double tmp1 = w20*(D_0 + D_1);
                                    const double tmp2 = w22*(D_0 + D_1 + D_2 + D_3);
                                    const double tmp3 = w21*(D_0 + D_1);
                                    const double tmp4 = w20*(D_2 + D_3);
                                    const double tmp5 = w22*(D_1 + D_2);
                                    const double tmp6 = w21*(D_0 + D_2);
                                    const double tmp7 = w20*(D_1 + D_3);
                                    const double tmp8 = w21*(D_1 + D_3);
                                    const double tmp9 = w20*(D_0 + D_2);
                                    const double tmp10 = w22*(D_0 + D_3);
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=D_0*w23 + D_3*w24 + tmp5;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0 + tmp1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp8 + tmp9;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp2;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp0 + tmp1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=D_1*w23 + D_2*w24 + tmp10;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp2;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp6 + tmp7;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp8 + tmp9;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp2;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=D_1*w24 + D_2*w23 + tmp10;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp3 + tmp4;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp2;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp6 + tmp7;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp3 + tmp4;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=D_0*w24 + D_3*w23 + tmp5;
                                }
                             }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double D_0 = D_p[INDEX2(k, m, numEq)];
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=16*D_0*w22;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=8*D_0*w22;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=8*D_0*w22;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=4*D_0*w22;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=8*D_0*w22;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=16*D_0*w22;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=4*D_0*w22;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=8*D_0*w22;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=8*D_0*w22;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=4*D_0*w22;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=16*D_0*w22;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=8*D_0*w22;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=4*D_0*w22;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=8*D_0*w22;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=8*D_0*w22;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=16*D_0*w22;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process X //
                    ///////////////
                    if (!X.isEmpty()) {
                        add_EM_F=true;
                        const double* X_p=const_cast<escript::Data*>(&X)->getSampleDataRO(e);
                        if (X.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                const double X_0_0 = X_p[INDEX3(k,0,0,numEq,2)];
                                const double X_1_0 = X_p[INDEX3(k,1,0,numEq,2)];
                                const double X_0_1 = X_p[INDEX3(k,0,1,numEq,2)];
                                const double X_1_1 = X_p[INDEX3(k,1,1,numEq,2)];
                                const double X_0_2 = X_p[INDEX3(k,0,2,numEq,2)];
                                const double X_1_2 = X_p[INDEX3(k,1,2,numEq,2)];
                                const double X_0_3 = X_p[INDEX3(k,0,3,numEq,2)];
                                const double X_1_3 = X_p[INDEX3(k,1,3,numEq,2)];
                                const double tmp0 = 6*w15*(X_0_2 + X_0_3);
                                const double tmp1 = 6*w10*(X_0_0 + X_0_1);
                                const double tmp2 = 6*w11*(X_1_0 + X_1_2);
                                const double tmp3 = 6*w14*(X_1_1 + X_1_3);
                                const double tmp4 = 6*w11*(X_1_1 + X_1_3);
                                const double tmp5 = w25*(X_0_0 + X_0_1);
                                const double tmp6 = w26*(X_0_2 + X_0_3);
                                const double tmp7 = 6*w14*(X_1_0 + X_1_2);
                                const double tmp8 = w27*(X_1_0 + X_1_2);
                                const double tmp9 = w28*(X_1_1 + X_1_3);
                                const double tmp10 = w25*(-X_0_2 - X_0_3);
                                const double tmp11 = w26*(-X_0_0 - X_0_1);
                                const double tmp12 = w27*(X_1_1 + X_1_3);
                                const double tmp13 = w28*(X_1_0 + X_1_2);
                                const double tmp14 = w25*(X_0_2 + X_0_3);
                                const double tmp15 = w26*(X_0_0 + X_0_1);
                                EM_F[INDEX2(k,0,numEq)]+=tmp0 + tmp1 + tmp2 + tmp3;
                                EM_F[INDEX2(k,1,numEq)]+=tmp4 + tmp5 + tmp6 + tmp7;
                                EM_F[INDEX2(k,2,numEq)]+=tmp10 + tmp11 + tmp8 + tmp9;
                                EM_F[INDEX2(k,3,numEq)]+=tmp12 + tmp13 + tmp14 + tmp15;
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                const double wX0 = X_p[INDEX2(k, 0, numEq)]*w18;
                                const double wX1 = X_p[INDEX2(k, 1, numEq)]*w19;
                                EM_F[INDEX2(k,0,numEq)]+= 6*wX0 + 6*wX1;
                                EM_F[INDEX2(k,1,numEq)]+=-6*wX0 + 6*wX1;
                                EM_F[INDEX2(k,2,numEq)]+= 6*wX0 - 6*wX1;
                                EM_F[INDEX2(k,3,numEq)]+=-6*wX0 - 6*wX1;
                            }
                        }
                    }
                    ///////////////
                    // process Y //
                    ///////////////
                    if (!Y.isEmpty()) {
                        add_EM_F=true;
                        const double* Y_p=const_cast<escript::Data*>(&Y)->getSampleDataRO(e);
                        if (Y.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                const double Y_0 = Y_p[INDEX2(k, 0, numEq)];
                                const double Y_1 = Y_p[INDEX2(k, 1, numEq)];
                                const double Y_2 = Y_p[INDEX2(k, 2, numEq)];
                                const double Y_3 = Y_p[INDEX2(k, 3, numEq)];
                                const double tmp0 = 6*w22*(Y_1 + Y_2);
                                const double tmp1 = 6*w22*(Y_0 + Y_3);
                                EM_F[INDEX2(k,0,numEq)]+=6*Y_0*w20 + 6*Y_3*w21 + tmp0;
                                EM_F[INDEX2(k,1,numEq)]+=6*Y_1*w20 + 6*Y_2*w21 + tmp1;
                                EM_F[INDEX2(k,2,numEq)]+=6*Y_1*w21 + 6*Y_2*w20 + tmp1;
                                EM_F[INDEX2(k,3,numEq)]+=6*Y_0*w21 + 6*Y_3*w20 + tmp0;
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                EM_F[INDEX2(k,0,numEq)]+=36*Y_p[k]*w22;
                                EM_F[INDEX2(k,1,numEq)]+=36*Y_p[k]*w22;
                                EM_F[INDEX2(k,2,numEq)]+=36*Y_p[k]*w22;
                                EM_F[INDEX2(k,3,numEq)]+=36*Y_p[k]*w22;
                            }
                        }
                    }

                    // add to matrix (if add_EM_S) and RHS (if add_EM_F)
                    const index_t firstNode=m_NN[0]*k1+k0;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F,
                            firstNode, numEq, numComp);
                } // end k0 loop
            } // end k1 loop
        } // end of colouring
    } // end of parallel region
}

//protected
void Rectangle::assemblePDESystemReduced(Paso_SystemMatrix* mat,
        escript::Data& rhs, const escript::Data& A, const escript::Data& B,
        const escript::Data& C, const escript::Data& D,
        const escript::Data& X, const escript::Data& Y) const
{
    dim_t numEq, numComp;
    if (!mat)
        numEq=numComp=(rhs.isEmpty() ? 1 : rhs.getDataPointSize());
    else {
        numEq=mat->logical_row_block_size;
        numComp=mat->logical_col_block_size;
    }

    const double w0 = 1./4;
    const double w1 = m_dx[0]/8;
    const double w2 = m_dx[1]/8;
    const double w3 = m_dx[0]*m_dx[1]/16;
    const double w4 = m_dx[0]/(4*m_dx[1]);
    const double w5 = m_dx[1]/(4*m_dx[0]);

    rhs.requireWrite();
#pragma omp parallel
    {
        for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
            for (index_t k1=k1_0; k1<m_NE[1]; k1+=2) {
                for (index_t k0=0; k0<m_NE[0]; ++k0)  {
                    bool add_EM_S=false;
                    bool add_EM_F=false;
                    vector<double> EM_S(4*4*numEq*numComp, 0);
                    vector<double> EM_F(4*numEq, 0);
                    const index_t e = k0 + m_NE[0]*k1;
                    ///////////////
                    // process A //
                    ///////////////
                    if (!A.isEmpty()) {
                        add_EM_S=true;
                        const double* A_p=const_cast<escript::Data*>(&A)->getSampleDataRO(e);
                        for (index_t k=0; k<numEq; k++) {
                            for (index_t m=0; m<numComp; m++) {
                                const double Aw00 = A_p[INDEX4(k,0,m,0, numEq,2, numComp)]*w5;
                                const double Aw01 = A_p[INDEX4(k,0,m,1, numEq,2, numComp)]*w0;
                                const double Aw10 = A_p[INDEX4(k,1,m,0, numEq,2, numComp)]*w0;
                                const double Aw11 = A_p[INDEX4(k,1,m,1, numEq,2, numComp)]*w4;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+= Aw00 + Aw01 + Aw10 + Aw11;
                                EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=-Aw00 - Aw01 + Aw10 + Aw11;
                                EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+= Aw00 + Aw01 - Aw10 - Aw11;
                                EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=-Aw00 - Aw01 - Aw10 - Aw11;
                                EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=-Aw00 + Aw01 - Aw10 + Aw11;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+= Aw00 - Aw01 - Aw10 + Aw11;
                                EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=-Aw00 + Aw01 + Aw10 - Aw11;
                                EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+= Aw00 - Aw01 + Aw10 - Aw11;
                                EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+= Aw00 - Aw01 + Aw10 - Aw11;
                                EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=-Aw00 + Aw01 + Aw10 - Aw11;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+= Aw00 - Aw01 - Aw10 + Aw11;
                                EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=-Aw00 + Aw01 - Aw10 + Aw11;
                                EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=-Aw00 - Aw01 - Aw10 - Aw11;
                                EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+= Aw00 + Aw01 - Aw10 - Aw11;
                                EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=-Aw00 - Aw01 + Aw10 + Aw11;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+= Aw00 + Aw01 + Aw10 + Aw11;
                            }
                        }
                    }
                    ///////////////
                    // process B //
                    ///////////////
                    if (!B.isEmpty()) {
                        add_EM_S=true;
                        const double* B_p=const_cast<escript::Data*>(&B)->getSampleDataRO(e);
                        for (index_t k=0; k<numEq; k++) {
                            for (index_t m=0; m<numComp; m++) {
                                const double wB0 = B_p[INDEX3(k,0,m, numEq, 2)]*w2;
                                const double wB1 = B_p[INDEX3(k,1,m, numEq, 2)]*w1;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=-wB0 - wB1;
                                EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=-wB0 - wB1;
                                EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=-wB0 - wB1;
                                EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=-wB0 - wB1;
                                EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+= wB0 - wB1;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+= wB0 - wB1;
                                EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+= wB0 - wB1;
                                EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+= wB0 - wB1;
                                EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=-wB0 + wB1;
                                EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=-wB0 + wB1;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=-wB0 + wB1;
                                EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=-wB0 + wB1;
                                EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+= wB0 + wB1;
                                EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+= wB0 + wB1;
                                EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+= wB0 + wB1;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+= wB0 + wB1;
                            }
                        }
                    }
                    ///////////////
                    // process C //
                    ///////////////
                    if (!C.isEmpty()) {
                        add_EM_S=true;
                        const double* C_p=const_cast<escript::Data*>(&C)->getSampleDataRO(e);
                        for (index_t k=0; k<numEq; k++) {
                            for (index_t m=0; m<numComp; m++) {
                                const double wC0 = C_p[INDEX3(k, m, 0, numEq, numComp)]*w2;
                                const double wC1 = C_p[INDEX3(k, m, 1, numEq, numComp)]*w1;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=-wC0 - wC1;
                                EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=-wC0 - wC1;
                                EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=-wC0 - wC1;
                                EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=-wC0 - wC1;
                                EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+= wC0 - wC1;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+= wC0 - wC1;
                                EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+= wC0 - wC1;
                                EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+= wC0 - wC1;
                                EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=-wC0 + wC1;
                                EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=-wC0 + wC1;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=-wC0 + wC1;
                                EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=-wC0 + wC1;
                                EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+= wC0 + wC1;
                                EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+= wC0 + wC1;
                                EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+= wC0 + wC1;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+= wC0 + wC1;
                            }
                        }
                    }
                    ///////////////
                    // process D //
                    ///////////////
                    if (!D.isEmpty()) {
                        add_EM_S=true;
                        const double* D_p=const_cast<escript::Data*>(&D)->getSampleDataRO(e);
                        for (index_t k=0; k<numEq; k++) {
                            for (index_t m=0; m<numComp; m++) {
                                const double wD0 = D_p[INDEX2(k, m, numEq)]*w3;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=wD0;
                            }
                        }
                    }
                    ///////////////
                    // process X //
                    ///////////////
                    if (!X.isEmpty()) {
                        add_EM_F=true;
                        const double* X_p=const_cast<escript::Data*>(&X)->getSampleDataRO(e);
                        for (index_t k=0; k<numEq; k++) {
                            const double wX0 = 4*X_p[INDEX2(k, 0, numEq)]*w2;
                            const double wX1 = 4*X_p[INDEX2(k, 1, numEq)]*w1;
                            EM_F[INDEX2(k,0,numEq)]+=-wX0 - wX1;
                            EM_F[INDEX2(k,1,numEq)]+= wX0 - wX1;
                            EM_F[INDEX2(k,2,numEq)]+=-wX0 + wX1;
                            EM_F[INDEX2(k,3,numEq)]+= wX0 + wX1;
                        }
                    }
                    ///////////////
                    // process Y //
                    ///////////////
                    if (!Y.isEmpty()) {
                        add_EM_F=true;
                        const double* Y_p=const_cast<escript::Data*>(&Y)->getSampleDataRO(e);
                        for (index_t k=0; k<numEq; k++) {
                            EM_F[INDEX2(k,0,numEq)]+=4*Y_p[k]*w3;
                            EM_F[INDEX2(k,1,numEq)]+=4*Y_p[k]*w3;
                            EM_F[INDEX2(k,2,numEq)]+=4*Y_p[k]*w3;
                            EM_F[INDEX2(k,3,numEq)]+=4*Y_p[k]*w3;
                        }
                    }

                    // add to matrix (if add_EM_S) and RHS (if add_EM_F)
                    const index_t firstNode=m_NN[0]*k1+k0;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F,
                            firstNode, numEq, numComp);
                } // end k0 loop
            } // end k1 loop
        } // end of colouring
    } // end of parallel region
}

//protected
void Rectangle::assemblePDEBoundarySingle(Paso_SystemMatrix* mat,
      escript::Data& rhs, const escript::Data& d, const escript::Data& y) const
{
    const double SQRT3 = 1.73205080756887719318;
    const double w5 = m_dx[0]/12;
    const double w6 = w5*(SQRT3 + 2);
    const double w7 = w5*(-SQRT3 + 2);
    const double w8 = w5*(SQRT3 + 3);
    const double w9 = w5*(-SQRT3 + 3);
    const double w2 = m_dx[1]/12;
    const double w0 = w2*(SQRT3 + 2);
    const double w1 = w2*(-SQRT3 + 2);
    const double w3 = w2*(SQRT3 + 3);
    const double w4 = w2*(-SQRT3 + 3);
    const bool add_EM_S=!d.isEmpty();
    const bool add_EM_F=!y.isEmpty();
    rhs.requireWrite();
#pragma omp parallel
    {
        if (m_faceOffset[0] > -1) {
            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
                for (index_t k1=k1_0; k1<m_NE[1]; k1+=2) {
                    vector<double> EM_S(4*4, 0);
                    vector<double> EM_F(4, 0);
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(k1);
                        if (d.actsExpanded()) {
                            const double d_0 = d_p[0];
                            const double d_1 = d_p[1];
                            const double tmp0 = w2*(d_0 + d_1);
                            EM_S[INDEX2(0,0,4)]+=d_0*w0 + d_1*w1;
                            EM_S[INDEX2(2,0,4)]+=tmp0;
                            EM_S[INDEX2(0,2,4)]+=tmp0;
                            EM_S[INDEX2(2,2,4)]+=d_0*w1 + d_1*w0;
                        } else { // constant data
                            EM_S[INDEX2(0,0,4)]+=4*d_p[0]*w2;
                            EM_S[INDEX2(2,0,4)]+=2*d_p[0]*w2;
                            EM_S[INDEX2(0,2,4)]+=2*d_p[0]*w2;
                            EM_S[INDEX2(2,2,4)]+=4*d_p[0]*w2;
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(k1);
                        if (y.actsExpanded()) {
                            EM_F[0]+=w3*y_p[0] + w4*y_p[1];
                            EM_F[2]+=w3*y_p[1] + w4*y_p[0];
                        } else { // constant data
                            EM_F[0]+=6*w2*y_p[0];
                            EM_F[2]+=6*w2*y_p[0];
                        }
                    }
                    const index_t firstNode=m_NN[0]*k1;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F, firstNode);
                }
            } // end colouring
        }

        if (m_faceOffset[1] > -1) {
            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring            
#pragma omp for
                for (index_t k1=k1_0; k1<m_NE[1]; k1+=2) {
                    vector<double> EM_S(4*4, 0);
                    vector<double> EM_F(4, 0);
                    const index_t e = m_faceOffset[1]+k1;
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(e);
                        if (d.actsExpanded()) {
                            const double d_0 = d_p[0];
                            const double d_1 = d_p[1];
                            const double tmp0 = w2*(d_0 + d_1);
                            EM_S[INDEX2(1,1,4)]+=d_0*w0 + d_1*w1;
                            EM_S[INDEX2(3,1,4)]+=tmp0;
                            EM_S[INDEX2(1,3,4)]+=tmp0;
                            EM_S[INDEX2(3,3,4)]+=d_0*w1 + d_1*w0;
                        } else { // constant data
                            EM_S[INDEX2(1,1,4)]+=4*d_p[0]*w2;
                            EM_S[INDEX2(3,1,4)]+=2*d_p[0]*w2;
                            EM_S[INDEX2(1,3,4)]+=2*d_p[0]*w2;
                            EM_S[INDEX2(3,3,4)]+=4*d_p[0]*w2;
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        if (y.actsExpanded()) {
                            EM_F[1]+=w3*y_p[0] + w4*y_p[1];
                            EM_F[3]+=w3*y_p[1] + w4*y_p[0];
                        } else { // constant data
                            EM_F[1]+=6*w2*y_p[0];
                            EM_F[3]+=6*w2*y_p[0];
                        }
                    }
                    const index_t firstNode=m_NN[0]*(k1+1)-2;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F, firstNode);
                }
            } // end colouring
        }

        if (m_faceOffset[2] > -1) {
            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < m_NE[0]; k0+=2) {
                    vector<double> EM_S(4*4, 0);
                    vector<double> EM_F(4, 0);
                    const index_t e = m_faceOffset[2]+k0;
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(e);
                        if (d.actsExpanded()) {
                            const double d_0 = d_p[0];
                            const double d_1 = d_p[1];
                            const double tmp0 = w5*(d_0 + d_1);
                            EM_S[INDEX2(0,0,4)]+=d_0*w6 + d_1*w7;
                            EM_S[INDEX2(1,0,4)]+=tmp0;
                            EM_S[INDEX2(0,1,4)]+=tmp0;
                            EM_S[INDEX2(1,1,4)]+=d_0*w7 + d_1*w6;
                        } else { // constant data
                            EM_S[INDEX2(0,0,4)]+=4*d_p[0]*w5;
                            EM_S[INDEX2(1,0,4)]+=2*d_p[0]*w5;
                            EM_S[INDEX2(0,1,4)]+=2*d_p[0]*w5;
                            EM_S[INDEX2(1,1,4)]+=4*d_p[0]*w5;
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        if (y.actsExpanded()) {
                            EM_F[0]+=w8*y_p[0] + w9*y_p[1];
                            EM_F[1]+=w8*y_p[1] + w9*y_p[0];
                        } else { // constant data
                            EM_F[0]+=6*w5*y_p[0];
                            EM_F[1]+=6*w5*y_p[0];
                        }
                    }
                    const index_t firstNode=k0;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F, firstNode);
                }
            } // end colouring
        }

        if (m_faceOffset[3] > -1) {
            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < m_NE[0]; k0+=2) {
                    const index_t e = m_faceOffset[3]+k0;
                    vector<double> EM_S(4*4, 0);
                    vector<double> EM_F(4, 0);
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(e);
                        if (d.actsExpanded()) {
                            const double d_0 = d_p[0];
                            const double d_1 = d_p[1];
                            const double tmp0 = w5*(d_0 + d_1);
                            EM_S[INDEX2(2,2,4)]+=d_0*w6 + d_1*w7;
                            EM_S[INDEX2(3,2,4)]+=tmp0;
                            EM_S[INDEX2(2,3,4)]+=tmp0;
                            EM_S[INDEX2(3,3,4)]+=d_0*w7 + d_1*w6;
                        } else { // constant data
                            EM_S[INDEX2(2,2,4)]+=4*d_p[0]*w5;
                            EM_S[INDEX2(3,2,4)]+=2*d_p[0]*w5;
                            EM_S[INDEX2(2,3,4)]+=2*d_p[0]*w5;
                            EM_S[INDEX2(3,3,4)]+=4*d_p[0]*w5;
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        if (y.actsExpanded()) {
                            EM_F[2]+=w8*y_p[0] + w9*y_p[1];
                            EM_F[3]+=w8*y_p[1] + w9*y_p[0];
                        } else { // constant data
                            EM_F[2]+=6*w5*y_p[0];
                            EM_F[3]+=6*w5*y_p[0];
                        }
                    }
                    const index_t firstNode=m_NN[0]*(m_NN[1]-2)+k0;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F, firstNode);
                }
            } // end colouring
        }
    } // end of parallel section
}

//protected
void Rectangle::assemblePDEBoundarySingleReduced(Paso_SystemMatrix* mat,
      escript::Data& rhs, const escript::Data& d, const escript::Data& y) const
{
    const double w0 = m_dx[0]/4;
    const double w1 = m_dx[1]/4;
    const bool add_EM_S=!d.isEmpty();
    const bool add_EM_F=!y.isEmpty();
    rhs.requireWrite();
#pragma omp parallel
    {
        if (m_faceOffset[0] > -1) {
            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
                for (index_t k1=k1_0; k1<m_NE[1]; k1+=2) {
                    vector<double> EM_S(4*4, 0);
                    vector<double> EM_F(4, 0);
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(k1);
                        EM_S[INDEX2(0,0,4)]+=d_p[0]*w1;
                        EM_S[INDEX2(2,0,4)]+=d_p[0]*w1;
                        EM_S[INDEX2(0,2,4)]+=d_p[0]*w1;
                        EM_S[INDEX2(2,2,4)]+=d_p[0]*w1;
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(k1);
                        EM_F[0]+=2*w1*y_p[0];
                        EM_F[2]+=2*w1*y_p[0];
                    }
                    const index_t firstNode=m_NN[0]*k1;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F, firstNode);
                }
            } // end colouring
        }

        if (m_faceOffset[1] > -1) {
            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring            
#pragma omp for
                for (index_t k1=k1_0; k1<m_NE[1]; k1+=2) {
                    vector<double> EM_S(4*4, 0);
                    vector<double> EM_F(4, 0);
                    const index_t e = m_faceOffset[1]+k1;
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(e);
                        EM_S[INDEX2(1,1,4)]+=d_p[0]*w1;
                        EM_S[INDEX2(3,1,4)]+=d_p[0]*w1;
                        EM_S[INDEX2(1,3,4)]+=d_p[0]*w1;
                        EM_S[INDEX2(3,3,4)]+=d_p[0]*w1;
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        EM_F[1]+=2*w1*y_p[0];
                        EM_F[3]+=2*w1*y_p[0];
                    }
                    const index_t firstNode=m_NN[0]*(k1+1)-2;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F, firstNode);
                }
            } // end colouring
        }

        if (m_faceOffset[2] > -1) {
            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < m_NE[0]; k0+=2) {
                    vector<double> EM_S(4*4, 0);
                    vector<double> EM_F(4, 0);
                    const index_t e = m_faceOffset[2]+k0;
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(e);
                        EM_S[INDEX2(0,0,4)]+=d_p[0]*w0;
                        EM_S[INDEX2(1,0,4)]+=d_p[0]*w0;
                        EM_S[INDEX2(0,1,4)]+=d_p[0]*w0;
                        EM_S[INDEX2(1,1,4)]+=d_p[0]*w0;
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        EM_F[0]+=2*w0*y_p[0];
                        EM_F[1]+=2*w0*y_p[0];
                    }
                    const index_t firstNode=k0;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F, firstNode);
                }
            } // end colouring
        }

        if (m_faceOffset[3] > -1) {
            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < m_NE[0]; k0+=2) {
                    vector<double> EM_S(4*4, 0);
                    vector<double> EM_F(4, 0);
                    const index_t e = m_faceOffset[3]+k0;
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(e);
                        EM_S[INDEX2(2,2,4)]+=d_p[0]*w0;
                        EM_S[INDEX2(3,2,4)]+=d_p[0]*w0;
                        EM_S[INDEX2(2,3,4)]+=d_p[0]*w0;
                        EM_S[INDEX2(3,3,4)]+=d_p[0]*w0;
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        EM_F[2]+=2*w0*y_p[0];
                        EM_F[3]+=2*w0*y_p[0];
                    }
                    const index_t firstNode=m_NN[0]*(m_NN[1]-2)+k0;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F, firstNode);
                }
            } // end colouring
        }
    } // end of parallel section
}

//protected
void Rectangle::assemblePDEBoundarySystem(Paso_SystemMatrix* mat,
      escript::Data& rhs, const escript::Data& d, const escript::Data& y) const
{
    dim_t numEq, numComp;
    if (!mat) {
        numEq=numComp=(rhs.isEmpty() ? 1 : rhs.getDataPointSize());
    } else {
        numEq=mat->logical_row_block_size;
        numComp=mat->logical_col_block_size;
    }
    const double SQRT3 = 1.73205080756887719318;
    const double w5 = m_dx[0]/12;
    const double w6 = w5*(SQRT3 + 2);
    const double w7 = w5*(-SQRT3 + 2);
    const double w8 = w5*(SQRT3 + 3);
    const double w9 = w5*(-SQRT3 + 3);
    const double w2 = m_dx[1]/12;
    const double w0 = w2*(SQRT3 + 2);
    const double w1 = w2*(-SQRT3 + 2);
    const double w3 = w2*(SQRT3 + 3);
    const double w4 = w2*(-SQRT3 + 3);
    const bool add_EM_S=!d.isEmpty();
    const bool add_EM_F=!y.isEmpty();
    rhs.requireWrite();
#pragma omp parallel
    {
        if (m_faceOffset[0] > -1) {
            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
                for (index_t k1=k1_0; k1<m_NE[1]; k1+=2) {
                    vector<double> EM_S(4*4*numEq*numComp, 0);
                    vector<double> EM_F(4*numEq, 0);
                    const index_t e = k1;
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(e);
                        if (d.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                    const double d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                    const double tmp0 = w2*(d_0 + d_1);
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=d_0*w0 + d_1*w1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp0;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp0;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=d_0*w1 + d_1*w0;
                                }
                             }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX2(k, m, numEq)];
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=4*d_0*w2;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=2*d_0*w2;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=2*d_0*w2;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=4*d_0*w2;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        if (y.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                const double y_0 = y_p[INDEX2(k, 0, numEq)];
                                const double y_1 = y_p[INDEX2(k, 1, numEq)];
                                EM_F[INDEX2(k,0,numEq)]+=w3*y_0 + w4*y_1;
                                EM_F[INDEX2(k,2,numEq)]+=w3*y_1 + w4*y_0;
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                EM_F[INDEX2(k,0,numEq)]+=6*w2*y_p[k];
                                EM_F[INDEX2(k,2,numEq)]+=6*w2*y_p[k];
                            }
                        }
                    }
                    const index_t firstNode=m_NN[0]*k1;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F,
                            firstNode, numEq, numComp);
                }
            } // end colouring
        }

        if (m_faceOffset[1] > -1) {
            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring            
#pragma omp for
                for (index_t k1=k1_0; k1<m_NE[1]; k1+=2) {
                    vector<double> EM_S(4*4*numEq*numComp, 0);
                    vector<double> EM_F(4*numEq, 0);
                    const index_t e = m_faceOffset[1]+k1;
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(e);
                        if (d.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                    const double d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                    const double tmp0 = w2*(d_0 + d_1);
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=d_0*w0 + d_1*w1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp0;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp0;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=d_0*w1 + d_1*w0;
                                }
                             }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX2(k, m, numEq)];
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=4*d_0*w2;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=2*d_0*w2;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=2*d_0*w2;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=4*d_0*w2;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        if (y.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                const double y_0 = y_p[INDEX2(k, 0, numEq)];
                                const double y_1 = y_p[INDEX2(k, 1, numEq)];
                                EM_F[INDEX2(k,1,numEq)]+=w3*y_0 + w4*y_1;
                                EM_F[INDEX2(k,3,numEq)]+=w3*y_1 + w4*y_0;
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                EM_F[INDEX2(k,1,numEq)]+=6*w2*y_p[k];
                                EM_F[INDEX2(k,3,numEq)]+=6*w2*y_p[k];
                            }
                        }
                    }
                    const index_t firstNode=m_NN[0]*(k1+1)-2;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F,
                            firstNode, numEq, numComp);
                }
            } // end colouring
        }

        if (m_faceOffset[2] > -1) {
            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < m_NE[0]; k0+=2) {
                    vector<double> EM_S(4*4*numEq*numComp, 0);
                    vector<double> EM_F(4*numEq, 0);
                    const index_t e = m_faceOffset[2]+k0;
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(e);
                        if (d.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                    const double d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                    const double tmp0 = w5*(d_0 + d_1);
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=d_0*w6 + d_1*w7;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp0;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=d_0*w7 + d_1*w6;
                                }
                             }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX2(k, m, numEq)];
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=4*d_0*w5;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=2*d_0*w5;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=2*d_0*w5;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=4*d_0*w5;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        if (y.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                const double y_0 = y_p[INDEX2(k, 0, numEq)];
                                const double y_1 = y_p[INDEX2(k, 1, numEq)];
                                EM_F[INDEX2(k,0,numEq)]+=w8*y_0 + w9*y_1;
                                EM_F[INDEX2(k,1,numEq)]+=w8*y_1 + w9*y_0;
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                EM_F[INDEX2(k,0,numEq)]+=6*w5*y_p[k];
                                EM_F[INDEX2(k,1,numEq)]+=6*w5*y_p[k];
                            }
                        }
                    }
                    const index_t firstNode=k0;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F,
                            firstNode, numEq, numComp);
                }
            } // end colouring
        }

        if (m_faceOffset[3] > -1) {
            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < m_NE[0]; k0+=2) {
                    vector<double> EM_S(4*4*numEq*numComp, 0);
                    vector<double> EM_F(4*numEq, 0);
                    const index_t e = m_faceOffset[3]+k0;
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(e);
                        if (d.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                    const double d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                    const double tmp0 = w5*(d_0 + d_1);
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=d_0*w6 + d_1*w7;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp0;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp0;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=d_0*w7 + d_1*w6;
                                }
                             }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX2(k, m, numEq)];
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=4*d_0*w5;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=2*d_0*w5;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=2*d_0*w5;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=4*d_0*w5;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        if (y.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                const double y_0 = y_p[INDEX2(k, 0, numEq)];
                                const double y_1 = y_p[INDEX2(k, 1, numEq)];
                                EM_F[INDEX2(k,2,numEq)]+=w8*y_0 + w9*y_1;
                                EM_F[INDEX2(k,3,numEq)]+=w8*y_1 + w9*y_0;
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                EM_F[INDEX2(k,2,numEq)]+=6*w5*y_p[k];
                                EM_F[INDEX2(k,3,numEq)]+=6*w5*y_p[k];
                            }
                        }
                    }
                    const index_t firstNode=m_NN[0]*(m_NN[1]-2)+k0;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F,
                            firstNode, numEq, numComp);
                }
            } // end colouring
        }
    } // end of parallel section
}

//protected
void Rectangle::assemblePDEBoundarySystemReduced(Paso_SystemMatrix* mat,
      escript::Data& rhs, const escript::Data& d, const escript::Data& y) const
{
    dim_t numEq, numComp;
    if (!mat)
        numEq=numComp=(rhs.isEmpty() ? 1 : rhs.getDataPointSize());
    else {
        numEq=mat->logical_row_block_size;
        numComp=mat->logical_col_block_size;
    }
    const double w0 = m_dx[0]/4;
    const double w1 = m_dx[1]/4;
    const bool add_EM_S=!d.isEmpty();
    const bool add_EM_F=!y.isEmpty();
    rhs.requireWrite();
#pragma omp parallel
    {
        if (m_faceOffset[0] > -1) {
            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
                for (index_t k1=k1_0; k1<m_NE[1]; k1+=2) {
                    vector<double> EM_S(4*4*numEq*numComp, 0);
                    vector<double> EM_F(4*numEq, 0);
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(k1);
                        for (index_t k=0; k<numEq; k++) {
                            for (index_t m=0; m<numComp; m++) {
                                const double tmp0 = d_p[INDEX2(k, m, numEq)]*w1;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp0;
                                EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp0;
                                EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp0;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp0;
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(k1);
                        for (index_t k=0; k<numEq; k++) {
                            EM_F[INDEX2(k,0,numEq)]+=2*w1*y_p[k];
                            EM_F[INDEX2(k,2,numEq)]+=2*w1*y_p[k];
                        }
                    }
                    const index_t firstNode=m_NN[0]*k1;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F,
                            firstNode, numEq, numComp);
                }
            } // end colouring
        }

        if (m_faceOffset[1] > -1) {
            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring            
#pragma omp for
                for (index_t k1=k1_0; k1<m_NE[1]; k1+=2) {
                    vector<double> EM_S(4*4*numEq*numComp, 0);
                    vector<double> EM_F(4*numEq, 0);
                    const index_t e = m_faceOffset[1]+k1;
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(e);
                        for (index_t k=0; k<numEq; k++) {
                            for (index_t m=0; m<numComp; m++) {
                                const double tmp0 = d_p[INDEX2(k, m, numEq)]*w1;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp0;
                                EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp0;
                                EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp0;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp0;
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        for (index_t k=0; k<numEq; k++) {
                            EM_F[INDEX2(k,1,numEq)]+=2*w1*y_p[k];
                            EM_F[INDEX2(k,3,numEq)]+=2*w1*y_p[k];
                        }
                    }
                    const index_t firstNode=m_NN[0]*(k1+1)-2;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F,
                            firstNode, numEq, numComp);
                }
            } // end colouring
        }

        if (m_faceOffset[2] > -1) {
            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < m_NE[0]; k0+=2) {
                    vector<double> EM_S(4*4*numEq*numComp, 0);
                    vector<double> EM_F(4*numEq, 0);
                    const index_t e = m_faceOffset[2]+k0;
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(e);
                        for (index_t k=0; k<numEq; k++) {
                            for (index_t m=0; m<numComp; m++) {
                                const double tmp0 = d_p[INDEX2(k, m, numEq)]*w0;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp0;
                                EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp0;
                                EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp0;
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        for (index_t k=0; k<numEq; k++) {
                            EM_F[INDEX2(k,0,numEq)]+=2*w0*y_p[k];
                            EM_F[INDEX2(k,1,numEq)]+=2*w0*y_p[k];
                        }
                    }
                    const index_t firstNode=k0;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F,
                            firstNode, numEq, numComp);
                }
            } // end colouring
        }

        if (m_faceOffset[3] > -1) {
            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < m_NE[0]; k0+=2) {
                    vector<double> EM_S(4*4*numEq*numComp, 0);
                    vector<double> EM_F(4*numEq, 0);
                    const index_t e = m_faceOffset[3]+k0;
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(e);
                        for (index_t k=0; k<numEq; k++) {
                            for (index_t m=0; m<numComp; m++) {
                                const double tmp0 = d_p[INDEX2(k, m, numEq)]*w0;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp0;
                                EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp0;
                                EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp0;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp0;
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        for (index_t k=0; k<numEq; k++) {
                            EM_F[INDEX2(k,2,numEq)]+=2*w0*y_p[k];
                            EM_F[INDEX2(k,3,numEq)]+=2*w0*y_p[k];
                        }
                    }
                    const index_t firstNode=m_NN[0]*(m_NN[1]-2)+k0;
                    addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F,
                            firstNode, numEq, numComp);
                }
            } // end colouring
        }
    } // end of parallel section
}

namespace
{
    // Calculates a guassian blur colvolution matrix for 2D
    double* get2DGauss(unsigned radius, double sigma)
    {
        double* arr=new double[(radius*2+1)*(radius*2+1)];
        double common=M_1_PI*0.5*1/(sigma*sigma);
	double total=0;
	int r=static_cast<int>(radius);
	for (int y=-r;y<=r;++y)
	{
	    for (int x=-r;x<=r;++x)
	    {	      
	        arr[(x+r)+(y+r)*(r*2+1)]=common*exp(-(x*x+y*y)/(2*sigma*sigma));
// cout << (x+y*(r*2+1)) << " " << arr[(x+r)+(y+r)*(r*2+1)] << endl;
	        total+=arr[(x+r)+(y+r)*(r*2+1)];
	    }
	}
	double invtotal=1/total;
//cout << "Inv total is "	 << invtotal << endl;
	for (size_t p=0;p<(radius*2+1)*(radius*2+1);++p)
	{
	    arr[p]*=invtotal; 
//cout << p << "->" << arr[p] << endl;	    
	}
	return arr;
    }
    
    // applies conv to source to get a point.
    // (xp, yp) are the coords in the source matrix not the destination matrix
    double Convolve2D(double* conv, double* source, size_t xp, size_t yp, unsigned radius, size_t width)
    {
        size_t bx=xp-radius, by=yp-radius;
	size_t sbase=bx+by*width;
	double result=0;
	for (int y=0;y<2*radius+1;++y)
	{	  
	    for (int x=0;x<2*radius+1;++x)
	    {
	        result+=conv[x+y*(2*radius+1)] * source[sbase + x+y*width];
	    }
	}
        return result;      
    }
}



escript::Data Rectangle::randomFill(long seed, const boost::python::tuple& filter) const
{
//     if (m_mpiInfo->size!=1)
//     {
//         throw RipleyException("This type of random does not support MPI yet.");
//     }
    if (m_numDim!=2)
    {
        throw RipleyException("Only 2D supported at this time.");
    }
    if (len(filter)!=3) {
        throw RipleyException("Unsupported random filter");
    }
    boost::python::extract<string> ex(filter[0]);
    if (!ex.check() || (ex()!="gaussian")) 
    {
        throw RipleyException("Unsupported random filter");
    }
    boost::python::extract<unsigned int> ex1(filter[1]);
    if (!ex1.check())
    {
        throw RipleyException("Radius of gaussian filter must be a positive integer.");
    }
    unsigned int radius=ex1();
#ifdef ESYS_MPI    

    // Need to check to see that radius would not cause the overlap to cover 
    // more than one cell (eg each rank holds 4 columns and the radius is 5).
    // Also need to take special care with narrow cells
    
    
    // In fact it needs to be stricter than this, if a rank has neighbours on both sides, the borders can't overlap.
    
#endif    
    double sigma=0.5;
    boost::python::extract<double> ex2(filter[2]);
    if (!ex2.check() || (sigma=ex2())<=0)
    {
        throw RipleyException("Sigma must be a postive floating point number.");
    }    
    size_t numpoints[2];
    numpoints[0]=m_ownNE[0]+1;
    numpoints[1]=m_ownNE[1]+1;
    size_t padding=max((unsigned)max((m_NE[0]-m_ownNE[0])/2, (m_NE[1]-m_ownNE[1])/2), radius);
    size_t width=(numpoints[0]+2*padding);  	// width of one row in points
    size_t height=(numpoints[1]+2*padding);	// height of one row in points
    size_t dsize=width * height; // size of padded source grid 
    
    double* src=new double[dsize];
    esysUtils::randomFillArray(seed, src, dsize);  
      
    // Now we need to copy the regions owned by other ranks over here  
#ifdef ESYS_MPI    
    
    
    dim_t X=m_mpiInfo->rank%m_NX[0];
    dim_t Y=m_mpiInfo->rank/m_NX[0];
    dim_t row=m_NX[0];

    MPI_Request reqs[10];
    MPI_Status stats[10];
    short rused=0;
    double* SWin=new double[radius*radius];  memset(SWin, 0, radius*radius*sizeof(double));
    double* SEin=new double[radius*radius];  memset(SEin, 0, radius*radius*sizeof(double));
    double* NWin=new double[radius*radius];  memset(NWin, 0, radius*radius*sizeof(double));
    double* Sin=new double[radius*numpoints[0]];  memset(Sin, 0, radius*numpoints[0]*sizeof(double));
    double* Win=new double[radius*numpoints[1]];  memset(Win, 0, radius*numpoints[1]*sizeof(double));

    double* NEout=new double[radius*radius];  memset(NEout, 0, radius*radius*sizeof(double));
    double* NWout=new double[radius*radius];  memset(NWout, 0, radius*radius*sizeof(double));
    double* SEout=new double[radius*radius];  memset(SEout, 0, radius*radius*sizeof(double));
    double* Nout=new double[radius*numpoints[0]];  memset(Nout, 0, radius*numpoints[0]*sizeof(double));
    double* Eout=new double[radius*numpoints[1]];  memset(Eout, 0, radius*numpoints[1]*sizeof(double));
    
//cerr << m_mpiInfo->rank << "[" << __LINE__ << "]: " << SWin << " " << SEin << " " << NWin << " "<< Sin << " "<< Win << " "<< NEout << " "<< NWout << " "
//<< SEout << " "
//<< Nout << " "
//<< Eout << " "
//<< endl;  	    
    

    int comserr=0;
    if (Y!=0)	// not on bottom row, 
    {
	if (X!=0)	// not on the left hand edge
	{
	    // recv bottomleft from SW
//cerr << m_mpiInfo->rank << "[" << __LINE__ << "]: " << "Recv SW (7) from " << (X-1)+(Y-1)*row << endl;
	    comserr|=MPI_Irecv(SWin, radius*radius, MPI_DOUBLE, (X-1)+(Y-1)*row, 7, m_mpiInfo->comm, reqs+(rused++));
//cerr << m_mpiInfo->rank << "[" << __LINE__ << "]: " << "Recv W (10) from " <<  X-1+Y*row << endl;  	    
	    comserr|=MPI_Irecv(Win, numpoints[1]*radius, MPI_DOUBLE, X-1+Y*row, 10, m_mpiInfo->comm, reqs+(rused++));
    
	}
	else
	{
//cerr << m_mpiInfo->rank << "[" << __LINE__ << "]: " << "Recv SW (7) from " << (Y-1)*row  << endl;  	    	  
	    comserr|=MPI_Irecv(SWin, radius*radius, MPI_DOUBLE, (Y-1)*row, 7, m_mpiInfo->comm, reqs+(rused++));
	}
//cerr << m_mpiInfo->rank << "[" << __LINE__ << "]: " << "Recv S (8) from " << X+(Y-1)*row  << endl;  		
	comserr|=MPI_Irecv(Sin, numpoints[0]*radius, MPI_DOUBLE, X+(Y-1)*row, 8, m_mpiInfo->comm, reqs+(rused++));
//cerr << m_mpiInfo->rank << "[" << __LINE__ << "]: " << "Recv SE (7) from " << X+(Y-1)*row << endl;  	
	comserr|=MPI_Irecv(SEin, radius*radius, MPI_DOUBLE, X+(Y-1)*row, 7, m_mpiInfo->comm, reqs+(rused++));

      
    }
    else		// on the bottom row
    {
	if (X!=0) 
	{
//cerr << m_mpiInfo->rank << "[" << __LINE__ << "]: " << "Recv W (10) from " << X-1+Y*row << endl;	  
	    comserr|=MPI_Irecv(Win, numpoints[1]*radius, MPI_DOUBLE, X-1+Y*row, 10, m_mpiInfo->comm, reqs+(rused++));
//cerr << m_mpiInfo->rank << "[" << __LINE__ << "]: " << "Recv NW (7) from " << X-1+Y*row << endl;	    
	    comserr|=MPI_Irecv(NWin, radius*radius, MPI_DOUBLE, X-1+Y*row, 7, m_mpiInfo->comm, reqs+(rused++));
	}
	if (X!=(row-1))
	{
//cerr << m_mpiInfo->rank << "[" << __LINE__ << "]: " << "Send SE (7) to " <<  X+1+(Y)*row << endl; 		  
	    comserr|=MPI_Isend(SEout, radius*radius, MPI_DOUBLE, X+1+(Y)*row, 7, m_mpiInfo->comm, reqs+(rused++));	
    
	}
    }
    
    if (Y!=(m_NX[1]-1))	// not on the top row
    {
//cerr << "Y=" << Y << "  (numpoints[1]-1)=" << (numpoints[1]-1) << endl; 
//cerr << m_mpiInfo->rank << "[" << __LINE__ << "]: " << "Send N (8) to " << X+(Y+1)*row  << endl;
 	comserr|=MPI_Isend(Nout, radius*numpoints[0], MPI_DOUBLE, X+(Y+1)*row, 8, m_mpiInfo->comm, reqs+(rused++));
//cerr << m_mpiInfo->rank << "[" << __LINE__ << "]: " << "Send NE (7) to " <<  X+(Y+1)*row << endl; 
	comserr|=MPI_Isend(NEout, radius*radius, MPI_DOUBLE, X+(Y+1)*row, 7, m_mpiInfo->comm, reqs+(rused++));
	if (X!=(row-1))	// not on right hand edge
	{
//cerr << m_mpiInfo->rank << "[" << __LINE__ << "]: " << "Send NE (7) to " <<  X+1+(Y+1)*row << endl; 	    	  
	    comserr|=MPI_Isend(NEout, radius*radius, MPI_DOUBLE, X+1+(Y+1)*row, 7, m_mpiInfo->comm, reqs+(rused++));
	}
	if (X==0)	// left hand edge
	{
//cerr << m_mpiInfo->rank << "[" << __LINE__ << "]: " << "Send NW (7) to " << (Y+1)*row << endl; 	  
	    comserr|=MPI_Isend(NWout, radius*radius, MPI_DOUBLE, (Y+1)*row,7, m_mpiInfo->comm, reqs+(rused++));	    
	}	
    }
    if (X!=(row-1))	// not on right hand edge
    {
//cerr << m_mpiInfo->rank << "[" << __LINE__ << "]: " << "Send NE (7) to " << X+1+(Y)*row << endl;       
	comserr|=MPI_Isend(NEout, radius*radius, MPI_DOUBLE, X+1+(Y)*row, 7, m_mpiInfo->comm, reqs+(rused++));
//cerr << m_mpiInfo->rank << "[" << __LINE__ << "]: " << "Send E (10) to " << X+1+(Y)*row << endl; 	
	comserr|=MPI_Isend(Eout, numpoints[1]*radius, MPI_DOUBLE, X+1+(Y)*row, 10, m_mpiInfo->comm, reqs+(rused++));
    }
    if (X!=0)
    {
//cerr << m_mpiInfo->rank << "[" << __LINE__ << "]: " << "Recv NW (7) from " << (X-1)+Y*row << endl;  	
      
	comserr|=MPI_Irecv(NWin, radius*radius, MPI_DOUBLE, (X-1)+Y*row, 7, m_mpiInfo->comm, reqs+(rused++));
      
      
    }
    
    if (!comserr)
    {
//cerr << rused << ": " <<   m_mpiInfo->rank << "[" << __LINE__ << "]\n";        
        comserr=MPI_Waitall(rused, reqs, stats);
    }

    if (comserr)
    {
	// Yes this is throwing an exception as a result of an MPI error.
	// and no we don't inform the other ranks that we are doing this.
	// however, we have no reason to believe coms work at this point anyway
        throw RipleyException("Error in coms for randomFill");      
    }
    
    
    delete[] SWin;
    delete[] SEin;
    delete[] NWin;
    delete[] Sin;
    delete[] Win;

    delete[] NEout;
    delete[] NWout;
    delete[] SEout;
    delete[] Nout;
    delete[] Eout;
    
    
    
    
#endif    
    // Lets call that done for now
    escript::FunctionSpace fs(getPtr(), getContinuousFunctionCode());
    escript::Data resdat(0, escript::DataTypes::scalarShape, fs , true);
    // don't need to check for exwrite because we just made it
    escript::DataVector& dv=resdat.getExpandedVectorReference();
    double* convolution=get2DGauss(radius, sigma);
    for (size_t y=0;y<(m_ownNE[1]+1);++y)    
    {
        for (size_t x=0;x<(m_ownNE[0]+1);++x)
	{	  
	    dv[x+y*(m_ownNE[0]+1)]=Convolve2D(convolution, src, x+radius, y+radius, radius, width);
	}
    }
    delete[] convolution;
    return resdat;
}




} // end of namespace ripley

