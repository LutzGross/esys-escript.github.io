
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
    const long *edges = var->edges();

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
                               const vector<int>& multiplier) const
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
    const int reqsize = numValues[0]*numValues[1]*numComp*sizeof(float);
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
    vector<float> values(num0*numComp);
    const int dpp = out.getNumDataPointsPerSample();

    for (index_t y=0; y<num1; y++) {
        const int fileofs = numComp*(idx0+(idx1+y)*numValues[0]);
        f.seekg(fileofs*sizeof(float));
        f.read((char*)&values[0], num0*numComp*sizeof(float));
        for (index_t x=0; x<num0; x++) {
            const int baseIndex = first0+x*multiplier[0]
                                    +(first1+y*multiplier[1])*myN0;
            for (index_t m1=0; m1<multiplier[1]; m1++) {
                for (index_t m0=0; m0<multiplier[0]; m0++) {
                    const int dataIndex = baseIndex+m0+m1*myN0;
                    double* dest = out.getSampleDataRW(dataIndex);
                    for (index_t c=0; c<numComp; c++) {
                        if (!std::isnan(values[x*numComp+c])) {
                            for (index_t q=0; q<dpp; q++) {
                                *dest++ = static_cast<double>(values[x*numComp+c]);
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
    FileWriter* fw = new FileWriter();
    fw->openFile(filename, fileSize);
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
        fw->writeAt(oss, fileofs);
    }
    fw->close();
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
    const double cx0 = -1./m_dx[0];
    const double cx1 = -.78867513459481288225/m_dx[0];
    const double cx2 = -.5/m_dx[0];
    const double cx3 = -.21132486540518711775/m_dx[0];
    const double cx4 = .21132486540518711775/m_dx[0];
    const double cx5 = .5/m_dx[0];
    const double cx6 = .78867513459481288225/m_dx[0];
    const double cx7 = 1./m_dx[0];
    const double cy0 = -1./m_dx[1];
    const double cy1 = -.78867513459481288225/m_dx[1];
    const double cy2 = -.5/m_dx[1];
    const double cy3 = -.21132486540518711775/m_dx[1];
    const double cy4 = .21132486540518711775/m_dx[1];
    const double cy5 = .5/m_dx[1];
    const double cy6 = .78867513459481288225/m_dx[1];
    const double cy7 = 1./m_dx[1];

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
                        o[INDEX3(i,0,0,numComp,2)] = f_00[i]*cx1 + f_01[i]*cx3 + f_10[i]*cx6 + f_11[i]*cx4;
                        o[INDEX3(i,1,0,numComp,2)] = f_00[i]*cy1 + f_01[i]*cy6 + f_10[i]*cy3 + f_11[i]*cy4;
                        o[INDEX3(i,0,1,numComp,2)] = f_00[i]*cx1 + f_01[i]*cx3 + f_10[i]*cx6 + f_11[i]*cx4;
                        o[INDEX3(i,1,1,numComp,2)] = f_00[i]*cy3 + f_01[i]*cy4 + f_10[i]*cy1 + f_11[i]*cy6;
                        o[INDEX3(i,0,2,numComp,2)] = f_00[i]*cx3 + f_01[i]*cx1 + f_10[i]*cx4 + f_11[i]*cx6;
                        o[INDEX3(i,1,2,numComp,2)] = f_00[i]*cy1 + f_01[i]*cy6 + f_10[i]*cy3 + f_11[i]*cy4;
                        o[INDEX3(i,0,3,numComp,2)] = f_00[i]*cx3 + f_01[i]*cx1 + f_10[i]*cx4 + f_11[i]*cx6;
                        o[INDEX3(i,1,3,numComp,2)] = f_00[i]*cy3 + f_01[i]*cy4 + f_10[i]*cy1 + f_11[i]*cy6;
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
                        o[INDEX3(i,0,0,numComp,2)] = cx5*(f_10[i] + f_11[i]) + cx2*(f_00[i] + f_01[i]);
                        o[INDEX3(i,1,0,numComp,2)] = cy2*(f_00[i] + f_10[i]) + cy5*(f_01[i] + f_11[i]);
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
                        o[INDEX3(i,0,0,numComp,2)] = f_00[i]*cx1 + f_01[i]*cx3 + f_10[i]*cx6 + f_11[i]*cx4;
                        o[INDEX3(i,1,0,numComp,2)] = f_00[i]*cy0 + f_01[i]*cy7;
                        o[INDEX3(i,0,1,numComp,2)] = f_00[i]*cx3 + f_01[i]*cx1 + f_10[i]*cx4 + f_11[i]*cx6;
                        o[INDEX3(i,1,1,numComp,2)] = f_00[i]*cy0 + f_01[i]*cy7;
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
                        o[INDEX3(i,0,0,numComp,2)] = f_00[i]*cx1 + f_01[i]*cx3 + f_10[i]*cx6 + f_11[i]*cx4;
                        o[INDEX3(i,1,0,numComp,2)] = f_10[i]*cy0 + f_11[i]*cy7;
                        o[INDEX3(i,0,1,numComp,2)] = f_00[i]*cx3 + f_01[i]*cx1 + f_10[i]*cx4 + f_11[i]*cx6;
                        o[INDEX3(i,1,1,numComp,2)] = f_10[i]*cy0 + f_11[i]*cy7;
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
                        o[INDEX3(i,0,0,numComp,2)] = f_00[i]*cx0 + f_10[i]*cx7;
                        o[INDEX3(i,1,0,numComp,2)] = f_00[i]*cy1 + f_01[i]*cy6 + f_10[i]*cy3 + f_11[i]*cy4;
                        o[INDEX3(i,0,1,numComp,2)] = f_00[i]*cx0 + f_10[i]*cx7;
                        o[INDEX3(i,1,1,numComp,2)] = f_00[i]*cy3 + f_01[i]*cy4 + f_10[i]*cy1 + f_11[i]*cy6;
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
                        o[INDEX3(i,0,0,numComp,2)] = f_01[i]*cx0 + f_11[i]*cx7;
                        o[INDEX3(i,1,0,numComp,2)] = f_00[i]*cy1 + f_01[i]*cy6 + f_10[i]*cy3 + f_11[i]*cy4;
                        o[INDEX3(i,0,1,numComp,2)] = f_01[i]*cx0 + f_11[i]*cx7;
                        o[INDEX3(i,1,1,numComp,2)] = f_00[i]*cy3 + f_01[i]*cy4 + f_10[i]*cy1 + f_11[i]*cy6;
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
                        o[INDEX3(i,0,0,numComp,2)] = cx5*(f_10[i] + f_11[i]) + cx2*(f_00[i] + f_01[i]);
                        o[INDEX3(i,1,0,numComp,2)] = f_00[i]*cy0 + f_01[i]*cy7;
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
                        o[INDEX3(i,0,0,numComp,2)] = cx5*(f_10[i] + f_11[i]) + cx2*(f_00[i] + f_01[i]);
                        o[INDEX3(i,1,0,numComp,2)] = f_10[i]*cy0 + f_11[i]*cy7;
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
                        o[INDEX3(i,0,0,numComp,2)] = f_00[i]*cx0 + f_10[i]*cx7;
                        o[INDEX3(i,1,0,numComp,2)] = cy2*(f_00[i] + f_10[i]) + cy5*(f_01[i] + f_11[i]);
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
                        o[INDEX3(i,0,0,numComp,2)] = f_01[i]*cx0 + f_11[i]*cx7;
                        o[INDEX3(i,1,0,numComp,2)] = cy5*(f_01[i] + f_11[i]) + cy2*(f_00[i] + f_10[i]);
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
        const double c0 = 0.5;
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
                        o[INDEX2(i,numComp,0)] = c0*(f_00[i] + f_01[i]);
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
                        o[INDEX2(i,numComp,0)] = c0*(f_10[i] + f_11[i]);
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
                        o[INDEX2(i,numComp,0)] = c0*(f_00[i] + f_10[i]);
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
                        o[INDEX2(i,numComp,0)] = c0*(f_01[i] + f_11[i]);
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
    const double w0 = -0.1555021169820365539*m_dx[1]/m_dx[0];
    const double w1 = 0.041666666666666666667;
    const double w2 = -0.15550211698203655390;
    const double w3 = 0.041666666666666666667*m_dx[0]/m_dx[1];
    const double w4 = 0.15550211698203655390;
    const double w5 = -0.041666666666666666667;
    const double w6 = -0.01116454968463011277*m_dx[1]/m_dx[0];
    const double w7 = 0.011164549684630112770;
    const double w8 = -0.011164549684630112770;
    const double w9 = -0.041666666666666666667*m_dx[1]/m_dx[0];
    const double w10 = -0.041666666666666666667*m_dx[0]/m_dx[1];
    const double w11 = 0.1555021169820365539*m_dx[1]/m_dx[0];
    const double w12 = 0.1555021169820365539*m_dx[0]/m_dx[1];
    const double w13 = 0.01116454968463011277*m_dx[0]/m_dx[1];
    const double w14 = 0.01116454968463011277*m_dx[1]/m_dx[0];
    const double w15 = 0.041666666666666666667*m_dx[1]/m_dx[0];
    const double w16 = -0.01116454968463011277*m_dx[0]/m_dx[1];
    const double w17 = -0.1555021169820365539*m_dx[0]/m_dx[1];
    const double w18 = -0.33333333333333333333*m_dx[1]/m_dx[0];
    const double w19 = 0.25;
    const double w20 = -0.25;
    const double w21 = 0.16666666666666666667*m_dx[0]/m_dx[1];
    const double w22 = -0.16666666666666666667*m_dx[1]/m_dx[0];
    const double w23 = -0.16666666666666666667*m_dx[0]/m_dx[1];
    const double w24 = 0.33333333333333333333*m_dx[1]/m_dx[0];
    const double w25 = 0.33333333333333333333*m_dx[0]/m_dx[1];
    const double w26 = 0.16666666666666666667*m_dx[1]/m_dx[0];
    const double w27 = -0.33333333333333333333*m_dx[0]/m_dx[1];
    const double w28 = -0.032861463941450536761*m_dx[1];
    const double w29 = -0.032861463941450536761*m_dx[0];
    const double w30 = -0.12264065304058601714*m_dx[1];
    const double w31 = -0.0023593469594139828636*m_dx[1];
    const double w32 = -0.008805202725216129906*m_dx[0];
    const double w33 = -0.008805202725216129906*m_dx[1];
    const double w34 = 0.032861463941450536761*m_dx[1];
    const double w35 = 0.008805202725216129906*m_dx[1];
    const double w36 = 0.008805202725216129906*m_dx[0];
    const double w37 = 0.0023593469594139828636*m_dx[1];
    const double w38 = 0.12264065304058601714*m_dx[1];
    const double w39 = 0.032861463941450536761*m_dx[0];
    const double w40 = -0.12264065304058601714*m_dx[0];
    const double w41 = -0.0023593469594139828636*m_dx[0];
    const double w42 = 0.0023593469594139828636*m_dx[0];
    const double w43 = 0.12264065304058601714*m_dx[0];
    const double w44 = -0.16666666666666666667*m_dx[1];
    const double w45 = -0.083333333333333333333*m_dx[0];
    const double w46 = 0.083333333333333333333*m_dx[1];
    const double w47 = 0.16666666666666666667*m_dx[1];
    const double w48 = 0.083333333333333333333*m_dx[0];
    const double w49 = -0.16666666666666666667*m_dx[0];
    const double w50 = 0.16666666666666666667*m_dx[0];
    const double w51 = -0.083333333333333333333*m_dx[1];
    const double w52 = 0.025917019497006092316*m_dx[0]*m_dx[1];
    const double w53 = 0.0018607582807716854616*m_dx[0]*m_dx[1];
    const double w54 = 0.0069444444444444444444*m_dx[0]*m_dx[1];
    const double w55 = 0.09672363354357992482*m_dx[0]*m_dx[1];
    const double w56 = 0.00049858867864229740201*m_dx[0]*m_dx[1];
    const double w57 = 0.055555555555555555556*m_dx[0]*m_dx[1];
    const double w58 = 0.027777777777777777778*m_dx[0]*m_dx[1];
    const double w59 = 0.11111111111111111111*m_dx[0]*m_dx[1];
    const double w60 = -0.19716878364870322056*m_dx[1];
    const double w61 = -0.19716878364870322056*m_dx[0];
    const double w62 = -0.052831216351296779436*m_dx[0];
    const double w63 = -0.052831216351296779436*m_dx[1];
    const double w64 = 0.19716878364870322056*m_dx[1];
    const double w65 = 0.052831216351296779436*m_dx[1];
    const double w66 = 0.19716878364870322056*m_dx[0];
    const double w67 = 0.052831216351296779436*m_dx[0];
    const double w68 = -0.5*m_dx[1];
    const double w69 = -0.5*m_dx[0];
    const double w70 = 0.5*m_dx[1];
    const double w71 = 0.5*m_dx[0];
    const double w72 = 0.1555021169820365539*m_dx[0]*m_dx[1];
    const double w73 = 0.041666666666666666667*m_dx[0]*m_dx[1];
    const double w74 = 0.01116454968463011277*m_dx[0]*m_dx[1];
    const double w75 = 0.25*m_dx[0]*m_dx[1];

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
                        if (A.actsExpanded()) {
                            const double A_00_0 = A_p[INDEX3(0,0,0,2,2)];
                            const double A_10_0 = A_p[INDEX3(1,0,0,2,2)];
                            const double A_01_0 = A_p[INDEX3(0,1,0,2,2)];
                            const double A_11_0 = A_p[INDEX3(1,1,0,2,2)];
                            const double A_00_1 = A_p[INDEX3(0,0,1,2,2)];
                            const double A_10_1 = A_p[INDEX3(1,0,1,2,2)];
                            const double A_01_1 = A_p[INDEX3(0,1,1,2,2)];
                            const double A_11_1 = A_p[INDEX3(1,1,1,2,2)];
                            const double A_00_2 = A_p[INDEX3(0,0,2,2,2)];
                            const double A_10_2 = A_p[INDEX3(1,0,2,2,2)];
                            const double A_01_2 = A_p[INDEX3(0,1,2,2,2)];
                            const double A_11_2 = A_p[INDEX3(1,1,2,2,2)];
                            const double A_00_3 = A_p[INDEX3(0,0,3,2,2)];
                            const double A_10_3 = A_p[INDEX3(1,0,3,2,2)];
                            const double A_01_3 = A_p[INDEX3(0,1,3,2,2)];
                            const double A_11_3 = A_p[INDEX3(1,1,3,2,2)];
                            const double tmp0_0 = A_01_0 + A_01_3;
                            const double tmp1_0 = A_00_0 + A_00_1;
                            const double tmp2_0 = A_11_0 + A_11_1 + A_11_2 + A_11_3;
                            const double tmp3_0 = A_00_2 + A_00_3;
                            const double tmp4_0 = A_10_1 + A_10_2;
                            const double tmp5_0 = A_00_0 + A_00_1 + A_00_2 + A_00_3;
                            const double tmp6_0 = A_01_3 + A_10_0;
                            const double tmp7_0 = A_01_0 + A_10_3;
                            const double tmp8_0 = A_01_1 + A_01_2 + A_10_1 + A_10_2;
                            const double tmp9_0 = A_01_0 + A_10_0;
                            const double tmp12_0 = A_11_0 + A_11_2;
                            const double tmp10_0 = A_01_3 + A_10_3;
                            const double tmp14_0 = A_01_0 + A_01_3 + A_10_0 + A_10_3;
                            const double tmp13_0 = A_01_2 + A_10_1;
                            const double tmp11_0 = A_11_1 + A_11_3;
                            const double tmp18_0 = A_01_1 + A_10_1;
                            const double tmp15_0 = A_01_1 + A_10_2;
                            const double tmp16_0 = A_10_0 + A_10_3;
                            const double tmp17_0 = A_01_1 + A_01_2;
                            const double tmp19_0 = A_01_2 + A_10_2;
                            const double tmp0_1 = A_10_3*w8;
                            const double tmp1_1 = tmp0_0*w1;
                            const double tmp2_1 = A_01_1*w4;
                            const double tmp3_1 = tmp1_0*w0;
                            const double tmp4_1 = A_01_2*w7;
                            const double tmp5_1 = tmp2_0*w3;
                            const double tmp6_1 = tmp3_0*w6;
                            const double tmp7_1 = A_10_0*w2;
                            const double tmp8_1 = tmp4_0*w5;
                            const double tmp9_1 = tmp2_0*w10;
                            const double tmp14_1 = A_10_0*w8;
                            const double tmp23_1 = tmp3_0*w14;
                            const double tmp35_1 = A_01_0*w8;
                            const double tmp54_1 = tmp13_0*w8;
                            const double tmp20_1 = tmp9_0*w4;
                            const double tmp25_1 = tmp12_0*w12;
                            const double tmp44_1 = tmp7_0*w7;
                            const double tmp26_1 = tmp10_0*w4;
                            const double tmp52_1 = tmp18_0*w8;
                            const double tmp48_1 = A_10_1*w7;
                            const double tmp46_1 = A_01_3*w8;
                            const double tmp50_1 = A_01_0*w2;
                            const double tmp56_1 = tmp19_0*w8;
                            const double tmp19_1 = A_10_3*w2;
                            const double tmp47_1 = A_10_2*w4;
                            const double tmp16_1 = tmp3_0*w0;
                            const double tmp18_1 = tmp1_0*w6;
                            const double tmp31_1 = tmp11_0*w12;
                            const double tmp55_1 = tmp15_0*w2;
                            const double tmp39_1 = A_10_2*w7;
                            const double tmp11_1 = tmp6_0*w7;
                            const double tmp40_1 = tmp11_0*w17;
                            const double tmp34_1 = tmp15_0*w8;
                            const double tmp33_1 = tmp14_0*w5;
                            const double tmp24_1 = tmp11_0*w13;
                            const double tmp43_1 = tmp17_0*w5;
                            const double tmp15_1 = A_01_2*w4;
                            const double tmp53_1 = tmp19_0*w2;
                            const double tmp27_1 = tmp3_0*w11;
                            const double tmp32_1 = tmp13_0*w2;
                            const double tmp10_1 = tmp5_0*w9;
                            const double tmp37_1 = A_10_1*w4;
                            const double tmp38_1 = tmp5_0*w15;
                            const double tmp17_1 = A_01_1*w7;
                            const double tmp12_1 = tmp7_0*w4;
                            const double tmp22_1 = tmp10_0*w7;
                            const double tmp57_1 = tmp18_0*w2;
                            const double tmp28_1 = tmp9_0*w7;
                            const double tmp29_1 = tmp1_0*w14;
                            const double tmp51_1 = tmp11_0*w16;
                            const double tmp42_1 = tmp12_0*w16;
                            const double tmp49_1 = tmp12_0*w17;
                            const double tmp21_1 = tmp1_0*w11;
                            const double tmp45_1 = tmp6_0*w4;
                            const double tmp13_1 = tmp8_0*w1;
                            const double tmp36_1 = tmp16_0*w1;
                            const double tmp41_1 = A_01_3*w2;
                            const double tmp30_1 = tmp12_0*w13;
                            EM_S[INDEX2(0,0,4)]+=tmp13_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp24_1 + tmp25_1;
                            EM_S[INDEX2(1,0,4)]+=tmp36_1 + tmp37_1 + tmp39_1 + tmp3_1 + tmp43_1 + tmp46_1 + tmp50_1 + tmp5_1 + tmp6_1;
                            EM_S[INDEX2(2,0,4)]+=tmp0_1 + tmp15_1 + tmp17_1 + tmp1_1 + tmp38_1 + tmp49_1 + tmp51_1 + tmp7_1 + tmp8_1;
                            EM_S[INDEX2(3,0,4)]+=tmp10_1 + tmp32_1 + tmp33_1 + tmp34_1 + tmp9_1;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1;
                            EM_S[INDEX2(1,1,4)]+=tmp21_1 + tmp23_1 + tmp30_1 + tmp31_1 + tmp33_1 + tmp56_1 + tmp57_1;
                            EM_S[INDEX2(2,1,4)]+=tmp10_1 + tmp13_1 + tmp44_1 + tmp45_1 + tmp9_1;
                            EM_S[INDEX2(3,1,4)]+=tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp43_1;
                            EM_S[INDEX2(0,2,4)]+=tmp36_1 + tmp38_1 + tmp43_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                            EM_S[INDEX2(1,2,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp9_1;
                            EM_S[INDEX2(2,2,4)]+=tmp24_1 + tmp25_1 + tmp27_1 + tmp29_1 + tmp33_1 + tmp52_1 + tmp53_1;
                            EM_S[INDEX2(3,2,4)]+=tmp14_1 + tmp15_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp1_1 + tmp5_1 + tmp8_1;
                            EM_S[INDEX2(0,3,4)]+=tmp10_1 + tmp33_1 + tmp54_1 + tmp55_1 + tmp9_1;
                            EM_S[INDEX2(1,3,4)]+=tmp14_1 + tmp19_1 + tmp1_1 + tmp2_1 + tmp38_1 + tmp40_1 + tmp42_1 + tmp4_1 + tmp8_1;
                            EM_S[INDEX2(2,3,4)]+=tmp16_1 + tmp18_1 + tmp35_1 + tmp36_1 + tmp41_1 + tmp43_1 + tmp47_1 + tmp48_1 + tmp5_1;
                            EM_S[INDEX2(3,3,4)]+=tmp13_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1;
                        } else { // constant data
                            const double A_00 = A_p[INDEX2(0,0,2)];
                            const double A_10 = A_p[INDEX2(1,0,2)];
                            const double A_01 = A_p[INDEX2(0,1,2)];
                            const double A_11 = A_p[INDEX2(1,1,2)];
                            const double tmp0_0 = A_01 + A_10;
                            const double tmp0_1 = A_00*w18;
                            const double tmp1_1 = A_01*w19;
                            const double tmp2_1 = A_10*w20;
                            const double tmp3_1 = A_11*w21;
                            const double tmp4_1 = A_00*w22;
                            const double tmp5_1 = tmp0_0*w19;
                            const double tmp6_1 = A_11*w23;
                            const double tmp7_1 = A_11*w25;
                            const double tmp8_1 = A_00*w24;
                            const double tmp9_1 = tmp0_0*w20;
                            const double tmp10_1 = A_01*w20;
                            const double tmp11_1 = A_11*w27;
                            const double tmp12_1 = A_00*w26;
                            const double tmp13_1 = A_10*w19;
                            EM_S[INDEX2(0,0,4)]+=tmp5_1 + tmp7_1 + tmp8_1;
                            EM_S[INDEX2(1,0,4)]+=tmp0_1 + tmp10_1 + tmp13_1 + tmp3_1;
                            EM_S[INDEX2(2,0,4)]+=tmp11_1 + tmp12_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(3,0,4)]+=tmp4_1 + tmp6_1 + tmp9_1;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(1,1,4)]+=tmp7_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(2,1,4)]+=tmp4_1 + tmp5_1 + tmp6_1;
                            EM_S[INDEX2(3,1,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1;
                            EM_S[INDEX2(0,2,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1;
                            EM_S[INDEX2(1,2,4)]+=tmp4_1 + tmp5_1 + tmp6_1;
                            EM_S[INDEX2(2,2,4)]+=tmp7_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(3,2,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(0,3,4)]+=tmp4_1 + tmp6_1 + tmp9_1;
                            EM_S[INDEX2(1,3,4)]+=tmp11_1 + tmp12_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(2,3,4)]+=tmp0_1 + tmp10_1 + tmp13_1 + tmp3_1;
                            EM_S[INDEX2(3,3,4)]+=tmp5_1 + tmp7_1 + tmp8_1;
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
                            const double tmp0_0 = B_1_0 + B_1_1;
                            const double tmp1_0 = B_1_2 + B_1_3;
                            const double tmp2_0 = B_0_1 + B_0_3;
                            const double tmp3_0 = B_0_0 + B_0_2;
                            const double tmp63_1 = B_1_1*w42;
                            const double tmp79_1 = B_1_1*w40;
                            const double tmp37_1 = tmp3_0*w35;
                            const double tmp8_1 = tmp0_0*w32;
                            const double tmp71_1 = B_0_1*w34;
                            const double tmp19_1 = B_0_3*w31;
                            const double tmp15_1 = B_0_3*w34;
                            const double tmp9_1 = tmp3_0*w34;
                            const double tmp35_1 = B_1_0*w36;
                            const double tmp66_1 = B_0_3*w28;
                            const double tmp28_1 = B_1_0*w42;
                            const double tmp22_1 = B_1_0*w40;
                            const double tmp16_1 = B_1_2*w29;
                            const double tmp6_1 = tmp2_0*w35;
                            const double tmp55_1 = B_1_3*w40;
                            const double tmp50_1 = B_1_3*w42;
                            const double tmp7_1 = tmp1_0*w29;
                            const double tmp1_1 = tmp1_0*w32;
                            const double tmp57_1 = B_0_3*w30;
                            const double tmp18_1 = B_1_1*w32;
                            const double tmp53_1 = B_1_0*w41;
                            const double tmp61_1 = B_1_3*w36;
                            const double tmp27_1 = B_0_3*w38;
                            const double tmp64_1 = B_0_2*w30;
                            const double tmp76_1 = B_0_1*w38;
                            const double tmp39_1 = tmp2_0*w34;
                            const double tmp62_1 = B_0_1*w31;
                            const double tmp56_1 = B_0_0*w31;
                            const double tmp49_1 = B_1_1*w36;
                            const double tmp2_1 = B_0_2*w31;
                            const double tmp23_1 = B_0_2*w33;
                            const double tmp38_1 = B_1_1*w43;
                            const double tmp74_1 = B_1_2*w41;
                            const double tmp43_1 = B_1_1*w41;
                            const double tmp58_1 = B_0_2*w28;
                            const double tmp67_1 = B_0_0*w33;
                            const double tmp33_1 = tmp0_0*w39;
                            const double tmp4_1 = B_0_0*w28;
                            const double tmp20_1 = B_0_0*w30;
                            const double tmp13_1 = B_0_2*w38;
                            const double tmp65_1 = B_1_2*w43;
                            const double tmp0_1 = tmp0_0*w29;
                            const double tmp41_1 = tmp3_0*w33;
                            const double tmp73_1 = B_0_2*w37;
                            const double tmp69_1 = B_0_0*w38;
                            const double tmp48_1 = B_1_2*w39;
                            const double tmp59_1 = B_0_1*w33;
                            const double tmp17_1 = B_1_3*w41;
                            const double tmp5_1 = B_0_3*w33;
                            const double tmp3_1 = B_0_1*w30;
                            const double tmp21_1 = B_0_1*w28;
                            const double tmp42_1 = B_1_0*w29;
                            const double tmp54_1 = B_1_2*w32;
                            const double tmp60_1 = B_1_0*w39;
                            const double tmp32_1 = tmp1_0*w36;
                            const double tmp10_1 = B_0_1*w37;
                            const double tmp14_1 = B_0_0*w35;
                            const double tmp29_1 = B_0_1*w35;
                            const double tmp26_1 = B_1_2*w36;
                            const double tmp30_1 = B_1_3*w43;
                            const double tmp70_1 = B_0_2*w35;
                            const double tmp34_1 = B_1_3*w39;
                            const double tmp51_1 = B_1_0*w43;
                            const double tmp31_1 = B_0_2*w34;
                            const double tmp45_1 = tmp3_0*w28;
                            const double tmp11_1 = tmp1_0*w39;
                            const double tmp52_1 = B_1_1*w29;
                            const double tmp44_1 = B_1_3*w32;
                            const double tmp25_1 = B_1_1*w39;
                            const double tmp47_1 = tmp2_0*w33;
                            const double tmp72_1 = B_1_3*w29;
                            const double tmp40_1 = tmp2_0*w28;
                            const double tmp46_1 = B_1_2*w40;
                            const double tmp36_1 = B_1_2*w42;
                            const double tmp24_1 = B_0_0*w37;
                            const double tmp77_1 = B_0_3*w35;
                            const double tmp68_1 = B_0_3*w37;
                            const double tmp78_1 = B_0_0*w34;
                            const double tmp12_1 = tmp0_0*w36;
                            const double tmp75_1 = B_1_0*w32;
                            EM_S[INDEX2(0,0,4)]+=tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1;
                            EM_S[INDEX2(1,0,4)]+=tmp0_1 + tmp1_1 + tmp68_1 + tmp69_1 + tmp70_1 + tmp71_1;
                            EM_S[INDEX2(2,0,4)]+=tmp45_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                            EM_S[INDEX2(3,0,4)]+=tmp32_1 + tmp33_1 + tmp6_1 + tmp9_1;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(1,1,4)]+=tmp72_1 + tmp73_1 + tmp74_1 + tmp75_1 + tmp76_1 + tmp77_1 + tmp78_1 + tmp79_1;
                            EM_S[INDEX2(2,1,4)]+=tmp32_1 + tmp33_1 + tmp40_1 + tmp41_1;
                            EM_S[INDEX2(3,1,4)]+=tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1;
                            EM_S[INDEX2(0,2,4)]+=tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1;
                            EM_S[INDEX2(1,2,4)]+=tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(2,2,4)]+=tmp60_1 + tmp61_1 + tmp62_1 + tmp63_1 + tmp64_1 + tmp65_1 + tmp66_1 + tmp67_1;
                            EM_S[INDEX2(3,2,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1;
                            EM_S[INDEX2(0,3,4)]+=tmp40_1 + tmp41_1 + tmp7_1 + tmp8_1;
                            EM_S[INDEX2(1,3,4)]+=tmp37_1 + tmp39_1 + tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1;
                            EM_S[INDEX2(2,3,4)]+=tmp11_1 + tmp12_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1;
                            EM_S[INDEX2(3,3,4)]+=tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1;
                        } else { // constant data
                            const double B_0 = B_p[0];
                            const double B_1 = B_p[1];
                            const double tmp0_1 = B_0*w44;
                            const double tmp1_1 = B_1*w45;
                            const double tmp2_1 = B_0*w46;
                            const double tmp3_1 = B_0*w47;
                            const double tmp4_1 = B_1*w48;
                            const double tmp5_1 = B_1*w49;
                            const double tmp6_1 = B_1*w50;
                            const double tmp7_1 = B_0*w51;
                            EM_S[INDEX2(0,0,4)]+=tmp0_1 + tmp5_1;
                            EM_S[INDEX2(1,0,4)]+=tmp1_1 + tmp3_1;
                            EM_S[INDEX2(2,0,4)]+=tmp6_1 + tmp7_1;
                            EM_S[INDEX2(3,0,4)]+=tmp2_1 + tmp4_1;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1;
                            EM_S[INDEX2(1,1,4)]+=tmp3_1 + tmp5_1;
                            EM_S[INDEX2(2,1,4)]+=tmp4_1 + tmp7_1;
                            EM_S[INDEX2(3,1,4)]+=tmp2_1 + tmp6_1;
                            EM_S[INDEX2(0,2,4)]+=tmp5_1 + tmp7_1;
                            EM_S[INDEX2(1,2,4)]+=tmp1_1 + tmp2_1;
                            EM_S[INDEX2(2,2,4)]+=tmp0_1 + tmp6_1;
                            EM_S[INDEX2(3,2,4)]+=tmp3_1 + tmp4_1;
                            EM_S[INDEX2(0,3,4)]+=tmp1_1 + tmp7_1;
                            EM_S[INDEX2(1,3,4)]+=tmp2_1 + tmp5_1;
                            EM_S[INDEX2(2,3,4)]+=tmp0_1 + tmp4_1;
                            EM_S[INDEX2(3,3,4)]+=tmp3_1 + tmp6_1;
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
                            const double tmp0_0 = C_1_0 + C_1_1;
                            const double tmp1_0 = C_1_2 + C_1_3;
                            const double tmp2_0 = C_0_1 + C_0_3;
                            const double tmp3_0 = C_0_0 + C_0_2;
                            const double tmp64_1 = C_0_2*w30;
                            const double tmp14_1 = C_0_2*w28;
                            const double tmp19_1 = C_0_3*w31;
                            const double tmp22_1 = C_1_0*w40;
                            const double tmp37_1 = tmp3_0*w35;
                            const double tmp29_1 = C_0_1*w35;
                            const double tmp73_1 = C_0_2*w37;
                            const double tmp74_1 = C_1_2*w41;
                            const double tmp52_1 = C_1_3*w39;
                            const double tmp25_1 = C_1_1*w39;
                            const double tmp62_1 = C_0_1*w31;
                            const double tmp79_1 = C_1_1*w40;
                            const double tmp43_1 = C_1_1*w36;
                            const double tmp27_1 = C_0_3*w38;
                            const double tmp28_1 = C_1_0*w42;
                            const double tmp63_1 = C_1_1*w42;
                            const double tmp59_1 = C_0_3*w34;
                            const double tmp72_1 = C_1_3*w29;
                            const double tmp40_1 = tmp2_0*w35;
                            const double tmp13_1 = C_0_3*w30;
                            const double tmp51_1 = C_1_2*w40;
                            const double tmp54_1 = C_1_2*w42;
                            const double tmp12_1 = C_0_0*w31;
                            const double tmp2_1 = tmp1_0*w32;
                            const double tmp68_1 = C_0_2*w31;
                            const double tmp75_1 = C_1_0*w32;
                            const double tmp49_1 = C_1_1*w41;
                            const double tmp4_1 = C_0_2*w35;
                            const double tmp66_1 = C_0_3*w28;
                            const double tmp56_1 = C_0_1*w37;
                            const double tmp5_1 = C_0_1*w34;
                            const double tmp38_1 = tmp2_0*w34;
                            const double tmp76_1 = C_0_1*w38;
                            const double tmp21_1 = C_0_1*w28;
                            const double tmp69_1 = C_0_1*w30;
                            const double tmp53_1 = C_1_0*w36;
                            const double tmp42_1 = C_1_2*w39;
                            const double tmp32_1 = tmp1_0*w29;
                            const double tmp45_1 = C_1_0*w43;
                            const double tmp33_1 = tmp0_0*w32;
                            const double tmp35_1 = C_1_0*w41;
                            const double tmp26_1 = C_1_2*w36;
                            const double tmp67_1 = C_0_0*w33;
                            const double tmp31_1 = C_0_2*w34;
                            const double tmp20_1 = C_0_0*w30;
                            const double tmp70_1 = C_0_0*w28;
                            const double tmp8_1 = tmp0_0*w39;
                            const double tmp30_1 = C_1_3*w43;
                            const double tmp0_1 = tmp0_0*w29;
                            const double tmp17_1 = C_1_3*w41;
                            const double tmp58_1 = C_0_0*w35;
                            const double tmp9_1 = tmp3_0*w33;
                            const double tmp61_1 = C_1_3*w36;
                            const double tmp41_1 = tmp3_0*w34;
                            const double tmp50_1 = C_1_3*w32;
                            const double tmp18_1 = C_1_1*w32;
                            const double tmp6_1 = tmp1_0*w36;
                            const double tmp3_1 = C_0_0*w38;
                            const double tmp34_1 = C_1_1*w29;
                            const double tmp77_1 = C_0_3*w35;
                            const double tmp65_1 = C_1_2*w43;
                            const double tmp71_1 = C_0_3*w33;
                            const double tmp55_1 = C_1_1*w43;
                            const double tmp46_1 = tmp3_0*w28;
                            const double tmp24_1 = C_0_0*w37;
                            const double tmp10_1 = tmp1_0*w39;
                            const double tmp48_1 = C_1_0*w29;
                            const double tmp15_1 = C_0_1*w33;
                            const double tmp36_1 = C_1_2*w32;
                            const double tmp60_1 = C_1_0*w39;
                            const double tmp47_1 = tmp2_0*w33;
                            const double tmp16_1 = C_1_2*w29;
                            const double tmp1_1 = C_0_3*w37;
                            const double tmp7_1 = tmp2_0*w28;
                            const double tmp39_1 = C_1_3*w40;
                            const double tmp44_1 = C_1_3*w42;
                            const double tmp57_1 = C_0_2*w38;
                            const double tmp78_1 = C_0_0*w34;
                            const double tmp11_1 = tmp0_0*w36;
                            const double tmp23_1 = C_0_2*w33;
                            EM_S[INDEX2(0,0,4)]+=tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1;
                            EM_S[INDEX2(1,0,4)]+=tmp0_1 + tmp2_1 + tmp68_1 + tmp69_1 + tmp70_1 + tmp71_1;
                            EM_S[INDEX2(2,0,4)]+=tmp46_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                            EM_S[INDEX2(3,0,4)]+=tmp32_1 + tmp33_1 + tmp7_1 + tmp9_1;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(1,1,4)]+=tmp72_1 + tmp73_1 + tmp74_1 + tmp75_1 + tmp76_1 + tmp77_1 + tmp78_1 + tmp79_1;
                            EM_S[INDEX2(2,1,4)]+=tmp32_1 + tmp33_1 + tmp40_1 + tmp41_1;
                            EM_S[INDEX2(3,1,4)]+=tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1;
                            EM_S[INDEX2(0,2,4)]+=tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1;
                            EM_S[INDEX2(1,2,4)]+=tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(2,2,4)]+=tmp60_1 + tmp61_1 + tmp62_1 + tmp63_1 + tmp64_1 + tmp65_1 + tmp66_1 + tmp67_1;
                            EM_S[INDEX2(3,2,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1;
                            EM_S[INDEX2(0,3,4)]+=tmp40_1 + tmp41_1 + tmp6_1 + tmp8_1;
                            EM_S[INDEX2(1,3,4)]+=tmp37_1 + tmp38_1 + tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1;
                            EM_S[INDEX2(2,3,4)]+=tmp10_1 + tmp11_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1;
                            EM_S[INDEX2(3,3,4)]+=tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1;
                        } else { // constant data
                            const double C_0 = C_p[0];
                            const double C_1 = C_p[1];
                            const double tmp0_1 = C_0*w47;
                            const double tmp1_1 = C_1*w45;
                            const double tmp2_1 = C_1*w48;
                            const double tmp3_1 = C_0*w51;
                            const double tmp4_1 = C_0*w44;
                            const double tmp5_1 = C_1*w49;
                            const double tmp6_1 = C_1*w50;
                            const double tmp7_1 = C_0*w46;
                            EM_S[INDEX2(0,0,4)]+=tmp4_1 + tmp5_1;
                            EM_S[INDEX2(1,0,4)]+=tmp1_1 + tmp4_1;
                            EM_S[INDEX2(2,0,4)]+=tmp3_1 + tmp5_1;
                            EM_S[INDEX2(3,0,4)]+=tmp1_1 + tmp3_1;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1;
                            EM_S[INDEX2(1,1,4)]+=tmp0_1 + tmp5_1;
                            EM_S[INDEX2(2,1,4)]+=tmp1_1 + tmp7_1;
                            EM_S[INDEX2(3,1,4)]+=tmp5_1 + tmp7_1;
                            EM_S[INDEX2(0,2,4)]+=tmp3_1 + tmp6_1;
                            EM_S[INDEX2(1,2,4)]+=tmp2_1 + tmp3_1;
                            EM_S[INDEX2(2,2,4)]+=tmp4_1 + tmp6_1;
                            EM_S[INDEX2(3,2,4)]+=tmp2_1 + tmp4_1;
                            EM_S[INDEX2(0,3,4)]+=tmp2_1 + tmp7_1;
                            EM_S[INDEX2(1,3,4)]+=tmp6_1 + tmp7_1;
                            EM_S[INDEX2(2,3,4)]+=tmp0_1 + tmp2_1;
                            EM_S[INDEX2(3,3,4)]+=tmp0_1 + tmp6_1;
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
                            const double tmp0_0 = D_0 + D_1;
                            const double tmp1_0 = D_2 + D_3;
                            const double tmp2_0 = D_0 + D_1 + D_2 + D_3;
                            const double tmp3_0 = D_1 + D_2;
                            const double tmp4_0 = D_1 + D_3;
                            const double tmp5_0 = D_0 + D_2;
                            const double tmp6_0 = D_0 + D_3;
                            const double tmp0_1 = tmp0_0*w52;
                            const double tmp1_1 = tmp1_0*w53;
                            const double tmp2_1 = tmp2_0*w54;
                            const double tmp3_1 = tmp1_0*w52;
                            const double tmp4_1 = tmp0_0*w53;
                            const double tmp5_1 = tmp3_0*w54;
                            const double tmp6_1 = D_0*w55;
                            const double tmp7_1 = D_3*w56;
                            const double tmp8_1 = D_3*w55;
                            const double tmp9_1 = D_0*w56;
                            const double tmp10_1 = tmp4_0*w52;
                            const double tmp11_1 = tmp5_0*w53;
                            const double tmp12_1 = tmp5_0*w52;
                            const double tmp13_1 = tmp4_0*w53;
                            const double tmp14_1 = tmp6_0*w54;
                            const double tmp15_1 = D_2*w55;
                            const double tmp16_1 = D_1*w56;
                            const double tmp17_1 = D_1*w55;
                            const double tmp18_1 = D_2*w56;
                            EM_S[INDEX2(0,0,4)]+=tmp5_1 + tmp6_1 + tmp7_1;
                            EM_S[INDEX2(1,0,4)]+=tmp0_1 + tmp1_1;
                            EM_S[INDEX2(2,0,4)]+=tmp12_1 + tmp13_1;
                            EM_S[INDEX2(3,0,4)]+=tmp2_1;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1;
                            EM_S[INDEX2(1,1,4)]+=tmp14_1 + tmp17_1 + tmp18_1;
                            EM_S[INDEX2(2,1,4)]+=tmp2_1;
                            EM_S[INDEX2(3,1,4)]+=tmp10_1 + tmp11_1;
                            EM_S[INDEX2(0,2,4)]+=tmp12_1 + tmp13_1;
                            EM_S[INDEX2(1,2,4)]+=tmp2_1;
                            EM_S[INDEX2(2,2,4)]+=tmp14_1 + tmp15_1 + tmp16_1;
                            EM_S[INDEX2(3,2,4)]+=tmp3_1 + tmp4_1;
                            EM_S[INDEX2(0,3,4)]+=tmp2_1;
                            EM_S[INDEX2(1,3,4)]+=tmp10_1 + tmp11_1;
                            EM_S[INDEX2(2,3,4)]+=tmp3_1 + tmp4_1;
                            EM_S[INDEX2(3,3,4)]+=tmp5_1 + tmp8_1 + tmp9_1;
                        } else { // constant data
                            const double tmp0_1 = D_p[0]*w57;
                            const double tmp1_1 = D_p[0]*w58;
                            const double tmp2_1 = D_p[0]*w59;
                            EM_S[INDEX2(0,0,4)]+=tmp2_1;
                            EM_S[INDEX2(1,0,4)]+=tmp0_1;
                            EM_S[INDEX2(2,0,4)]+=tmp0_1;
                            EM_S[INDEX2(3,0,4)]+=tmp1_1;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1;
                            EM_S[INDEX2(1,1,4)]+=tmp2_1;
                            EM_S[INDEX2(2,1,4)]+=tmp1_1;
                            EM_S[INDEX2(3,1,4)]+=tmp0_1;
                            EM_S[INDEX2(0,2,4)]+=tmp0_1;
                            EM_S[INDEX2(1,2,4)]+=tmp1_1;
                            EM_S[INDEX2(2,2,4)]+=tmp2_1;
                            EM_S[INDEX2(3,2,4)]+=tmp0_1;
                            EM_S[INDEX2(0,3,4)]+=tmp1_1;
                            EM_S[INDEX2(1,3,4)]+=tmp0_1;
                            EM_S[INDEX2(2,3,4)]+=tmp0_1;
                            EM_S[INDEX2(3,3,4)]+=tmp2_1;
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
                            const double tmp0_0 = X_0_2 + X_0_3;
                            const double tmp1_0 = X_1_1 + X_1_3;
                            const double tmp2_0 = X_1_0 + X_1_2;
                            const double tmp3_0 = X_0_0 + X_0_1;
                            const double tmp0_1 = tmp0_0*w63;
                            const double tmp1_1 = tmp1_0*w62;
                            const double tmp2_1 = tmp2_0*w61;
                            const double tmp3_1 = tmp3_0*w60;
                            const double tmp4_1 = tmp0_0*w65;
                            const double tmp5_1 = tmp3_0*w64;
                            const double tmp6_1 = tmp2_0*w62;
                            const double tmp7_1 = tmp1_0*w61;
                            const double tmp8_1 = tmp2_0*w66;
                            const double tmp9_1 = tmp3_0*w63;
                            const double tmp10_1 = tmp0_0*w60;
                            const double tmp11_1 = tmp1_0*w67;
                            const double tmp12_1 = tmp1_0*w66;
                            const double tmp13_1 = tmp3_0*w65;
                            const double tmp14_1 = tmp0_0*w64;
                            const double tmp15_1 = tmp2_0*w67;
                            EM_F[0]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                            EM_F[1]+=tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1;
                            EM_F[2]+=tmp10_1 + tmp11_1 + tmp8_1 + tmp9_1;
                            EM_F[3]+=tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1;
                        } else { // constant data
                            const double tmp0_1 = X_p[1]*w69;
                            const double tmp1_1 = X_p[0]*w68;
                            const double tmp2_1 = X_p[0]*w70;
                            const double tmp3_1 = X_p[1]*w71;
                            EM_F[0]+=tmp0_1 + tmp1_1;
                            EM_F[1]+=tmp0_1 + tmp2_1;
                            EM_F[2]+=tmp1_1 + tmp3_1;
                            EM_F[3]+=tmp2_1 + tmp3_1;
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
                            const double tmp0_0 = Y_1 + Y_2;
                            const double tmp1_0 = Y_0 + Y_3;
                            const double tmp0_1 = Y_0*w72;
                            const double tmp1_1 = tmp0_0*w73;
                            const double tmp2_1 = Y_3*w74;
                            const double tmp3_1 = Y_1*w72;
                            const double tmp4_1 = tmp1_0*w73;
                            const double tmp5_1 = Y_2*w74;
                            const double tmp6_1 = Y_2*w72;
                            const double tmp7_1 = Y_1*w74;
                            const double tmp8_1 = Y_3*w72;
                            const double tmp9_1 = Y_0*w74;
                            EM_F[0]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_F[1]+=tmp3_1 + tmp4_1 + tmp5_1;
                            EM_F[2]+=tmp4_1 + tmp6_1 + tmp7_1;
                            EM_F[3]+=tmp1_1 + tmp8_1 + tmp9_1;
                        } else { // constant data
                            const double tmp0_1 = Y_p[0]*w75;
                            EM_F[0]+=tmp0_1;
                            EM_F[1]+=tmp0_1;
                            EM_F[2]+=tmp0_1;
                            EM_F[3]+=tmp0_1;
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
    const double w0 = -.25*m_dx[1]/m_dx[0];
    const double w1 = .25;
    const double w2 = -.25;
    const double w3 = .25*m_dx[0]/m_dx[1];
    const double w4 = -.25*m_dx[0]/m_dx[1];
    const double w5 = .25*m_dx[1]/m_dx[0];
    const double w6 = -.125*m_dx[1];
    const double w7 = -.125*m_dx[0];
    const double w8 = .125*m_dx[1];
    const double w9 = .125*m_dx[0];
    const double w10 = .0625*m_dx[0]*m_dx[1];
    const double w11 = -.5*m_dx[1];
    const double w12 = -.5*m_dx[0];
    const double w13 = .5*m_dx[1];
    const double w14 = .5*m_dx[0];
    const double w15 = .25*m_dx[0]*m_dx[1];

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
                        const double tmp0_0 = A_01 + A_10;
                        const double tmp0_1 = A_11*w3;
                        const double tmp1_1 = A_00*w0;
                        const double tmp2_1 = A_01*w1;
                        const double tmp3_1 = A_10*w2;
                        const double tmp4_1 = tmp0_0*w1;
                        const double tmp5_1 = A_11*w4;
                        const double tmp6_1 = A_00*w5;
                        const double tmp7_1 = tmp0_0*w2;
                        const double tmp8_1 = A_10*w1;
                        const double tmp9_1 = A_01*w2;
                        EM_S[INDEX2(0,0,4)]+=tmp0_1 + tmp4_1 + tmp6_1;
                        EM_S[INDEX2(1,0,4)]+=tmp0_1 + tmp1_1 + tmp8_1 + tmp9_1;
                        EM_S[INDEX2(2,0,4)]+=tmp2_1 + tmp3_1 + tmp5_1 + tmp6_1;
                        EM_S[INDEX2(3,0,4)]+=tmp1_1 + tmp5_1 + tmp7_1;
                        EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                        EM_S[INDEX2(1,1,4)]+=tmp0_1 + tmp6_1 + tmp7_1;
                        EM_S[INDEX2(2,1,4)]+=tmp1_1 + tmp4_1 + tmp5_1;
                        EM_S[INDEX2(3,1,4)]+=tmp5_1 + tmp6_1 + tmp8_1 + tmp9_1;
                        EM_S[INDEX2(0,2,4)]+=tmp5_1 + tmp6_1 + tmp8_1 + tmp9_1;
                        EM_S[INDEX2(1,2,4)]+=tmp1_1 + tmp4_1 + tmp5_1;
                        EM_S[INDEX2(2,2,4)]+=tmp0_1 + tmp6_1 + tmp7_1;
                        EM_S[INDEX2(3,2,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                        EM_S[INDEX2(0,3,4)]+=tmp1_1 + tmp5_1 + tmp7_1;
                        EM_S[INDEX2(1,3,4)]+=tmp2_1 + tmp3_1 + tmp5_1 + tmp6_1;
                        EM_S[INDEX2(2,3,4)]+=tmp0_1 + tmp1_1 + tmp8_1 + tmp9_1;
                        EM_S[INDEX2(3,3,4)]+=tmp0_1 + tmp4_1 + tmp6_1;
                    }
                    ///////////////
                    // process B //
                    ///////////////
                    if (!B.isEmpty()) {
                        add_EM_S=true;
                        const double* B_p=const_cast<escript::Data*>(&B)->getSampleDataRO(e);
                        const double tmp2_1 = B_p[0]*w8;
                        const double tmp0_1 = B_p[0]*w6;
                        const double tmp3_1 = B_p[1]*w9;
                        const double tmp1_1 = B_p[1]*w7;
                        EM_S[INDEX2(0,0,4)]+=tmp0_1 + tmp1_1;
                        EM_S[INDEX2(1,0,4)]+=tmp1_1 + tmp2_1;
                        EM_S[INDEX2(2,0,4)]+=tmp0_1 + tmp3_1;
                        EM_S[INDEX2(3,0,4)]+=tmp2_1 + tmp3_1;
                        EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1;
                        EM_S[INDEX2(1,1,4)]+=tmp1_1 + tmp2_1;
                        EM_S[INDEX2(2,1,4)]+=tmp0_1 + tmp3_1;
                        EM_S[INDEX2(3,1,4)]+=tmp2_1 + tmp3_1;
                        EM_S[INDEX2(0,2,4)]+=tmp0_1 + tmp1_1;
                        EM_S[INDEX2(1,2,4)]+=tmp1_1 + tmp2_1;
                        EM_S[INDEX2(2,2,4)]+=tmp0_1 + tmp3_1;
                        EM_S[INDEX2(3,2,4)]+=tmp2_1 + tmp3_1;
                        EM_S[INDEX2(0,3,4)]+=tmp0_1 + tmp1_1;
                        EM_S[INDEX2(1,3,4)]+=tmp1_1 + tmp2_1;
                        EM_S[INDEX2(2,3,4)]+=tmp0_1 + tmp3_1;
                        EM_S[INDEX2(3,3,4)]+=tmp2_1 + tmp3_1;
                    }
                    ///////////////
                    // process C //
                    ///////////////
                    if (!C.isEmpty()) {
                        add_EM_S=true;
                        const double* C_p=const_cast<escript::Data*>(&C)->getSampleDataRO(e);
                        const double tmp1_1 = C_p[1]*w7;
                        const double tmp0_1 = C_p[0]*w8;
                        const double tmp3_1 = C_p[0]*w6;
                        const double tmp2_1 = C_p[1]*w9;
                        EM_S[INDEX2(0,0,4)]+=tmp1_1 + tmp3_1;
                        EM_S[INDEX2(1,0,4)]+=tmp1_1 + tmp3_1;
                        EM_S[INDEX2(2,0,4)]+=tmp1_1 + tmp3_1;
                        EM_S[INDEX2(3,0,4)]+=tmp1_1 + tmp3_1;
                        EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1;
                        EM_S[INDEX2(1,1,4)]+=tmp0_1 + tmp1_1;
                        EM_S[INDEX2(2,1,4)]+=tmp0_1 + tmp1_1;
                        EM_S[INDEX2(3,1,4)]+=tmp0_1 + tmp1_1;
                        EM_S[INDEX2(0,2,4)]+=tmp2_1 + tmp3_1;
                        EM_S[INDEX2(1,2,4)]+=tmp2_1 + tmp3_1;
                        EM_S[INDEX2(2,2,4)]+=tmp2_1 + tmp3_1;
                        EM_S[INDEX2(3,2,4)]+=tmp2_1 + tmp3_1;
                        EM_S[INDEX2(0,3,4)]+=tmp0_1 + tmp2_1;
                        EM_S[INDEX2(1,3,4)]+=tmp0_1 + tmp2_1;
                        EM_S[INDEX2(2,3,4)]+=tmp0_1 + tmp2_1;
                        EM_S[INDEX2(3,3,4)]+=tmp0_1 + tmp2_1;
                    }
                    ///////////////
                    // process D //
                    ///////////////
                    if (!D.isEmpty()) {
                        add_EM_S=true;
                        const double* D_p=const_cast<escript::Data*>(&D)->getSampleDataRO(e);
                        const double tmp0_1 = D_p[0]*w10;
                        EM_S[INDEX2(0,0,4)]+=tmp0_1;
                        EM_S[INDEX2(1,0,4)]+=tmp0_1;
                        EM_S[INDEX2(2,0,4)]+=tmp0_1;
                        EM_S[INDEX2(3,0,4)]+=tmp0_1;
                        EM_S[INDEX2(0,1,4)]+=tmp0_1;
                        EM_S[INDEX2(1,1,4)]+=tmp0_1;
                        EM_S[INDEX2(2,1,4)]+=tmp0_1;
                        EM_S[INDEX2(3,1,4)]+=tmp0_1;
                        EM_S[INDEX2(0,2,4)]+=tmp0_1;
                        EM_S[INDEX2(1,2,4)]+=tmp0_1;
                        EM_S[INDEX2(2,2,4)]+=tmp0_1;
                        EM_S[INDEX2(3,2,4)]+=tmp0_1;
                        EM_S[INDEX2(0,3,4)]+=tmp0_1;
                        EM_S[INDEX2(1,3,4)]+=tmp0_1;
                        EM_S[INDEX2(2,3,4)]+=tmp0_1;
                        EM_S[INDEX2(3,3,4)]+=tmp0_1;
                    }
                    ///////////////
                    // process X //
                    ///////////////
                    if (!X.isEmpty()) {
                        add_EM_F=true;
                        const double* X_p=const_cast<escript::Data*>(&X)->getSampleDataRO(e);
                        const double tmp0_1 = X_p[0]*w11;
                        const double tmp2_1 = X_p[0]*w13;
                        const double tmp1_1 = X_p[1]*w12;
                        const double tmp3_1 = X_p[1]*w14;
                        EM_F[0]+=tmp0_1 + tmp1_1;
                        EM_F[1]+=tmp1_1 + tmp2_1;
                        EM_F[2]+=tmp0_1 + tmp3_1;
                        EM_F[3]+=tmp2_1 + tmp3_1;
                    }
                    ///////////////
                    // process Y //
                    ///////////////
                    if (!Y.isEmpty()) {
                        add_EM_F=true;
                        const double* Y_p=const_cast<escript::Data*>(&Y)->getSampleDataRO(e);
                        const double tmp0_1 = Y_p[0]*w15;
                        EM_F[0]+=tmp0_1;
                        EM_F[1]+=tmp0_1;
                        EM_F[2]+=tmp0_1;
                        EM_F[3]+=tmp0_1;
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

    const double w0 = -0.1555021169820365539*m_dx[1]/m_dx[0];
    const double w1 = 0.041666666666666666667;
    const double w2 = -0.15550211698203655390;
    const double w3 = 0.041666666666666666667*m_dx[0]/m_dx[1];
    const double w4 = 0.15550211698203655390;
    const double w5 = -0.041666666666666666667;
    const double w6 = -0.01116454968463011277*m_dx[1]/m_dx[0];
    const double w7 = 0.011164549684630112770;
    const double w8 = -0.011164549684630112770;
    const double w9 = -0.041666666666666666667*m_dx[1]/m_dx[0];
    const double w10 = -0.041666666666666666667*m_dx[0]/m_dx[1];
    const double w11 = 0.1555021169820365539*m_dx[1]/m_dx[0];
    const double w12 = 0.1555021169820365539*m_dx[0]/m_dx[1];
    const double w13 = 0.01116454968463011277*m_dx[0]/m_dx[1];
    const double w14 = 0.01116454968463011277*m_dx[1]/m_dx[0];
    const double w15 = 0.041666666666666666667*m_dx[1]/m_dx[0];
    const double w16 = -0.01116454968463011277*m_dx[0]/m_dx[1];
    const double w17 = -0.1555021169820365539*m_dx[0]/m_dx[1];
    const double w18 = -0.33333333333333333333*m_dx[1]/m_dx[0];
    const double w19 = 0.25000000000000000000;
    const double w20 = -0.25000000000000000000;
    const double w21 = 0.16666666666666666667*m_dx[0]/m_dx[1];
    const double w22 = -0.16666666666666666667*m_dx[1]/m_dx[0];
    const double w23 = -0.16666666666666666667*m_dx[0]/m_dx[1];
    const double w24 = 0.33333333333333333333*m_dx[1]/m_dx[0];
    const double w25 = 0.33333333333333333333*m_dx[0]/m_dx[1];
    const double w26 = 0.16666666666666666667*m_dx[1]/m_dx[0];
    const double w27 = -0.33333333333333333333*m_dx[0]/m_dx[1];
    const double w28 = -0.032861463941450536761*m_dx[1];
    const double w29 = -0.032861463941450536761*m_dx[0];
    const double w30 = -0.12264065304058601714*m_dx[1];
    const double w31 = -0.0023593469594139828636*m_dx[1];
    const double w32 = -0.008805202725216129906*m_dx[0];
    const double w33 = -0.008805202725216129906*m_dx[1];
    const double w34 = 0.032861463941450536761*m_dx[1];
    const double w35 = 0.008805202725216129906*m_dx[1];
    const double w36 = 0.008805202725216129906*m_dx[0];
    const double w37 = 0.0023593469594139828636*m_dx[1];
    const double w38 = 0.12264065304058601714*m_dx[1];
    const double w39 = 0.032861463941450536761*m_dx[0];
    const double w40 = -0.12264065304058601714*m_dx[0];
    const double w41 = -0.0023593469594139828636*m_dx[0];
    const double w42 = 0.0023593469594139828636*m_dx[0];
    const double w43 = 0.12264065304058601714*m_dx[0];
    const double w44 = -0.16666666666666666667*m_dx[1];
    const double w45 = -0.083333333333333333333*m_dx[0];
    const double w46 = 0.083333333333333333333*m_dx[1];
    const double w47 = 0.16666666666666666667*m_dx[1];
    const double w48 = 0.083333333333333333333*m_dx[0];
    const double w49 = -0.16666666666666666667*m_dx[0];
    const double w50 = 0.16666666666666666667*m_dx[0];
    const double w51 = -0.083333333333333333333*m_dx[1];
    const double w52 = 0.025917019497006092316*m_dx[0]*m_dx[1];
    const double w53 = 0.0018607582807716854616*m_dx[0]*m_dx[1];
    const double w54 = 0.0069444444444444444444*m_dx[0]*m_dx[1];
    const double w55 = 0.09672363354357992482*m_dx[0]*m_dx[1];
    const double w56 = 0.00049858867864229740201*m_dx[0]*m_dx[1];
    const double w57 = 0.055555555555555555556*m_dx[0]*m_dx[1];
    const double w58 = 0.027777777777777777778*m_dx[0]*m_dx[1];
    const double w59 = 0.11111111111111111111*m_dx[0]*m_dx[1];
    const double w60 = -0.19716878364870322056*m_dx[1];
    const double w61 = -0.19716878364870322056*m_dx[0];
    const double w62 = -0.052831216351296779436*m_dx[0];
    const double w63 = -0.052831216351296779436*m_dx[1];
    const double w64 = 0.19716878364870322056*m_dx[1];
    const double w65 = 0.052831216351296779436*m_dx[1];
    const double w66 = 0.19716878364870322056*m_dx[0];
    const double w67 = 0.052831216351296779436*m_dx[0];
    const double w68 = -0.5*m_dx[1];
    const double w69 = -0.5*m_dx[0];
    const double w70 = 0.5*m_dx[1];
    const double w71 = 0.5*m_dx[0];
    const double w72 = 0.1555021169820365539*m_dx[0]*m_dx[1];
    const double w73 = 0.041666666666666666667*m_dx[0]*m_dx[1];
    const double w74 = 0.01116454968463011277*m_dx[0]*m_dx[1];
    const double w75 = 0.25*m_dx[0]*m_dx[1];

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
                        if (A.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double A_00_0 = A_p[INDEX5(k,0,m,0,0, numEq,2,numComp,2)];
                                    const double A_01_0 = A_p[INDEX5(k,0,m,1,0, numEq,2,numComp,2)];
                                    const double A_10_0 = A_p[INDEX5(k,1,m,0,0, numEq,2,numComp,2)];
                                    const double A_11_0 = A_p[INDEX5(k,1,m,1,0, numEq,2,numComp,2)];
                                    const double A_00_1 = A_p[INDEX5(k,0,m,0,1, numEq,2,numComp,2)];
                                    const double A_01_1 = A_p[INDEX5(k,0,m,1,1, numEq,2,numComp,2)];
                                    const double A_10_1 = A_p[INDEX5(k,1,m,0,1, numEq,2,numComp,2)];
                                    const double A_11_1 = A_p[INDEX5(k,1,m,1,1, numEq,2,numComp,2)];
                                    const double A_00_2 = A_p[INDEX5(k,0,m,0,2, numEq,2,numComp,2)];
                                    const double A_01_2 = A_p[INDEX5(k,0,m,1,2, numEq,2,numComp,2)];
                                    const double A_10_2 = A_p[INDEX5(k,1,m,0,2, numEq,2,numComp,2)];
                                    const double A_11_2 = A_p[INDEX5(k,1,m,1,2, numEq,2,numComp,2)];
                                    const double A_00_3 = A_p[INDEX5(k,0,m,0,3, numEq,2,numComp,2)];
                                    const double A_01_3 = A_p[INDEX5(k,0,m,1,3, numEq,2,numComp,2)];
                                    const double A_10_3 = A_p[INDEX5(k,1,m,0,3, numEq,2,numComp,2)];
                                    const double A_11_3 = A_p[INDEX5(k,1,m,1,3, numEq,2,numComp,2)];
                                    const double tmp0_0 = A_01_0 + A_01_3;
                                    const double tmp1_0 = A_00_0 + A_00_1;
                                    const double tmp2_0 = A_11_0 + A_11_1 + A_11_2 + A_11_3;
                                    const double tmp3_0 = A_00_2 + A_00_3;
                                    const double tmp4_0 = A_10_1 + A_10_2;
                                    const double tmp5_0 = A_00_0 + A_00_1 + A_00_2 + A_00_3;
                                    const double tmp6_0 = A_01_3 + A_10_0;
                                    const double tmp7_0 = A_01_0 + A_10_3;
                                    const double tmp8_0 = A_01_1 + A_01_2 + A_10_1 + A_10_2;
                                    const double tmp9_0 = A_01_0 + A_10_0;
                                    const double tmp10_0 = A_01_3 + A_10_3;
                                    const double tmp11_0 = A_11_1 + A_11_3;
                                    const double tmp12_0 = A_11_0 + A_11_2;
                                    const double tmp13_0 = A_01_2 + A_10_1;
                                    const double tmp14_0 = A_01_0 + A_01_3 + A_10_0 + A_10_3;
                                    const double tmp15_0 = A_01_1 + A_10_2;
                                    const double tmp16_0 = A_10_0 + A_10_3;
                                    const double tmp17_0 = A_01_1 + A_01_2;
                                    const double tmp18_0 = A_01_1 + A_10_1;
                                    const double tmp19_0 = A_01_2 + A_10_2;
                                    const double tmp0_1 = A_10_3*w8;
                                    const double tmp1_1 = tmp0_0*w1;
                                    const double tmp2_1 = A_01_1*w4;
                                    const double tmp3_1 = tmp1_0*w0;
                                    const double tmp4_1 = A_01_2*w7;
                                    const double tmp5_1 = tmp2_0*w3;
                                    const double tmp6_1 = tmp3_0*w6;
                                    const double tmp7_1 = A_10_0*w2;
                                    const double tmp8_1 = tmp4_0*w5;
                                    const double tmp9_1 = tmp2_0*w10;
                                    const double tmp10_1 = tmp5_0*w9;
                                    const double tmp11_1 = tmp6_0*w7;
                                    const double tmp12_1 = tmp7_0*w4;
                                    const double tmp13_1 = tmp8_0*w1;
                                    const double tmp14_1 = A_10_0*w8;
                                    const double tmp15_1 = A_01_2*w4;
                                    const double tmp16_1 = tmp3_0*w0;
                                    const double tmp17_1 = A_01_1*w7;
                                    const double tmp18_1 = tmp1_0*w6;
                                    const double tmp19_1 = A_10_3*w2;
                                    const double tmp20_1 = tmp9_0*w4;
                                    const double tmp21_1 = tmp1_0*w11;
                                    const double tmp22_1 = tmp10_0*w7;
                                    const double tmp23_1 = tmp3_0*w14;
                                    const double tmp24_1 = tmp11_0*w13;
                                    const double tmp25_1 = tmp12_0*w12;
                                    const double tmp26_1 = tmp10_0*w4;
                                    const double tmp27_1 = tmp3_0*w11;
                                    const double tmp28_1 = tmp9_0*w7;
                                    const double tmp29_1 = tmp1_0*w14;
                                    const double tmp30_1 = tmp12_0*w13;
                                    const double tmp31_1 = tmp11_0*w12;
                                    const double tmp32_1 = tmp13_0*w2;
                                    const double tmp33_1 = tmp14_0*w5;
                                    const double tmp34_1 = tmp15_0*w8;
                                    const double tmp35_1 = A_01_0*w8;
                                    const double tmp36_1 = tmp16_0*w1;
                                    const double tmp37_1 = A_10_1*w4;
                                    const double tmp38_1 = tmp5_0*w15;
                                    const double tmp39_1 = A_10_2*w7;
                                    const double tmp40_1 = tmp11_0*w17;
                                    const double tmp41_1 = A_01_3*w2;
                                    const double tmp42_1 = tmp12_0*w16;
                                    const double tmp43_1 = tmp17_0*w5;
                                    const double tmp44_1 = tmp7_0*w7;
                                    const double tmp45_1 = tmp6_0*w4;
                                    const double tmp46_1 = A_01_3*w8;
                                    const double tmp47_1 = A_10_2*w4;
                                    const double tmp48_1 = A_10_1*w7;
                                    const double tmp49_1 = tmp12_0*w17;
                                    const double tmp50_1 = A_01_0*w2;
                                    const double tmp51_1 = tmp11_0*w16;
                                    const double tmp52_1 = tmp18_0*w8;
                                    const double tmp53_1 = tmp19_0*w2;
                                    const double tmp54_1 = tmp13_0*w8;
                                    const double tmp55_1 = tmp15_0*w2;
                                    const double tmp56_1 = tmp19_0*w8;
                                    const double tmp57_1 = tmp18_0*w2;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp13_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp24_1 + tmp25_1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp36_1 + tmp37_1 + tmp39_1 + tmp3_1 + tmp43_1 + tmp46_1 + tmp50_1 + tmp5_1 + tmp6_1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp0_1 + tmp15_1 + tmp17_1 + tmp1_1 + tmp38_1 + tmp49_1 + tmp51_1 + tmp7_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp10_1 + tmp32_1 + tmp33_1 + tmp34_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp21_1 + tmp23_1 + tmp30_1 + tmp31_1 + tmp33_1 + tmp56_1 + tmp57_1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp10_1 + tmp13_1 + tmp44_1 + tmp45_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp43_1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp36_1 + tmp38_1 + tmp43_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp24_1 + tmp25_1 + tmp27_1 + tmp29_1 + tmp33_1 + tmp52_1 + tmp53_1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp14_1 + tmp15_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp1_1 + tmp5_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp10_1 + tmp33_1 + tmp54_1 + tmp55_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp14_1 + tmp19_1 + tmp1_1 + tmp2_1 + tmp38_1 + tmp40_1 + tmp42_1 + tmp4_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp16_1 + tmp18_1 + tmp35_1 + tmp36_1 + tmp41_1 + tmp43_1 + tmp47_1 + tmp48_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp13_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1;
                                }
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double A_00 = A_p[INDEX4(k,0,m,0, numEq,2, numComp)];
                                    const double A_01 = A_p[INDEX4(k,0,m,1, numEq,2, numComp)];
                                    const double A_10 = A_p[INDEX4(k,1,m,0, numEq,2, numComp)];
                                    const double A_11 = A_p[INDEX4(k,1,m,1, numEq,2, numComp)];
                                    const double tmp0_0 = A_01 + A_10;
                                    const double tmp0_1 = A_00*w18;
                                    const double tmp1_1 = A_01*w19;
                                    const double tmp2_1 = A_10*w20;
                                    const double tmp3_1 = A_11*w21;
                                    const double tmp4_1 = A_00*w22;
                                    const double tmp5_1 = tmp0_0*w19;
                                    const double tmp6_1 = A_11*w23;
                                    const double tmp7_1 = A_11*w25;
                                    const double tmp8_1 = A_00*w24;
                                    const double tmp9_1 = tmp0_0*w20;
                                    const double tmp10_1 = A_01*w20;
                                    const double tmp11_1 = A_11*w27;
                                    const double tmp12_1 = A_00*w26;
                                    const double tmp13_1 = A_10*w19;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp5_1 + tmp7_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp0_1 + tmp10_1 + tmp13_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp11_1 + tmp12_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp4_1 + tmp6_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp7_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp4_1 + tmp5_1 + tmp6_1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp4_1 + tmp5_1 + tmp6_1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp7_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp4_1 + tmp6_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp11_1 + tmp12_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp0_1 + tmp10_1 + tmp13_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp5_1 + tmp7_1 + tmp8_1;
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
                                    const double tmp0_0 = B_1_0 + B_1_1;
                                    const double tmp1_0 = B_1_2 + B_1_3;
                                    const double tmp2_0 = B_0_1 + B_0_3;
                                    const double tmp3_0 = B_0_0 + B_0_2;
                                    const double tmp63_1 = B_1_1*w42;
                                    const double tmp79_1 = B_1_1*w40;
                                    const double tmp37_1 = tmp3_0*w35;
                                    const double tmp8_1 = tmp0_0*w32;
                                    const double tmp71_1 = B_0_1*w34;
                                    const double tmp19_1 = B_0_3*w31;
                                    const double tmp15_1 = B_0_3*w34;
                                    const double tmp9_1 = tmp3_0*w34;
                                    const double tmp35_1 = B_1_0*w36;
                                    const double tmp66_1 = B_0_3*w28;
                                    const double tmp28_1 = B_1_0*w42;
                                    const double tmp22_1 = B_1_0*w40;
                                    const double tmp16_1 = B_1_2*w29;
                                    const double tmp6_1 = tmp2_0*w35;
                                    const double tmp55_1 = B_1_3*w40;
                                    const double tmp50_1 = B_1_3*w42;
                                    const double tmp7_1 = tmp1_0*w29;
                                    const double tmp1_1 = tmp1_0*w32;
                                    const double tmp57_1 = B_0_3*w30;
                                    const double tmp18_1 = B_1_1*w32;
                                    const double tmp53_1 = B_1_0*w41;
                                    const double tmp61_1 = B_1_3*w36;
                                    const double tmp27_1 = B_0_3*w38;
                                    const double tmp64_1 = B_0_2*w30;
                                    const double tmp76_1 = B_0_1*w38;
                                    const double tmp39_1 = tmp2_0*w34;
                                    const double tmp62_1 = B_0_1*w31;
                                    const double tmp56_1 = B_0_0*w31;
                                    const double tmp49_1 = B_1_1*w36;
                                    const double tmp2_1 = B_0_2*w31;
                                    const double tmp23_1 = B_0_2*w33;
                                    const double tmp38_1 = B_1_1*w43;
                                    const double tmp74_1 = B_1_2*w41;
                                    const double tmp43_1 = B_1_1*w41;
                                    const double tmp58_1 = B_0_2*w28;
                                    const double tmp67_1 = B_0_0*w33;
                                    const double tmp33_1 = tmp0_0*w39;
                                    const double tmp4_1 = B_0_0*w28;
                                    const double tmp20_1 = B_0_0*w30;
                                    const double tmp13_1 = B_0_2*w38;
                                    const double tmp65_1 = B_1_2*w43;
                                    const double tmp0_1 = tmp0_0*w29;
                                    const double tmp41_1 = tmp3_0*w33;
                                    const double tmp73_1 = B_0_2*w37;
                                    const double tmp69_1 = B_0_0*w38;
                                    const double tmp48_1 = B_1_2*w39;
                                    const double tmp59_1 = B_0_1*w33;
                                    const double tmp17_1 = B_1_3*w41;
                                    const double tmp5_1 = B_0_3*w33;
                                    const double tmp3_1 = B_0_1*w30;
                                    const double tmp21_1 = B_0_1*w28;
                                    const double tmp42_1 = B_1_0*w29;
                                    const double tmp54_1 = B_1_2*w32;
                                    const double tmp60_1 = B_1_0*w39;
                                    const double tmp32_1 = tmp1_0*w36;
                                    const double tmp10_1 = B_0_1*w37;
                                    const double tmp14_1 = B_0_0*w35;
                                    const double tmp29_1 = B_0_1*w35;
                                    const double tmp26_1 = B_1_2*w36;
                                    const double tmp30_1 = B_1_3*w43;
                                    const double tmp70_1 = B_0_2*w35;
                                    const double tmp34_1 = B_1_3*w39;
                                    const double tmp51_1 = B_1_0*w43;
                                    const double tmp31_1 = B_0_2*w34;
                                    const double tmp45_1 = tmp3_0*w28;
                                    const double tmp11_1 = tmp1_0*w39;
                                    const double tmp52_1 = B_1_1*w29;
                                    const double tmp44_1 = B_1_3*w32;
                                    const double tmp25_1 = B_1_1*w39;
                                    const double tmp47_1 = tmp2_0*w33;
                                    const double tmp72_1 = B_1_3*w29;
                                    const double tmp40_1 = tmp2_0*w28;
                                    const double tmp46_1 = B_1_2*w40;
                                    const double tmp36_1 = B_1_2*w42;
                                    const double tmp24_1 = B_0_0*w37;
                                    const double tmp77_1 = B_0_3*w35;
                                    const double tmp68_1 = B_0_3*w37;
                                    const double tmp78_1 = B_0_0*w34;
                                    const double tmp12_1 = tmp0_0*w36;
                                    const double tmp75_1 = B_1_0*w32;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp0_1 + tmp1_1 + tmp68_1 + tmp69_1 + tmp70_1 + tmp71_1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp45_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp32_1 + tmp33_1 + tmp6_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp72_1 + tmp73_1 + tmp74_1 + tmp75_1 + tmp76_1 + tmp77_1 + tmp78_1 + tmp79_1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp32_1 + tmp33_1 + tmp40_1 + tmp41_1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp60_1 + tmp61_1 + tmp62_1 + tmp63_1 + tmp64_1 + tmp65_1 + tmp66_1 + tmp67_1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp40_1 + tmp41_1 + tmp7_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp37_1 + tmp39_1 + tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp11_1 + tmp12_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1;
                                }
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double B_0 = B_p[INDEX3(k,0,m, numEq, 2)];
                                    const double B_1 = B_p[INDEX3(k,1,m, numEq, 2)];
                                    const double tmp6_1 = B_1*w50;
                                    const double tmp1_1 = B_1*w45;
                                    const double tmp5_1 = B_1*w49;
                                    const double tmp4_1 = B_1*w48;
                                    const double tmp0_1 = B_0*w44;
                                    const double tmp2_1 = B_0*w46;
                                    const double tmp7_1 = B_0*w51;
                                    const double tmp3_1 = B_0*w47;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp0_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp1_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp6_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp2_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0_1 + tmp1_1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp4_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp2_1 + tmp6_1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp5_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp0_1 + tmp6_1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp1_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp2_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp0_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp3_1 + tmp6_1;
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
                                    const double tmp2_0 = C_0_1 + C_0_3;
                                    const double tmp1_0 = C_1_2 + C_1_3;
                                    const double tmp3_0 = C_0_0 + C_0_2;
                                    const double tmp0_0 = C_1_0 + C_1_1;
                                    const double tmp64_1 = C_0_2*w30;
                                    const double tmp14_1 = C_0_2*w28;
                                    const double tmp19_1 = C_0_3*w31;
                                    const double tmp22_1 = C_1_0*w40;
                                    const double tmp37_1 = tmp3_0*w35;
                                    const double tmp29_1 = C_0_1*w35;
                                    const double tmp73_1 = C_0_2*w37;
                                    const double tmp74_1 = C_1_2*w41;
                                    const double tmp52_1 = C_1_3*w39;
                                    const double tmp25_1 = C_1_1*w39;
                                    const double tmp62_1 = C_0_1*w31;
                                    const double tmp79_1 = C_1_1*w40;
                                    const double tmp43_1 = C_1_1*w36;
                                    const double tmp27_1 = C_0_3*w38;
                                    const double tmp28_1 = C_1_0*w42;
                                    const double tmp63_1 = C_1_1*w42;
                                    const double tmp59_1 = C_0_3*w34;
                                    const double tmp72_1 = C_1_3*w29;
                                    const double tmp40_1 = tmp2_0*w35;
                                    const double tmp13_1 = C_0_3*w30;
                                    const double tmp51_1 = C_1_2*w40;
                                    const double tmp54_1 = C_1_2*w42;
                                    const double tmp12_1 = C_0_0*w31;
                                    const double tmp2_1 = tmp1_0*w32;
                                    const double tmp68_1 = C_0_2*w31;
                                    const double tmp75_1 = C_1_0*w32;
                                    const double tmp49_1 = C_1_1*w41;
                                    const double tmp4_1 = C_0_2*w35;
                                    const double tmp66_1 = C_0_3*w28;
                                    const double tmp56_1 = C_0_1*w37;
                                    const double tmp5_1 = C_0_1*w34;
                                    const double tmp38_1 = tmp2_0*w34;
                                    const double tmp76_1 = C_0_1*w38;
                                    const double tmp21_1 = C_0_1*w28;
                                    const double tmp69_1 = C_0_1*w30;
                                    const double tmp53_1 = C_1_0*w36;
                                    const double tmp42_1 = C_1_2*w39;
                                    const double tmp32_1 = tmp1_0*w29;
                                    const double tmp45_1 = C_1_0*w43;
                                    const double tmp33_1 = tmp0_0*w32;
                                    const double tmp35_1 = C_1_0*w41;
                                    const double tmp26_1 = C_1_2*w36;
                                    const double tmp67_1 = C_0_0*w33;
                                    const double tmp31_1 = C_0_2*w34;
                                    const double tmp20_1 = C_0_0*w30;
                                    const double tmp70_1 = C_0_0*w28;
                                    const double tmp8_1 = tmp0_0*w39;
                                    const double tmp30_1 = C_1_3*w43;
                                    const double tmp0_1 = tmp0_0*w29;
                                    const double tmp17_1 = C_1_3*w41;
                                    const double tmp58_1 = C_0_0*w35;
                                    const double tmp9_1 = tmp3_0*w33;
                                    const double tmp61_1 = C_1_3*w36;
                                    const double tmp41_1 = tmp3_0*w34;
                                    const double tmp50_1 = C_1_3*w32;
                                    const double tmp18_1 = C_1_1*w32;
                                    const double tmp6_1 = tmp1_0*w36;
                                    const double tmp3_1 = C_0_0*w38;
                                    const double tmp34_1 = C_1_1*w29;
                                    const double tmp77_1 = C_0_3*w35;
                                    const double tmp65_1 = C_1_2*w43;
                                    const double tmp71_1 = C_0_3*w33;
                                    const double tmp55_1 = C_1_1*w43;
                                    const double tmp46_1 = tmp3_0*w28;
                                    const double tmp24_1 = C_0_0*w37;
                                    const double tmp10_1 = tmp1_0*w39;
                                    const double tmp48_1 = C_1_0*w29;
                                    const double tmp15_1 = C_0_1*w33;
                                    const double tmp36_1 = C_1_2*w32;
                                    const double tmp60_1 = C_1_0*w39;
                                    const double tmp47_1 = tmp2_0*w33;
                                    const double tmp16_1 = C_1_2*w29;
                                    const double tmp1_1 = C_0_3*w37;
                                    const double tmp7_1 = tmp2_0*w28;
                                    const double tmp39_1 = C_1_3*w40;
                                    const double tmp44_1 = C_1_3*w42;
                                    const double tmp57_1 = C_0_2*w38;
                                    const double tmp78_1 = C_0_0*w34;
                                    const double tmp11_1 = tmp0_0*w36;
                                    const double tmp23_1 = C_0_2*w33;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp32_1 + tmp33_1 + tmp7_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp32_1 + tmp33_1 + tmp40_1 + tmp41_1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp46_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp37_1 + tmp38_1 + tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp10_1 + tmp11_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp60_1 + tmp61_1 + tmp62_1 + tmp63_1 + tmp64_1 + tmp65_1 + tmp66_1 + tmp67_1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp0_1 + tmp2_1 + tmp68_1 + tmp69_1 + tmp70_1 + tmp71_1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp40_1 + tmp41_1 + tmp6_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp72_1 + tmp73_1 + tmp74_1 + tmp75_1 + tmp76_1 + tmp77_1 + tmp78_1 + tmp79_1;
                                }
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double C_0 = C_p[INDEX3(k, m, 0, numEq, numComp)];
                                    const double C_1 = C_p[INDEX3(k, m, 1, numEq, numComp)];
                                    const double tmp1_1 = C_1*w45;
                                    const double tmp3_1 = C_0*w51;
                                    const double tmp4_1 = C_0*w44;
                                    const double tmp7_1 = C_0*w46;
                                    const double tmp5_1 = C_1*w49;
                                    const double tmp2_1 = C_1*w48;
                                    const double tmp0_1 = C_0*w47;
                                    const double tmp6_1 = C_1*w50;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0_1 + tmp1_1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp2_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp0_1 + tmp6_1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp1_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp5_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp1_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp3_1 + tmp6_1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp6_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp0_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp4_1 + tmp6_1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp1_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp2_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp0_1 + tmp5_1;
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
                                    const double D_0 = D_p[INDEX3(k, m, 0, numEq, numComp)];
                                    const double D_1 = D_p[INDEX3(k, m, 1, numEq, numComp)];
                                    const double D_2 = D_p[INDEX3(k, m, 2, numEq, numComp)];
                                    const double D_3 = D_p[INDEX3(k, m, 3, numEq, numComp)];
                                    const double tmp4_0 = D_1 + D_3;
                                    const double tmp2_0 = D_0 + D_1 + D_2 + D_3;
                                    const double tmp5_0 = D_0 + D_2;
                                    const double tmp0_0 = D_0 + D_1;
                                    const double tmp6_0 = D_0 + D_3;
                                    const double tmp1_0 = D_2 + D_3;
                                    const double tmp3_0 = D_1 + D_2;
                                    const double tmp16_1 = D_1*w56;
                                    const double tmp14_1 = tmp6_0*w54;
                                    const double tmp8_1 = D_3*w55;
                                    const double tmp2_1 = tmp2_0*w54;
                                    const double tmp12_1 = tmp5_0*w52;
                                    const double tmp4_1 = tmp0_0*w53;
                                    const double tmp3_1 = tmp1_0*w52;
                                    const double tmp13_1 = tmp4_0*w53;
                                    const double tmp10_1 = tmp4_0*w52;
                                    const double tmp1_1 = tmp1_0*w53;
                                    const double tmp7_1 = D_3*w56;
                                    const double tmp0_1 = tmp0_0*w52;
                                    const double tmp11_1 = tmp5_0*w53;
                                    const double tmp9_1 = D_0*w56;
                                    const double tmp5_1 = tmp3_0*w54;
                                    const double tmp18_1 = D_2*w56;
                                    const double tmp17_1 = D_1*w55;
                                    const double tmp6_1 = D_0*w55;
                                    const double tmp15_1 = D_2*w55;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0_1 + tmp1_1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp2_1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp5_1 + tmp6_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp5_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp2_1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp10_1 + tmp11_1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp2_1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp12_1 + tmp13_1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp12_1 + tmp13_1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp10_1 + tmp11_1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp14_1 + tmp15_1 + tmp16_1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp0_1 + tmp1_1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp2_1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp14_1 + tmp17_1 + tmp18_1;
                                }
                             }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double D_0 = D_p[INDEX2(k, m, numEq)];
                                    const double tmp2_1 = D_0*w59;
                                    const double tmp1_1 = D_0*w58;
                                    const double tmp0_1 = D_0*w57;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp1_1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp2_1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp2_1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp1_1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp1_1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp2_1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp1_1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp2_1;
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
                                const double X_0_0 = X_p[INDEX3(k, 0, 0, numEq, 2)];
                                const double X_1_0 = X_p[INDEX3(k, 1, 0, numEq, 2)];
                                const double X_0_1 = X_p[INDEX3(k, 0, 1, numEq, 2)];
                                const double X_1_1 = X_p[INDEX3(k, 1, 1, numEq, 2)];
                                const double X_0_2 = X_p[INDEX3(k, 0, 2, numEq, 2)];
                                const double X_1_2 = X_p[INDEX3(k, 1, 2, numEq, 2)];
                                const double X_0_3 = X_p[INDEX3(k, 0, 3, numEq, 2)];
                                const double X_1_3 = X_p[INDEX3(k, 1, 3, numEq, 2)];
                                const double tmp1_0 = X_1_1 + X_1_3;
                                const double tmp3_0 = X_0_0 + X_0_1;
                                const double tmp2_0 = X_1_0 + X_1_2;
                                const double tmp0_0 = X_0_2 + X_0_3;
                                const double tmp8_1 = tmp2_0*w66;
                                const double tmp5_1 = tmp3_0*w64;
                                const double tmp14_1 = tmp0_0*w64;
                                const double tmp3_1 = tmp3_0*w60;
                                const double tmp9_1 = tmp3_0*w63;
                                const double tmp13_1 = tmp3_0*w65;
                                const double tmp12_1 = tmp1_0*w66;
                                const double tmp10_1 = tmp0_0*w60;
                                const double tmp2_1 = tmp2_0*w61;
                                const double tmp6_1 = tmp2_0*w62;
                                const double tmp4_1 = tmp0_0*w65;
                                const double tmp11_1 = tmp1_0*w67;
                                const double tmp1_1 = tmp1_0*w62;
                                const double tmp7_1 = tmp1_0*w61;
                                const double tmp0_1 = tmp0_0*w63;
                                const double tmp15_1 = tmp2_0*w67;
                                EM_F[INDEX2(k,0,numEq)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                                EM_F[INDEX2(k,1,numEq)]+=tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1;
                                EM_F[INDEX2(k,2,numEq)]+=tmp10_1 + tmp11_1 + tmp8_1 + tmp9_1;
                                EM_F[INDEX2(k,3,numEq)]+=tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1;
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                const double X_0 = X_p[INDEX2(k, 0, numEq)];
                                const double X_1 = X_p[INDEX2(k, 1, numEq)];
                                const double tmp0_1 = X_1*w69;
                                const double tmp1_1 = X_0*w68;
                                const double tmp2_1 = X_0*w70;
                                const double tmp3_1 = X_1*w71;
                                EM_F[INDEX2(k,0,numEq)]+=tmp0_1 + tmp1_1;
                                EM_F[INDEX2(k,1,numEq)]+=tmp0_1 + tmp2_1;
                                EM_F[INDEX2(k,2,numEq)]+=tmp1_1 + tmp3_1;
                                EM_F[INDEX2(k,3,numEq)]+=tmp2_1 + tmp3_1;
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
                                const double tmp0_0 = Y_1 + Y_2;
                                const double tmp1_0 = Y_0 + Y_3;
                                const double tmp0_1 = Y_0*w72;
                                const double tmp1_1 = tmp0_0*w73;
                                const double tmp2_1 = Y_3*w74;
                                const double tmp3_1 = Y_1*w72;
                                const double tmp4_1 = tmp1_0*w73;
                                const double tmp5_1 = Y_2*w74;
                                const double tmp6_1 = Y_2*w72;
                                const double tmp7_1 = Y_1*w74;
                                const double tmp8_1 = Y_3*w72;
                                const double tmp9_1 = Y_0*w74;
                                EM_F[INDEX2(k,0,numEq)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                EM_F[INDEX2(k,1,numEq)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                EM_F[INDEX2(k,2,numEq)]+=tmp4_1 + tmp6_1 + tmp7_1;
                                EM_F[INDEX2(k,3,numEq)]+=tmp1_1 + tmp8_1 + tmp9_1;
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                const double tmp0_1 = Y_p[k]*w75;
                                EM_F[INDEX2(k,0,numEq)]+=tmp0_1;
                                EM_F[INDEX2(k,1,numEq)]+=tmp0_1;
                                EM_F[INDEX2(k,2,numEq)]+=tmp0_1;
                                EM_F[INDEX2(k,3,numEq)]+=tmp0_1;
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

    const double w0 = -.25*m_dx[1]/m_dx[0];
    const double w1 = .25;
    const double w2 = -.25;
    const double w3 = .25*m_dx[0]/m_dx[1];
    const double w4 = -.25*m_dx[0]/m_dx[1];
    const double w5 = .25*m_dx[1]/m_dx[0];
    const double w6 = -.125*m_dx[1];
    const double w7 = -.125*m_dx[0];
    const double w8 = .125*m_dx[1];
    const double w9 = .125*m_dx[0];
    const double w10 = .0625*m_dx[0]*m_dx[1];
    const double w11 = -.5*m_dx[1];
    const double w12 = -.5*m_dx[0];
    const double w13 = .5*m_dx[1];
    const double w14 = .5*m_dx[0];
    const double w15 = .25*m_dx[0]*m_dx[1];

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
                                const double A_00 = A_p[INDEX4(k,0,m,0, numEq,2, numComp)];
                                const double A_01 = A_p[INDEX4(k,0,m,1, numEq,2, numComp)];
                                const double A_10 = A_p[INDEX4(k,1,m,0, numEq,2, numComp)];
                                const double A_11 = A_p[INDEX4(k,1,m,1, numEq,2, numComp)];
                                const double tmp0_0 = A_01 + A_10;
                                const double tmp0_1 = A_11*w3;
                                const double tmp1_1 = A_00*w0;
                                const double tmp2_1 = A_01*w1;
                                const double tmp3_1 = A_10*w2;
                                const double tmp4_1 = tmp0_0*w1;
                                const double tmp5_1 = A_11*w4;
                                const double tmp6_1 = A_00*w5;
                                const double tmp7_1 = tmp0_0*w2;
                                const double tmp8_1 = A_10*w1;
                                const double tmp9_1 = A_01*w2;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp0_1 + tmp4_1 + tmp6_1;
                                EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp0_1 + tmp1_1 + tmp8_1 + tmp9_1;
                                EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp2_1 + tmp3_1 + tmp5_1 + tmp6_1;
                                EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp1_1 + tmp5_1 + tmp7_1;
                                EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp0_1 + tmp6_1 + tmp7_1;
                                EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp5_1 + tmp6_1 + tmp8_1 + tmp9_1;
                                EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp5_1 + tmp6_1 + tmp8_1 + tmp9_1;
                                EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp0_1 + tmp6_1 + tmp7_1;
                                EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                                EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp1_1 + tmp5_1 + tmp7_1;
                                EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp2_1 + tmp3_1 + tmp5_1 + tmp6_1;
                                EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp0_1 + tmp1_1 + tmp8_1 + tmp9_1;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp0_1 + tmp4_1 + tmp6_1;
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
                                const double B_0 = B_p[INDEX3(k,0,m, numEq, 2)];
                                const double B_1 = B_p[INDEX3(k,1,m, numEq, 2)];
                                const double tmp0_1 = B_0*w6;
                                const double tmp1_1 = B_1*w7;
                                const double tmp2_1 = B_0*w8;
                                const double tmp3_1 = B_1*w9;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp0_1 + tmp1_1;
                                EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp1_1 + tmp2_1;
                                EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp0_1 + tmp3_1;
                                EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp2_1 + tmp3_1;
                                EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0_1 + tmp1_1;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp1_1 + tmp2_1;
                                EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp0_1 + tmp3_1;
                                EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp2_1 + tmp3_1;
                                EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp0_1 + tmp1_1;
                                EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp1_1 + tmp2_1;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp0_1 + tmp3_1;
                                EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp2_1 + tmp3_1;
                                EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp0_1 + tmp1_1;
                                EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp1_1 + tmp2_1;
                                EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp0_1 + tmp3_1;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp2_1 + tmp3_1;
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
                                const double C_0 = C_p[INDEX3(k, m, 0, numEq, numComp)];
                                const double C_1 = C_p[INDEX3(k, m, 1, numEq, numComp)];
                                const double tmp0_1 = C_0*w8;
                                const double tmp1_1 = C_1*w7;
                                const double tmp2_1 = C_1*w9;
                                const double tmp3_1 = C_0*w6;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp1_1 + tmp3_1;
                                EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp1_1 + tmp3_1;
                                EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp1_1 + tmp3_1;
                                EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp1_1 + tmp3_1;
                                EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0_1 + tmp1_1;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp0_1 + tmp1_1;
                                EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp0_1 + tmp1_1;
                                EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp0_1 + tmp1_1;
                                EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp2_1 + tmp3_1;
                                EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp2_1 + tmp3_1;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp2_1 + tmp3_1;
                                EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp2_1 + tmp3_1;
                                EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp0_1 + tmp2_1;
                                EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp0_1 + tmp2_1;
                                EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp0_1 + tmp2_1;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp0_1 + tmp2_1;
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
                                const double tmp0_1 = D_p[INDEX2(k, m, numEq)]*w10;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp0_1;
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
                            const double X_0 = X_p[INDEX2(k, 0, numEq)];
                            const double X_1 = X_p[INDEX2(k, 1, numEq)];
                            const double tmp0_1 = X_0*w11;
                            const double tmp1_1 = X_1*w12;
                            const double tmp2_1 = X_0*w13;
                            const double tmp3_1 = X_1*w14;
                            EM_F[INDEX2(k,0,numEq)]+=tmp0_1 + tmp1_1;
                            EM_F[INDEX2(k,1,numEq)]+=tmp1_1 + tmp2_1;
                            EM_F[INDEX2(k,2,numEq)]+=tmp0_1 + tmp3_1;
                            EM_F[INDEX2(k,3,numEq)]+=tmp2_1 + tmp3_1;
                        }
                    }
                    ///////////////
                    // process Y //
                    ///////////////
                    if (!Y.isEmpty()) {
                        add_EM_F=true;
                        const double* Y_p=const_cast<escript::Data*>(&Y)->getSampleDataRO(e);
                        for (index_t k=0; k<numEq; k++) {
                            const double tmp0_1 = Y_p[k]*w15;
                            EM_F[INDEX2(k,0,numEq)]+=tmp0_1;
                            EM_F[INDEX2(k,1,numEq)]+=tmp0_1;
                            EM_F[INDEX2(k,2,numEq)]+=tmp0_1;
                            EM_F[INDEX2(k,3,numEq)]+=tmp0_1;
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
    const double w0 = 0.31100423396407310779*m_dx[1];
    const double w1 = 0.022329099369260225539*m_dx[1];
    const double w10 = 0.022329099369260225539*m_dx[0];
    const double w11 = 0.16666666666666666667*m_dx[0];
    const double w12 = 0.33333333333333333333*m_dx[0];
    const double w13 = 0.39433756729740644113*m_dx[0];
    const double w14 = 0.10566243270259355887*m_dx[0];
    const double w15 = 0.5*m_dx[0];
    const double w2 = 0.083333333333333333333*m_dx[1];
    const double w3 = 0.33333333333333333333*m_dx[1];
    const double w4 = 0.16666666666666666667*m_dx[1];
    const double w5 = 0.39433756729740644113*m_dx[1];
    const double w6 = 0.10566243270259355887*m_dx[1];
    const double w7 = 0.5*m_dx[1];
    const double w8 = 0.083333333333333333333*m_dx[0];
    const double w9 = 0.31100423396407310779*m_dx[0];
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
                    const index_t e = k1;
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(e);
                        if (d.actsExpanded()) {
                            const double d_0 = d_p[0];
                            const double d_1 = d_p[1];
                            const double tmp0_0 = d_0 + d_1;
                            const double tmp1_1 = d_1*w1;
                            const double tmp4_1 = d_0*w1;
                            const double tmp0_1 = d_0*w0;
                            const double tmp3_1 = d_1*w0;
                            const double tmp2_1 = tmp0_0*w2;
                            EM_S[INDEX2(0,0,4)]+=tmp0_1 + tmp1_1;
                            EM_S[INDEX2(2,0,4)]+=tmp2_1;
                            EM_S[INDEX2(0,2,4)]+=tmp2_1;
                            EM_S[INDEX2(2,2,4)]+=tmp3_1 + tmp4_1;
                        } else { /* constant data */
                            const double d_0 = d_p[0];
                            const double tmp1_1 = d_0*w4;
                            const double tmp0_1 = d_0*w3;
                            EM_S[INDEX2(0,0,4)]+=tmp0_1;
                            EM_S[INDEX2(2,0,4)]+=tmp1_1;
                            EM_S[INDEX2(0,2,4)]+=tmp1_1;
                            EM_S[INDEX2(2,2,4)]+=tmp0_1;
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        if (y.actsExpanded()) {
                            const double y_0 = y_p[0];
                            const double y_1 = y_p[1];
                            const double tmp3_1 = w5*y_1;
                            const double tmp2_1 = w6*y_0;
                            const double tmp0_1 = w6*y_1;
                            const double tmp1_1 = w5*y_0;
                            EM_F[0]+=tmp0_1 + tmp1_1;
                            EM_F[2]+=tmp2_1 + tmp3_1;
                        } else { /* constant data */
                            const double y_0 = y_p[0];
                            const double tmp0_1 = w7*y_0;
                            EM_F[0]+=tmp0_1;
                            EM_F[2]+=tmp0_1;
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
                            const double tmp0_0 = d_0 + d_1;
                            const double tmp4_1 = d_1*w1;
                            const double tmp1_1 = d_0*w1;
                            const double tmp3_1 = d_0*w0;
                            const double tmp0_1 = d_1*w0;
                            const double tmp2_1 = tmp0_0*w2;
                            EM_S[INDEX2(1,1,4)]+=tmp3_1 + tmp4_1;
                            EM_S[INDEX2(3,1,4)]+=tmp2_1;
                            EM_S[INDEX2(1,3,4)]+=tmp2_1;
                            EM_S[INDEX2(3,3,4)]+=tmp0_1 + tmp1_1;
                        } else { /* constant data */
                            const double d_0 = d_p[0];
                            const double tmp1_1 = d_0*w4;
                            const double tmp0_1 = d_0*w3;
                            EM_S[INDEX2(1,1,4)]+=tmp0_1;
                            EM_S[INDEX2(3,1,4)]+=tmp1_1;
                            EM_S[INDEX2(1,3,4)]+=tmp1_1;
                            EM_S[INDEX2(3,3,4)]+=tmp0_1;
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        if (y.actsExpanded()) {
                            const double y_0 = y_p[0];
                            const double y_1 = y_p[1];
                            const double tmp3_1 = w5*y_1;
                            const double tmp2_1 = w6*y_0;
                            const double tmp0_1 = w6*y_1;
                            const double tmp1_1 = w5*y_0;
                            EM_F[1]+=tmp0_1 + tmp1_1;
                            EM_F[3]+=tmp2_1 + tmp3_1;
                        } else { /* constant data */
                            const double y_0 = y_p[0];
                            const double tmp0_1 = w7*y_0;
                            EM_F[1]+=tmp0_1;
                            EM_F[3]+=tmp0_1;
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
                            const double tmp0_0 = d_0 + d_1;
                            const double tmp4_1 = d_1*w9;
                            const double tmp2_1 = d_0*w9;
                            const double tmp0_1 = tmp0_0*w8;
                            const double tmp1_1 = d_1*w10;
                            const double tmp3_1 = d_0*w10;
                            EM_S[INDEX2(0,0,4)]+=tmp1_1 + tmp2_1;
                            EM_S[INDEX2(1,0,4)]+=tmp0_1;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1;
                            EM_S[INDEX2(1,1,4)]+=tmp3_1 + tmp4_1;
                        } else { /* constant data */
                            const double d_0 = d_p[0];
                            const double tmp0_1 = d_0*w11;
                            const double tmp1_1 = d_0*w12;
                            EM_S[INDEX2(0,0,4)]+=tmp1_1;
                            EM_S[INDEX2(1,0,4)]+=tmp0_1;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1;
                            EM_S[INDEX2(1,1,4)]+=tmp1_1;
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        if (y.actsExpanded()) {
                            const double y_0 = y_p[0];
                            const double y_1 = y_p[1];
                            const double tmp2_1 = w13*y_1;
                            const double tmp1_1 = w14*y_1;
                            const double tmp3_1 = w14*y_0;
                            const double tmp0_1 = w13*y_0;
                            EM_F[0]+=tmp0_1 + tmp1_1;
                            EM_F[1]+=tmp2_1 + tmp3_1;
                        } else { /* constant data */
                            const double y_0 = y_p[0];
                            const double tmp0_1 = w15*y_0;
                            EM_F[0]+=tmp0_1;
                            EM_F[1]+=tmp0_1;
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
                            const double tmp0_0 = d_0 + d_1;
                            const double tmp2_1 = d_1*w9;
                            const double tmp4_1 = d_0*w9;
                            const double tmp0_1 = tmp0_0*w8;
                            const double tmp3_1 = d_1*w10;
                            const double tmp1_1 = d_0*w10;
                            EM_S[INDEX2(2,2,4)]+=tmp3_1 + tmp4_1;
                            EM_S[INDEX2(3,2,4)]+=tmp0_1;
                            EM_S[INDEX2(2,3,4)]+=tmp0_1;
                            EM_S[INDEX2(3,3,4)]+=tmp1_1 + tmp2_1;
                        } else { /* constant data */
                            const double d_0 = d_p[0];
                            const double tmp0_1 = d_0*w11;
                            const double tmp1_1 = d_0*w12;
                            EM_S[INDEX2(2,2,4)]+=tmp1_1;
                            EM_S[INDEX2(3,2,4)]+=tmp0_1;
                            EM_S[INDEX2(2,3,4)]+=tmp0_1;
                            EM_S[INDEX2(3,3,4)]+=tmp1_1;
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        if (y.actsExpanded()) {
                            const double y_0 = y_p[0];
                            const double y_1 = y_p[1];
                            const double tmp2_1 = w13*y_1;
                            const double tmp1_1 = w14*y_1;
                            const double tmp3_1 = w14*y_0;
                            const double tmp0_1 = w13*y_0;
                            EM_F[2]+=tmp0_1 + tmp1_1;
                            EM_F[3]+=tmp2_1 + tmp3_1;
                        } else { /* constant data */
                            const double y_0 = y_p[0];
                            const double tmp0_1 = w15*y_0;
                            EM_F[2]+=tmp0_1;
                            EM_F[3]+=tmp0_1;
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
    const double w0 = 0.25*m_dx[1];
    const double w1 = 0.5*m_dx[1];
    const double w2 = 0.25*m_dx[0];
    const double w3 = 0.5*m_dx[0];
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
                    const index_t e = k1;
                    ///////////////
                    // process d //
                    ///////////////
                    if (add_EM_S) {
                        const double* d_p=const_cast<escript::Data*>(&d)->getSampleDataRO(e);
                        const double d_0 = d_p[0];
                        const double tmp0_1 = d_0*w0;
                        EM_S[INDEX2(0,0,4)]+=tmp0_1;
                        EM_S[INDEX2(2,0,4)]+=tmp0_1;
                        EM_S[INDEX2(0,2,4)]+=tmp0_1;
                        EM_S[INDEX2(2,2,4)]+=tmp0_1;
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        const double tmp0_1 = w1*y_p[0];
                        EM_F[0]+=tmp0_1;
                        EM_F[2]+=tmp0_1;
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
                        const double d_0 = d_p[0];
                        const double tmp0_1 = d_0*w0;
                        EM_S[INDEX2(1,1,4)]+=tmp0_1;
                        EM_S[INDEX2(3,1,4)]+=tmp0_1;
                        EM_S[INDEX2(1,3,4)]+=tmp0_1;
                        EM_S[INDEX2(3,3,4)]+=tmp0_1;
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        const double tmp0_1 = w1*y_p[0];
                        EM_F[1]+=tmp0_1;
                        EM_F[3]+=tmp0_1;
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
                        const double tmp0_1 = d_p[0]*w2;
                        EM_S[INDEX2(0,0,4)]+=tmp0_1;
                        EM_S[INDEX2(1,0,4)]+=tmp0_1;
                        EM_S[INDEX2(0,1,4)]+=tmp0_1;
                        EM_S[INDEX2(1,1,4)]+=tmp0_1;
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        const double tmp0_1 = w3*y_p[0];
                        EM_F[0]+=tmp0_1;
                        EM_F[1]+=tmp0_1;
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
                        const double tmp0_1 = d_p[0]*w2;
                        EM_S[INDEX2(2,2,4)]+=tmp0_1;
                        EM_S[INDEX2(3,2,4)]+=tmp0_1;
                        EM_S[INDEX2(2,3,4)]+=tmp0_1;
                        EM_S[INDEX2(3,3,4)]+=tmp0_1;
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        const double tmp0_1 = w3*y_p[0];
                        EM_F[2]+=tmp0_1;
                        EM_F[3]+=tmp0_1;
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
    const double w0 = 0.31100423396407310779*m_dx[1];
    const double w1 = 0.022329099369260225539*m_dx[1];
    const double w10 = 0.022329099369260225539*m_dx[0];
    const double w11 = 0.16666666666666666667*m_dx[0];
    const double w12 = 0.33333333333333333333*m_dx[0];
    const double w13 = 0.39433756729740644113*m_dx[0];
    const double w14 = 0.10566243270259355887*m_dx[0];
    const double w15 = 0.5*m_dx[0];
    const double w2 = 0.083333333333333333333*m_dx[1];
    const double w3 = 0.33333333333333333333*m_dx[1];
    const double w4 = 0.16666666666666666667*m_dx[1];
    const double w5 = 0.39433756729740644113*m_dx[1];
    const double w6 = 0.10566243270259355887*m_dx[1];
    const double w7 = 0.5*m_dx[1];
    const double w8 = 0.083333333333333333333*m_dx[0];
    const double w9 = 0.31100423396407310779*m_dx[0];
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
                                    const double d_0 = d_p[INDEX3(k, m, 0, numEq, numComp)];
                                    const double d_1 = d_p[INDEX3(k, m, 1, numEq, numComp)];
                                    const double tmp0_0 = d_0 + d_1;
                                    const double tmp0_1 = d_0*w0;
                                    const double tmp1_1 = d_1*w1;
                                    const double tmp2_1 = tmp0_0*w2;
                                    const double tmp3_1 = d_1*w0;
                                    const double tmp4_1 = d_0*w1;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp0_1 + tmp1_1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp2_1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp2_1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp3_1 + tmp4_1;
                                }
                             }
                        } else { /* constant data */
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX2(k, m, numEq)];
                                    const double tmp0_1 = d_0*w3;
                                    const double tmp1_1 = d_0*w4;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp1_1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp1_1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp0_1;
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
                                const double tmp3_1 = w5*y_1;
                                const double tmp2_1 = w6*y_0;
                                const double tmp0_1 = w6*y_1;
                                const double tmp1_1 = w5*y_0;
                                EM_F[INDEX2(k,0,numEq)]+=tmp0_1 + tmp1_1;
                                EM_F[INDEX2(k,2,numEq)]+=tmp2_1 + tmp3_1;
                            }
                        } else { /* constant data */
                            for (index_t k=0; k<numEq; k++) {
                                const double y_0 = y_p[k];
                                const double tmp0_1 = w7*y_0;
                                EM_F[INDEX2(k,0,numEq)]+=tmp0_1;
                                EM_F[INDEX2(k,2,numEq)]+=tmp0_1;
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
                                    const double d_0 = d_p[INDEX3(k, m, 0, numEq, numComp)];
                                    const double d_1 = d_p[INDEX3(k, m, 1, numEq, numComp)];
                                    const double tmp0_0 = d_0 + d_1;
                                    const double tmp4_1 = d_1*w1;
                                    const double tmp1_1 = d_0*w1;
                                    const double tmp3_1 = d_0*w0;
                                    const double tmp0_1 = d_1*w0;
                                    const double tmp2_1 = tmp0_0*w2;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp2_1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp2_1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp0_1 + tmp1_1;
                                }
                             }
                        } else { /* constant data */
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX2(k, m, numEq)];
                                    const double tmp1_1 = d_0*w4;
                                    const double tmp0_1 = d_0*w3;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp1_1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp1_1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp0_1;
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
                                const double tmp3_1 = w5*y_1;
                                const double tmp2_1 = w6*y_0;
                                const double tmp0_1 = w6*y_1;
                                const double tmp1_1 = w5*y_0;
                                EM_F[INDEX2(k,1,numEq)]+=tmp0_1 + tmp1_1;
                                EM_F[INDEX2(k,3,numEq)]+=tmp2_1 + tmp3_1;
                            }
                        } else { /* constant data */
                            for (index_t k=0; k<numEq; k++) {
                                const double y_0 = y_p[k];
                                const double tmp0_1 = w7*y_0;
                                EM_F[INDEX2(k,1,numEq)]+=tmp0_1;
                                EM_F[INDEX2(k,3,numEq)]+=tmp0_1;
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
                                    const double d_0 = d_p[INDEX3(k, m, 0, numEq, numComp)];
                                    const double d_1 = d_p[INDEX3(k, m, 1, numEq, numComp)];
                                    const double tmp0_0 = d_0 + d_1;
                                    const double tmp4_1 = d_1*w9;
                                    const double tmp2_1 = d_0*w9;
                                    const double tmp0_1 = tmp0_0*w8;
                                    const double tmp1_1 = d_1*w10;
                                    const double tmp3_1 = d_0*w10;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp3_1 + tmp4_1;
                                }
                             }
                        } else { /* constant data */
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX2(k, m, numEq)];
                                    const double tmp0_1 = d_0*w11;
                                    const double tmp1_1 = d_0*w12;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp1_1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp1_1;
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
                                const double tmp2_1 = w13*y_1;
                                const double tmp1_1 = w14*y_1;
                                const double tmp3_1 = w14*y_0;
                                const double tmp0_1 = w13*y_0;
                                EM_F[INDEX2(k,0,numEq)]+=tmp0_1 + tmp1_1;
                                EM_F[INDEX2(k,1,numEq)]+=tmp2_1 + tmp3_1;
                            }
                        } else { /* constant data */
                            for (index_t k=0; k<numEq; k++) {
                                const double y_0 = y_p[k];
                                const double tmp0_1 = w15*y_0;
                                EM_F[INDEX2(k,0,numEq)]+=tmp0_1;
                                EM_F[INDEX2(k,1,numEq)]+=tmp0_1;
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
                                    const double d_0 = d_p[INDEX3(k, m, 0, numEq, numComp)];
                                    const double d_1 = d_p[INDEX3(k, m, 1, numEq, numComp)];
                                    const double tmp0_0 = d_0 + d_1;
                                    const double tmp2_1 = d_1*w9;
                                    const double tmp4_1 = d_0*w9;
                                    const double tmp0_1 = tmp0_0*w8;
                                    const double tmp3_1 = d_1*w10;
                                    const double tmp1_1 = d_0*w10;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp1_1 + tmp2_1;
                                }
                             }
                        } else { /* constant data */
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX2(k, m, numEq)];
                                    const double tmp0_1 = d_0*w11;
                                    const double tmp1_1 = d_0*w12;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp1_1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp1_1;
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
                                const double tmp2_1 = w13*y_1;
                                const double tmp1_1 = w14*y_1;
                                const double tmp3_1 = w14*y_0;
                                const double tmp0_1 = w13*y_0;
                                EM_F[INDEX2(k,2,numEq)]+=tmp0_1 + tmp1_1;
                                EM_F[INDEX2(k,3,numEq)]+=tmp2_1 + tmp3_1;
                            }
                        } else { /* constant data */
                            for (index_t k=0; k<numEq; k++) {
                                const double y_0 = y_p[k];
                                const double tmp0_1 = w15*y_0;
                                EM_F[INDEX2(k,2,numEq)]+=tmp0_1;
                                EM_F[INDEX2(k,3,numEq)]+=tmp0_1;
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
    const double w0 = 0.25*m_dx[1];
    const double w1 = 0.5*m_dx[1];
    const double w2 = 0.25*m_dx[0];
    const double w3 = 0.5*m_dx[0];
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
                        for (index_t k=0; k<numEq; k++) {
                            for (index_t m=0; m<numComp; m++) {
                                const double d_0 = d_p[INDEX2(k, m, numEq)];
                                const double tmp0_1 = d_0*w0;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp0_1;
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        for (index_t k=0; k<numEq; k++) {
                            const double tmp0_1 = w1*y_p[k];
                            EM_F[INDEX2(k,0,numEq)]+=tmp0_1;
                            EM_F[INDEX2(k,2,numEq)]+=tmp0_1;
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
                                const double d_0 = d_p[INDEX2(k, m, numEq)];
                                const double tmp0_1 = d_0*w0;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp0_1;
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        for (index_t k=0; k<numEq; k++) {
                            const double tmp0_1 = w1*y_p[k];
                            EM_F[INDEX2(k,1,numEq)]+=tmp0_1;
                            EM_F[INDEX2(k,3,numEq)]+=tmp0_1;
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
                                const double d_0 = d_p[INDEX2(k, m, numEq)];
                                const double tmp0_1 = d_0*w2;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp0_1;
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        for (index_t k=0; k<numEq; k++) {
                            const double tmp0_1 = w3*y_p[k];
                            EM_F[INDEX2(k,0,numEq)]+=tmp0_1;
                            EM_F[INDEX2(k,1,numEq)]+=tmp0_1;
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
                                const double d_0 = d_p[INDEX2(k, m, numEq)];
                                const double tmp0_1 = d_0*w2;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp0_1;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp0_1;
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (add_EM_F) {
                        const double* y_p=const_cast<escript::Data*>(&y)->getSampleDataRO(e);
                        for (index_t k=0; k<numEq; k++) {
                            const double tmp0_1 = w3*y_p[k];
                            EM_F[INDEX2(k,2,numEq)]+=tmp0_1;
                            EM_F[INDEX2(k,3,numEq)]+=tmp0_1;
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


} // end of namespace ripley

