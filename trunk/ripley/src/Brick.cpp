
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

#include <ripley/Brick.h>
#include <paso/SystemMatrix.h>
#include <esysUtils/esysFileWriter.h>
#include <ripley/DefaultAssembler3D.h>
#include <ripley/WaveAssembler3D.h>
#include <boost/scoped_array.hpp>

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

#include "esysUtils/EsysRandom.h"
#include "blocktools.h"


using namespace std;
using esysUtils::FileWriter;

namespace ripley {

Brick::Brick(int n0, int n1, int n2, double x0, double y0, double z0,
             double x1, double y1, double z1, int d0, int d1, int d2,
             const std::vector<double>& points, const std::vector<int>& tags,
             const simap_t& tagnamestonums) :
    RipleyDomain(3)
{
    // ignore subdivision parameters for serial run
    if (m_mpiInfo->size == 1) {
        d0=1;
        d1=1;
        d2=1;
    }
    bool warn=false;
    // if number of subdivisions is non-positive, try to subdivide by the same
    // ratio as the number of elements
    if (d0<=0 && d1<=0 && d2<=0) {
        warn=true;
        d0=(int)(pow(m_mpiInfo->size*(n0+1)*(n0+1)/((float)(n1+1)*(n2+1)), 1.f/3));
        d0=max(1, d0);
        d1=max(1, (int)(d0*n1/(float)n0));
        d2=m_mpiInfo->size/(d0*d1);
        if (d0*d1*d2 != m_mpiInfo->size) {
            // ratios not the same so leave "smallest" side undivided and try
            // dividing 2 sides only
            if (n0>=n1) {
                if (n1>=n2) {
                    d0=d1=0;
                    d2=1;
                } else {
                    d0=d2=0;
                    d1=1;
                }
            } else {
                if (n0>=n2) {
                    d0=d1=0;
                    d2=1;
                } else {
                    d0=1;
                    d1=d2=0;
                }
            }
        }
    }
    if (d0<=0 && d1<=0) {
        warn=true;
        d0=max(1, int(sqrt(m_mpiInfo->size*(n0+1)/(float)(n1+1))));
        d1=m_mpiInfo->size/d0;
        if (d0*d1*d2 != m_mpiInfo->size) {
            // ratios not the same so subdivide side with more elements only
            if (n0>n1) {
                d0=0;
                d1=1;
            } else {
                d0=1;
                d1=0;
            }
        }
    } else if (d0<=0 && d2<=0) {
        warn=true;
        d0=max(1, int(sqrt(m_mpiInfo->size*(n0+1)/(float)(n2+1))));
        d2=m_mpiInfo->size/d0;
        if (d0*d1*d2 != m_mpiInfo->size) {
            // ratios not the same so subdivide side with more elements only
            if (n0>n2) {
                d0=0;
                d2=1;
            } else {
                d0=1;
                d2=0;
            }
        }
    } else if (d1<=0 && d2<=0) {
        warn=true;
        d1=max(1, int(sqrt(m_mpiInfo->size*(n1+1)/(float)(n2+1))));
        d2=m_mpiInfo->size/d1;
        if (d0*d1*d2 != m_mpiInfo->size) {
            // ratios not the same so subdivide side with more elements only
            if (n1>n2) {
                d1=0;
                d2=1;
            } else {
                d1=1;
                d2=0;
            }
        }
    }
    if (d0<=0) {
        // d1,d2 are preset, determine d0
        d0=m_mpiInfo->size/(d1*d2);
    } else if (d1<=0) {
        // d0,d2 are preset, determine d1
        d1=m_mpiInfo->size/(d0*d2);
    } else if (d2<=0) {
        // d0,d1 are preset, determine d2
        d2=m_mpiInfo->size/(d0*d1);
    }

    // ensure number of subdivisions is valid and nodes can be distributed
    // among number of ranks
    if (d0*d1*d2 != m_mpiInfo->size)
        throw RipleyException("Invalid number of spatial subdivisions");

    if (warn) {
        cout << "Warning: Automatic domain subdivision (d0=" << d0 << ", d1="
            << d1 << ", d2=" << d2 << "). This may not be optimal!" << endl;
    }

    double l0 = x1-x0;
    double l1 = y1-y0;
    double l2 = z1-z0;
    m_dx[0] = l0/n0;
    m_dx[1] = l1/n1;
    m_dx[2] = l2/n2;

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
    if ((n2+1)%d2 > 0) {
        n2=(int)round((float)(n2+1)/d2+0.5)*d2-1;
        l2=m_dx[2]*n2;
        cout << "Warning: Adjusted number of elements and length. N2="
            << n2 << ", l2=" << l2 << endl;
    }

    if ((d0 > 1 && (n0+1)/d0<2) || (d1 > 1 && (n1+1)/d1<2) || (d2 > 1 && (n2+1)/d2<2))
        throw RipleyException("Too few elements for the number of ranks");

    m_gNE[0] = n0;
    m_gNE[1] = n1;
    m_gNE[2] = n2;
    m_origin[0] = x0;
    m_origin[1] = y0;
    m_origin[2] = z0;
    m_length[0] = l0;
    m_length[1] = l1;
    m_length[2] = l2;
    m_NX[0] = d0;
    m_NX[1] = d1;
    m_NX[2] = d2;

    // local number of elements (including overlap)
    m_NE[0] = m_ownNE[0] = (d0>1 ? (n0+1)/d0 : n0);
    if (m_mpiInfo->rank%d0>0 && m_mpiInfo->rank%d0<d0-1)
        m_NE[0]++;
    else if (d0>1 && m_mpiInfo->rank%d0==d0-1)
        m_ownNE[0]--;

    m_NE[1] = m_ownNE[1] = (d1>1 ? (n1+1)/d1 : n1);
    if (m_mpiInfo->rank%(d0*d1)/d0>0 && m_mpiInfo->rank%(d0*d1)/d0<d1-1)
        m_NE[1]++;
    else if (d1>1 && m_mpiInfo->rank%(d0*d1)/d0==d1-1)
        m_ownNE[1]--;

    m_NE[2] = m_ownNE[2] = (d2>1 ? (n2+1)/d2 : n2);
    if (m_mpiInfo->rank/(d0*d1)>0 && m_mpiInfo->rank/(d0*d1)<d2-1)
        m_NE[2]++;
    else if (d2>1 && m_mpiInfo->rank/(d0*d1)==d2-1)
        m_ownNE[2]--;

    // local number of nodes
    m_NN[0] = m_NE[0]+1;
    m_NN[1] = m_NE[1]+1;
    m_NN[2] = m_NE[2]+1;

    // bottom-left-front node is at (offset0,offset1,offset2) in global mesh
    m_offset[0] = (n0+1)/d0*(m_mpiInfo->rank%d0);
    if (m_offset[0] > 0)
        m_offset[0]--;
    m_offset[1] = (n1+1)/d1*(m_mpiInfo->rank%(d0*d1)/d0);
    if (m_offset[1] > 0)
        m_offset[1]--;
    m_offset[2] = (n2+1)/d2*(m_mpiInfo->rank/(d0*d1));
    if (m_offset[2] > 0)
        m_offset[2]--;

    populateSampleIds();
    createPattern();
    
    assembler = new DefaultAssembler3D(this, m_dx, m_NX, m_NE, m_NN);
    for (map<string, int>::const_iterator i = tagnamestonums.begin();
            i != tagnamestonums.end(); i++) {
        setTagMap(i->first, i->second);
    }
    addPoints(tags.size(), &points[0], &tags[0]);
}


Brick::~Brick()
{
    Paso_SystemMatrixPattern_free(m_pattern);
    Paso_Connector_free(m_connector);
    delete assembler;
}

string Brick::getDescription() const
{
    return "ripley::Brick";
}

bool Brick::operator==(const AbstractDomain& other) const
{
    const Brick* o=dynamic_cast<const Brick*>(&other);
    if (o) {
        return (RipleyDomain::operator==(other) &&
                m_gNE[0]==o->m_gNE[0] && m_gNE[1]==o->m_gNE[1] && m_gNE[2]==o->m_gNE[2]
                && m_origin[0]==o->m_origin[0] && m_origin[1]==o->m_origin[1] && m_origin[2]==o->m_origin[2]
                && m_length[0]==o->m_length[0] && m_length[1]==o->m_length[1] && m_length[2]==o->m_length[2]
                && m_NX[0]==o->m_NX[0] && m_NX[1]==o->m_NX[1] && m_NX[2]==o->m_NX[2]);
    }

    return false;
}

void Brick::readNcGrid(escript::Data& out, string filename, string varname,
            const ReaderParameters& params) const
{
#ifdef USE_NETCDF
    // check destination function space
    int myN0, myN1, myN2;
    if (out.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NN[0];
        myN1 = m_NN[1];
        myN2 = m_NN[2];
    } else if (out.getFunctionSpace().getTypeCode() == Elements ||
                out.getFunctionSpace().getTypeCode() == ReducedElements) {
        myN0 = m_NE[0];
        myN1 = m_NE[1];
        myN2 = m_NE[2];
    } else
        throw RipleyException("readNcGrid(): invalid function space for output data object");

    if (params.first.size() != 3)
        throw RipleyException("readNcGrid(): argument 'first' must have 3 entries");

    if (params.numValues.size() != 3)
        throw RipleyException("readNcGrid(): argument 'numValues' must have 3 entries");

    if (params.multiplier.size() != 3)
        throw RipleyException("readNcGrid(): argument 'multiplier' must have 3 entries");
    for (size_t i=0; i<params.multiplier.size(); i++)
        if (params.multiplier[i]<1)
            throw RipleyException("readNcGrid(): all multipliers must be positive");

    // check file existence and size
    NcFile f(filename.c_str(), NcFile::ReadOnly);
    if (!f.is_valid())
        throw RipleyException("readNcGrid(): cannot open file");

    NcVar* var = f.get_var(varname.c_str());
    if (!var)
        throw RipleyException("readNcGrid(): invalid variable name");

    // TODO: rank>0 data support
    const int numComp = out.getDataPointSize();
    if (numComp > 1)
        throw RipleyException("readNcGrid(): only scalar data supported");

    const int dims = var->num_dims();
    boost::scoped_array<long> edges(var->edges());

    // is this a slice of the data object (dims!=3)?
    // note the expected ordering of edges (as in numpy: z,y,x)
    if ( (dims==3 && (params.numValues[2] > edges[0] ||
                      params.numValues[1] > edges[1] ||
                      params.numValues[0] > edges[2]))
            || (dims==2 && params.numValues[2]>1)
            || (dims==1 && (params.numValues[2]>1 || params.numValues[1]>1)) ) {
        throw RipleyException("readNcGrid(): not enough data in file");
    }

    // check if this rank contributes anything
    if (params.first[0] >= m_offset[0]+myN0 ||
            params.first[0]+params.numValues[0]*params.multiplier[0] <= m_offset[0] ||
            params.first[1] >= m_offset[1]+myN1 ||
            params.first[1]+params.numValues[1]*params.multiplier[1] <= m_offset[1] ||
            params.first[2] >= m_offset[2]+myN2 ||
            params.first[2]+params.numValues[2]*params.multiplier[2] <= m_offset[2]) {
        return;
    }

    // now determine how much this rank has to write

    // first coordinates in data object to write to
    const int first0 = max(0, params.first[0]-m_offset[0]);
    const int first1 = max(0, params.first[1]-m_offset[1]);
    const int first2 = max(0, params.first[2]-m_offset[2]);
    // indices to first value in file (not accounting for reverse yet)
    int idx0 = max(0, m_offset[0]-params.first[0]);
    int idx1 = max(0, m_offset[1]-params.first[1]);
    int idx2 = max(0, m_offset[2]-params.first[2]);
    // number of values to read
    const int num0 = min(params.numValues[0]-idx0, myN0-first0);
    const int num1 = min(params.numValues[1]-idx1, myN1-first1);
    const int num2 = min(params.numValues[2]-idx2, myN2-first2);

    // make sure we read the right block if going backwards through file
    if (params.reverse[0])
        idx0 = edges[dims-1]-num0-idx0;
    if (dims>1 && params.reverse[1])
        idx1 = edges[dims-2]-num1-idx1;
    if (dims>2 && params.reverse[2])
        idx2 = edges[dims-3]-num2-idx2;


    vector<double> values(num0*num1*num2);
    if (dims==3) {
        var->set_cur(idx2, idx1, idx0);
        var->get(&values[0], num2, num1, num0);
    } else if (dims==2) {
        var->set_cur(idx1, idx0);
        var->get(&values[0], num1, num0);
    } else {
        var->set_cur(idx0);
        var->get(&values[0], num0);
    }

    const int dpp = out.getNumDataPointsPerSample();
    out.requireWrite();

    // helpers for reversing
    const int x0 = (params.reverse[0] ? num0-1 : 0);
    const int x_mult = (params.reverse[0] ? -1 : 1);
    const int y0 = (params.reverse[1] ? num1-1 : 0);
    const int y_mult = (params.reverse[1] ? -1 : 1);
    const int z0 = (params.reverse[2] ? num2-1 : 0);
    const int z_mult = (params.reverse[2] ? -1 : 1);

    for (index_t z=0; z<num2; z++) {
        for (index_t y=0; y<num1; y++) {
#pragma omp parallel for
            for (index_t x=0; x<num0; x++) {
                const int baseIndex = first0+x*params.multiplier[0]
                                     +(first1+y*params.multiplier[1])*myN0
                                     +(first2+z*params.multiplier[2])*myN0*myN1;
                const int srcIndex=(z0+z_mult*z)*num1*num0
                                  +(y0+y_mult*y)*num0
                                  +(x0+x_mult*x);
                if (!isnan(values[srcIndex])) {
                    for (index_t m2=0; m2<params.multiplier[2]; m2++) {
                        for (index_t m1=0; m1<params.multiplier[1]; m1++) {
                            for (index_t m0=0; m0<params.multiplier[0]; m0++) {
                                const int dataIndex = baseIndex+m0
                                               +m1*myN0
                                               +m2*myN0*myN1;
                                double* dest = out.getSampleDataRW(dataIndex);
                                for (index_t q=0; q<dpp; q++) {
                                    *dest++ = values[srcIndex];
                                }
                            }
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

void Brick::readBinaryGrid(escript::Data& out, string filename,
                           const ReaderParameters& params) const
{
    // the mapping is not universally correct but should work on our
    // supported platforms
    switch (params.dataType) {
        case DATATYPE_INT32:
            readBinaryGridImpl<int>(out, filename, params);
            break;
        case DATATYPE_FLOAT32:
            readBinaryGridImpl<float>(out, filename, params);
            break;
        case DATATYPE_FLOAT64:
            readBinaryGridImpl<double>(out, filename, params);
            break;
        default:
            throw RipleyException("readBinaryGrid(): invalid or unsupported datatype");
    }
}

template<typename ValueType>
void Brick::readBinaryGridImpl(escript::Data& out, const string& filename,
                               const ReaderParameters& params) const
{
    // check destination function space
    int myN0, myN1, myN2;
    if (out.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NN[0];
        myN1 = m_NN[1];
        myN2 = m_NN[2];
    } else if (out.getFunctionSpace().getTypeCode() == Elements ||
                out.getFunctionSpace().getTypeCode() == ReducedElements) {
        myN0 = m_NE[0];
        myN1 = m_NE[1];
        myN2 = m_NE[2];
    } else
        throw RipleyException("readBinaryGrid(): invalid function space for output data object");

    if (params.first.size() != 3)
        throw RipleyException("readBinaryGrid(): argument 'first' must have 3 entries");

    if (params.numValues.size() != 3)
        throw RipleyException("readBinaryGrid(): argument 'numValues' must have 3 entries");

    if (params.multiplier.size() != 3)
        throw RipleyException("readBinaryGrid(): argument 'multiplier' must have 3 entries");
    for (size_t i=0; i<params.multiplier.size(); i++)
        if (params.multiplier[i]<1)
            throw RipleyException("readBinaryGrid(): all multipliers must be positive");

    // check file existence and size
    ifstream f(filename.c_str(), ifstream::binary);
    if (f.fail()) {
        throw RipleyException("readBinaryGrid(): cannot open file");
    }
    f.seekg(0, ios::end);
    const int numComp = out.getDataPointSize();
    const int filesize = f.tellg();
    const int reqsize = params.numValues[0]*params.numValues[1]*params.numValues[2]*numComp*sizeof(ValueType);
    if (filesize < reqsize) {
        f.close();
        throw RipleyException("readBinaryGrid(): not enough data in file");
    }

    // check if this rank contributes anything
    if (params.first[0] >= m_offset[0]+myN0 ||
            params.first[0]+params.numValues[0]*params.multiplier[0] <= m_offset[0] ||
            params.first[1] >= m_offset[1]+myN1 ||
            params.first[1]+params.numValues[1]*params.multiplier[1] <= m_offset[1] ||
            params.first[2] >= m_offset[2]+myN2 ||
            params.first[2]+params.numValues[2]*params.multiplier[2] <= m_offset[2]) {
        f.close();
        return;
    }

    // now determine how much this rank has to write

    // first coordinates in data object to write to
    const int first0 = max(0, params.first[0]-m_offset[0]);
    const int first1 = max(0, params.first[1]-m_offset[1]);
    const int first2 = max(0, params.first[2]-m_offset[2]);
    // indices to first value in file
    const int idx0 = max(0, m_offset[0]-params.first[0]);
    const int idx1 = max(0, m_offset[1]-params.first[1]);
    const int idx2 = max(0, m_offset[2]-params.first[2]);
    // number of values to read
    const int num0 = min(params.numValues[0]-idx0, myN0-first0);
    const int num1 = min(params.numValues[1]-idx1, myN1-first1);
    const int num2 = min(params.numValues[2]-idx2, myN2-first2);

    out.requireWrite();
    vector<ValueType> values(num0*numComp);
    const int dpp = out.getNumDataPointsPerSample();

    for (int z=0; z<num2; z++) {
        for (int y=0; y<num1; y++) {
            const int fileofs = numComp*(idx0+(idx1+y)*params.numValues[0]
                             +(idx2+z)*params.numValues[0]*params.numValues[1]);
            f.seekg(fileofs*sizeof(ValueType));
            f.read((char*)&values[0], num0*numComp*sizeof(ValueType));

            for (int x=0; x<num0; x++) {
                const int baseIndex = first0+x*params.multiplier[0]
                                     +(first1+y*params.multiplier[1])*myN0
                                     +(first2+z*params.multiplier[2])*myN0*myN1;
                for (int m2=0; m2<params.multiplier[2]; m2++) {
                    for (int m1=0; m1<params.multiplier[1]; m1++) {
                        for (int m0=0; m0<params.multiplier[0]; m0++) {
                            const int dataIndex = baseIndex+m0
                                           +m1*myN0
                                           +m2*myN0*myN1;
                            double* dest = out.getSampleDataRW(dataIndex);
                            for (int c=0; c<numComp; c++) {
                                ValueType val = values[x*numComp+c];

                                if (params.byteOrder != BYTEORDER_NATIVE) {
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
        }
    }

    f.close();
}

void Brick::writeBinaryGrid(const escript::Data& in, string filename,
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
void Brick::writeBinaryGridImpl(const escript::Data& in,
                                const string& filename, int byteOrder) const
{
    // check function space and determine number of points
    int myN0, myN1, myN2;
    int totalN0, totalN1, totalN2;
    if (in.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NN[0];
        myN1 = m_NN[1];
        myN2 = m_NN[2];
        totalN0 = m_gNE[0]+1;
        totalN1 = m_gNE[1]+1;
        totalN2 = m_gNE[2]+1;
    } else if (in.getFunctionSpace().getTypeCode() == Elements ||
                in.getFunctionSpace().getTypeCode() == ReducedElements) {
        myN0 = m_NE[0];
        myN1 = m_NE[1];
        myN2 = m_NE[2];
        totalN0 = m_gNE[0];
        totalN1 = m_gNE[1];
        totalN2 = m_gNE[2];
    } else
        throw RipleyException("writeBinaryGrid(): invalid function space of data object");

    const int numComp = in.getDataPointSize();
    const int dpp = in.getNumDataPointsPerSample();
    const int fileSize = sizeof(ValueType)*numComp*dpp*totalN0*totalN1*totalN2;

    if (numComp > 1 || dpp > 1)
        throw RipleyException("writeBinaryGrid(): only scalar, single-value data supported");

    // from here on we know that each sample consists of one value
    FileWriter fw;
    fw.openFile(filename, fileSize);
    MPIBarrier();

    for (index_t z=0; z<myN2; z++) {
        for (index_t y=0; y<myN1; y++) {
            const int fileofs = (m_offset[0]+(m_offset[1]+y)*totalN0
                                +(m_offset[2]+z)*totalN0*totalN1)*sizeof(ValueType);
            ostringstream oss;

            for (index_t x=0; x<myN0; x++) {
                const double* sample = in.getSampleDataRO(z*myN0*myN1+y*myN0+x);
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
    }
    fw.close();
}

void Brick::dump(const string& fileName) const
{
#if USE_SILO
    string fn(fileName);
    if (fileName.length() < 6 || fileName.compare(fileName.length()-5, 5, ".silo") != 0) {
        fn+=".silo";
    }

    int driver=DB_HDF5;    
    string siloPath;
    DBfile* dbfile = NULL;

#ifdef ESYS_MPI
    PMPIO_baton_t* baton = NULL;
    const int NUM_SILO_FILES = 1;
    const char* blockDirFmt = "/block%04d";
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
            char str[64];
            snprintf(str, 64, blockDirFmt, PMPIO_RankInGroup(baton, m_mpiInfo->rank));
            siloPath = str;
            dbfile = (DBfile*) PMPIO_WaitForBaton(baton, fn.c_str(), siloPath.c_str());
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
    boost::scoped_ptr<double> z(new double[m_NN[2]]);
    double* coords[3] = { x.get(), y.get(), z.get() };
#pragma omp parallel
    {
#pragma omp for
        for (dim_t i0 = 0; i0 < m_NN[0]; i0++) {
            coords[0][i0]=getLocalCoordinate(i0, 0);
        }
#pragma omp for
        for (dim_t i1 = 0; i1 < m_NN[1]; i1++) {
            coords[1][i1]=getLocalCoordinate(i1, 1);
        }
#pragma omp for
        for (dim_t i2 = 0; i2 < m_NN[2]; i2++) {
            coords[2][i2]=getLocalCoordinate(i2, 2);
        }
    }
    int* dims = const_cast<int*>(getNumNodesPerDim());

    // write mesh
    DBPutQuadmesh(dbfile, "mesh", NULL, coords, dims, 3, DB_DOUBLE,
            DB_COLLINEAR, NULL);

    // write node ids
    DBPutQuadvar1(dbfile, "nodeId", "mesh", (void*)&m_nodeId[0], dims, 3,
            NULL, 0, DB_INT, DB_NODECENT, NULL);

    // write element ids
    dims = const_cast<int*>(getNumElementsPerDim());
    DBPutQuadvar1(dbfile, "elementId", "mesh", (void*)&m_elementId[0],
            dims, 3, NULL, 0, DB_INT, DB_ZONECENT, NULL);

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

const int* Brick::borrowSampleReferenceIDs(int fsType) const
{
    switch (fsType) {
        case Nodes:
        case ReducedNodes: //FIXME: reduced
            return &m_nodeId[0];
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom: //FIXME: reduced
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
    msg << "borrowSampleReferenceIDs: invalid function space type "<<fsType;
    throw RipleyException(msg.str());
}

bool Brick::ownSample(int fsType, index_t id) const
{
    if (getMPISize()==1)
        return true;

    switch (fsType) {
        case Nodes:
        case ReducedNodes: //FIXME: reduced
            return (m_dofMap[id] < getNumDOF());
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
            return true;
        case Elements:
        case ReducedElements:
            {
                // check ownership of element's _last_ node
                const index_t x=id%m_NE[0] + 1;
                const index_t y=id%(m_NE[0]*m_NE[1])/m_NE[0] + 1;
                const index_t z=id/(m_NE[0]*m_NE[1]) + 1;
                return (m_dofMap[x + m_NN[0]*y + m_NN[0]*m_NN[1]*z] < getNumDOF());
            }
        case FaceElements:
        case ReducedFaceElements:
            {
                // check ownership of face element's last node
                dim_t n=0;
                for (size_t i=0; i<6; i++) {
                    n+=m_faceCount[i];
                    if (id<n) {
                        const index_t j=id-n+m_faceCount[i];
                        if (i>=4) { // front or back
                            const index_t first=(i==4 ? 0 : m_NN[0]*m_NN[1]*(m_NN[2]-1));
                            return (m_dofMap[first+j%m_NE[0]+1+(j/m_NE[0]+1)*m_NN[0]] < getNumDOF());
                        } else if (i>=2) { // bottom or top
                            const index_t first=(i==2 ? 0 : m_NN[0]*(m_NN[1]-1));
                            return (m_dofMap[first+j%m_NE[0]+1+(j/m_NE[0]+1)*m_NN[0]*m_NN[1]] < getNumDOF());
                        } else { // left or right
                            const index_t first=(i==0 ? 0 : m_NN[0]-1);
                            return (m_dofMap[first+(j%m_NE[1]+1)*m_NN[0]+(j/m_NE[1]+1)*m_NN[0]*m_NN[1]] < getNumDOF());
                        }
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

void Brick::setToNormal(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == FaceElements) {
        out.requireWrite();
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE[2]; ++k2) {
                    for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                        double* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]));
                        // set vector at four quadrature points
                        *o++ = -1.; *o++ = 0.; *o++ = 0.;
                        *o++ = -1.; *o++ = 0.; *o++ = 0.;
                        *o++ = -1.; *o++ = 0.; *o++ = 0.;
                        *o++ = -1.; *o++ = 0.; *o = 0.;
                    }
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE[2]; ++k2) {
                    for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                        double* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]));
                        // set vector at four quadrature points
                        *o++ = 1.; *o++ = 0.; *o++ = 0.;
                        *o++ = 1.; *o++ = 0.; *o++ = 0.;
                        *o++ = 1.; *o++ = 0.; *o++ = 0.;
                        *o++ = 1.; *o++ = 0.; *o = 0.;
                    }
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE[2]; ++k2) {
                    for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]));
                        // set vector at four quadrature points
                        *o++ = 0.; *o++ = -1.; *o++ = 0.;
                        *o++ = 0.; *o++ = -1.; *o++ = 0.;
                        *o++ = 0.; *o++ = -1.; *o++ = 0.;
                        *o++ = 0.; *o++ = -1.; *o = 0.;
                    }
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE[2]; ++k2) {
                    for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]));
                        // set vector at four quadrature points
                        *o++ = 0.; *o++ = 1.; *o++ = 0.;
                        *o++ = 0.; *o++ = 1.; *o++ = 0.;
                        *o++ = 0.; *o++ = 1.; *o++ = 0.;
                        *o++ = 0.; *o++ = 1.; *o = 0.;
                    }
                }
            }

            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]));
                        // set vector at four quadrature points
                        *o++ = 0.; *o++ = 0.; *o++ = -1.;
                        *o++ = 0.; *o++ = 0.; *o++ = -1.;
                        *o++ = 0.; *o++ = 0.; *o++ = -1.;
                        *o++ = 0.; *o++ = 0.; *o = -1.;
                    }
                }
            }

            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]));
                        // set vector at four quadrature points
                        *o++ = 0.; *o++ = 0.; *o++ = 1.;
                        *o++ = 0.; *o++ = 0.; *o++ = 1.;
                        *o++ = 0.; *o++ = 0.; *o++ = 1.;
                        *o++ = 0.; *o++ = 0.; *o = 1.;
                    }
                }
            }
        } // end of parallel section
    } else if (out.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        out.requireWrite();
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE[2]; ++k2) {
                    for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                        double* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]));
                        *o++ = -1.;
                        *o++ = 0.;
                        *o = 0.;
                    }
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE[2]; ++k2) {
                    for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                        double* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]));
                        *o++ = 1.;
                        *o++ = 0.;
                        *o = 0.;
                    }
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE[2]; ++k2) {
                    for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]));
                        *o++ = 0.;
                        *o++ = -1.;
                        *o = 0.;
                    }
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE[2]; ++k2) {
                    for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]));
                        *o++ = 0.;
                        *o++ = 1.;
                        *o = 0.;
                    }
                }
            }

            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]));
                        *o++ = 0.;
                        *o++ = 0.;
                        *o = -1.;
                    }
                }
            }

            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]));
                        *o++ = 0.;
                        *o++ = 0.;
                        *o = 1.;
                    }
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

void Brick::setToSize(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == Elements
            || out.getFunctionSpace().getTypeCode() == ReducedElements) {
        out.requireWrite();
        const dim_t numQuad=out.getNumDataPointsPerSample();
        const double size=sqrt(m_dx[0]*m_dx[0]+m_dx[1]*m_dx[1]+m_dx[2]*m_dx[2]);
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
                const double size=min(m_dx[1],m_dx[2]);
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE[2]; ++k2) {
                    for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                        double* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]));
                        fill(o, o+numQuad, size);
                    }
                }
            }

            if (m_faceOffset[1] > -1) {
                const double size=min(m_dx[1],m_dx[2]);
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE[2]; ++k2) {
                    for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                        double* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]));
                        fill(o, o+numQuad, size);
                    }
                }
            }

            if (m_faceOffset[2] > -1) {
                const double size=min(m_dx[0],m_dx[2]);
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE[2]; ++k2) {
                    for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]));
                        fill(o, o+numQuad, size);
                    }
                }
            }

            if (m_faceOffset[3] > -1) {
                const double size=min(m_dx[0],m_dx[2]);
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE[2]; ++k2) {
                    for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]));
                        fill(o, o+numQuad, size);
                    }
                }
            }

            if (m_faceOffset[4] > -1) {
                const double size=min(m_dx[0],m_dx[1]);
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]));
                        fill(o, o+numQuad, size);
                    }
                }
            }

            if (m_faceOffset[5] > -1) {
                const double size=min(m_dx[0],m_dx[1]);
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0 = 0; k0 < m_NE[0]; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]));
                        fill(o, o+numQuad, size);
                    }
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

void Brick::Print_Mesh_Info(const bool full) const
{
    RipleyDomain::Print_Mesh_Info(full);
    if (full) {
        cout << "     Id  Coordinates" << endl;
        cout.precision(15);
        cout.setf(ios::scientific, ios::floatfield);
        for (index_t i=0; i < getNumNodes(); i++) {
            cout << "  " << setw(5) << m_nodeId[i]
                << "  " << getLocalCoordinate(i%m_NN[0], 0)
                << "  " << getLocalCoordinate(i%(m_NN[0]*m_NN[1])/m_NN[0], 1)
                << "  " << getLocalCoordinate(i/(m_NN[0]*m_NN[1]), 2) << endl;
        }
    }
}


//protected
void Brick::assembleCoordinates(escript::Data& arg) const
{
    escriptDataC x = arg.getDataC();
    int numDim = m_numDim;
    if (!isDataPointShapeEqual(&x, 1, &numDim))
        throw RipleyException("setToX: Invalid Data object shape");
    if (!numSamplesEqual(&x, 1, getNumNodes()))
        throw RipleyException("setToX: Illegal number of samples in Data object");

    arg.requireWrite();
#pragma omp parallel for
    for (dim_t i2 = 0; i2 < m_NN[2]; i2++) {
        for (dim_t i1 = 0; i1 < m_NN[1]; i1++) {
            for (dim_t i0 = 0; i0 < m_NN[0]; i0++) {
                double* point = arg.getSampleDataRW(i0+m_NN[0]*i1+m_NN[0]*m_NN[1]*i2);
                point[0] = getLocalCoordinate(i0, 0);
                point[1] = getLocalCoordinate(i1, 1);
                point[2] = getLocalCoordinate(i2, 2);
            }
        }
    }
}

//protected
void Brick::assembleGradient(escript::Data& out, const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    const double C0 = .044658198738520451079;
    const double C1 = .16666666666666666667;
    const double C2 = .21132486540518711775;
    const double C3 = .25;
    const double C4 = .5;
    const double C5 = .62200846792814621559;
    const double C6 = .78867513459481288225;

    if (out.getFunctionSpace().getTypeCode() == Elements) {
        out.requireWrite();
#pragma omp parallel
        {
            vector<double> f_000(numComp);
            vector<double> f_001(numComp);
            vector<double> f_010(numComp);
            vector<double> f_011(numComp);
            vector<double> f_100(numComp);
            vector<double> f_101(numComp);
            vector<double> f_110(numComp);
            vector<double> f_111(numComp);
#pragma omp for
            for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(INDEX3(k0,k1,k2,m_NE[0],m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            const double V0=((f_100[i]-f_000[i])*C5 + (f_111[i]-f_011[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / m_dx[0];
                            const double V1=((f_110[i]-f_010[i])*C5 + (f_101[i]-f_001[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / m_dx[0];
                            const double V2=((f_101[i]-f_001[i])*C5 + (f_110[i]-f_010[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / m_dx[0];
                            const double V3=((f_111[i]-f_011[i])*C5 + (f_100[i]-f_000[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / m_dx[0];
                            const double V4=((f_010[i]-f_000[i])*C5 + (f_111[i]-f_101[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / m_dx[1];
                            const double V5=((f_110[i]-f_100[i])*C5 + (f_011[i]-f_001[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / m_dx[1];
                            const double V6=((f_011[i]-f_001[i])*C5 + (f_110[i]-f_100[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / m_dx[1];
                            const double V7=((f_111[i]-f_101[i])*C5 + (f_010[i]-f_000[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / m_dx[1];
                            const double V8=((f_001[i]-f_000[i])*C5 + (f_111[i]-f_110[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / m_dx[2];
                            const double V9=((f_101[i]-f_100[i])*C5 + (f_011[i]-f_010[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / m_dx[2];
                            const double V10=((f_011[i]-f_010[i])*C5 + (f_101[i]-f_100[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / m_dx[2];
                            const double V11=((f_111[i]-f_110[i])*C5 + (f_001[i]-f_000[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / m_dx[2];
                            o[INDEX3(i,0,0,numComp,3)] = V0;
                            o[INDEX3(i,1,0,numComp,3)] = V4;
                            o[INDEX3(i,2,0,numComp,3)] = V8;
                            o[INDEX3(i,0,1,numComp,3)] = V0;
                            o[INDEX3(i,1,1,numComp,3)] = V5;
                            o[INDEX3(i,2,1,numComp,3)] = V9;
                            o[INDEX3(i,0,2,numComp,3)] = V1;
                            o[INDEX3(i,1,2,numComp,3)] = V4;
                            o[INDEX3(i,2,2,numComp,3)] = V10;
                            o[INDEX3(i,0,3,numComp,3)] = V1;
                            o[INDEX3(i,1,3,numComp,3)] = V5;
                            o[INDEX3(i,2,3,numComp,3)] = V11;
                            o[INDEX3(i,0,4,numComp,3)] = V2;
                            o[INDEX3(i,1,4,numComp,3)] = V6;
                            o[INDEX3(i,2,4,numComp,3)] = V8;
                            o[INDEX3(i,0,5,numComp,3)] = V2;
                            o[INDEX3(i,1,5,numComp,3)] = V7;
                            o[INDEX3(i,2,5,numComp,3)] = V9;
                            o[INDEX3(i,0,6,numComp,3)] = V3;
                            o[INDEX3(i,1,6,numComp,3)] = V6;
                            o[INDEX3(i,2,6,numComp,3)] = V10;
                            o[INDEX3(i,0,7,numComp,3)] = V3;
                            o[INDEX3(i,1,7,numComp,3)] = V7;
                            o[INDEX3(i,2,7,numComp,3)] = V11;
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of k2 loop
        } // end of parallel section
    } else if (out.getFunctionSpace().getTypeCode() == ReducedElements) {
        out.requireWrite();
#pragma omp parallel
        {
            vector<double> f_000(numComp);
            vector<double> f_001(numComp);
            vector<double> f_010(numComp);
            vector<double> f_011(numComp);
            vector<double> f_100(numComp);
            vector<double> f_101(numComp);
            vector<double> f_110(numComp);
            vector<double> f_111(numComp);
#pragma omp for
            for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(INDEX3(k0,k1,k2,m_NE[0],m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_101[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_010[i]-f_011[i])*C3 / m_dx[0];
                            o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_011[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_100[i]-f_101[i])*C3 / m_dx[1];
                            o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_011[i]+f_101[i]+f_111[i]-f_000[i]-f_010[i]-f_100[i]-f_110[i])*C3 / m_dx[2];
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of k2 loop
        } // end of parallel section
    } else if (out.getFunctionSpace().getTypeCode() == FaceElements) {
        out.requireWrite();
#pragma omp parallel
        {
            vector<double> f_000(numComp);
            vector<double> f_001(numComp);
            vector<double> f_010(numComp);
            vector<double> f_011(numComp);
            vector<double> f_100(numComp);
            vector<double> f_101(numComp);
            vector<double> f_110(numComp);
            vector<double> f_111(numComp);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(0,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(0,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(0,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(0,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(1,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(1,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(1,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(1,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            const double V0=((f_010[i]-f_000[i])*C6 + (f_011[i]-f_001[i])*C2) / m_dx[1];
                            const double V1=((f_010[i]-f_000[i])*C2 + (f_011[i]-f_001[i])*C6) / m_dx[1];
                            const double V2=((f_001[i]-f_000[i])*C6 + (f_010[i]-f_011[i])*C2) / m_dx[2];
                            const double V3=((f_001[i]-f_000[i])*C2 + (f_011[i]-f_010[i])*C6) / m_dx[2];
                            o[INDEX3(i,0,0,numComp,3)] = ((f_100[i]-f_000[i])*C5 + (f_111[i]-f_011[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / m_dx[0];
                            o[INDEX3(i,1,0,numComp,3)] = V0;
                            o[INDEX3(i,2,0,numComp,3)] = V2;
                            o[INDEX3(i,0,1,numComp,3)] = ((f_110[i]-f_010[i])*C5 + (f_101[i]-f_001[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / m_dx[0];
                            o[INDEX3(i,1,1,numComp,3)] = V0;
                            o[INDEX3(i,2,1,numComp,3)] = V3;
                            o[INDEX3(i,0,2,numComp,3)] = ((f_101[i]-f_001[i])*C5 + (f_110[i]-f_010[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / m_dx[0];
                            o[INDEX3(i,1,2,numComp,3)] = V1;
                            o[INDEX3(i,2,2,numComp,3)] = V2;
                            o[INDEX3(i,0,3,numComp,3)] = ((f_111[i]-f_011[i])*C5 + (f_100[i]-f_000[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / m_dx[0];
                            o[INDEX3(i,1,3,numComp,3)] = V1;
                            o[INDEX3(i,2,3,numComp,3)] = V3;
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 0
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            const double V0=((f_110[i]-f_100[i])*C6 + (f_111[i]-f_101[i])*C2) / m_dx[1];
                            const double V1=((f_110[i]-f_100[i])*C2 + (f_111[i]-f_101[i])*C6) / m_dx[1];
                            const double V2=((f_101[i]-f_100[i])*C6 + (f_111[i]-f_110[i])*C2) / m_dx[2];
                            const double V3=((f_101[i]-f_100[i])*C2 + (f_111[i]-f_110[i])*C6) / m_dx[2];
                            o[INDEX3(i,0,0,numComp,3)] = ((f_100[i]-f_000[i])*C5 + (f_111[i]-f_011[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / m_dx[0];
                            o[INDEX3(i,1,0,numComp,3)] = V0;
                            o[INDEX3(i,2,0,numComp,3)] = V2;
                            o[INDEX3(i,0,1,numComp,3)] = ((f_110[i]-f_010[i])*C5 + (f_101[i]-f_001[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / m_dx[0];
                            o[INDEX3(i,1,1,numComp,3)] = V0;
                            o[INDEX3(i,2,1,numComp,3)] = V3;
                            o[INDEX3(i,0,2,numComp,3)] = ((f_101[i]-f_001[i])*C5 + (f_110[i]-f_010[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / m_dx[0];
                            o[INDEX3(i,1,2,numComp,3)] = V1;
                            o[INDEX3(i,2,2,numComp,3)] = V2;
                            o[INDEX3(i,0,3,numComp,3)] = ((f_111[i]-f_011[i])*C5 + (f_100[i]-f_000[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / m_dx[0];
                            o[INDEX3(i,1,3,numComp,3)] = V1;
                            o[INDEX3(i,2,3,numComp,3)] = V3;
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 1
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,0,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,0,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,0,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,0,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            const double V0=((f_100[i]-f_000[i])*C6 + (f_101[i]-f_001[i])*C2) / m_dx[0];
                            const double V1=((f_001[i]-f_000[i])*C6 + (f_101[i]-f_100[i])*C2) / m_dx[2];
                            const double V2=((f_001[i]-f_000[i])*C2 + (f_101[i]-f_100[i])*C6) / m_dx[2];
                            o[INDEX3(i,0,0,numComp,3)] = V0;
                            o[INDEX3(i,1,0,numComp,3)] = ((f_010[i]-f_000[i])*C5 + (f_111[i]-f_101[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / m_dx[1];
                            o[INDEX3(i,2,0,numComp,3)] = V1;
                            o[INDEX3(i,0,1,numComp,3)] = V0;
                            o[INDEX3(i,1,1,numComp,3)] = ((f_110[i]-f_100[i])*C5 + (f_011[i]-f_001[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / m_dx[1];
                            o[INDEX3(i,2,1,numComp,3)] = V2;
                            o[INDEX3(i,0,2,numComp,3)] = V0;
                            o[INDEX3(i,1,2,numComp,3)] = ((f_011[i]-f_001[i])*C5 + (f_110[i]-f_100[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / m_dx[1];
                            o[INDEX3(i,2,2,numComp,3)] = V1;
                            o[INDEX3(i,0,3,numComp,3)] = V0;
                            o[INDEX3(i,1,3,numComp,3)] = ((f_111[i]-f_101[i])*C5 + (f_010[i]-f_000[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / m_dx[1];
                            o[INDEX3(i,2,3,numComp,3)] = V2;
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 2
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-2,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-2,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-2,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-2,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            const double V0=((f_110[i]-f_010[i])*C6 + (f_111[i]-f_011[i])*C2) / m_dx[0];
                            const double V1=((f_110[i]-f_010[i])*C2 + (f_111[i]-f_011[i])*C6) / m_dx[0];
                            const double V2=((f_011[i]-f_010[i])*C6 + (f_111[i]-f_110[i])*C2) / m_dx[2];
                            const double V3=((f_011[i]-f_010[i])*C2 + (f_111[i]-f_110[i])*C6) / m_dx[2];
                            o[INDEX3(i,0,0,numComp,3)] = V0;
                            o[INDEX3(i,1,0,numComp,3)] = ((f_010[i]-f_000[i])*C5 + (f_111[i]-f_101[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / m_dx[1];
                            o[INDEX3(i,2,0,numComp,3)] = V2;
                            o[INDEX3(i,0,1,numComp,3)] = V0;
                            o[INDEX3(i,1,1,numComp,3)] = ((f_110[i]-f_100[i])*C5 + (f_011[i]-f_001[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / m_dx[1];
                            o[INDEX3(i,2,1,numComp,3)] = V3;
                            o[INDEX3(i,0,2,numComp,3)] = V1;
                            o[INDEX3(i,1,2,numComp,3)] = ((f_011[i]-f_001[i])*C5 + (f_110[i]-f_100[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / m_dx[1];
                            o[INDEX3(i,2,2,numComp,3)] = V2;
                            o[INDEX3(i,0,3,numComp,3)] = V1;
                            o[INDEX3(i,1,3,numComp,3)] = ((f_111[i]-f_101[i])*C5 + (f_010[i]-f_000[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / m_dx[1];
                            o[INDEX3(i,2,3,numComp,3)] = V3;
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 3
            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,0, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,0, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,0, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,0, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            const double V0=((f_100[i]-f_000[i])*C6 + (f_110[i]-f_010[i])*C2) / m_dx[0];
                            const double V1=((f_100[i]-f_000[i])*C2 + (f_110[i]-f_010[i])*C6) / m_dx[0];
                            const double V2=((f_010[i]-f_000[i])*C6 + (f_110[i]-f_100[i])*C2) / m_dx[1];
                            const double V3=((f_010[i]-f_000[i])*C2 + (f_110[i]-f_100[i])*C6) / m_dx[1];
                            o[INDEX3(i,0,0,numComp,3)] = V0;
                            o[INDEX3(i,1,0,numComp,3)] = V2;
                            o[INDEX3(i,2,0,numComp,3)] = ((f_001[i]-f_000[i])*C5 + (f_111[i]-f_110[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / m_dx[2];
                            o[INDEX3(i,0,1,numComp,3)] = V0;
                            o[INDEX3(i,1,1,numComp,3)] = V3;
                            o[INDEX3(i,2,1,numComp,3)] = ((f_101[i]-f_100[i])*C5 + (f_011[i]-f_010[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / m_dx[2];
                            o[INDEX3(i,0,2,numComp,3)] = V1;
                            o[INDEX3(i,1,2,numComp,3)] = V2;
                            o[INDEX3(i,2,2,numComp,3)] = ((f_011[i]-f_010[i])*C5 + (f_101[i]-f_100[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / m_dx[2];
                            o[INDEX3(i,0,3,numComp,3)] = V1;
                            o[INDEX3(i,1,3,numComp,3)] = V3;
                            o[INDEX3(i,2,3,numComp,3)] = ((f_111[i]-f_110[i])*C5 + (f_001[i]-f_000[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / m_dx[2];
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of face 4
            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            const double V0=((f_101[i]-f_001[i])*C6 + (f_111[i]-f_011[i])*C2) / m_dx[0];
                            const double V1=((f_101[i]-f_001[i])*C2 + (f_111[i]-f_011[i])*C6) / m_dx[0];
                            const double V2=((f_011[i]-f_001[i])*C6 + (f_111[i]-f_101[i])*C2) / m_dx[1];
                            const double V3=((f_011[i]-f_001[i])*C2 + (f_111[i]-f_101[i])*C6) / m_dx[1];
                            o[INDEX3(i,0,0,numComp,3)] = V0;
                            o[INDEX3(i,1,0,numComp,3)] = V2;
                            o[INDEX3(i,2,0,numComp,3)] = ((f_001[i]-f_000[i])*C5 + (f_111[i]-f_110[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / m_dx[2];
                            o[INDEX3(i,0,1,numComp,3)] = V0;
                            o[INDEX3(i,1,1,numComp,3)] = V3;
                            o[INDEX3(i,2,1,numComp,3)] = ((f_011[i]-f_010[i])*C0 + (f_101[i]-f_100[i])*C5 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / m_dx[2];
                            o[INDEX3(i,0,2,numComp,3)] = V1;
                            o[INDEX3(i,1,2,numComp,3)] = V2;
                            o[INDEX3(i,2,2,numComp,3)] = ((f_011[i]-f_010[i])*C5 + (f_101[i]-f_100[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / m_dx[2];
                            o[INDEX3(i,0,3,numComp,3)] = V1;
                            o[INDEX3(i,1,3,numComp,3)] = V3;
                            o[INDEX3(i,2,3,numComp,3)] = ((f_001[i]-f_000[i])*C0 + (f_111[i]-f_110[i])*C5 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / m_dx[2];
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of face 5
        } // end of parallel section
    } else if (out.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        out.requireWrite();
#pragma omp parallel
        {
            vector<double> f_000(numComp);
            vector<double> f_001(numComp);
            vector<double> f_010(numComp);
            vector<double> f_011(numComp);
            vector<double> f_100(numComp);
            vector<double> f_101(numComp);
            vector<double> f_110(numComp);
            vector<double> f_111(numComp);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(0,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(0,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(0,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(0,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(1,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(1,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(1,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(1,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_101[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_010[i]-f_011[i])*C3 / m_dx[0];
                            o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_011[i]-f_000[i]-f_001[i])*C4 / m_dx[1];
                            o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_011[i]-f_000[i]-f_010[i])*C4 / m_dx[2];
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 0
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_101[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_010[i]-f_011[i])*C3 / m_dx[0];
                            o[INDEX3(i,1,0,numComp,3)] = (f_110[i]+f_111[i]-f_100[i]-f_101[i])*C4 / m_dx[1];
                            o[INDEX3(i,2,0,numComp,3)] = (f_101[i]+f_111[i]-f_100[i]-f_110[i])*C4 / m_dx[2];
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 1
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,0,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,0,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,0,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,0,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_101[i]-f_000[i]-f_001[i])*C4 / m_dx[0];
                            o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_011[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_100[i]-f_101[i])*C3 / m_dx[1];
                            o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_101[i]-f_000[i]-f_100[i])*C4 / m_dx[2];
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 2
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-2,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-2,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-2,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-2,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,3)] = (f_110[i]+f_111[i]-f_010[i]-f_011[i])*C4 / m_dx[0];
                            o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_011[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_100[i]-f_101[i])*C3 / m_dx[1];
                            o[INDEX3(i,2,0,numComp,3)] = (f_011[i]+f_111[i]-f_010[i]-f_110[i])*C4 / m_dx[2];
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 3
            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,0, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,0, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,0, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,0, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_110[i]-f_000[i]-f_010[i])*C4 / m_dx[0];
                            o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_110[i]-f_000[i]-f_100[i])*C4 / m_dx[1];
                            o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_011[i]+f_101[i]+f_111[i]-f_000[i]-f_010[i]-f_100[i]-f_110[i])*C4 / m_dx[2];
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of face 4
            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,3)] = (f_101[i]+f_111[i]-f_001[i]-f_011[i])*C4 / m_dx[0];
                            o[INDEX3(i,1,0,numComp,3)] = (f_011[i]+f_111[i]-f_001[i]-f_101[i])*C4 / m_dx[1];
                            o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_011[i]+f_101[i]+f_111[i]-f_000[i]-f_010[i]-f_100[i]-f_110[i])*C3 / m_dx[2];
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of face 5
        } // end of parallel section
    }
}

//protected
void Brick::assembleIntegrate(vector<double>& integrals, const escript::Data& arg) const
{
    const dim_t numComp = arg.getDataPointSize();
    const index_t left = (m_offset[0]==0 ? 0 : 1);
    const index_t bottom = (m_offset[1]==0 ? 0 : 1);
    const index_t front = (m_offset[2]==0 ? 0 : 1);
    const int fs = arg.getFunctionSpace().getTypeCode();
    if (fs == Elements && arg.actsExpanded()) {
        const double w_0 = m_dx[0]*m_dx[1]*m_dx[2]/8.;
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
#pragma omp for nowait
            for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const double* f = arg.getSampleDataRO(INDEX3(k0, k1, k2, m_NE[0], m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            const double f_0 = f[INDEX2(i,0,numComp)];
                            const double f_1 = f[INDEX2(i,1,numComp)];
                            const double f_2 = f[INDEX2(i,2,numComp)];
                            const double f_3 = f[INDEX2(i,3,numComp)];
                            const double f_4 = f[INDEX2(i,4,numComp)];
                            const double f_5 = f[INDEX2(i,5,numComp)];
                            const double f_6 = f[INDEX2(i,6,numComp)];
                            const double f_7 = f[INDEX2(i,7,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3+f_4+f_5+f_6+f_7)*w_0;
                        }  // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of k2 loop

#pragma omp critical
            for (index_t i=0; i<numComp; i++)
                integrals[i]+=int_local[i];
        } // end of parallel section

    } else if (fs==ReducedElements || (fs==Elements && !arg.actsExpanded())) {
        const double w_0 = m_dx[0]*m_dx[1]*m_dx[2];
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
#pragma omp for nowait
            for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const double* f = arg.getSampleDataRO(INDEX3(k0, k1, k2, m_NE[0], m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_0;
                        }  // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of k2 loop

#pragma omp critical
            for (index_t i=0; i<numComp; i++)
                integrals[i]+=int_local[i];
        } // end of parallel section

    } else if (fs == FaceElements && arg.actsExpanded()) {
        const double w_0 = m_dx[1]*m_dx[2]/4.;
        const double w_1 = m_dx[0]*m_dx[2]/4.;
        const double w_2 = m_dx[0]*m_dx[1]/4.;
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                    for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            const double f_0 = f[INDEX2(i,0,numComp)];
                            const double f_1 = f[INDEX2(i,1,numComp)];
                            const double f_2 = f[INDEX2(i,2,numComp)];
                            const double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_0;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                    for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            const double f_0 = f[INDEX2(i,0,numComp)];
                            const double f_1 = f[INDEX2(i,1,numComp)];
                            const double f_2 = f[INDEX2(i,2,numComp)];
                            const double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_0;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            const double f_0 = f[INDEX2(i,0,numComp)];
                            const double f_1 = f[INDEX2(i,1,numComp)];
                            const double f_2 = f[INDEX2(i,2,numComp)];
                            const double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_1;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            const double f_0 = f[INDEX2(i,0,numComp)];
                            const double f_1 = f[INDEX2(i,1,numComp)];
                            const double f_2 = f[INDEX2(i,2,numComp)];
                            const double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_1;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            const double f_0 = f[INDEX2(i,0,numComp)];
                            const double f_1 = f[INDEX2(i,1,numComp)];
                            const double f_2 = f[INDEX2(i,2,numComp)];
                            const double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_2;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            const double f_0 = f[INDEX2(i,0,numComp)];
                            const double f_1 = f[INDEX2(i,1,numComp)];
                            const double f_2 = f[INDEX2(i,2,numComp)];
                            const double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_2;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

#pragma omp critical
            for (index_t i=0; i<numComp; i++)
                integrals[i]+=int_local[i];
        } // end of parallel section

    } else if (fs==ReducedFaceElements || (fs==FaceElements && !arg.actsExpanded())) {
        const double w_0 = m_dx[1]*m_dx[2];
        const double w_1 = m_dx[0]*m_dx[2];
        const double w_2 = m_dx[0]*m_dx[1];
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                    for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_0;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                    for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_0;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_1;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_1;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_2;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_2;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

#pragma omp critical
            for (index_t i=0; i<numComp; i++)
                integrals[i]+=int_local[i];
        } // end of parallel section
    } // function space selector
}

//protected
dim_t Brick::insertNeighbourNodes(IndexVector& index, index_t node) const
{
    const dim_t nDOF0 = (m_gNE[0]+1)/m_NX[0];
    const dim_t nDOF1 = (m_gNE[1]+1)/m_NX[1];
    const dim_t nDOF2 = (m_gNE[2]+1)/m_NX[2];
    const int x=node%nDOF0;
    const int y=node%(nDOF0*nDOF1)/nDOF0;
    const int z=node/(nDOF0*nDOF1);
    int num=0;
    // loop through potential neighbours and add to index if positions are
    // within bounds
    for (int i2=-1; i2<2; i2++) {
        for (int i1=-1; i1<2; i1++) {
            for (int i0=-1; i0<2; i0++) {
                // skip node itself
                if (i0==0 && i1==0 && i2==0)
                    continue;
                // location of neighbour node
                const int nx=x+i0;
                const int ny=y+i1;
                const int nz=z+i2;
                if (nx>=0 && ny>=0 && nz>=0
                        && nx<nDOF0 && ny<nDOF1 && nz<nDOF2) {
                    index.push_back(nz*nDOF0*nDOF1+ny*nDOF0+nx);
                    num++;
                }
            }
        }
    }

    return num;
}

//protected
void Brick::nodesToDOF(escript::Data& out, const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    out.requireWrite();

    const index_t left = (m_offset[0]==0 ? 0 : 1);
    const index_t bottom = (m_offset[1]==0 ? 0 : 1);
    const index_t front = (m_offset[2]==0 ? 0 : 1);
    const dim_t nDOF0 = (m_gNE[0]+1)/m_NX[0];
    const dim_t nDOF1 = (m_gNE[1]+1)/m_NX[1];
    const dim_t nDOF2 = (m_gNE[2]+1)/m_NX[2];
#pragma omp parallel for
    for (index_t i=0; i<nDOF2; i++) {
        for (index_t j=0; j<nDOF1; j++) {
            for (index_t k=0; k<nDOF0; k++) {
                const index_t n=k+left+(j+bottom)*m_NN[0]+(i+front)*m_NN[0]*m_NN[1];
                const double* src=in.getSampleDataRO(n);
                copy(src, src+numComp, out.getSampleDataRW(k+j*nDOF0+i*nDOF0*nDOF1));
            }
        }
    }
}

//protected
void Brick::dofToNodes(escript::Data& out, const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    Paso_Coupler* coupler = Paso_Coupler_alloc(m_connector, numComp);
    // expand data object if necessary to be able to grab the whole data
    const_cast<escript::Data*>(&in)->expand();
    Paso_Coupler_startCollect(coupler, in.getSampleDataRO(0));

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
void Brick::populateSampleIds()
{
    // degrees of freedom are numbered from left to right, bottom to top, front
    // to back in each rank, continuing on the next rank (ranks also go
    // left-right, bottom-top, front-back).
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
        m_faceCount[0]=m_NE[1]*m_NE[2];
    else
        m_faceCount[0]=0;
    //right
    if (m_mpiInfo->rank%m_NX[0]==m_NX[0]-1)
        m_faceCount[1]=m_NE[1]*m_NE[2];
    else
        m_faceCount[1]=0;
    //bottom
    if (m_offset[1]==0)
        m_faceCount[2]=m_NE[0]*m_NE[2];
    else
        m_faceCount[2]=0;
    //top
    if (m_mpiInfo->rank%(m_NX[0]*m_NX[1])/m_NX[0]==m_NX[1]-1)
        m_faceCount[3]=m_NE[0]*m_NE[2];
    else
        m_faceCount[3]=0;
    //front
    if (m_offset[2]==0)
        m_faceCount[4]=m_NE[0]*m_NE[1];
    else
        m_faceCount[4]=0;
    //back
    if (m_mpiInfo->rank/(m_NX[0]*m_NX[1])==m_NX[2]-1)
        m_faceCount[5]=m_NE[0]*m_NE[1];
    else
        m_faceCount[5]=0;

    m_faceId.resize(getNumFaceElements());

    const index_t left = (m_offset[0]==0 ? 0 : 1);
    const index_t bottom = (m_offset[1]==0 ? 0 : 1);
    const index_t front = (m_offset[2]==0 ? 0 : 1);
    const dim_t nDOF0 = (m_gNE[0]+1)/m_NX[0];
    const dim_t nDOF1 = (m_gNE[1]+1)/m_NX[1];
    const dim_t nDOF2 = (m_gNE[2]+1)/m_NX[2];

    // the following is a compromise between efficiency and code length to
    // set the node id's according to the order mentioned above.
    // First we set all the edge and corner id's in a rather slow way since
    // they might or might not be owned by this rank. Next come the own
    // node id's which are identical to the DOF id's (simple loop), and finally
    // the 6 faces are set but only if required...

#define globalNodeId(x,y,z) \
    ((m_offset[0]+x)/nDOF0)*nDOF0*nDOF1*nDOF2+(m_offset[0]+x)%nDOF0\
    + ((m_offset[1]+y)/nDOF1)*nDOF0*nDOF1*nDOF2*m_NX[0]+((m_offset[1]+y)%nDOF1)*nDOF0\
    + ((m_offset[2]+z)/nDOF2)*nDOF0*nDOF1*nDOF2*m_NX[0]*m_NX[1]+((m_offset[2]+z)%nDOF2)*nDOF0*nDOF1

#pragma omp parallel
    {
        // set edge id's
        // edges in x-direction, including corners
#pragma omp for nowait
        for (dim_t i=0; i<m_NN[0]; i++) {
            m_nodeId[i] = globalNodeId(i, 0, 0); // LF
            m_nodeId[m_NN[0]*(m_NN[1]-1)+i] = globalNodeId(i, m_NN[1]-1, 0); // UF
            m_nodeId[m_NN[0]*m_NN[1]*(m_NN[2]-1)+i] = globalNodeId(i, 0, m_NN[2]-1); // LB
            m_nodeId[m_NN[0]*m_NN[1]*m_NN[2]-m_NN[0]+i] = globalNodeId(i, m_NN[1]-1, m_NN[2]-1); // UB
        }
        // edges in y-direction, without corners
#pragma omp for nowait
        for (dim_t i=1; i<m_NN[1]-1; i++) {
            m_nodeId[m_NN[0]*i] = globalNodeId(0, i, 0); // FL
            m_nodeId[m_NN[0]*(i+1)-1] = globalNodeId(m_NN[0]-1, i, 0); // FR
            m_nodeId[m_NN[0]*m_NN[1]*(m_NN[2]-1)+m_NN[0]*i] = globalNodeId(0, i, m_NN[2]-1); // BL
            m_nodeId[m_NN[0]*m_NN[1]*(m_NN[2]-1)+m_NN[0]*(i+1)-1] = globalNodeId(m_NN[0]-1, i, m_NN[2]-1); // BR
        }
        // edges in z-direction, without corners
#pragma omp for
        for (dim_t i=1; i<m_NN[2]-1; i++) {
            m_nodeId[m_NN[0]*m_NN[1]*i] = globalNodeId(0, 0, i); // LL
            m_nodeId[m_NN[0]*m_NN[1]*i+m_NN[0]-1] = globalNodeId(m_NN[0]-1, 0, i); // LR
            m_nodeId[m_NN[0]*m_NN[1]*(i+1)-m_NN[0]] = globalNodeId(0, m_NN[1]-1, i); // UL
            m_nodeId[m_NN[0]*m_NN[1]*(i+1)-1] = globalNodeId(m_NN[0]-1, m_NN[1]-1, i); // UR
        }
        // implicit barrier here because some node IDs will be overwritten
        // below

        // populate degrees of freedom and own nodes (identical id)
#pragma omp for nowait
        for (dim_t i=0; i<nDOF2; i++) {
            for (dim_t j=0; j<nDOF1; j++) {
                for (dim_t k=0; k<nDOF0; k++) {
                    const index_t nodeIdx=k+left+(j+bottom)*m_NN[0]+(i+front)*m_NN[0]*m_NN[1];
                    const index_t dofIdx=k+j*nDOF0+i*nDOF0*nDOF1;
                    m_dofId[dofIdx] = m_nodeId[nodeIdx]
                        = m_nodeDistribution[m_mpiInfo->rank]+dofIdx;
                }
            }
        }

        // populate the rest of the nodes (shared with other ranks)
        if (m_faceCount[0]==0) { // left plane
#pragma omp for nowait
            for (dim_t i=0; i<nDOF2; i++) {
                for (dim_t j=0; j<nDOF1; j++) {
                    const index_t nodeIdx=(j+bottom)*m_NN[0]+(i+front)*m_NN[0]*m_NN[1];
                    const index_t dofId=(j+1)*nDOF0-1+i*nDOF0*nDOF1;
                    m_nodeId[nodeIdx]
                        = m_nodeDistribution[m_mpiInfo->rank-1]+dofId;
                }
            }
        }
        if (m_faceCount[1]==0) { // right plane
#pragma omp for nowait
            for (dim_t i=0; i<nDOF2; i++) {
                for (dim_t j=0; j<nDOF1; j++) {
                    const index_t nodeIdx=(j+bottom+1)*m_NN[0]-1+(i+front)*m_NN[0]*m_NN[1];
                    const index_t dofId=j*nDOF0+i*nDOF0*nDOF1;
                    m_nodeId[nodeIdx]
                        = m_nodeDistribution[m_mpiInfo->rank+1]+dofId;
                }
            }
        }
        if (m_faceCount[2]==0) { // bottom plane
#pragma omp for nowait
            for (dim_t i=0; i<nDOF2; i++) {
                for (dim_t k=0; k<nDOF0; k++) {
                    const index_t nodeIdx=k+left+(i+front)*m_NN[0]*m_NN[1];
                    const index_t dofId=nDOF0*(nDOF1-1)+k+i*nDOF0*nDOF1;
                    m_nodeId[nodeIdx]
                        = m_nodeDistribution[m_mpiInfo->rank-m_NX[0]]+dofId;
                }
            }
        }
        if (m_faceCount[3]==0) { // top plane
#pragma omp for nowait
            for (dim_t i=0; i<nDOF2; i++) {
                for (dim_t k=0; k<nDOF0; k++) {
                    const index_t nodeIdx=k+left+(i+front)*m_NN[0]*m_NN[1]+m_NN[0]*(m_NN[1]-1);
                    const index_t dofId=k+i*nDOF0*nDOF1;
                    m_nodeId[nodeIdx]
                        = m_nodeDistribution[m_mpiInfo->rank+m_NX[0]]+dofId;
                }
            }
        }
        if (m_faceCount[4]==0) { // front plane
#pragma omp for nowait
            for (dim_t j=0; j<nDOF1; j++) {
                for (dim_t k=0; k<nDOF0; k++) {
                    const index_t nodeIdx=k+left+(j+bottom)*m_NN[0];
                    const index_t dofId=k+j*nDOF0+nDOF0*nDOF1*(nDOF2-1);
                    m_nodeId[nodeIdx]
                        = m_nodeDistribution[m_mpiInfo->rank-m_NX[0]*m_NX[1]]+dofId;
                }
            }
        }
        if (m_faceCount[5]==0) { // back plane
#pragma omp for nowait
            for (dim_t j=0; j<nDOF1; j++) {
                for (dim_t k=0; k<nDOF0; k++) {
                    const index_t nodeIdx=k+left+(j+bottom)*m_NN[0]+m_NN[0]*m_NN[1]*(m_NN[2]-1);
                    const index_t dofId=k+j*nDOF0;
                    m_nodeId[nodeIdx]
                        = m_nodeDistribution[m_mpiInfo->rank+m_NX[0]*m_NX[1]]+dofId;
                }
            }
        }

        // populate element id's
#pragma omp for nowait
        for (dim_t i2=0; i2<m_NE[2]; i2++) {
            for (dim_t i1=0; i1<m_NE[1]; i1++) {
                for (dim_t i0=0; i0<m_NE[0]; i0++) {
                    m_elementId[i0+i1*m_NE[0]+i2*m_NE[0]*m_NE[1]] =
                        (m_offset[2]+i2)*m_gNE[0]*m_gNE[1]
                        +(m_offset[1]+i1)*m_gNE[0]
                        +m_offset[0]+i0;
                }
            }
        }

        // face elements
#pragma omp for
        for (dim_t k=0; k<getNumFaceElements(); k++)
            m_faceId[k]=k;
    } // end parallel section

#undef globalNodeId

    m_nodeTags.assign(getNumNodes(), 0);
    updateTagsInUse(Nodes);

    m_elementTags.assign(getNumElements(), 0);
    updateTagsInUse(Elements);

    // generate face offset vector and set face tags
    const index_t LEFT=1, RIGHT=2, BOTTOM=10, TOP=20, FRONT=100, BACK=200;
    const index_t faceTag[] = { LEFT, RIGHT, BOTTOM, TOP, FRONT, BACK };
    m_faceOffset.assign(6, -1);
    m_faceTags.clear();
    index_t offset=0;
    for (size_t i=0; i<6; i++) {
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
    setTagMap("front", FRONT);
    setTagMap("back", BACK);
    updateTagsInUse(FaceElements);
}

//private
void Brick::createPattern()
{
    const dim_t nDOF0 = (m_gNE[0]+1)/m_NX[0];
    const dim_t nDOF1 = (m_gNE[1]+1)/m_NX[1];
    const dim_t nDOF2 = (m_gNE[2]+1)/m_NX[2];
    const index_t left = (m_offset[0]==0 ? 0 : 1);
    const index_t bottom = (m_offset[1]==0 ? 0 : 1);
    const index_t front = (m_offset[2]==0 ? 0 : 1);

    // populate node->DOF mapping with own degrees of freedom.
    // The rest is assigned in the loop further down
    m_dofMap.assign(getNumNodes(), 0);
#pragma omp parallel for
    for (index_t i=front; i<front+nDOF2; i++) {
        for (index_t j=bottom; j<bottom+nDOF1; j++) {
            for (index_t k=left; k<left+nDOF0; k++) {
                m_dofMap[i*m_NN[0]*m_NN[1]+j*m_NN[0]+k]=(i-front)*nDOF0*nDOF1+(j-bottom)*nDOF0+k-left;
            }
        }
    }

    // build list of shared components and neighbours by looping through
    // all potential neighbouring ranks and checking if positions are
    // within bounds
    const dim_t numDOF=nDOF0*nDOF1*nDOF2;
    vector<IndexVector> colIndices(numDOF); // for the couple blocks
    RankVector neighbour;
    IndexVector offsetInShared(1,0);
    IndexVector sendShared, recvShared;
    int numShared=0;
    const int x=m_mpiInfo->rank%m_NX[0];
    const int y=m_mpiInfo->rank%(m_NX[0]*m_NX[1])/m_NX[0];
    const int z=m_mpiInfo->rank/(m_NX[0]*m_NX[1]);
    for (int i2=-1; i2<2; i2++) {
        for (int i1=-1; i1<2; i1++) {
            for (int i0=-1; i0<2; i0++) {
                // skip this rank
                if (i0==0 && i1==0 && i2==0)
                    continue;
                // location of neighbour rank
                const int nx=x+i0;
                const int ny=y+i1;
                const int nz=z+i2;
                if (nx>=0 && ny>=0 && nz>=0 && nx<m_NX[0] && ny<m_NX[1] && nz<m_NX[2]) {
                    neighbour.push_back(nz*m_NX[0]*m_NX[1]+ny*m_NX[0]+nx);
                    if (i0==0 && i1==0) {
                        // sharing front or back plane
                        offsetInShared.push_back(offsetInShared.back()+nDOF0*nDOF1);
                        for (dim_t i=0; i<nDOF1; i++) {
                            const int firstDOF=(i2==-1 ? i*nDOF0
                                    : i*nDOF0 + nDOF0*nDOF1*(nDOF2-1));
                            const int firstNode=(i2==-1 ? left+(i+bottom)*m_NN[0]
                                    : left+(i+bottom)*m_NN[0]+m_NN[0]*m_NN[1]*(m_NN[2]-1));
                            for (dim_t j=0; j<nDOF0; j++, numShared++) {
                                sendShared.push_back(firstDOF+j);
                                recvShared.push_back(numDOF+numShared);
                                if (j>0) {
                                    if (i>0)
                                        colIndices[firstDOF+j-1-nDOF0].push_back(numShared);
                                    colIndices[firstDOF+j-1].push_back(numShared);
                                    if (i<nDOF1-1)
                                        colIndices[firstDOF+j-1+nDOF0].push_back(numShared);
                                }
                                if (i>0)
                                    colIndices[firstDOF+j-nDOF0].push_back(numShared);
                                colIndices[firstDOF+j].push_back(numShared);
                                if (i<nDOF1-1)
                                    colIndices[firstDOF+j+nDOF0].push_back(numShared);
                                if (j<nDOF0-1) {
                                    if (i>0)
                                        colIndices[firstDOF+j+1-nDOF0].push_back(numShared);
                                    colIndices[firstDOF+j+1].push_back(numShared);
                                    if (i<nDOF1-1)
                                        colIndices[firstDOF+j+1+nDOF0].push_back(numShared);
                                }
                                m_dofMap[firstNode+j]=numDOF+numShared;
                            }
                        }
                    } else if (i0==0 && i2==0) {
                        // sharing top or bottom plane
                        offsetInShared.push_back(offsetInShared.back()+nDOF0*nDOF2);
                        for (dim_t i=0; i<nDOF2; i++) {
                            const int firstDOF=(i1==-1 ? i*nDOF0*nDOF1
                                    : nDOF0*((i+1)*nDOF1-1));
                            const int firstNode=(i1==-1 ?
                                    left+(i+front)*m_NN[0]*m_NN[1]
                                    : left+m_NN[0]*((i+1+front)*m_NN[1]-1));
                            for (dim_t j=0; j<nDOF0; j++, numShared++) {
                                sendShared.push_back(firstDOF+j);
                                recvShared.push_back(numDOF+numShared);
                                if (j>0) {
                                    if (i>0)
                                        colIndices[firstDOF+j-1-nDOF0*nDOF1].push_back(numShared);
                                    colIndices[firstDOF+j-1].push_back(numShared);
                                    if (i<nDOF2-1)
                                        colIndices[firstDOF+j-1+nDOF0*nDOF1].push_back(numShared);
                                }
                                if (i>0)
                                    colIndices[firstDOF+j-nDOF0*nDOF1].push_back(numShared);
                                colIndices[firstDOF+j].push_back(numShared);
                                if (i<nDOF2-1)
                                    colIndices[firstDOF+j+nDOF0*nDOF1].push_back(numShared);
                                if (j<nDOF0-1) {
                                    if (i>0)
                                        colIndices[firstDOF+j+1-nDOF0*nDOF1].push_back(numShared);
                                    colIndices[firstDOF+j+1].push_back(numShared);
                                    if (i<nDOF2-1)
                                        colIndices[firstDOF+j+1+nDOF0*nDOF1].push_back(numShared);
                                }
                                m_dofMap[firstNode+j]=numDOF+numShared;
                            }
                        }
                    } else if (i1==0 && i2==0) {
                        // sharing left or right plane
                        offsetInShared.push_back(offsetInShared.back()+nDOF1*nDOF2);
                        for (dim_t i=0; i<nDOF2; i++) {
                            const int firstDOF=(i0==-1 ? i*nDOF0*nDOF1
                                    : nDOF0*(1+i*nDOF1)-1);
                            const int firstNode=(i0==-1 ?
                                    (bottom+(i+front)*m_NN[1])*m_NN[0]
                                    : (bottom+1+(i+front)*m_NN[1])*m_NN[0]-1);
                            for (dim_t j=0; j<nDOF1; j++, numShared++) {
                                sendShared.push_back(firstDOF+j*nDOF0);
                                recvShared.push_back(numDOF+numShared);
                                if (j>0) {
                                    if (i>0)
                                        colIndices[firstDOF+(j-1)*nDOF0-nDOF0*nDOF1].push_back(numShared);
                                    colIndices[firstDOF+(j-1)*nDOF0].push_back(numShared);
                                    if (i<nDOF2-1)
                                        colIndices[firstDOF+(j-1)*nDOF0+nDOF0*nDOF1].push_back(numShared);
                                }
                                if (i>0)
                                    colIndices[firstDOF+j*nDOF0-nDOF0*nDOF1].push_back(numShared);
                                colIndices[firstDOF+j*nDOF0].push_back(numShared);
                                if (i<nDOF2-1)
                                    colIndices[firstDOF+j*nDOF0+nDOF0*nDOF1].push_back(numShared);
                                if (j<nDOF1-1) {
                                    if (i>0)
                                        colIndices[firstDOF+(j+1)*nDOF0-nDOF0*nDOF1].push_back(numShared);
                                    colIndices[firstDOF+(j+1)*nDOF0].push_back(numShared);
                                    if (i<nDOF2-1)
                                        colIndices[firstDOF+(j+1)*nDOF0+nDOF0*nDOF1].push_back(numShared);
                                }
                                m_dofMap[firstNode+j*m_NN[0]]=numDOF+numShared;
                            }
                        }
                    } else if (i0==0) {
                        // sharing an edge in x direction
                        offsetInShared.push_back(offsetInShared.back()+nDOF0);
                        const int firstDOF=(i1+1)/2*nDOF0*(nDOF1-1)
                                           +(i2+1)/2*nDOF0*nDOF1*(nDOF2-1);
                        const int firstNode=left+(i1+1)/2*m_NN[0]*(m_NN[1]-1)
                                            +(i2+1)/2*m_NN[0]*m_NN[1]*(m_NN[2]-1);
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
                        // sharing an edge in y direction
                        offsetInShared.push_back(offsetInShared.back()+nDOF1);
                        const int firstDOF=(i0+1)/2*(nDOF0-1)
                                           +(i2+1)/2*nDOF0*nDOF1*(nDOF2-1);
                        const int firstNode=bottom*m_NN[0]
                                            +(i0+1)/2*(m_NN[0]-1)
                                            +(i2+1)/2*m_NN[0]*m_NN[1]*(m_NN[2]-1);
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
                    } else if (i2==0) {
                        // sharing an edge in z direction
                        offsetInShared.push_back(offsetInShared.back()+nDOF2);
                        const int firstDOF=(i0+1)/2*(nDOF0-1)
                                           +(i1+1)/2*nDOF0*(nDOF1-1);
                        const int firstNode=front*m_NN[0]*m_NN[1]
                                            +(i0+1)/2*(m_NN[0]-1)
                                            +(i1+1)/2*m_NN[0]*(m_NN[1]-1);
                        for (dim_t i=0; i<nDOF2; i++, numShared++) {
                            sendShared.push_back(firstDOF+i*nDOF0*nDOF1);
                            recvShared.push_back(numDOF+numShared);
                            if (i>0)
                                colIndices[firstDOF+(i-1)*nDOF0*nDOF1].push_back(numShared);
                            colIndices[firstDOF+i*nDOF0*nDOF1].push_back(numShared);
                            if (i<nDOF2-1)
                                colIndices[firstDOF+(i+1)*nDOF0*nDOF1].push_back(numShared);
                            m_dofMap[firstNode+i*m_NN[0]*m_NN[1]]=numDOF+numShared;
                        }
                    } else {
                        // sharing a node
                        const int dof=(i0+1)/2*(nDOF0-1)
                                      +(i1+1)/2*nDOF0*(nDOF1-1)
                                      +(i2+1)/2*nDOF0*nDOF1*(nDOF2-1);
                        const int node=(i0+1)/2*(m_NN[0]-1)
                                       +(i1+1)/2*m_NN[0]*(m_NN[1]-1)
                                       +(i2+1)/2*m_NN[0]*m_NN[1]*(m_NN[2]-1);
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
void Brick::addToMatrixAndRHS(Paso_SystemMatrix* S, escript::Data& F,
         const vector<double>& EM_S, const vector<double>& EM_F, bool addS,
         bool addF, index_t firstNode, dim_t nEq, dim_t nComp) const
{
    IndexVector rowIndex;
    rowIndex.push_back(m_dofMap[firstNode]);
    rowIndex.push_back(m_dofMap[firstNode+1]);
    rowIndex.push_back(m_dofMap[firstNode+m_NN[0]]);
    rowIndex.push_back(m_dofMap[firstNode+m_NN[0]+1]);
    rowIndex.push_back(m_dofMap[firstNode+m_NN[0]*m_NN[1]]);
    rowIndex.push_back(m_dofMap[firstNode+m_NN[0]*m_NN[1]+1]);
    rowIndex.push_back(m_dofMap[firstNode+m_NN[0]*(m_NN[1]+1)]);
    rowIndex.push_back(m_dofMap[firstNode+m_NN[0]*(m_NN[1]+1)+1]);
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
void Brick::interpolateNodesOnElements(escript::Data& out,
                                       const escript::Data& in,
                                       bool reduced) const
{
    const dim_t numComp = in.getDataPointSize();
    if (reduced) {
        out.requireWrite();
#pragma omp parallel
        {
            vector<double> f_000(numComp);
            vector<double> f_001(numComp);
            vector<double> f_010(numComp);
            vector<double> f_011(numComp);
            vector<double> f_100(numComp);
            vector<double> f_101(numComp);
            vector<double> f_110(numComp);
            vector<double> f_111(numComp);
#pragma omp for
            for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(INDEX3(k0,k1,k2,m_NE[0],m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = (f_000[i] + f_001[i] + f_010[i] + f_011[i] + f_100[i] + f_101[i] + f_110[i] + f_111[i])/8;
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of k2 loop
        } // end of parallel section
    } else {
        out.requireWrite();
        const double c0 = .0094373878376559314545;
        const double c1 = .035220810900864519624;
        const double c2 = .13144585576580214704;
        const double c3 = .49056261216234406855;
#pragma omp parallel
        {
            vector<double> f_000(numComp);
            vector<double> f_001(numComp);
            vector<double> f_010(numComp);
            vector<double> f_011(numComp);
            vector<double> f_100(numComp);
            vector<double> f_101(numComp);
            vector<double> f_110(numComp);
            vector<double> f_111(numComp);
#pragma omp for
            for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(INDEX3(k0,k1,k2,m_NE[0],m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = f_000[i]*c3 + f_111[i]*c0 + c2*(f_001[i] + f_010[i] + f_100[i]) + c1*(f_011[i] + f_101[i] + f_110[i]);
                            o[INDEX2(i,numComp,1)] = f_011[i]*c0 + f_100[i]*c3 + c2*(f_000[i] + f_101[i] + f_110[i]) + c1*(f_001[i] + f_010[i] + f_111[i]);
                            o[INDEX2(i,numComp,2)] = f_010[i]*c3 + f_101[i]*c0 + c2*(f_000[i] + f_011[i] + f_110[i]) + c1*(f_001[i] + f_100[i] + f_111[i]);
                            o[INDEX2(i,numComp,3)] = f_001[i]*c0 + f_110[i]*c3 + c2*(f_010[i] + f_100[i] + f_111[i]) + c1*(f_000[i] + f_011[i] + f_101[i]);
                            o[INDEX2(i,numComp,4)] = f_001[i]*c3 + f_110[i]*c0 + c2*(f_000[i] + f_011[i] + f_101[i]) + c1*(f_010[i] + f_100[i] + f_111[i]);
                            o[INDEX2(i,numComp,5)] = f_010[i]*c0 + f_101[i]*c3 + c2*(f_001[i] + f_100[i] + f_111[i]) + c1*(f_000[i] + f_011[i] + f_110[i]);
                            o[INDEX2(i,numComp,6)] = f_011[i]*c3 + f_100[i]*c0 + c2*(f_001[i] + f_010[i] + f_111[i]) + c1*(f_000[i] + f_101[i] + f_110[i]);
                            o[INDEX2(i,numComp,7)] = f_000[i]*c0 + f_111[i]*c3 + c2*(f_011[i] + f_101[i] + f_110[i]) + c1*(f_001[i] + f_010[i] + f_100[i]);
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of k2 loop
        } // end of parallel section
    }
}

//protected
void Brick::interpolateNodesOnFaces(escript::Data& out, const escript::Data& in,
                                    bool reduced) const
{
    const dim_t numComp = in.getDataPointSize();
    if (reduced) {
        out.requireWrite();
#pragma omp parallel
        {
            vector<double> f_000(numComp);
            vector<double> f_001(numComp);
            vector<double> f_010(numComp);
            vector<double> f_011(numComp);
            vector<double> f_100(numComp);
            vector<double> f_101(numComp);
            vector<double> f_110(numComp);
            vector<double> f_111(numComp);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(0,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(0,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(0,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(0,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = (f_000[i] + f_001[i] + f_010[i] + f_011[i])/4;
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 0
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = (f_100[i] + f_101[i] + f_110[i] + f_111[i])/4;
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 1
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,0,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,0,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,0,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,0,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = (f_000[i] + f_001[i] + f_100[i] + f_101[i])/4;
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 2
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = (f_010[i] + f_011[i] + f_110[i] + f_111[i])/4;
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 3
            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,0, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,0, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,0, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,0, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = (f_000[i] + f_010[i] + f_100[i] + f_110[i])/4;
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of face 4
            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = (f_001[i] + f_011[i] + f_101[i] + f_111[i])/4;
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of face 5
        } // end of parallel section
    } else {
        out.requireWrite();
        const double c0 = 0.044658198738520451079;
        const double c1 = 0.16666666666666666667;
        const double c2 = 0.62200846792814621559;
#pragma omp parallel
        {
            vector<double> f_000(numComp);
            vector<double> f_001(numComp);
            vector<double> f_010(numComp);
            vector<double> f_011(numComp);
            vector<double> f_100(numComp);
            vector<double> f_101(numComp);
            vector<double> f_110(numComp);
            vector<double> f_111(numComp);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(0,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(0,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(0,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(0,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = f_000[i]*c2 + f_011[i]*c0 + c1*(f_001[i] + f_010[i]);
                            o[INDEX2(i,numComp,1)] = f_001[i]*c0 + f_010[i]*c2 + c1*(f_000[i] + f_011[i]);
                            o[INDEX2(i,numComp,2)] = f_001[i]*c2 + f_010[i]*c0 + c1*(f_000[i] + f_011[i]);
                            o[INDEX2(i,numComp,3)] = f_000[i]*c0 + f_011[i]*c2 + c1*(f_001[i] + f_010[i]);
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 0
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = f_100[i]*c2 + f_111[i]*c0 + c1*(f_101[i] + f_110[i]);
                            o[INDEX2(i,numComp,1)] = f_101[i]*c0 + f_110[i]*c2 + c1*(f_100[i] + f_111[i]);
                            o[INDEX2(i,numComp,2)] = f_101[i]*c2 + f_110[i]*c0 + c1*(f_100[i] + f_111[i]);
                            o[INDEX2(i,numComp,3)] = f_100[i]*c0 + f_111[i]*c2 + c1*(f_101[i] + f_110[i]);
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 1
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,0,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,0,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,0,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,0,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = f_000[i]*c2 + f_101[i]*c0 + c1*(f_001[i] + f_100[i]);
                            o[INDEX2(i,numComp,1)] = f_001[i]*c0 + f_100[i]*c2 + c1*(f_000[i] + f_101[i]);
                            o[INDEX2(i,numComp,2)] = f_001[i]*c2 + f_100[i]*c0 + c1*(f_000[i] + f_101[i]);
                            o[INDEX2(i,numComp,3)] = f_000[i]*c0 + f_101[i]*c2 + c1*(f_001[i] + f_100[i]);
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 2
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2+1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = f_010[i]*c2 + f_111[i]*c0 + c1*(f_011[i] + f_110[i]);
                            o[INDEX2(i,numComp,1)] = f_011[i]*c0 + f_110[i]*c2 + c1*(f_010[i] + f_111[i]);
                            o[INDEX2(i,numComp,2)] = f_011[i]*c2 + f_110[i]*c0 + c1*(f_010[i] + f_111[i]);
                            o[INDEX2(i,numComp,3)] = f_010[i]*c0 + f_111[i]*c2 + c1*(f_011[i] + f_110[i]);
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 3
            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,0, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,0, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,0, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,0, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = f_000[i]*c2 + f_110[i]*c0 + c1*(f_010[i] + f_100[i]);
                            o[INDEX2(i,numComp,1)] = f_010[i]*c0 + f_100[i]*c2 + c1*(f_000[i] + f_110[i]);
                            o[INDEX2(i,numComp,2)] = f_010[i]*c2 + f_100[i]*c0 + c1*(f_000[i] + f_110[i]);
                            o[INDEX2(i,numComp,3)] = f_000[i]*c0 + f_110[i]*c2 + c1*(f_010[i] + f_100[i]);
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of face 4
            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-1, m_NN[0],m_NN[1])), numComp*sizeof(double));
                        double* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = f_001[i]*c2 + f_111[i]*c0 + c1*(f_011[i] + f_101[i]);
                            o[INDEX2(i,numComp,1)] = f_011[i]*c0 + f_101[i]*c2 + c1*(f_001[i] + f_111[i]);
                            o[INDEX2(i,numComp,2)] = f_011[i]*c2 + f_101[i]*c0 + c1*(f_001[i] + f_111[i]);
                            o[INDEX2(i,numComp,3)] = f_001[i]*c0 + f_111[i]*c2 + c1*(f_011[i] + f_101[i]);
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of face 5
        } // end of parallel section
    }
}

namespace
{
    // Calculates a guassian blur colvolution matrix for 3D
    // See wiki article on the subject
    double* get3DGauss(unsigned radius, double sigma)
    {
        double* arr=new double[(radius*2+1)*(radius*2+1)*(radius*2+1)];
        double common=pow(M_1_PI*0.5*1/(sigma*sigma), 3./2);
	double total=0;
	int r=static_cast<int>(radius);
	for (int z=-r;z<=r;++z)
	{
	    for (int y=-r;y<=r;++y)
	    {
		for (int x=-r;x<=r;++x)
		{	      
		    arr[(x+r)+(y+r)*(r*2+1)]=common*exp(-(x*x+y*y+z*z)/(2*sigma*sigma));
		    total+=arr[(x+r)+(y+r)*(r*2+1)];
		}
	    }
	}
	double invtotal=1/total;
	for (size_t p=0;p<(radius*2+1)*(radius*2+1);++p)
	{
	    arr[p]*=invtotal; 
	}
	return arr;
    }
    
    // applies conv to source to get a point.
    // (xp, yp) are the coords in the source matrix not the destination matrix
    double Convolve3D(double* conv, double* source, size_t xp, size_t yp, size_t zp, unsigned radius, size_t width, size_t height)
    {
        size_t bx=xp-radius, by=yp-radius, bz=zp-radius;
	size_t sbase=bx+by*width+bz*width*height;
	double result=0;
	for (int z=0;z<2*radius+1;++z)
	{
	    for (int y=0;y<2*radius+1;++y)
	    {	  
		for (int x=0;x<2*radius+1;++x)
		{
		    result+=conv[x+y*(2*radius+1)+z*(2*radius+1)*(2*radius+1)] * source[sbase + x+y*width+z*width*height];
		}
	    }
	}
        return result;      
    }
}


/* This routine produces a Data object filled with smoothed random data.
The dimensions of the rectangle being filled are internal[0] x internal[1] x internal[2] points.
A parameter radius  gives the size of the stencil used for the smoothing.
A point on the left hand edge for example, will still require `radius` extra points to the left
in order to complete the stencil.

All local calculation is done on an array called `src`, with 
dimensions = ext[0] * ext[1] *ext[2].
Where ext[i]= internal[i]+2*radius.

Now for MPI there is overlap to deal with. We need to share both the overlapping 
values themselves but also the external region.

In a hypothetical 1-D case:


1234567
would be split into two ranks thus:
123(4)  (4)567     [4 being a shared element]

If the radius is 2.   There will be padding elements on the outside:

pp123(4)  (4)567pp

To ensure that 4 can be correctly computed on both ranks, values from the other rank
need to be known.

pp123(4)56   23(4)567pp

Now in our case, we wout set all the values 23456 on the left rank and send them to the 
right hand rank.

So the edges _may_ need to be shared at a distance `inset` from all boundaries.
inset=2*radius+1    
This is to ensure that values at distance `radius` from the shared/overlapped element
that ripley has.
*/
escript::Data Brick::randomFill(long seed, const boost::python::tuple& filter) const
{
    if (m_numDim!=3)
    {
        throw RipleyException("Brick must be 3D.");
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
    double sigma=0.5;
    boost::python::extract<double> ex2(filter[2]);
    if (!ex2.check() || (sigma=ex2())<=0)
    {
        throw RipleyException("Sigma must be a postive floating point number.");
    }    
    
    size_t internal[3];
    internal[0]=m_NE[0]+1;	// number of points in the internal region
    internal[1]=m_NE[1]+1;	// that is, the ones we need smoothed versions of
    internal[2]=m_NE[2]+1;	// that is, the ones we need smoothed versions of
    size_t ext[3];
    ext[0]=internal[0]+2*radius;	// includes points we need as input
    ext[1]=internal[1]+2*radius;	// for smoothing
    ext[2]=internal[2]+2*radius;	// for smoothing
    
    // now we check to see if the radius is acceptable 
    // That is, would not cross multiple ranks in MPI

    if ((2*radius>=internal[0]) || (2*radius>=internal[1]) || (2*radius>=internal[2]))
    {
        throw RipleyException("Radius of gaussian filter must be less than half the width/height of a rank");
    }
    
   
    double* src=new double[ext[0]*ext[1]*ext[2]];
    esysUtils::randomFillArray(seed, src, ext[0]*ext[1]*ext[2]);  
    
   
#ifdef ESYS_MPI
    
    dim_t X=m_mpiInfo->rank%m_NX[0];
    dim_t Y=m_mpiInfo->rank%(m_NX[0]*m_NX[1])/m_NX[0];
    dim_t Z=m_mpiInfo->rank/(m_NX[0]*m_NX[1]);
    
    BlockGrid grid(m_NX[0]-1, m_NX[1]-1, m_NX[2]-1);
    size_t inset=2*radius+1;	
    
    size_t xmidlen=ext[0]-2*inset;	// how wide is the x-dimension between the two insets
    size_t ymidlen=ext[1]-2*inset;	
    size_t zmidlen=ext[2]-2*inset;
    
    Block block(ext[0], ext[1], ext[2], inset, xmidlen, ymidlen, zmidlen);    
    
    MPI_Request reqs[50];		// a non-tight upper bound on how many we need
    MPI_Status stats[50];
    short rused=0;
    
    messvec incoms;
    messvec outcoms;
    
    grid.generateInNeighbours(X, Y, Z ,incoms);
    grid.generateOutNeighbours(X, Y, Z, outcoms);
    
    
    block.copyUsedFromBuffer(src);
    
    
    int comserr=0;    
    for (size_t i=0;i<incoms.size();++i)
    {
	message& m=incoms[i];
	comserr|=MPI_Irecv(block.getInBuffer(m.destbuffid), block.getBuffSize(m.destbuffid) , MPI_DOUBLE, m.sourceID, m.tag, m_mpiInfo->comm, reqs+(rused++));
	block.setUsed(m.destbuffid);
    }

    for (size_t i=0;i<outcoms.size();++i)
    {
	message& m=outcoms[i];
	comserr|=MPI_Isend(block.getOutBuffer(m.srcbuffid), block.getBuffSize(m.srcbuffid) , MPI_DOUBLE, m.destID, m.tag, m_mpiInfo->comm, reqs+(rused++));
    }    
    
    if (!comserr)
    {     
        comserr=MPI_Waitall(rused, reqs, stats);
    }

    if (comserr)
    {
	// Yes this is throwing an exception as a result of an MPI error.
	// and no we don't inform the other ranks that we are doing this.
	// however, we have no reason to believe coms work at this point anyway
        throw RipleyException("Error in coms for randomFill");      
    }
    
    block.copyUsedFromBuffer(src);

#endif    
    escript::FunctionSpace fs(getPtr(), getContinuousFunctionCode());
    escript::Data resdat(0, escript::DataTypes::scalarShape, fs , true);
    // don't need to check for exwrite because we just made it
    escript::DataVector& dv=resdat.getExpandedVectorReference();
    double* convolution=get3DGauss(radius, sigma);
    for (size_t z=0;z<(internal[2]);++z)
    {
	for (size_t y=0;y<(internal[1]);++y)    
	{
	    for (size_t x=0;x<(internal[0]);++x)
	    {	  
		dv[x+y*(internal[0])+z*internal[0]*internal[1]]=Convolve3D(convolution, src, x+radius, y+radius, z+radius, radius, ext[0], ext[1]);
		
	    }
	}
    }
    delete[] convolution;
    delete[] src;
    return resdat;
}






int Brick::findNode(const double *coords) const {
    const int NOT_MINE = -1;
    //is the found element even owned by this rank
    for (int dim = 0; dim < m_numDim; dim++) {
        if (m_origin[dim] + m_offset[dim] > coords[dim]  || m_origin[dim] 
                + m_offset[dim] + m_dx[dim]*m_ownNE[dim] < coords[dim]) {
            return NOT_MINE;
        }
    }
    // get distance from origin
    double x = coords[0] - m_origin[0];
    double y = coords[1] - m_origin[1];
    double z = coords[2] - m_origin[2];
    // distance in elements
    int ex = (int) floor(x / m_dx[0]);
    int ey = (int) floor(y / m_dx[1]);
    int ez = (int) floor(z / m_dx[2]);
    // set the min distance high enough to be outside the element plus a bit
    int closest = NOT_MINE;
    double minDist = 1;
    for (int dim = 0; dim < m_numDim; dim++) {
        minDist += m_dx[dim]*m_dx[dim];
    }
    //find the closest node
    for (int dx = 0; dx < 2; dx++) {
        double xdist = x - (ex + dx)*m_dx[0];
        for (int dy = 0; dy < 2; dy++) {
            double ydist = y - (ey + dy)*m_dx[1];
            for (int dz = 0; dz < 2; dz++) {
                double zdist = z - (ez + dz)*m_dx[2];
                double total = xdist*xdist + ydist*ydist + zdist*zdist;
                if (total < minDist) {
                    closest = INDEX3(ex+dy, ey+dy, ez+dz, m_NE[0]+1, m_NE[1]+1);
                    minDist = total;
                }
            }
        }
    }
    return closest;
}

void Brick::setAssembler(std::string type, std::map<std::string,
        escript::Data> constants) {
    if (type.compare("WaveAssembler") == 0) {
        delete assembler;
        assembler = new WaveAssembler3D(this, m_dx, m_NX, m_NE, m_NN, constants);
    } else { //else ifs would go before this for other types
        throw RipleyException("Ripley::Rectangle does not support the"
                                " requested assembler");
    }
}

} // end of namespace ripley

