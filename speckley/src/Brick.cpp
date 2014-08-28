
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include <speckley/Brick.h>
#include <speckley/DefaultAssembler3D.h>
#include <speckley/blocktools.h>
#include <esysUtils/esysFileWriter.h>
#include <esysUtils/EsysRandom.h>
#include <esysUtils/index.h>
#include <escript/FunctionSpaceFactory.h>

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
#include <limits>

using namespace std;
using esysUtils::FileWriter;

namespace speckley {

int indexOfMax(int a, int b, int c) {
    if (a > b) {
        if (c > a) {
            return 2;
        }
        return 0;
    } else if (b > c) {
        return 1;
    }
    return 2;
}

Brick::Brick(int order, int n0, int n1, int n2, double x0, double y0, double z0,
             double x1, double y1, double z1, int d0, int d1, int d2,
             const std::vector<double>& points, const std::vector<int>& tags,
             const simap_t& tagnamestonums,
             escript::SubWorld_ptr w) :
    SpeckleyDomain(3, order, w)
{
    if (static_cast<long>(n0 + 1) * static_cast<long>(n1 + 1) 
            * static_cast<long>(n2 + 1) > std::numeric_limits<int>::max())
        throw SpeckleyException("The number of elements has overflowed, this "
                "limit may be raised in future releases.");

    if (n0 <= 0 || n1 <= 0 || n2 <= 0)
        throw SpeckleyException("Number of elements in each spatial dimension "
                "must be positive");

    // ignore subdivision parameters for serial run
    if (m_mpiInfo->size == 1) {
        d0=1;
        d1=1;
        d2=1;
    } else {
        throw SpeckleyException("MPI not currently supported by Brick");
    }
    bool warn=false;

    std::vector<int> factors;
    int ranks = m_mpiInfo->size;
    int epr[3] = {n0,n1,n2};
    int d[3] = {d0,d1,d2};
    if (d0<=0 || d1<=0 || d2<=0) {
        for (int i = 0; i < 3; i++) {
            if (d[i] < 1) {
                d[i] = 1;
                continue;
            }
            epr[i] = -1; // can no longer be max
            if (ranks % d[i] != 0) {
                throw SpeckleyException("Invalid number of spatial subdivisions");
            }
            ranks /= d[i];
        }
        factorise(factors, ranks);
        if (factors.size() != 0) {
            warn = true;
        }
    } 
    while (factors.size() > 0) {
        int i = indexOfMax(epr[0],epr[1],epr[2]);
        int f = factors.back();
        factors.pop_back();
        d[i] *= f;
        epr[i] /= f;
    }
    d0 = d[0]; d1 = d[1]; d2 = d[2];

    // ensure number of subdivisions is valid and nodes can be distributed
    // among number of ranks
    if (d0*d1*d2 != m_mpiInfo->size)
        throw SpeckleyException("Invalid number of spatial subdivisions");

    if (warn) {
        std::cout << "Warning: Automatic domain subdivision (d0=" << d0 << ", d1="
            << d1 << ", d2=" << d2 << "). This may not be optimal!" << std::endl;
    }

    double l0 = x1-x0;
    double l1 = y1-y0;
    double l2 = z1-z0;
    m_dx[0] = l0/n0;
    m_dx[1] = l1/n1;
    m_dx[2] = l2/n2;

    if (n0 % d0 > 0) {
        n0 += d0 - (n0 % d0);
        l0 = m_dx[0]*n0;
        std::cout << "Warning: Adjusted number of elements and length. N0="
            << n0 << ", l0=" << l0 << std::endl;
    }
    if (n1 % d1 > 0) {
        n1 += d1 - (n1 % d1);
        l1 = m_dx[1]*n1;
        std::cout << "Warning: Adjusted number of elements and length. N1="
            << n1 << ", l1=" << l1 << std::endl;
    }
    if (n2 % d2 > 0) {
        n2 += d2 - (n2 % d2);
        l2 = m_dx[2]*n2;
        std::cout << "Warning: Adjusted number of elements and length. N2="
            << n2 << ", l2=" << l2 << std::endl;
    }

    if (n0/d0 < 1 || n1/d1 < 1 || n2/d2 < 1)
        throw SpeckleyException("Too few elements for the number of ranks");

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
    m_NE[0] = n0 / d0;
    m_NE[1] = n1 / d1;
    m_NE[2] = n2 / d2;

    // local number of nodes
    m_NN[0] = m_NE[0]*m_order+1;
    m_NN[1] = m_NE[1]*m_order+1;
    m_NN[2] = m_NE[2]*m_order+1;

    // bottom-left-front node is at (offset0,offset1,offset2) in global mesh
    m_offset[0] = n0/d0*(m_mpiInfo->rank%d0);
    m_offset[1] = n1/d1*(m_mpiInfo->rank%(d0*d1)/d0);
    m_offset[2] = n2/d2*(m_mpiInfo->rank/(d0*d1));

    populateSampleIds();
    for (std::map<std::string, int>::const_iterator i = tagnamestonums.begin();
            i != tagnamestonums.end(); i++) {
        setTagMap(i->first, i->second);
    }
    addPoints(points, tags);
}


Brick::~Brick()
{
}

string Brick::getDescription() const
{
    return "speckley::Brick";
}

bool Brick::operator==(const AbstractDomain& other) const
{
    const Brick* o=dynamic_cast<const Brick*>(&other);
    if (o) {
        return (SpeckleyDomain::operator==(other) &&
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
        throw SpeckleyException("readNcGrid(): invalid function space for output data object");

    if (params.first.size() != 3)
        throw SpeckleyException("readNcGrid(): argument 'first' must have 3 entries");

    if (params.numValues.size() != 3)
        throw SpeckleyException("readNcGrid(): argument 'numValues' must have 3 entries");

    if (params.multiplier.size() != 3)
        throw SpeckleyException("readNcGrid(): argument 'multiplier' must have 3 entries");
    for (size_t i=0; i<params.multiplier.size(); i++)
        if (params.multiplier[i]<1)
            throw SpeckleyException("readNcGrid(): all multipliers must be positive");

    // check file existence and size
    NcFile f(filename.c_str(), NcFile::ReadOnly);
    if (!f.is_valid())
        throw SpeckleyException("readNcGrid(): cannot open file");

    NcVar* var = f.get_var(varname.c_str());
    if (!var)
        throw SpeckleyException("readNcGrid(): invalid variable name");

    // TODO: rank>0 data support
    const int numComp = out.getDataPointSize();
    if (numComp > 1)
        throw SpeckleyException("readNcGrid(): only scalar data supported");

    const int dims = var->num_dims();
    boost::scoped_array<long> edges(var->edges());

    // is this a slice of the data object (dims!=3)?
    // note the expected ordering of edges (as in numpy: z,y,x)
    if ( (dims==3 && (params.numValues[2] > edges[0] ||
                      params.numValues[1] > edges[1] ||
                      params.numValues[0] > edges[2]))
            || (dims==2 && params.numValues[2]>1)
            || (dims==1 && (params.numValues[2]>1 || params.numValues[1]>1)) ) {
        throw SpeckleyException("readNcGrid(): not enough data in file");
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
    throw SpeckleyException("readNcGrid(): not compiled with netCDF support");
#endif
}

#ifdef USE_BOOSTIO
void Brick::readBinaryGridFromZipped(escript::Data& out, string filename,
                           const ReaderParameters& params) const
{
    // the mapping is not universally correct but should work on our
    // supported platforms
    switch (params.dataType) {
        case DATATYPE_INT32:
            readBinaryGridZippedImpl<int>(out, filename, params);
            break;
        case DATATYPE_FLOAT32:
            readBinaryGridZippedImpl<float>(out, filename, params);
            break;
        case DATATYPE_FLOAT64:
            readBinaryGridZippedImpl<double>(out, filename, params);
            break;
        default:
            throw SpeckleyException("readBinaryGrid(): invalid or unsupported datatype");
    }
}
#endif

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
            throw SpeckleyException("readBinaryGrid(): invalid or unsupported datatype");
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
        throw SpeckleyException("readBinaryGrid(): invalid function space for output data object");

    if (params.first.size() != 3)
        throw SpeckleyException("readBinaryGrid(): argument 'first' must have 3 entries");

    if (params.numValues.size() != 3)
        throw SpeckleyException("readBinaryGrid(): argument 'numValues' must have 3 entries");

    if (params.multiplier.size() != 3)
        throw SpeckleyException("readBinaryGrid(): argument 'multiplier' must have 3 entries");
    for (size_t i=0; i<params.multiplier.size(); i++)
        if (params.multiplier[i]<1)
            throw SpeckleyException("readBinaryGrid(): all multipliers must be positive");
    if (params.reverse[0] != 0 || params.reverse[1] != 0)
        throw SpeckleyException("readBinaryGrid(): reversing only supported in Z-direction currently");

    // check file existence and size
    ifstream f(filename.c_str(), ifstream::binary);
    if (f.fail()) {
        throw SpeckleyException("readBinaryGrid(): cannot open file");
    }
    f.seekg(0, ios::end);
    const int numComp = out.getDataPointSize();
    const int filesize = f.tellg();
    const int reqsize = params.numValues[0]*params.numValues[1]*params.numValues[2]*numComp*sizeof(ValueType);
    if (filesize < reqsize) {
        f.close();
        throw SpeckleyException("readBinaryGrid(): not enough data in file");
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
    // indices to first value in file (not accounting for reverse yet)
    int idx0 = max(0, m_offset[0]-params.first[0]);
    int idx1 = max(0, m_offset[1]-params.first[1]);
    int idx2 = max(0, m_offset[2]-params.first[2]);
    // number of values to read
    const int num0 = min(params.numValues[0]-idx0, myN0-first0);
    const int num1 = min(params.numValues[1]-idx1, myN1-first1);
    const int num2 = min(params.numValues[2]-idx2, myN2-first2);

    // make sure we read the right block if going backwards through file
    if (params.reverse[2])
        idx2 = params.numValues[2]-idx2-1;

    // helpers for reversing
    const int z_mult = (params.reverse[2] ? -1 : 1);

    out.requireWrite();
    vector<ValueType> values(num0*numComp);
    const int dpp = out.getNumDataPointsPerSample();

    for (int z=0; z<num2; z++) {
        for (int y=0; y<num1; y++) {
            const int fileofs = numComp*(idx0 +
                                (idx1+y)*params.numValues[0] +
                                (idx2+z_mult*z)*params.numValues[0]*params.numValues[1]);
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
                                if (!isnan(val)) {
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

#ifdef USE_BOOSTIO
template<typename ValueType>
void Brick::readBinaryGridZippedImpl(escript::Data& out, const string& filename,
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
        throw SpeckleyException("readBinaryGridFromZipped(): invalid function space for output data object");

    if (params.first.size() != 3)
        throw SpeckleyException("readBinaryGridFromZipped(): argument 'first' must have 3 entries");

    if (params.numValues.size() != 3)
        throw SpeckleyException("readBinaryGridFromZipped(): argument 'numValues' must have 3 entries");

    if (params.multiplier.size() != 3)
        throw SpeckleyException("readBinaryGridFromZipped(): argument 'multiplier' must have 3 entries");
    for (size_t i=0; i<params.multiplier.size(); i++)
        if (params.multiplier[i]<1)
            throw SpeckleyException("readBinaryGridFromZipped(): all multipliers must be positive");

    // check file existence and size
    ifstream f(filename.c_str(), ifstream::binary);
    if (f.fail()) {
        throw SpeckleyException("readBinaryGridFromZipped(): cannot open file");
    }
    f.seekg(0, ios::end);
    const int numComp = out.getDataPointSize();
    int filesize = f.tellg();
    f.seekg(0, ios::beg);
    std::vector<char> compressed(filesize);
    f.read((char*)&compressed[0], filesize);
    f.close();
    std::vector<char> decompressed = unzip(compressed);
    filesize = decompressed.size();
    const int reqsize = params.numValues[0]*params.numValues[1]*params.numValues[2]*numComp*sizeof(ValueType);
    if (filesize < reqsize) {
        throw SpeckleyException("readBinaryGridFromZipped(): not enough data in file");
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
            memcpy((char*)&values[0], (char*)&decompressed[fileofs*sizeof(ValueType)], num0*numComp*sizeof(ValueType));
            
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
                                if (!isnan(val)) {
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
}
#endif

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
            throw SpeckleyException("writeBinaryGrid(): invalid or unsupported datatype");
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
        throw SpeckleyException("writeBinaryGrid(): invalid function space of data object");

    const int numComp = in.getDataPointSize();
    const int dpp = in.getNumDataPointsPerSample();
    const int fileSize = sizeof(ValueType)*numComp*dpp*totalN0*totalN1*totalN2;

    if (numComp > 1 || dpp > 1)
        throw SpeckleyException("writeBinaryGrid(): only scalar, single-value data supported");

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
                    if (sizeof(fvalue)>4) {
                        byte_swap64(value);
                    } else {
                        byte_swap32(value);
                    }
                    oss.write(value, sizeof(fvalue));
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
        throw SpeckleyException("dump: Could not create Silo file");

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
    const int NN0 = m_NN[0];
    const int NN1 = m_NN[1];
    const int NN2 = m_NN[2];

#pragma omp parallel
    {
#pragma omp for
        for (dim_t i0 = 0; i0 < NN0; i0++) {
            coords[0][i0]=getLocalCoordinate(i0, 0);
        }
#pragma omp for
        for (dim_t i1 = 0; i1 < NN1; i1++) {
            coords[1][i1]=getLocalCoordinate(i1, 1);
        }
#pragma omp for
        for (dim_t i2 = 0; i2 < NN2; i2++) {
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
    throw SpeckleyException("dump: no Silo support");
#endif
}

const int* Brick::borrowSampleReferenceIDs(int fsType) const
{
    switch (fsType) {
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom: // FIXME: reduced
//            return &m_dofId[0];
        case Nodes:
        case ReducedNodes: // FIXME: reduced
            return &m_nodeId[0];
        case Elements:
        case ReducedElements:
            return &m_elementId[0];
        case FaceElements:
        case ReducedFaceElements:
            return &m_faceId[0];
        case Points:
            return &m_diracPointNodeIDs[0];
        default:
            break;
    }

    std::stringstream msg;
    msg << "borrowSampleReferenceIDs: invalid function space type " << fsType;
    throw SpeckleyException(msg.str());
}

bool Brick::ownSample(int fsType, index_t id) const
{
    throw SpeckleyException("ownSample not implemented");
}

void Brick::setToNormal(escript::Data& out) const
{
    throw SpeckleyException("setToNormal not implemented");
}

void Brick::setToSize(escript::Data& out) const
{
    std::cerr << "Brick::setToSize not updated\n";
    if (out.getFunctionSpace().getTypeCode() == Elements
            || out.getFunctionSpace().getTypeCode() == ReducedElements) {
        out.requireWrite();
        const dim_t numQuad = out.getNumDataPointsPerSample();
        const double size=sqrt(m_dx[0]*m_dx[0]+m_dx[1]*m_dx[1]+m_dx[2]*m_dx[2]);
        const int numElements = getNumElements();
#pragma omp parallel for
        for (index_t k = 0; k < numElements; ++k) {
            double* o = out.getSampleDataRW(k);
            std::fill(o, o+numQuad, size);
        }
    } else if (out.getFunctionSpace().getTypeCode() == FaceElements
            || out.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        out.requireWrite();
        const dim_t numQuad=out.getNumDataPointsPerSample();
        const int NE0 = m_NE[0];
        const int NE1 = m_NE[1];
        const int NE2 = m_NE[2];
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
                const double size=min(m_dx[1],m_dx[2]);
#pragma omp for nowait
                for (index_t k2 = 0; k2 < NE2; ++k2) {
                    for (index_t k1 = 0; k1 < NE1; ++k1) {
                        double* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]));
                        fill(o, o+numQuad, size);
                    }
                }
            }

            if (m_faceOffset[1] > -1) {
                const double size=min(m_dx[1],m_dx[2]);
#pragma omp for nowait
                for (index_t k2 = 0; k2 < NE2; ++k2) {
                    for (index_t k1 = 0; k1 < NE1; ++k1) {
                        double* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]));
                        fill(o, o+numQuad, size);
                    }
                }
            }

            if (m_faceOffset[2] > -1) {
                const double size=min(m_dx[0],m_dx[2]);
#pragma omp for nowait
                for (index_t k2 = 0; k2 < NE2; ++k2) {
                    for (index_t k0 = 0; k0 < NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]));
                        fill(o, o+numQuad, size);
                    }
                }
            }

            if (m_faceOffset[3] > -1) {
                const double size=min(m_dx[0],m_dx[2]);
#pragma omp for nowait
                for (index_t k2 = 0; k2 < NE2; ++k2) {
                    for (index_t k0 = 0; k0 < NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]));
                        fill(o, o+numQuad, size);
                    }
                }
            }

            if (m_faceOffset[4] > -1) {
                const double size=min(m_dx[0],m_dx[1]);
#pragma omp for nowait
                for (index_t k1 = 0; k1 < NE1; ++k1) {
                    for (index_t k0 = 0; k0 < NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]));
                        fill(o, o+numQuad, size);
                    }
                }
            }

            if (m_faceOffset[5] > -1) {
                const double size=min(m_dx[0],m_dx[1]);
#pragma omp for nowait
                for (index_t k1 = 0; k1 < NE1; ++k1) {
                    for (index_t k0 = 0; k0 < NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]));
                        fill(o, o+numQuad, size);
                    }
                }
            }
        } // end of parallel section

    } else {
        std::stringstream msg;
        msg << "setToSize: invalid function space type "
            << out.getFunctionSpace().getTypeCode();
        throw SpeckleyException(msg.str());
    }
}

void Brick::Print_Mesh_Info(const bool full) const
{
    SpeckleyDomain::Print_Mesh_Info(full);
    if (full) {
        std::cout << "     Id  Coordinates" << std::endl;
        std::cout.precision(15);
        std::cout.setf(std::ios::scientific, std::ios::floatfield);
        for (index_t i=0; i < getNumNodes(); i++) {
            std::cout << "  " << std::setw(5) << m_nodeId[i]
                << "  " << getLocalCoordinate(i%m_NN[0], 0)
                << "  " << getLocalCoordinate(i%(m_NN[0]*m_NN[1])/m_NN[0], 1)
                << "  " << getLocalCoordinate(i/(m_NN[0]*m_NN[1]), 2) << std::endl;
        }
    }
}


//protected
void Brick::assembleCoordinates(escript::Data& arg) const
{
    const int index = m_order-2;
    const int max_x = m_order*m_NE[0] + 1;
    const int max_y = m_order*m_NE[1] + 1;
    const int max_z = m_order*m_NE[2] + 1;

    escriptDataC dc = arg.getDataC();
    int numDim = m_numDim;
    if (!isDataPointShapeEqual(&dc, 1, &numDim))
        throw SpeckleyException("setToX: Invalid Data object shape");
    if (!numSamplesEqual(&dc, 1, getNumNodes()))
        throw SpeckleyException("setToX: Illegal number of samples in Data object");

    arg.requireWrite();
#pragma omp parallel for
    for (dim_t z = 0; z < max_z; z++) {
        for (dim_t y = 0; y < max_y; y++) {
            for (dim_t x = 0; x < max_x; x++) {
                double *point = arg.getSampleDataRW(INDEX3(x,y,z,max_x,max_y));
                point[0] = m_origin[0] + m_dx[0]*(m_offset[0] + x/m_order + point_locations[index][x%m_order]);
                point[1] = m_origin[1] + m_dx[1]*(m_offset[1] + y/m_order + point_locations[index][y%m_order]);
                point[2] = m_origin[2] + m_dx[2]*(m_offset[2] + z/m_order + point_locations[index][z%m_order]);
            }
        }
    }
}

//protected
void Brick::assembleGradient(escript::Data& out, const escript::Data& in) const
{
    escript::Data converted;

    if (in.getFunctionSpace().getTypeCode() != Elements) {
        converted = escript::Data(in, escript::function(*this));
    } else {
        converted = in;
    }

    if (m_order == 2) {
        gradient_order2(out,converted);
    } else if (m_order == 3) {
        gradient_order3(out,converted);
    } else if (m_order == 4) {
        gradient_order4(out,converted);
    } else if (m_order == 5) {
        gradient_order5(out,converted);
    } else if (m_order == 6) {
        gradient_order6(out,converted);
    } else if (m_order == 7) {
        gradient_order7(out,converted);
    } else if (m_order == 8) {
        gradient_order8(out,converted);
    } else if (m_order == 9) {
        gradient_order9(out,converted);
    } else if (m_order == 10) {
        gradient_order10(out,converted);
    }
}

//protected
void Brick::assembleIntegrate(vector<double>& integrals, const escript::Data& arg) const
{
    const int fs = arg.getFunctionSpace().getTypeCode();
    if (fs != Elements)
        throw new SpeckleyException("Speckley doesn't currently support integrals of non-Element functionspaces");
    if (!arg.actsExpanded())
        throw new SpeckleyException("Speckley doesn't currently support unexpanded data");
    
    if (m_order == 2) {
        integral_order2(integrals, arg);
    } else if (m_order == 3) {
        integral_order3(integrals, arg);
    } else if (m_order == 4) {
        integral_order4(integrals, arg);
    } else if (m_order == 5) {
        integral_order5(integrals, arg);
    } else if (m_order == 6) {
        integral_order6(integrals, arg);
    } else if (m_order == 7) {
        integral_order7(integrals, arg);
    } else if (m_order == 8) {
        integral_order8(integrals, arg);
    } else if (m_order == 9) {
        integral_order9(integrals, arg);
    } else if (m_order == 10) {
        integral_order10(integrals, arg);
    }
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

    for (dim_t k = 0; k < m_mpiInfo->size - 1; k++) {
        index_t rank_left = k % m_NX[0] == 0 ? 0 : 1;
        index_t rank_bottom = k % (m_NX[0]*m_NX[1])/m_NX[0] == 0 ? 0 : 1;
        index_t rank_front = k / (m_NX[0]*m_NX[1]) == 0 ? 0 : 1;
        m_nodeDistribution[k+1] = m_nodeDistribution[k] 
                                + (m_NN[0]-rank_left) * (m_NN[1]-rank_bottom)
                                * (m_NN[2]-rank_front);
    }

    m_nodeDistribution[m_mpiInfo->size]=getNumDataPointsGlobal();
    
    try {
        m_nodeId.resize(getNumNodes());
        m_elementId.resize(getNumElements());
    } catch (const std::length_error& le) {
        throw SpeckleyException("The system does not have sufficient memory for a domain of this size.");
    }

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

    const dim_t NFE = getNumFaceElements();
    m_faceId.resize(NFE);

//    const dim_t NN0 = m_NN[0];
//    const dim_t NN1 = m_NN[1];
//    const dim_t NN2 = m_NN[2];

    const int left = (m_offset[0]==0 ? 0 : 1);
    const int front = (m_offset[1]==0 ? 0 : 1);
    const int bottom = (m_offset[2]==0 ? 0 : 1);

    const int rank = m_mpiInfo->rank;

    //re-use bottom-left-front corner node
    if (bottom && left && front) {
        //get lower-left-front node
        int rank_wanted = rank 
                    - m_NX[0]*m_NX[1]   //down a layer
                    - m_NX[0]           //towards the front
                    - 1;                //to the left
        m_nodeId[0] = m_nodeDistribution[rank_wanted] - 1;
    }
    //re-use nodes on front border
    if (front) {
        int neighbour = rank - m_NX[0];
        index_t begin = m_nodeDistribution[neighbour] + (m_NN[0] - 1)*m_NN[1];
#pragma omp parallel for
        for (index_t z = bottom; z < m_NN[2]; z++) {
            for (index_t x = left; x < m_NN[0]; x++) {
                index_t i = x + z*m_NN[0]*m_NN[1];
                m_nodeId[i] = begin + x + z*m_NN[0]*(m_NN[1]);
            }
        }
    }
    //re-use nodes on left border
    if (left) {
        //is the rank to the left itself right of another rank
        index_t rank_left = (rank - 1) % m_NX[0] == 0 ? 0 : 1;
        //end of first owned row of neighbouring rank
        index_t end = m_nodeDistribution[m_mpiInfo->rank - 1]
                    + m_NN[0] - rank_left - 1;
#pragma omp parallel for
        for (index_t z = bottom; z < m_NN[2]; z++) {
            for (index_t y = front; y < m_NN[1]; y++) {
                index_t i = INDEX3(0,y,z,m_NN[0],m_NN[1]);
                m_nodeId[i] = end + INDEX3(0,(y-bottom),(z-front),(m_NN[0]-rank_left),(m_NN[2]-bottom));
            }
        }
    }
    //re-use nodes on bottom border
    if (bottom) {
        int neighbour = rank - m_NX[0]*m_NX[1];
        //DOF size of the neighbouring rank 
        index_t rankDOF = m_nodeDistribution[neighbour + 1] 
                        - m_nodeDistribution[neighbour];
        //beginning of first row
        index_t begin = m_nodeDistribution[neighbour] + rankDOF - m_NN[0]*m_NN[1];
#pragma omp parallel for
        for (index_t y = front; y < m_NN[1]; y++) {
            for (index_t x = left; x < m_NN[0]; x++) {
                index_t i = INDEX2(x,y,m_NN[0]);
                m_nodeId[i] = begin + INDEX2(x,y,m_NN[0]);
//                fprintf(stderr, "r%d: B %d -> %d\n", rank, i, m_nodeId[i]);
            }
        }
    }
    
    //now the new nodes
    const index_t start = m_nodeDistribution[m_mpiInfo->rank];
#pragma omp parallel for
    for (index_t z = front; z < m_NN[2]; z++) {
        index_t z_chunk = (z-front)*(m_NN[0]-left)*(m_NN[1]-bottom);
        for (index_t y = bottom; y < m_NN[1]; y++) {
            index_t y_chunk = (y-bottom)*(m_NN[0]-left);
            for (index_t x = left; x < m_NN[0]; x++) {
                m_nodeId[INDEX3(x, y, z, m_NN[0], m_NN[1])] 
                        =  start + z_chunk + y_chunk + (x-left);
            }
        }
    }
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
void Brick::addToMatrixAndRHS(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<double>& EM_S, const std::vector<double>& EM_F, bool addS,
         bool addF, index_t firstNode, dim_t nEq, dim_t nComp) const
{
    throw SpeckleyException("Rectangle::addToMatrixAndRHS not implemented");
}

void Brick::interpolateElementsOnNodes(escript::Data& out,
                                  const escript::Data& in, bool reduced) const {
    if (reduced)
        throw SpeckleyException("Speckley does not support "
                "reduced functionspaces");
    const dim_t numComp = in.getDataPointSize();
    const int NE0 = m_NE[0];
    const int NE1 = m_NE[1];
    const int NE2 = m_NE[2];
    const int quads = m_order + 1;
    const int max_x = (m_order*NE0) + 1;
    const int max_y = (m_order*NE1) + 1;
    const int max_z = (m_order*NE2) + 1;
    out.requireWrite();
    //init to zero so we can do some sums without undefined, may not be required
    memset(out.getSampleDataRW(0), 0, sizeof(double)*quads*quads*numComp);
    // the summation portion
    for (dim_t colouring = 0; colouring < 2; colouring++) {
#pragma omp parallel for
        for (dim_t ez = colouring; ez < NE2; ez += 2) {
            for (dim_t ey = 0; ey < NE1; ey++) {
                for (dim_t ex = 0; ex < NE0; ex++) {
                    int start = m_order * (INDEX3(ex, ey, ez, max_x, max_y));
                    const double *e_in = in.getSampleDataRO(INDEX3(ex,ey,ez,NE0,NE1));
                    for (int qz = 0; qz < quads; qz++) {
                        for (int qy = 0; qy < quads; qy++) {
                            for (int qx = 0; qx < quads; qx++) {
                                double *n_out = out.getSampleDataRW(start + INDEX3(qx, qy, qz, max_x, max_y));
                                for (dim_t comp = 0; comp < numComp; comp++) {
                                    n_out[comp] += e_in[INDEX4(comp, qx, qy, qz, numComp, quads, quads)];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    /* the averaging out (each point divided by the number of additions in
       the summation step). By doing each edge along an axis, those points
       requiring division by 2, 4 and 8 will be divided by 2 the right 
       number of times
       
       border edges are skipped because they aren't shared but for the 
       points that lie along a different axes non-border edge
    */
    // for every non-border edge in x
#pragma omp parallel for
    for (int qz = 0; qz < max_z; qz++) {
        for (int qy = 0; qy < max_y; qy++) {
            for (int qx = m_order; qx < max_x - m_order; qx += m_order) {
                double *n_out = out.getSampleDataRW(INDEX3(qx, qy, qz, max_x, max_y));
                for (int comp = 0; comp < numComp; comp++) {
                    n_out[comp] /= 2.;
                }
            }
        }
    }
    // for every non-border edge in y
#pragma omp parallel for
    for (int qz = 0; qz < max_z; qz++) {
        for (int qy = m_order; qy < max_y - m_order; qy += m_order) {
            for (int qx = 0; qx < max_x; qx ++) {
                double *n_out = out.getSampleDataRW(INDEX3(qx, qy, qz, max_x, max_y));
                for (int comp = 0; comp < numComp; comp++) {
                    n_out[comp] /= 2.;
                }
            }
        }
    }
    // for every non-border edge in z
#pragma omp parallel for
    for (int qz = m_order; qz < max_z - m_order; qz += m_order) {
        for (int qy = 0; qy < max_y; qy++) {
            for (int qx = 0; qx < max_x; qx++) {
                double *n_out = out.getSampleDataRW(INDEX3(qx, qy, qz, max_x, max_y));
                for (int comp = 0; comp < numComp; comp++) {
                    n_out[comp] /= 2.;
                }
            }
        }
    }
}

//protected
void Brick::interpolateNodesOnElements(escript::Data& out,
                                       const escript::Data& in,
                                       bool reduced) const
{
    const dim_t numComp = in.getDataPointSize();
    const int NE0 = m_NE[0];
    const int NE1 = m_NE[1];
    const int NE2 = m_NE[2];
    const int quads = m_order + 1;
    const int max_x = (m_order*NE0) + 1;
    const int max_y = (m_order*NE1) + 1;
    out.requireWrite();
#pragma omp parallel for
    for (dim_t ez = 0; ez < NE2; ez++) {
        for (dim_t ey = 0; ey < NE1; ey++) {
            for (dim_t ex = 0; ex < NE0; ex++) {
                double *e_out = out.getSampleDataRW(INDEX3(ex, ey, ez, NE0, NE1));
                int start = m_order * INDEX3(ex, ey, ez, max_x, max_y);
                int quad = 0;
                for (int qz = 0; qz < quads; qz++) {
                    for (int qy = 0; qy < quads; qy++) {
                        for (int qx = 0; qx < quads; qx++, quad++) {
                            const double *n_in = in.getSampleDataRO(start + INDEX3(qx,qy,qz,max_x,max_y));
                            for (int comp = 0; comp < numComp; comp++) {
                                e_out[INDEX4(comp, qx, qy, qz, numComp, quads, quads)] = n_in[comp];
                            }
                        }
                    }
                }
            }
        }
    }
}

//protected
void Brick::interpolateNodesOnFaces(escript::Data& out, const escript::Data& in,
                                    bool reduced) const
{
    throw SpeckleyException("Brick::interpolateNodesOnFaces not implemented");
}

escript::Data Brick::randomFill(const escript::DataTypes::ShapeType& shape,
       const escript::FunctionSpace& fs,
       long seed, const boost::python::tuple& filter) const
{
    int numvals=escript::DataTypes::noValues(shape);
    int per_element = (m_order+1)*(m_order+1)*(m_order+1)*numvals;
    if (len(filter)>0) {
        throw SpeckleyException("Speckley does not support filters.");
    }

    double* src=new double[m_NE[0]*m_NE[1]*m_NE[2]*per_element*numvals];
    esysUtils::randomFillArray(seed, src, m_NE[0]*m_NE[1]*m_NE[2]*per_element);
    escript::Data res(0, shape, escript::function(*this), true);
    int current = 0;
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                double *e = res.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                memcpy(e, &src[current], sizeof(double)*per_element);
                current += per_element;
            }
        }
    }
    if (res.getFunctionSpace() != fs) {
        return escript::Data(res, fs);
    }
    return res;
}

escript::Data Brick::randomFillWorker(const escript::DataTypes::ShapeType& shape,
        long seed, const boost::python::tuple& filter) const {
    throw SpeckleyException("Brick::randomFillWorker not yet implemented");
}

int Brick::findNode(const double *coords) const {
    const int NOT_MINE = -1;
    //is the found element even owned by this rank
    for (int dim = 0; dim < m_numDim; dim++) {
        double min = m_origin[dim] + m_offset[dim]* m_dx[dim]
                - m_dx[dim]/2.; //allows for point outside mapping onto node
        double max = m_origin[dim] + (m_offset[dim] + m_NE[dim])*m_dx[dim]
                + m_dx[dim]/2.;
        if (min > coords[dim] || max < coords[dim]) {
            return NOT_MINE;
        }
    }
    // get distance from origin
    double x = coords[0] - m_origin[0];
    double y = coords[1] - m_origin[1];
    double z = coords[2] - m_origin[2];
    
    //check if the point is even inside the domain
    if (x < 0 || y < 0 || z < 0 
            || x > m_length[0] || y > m_length[1] || z > m_length[2])
        return NOT_MINE;
        
    // distance in elements
    int ex = (int) floor(x / m_dx[0] + 0.01*m_dx[0]);
    int ey = (int) floor(y / m_dx[1]+ 0.01*m_dx[1]);
    int ez = (int) floor(z / m_dx[2]+ 0.01*m_dx[2]);
    int start = m_order*(INDEX3(ex,ey,ez,m_NN[0],m_NN[1]));
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
                    closest = start + INDEX3(m_order*dx,dy,dz,m_NN[0],m_NN[1]);
                    minDist = total;
                }
            }
        }
    }
    //if this happens, we've let a dirac point slip through, which is awful
    if (closest == NOT_MINE) {
        throw SpeckleyException("Unable to map appropriate dirac point to a "
                "node, implementation problem in Brick::findNode()");
    }
    return closest;
}

Assembler_ptr Brick::createAssembler(std::string type,
        const DataMap& options) const {
    if (type.compare("DefaultAssembler") == 0) {
        return Assembler_ptr(new DefaultAssembler3D(shared_from_this(), m_dx,
                m_NX, m_NE, m_NN));
    } else { //else ifs would go before this for other types
        throw SpeckleyException("Speckley::Brick does not support the"
                                " requested assembler");
    }
}

} // end of namespace speckley

