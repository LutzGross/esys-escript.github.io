
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
#include <esysUtils/esysFileWriter.h>
#include <esysUtils/EsysRandom.h>
#include <esysUtils/index.h>
#include <escript/FunctionSpaceFactory.h>

#include <boost/scoped_array.hpp>

#ifdef USE_RIPLEY
#include <speckley/CrossDomainCoupler.h>
#endif

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

inline int indexOfMax(dim_t a, dim_t b, dim_t c)
{
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

Brick::Brick(int order, dim_t n0, dim_t n1, dim_t n2, double x0, double y0, double z0,
             double x1, double y1, double z1, int d0, int d1, int d2,
             const std::vector<double>& points, const std::vector<int>& tags,
             const TagMap& tagnamestonums,
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
    }
    bool warn=false;

    std::vector<int> factors;
    int ranks = m_mpiInfo->size;
    dim_t epr[3] = {n0,n1,n2};
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

    for (TagMap::const_iterator i = tagnamestonums.begin();
            i != tagnamestonums.end(); i++) {
        setTagMap(i->first, i->second);
    }
    addPoints(points, tags);

#ifdef ESYS_MPI
    setCornerNeighbours();
#endif

#ifdef USE_RIPLEY
    coupler = NULL;
#endif
}


Brick::~Brick()
{
#ifdef USE_RIPLEY
    if (coupler != NULL)
        delete coupler;
#endif
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
    dim_t myN0, myN1, myN2;
    if (out.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NN[0];
        myN1 = m_NN[1];
        myN2 = m_NN[2];
    } else if (out.getFunctionSpace().getTypeCode() == Elements) {
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
    const dim_t first0 = max(dim_t(0), params.first[0]-m_offset[0]);
    const dim_t first1 = max(dim_t(0), params.first[1]-m_offset[1]);
    const dim_t first2 = max(dim_t(0), params.first[2]-m_offset[2]);
    // indices to first value in file (not accounting for reverse yet)
    dim_t idx0 = max(dim_t(0), m_offset[0]-params.first[0]);
    dim_t idx1 = max(dim_t(0), m_offset[1]-params.first[1]);
    dim_t idx2 = max(dim_t(0), m_offset[2]-params.first[2]);
    // number of values to read
    const dim_t num0 = min(params.numValues[0]-idx0, myN0-first0);
    const dim_t num1 = min(params.numValues[1]-idx1, myN1-first1);
    const dim_t num2 = min(params.numValues[2]-idx2, myN2-first2);

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
    const dim_t x0 = (params.reverse[0] ? num0-1 : 0);
    const int x_mult = (params.reverse[0] ? -1 : 1);
    const dim_t y0 = (params.reverse[1] ? num1-1 : 0);
    const int y_mult = (params.reverse[1] ? -1 : 1);
    const dim_t z0 = (params.reverse[2] ? num2-1 : 0);
    const int z_mult = (params.reverse[2] ? -1 : 1);

    for (index_t z=0; z<num2; z++) {
        for (index_t y=0; y<num1; y++) {
#pragma omp parallel for
            for (index_t x=0; x<num0; x++) {
                const dim_t baseIndex = first0+x*params.multiplier[0]
                                     +(first1+y*params.multiplier[1])*myN0
                                     +(first2+z*params.multiplier[2])*myN0*myN1;
                const dim_t srcIndex=(z0+z_mult*z)*num1*num0
                                  +(y0+y_mult*y)*num0
                                  +(x0+x_mult*x);
                if (!isnan(values[srcIndex])) {
                    for (index_t m2=0; m2<params.multiplier[2]; m2++) {
                        for (index_t m1=0; m1<params.multiplier[1]; m1++) {
                            for (index_t m0=0; m0<params.multiplier[0]; m0++) {
                                const dim_t dataIndex = baseIndex+m0
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
    dim_t myN0, myN1, myN2;
    if (out.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NE[0] + 1;
        myN1 = m_NE[1] + 1;
        myN2 = m_NE[2] + 1;
//    } else if (out.getFunctionSpace().getTypeCode() == Elements) {
//        myN0 = m_NE[0];
//        myN1 = m_NE[1];
//        myN2 = m_NE[2];
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
    const dim_t filesize = f.tellg();
    const dim_t reqsize = params.numValues[0]*params.numValues[1]*params.numValues[2]*numComp*sizeof(ValueType);
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
    const dim_t first0 = max(dim_t(0), params.first[0]-m_offset[0]);
    const dim_t first1 = max(dim_t(0), params.first[1]-m_offset[1]);
    const dim_t first2 = max(dim_t(0), params.first[2]-m_offset[2]);
    // indices to first value in file (not accounting for reverse yet)
    dim_t idx0 = max(dim_t(0), (m_offset[0]/params.multiplier[0])-params.first[0]);
    dim_t idx1 = max(dim_t(0), (m_offset[1]/params.multiplier[1])-params.first[1]);
    dim_t idx2 = max(dim_t(0), (m_offset[2]/params.multiplier[2])-params.first[2]);
    // if restX > 0 the first value in the respective dimension has been
    // written restX times already in a previous rank so this rank only
    // contributes (multiplier-rank) copies of that value
    const dim_t rest0 = m_offset[0]%params.multiplier[0];
    const dim_t rest1 = m_offset[1]%params.multiplier[1];
    const dim_t rest2 = m_offset[2]%params.multiplier[2];

    // number of values to read
    const dim_t num0 = min(params.numValues[0]-idx0, myN0-first0);
    const dim_t num1 = min(params.numValues[1]-idx1, myN1-first1);
    const dim_t num2 = min(params.numValues[2]-idx2, myN2-first2);

    // make sure we read the right block if going backwards through file
    if (params.reverse[2])
        idx2 = params.numValues[2]-idx2-1;

    // helpers for reversing
    const int z_mult = (params.reverse[2] ? -1 : 1);

    out.requireWrite();
    vector<ValueType> values(num0*numComp);
    const dim_t dpp = out.getNumDataPointsPerSample();

    for (dim_t z=0; z<num2; z++) {
        const dim_t m2limit = (z==0 ? params.multiplier[2]-rest2 : params.multiplier[2]);
        dim_t dataZbase = first2 + z*params.multiplier[2];
        if (z>0)
            dataZbase -= rest2;

        for (dim_t y=0; y<num1; y++) {
            const dim_t fileofs = numComp*(idx0 +
                                (idx1+y)*params.numValues[0] +
                                (idx2+z_mult*z)*params.numValues[0]*params.numValues[1]);
            f.seekg(fileofs*sizeof(ValueType));
            f.read((char*)&values[0], num0*numComp*sizeof(ValueType));
            const dim_t m1limit = (y==0 ? params.multiplier[1]-rest1 : params.multiplier[1]);
            dim_t dataYbase = first1 + y*params.multiplier[1];
            if (y>0)
                dataYbase -= rest1;

            for (dim_t x=0; x<num0; x++) {
                const dim_t m0limit = (x==0 ? params.multiplier[0]-rest0 : params.multiplier[0]);
                dim_t dataXbase = first0 + x*params.multiplier[0];
                if (x>0)
                    dataXbase -= rest0;
                // write a block of mult0 x mult1 x mult2 identical values into
                // Data object
                for (dim_t m2=0; m2 < m2limit; m2++) {
                    const dim_t dataZ = dataZbase + m2;
                    if (dataZ >= myN2)
                        break;
                    for (dim_t m1=0; m1 < m1limit; m1++) {
                        const dim_t dataY = dataYbase + m1;
                        if (dataY >= myN1)
                            break;
                        for (dim_t m0=0; m0 < m0limit; m0++) {
                            const dim_t dataX = dataXbase + m0;
                            if (dataX >= myN0)
                                break;
                            const dim_t dataIndex = INDEX3(dataX, dataY, dataZ, m_NN[0],m_NN[1]);
                            double* dest = out.getSampleDataRW(dataIndex*m_order);
                            for (int c=0; c<numComp; c++) {
                                ValueType val = values[x*numComp+c];

                                if (params.byteOrder != BYTEORDER_NATIVE) {
                                    char* cval = reinterpret_cast<char*>(&val);
                                    // this will alter val!!
                                    if (sizeof(ValueType)>4) {
                                        byte_swap64(cval);
                                    } else {
                                        byte_swap32(cval);
                                    }
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

    interpolateFromCorners(out);
}

#ifdef USE_BOOSTIO
template<typename ValueType>
void Brick::readBinaryGridZippedImpl(escript::Data& out, const string& filename,
                               const ReaderParameters& params) const
{
    // check destination function space
    dim_t myN0, myN1, myN2;
    if (out.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NE[0] + 1;
        myN1 = m_NE[1] + 1;
        myN2 = m_NE[2] + 1;
//    } else if (out.getFunctionSpace().getTypeCode() == Elements) {
//        myN0 = m_NE[0];
//        myN1 = m_NE[1];
//        myN2 = m_NE[2];
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
    dim_t filesize = f.tellg();
    f.seekg(0, ios::beg);
    std::vector<char> compressed(filesize);
    f.read((char*)&compressed[0], filesize);
    f.close();
    std::vector<char> decompressed = unzip(compressed);
    filesize = decompressed.size();
    const dim_t reqsize = params.numValues[0]*params.numValues[1]*params.numValues[2]*numComp*sizeof(ValueType);
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
    const dim_t first0 = max(dim_t(0), params.first[0]-m_offset[0]);
    const dim_t first1 = max(dim_t(0), params.first[1]-m_offset[1]);
    const dim_t first2 = max(dim_t(0), params.first[2]-m_offset[2]);
    // indices to first value in file (not accounting for reverse yet)
    dim_t idx0 = max(dim_t(0), (m_offset[0]/params.multiplier[0])-params.first[0]);
    dim_t idx1 = max(dim_t(0), (m_offset[1]/params.multiplier[1])-params.first[1]);
    dim_t idx2 = max(dim_t(0), (m_offset[2]/params.multiplier[2])-params.first[2]);
    // if restX > 0 the first value in the respective dimension has been
    // written restX times already in a previous rank so this rank only
    // contributes (multiplier-rank) copies of that value
    const dim_t rest0 = m_offset[0]%params.multiplier[0];
    const dim_t rest1 = m_offset[1]%params.multiplier[1];
    const dim_t rest2 = m_offset[2]%params.multiplier[2];

    // number of values to read
    const dim_t num0 = min(params.numValues[0]-idx0, myN0-first0);
    const dim_t num1 = min(params.numValues[1]-idx1, myN1-first1);
    const dim_t num2 = min(params.numValues[2]-idx2, myN2-first2);

    // make sure we read the right block if going backwards through file
    if (params.reverse[2])
        idx2 = params.numValues[2]-idx2-1;

    // helpers for reversing
    const int z_mult = (params.reverse[2] ? -1 : 1);

    out.requireWrite();
    vector<ValueType> values(num0*numComp);
    const int dpp = out.getNumDataPointsPerSample();

    for (dim_t z=0; z<num2; z++) {
        const dim_t m2limit = (z==0 ? params.multiplier[2]-rest2 : params.multiplier[2]);
        dim_t dataZbase = first2 + z*params.multiplier[2];
        if (z>0)
            dataZbase -= rest2;

        for (dim_t y=0; y<num1; y++) {
            const dim_t fileofs = numComp*(idx0 +
                                (idx1+y)*params.numValues[0] +
                                (idx2+z_mult*z)*params.numValues[0]*params.numValues[1]);
            memcpy((char*)&values[0], (char*)&decompressed[fileofs*sizeof(ValueType)], num0*numComp*sizeof(ValueType));
            const dim_t m1limit = (y==0 ? params.multiplier[1]-rest1 : params.multiplier[1]);
            dim_t dataYbase = first1 + y*params.multiplier[1];
            if (y>0)
                dataYbase -= rest1;

            for (dim_t x=0; x<num0; x++) {
                const dim_t m0limit = (x==0 ? params.multiplier[0]-rest0 : params.multiplier[0]);
                dim_t dataXbase = first0 + x*params.multiplier[0];
                if (x>0)
                    dataXbase -= rest0;
                // write a block of mult0 x mult1 x mult2 identical values into
                // Data object
                for (dim_t m2=0; m2 < m2limit; m2++) {
                    const dim_t dataZ = dataZbase + m2;
                    if (dataZ >= myN2)
                        break;
                    for (dim_t m1=0; m1 < m1limit; m1++) {
                        const dim_t dataY = dataYbase + m1;
                        if (dataY >= myN1)
                            break;
                        for (dim_t m0=0; m0 < m0limit; m0++) {
                            const dim_t dataX = dataXbase + m0;
                            if (dataX >= myN0)
                                break;
                            const dim_t dataIndex = INDEX3(dataX, dataY, dataZ, m_NN[0],m_NN[1]);
                            double* dest = out.getSampleDataRW(dataIndex*m_order);
                            for (int c=0; c<numComp; c++) {
                                ValueType val = values[x*numComp+c];

                                if (params.byteOrder != BYTEORDER_NATIVE) {
                                    char* cval = reinterpret_cast<char*>(&val);
                                    // this will alter val!!
                                    if (sizeof(ValueType)>4) {
                                        byte_swap64(cval);
                                    } else {
                                        byte_swap32(cval);
                                    }
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
    interpolateFromCorners(out);
}
#endif

void Brick::interpolateFromCorners(escript::Data &out) const
{
    const int numComp = out.getDataPointSize();
    //interpolate the missing portions
#pragma omp parallel for
    for (dim_t z = 0; z < m_NN[2]; z++) {
        const double pz = point_locations[m_order-2][z%m_order];
        for (dim_t y = 0; y < m_NN[1]; y++) {
            const double py = point_locations[m_order-2][y%m_order];
            for (dim_t x = 0; x < m_NN[0]; x++) {
                //skip the points we have values for
                if (y % m_order == 0 && x % m_order == 0 && z % m_order == 0)
                    continue;
                //point location in element: x,y
                const double px = point_locations[m_order-2][x%m_order];

                //the point we're interpolating a value for
                double *point = out.getSampleDataRW(
                        INDEX3(x, y, z, m_NN[0], m_NN[1]));

                const dim_t left = x - x%m_order;
                const dim_t right = left < m_NN[0] - 1 ? left + m_order : left;
                const dim_t front = y - y%m_order;
                const dim_t back = front < m_NN[1] - 1 ? front + m_order : front;
                const dim_t down = z - z%m_order;
                const dim_t up = down < m_NN[2] - 1 ? down + m_order : down;
         
                //corner values
                const double *dlf = out.getSampleDataRO(
                        INDEX3(left, front, down, m_NN[0], m_NN[1]));
                const double *dlb = out.getSampleDataRO(
                        INDEX3(left, back, down, m_NN[0], m_NN[1]));
                const double *drf = out.getSampleDataRO(
                        INDEX3(right, front, down, m_NN[0], m_NN[1]));
                const double *drb = out.getSampleDataRO(
                        INDEX3(right, back, down, m_NN[0], m_NN[1]));
                const double *ulf = out.getSampleDataRO(
                        INDEX3(left, front, up, m_NN[0], m_NN[1]));
                const double *ulb = out.getSampleDataRO(
                        INDEX3(left, back, up, m_NN[0], m_NN[1]));
                const double *urf = out.getSampleDataRO(
                        INDEX3(right, front, up, m_NN[0], m_NN[1]));
                const double *urb = out.getSampleDataRO(
                        INDEX3(right, back, up, m_NN[0], m_NN[1]));

                //the interpolation itself
                for (int comp = 0; comp < numComp; comp++) {
                    point[comp] = urb[comp]*px    *py    *pz
                                + drb[comp]*px    *py    *(1-pz)
                                + urf[comp]*px    *(1-py)*pz
                                + drf[comp]*px    *(1-py)*(1-pz)
                                + ulb[comp]*(1-px)*py    *pz
                                + dlb[comp]*(1-px)*py    *(1-pz)
                                + ulf[comp]*(1-px)*(1-py)*pz
                                + dlf[comp]*(1-px)*(1-py)*(1-pz);
                }
            }
        }
    }
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
            throw SpeckleyException("writeBinaryGrid(): invalid or unsupported datatype");
    }
}

template<typename ValueType>
void Brick::writeBinaryGridImpl(const escript::Data& in,
                                const string& filename, int byteOrder) const
{
    // check function space and determine number of points
    dim_t myN0, myN1, myN2;
    dim_t totalN0, totalN1, totalN2;
    if (in.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NE[0] + 1;
        myN1 = m_NE[1] + 1;
        myN2 = m_NE[2] + 1;
        totalN0 = m_gNE[0]+1;
        totalN1 = m_gNE[1]+1;
        totalN2 = m_gNE[2]+1;
    } else if (in.getFunctionSpace().getTypeCode() == Elements) {
        myN0 = m_NE[0];
        myN1 = m_NE[1];
        myN2 = m_NE[2];
        totalN0 = m_gNE[0];
        totalN1 = m_gNE[1];
        totalN2 = m_gNE[2];
    } else
        throw SpeckleyException("writeBinaryGrid(): invalid function space of data object");

    const int numComp = in.getDataPointSize();
    const dim_t dpp = in.getNumDataPointsPerSample();
    const dim_t fileSize = sizeof(ValueType)*numComp*dpp*totalN0*totalN1*totalN2;

    if (numComp > 1 || dpp > 1)
        throw SpeckleyException("writeBinaryGrid(): only scalar, single-value data supported");

    // from here on we know that each sample consists of one value
    FileWriter fw;
    fw.openFile(filename, fileSize);
    MPIBarrier();

    for (index_t z=0; z<myN2; z++) {
        for (index_t y=0; y<myN1; y++) {
            const dim_t fileofs = (m_offset[0]+(m_offset[1]+y)*totalN0
                                +(m_offset[2]+z)*totalN0*totalN1)*sizeof(ValueType);
            ostringstream oss;

            for (index_t x=0; x<myN0; x++) {
                const double* sample = in.getSampleDataRO(
                                INDEX3(x,y,z,m_NN[0],m_NN[1])*m_order);
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

void Brick::write(const std::string& filename) const
{
    throw SpeckleyException("write: not supported");
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
    const dim_t NN0 = m_NN[0];
    const dim_t NN1 = m_NN[1];
    const dim_t NN2 = m_NN[2];

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
    std::vector<int> dims(m_NN, m_NN+3);

    // write mesh
    DBPutQuadmesh(dbfile, "mesh", NULL, coords, &dims[0], 3, DB_DOUBLE,
            DB_COLLINEAR, NULL);

    // write node ids
    DBPutQuadvar1(dbfile, "nodeId", "mesh", (void*)&m_nodeId[0], &dims[0], 3,
            NULL, 0, DB_INT, DB_NODECENT, NULL);

    // write element ids
    dims.assign(m_NE, m_NE+3);
    DBPutQuadvar1(dbfile, "elementId", "mesh", (void*)&m_elementId[0],
            &dims[0], 3, NULL, 0, DB_INT, DB_ZONECENT, NULL);

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

const dim_t* Brick::borrowSampleReferenceIDs(int fsType) const
{
    switch (fsType) {
        case DegreesOfFreedom:
        case Nodes:
            return &m_nodeId[0];
        case Elements:
            return &m_elementId[0];
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
    if (out.getFunctionSpace().getTypeCode() == Elements) {
        out.requireWrite();
        const dim_t numQuad = out.getNumDataPointsPerSample();
        const double size=sqrt(m_dx[0]*m_dx[0]+m_dx[1]*m_dx[1]+m_dx[2]*m_dx[2]);
        const dim_t numElements = getNumElements();
#pragma omp parallel for
        for (index_t k = 0; k < numElements; ++k) {
            double* o = out.getSampleDataRW(k);
            std::fill(o, o+numQuad, size);
        }
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
    escriptDataC dc = arg.getDataC();
    int numDim = m_numDim;
    if (!isDataPointShapeEqual(&dc, 1, &numDim))
        throw SpeckleyException("setToX: Invalid Data object shape");
    if (!numSamplesEqual(&dc, 1, getNumNodes()))
        throw SpeckleyException("setToX: Illegal number of samples in Data object");

    const dim_t NN0 = m_NN[0];
    const dim_t NN1 = m_NN[1];
    const dim_t NN2 = m_NN[2];
    arg.requireWrite();
#pragma omp parallel for
    for (dim_t i2 = 0; i2 < NN2; i2++) {
        for (dim_t i1 = 0; i1 < NN1; i1++) {
            for (dim_t i0 = 0; i0 < NN0; i0++) {
                double* point = arg.getSampleDataRW(i0+NN0*i1+NN0*NN1*i2);
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

#define RANK_LEFT(__rank__) ((__rank__) % m_NX[0] == 0 ? 0 : 1)
#define RANK_FRONT(__rank__) ((__rank__) % (m_NX[0]*m_NX[1])/m_NX[0] == 0 ? 0 : 1)
#define RANK_BELOW(__rank__) ((__rank__) / (m_NX[0]*m_NX[1]) == 0 ? 0 : 1)

    m_nodeDistribution.assign(m_mpiInfo->size+1, 0);

    for (dim_t k = 0; k < m_mpiInfo->size - 1; k++) {
        m_nodeDistribution[k+1] = m_nodeDistribution[k]
                                + (m_NN[0]-RANK_LEFT(k))
                                * (m_NN[1]-RANK_FRONT(k))
                                * (m_NN[2]-RANK_BELOW(k));
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

    const int rank = m_mpiInfo->rank;

    const index_t left = RANK_LEFT(rank);
    const index_t front = RANK_FRONT(rank);
    const index_t bottom = RANK_BELOW(rank);


    if (left && front) {
        //re-use bottom-left-front corner node
        if (bottom) {
            //get lower-left-front node
            int rank_wanted = rank
                        - m_NX[0]*m_NX[1]   //down a layer
                        - m_NX[0];          //towards the front
            m_nodeId[0] = m_nodeDistribution[rank_wanted] - 1; //end of the left
        }
        //front left edge
        int neighbour = rank - m_NX[0] - 1;
        int neighboursLeft = RANK_LEFT(neighbour);
        int neighboursFront = RANK_FRONT(neighbour);
        index_t begin = m_nodeDistribution[neighbour] 
                    + (m_NN[0]-neighboursLeft)*(m_NN[1]-neighboursFront) - 1;
#pragma omp parallel for
        for (index_t z = bottom; z < m_NN[2]; z++) {
            index_t theirs = z*(m_NN[0]-neighboursLeft)*(m_NN[1]-neighboursFront);
            index_t mine = z*m_NN[0]*m_NN[1];
            m_nodeId[mine] = begin + theirs;
        }
    }
    //re-use nodes on bottom border
    if (bottom) {
        int neighbour = rank - m_NX[0]*m_NX[1];
        //beginning, top left front of rank underneath
        index_t begin = m_nodeDistribution[neighbour + 1] - m_NN[0]*m_NN[1];
#pragma omp parallel for
        for (index_t y = front; y < m_NN[1]; y++) {
            for (index_t x = left; x < m_NN[0]; x++) {
                index_t i = INDEX2(x,y,m_NN[0]);
                m_nodeId[i] = begin + INDEX2(x,y,m_NN[0]);
            }
        }
    }

    //re-use nodes on front border
    if (front) {
        int neighbour = rank - m_NX[0];
        index_t begin = m_nodeDistribution[neighbour] 
                        + (m_NN[0]-RANK_LEFT(neighbour))
                        *(m_NN[1]-RANK_FRONT(neighbour) - 1);
#pragma omp parallel for
        for (index_t z = bottom; z < m_NN[2]; z++) {
            index_t mine = z*m_NN[0]*m_NN[1] + left;
            index_t theirs = z*(m_NN[0]-RANK_LEFT(neighbour))*(m_NN[1]-RANK_FRONT(neighbour));
            for (index_t x = left; x < m_NN[0]; x++, theirs++, mine++) {
                m_nodeId[mine] = begin + theirs;
            }
        }
    }
    //re-use nodes on left border
    if (left) {
        //is the rank to the left itself right of another rank
        index_t neighboursLeft = RANK_LEFT(rank - 1);
        index_t neighboursFront = RANK_FRONT(rank - 1);
        index_t neighboursBelow = RANK_BELOW(rank - 1);
        //end of first owned row of neighbouring rank
        index_t end = m_nodeDistribution[rank - 1] + m_NN[0]-neighboursLeft - 1;
#pragma omp parallel for
        for (index_t z = bottom; z < m_NN[2]; z++) {
            for (index_t y = front; y < m_NN[1]; y++) {
                index_t mine = INDEX3(0,y,z,m_NN[0],m_NN[1]);
                index_t theirs = INDEX3(0,(y-neighboursFront),(z-neighboursBelow),(m_NN[0]-neighboursLeft),(m_NN[1]-neighboursFront));
                m_nodeId[mine] = end + theirs;
            }
        }
    }

    //now the new nodes
    const index_t start = m_nodeDistribution[m_mpiInfo->rank];
#pragma omp parallel for
    for (index_t z = bottom; z < m_NN[2]; z++) {
        index_t z_chunk = (z-bottom)*(m_NN[0]-left)*(m_NN[1]-front);
        for (index_t y = front; y < m_NN[1]; y++) {
            index_t y_chunk = (y-front)*(m_NN[0]-left);
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
#undef RANK_TO_LEFT
#undef RANK_TO_RIGHT
#undef RANK_TO_FRONT
}

//private
void Brick::addToMatrixAndRHS(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<double>& EM_S, const std::vector<double>& EM_F, bool addS,
         bool addF, index_t firstNode, int nEq, int nComp) const
{
    throw SpeckleyException("Rectangle::addToMatrixAndRHS not implemented");
}

void Brick::interpolateElementsOnNodes(escript::Data& out,
                                  const escript::Data& in) const {
    const dim_t numComp = in.getDataPointSize();
    const dim_t NE0 = m_NE[0];
    const dim_t NE1 = m_NE[1];
    const dim_t NE2 = m_NE[2];
    const int quads = m_order + 1;
    const dim_t max_x = m_NN[0];
    const dim_t max_y = m_NN[1];
    const dim_t max_z = m_NN[2];
    out.requireWrite();
    //init to zero so we can do some sums without undefined, may not be required
    memset(out.getSampleDataRW(0), 0, sizeof(double)*quads*quads*numComp);
    // the summation portion
    for (dim_t colouring = 0; colouring < 2; colouring++) {
#pragma omp parallel for
        for (dim_t ez = colouring; ez < NE2; ez += 2) {
            for (dim_t ey = 0; ey < NE1; ey++) {
                for (dim_t ex = 0; ex < NE0; ex++) {
                    dim_t start = m_order * (INDEX3(ex, ey, ez, max_x, max_y));
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
#ifdef ESYS_MPI
    //sum and average neighbours before we average out our internal structure
    balanceNeighbours(out, true);
#endif

    /* the averaging out (each point divided by the number of additions in
       the summation step). By doing each edge along an axis, those points
       requiring division by 2, 4 and 8 will be divided by 2 the right
       number of times

       border edges are skipped because they aren't shared but for the
       points that lie along a different axes non-border edge
    */
    // for every non-border edge in x
#pragma omp parallel for
    for (index_t qz = 0; qz < max_z; qz++) {
        for (index_t qy = 0; qy < max_y; qy++) {
            for (index_t qx = m_order; qx < max_x - m_order; qx += m_order) {
                double *n_out = out.getSampleDataRW(INDEX3(qx, qy, qz, max_x, max_y));
                for (int comp = 0; comp < numComp; comp++) {
                    n_out[comp] /= 2.;
                }
            }
        }
    }
    // for every non-border edge in y
#pragma omp parallel for
    for (index_t qz = 0; qz < max_z; qz++) {
        for (index_t qy = m_order; qy < max_y - m_order; qy += m_order) {
            for (index_t qx = 0; qx < max_x; qx ++) {
                double *n_out = out.getSampleDataRW(INDEX3(qx, qy, qz, max_x, max_y));
                for (int comp = 0; comp < numComp; comp++) {
                    n_out[comp] /= 2.;
                }
            }
        }
    }
    // for every non-border edge in z
#pragma omp parallel for
    for (index_t qz = m_order; qz < max_z - m_order; qz += m_order) {
        for (index_t qy = 0; qy < max_y; qy++) {
            for (index_t qx = 0; qx < max_x; qx++) {
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
                                       const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    const dim_t NE0 = m_NE[0];
    const dim_t NE1 = m_NE[1];
    const dim_t NE2 = m_NE[2];
    const int quads = m_order + 1;
    const dim_t max_x = m_NN[0];
    const dim_t max_y = m_NN[1];
    out.requireWrite();
#pragma omp parallel for
    for (dim_t ez = 0; ez < NE2; ez++) {
        for (dim_t ey = 0; ey < NE1; ey++) {
            for (dim_t ex = 0; ex < NE0; ex++) {
                double *e_out = out.getSampleDataRW(INDEX3(ex, ey, ez, NE0, NE1));
                dim_t start = m_order * INDEX3(ex, ey, ez, max_x, max_y);
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

#ifdef ESYS_MPI

//protected
void Brick::balanceNeighbours(escript::Data& out, bool average) const
{
    // skip all this if we aren't even subdividing the domain
    if (m_NX[0] * m_NX[1] * m_NX[2] == 1) {
        return;
    }
    const int numComp = out.getDataPointSize();
    int rx = m_mpiInfo->rank % m_NX[0];
    int ry = m_mpiInfo->rank % (m_NX[0]*m_NX[1]) / m_NX[0];
    int rz = m_mpiInfo->rank / (m_NX[0]*m_NX[1]);
    //include bordering ranks in summation
    //averaging waits til after all sharing
    shareFaces(out, rx, ry, rz);
    if ((m_NX[0] > 1 && m_NX[1] > 1)
            || (m_NX[1] > 1 && m_NX[2] > 1)
            || (m_NX[0] > 1 && m_NX[2] > 1))
        shareEdges(out, rx, ry, rz);
    if (m_NX[0] != 1 && m_NX[1] != 1 && m_NX[2] != 1) {
        shareCorners(out);
        if (!average)
            return;
        //averaging out corners now that all sharing done
        for (int z = 0; z < 2; z++) {
            for (int y = 0; y < 2; y++) {
                for (int x = 0; x < 2; x++) {
                    if (!neighbour_exists[INDEX3(x,y,z,2,2)])
                        continue;
                    double *values = out.getSampleDataRW(
                            INDEX3(x*(m_NN[0]-1), y*(m_NN[1]-1), z*(m_NN[2]-1),
                                            m_NN[0], m_NN[1]));
                    for (int comp = 0; comp < numComp; comp++) {
                        values[comp] /= 2;
                    }
                }
            }
        }
    }
    if (!average)
        return;
    // average shared-edges

    const bool left = rx;
    const bool right = rx < m_NX[0] - 1;
    const bool front = ry;
    const bool back = ry < m_NX[1] - 1;
    const bool bottom = rz;
    const bool top = rz < m_NX[2] - 1;
    if (left) {
        if (front) { //average Z lines
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[2]; i++) {
                double *data = out.getSampleDataRW(
                                INDEX3(0, 0, i, m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;
                }
            }
        }
        if (back) {//average Z lines
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[2]; i++) {
                double *data = out.getSampleDataRW(
                                INDEX3(0, m_NN[1]-1, i, m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;
                }
            }
        }
        if (top) {//average Y lines
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[1]; i++) {
                double *data = out.getSampleDataRW(
                                INDEX3(0, i, m_NN[2]-1, m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;;
                }
            }
        }
        if (bottom) {//average Y lines
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[1]; i++) {
                double *data = out.getSampleDataRW(
                                INDEX3(0,i,0,m_NN[0],m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;
                }
            }
        }
    }
    if (right) {
        if (front) { //average Z lines
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[2]; i++) {
                double *data = out.getSampleDataRW(
                                INDEX3(m_NN[0]-1, 0, i, m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;
                }
            }
        }
        if (back) {//average Z lines
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[2]; i++) {
                double *data = out.getSampleDataRW(
                                INDEX3(m_NN[0]-1,m_NN[1]-1,i,m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;
                }
            }
        }
        if (top) {//average Y lines
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[1]; i++) {
                double *data = out.getSampleDataRW(
                                INDEX3(m_NN[0]-1,i,m_NN[2]-1, m_NN[0],m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;
                }
            }
        }
        if (bottom) {//average Y lines
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[1]; i++) {
                double *data = out.getSampleDataRW(
                                INDEX3(m_NN[0]-1, i, 0, m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;
                }
            }
        }
    }

    if (top) {
        if (front) {//average X lines
#pragma omp parallel for
            for (dim_t x = 0; x < m_NN[0]; x++) {
                double *data = out.getSampleDataRW(
                                INDEX3(x, 0, m_NN[2]-1, m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;
                }
            }
        }
        if (back) {//average X lines
#pragma omp parallel for
            for (dim_t x = 0; x < m_NN[0]; x++) {
                double *data = out.getSampleDataRW(
                                INDEX3(x,m_NN[1]-1,m_NN[2]-1,m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;
                }
            }
        }
    }
    if (bottom) {
        if (front) {//average X lines
#pragma omp parallel for
            for (dim_t x = 0; x < m_NN[0]; x++) {
                double *data = out.getSampleDataRW(x);
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;
                }
            }
        }
        if (back) {//average X lines
#pragma omp parallel for
            for (dim_t x = 0; x < m_NN[0]; x++) {
                double *data = out.getSampleDataRW(
                                INDEX3(x, m_NN[1]-1, 0, m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;
                }
            }
        }
    }

    //average shared faces
    //up and down
    if (top) {
#pragma omp parallel for
        for (dim_t y = 0; y <m_NN[1]; y++) {
            for (dim_t x = 0; x <m_NN[0]; x++) {
                double *data = out.getSampleDataRW(
                                INDEX3(x, y, m_NN[2]-1, m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;
                }
            }
        }
    }
    if (bottom) {
#pragma omp parallel for
        for (dim_t y = 0; y <m_NN[1]; y++) {
            for (dim_t x = 0; x <m_NN[0]; x++) {
                double *data = out.getSampleDataRW(INDEX2(x, y, m_NN[0]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;
                }
            }
        }
    }
    //left and right
    if (right) {
#pragma omp parallel for
        for (dim_t z = 0; z <m_NN[2]; z++) {
            for (dim_t y = 0; y <m_NN[1]; y++) {
                double *data = out.getSampleDataRW(
                                INDEX3(m_NN[0]-1, y, z, m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;
                }
            }
        }
    }
    if (left) {
#pragma omp parallel for
        for (dim_t z = 0; z <m_NN[2]; z++) {
            for (dim_t y = 0; y <m_NN[1]; y++) {
                double *data = out.getSampleDataRW(
                                INDEX3(0, y, z, m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;
                }
            }
        }
    }
    //front and back
    if (back) {
#pragma omp parallel for
        for (dim_t z = 0; z <m_NN[2]; z++) {
            for (dim_t x = 0; x <m_NN[0]; x++) {
                double *data = out.getSampleDataRW(
                                INDEX3(x, m_NN[1]-1, z, m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;
                }
            }
        }
    }
    if (front) {
#pragma omp parallel for
        for (dim_t z = 0; z <m_NN[2]; z++) {
            for (dim_t x = 0; x <m_NN[0]; x++) {
                double *data = out.getSampleDataRW(
                                INDEX3(x, 0, z, m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] /= 2;
                }
            }
        }
    }
}

//private
void Brick::setCornerNeighbours()
{
    const int rank = m_mpiInfo->rank;

    const int rx = rank % m_NX[0];
    const int ry = rank % (m_NX[0]*m_NX[1])/m_NX[0];
    const int rz = rank / (m_NX[0]*m_NX[1]);

    const bool left = rx;
    const bool right = rx < m_NX[0] - 1;
    const bool front = ry;
    const bool back = ry < m_NX[1] - 1;
    const bool bottom = rz;
    const bool top = rz < m_NX[2] - 1;

    neighbour_exists[0] = left && front && bottom;
    neighbour_exists[1] = right && front && bottom;
    neighbour_exists[2] = left && back && bottom;
    neighbour_exists[3] = right && back && bottom;
    neighbour_exists[4] = left && front && top;
    neighbour_exists[5] = right && front && top;
    neighbour_exists[6] = left && back && top;
    neighbour_exists[7] = right && back && top;

    neighbour_ranks[0] = rank - m_NX[0]*m_NX[1] - m_NX[0] - 1;
    neighbour_ranks[1] = rank - m_NX[0]*m_NX[1] - m_NX[0] + 1;
    neighbour_ranks[2] = rank - m_NX[0]*m_NX[1] + m_NX[0] - 1;
    neighbour_ranks[3] = rank - m_NX[0]*m_NX[1] + m_NX[0] + 1;
    neighbour_ranks[4] = rank + m_NX[0]*m_NX[1] - m_NX[0] - 1;
    neighbour_ranks[5] = rank + m_NX[0]*m_NX[1] - m_NX[0] + 1;
    neighbour_ranks[6] = rank + m_NX[0]*m_NX[1] + m_NX[0] - 1;
    neighbour_ranks[7] = rank + m_NX[0]*m_NX[1] + m_NX[0] + 1;
}

//private
void Brick::shareCorners(escript::Data& out) const
{
    //setup
    const int tag = 0;
    MPI_Status status;
    const int numComp = out.getDataPointSize();
    const int count = numComp;
    std::vector<double> inbuf(count, 0);

    //share
    for (int z = 0; z < 2; z++) {
        for (int y = 0; y < 2; y++) {
            for (int x = 0; x < 2; x++) {
                int i = INDEX3(x,y,z,2,2);
                if (neighbour_exists[i]) {
                    double *data = out.getSampleDataRW(
                                           x*(m_NN[0]-1)
                                         + y*(m_NN[1]-1)*m_NN[0]
                                         + z*(m_NN[2]-1)*m_NN[0]*m_NN[1]
                                        );

                    MPI_Sendrecv(data, numComp, MPI_DOUBLE, neighbour_ranks[i], tag,
                            &inbuf[0], numComp, MPI_DOUBLE, neighbour_ranks[i],
                            tag, m_mpiInfo->comm, &status);
                    //unpack
                    for (int comp = 0; comp < numComp; comp++) {
                        data[comp] += inbuf[comp];
                    }
                }
            }
        }
    }
}

//private
void Brick::shareEdges(escript::Data& out, int rx, int ry, int rz) const
{
    const int rank = m_mpiInfo->rank;

    const int tag = 0;
    MPI_Status status;
    const int numComp = out.getDataPointSize();

    const bool left = rx;
    const bool right = rx < m_NX[0] - 1;
    const bool front = ry;
    const bool back = ry < m_NX[1] - 1;
    const bool bottom = rz;
    const bool top = rz < m_NX[2] - 1;


    if (left) {
        if (front) { //share Z lines
            int neighbour = rank - m_NX[0] - 1;
            const double count = m_NN[2]*numComp;
            std::vector<double> buf(count);
            std::vector<double> outbuf(count);
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[2]; i++) {
                double *data = out.getSampleDataRW(INDEX3(0, 0, i, m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    outbuf[i*numComp + comp] = data[comp];
                }
            }
            MPI_Sendrecv(&outbuf[0], count, MPI_DOUBLE, neighbour, tag,
                            &buf[0], count, MPI_DOUBLE, neighbour, tag,
                            m_mpiInfo->comm, &status);
            for (dim_t i = 0; i < m_NN[2]; i++) {
                double *data = out.getSampleDataRW(INDEX3(0, 0, i, m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += buf[i*numComp + comp];
                }
            }
        }
        if (back) {//share Z lines
            int neighbour = rank + m_NX[0] - 1;
            const double count = m_NN[2]*numComp;
            std::vector<double> buf(count);
            std::vector<double> outbuf(count);
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[2]; i++) {
                double *data = out.getSampleDataRW(INDEX3(0,m_NN[1]-1,i,m_NN[0],m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    outbuf[i*numComp + comp] = data[comp];
                }
            }
            MPI_Sendrecv(&outbuf[0], count, MPI_DOUBLE, neighbour, tag,
                            &buf[0], count, MPI_DOUBLE, neighbour, tag,
                            m_mpiInfo->comm, &status);
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[2]; i++) {
                double *data = out.getSampleDataRW(INDEX3(0,m_NN[1]-1,i,m_NN[0],m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += buf[i*numComp + comp];
                }
            }
        }
        if (top) {//share Y lines
            int neighbour = rank + m_NX[0]*m_NX[1] - 1;
            const double count = m_NN[1]*numComp;
            std::vector<double> buf(count);
            std::vector<double> outbuf(count);
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[1]; i++) {
                double *data = out.getSampleDataRW(INDEX3(0,i,m_NN[2]-1,m_NN[0],m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    outbuf[i*numComp + comp] = data[comp];
                }
            }
            MPI_Sendrecv(&outbuf[0], count, MPI_DOUBLE, neighbour, tag,
                            &buf[0], count, MPI_DOUBLE, neighbour, tag,
                            m_mpiInfo->comm, &status);
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[1]; i++) {
                double *data = out.getSampleDataRW(INDEX3(0,i,m_NN[2]-1,m_NN[0],m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += buf[i*numComp + comp];
                }
            }
        }
        if (bottom) {//share Y lines
            int neighbour = rank - m_NX[0]*m_NX[1] - 1;
            const double count = m_NN[1]*numComp;
            std::vector<double> buf(count);
            std::vector<double> outbuf(count);
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[1]; i++) {
                double *data = out.getSampleDataRW(INDEX3(0,i,0,m_NN[0],m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    outbuf[i*numComp + comp] = data[comp];
                }
            }
            MPI_Sendrecv(&outbuf[0], count, MPI_DOUBLE, neighbour, tag,
                            &buf[0], count, MPI_DOUBLE, neighbour, tag,
                            m_mpiInfo->comm, &status);
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[1]; i++) {
                double *data = out.getSampleDataRW(INDEX3(0,i,0,m_NN[0],m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += buf[i*numComp + comp];
                }
            }
        }
    }
    if (right) {
        if (front) { //share Z lines
            int neighbour = rank - m_NX[0] + 1;
            const double count = m_NN[2]*numComp;
            std::vector<double> buf(count);
            std::vector<double> outbuf(count);
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[2]; i++) {
                double *data = out.getSampleDataRW(
                            INDEX3(m_NN[0]-1, 0, i, m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    outbuf[i*numComp + comp] = data[comp];
                }
            }
            MPI_Sendrecv(&outbuf[0], count, MPI_DOUBLE, neighbour, tag,
                            &buf[0], count, MPI_DOUBLE, neighbour, tag,
                            m_mpiInfo->comm, &status);
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[2]; i++) {
                double *data = out.getSampleDataRW(
                            INDEX3(m_NN[0]-1, 0, i, m_NN[0], m_NN[1]));     
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += buf[i*numComp + comp];
                }
            }
        }
        if (back) {//share Z lines
            int neighbour = rank + m_NX[0] + 1;
            const double count = m_NN[2]*numComp;
            std::vector<double> buf(count);
            std::vector<double> outbuf(count);
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[2]; i++) {
                double *data = out.getSampleDataRW(
                            INDEX3(m_NN[0]-1, m_NN[1]-1, i, m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    outbuf[i*numComp + comp] = data[comp];
                }
            }
            MPI_Sendrecv(&outbuf[0], count, MPI_DOUBLE, neighbour, tag,
                            &buf[0], count, MPI_DOUBLE, neighbour, tag,
                            m_mpiInfo->comm, &status);
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[2]; i++) {
                double *data = out.getSampleDataRW(
                            INDEX3(m_NN[0]-1, m_NN[1]-1, i, m_NN[0], m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += buf[i*numComp + comp];
                }
            }
        }
        if (top) {//share Y lines
            int neighbour = rank + m_NX[0]*m_NX[1] + 1;
            const double count = m_NN[1]*numComp;
            std::vector<double> buf(count);
            std::vector<double> outbuf(count);
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[1]; i++) {
                double *data = out.getSampleDataRW(INDEX3(m_NN[0]-1,i,m_NN[2]-1,m_NN[0],m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    outbuf[i*numComp + comp] = data[comp];
                }
            }
            MPI_Sendrecv(&outbuf[0], count, MPI_DOUBLE, neighbour, tag,
                            &buf[0], count, MPI_DOUBLE, neighbour, tag,
                            m_mpiInfo->comm, &status);
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[1]; i++) {
                double *data = out.getSampleDataRW(INDEX3(m_NN[0]-1,i,m_NN[2]-1,m_NN[0],m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += buf[i*numComp + comp];
                }
            }
        }
        if (bottom) {//share Y lines
            int neighbour = rank - m_NX[0]*m_NX[1] + 1;
            const double count = m_NN[1]*numComp;
            std::vector<double> buf(count);
            std::vector<double> outbuf(count);
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[1]; i++) {
                double *data = out.getSampleDataRW(INDEX3(m_NN[0]-1,i,0,m_NN[0],m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    outbuf[i*numComp + comp] = data[comp];
                }
            }
            MPI_Sendrecv(&outbuf[0], count, MPI_DOUBLE, neighbour, tag,
                            &buf[0], count, MPI_DOUBLE, neighbour, tag,
                            m_mpiInfo->comm, &status);
#pragma omp parallel for
            for (dim_t i = 0; i < m_NN[1]; i++) {
                double *data = out.getSampleDataRW(INDEX3(m_NN[0]-1,i,0,m_NN[0],m_NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += buf[i*numComp + comp];
                }
            }
        }
    }

    if (top) {
        const double count = m_NN[0]*numComp;
        std::vector<double> buf(count);
        if (front) {//share X lines
            int neighbour = rank + m_NX[0]*m_NX[1] - m_NX[0];
            double *data = out.getSampleDataRW(m_NN[0]*m_NN[1]*(m_NN[2]-1));
                MPI_Sendrecv(data, count, MPI_DOUBLE, neighbour, tag,
                                &buf[0], count, MPI_DOUBLE, neighbour, tag,
                                m_mpiInfo->comm, &status);
            for (dim_t i = 0; i < count; i++) {
                data[i] += buf[i];
            }
        }
        if (back) {//share X lines
            int neighbour = rank + m_NX[0]*m_NX[1] + m_NX[0];
            double *data = out.getSampleDataRW(m_NN[0]*m_NN[1]*(m_NN[2]-1)
                                             + m_NN[0]*(m_NN[1]-1));
            MPI_Sendrecv(data, count, MPI_DOUBLE, neighbour, tag,
                            &buf[0], count, MPI_DOUBLE, neighbour, tag,
                            m_mpiInfo->comm, &status);
            for (dim_t i = 0; i < count; i++) {
                data[i] += buf[i];
            }
        }
    }
    if (bottom) {
        const double count = m_NN[0]*numComp;
        std::vector<double> buf(count);
        if (front) {//share X lines
            int neighbour = rank - m_NX[0]*m_NX[1] - m_NX[0];
            double *data = out.getSampleDataRW(0);
            MPI_Sendrecv(data, count, MPI_DOUBLE, neighbour, tag,
                            &buf[0], count, MPI_DOUBLE, neighbour, tag,
                            m_mpiInfo->comm, &status);
            for (dim_t i = 0; i < count; i++) {
                data[i] += buf[i];
            }
        }
        if (back) {//share X lines
            int neighbour = rank - m_NX[0]*m_NX[1] + m_NX[0];
            double *data = out.getSampleDataRW(m_NN[0]*(m_NN[1]-1));
            MPI_Sendrecv(data, count, MPI_DOUBLE, neighbour, tag,
                            &buf[0], count, MPI_DOUBLE, neighbour, tag,
                            m_mpiInfo->comm, &status);
            for (dim_t i = 0; i < count; i++) {
                data[i] += buf[i];
            }
        }
    }
}


void frontAndBack(escript::Data& out, int ry, const int numComp, int rank,
                    const dim_t NN[3], const int NX[3], MPI_Comm& comm) {
    const int tag = 0;
    MPI_Status status;
    const int front_neighbour = rank - NX[0];
    const int back_neighbour = rank + NX[0];
    const dim_t count = NN[0]*NN[2]*numComp;
    std::vector<double> front(count);
    std::vector<double> back(count);
    std::vector<double> recv(count);
#pragma omp parallel for
    for (dim_t z = 0; z < NN[2]; z++) {
        for (dim_t x = 0; x < NN[0]; x++) {
            const dim_t index = INDEX2(x,z,NN[0])*numComp;
            const double *frontData = out.getSampleDataRO(INDEX3(x, 0, z, NN[0], NN[1]));
            std::copy(frontData, frontData + numComp, &front[index]);

            const double *backData = out.getSampleDataRO(INDEX3(x, NN[1]-1, z, NN[0], NN[1]));
            std::copy(backData, backData + numComp, &back[index]);
        }
    }
    //front
    if (ry) {
        MPI_Sendrecv(&front[0], count, MPI_DOUBLE, front_neighbour, tag,
                    &recv[0], count, MPI_DOUBLE, front_neighbour, tag,
                    comm, &status);
        //unpack front
#pragma omp parallel for
        for (dim_t z = 0; z < NN[2]; z++) {
            for (dim_t x = 0; x < NN[0]; x++) {
                const dim_t index = INDEX2(x,z,NN[0])*numComp;
                double *data = out.getSampleDataRW(INDEX3(x, 0, z, NN[0], NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += recv[index+comp];
                }
            }
        }
    }

    //back
    if (ry < NX[1] - 1) {
        MPI_Sendrecv(&back[0], count, MPI_DOUBLE, back_neighbour, tag,
                &recv[0], count, MPI_DOUBLE, back_neighbour, tag,
                comm, &status);
        //unpack back
#pragma omp parallel for
        for (dim_t z = 0; z < NN[2]; z++) {
            for (dim_t x = 0; x < NN[0]; x++) {
                const dim_t index = INDEX2(x,z,NN[0])*numComp;
                double *data = out.getSampleDataRW(INDEX3(x, NN[1] - 1, z, NN[0], NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += recv[index+comp];
                }
            }
        }
    }
}

void topAndBottom(escript::Data& out, int rz, int numComp, int rank,
                    const dim_t NN[3], const int NX[3], MPI_Comm& comm) {
    const int tag = 0;
    MPI_Status status;
    const int top_neighbour = rank + NX[0]*NX[1];
    const int bottom_neighbour = rank - NX[0]*NX[1];
    const dim_t count = NN[0]*NN[1]*numComp;
    std::vector<double> top(count);
    std::vector<double> bottom(count);
    std::vector<double> recv(count);
#pragma omp parallel for
    for (dim_t y = 0; y < NN[1]; y++) {
        for (dim_t x = 0; x < NN[0]; x++) {
            const dim_t index = INDEX2(x,y,NN[0])*numComp;
            const double *bottomData = out.getSampleDataRO(INDEX2(x, y, NN[0]));
            std::copy(bottomData, bottomData + numComp, &bottom[index]);

            const double *topData = out.getSampleDataRO(INDEX3(x, y, NN[2]-1, NN[0], NN[1]));
            std::copy(topData, topData + numComp, &top[index]);
        }
    }

    //bottom
    if (rz) {
        MPI_Sendrecv(&bottom[0], count, MPI_DOUBLE, bottom_neighbour, tag,
                    &recv[0], count, MPI_DOUBLE, bottom_neighbour, tag,
                    comm, &status);
        //unpack to bottom
#pragma omp parallel for
        for (dim_t y = 0; y < NN[1]; y++) {
            for (dim_t x = 0; x < NN[0]; x++) {
                const dim_t index = INDEX2(x,y,NN[0])*numComp;
                double *data = out.getSampleDataRW(INDEX2(x, y, NN[0]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += recv[index+comp];
                }
            }
        }
    }

    //top
    if (rz < NX[2] - 1) {
        MPI_Sendrecv(&top[0], count, MPI_DOUBLE, top_neighbour, tag,
                    &recv[0], count, MPI_DOUBLE, top_neighbour, tag,
                    comm, &status);
        //unpack to top
#pragma omp parallel for
        for (dim_t y = 0; y < NN[1]; y++) {
            for (dim_t x = 0; x < NN[0]; x++) {
                const dim_t index = INDEX2(x,y,NN[0])*numComp;
                double *data = out.getSampleDataRW(INDEX3(x, y, NN[2]-1, NN[0], NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += recv[index+comp];
                }
            }
        }
    }
}


void leftAndRight(escript::Data& out, int rx, int numComp, int rank,
                    const dim_t NN[3], const int NX[3], MPI_Comm& comm) {
    const int tag = 0;
    MPI_Status status;
    const int left_neighbour = rank - 1;
    const int right_neighbour = rank + 1;
    const dim_t count = NN[2]*NN[1]*numComp;
    std::vector<double> left(count);
    std::vector<double> right(count);
    std::vector<double> recv(count);
#pragma omp parallel for
    for (dim_t z = 0; z < NN[2]; z++) {
        for (dim_t y = 0; y < NN[1]; y++) {
            const dim_t index = INDEX2(y,z,NN[1])*numComp;
            const double *leftData = out.getSampleDataRO(INDEX3(0, y, z, NN[0], NN[1]));
            std::copy(leftData, leftData + numComp, &left[index]);

            const double *rightData = out.getSampleDataRO(INDEX3(NN[0]-1, y, z, NN[0], NN[1]));
            std::copy(rightData, rightData + numComp, &right[index]);
        }

    }

    //right
    if (rx < NX[0] - 1) {
        MPI_Sendrecv(&right[0], count, MPI_DOUBLE, right_neighbour, tag,
                &recv[0], count, MPI_DOUBLE, right_neighbour, tag,
                comm, &status);
        //unpack to right
#pragma omp parallel for
        for (dim_t z = 0; z < NN[2]; z++) {
            for (dim_t y = 0; y < NN[1]; y++) {
                const dim_t index = INDEX2(y,z,NN[1])*numComp;
                double *data = out.getSampleDataRW(INDEX3(NN[0]-1, y, z, NN[0], NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += recv[index+comp];
                }
            }
        }
    }
    //left
    if (rx) {
        MPI_Sendrecv(&left[0], count, MPI_DOUBLE, left_neighbour, tag,
                    &recv[0], count, MPI_DOUBLE, left_neighbour, tag,
                    comm, &status);
        //unpack to left
#pragma omp parallel for
        for (dim_t z = 0; z < NN[2]; z++) {
            for (dim_t y = 0; y < NN[1]; y++) {
                const dim_t index = INDEX2(y,z,NN[1])*numComp;
                double *data = out.getSampleDataRW(INDEX3(0, y, z, NN[0], NN[1]));
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += recv[index+comp];
                }
            }
        }
    }
}


//private
void Brick::shareFaces(escript::Data& out, int rx, int ry, int rz) const
{
    const int numComp = out.getDataPointSize();
    if (m_NX[0] != 1)
        leftAndRight(out, rx, numComp, m_mpiInfo->rank, m_NN, m_NX, m_mpiInfo->comm);
    if (m_NX[1] != 1)
        frontAndBack(out, ry, numComp, m_mpiInfo->rank, m_NN, m_NX, m_mpiInfo->comm);
    if (m_NX[2] != 1)
        topAndBottom(out, rz, numComp, m_mpiInfo->rank, m_NN, m_NX, m_mpiInfo->comm);
}
#endif //#ifdef ESYS_MPI


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

index_t Brick::findNode(const double *coords) const {
    const index_t NOT_MINE = -1;
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
    // get distance from subdivision origin
    double x = coords[0] - m_origin[0] - m_offset[0]*m_dx[0];
    double y = coords[1] - m_origin[1] - m_offset[1]*m_dx[1];
    double z = coords[2] - m_origin[2] - m_offset[2]*m_dx[2];

    // distance in elements
    dim_t ex = (dim_t) floor((x + 0.01*m_dx[0]) / m_dx[0]);
    dim_t ey = (dim_t) floor((y + 0.01*m_dx[1]) / m_dx[1]);
    dim_t ez = (dim_t) floor((z + 0.01*m_dx[2]) / m_dx[2]);
    dim_t start = m_order*(INDEX3(ex,ey,ez,m_NN[0],m_NN[1]));
    // set the min distance high enough to be outside the element plus a bit
    index_t closest = NOT_MINE;
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
                m_NE, m_NN));
    } else { //else ifs would go before this for other types
        throw SpeckleyException("Speckley::Brick does not support the"
                                " requested assembler");
    }
}

bool Brick::probeInterpolationAcross(int fsType_source,
        const escript::AbstractDomain& domain, int fsType_target) const
{
#ifdef USE_RIPLEY
    return speckley::probeInterpolationAcross(fsType_source, domain,
            fsType_target, 3);
#else
    return false;
#endif
}

void Brick::interpolateAcross(escript::Data& target, const escript::Data& source) const
{
#ifdef USE_RIPLEY
    if (coupler == NULL) {
        coupler = new RipleyCoupler(this, m_dx, m_mpiInfo->rank);
    }
    coupler->interpolate(target, source);
#else
    throw SpeckleyException("Speckley::Brick interpolation to unsupported domain");
#endif
}

} // end of namespace speckley

