
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

#include <speckley/Rectangle.h>
#include <speckley/Speckley.h>
#include <speckley/DefaultAssembler2D.h>
#include <speckley/WaveAssembler2D.h>
#ifdef USE_RIPLEY
#include <speckley/CrossDomainCoupler.h>
#endif

#include <escript/index.h>
#include <escript/FileWriter.h>
#include <escript/FunctionSpaceFactory.h>
#include <escript/Random.h>
#include <escript/Utils.h>

#include <boost/scoped_array.hpp>
#include <boost/math/special_functions/fpclassify.hpp> // for isnan
#include <vector>

#ifdef ESYS_HAVE_SILO
#include <silo.h>
#ifdef ESYS_MPI
#include <pmpio.h>
#endif
#endif

#include <algorithm>
#include <iomanip>
#include <limits>

namespace bm = boost::math;
namespace bp = boost::python;
using escript::FileWriter;

#ifdef ESYS_HAVE_NETCDF4
#include <ncFile.h>
#include <ncVar.h>
#include <ncDim.h>
using namespace netCDF;
#endif

namespace speckley {

Rectangle::Rectangle(escript::JMPI jmpi, int order, dim_t n0, dim_t n1, double x0, double y0, double x1,
                     double y1, int d0, int d1,
                     const std::vector<double>& points,
                     const std::vector<int>& tags,
                     const TagMap& tagnamestonums) :
    SpeckleyDomain(2, order, jmpi)
{
    if (static_cast<long>(n0 + 1) * static_cast<long>(n1 + 1)
            > std::numeric_limits<dim_t>::max())
        throw SpeckleyException("The number of elements has overflowed, this "
                "limit may be raised in future releases.");

    if (n0 <= 0 || n1 <= 0)
        throw SpeckleyException("Number of elements in each spatial dimension "
                "must be positive");

    // ignore subdivision parameters for serial run
    if (m_mpiInfo->size == 1) {
        d0=1;
        d1=1;
    }

    bool warn=false;
    std::vector<int> factors;
    int ranks = m_mpiInfo->size;
    dim_t epr[2] = {n0,n1};
    int d[2] = {d0,d1};
    if (d0<=0 || d1<=0) {
        for (int i = 0; i < 2; i++) {
            if (d[i] < 1) {
                d[i] = 1;
                continue;
            }
            epr[i] = -1; // can no longer be max
            //remove
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
        int i = epr[0] > epr[1] ? 0 : 1;
        int f = factors.back();
        factors.pop_back();
        d[i] *= f;
        epr[i] /= f;
    }
    d0 = d[0]; d1 = d[1];

    // ensure number of subdivisions is valid and nodes can be distributed
    // among number of ranks
    if (d0*d1 != m_mpiInfo->size)
        throw SpeckleyException("Invalid number of spatial subdivisions");

    if (warn) {
        std::cout << "speckley.Rectangle: Warning: Automatic domain subdivision (d0=" << d0 << ", d1="
            << d1 << "). This may not be optimal!" << std::endl;
    }

    double l0 = x1-x0;
    double l1 = y1-y0;
    m_dx[0] = l0/n0;
    m_dx[1] = l1/n1;

    if (n0 % d0 > 0) {
        n0 += d0 - (n0 % d0);
        l0 = m_dx[0]*n0;
        std::cout << "speckley.Rectangle: Warning: Adjusted number of elements and length. N0="
            << n0 << ", l0=" << l0 << std::endl;
    }
    if (n1 % d1 > 0) {
        n1 += d1 - (n1 % d1);
        l1 = m_dx[1]*n1;
        std::cout << "speckley.Rectangle: Warning: Adjusted number of elements and length. N1="
            << n1 << ", l1=" << l1 << std::endl;
    }

    if (n0/d0 < 2 || n1/d1 < 2)
        throw SpeckleyException("Too few elements for the number of ranks");

    m_gNE[0] = n0;
    m_gNE[1] = n1;
    m_origin[0] = x0;
    m_origin[1] = y0;
    m_length[0] = l0;
    m_length[1] = l1;
    m_NX[0] = d0;
    m_NX[1] = d1;

    // local number of elements
    m_NE[0] = m_gNE[0] / d0;
    m_NE[1] = m_gNE[1] / d1;

    // local number of nodes
    m_NN[0] = m_NE[0]*m_order+1;
    m_NN[1] = m_NE[1]*m_order+1;

    // bottom-left node is at (offset0,offset1) in global mesh
    m_offset[0] = n0/d0*(m_mpiInfo->rank%d0);
    m_offset[1] = n1/d1*(m_mpiInfo->rank/d0);

    populateSampleIds();

    for (TagMap::const_iterator i = tagnamestonums.begin();
            i != tagnamestonums.end(); i++) {
        setTagMap(i->first, i->second);
    }
    addPoints(points, tags);


#ifdef USE_RIPLEY
    coupler = NULL;
#endif
}


Rectangle::~Rectangle()
{
#ifdef USE_RIPLEY
    if (coupler != NULL)
        delete coupler;
#endif
}

std::string Rectangle::getDescription() const
{
    return "speckley::Rectangle";
}

bool Rectangle::operator==(const escript::AbstractDomain& other) const
{
    const Rectangle* o=dynamic_cast<const Rectangle*>(&other);
    if (o) {
        return (SpeckleyDomain::operator==(other) &&
                m_order == o->m_order &&
                m_gNE[0]==o->m_gNE[0] && m_gNE[1]==o->m_gNE[1]
                && m_origin[0]==o->m_origin[0] && m_origin[1]==o->m_origin[1]
                && m_length[0]==o->m_length[0] && m_length[1]==o->m_length[1]
                && m_NX[0]==o->m_NX[0] && m_NX[1]==o->m_NX[1]);
    }

    return false;
}

void Rectangle::readNcGrid(escript::Data& out, std::string filename,
        std::string varname, const ReaderParameters& params) const
{
#ifdef ESYS_HAVE_NETCDF4
    // check destination function space
    dim_t myN0, myN1;
    if (out.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NE[0] + 1;
        myN1 = m_NE[1] + 1;
//    } else if (out.getFunctionSpace().getTypeCode() == Elements) {
//        myN0 = m_NE[0];
//        myN1 = m_NE[1];
    } else
        throw SpeckleyException("readNcGrid(): invalid function space for output data object");

    if (params.first.size() != 2)
        throw SpeckleyException("readNcGrid(): argument 'first' must have 2 entries");

    if (params.numValues.size() != 2)
        throw SpeckleyException("readNcGrid(): argument 'numValues' must have 2 entries");

    if (params.multiplier.size() != 2)
        throw SpeckleyException("readNcGrid(): argument 'multiplier' must have 2 entries");
    for (size_t i=0; i<params.multiplier.size(); i++)
        if (params.multiplier[i]<1)
            throw SpeckleyException("readNcGrid(): all multipliers must be positive");
    if (params.reverse.size() != 2)
        throw SpeckleyException("readNcGrid(): argument 'reverse' must have 2 entries");

    // check file existence and size
    NcFile f;
    if (!escript::openNcFile(f, filename))
    {
        throw SpeckleyException("readNcGrid(): cannot open file");
    }

    NcVar var = f.getVar(varname.c_str());
    if (var.isNull())
        throw SpeckleyException("readNcGrid(): invalid variable");

    // TODO: rank>0 data support
    const int numComp = out.getDataPointSize();
    if (numComp > 1)
        throw SpeckleyException("readNcGrid(): only scalar data supported");

    const int dims = var.getDimCount();
    std::vector<long> edges(dims);
    std::vector< NcDim > vard=var.getDims();
    for (size_t i=0;i<vard.size();++i)
    {
        edges[i]=vard[i].getSize();
    }

    // is this a slice of the data object (dims!=2)?
    // note the expected ordering of edges (as in numpy: y,x)
    if ( (dims==2 && (params.numValues[1] > edges[0] || params.numValues[0] > edges[1]))
            || (dims==1 && params.numValues[1]>1) ) {
        throw SpeckleyException("readNcGrid(): not enough data in file");
    }

    // check if this rank contributes anything
    if (params.first[0] >= m_offset[0]+myN0 ||
            params.first[0]+params.numValues[0]*params.multiplier[0] <= m_offset[0] ||
            params.first[1] >= m_offset[1]+myN1 ||
            params.first[1]+params.numValues[1]*params.multiplier[1] <= m_offset[1])
        return;

    // now determine how much this rank has to write

    // first coordinates in data object to write to
    const dim_t first0 = std::max(dim_t(0), params.first[0]-m_offset[0]);
    const dim_t first1 = std::max(dim_t(0), params.first[1]-m_offset[1]);
    // indices to first value in file (not accounting for reverse yet)
    dim_t idx0 = std::max(dim_t(0), m_offset[0]-params.first[0]);
    dim_t idx1 = std::max(dim_t(0), m_offset[1]-params.first[1]);
    // number of values to read
    const dim_t num0 = std::min(params.numValues[0]-idx0, myN0-first0);
    const dim_t num1 = std::min(params.numValues[1]-idx1, myN1-first1);

    // make sure we read the right block if going backwards through file
    if (params.reverse[0])
        idx0 = edges[dims-1]-num0-idx0;
    if (dims>1 && params.reverse[1])
        idx1 = edges[dims-2]-num1-idx1;

    std::vector<double> values(num0*num1);
    std::vector<size_t> startindex;
    std::vector<size_t> counts;
    if (dims==2) {
        // var->set_cur(idx1, idx0);                // from old API
        startindex.push_back(idx1);
        startindex.push_back(idx0);
        counts.push_back(num1);
        counts.push_back(num0);
        var.getVar(startindex, counts, &values[0]);
    } else {
        //var->set_cur(idx0);
        //var->get(&values[0], num0);
        startindex.push_back(idx0);
        counts.push_back(num0);
        var.getVar(startindex, counts, &values[0]);
    }

    const int dpp = out.getNumDataPointsPerSample();
    out.requireWrite();

    // helpers for reversing
    const dim_t x0 = (params.reverse[0] ? num0-1 : 0);
    const int x_mult = (params.reverse[0] ? -1 : 1);
    const dim_t y0 = (params.reverse[1] ? num1-1 : 0);
    const int y_mult = (params.reverse[1] ? -1 : 1);

    for (index_t y=0; y<num1; y++) {
#pragma omp parallel for
        for (index_t x=0; x<num0; x++) {
            const dim_t baseIndex = first0+x*params.multiplier[0]
                                  +(first1+y*params.multiplier[1])*myN0;
            const dim_t srcIndex = (y0+y_mult*y)*num0+(x0+x_mult*x);
            if (!bm::isnan(values[srcIndex])) {
                for (index_t m1=0; m1<params.multiplier[1]; m1++) {
                    for (index_t m0=0; m0<params.multiplier[0]; m0++) {
                        const dim_t dataIndex = baseIndex+m0+m1*myN0;
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
    throw SpeckleyException("readNcGrid(): not compiled with netCDF support");
#endif
}

void Rectangle::readBinaryGrid(escript::Data& out, std::string filename,
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

void Rectangle::readBinaryGridFromZipped(escript::Data& out, std::string filename,
                               const ReaderParameters& params) const
{
#ifdef ESYS_HAVE_BOOST_IO
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
            throw SpeckleyException("readBinaryGridFromZipped(): invalid or unsupported datatype");
    }
#else
    throw SpeckleyException("readBinaryGridFromZipped(): not built with zip support");
#endif
}

template<typename ValueType>
void Rectangle::readBinaryGridImpl(escript::Data& out, const std::string& filename,
                                   const ReaderParameters& params) const
{
    // check destination function space
    dim_t myN0, myN1;
    if (out.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NE[0] + 1;
        myN1 = m_NE[1] + 1;
//    } else if (out.getFunctionSpace().getTypeCode() == Elements) {
//        myN0 = m_NE[0];
//        myN1 = m_NE[1];
    } else
        throw SpeckleyException("readBinaryGrid(): invalid function space for output data object");

    if (params.first.size() != 2)
        throw SpeckleyException("readBinaryGrid(): argument 'first' must have 2 entries");

    if (params.numValues.size() != 2)
        throw SpeckleyException("readBinaryGrid(): argument 'numValues' must have 2 entries");

    if (params.multiplier.size() != 2)
        throw SpeckleyException("readBinaryGrid(): argument 'multiplier' must have 2 entries");
    for (size_t i=0; i<params.multiplier.size(); i++)
        if (params.multiplier[i]<1)
            throw SpeckleyException("readBinaryGrid(): all multipliers must be positive");
    if (params.reverse[0] != 0 || params.reverse[1] != 0)
        throw SpeckleyException("readBinaryGrid(): reversing not supported yet");

    // check file existence and size
    std::ifstream f(filename.c_str(), std::ifstream::binary);
    if (f.fail()) {
        throw SpeckleyException("readBinaryGrid(): cannot open file " + filename);
    }
    f.seekg(0, std::ios::end);
    const int numComp = out.getDataPointSize();
    const dim_t filesize = f.tellg();
    const dim_t reqsize = params.numValues[0]*params.numValues[1]*numComp*sizeof(ValueType);
    if (filesize < reqsize) {
        f.close();
        throw SpeckleyException("readBinaryGrid(): not enough data in file");
    }

    // check if this rank contributes anything
    if (params.first[0] >= m_offset[0]+myN0 ||
            params.first[0]+(params.numValues[0]*params.multiplier[0]) <= m_offset[0] ||
            params.first[1] >= m_offset[1]+myN1 ||
            params.first[1]+(params.numValues[1]*params.multiplier[1]) <= m_offset[1]) {
        f.close();
        return;
    }

    // now determine how much this rank has to write

    // first coordinates in data object to write to
    const dim_t first0 = std::max(dim_t(0), params.first[0]-m_offset[0]);
    const dim_t first1 = std::max(dim_t(0), params.first[1]-m_offset[1]);
    // indices to first value in file
    const dim_t idx0 = std::max(dim_t(0), (m_offset[0]/params.multiplier[0])-params.first[0]);
    const dim_t idx1 = std::max(dim_t(0), (m_offset[1]/params.multiplier[1])-params.first[1]);
    // if restX > 0 the first value in the respective dimension has been
    // written restX times already in a previous rank so this rank only
    // contributes (multiplier-rank) copies of that value
    const dim_t rest0 = m_offset[0]%params.multiplier[0];
    const dim_t rest1 = m_offset[1]%params.multiplier[1];
    // number of values to read
    const dim_t num0 = std::min(params.numValues[0]-idx0, myN0-first0);
    const dim_t num1 = std::min(params.numValues[1]-idx1, myN1-first1);

    out.requireWrite();
    std::vector<ValueType> values(num0*numComp);
    const int dpp = out.getNumDataPointsPerSample();

    for (dim_t y=0; y<num1; y++) {
        const dim_t fileofs = numComp*(idx0+(idx1+y)*params.numValues[0]);
        f.seekg(fileofs*sizeof(ValueType));
        f.read((char*)&values[0], num0*numComp*sizeof(ValueType));
        const dim_t m1limit = (y==0 ? params.multiplier[1]-rest1 : params.multiplier[1]);
        dim_t dataYbase = first1+y*params.multiplier[1];
        if (y>0)
            dataYbase -= rest1;
        for (dim_t x=0; x<num0; x++) {
            const dim_t m0limit = (x==0 ? params.multiplier[0]-rest0 : params.multiplier[0]);
            dim_t dataXbase = first0+x*params.multiplier[0];
            if (x>0)
                dataXbase -= rest0;
            // write a block of mult0 x mult1 identical values into Data object
            for (dim_t m1=0; m1<m1limit; m1++) {
                const dim_t dataY = dataYbase+m1;
                if (dataY >= myN1)
                    break;
                for (dim_t m0=0; m0<m0limit; m0++) {
                    const dim_t dataX = dataXbase+m0;
                    if (dataX >= myN0)
                        break;
                    const dim_t dataIndex = dataX+dataY*m_NN[0];
                    double* dest = out.getSampleDataRW(dataIndex*m_order);
                    for (int c=0; c<numComp; c++) {
                        ValueType val = values[x*numComp+c];
                        if (params.byteOrder != BYTEORDER_NATIVE) {
                            char* cval = reinterpret_cast<char*>(&val);
                            // this will alter val!!
                            if (sizeof(ValueType) > 4) {
                                byte_swap64(cval);
                            } else {
                                byte_swap32(cval);
                            }
                        }
                        if (!bm::isnan(val)) {
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
    interpolateFromCorners(out);
}

#ifdef ESYS_HAVE_BOOST_IO
template<typename ValueType>
void Rectangle::readBinaryGridZippedImpl(escript::Data& out, const std::string& filename,
                                   const ReaderParameters& params) const
{
    // check destination function space
    dim_t myN0, myN1;
    if (out.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NE[0] + 1;
        myN1 = m_NE[1] + 1;
//    } else if (out.getFunctionSpace().getTypeCode() == Elements) {
//        myN0 = m_NE[0];
//        myN1 = m_NE[1];
    } else
        throw SpeckleyException("readBinaryGrid(): invalid function space for output data object");

    // check file existence and size
    std::ifstream f(filename.c_str(), std::ifstream::binary);
    if (f.fail()) {
        throw SpeckleyException(strerror(errno));//"readBinaryGridFromZipped(): cannot open file");
    }
    f.seekg(0, std::ios::end);
    const int numComp = out.getDataPointSize();
    dim_t filesize = f.tellg();
    f.seekg(0, std::ios::beg);
    std::vector<char> compressed(filesize);
    f.read((char*)&compressed[0], filesize);
    f.close();
    std::vector<char> decompressed = unzip(compressed);
    filesize = decompressed.size();
    const dim_t reqsize = params.numValues[0]*params.numValues[1]*numComp*sizeof(ValueType);
    if (filesize < reqsize) {
        throw SpeckleyException("readBinaryGridFromZipped(): not enough data in file");
    }

    // check if this rank contributes anything
    if (params.first[0] >= m_offset[0]+myN0 ||
            params.first[0]+(params.numValues[0]*params.multiplier[0]) <= m_offset[0] ||
            params.first[1] >= m_offset[1]+myN1 ||
            params.first[1]+(params.numValues[1]*params.multiplier[1]) <= m_offset[1]) {
        return;
    }

    // now determine how much this rank has to write

    // first coordinates in data object to write to
    const dim_t first0 = std::max(dim_t(0), params.first[0]-m_offset[0]);
    const dim_t first1 = std::max(dim_t(0), params.first[1]-m_offset[1]);
    // indices to first value in file
    const dim_t idx0 = std::max(dim_t(0), (m_offset[0]/params.multiplier[0])-params.first[0]);
    const dim_t idx1 = std::max(dim_t(0), (m_offset[1]/params.multiplier[1])-params.first[1]);
    // if restX > 0 the first value in the respective dimension has been
    // written restX times already in a previous rank so this rank only
    // contributes (multiplier-rank) copies of that value
    const dim_t rest0 = m_offset[0]%params.multiplier[0];
    const dim_t rest1 = m_offset[1]%params.multiplier[1];
    // number of values to read
    const dim_t num0 = std::min(params.numValues[0]-idx0, myN0-first0);
    const dim_t num1 = std::min(params.numValues[1]-idx1, myN1-first1);

    out.requireWrite();
    std::vector<ValueType> values(num0*numComp);
    const int dpp = out.getNumDataPointsPerSample();

    for (dim_t y=0; y<num1; y++) {
        const dim_t fileofs = numComp*(idx0+(idx1+y)*params.numValues[0]);
        memcpy((char*)&values[0], (char*)&decompressed[fileofs*sizeof(ValueType)], num0*numComp*sizeof(ValueType));
        const dim_t m1limit = (y==0 ? params.multiplier[1]-rest1 : params.multiplier[1]);
        dim_t dataYbase = first1+y*params.multiplier[1];
        if (y>0)
            dataYbase -= rest1;
        for (dim_t x=0; x<num0; x++) {
            const dim_t m0limit = (x==0 ? params.multiplier[0]-rest0 : params.multiplier[0]);
            dim_t dataXbase = first0+x*params.multiplier[0];
            if (x>0)
                dataXbase -= rest0;
            // write a block of mult0 x mult1 identical values into Data object
            for (dim_t m1=0; m1<m1limit; m1++) {
                const dim_t dataY = dataYbase+m1;
                if (dataY >= myN1)
                    break;
                for (dim_t m0=0; m0<m0limit; m0++) {
                    const dim_t dataX = dataXbase+m0;
                    if (dataX >= myN0)
                        break;
                    const dim_t dataIndex = dataX+dataY*m_NN[0];
                    double* dest = out.getSampleDataRW(dataIndex*m_order);
                    for (int c=0; c<numComp; c++) {
                        ValueType val = values[x*numComp+c];

                        if (params.byteOrder != BYTEORDER_NATIVE) {
                            char* cval = reinterpret_cast<char*>(&val);
                            // this will alter val!!
                            if (sizeof(ValueType) > 4) {
                                byte_swap64(cval);
                            } else {
                                byte_swap32(cval);
                            }
                        }
                        if (!bm::isnan(val)) {
                            for (int q=0; q<dpp; q++) {
                                *dest++ = static_cast<double>(val);
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

void Rectangle::interpolateFromCorners(escript::Data &out) const
{
    const int numComp = out.getDataPointSize();
    //interpolate the missing portions
#pragma omp parallel for
    for (dim_t y = 0; y < m_NN[1]; y++) {
        for (dim_t x = 0; x < m_NN[0]; x++) {
            //skip the points we have values for
            if (y % m_order == 0 && x % m_order == 0)
                continue;
            //point location in element: x,y
            const double px = point_locations[m_order-2][x%m_order];
            const double py = point_locations[m_order-2][y%m_order];

            double *point = out.getSampleDataRW(INDEX2(x, y, m_NN[0]));

            const dim_t left = x - x%m_order;
            const dim_t right = left < m_NN[0] - 1 ? left + m_order : left;
            const dim_t front = y - y%m_order;
            const dim_t back = front < m_NN[1] - 1 ? front + m_order : front;


            const double *lowleft = out.getSampleDataRO(
                    INDEX2(left, front, m_NN[0]));
            const double *lowright = out.getSampleDataRO(
                    INDEX2(right, front, m_NN[0]));
            const double *highleft = out.getSampleDataRO(
                    INDEX2(left, back, m_NN[0]));
            const double *highright = out.getSampleDataRO(
                    INDEX2(right, back, m_NN[0]));

            for (int comp = 0; comp < numComp; comp++) {
                point[comp] = highright[comp]*px*py
                            + highleft[comp]*(1-px)*py
                            + lowright[comp]*px*(1-py)
                            + lowleft[comp]*(1-px)*(1-py);
            }
        }
    }
}

void Rectangle::writeBinaryGrid(const escript::Data& in, std::string filename,
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
void Rectangle::writeBinaryGridImpl(const escript::Data& in,
                                    const std::string& filename, int byteOrder) const
{
    // check function space and determine number of points
    dim_t myN0, myN1;
    dim_t totalN0, totalN1;
    if (in.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NE[0]+1;
        myN1 = m_NE[1]+1;
        totalN0 = m_gNE[0]+1;
        totalN1 = m_gNE[1]+1;
//    } else if (in.getFunctionSpace().getTypeCode() == Elements) {
//        myN0 = m_NE[0];
//        myN1 = m_NE[1];
//        totalN0 = m_gNE[0];
//        totalN1 = m_gNE[1];
    } else
        throw SpeckleyException("writeBinaryGrid(): invalid function space of data object");

    const int numComp = in.getDataPointSize();
    const int dpp = in.getNumDataPointsPerSample();

    if (numComp > 1 || dpp > 1)
        throw SpeckleyException("writeBinaryGrid(): only scalar, single-value data supported");

    const dim_t fileSize = sizeof(ValueType)*numComp*dpp*totalN0*totalN1;

    // from here on we know that each sample consists of one value
    FileWriter fw;
#ifdef _WIN32
    fw.openFile(filename, fileSize, true); // open in binary mode to allow seek
#else
    fw.openFile(filename, fileSize);
#endif
    MPIBarrier();

    for (index_t y=0; y<myN1; y++) {
        const dim_t fileofs = (m_offset[0]+(m_offset[1]+y)*totalN0)*sizeof(ValueType);
        std::ostringstream oss;

        for (index_t x=0; x<myN0; x++) {
            const double* sample = in.getSampleDataRO((y*m_NN[0]+x)*m_order);
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
    fw.close();
}

void Rectangle::write(const std::string& filename) const
{
    throw SpeckleyException("write: not supported");
}

void Rectangle::dump(const std::string& fileName) const
{
#ifdef ESYS_HAVE_SILO
    std::string fn(fileName);
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
        throw SpeckleyException("dump: Could not create Silo file");

    /*
    if (driver==DB_HDF5) {
        // gzip level 1 already provides good compression with minimal
        // performance penalty. Some tests showed that gzip levels >3 performed
        // rather badly on escript data both in terms of time and space
        DBSetCompression("ERRMODE=FALLBACK METHOD=GZIP LEVEL=1");
    }
    */

    const dim_t NN0 = m_NN[0];
    const dim_t NN1 = m_NN[1];
    std::unique_ptr<double[]> x(new double[NN0]);
    std::unique_ptr<double[]> y(new double[NN1]);
    double* coords[2] = { x.get(), y.get() };
#pragma omp parallel
    {
#pragma omp for nowait
        for (dim_t i0 = 0; i0 < NN0; i0++) {
            coords[0][i0]=getLocalCoordinate(i0, 0);
        }
#pragma omp for nowait
        for (dim_t i1 = 0; i1 < NN1; i1++) {
            coords[1][i1]=getLocalCoordinate(i1, 1);
        }
    }
    std::vector<int> dims(m_NN, m_NN+2);

    // write mesh
    DBPutQuadmesh(dbfile, "mesh", NULL, coords, &dims[0], 2, DB_DOUBLE,
            DB_COLLINEAR, NULL);

    // write node ids
    DBPutQuadvar1(dbfile, "nodeId", "mesh", (void*)&m_nodeId[0], &dims[0], 2,
            NULL, 0, DB_INT, DB_NODECENT, NULL);

    // write element ids
    dims.assign(m_NE, m_NE+2);
    DBPutQuadvar1(dbfile, "elementId", "mesh", (void*)&m_elementId[0],
            &dims[0], 2, NULL, 0, DB_INT, DB_ZONECENT, NULL);

    // rank 0 writes multimesh and multivar
    if (m_mpiInfo->rank == 0) {
        std::vector<std::string> tempstrings;
        std::vector<char*> names;
        for (dim_t i=0; i<m_mpiInfo->size; i++) {
            std::stringstream path;
            path << "/block" << std::setw(4) << std::setfill('0') << std::right << i << "/mesh";
            tempstrings.push_back(path.str());
            names.push_back((char*)tempstrings.back().c_str());
        }
        std::vector<int> types(m_mpiInfo->size, DB_QUAD_RECT);
        DBSetDir(dbfile, "/");
        DBPutMultimesh(dbfile, "multimesh", m_mpiInfo->size, &names[0],
               &types[0], NULL);
        tempstrings.clear();
        names.clear();
        for (dim_t i=0; i<m_mpiInfo->size; i++) {
            std::stringstream path;
            path << "/block" << std::setw(4) << std::setfill('0') << std::right << i << "/nodeId";
            tempstrings.push_back(path.str());
            names.push_back((char*)tempstrings.back().c_str());
        }
        types.assign(m_mpiInfo->size, DB_QUADVAR);
        DBPutMultivar(dbfile, "nodeId", m_mpiInfo->size, &names[0],
               &types[0], NULL);
        tempstrings.clear();
        names.clear();
        for (dim_t i=0; i<m_mpiInfo->size; i++) {
            std::stringstream path;
            path << "/block" << std::setw(4) << std::setfill('0') << std::right << i << "/elementId";
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

#else // ESYS_HAVE_SILO
    throw SpeckleyException("dump: no Silo support");
#endif
}

const dim_t* Rectangle::borrowSampleReferenceIDs(int fsType) const
{
    switch (fsType) {
        case DegreesOfFreedom:
        case Nodes:
            return &m_nodeId[0];
        case Elements:
        case ReducedElements:
            return &m_elementId[0];
        case Points:
            return &m_diracPointNodeIDs[0];
        default:
            break;
    }

    std::stringstream msg;
    msg << "borrowSampleReferenceIDs: invalid function space type" << fsType;
    throw SpeckleyException(msg.str());
}

bool Rectangle::ownSample(int fsType, index_t id) const
{
#ifdef ESYS_MPI
    if (getMPISize() > 1) {
        if (fsType == Nodes || fsType == Elements) {
            const index_t myFirstNode = m_nodeDistribution[getMPIRank()];
            const index_t myLastNode = m_nodeDistribution[getMPIRank()+1];
            const index_t k = m_nodeId[id];
            return (myFirstNode <= k && k < myLastNode);
        } else {
            throw SpeckleyException("ownSample: unsupported function space type");
        }
    }
#endif
    return true;
}

void Rectangle::setToNormal(escript::Data& out) const
{
    throw SpeckleyException("setToNormal not implemented");
}

void Rectangle::setToSize(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == Elements) {
        out.requireWrite();
        const dim_t numQuad = m_order + 1;
        const dim_t numElements = getNumElements();
        const double *quad_locs = point_locations[m_order-2];
        //since elements are uniform, calc the first and copy to others
        double* first_element = out.getSampleDataRW(0);
#pragma omp parallel for
        for (short qy = 0; qy < m_order; qy++) {
            const double y = quad_locs[qy+1] - quad_locs[qy];
            for (short qx = 0; qx < m_order; qx++) {
                const double x = quad_locs[qx+1] - quad_locs[qx];
                first_element[qx + qy*numQuad] = sqrt(x*x + y*y);
            }
        }
        const short top_start = (numQuad - 1)*numQuad;
        for (short edge = 0; edge < m_order; edge++) {
            //right edge = left edge
            first_element[(edge+1)*numQuad - 1] = first_element[edge*numQuad];
            //top edge = bottom edge
            first_element[top_start + edge] = first_element[edge];
        }
        //top-right corner
        first_element[numQuad*numQuad - 1] = first_element[0];
        const size_t size = numQuad*numQuad*sizeof(double);
#pragma omp parallel for
        for (index_t k = 0; k < numElements; ++k) {
            double* o = out.getSampleDataRW(k);
            memcpy(o, first_element, size);
        }
    } else {
        std::stringstream msg;
        msg << "setToSize: invalid function space type "
            << out.getFunctionSpace().getTypeCode();
        throw SpeckleyException(msg.str());
    }
}

void Rectangle::Print_Mesh_Info(const bool full) const
{
    SpeckleyDomain::Print_Mesh_Info(full);
    if (full) {
        std::cout << "     Id  Coordinates" << std::endl;
        std::cout.precision(15);
        std::cout.setf(std::ios::scientific, std::ios::floatfield);
        for (index_t i=0; i < getNumNodes(); i++) {
            std::cout << "  " << std::setw(5) << m_nodeId[i]
                << "  " << getLocalCoordinate(i%m_NN[0], 0)
                << "  " << getLocalCoordinate(i/m_NN[0], 1) << std::endl;
        }
    }
}


//protected
void Rectangle::assembleCoordinates(escript::Data& arg) const
{
    int numDim = m_numDim;
    if (!arg.isDataPointShapeEqual(1, &numDim))
        throw SpeckleyException("setToX: Invalid Data object shape");
    if (!arg.numSamplesEqual(1, getNumNodes()))
        throw SpeckleyException("setToX: Illegal number of samples in Data object");

    const dim_t NN0 = m_NN[0];
    const dim_t NN1 = m_NN[1];
    arg.requireWrite();
#pragma omp parallel for
    for (dim_t y = 0; y < NN1; y++) {
        for (dim_t x = 0; x < NN0; x++) {
            double *point = arg.getSampleDataRW(y*NN0 + x);
            point[0] = getLocalCoordinate(x, 0);
            point[1] = getLocalCoordinate(y, 1);
        }
    }
}

//protected
void Rectangle::assembleGradient(escript::Data& out, const escript::Data& in) const
{
    escript::Data converted;

    if (in.getFunctionSpace().getTypeCode() != Elements) {
        converted = escript::Data(in, escript::function(*this));
    } else {
        converted = in;
    }

    if (m_order == 2) {
        if (in.isComplex())
            gradient_order2<cplx_t>(out,converted);
        else
            gradient_order2<real_t>(out,converted);
    } else if (m_order == 3) {
        if (in.isComplex())
            gradient_order3<cplx_t>(out,converted);
        else
            gradient_order3<real_t>(out,converted);
    } else if (m_order == 4) {
        if (in.isComplex())
            gradient_order4<cplx_t>(out,converted);
        else
            gradient_order4<real_t>(out,converted);
    } else if (m_order == 5) {
        if (in.isComplex())
            gradient_order5<cplx_t>(out,converted);
        else
            gradient_order5<real_t>(out,converted);
    } else if (m_order == 6) {
        if (in.isComplex())
            gradient_order6<cplx_t>(out,converted);
        else
            gradient_order6<real_t>(out,converted);
    } else if (m_order == 7) {
        if (in.isComplex())
            gradient_order7<cplx_t>(out,converted);
        else
            gradient_order7<real_t>(out,converted);
    } else if (m_order == 8) {
        if (in.isComplex())
            gradient_order8<cplx_t>(out,converted);
        else
            gradient_order8<real_t>(out,converted);
    } else if (m_order == 9) {
        if (in.isComplex())
            gradient_order9<cplx_t>(out,converted);
        else
            gradient_order9<real_t>(out,converted);
    } else if (m_order == 10) {
        if (in.isComplex())
            gradient_order10<cplx_t>(out,converted);
        else
            gradient_order10<real_t>(out,converted);
    }
}

//protected
void Rectangle::assembleIntegrate(std::vector<real_t>& integrals,
                                  const escript::Data& arg) const
{
    assembleIntegrateWorker<real_t>(integrals, arg);
}

//protected
void Rectangle::assembleIntegrate(std::vector<cplx_t>& integrals,
                                  const escript::Data& arg) const
{
    assembleIntegrateWorker<cplx_t>(integrals, arg);
}

//private
template<typename Scalar>
void Rectangle::assembleIntegrateWorker(std::vector<Scalar>& integrals, const escript::Data& arg) const
{
    const int fs = arg.getFunctionSpace().getTypeCode();
    if (fs != Elements && fs != Points)
        throw new SpeckleyException("Speckley doesn't currently support integrals of non-Element functionspaces");
    if (!arg.actsExpanded() && fs != Points)
        throw new SpeckleyException("Speckley doesn't currently support unexpanded data");

    if(fs == Points){
#ifdef ESYS_MPI
        if(getMPIRank() == 0)
#endif
        integrals[0] += arg.getNumberOfTaggedValues();
    } else if (m_order == 2) {
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

/* This is a wrapper for filtered (and non-filtered) randoms
 * For detailed doco see randomFillWorker
*/
escript::Data Rectangle::randomFill(const escript::DataTypes::ShapeType& shape,
                                    const escript::FunctionSpace& fs,
                                    long seed, const bp::tuple& filter) const
{
    const int numvals = escript::DataTypes::noValues(shape);
    const int per_element = (m_order+1) * (m_order+1) * numvals;
    if (len(filter) > 0) {
        throw SpeckleyException("Speckley does not support filters.");
    }

    double* src = new double[m_NE[0] * m_NE[1] * per_element * numvals];
    escript::randomFillArray(seed, src, m_NE[0]*m_NE[1]*per_element, m_mpiInfo);
    escript::Data res(0, shape, escript::function(*this), true);
    int current = 0;
    for (index_t ei = 0; ei < m_NE[1]; ++ei) {
        for (index_t ej = 0; ej < m_NE[0]; ++ej) {
            double* e = res.getSampleDataRW(INDEX2(ej, ei, m_NE[0]));
            memcpy(e, &src[current], sizeof(double)*per_element);
            current += per_element;
        }
    }
    delete[] src;

    if (res.getFunctionSpace() != fs) {
        return escript::Data(res, fs);
    }
    return res;
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
    const index_t left = (m_offset[0]==0 ? 0 : 1);
    const index_t bottom = (m_offset[1]==0 ? 0 : 1);
    for (dim_t k=1; k<m_mpiInfo->size; k++) {
        index_t rank_left = (k-1)%m_NX[0] == 0 ? 0 : 1;
        index_t rank_bottom = (k-1)/m_NX[0] == 0 ? 0 : 1;
        m_nodeDistribution[k] = m_nodeDistribution[k-1]
                                + (m_NN[0]-rank_left)*(m_NN[1]-rank_bottom);
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


    if (bottom && left) {
        //get lower-left node
        m_nodeId[0] = m_nodeDistribution[m_mpiInfo->rank - m_NX[0]] - 1;
    }
    if (bottom) {
        //DOF size of the neighbouring rank
        index_t rankDOF = m_nodeDistribution[m_mpiInfo->rank - m_NX[0] + 1]
                        - m_nodeDistribution[m_mpiInfo->rank - m_NX[0]];
        //beginning of last row of neighbouring rank
        index_t begin = m_nodeDistribution[m_mpiInfo->rank - m_NX[0]]
                      + rankDOF - m_NN[0];
        for (index_t x = left; x < m_NN[0]; x++) {
            m_nodeId[x] = begin + x;
        }
    }
    if (left) {
        //is the rank to the left itself right of another rank
        index_t rank_left = (m_mpiInfo->rank - 1) % m_NX[0] == 0 ? 0 : 1;
        //end of first owned row of neighbouring rank
        index_t end = m_nodeDistribution[m_mpiInfo->rank - 1]
                    + m_NN[0] - rank_left - 1;
        for (index_t y = bottom; y < m_NN[1]; y++) {
            m_nodeId[y*m_NN[0]] = end + (y-bottom)*(m_NN[0]-rank_left);
        }
    }

#pragma omp parallel for
    for (index_t y = bottom; y < m_NN[1]; y++) {
        for (index_t x = left; x < m_NN[0]; x++) {
            m_nodeId[y*m_NN[0]+x] = m_nodeDistribution[m_mpiInfo->rank] + (y-bottom)*(m_NN[0]-left) + (x-left);
        }
    }

    m_nodeTags.assign(getNumNodes(), 0);
    updateTagsInUse(Nodes);

    m_elementTags.assign(getNumElements(), 0);
    updateTagsInUse(Elements);
}

void Rectangle::reduceElements(escript::Data& out, const escript::Data& in) const
{
    if (m_order == 2) {
        if (in.isComplex())
            reduction_order2<cplx_t>(in, out);
        else
            reduction_order2<real_t>(in, out);
    } else if (m_order == 3) {
        if (in.isComplex())
            reduction_order3<cplx_t>(in, out);
        else
            reduction_order3<real_t>(in, out);
    } else if (m_order == 4) {
        if (in.isComplex())
            reduction_order4<cplx_t>(in, out);
        else
            reduction_order4<real_t>(in, out);
    } else if (m_order == 5) {
        if (in.isComplex())
            reduction_order5<cplx_t>(in, out);
        else
            reduction_order5<real_t>(in, out);
    } else if (m_order == 6) {
        if (in.isComplex())
            reduction_order6<cplx_t>(in, out);
        else
            reduction_order6<real_t>(in, out);
    } else if (m_order == 7) {
        if (in.isComplex())
            reduction_order7<cplx_t>(in, out);
        else
            reduction_order7<real_t>(in, out);
    } else if (m_order == 8) {
        if (in.isComplex())
            reduction_order8<cplx_t>(in, out);
        else
            reduction_order8<real_t>(in, out);
    } else if (m_order == 9) {
        if (in.isComplex())
            reduction_order9<cplx_t>(in, out);
        else
            reduction_order9<real_t>(in, out);
    } else if (m_order == 10) {
        if (in.isComplex())
            reduction_order10<cplx_t>(in, out);
        else
            reduction_order10<real_t>(in, out);
    }
}

//protected
void Rectangle::interpolateNodesOnElements(escript::Data& out,
                                           const escript::Data& in,
                                           bool reduced) const
{
    if (in.isComplex())
        interpolateNodesOnElementsWorker<cplx_t>(out, in, reduced);
    else
        interpolateNodesOnElementsWorker<real_t>(out, in, reduced);
}

//private
template<typename Scalar>
void Rectangle::interpolateNodesOnElementsWorker(escript::Data& out,
                                                 const escript::Data& in,
                                                 bool reduced) const
{
    const dim_t numComp = in.getDataPointSize();
    const dim_t NE0 = m_NE[0];
    const dim_t NE1 = m_NE[1];
    const int quads = m_order + 1;
    const int max_x = m_NN[0];
    const Scalar zero = static_cast<Scalar>(0);

    out.requireWrite();
    if (reduced) { //going to ReducedElements
        escript::Data funcIn(in, escript::function(*this));
        reduceElements(out, funcIn);
        return;
    }
#pragma omp parallel for
    for (dim_t ey = 0; ey < NE1; ey++) {
        for (dim_t ex = 0; ex < NE0; ex++) {
            Scalar* e_out = out.getSampleDataRW(ex + ey*NE0, zero);
            dim_t start = ex*m_order + ey*max_x*m_order;
            int quad = 0;
            for (int qy = 0; qy < quads; qy++) {
                for (int qx = 0; qx < quads; qx++, quad++) {
                    const Scalar* n_in = in.getSampleDataRO(start + max_x*qy + qx, zero);
                    memcpy(e_out+quad*numComp, n_in, sizeof(Scalar) * numComp);
                }
            }
        }
    }
}

#ifdef ESYS_MPI
//protected
void Rectangle::balanceNeighbours(escript::Data& data, bool average) const
{
    if (data.isComplex())
        balanceNeighboursWorker<cplx_t>(data, average);
    else
        balanceNeighboursWorker<real_t>(data, average);
}

//private
template<typename Scalar>
void Rectangle::balanceNeighboursWorker(escript::Data& data, bool average) const
{
    if (m_NX[0] * m_NX[1] == 1) {
        return;
    }
    const int numComp = data.getDataPointSize();
    const int rx = m_mpiInfo->rank % m_NX[0];
    const int ry = m_mpiInfo->rank / m_NX[0];
    const Scalar zero = static_cast<Scalar>(0);

    //include bordering ranks in summation
    if (m_NX[1] != 1)
        shareVertical<Scalar>(data, rx, ry);
    if (m_NX[0] != 1)
        shareSides<Scalar>(data, rx, ry);
    if (m_NX[0] != 1 && m_NX[1] != 1) {
        shareCorners<Scalar>(data, rx, ry);
        if (!average)
            return;
        //averaging out corners
        // bottom left
        if (rx && ry) {
            Scalar* values = data.getSampleDataRW(0, zero);
            for (int comp = 0; comp < numComp; comp++) {
                values[comp] /= 2;
            }
        }
        // bottom right
        if (rx < (m_NX[0] - 1) && ry) {
            Scalar* values = data.getSampleDataRW(m_NN[0]-1, zero);
            for (int comp = 0; comp < numComp; comp++) {
                values[comp] /= 2;
            }
        }
        // top left
        if (rx && ry < (m_NX[0] - 1)) {
            Scalar* values = data.getSampleDataRW((m_NN[1]-1)*m_NN[0], zero);
            for (int comp = 0; comp < numComp; comp++) {
                values[comp] /= 2;
            }
        }
        // top right
        if (rx < (m_NX[0] - 1) && ry < (m_NX[0] - 1)) {
            Scalar* values = data.getSampleDataRW(m_NN[1]*m_NN[0] - 1, zero);
            for (int comp = 0; comp < numComp; comp++) {
                values[comp] /= 2;
            }
        }
    }
    if (!average)
        return;
    // average shared-edges in x and y
    //left
    if (rx) {
#pragma omp parallel for
        for (dim_t qy = 0; qy < m_NN[1]; qy++) {
            Scalar* values = data.getSampleDataRW(qy*m_NN[0], zero);
            for (int comp = 0; comp < numComp; comp++) {
                values[comp] /= 2;
            }
        }
    }
    //right
    if (rx < m_NX[0] - 1) {
#pragma omp parallel for
        for (dim_t qy = 0; qy < m_NN[1]; qy++) {
            Scalar* values = data.getSampleDataRW(qy*m_NN[0] + m_NN[0] - 1, zero);
            for (int comp = 0; comp < numComp; comp++) {
                values[comp] /= 2;
            }
        }
    }
    //bottom
    if (ry) {
#pragma omp parallel for
        for (dim_t qx = 0; qx < m_NN[0]; qx++) {
            Scalar* values = data.getSampleDataRW(qx, zero);
            for (int comp = 0; comp < numComp; comp++) {
                values[comp] /= 2;
            }
        }
    }
    //top
    if (ry < m_NX[1] - 1) {
        const dim_t start = (m_NN[1]-1)*m_NN[0];
#pragma omp parallel for
        for (dim_t qx = 0; qx < m_NN[0]; qx++) {
            Scalar* values = data.getSampleDataRW(start + qx, zero);
            for (int comp = 0; comp < numComp; comp++) {
                values[comp] /= 2;
            }
        }
    }
}
#endif //#ifdef ESYS_MPI

//protected
void Rectangle::interpolateElementsOnNodes(escript::Data& out,
                                           const escript::Data& in) const
{
    if (in.isComplex())
        interpolateElementsOnNodesWorker<cplx_t>(out, in);
    else
        interpolateElementsOnNodesWorker<real_t>(out, in);
}

//private
template<typename Scalar>
void Rectangle::interpolateElementsOnNodesWorker(escript::Data& out,
                                                 const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    const dim_t NE0 = m_NE[0];
    const dim_t NE1 = m_NE[1];
    const int quads = m_order + 1;
    const dim_t max_x = (m_order*NE0) + 1;
    const dim_t max_y = (m_order*NE1) + 1;
    const int inFS = in.getFunctionSpace().getTypeCode();
    const Scalar zero = static_cast<Scalar>(0);
    out.requireWrite();

    // the summation portion
    if (inFS == ReducedElements) {
        for (dim_t colouring = 0; colouring < 2; colouring++) {
#pragma omp parallel for
            for (dim_t ey = colouring; ey < NE1; ey += 2) {
                for (dim_t ex = 0; ex < NE0; ex++) {
                    dim_t start = ex*m_order + ey*max_x*m_order;
                    const Scalar* e_in = in.getSampleDataRO(ex + ey*NE0, zero);
                    for (int qy = 0; qy < quads; qy++) {
                        for (int qx = 0; qx < quads; qx++) {
                            Scalar* n_out = out.getSampleDataRW(start + max_x*qy + qx, zero);
                            for (int comp = 0; comp < numComp; comp++) {
                                n_out[comp] += e_in[comp];
                            }
                        }
                    }
                }
            }
        }
    } else { //inFS == Elements
        for (dim_t colouring = 0; colouring < 2; colouring++) {
#pragma omp parallel for
            for (dim_t ey = colouring; ey < NE1; ey += 2) {
                for (dim_t ex = 0; ex < NE0; ex++) {
                    dim_t start = ex*m_order + ey*max_x*m_order;
                    const Scalar* e_in = in.getSampleDataRO(ex + ey*NE0, zero);
                    for (int qy = 0; qy < quads; qy++) {
                        for (int qx = 0; qx < quads; qx++) {
                            Scalar* n_out = out.getSampleDataRW(start + max_x*qy + qx, zero);
                            for (int comp = 0; comp < numComp; comp++) {
                                n_out[comp] += e_in[INDEX3(comp, qx, qy, numComp, quads)];
                            }
                        }
                    }
                }
            }
        }
    }
#ifdef ESYS_MPI
    //share and average shared edges/corners
    balanceNeighbours(out, true);
#endif
    // for every internal edge in x
#pragma omp parallel for
    for (dim_t qy = 0; qy < max_y; qy++) {
        for (dim_t qx = m_order; qx < max_x - m_order; qx += m_order) {
            Scalar* n_out = out.getSampleDataRW(qx + qy*max_x, zero);
            for (int comp = 0; comp < numComp; comp++) {
                n_out[comp] /= 2;
            }
        }
    }

    // for every internal edge in y
    const dim_t order = m_order;
#pragma omp parallel for
    for (dim_t qy = order; qy < max_y - order; qy += order) {
        for (dim_t qx = 0; qx < max_x; qx ++) {
            Scalar* n_out = out.getSampleDataRW(qx + qy*max_x, zero);
            for (int comp = 0; comp < numComp; comp++) {
                n_out[comp] /= 2;
            }
        }
    }
}

#ifdef ESYS_MPI
//private
template<typename Scalar>
void Rectangle::shareCorners(escript::Data& out, int rx, int ry) const
{
    //setup
    const int tag = 0;
    MPI_Status status;
    MPI_Request request[4];
    const int numComp = out.getDataPointSize();
    const int count = 4 * numComp;
    std::vector<Scalar> outbuf(count, 0);
    std::vector<Scalar> inbuf(count, 0);
    const int rank = m_mpiInfo->rank;
    const Scalar zero = static_cast<Scalar>(0);
    //precalc bounds so we can loop nicely, can probably be cleaned up
    const bool conds[4] = {rx && ry,
                           rx < (m_NX[0] - 1) && ry,
                           rx && ry < (m_NX[1] - 1),
                           rx < (m_NX[0] - 1) && ry < (m_NX[1] - 1)};
    const int ranks[4] = {rank-m_NX[0]-1,
                          rank-m_NX[0]+1,
                          rank+m_NX[0]-1,
                          rank+m_NX[0]+1};
    //fill everything, regardless of whether we're sharing that corner or not
    for (int y = 0; y < 2; y++) {
        for (int x = 0; x < 2; x++) {
            const Scalar* data = out.getSampleDataRO(x*(m_NN[0]-1) + y*(m_NN[1]-1)*m_NN[0], zero);
            std::copy(data, data + numComp, &outbuf[(x + 2*y)*numComp]);
        }
    }

    //share
    for (int i = 0; i < 4; i++) {
        if (conds[i]) {
            if (sizeof(Scalar) == sizeof(double)) {
                MPI_Isend(&outbuf[i], numComp, MPI_DOUBLE, ranks[i], tag,
                          m_mpiInfo->comm, &request[i]);
            } else {
                MPI_Isend(&outbuf[i], numComp, MPI_DOUBLE_COMPLEX, ranks[i],
                          tag, m_mpiInfo->comm, &request[i]);
            }
        }
    }

    //unpack
    for (int y = 0; y < 2; y++) {
        for (int x = 0; x < 2; x++) {
            int i = 2*y+x;
            if (conds[i]) {
                if (sizeof(Scalar) == sizeof(double)) {
                    MPI_Recv(&inbuf[i], numComp, MPI_DOUBLE, ranks[i], tag,
                             m_mpiInfo->comm, &status);
                } else {
                    MPI_Recv(&inbuf[i], numComp, MPI_DOUBLE_COMPLEX, ranks[i],
                             tag, m_mpiInfo->comm, &status);
                }
                Scalar* data = out.getSampleDataRW(x*(m_NN[0]-1) + y*(m_NN[1]-1)*m_NN[0], zero);
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += inbuf[i*numComp + comp];
                }
            }
        }
    }
    for (int i = 0; i < 4; i++) {
        if (conds[i]) {
            MPI_Wait(request + i, &status);
        }
    }
}

//private
template<typename Scalar>
void Rectangle::shareVertical(escript::Data& out, int rx, int ry) const
{
    const int tag = 0;
    MPI_Status status;
    const int numComp = out.getDataPointSize();
    const dim_t count = m_NN[0]*numComp;
    const int up_neighbour = m_mpiInfo->rank + m_NX[0];
    const int down_neighbour = m_mpiInfo->rank - m_NX[0];
    //allocate some space to recieve
    std::vector<Scalar> recv(count);
    //get our sources
    const Scalar zero = static_cast<Scalar>(0);
    Scalar* top = out.getSampleDataRW((m_NN[1]-1) * m_NN[0], zero);
    Scalar* bottom = out.getSampleDataRW(0, zero);
    MPI_Datatype mpiType = (sizeof(Scalar) == sizeof(double) ? MPI_DOUBLE : MPI_DOUBLE_COMPLEX);

    MPI_Request request[2];

    if (ry) {
        MPI_Isend(bottom, count, mpiType, down_neighbour, tag,
                m_mpiInfo->comm, request);
    }

    if (ry < m_NX[1] - 1) {
        MPI_Isend(top, count, mpiType, up_neighbour, tag,
            m_mpiInfo->comm, request+1);
    }

    //read down
    if (ry) {
        MPI_Recv(&recv[0], count, mpiType, down_neighbour, tag,
                m_mpiInfo->comm, &status);

        //unpack bottom
        for (dim_t i = 0; i < count; i++) {
            bottom[i] += recv[i];
        }
    }

    //read up, send up
    if (ry < m_NX[1] - 1) {
        MPI_Recv(&recv[0], count, mpiType, up_neighbour, tag,
                m_mpiInfo->comm, &status);

        //unpack up
        for (dim_t i = 0; i < count; i++) {
            top[i] += recv[i];
        }
    }

    if (ry) {
        MPI_Wait(request, &status);
    }
    if (ry < m_NX[1] - 1) {
        MPI_Wait(request+1, &status);
    }
}

//private
template<typename Scalar>
void Rectangle::shareSides(escript::Data& out, int rx, int ry) const
{
    const int tag = 0;
    MPI_Status status;
    const int numComp = out.getDataPointSize();
    const dim_t count = m_NN[1]*numComp;
    const int left_neighbour = m_mpiInfo->rank - 1;
    const int right_neighbour = m_mpiInfo->rank + 1;
    const Scalar zero = static_cast<Scalar>(0);
    //allocate some space
    std::vector<Scalar> left(count, 170);
    std::vector<Scalar> right(count, 17000);
    std::vector<Scalar> recv(count, 1700);
    MPI_Datatype mpiType = (sizeof(Scalar) == sizeof(double) ? MPI_DOUBLE : MPI_DOUBLE_COMPLEX);

    MPI_Request request[2];

    if (rx) {
        for (dim_t n = 0; n < m_NN[1]; n++) {
            index_t index = n*m_NN[0];
            const Scalar* leftData = out.getSampleDataRO(index, zero);
            std::copy(leftData, leftData + numComp, &left[n*numComp]);
        }
        MPI_Isend(&left[0], count, mpiType, left_neighbour, tag,
                m_mpiInfo->comm, request);
    }

    if (rx < m_NX[0] - 1) {
        for (dim_t n = 0; n < m_NN[1]; n++) {
            index_t index = n*m_NN[0];
            const Scalar* rightData = out.getSampleDataRO(index+m_NN[0]-1, zero);
            std::copy(rightData, rightData + numComp, &right[n*numComp]);
        }
        MPI_Isend(&right[0], count, mpiType, right_neighbour, tag,
                m_mpiInfo->comm, request+1);
    }

    //read left
    if (rx) {
        MPI_Recv(&recv[0], count, mpiType, left_neighbour, tag,
                m_mpiInfo->comm, &status);
        //unpack to left
        for (dim_t i = 0; i < m_NN[1]; i++) {
            Scalar* data = out.getSampleDataRW(i*m_NN[0], zero);
            for (int comp = 0; comp < numComp; comp++) {
                data[comp] += recv[i*numComp+comp];
            }
        }
    }

    //read right
    if (rx < m_NX[0] - 1) {
        MPI_Recv(&recv[0], count, mpiType, right_neighbour, tag,
                m_mpiInfo->comm, &status);
        //unpack to right
        for (dim_t i = 0; i < m_NN[1]; i++) {
            Scalar* data = out.getSampleDataRW((i + 1) * m_NN[0] - 1, zero);
            for (int comp = 0; comp < numComp; comp++) {
                data[comp] += recv[i*numComp+comp];
            }
        }
    }
    if (rx) {
        MPI_Wait(request, &status);
    }
    if (rx < m_NX[0] - 1) {
        MPI_Wait(request+1, &status);
    }
}
#endif //#ifdef ESYS_MPI

dim_t Rectangle::findNode(const double *coords) const
{
    const dim_t NOT_MINE = -1;
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

    //check if the point is even inside the domain
    if (x < 0 || y < 0 || x > m_length[0] || y > m_length[1])
        return NOT_MINE;

    // trim to rank reference point
    x -= m_offset[0] * m_dx[0];
    y -= m_offset[1] * m_dx[1];

    // distance in elements
    dim_t ex = (dim_t) floor((x + 0.01*m_dx[0]) / m_dx[0]);
    dim_t ey = (dim_t) floor((y + 0.01*m_dx[1]) / m_dx[1]);
    // set the min distance high enough to be outside the element plus a bit
    dim_t closest = NOT_MINE;
    double minDist = 1;
    for (int dim = 0; dim < m_numDim; dim++) {
        minDist += m_dx[dim]*m_dx[dim];
    }
    //find the closest node
    for (int dx = 0; dx < 2; dx++) {
        double xdist = x - (ex + dx)*m_dx[0];
        for (int dy = 0; dy < 2; dy++) {
            double ydist = y - (ey + dy)*m_dx[1];
            double total = xdist*xdist + ydist*ydist;
            if (total < minDist) {
                closest = (ex+dx)*m_order + (ey+dy)*m_order*m_NN[0];
                minDist = total;
            }
        }
    }
    //if this happens, we've let a dirac point slip through, which is awful
    if (closest == NOT_MINE) {
        throw SpeckleyException("Unable to map appropriate dirac point to a node,"
                " implementation problem in Rectangle::findNode()");
    }
    return closest;
}

Assembler_ptr Rectangle::createAssembler(std::string type,
        const DataMap& options) const
{
    if (type.compare("DefaultAssembler") == 0) {
        return Assembler_ptr(new DefaultAssembler2D(shared_from_this(), m_dx, m_NE, m_NN));
    } else if (type.compare("WaveAssembler") == 0) {
        return Assembler_ptr(new WaveAssembler2D(shared_from_this(), m_dx, m_NE, m_NN, options));
    }
    throw SpeckleyException("Speckley::Rectangle does not support the"
            " requested assembler");
}

bool Rectangle::probeInterpolationAcross(int fsType_source,
        const escript::AbstractDomain& domain, int fsType_target) const
{
#ifdef USE_RIPLEY
    return speckley::probeInterpolationAcross(fsType_source, domain,
            fsType_target, 2);
#else
    return false;
#endif
}

void Rectangle::interpolateAcross(escript::Data& target, const escript::Data& source) const
{
#ifdef USE_RIPLEY
    if (coupler == NULL) {
        coupler = new RipleyCoupler(this, m_dx, m_mpiInfo->rank);
    }
    coupler->interpolate(target, source);
#else
    throw SpeckleyException("Speckley::Rectangle interpolation to unsupported domain");
#endif
}

} // end of namespace speckley
