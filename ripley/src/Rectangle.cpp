
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <ripley/Rectangle.h>
#include <ripley/DefaultAssembler2D.h>
#include <ripley/LameAssembler2D.h>
#include <ripley/WaveAssembler2D.h>
#include <ripley/blocktools.h>
#include <ripley/domainhelpers.h>

#include <escript/FileWriter.h>
#include <escript/index.h>
#include <escript/Random.h>
#include <escript/Utils.h>

#ifdef ESYS_HAVE_PASO
#include <paso/SystemMatrix.h>
#endif

#ifdef ESYS_HAVE_SILO
#include <silo.h>
#ifdef ESYS_MPI
#include <pmpio.h>
#endif
#endif

#include <boost/math/special_functions/fpclassify.hpp>	// for isnan
#include <boost/scoped_array.hpp>

#include <iomanip>
#include <limits>

namespace bp = boost::python;
namespace bm = boost::math;
using escript::AbstractSystemMatrix;
using escript::FileWriter;
using escript::IOError;
using escript::NotImplementedError;
using escript::ValueError;
using std::vector;
using std::string;
using std::min;
using std::max;
using std::copy;
using std::ios;
using std::fill;

#ifdef ESYS_HAVE_NETCDF4
#include <ncFile.h>
#include <ncVar.h>
#include <ncDim.h>
using namespace netCDF;
#endif

namespace ripley {

Rectangle::Rectangle(dim_t n0, dim_t n1, double x0, double y0, double x1,
                     double y1, int d0, int d1,
                     const vector<double>& points,
                     const vector<int>& tags,
                     const TagMap& tagnamestonums) :
    RipleyDomain(2)
{
    if (static_cast<long>(n0 + 1) * static_cast<long>(n1 + 1)
            > std::numeric_limits<dim_t>::max())
        throw RipleyException("The number of elements has overflowed, this "
                "limit may be raised in future releases.");

    if (n0 <= 0 || n1 <= 0)
        throw ValueError("Number of elements in each spatial dimension "
                "must be positive");

    // ignore subdivision parameters for serial run
    if (m_mpiInfo->size == 1) {
        d0=1;
        d1=1;
    }

    bool warn = false;
    vector<int> factors;
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
                throw ValueError("Invalid number of spatial subdivisions");
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
        throw ValueError("Invalid number of spatial subdivisions");

    if (warn) {
        std::cout << "ripley.Rectangle: Warning: Automatic domain subdivision (d0=" << d0 << ", d1="
            << d1 << "). This may not be optimal!" << std::endl;
    }

    double l0 = x1-x0;
    double l1 = y1-y0;
    m_dx[0] = l0/n0;
    m_dx[1] = l1/n1;

    warn = false;
    if ((n0+1)%d0 > 0) {
        switch (getDecompositionPolicy()) {
            case DECOMP_EXPAND:
                l0 = m_dx[0]*n0; // fall through
            case DECOMP_ADD_ELEMENTS:
                n0 = (dim_t)round((float)(n0+1)/d0+0.5)*d0-1; // fall through
            case DECOMP_STRICT:
                warn = true;
                break;
        }
        // reset spacing
        m_dx[0] = l0/n0;
    }
    if ((n1+1)%d1 > 0) {
        switch (getDecompositionPolicy()) {
            case DECOMP_EXPAND:
                l1 = m_dx[1]*n1; // fall through
            case DECOMP_ADD_ELEMENTS:
                n1 = (dim_t)round((float)(n1+1)/d1+0.5)*d1-1; // fall through
            case DECOMP_STRICT:
                warn = true;
                break;
        }
        // reset spacing
        m_dx[1] = l1/n1;
    }

    if ((d0 > 1 && (n0+1)/d0<2) || (d1 > 1 && (n1+1)/d1<2))
        throw ValueError("Too few elements for the number of ranks");

    if (warn) {
        if (getDecompositionPolicy() == DECOMP_STRICT) {
            throw ValueError("Unable to decompose domain to the number of "
                    "MPI ranks without adding elements and the policy "
                    "is set to STRICT. Use setDecompositionPolicy() "
                    "to allow adding elements.");
        } else {
            std::cout << "ripley.Rectangle: Warning: Domain setup has been adjusted as follows "
                    "to allow decomposition into " << m_mpiInfo->size
                    << " MPI ranks:" << std::endl
                    << "    N0=" << n0 << ", l0=" << l0 << std::endl
                    << "    N1=" << n1 << ", l1=" << l1 << std::endl;
        }
    }
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

    for (TagMap::const_iterator i = tagnamestonums.begin();
            i != tagnamestonums.end(); i++) {
        setTagMap(i->first, i->second);
    }
    addPoints(points, tags);
}

Rectangle::Rectangle(escript::JMPI jmpi, dim_t n0, dim_t n1, double x0, double y0,
                     double x1, double y1, int d0, int d1,
                     const vector<double>& points,
                     const vector<int>& tags,
                     const TagMap& tagnamestonums) :
    RipleyDomain(2, jmpi)
{
    // This constructor delegates most work to the existing constructor logic
    // but uses a custom MPI communicator via the base class

    if (static_cast<long>(n0 + 1) * static_cast<long>(n1 + 1)
            > std::numeric_limits<dim_t>::max())
        throw RipleyException("The number of elements has overflowed, this "
                "limit may be raised in future releases.");

    if (n0 <= 0 || n1 <= 0)
        throw ValueError("Number of elements in each spatial dimension "
                "must be positive");

    // ignore subdivision parameters for serial run
    if (m_mpiInfo->size == 1) {
        d0=1;
        d1=1;
    }

    bool warn = false;
    vector<int> factors;
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
                throw ValueError("Invalid number of spatial subdivisions");
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
        throw ValueError("Invalid number of spatial subdivisions");

    if (warn) {
        std::cout << "ripley.Rectangle: Warning: Automatic domain subdivision (d0=" << d0 << ", d1="
            << d1 << "). This may not be optimal!" << std::endl;
    }

    double l0 = x1-x0;
    double l1 = y1-y0;
    m_dx[0] = l0/n0;
    m_dx[1] = l1/n1;

    warn = false;
    if ((n0+1)%d0 > 0) {
        switch (getDecompositionPolicy()) {
            case DECOMP_EXPAND:
                l0 = m_dx[0]*n0; // fall through
            case DECOMP_ADD_ELEMENTS:
                n0 = (dim_t)round((float)(n0+1)/d0+0.5)*d0-1; // fall through
            case DECOMP_STRICT:
                warn = true;
                break;
        }
        // reset spacing
        m_dx[0] = l0/n0;
    }
    if ((n1+1)%d1 > 0) {
        switch (getDecompositionPolicy()) {
            case DECOMP_EXPAND:
                l1 = m_dx[1]*n1; // fall through
            case DECOMP_ADD_ELEMENTS:
                n1 = (dim_t)round((float)(n1+1)/d1+0.5)*d1-1; // fall through
            case DECOMP_STRICT:
                warn = true;
                break;
        }
        // reset spacing
        m_dx[1] = l1/n1;
    }

    if ((d0 > 1 && (n0+1)/d0<2) || (d1 > 1 && (n1+1)/d1<2))
        throw ValueError("Too few elements for the number of ranks");

    if (warn) {
        if (getDecompositionPolicy() == DECOMP_STRICT) {
            throw ValueError("Unable to decompose domain to the number of "
                    "MPI ranks without adding elements and the policy "
                    "is set to STRICT. Use setDecompositionPolicy() "
                    "to allow adding elements.");
        } else {
            std::cout << "ripley.Rectangle: Warning: Domain setup has been adjusted as follows "
                    "to allow decomposition into " << m_mpiInfo->size
                    << " MPI ranks:" << std::endl
                    << "    N0=" << n0 << ", l0=" << l0 << std::endl
                    << "    N1=" << n1 << ", l1=" << l1 << std::endl;
        }
    }
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

    for (TagMap::const_iterator i = tagnamestonums.begin();
            i != tagnamestonums.end(); i++) {
        setTagMap(i->first, i->second);
    }
    addPoints(points, tags);
}

Rectangle::~Rectangle()
{
}

string Rectangle::getDescription() const
{
    return "ripley::Rectangle";
}

bool Rectangle::operator==(const escript::AbstractDomain& other) const
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
            const ReaderParameters& params) const
{
#ifdef ESYS_HAVE_NETCDF4
    // check destination function space
    dim_t myN0, myN1;
    if (out.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NN[0];
        myN1 = m_NN[1];
    } else if (out.getFunctionSpace().getTypeCode() == Elements ||
                out.getFunctionSpace().getTypeCode() == ReducedElements) {
        myN0 = m_NE[0];
        myN1 = m_NE[1];
    } else
        throw ValueError("readNcGrid(): invalid function space for output data object");

    if (params.first.size() != 2)
        throw ValueError("readNcGrid(): argument 'first' must have 2 entries");

    if (params.numValues.size() != 2)
        throw ValueError("readNcGrid(): argument 'numValues' must have 2 entries");

    if (params.multiplier.size() != 2)
        throw ValueError("readNcGrid(): argument 'multiplier' must have 2 entries");
    for (size_t i=0; i<params.multiplier.size(); i++)
        if (params.multiplier[i]<1)
            throw ValueError("readNcGrid(): all multipliers must be positive");
    if (params.reverse.size() != 2)
        throw ValueError("readNcGrid(): argument 'reverse' must have 2 entries");

    NcFile f;
    if (!escript::openNcFile(f, filename))
    {
        throw RipleyException("readNcGrid(): cannot open file");
    }

    NcVar var = f.getVar(varname.c_str());
    if (var.isNull())
        throw RipleyException("readNcGrid(): invalid variable name");

    // TODO: rank>0 data support
    const int numComp = out.getDataPointSize();
    if (numComp > 1)
        throw NotImplementedError("readNcGrid(): only scalar data supported");

    const int dims = var.getDimCount();
    vector<long> edges(dims);
    std::vector< NcDim > vard=var.getDims();
    for (size_t i=0;i<vard.size();++i)
    {
        edges[i]=vard[i].getSize();
    }

    // is this a slice of the data object (dims!=2)?
    // note the expected ordering of edges (as in numpy: y,x)
    if ( (dims==2 && (params.numValues[1] > edges[0] || params.numValues[0] > edges[1]))
            || (dims==1 && params.numValues[1]>1) ) {
        throw IOError("readNcGrid(): not enough data in file");
    }

    // check if this rank contributes anything
    if (params.first[0] >= m_offset[0]+myN0 ||
            params.first[0]+params.numValues[0]*params.multiplier[0] <= m_offset[0] ||
            params.first[1] >= m_offset[1]+myN1 ||
            params.first[1]+params.numValues[1]*params.multiplier[1] <= m_offset[1])
        return;

    // now determine how much this rank has to write

    // first coordinates in data object to write to
    const dim_t first0 = max(dim_t(0), params.first[0]-m_offset[0]);
    const dim_t first1 = max(dim_t(0), params.first[1]-m_offset[1]);
    // indices to first value in file (not accounting for reverse yet)
    dim_t idx0 = max(dim_t(0), m_offset[0]-params.first[0]);
    dim_t idx1 = max(dim_t(0), m_offset[1]-params.first[1]);
    // number of values to read
    const dim_t num0 = min(params.numValues[0]-idx0, myN0-first0);
    const dim_t num1 = min(params.numValues[1]-idx1, myN1-first1);

    // make sure we read the right block if going backwards through file
    if (params.reverse[0])
        idx0 = edges[dims-1]-num0-idx0;
    if (dims>1 && params.reverse[1])
        idx1 = edges[dims-2]-num1-idx1;

    vector<double> values(num0*num1);
    vector<size_t> startindex;
    vector<size_t> counts;
    if (dims==2) {
//         var->set_cur(idx1, idx0);
//         var->get(&values[0], num1, num0);
        startindex.push_back(idx1);
        startindex.push_back(idx0);
        counts.push_back(num1);
        counts.push_back(num0);
        var.getVar(startindex, counts, &values[0]);
    } else {
//         var->set_cur(idx0);
//         var->get(&values[0], num0);
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
    throw RipleyException("readNcGrid(): not compiled with netCDF support");
#endif
}

void Rectangle::readBinaryGrid(escript::Data& out, string filename,
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
            throw ValueError("readBinaryGrid(): invalid or unsupported datatype");
    }
}

void Rectangle::readBinaryGridFromZipped(escript::Data& out, string filename,
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
            throw ValueError("readBinaryGridFromZipped(): invalid or unsupported datatype");
    }
#else
    throw ValueError("readBinaryGridFromZipped(): not compiled with zip support");
#endif
}

template<typename ValueType>
void Rectangle::readBinaryGridImpl(escript::Data& out, const string& filename,
                                   const ReaderParameters& params) const
{
    // check destination function space
    dim_t myN0, myN1;
    if (out.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NN[0];
        myN1 = m_NN[1];
    } else if (out.getFunctionSpace().getTypeCode() == Elements ||
                out.getFunctionSpace().getTypeCode() == ReducedElements) {
        myN0 = m_NE[0];
        myN1 = m_NE[1];
    } else
        throw ValueError("readBinaryGrid(): invalid function space for output data object");

    if (params.first.size() != 2)
        throw ValueError("readBinaryGrid(): argument 'first' must have 2 entries");

    if (params.numValues.size() != 2)
        throw ValueError("readBinaryGrid(): argument 'numValues' must have 2 entries");

    if (params.multiplier.size() != 2)
        throw ValueError("readBinaryGrid(): argument 'multiplier' must have 2 entries");
    for (size_t i=0; i<params.multiplier.size(); i++)
        if (params.multiplier[i]<1)
            throw ValueError("readBinaryGrid(): all multipliers must be positive");
    if (params.reverse[0] != 0 || params.reverse[1] != 0)
        throw NotImplementedError("readBinaryGrid(): reversing not supported yet");

    // check file existence and size
    std::ifstream f(filename.c_str(), std::ifstream::binary);
    if (f.fail()) {
        throw IOError("readBinaryGrid(): cannot open file " + filename);
    }
    f.seekg(0, ios::end);
    const int numComp = out.getDataPointSize();
    const dim_t filesize = f.tellg();
    const dim_t reqsize = params.numValues[0]*params.numValues[1]*numComp*sizeof(ValueType);
    if (filesize < reqsize) {
        f.close();
        throw IOError("readBinaryGrid(): not enough data in file");
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
    const dim_t first0 = max(dim_t(0), params.first[0]-m_offset[0]);
    const dim_t first1 = max(dim_t(0), params.first[1]-m_offset[1]);
    // indices to first value in file
    const dim_t idx0 = max(dim_t(0), (m_offset[0]/params.multiplier[0])-params.first[0]);
    const dim_t idx1 = max(dim_t(0), (m_offset[1]/params.multiplier[1])-params.first[1]);
    // if restX > 0 the first value in the respective dimension has been
    // written restX times already in a previous rank so this rank only
    // contributes (multiplier-rank) copies of that value
    const dim_t rest0 = m_offset[0]%params.multiplier[0];
    const dim_t rest1 = m_offset[1]%params.multiplier[1];
    // number of values to read
    const dim_t num0 = min(params.numValues[0]-idx0, myN0-first0);
    const dim_t num1 = min(params.numValues[1]-idx1, myN1-first1);

    out.requireWrite();
    vector<ValueType> values(num0*numComp);
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
                    const dim_t dataIndex = dataX+dataY*myN0;
                    double* dest = out.getSampleDataRW(dataIndex);
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
}

#ifdef ESYS_HAVE_BOOST_IO
template<typename ValueType>
void Rectangle::readBinaryGridZippedImpl(escript::Data& out, const string& filename,
                                   const ReaderParameters& params) const
{
    // check destination function space
    dim_t myN0, myN1;
    if (out.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NN[0];
        myN1 = m_NN[1];
    } else if (out.getFunctionSpace().getTypeCode() == Elements ||
                out.getFunctionSpace().getTypeCode() == ReducedElements) {
        myN0 = m_NE[0];
        myN1 = m_NE[1];
    } else
        throw ValueError("readBinaryGrid(): invalid function space for output data object");

    // check file existence and size
    std::ifstream f(filename.c_str(), std::ifstream::binary);
    if (f.fail()) {
        throw IOError("readBinaryGridFromZipped(): cannot open file" + filename);
    }
    f.seekg(0, ios::end);
    const int numComp = out.getDataPointSize();
    dim_t filesize = f.tellg();
    f.seekg(0, ios::beg);
    vector<char> compressed(filesize);
    f.read((char*)&compressed[0], filesize);
    f.close();
    vector<char> decompressed = unzip(compressed);
    filesize = decompressed.size();
    const dim_t reqsize = params.numValues[0]*params.numValues[1]*numComp*sizeof(ValueType);
    if (filesize < reqsize) {
        throw IOError("readBinaryGridFromZipped(): not enough data in file");
    }

    // check if this rank contributes anything
    if (params.first[0] >= m_offset[0]+myN0 ||
            params.first[0]+params.numValues[0] <= m_offset[0] ||
            params.first[1] >= m_offset[1]+myN1 ||
            params.first[1]+params.numValues[1] <= m_offset[1]) {
        f.close();
        return;
    }

    // now determine how much this rank has to write

    // first coordinates in data object to write to
    const dim_t first0 = max(dim_t(0), params.first[0]-m_offset[0]);
    const dim_t first1 = max(dim_t(0), params.first[1]-m_offset[1]);
    // indices to first value in file
    const dim_t idx0 = max(dim_t(0), m_offset[0]-params.first[0]);
    const dim_t idx1 = max(dim_t(0), m_offset[1]-params.first[1]);
    // number of values to read
    const dim_t num0 = min(params.numValues[0]-idx0, myN0-first0);
    const dim_t num1 = min(params.numValues[1]-idx1, myN1-first1);

    out.requireWrite();
    vector<ValueType> values(num0*numComp);
    const int dpp = out.getNumDataPointsPerSample();

    for (dim_t y=0; y<num1; y++) {
        const dim_t fileofs = numComp*(idx0+(idx1+y)*params.numValues[0]);
            memcpy((char*)&values[0], (char*)&decompressed[fileofs*sizeof(ValueType)], num0*numComp*sizeof(ValueType));
        for (dim_t x=0; x<num0; x++) {
            const dim_t baseIndex = first0+x*params.multiplier[0]
                                    +(first1+y*params.multiplier[1])*myN0;
            for (int m1=0; m1<params.multiplier[1]; m1++) {
                for (int m0=0; m0<params.multiplier[0]; m0++) {
                    const dim_t dataIndex = baseIndex+m0+m1*myN0;
                    double* dest = out.getSampleDataRW(dataIndex);
                    for (int c=0; c<numComp; c++) {
                        ValueType val = values[x*numComp+c];

                        if (params.byteOrder != BYTEORDER_NATIVE) {
                            char* cval = reinterpret_cast<char*>(&val);
                            // this will alter val!!
                            byte_swap32(cval);
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
}
#endif

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
            throw ValueError("writeBinaryGrid(): invalid or unsupported datatype");
    }
}

template<typename ValueType>
void Rectangle::writeBinaryGridImpl(const escript::Data& in,
                                    const string& filename, int byteOrder) const
{
    // check function space and determine number of points
    dim_t myN0, myN1;
    dim_t totalN0, totalN1;
    dim_t offset0, offset1;
    if (in.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NN[0];
        myN1 = m_NN[1];
        totalN0 = m_gNE[0]+1;
        totalN1 = m_gNE[1]+1;
        offset0 = m_offset[0];
        offset1 = m_offset[1];
    } else if (in.getFunctionSpace().getTypeCode() == DegreesOfFreedom ||
            in.getFunctionSpace().getTypeCode() == ReducedDegreesOfFreedom) {
        myN0 = (m_gNE[0]+1)/m_NX[0];
        myN1 = (m_gNE[1]+1)/m_NX[1];
        totalN0 = m_gNE[0]+1;
        totalN1 = m_gNE[1]+1;
        offset0 = (m_offset[0]>0 ? m_offset[0]+1 : 0);
        offset1 = (m_offset[1]>0 ? m_offset[1]+1 : 0);
    } else if (in.getFunctionSpace().getTypeCode() == Elements ||
                in.getFunctionSpace().getTypeCode() == ReducedElements) {
        myN0 = m_NE[0];
        myN1 = m_NE[1];
        totalN0 = m_gNE[0];
        totalN1 = m_gNE[1];
        offset0 = m_offset[0];
        offset1 = m_offset[1];
    } else
        throw ValueError("writeBinaryGrid(): unsupported function space");

    const int numComp = in.getDataPointSize();
    const int dpp = in.getNumDataPointsPerSample();

    if (numComp > 1 || dpp > 1)
        throw NotImplementedError("writeBinaryGrid(): only scalar, single-value data supported");

    const dim_t fileSize = sizeof(ValueType)*numComp*dpp*totalN0*totalN1;

    // from here on we know that each sample consists of one value
    FileWriter fw(m_mpiInfo->comm);
    fw.openFile(filename, fileSize, true);
    MPIBarrier();

    for (index_t y=0; y<myN1; y++) {
        const dim_t fileofs = (offset0+(offset1+y)*totalN0)*sizeof(ValueType);
        std::ostringstream oss;

        for (index_t x=0; x<myN0; x++) {
            const double* sample = in.getSampleDataRO(y*myN0+x);
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
    throw NotImplementedError("write: not supported");
}

void Rectangle::dump(const string& fileName) const
{
#ifdef ESYS_HAVE_SILO
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
        throw IOError("dump: Could not create Silo file");

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
    boost::scoped_ptr<double> x(new double[NN0]);
    boost::scoped_ptr<double> y(new double[NN1]);
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

    // casting to int!!
    vector<int> dims(m_NN, m_NN+2);

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
        vector<string> tempstrings;
        vector<char*> names;
        for (dim_t i=0; i<m_mpiInfo->size; i++) {
            std::stringstream path;
            path << "/block" << std::setw(4) << std::setfill('0') << std::right << i << "/mesh";
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
    throw RipleyException("dump: no Silo support");
#endif
}

const dim_t* Rectangle::borrowSampleReferenceIDs(int fsType) const
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
        case Points:
            return &m_diracPointNodeIDs[0];
        default:
            break;
    }

    std::stringstream msg;
    msg << "borrowSampleReferenceIDs: invalid function space type " << fsType;
    throw ValueError(msg.str());
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

    std::stringstream msg;
    msg << "ownSample: invalid function space type " << fsType;
    throw ValueError(msg.str());
}

RankVector Rectangle::getOwnerVector(int fsType) const
{
    RankVector owner;
    const int rank = m_mpiInfo->rank;

    if (fsType == Elements || fsType == ReducedElements) {
        owner.assign(getNumElements(), rank);
        if (m_faceCount[0] == 0) {
            owner[0]=(m_faceCount[2]==0 ? rank-m_NX[0]-1 : rank-1);
            for (dim_t i=1; i<m_NE[1]; i++)
                owner[i*m_NE[0]] = rank-1;
        }
        if (m_faceCount[2]==0) {
            const int first=(m_faceCount[0]==0 ? 1 : 0);
            for (dim_t i=first; i<m_NE[0]; i++)
                owner[i] = rank-m_NX[0];
        }

    } else if (fsType == FaceElements || fsType == ReducedFaceElements) {
        owner.assign(getNumFaceElements(), rank);
        if (m_faceCount[0] == 0) {
            if (m_faceCount[2] > 0)
                owner[m_faceCount[1]] = rank-1;
            if (m_faceCount[3] > 0)
                owner[m_faceCount[1]+m_faceCount[2]] = rank-1;
        }
        if (m_faceCount[2] == 0) {
            if (m_faceCount[0] > 0)
                owner[0] = rank-m_NX[0];
            if (m_faceCount[1] > 0)
                owner[m_faceCount[0]] = rank-m_NX[0];
        }

    } else {
        throw ValueError("getOwnerVector: only valid for element types");
    }

    return owner;
}

void Rectangle::setToNormal(escript::Data& out) const
{
    const dim_t NE0=m_NE[0];
    const dim_t NE1=m_NE[1];
    if (out.getFunctionSpace().getTypeCode() == FaceElements) {
        out.requireWrite();
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < NE1; ++k1) {
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
                for (index_t k1 = 0; k1 < NE1; ++k1) {
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
                for (index_t k0 = 0; k0 < NE0; ++k0) {
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
                for (index_t k0 = 0; k0 < NE0; ++k0) {
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
                for (index_t k1 = 0; k1 < NE1; ++k1) {
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k1);
                    *o++ = -1.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < NE1; ++k1) {
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k1);
                    *o++ = 1.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < NE0; ++k0) {
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k0);
                    *o++ = 0.;
                    *o = -1.;
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < NE0; ++k0) {
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k0);
                    *o++ = 0.;
                    *o = 1.;
                }
            }
        } // end of parallel section

    } else {
        std::stringstream msg;
        msg << "setToNormal: invalid function space type "
            << out.getFunctionSpace().getTypeCode();
        throw ValueError(msg.str());
    }
}

void Rectangle::setToSize(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == Elements
            || out.getFunctionSpace().getTypeCode() == ReducedElements) {
        out.requireWrite();
        const dim_t numQuad = out.getNumDataPointsPerSample();
        const dim_t numElements = getNumElements();
        const double size=sqrt(m_dx[0]*m_dx[0]+m_dx[1]*m_dx[1]);
#pragma omp parallel for
        for (index_t k = 0; k < numElements; ++k) {
            double* o = out.getSampleDataRW(k);
            fill(o, o+numQuad, size);
        }
    } else if (out.getFunctionSpace().getTypeCode() == FaceElements
            || out.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        out.requireWrite();
        const dim_t numQuad=out.getNumDataPointsPerSample();
        const dim_t NE0 = m_NE[0];
        const dim_t NE1 = m_NE[1];
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < NE1; ++k1) {
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k1);
                    fill(o, o+numQuad, m_dx[1]);
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < NE1; ++k1) {
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k1);
                    fill(o, o+numQuad, m_dx[1]);
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < NE0; ++k0) {
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k0);
                    fill(o, o+numQuad, m_dx[0]);
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < NE0; ++k0) {
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k0);
                    fill(o, o+numQuad, m_dx[0]);
                }
            }
        } // end of parallel section

    } else {
        std::stringstream msg;
        msg << "setToSize: invalid function space type "
            << out.getFunctionSpace().getTypeCode();
        throw ValueError(msg.str());
    }
}

void Rectangle::Print_Mesh_Info(const bool full) const
{
    RipleyDomain::Print_Mesh_Info(full);
    if (full) {
        std::cout << "     Id  Coordinates" << std::endl;
        std::cout.precision(15);
        std::cout.setf(ios::scientific, ios::floatfield);
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
        throw ValueError("setToX: Invalid Data object shape");
    if (!arg.numSamplesEqual(1, getNumNodes()))
        throw ValueError("setToX: Illegal number of samples in Data object");

    const dim_t NN0 = m_NN[0];
    const dim_t NN1 = m_NN[1];
    arg.requireWrite();
#pragma omp parallel for
    for (dim_t i1 = 0; i1 < NN1; i1++) {
        for (dim_t i0 = 0; i0 < NN0; i0++) {
            double* point = arg.getSampleDataRW(i0+NN0*i1);
            point[0] = getLocalCoordinate(i0, 0);
            point[1] = getLocalCoordinate(i1, 1);
        }
    }
}

//protected
void Rectangle::assembleGradient(escript::Data& out,
                                 const escript::Data& in) const
{
    if (out.isComplex() && in.isComplex())
        assembleGradientImpl<cplx_t>(out, in);
    else if (!out.isComplex() && !in.isComplex())
        assembleGradientImpl<real_t>(out, in);
    else
        throw ValueError("Gradient: input & output complexity must match.");
}

//protected
template<typename Scalar>
void Rectangle::assembleGradientImpl(escript::Data& out,
                                     const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    const dim_t NE0 = m_NE[0];
    const dim_t NE1 = m_NE[1];
    const double cx0 = .21132486540518711775/m_dx[0];
    const double cx1 = .78867513459481288225/m_dx[0];
    const double cx2 = 1./m_dx[0];
    const double cy0 = .21132486540518711775/m_dx[1];
    const double cy1 = .78867513459481288225/m_dx[1];
    const double cy2 = 1./m_dx[1];
    const Scalar zero = static_cast<Scalar>(0);

    if (out.getFunctionSpace().getTypeCode() == Elements) {
        out.requireWrite();
#pragma omp parallel
        {
            vector<Scalar> f_00(numComp, zero);
            vector<Scalar> f_01(numComp, zero);
            vector<Scalar> f_10(numComp, zero);
            vector<Scalar> f_11(numComp, zero);
#pragma omp for
            for (index_t k1 = 0; k1 < NE1; ++k1) {
                for (index_t k0 = 0; k0 < NE0; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,k1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,k1+1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,k1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,k1+1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(INDEX2(k0,k1,m_NE[0]), zero);
                    for (index_t i = 0; i < numComp; ++i) {
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
            vector<Scalar> f_00(numComp, zero);
            vector<Scalar> f_01(numComp, zero);
            vector<Scalar> f_10(numComp, zero);
            vector<Scalar> f_11(numComp, zero);
#pragma omp for
            for (index_t k1 = 0; k1 < NE1; ++k1) {
                for (index_t k0 = 0; k0 < NE0; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,k1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,k1+1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,k1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,k1+1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(INDEX2(k0,k1,m_NE[0]), zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i] + f_11[i] - f_00[i] - f_01[i])*cx2 * 0.5;
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i] + f_11[i] - f_00[i] - f_10[i])*cy2 * 0.5;
                    } // end of component loop i
                } // end of k0 loop
            } // end of k1 loop
        } // end of parallel section
    } else if (out.getFunctionSpace().getTypeCode() == FaceElements) {
        out.requireWrite();
#pragma omp parallel
        {
            vector<Scalar> f_00(numComp, zero);
            vector<Scalar> f_01(numComp, zero);
            vector<Scalar> f_10(numComp, zero);
            vector<Scalar> f_11(numComp, zero);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < NE1; ++k1) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(0,k1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(0,k1+1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(1,k1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(1,k1+1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(m_faceOffset[0]+k1, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx1 + (f_11[i]-f_01[i])*cx0;
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy2;
                        o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx0 + (f_11[i]-f_01[i])*cx1;
                        o[INDEX3(i,1,1,numComp,2)] = (f_01[i]-f_00[i])*cy2;
                    } // end of component loop i
                } // end of k1 loop
            } // end of face 0
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < NE1; ++k1) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(m_NN[0]-2,k1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(m_NN[0]-2,k1+1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1+1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(m_faceOffset[1]+k1, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx1 + (f_11[i]-f_01[i])*cx0;
                        o[INDEX3(i,1,0,numComp,2)] = (f_11[i]-f_10[i])*cy2;
                        o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx0 + (f_11[i]-f_01[i])*cx1;
                        o[INDEX3(i,1,1,numComp,2)] = (f_11[i]-f_10[i])*cy2;
                    } // end of component loop i
                } // end of k1 loop
            } // end of face 1
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < NE0; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,0, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,0, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(m_faceOffset[2]+k0, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx2;
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy1 + (f_11[i]-f_10[i])*cy0;
                        o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx2;
                        o[INDEX3(i,1,1,numComp,2)] = (f_01[i]-f_00[i])*cy0 + (f_11[i]-f_10[i])*cy1;
                    } // end of component loop i
                } // end of k0 loop
            } // end of face 2
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < NE0; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,m_NN[1]-2, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,m_NN[1]-1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,m_NN[1]-2, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,m_NN[1]-1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(m_faceOffset[3]+k0, zero);
                    for (index_t i = 0; i < numComp; ++i) {
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
            vector<Scalar> f_00(numComp, zero);
            vector<Scalar> f_01(numComp, zero);
            vector<Scalar> f_10(numComp, zero);
            vector<Scalar> f_11(numComp, zero);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < NE1; ++k1) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(0,k1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(0,k1+1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(1,k1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(1,k1+1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(m_faceOffset[0]+k1, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i] + f_11[i] - f_00[i] - f_01[i])*cx2 * 0.5;
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy2;
                    } // end of component loop i
                } // end of k1 loop
            } // end of face 0
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < NE1; ++k1) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(m_NN[0]-2,k1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(m_NN[0]-2,k1+1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1+1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(m_faceOffset[1]+k1, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i] + f_11[i] - f_00[i] - f_01[i])*cx2 * 0.5;
                        o[INDEX3(i,1,0,numComp,2)] = (f_11[i]-f_10[i])*cy2;
                    } // end of component loop i
                } // end of k1 loop
            } // end of face 1
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < NE0; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,0, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,0, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(m_faceOffset[2]+k0, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx2;
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i] + f_11[i] - f_00[i] - f_10[i])*cy2 * 0.5;
                    } // end of component loop i
                } // end of k0 loop
            } // end of face 2
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < NE0; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,m_NN[1]-2, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,m_NN[1]-1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,m_NN[1]-2, m_NN[0]), zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,m_NN[1]-1, m_NN[0]), zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(m_faceOffset[3]+k0, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_11[i]-f_01[i])*cx2;
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i] + f_11[i] - f_00[i] - f_10[i])*cy2 * 0.5;
                    } // end of component loop i
                } // end of k0 loop
            } // end of face 3
        } // end of parallel section
    }
}

// instantiate our two supported versions
template
void Rectangle::assembleGradientImpl<real_t>(escript::Data& out,
                                             const escript::Data& in) const;

template
void Rectangle::assembleGradientImpl<cplx_t>(escript::Data& out,
                                             const escript::Data& in) const;

//protected
void Rectangle::assembleIntegrate(vector<real_t>& integrals,
                                  const escript::Data& arg) const
{
    assembleIntegrateImpl<real_t>(integrals, arg);
}

//protected
void Rectangle::assembleIntegrate(vector<cplx_t>& integrals,
                                  const escript::Data& arg) const
{
    assembleIntegrateImpl<cplx_t>(integrals, arg);
}

//private
template<typename Scalar>
void Rectangle::assembleIntegrateImpl(vector<Scalar>& integrals,
                                      const escript::Data& arg) const
{
    const dim_t numComp = arg.getDataPointSize();
    const index_t left = getFirstInDim(0);
    const index_t bottom = getFirstInDim(1);
    const int fs = arg.getFunctionSpace().getTypeCode();
    const Scalar zero = static_cast<Scalar>(0);

    if(fs == Points ) {
        for (index_t k1 = 0; k1 < m_diracPoints.size(); k1++) { //only for this rank
            const Scalar* f  = arg.getSampleDataRO(k1, zero);
            for (index_t i = 0; i < numComp; ++i) {
                integrals[i]+=f[i];
            }
        }

//     bool HavePointData = arg.getFunctionSpace().getTypeCode() == Points;

// #ifdef ESYS_MPI
//     if(HavePointData && escript::getMPIRankWorld() == 0) {
// #else
//     if(HavePointData) {
// #endif
//         integrals[0] += arg.getNumberOfTaggedValues();

    } else if (fs == Elements && arg.actsExpanded()) {
#pragma omp parallel
        {
            vector<Scalar> int_local(numComp, zero);
            const real_t w = m_dx[0]*m_dx[1]/4.;
#pragma omp for nowait
            for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                    const Scalar* f = arg.getSampleDataRO(INDEX2(k0, k1, m_NE[0]), zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        const Scalar f0 = f[INDEX2(i,0,numComp)];
                        const Scalar f1 = f[INDEX2(i,1,numComp)];
                        const Scalar f2 = f[INDEX2(i,2,numComp)];
                        const Scalar f3 = f[INDEX2(i,3,numComp)];
                        int_local[i] += (f0+f1+f2+f3)*w;
                    }  // end of component loop i
                } // end of k0 loop
            } // end of k1 loop
#pragma omp critical
            for (index_t i=0; i<numComp; i++)
                integrals[i] += int_local[i];
        } // end of parallel section

    } else if (fs==ReducedElements || (fs==Elements && !arg.actsExpanded())) {
        const real_t w = m_dx[0]*m_dx[1];
#pragma omp parallel
        {
            vector<Scalar> int_local(numComp, 0);
#pragma omp for nowait
            for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                    const Scalar* f = arg.getSampleDataRO(INDEX2(k0, k1, m_NE[0]), zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        int_local[i] += f[i]*w;
                    }
                }
            }
#pragma omp critical
            for (index_t i = 0; i < numComp; i++)
                integrals[i] += int_local[i];
        } // end of parallel section

    } else if (fs == FaceElements && arg.actsExpanded()) {
#pragma omp parallel
        {
            vector<Scalar> int_local(numComp, zero);
            const real_t w0 = m_dx[0]/2.;
            const real_t w1 = m_dx[1]/2.;
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    const Scalar* f = arg.getSampleDataRO(m_faceOffset[0]+k1, zero);
                    for (index_t i=0; i < numComp; ++i) {
                        const Scalar f0 = f[INDEX2(i,0,numComp)];
                        const Scalar f1 = f[INDEX2(i,1,numComp)];
                        int_local[i] += (f0+f1)*w1;
                    }  // end of component loop i
                } // end of k1 loop
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    const Scalar* f = arg.getSampleDataRO(m_faceOffset[1]+k1, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        const Scalar f0 = f[INDEX2(i,0,numComp)];
                        const Scalar f1 = f[INDEX2(i,1,numComp)];
                        int_local[i] += (f0+f1)*w1;
                    }  // end of component loop i
                } // end of k1 loop
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                    const Scalar* f = arg.getSampleDataRO(m_faceOffset[2]+k0, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        const Scalar f0 = f[INDEX2(i,0,numComp)];
                        const Scalar f1 = f[INDEX2(i,1,numComp)];
                        int_local[i] += (f0+f1)*w0;
                    }  // end of component loop i
                } // end of k0 loop
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                    const Scalar* f = arg.getSampleDataRO(m_faceOffset[3]+k0, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        const Scalar f0 = f[INDEX2(i,0,numComp)];
                        const Scalar f1 = f[INDEX2(i,1,numComp)];
                        int_local[i] += (f0+f1)*w0;
                    }  // end of component loop i
                } // end of k0 loop
            }
#pragma omp critical
            for (index_t i = 0; i < numComp; i++)
                integrals[i] += int_local[i];
        } // end of parallel section

    } else if (fs==ReducedFaceElements || (fs==FaceElements && !arg.actsExpanded())) {
#pragma omp parallel
        {
            vector<Scalar> int_local(numComp, 0);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    const Scalar* f = arg.getSampleDataRO(m_faceOffset[0]+k1, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        int_local[i] += f[i]*m_dx[1];
                    }
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    const Scalar* f = arg.getSampleDataRO(m_faceOffset[1]+k1, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        int_local[i] += f[i]*m_dx[1];
                    }
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                    const Scalar* f = arg.getSampleDataRO(m_faceOffset[2]+k0, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        int_local[i] += f[i]*m_dx[0];
                    }
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                    const Scalar* f = arg.getSampleDataRO(m_faceOffset[3]+k0, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        int_local[i] += f[i]*m_dx[0];
                    }
                }
            }

#pragma omp critical
            for (index_t i = 0; i < numComp; i++)
                integrals[i] += int_local[i];
        } // end of parallel section
    } // function space selector
}

//protected
IndexVector Rectangle::getDiagonalIndices(bool upperOnly) const
{
    IndexVector ret;
    // only store non-negative indices if requested
    if (upperOnly)
        ret.resize(5);
    else
        ret.resize(9);
    const dim_t nDOF0 = getNumDOFInAxis(0);
    size_t idx = 0;
    for (int i1=-1; i1<2; i1++) {
        for (int i0=-1; i0<2; i0++) {
            const int index = i1*nDOF0 + i0;
            if (!upperOnly || index >= 0)
                ret[idx++] = index;
        }
    }

    return ret;
}

//protected
void Rectangle::nodesToDOF(escript::Data& out, const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    out.requireWrite();

    const index_t left = getFirstInDim(0);
    const index_t bottom = getFirstInDim(1);
    const dim_t nDOF0 = getNumDOFInAxis(0);
    const dim_t nDOF1 = getNumDOFInAxis(1);
#pragma omp parallel for
    for (index_t i=0; i<nDOF1; i++) {
        for (index_t j=0; j<nDOF0; j++) {
            const index_t n=j+left+(i+bottom)*m_NN[0];
            const double* src=in.getSampleDataRO(n);
            copy(src, src+numComp, out.getSampleDataRW(j+i*nDOF0));
        }
    }
}

#ifdef ESYS_HAVE_TRILINOS
//protected
esys_trilinos::const_TrilinosGraph_ptr Rectangle::getTrilinosGraph() const
{
    if (m_graph.is_null()) {
        m_graph = createTrilinosGraph(m_dofId, m_nodeId);
    }
    return m_graph;
}
#endif

#ifdef ESYS_HAVE_PASO
//protected
paso::SystemMatrixPattern_ptr Rectangle::getPasoMatrixPattern(
                                                    bool reducedRowOrder,
                                                    bool reducedColOrder) const
{
    if (m_pattern.get())
        return m_pattern;

    // first call - create pattern, then return
    paso::Connector_ptr conn(getPasoConnector());
    const dim_t numDOF = getNumDOF();
    const dim_t numShared = conn->send->numSharedComponents;
    const dim_t numNeighbours = conn->send->neighbour.size();
    const std::vector<index_t>& offsetInShared(conn->send->offsetInShared);
    const index_t* sendShared = conn->send->shared;

    // these are for the couple blocks
    vector<IndexVector> colIndices(numDOF);
    vector<IndexVector> rowIndices(numShared);

    for (dim_t i=0; i<numNeighbours; i++) {
        const dim_t start = offsetInShared[i];
        const dim_t end = offsetInShared[i+1];
        for (dim_t j = start; j < end; j++) {
            if (j > start)
                doublyLink(colIndices, rowIndices, sendShared[j-1], j);
            doublyLink(colIndices, rowIndices, sendShared[j], j);
            if (j < end-1)
                doublyLink(colIndices, rowIndices, sendShared[j+1], j);
        }
    }
#pragma omp parallel for
    for (dim_t i = 0; i < numShared; i++) {
        sort(rowIndices[i].begin(), rowIndices[i].end());
    }

    // create main and couple blocks
    paso::Pattern_ptr mainPattern = createPasoPattern(getConnections(), numDOF);
    paso::Pattern_ptr colPattern = createPasoPattern(colIndices, numShared);
    paso::Pattern_ptr rowPattern = createPasoPattern(rowIndices, numDOF);

    // allocate paso distribution
    escript::Distribution_ptr distribution(new escript::Distribution(
                                               m_mpiInfo, m_nodeDistribution));

    // finally create the system matrix pattern
    m_pattern.reset(new paso::SystemMatrixPattern(MATRIX_FORMAT_DEFAULT,
            distribution, distribution, mainPattern, colPattern, rowPattern,
            conn, conn));
    return m_pattern;
}
#endif // ESYS_HAVE_PASO

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
    try {
        m_nodeId.resize(getNumNodes());
        m_dofId.resize(numDOF);
        m_elementId.resize(getNumElements());
    } catch (const std::length_error& le) {
        throw RipleyException("The system does not have sufficient memory for a domain of this size.");
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

    const dim_t NFE = getNumFaceElements();
    m_faceId.resize(NFE);

    const index_t left = getFirstInDim(0);
    const index_t bottom = getFirstInDim(1);
    const dim_t nDOF0 = getNumDOFInAxis(0);
    const dim_t nDOF1 = getNumDOFInAxis(1);
    const dim_t NE0 = m_NE[0];
    const dim_t NE1 = m_NE[1];

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
        for (dim_t i1=0; i1<NE1; i1++) {
            for (dim_t i0=0; i0<NE0; i0++) {
                m_elementId[i0+i1*NE0]=(m_offset[1]+i1)*m_gNE[0]+m_offset[0]+i0;
            }
        }

        // face elements
#pragma omp for
        for (dim_t k=0; k<NFE; k++)
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

    populateDofMap();
}

//private
vector<IndexVector> Rectangle::getConnections(bool includeShared) const
{
    // returns a vector v of size numDOF where v[i] is a vector with indices
    // of DOFs connected to i (up to 9 in 2D).
    // In other words this method returns the occupied (local) matrix columns
    // for all (local) matrix rows.
    // If includeShared==true then connections to non-owned DOFs are also
    // returned (i.e. indices of the column couplings)
    const dim_t nDOF0 = getNumDOFInAxis(0);
    const dim_t nDOF1 = getNumDOFInAxis(1);
    const dim_t numMatrixRows = nDOF0*nDOF1;
    vector<IndexVector> indices(numMatrixRows);

    if (includeShared) {
        const index_t left = getFirstInDim(0);
        const index_t bottom = getFirstInDim(1);
        const dim_t NN0 = m_NN[0];
        const dim_t NN1 = m_NN[1];
#pragma omp parallel for
        for (index_t i=0; i < numMatrixRows; i++) {
            const index_t x = left + i % nDOF0;
            const index_t y = bottom + i / nDOF0;
            // loop through potential neighbours and add to index if positions
            // are within bounds
            for (dim_t i1=y-1; i1<y+2; i1++) {
                for (dim_t i0=x-1; i0<x+2; i0++) {
                    if (i0>=0 && i1>=0 && i0<NN0 && i1<NN1) {
                        indices[i].push_back(m_dofMap[i1*NN0 + i0]);
                    }
                }
            }
            sort(indices[i].begin(), indices[i].end());
        }
    } else {
#pragma omp parallel for
        for (index_t i=0; i < numMatrixRows; i++) {
            const index_t x = i % nDOF0;
            const index_t y = i / nDOF0;
            // loop through potential neighbours and add to index if positions
            // are within bounds
            for (dim_t i1=y-1; i1<y+2; i1++) {
                for (dim_t i0=x-1; i0<x+2; i0++) {
                    if (i0>=0 && i1>=0 && i0<nDOF0 && i1<nDOF1) {
                        indices[i].push_back(i1*nDOF0 + i0);
                    }
                }
            }
        }
    }
    return indices;
}

//private
void Rectangle::populateDofMap()
{
    const dim_t nDOF0 = getNumDOFInAxis(0);
    const dim_t nDOF1 = getNumDOFInAxis(1);
    const index_t left = getFirstInDim(0);
    const index_t bottom = getFirstInDim(1);

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
    RankVector neighbour;
    IndexVector offsetInShared(1,0);
    IndexVector sendShared, recvShared;
    const int x=m_mpiInfo->rank%m_NX[0];
    const int y=m_mpiInfo->rank/m_NX[0];
    // numShared will contain the number of shared DOFs after the following
    // blocks
    dim_t numShared=0;
    // sharing bottom edge
    if (y > 0) {
        neighbour.push_back((y-1)*m_NX[0] + x);
        const dim_t num = nDOF0;
        offsetInShared.push_back(offsetInShared.back()+num);
        for (dim_t i=0; i<num; i++, numShared++) {
            sendShared.push_back(i);
            recvShared.push_back(numDOF+numShared);
            m_dofMap[left+i]=numDOF+numShared;
        }
    }
    // sharing top edge
    if (y < m_NX[1] - 1) {
        neighbour.push_back((y+1)*m_NX[0] + x);
        const dim_t num = nDOF0;
        offsetInShared.push_back(offsetInShared.back()+num);
        for (dim_t i=0; i<num; i++, numShared++) {
            sendShared.push_back(numDOF-num+i);
            recvShared.push_back(numDOF+numShared);
            m_dofMap[m_NN[0]*(m_NN[1]-1)+left+i]=numDOF+numShared;
        }
    }
    // sharing left edge
    if (x > 0) {
        neighbour.push_back(y*m_NX[0] + x-1);
        const dim_t num = nDOF1;
        offsetInShared.push_back(offsetInShared.back()+num);
        for (dim_t i=0; i<num; i++, numShared++) {
            sendShared.push_back(i*nDOF0);
            recvShared.push_back(numDOF+numShared);
            m_dofMap[(bottom+i)*m_NN[0]]=numDOF+numShared;
        }
    }
    // sharing right edge
    if (x < m_NX[0] - 1) {
        neighbour.push_back(y*m_NX[0] + x+1);
        const dim_t num = nDOF1;
        offsetInShared.push_back(offsetInShared.back()+num);
        for (dim_t i=0; i<num; i++, numShared++) {
            sendShared.push_back((i+1)*nDOF0-1);
            recvShared.push_back(numDOF+numShared);
            m_dofMap[(bottom+1+i)*m_NN[0]-1]=numDOF+numShared;
        }
    }
    // sharing bottom-left node
    if (x > 0 && y > 0) {
        neighbour.push_back((y-1)*m_NX[0] + x-1);
        // sharing a node
        offsetInShared.push_back(offsetInShared.back()+1);
        sendShared.push_back(0);
        recvShared.push_back(numDOF+numShared);
        m_dofMap[0]=numDOF+numShared;
        ++numShared;
    }
    // sharing top-left node
    if (x > 0 && y < m_NX[1]-1) {
        neighbour.push_back((y+1)*m_NX[0] + x-1);
        offsetInShared.push_back(offsetInShared.back()+1);
        sendShared.push_back(numDOF-nDOF0);
        recvShared.push_back(numDOF+numShared);
        m_dofMap[m_NN[0]*(m_NN[1]-1)]=numDOF+numShared;
        ++numShared;
    }
    // sharing bottom-right node
    if (x < m_NX[0]-1 && y > 0) {
        neighbour.push_back((y-1)*m_NX[0] + x+1);
        offsetInShared.push_back(offsetInShared.back()+1);
        sendShared.push_back(nDOF0-1);
        recvShared.push_back(numDOF+numShared);
        m_dofMap[m_NN[0]-1]=numDOF+numShared;
        ++numShared;
    }
    // sharing top-right node
    if (x < m_NX[0]-1 && y < m_NX[1]-1) {
        neighbour.push_back((y+1)*m_NX[0] + x+1);
        offsetInShared.push_back(offsetInShared.back()+1);
        sendShared.push_back(numDOF-1);
        recvShared.push_back(numDOF+numShared);
        m_dofMap[m_NN[0]*m_NN[1]-1]=numDOF+numShared;
        ++numShared;
    }

#ifdef ESYS_HAVE_PASO
    createPasoConnector(neighbour, offsetInShared, offsetInShared, sendShared,
                        recvShared);
#endif

    // useful debug output
    /*
    std::cout << "--- rcv_shcomp ---" << std::endl;
    std::cout << "numDOF=" << numDOF << ", numNeighbors=" << neighbour.size() << std::endl;
    for (size_t i=0; i<neighbour.size(); i++) {
        std::cout << "neighbor[" << i << "]=" << neighbour[i]
            << " offsetInShared[" << i+1 << "]=" << offsetInShared[i+1] << std::endl;
    }
    for (size_t i=0; i<recvShared.size(); i++) {
        std::cout << "shared[" << i << "]=" << recvShared[i] << std::endl;
    }
    std::cout << "--- snd_shcomp ---" << std::endl;
    for (size_t i=0; i<sendShared.size(); i++) {
        std::cout << "shared[" << i << "]=" << sendShared[i] << std::endl;
    }
    std::cout << "--- dofMap ---" << std::endl;
    for (size_t i=0; i<m_dofMap.size(); i++) {
        std::cout << "m_dofMap[" << i << "]=" << m_dofMap[i] << std::endl;
    }
    */
}

//private
template<typename Scalar>
void Rectangle::addToMatrixAndRHS(AbstractSystemMatrix* S, escript::Data& F,
         const vector<Scalar>& EM_S, const vector<Scalar>& EM_F, bool addS,
         bool addF, index_t firstNode, int nEq, int nComp) const
{
    IndexVector rowIndex(4);
    rowIndex[0] = m_dofMap[firstNode];
    rowIndex[1] = m_dofMap[firstNode+1];
    rowIndex[2] = m_dofMap[firstNode+m_NN[0]];
    rowIndex[3] = m_dofMap[firstNode+m_NN[0]+1];
    if (addF) {
        Scalar* F_p = F.getSampleDataRW(0, static_cast<Scalar>(0));
        for (index_t i=0; i<rowIndex.size(); i++) {
            if (rowIndex[i]<getNumDOF()) {
                for (int eq=0; eq<nEq; eq++) {
                    F_p[INDEX2(eq, rowIndex[i], nEq)]+=EM_F[INDEX2(eq,i,nEq)];
                }
            }
        }
    }
    if (addS) {
        addToSystemMatrix<Scalar>(S, rowIndex, nEq, EM_S);
    }
}

template
void Rectangle::addToMatrixAndRHS<real_t>(AbstractSystemMatrix* S, escript::Data& F,
         const vector<real_t>& EM_S, const vector<real_t>& EM_F, bool addS,
         bool addF, index_t firstNode, int nEq, int nComp) const;

template
void Rectangle::addToMatrixAndRHS<cplx_t>(AbstractSystemMatrix* S, escript::Data& F,
         const vector<cplx_t>& EM_S, const vector<cplx_t>& EM_F, bool addS,
         bool addF, index_t firstNode, int nEq, int nComp) const;

//protected
void Rectangle::interpolateNodesOnElements(escript::Data& out,
                                           const escript::Data& in,
                                           bool reduced) const
{
    if (in.isComplex()!=out.isComplex())
    {
        throw RipleyException("Programmer Error: in and out parameters do not have the same complexity.");
    }
    if (in.isComplex())
    {
        interpolateNodesOnElementsWorker(out, in, reduced, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateNodesOnElementsWorker(out, in, reduced, escript::DataTypes::real_t(0));
    }
}

//protected
void Rectangle::interpolateNodesOnFaces(escript::Data& out,
                                           const escript::Data& in,
                                           bool reduced) const
{
    if (in.isComplex()!=out.isComplex())
    {
        throw RipleyException("Programmer Error: in and out parameters do not have the same complexity.");
    }
    if (in.isComplex())
    {
        interpolateNodesOnFacesWorker(out, in, reduced, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateNodesOnFacesWorker(out, in, reduced, escript::DataTypes::real_t(0));
    }
}

// private
template <typename S>
void Rectangle::interpolateNodesOnElementsWorker(escript::Data& out,
                                           const escript::Data& in,
                                           bool reduced, S sentinel) const
{
    const dim_t numComp = in.getDataPointSize();
    const dim_t NE0 = m_NE[0];
    const dim_t NE1 = m_NE[1];
    if (reduced) {
        out.requireWrite();
        const S c0 = 0.25;
#pragma omp parallel
        {
            vector<S> f_00(numComp);
            vector<S> f_01(numComp);
            vector<S> f_10(numComp);
            vector<S> f_11(numComp);
#pragma omp for
            for (index_t k1=0; k1 < NE1; ++k1) {
                for (index_t k0=0; k0 < NE0; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,k1, m_NN[0]), sentinel), numComp*sizeof(S));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,k1, m_NN[0]), sentinel), numComp*sizeof(S));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));
                    S* o = out.getSampleDataRW(INDEX2(k0,k1,m_NE[0]), sentinel);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = c0*(f_00[i] + f_01[i] + f_10[i] + f_11[i]);
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of k1 loop */
        } /* end of parallel section */
    } else {
        out.requireWrite();
        const S c0 = 0.16666666666666666667;
        const S c1 = 0.044658198738520451079;
        const S c2 = 0.62200846792814621559;
#pragma omp parallel
        {
            vector<S> f_00(numComp);
            vector<S> f_01(numComp);
            vector<S> f_10(numComp);
            vector<S> f_11(numComp);
#pragma omp for
            for (index_t k1=0; k1 < NE1; ++k1) {
                for (index_t k0=0; k0 < NE0; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,k1, m_NN[0]), sentinel), numComp*sizeof(S));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,k1, m_NN[0]), sentinel), numComp*sizeof(S));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));
                    S* o = out.getSampleDataRW(INDEX2(k0,k1,m_NE[0]), sentinel);
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

//private
template <typename S>
void Rectangle::interpolateNodesOnFacesWorker(escript::Data& out,
                                        const escript::Data& in,
                                        bool reduced, S sentinel) const
{
    const dim_t numComp = in.getDataPointSize();
    const dim_t NE0 = m_NE[0];
    const dim_t NE1 = m_NE[1];
    if (reduced) {
        out.requireWrite();
#pragma omp parallel
        {
            vector<S> f_00(numComp);
            vector<S> f_01(numComp);
            vector<S> f_10(numComp);
            vector<S> f_11(numComp);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < NE1; ++k1) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(0,k1, m_NN[0]), sentinel), numComp*sizeof(S));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(0,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));
                    S* o = out.getSampleDataRW(m_faceOffset[0]+k1, sentinel);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = (f_00[i] + f_01[i])/static_cast<S>(2);
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of face 0 */
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < NE1; ++k1) {
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1, m_NN[0]), sentinel), numComp*sizeof(S));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));
                    S* o = out.getSampleDataRW(m_faceOffset[1]+k1, sentinel);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = (f_10[i] + f_11[i])/static_cast<S>(2);
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of face 1 */
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < NE0; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,0, m_NN[0]), sentinel), numComp*sizeof(S));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,0, m_NN[0]), sentinel), numComp*sizeof(S));
                    S* o = out.getSampleDataRW(m_faceOffset[2]+k0, sentinel);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = (f_00[i] + f_10[i])/static_cast<S>(2);
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of face 2 */
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < NE0; ++k0) {
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,m_NN[1]-1, m_NN[0]), sentinel), numComp*sizeof(S));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,m_NN[1]-1, m_NN[0]), sentinel), numComp*sizeof(S));
                    S* o = out.getSampleDataRW(m_faceOffset[3]+k0, sentinel);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = (f_01[i] + f_11[i])/static_cast<S>(2);
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of face 3 */
        } /* end of parallel section */
    } else {
        out.requireWrite();
        const S c0 = 0.21132486540518711775;
        const S c1 = 0.78867513459481288225;
#pragma omp parallel
        {
            vector<S> f_00(numComp);
            vector<S> f_01(numComp);
            vector<S> f_10(numComp);
            vector<S> f_11(numComp);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < NE1; ++k1) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(0,k1, m_NN[0]), sentinel), numComp*sizeof(S));
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(0,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));
                    S* o = out.getSampleDataRW(m_faceOffset[0]+k1, sentinel);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = c0*f_01[i] + c1*f_00[i];
                        o[INDEX2(i,numComp,1)] = c0*f_00[i] + c1*f_01[i];
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of face 0 */
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < NE1; ++k1) {
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1, m_NN[0]), sentinel), numComp*sizeof(S));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));
                    S* o = out.getSampleDataRW(m_faceOffset[1]+k1, sentinel);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = c1*f_10[i] + c0*f_11[i];
                        o[INDEX2(i,numComp,1)] = c1*f_11[i] + c0*f_10[i];
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of face 1 */
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < NE0; ++k0) {
                    memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,0, m_NN[0]), sentinel), numComp*sizeof(S));
                    memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,0, m_NN[0]), sentinel), numComp*sizeof(S));
                    S* o = out.getSampleDataRW(m_faceOffset[2]+k0, sentinel);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = c0*f_10[i] + c1*f_00[i];
                        o[INDEX2(i,numComp,1)] = c0*f_00[i] + c1*f_10[i];
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of face 2 */
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < NE0; ++k0) {
                    memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,m_NN[1]-1, m_NN[0]), sentinel), numComp*sizeof(S));
                    memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,m_NN[1]-1, m_NN[0]), sentinel), numComp*sizeof(S));
                    S* o = out.getSampleDataRW(m_faceOffset[3]+k0, sentinel);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = c0*f_11[i] + c1*f_01[i];
                        o[INDEX2(i,numComp,1)] = c0*f_01[i] + c1*f_11[i];
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of face 3 */
        } /* end of parallel section */
    }
}



namespace
{
    // Calculates a Gaussian blur convolution matrix for 2D
    // See wiki article on the subject
    double* get2DGauss(unsigned radius, double sigma)
    {
        double* arr = new double[(radius*2+1)*(radius*2+1)];
        const double common = M_1_PI * 0.5 / (sigma*sigma);
        const int r = static_cast<int>(radius);
        double total = 0;
        for (int y = -r; y <= r; ++y) {
            for (int x = -r; x <= r; ++x) {
                arr[(x+r)+(y+r)*(r*2+1)]=common*exp(-(x*x+y*y)/(2*sigma*sigma));
                total+=arr[(x+r)+(y+r)*(r*2+1)];
            }
        }
        const double invtotal = 1/total;
        for (size_t p=0; p<(radius*2+1)*(radius*2+1); ++p) {
            arr[p] *= invtotal;
        }
        return arr;
    }

    // applies conv to source to get a point.
    // (xp, yp) are the coords in the source matrix not the destination matrix
    double Convolve2D(double* conv, double* source, size_t xp, size_t yp,
                      unsigned radius, size_t width)
    {
        const size_t bx = xp-radius, by=yp-radius;
        const size_t sbase = bx+by*width;
        double result = 0;
        for (int y=0; y<2*radius+1; ++y) {
            for (int x=0; x<2*radius+1; ++x) {
                result += conv[x+y*(2*radius+1)] * source[sbase + x+y*width];
            }
        }
        return result;
    }
}

/* This is a wrapper for filtered (and non-filtered) randoms
 * For detailed doco see randomFillWorker
 */
escript::Data Rectangle::randomFill(const escript::DataTypes::ShapeType& shape,
                                const escript::FunctionSpace& what, long seed,
                                const bp::tuple& filter) const
{
    int numvals=escript::DataTypes::noValues(shape);
    if (len(filter) > 0 && numvals != 1)
        throw NotImplementedError("Ripley only supports filters for scalar data.");

    escript::Data res = randomFillWorker(shape, seed, filter);
    if (res.getFunctionSpace() != what) {
        escript::Data r(res, what);
        return r;
    }
    return res;
}


/* This routine produces a Data object filled with smoothed random data.
 * The dimensions of the rectangle being filled are internal[0] x internal[1]
 * points. A parameter radius gives the size of the stencil used for the
 * smoothing.  A point on the left hand edge for example, will still require
 * `radius` extra points to the left in order to complete the stencil.
 *
 * All local calculation is done on an array called `src`, with
 * dimensions = ext[0] * ext[1], where ext[i]= internal[i]+2*radius.
 *
 * Now for MPI there is overlap to deal with. We need to share both the
 * overlapping values themselves but also the external region.
 *
 * In a hypothetical 1-D case:
 *
 * 1234567 would be split into two ranks thus:
 * 123(4)  (4)567     [4 being a shared element]
 *
 * If the radius is 2. There will be padding elements on the outside:
 * pp123(4)  (4)567pp
 *
 * To ensure that 4 can be correctly computed on both ranks, values from the
 * other rank need to be known.
 *
 * pp123(4)56   23(4)567pp
 *
 * Now in our case, we set all the values 23456 on the left rank and send them
 * to the right hand rank.
 *
 * So the edges _may_ need to be shared at a distance `inset` from all
 * boundaries.
 *
 * inset=2*radius+1
 * This is to ensure that values at distance `radius` from the
 * shared/overlapped element that ripley has.
 */
escript::Data Rectangle::randomFillWorker(
                        const escript::DataTypes::ShapeType& shape, long seed,
                        const bp::tuple& filter) const
{
    unsigned int radius=0;  // these are only used by gaussian
    double sigma=0.5;

    unsigned int numvals=escript::DataTypes::noValues(shape);

    if (len(filter) == 0) {
        // nothing special required here yet
    } else if (len(filter) == 3) {
        bp::extract<string> ex(filter[0]);
        if (!ex.check() || (ex()!="gaussian")) {
            throw ValueError("Unsupported random filter");
        }
        bp::extract<unsigned int> ex1(filter[1]);
        if (!ex1.check()) {
            throw ValueError("Radius of Gaussian filter must be a positive integer.");
        }
        radius = ex1();
        sigma = 0.5;
        bp::extract<double> ex2(filter[2]);
        if (!ex2.check() || (sigma=ex2()) <= 0) {
            throw ValueError("Sigma must be a positive floating point number.");
        }
    } else {
        throw ValueError("Unsupported random filter for Rectangle.");
    }

    // number of points in the internal region
    // that is, the ones we need smoothed versions of
    const dim_t internal[2] = { m_NN[0], m_NN[1] };
    size_t ext[2];
    ext[0]=(size_t)internal[0]+2*radius; // includes points we need as input
    ext[1]=(size_t)internal[1]+2*radius; // for smoothing

    // now we check to see if the radius is acceptable
    // That is, would not cross multiple ranks in MPI

    if (2*radius >= internal[0]-4) {
        throw ValueError("Radius of gaussian filter is too large for X dimension of a rank");
    }
    if (2*radius >= internal[1]-4) {
        throw ValueError("Radius of gaussian filter is too large for Y dimension of a rank");
    }

    double* src = new double[ext[0]*ext[1]*numvals];
    escript::randomFillArray(seed, src, ext[0]*ext[1]*numvals, m_mpiInfo);

#ifdef ESYS_MPI
    if ((internal[0] < 5) || (internal[1] < 5)) {
        // since the dimensions are equal for all ranks, this exception
        // will be thrown on all ranks
        throw RipleyException("Random Data in Ripley requires at least five elements per side per rank.");
    }
    dim_t X = m_mpiInfo->rank%m_NX[0];
    dim_t Y = m_mpiInfo->rank/m_NX[0];
#endif

/*
    // if we wanted to test a repeating pattern
    size_t basex=0;
    size_t basey=0;
#ifdef ESYS_MPI
    basex=X*m_gNE[0]/m_NX[0];
    basey=Y*m_gNE[1]/m_NX[1];
#endif

    escript::patternFillArray2D(ext[0], ext[1], src, 4, basex, basey, numvals);
*/

#ifdef ESYS_MPI
    BlockGrid2 grid(m_NX[0]-1, m_NX[1]-1);
    // it's +2 not +1 because a whole element is shared (and hence there is
    // an overlap of two points both of which need to have "radius" points on
    // either side.
    size_t inset=2*radius+2;

    // how wide is the x-dimension between the two insets
    size_t xmidlen=ext[0]-2*inset;
    size_t ymidlen=ext[1]-2*inset;

    Block2 block(ext[0], ext[1], inset, xmidlen, ymidlen, numvals);

    // a non-tight upper bound on how many we need
    MPI_Request reqs[40];
    MPI_Status stats[40];
    short rused=0;

    messvec incoms;
    messvec outcoms;

    grid.generateInNeighbours(X, Y, incoms);
    grid.generateOutNeighbours(X, Y, outcoms);

    block.copyAllToBuffer(src);

    int comserr = 0;
    for (size_t i=0; i < incoms.size(); ++i) {
        message& m = incoms[i];
        comserr |= MPI_Irecv(block.getInBuffer(m.destbuffid),
                             block.getBuffSize(m.destbuffid), MPI_DOUBLE,
                             m.sourceID, m.tag, m_mpiInfo->comm,
                             reqs+(rused++));
        block.setUsed(m.destbuffid);
    }

    for (size_t i=0; i < outcoms.size(); ++i) {
        message& m = outcoms[i];
        comserr |= MPI_Isend(block.getOutBuffer(m.srcbuffid),
                             block.getBuffSize(m.srcbuffid), MPI_DOUBLE,
                             m.destID, m.tag, m_mpiInfo->comm, reqs+(rused++));
    }

    if (!comserr) {
        comserr = MPI_Waitall(rused, reqs, stats);
    }

    if (comserr) {
        // Yes this is throwing an exception as a result of an MPI error
        // and no we don't inform the other ranks that we are doing this.
        // However, we have no reason to believe coms work at this point anyway
        throw RipleyException("Error in coms for randomFill");
    }

    block.copyUsedFromBuffer(src);
#endif

    // the truth of either should imply the truth of the other but let's be safe
    if (radius==0 || numvals > 1) {
        escript::FunctionSpace fs(getPtr(), getContinuousFunctionCode());
        escript::Data resdat(0, shape, fs, true);
        // don't need to check for exwrite because we just made it
        escript::DataTypes::RealVectorType& dv = resdat.getExpandedVectorReference();

        // now we need to copy values over
        for (size_t y=0; y < internal[1]; ++y) {
            for (size_t x=0; x < internal[0]; ++x) {
                for (unsigned int i=0; i < numvals; ++i) {
                    dv[i+(x+y*(internal[0]))*numvals]=src[i+(x+y*ext[0])*numvals];
                }
            }
        }
        delete[] src;
        return resdat;
    } else { // filter enabled
        escript::FunctionSpace fs(getPtr(), getContinuousFunctionCode());
        escript::Data resdat(0, escript::DataTypes::scalarShape, fs, true);
        // don't need to check for exwrite because we just made it
        escript::DataTypes::RealVectorType& dv=resdat.getExpandedVectorReference();
        double* convolution=get2DGauss(radius, sigma);
        for (size_t y=0; y < internal[1]; ++y) {
            for (size_t x=0; x < internal[0]; ++x) {
                dv[x+y*(internal[0])] = Convolve2D(convolution, src, x+radius, y+radius, radius, ext[0]);
            }
        }
        delete[] convolution;
        delete[] src;
        return resdat;
    }
}

dim_t Rectangle::findNode(const double *coords) const
{
    const dim_t NOT_MINE = -1;
    //is the found element even owned by this rank
    // (inside owned or shared elements but will map to an owned element)
    for (int dim = 0; dim < m_numDim; dim++) {
        //allows for point outside mapping onto node
        double min = m_origin[dim] + m_offset[dim]* m_dx[dim]
                - m_dx[dim]/2. + escript::DataTypes::real_t_eps();
        double max = m_origin[dim] + (m_offset[dim] + m_NE[dim])*m_dx[dim]
                + m_dx[dim]/2. - escript::DataTypes::real_t_eps();
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
    for (int dx = 0; dx < 1; dx++) {
        double xdist = (x - (ex + dx)*m_dx[0]);
        for (int dy = 0; dy < 1; dy++) {
            double ydist = (y - (ey + dy)*m_dx[1]);
            double total = xdist*xdist + ydist*ydist;
            if (total < minDist) {
                closest = INDEX2(ex+dx-m_offset[0], ey+dy-m_offset[1], m_NN[0]);
                minDist = total;
            }
        }
    }
    //if this happens, we've let a dirac point slip through, which is awful
    if (closest == NOT_MINE) {
        throw RipleyException("Unable to map appropriate dirac point to a node,"
                " implementation problem in Rectangle::findNode()");
    }
    return closest;
}

Assembler_ptr Rectangle::createAssembler(string type,
                                         const DataMap& constants) const
{
    bool isComplex = false;
    DataMap::const_iterator it;
    for (it = constants.begin(); it != constants.end(); it++) {
        if (!it->second.isEmpty() && it->second.isComplex()) {
            isComplex = true;
            break;
        }
    }

    if (type.compare("DefaultAssembler") == 0) {
        if (isComplex) {
            return Assembler_ptr(new DefaultAssembler2D<cplx_t>(shared_from_this(), m_dx, m_NE, m_NN));
        } else {
            return Assembler_ptr(new DefaultAssembler2D<real_t>(shared_from_this(), m_dx, m_NE, m_NN));
        }
    } else if (type.compare("WaveAssembler") == 0) {
        return Assembler_ptr(new WaveAssembler2D(shared_from_this(), m_dx, m_NE, m_NN, constants));
    } else if (type.compare("LameAssembler") == 0) {
        return Assembler_ptr(new LameAssembler2D(shared_from_this(), m_dx, m_NE, m_NN));
    }
    throw NotImplementedError("Ripley::Rectangle does not support the"
                              " requested assembler");
}

} // end of namespace ripley
