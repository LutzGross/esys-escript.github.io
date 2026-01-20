
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <ripley/Brick.h>
#include <ripley/DefaultAssembler3D.h>
#include <ripley/LameAssembler3D.h>
#include <ripley/WaveAssembler3D.h>
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
using escript::NotImplementedError;
using escript::ValueError;
using std::vector;
using std::string;
using std::min;
using std::max;
using std::ios;
using std::fill;

#ifdef ESYS_HAVE_NETCDF4
#include <ncFile.h>
#include <ncVar.h>
#include <ncDim.h>
using namespace netCDF;
#endif

namespace ripley {

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

Brick::Brick(dim_t n0, dim_t n1, dim_t n2, double x0, double y0, double z0,
             double x1, double y1, double z1, int d0, int d1, int d2,
             const vector<double>& points, const vector<int>& tags,
             const TagMap& tagnamestonums) :
    RipleyDomain(3)
{
    if (static_cast<long>(n0 + 1) * static_cast<long>(n1 + 1)
            * static_cast<long>(n2 + 1) > std::numeric_limits<dim_t>::max())
        throw RipleyException("The number of elements has overflowed, this "
                "limit may be raised in future releases.");

    if (n0 <= 0 || n1 <= 0 || n2 <= 0)
        throw ValueError("Number of elements in each spatial dimension "
                "must be positive");

    // ignore subdivision parameters for serial run
    if (m_mpiInfo->size == 1) {
        d0=1;
        d1=1;
        d2=1;
    }
    bool warn=false;

    vector<int> factors;
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
                throw ValueError("Invalid number of spatial subdivisions");
            }
            //remove
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
    if (d0*d1*d2 != m_mpiInfo->size){
        throw ValueError("Invalid number of spatial subdivisions");
    }
    if (warn) {
        std::cout << "ripley.Brick: Warning: Automatic domain subdivision (d0=" << d0 << ", d1="
            << d1 << ", d2=" << d2 << "). This may not be optimal!" << std::endl;
    }

    double l0 = x1-x0;
    double l1 = y1-y0;
    double l2 = z1-z0;
    m_dx[0] = l0/n0;
    m_dx[1] = l1/n1;
    m_dx[2] = l2/n2;

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
    if ((n2+1)%d2 > 0) {
        switch (getDecompositionPolicy()) {
            case DECOMP_EXPAND:
                l2 = m_dx[2]*n2; // fall through
            case DECOMP_ADD_ELEMENTS:
                n2 = (dim_t)round((float)(n2+1)/d2+0.5)*d2-1; // fall through
            case DECOMP_STRICT:
                warn = true;
                break;
        }
        // reset spacing
        m_dx[2] = l2/n2;
    }

    if ((d0 > 1 && (n0+1)/d0<2) || (d1 > 1 && (n1+1)/d1<2) || (d2 > 1 && (n2+1)/d2<2))
        throw ValueError("Too few elements for the number of ranks");

    if (warn) {
        if (getDecompositionPolicy() == DECOMP_STRICT) {
            throw ValueError("Unable to decompose domain to the number of "
                    "MPI ranks without adding elements and the policy "
                    "is set to STRICT. Use setDecompositionPolicy() "
                    "to allow adding elements.");
        } else {
            std::cout << "ripley.Brick: Warning: Domain setup has been adjusted as follows "
                    "to allow decomposition into " << m_mpiInfo->size
                    << " MPI ranks:" << std::endl
                    << "    N0=" << n0 << ", l0=" << l0 << std::endl
                    << "    N1=" << n1 << ", l1=" << l1 << std::endl
                    << "    N2=" << n2 << ", l2=" << l1 << std::endl;
        }
    }
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

    for (TagMap::const_iterator i = tagnamestonums.begin();
            i != tagnamestonums.end(); i++) {
        setTagMap(i->first, i->second);
    }
    addPoints(points, tags);
}

Brick::Brick(escript::JMPI jmpi, dim_t n0, dim_t n1, dim_t n2, double x0, double y0, double z0,
             double x1, double y1, double z1, int d0, int d1, int d2,
             const vector<double>& points, const vector<int>& tags,
             const TagMap& tagnamestonums) :
    RipleyDomain(3, jmpi)
{
    if (static_cast<long>(n0 + 1) * static_cast<long>(n1 + 1)
            * static_cast<long>(n2 + 1) > std::numeric_limits<dim_t>::max())
        throw RipleyException("The number of elements has overflowed, this "
                "limit may be raised in future releases.");

    if (n0 <= 0 || n1 <= 0 || n2 <= 0)
        throw ValueError("Number of elements in each spatial dimension "
                "must be positive");

    // ignore subdivision parameters for serial run
    if (m_mpiInfo->size == 1) {
        d0=1;
        d1=1;
        d2=1;
    }
    bool warn=false;

    vector<int> factors;
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
                throw ValueError("Invalid number of spatial subdivisions");
            }
            //remove
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
    if (d0*d1*d2 != m_mpiInfo->size){
        throw ValueError("Invalid number of spatial subdivisions");
    }
    if (warn) {
        std::cout << "ripley.Brick: Warning: Automatic domain subdivision (d0=" << d0 << ", d1="
            << d1 << ", d2=" << d2 << "). This may not be optimal!" << std::endl;
    }

    double l0 = x1-x0;
    double l1 = y1-y0;
    double l2 = z1-z0;
    m_dx[0] = l0/n0;
    m_dx[1] = l1/n1;
    m_dx[2] = l2/n2;

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
    if ((n2+1)%d2 > 0) {
        switch (getDecompositionPolicy()) {
            case DECOMP_EXPAND:
                l2 = m_dx[2]*n2; // fall through
            case DECOMP_ADD_ELEMENTS:
                n2 = (dim_t)round((float)(n2+1)/d2+0.5)*d2-1; // fall through
            case DECOMP_STRICT:
                warn = true;
                break;
        }
        // reset spacing
        m_dx[2] = l2/n2;
    }

    if ((d0 > 1 && (n0+1)/d0<2) || (d1 > 1 && (n1+1)/d1<2) || (d2 > 1 && (n2+1)/d2<2))
        throw ValueError("Too few elements for the number of ranks");

    if (warn) {
        if (getDecompositionPolicy() == DECOMP_STRICT) {
            throw ValueError("Unable to decompose domain to the number of "
                    "MPI ranks without adding elements and the policy "
                    "is set to STRICT. Use setDecompositionPolicy() "
                    "to allow adding elements.");
        } else {
            std::cout << "ripley.Brick: Warning: Domain setup has been adjusted as follows "
                    "to allow decomposition into " << m_mpiInfo->size
                    << " MPI ranks:" << std::endl
                    << "    N0=" << n0 << ", l0=" << l0 << std::endl
                    << "    N1=" << n1 << ", l1=" << l1 << std::endl
                    << "    N2=" << n2 << ", l2=" << l1 << std::endl;
        }
    }
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

    for (TagMap::const_iterator i = tagnamestonums.begin();
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
    return "ripley::Brick";
}

bool Brick::operator==(const escript::AbstractDomain& other) const
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
#ifdef ESYS_HAVE_NETCDF4
    // check destination function space
    dim_t myN0, myN1, myN2;
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
        throw ValueError("readNcGrid(): invalid function space for output data object");

    if (params.first.size() != 3)
        throw ValueError("readNcGrid(): argument 'first' must have 3 entries");

    if (params.numValues.size() != 3)
        throw ValueError("readNcGrid(): argument 'numValues' must have 3 entries");

    if (params.multiplier.size() != 3)
        throw ValueError("readNcGrid(): argument 'multiplier' must have 3 entries");
    for (size_t i=0; i<params.multiplier.size(); i++)
        if (params.multiplier[i]<1)
            throw ValueError("readNcGrid(): all multipliers must be positive");

    // check file existence and size

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
        throw RipleyException("readNcGrid(): only scalar data supported");

    const int dims = var.getDimCount();
    vector<long> edges(dims);
    std::vector< NcDim > vard=var.getDims();
    for (size_t i=0;i<vard.size();++i)
    {
        edges[i]=vard[i].getSize();
    }

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
    vector<size_t> startindex;
    vector<size_t> counts;
    if (dims==3) {
        //var->set_cur(idx2, idx1, idx0);             // from old API
        startindex.push_back(idx2);
        startindex.push_back(idx1);
        startindex.push_back(idx0);
        counts.push_back(num2);
        counts.push_back(num1);
        counts.push_back(num0);
        var.getVar(startindex, counts, &values[0]);
    } else if (dims==2) {
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
                if (!bm::isnan(values[srcIndex])) {
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
    throw RipleyException("readNcGrid(): not compiled with netCDF support");
#endif
}

void Brick::readBinaryGridFromZipped(escript::Data& out, string filename,
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
            throw ValueError("readBinaryGridZipped(): invalid or unsupported datatype");
    }
#else
    throw RipleyException("readBinaryGridZipped(): not compiled with zip support");
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
            throw ValueError("readBinaryGrid(): invalid or unsupported datatype");
    }
}

template<typename ValueType>
void Brick::readBinaryGridImpl(escript::Data& out, const string& filename,
                               const ReaderParameters& params) const
{
    // check destination function space
    dim_t myN0, myN1, myN2;
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
        throw ValueError("readBinaryGrid(): invalid function space for output data object");

    if (params.first.size() != 3)
        throw ValueError("readBinaryGrid(): argument 'first' must have 3 entries");

    if (params.numValues.size() != 3)
        throw ValueError("readBinaryGrid(): argument 'numValues' must have 3 entries");

    if (params.multiplier.size() != 3)
        throw ValueError("readBinaryGrid(): argument 'multiplier' must have 3 entries");
    for (size_t i=0; i<params.multiplier.size(); i++)
        if (params.multiplier[i]<1)
            throw ValueError("readBinaryGrid(): all multipliers must be positive");
    if (params.reverse[0] != 0 || params.reverse[1] != 0)
        throw RipleyException("readBinaryGrid(): reversing only supported in Z-direction currently");

    // check file existence and size
    std::ifstream f(filename.c_str(), std::ifstream::binary);
    if (f.fail()) {
        throw RipleyException("readBinaryGrid(): cannot open file " + filename);
    }
    f.seekg(0, std::ios::end);
    const int numComp = out.getDataPointSize();
    const dim_t filesize = f.tellg();
    const dim_t reqsize = params.numValues[0]*params.numValues[1]*params.numValues[2]*numComp*sizeof(ValueType);
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
                            const dim_t dataIndex = dataX + dataY*myN0 + dataZ*myN0*myN1;
                            double* dest = out.getSampleDataRW(dataIndex);
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
        }
    }

    f.close();
}

#ifdef ESYS_HAVE_BOOST_IO
template<typename ValueType>
void Brick::readBinaryGridZippedImpl(escript::Data& out, const string& filename,
                               const ReaderParameters& params) const
{
    // check destination function space
    dim_t myN0, myN1, myN2;
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
        throw ValueError("readBinaryGridFromZipped(): invalid function space for output data object");

    if (params.first.size() != 3)
        throw ValueError("readBinaryGridFromZipped(): argument 'first' must have 3 entries");

    if (params.numValues.size() != 3)
        throw ValueError("readBinaryGridFromZipped(): argument 'numValues' must have 3 entries");

    if (params.multiplier.size() != 3)
        throw ValueError("readBinaryGridFromZipped(): argument 'multiplier' must have 3 entries");
    for (size_t i=0; i<params.multiplier.size(); i++)
        if (params.multiplier[i]<1)
            throw ValueError("readBinaryGridFromZipped(): all multipliers must be positive");

    // check file existence and size
    std::ifstream f(filename.c_str(), std::ifstream::binary);
    if (f.fail()) {
        throw RipleyException("readBinaryGridFromZipped(): cannot open file " + filename);
    }
    f.seekg(0, std::ios::end);
    const int numComp = out.getDataPointSize();
    dim_t filesize = f.tellg();
    f.seekg(0, ios::beg);
    vector<char> compressed(filesize);
    f.read((char*)&compressed[0], filesize);
    f.close();
    vector<char> decompressed = unzip(compressed);
    filesize = decompressed.size();
    const dim_t reqsize = params.numValues[0]*params.numValues[1]*params.numValues[2]*numComp*sizeof(ValueType);
    if (filesize < reqsize) {
        throw RipleyException("readBinaryGridFromZipped(): not enough data in file");
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
    // indices to first value in file
    const dim_t idx0 = max(dim_t(0), m_offset[0]-params.first[0]);
    const dim_t idx1 = max(dim_t(0), m_offset[1]-params.first[1]);
    const dim_t idx2 = max(dim_t(0), m_offset[2]-params.first[2]);
    // number of values to read
    const dim_t num0 = min(params.numValues[0]-idx0, myN0-first0);
    const dim_t num1 = min(params.numValues[1]-idx1, myN1-first1);
    const dim_t num2 = min(params.numValues[2]-idx2, myN2-first2);

    out.requireWrite();
    vector<ValueType> values(num0*numComp);
    const int dpp = out.getNumDataPointsPerSample();

    for (dim_t z=0; z<num2; z++) {
        for (dim_t y=0; y<num1; y++) {
            const dim_t fileofs = numComp*(idx0+(idx1+y)*params.numValues[0]
                             +(idx2+z)*params.numValues[0]*params.numValues[1]);
            memcpy((char*)&values[0], (char*)&decompressed[fileofs*sizeof(ValueType)], num0*numComp*sizeof(ValueType));

            for (dim_t x=0; x<num0; x++) {
                const dim_t baseIndex = first0+x*params.multiplier[0]
                                     +(first1+y*params.multiplier[1])*myN0
                                     +(first2+z*params.multiplier[2])*myN0*myN1;
                for (dim_t m2=0; m2<params.multiplier[2]; m2++) {
                    for (dim_t m1=0; m1<params.multiplier[1]; m1++) {
                        for (dim_t m0=0; m0<params.multiplier[0]; m0++) {
                            const dim_t dataIndex = baseIndex+m0
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
            throw ValueError("writeBinaryGrid(): invalid or unsupported datatype");
    }
}

template<typename ValueType>
void Brick::writeBinaryGridImpl(const escript::Data& in,
                                const string& filename, int byteOrder) const
{
    // check function space and determine number of points
    dim_t myN0, myN1, myN2;
    dim_t totalN0, totalN1, totalN2;
    dim_t offset0, offset1, offset2;
    if (in.getFunctionSpace().getTypeCode() == Nodes) {
        myN0 = m_NN[0];
        myN1 = m_NN[1];
        myN2 = m_NN[2];
        totalN0 = m_gNE[0]+1;
        totalN1 = m_gNE[1]+1;
        totalN2 = m_gNE[2]+1;
        offset0 = m_offset[0];
        offset1 = m_offset[1];
        offset2 = m_offset[2];
    } else if (in.getFunctionSpace().getTypeCode() == DegreesOfFreedom ||
            in.getFunctionSpace().getTypeCode() == ReducedDegreesOfFreedom) {
        myN0 = (m_gNE[0]+1)/m_NX[0];
        myN1 = (m_gNE[1]+1)/m_NX[1];
        myN2 = (m_gNE[2]+1)/m_NX[2];
        totalN0 = m_gNE[0]+1;
        totalN1 = m_gNE[1]+1;
        totalN2 = m_gNE[2]+1;
        offset0 = (m_offset[0]>0 ? m_offset[0]+1 : 0);
        offset1 = (m_offset[1]>0 ? m_offset[1]+1 : 0);
        offset2 = (m_offset[2]>0 ? m_offset[2]+1 : 0);
    } else if (in.getFunctionSpace().getTypeCode() == Elements ||
                in.getFunctionSpace().getTypeCode() == ReducedElements) {
        myN0 = m_NE[0];
        myN1 = m_NE[1];
        myN2 = m_NE[2];
        totalN0 = m_gNE[0];
        totalN1 = m_gNE[1];
        totalN2 = m_gNE[2];
        offset0 = m_offset[0];
        offset1 = m_offset[1];
        offset2 = m_offset[2];
    } else
        throw RipleyException("writeBinaryGrid(): unsupported function space");

    const int numComp = in.getDataPointSize();
    const int dpp = in.getNumDataPointsPerSample();
    const dim_t fileSize = sizeof(ValueType)*numComp*dpp*totalN0*totalN1*totalN2;

    if (numComp > 1 || dpp > 1)
        throw RipleyException("writeBinaryGrid(): only scalar, single-value data supported");

    // from here on we know that each sample consists of one value
    FileWriter fw(m_mpiInfo->comm);
    fw.openFile(filename, fileSize, true);
    MPIBarrier();

    for (index_t z=0; z<myN2; z++) {
        for (index_t y=0; y<myN1; y++) {
            const dim_t fileofs = (offset0+(offset1+y)*totalN0
                                +(offset2+z)*totalN0*totalN1)*sizeof(ValueType);
            std::ostringstream oss;

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

void Brick::write(const std::string& filename) const
{
    throw RipleyException("write: not supported");
}

void Brick::dump(const string& fileName) const
{
#ifdef ESYS_HAVE_SILO
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
    // converting to int!
    vector<int> dims(m_NN, m_NN+3);

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

const dim_t* Brick::borrowSampleReferenceIDs(int fsType) const
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
        case Points:
            return &m_diracPointNodeIDs[0];
        default:
            break;
    }

    std::stringstream msg;
    msg << "borrowSampleReferenceIDs: invalid function space type "<<fsType;
    throw ValueError(msg.str());
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

    std::stringstream msg;
    msg << "ownSample: invalid function space type " << fsType;
    throw ValueError(msg.str());
}

RankVector Brick::getOwnerVector(int fsType) const
{
    RankVector owner;
    const int rank = m_mpiInfo->rank;

    if (fsType == Elements || fsType == ReducedElements) {
        owner.assign(getNumElements(), rank);
        // shared plane in Y-Z
        if (m_faceCount[0]==0) {
            for (dim_t k2=0; k2<m_NE[2]; k2++) {
                for (dim_t k1=0; k1<m_NE[1]; k1++) {
                    const dim_t e=k2*m_NE[0]*m_NE[1]+k1*m_NE[0];
                    owner[e] = rank-1;
                }
            }
        }
        // shared plane in X-Z
        if (m_faceCount[2]==0) {
            for (dim_t k2=0; k2<m_NE[2]; k2++) {
                for (dim_t k0=0; k0<m_NE[0]; k0++) {
                    const dim_t e=k2*m_NE[0]*m_NE[1]+k0;
                    owner[e] = rank-m_NX[0];
                }
            }
        }
        // shared plane in X-Y
        if (m_faceCount[4]==0) {
            for (dim_t k1=0; k1<m_NE[1]; k1++) {
                for (dim_t k0=0; k0<m_NE[0]; k0++) {
                    const dim_t e=k1*m_NE[0]+k0;
                    owner[e] = rank-m_NX[0]*m_NX[1];
                }
            }
        }

    } else if (fsType == FaceElements || fsType == ReducedFaceElements) {
        owner.assign(getNumFaceElements(), rank);
        index_t offset=0;
        if (m_faceCount[0] > 0) {
            if (m_faceCount[2]==0) {
                for (dim_t k2=0; k2<m_NE[2]; k2++)
                    owner[k2*m_NE[1]] = rank-m_NX[0];
            }
            if (m_faceCount[4]==0) {
                for (dim_t k1=0; k1<m_NE[1]; k1++)
                    owner[k1] = rank-m_NX[0]*m_NX[1];
            }
            offset += m_faceCount[0];
        }
        if (m_faceCount[1] > 0) {
            if (m_faceCount[2]==0) {
                for (dim_t k2=0; k2<m_NE[2]; k2++)
                    owner[offset+k2*m_NE[1]] = rank-m_NX[0];
            }
            if (m_faceCount[4]==0) {
                for (dim_t k1=0; k1<m_NE[1]; k1++)
                    owner[offset+k1] = rank-m_NX[0]*m_NX[1];
            }
            offset += m_faceCount[1];
        }
        if (m_faceCount[2] > 0) {
            if (m_faceCount[0]==0) {
                for (dim_t k2=0; k2<m_NE[2]; k2++)
                    owner[offset+k2*m_NE[0]] = rank-1;
            }
            if (m_faceCount[4]==0) {
                for (dim_t k0=0; k0<m_NE[0]; k0++)
                    owner[offset+k0] = rank-m_NX[0]*m_NX[1];
            }
            offset += m_faceCount[2];
        }
        if (m_faceCount[3] > 0) {
            if (m_faceCount[0]==0) {
                for (dim_t k2=0; k2<m_NE[2]; k2++)
                    owner[offset+k2*m_NE[0]] = rank-1;
            }
            if (m_faceCount[4]==0) {
                for (dim_t k0=0; k0<m_NE[0]; k0++)
                    owner[offset+k0] = rank-m_NX[0]*m_NX[1];
            }
            offset += m_faceCount[3];
        }
        if (m_faceCount[4] > 0) {
            if (m_faceCount[0]==0) {
                for (dim_t k1=0; k1<m_NE[1]; k1++)
                    owner[offset+k1*m_NE[0]] = rank-1;
            }
            if (m_faceCount[2]==0) {
                for (dim_t k0=0; k0<m_NE[0]; k0++)
                    owner[offset+k0] = rank-m_NX[0];
            }
            offset += m_faceCount[4];
        }
        if (m_faceCount[5] > 0) {
            if (m_faceCount[0]==0) {
                for (dim_t k1=0; k1<m_NE[1]; k1++)
                    owner[offset+k1*m_NE[0]] = rank-1;
            }
            if (m_faceCount[2]==0) {
                for (dim_t k0=0; k0<m_NE[0]; k0++)
                    owner[offset+k0] = rank-m_NX[0];
            }
        }
    } else {
        throw ValueError("getOwnerVector: only valid for element types");
    }

    return owner;
}

void Brick::setToNormal(escript::Data& out) const
{
    const dim_t NE0 = m_NE[0];
    const dim_t NE1 = m_NE[1];
    const dim_t NE2 = m_NE[2];

    if (out.getFunctionSpace().getTypeCode() == FaceElements) {
        out.requireWrite();
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2 = 0; k2 < NE2; ++k2) {
                    for (index_t k1 = 0; k1 < NE1; ++k1) {
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
                for (index_t k2 = 0; k2 < NE2; ++k2) {
                    for (index_t k1 = 0; k1 < NE1; ++k1) {
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
                for (index_t k2 = 0; k2 < NE2; ++k2) {
                    for (index_t k0 = 0; k0 < NE0; ++k0) {
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
                for (index_t k2 = 0; k2 < NE2; ++k2) {
                    for (index_t k0 = 0; k0 < NE0; ++k0) {
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
                for (index_t k1 = 0; k1 < NE1; ++k1) {
                    for (index_t k0 = 0; k0 < NE0; ++k0) {
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
                for (index_t k1 = 0; k1 < NE1; ++k1) {
                    for (index_t k0 = 0; k0 < NE0; ++k0) {
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
                for (index_t k2 = 0; k2 < NE2; ++k2) {
                    for (index_t k1 = 0; k1 < NE1; ++k1) {
                        double* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]));
                        *o++ = -1.;
                        *o++ = 0.;
                        *o = 0.;
                    }
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2 = 0; k2 < NE2; ++k2) {
                    for (index_t k1 = 0; k1 < NE1; ++k1) {
                        double* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]));
                        *o++ = 1.;
                        *o++ = 0.;
                        *o = 0.;
                    }
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2 = 0; k2 < NE2; ++k2) {
                    for (index_t k0 = 0; k0 < NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]));
                        *o++ = 0.;
                        *o++ = -1.;
                        *o = 0.;
                    }
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2 = 0; k2 < NE2; ++k2) {
                    for (index_t k0 = 0; k0 < NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]));
                        *o++ = 0.;
                        *o++ = 1.;
                        *o = 0.;
                    }
                }
            }

            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < NE1; ++k1) {
                    for (index_t k0 = 0; k0 < NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]));
                        *o++ = 0.;
                        *o++ = 0.;
                        *o = -1.;
                    }
                }
            }

            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < NE1; ++k1) {
                    for (index_t k0 = 0; k0 < NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]));
                        *o++ = 0.;
                        *o++ = 0.;
                        *o = 1.;
                    }
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

void Brick::setToSize(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == Elements
            || out.getFunctionSpace().getTypeCode() == ReducedElements) {
        out.requireWrite();
        const dim_t numQuad=out.getNumDataPointsPerSample();
        const double size=sqrt(m_dx[0]*m_dx[0]+m_dx[1]*m_dx[1]+m_dx[2]*m_dx[2]);
        const dim_t NE = getNumElements();
#pragma omp parallel for
        for (index_t k = 0; k < NE; ++k) {
            double* o = out.getSampleDataRW(k);
            fill(o, o+numQuad, size);
        }
    } else if (out.getFunctionSpace().getTypeCode() == FaceElements
            || out.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        out.requireWrite();
        const dim_t numQuad=out.getNumDataPointsPerSample();
        const dim_t NE0 = m_NE[0];
        const dim_t NE1 = m_NE[1];
        const dim_t NE2 = m_NE[2];
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
        throw ValueError(msg.str());
    }
}

void Brick::Print_Mesh_Info(const bool full) const
{
    RipleyDomain::Print_Mesh_Info(full);
    if (full) {
        std::cout << "     Id  Coordinates" << std::endl;
        std::cout.precision(15);
        std::cout.setf(ios::scientific, std::ios::floatfield);
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
    int numDim = m_numDim;
    if (!arg.isDataPointShapeEqual(1, &numDim))
        throw ValueError("setToX: Invalid Data object shape");
    if (!arg.numSamplesEqual(1, getNumNodes()))
        throw ValueError("setToX: Illegal number of samples in Data object");

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
void Brick::assembleGradient(escript::Data& out,
                             const escript::Data& in) const
{
    if (out.isComplex() != in.isComplex())
        throw ValueError("Gradient: input & output complexity must match.");
    else if (in.isComplex())
        assembleGradientImpl<cplx_t>(out, in);
    else
        assembleGradientImpl<real_t>(out, in);
}

//protected
template<typename Scalar>
void Brick::assembleGradientImpl(escript::Data& out,
                                 const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    const double C0 = .044658198738520451079;
    const double C1 = .16666666666666666667;
    const double C2 = .21132486540518711775;
    const double C3 = .25;
    const double C4 = .5;
    const double C5 = .62200846792814621559;
    const double C6 = .78867513459481288225;
    const dim_t NE0 = m_NE[0];
    const dim_t NE1 = m_NE[1];
    const dim_t NE2 = m_NE[2];
    const Scalar zero = static_cast<Scalar>(0);

    if (out.getFunctionSpace().getTypeCode() == Elements) {
        out.requireWrite();
#pragma omp parallel
        {
            vector<Scalar> f_000(numComp, zero);
            vector<Scalar> f_001(numComp, zero);
            vector<Scalar> f_010(numComp, zero);
            vector<Scalar> f_011(numComp, zero);
            vector<Scalar> f_100(numComp, zero);
            vector<Scalar> f_101(numComp, zero);
            vector<Scalar> f_110(numComp, zero);
            vector<Scalar> f_111(numComp, zero);
#pragma omp for
            for (index_t k2=0; k2 < NE2; ++k2) {
                for (index_t k1=0; k1 < NE1; ++k1) {
                    for (index_t k0=0; k0 < NE0; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        Scalar* o = out.getSampleDataRW(INDEX3(k0,k1,k2,NE0,NE1), zero);
                        for (index_t i=0; i < numComp; ++i) {
                            const Scalar V0=((f_100[i]-f_000[i])*C5 + (f_111[i]-f_011[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / m_dx[0];
                            const Scalar V1=((f_110[i]-f_010[i])*C5 + (f_101[i]-f_001[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / m_dx[0];
                            const Scalar V2=((f_101[i]-f_001[i])*C5 + (f_110[i]-f_010[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / m_dx[0];
                            const Scalar V3=((f_111[i]-f_011[i])*C5 + (f_100[i]-f_000[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / m_dx[0];
                            const Scalar V4=((f_010[i]-f_000[i])*C5 + (f_111[i]-f_101[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / m_dx[1];
                            const Scalar V5=((f_110[i]-f_100[i])*C5 + (f_011[i]-f_001[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / m_dx[1];
                            const Scalar V6=((f_011[i]-f_001[i])*C5 + (f_110[i]-f_100[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / m_dx[1];
                            const Scalar V7=((f_111[i]-f_101[i])*C5 + (f_010[i]-f_000[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / m_dx[1];
                            const Scalar V8=((f_001[i]-f_000[i])*C5 + (f_111[i]-f_110[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / m_dx[2];
                            const Scalar V9=((f_101[i]-f_100[i])*C5 + (f_011[i]-f_010[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / m_dx[2];
                            const Scalar V10=((f_011[i]-f_010[i])*C5 + (f_101[i]-f_100[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / m_dx[2];
                            const Scalar V11=((f_111[i]-f_110[i])*C5 + (f_001[i]-f_000[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / m_dx[2];
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
            vector<Scalar> f_000(numComp, zero);
            vector<Scalar> f_001(numComp, zero);
            vector<Scalar> f_010(numComp, zero);
            vector<Scalar> f_011(numComp, zero);
            vector<Scalar> f_100(numComp, zero);
            vector<Scalar> f_101(numComp, zero);
            vector<Scalar> f_110(numComp, zero);
            vector<Scalar> f_111(numComp, zero);
#pragma omp for
            for (index_t k2=0; k2 < NE2; ++k2) {
                for (index_t k1=0; k1 < NE1; ++k1) {
                    for (index_t k0=0; k0 < NE0; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        Scalar* o = out.getSampleDataRW(INDEX3(k0,k1,k2,NE0,NE1), zero);
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
            vector<Scalar> f_000(numComp, zero);
            vector<Scalar> f_001(numComp, zero);
            vector<Scalar> f_010(numComp, zero);
            vector<Scalar> f_011(numComp, zero);
            vector<Scalar> f_100(numComp, zero);
            vector<Scalar> f_101(numComp, zero);
            vector<Scalar> f_110(numComp, zero);
            vector<Scalar> f_111(numComp, zero);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < NE2; ++k2) {
                    for (index_t k1=0; k1 < NE1; ++k1) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(0,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(0,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(0,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(0,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(1,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(1,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(1,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(1,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        Scalar* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,NE1), zero);
                        for (index_t i=0; i < numComp; ++i) {
                            const Scalar V0=((f_010[i]-f_000[i])*C6 + (f_011[i]-f_001[i])*C2) / m_dx[1];
                            const Scalar V1=((f_010[i]-f_000[i])*C2 + (f_011[i]-f_001[i])*C6) / m_dx[1];
                            const Scalar V2=((f_001[i]-f_000[i])*C6 + (f_010[i]-f_011[i])*C2) / m_dx[2];
                            const Scalar V3=((f_001[i]-f_000[i])*C2 + (f_011[i]-f_010[i])*C6) / m_dx[2];
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
                for (index_t k2=0; k2 < NE2; ++k2) {
                    for (index_t k1=0; k1 < NE1; ++k1) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        Scalar* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,NE1), zero);
                        for (index_t i=0; i < numComp; ++i) {
                            const Scalar V0=((f_110[i]-f_100[i])*C6 + (f_111[i]-f_101[i])*C2) / m_dx[1];
                            const Scalar V1=((f_110[i]-f_100[i])*C2 + (f_111[i]-f_101[i])*C6) / m_dx[1];
                            const Scalar V2=((f_101[i]-f_100[i])*C6 + (f_111[i]-f_110[i])*C2) / m_dx[2];
                            const Scalar V3=((f_101[i]-f_100[i])*C2 + (f_111[i]-f_110[i])*C6) / m_dx[2];
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
                for (index_t k2=0; k2 < NE2; ++k2) {
                    for (index_t k0=0; k0 < NE0; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,0,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,0,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,0,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,0,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        Scalar* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,NE0), zero);
                        for (index_t i=0; i < numComp; ++i) {
                            const Scalar V0=((f_100[i]-f_000[i])*C6 + (f_101[i]-f_001[i])*C2) / m_dx[0];
                            const Scalar V1=((f_001[i]-f_000[i])*C6 + (f_101[i]-f_100[i])*C2) / m_dx[2];
                            const Scalar V2=((f_001[i]-f_000[i])*C2 + (f_101[i]-f_100[i])*C6) / m_dx[2];
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
                for (index_t k2=0; k2 < NE2; ++k2) {
                    for (index_t k0=0; k0 < NE0; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-2,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-2,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-2,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-2,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        Scalar* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,NE0), zero);
                        for (index_t i=0; i < numComp; ++i) {
                            const Scalar V0=((f_110[i]-f_010[i])*C6 + (f_111[i]-f_011[i])*C2) / m_dx[0];
                            const Scalar V1=((f_110[i]-f_010[i])*C2 + (f_111[i]-f_011[i])*C6) / m_dx[0];
                            const Scalar V2=((f_011[i]-f_010[i])*C6 + (f_111[i]-f_110[i])*C2) / m_dx[2];
                            const Scalar V3=((f_011[i]-f_010[i])*C2 + (f_111[i]-f_110[i])*C6) / m_dx[2];
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
                for (index_t k1=0; k1 < NE1; ++k1) {
                    for (index_t k0=0; k0 < NE0; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,0, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,0, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,0, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,0, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        Scalar* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,NE0), zero);
                        for (index_t i=0; i < numComp; ++i) {
                            const Scalar V0=((f_100[i]-f_000[i])*C6 + (f_110[i]-f_010[i])*C2) / m_dx[0];
                            const Scalar V1=((f_100[i]-f_000[i])*C2 + (f_110[i]-f_010[i])*C6) / m_dx[0];
                            const Scalar V2=((f_010[i]-f_000[i])*C6 + (f_110[i]-f_100[i])*C2) / m_dx[1];
                            const Scalar V3=((f_010[i]-f_000[i])*C2 + (f_110[i]-f_100[i])*C6) / m_dx[1];
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
                for (index_t k1=0; k1 < NE1; ++k1) {
                    for (index_t k0=0; k0 < NE0; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        Scalar* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,NE0), zero);
                        for (index_t i=0; i < numComp; ++i) {
                            const Scalar V0=((f_101[i]-f_001[i])*C6 + (f_111[i]-f_011[i])*C2) / m_dx[0];
                            const Scalar V1=((f_101[i]-f_001[i])*C2 + (f_111[i]-f_011[i])*C6) / m_dx[0];
                            const Scalar V2=((f_011[i]-f_001[i])*C6 + (f_111[i]-f_101[i])*C2) / m_dx[1];
                            const Scalar V3=((f_011[i]-f_001[i])*C2 + (f_111[i]-f_101[i])*C6) / m_dx[1];
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
            vector<Scalar> f_000(numComp, zero);
            vector<Scalar> f_001(numComp, zero);
            vector<Scalar> f_010(numComp, zero);
            vector<Scalar> f_011(numComp, zero);
            vector<Scalar> f_100(numComp, zero);
            vector<Scalar> f_101(numComp, zero);
            vector<Scalar> f_110(numComp, zero);
            vector<Scalar> f_111(numComp, zero);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < NE2; ++k2) {
                    for (index_t k1=0; k1 < NE1; ++k1) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(0,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(0,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(0,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(0,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(1,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(1,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(1,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(1,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        Scalar* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,NE1), zero);
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
                for (index_t k2=0; k2 < NE2; ++k2) {
                    for (index_t k1=0; k1 < NE1; ++k1) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        Scalar* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,NE1), zero);
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
                for (index_t k2=0; k2 < NE2; ++k2) {
                    for (index_t k0=0; k0 < NE0; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,0,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,0,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,0,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,0,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        Scalar* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,NE0), zero);
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
                for (index_t k2=0; k2 < NE2; ++k2) {
                    for (index_t k0=0; k0 < NE0; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-2,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-2,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-2,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-2,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        Scalar* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,NE0), zero);
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
                for (index_t k1=0; k1 < NE1; ++k1) {
                    for (index_t k0=0; k0 < NE0; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,0, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,0, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,0, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,0, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        Scalar* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,NE0), zero);
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
                for (index_t k1=0; k1 < NE1; ++k1) {
                    for (index_t k0=0; k0 < NE0; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
                        Scalar* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,NE0), zero);
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

// instantiate our two supported versions
template
void Brick::assembleGradientImpl<real_t>(escript::Data& out,
                                         const escript::Data& in) const;

template
void Brick::assembleGradientImpl<cplx_t>(escript::Data& out,
                                         const escript::Data& in) const;

//protected
void Brick::assembleIntegrate(vector<real_t>& integrals, const escript::Data& arg) const
{
    assembleIntegrateImpl<real_t>(integrals, arg);
}

//protected
void Brick::assembleIntegrate(vector<cplx_t>& integrals, const escript::Data& arg) const
{
    assembleIntegrateImpl<cplx_t>(integrals, arg);
}

//private
template<typename Scalar>
void Brick::assembleIntegrateImpl(vector<Scalar>& integrals, const escript::Data& arg) const
{
    const dim_t numComp = arg.getDataPointSize();
    const index_t left = (m_offset[0]==0 ? 0 : 1);
    const index_t bottom = (m_offset[1]==0 ? 0 : 1);
    const index_t front = (m_offset[2]==0 ? 0 : 1);
    const int fs = arg.getFunctionSpace().getTypeCode();
    const Scalar zero = static_cast<Scalar>(0);

    if(fs == Points ) {
        for (index_t k1 = 0; k1 < m_diracPoints.size(); k1++) { //only for this rank
            const Scalar* f  = arg.getSampleDataRO(k1, zero);
            for (index_t i = 0; i < numComp; ++i) {
                integrals[i]+=f[i];
            }
        }
    } else if (fs == Elements && arg.actsExpanded()) {
        const real_t w_0 = m_dx[0]*m_dx[1]*m_dx[2]/8.;
#pragma omp parallel
        {
            vector<Scalar> int_local(numComp, zero);
#pragma omp for nowait
            for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const Scalar* f = arg.getSampleDataRO(INDEX3(k0, k1, k2, m_NE[0], m_NE[1]), zero);
                        for (index_t i = 0; i < numComp; ++i) {
                            const Scalar f_0 = f[INDEX2(i,0,numComp)];
                            const Scalar f_1 = f[INDEX2(i,1,numComp)];
                            const Scalar f_2 = f[INDEX2(i,2,numComp)];
                            const Scalar f_3 = f[INDEX2(i,3,numComp)];
                            const Scalar f_4 = f[INDEX2(i,4,numComp)];
                            const Scalar f_5 = f[INDEX2(i,5,numComp)];
                            const Scalar f_6 = f[INDEX2(i,6,numComp)];
                            const Scalar f_7 = f[INDEX2(i,7,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3+f_4+f_5+f_6+f_7)*w_0;
                        }  // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of k2 loop

#pragma omp critical
            for (index_t i = 0; i < numComp; i++)
                integrals[i] += int_local[i];
        } // end of parallel section

    } else if (fs==ReducedElements || (fs==Elements && !arg.actsExpanded())) {
        const real_t w_0 = m_dx[0]*m_dx[1]*m_dx[2];
#pragma omp parallel
        {
            vector<Scalar> int_local(numComp, zero);
#pragma omp for nowait
            for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const Scalar* f = arg.getSampleDataRO(INDEX3(k0, k1, k2, m_NE[0], m_NE[1]), zero);
                        for (index_t i = 0; i < numComp; ++i) {
                            int_local[i] += f[i]*w_0;
                        }  // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of k2 loop

#pragma omp critical
            for (index_t i = 0; i < numComp; i++)
                integrals[i] += int_local[i];
        } // end of parallel section

    } else if (fs == FaceElements && arg.actsExpanded()) {
        const real_t w_0 = m_dx[1]*m_dx[2]/4.;
        const real_t w_1 = m_dx[0]*m_dx[2]/4.;
        const real_t w_2 = m_dx[0]*m_dx[1]/4.;
#pragma omp parallel
        {
            vector<Scalar> int_local(numComp, zero);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                    for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                        const Scalar* f = arg.getSampleDataRO(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]), zero);
                        for (index_t i = 0; i < numComp; ++i) {
                            const Scalar f_0 = f[INDEX2(i,0,numComp)];
                            const Scalar f_1 = f[INDEX2(i,1,numComp)];
                            const Scalar f_2 = f[INDEX2(i,2,numComp)];
                            const Scalar f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i] += (f_0+f_1+f_2+f_3)*w_0;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                    for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                        const Scalar* f = arg.getSampleDataRO(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]), zero);
                        for (index_t i = 0; i < numComp; ++i) {
                            const Scalar f_0 = f[INDEX2(i,0,numComp)];
                            const Scalar f_1 = f[INDEX2(i,1,numComp)];
                            const Scalar f_2 = f[INDEX2(i,2,numComp)];
                            const Scalar f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_0;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const Scalar* f = arg.getSampleDataRO(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]), zero);
                        for (index_t i = 0; i < numComp; ++i) {
                            const Scalar f_0 = f[INDEX2(i,0,numComp)];
                            const Scalar f_1 = f[INDEX2(i,1,numComp)];
                            const Scalar f_2 = f[INDEX2(i,2,numComp)];
                            const Scalar f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_1;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const Scalar* f = arg.getSampleDataRO(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]), zero);
                        for (index_t i = 0; i < numComp; ++i) {
                            const Scalar f_0 = f[INDEX2(i,0,numComp)];
                            const Scalar f_1 = f[INDEX2(i,1,numComp)];
                            const Scalar f_2 = f[INDEX2(i,2,numComp)];
                            const Scalar f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i] += (f_0+f_1+f_2+f_3)*w_1;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const Scalar* f = arg.getSampleDataRO(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]), zero);
                        for (index_t i = 0; i < numComp; ++i) {
                            const Scalar f_0 = f[INDEX2(i,0,numComp)];
                            const Scalar f_1 = f[INDEX2(i,1,numComp)];
                            const Scalar f_2 = f[INDEX2(i,2,numComp)];
                            const Scalar f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i] += (f_0+f_1+f_2+f_3)*w_2;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const Scalar* f = arg.getSampleDataRO(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]), zero);
                        for (index_t i = 0; i < numComp; ++i) {
                            const Scalar f_0 = f[INDEX2(i,0,numComp)];
                            const Scalar f_1 = f[INDEX2(i,1,numComp)];
                            const Scalar f_2 = f[INDEX2(i,2,numComp)];
                            const Scalar f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_2;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

#pragma omp critical
            for (index_t i = 0; i < numComp; i++)
                integrals[i] += int_local[i];
        } // end of parallel section

    } else if (fs==ReducedFaceElements || (fs==FaceElements && !arg.actsExpanded())) {
        const real_t w_0 = m_dx[1]*m_dx[2];
        const real_t w_1 = m_dx[0]*m_dx[2];
        const real_t w_2 = m_dx[0]*m_dx[1];
#pragma omp parallel
        {
            vector<Scalar> int_local(numComp, zero);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                    for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                        const Scalar* f = arg.getSampleDataRO(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]), zero);
                        for (index_t i = 0; i < numComp; ++i) {
                            int_local[i] += f[i]*w_0;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                    for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                        const Scalar* f = arg.getSampleDataRO(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]), zero);
                        for (index_t i = 0; i < numComp; ++i) {
                            int_local[i] += f[i]*w_0;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const Scalar* f = arg.getSampleDataRO(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]), zero);
                        for (index_t i = 0; i < numComp; ++i) {
                            int_local[i] += f[i]*w_1;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const Scalar* f = arg.getSampleDataRO(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]), zero);
                        for (index_t i = 0; i < numComp; ++i) {
                            int_local[i] += f[i]*w_1;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const Scalar* f = arg.getSampleDataRO(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]), zero);
                        for (index_t i = 0; i < numComp; ++i) {
                            int_local[i] += f[i]*w_2;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
                        const Scalar* f = arg.getSampleDataRO(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]), zero);
                        for (index_t i = 0; i < numComp; ++i) {
                            int_local[i] += f[i]*w_2;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

#pragma omp critical
            for (index_t i = 0; i < numComp; i++)
                integrals[i] += int_local[i];
        } // end of parallel section
    } // function space selector
}

//protected
IndexVector Brick::getDiagonalIndices(bool upperOnly) const
{
    IndexVector ret;
    // only store non-negative indices if requested
    if (upperOnly)
        ret.resize(14);
    else
        ret.resize(27);

    const dim_t nDOF0 = (m_gNE[0]+1)/m_NX[0];
    const dim_t nDOF1 = (m_gNE[1]+1)/m_NX[1];
    size_t idx = 0;
    for (int i2=-1; i2<2; i2++) {
        for (int i1=-1; i1<2; i1++) {
            for (int i0=-1; i0<2; i0++) {
                const int index = i2*nDOF0*nDOF1 + i1*nDOF0 + i0;
                if (!upperOnly || index >= 0)
                    ret[idx++] = index;
            }
        }
    }

    return ret;
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
                std::copy(src, src+numComp, out.getSampleDataRW(k+j*nDOF0+i*nDOF0*nDOF1));
            }
        }
    }
}

#ifdef ESYS_HAVE_TRILINOS
//protected
esys_trilinos::const_TrilinosGraph_ptr Brick::getTrilinosGraph() const
{
    if (m_graph.is_null()) {
        m_graph = createTrilinosGraph(m_dofId, m_nodeId);
    }
    return m_graph;
}
#endif

#ifdef ESYS_HAVE_PASO
//protected
paso::SystemMatrixPattern_ptr Brick::getPasoMatrixPattern(
                                                    bool reducedRowOrder,
                                                    bool reducedColOrder) const
{
    if (m_pattern.get())
        return m_pattern;

    // first call to this method -> create the pattern, then return it
    paso::Connector_ptr conn(getPasoConnector());
    const dim_t nDOF0 = (m_gNE[0]+1)/m_NX[0];
    const dim_t nDOF1 = (m_gNE[1]+1)/m_NX[1];
    const dim_t nDOF2 = (m_gNE[2]+1)/m_NX[2];
    const dim_t numDOF = nDOF0*nDOF1*nDOF2;
    const dim_t numShared = conn->send->numSharedComponents;
    const index_t* sendShared = conn->send->shared;
    const int x = m_mpiInfo->rank%m_NX[0];
    const int y = m_mpiInfo->rank%(m_NX[0]*m_NX[1])/m_NX[0];
    const int z = m_mpiInfo->rank/(m_NX[0]*m_NX[1]);

    // these are for the couple blocks
    vector<IndexVector> colIndices(numDOF);
    vector<IndexVector> rowIndices(numShared);

    for (dim_t i=0; i < conn->send->neighbour.size(); i++) {
        const dim_t start = conn->send->offsetInShared[i];
        const dim_t end = conn->send->offsetInShared[i+1];
        // location of neighbour rank relative to this rank
        const int xDiff = conn->send->neighbour[i]%m_NX[0] - x;
        const int yDiff = conn->send->neighbour[i]%(m_NX[0]*m_NX[1])/m_NX[0] - y;
        const int zDiff = conn->send->neighbour[i]/(m_NX[0]*m_NX[1]) - z;

        if (xDiff==0 && yDiff==0) {
            // sharing front or back plane
            for (dim_t j = start; j < end; j++) {
                const dim_t i0 = (j-start)%nDOF0;
                const dim_t i1 = (j-start)/nDOF0;
                if (i0 > 0) {
                    if (i1 > 0)
                        doublyLink(colIndices, rowIndices, sendShared[j-1-nDOF0], j);
                    doublyLink(colIndices, rowIndices, sendShared[j-1], j);
                    if (i1 < nDOF1-1)
                        doublyLink(colIndices, rowIndices, sendShared[j-1+nDOF0], j);
                }
                if (i1 > 0)
                    doublyLink(colIndices, rowIndices, sendShared[j-nDOF0], j);
                doublyLink(colIndices, rowIndices, sendShared[j], j);
                if (i1 < nDOF1-1)
                    doublyLink(colIndices, rowIndices, sendShared[j+nDOF0], j);
                if (i0 < nDOF0-1) {
                    if (i1 > 0)
                        doublyLink(colIndices, rowIndices, sendShared[j+1-nDOF0], j);
                    doublyLink(colIndices, rowIndices, sendShared[j+1], j);
                    if (i1 < nDOF1-1)
                        doublyLink(colIndices, rowIndices, sendShared[j+1+nDOF0], j);
                }
            }
        } else if (xDiff==0 && zDiff==0) {
            // sharing top or bottom plane
            for (dim_t j = start; j < end; j++) {
                const dim_t i0 = (j-start)%nDOF0;
                const dim_t i1 = (j-start)/nDOF0;
                if (i0 > 0) {
                    if (i1 > 0)
                        doublyLink(colIndices, rowIndices, sendShared[j]-1-nDOF0*nDOF1, j);
                    doublyLink(colIndices, rowIndices, sendShared[j]-1, j);
                    if (i1 < nDOF2-1)
                        doublyLink(colIndices, rowIndices, sendShared[j]-1+nDOF0*nDOF1, j);
                }
                if (i1 > 0)
                    doublyLink(colIndices, rowIndices, sendShared[j]-nDOF0*nDOF1, j);
                doublyLink(colIndices, rowIndices, sendShared[j], j);
                if (i1 < nDOF2-1)
                    doublyLink(colIndices, rowIndices, sendShared[j]+nDOF0*nDOF1, j);
                if (i0 < nDOF0-1) {
                    if (i1 > 0)
                        doublyLink(colIndices, rowIndices, sendShared[j]+1-nDOF0*nDOF1, j);
                    doublyLink(colIndices, rowIndices, sendShared[j]+1, j);
                    if (i1 < nDOF2-1)
                        doublyLink(colIndices, rowIndices, sendShared[j]+1+nDOF0*nDOF1, j);
                }
            }
        } else if (yDiff==0 && zDiff==0) {
            // sharing left or right plane
            for (dim_t j = start; j < end; j++) {
                const dim_t i0 = (j-start)%nDOF1;
                const dim_t i1 = (j-start)/nDOF1;
                if (i0 > 0) {
                    if (i1 > 0)
                        doublyLink(colIndices, rowIndices, sendShared[j]-nDOF0-nDOF0*nDOF1, j);
                    doublyLink(colIndices, rowIndices, sendShared[j]-nDOF0, j);
                    if (i1 < nDOF2-1)
                        doublyLink(colIndices, rowIndices, sendShared[j]-nDOF0+nDOF0*nDOF1, j);
                }
                if (i1 > 0)
                    doublyLink(colIndices, rowIndices, sendShared[j]-nDOF0*nDOF1, j);
                doublyLink(colIndices, rowIndices, sendShared[j], j);
                if (i1 < nDOF2-1)
                    doublyLink(colIndices, rowIndices, sendShared[j]+nDOF0*nDOF1, j);
                if (i0 < nDOF1-1) {
                    if (i1 > 0)
                        doublyLink(colIndices, rowIndices, sendShared[j]+nDOF0-nDOF0*nDOF1, j);
                    doublyLink(colIndices, rowIndices, sendShared[j]+nDOF0, j);
                    if (i1 < nDOF2-1)
                        doublyLink(colIndices, rowIndices, sendShared[j]+nDOF0+nDOF0*nDOF1, j);
                }
            }
        } else if (xDiff == 0) {
            // sharing an edge in x direction
            for (dim_t j = start; j < end; j++) {
                if (j > start)
                    doublyLink(colIndices, rowIndices, sendShared[j]-1, j);
                doublyLink(colIndices, rowIndices, sendShared[j], j);
                if (j < end-1)
                    doublyLink(colIndices, rowIndices, sendShared[j]+1, j);
            }
        } else if (yDiff == 0) {
            // sharing an edge in y direction
            for (dim_t j = start; j < end; j++) {
                if (j > start)
                    doublyLink(colIndices, rowIndices, sendShared[j]-nDOF0, j);
                doublyLink(colIndices, rowIndices, sendShared[j], j);
                if (j < end-1)
                    doublyLink(colIndices, rowIndices, sendShared[j]+nDOF0, j);
            }
        } else if (zDiff == 0) {
            // sharing an edge in z direction
            for (dim_t j = start; j < end; j++) {
                if (j > start)
                    doublyLink(colIndices, rowIndices, sendShared[j]-nDOF0*nDOF1, j);
                doublyLink(colIndices, rowIndices, sendShared[j], j);
                if (j < end-1)
                    doublyLink(colIndices, rowIndices, sendShared[j]+nDOF0*nDOF1, j);
            }
        } else {
            // sharing a node
            doublyLink(colIndices, rowIndices, sendShared[start], start);
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

    // useful debug output
    /*
    std::cout << "--- colIndices ---" << std::endl;
    for (size_t i=0; i<colIndices.size(); i++) {
        std::cout << "colIndices[" << i << "].size()=" << colIndices[i].size() << std::endl;
    }
    std::cout << "--- rowIndices ---" << std::endl;
    for (size_t i=0; i<rowIndices.size(); i++) {
        std::cout << "rowIndices[" << i << "].size()=" << rowIndices[i].size() << std::endl;
    }
    */
    /*
    std::cout << "--- main_pattern ---" << std::endl;
    std::cout << "M=" << mainPattern->numOutput << ", N=" << mainPattern->numInput << std::endl;
    for (size_t i=0; i<mainPattern->numOutput+1; i++) {
        std::cout << "ptr[" << i << "]=" << mainPattern->ptr[i] << std::endl;
    }
    for (size_t i=0; i<mainPattern->ptr[mainPattern->numOutput]; i++) {
        std::cout << "index[" << i << "]=" << mainPattern->index[i] << std::endl;
    }
    */
    /*
    std::cout << "--- colCouple_pattern ---" << std::endl;
    std::cout << "M=" << colPattern->numOutput << ", N=" << colPattern->numInput << std::endl;
    for (size_t i=0; i<colPattern->numOutput+1; i++) {
        std::cout << "ptr[" << i << "]=" << colPattern->ptr[i] << std::endl;
    }
    for (size_t i=0; i<colPattern->ptr[colPattern->numOutput]; i++) {
        std::cout << "index[" << i << "]=" << colPattern->index[i] << std::endl;
    }
    */
    /*
    std::cout << "--- rowCouple_pattern ---" << std::endl;
    std::cout << "M=" << rowPattern->numOutput << ", N=" << rowPattern->numInput << std::endl;
    for (size_t i=0; i<rowPattern->numOutput+1; i++) {
        std::cout << "ptr[" << i << "]=" << rowPattern->ptr[i] << std::endl;
    }
    for (size_t i=0; i<rowPattern->ptr[rowPattern->numOutput]; i++) {
        std::cout << "index[" << i << "]=" << rowPattern->index[i] << std::endl;
    }
    */

    return m_pattern;
}
#endif // ESYS_HAVE_PASO

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

    const index_t left = (m_offset[0]==0 ? 0 : 1);
    const index_t bottom = (m_offset[1]==0 ? 0 : 1);
    const index_t front = (m_offset[2]==0 ? 0 : 1);
    const dim_t nDOF0 = (m_gNE[0]+1)/m_NX[0];
    const dim_t nDOF1 = (m_gNE[1]+1)/m_NX[1];
    const dim_t nDOF2 = (m_gNE[2]+1)/m_NX[2];
    const dim_t NN0 = m_NN[0];
    const dim_t NN1 = m_NN[1];
    const dim_t NN2 = m_NN[2];
    const dim_t NE0 = m_NE[0];
    const dim_t NE1 = m_NE[1];
    const dim_t NE2 = m_NE[2];

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
        for (dim_t i=0; i<NN0; i++) {
            m_nodeId[i] = globalNodeId(i, 0, 0); // LF
            m_nodeId[NN0*(NN1-1)+i] = globalNodeId(i, NN1-1, 0); // UF
            m_nodeId[NN0*NN1*(NN2-1)+i] = globalNodeId(i, 0, NN2-1); // LB
            m_nodeId[NN0*NN1*NN2-NN0+i] = globalNodeId(i, NN1-1, NN2-1); // UB
        }
        // edges in y-direction, without corners
#pragma omp for nowait
        for (dim_t i=1; i<NN1-1; i++) {
            m_nodeId[NN0*i] = globalNodeId(0, i, 0); // FL
            m_nodeId[NN0*(i+1)-1] = globalNodeId(NN0-1, i, 0); // FR
            m_nodeId[NN0*NN1*(NN2-1)+NN0*i] = globalNodeId(0, i, NN2-1); // BL
            m_nodeId[NN0*NN1*(NN2-1)+NN0*(i+1)-1] = globalNodeId(NN0-1, i, NN2-1); // BR
        }
        // edges in z-direction, without corners
#pragma omp for
        for (dim_t i=1; i<NN2-1; i++) {
            m_nodeId[NN0*NN1*i] = globalNodeId(0, 0, i); // LL
            m_nodeId[NN0*NN1*i+NN0-1] = globalNodeId(NN0-1, 0, i); // LR
            m_nodeId[NN0*NN1*(i+1)-NN0] = globalNodeId(0, NN1-1, i); // UL
            m_nodeId[NN0*NN1*(i+1)-1] = globalNodeId(NN0-1, NN1-1, i); // UR
        }
        // implicit barrier here because some node IDs will be overwritten
        // below

        // populate degrees of freedom and own nodes (identical id)
#pragma omp for nowait
        for (dim_t i=0; i<nDOF2; i++) {
            for (dim_t j=0; j<nDOF1; j++) {
                for (dim_t k=0; k<nDOF0; k++) {
                    const index_t nodeIdx=k+left+(j+bottom)*NN0+(i+front)*NN0*NN1;
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
                    const index_t nodeIdx=(j+bottom)*NN0+(i+front)*NN0*NN1;
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
                    const index_t nodeIdx=(j+bottom+1)*NN0-1+(i+front)*NN0*NN1;
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
                    const index_t nodeIdx=k+left+(i+front)*NN0*NN1;
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
                    const index_t nodeIdx=k+left+(i+front)*NN0*NN1+NN0*(NN1-1);
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
                    const index_t nodeIdx=k+left+(j+bottom)*NN0;
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
                    const index_t nodeIdx=k+left+(j+bottom)*NN0+NN0*NN1*(NN2-1);
                    const index_t dofId=k+j*nDOF0;
                    m_nodeId[nodeIdx]
                        = m_nodeDistribution[m_mpiInfo->rank+m_NX[0]*m_NX[1]]+dofId;
                }
            }
        }

        // populate element id's
#pragma omp for nowait
        for (dim_t i2=0; i2<NE2; i2++) {
            for (dim_t i1=0; i1<NE1; i1++) {
                for (dim_t i0=0; i0<NE0; i0++) {
                    m_elementId[i0+i1*NE0+i2*NE0*NE1] =
                        (m_offset[2]+i2)*m_gNE[0]*m_gNE[1]
                        +(m_offset[1]+i1)*m_gNE[0]
                        +m_offset[0]+i0;
                }
            }
        }

        // face elements
#pragma omp for
        for (dim_t k=0; k<NFE; k++)
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

    populateDofMap();
}

//protected
vector<IndexVector> Brick::getConnections(bool includeShared) const
{
    // returns a vector v of size numDOF where v[i] is a vector with indices
    // of DOFs connected to i (up to 27 in 3D).
    // In other words this method returns the occupied (local) matrix columns
    // for all (local) matrix rows.
    // If includeShared==true then connections to non-owned DOFs are also
    // returned (i.e. indices of the column couplings)
    const dim_t nDOF0 = getNumDOFInAxis(0);
    const dim_t nDOF1 = getNumDOFInAxis(1);
    const dim_t nDOF2 = getNumDOFInAxis(2);
    const dim_t numMatrixRows = nDOF0*nDOF1*nDOF2;
    vector<IndexVector> indices(numMatrixRows);

    if (includeShared) {
        const index_t left = getFirstInDim(0);
        const index_t bottom = getFirstInDim(1);
        const index_t front = getFirstInDim(2);
        const dim_t NN0 = m_NN[0];
        const dim_t NN1 = m_NN[1];
        const dim_t NN2 = m_NN[2];
#pragma omp parallel for
        for (index_t i=0; i < numMatrixRows; i++) {
            const index_t x = left + i % nDOF0;
            const index_t y = bottom + i % (nDOF0*nDOF1)/nDOF0;
            const index_t z = front + i / (nDOF0*nDOF1);
            // loop through potential neighbours and add to index if positions
            // are within bounds
            for (int i2=z-1; i2<z+2; i2++) {
                for (int i1=y-1; i1<y+2; i1++) {
                    for (int i0=x-1; i0<x+2; i0++) {
                        if (i0>=0 && i1>=0 && i2>=0
                                && i0<NN0 && i1<NN1 && i2<NN2) {
                            indices[i].push_back(m_dofMap[i2*NN0*NN1+i1*NN0+i0]);
                        }
                    }
                }
            }
            sort(indices[i].begin(), indices[i].end());
        }
    } else {
#pragma omp parallel for
        for (index_t i=0; i < numMatrixRows; i++) {
            const index_t x = i % nDOF0;
            const index_t y = i % (nDOF0*nDOF1)/nDOF0;
            const index_t z = i / (nDOF0*nDOF1);
            // loop through potential neighbours and add to index if positions
            // are within bounds
            for (int i2=z-1; i2<z+2; i2++) {
                for (int i1=y-1; i1<y+2; i1++) {
                    for (int i0=x-1; i0<x+2; i0++) {
                        if (i0>=0 && i1>=0 && i2>=0
                                && i0<nDOF0 && i1<nDOF1 && i2<nDOF2) {
                            indices[i].push_back(i2*nDOF0*nDOF1+i1*nDOF0+i0);
                        }
                    }
                }
            }
        }
    }
    return indices;
}

//private
void Brick::populateDofMap()
{
    const dim_t nDOF0 = getNumDOFInAxis(0);
    const dim_t nDOF1 = getNumDOFInAxis(1);
    const dim_t nDOF2 = getNumDOFInAxis(2);
    const index_t left = getFirstInDim(0);
    const index_t bottom = getFirstInDim(1);
    const index_t front = getFirstInDim(2);

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

    const dim_t numDOF=nDOF0*nDOF1*nDOF2;
    RankVector neighbour;
    IndexVector offsetInShared(1,0);
    IndexVector sendShared, recvShared;
    dim_t numShared=0;
    const int x=m_mpiInfo->rank%m_NX[0];
    const int y=m_mpiInfo->rank%(m_NX[0]*m_NX[1])/m_NX[0];
    const int z=m_mpiInfo->rank/(m_NX[0]*m_NX[1]);

    // build list of shared components and neighbours by looping through
    // all potential neighbouring ranks and checking if positions are
    // within bounds
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
                            const dim_t firstDOF=(i2==-1 ? i*nDOF0
                                    : i*nDOF0 + nDOF0*nDOF1*(nDOF2-1));
                            const dim_t firstNode=(i2==-1 ? left+(i+bottom)*m_NN[0]
                                    : left+(i+bottom)*m_NN[0]+m_NN[0]*m_NN[1]*(m_NN[2]-1));
                            for (dim_t j=0; j<nDOF0; j++, numShared++) {
                                sendShared.push_back(firstDOF+j);
                                recvShared.push_back(numDOF+numShared);
                                m_dofMap[firstNode+j]=numDOF+numShared;
                            }
                        }
                    } else if (i0==0 && i2==0) {
                        // sharing top or bottom plane
                        offsetInShared.push_back(offsetInShared.back()+nDOF0*nDOF2);
                        for (dim_t i=0; i<nDOF2; i++) {
                            const dim_t firstDOF=(i1==-1 ? i*nDOF0*nDOF1
                                    : nDOF0*((i+1)*nDOF1-1));
                            const dim_t firstNode=(i1==-1 ?
                                    left+(i+front)*m_NN[0]*m_NN[1]
                                    : left+m_NN[0]*((i+1+front)*m_NN[1]-1));
                            for (dim_t j=0; j<nDOF0; j++, numShared++) {
                                sendShared.push_back(firstDOF+j);
                                recvShared.push_back(numDOF+numShared);
                                m_dofMap[firstNode+j]=numDOF+numShared;
                            }
                        }
                    } else if (i1==0 && i2==0) {
                        // sharing left or right plane
                        offsetInShared.push_back(offsetInShared.back()+nDOF1*nDOF2);
                        for (dim_t i=0; i<nDOF2; i++) {
                            const dim_t firstDOF=(i0==-1 ? i*nDOF0*nDOF1
                                    : nDOF0*(1+i*nDOF1)-1);
                            const dim_t firstNode=(i0==-1 ?
                                    (bottom+(i+front)*m_NN[1])*m_NN[0]
                                    : (bottom+1+(i+front)*m_NN[1])*m_NN[0]-1);
                            for (dim_t j=0; j<nDOF1; j++, numShared++) {
                                sendShared.push_back(firstDOF+j*nDOF0);
                                recvShared.push_back(numDOF+numShared);
                                m_dofMap[firstNode+j*m_NN[0]]=numDOF+numShared;
                            }
                        }
                    } else if (i0==0) {
                        // sharing an edge in x direction
                        offsetInShared.push_back(offsetInShared.back()+nDOF0);
                        const dim_t firstDOF=(i1+1)/2*nDOF0*(nDOF1-1)
                                           +(i2+1)/2*nDOF0*nDOF1*(nDOF2-1);
                        const dim_t firstNode=left+(i1+1)/2*m_NN[0]*(m_NN[1]-1)
                                            +(i2+1)/2*m_NN[0]*m_NN[1]*(m_NN[2]-1);
                        for (dim_t i=0; i<nDOF0; i++, numShared++) {
                            sendShared.push_back(firstDOF+i);
                            recvShared.push_back(numDOF+numShared);
                            m_dofMap[firstNode+i]=numDOF+numShared;
                        }
                    } else if (i1==0) {
                        // sharing an edge in y direction
                        offsetInShared.push_back(offsetInShared.back()+nDOF1);
                        const dim_t firstDOF=(i0+1)/2*(nDOF0-1)
                                           +(i2+1)/2*nDOF0*nDOF1*(nDOF2-1);
                        const dim_t firstNode=bottom*m_NN[0]
                                            +(i0+1)/2*(m_NN[0]-1)
                                            +(i2+1)/2*m_NN[0]*m_NN[1]*(m_NN[2]-1);
                        for (dim_t i=0; i<nDOF1; i++, numShared++) {
                            sendShared.push_back(firstDOF+i*nDOF0);
                            recvShared.push_back(numDOF+numShared);
                            m_dofMap[firstNode+i*m_NN[0]]=numDOF+numShared;
                        }
                    } else if (i2==0) {
                        // sharing an edge in z direction
                        offsetInShared.push_back(offsetInShared.back()+nDOF2);
                        const dim_t firstDOF=(i0+1)/2*(nDOF0-1)
                                           +(i1+1)/2*nDOF0*(nDOF1-1);
                        const dim_t firstNode=front*m_NN[0]*m_NN[1]
                                            +(i0+1)/2*(m_NN[0]-1)
                                            +(i1+1)/2*m_NN[0]*(m_NN[1]-1);
                        for (dim_t i=0; i<nDOF2; i++, numShared++) {
                            sendShared.push_back(firstDOF+i*nDOF0*nDOF1);
                            recvShared.push_back(numDOF+numShared);
                            m_dofMap[firstNode+i*m_NN[0]*m_NN[1]]=numDOF+numShared;
                        }
                    } else {
                        // sharing a node
                        const dim_t dof = (i0+1)/2*(nDOF0-1)
                                       +(i1+1)/2*nDOF0*(nDOF1-1)
                                       +(i2+1)/2*nDOF0*nDOF1*(nDOF2-1);
                        const dim_t node = (i0+1)/2*(m_NN[0]-1)
                                        +(i1+1)/2*m_NN[0]*(m_NN[1]-1)
                                        +(i2+1)/2*m_NN[0]*m_NN[1]*(m_NN[2]-1);
                        offsetInShared.push_back(offsetInShared.back()+1);
                        sendShared.push_back(dof);
                        recvShared.push_back(numDOF+numShared);
                        m_dofMap[node] = numDOF+numShared;
                        ++numShared;
                    }
                }
            }
        }
    }

#ifdef ESYS_HAVE_PASO
    createPasoConnector(neighbour, offsetInShared, offsetInShared, sendShared,
                        recvShared);
#endif

    // useful debug output
    /*
    std::cout << "--- rcv_shcomp ---" << std::endl;
    std::cout << "numDOF=" << numDOF << ", numNeighbours=" << neighbour.size() << std::endl;
    for (size_t i=0; i<neighbour.size(); i++) {
        std::cout << "neighbour[" << i << "]=" << neighbour[i]
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
void Brick::addToMatrixAndRHS(AbstractSystemMatrix* S, escript::Data& F,
         const vector<Scalar>& EM_S, const vector<Scalar>& EM_F, bool addS,
         bool addF, index_t firstNode, int nEq, int nComp) const
{
    IndexVector rowIndex(8);
    rowIndex[0] = m_dofMap[firstNode];
    rowIndex[1] = m_dofMap[firstNode+1];
    rowIndex[2] = m_dofMap[firstNode+m_NN[0]];
    rowIndex[3] = m_dofMap[firstNode+m_NN[0]+1];
    rowIndex[4] = m_dofMap[firstNode+m_NN[0]*m_NN[1]];
    rowIndex[5] = m_dofMap[firstNode+m_NN[0]*m_NN[1]+1];
    rowIndex[6] = m_dofMap[firstNode+m_NN[0]*(m_NN[1]+1)];
    rowIndex[7] = m_dofMap[firstNode+m_NN[0]*(m_NN[1]+1)+1];
    if (addF) {
        Scalar* F_p = F.getSampleDataRW(0, static_cast<Scalar>(0));
        for (index_t i = 0; i < rowIndex.size(); i++) {
            if (rowIndex[i] < getNumDOF()) {
                for (int eq = 0; eq < nEq; eq++) {
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
void Brick::addToMatrixAndRHS<real_t>(AbstractSystemMatrix* S, escript::Data& F,
         const vector<real_t>& EM_S, const vector<real_t>& EM_F, bool addS,
         bool addF, index_t firstNode, int nEq, int nComp) const;

template
void Brick::addToMatrixAndRHS<cplx_t>(AbstractSystemMatrix* S, escript::Data& F,
         const vector<cplx_t>& EM_S, const vector<cplx_t>& EM_F, bool addS,
         bool addF, index_t firstNode, int nEq, int nComp) const;

//protected
void Brick::interpolateNodesOnElements(escript::Data& out,
                                       const escript::Data& in,
                                       bool reduced) const
{
    if (out.isComplex()!=in.isComplex())
    {
        throw RipleyException("Programmer Error: in and out parameters do not have the same complexity.");
    }
    if (out.isComplex())
    {
        interpolateNodesOnElementsWorker(out, in, reduced, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateNodesOnElementsWorker(out, in, reduced, escript::DataTypes::real_t(0));
    }
}
//protected
void Brick::interpolateNodesOnFaces(escript::Data& out, const escript::Data& in,
                                    bool reduced) const
{
    if (out.isComplex()!=in.isComplex())
    {
        throw RipleyException("Programmer Error: in and out parameters do not have the same complexity.");
    }
    if (out.isComplex())
    {
        interpolateNodesOnFacesWorker(out, in, reduced, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateNodesOnFacesWorker(out, in, reduced, escript::DataTypes::real_t(0));
    }
}


//private
template <typename S>
void Brick::interpolateNodesOnElementsWorker(escript::Data& out,
                                       const escript::Data& in,
                                       bool reduced, S sentinel) const
{
    const dim_t numComp = in.getDataPointSize();
    if (reduced) {
        out.requireWrite();
#pragma omp parallel
        {
            vector<S> f_000(numComp);
            vector<S> f_001(numComp);
            vector<S> f_010(numComp);
            vector<S> f_011(numComp);
            vector<S> f_100(numComp);
            vector<S> f_101(numComp);
            vector<S> f_110(numComp);
            vector<S> f_111(numComp);
#pragma omp for
            for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        S* o = out.getSampleDataRW(INDEX3(k0,k1,k2,m_NE[0],m_NE[1]), sentinel);
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = (f_000[i] + f_001[i] + f_010[i] + f_011[i] + f_100[i] + f_101[i] + f_110[i] + f_111[i])/static_cast<S>(8);
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of k2 loop
        } // end of parallel section
    } else {
        out.requireWrite();
        const S c0 = .0094373878376559314545;
        const S c1 = .035220810900864519624;
        const S c2 = .13144585576580214704;
        const S c3 = .49056261216234406855;
#pragma omp parallel
        {
            vector<S> f_000(numComp);
            vector<S> f_001(numComp);
            vector<S> f_010(numComp);
            vector<S> f_011(numComp);
            vector<S> f_100(numComp);
            vector<S> f_101(numComp);
            vector<S> f_110(numComp);
            vector<S> f_111(numComp);
#pragma omp for
            for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        S* o = out.getSampleDataRW(INDEX3(k0,k1,k2,m_NE[0],m_NE[1]), sentinel);
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

//private
template <typename S>
void Brick::interpolateNodesOnFacesWorker(escript::Data& out, const escript::Data& in,
                                    bool reduced, S sentinel) const
{
    const dim_t numComp = in.getDataPointSize();
    if (reduced) {
        out.requireWrite();
#pragma omp parallel
        {
            vector<S> f_000(numComp);
            vector<S> f_001(numComp);
            vector<S> f_010(numComp);
            vector<S> f_011(numComp);
            vector<S> f_100(numComp);
            vector<S> f_101(numComp);
            vector<S> f_110(numComp);
            vector<S> f_111(numComp);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(0,k1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(0,k1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(0,k1+1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(0,k1+1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        S* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]), sentinel);
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = (f_000[i] + f_001[i] + f_010[i] + f_011[i])/static_cast<S>(4);
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 0
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        S* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]), sentinel);
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = (f_100[i] + f_101[i] + f_110[i] + f_111[i])/static_cast<S>(4);
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 1
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,0,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,0,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,0,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,0,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        S* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]), sentinel);
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = (f_000[i] + f_001[i] + f_100[i] + f_101[i])/static_cast<S>(4);
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 2
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        S* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]), sentinel);
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = (f_010[i] + f_011[i] + f_110[i] + f_111[i])/static_cast<S>(4);
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 3
            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,0, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,0, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,0, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,0, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        S* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]), sentinel);
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = (f_000[i] + f_010[i] + f_100[i] + f_110[i])/static_cast<S>(4);
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of face 4
            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                    for (index_t k0=0; k0 < m_NE[0]; ++k0) {
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        S* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]), sentinel);
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = (f_001[i] + f_011[i] + f_101[i] + f_111[i])/static_cast<S>(4);
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of face 5
        } // end of parallel section
    } else {
        out.requireWrite();
        const S c0 = 0.044658198738520451079;
        const S c1 = 0.16666666666666666667;
        const S c2 = 0.62200846792814621559;
#pragma omp parallel
        {
            vector<S> f_000(numComp);
            vector<S> f_001(numComp);
            vector<S> f_010(numComp);
            vector<S> f_011(numComp);
            vector<S> f_100(numComp);
            vector<S> f_101(numComp);
            vector<S> f_110(numComp);
            vector<S> f_111(numComp);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE[2]; ++k2) {
                    for (index_t k1=0; k1 < m_NE[1]; ++k1) {
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(0,k1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(0,k1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(0,k1+1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(0,k1+1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        S* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]), sentinel);
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
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        S* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]), sentinel);
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
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,0,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,0,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,0,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,0,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        S* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]), sentinel);
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
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        S* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]), sentinel);
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
                        memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,0, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,0, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,0, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,0, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        S* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]), sentinel);
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
                        memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
                        S* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]), sentinel);
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
    // Calculates a gaussian blur convolution matrix for 3D
    // See wiki article on the subject
    double* get3DGauss(unsigned radius, double sigma)
    {
        double* arr=new double[(radius*2+1)*(radius*2+1)*(radius*2+1)];
        double common = pow(M_1_PI * 0.5 / (sigma*sigma), 3./2);
        double total=0;
        int r=static_cast<int>(radius);
        for (int z=-r;z<=r;++z) {
            for (int y=-r;y<=r;++y) {
                for (int x=-r;x<=r;++x) {
                    arr[(x+r)+(y+r)*(r*2+1)+(z+r)*(r*2+1)*(r*2+1)]=common*exp(-(x*x+y*y+z*z)/(2*sigma*sigma));
                    total+=arr[(x+r)+(y+r)*(r*2+1)+(z+r)*(r*2+1)*(r*2+1)];
                }
            }
        }
        double invtotal=1/total;
        for (size_t p=0;p<(radius*2+1)*(radius*2+1)*(radius*2+1);++p) {
            arr[p]*=invtotal;
        }
        return arr;
    }

    // applies conv to source to get a point.
    // (xp, yp) are the coords in the source matrix not the destination matrix
    double Convolve3D(double* conv, double* source, size_t xp, size_t yp,
                      size_t zp, unsigned radius, size_t width, size_t height)
    {
        size_t bx=xp-radius, by=yp-radius, bz=zp-radius;
        size_t sbase=bx+by*width+bz*width*height;
        double result=0;
        for (int z=0; z<2*radius+1; ++z) {
            for (int y=0; y<2*radius+1; ++y) {
                for (int x=0; x<2*radius+1; ++x) {
                    result += conv[x+y*(2*radius+1)+z*(2*radius+1)*(2*radius+1)] * source[sbase + x+y*width+z*width*height];
                }
            }
        }
        // use this line for "pass-though" (return the centre point value)
//      return source[sbase+(radius)+(radius)*width+(radius)*width*height];
        return result;
    }
}

/* This is a wrapper for filtered (and non-filtered) randoms
 * For detailed doco see randomFillWorker
 */
escript::Data Brick::randomFill(const escript::DataTypes::ShapeType& shape,
                                const escript::FunctionSpace& what,
                                long seed, const bp::tuple& filter) const
{
    int numvals=escript::DataTypes::noValues(shape);
    if (len(filter) > 0 && numvals != 1) {
        throw NotImplementedError("Ripley only supports filters for scalar data.");
    }
    escript::Data res = randomFillWorker(shape, seed, filter);
    if (res.getFunctionSpace()!=what) {
        escript::Data r(res, what);
        return r;
    }
    return res;
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
escript::Data Brick::randomFillWorker(
                        const escript::DataTypes::ShapeType& shape, long seed,
                        const bp::tuple& filter) const
{
    unsigned int radius=0;  // these are only used by gaussian
    double sigma=0.5;

    unsigned int numvals=escript::DataTypes::noValues(shape);

    if (len(filter)==0) {
    // nothing special required here yet
    } else if (len(filter) == 3) {
        bp::extract<string> ex(filter[0]);
        if (!ex.check() || (ex() != "gaussian")) {
            throw ValueError("Unsupported random filter for Brick.");
        }
        bp::extract<unsigned int> ex1(filter[1]);
        if (!ex1.check()) {
            throw ValueError("Radius of gaussian filter must be a positive integer.");
        }
        radius=ex1();
        sigma=0.5;
        bp::extract<double> ex2(filter[2]);
        if (!ex2.check() || (sigma=ex2()) <= 0) {
            throw ValueError("Sigma must be a positive floating point number.");
        }
    } else {
        throw ValueError("Unsupported random filter");
    }

    // number of points in the internal region
    // that is, the ones we need smoothed versions of
    const dim_t internal[3] = { m_NN[0], m_NN[1], m_NN[2] };
    size_t ext[3];
    ext[0]=(size_t)internal[0]+2*radius;  // includes points we need as input
    ext[1]=(size_t)internal[1]+2*radius;  // for smoothing
    ext[2]=(size_t)internal[2]+2*radius;  // for smoothing

    // now we check to see if the radius is acceptable
    // That is, would not cross multiple ranks in MPI

    if (2*radius>=internal[0]-4) {
        throw ValueError("Radius of gaussian filter is too large for X dimension of a rank");
    }
    if (2*radius>=internal[1]-4) {
        throw ValueError("Radius of gaussian filter is too large for Y dimension of a rank");
    }
    if (2*radius>=internal[2]-4) {
        throw ValueError("Radius of gaussian filter is too large for Z dimension of a rank");
    }

    double* src=new double[ext[0]*ext[1]*ext[2]*numvals];
    escript::randomFillArray(seed, src, ext[0]*ext[1]*ext[2]*numvals, m_mpiInfo);

#ifdef ESYS_MPI
    if ((internal[0]<5) || (internal[1]<5) || (internal[2]<5)) {
        // since the dimensions are equal for all ranks, this exception
        // will be thrown on all ranks
        throw ValueError("Random Data in Ripley requires at least five elements per side per rank.");
    }
    dim_t X=m_mpiInfo->rank%m_NX[0];
    dim_t Y=m_mpiInfo->rank%(m_NX[0]*m_NX[1])/m_NX[0];
    dim_t Z=m_mpiInfo->rank/(m_NX[0]*m_NX[1]);
#endif

/*
    // if we wanted to test a repeating pattern
    size_t basex=0;
    size_t basey=0;
    size_t basez=0;
#ifdef ESYS_MPI
    basex=X*m_gNE[0]/m_NX[0];
    basey=Y*m_gNE[1]/m_NX[1];
    basez=Z*m_gNE[2]/m_NX[2];
    std::cout << "basex=" << basex << " basey=" << basey << " basez=" << basez << std::endl;
#endif
    escript::patternFillArray(1, ext[0],ext[1],ext[2], src, 4, basex, basey, basez, numvals);
*/

#ifdef ESYS_MPI
    BlockGrid grid(m_NX[0]-1, m_NX[1]-1, m_NX[2]-1);
    // it's +2 not +1 because a whole element is shared (and hence there is
    // an overlap of two points both of which need to have "radius" points on
    // either side.
    size_t inset=2*radius+2;

    // how wide is the x-dimension between the two insets
    size_t xmidlen=ext[0]-2*inset;
    size_t ymidlen=ext[1]-2*inset;
    size_t zmidlen=ext[2]-2*inset;

    Block block(ext[0], ext[1], ext[2], inset, xmidlen, ymidlen, zmidlen, numvals);

    MPI_Request reqs[50]; // a non-tight upper bound on how many we need
    MPI_Status stats[50];
    short rused=0;

    messvec incoms;
    messvec outcoms;

    grid.generateInNeighbours(X, Y, Z ,incoms);
    grid.generateOutNeighbours(X, Y, Z, outcoms);

    block.copyAllToBuffer(src);

    int comserr=0;
    for (size_t i=0; i < incoms.size(); ++i) {
        message& m = incoms[i];
        comserr |= MPI_Irecv(block.getInBuffer(m.destbuffid),
                           block.getBuffSize(m.destbuffid), MPI_DOUBLE,
                           m.sourceID, m.tag, m_mpiInfo->comm, reqs+(rused++));
        block.setUsed(m.destbuffid);
    }

    for (size_t i=0; i<outcoms.size(); ++i) {
        message& m=outcoms[i];
        comserr |= MPI_Isend(block.getOutBuffer(m.srcbuffid),
                             block.getBuffSize(m.srcbuffid), MPI_DOUBLE,
                             m.destID, m.tag, m_mpiInfo->comm, reqs+(rused++));
    }

    if (!comserr) {
        comserr=MPI_Waitall(rused, reqs, stats);
    }

    if (comserr) {
        // Yes this is throwing an exception as a result of an MPI error.
        // and no we don't inform the other ranks that we are doing this.
        // however, we have no reason to believe coms work at this point anyway
        throw RipleyException("Error in coms for randomFill");
    }

    block.copyUsedFromBuffer(src);
#endif // ESYS_MPI

    // the truth of either should imply the truth of the other but let's be safe
    if (radius==0 || numvals>1) {
        escript::FunctionSpace fs(getPtr(), getContinuousFunctionCode());
        escript::Data resdat(0, shape, fs , true);
        // don't need to check for exwrite because we just made it
        escript::DataTypes::RealVectorType& dv=resdat.getExpandedVectorReference();

        // now we need to copy values over
        for (size_t z=0; z < internal[2]; ++z) {
            for (size_t y=0; y < internal[1]; ++y) {
                for (size_t x=0; x < internal[0]; ++x) {
                    for (unsigned int i=0; i < numvals; ++i) {
                        dv[i+(x+y*(internal[0])+z*internal[0]*internal[1])*numvals]=src[i+(x+y*ext[0]+z*ext[0]*ext[1])*numvals];
                    }
                }
            }
        }
        delete[] src;
        return resdat;
    } else { // filter enabled
        escript::FunctionSpace fs(getPtr(), getContinuousFunctionCode());
        escript::Data resdat(0, escript::DataTypes::scalarShape, fs , true);
        // don't need to check for exwrite because we just made it
        escript::DataTypes::RealVectorType& dv=resdat.getExpandedVectorReference();
        double* convolution=get3DGauss(radius, sigma);

        for (size_t z=0;z<(internal[2]);++z) {
            for (size_t y=0;y<(internal[1]);++y) {
                for (size_t x=0;x<(internal[0]);++x) {
                    dv[x+y*(internal[0])+z*internal[0]*internal[1]]=Convolve3D(convolution, src, x+radius, y+radius, z+radius, radius, ext[0], ext[1]);

                }
            }
        }

        delete[] convolution;
        delete[] src;
        return resdat;
    }
}

dim_t Brick::findNode(const double *coords) const
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
    double z = coords[2] - m_origin[2];

    //check if the point is even inside the domain
    if (x < 0 || y < 0 || z < 0
            || x > m_length[0] || y > m_length[1] || z > m_length[2])
        return NOT_MINE;

    // distance in elements
    dim_t ex = (dim_t) floor(x / m_dx[0]);
    dim_t ey = (dim_t) floor(y / m_dx[1]);
    dim_t ez = (dim_t) floor(z / m_dx[2]);
    // set the min distance high enough to be outside the element plus a bit
    dim_t closest = NOT_MINE;
    double minDist = 1;
    for (int dim = 0; dim < m_numDim; dim++) {
        minDist += m_dx[dim]*m_dx[dim];
    }
    //find the closest node
    for (int dx = 0; dx < 1; dx++) {
        double xdist = x - (ex + dx)*m_dx[0];
        for (int dy = 0; dy < 1; dy++) {
            double ydist = y - (ey + dy)*m_dx[1];
            for (int dz = 0; dz < 1; dz++) {
                double zdist = z - (ez + dz)*m_dx[2];
                double total = xdist*xdist + ydist*ydist + zdist*zdist;
                if (total < minDist) {
                    closest = INDEX3(ex+dy-m_offset[0], ey+dy-m_offset[1],
                            ez+dz-m_offset[2], m_NE[0]+1, m_NE[1]+1);
                    minDist = total;
                }
            }
        }
    }
    if (closest == NOT_MINE) {
        throw RipleyException("Unable to map appropriate dirac point to a "
                         "node, implementation problem in Brick::findNode()");
    }
    return closest;
}

Assembler_ptr Brick::createAssembler(string type, const DataMap& constants) const
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
            return Assembler_ptr(new DefaultAssembler3D<cplx_t>(shared_from_this(), m_dx, m_NE, m_NN));
        } else {
            return Assembler_ptr(new DefaultAssembler3D<real_t>(shared_from_this(), m_dx, m_NE, m_NN));
        }
    } else if (type.compare("WaveAssembler") == 0) {
        return Assembler_ptr(new WaveAssembler3D(shared_from_this(), m_dx, m_NE, m_NN, constants));
    } else if (type.compare("LameAssembler") == 0) {
        return Assembler_ptr(new LameAssembler3D(shared_from_this(), m_dx, m_NE, m_NN));
    } else { //else ifs would go before this for other types
        throw RipleyException("Ripley::Brick does not support the requested "
                              "assembler");
    }
}

} // end of namespace ripley
