
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include <ripley/Rectangle.h>
extern "C" {
#include "paso/SystemMatrixPattern.h"
}

#if USE_SILO
#include <silo.h>
#ifdef ESYS_MPI
#include <pmpio.h>
#endif
#endif

#include <iomanip>

using namespace std;

namespace ripley {

Rectangle::Rectangle(int n0, int n1, double l0, double l1, int d0, int d1) :
    RipleyDomain(2),
    m_gNE0(n0),
    m_gNE1(n1),
    m_l0(l0),
    m_l1(l1),
    m_NX(d0),
    m_NY(d1)
{
    // ensure number of subdivisions is valid and nodes can be distributed
    // among number of ranks
    if (m_NX*m_NY != m_mpiInfo->size)
        throw RipleyException("Invalid number of spatial subdivisions");

    if (n0%m_NX > 0 || n1%m_NY > 0)
        throw RipleyException("Number of elements must be separable into number of ranks in each dimension");

    // local number of elements
    m_NE0 = n0/m_NX;
    m_NE1 = n1/m_NY;
    // local number of nodes (not necessarily owned)
    m_N0 = m_NE0+1;
    m_N1 = m_NE1+1;
    // bottom-left node is at (offset0,offset1) in global mesh
    m_offset0 = m_NE0*(m_mpiInfo->rank%m_NX);
    m_offset1 = m_NE1*(m_mpiInfo->rank/m_NX);
    populateSampleIds();
}

Rectangle::~Rectangle()
{
}

string Rectangle::getDescription() const
{
    return "ripley::Rectangle";
}

bool Rectangle::operator==(const AbstractDomain& other) const
{
    if (dynamic_cast<const Rectangle*>(&other))
        return this==&other;

    return false;
}

void Rectangle::dump(const string& fileName) const
{
#if USE_SILO
    string fn(fileName);
    if (fileName.length() < 6 || fileName.compare(fileName.length()-5, 5, ".silo") != 0) {
        fn+=".silo";
    }

    const int NUM_SILO_FILES = 1;
    const char* blockDirFmt = "/block%04d";
    int driver=DB_HDF5;    
    string siloPath;
    DBfile* dbfile = NULL;

#ifdef ESYS_MPI
    PMPIO_baton_t* baton = NULL;
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

    boost::scoped_ptr<double> x(new double[m_N0]);
    boost::scoped_ptr<double> y(new double[m_N1]);
    double* coords[2] = { x.get(), y.get() };
    pair<double,double> xdx = getFirstCoordAndSpacing(0);
    pair<double,double> ydy = getFirstCoordAndSpacing(1);
#pragma omp parallel
    {
#pragma omp for
        for (dim_t i0 = 0; i0 < m_N0; i0++) {
            coords[0][i0]=xdx.first+i0*xdx.second;
        }
#pragma omp for
        for (dim_t i1 = 0; i1 < m_N1; i1++) {
            coords[1][i1]=ydy.first+i1*ydy.second;
        }
    }
    IndexVector dims = getNumNodesPerDim();

    // write mesh
    DBPutQuadmesh(dbfile, "mesh", NULL, coords, &dims[0], 2, DB_DOUBLE,
            DB_COLLINEAR, NULL);

    // write node ids
    DBPutQuadvar1(dbfile, "nodeId", "mesh", (void*)&m_nodeId[0], &dims[0], 2,
            NULL, 0, DB_INT, DB_NODECENT, NULL);

    // write element ids
    dims = getNumElementsPerDim();
    DBPutQuadvar1(dbfile, "elementId", "mesh", (void*)&m_elementId[0],
            &dims[0], 2, NULL, 0, DB_INT, DB_ZONECENT, NULL);

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
    throw RipleyException("dump(): no Silo support");
#endif
}

const int* Rectangle::borrowSampleReferenceIDs(int fsType) const
{
    switch (fsType) {
        case Nodes:
            return &m_nodeId[0];
        case Elements:
            return &m_elementId[0];
        case FaceElements:
            return &m_faceId[0];
        default:
            break;
    }

    stringstream msg;
    msg << "borrowSampleReferenceIDs() not implemented for function space type "
        << fsType;
    throw RipleyException(msg.str());
}

bool Rectangle::ownSample(int fsCode, index_t id) const
{
#ifdef ESYS_MPI
    if (fsCode == Nodes) {
        const index_t myFirst=getNumNodes()*m_mpiInfo->rank;
        const index_t myLast=getNumNodes()*(m_mpiInfo->rank+1)-1;
        return (m_nodeId[id]>=myFirst && m_nodeId[id]<=myLast);
    } else
        throw RipleyException("ownSample() only implemented for Nodes");
#else
    return true;
#endif
}

void Rectangle::interpolateOnDomain(escript::Data& target,
                                    const escript::Data& in) const
{
    const Rectangle& inDomain=dynamic_cast<const Rectangle&>(*(in.getFunctionSpace().getDomain()));
    const Rectangle& targetDomain=dynamic_cast<const Rectangle&>(*(target.getFunctionSpace().getDomain()));
    if (inDomain != *this)
        throw RipleyException("Illegal domain of interpolant");
    if (targetDomain != *this)
        throw RipleyException("Illegal domain of interpolation target");

    stringstream msg;
    msg << "interpolateOnDomain() not implemented for function space "
        << in.getFunctionSpace().getTypeCode() << " -> "
        << target.getFunctionSpace().getTypeCode();

    switch (in.getFunctionSpace().getTypeCode()) {
        case Nodes:
        case DegreesOfFreedom:
            switch (target.getFunctionSpace().getTypeCode()) {
                case Elements:
                {
                    const double tmp0_2 = 0.62200846792814621559;
                    const double tmp0_1 = 0.044658198738520451079;
                    const double tmp0_0 = 0.16666666666666666667;
                    const dim_t numComp = in.getDataPointSize();
                    escript::Data* inn=const_cast<escript::Data*>(&in);
#pragma omp parallel for
                    for (index_t k1=0; k1 < m_NE1; ++k1) {
                        for (index_t k0=0; k0 < m_NE0; ++k0) {
                            const register double* f_10 = inn->getSampleDataRO(INDEX2(k0+1,k1, m_N0));
                            const register double* f_11 = inn->getSampleDataRO(INDEX2(k0+1,k1+1, m_N0));
                            const register double* f_01 = inn->getSampleDataRO(INDEX2(k0,k1+1, m_N0));
                            const register double* f_00 = inn->getSampleDataRO(INDEX2(k0,k1, m_N0));
                            double* o = target.getSampleDataRW(INDEX2(k0,k1,m_NE0));
                            for (index_t i=0; i < numComp; ++i) {
                                o[INDEX2(i,numComp,0)] = f_00[i]*tmp0_2 + f_11[i]*tmp0_1 + tmp0_0*(f_01[i] + f_10[i]);
                                o[INDEX2(i,numComp,1)] = f_01[i]*tmp0_1 + f_10[i]*tmp0_2 + tmp0_0*(f_00[i] + f_11[i]);
                                o[INDEX2(i,numComp,2)] = f_01[i]*tmp0_2 + f_10[i]*tmp0_1 + tmp0_0*(f_00[i] + f_11[i]);
                                o[INDEX2(i,numComp,3)] = f_00[i]*tmp0_1 + f_11[i]*tmp0_2 + tmp0_0*(f_01[i] + f_10[i]);
                            } // close component loop i
                        } // close k0 loop
                    } // close k1 loop
                    break;
                }

                case DegreesOfFreedom:
                    target=in;
                    break;
                    
                default:
                    throw RipleyException(msg.str());
            }
        default:
            throw RipleyException(msg.str());
    }
}

Paso_SystemMatrixPattern* Rectangle::getPattern(bool reducedRowOrder,
                                                bool reducedColOrder) const
{
    if (reducedRowOrder || reducedColOrder)
        throw RipleyException("getPattern() not implemented for reduced order");

/*
    // create distribution
    IndexVector dist;
    for (index_t i=0; i<m_mpiInfo->size+1; i++)
        dist.push_back(i*getNumNodes());
    Paso_Distribution* distribution = Paso_Distribution_alloc(m_mpiInfo,
            &dist[0], 1, 0);

    // connectors
    dim_t numNeighbours = 0;
    RankVector neighbour(numNeighbours);
    IndexVector offsetInShared(numNeighbours+1);
    IndexVector shared(offsetInShared[numNeighbours]);

    Paso_SharedComponents *snd_shcomp = Paso_SharedComponents_alloc(
            getNumNodes(), numNeighbours, &neighbour[0], &shared[0],
            &offsetInShared[0], 1, 0, m_mpiInfo);
    Paso_SharedComponents *rcv_shcomp = Paso_SharedComponents_alloc(
            getNumNodes(), numNeighbours, &neighbour[0], &shared[0],
            &offsetInShared[0], 1, 0, m_mpiInfo);
    Paso_Connector* connector = Paso_Connector_alloc(snd_shcomp, rcv_shcomp);

    // create patterns
    dim_t M=0, N=0;
    int* ptr=NULL;
    index_t* index=NULL;
    Paso_Pattern *mainPattern = Paso_Pattern_alloc(MATRIX_FORMAT_DEFAULT,
            M, N, ptr, index);
    Paso_Pattern *colCouplePattern = Paso_Pattern_alloc(MATRIX_FORMAT_DEFAULT,
            M, N, ptr, index);
    Paso_Pattern *rowCouplePattern = Paso_Pattern_alloc(MATRIX_FORMAT_DEFAULT,
            M, N, ptr, index);

    Paso_SystemMatrixPattern* pattern = Paso_SystemMatrixPattern_alloc(
            MATRIX_FORMAT_DEFAULT, distribution, distribution,
            mainPattern, colCouplePattern, rowCouplePattern,
            connector, connector);
    Paso_Pattern_free(mainPattern);
    Paso_Pattern_free(colCouplePattern);
    Paso_Pattern_free(rowCouplePattern);
    Paso_Distribution_free(distribution);
    Paso_SharedComponents_free(snd_shcomp);
    Paso_SharedComponents_free(rcv_shcomp);
    return pattern;
*/
    throw RipleyException("getPattern() not implemented");
}

void Rectangle::Print_Mesh_Info(const bool full) const
{
    RipleyDomain::Print_Mesh_Info(full);
    if (full) {
        cout << "     Id  Coordinates" << endl;
        cout.precision(15);
        cout.setf(ios::scientific, ios::floatfield);
        pair<double,double> xdx = getFirstCoordAndSpacing(0);
        pair<double,double> ydy = getFirstCoordAndSpacing(1);
        for (index_t i=0; i < getNumNodes(); i++) {
            cout << "  " << setw(5) << m_nodeId[i]
                << "  " << xdx.first+(i%m_N0)*xdx.second
                << "  " << ydy.first+(i/m_N0)*ydy.second << endl;
        }
    }
}

IndexVector Rectangle::getNumNodesPerDim() const
{
    IndexVector ret;
    ret.push_back(m_N0);
    ret.push_back(m_N1);
    return ret;
}

IndexVector Rectangle::getNumElementsPerDim() const
{
    IndexVector ret;
    ret.push_back(m_NE0);
    ret.push_back(m_NE1);
    return ret;
}

IndexVector Rectangle::getNumFacesPerBoundary() const
{
    IndexVector ret(4, 0);
    //left
    if (m_offset0==0)
        ret[0]=m_NE1;
    //right
    if (m_mpiInfo->rank%m_NX==m_NX-1)
        ret[1]=m_NE1;
    //bottom
    if (m_offset1==0)
        ret[2]=m_NE0;
    //top
    if (m_mpiInfo->rank/m_NX==m_NY-1)
        ret[3]=m_NE0;
    return ret;
}

pair<double,double> Rectangle::getFirstCoordAndSpacing(dim_t dim) const
{
    if (dim==0) {
        return pair<double,double>((m_l0*m_offset0)/m_gNE0, m_l0/m_gNE0);
    } else if (dim==1) {
        return pair<double,double>((m_l1*m_offset1)/m_gNE1, m_l1/m_gNE1);
    }
    throw RipleyException("getFirstCoordAndSpacing(): invalid argument");
}

//protected
dim_t Rectangle::getNumFaceElements() const
{
    dim_t n=0;
    //left
    if (m_offset0==0)
        n+=m_NE1;
    //right
    if (m_mpiInfo->rank%m_NX==m_NX-1)
        n+=m_NE1;
    //bottom
    if (m_offset1==0)
        n+=m_NE0;
    //top
    if (m_mpiInfo->rank/m_NX==m_NY-1)
        n+=m_NE0;

    return n;
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

    pair<double,double> xdx = getFirstCoordAndSpacing(0);
    pair<double,double> ydy = getFirstCoordAndSpacing(1);
    arg.requireWrite();
#pragma omp parallel for
    for (dim_t i1 = 0; i1 < m_N1; i1++) {
        for (dim_t i0 = 0; i0 < m_N0; i0++) {
            double* point = arg.getSampleDataRW(i0+m_N0*i1);
            point[0] = xdx.first+i0*xdx.second;
            point[1] = ydy.first+i1*ydy.second;
        }
    }
}

//private
void Rectangle::populateSampleIds()
{
    // identifiers are ordered from left to right, bottom to top on each rank,
    // except for the shared nodes which are owned by the rank below / to the
    // left of the current rank

    // build node distribution vector first.
    // m_nodeDistribution[i] is the first node id on rank i, that is
    // rank i owns m_nodeDistribution[i+1]-nodeDistribution[i] nodes
    m_nodeDistribution.assign(m_mpiInfo->size+1, 0);
    m_nodeDistribution[1]=getNumNodes();
    for (dim_t k=1; k<m_mpiInfo->size-1; k++) {
        const index_t x=k%m_NX;
        const index_t y=k/m_NX;
        index_t numNodes=getNumNodes();
        if (x>0)
            numNodes-=m_N1;
        if (y>0)
            numNodes-=m_N0;
        if (x>0 && y>0)
            numNodes++; // subtracted corner twice -> fix that
        m_nodeDistribution[k+1]=m_nodeDistribution[k]+numNodes;
    }
    m_nodeDistribution[m_mpiInfo->size]=getNumDataPointsGlobal();

    m_nodeId.resize(getNumNodes());

    // the bottom row and left column are not owned by this rank so the
    // identifiers need to be computed accordingly
    const index_t left = (m_offset0==0 ? 0 : 1);
    const index_t bottom = (m_offset1==0 ? 0 : 1);
    if (left>0) {
        const int neighbour=m_mpiInfo->rank-1;
        const index_t leftN0=(neighbour%m_NX == 0 ? m_N0 : m_N0-1);
#pragma omp parallel for
        for (dim_t i1=bottom; i1<m_N1; i1++) {
            m_nodeId[i1*m_N0]=m_nodeDistribution[neighbour]
                + (i1-bottom+1)*leftN0-1;
        }
    }
    if (bottom>0) {
        const int neighbour=m_mpiInfo->rank-m_NX;
        const index_t bottomN0=(neighbour%m_NX == 0 ? m_N0 : m_N0-1);
        const index_t bottomN1=(neighbour/m_NX == 0 ? m_N1 : m_N1-1);
#pragma omp parallel for
        for (dim_t i0=left; i0<m_N0; i0++) {
            m_nodeId[i0]=m_nodeDistribution[neighbour]
                + (bottomN1-1)*bottomN0 + i0 - left;
        }
    }
    if (left>0 && bottom>0) {
        const int neighbour=m_mpiInfo->rank-m_NX-1;
        const index_t bottomN0=(neighbour%m_NX == 0 ? m_N0 : m_N0-1);
        const index_t bottomN1=(neighbour/m_NX == 0 ? m_N1 : m_N1-1);
        m_nodeId[0]=m_nodeDistribution[neighbour]+bottomN1*bottomN0-1;
    }

    // the rest of the id's are contiguous
    const index_t firstId=m_nodeDistribution[m_mpiInfo->rank];
#pragma omp parallel for
    for (dim_t i1=bottom; i1<m_N1; i1++) {
        for (dim_t i0=left; i0<m_N0; i0++) {
            m_nodeId[i0+i1*m_N0] = firstId+i0-left+(i1-bottom)*(m_N0-left);
        }
    }

    for (dim_t k=0; k<m_nodeDistribution.size(); k++) {
        cout << m_nodeDistribution[k] << endl;
    }
    cout<<endl;
    for (dim_t k=0; k<getNumNodes(); k++) {
        cout << m_nodeId[k] << endl;
    }

    // elements
    m_elementId.resize(getNumElements());
#pragma omp parallel for
    for (dim_t k=0; k<getNumElements(); k++) {
        m_elementId[k]=k;
    }

    // face elements
    m_faceId.resize(getNumFaceElements());
#pragma omp parallel for
    for (dim_t k=0; k<getNumFaceElements(); k++) {
        m_faceId[k]=k;
    }
}

} // end of namespace ripley

