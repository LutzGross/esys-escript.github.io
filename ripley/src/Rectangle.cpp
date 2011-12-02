
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
        << functionSpaceTypeAsString(fsType);
    throw RipleyException(msg.str());
}

bool Rectangle::ownSample(int fsCode, index_t id) const
{
#ifdef ESYS_MPI
    if (fsCode == Nodes) {
        const index_t myFirst=m_nodeDistribution[m_mpiInfo->rank];
        const index_t myLast=m_nodeDistribution[m_mpiInfo->rank+1]-1;
        return (m_nodeId[id]>=myFirst && m_nodeId[id]<=myLast);
    } else
        throw RipleyException("ownSample() only implemented for Nodes");
#else
    return true;
#endif
}

void Rectangle::interpolateNodesOnElements(escript::Data& out, escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    /* GENERATOR SNIP_INTERPOLATE_ELEMENTS TOP */
    const double tmp0_2 = 0.62200846792814621559;
    const double tmp0_1 = 0.044658198738520451079;
    const double tmp0_0 = 0.16666666666666666667;
#pragma omp parallel for
    for (index_t k1=0; k1 < m_NE1; ++k1) {
        for (index_t k0=0; k0 < m_NE0; ++k0) {
            const register double* f_10 = in.getSampleDataRO(INDEX2(k0+1,k1, m_N0));
            const register double* f_11 = in.getSampleDataRO(INDEX2(k0+1,k1+1, m_N0));
            const register double* f_01 = in.getSampleDataRO(INDEX2(k0,k1+1, m_N0));
            const register double* f_00 = in.getSampleDataRO(INDEX2(k0,k1, m_N0));
            double* o = out.getSampleDataRW(INDEX2(k0,k1,m_NE0));
            for (index_t i=0; i < numComp; ++i) {
                o[INDEX2(i,numComp,0)] = f_00[i]*tmp0_2 + f_11[i]*tmp0_1 + tmp0_0*(f_01[i] + f_10[i]);
                o[INDEX2(i,numComp,1)] = f_01[i]*tmp0_1 + f_10[i]*tmp0_2 + tmp0_0*(f_00[i] + f_11[i]);
                o[INDEX2(i,numComp,2)] = f_01[i]*tmp0_2 + f_10[i]*tmp0_1 + tmp0_0*(f_00[i] + f_11[i]);
                o[INDEX2(i,numComp,3)] = f_00[i]*tmp0_1 + f_11[i]*tmp0_2 + tmp0_0*(f_01[i] + f_10[i]);
            } /* end of component loop i */
        } /* end of k0 loop */
    } /* end of k1 loop */
    /* GENERATOR SNIP_INTERPOLATE_ELEMENTS BOTTOM */
}

void Rectangle::interpolateNodesOnFaces(escript::Data& out, escript::Data& in) const
{
    throw RipleyException("interpolateNodesOnFaces() not implemented");
}

void Rectangle::addPDEToSystem(escript::AbstractSystemMatrix& mat,
        escript::Data& rhs, const escript::Data& A, const escript::Data& B,
        const escript::Data& C, const escript::Data& D,
        const escript::Data& X, const escript::Data& Y,
        const escript::Data& d, const escript::Data& y,
        const escript::Data& d_contact, const escript::Data& y_contact,
        const escript::Data& d_dirac, const escript::Data& y_dirac) const
{
    throw RipleyException("addPDEToSystem() not implemented");
}

Paso_SystemMatrixPattern* Rectangle::getPattern(bool reducedRowOrder,
                                                bool reducedColOrder) const
{
    if (reducedRowOrder || reducedColOrder)
        throw RipleyException("getPattern() not implemented for reduced order");

    // connector
    RankVector neighbour;
    IndexVector offsetInShared(1,0);
    IndexVector sendShared, recvShared;
    const IndexVector faces=getNumFacesPerBoundary();
    const index_t left = (m_offset0==0 ? 0 : 1);
    const index_t bottom = (m_offset1==0 ? 0 : 1);
    // corner node from bottom-left
    if (faces[0] == 0 && faces[2] == 0) {
        neighbour.push_back(m_mpiInfo->rank-m_NX-1);
        offsetInShared.push_back(offsetInShared.back()+1);
        sendShared.push_back(m_nodeId[m_N0+left]);
        recvShared.push_back(m_nodeId[0]);
    }
    // bottom edge
    if (faces[2] == 0) {
        neighbour.push_back(m_mpiInfo->rank-m_NX);
        offsetInShared.push_back(offsetInShared.back()+m_N0-left);
        for (dim_t i=left; i<m_N0; i++) {
            // easy case, we know the neighbour id's
            sendShared.push_back(m_nodeId[m_N0+i]);
            recvShared.push_back(m_nodeId[i]);
        }
    }
    // corner node from bottom-right
    if (faces[1] == 0 && faces[2] == 0) {
        neighbour.push_back(m_mpiInfo->rank-m_NX+1);
        const index_t N0=(neighbour.back()%m_NX == 0 ? m_N0 : m_N0-1);
        const index_t N1=(neighbour.back()/m_NX == 0 ? m_N1 : m_N1-1);
        const index_t first=m_nodeDistribution[neighbour.back()];
        offsetInShared.push_back(offsetInShared.back()+1);
        sendShared.push_back(m_nodeId[(bottom+1)*m_N0-1]);
        recvShared.push_back(first+N0*(N1-1));
    }
    // left edge
    if (faces[0] == 0) {
        neighbour.push_back(m_mpiInfo->rank-1);
        offsetInShared.push_back(offsetInShared.back()+m_N1-bottom);
        for (dim_t i=bottom; i<m_N1; i++) {
            // easy case, we know the neighbour id's
            sendShared.push_back(m_nodeId[i*m_N0+1]);
            recvShared.push_back(m_nodeId[i*m_N0]);
        }
    }
    // right edge
    if (faces[1] == 0) {
        neighbour.push_back(m_mpiInfo->rank+1);
        const index_t rightN0=(neighbour.back()%m_NX == 0 ? m_N0 : m_N0-1);
        const index_t first=m_nodeDistribution[neighbour.back()];
        offsetInShared.push_back(offsetInShared.back()+m_N1-bottom);
        for (dim_t i=bottom; i<m_N1; i++) {
            sendShared.push_back(m_nodeId[(i+1)*m_N0-1]);
            recvShared.push_back(first+rightN0*(i-bottom));
        }
    }
    // corner node from top-left
    if (faces[0] == 0 && faces[3] == 0) {
        neighbour.push_back(m_mpiInfo->rank+m_NX-1);
        const index_t N0=(neighbour.back()%m_NX == 0 ? m_N0 : m_N0-1);
        const index_t first=m_nodeDistribution[neighbour.back()];
        offsetInShared.push_back(offsetInShared.back()+1);
        sendShared.push_back(m_nodeId[m_N0*(m_N1-1)+left]);
        recvShared.push_back(first+N0-1);
    }
    // top edge
    if (faces[3] == 0) {
        neighbour.push_back(m_mpiInfo->rank+m_NX);
        const index_t first=m_nodeDistribution[neighbour.back()];
        offsetInShared.push_back(offsetInShared.back()+m_N0-left);
        for (dim_t i=left; i<m_N0; i++) {
            sendShared.push_back(m_nodeId[m_N0*(m_N1-1)+i]);
            recvShared.push_back(first+i-left);
        }
    }
    // corner node from top-right
    if (faces[1] == 0 && faces[3] == 0) {
        neighbour.push_back(m_mpiInfo->rank+m_NX+1);
        const index_t first=m_nodeDistribution[neighbour.back()];
        offsetInShared.push_back(offsetInShared.back()+1);
        sendShared.push_back(m_nodeId[m_N0*m_N1-1]);
        recvShared.push_back(first);
    }
    const int numDOF=m_nodeDistribution[m_mpiInfo->rank+1]-m_nodeDistribution[m_mpiInfo->rank];
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

    Paso_SharedComponents *snd_shcomp = Paso_SharedComponents_alloc(
            numDOF, neighbour.size(), &neighbour[0], &sendShared[0],
            &offsetInShared[0], 1, 0, m_mpiInfo);
    Paso_SharedComponents *rcv_shcomp = Paso_SharedComponents_alloc(
            numDOF, neighbour.size(), &neighbour[0], &recvShared[0],
            &offsetInShared[0], 1, 0, m_mpiInfo);
    Paso_Connector* connector = Paso_Connector_alloc(snd_shcomp, rcv_shcomp);
    Paso_SharedComponents_free(snd_shcomp);
    Paso_SharedComponents_free(rcv_shcomp);

    // create patterns
    dim_t M, N;
    IndexVector ptr(1,0);
    IndexVector index;

    // main pattern
    for (index_t i=0; i<numDOF; i++) {
        // always add the node itself
        index.push_back(i);
        int num=insertNeighbours(index, i);
        ptr.push_back(ptr.back()+num+1);
    }
    M=N=ptr.size()-1;
    // paso will manage the memory
    index_t* indexC = MEMALLOC(index.size(),index_t);
    index_t* ptrC = MEMALLOC(ptr.size(), index_t);
    copy(index.begin(), index.end(), indexC);
    copy(ptr.begin(), ptr.end(), ptrC);
    Paso_Pattern *mainPattern = Paso_Pattern_alloc(MATRIX_FORMAT_DEFAULT,
            M, N, ptrC, indexC);

    cout << "--- main_pattern ---" << endl;
    cout << "M=" << M << ", N=" << N << endl;
    for (size_t i=0; i<ptr.size(); i++) {
        cout << "ptr[" << i << "]=" << ptr[i] << endl;
    }
    for (size_t i=0; i<index.size(); i++) {
        cout << "index[" << i << "]=" << index[i] << endl;
    }

    ptr.clear();
    index.clear();

    // column & row couple patterns
    Paso_Pattern *colCouplePattern, *rowCouplePattern;
    generateCouplePatterns(&colCouplePattern, &rowCouplePattern);

    // allocate paso distribution
    Paso_Distribution* distribution = Paso_Distribution_alloc(m_mpiInfo,
            const_cast<index_t*>(&m_nodeDistribution[0]), 1, 0);

    Paso_SystemMatrixPattern* pattern = Paso_SystemMatrixPattern_alloc(
            MATRIX_FORMAT_DEFAULT, distribution, distribution,
            mainPattern, colCouplePattern, rowCouplePattern,
            connector, connector);
    Paso_Pattern_free(mainPattern);
    Paso_Pattern_free(colCouplePattern);
    Paso_Pattern_free(rowCouplePattern);
    Paso_Distribution_free(distribution);
    return pattern;
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
    const IndexVector faces = getNumFacesPerBoundary();
    dim_t n=0;
    for (size_t i=0; i<faces.size(); i++)
        n+=faces[i];
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
        const index_t N0=(neighbour%m_NX == 0 ? m_N0 : m_N0-1);
        const index_t N1=(neighbour/m_NX == 0 ? m_N1 : m_N1-1);
        m_nodeId[0]=m_nodeDistribution[neighbour]+N0*N1-1;
    }

    // the rest of the id's are contiguous
    const index_t firstId=m_nodeDistribution[m_mpiInfo->rank];
#pragma omp parallel for
    for (dim_t i1=bottom; i1<m_N1; i1++) {
        for (dim_t i0=left; i0<m_N0; i0++) {
            m_nodeId[i0+i1*m_N0] = firstId+i0-left+(i1-bottom)*(m_N0-left);
        }
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

//private
int Rectangle::insertNeighbours(IndexVector& index, index_t node) const
{
    const dim_t myN0 = (m_offset0==0 ? m_N0 : m_N0-1);
    const dim_t myN1 = (m_offset1==0 ? m_N1 : m_N1-1);
    const int x=node%myN0;
    const int y=node/myN0;
    int num=0;
    if (y>0) {
        if (x>0) {
            // bottom-left
            index.push_back(node-myN0-1);
            num++;
        }
        // bottom
        index.push_back(node-myN0);
        num++;
        if (x<myN0-1) {
            // bottom-right
            index.push_back(node-myN0+1);
            num++;
        }
    }
    if (x<myN0-1) {
        // right
        index.push_back(node+1);
        num++;
        if (y<myN1-1) {
            // top-right
            index.push_back(node+myN0+1);
            num++;
        }
    }
    if (y<myN1-1) {
        // top
        index.push_back(node+myN0);
        num++;
        if (x>0) {
            // top-left
            index.push_back(node+myN0-1);
            num++;
        }
    }
    if (x>0) {
        // left
        index.push_back(node-1);
        num++;
    }

    return num;
}

//private
void Rectangle::generateCouplePatterns(Paso_Pattern** colPattern, Paso_Pattern** rowPattern) const
{
    IndexVector ptr(1,0);
    IndexVector index;
    const dim_t myN0 = (m_offset0==0 ? m_N0 : m_N0-1);
    const dim_t myN1 = (m_offset1==0 ? m_N1 : m_N1-1);
    const IndexVector faces=getNumFacesPerBoundary();

    // bottom edge
    dim_t n=0;
    if (faces[0] == 0) {
        index.push_back(2*(myN0+myN1+1));
        index.push_back(2*(myN0+myN1+1)+1);
        n+=2;
        if (faces[2] == 0) {
            index.push_back(0);
            index.push_back(1);
            index.push_back(2);
            n+=3;
        }
    } else if (faces[2] == 0) {
        index.push_back(1);
        index.push_back(2);
        n+=2;
    }
    // n=neighbours of bottom-left corner node
    ptr.push_back(ptr.back()+n);
    n=0;
    if (faces[2] == 0) {
        for (dim_t i=1; i<myN0-1; i++) {
            index.push_back(i);
            index.push_back(i+1);
            index.push_back(i+2);
            ptr.push_back(ptr.back()+3);
        }
        index.push_back(myN0-1);
        index.push_back(myN0);
        n+=2;
        if (faces[1] == 0) {
            index.push_back(myN0+1);
            index.push_back(myN0+2);
            index.push_back(myN0+3);
            n+=3;
        }
    } else {
        for (dim_t i=1; i<myN0-1; i++) {
            ptr.push_back(ptr.back());
        }
        if (faces[1] == 0) {
            index.push_back(myN0+2);
            index.push_back(myN0+3);
            n+=2;
        }
    }
    // n=neighbours of bottom-right corner node
    ptr.push_back(ptr.back()+n);

    // 2nd row to 2nd last row
    for (dim_t i=1; i<myN1-1; i++) {
        // left edge
        if (faces[0] == 0) {
            index.push_back(2*(myN0+myN1+2)-i);
            index.push_back(2*(myN0+myN1+2)-i-1);
            index.push_back(2*(myN0+myN1+2)-i-2);
            ptr.push_back(ptr.back()+3);
        } else {
            ptr.push_back(ptr.back());
        }
        for (dim_t j=1; j<myN0-1; j++) {
            ptr.push_back(ptr.back());
        }
        // right edge
        if (faces[1] == 0) {
            index.push_back(myN0+i+1);
            index.push_back(myN0+i+2);
            index.push_back(myN0+i+3);
            ptr.push_back(ptr.back()+3);
        } else {
            ptr.push_back(ptr.back());
        }
    }

    // top edge
    n=0;
    if (faces[0] == 0) {
        index.push_back(2*myN0+myN1+5);
        index.push_back(2*myN0+myN1+4);
        n+=2;
        if (faces[3] == 0) {
            index.push_back(2*myN0+myN1+3);
            index.push_back(2*myN0+myN1+2);
            index.push_back(2*myN0+myN1+1);
            n+=3;
        }
    } else if (faces[3] == 0) {
        index.push_back(2*myN0+myN1+2);
        index.push_back(2*myN0+myN1+1);
        n+=2;
    }
    // n=neighbours of top-left corner node
    ptr.push_back(ptr.back()+n);
    n=0;
    if (faces[3] == 0) {
        for (dim_t i=1; i<myN0-1; i++) {
            index.push_back(2*myN0+myN1+i+1);
            index.push_back(2*myN0+myN1+i);
            index.push_back(2*myN0+myN1+i-1);
            ptr.push_back(ptr.back()+3);
        }
        index.push_back(myN0+myN1+4);
        index.push_back(myN0+myN1+3);
        n+=2;
        if (faces[1] == 0) {
            index.push_back(myN0+myN1+2);
            index.push_back(myN0+myN1+1);
            index.push_back(myN0+myN1);
            n+=3;
        }
    } else {
        for (dim_t i=1; i<myN0-1; i++) {
            ptr.push_back(ptr.back());
        }
        if (faces[1] == 0) {
            index.push_back(myN0+myN1+1);
            index.push_back(myN0+myN1);
            n+=2;
        }
    }
    // n=neighbours of top-right corner node
    ptr.push_back(ptr.back()+n);

    dim_t M=ptr.size()-1;
    map<index_t,index_t> idMap;
    dim_t N=0;
    for (IndexVector::iterator it=index.begin(); it!=index.end(); it++) {
        if (idMap.find(*it)==idMap.end()) {
            idMap[*it]=N;
            *it=N++;
        } else {
            *it=idMap[*it];
        }
    }

    cout << "--- colCouple_pattern ---" << endl;
    cout << "M=" << M << ", N=" << N << endl;
    for (size_t i=0; i<ptr.size(); i++) {
        cout << "ptr[" << i << "]=" << ptr[i] << endl;
    }
    for (size_t i=0; i<index.size(); i++) {
        cout << "index[" << i << "]=" << index[i] << endl;
    }

    // now build the row couple pattern
    IndexVector ptr2(1,0);
    IndexVector index2;
    for (dim_t id=0; id<N; id++) {
        n=0;
        for (dim_t i=0; i<M; i++) {
            for (dim_t j=ptr[i]; j<ptr[i+1]; j++) {
                if (index[j]==id) {
                    index2.push_back(i);
                    n++;
                    break;
                }
            }
        }
        ptr2.push_back(ptr2.back()+n);
    }

    cout << "--- rowCouple_pattern ---" << endl;
    cout << "M=" << N << ", N=" << M << endl;
    for (size_t i=0; i<ptr2.size(); i++) {
        cout << "ptr[" << i << "]=" << ptr2[i] << endl;
    }
    for (size_t i=0; i<index2.size(); i++) {
        cout << "index[" << i << "]=" << index2[i] << endl;
    }

    // paso will manage the memory
    index_t* indexC = MEMALLOC(index.size(), index_t);
    index_t* ptrC = MEMALLOC(ptr.size(), index_t);
    copy(index.begin(), index.end(), indexC);
    copy(ptr.begin(), ptr.end(), ptrC);
    *colPattern=Paso_Pattern_alloc(MATRIX_FORMAT_DEFAULT, M, N, ptrC, indexC);

    // paso will manage the memory
    indexC = MEMALLOC(index2.size(), index_t);
    ptrC = MEMALLOC(ptr2.size(), index_t);
    copy(index2.begin(), index2.end(), indexC);
    copy(ptr2.begin(), ptr2.end(), ptrC);
    *rowPattern=Paso_Pattern_alloc(MATRIX_FORMAT_DEFAULT, N, M, ptrC, indexC);
}

} // end of namespace ripley

