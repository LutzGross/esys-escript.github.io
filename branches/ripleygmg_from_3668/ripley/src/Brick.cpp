
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

#include <ripley/Brick.h>
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

Brick::Brick(int n0, int n1, int n2, double l0, double l1, double l2, int d0,
             int d1, int d2) :
    RipleyDomain(3),
    m_gNE0(n0),
    m_gNE1(n1),
    m_gNE2(n2),
    m_l0(l0),
    m_l1(l1),
    m_l2(l2),
    m_NX(d0),
    m_NY(d1),
    m_NZ(d2)
{
    // ensure number of subdivisions is valid and nodes can be distributed
    // among number of ranks
    if (m_NX*m_NY*m_NZ != m_mpiInfo->size)
        throw RipleyException("Invalid number of spatial subdivisions");

    if (n0%m_NX > 0 || n1%m_NY > 0 || n2%m_NZ > 0)
        throw RipleyException("Number of elements must be separable into number of ranks in each dimension");

    // local number of elements
    m_NE0 = n0/m_NX;
    m_NE1 = n1/m_NY;
    m_NE2 = n2/m_NZ;
    // local number of nodes (not necessarily owned)
    m_N0 = m_NE0+1;
    m_N1 = m_NE1+1;
    m_N2 = m_NE2+1;
    // bottom-left-front node is at (offset0,offset1,offset2) in global mesh
    m_offset0 = m_NE0*(m_mpiInfo->rank%m_NX);
    m_offset1 = m_NE1*(m_mpiInfo->rank%(m_NX*m_NY)/m_NX);
    m_offset2 = m_NE2*(m_mpiInfo->rank/(m_NX*m_NY));
    populateSampleIds();
}


Brick::~Brick()
{
}

string Brick::getDescription() const
{
    return "ripley::Brick";
}

bool Brick::operator==(const AbstractDomain& other) const
{
    if (dynamic_cast<const Brick*>(&other))
        return this==&other;

    return false;
}

void Brick::dump(const string& fileName) const
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
    boost::scoped_ptr<double> z(new double[m_N2]);
    double* coords[3] = { x.get(), y.get(), z.get() };
    pair<double,double> xdx = getFirstCoordAndSpacing(0);
    pair<double,double> ydy = getFirstCoordAndSpacing(1);
    pair<double,double> zdz = getFirstCoordAndSpacing(2);
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
#pragma omp for
        for (dim_t i2 = 0; i2 < m_N2; i2++) {
            coords[2][i2]=zdz.first+i2*zdz.second;
        }
    }
    IndexVector dims = getNumNodesPerDim();
    DBPutQuadmesh(dbfile, "mesh", NULL, coords, &dims[0], 3, DB_DOUBLE,
            DB_COLLINEAR, NULL);

    DBPutQuadvar1(dbfile, "nodeId", "mesh", (void*)&m_nodeId[0], &dims[0], 3,
            NULL, 0, DB_INT, DB_NODECENT, NULL);

    // write element ids
    dims = getNumElementsPerDim();
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
    throw RipleyException("dump(): no Silo support");
#endif
}

const int* Brick::borrowSampleReferenceIDs(int fsType) const
{
    switch (fsType) {
        case Nodes:
        case DegreesOfFreedom: //FIXME
        case ReducedDegreesOfFreedom: //FIXME
            return &m_nodeId[0];
        case Elements:
        case ReducedElements:
            return &m_elementId[0];
        case ReducedFaceElements:
            return &m_faceId[0];
        default:
            break;
    }

    stringstream msg;
    msg << "borrowSampleReferenceIDs() not implemented for function space type "
        << fsType;
    throw RipleyException(msg.str());
}

bool Brick::ownSample(int fsCode, index_t id) const
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

void Brick::setToGradient(escript::Data& out, const escript::Data& cIn) const
{
    escript::Data& in = *const_cast<escript::Data*>(&cIn);
    const dim_t numComp = in.getDataPointSize();
    const double h0 = m_l0/m_gNE0;
    const double h1 = m_l1/m_gNE1;
    const double h2 = m_l1/m_gNE2;
    const double C0 = .044658198738520451079;
    const double C1 = .16666666666666666667;
    const double C2 = .21132486540518711775;
    const double C3 = .25;
    const double C4 = .5;
    const double C5 = .62200846792814621559;
    const double C6 = .78867513459481288225;

    if (out.getFunctionSpace().getTypeCode() == Elements) {
        /*** GENERATOR SNIP_GRAD_ELEMENTS TOP */
#pragma omp parallel for
        for (index_t k2=0; k2 < m_NE2; ++k2) {
            for (index_t k1=0; k1 < m_NE1; ++k1) {
                for (index_t k0=0; k0 < m_NE0; ++k0) {
                    const register double* f_000 = in.getSampleDataRO(INDEX3(k0,k1,k2, m_N0,m_N1));
                    const register double* f_001 = in.getSampleDataRO(INDEX3(k0,k1,k2+1, m_N0,m_N1));
                    const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,k1,k2+1, m_N0,m_N1));
                    const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,k1+1,k2+1, m_N0,m_N1));
                    const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,k1+1,k2, m_N0,m_N1));
                    const register double* f_011 = in.getSampleDataRO(INDEX3(k0,k1+1,k2+1, m_N0,m_N1));
                    const register double* f_010 = in.getSampleDataRO(INDEX3(k0,k1+1,k2, m_N0,m_N1));
                    const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,k1,k2, m_N0,m_N1));
                    double* o = out.getSampleDataRW(INDEX3(k0,k1,k2,m_NE0,m_NE1));
                    for (index_t i=0; i < numComp; ++i) {
                        const double V0=((f_100[i]-f_000[i])*C5 + (f_111[i]-f_011[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / h0;
                        const double V1=((f_110[i]-f_010[i])*C5 + (f_101[i]-f_001[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / h0;
                        const double V2=((f_101[i]-f_001[i])*C5 + (f_110[i]-f_010[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / h0;
                        const double V3=((f_111[i]-f_011[i])*C5 + (f_100[i]-f_000[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / h0;
                        const double V4=((f_010[i]-f_000[i])*C5 + (f_111[i]-f_101[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / h1;
                        const double V5=((f_110[i]-f_100[i])*C5 + (f_011[i]-f_001[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / h1;
                        const double V6=((f_011[i]-f_001[i])*C5 + (f_110[i]-f_100[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / h1;
                        const double V7=((f_111[i]-f_101[i])*C5 + (f_010[i]-f_000[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / h1;
                        const double V8=((f_001[i]-f_000[i])*C5 + (f_111[i]-f_110[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / h2;
                        const double V9=((f_101[i]-f_100[i])*C5 + (f_011[i]-f_010[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / h2;
                        const double V10=((f_011[i]-f_010[i])*C5 + (f_101[i]-f_100[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / h2;
                        const double V11=((f_111[i]-f_110[i])*C5 + (f_001[i]-f_000[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / h2;
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
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of k1 loop */
        } /* end of k2 loop */
        /* GENERATOR SNIP_GRAD_ELEMENTS BOTTOM */
    } else if (out.getFunctionSpace().getTypeCode() == ReducedElements) {
        /*** GENERATOR SNIP_GRAD_REDUCED_ELEMENTS TOP */
#pragma omp parallel for
        for (index_t k2=0; k2 < m_NE2; ++k2) {
            for (index_t k1=0; k1 < m_NE1; ++k1) {
                for (index_t k0=0; k0 < m_NE0; ++k0) {
                    const register double* f_000 = in.getSampleDataRO(INDEX3(k0,k1,k2, m_N0,m_N1));
                    const register double* f_001 = in.getSampleDataRO(INDEX3(k0,k1,k2+1, m_N0,m_N1));
                    const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,k1,k2+1, m_N0,m_N1));
                    const register double* f_011 = in.getSampleDataRO(INDEX3(k0,k1+1,k2+1, m_N0,m_N1));
                    const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,k1,k2, m_N0,m_N1));
                    const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,k1+1,k2, m_N0,m_N1));
                    const register double* f_010 = in.getSampleDataRO(INDEX3(k0,k1+1,k2, m_N0,m_N1));
                    const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,k1+1,k2+1, m_N0,m_N1));
                    double* o = out.getSampleDataRW(INDEX3(k0,k1,k2,m_NE0,m_NE1));
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_101[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_010[i]-f_011[i])*C3 / h0;
                        o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_011[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_100[i]-f_101[i])*C3 / h1;
                        o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_011[i]+f_101[i]+f_111[i]-f_000[i]-f_010[i]-f_100[i]-f_110[i])*C3 / h2;
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of k1 loop */
        } /* end of k2 loop */
        /* GENERATOR SNIP_GRAD_REDUCED_ELEMENTS BOTTOM */
    } else if (out.getFunctionSpace().getTypeCode() == FaceElements) {
        /*** GENERATOR SNIP_GRAD_FACES TOP */
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k1=0; k1 < m_NE1; ++k1) {
                        const register double* f_000 = in.getSampleDataRO(INDEX3(0,k1,k2, m_N0,m_N1));
                        const register double* f_001 = in.getSampleDataRO(INDEX3(0,k1,k2+1, m_N0,m_N1));
                        const register double* f_101 = in.getSampleDataRO(INDEX3(1,k1,k2+1, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(1,k1+1,k2+1, m_N0,m_N1));
                        const register double* f_110 = in.getSampleDataRO(INDEX3(1,k1+1,k2, m_N0,m_N1));
                        const register double* f_011 = in.getSampleDataRO(INDEX3(0,k1+1,k2+1, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(0,k1+1,k2, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(1,k1,k2, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            const double V0=((f_010[i]-f_000[i])*C6 + (f_011[i]-f_001[i])*C2) / h1;
                            const double V1=((f_010[i]-f_000[i])*C2 + (f_011[i]-f_001[i])*C6) / h1;
                            const double V2=((f_001[i]-f_000[i])*C6 + (f_010[i]-f_011[i])*C2) / h2;
                            const double V3=((f_001[i]-f_000[i])*C2 + (f_011[i]-f_010[i])*C6) / h2;
                            o[INDEX3(i,0,0,numComp,3)] = ((f_100[i]-f_000[i])*C5 + (f_111[i]-f_011[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / h0;
                            o[INDEX3(i,1,0,numComp,3)] = V0;
                            o[INDEX3(i,2,0,numComp,3)] = V2;
                            o[INDEX3(i,0,1,numComp,3)] = ((f_110[i]-f_010[i])*C5 + (f_101[i]-f_001[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / h0;
                            o[INDEX3(i,1,1,numComp,3)] = V0;
                            o[INDEX3(i,2,1,numComp,3)] = V3;
                            o[INDEX3(i,0,2,numComp,3)] = ((f_101[i]-f_001[i])*C5 + (f_110[i]-f_010[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / h0;
                            o[INDEX3(i,1,2,numComp,3)] = V1;
                            o[INDEX3(i,2,2,numComp,3)] = V2;
                            o[INDEX3(i,0,3,numComp,3)] = ((f_111[i]-f_011[i])*C5 + (f_100[i]-f_000[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / h0;
                            o[INDEX3(i,1,3,numComp,3)] = V1;
                            o[INDEX3(i,2,3,numComp,3)] = V3;
                        } /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            } /* end of face 0 */
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k1=0; k1 < m_NE1; ++k1) {
                        const register double* f_000 = in.getSampleDataRO(INDEX3(m_N0-2,k1,k2, m_N0,m_N1));
                        const register double* f_001 = in.getSampleDataRO(INDEX3(m_N0-2,k1,k2+1, m_N0,m_N1));
                        const register double* f_101 = in.getSampleDataRO(INDEX3(m_N0-1,k1,k2+1, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(m_N0-1,k1+1,k2+1, m_N0,m_N1));
                        const register double* f_110 = in.getSampleDataRO(INDEX3(m_N0-1,k1+1,k2, m_N0,m_N1));
                        const register double* f_011 = in.getSampleDataRO(INDEX3(m_N0-2,k1+1,k2+1, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(m_N0-2,k1+1,k2, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(m_N0-1,k1,k2, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            const double V0=((f_110[i]-f_100[i])*C6 + (f_111[i]-f_101[i])*C2) / h1;
                            const double V1=((f_110[i]-f_100[i])*C2 + (f_111[i]-f_101[i])*C6) / h1;
                            const double V2=((f_101[i]-f_100[i])*C6 + (f_111[i]-f_110[i])*C2) / h2;
                            const double V3=((f_101[i]-f_100[i])*C2 + (f_111[i]-f_110[i])*C6) / h2;
                            o[INDEX3(i,0,0,numComp,3)] = ((f_100[i]-f_000[i])*C5 + (f_111[i]-f_011[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / h0;
                            o[INDEX3(i,1,0,numComp,3)] = V0;
                            o[INDEX3(i,2,0,numComp,3)] = V2;
                            o[INDEX3(i,0,1,numComp,3)] = ((f_110[i]-f_010[i])*C5 + (f_101[i]-f_001[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / h0;
                            o[INDEX3(i,1,1,numComp,3)] = V0;
                            o[INDEX3(i,2,1,numComp,3)] = V3;
                            o[INDEX3(i,0,2,numComp,3)] = ((f_101[i]-f_001[i])*C5 + (f_110[i]-f_010[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / h0;
                            o[INDEX3(i,1,2,numComp,3)] = V1;
                            o[INDEX3(i,2,2,numComp,3)] = V2;
                            o[INDEX3(i,0,3,numComp,3)] = ((f_111[i]-f_011[i])*C5 + (f_100[i]-f_000[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / h0;
                            o[INDEX3(i,1,3,numComp,3)] = V1;
                            o[INDEX3(i,2,3,numComp,3)] = V3;
                        } /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            } /* end of face 1 */
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const register double* f_000 = in.getSampleDataRO(INDEX3(k0,0,k2, m_N0,m_N1));
                        const register double* f_001 = in.getSampleDataRO(INDEX3(k0,0,k2+1, m_N0,m_N1));
                        const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,0,k2+1, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,0,k2, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,1,k2+1, m_N0,m_N1));
                        const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,1,k2, m_N0,m_N1));
                        const register double* f_011 = in.getSampleDataRO(INDEX3(k0,1,k2+1, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(k0,1,k2, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            const double V0=((f_100[i]-f_000[i])*C6 + (f_101[i]-f_001[i])*C2) / h0;
                            const double V1=((f_001[i]-f_000[i])*C6 + (f_101[i]-f_100[i])*C2) / h2;
                            const double V2=((f_001[i]-f_000[i])*C2 + (f_101[i]-f_100[i])*C6) / h2;
                            o[INDEX3(i,0,0,numComp,3)] = V0;
                            o[INDEX3(i,1,0,numComp,3)] = ((f_010[i]-f_000[i])*C5 + (f_111[i]-f_101[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / h1;
                            o[INDEX3(i,2,0,numComp,3)] = V1;
                            o[INDEX3(i,0,1,numComp,3)] = V0;
                            o[INDEX3(i,1,1,numComp,3)] = ((f_110[i]-f_100[i])*C5 + (f_011[i]-f_001[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / h1;
                            o[INDEX3(i,2,1,numComp,3)] = V2;
                            o[INDEX3(i,0,2,numComp,3)] = V0;
                            o[INDEX3(i,1,2,numComp,3)] = ((f_011[i]-f_001[i])*C5 + (f_110[i]-f_100[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / h1;
                            o[INDEX3(i,2,2,numComp,3)] = V1;
                            o[INDEX3(i,0,3,numComp,3)] = V0;
                            o[INDEX3(i,1,3,numComp,3)] = ((f_111[i]-f_101[i])*C5 + (f_010[i]-f_000[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / h1;
                            o[INDEX3(i,2,3,numComp,3)] = V2;
                        } /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k2 loop */
            } /* end of face 2 */
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const register double* f_011 = in.getSampleDataRO(INDEX3(k0,m_N1-1,k2+1, m_N0,m_N1));
                        const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,m_N1-1,k2, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(k0,m_N1-1,k2, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,m_N1-1,k2+1, m_N0,m_N1));
                        const register double* f_000 = in.getSampleDataRO(INDEX3(k0,m_N1-2,k2, m_N0,m_N1));
                        const register double* f_001 = in.getSampleDataRO(INDEX3(k0,m_N1-2,k2+1, m_N0,m_N1));
                        const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,m_N1-2,k2+1, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,m_N1-2,k2, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            const double V0=((f_110[i]-f_010[i])*C6 + (f_111[i]-f_011[i])*C2) / h0;
                            const double V1=((f_110[i]-f_010[i])*C2 + (f_111[i]-f_011[i])*C6) / h0;
                            const double V2=((f_011[i]-f_010[i])*C6 + (f_111[i]-f_110[i])*C2) / h2;
                            const double V3=((f_011[i]-f_010[i])*C2 + (f_111[i]-f_110[i])*C6) / h2;
                            o[INDEX3(i,0,0,numComp,3)] = V0;
                            o[INDEX3(i,1,0,numComp,3)] = ((f_010[i]-f_000[i])*C5 + (f_111[i]-f_101[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / h1;
                            o[INDEX3(i,2,0,numComp,3)] = V2;
                            o[INDEX3(i,0,1,numComp,3)] = V0;
                            o[INDEX3(i,1,1,numComp,3)] = ((f_110[i]-f_100[i])*C5 + (f_011[i]-f_001[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / h1;
                            o[INDEX3(i,2,1,numComp,3)] = V3;
                            o[INDEX3(i,0,2,numComp,3)] = V1;
                            o[INDEX3(i,1,2,numComp,3)] = ((f_011[i]-f_001[i])*C5 + (f_110[i]-f_100[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / h1;
                            o[INDEX3(i,2,2,numComp,3)] = V2;
                            o[INDEX3(i,0,3,numComp,3)] = V1;
                            o[INDEX3(i,1,3,numComp,3)] = ((f_111[i]-f_101[i])*C5 + (f_010[i]-f_000[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / h1;
                            o[INDEX3(i,2,3,numComp,3)] = V3;
                        } /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k2 loop */
            } /* end of face 3 */
            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const register double* f_000 = in.getSampleDataRO(INDEX3(k0,k1,0, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,k1,0, m_N0,m_N1));
                        const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,k1+1,0, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(k0,k1+1,0, m_N0,m_N1));
                        const register double* f_001 = in.getSampleDataRO(INDEX3(k0,k1,1, m_N0,m_N1));
                        const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,k1,1, m_N0,m_N1));
                        const register double* f_011 = in.getSampleDataRO(INDEX3(k0,k1+1,1, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,k1+1,1, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            const double V0=((f_100[i]-f_000[i])*C6 + (f_110[i]-f_010[i])*C2) / h0;
                            const double V1=((f_100[i]-f_000[i])*C2 + (f_110[i]-f_010[i])*C6) / h0;
                            const double V2=((f_010[i]-f_000[i])*C6 + (f_110[i]-f_100[i])*C2) / h1;
                            const double V3=((f_010[i]-f_000[i])*C2 + (f_110[i]-f_100[i])*C6) / h1;
                            o[INDEX3(i,0,0,numComp,3)] = V0;
                            o[INDEX3(i,1,0,numComp,3)] = V2;
                            o[INDEX3(i,2,0,numComp,3)] = ((f_001[i]-f_000[i])*C5 + (f_111[i]-f_110[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / h2;
                            o[INDEX3(i,0,1,numComp,3)] = V0;
                            o[INDEX3(i,1,1,numComp,3)] = V3;
                            o[INDEX3(i,2,1,numComp,3)] = ((f_101[i]-f_100[i])*C5 + (f_011[i]-f_010[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / h2;
                            o[INDEX3(i,0,2,numComp,3)] = V1;
                            o[INDEX3(i,1,2,numComp,3)] = V2;
                            o[INDEX3(i,2,2,numComp,3)] = ((f_011[i]-f_010[i])*C5 + (f_101[i]-f_100[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / h2;
                            o[INDEX3(i,0,3,numComp,3)] = V1;
                            o[INDEX3(i,1,3,numComp,3)] = V3;
                            o[INDEX3(i,2,3,numComp,3)] = ((f_111[i]-f_110[i])*C5 + (f_001[i]-f_000[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / h2;
                        } /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k1 loop */
            } /* end of face 4 */
            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const register double* f_001 = in.getSampleDataRO(INDEX3(k0,k1,m_N2-1, m_N0,m_N1));
                        const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,k1,m_N2-1, m_N0,m_N1));
                        const register double* f_011 = in.getSampleDataRO(INDEX3(k0,k1+1,m_N2-1, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,k1+1,m_N2-1, m_N0,m_N1));
                        const register double* f_000 = in.getSampleDataRO(INDEX3(k0,k1,m_N2-2, m_N0,m_N1));
                        const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,k1+1,m_N2-2, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(k0,k1+1,m_N2-2, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,k1,m_N2-2, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            const double V0=((f_101[i]-f_001[i])*C6 + (f_111[i]-f_011[i])*C2) / h0;
                            const double V1=((f_101[i]-f_001[i])*C2 + (f_111[i]-f_011[i])*C6) / h0;
                            const double V2=((f_011[i]-f_001[i])*C6 + (f_111[i]-f_101[i])*C2) / h1;
                            const double V3=((f_011[i]-f_001[i])*C2 + (f_111[i]-f_101[i])*C6) / h1;
                            o[INDEX3(i,0,0,numComp,3)] = V0;
                            o[INDEX3(i,1,0,numComp,3)] = V2;
                            o[INDEX3(i,2,0,numComp,3)] = ((f_001[i]-f_000[i])*C5 + (f_111[i]-f_110[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / h2;
                            o[INDEX3(i,0,1,numComp,3)] = V0;
                            o[INDEX3(i,1,1,numComp,3)] = V3;
                            o[INDEX3(i,2,1,numComp,3)] = ((f_011[i]-f_010[i])*C0 + (f_101[i]-f_100[i])*C5 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / h2;
                            o[INDEX3(i,0,2,numComp,3)] = V1;
                            o[INDEX3(i,1,2,numComp,3)] = V2;
                            o[INDEX3(i,2,2,numComp,3)] = ((f_011[i]-f_010[i])*C5 + (f_101[i]-f_100[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / h2;
                            o[INDEX3(i,0,3,numComp,3)] = V1;
                            o[INDEX3(i,1,3,numComp,3)] = V3;
                            o[INDEX3(i,2,3,numComp,3)] = ((f_001[i]-f_000[i])*C0 + (f_111[i]-f_110[i])*C5 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / h2;
                        } /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k1 loop */
            } /* end of face 5 */
        } // end of parallel section
        /* GENERATOR SNIP_GRAD_FACES BOTTOM */
    } else if (out.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        /*** GENERATOR SNIP_GRAD_REDUCED_FACES TOP */
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k1=0; k1 < m_NE1; ++k1) {
                        const register double* f_000 = in.getSampleDataRO(INDEX3(0,k1,k2, m_N0,m_N1));
                        const register double* f_001 = in.getSampleDataRO(INDEX3(0,k1,k2+1, m_N0,m_N1));
                        const register double* f_101 = in.getSampleDataRO(INDEX3(1,k1,k2+1, m_N0,m_N1));
                        const register double* f_011 = in.getSampleDataRO(INDEX3(0,k1+1,k2+1, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(1,k1,k2, m_N0,m_N1));
                        const register double* f_110 = in.getSampleDataRO(INDEX3(1,k1+1,k2, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(0,k1+1,k2, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(1,k1+1,k2+1, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_101[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_010[i]-f_011[i])*C3 / h0;
                            o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_011[i]-f_000[i]-f_001[i])*C4 / h1;
                            o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_011[i]-f_000[i]-f_010[i])*C4 / h2;
                        } /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            } /* end of face 0 */
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k1=0; k1 < m_NE1; ++k1) {
                        const register double* f_000 = in.getSampleDataRO(INDEX3(m_N0-2,k1,k2, m_N0,m_N1));
                        const register double* f_001 = in.getSampleDataRO(INDEX3(m_N0-2,k1,k2+1, m_N0,m_N1));
                        const register double* f_101 = in.getSampleDataRO(INDEX3(m_N0-1,k1,k2+1, m_N0,m_N1));
                        const register double* f_011 = in.getSampleDataRO(INDEX3(m_N0-2,k1+1,k2+1, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(m_N0-1,k1,k2, m_N0,m_N1));
                        const register double* f_110 = in.getSampleDataRO(INDEX3(m_N0-1,k1+1,k2, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(m_N0-2,k1+1,k2, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(m_N0-1,k1+1,k2+1, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_101[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_010[i]-f_011[i])*C3 / h0;
                            o[INDEX3(i,1,0,numComp,3)] = (f_110[i]+f_111[i]-f_100[i]-f_101[i])*C4 / h1;
                            o[INDEX3(i,2,0,numComp,3)] = (f_101[i]+f_111[i]-f_100[i]-f_110[i])*C4 / h2;
                        } /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            } /* end of face 1 */
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const register double* f_000 = in.getSampleDataRO(INDEX3(k0,0,k2, m_N0,m_N1));
                        const register double* f_001 = in.getSampleDataRO(INDEX3(k0,0,k2+1, m_N0,m_N1));
                        const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,0,k2+1, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,0,k2, m_N0,m_N1));
                        const register double* f_011 = in.getSampleDataRO(INDEX3(k0,1,k2+1, m_N0,m_N1));
                        const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,1,k2, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(k0,1,k2, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,1,k2+1, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_101[i]-f_000[i]-f_001[i])*C4 / h0;
                            o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_011[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_100[i]-f_101[i])*C3 / h1;
                            o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_101[i]-f_000[i]-f_100[i])*C4 / h2;
                        } /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k2 loop */
            } /* end of face 2 */
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const register double* f_011 = in.getSampleDataRO(INDEX3(k0,m_N1-1,k2+1, m_N0,m_N1));
                        const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,m_N1-1,k2, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(k0,m_N1-1,k2, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,m_N1-1,k2+1, m_N0,m_N1));
                        const register double* f_000 = in.getSampleDataRO(INDEX3(k0,m_N1-2,k2, m_N0,m_N1));
                        const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,m_N1-2,k2+1, m_N0,m_N1));
                        const register double* f_001 = in.getSampleDataRO(INDEX3(k0,m_N1-2,k2+1, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,m_N1-2,k2, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,3)] = (f_110[i]+f_111[i]-f_010[i]-f_011[i])*C4 / h0;
                            o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_011[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_100[i]-f_101[i])*C3 / h1;
                            o[INDEX3(i,2,0,numComp,3)] = (f_011[i]+f_111[i]-f_010[i]-f_110[i])*C4 / h2;
                        } /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k2 loop */
            } /* end of face 3 */
            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const register double* f_000 = in.getSampleDataRO(INDEX3(k0,k1,0, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,k1,0, m_N0,m_N1));
                        const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,k1+1,0, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(k0,k1+1,0, m_N0,m_N1));
                        const register double* f_001 = in.getSampleDataRO(INDEX3(k0,k1,1, m_N0,m_N1));
                        const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,k1,1, m_N0,m_N1));
                        const register double* f_011 = in.getSampleDataRO(INDEX3(k0,k1+1,1, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,k1+1,1, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_110[i]-f_000[i]-f_010[i])*C4 / h0;
                            o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_110[i]-f_000[i]-f_100[i])*C4 / h1;
                            o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_011[i]+f_101[i]+f_111[i]-f_000[i]-f_010[i]-f_100[i]-f_110[i])*C4 / h2;
                        } /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k1 loop */
            } /* end of face 4 */
            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const register double* f_001 = in.getSampleDataRO(INDEX3(k0,k1,m_N2-1, m_N0,m_N1));
                        const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,k1,m_N2-1, m_N0,m_N1));
                        const register double* f_011 = in.getSampleDataRO(INDEX3(k0,k1+1,m_N2-1, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,k1+1,m_N2-1, m_N0,m_N1));
                        const register double* f_000 = in.getSampleDataRO(INDEX3(k0,k1,m_N2-2, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,k1,m_N2-2, m_N0,m_N1));
                        const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,k1+1,m_N2-2, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(k0,k1+1,m_N2-2, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,3)] = (f_101[i]+f_111[i]-f_001[i]-f_011[i])*C4 / h0;
                            o[INDEX3(i,1,0,numComp,3)] = (f_011[i]+f_111[i]-f_001[i]-f_101[i])*C4 / h1;
                            o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_011[i]+f_101[i]+f_111[i]-f_000[i]-f_010[i]-f_100[i]-f_110[i])*C3 / h2;
                        } /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k1 loop */
            } /* end of face 5 */
        } // end of parallel section
        /* GENERATOR SNIP_GRAD_REDUCED_FACES BOTTOM */
    } else {
        stringstream msg;
        msg << "setToGradient() not implemented for "
            << functionSpaceTypeAsString(out.getFunctionSpace().getTypeCode());
        throw RipleyException(msg.str());
    }
}

void Brick::setToIntegrals(vector<double>& integrals, const escript::Data& arg) const
{
    escript::Data& in = *const_cast<escript::Data*>(&arg);
    const dim_t numComp = in.getDataPointSize();
    const double h0 = m_l0/m_gNE0;
    const double h1 = m_l1/m_gNE1;
    const double h2 = m_l2/m_gNE2;
    if (arg.getFunctionSpace().getTypeCode() == Elements) {
        const double w_0 = h0*h1*h2/8.;
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
#pragma omp for nowait
            for (index_t k2 = 0; k2 < m_NE2; ++k2) {
                for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                    for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                        const double* f = in.getSampleDataRO(INDEX3(k0, k1, k2, m_NE0, m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            const register double f_0 = f[INDEX2(i,0,numComp)];
                            const register double f_1 = f[INDEX2(i,1,numComp)];
                            const register double f_2 = f[INDEX2(i,2,numComp)];
                            const register double f_3 = f[INDEX2(i,3,numComp)];
                            const register double f_4 = f[INDEX2(i,4,numComp)];
                            const register double f_5 = f[INDEX2(i,5,numComp)];
                            const register double f_6 = f[INDEX2(i,6,numComp)];
                            const register double f_7 = f[INDEX2(i,7,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3+f_4+f_5+f_6+f_7)*w_0;
                        }  /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k1 loop */
            } /* end of k2 loop */

#pragma omp critical
            for (index_t i=0; i<numComp; i++)
                integrals[i]+=int_local[i];
        } // end of parallel section
    } else if (arg.getFunctionSpace().getTypeCode() == ReducedElements) {
        const double w_0 = h0*h1*h2;
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
#pragma omp for nowait
            for (index_t k2 = 0; k2 < m_NE2; ++k2) {
                for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                    for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                        const double* f = in.getSampleDataRO(INDEX3(k0, k1, k2, m_NE0, m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_0;
                        }  /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k1 loop */
            } /* end of k2 loop */

#pragma omp critical
            for (index_t i=0; i<numComp; i++)
                integrals[i]+=int_local[i];
        } // end of parallel section
    } else if (arg.getFunctionSpace().getTypeCode() == FaceElements) {
        const double w_0 = h1*h2/4.;
        const double w_1 = h0*h2/4.;
        const double w_2 = h0*h1/4.;
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k1=0; k1 < m_NE1; ++k1) {
                        const double* f = in.getSampleDataRO(m_faceOffset[0]+INDEX2(k1,k2,m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            const register double f_0 = f[INDEX2(i,0,numComp)];
                            const register double f_1 = f[INDEX2(i,1,numComp)];
                            const register double f_2 = f[INDEX2(i,2,numComp)];
                            const register double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_0;
                        }  /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k1=0; k1 < m_NE1; ++k1) {
                        const double* f = in.getSampleDataRO(m_faceOffset[1]+INDEX2(k1,k2,m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            const register double f_0 = f[INDEX2(i,0,numComp)];
                            const register double f_1 = f[INDEX2(i,1,numComp)];
                            const register double f_2 = f[INDEX2(i,2,numComp)];
                            const register double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_0;
                        }  /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const double* f = in.getSampleDataRO(m_faceOffset[2]+INDEX2(k0,k2,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            const register double f_0 = f[INDEX2(i,0,numComp)];
                            const register double f_1 = f[INDEX2(i,1,numComp)];
                            const register double f_2 = f[INDEX2(i,2,numComp)];
                            const register double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_1;
                        }  /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const double* f = in.getSampleDataRO(m_faceOffset[3]+INDEX2(k0,k2,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            const register double f_0 = f[INDEX2(i,0,numComp)];
                            const register double f_1 = f[INDEX2(i,1,numComp)];
                            const register double f_2 = f[INDEX2(i,2,numComp)];
                            const register double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_1;
                        }  /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            }

            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const double* f = in.getSampleDataRO(m_faceOffset[4]+INDEX2(k0,k1,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            const register double f_0 = f[INDEX2(i,0,numComp)];
                            const register double f_1 = f[INDEX2(i,1,numComp)];
                            const register double f_2 = f[INDEX2(i,2,numComp)];
                            const register double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_2;
                        }  /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            }

            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const double* f = in.getSampleDataRO(m_faceOffset[5]+INDEX2(k0,k1,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            const register double f_0 = f[INDEX2(i,0,numComp)];
                            const register double f_1 = f[INDEX2(i,1,numComp)];
                            const register double f_2 = f[INDEX2(i,2,numComp)];
                            const register double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_2;
                        }  /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            }

#pragma omp critical
            for (index_t i=0; i<numComp; i++)
                integrals[i]+=int_local[i];
        } // end of parallel section

    } else if (arg.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        const double w_0 = h1*h2;
        const double w_1 = h0*h2;
        const double w_2 = h0*h1;
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k1=0; k1 < m_NE1; ++k1) {
                        const double* f = in.getSampleDataRO(m_faceOffset[0]+INDEX2(k1,k2,m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_0;
                        }  /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k1=0; k1 < m_NE1; ++k1) {
                        const double* f = in.getSampleDataRO(m_faceOffset[1]+INDEX2(k1,k2,m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_0;
                        }  /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const double* f = in.getSampleDataRO(m_faceOffset[2]+INDEX2(k0,k2,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_1;
                        }  /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const double* f = in.getSampleDataRO(m_faceOffset[3]+INDEX2(k0,k2,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_1;
                        }  /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            }

            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const double* f = in.getSampleDataRO(m_faceOffset[4]+INDEX2(k0,k1,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_2;
                        }  /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            }

            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const double* f = in.getSampleDataRO(m_faceOffset[5]+INDEX2(k0,k1,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_2;
                        }  /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            }

#pragma omp critical
            for (index_t i=0; i<numComp; i++)
                integrals[i]+=int_local[i];
        } // end of parallel section

    } else {
        stringstream msg;
        msg << "setToIntegrals() not implemented for "
            << functionSpaceTypeAsString(arg.getFunctionSpace().getTypeCode());
        throw RipleyException(msg.str());
    }
}

void Brick::setToNormal(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == FaceElements) {
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE2; ++k2) {
                    for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                        double* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE1));
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
                for (index_t k2 = 0; k2 < m_NE2; ++k2) {
                    for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                        double* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE1));
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
                for (index_t k2 = 0; k2 < m_NE2; ++k2) {
                    for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE0));
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
                for (index_t k2 = 0; k2 < m_NE2; ++k2) {
                    for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE0));
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
                for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                    for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE0));
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
                for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                    for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE0));
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
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE2; ++k2) {
                    for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                        double* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE1));
                        *o++ = -1.;
                        *o++ = 0.;
                        *o = 0.;
                    }
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE2; ++k2) {
                    for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                        double* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE1));
                        *o++ = 1.;
                        *o++ = 0.;
                        *o = 0.;
                    }
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE2; ++k2) {
                    for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE0));
                        *o++ = 0.;
                        *o++ = -1.;
                        *o = 0.;
                    }
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE2; ++k2) {
                    for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE0));
                        *o++ = 0.;
                        *o++ = 1.;
                        *o = 0.;
                    }
                }
            }

            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                    for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE0));
                        *o++ = 0.;
                        *o++ = 0.;
                        *o = -1.;
                    }
                }
            }

            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                    for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE0));
                        *o++ = 0.;
                        *o++ = 0.;
                        *o = 1.;
                    }
                }
            }
        } // end of parallel section

    } else {
        stringstream msg;
        msg << "setToNormal() not implemented for "
            << functionSpaceTypeAsString(out.getFunctionSpace().getTypeCode());
        throw RipleyException(msg.str());
    }
}

Paso_SystemMatrixPattern* Brick::getPattern(bool reducedRowOrder,
                                            bool reducedColOrder) const
{
    if (reducedRowOrder || reducedColOrder)
        throw RipleyException("getPattern() not implemented for reduced order");

    throw RipleyException("getPattern() not implemented");
}

void Brick::Print_Mesh_Info(const bool full) const
{
    RipleyDomain::Print_Mesh_Info(full);
    if (full) {
        cout << "     Id  Coordinates" << endl;
        cout.precision(15);
        cout.setf(ios::scientific, ios::floatfield);
        pair<double,double> xdx = getFirstCoordAndSpacing(0);
        pair<double,double> ydy = getFirstCoordAndSpacing(1);
        pair<double,double> zdz = getFirstCoordAndSpacing(2);
        for (index_t i=0; i < getNumNodes(); i++) {
            cout << "  " << setw(5) << m_nodeId[i]
                << "  " << xdx.first+(i%m_N0)*xdx.second
                << "  " << ydy.first+(i%(m_N0*m_N1)/m_N0)*ydy.second
                << "  " << zdz.first+(i/(m_N0*m_N1))*zdz.second << endl;
        }
    }
}

IndexVector Brick::getNumNodesPerDim() const
{
    IndexVector ret;
    ret.push_back(m_N0);
    ret.push_back(m_N1);
    ret.push_back(m_N2);
    return ret;
}

IndexVector Brick::getNumElementsPerDim() const
{
    IndexVector ret;
    ret.push_back(m_NE0);
    ret.push_back(m_NE1);
    ret.push_back(m_NE2);
    return ret;
}

IndexVector Brick::getNumFacesPerBoundary() const
{
    IndexVector ret(6, 0);
    //left
    if (m_offset0==0)
        ret[0]=m_NE1*m_NE2;
    //right
    if (m_mpiInfo->rank%m_NX==m_NX-1)
        ret[1]=m_NE1*m_NE2;
    //bottom
    if (m_offset1==0)
        ret[2]=m_NE0*m_NE2;
    //top
    if (m_mpiInfo->rank%(m_NX*m_NY)/m_NX==m_NY-1)
        ret[3]=m_NE0*m_NE2;
    //front
    if (m_offset2==0)
        ret[4]=m_NE0*m_NE1;
    //back
    if (m_mpiInfo->rank/(m_NX*m_NY)==m_NZ-1)
        ret[5]=m_NE0*m_NE1;
    return ret;
}

pair<double,double> Brick::getFirstCoordAndSpacing(dim_t dim) const
{
    if (dim==0)
        return pair<double,double>((m_l0*m_offset0)/m_gNE0, m_l0/m_gNE0);
    else if (dim==1)
        return pair<double,double>((m_l1*m_offset1)/m_gNE1, m_l1/m_gNE1);
    else if (dim==2)
        return pair<double,double>((m_l2*m_offset2)/m_gNE2, m_l2/m_gNE2);

    throw RipleyException("getFirstCoordAndSpacing(): invalid argument");
}


//protected
dim_t Brick::getNumFaceElements() const
{
    dim_t n=0;
    //left
    if (m_offset0==0)
        n+=m_NE1*m_NE2;
    //right
    if (m_mpiInfo->rank%m_NX==m_NX-1)
        n+=m_NE1*m_NE2;
    //bottom
    if (m_offset1==0)
        n+=m_NE0*m_NE2;
    //top
    if (m_mpiInfo->rank%(m_NX*m_NY)/m_NX==m_NY-1)
        n+=m_NE0*m_NE2;
    //front
    if (m_offset2==0)
        n+=m_NE0*m_NE1;
    //back
    if (m_mpiInfo->rank/(m_NX*m_NY)==m_NZ-1)
        n+=m_NE0*m_NE1;

    return n;
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

    pair<double,double> xdx = getFirstCoordAndSpacing(0);
    pair<double,double> ydy = getFirstCoordAndSpacing(1);
    pair<double,double> zdz = getFirstCoordAndSpacing(2);
    arg.requireWrite();
#pragma omp parallel for
    for (dim_t i2 = 0; i2 < m_N2; i2++) {
        for (dim_t i1 = 0; i1 < m_N1; i1++) {
            for (dim_t i0 = 0; i0 < m_N0; i0++) {
                double* point = arg.getSampleDataRW(i0+m_N0*i1+m_N0*m_N1*i2);
                point[0] = xdx.first+i0*xdx.second;
                point[1] = ydy.first+i1*ydy.second;
                point[2] = zdz.first+i2*zdz.second;
            }
        }
    }
}

//private
void Brick::populateSampleIds()
{
    // identifiers are ordered from left to right, bottom to top, front to back
    // on each rank, except for the shared nodes which are owned by the rank
    // below / to the left / to the front of the current rank

    // build node distribution vector first.
    // m_nodeDistribution[i] is the first node id on rank i, that is
    // rank i owns m_nodeDistribution[i+1]-nodeDistribution[i] nodes
    m_nodeDistribution.assign(m_mpiInfo->size+1, 0);
    m_nodeDistribution[1]=getNumNodes();
    for (dim_t k=1; k<m_mpiInfo->size-1; k++) {
        const index_t x = k%m_NX;
        const index_t y = k%(m_NX*m_NY)/m_NX;
        const index_t z = k/(m_NX*m_NY);
        index_t numNodes=getNumNodes();
        if (x>0)
            numNodes-=m_N1*m_N2;
        if (y>0)
            numNodes-=m_N0*m_N2;
        if (z>0)
            numNodes-=m_N0*m_N1;
        // if an edge was subtracted twice add it back
        if (x>0 && y>0)
            numNodes+=m_N2; 
        if (x>0 && z>0)
            numNodes+=m_N1;
        if (y>0 && z>0)
            numNodes+=m_N0;
        // the corner node was removed 3x and added back 3x, so subtract it
        if (x>0 && y>0 && z>0)
            numNodes--;
        m_nodeDistribution[k+1]=m_nodeDistribution[k]+numNodes;
    }
    m_nodeDistribution[m_mpiInfo->size]=getNumDataPointsGlobal();

    m_nodeId.resize(getNumNodes());

    // the bottom, left and front planes are not owned by this rank so the
    // identifiers need to be computed accordingly
    const index_t left = (m_offset0==0 ? 0 : 1);
    const index_t bottom = (m_offset1==0 ? 0 : 1);
    const index_t front = (m_offset2==0 ? 0 : 1);

    // case 1: all nodes on left plane are owned by rank on the left
    if (left>0) {
        const int neighbour=m_mpiInfo->rank-1;
        const index_t leftN0=(neighbour%m_NX == 0 ? m_N0 : m_N0-1);
        const index_t leftN1=(neighbour%(m_NX*m_NY)/m_NX==0 ? m_N1 : m_N1-1);
#pragma omp parallel for
        for (dim_t i2=front; i2<m_N2; i2++) {
            for (dim_t i1=bottom; i1<m_N1; i1++) {
                m_nodeId[i1*m_N0+i2*m_N0*m_N1]=m_nodeDistribution[neighbour]
                    + (i1-bottom+1)*leftN0
                    + (i2-front)*leftN0*leftN1 - 1;
            }
        }
    }
    // case 2: all nodes on bottom plane are owned by rank below
    if (bottom>0) {
        const int neighbour=m_mpiInfo->rank-m_NX;
        const index_t bottomN0=(neighbour%m_NX == 0 ? m_N0 : m_N0-1);
        const index_t bottomN1=(neighbour%(m_NX*m_NY)/m_NX==0 ? m_N1 : m_N1-1);
#pragma omp parallel for
        for (dim_t i2=front; i2<m_N2; i2++) {
            for (dim_t i0=left; i0<m_N0; i0++) {
                m_nodeId[i0+i2*m_N0*m_N1]=m_nodeDistribution[neighbour]
                    + bottomN0*(bottomN1-1)
                    + (i2-front)*bottomN0*bottomN1 + i0-left;
            }
        }
    }
    // case 3: all nodes on front plane are owned by rank in front
    if (front>0) {
        const int neighbour=m_mpiInfo->rank-m_NX*m_NY;
        const index_t N0=(neighbour%m_NX == 0 ? m_N0 : m_N0-1);
        const index_t N1=(neighbour%(m_NX*m_NY)/m_NX==0 ? m_N1 : m_N1-1);
        const index_t N2=(neighbour/(m_NX*m_NY)==0 ? m_N2 : m_N2-1);
#pragma omp parallel for
        for (dim_t i1=bottom; i1<m_N1; i1++) {
            for (dim_t i0=left; i0<m_N0; i0++) {
                m_nodeId[i0+i1*m_N0]=m_nodeDistribution[neighbour]
                    + N0*N1*(N2-1)+(i1-bottom)*N0 + i0-left;
            }
        }
    }
    // case 4: nodes on front bottom edge are owned by the corresponding rank
    if (front>0 && bottom>0) {
        const int neighbour=m_mpiInfo->rank-m_NX*(m_NY+1);
        const index_t N0=(neighbour%m_NX == 0 ? m_N0 : m_N0-1);
        const index_t N1=(neighbour%(m_NX*m_NY)/m_NX==0 ? m_N1 : m_N1-1);
        const index_t N2=(neighbour/(m_NX*m_NY)==0 ? m_N2 : m_N2-1);
#pragma omp parallel for
        for (dim_t i0=left; i0<m_N0; i0++) {
            m_nodeId[i0]=m_nodeDistribution[neighbour]
                + N0*N1*(N2-1)+(N1-1)*N0 + i0-left;
        }
    }
    // case 5: nodes on left bottom edge are owned by the corresponding rank
    if (left>0 && bottom>0) {
        const int neighbour=m_mpiInfo->rank-m_NX-1;
        const index_t N0=(neighbour%m_NX == 0 ? m_N0 : m_N0-1);
        const index_t N1=(neighbour%(m_NX*m_NY)/m_NX==0 ? m_N1 : m_N1-1);
#pragma omp parallel for
        for (dim_t i2=front; i2<m_N2; i2++) {
            m_nodeId[i2*m_N0*m_N1]=m_nodeDistribution[neighbour]
                + (1+i2-front)*N0*N1-1;
        }
    }
    // case 6: nodes on left front edge are owned by the corresponding rank
    if (left>0 && front>0) {
        const int neighbour=m_mpiInfo->rank-m_NX*m_NY-1;
        const index_t N0=(neighbour%m_NX == 0 ? m_N0 : m_N0-1);
        const index_t N1=(neighbour%(m_NX*m_NY)/m_NX==0 ? m_N1 : m_N1-1);
        const index_t N2=(neighbour/(m_NX*m_NY)==0 ? m_N2 : m_N2-1);
#pragma omp parallel for
        for (dim_t i1=bottom; i1<m_N1; i1++) {
            m_nodeId[i1*m_N0]=m_nodeDistribution[neighbour]
                + N0*N1*(N2-1)+N0-1+(i1-bottom)*N0;
        }
    }
    // case 7: bottom-left-front corner node owned by corresponding rank
    if (left>0 && bottom>0 && front>0) {
        const int neighbour=m_mpiInfo->rank-m_NX*(m_NY+1)-1;
        const index_t N0=(neighbour%m_NX == 0 ? m_N0 : m_N0-1);
        const index_t N1=(neighbour%(m_NX*m_NY)/m_NX==0 ? m_N1 : m_N1-1);
        const index_t N2=(neighbour/(m_NX*m_NY) == 0 ? m_N2 : m_N2-1);
        m_nodeId[0]=m_nodeDistribution[neighbour]+N0*N1*N2-1;
    }

    // the rest of the id's are contiguous
    const index_t firstId=m_nodeDistribution[m_mpiInfo->rank];
#pragma omp parallel for
    for (dim_t i2=front; i2<m_N2; i2++) {
        for (dim_t i1=bottom; i1<m_N1; i1++) {
            for (dim_t i0=left; i0<m_N0; i0++) {
                m_nodeId[i0+i1*m_N0+i2*m_N0*m_N1] = firstId+i0-left
                    +(i1-bottom)*(m_N0-left)
                    +(i2-front)*(m_N0-left)*(m_N1-bottom);
            }
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

    // generate face offset vector and set face tags
    const IndexVector facesPerEdge = getNumFacesPerBoundary();
    const index_t LEFT=1, RIGHT=2, BOTTOM=10, TOP=20, FRONT=100, BACK=200;
    const index_t faceTag[] = { LEFT, RIGHT, BOTTOM, TOP, FRONT, BACK };
    m_faceOffset.assign(facesPerEdge.size(), -1);
    m_faceTags.clear();
    index_t offset=0;
    for (size_t i=0; i<facesPerEdge.size(); i++) {
        if (facesPerEdge[i]>0) {
            m_faceOffset[i]=offset;
            offset+=facesPerEdge[i];
            m_faceTags.insert(m_faceTags.end(), facesPerEdge[i], faceTag[i]);
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

//protected
void Brick::interpolateNodesOnElements(escript::Data& out, escript::Data& in,
                                       bool reduced) const
{
    const dim_t numComp = in.getDataPointSize();
    if (reduced) {
        /*** GENERATOR SNIP_INTERPOLATE_REDUCED_ELEMENTS TOP */
        const double c0 = .125;
#pragma omp parallel for
        for (index_t k2=0; k2 < m_NE2; ++k2) {
            for (index_t k1=0; k1 < m_NE1; ++k1) {
                for (index_t k0=0; k0 < m_NE0; ++k0) {
                    const register double* f_000 = in.getSampleDataRO(INDEX3(k0,k1,k2, m_N0,m_N1));
                    const register double* f_001 = in.getSampleDataRO(INDEX3(k0,k1,k2+1, m_N0,m_N1));
                    const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,k1,k2+1, m_N0,m_N1));
                    const register double* f_011 = in.getSampleDataRO(INDEX3(k0,k1+1,k2+1, m_N0,m_N1));
                    const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,k1,k2, m_N0,m_N1));
                    const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,k1+1,k2, m_N0,m_N1));
                    const register double* f_010 = in.getSampleDataRO(INDEX3(k0,k1+1,k2, m_N0,m_N1));
                    const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,k1+1,k2+1, m_N0,m_N1));
                    double* o = out.getSampleDataRW(INDEX3(k0,k1,k2,m_NE0,m_NE1));
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = c0*(f_000[i] + f_001[i] + f_010[i] + f_011[i] + f_100[i] + f_101[i] + f_110[i] + f_111[i]);
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of k1 loop */
        } /* end of k2 loop */
        /* GENERATOR SNIP_INTERPOLATE_REDUCED_ELEMENTS BOTTOM */
    } else {
        /*** GENERATOR SNIP_INTERPOLATE_ELEMENTS TOP */
        const double c0 = .0094373878376559314545;
        const double c1 = .035220810900864519624;
        const double c2 = .13144585576580214704;
        const double c3 = .49056261216234406855;
#pragma omp parallel for
        for (index_t k2=0; k2 < m_NE2; ++k2) {
            for (index_t k1=0; k1 < m_NE1; ++k1) {
                for (index_t k0=0; k0 < m_NE0; ++k0) {
                    const register double* f_000 = in.getSampleDataRO(INDEX3(k0,k1,k2, m_N0,m_N1));
                    const register double* f_001 = in.getSampleDataRO(INDEX3(k0,k1,k2+1, m_N0,m_N1));
                    const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,k1,k2+1, m_N0,m_N1));
                    const register double* f_011 = in.getSampleDataRO(INDEX3(k0,k1+1,k2+1, m_N0,m_N1));
                    const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,k1+1,k2, m_N0,m_N1));
                    const register double* f_010 = in.getSampleDataRO(INDEX3(k0,k1+1,k2, m_N0,m_N1));
                    const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,k1,k2, m_N0,m_N1));
                    const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,k1+1,k2+1, m_N0,m_N1));
                    double* o = out.getSampleDataRW(INDEX3(k0,k1,k2,m_NE0,m_NE1));
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = f_000[i]*c3 + f_111[i]*c0 + c2*(f_001[i] + f_010[i] + f_100[i]) + c1*(f_011[i] + f_101[i] + f_110[i]);
                        o[INDEX2(i,numComp,1)] = f_011[i]*c0 + f_100[i]*c3 + c2*(f_000[i] + f_101[i] + f_110[i]) + c1*(f_001[i] + f_010[i] + f_111[i]);
                        o[INDEX2(i,numComp,2)] = f_010[i]*c3 + f_101[i]*c0 + c2*(f_000[i] + f_011[i] + f_110[i]) + c1*(f_001[i] + f_100[i] + f_111[i]);
                        o[INDEX2(i,numComp,3)] = f_001[i]*c0 + f_110[i]*c3 + c2*(f_010[i] + f_100[i] + f_111[i]) + c1*(f_000[i] + f_011[i] + f_101[i]);
                        o[INDEX2(i,numComp,4)] = f_001[i]*c3 + f_110[i]*c0 + c2*(f_000[i] + f_011[i] + f_101[i]) + c1*(f_010[i] + f_100[i] + f_111[i]);
                        o[INDEX2(i,numComp,5)] = f_010[i]*c0 + f_101[i]*c3 + c2*(f_001[i] + f_100[i] + f_111[i]) + c1*(f_000[i] + f_011[i] + f_110[i]);
                        o[INDEX2(i,numComp,6)] = f_011[i]*c3 + f_100[i]*c0 + c2*(f_001[i] + f_010[i] + f_111[i]) + c1*(f_000[i] + f_101[i] + f_110[i]);
                        o[INDEX2(i,numComp,7)] = f_000[i]*c0 + f_111[i]*c3 + c2*(f_011[i] + f_101[i] + f_110[i]) + c1*(f_001[i] + f_010[i] + f_100[i]);
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of k1 loop */
        } /* end of k2 loop */
        /* GENERATOR SNIP_INTERPOLATE_ELEMENTS BOTTOM */
    }
}

//protected
void Brick::interpolateNodesOnFaces(escript::Data& out, escript::Data& in,
                                    bool reduced) const
{
    const dim_t numComp = in.getDataPointSize();
    if (reduced) {
        const double c0 = .25;
#pragma omp parallel
        {
            /*** GENERATOR SNIP_INTERPOLATE_REDUCED_FACES TOP */
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k1=0; k1 < m_NE1; ++k1) {
                        const register double* f_011 = in.getSampleDataRO(INDEX3(0,k1+1,k2+1, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(0,k1+1,k2, m_N0,m_N1));
                        const register double* f_001 = in.getSampleDataRO(INDEX3(0,k1,k2+1, m_N0,m_N1));
                        const register double* f_000 = in.getSampleDataRO(INDEX3(0,k1,k2, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = c0*(f_000[i] + f_001[i] + f_010[i] + f_011[i]);
                        } /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            } /* end of face 0 */
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k1=0; k1 < m_NE1; ++k1) {
                        const register double* f_110 = in.getSampleDataRO(INDEX3(m_N0-1,k1+1,k2, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(m_N0-1,k1,k2, m_N0,m_N1));
                        const register double* f_101 = in.getSampleDataRO(INDEX3(m_N0-1,k1,k2+1, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(m_N0-1,k1+1,k2+1, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = c0*(f_100[i] + f_101[i] + f_110[i] + f_111[i]);
                        } /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            } /* end of face 1 */
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const register double* f_000 = in.getSampleDataRO(INDEX3(k0,0,k2, m_N0,m_N1));
                        const register double* f_001 = in.getSampleDataRO(INDEX3(k0,0,k2+1, m_N0,m_N1));
                        const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,0,k2+1, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,0,k2, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = c0*(f_000[i] + f_001[i] + f_100[i] + f_101[i]);
                        } /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k2 loop */
            } /* end of face 2 */
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,m_N1-1,k2, m_N0,m_N1));
                        const register double* f_011 = in.getSampleDataRO(INDEX3(k0,m_N1-1,k2+1, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(k0,m_N1-1,k2, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,m_N1-1,k2+1, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = c0*(f_010[i] + f_011[i] + f_110[i] + f_111[i]);
                        } /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k2 loop */
            } /* end of face 3 */
            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,k1+1,0, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(k0,k1+1,0, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,k1,0, m_N0,m_N1));
                        const register double* f_000 = in.getSampleDataRO(INDEX3(k0,k1,0, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = c0*(f_000[i] + f_010[i] + f_100[i] + f_110[i]);
                        } /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k1 loop */
            } /* end of face 4 */
            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const register double* f_011 = in.getSampleDataRO(INDEX3(k0,k1+1,m_N2-1, m_N0,m_N1));
                        const register double* f_001 = in.getSampleDataRO(INDEX3(k0,k1,m_N2-1, m_N0,m_N1));
                        const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,k1,m_N2-1, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,k1+1,m_N2-1, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = c0*(f_001[i] + f_011[i] + f_101[i] + f_111[i]);
                        } /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k1 loop */
            } /* end of face 5 */
            /* GENERATOR SNIP_INTERPOLATE_REDUCED_FACES BOTTOM */
        } // end of parallel section
    } else {
        const double c0 = 0.044658198738520451079;
        const double c1 = 0.16666666666666666667;
        const double c2 = 0.62200846792814621559;
#pragma omp parallel
        {
            /*** GENERATOR SNIP_INTERPOLATE_FACES TOP */
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k1=0; k1 < m_NE1; ++k1) {
                        const register double* f_000 = in.getSampleDataRO(INDEX3(0,k1,k2, m_N0,m_N1));
                        const register double* f_001 = in.getSampleDataRO(INDEX3(0,k1,k2+1, m_N0,m_N1));
                        const register double* f_011 = in.getSampleDataRO(INDEX3(0,k1+1,k2+1, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(0,k1+1,k2, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = f_000[i]*c2 + f_011[i]*c0 + c1*(f_001[i] + f_010[i]);
                            o[INDEX2(i,numComp,1)] = f_001[i]*c0 + f_010[i]*c2 + c1*(f_000[i] + f_011[i]);
                            o[INDEX2(i,numComp,2)] = f_001[i]*c2 + f_010[i]*c0 + c1*(f_000[i] + f_011[i]);
                            o[INDEX2(i,numComp,3)] = f_000[i]*c0 + f_011[i]*c2 + c1*(f_001[i] + f_010[i]);
                        } /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            } /* end of face 0 */
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k1=0; k1 < m_NE1; ++k1) {
                        const register double* f_101 = in.getSampleDataRO(INDEX3(m_N0-1,k1,k2+1, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(m_N0-1,k1,k2, m_N0,m_N1));
                        const register double* f_110 = in.getSampleDataRO(INDEX3(m_N0-1,k1+1,k2, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(m_N0-1,k1+1,k2+1, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = f_100[i]*c2 + f_111[i]*c0 + c1*(f_101[i] + f_110[i]);
                            o[INDEX2(i,numComp,1)] = f_101[i]*c0 + f_110[i]*c2 + c1*(f_100[i] + f_111[i]);
                            o[INDEX2(i,numComp,2)] = f_101[i]*c2 + f_110[i]*c0 + c1*(f_100[i] + f_111[i]);
                            o[INDEX2(i,numComp,3)] = f_100[i]*c0 + f_111[i]*c2 + c1*(f_101[i] + f_110[i]);
                        } /* end of component loop i */
                    } /* end of k1 loop */
                } /* end of k2 loop */
            } /* end of face 1 */
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const register double* f_000 = in.getSampleDataRO(INDEX3(k0,0,k2, m_N0,m_N1));
                        const register double* f_001 = in.getSampleDataRO(INDEX3(k0,0,k2+1, m_N0,m_N1));
                        const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,0,k2+1, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,0,k2, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = f_000[i]*c2 + f_101[i]*c0 + c1*(f_001[i] + f_100[i]);
                            o[INDEX2(i,numComp,1)] = f_001[i]*c0 + f_100[i]*c2 + c1*(f_000[i] + f_101[i]);
                            o[INDEX2(i,numComp,2)] = f_001[i]*c2 + f_100[i]*c0 + c1*(f_000[i] + f_101[i]);
                            o[INDEX2(i,numComp,3)] = f_000[i]*c0 + f_101[i]*c2 + c1*(f_001[i] + f_100[i]);
                        } /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k2 loop */
            } /* end of face 2 */
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2=0; k2 < m_NE2; ++k2) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,m_N1-1,k2, m_N0,m_N1));
                        const register double* f_011 = in.getSampleDataRO(INDEX3(k0,m_N1-1,k2+1, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(k0,m_N1-1,k2, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,m_N1-1,k2+1, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = f_010[i]*c2 + f_111[i]*c0 + c1*(f_011[i] + f_110[i]);
                            o[INDEX2(i,numComp,1)] = f_011[i]*c0 + f_110[i]*c2 + c1*(f_010[i] + f_111[i]);
                            o[INDEX2(i,numComp,2)] = f_011[i]*c2 + f_110[i]*c0 + c1*(f_010[i] + f_111[i]);
                            o[INDEX2(i,numComp,3)] = f_010[i]*c0 + f_111[i]*c2 + c1*(f_011[i] + f_110[i]);
                        } /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k2 loop */
            } /* end of face 3 */
            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const register double* f_000 = in.getSampleDataRO(INDEX3(k0,k1,0, m_N0,m_N1));
                        const register double* f_100 = in.getSampleDataRO(INDEX3(k0+1,k1,0, m_N0,m_N1));
                        const register double* f_110 = in.getSampleDataRO(INDEX3(k0+1,k1+1,0, m_N0,m_N1));
                        const register double* f_010 = in.getSampleDataRO(INDEX3(k0,k1+1,0, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = f_000[i]*c2 + f_110[i]*c0 + c1*(f_010[i] + f_100[i]);
                            o[INDEX2(i,numComp,1)] = f_010[i]*c0 + f_100[i]*c2 + c1*(f_000[i] + f_110[i]);
                            o[INDEX2(i,numComp,2)] = f_010[i]*c2 + f_100[i]*c0 + c1*(f_000[i] + f_110[i]);
                            o[INDEX2(i,numComp,3)] = f_000[i]*c0 + f_110[i]*c2 + c1*(f_010[i] + f_100[i]);
                        } /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k1 loop */
            } /* end of face 4 */
            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    for (index_t k0=0; k0 < m_NE0; ++k0) {
                        const register double* f_001 = in.getSampleDataRO(INDEX3(k0,k1,m_N2-1, m_N0,m_N1));
                        const register double* f_101 = in.getSampleDataRO(INDEX3(k0+1,k1,m_N2-1, m_N0,m_N1));
                        const register double* f_011 = in.getSampleDataRO(INDEX3(k0,k1+1,m_N2-1, m_N0,m_N1));
                        const register double* f_111 = in.getSampleDataRO(INDEX3(k0+1,k1+1,m_N2-1, m_N0,m_N1));
                        double* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            o[INDEX2(i,numComp,0)] = f_001[i]*c2 + f_111[i]*c0 + c1*(f_011[i] + f_101[i]);
                            o[INDEX2(i,numComp,1)] = f_011[i]*c0 + f_101[i]*c2 + c1*(f_001[i] + f_111[i]);
                            o[INDEX2(i,numComp,2)] = f_011[i]*c2 + f_101[i]*c0 + c1*(f_001[i] + f_111[i]);
                            o[INDEX2(i,numComp,3)] = f_001[i]*c0 + f_111[i]*c2 + c1*(f_011[i] + f_101[i]);
                        } /* end of component loop i */
                    } /* end of k0 loop */
                } /* end of k1 loop */
            } /* end of face 5 */
            /* GENERATOR SNIP_INTERPOLATE_FACES BOTTOM */
        } // end of parallel section
    }
}

} // end of namespace ripley

