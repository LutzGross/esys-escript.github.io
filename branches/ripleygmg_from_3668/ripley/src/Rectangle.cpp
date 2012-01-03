
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
#include "paso/SystemMatrix.h"
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

    if ((n0+1)%m_NX > 0 || (n1+1)%m_NY > 0)
        throw RipleyException("Number of elements+1 must be separable into number of ranks in each dimension");

    if ((m_NX > 1 && (n0+1)/m_NX<2) || (m_NY > 1 && (n1+1)/m_NY<2))
        throw RipleyException("Too few elements for the number of ranks");

    // local number of elements (including overlap)
    m_NE0 = (m_NX>1 ? (n0+1)/m_NX : n0);
    if (m_mpiInfo->rank%m_NX>0 && m_mpiInfo->rank%m_NX<m_NX-1)
        m_NE0++;
    m_NE1 = (m_NY>1 ? (n1+1)/m_NY : n1);
    if (m_mpiInfo->rank/m_NX>0 && m_mpiInfo->rank/m_NX<m_NY-1)
        m_NE1++;

    // local number of nodes
    m_N0 = m_NE0+1;
    m_N1 = m_NE1+1;

    // bottom-left node is at (offset0,offset1) in global mesh
    m_offset0 = (n0+1)/m_NX*(m_mpiInfo->rank%m_NX);
    if (m_offset0 > 0)
        m_offset0--;
    m_offset1 = (n1+1)/m_NY*(m_mpiInfo->rank/m_NX);
    if (m_offset1 > 0)
        m_offset1--;

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
    const Rectangle* o=dynamic_cast<const Rectangle*>(&other);
    if (o) {
        return (RipleyDomain::operator==(other) &&
                m_gNE0==o->m_gNE0 && m_gNE1==o->m_gNE1
                && m_l0==o->m_l0 && m_l1==o->m_l1
                && m_NX==o->m_NX && m_NY==o->m_NY);
    }

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
#pragma omp for nowait
        for (dim_t i0 = 0; i0 < m_N0; i0++) {
            coords[0][i0]=xdx.first+i0*xdx.second;
        }
#pragma omp for nowait
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
    msg << "borrowSampleReferenceIDs() not implemented for function space type "
        << functionSpaceTypeAsString(fsType);
    throw RipleyException(msg.str());
}

bool Rectangle::ownSample(int fsCode, index_t id) const
{
#ifdef ESYS_MPI
    if (fsCode == Nodes) {
        const index_t left = (m_offset0==0 ? 0 : 1);
        const index_t bottom = (m_offset1==0 ? 0 : 1);
        const index_t right = (m_mpiInfo->rank%m_NX==m_NX-1 ? m_N0 : m_N0-1);
        const index_t top = (m_mpiInfo->rank/m_NX==m_NY-1 ? m_N1 : m_N1-1);
        const index_t x=id%m_N0;
        const index_t y=id/m_N0;
        return (x>=left && x<right && y>=bottom && y<top);
    } else {
        stringstream msg;
        msg << "ownSample() not implemented for "
            << functionSpaceTypeAsString(fsCode);
        throw RipleyException(msg.str());
    }
#else
    return true;
#endif
}

void Rectangle::setToGradient(escript::Data& out, const escript::Data& cIn) const
{
    escript::Data& in = *const_cast<escript::Data*>(&cIn);
    const dim_t numComp = in.getDataPointSize();
    const double h0 = m_l0/m_gNE0;
    const double h1 = m_l1/m_gNE1;
    const double cx0 = -1./h0;
    const double cx1 = -.78867513459481288225/h0;
    const double cx2 = -.5/h0;
    const double cx3 = -.21132486540518711775/h0;
    const double cx4 = .21132486540518711775/h0;
    const double cx5 = .5/h0;
    const double cx6 = .78867513459481288225/h0;
    const double cx7 = 1./h0;
    const double cy0 = -1./h1;
    const double cy1 = -.78867513459481288225/h1;
    const double cy2 = -.5/h1;
    const double cy3 = -.21132486540518711775/h1;
    const double cy4 = .21132486540518711775/h1;
    const double cy5 = .5/h1;
    const double cy6 = .78867513459481288225/h1;
    const double cy7 = 1./h1;

    if (out.getFunctionSpace().getTypeCode() == Elements) {
        /*** GENERATOR SNIP_GRAD_ELEMENTS TOP */
#pragma omp parallel for
        for (index_t k1=0; k1 < m_NE1; ++k1) {
            for (index_t k0=0; k0 < m_NE0; ++k0) {
                const register double* f_10 = in.getSampleDataRO(INDEX2(k0+1,k1, m_N0));
                const register double* f_11 = in.getSampleDataRO(INDEX2(k0+1,k1+1, m_N0));
                const register double* f_01 = in.getSampleDataRO(INDEX2(k0,k1+1, m_N0));
                const register double* f_00 = in.getSampleDataRO(INDEX2(k0,k1, m_N0));
                double* o = out.getSampleDataRW(INDEX2(k0,k1,m_NE0));
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX3(i,0,0,numComp,2)] = f_00[i]*cx1 + f_01[i]*cx3 + f_10[i]*cx6 + f_11[i]*cx4;
                    o[INDEX3(i,1,0,numComp,2)] = f_00[i]*cy1 + f_01[i]*cy6 + f_10[i]*cy3 + f_11[i]*cy4;
                    o[INDEX3(i,0,1,numComp,2)] = f_00[i]*cx1 + f_01[i]*cx3 + f_10[i]*cx6 + f_11[i]*cx4;
                    o[INDEX3(i,1,1,numComp,2)] = f_00[i]*cy3 + f_01[i]*cy4 + f_10[i]*cy1 + f_11[i]*cy6;
                    o[INDEX3(i,0,2,numComp,2)] = f_00[i]*cx3 + f_01[i]*cx1 + f_10[i]*cx4 + f_11[i]*cx6;
                    o[INDEX3(i,1,2,numComp,2)] = f_00[i]*cy1 + f_01[i]*cy6 + f_10[i]*cy3 + f_11[i]*cy4;
                    o[INDEX3(i,0,3,numComp,2)] = f_00[i]*cx3 + f_01[i]*cx1 + f_10[i]*cx4 + f_11[i]*cx6;
                    o[INDEX3(i,1,3,numComp,2)] = f_00[i]*cy3 + f_01[i]*cy4 + f_10[i]*cy1 + f_11[i]*cy6;
                } /* end of component loop i */
            } /* end of k0 loop */
        } /* end of k1 loop */
        /* GENERATOR SNIP_GRAD_ELEMENTS BOTTOM */
    } else if (out.getFunctionSpace().getTypeCode() == ReducedElements) {
        /*** GENERATOR SNIP_GRAD_REDUCED_ELEMENTS TOP */
#pragma omp parallel for
        for (index_t k1=0; k1 < m_NE1; ++k1) {
            for (index_t k0=0; k0 < m_NE0; ++k0) {
                const register double* f_10 = in.getSampleDataRO(INDEX2(k0+1,k1, m_N0));
                const register double* f_11 = in.getSampleDataRO(INDEX2(k0+1,k1+1, m_N0));
                const register double* f_01 = in.getSampleDataRO(INDEX2(k0,k1+1, m_N0));
                const register double* f_00 = in.getSampleDataRO(INDEX2(k0,k1, m_N0));
                double* o = out.getSampleDataRW(INDEX2(k0,k1,m_NE0));
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX3(i,0,0,numComp,2)] = cx5*(f_10[i] + f_11[i]) + cx2*(f_00[i] + f_01[i]);
                    o[INDEX3(i,1,0,numComp,2)] = cy2*(f_00[i] + f_10[i]) + cy5*(f_01[i] + f_11[i]);
                } /* end of component loop i */
            } /* end of k0 loop */
        } /* end of k1 loop */
        /* GENERATOR SNIP_GRAD_REDUCED_ELEMENTS BOTTOM */
    } else if (out.getFunctionSpace().getTypeCode() == FaceElements) {
#pragma omp parallel
        {
            /*** GENERATOR SNIP_GRAD_FACES TOP */
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    const register double* f_10 = in.getSampleDataRO(INDEX2(1,k1, m_N0));
                    const register double* f_11 = in.getSampleDataRO(INDEX2(1,k1+1, m_N0));
                    const register double* f_01 = in.getSampleDataRO(INDEX2(0,k1+1, m_N0));
                    const register double* f_00 = in.getSampleDataRO(INDEX2(0,k1, m_N0));
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = f_00[i]*cx1 + f_01[i]*cx3 + f_10[i]*cx6 + f_11[i]*cx4;
                        o[INDEX3(i,1,0,numComp,2)] = f_00[i]*cy0 + f_01[i]*cy7;
                        o[INDEX3(i,0,1,numComp,2)] = f_00[i]*cx3 + f_01[i]*cx1 + f_10[i]*cx4 + f_11[i]*cx6;
                        o[INDEX3(i,1,1,numComp,2)] = f_00[i]*cy0 + f_01[i]*cy7;
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of face 0 */
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    const register double* f_10 = in.getSampleDataRO(INDEX2(m_N0-1,k1, m_N0));
                    const register double* f_11 = in.getSampleDataRO(INDEX2(m_N0-1,k1+1, m_N0));
                    const register double* f_01 = in.getSampleDataRO(INDEX2(m_N0-2,k1+1, m_N0));
                    const register double* f_00 = in.getSampleDataRO(INDEX2(m_N0-2,k1, m_N0));
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = f_00[i]*cx1 + f_01[i]*cx3 + f_10[i]*cx6 + f_11[i]*cx4;
                        o[INDEX3(i,1,0,numComp,2)] = f_10[i]*cy0 + f_11[i]*cy7;
                        o[INDEX3(i,0,1,numComp,2)] = f_00[i]*cx3 + f_01[i]*cx1 + f_10[i]*cx4 + f_11[i]*cx6;
                        o[INDEX3(i,1,1,numComp,2)] = f_10[i]*cy0 + f_11[i]*cy7;
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of face 1 */
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < m_NE0; ++k0) {
                    const register double* f_00 = in.getSampleDataRO(INDEX2(k0,0, m_N0));
                    const register double* f_10 = in.getSampleDataRO(INDEX2(k0+1,0, m_N0));
                    const register double* f_11 = in.getSampleDataRO(INDEX2(k0+1,1, m_N0));
                    const register double* f_01 = in.getSampleDataRO(INDEX2(k0,1, m_N0));
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = f_00[i]*cx0 + f_10[i]*cx7;
                        o[INDEX3(i,1,0,numComp,2)] = f_00[i]*cy1 + f_01[i]*cy6 + f_10[i]*cy3 + f_11[i]*cy4;
                        o[INDEX3(i,0,1,numComp,2)] = f_00[i]*cx0 + f_10[i]*cx7;
                        o[INDEX3(i,1,1,numComp,2)] = f_00[i]*cy3 + f_01[i]*cy4 + f_10[i]*cy1 + f_11[i]*cy6;
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of face 2 */
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < m_NE0; ++k0) {
                    const register double* f_11 = in.getSampleDataRO(INDEX2(k0+1,m_N1-1, m_N0));
                    const register double* f_01 = in.getSampleDataRO(INDEX2(k0,m_N1-1, m_N0));
                    const register double* f_10 = in.getSampleDataRO(INDEX2(k0+1,m_N1-2, m_N0));
                    const register double* f_00 = in.getSampleDataRO(INDEX2(k0,m_N1-2, m_N0));
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = f_01[i]*cx0 + f_11[i]*cx7;
                        o[INDEX3(i,1,0,numComp,2)] = f_00[i]*cy1 + f_01[i]*cy6 + f_10[i]*cy3 + f_11[i]*cy4;
                        o[INDEX3(i,0,1,numComp,2)] = f_01[i]*cx0 + f_11[i]*cx7;
                        o[INDEX3(i,1,1,numComp,2)] = f_00[i]*cy3 + f_01[i]*cy4 + f_10[i]*cy1 + f_11[i]*cy6;
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of face 3 */
            /* GENERATOR SNIP_GRAD_FACES BOTTOM */
        } // end of parallel section
    } else if (out.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
#pragma omp parallel
        {
            /*** GENERATOR SNIP_GRAD_REDUCED_FACES TOP */
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    const register double* f_10 = in.getSampleDataRO(INDEX2(1,k1, m_N0));
                    const register double* f_11 = in.getSampleDataRO(INDEX2(1,k1+1, m_N0));
                    const register double* f_01 = in.getSampleDataRO(INDEX2(0,k1+1, m_N0));
                    const register double* f_00 = in.getSampleDataRO(INDEX2(0,k1, m_N0));
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = cx5*(f_10[i] + f_11[i]) + cx2*(f_00[i] + f_01[i]);
                        o[INDEX3(i,1,0,numComp,2)] = f_00[i]*cy0 + f_01[i]*cy7;
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of face 0 */
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    const register double* f_10 = in.getSampleDataRO(INDEX2(m_N0-1,k1, m_N0));
                    const register double* f_11 = in.getSampleDataRO(INDEX2(m_N0-1,k1+1, m_N0));
                    const register double* f_01 = in.getSampleDataRO(INDEX2(m_N0-2,k1+1, m_N0));
                    const register double* f_00 = in.getSampleDataRO(INDEX2(m_N0-2,k1, m_N0));
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = cx5*(f_10[i] + f_11[i]) + cx2*(f_00[i] + f_01[i]);
                        o[INDEX3(i,1,0,numComp,2)] = f_10[i]*cy0 + f_11[i]*cy7;
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of face 1 */
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < m_NE0; ++k0) {
                    const register double* f_00 = in.getSampleDataRO(INDEX2(k0,0, m_N0));
                    const register double* f_10 = in.getSampleDataRO(INDEX2(k0+1,0, m_N0));
                    const register double* f_11 = in.getSampleDataRO(INDEX2(k0+1,1, m_N0));
                    const register double* f_01 = in.getSampleDataRO(INDEX2(k0,1, m_N0));
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = f_00[i]*cx0 + f_10[i]*cx7;
                        o[INDEX3(i,1,0,numComp,2)] = cy2*(f_00[i] + f_10[i]) + cy5*(f_01[i] + f_11[i]);
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of face 2 */
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < m_NE0; ++k0) {
                    const register double* f_11 = in.getSampleDataRO(INDEX2(k0+1,m_N1-1, m_N0));
                    const register double* f_01 = in.getSampleDataRO(INDEX2(k0,m_N1-1, m_N0));
                    const register double* f_10 = in.getSampleDataRO(INDEX2(k0+1,m_N1-2, m_N0));
                    const register double* f_00 = in.getSampleDataRO(INDEX2(k0,m_N1-2, m_N0));
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = f_01[i]*cx0 + f_11[i]*cx7;
                        o[INDEX3(i,1,0,numComp,2)] = cy5*(f_01[i] + f_11[i]) + cy2*(f_00[i] + f_10[i]);
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of face 3 */
            /* GENERATOR SNIP_GRAD_REDUCED_FACES BOTTOM */
        } // end of parallel section
    } else {
        stringstream msg;
        msg << "setToGradient() not implemented for "
            << functionSpaceTypeAsString(out.getFunctionSpace().getTypeCode());
        throw RipleyException(msg.str());
    }
}

void Rectangle::setToIntegrals(vector<double>& integrals, const escript::Data& arg) const
{
    escript::Data& in = *const_cast<escript::Data*>(&arg);
    const dim_t numComp = in.getDataPointSize();
    const double h0 = m_l0/m_gNE0;
    const double h1 = m_l1/m_gNE1;
    if (arg.getFunctionSpace().getTypeCode() == Elements) {
        const double w_0 = h0*h1/4.;
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
#pragma omp for nowait
            for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                    const double* f = in.getSampleDataRO(INDEX2(k0, k1, m_NE0));
                    for (index_t i=0; i < numComp; ++i) {
                        const register double f_0 = f[INDEX2(i,0,numComp)];
                        const register double f_1 = f[INDEX2(i,1,numComp)];
                        const register double f_2 = f[INDEX2(i,2,numComp)];
                        const register double f_3 = f[INDEX2(i,3,numComp)];
                        int_local[i]+=(f_0+f_1+f_2+f_3)*w_0;
                    }  /* end of component loop i */
                } /* end of k0 loop */
            } /* end of k1 loop */

#pragma omp critical
            for (index_t i=0; i<numComp; i++)
                integrals[i]+=int_local[i];
        } // end of parallel section
    } else if (arg.getFunctionSpace().getTypeCode() == ReducedElements) {
        const double w_0 = h0*h1;
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
#pragma omp for nowait
            for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                    const double* f = in.getSampleDataRO(INDEX2(k0, k1, m_NE0));
                    for (index_t i=0; i < numComp; ++i) {
                        int_local[i]+=f[i]*w_0;
                    }  /* end of component loop i */
                } /* end of k0 loop */
            } /* end of k1 loop */

#pragma omp critical
            for (index_t i=0; i<numComp; i++)
                integrals[i]+=int_local[i];
        } // end of parallel section
    } else if (arg.getFunctionSpace().getTypeCode() == FaceElements) {
        const double w_0 = h0/2.;
        const double w_1 = h1/2.;
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                    const double* f = in.getSampleDataRO(m_faceOffset[0]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        const register double f_0 = f[INDEX2(i,0,numComp)];
                        const register double f_1 = f[INDEX2(i,1,numComp)];
                        int_local[i]+=(f_0+f_1)*w_1;
                    }  /* end of component loop i */
                } /* end of k1 loop */
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                    const double* f = in.getSampleDataRO(m_faceOffset[1]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        const register double f_0 = f[INDEX2(i,0,numComp)];
                        const register double f_1 = f[INDEX2(i,1,numComp)];
                        int_local[i]+=(f_0+f_1)*w_1;
                    }  /* end of component loop i */
                } /* end of k1 loop */
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                    const double* f = in.getSampleDataRO(m_faceOffset[2]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        const register double f_0 = f[INDEX2(i,0,numComp)];
                        const register double f_1 = f[INDEX2(i,1,numComp)];
                        int_local[i]+=(f_0+f_1)*w_0;
                    }  /* end of component loop i */
                } /* end of k0 loop */
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                    const double* f = in.getSampleDataRO(m_faceOffset[3]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        const register double f_0 = f[INDEX2(i,0,numComp)];
                        const register double f_1 = f[INDEX2(i,1,numComp)];
                        int_local[i]+=(f_0+f_1)*w_0;
                    }  /* end of component loop i */
                } /* end of k0 loop */
            }

#pragma omp critical
            for (index_t i=0; i<numComp; i++)
                integrals[i]+=int_local[i];
        } // end of parallel section
    } else if (arg.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                    const double* f = in.getSampleDataRO(m_faceOffset[0]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        int_local[i]+=f[i]*h1;
                    }  /* end of component loop i */
                } /* end of k1 loop */
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                    const double* f = in.getSampleDataRO(m_faceOffset[1]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        int_local[i]+=f[i]*h1;
                    }  /* end of component loop i */
                } /* end of k1 loop */
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                    const double* f = in.getSampleDataRO(m_faceOffset[2]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        int_local[i]+=f[i]*h0;
                    }  /* end of component loop i */
                } /* end of k0 loop */
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                    const double* f = in.getSampleDataRO(m_faceOffset[3]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        int_local[i]+=f[i]*h0;
                    }  /* end of component loop i */
                } /* end of k0 loop */
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

void Rectangle::setToNormal(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == FaceElements) {
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE1; ++k1) {
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
                for (index_t k1 = 0; k1 < m_NE1; ++k1) {
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
                for (index_t k0 = 0; k0 < m_NE0; ++k0) {
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
                for (index_t k0 = 0; k0 < m_NE0; ++k0) {
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
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k1);
                    *o++ = -1.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k1);
                    *o++ = 1.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k0);
                    *o++ = 0.;
                    *o = -1.;
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k0);
                    *o++ = 0.;
                    *o = 1.;
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
    const index_t nDOF0 = (m_gNE0+1)/m_NX;
    const index_t nDOF1 = (m_gNE1+1)/m_NY;
    const int numDOF=nDOF0*nDOF1;
    const index_t left = (m_offset0==0 ? 0 : 1);
    const index_t bottom = (m_offset1==0 ? 0 : 1);
    vector<IndexVector> colIndices(numDOF); // for the couple blocks
    int numShared=0;

    m_dofMap.assign(getNumNodes(), 0);
#pragma omp parallel for
    for (index_t i=bottom; i<m_N1; i++) {
        for (index_t j=left; j<m_N0; j++) {
            m_dofMap[i*m_N0+j]=(i-bottom)*nDOF0+j-left;
        }
    }

    // corner node from bottom-left
    if (faces[0] == 0 && faces[2] == 0) {
        neighbour.push_back(m_mpiInfo->rank-m_NX-1);
        offsetInShared.push_back(offsetInShared.back()+1);
        sendShared.push_back(0);
        recvShared.push_back(numDOF+numShared);
        colIndices[0].push_back(numShared);
        m_dofMap[0]=numDOF+numShared;
        ++numShared;
    }
    // bottom edge
    if (faces[2] == 0) {
        neighbour.push_back(m_mpiInfo->rank-m_NX);
        offsetInShared.push_back(offsetInShared.back()+nDOF0);
        for (dim_t i=0; i<nDOF0; i++, numShared++) {
            sendShared.push_back(i);
            recvShared.push_back(numDOF+numShared);
            if (i>0)
                colIndices[i-1].push_back(numShared);
            colIndices[i].push_back(numShared);
            if (i<nDOF0-1)
                colIndices[i+1].push_back(numShared);
            m_dofMap[i+left]=numDOF+numShared;
        }
    }
    // corner node from bottom-right
    if (faces[1] == 0 && faces[2] == 0) {
        neighbour.push_back(m_mpiInfo->rank-m_NX+1);
        offsetInShared.push_back(offsetInShared.back()+1);
        sendShared.push_back(nDOF0-1);
        recvShared.push_back(numDOF+numShared);
        colIndices[nDOF0-2].push_back(numShared);
        colIndices[nDOF0-1].push_back(numShared);
        m_dofMap[m_N0-1]=numDOF+numShared;
        ++numShared;
    }
    // right edge
    if (faces[1] == 0) {
        neighbour.push_back(m_mpiInfo->rank+1);
        offsetInShared.push_back(offsetInShared.back()+nDOF1);
        for (dim_t i=0; i<nDOF1; i++, numShared++) {
            sendShared.push_back((i+1)*(nDOF0)-1);
            recvShared.push_back(numDOF+numShared);
            if (i>0)
                colIndices[i*(nDOF0)-1].push_back(numShared);
            colIndices[(i+1)*(nDOF0)-1].push_back(numShared);
            if (i<nDOF1-1)
                colIndices[(i+2)*(nDOF0)-1].push_back(numShared);
            m_dofMap[(i+bottom+1)*m_N0-1]=numDOF+numShared;
        }
    }
    // corner node from top-right
    if (faces[1] == 0 && faces[3] == 0) {
        neighbour.push_back(m_mpiInfo->rank+m_NX+1);
        offsetInShared.push_back(offsetInShared.back()+1);
        sendShared.push_back(numDOF-1);
        recvShared.push_back(numDOF+numShared);
        colIndices[numDOF-1].push_back(numShared);
        m_dofMap[m_N0*m_N1-1]=numDOF+numShared;
        ++numShared;
    }
    // top edge
    if (faces[3] == 0) {
        neighbour.push_back(m_mpiInfo->rank+m_NX);
        offsetInShared.push_back(offsetInShared.back()+nDOF0);
        for (dim_t i=0; i<nDOF0; i++, numShared++) {
            sendShared.push_back(numDOF-nDOF0+i);
            recvShared.push_back(numDOF+numShared);
            if (i>0)
                colIndices[numDOF-nDOF0+i-1].push_back(numShared);
            colIndices[numDOF-nDOF0+i].push_back(numShared);
            if (i<nDOF0-1)
                colIndices[numDOF-nDOF0+i+1].push_back(numShared);
            m_dofMap[m_N0*(m_N1-1)+i+left]=numDOF+numShared;
        }
    }
    // corner node from top-left
    if (faces[0] == 0 && faces[3] == 0) {
        neighbour.push_back(m_mpiInfo->rank+m_NX-1);
        offsetInShared.push_back(offsetInShared.back()+1);
        sendShared.push_back(numDOF-nDOF0);
        recvShared.push_back(numDOF+numShared);
        colIndices[numDOF-nDOF0].push_back(numShared);
        m_dofMap[m_N0*(m_N1-1)]=numDOF+numShared;
        ++numShared;
    }
    // left edge
    if (faces[0] == 0) {
        neighbour.push_back(m_mpiInfo->rank-1);
        offsetInShared.push_back(offsetInShared.back()+nDOF1);
        for (dim_t i=0; i<nDOF1; i++, numShared++) {
            sendShared.push_back(i*nDOF0);
            recvShared.push_back(numDOF+numShared);
            if (i>0)
                colIndices[(i-1)*nDOF0].push_back(numShared);
            colIndices[i*nDOF0].push_back(numShared);
            if (i<nDOF1-1)
                colIndices[(i+1)*nDOF0].push_back(numShared);
            m_dofMap[(i+bottom)*m_N0]=numDOF+numShared;
        }
    }

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
    */

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
        const int num=insertNeighbours(index, i);
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

    /*
    cout << "--- main_pattern ---" << endl;
    cout << "M=" << M << ", N=" << N << endl;
    for (size_t i=0; i<ptr.size(); i++) {
        cout << "ptr[" << i << "]=" << ptr[i] << endl;
    }
    for (size_t i=0; i<index.size(); i++) {
        cout << "index[" << i << "]=" << index[i] << endl;
    }
    */

    // column & row couple patterns
    ptr.assign(1, 0);
    index.clear();

    for (index_t i=0; i<numDOF; i++) {
        index.insert(index.end(), colIndices[i].begin(), colIndices[i].end());
        ptr.push_back(ptr.back()+colIndices[i].size());
    }

    // paso will manage the memory
    indexC = MEMALLOC(index.size(), index_t);
    ptrC = MEMALLOC(ptr.size(), index_t);
    copy(index.begin(), index.end(), indexC);
    copy(ptr.begin(), ptr.end(), ptrC);
    M=ptr.size()-1;
    N=numShared;
    Paso_Pattern *colCouplePattern=Paso_Pattern_alloc(MATRIX_FORMAT_DEFAULT,
            M, N, ptrC, indexC);

    /*
    cout << "--- colCouple_pattern ---" << endl;
    cout << "M=" << M << ", N=" << N << endl;
    for (size_t i=0; i<ptr.size(); i++) {
        cout << "ptr[" << i << "]=" << ptr[i] << endl;
    }
    for (size_t i=0; i<index.size(); i++) {
        cout << "index[" << i << "]=" << index[i] << endl;
    }
    */

    // now build the row couple pattern
    IndexVector ptr2(1,0);
    IndexVector index2;
    for (dim_t id=0; id<N; id++) {
        dim_t n=0;
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

    // paso will manage the memory
    indexC = MEMALLOC(index2.size(), index_t);
    ptrC = MEMALLOC(ptr2.size(), index_t);
    copy(index2.begin(), index2.end(), indexC);
    copy(ptr2.begin(), ptr2.end(), ptrC);
    Paso_Pattern *rowCouplePattern=Paso_Pattern_alloc(MATRIX_FORMAT_DEFAULT,
            N, M, ptrC, indexC);

    /*
    cout << "--- rowCouple_pattern ---" << endl;
    cout << "M=" << N << ", N=" << M << endl;
    for (size_t i=0; i<ptr2.size(); i++) {
        cout << "ptr[" << i << "]=" << ptr2[i] << endl;
    }
    for (size_t i=0; i<index2.size(); i++) {
        cout << "index[" << i << "]=" << index2[i] << endl;
    }
    */

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
dim_t Rectangle::getNumDOF() const
{
    return (m_gNE0+1)/m_NX*(m_gNE1+1)/m_NY;
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
    // identifiers are ordered from left to right, bottom to top globablly.

    // build node distribution vector first.
    // rank i owns m_nodeDistribution[i+1]-nodeDistribution[i] nodes
    m_nodeDistribution.assign(m_mpiInfo->size+1, 0);
    const dim_t numDOF=getNumDOF();
    for (dim_t k=1; k<m_mpiInfo->size; k++) {
        m_nodeDistribution[k]=k*numDOF;
    }
    m_nodeDistribution[m_mpiInfo->size]=getNumDataPointsGlobal();
    m_nodeId.resize(getNumNodes());
    m_dofId.resize(numDOF);
    m_elementId.resize(getNumElements());
    m_faceId.resize(getNumFaceElements());

#pragma omp parallel
    {
        // nodes
#pragma omp for nowait
        for (dim_t i1=0; i1<m_N1; i1++) {
            for (dim_t i0=0; i0<m_N0; i0++) {
                m_nodeId[i0+i1*m_N0] = (m_offset1+i1)*(m_gNE0+1)+m_offset0+i0;
            }
        }

        // degrees of freedom
#pragma omp for nowait
        for (dim_t k=0; k<numDOF; k++)
            m_dofId[k] = m_nodeDistribution[m_mpiInfo->rank]+k;

        // elements
#pragma omp for nowait
        for (dim_t k=0; k<getNumElements(); k++)
            m_elementId[k]=k;

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
    const IndexVector facesPerEdge = getNumFacesPerBoundary();
    const index_t LEFT=1, RIGHT=2, BOTTOM=10, TOP=20;
    const index_t faceTag[] = { LEFT, RIGHT, BOTTOM, TOP };
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
    updateTagsInUse(FaceElements);
}

//private
int Rectangle::insertNeighbours(IndexVector& index, index_t node) const
{
    const index_t myN0 = (m_gNE0+1)/m_NX;
    const index_t myN1 = (m_gNE1+1)/m_NY;
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

//protected
void Rectangle::interpolateNodesOnElements(escript::Data& out,
                                        escript::Data& in, bool reduced) const
{
    const dim_t numComp = in.getDataPointSize();
    if (reduced) {
        /*** GENERATOR SNIP_INTERPOLATE_REDUCED_ELEMENTS TOP */
        const double c0 = .25;
#pragma omp parallel for
        for (index_t k1=0; k1 < m_NE1; ++k1) {
            for (index_t k0=0; k0 < m_NE0; ++k0) {
                const register double* f_11 = in.getSampleDataRO(INDEX2(k0+1,k1+1, m_N0));
                const register double* f_10 = in.getSampleDataRO(INDEX2(k0+1,k1, m_N0));
                const register double* f_00 = in.getSampleDataRO(INDEX2(k0,k1, m_N0));
                const register double* f_01 = in.getSampleDataRO(INDEX2(k0,k1+1, m_N0));
                double* o = out.getSampleDataRW(INDEX2(k0,k1,m_NE0));
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = c0*(f_00[i] + f_01[i] + f_10[i] + f_11[i]);
                } /* end of component loop i */
            } /* end of k0 loop */
        } /* end of k1 loop */
        /* GENERATOR SNIP_INTERPOLATE_REDUCED_ELEMENTS BOTTOM */
    } else {
        /*** GENERATOR SNIP_INTERPOLATE_ELEMENTS TOP */
        const double c0 = .16666666666666666667;
        const double c1 = .044658198738520451079;
        const double c2 = .62200846792814621559;
#pragma omp parallel for
        for (index_t k1=0; k1 < m_NE1; ++k1) {
            for (index_t k0=0; k0 < m_NE0; ++k0) {
                const register double* f_10 = in.getSampleDataRO(INDEX2(k0+1,k1, m_N0));
                const register double* f_11 = in.getSampleDataRO(INDEX2(k0+1,k1+1, m_N0));
                const register double* f_01 = in.getSampleDataRO(INDEX2(k0,k1+1, m_N0));
                const register double* f_00 = in.getSampleDataRO(INDEX2(k0,k1, m_N0));
                double* o = out.getSampleDataRW(INDEX2(k0,k1,m_NE0));
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = f_00[i]*c2 + f_11[i]*c1 + c0*(f_01[i] + f_10[i]);
                    o[INDEX2(i,numComp,1)] = f_01[i]*c1 + f_10[i]*c2 + c0*(f_00[i] + f_11[i]);
                    o[INDEX2(i,numComp,2)] = f_01[i]*c2 + f_10[i]*c1 + c0*(f_00[i] + f_11[i]);
                    o[INDEX2(i,numComp,3)] = f_00[i]*c1 + f_11[i]*c2 + c0*(f_01[i] + f_10[i]);
                } /* end of component loop i */
            } /* end of k0 loop */
        } /* end of k1 loop */
        /* GENERATOR SNIP_INTERPOLATE_ELEMENTS BOTTOM */
    }
}

//protected
void Rectangle::interpolateNodesOnFaces(escript::Data& out, escript::Data& in,
                                        bool reduced) const
{
    const dim_t numComp = in.getDataPointSize();
    if (reduced) {
        const double c0 = .5;
#pragma omp parallel
        {
            /*** GENERATOR SNIP_INTERPOLATE_REDUCED_FACES TOP */
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    const register double* f_00 = in.getSampleDataRO(INDEX2(0,k1, m_N0));
                    const register double* f_01 = in.getSampleDataRO(INDEX2(0,k1+1, m_N0));
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = c0*(f_00[i] + f_01[i]);
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of face 0 */
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    const register double* f_10 = in.getSampleDataRO(INDEX2(m_N0-1,k1, m_N0));
                    const register double* f_11 = in.getSampleDataRO(INDEX2(m_N0-1,k1+1, m_N0));
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = c0*(f_10[i] + f_11[i]);
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of face 1 */
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < m_NE0; ++k0) {
                    const register double* f_10 = in.getSampleDataRO(INDEX2(k0+1,0, m_N0));
                    const register double* f_00 = in.getSampleDataRO(INDEX2(k0,0, m_N0));
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = c0*(f_00[i] + f_10[i]);
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of face 2 */
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < m_NE0; ++k0) {
                    const register double* f_01 = in.getSampleDataRO(INDEX2(k0,m_N1-1, m_N0));
                    const register double* f_11 = in.getSampleDataRO(INDEX2(k0+1,m_N1-1, m_N0));
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = c0*(f_01[i] + f_11[i]);
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of face 3 */
            /* GENERATOR SNIP_INTERPOLATE_REDUCED_FACES BOTTOM */
        } // end of parallel section
    } else {
        const double c0 = 0.21132486540518711775;
        const double c1 = 0.78867513459481288225;
#pragma omp parallel
        {
            /*** GENERATOR SNIP_INTERPOLATE_FACES TOP */
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    const register double* f_01 = in.getSampleDataRO(INDEX2(0,k1+1, m_N0));
                    const register double* f_00 = in.getSampleDataRO(INDEX2(0,k1, m_N0));
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = f_00[i]*c1 + f_01[i]*c0;
                        o[INDEX2(i,numComp,1)] = f_00[i]*c0 + f_01[i]*c1;
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of face 0 */
            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1=0; k1 < m_NE1; ++k1) {
                    const register double* f_10 = in.getSampleDataRO(INDEX2(m_N0-1,k1, m_N0));
                    const register double* f_11 = in.getSampleDataRO(INDEX2(m_N0-1,k1+1, m_N0));
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k1);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = f_10[i]*c1 + f_11[i]*c0;
                        o[INDEX2(i,numComp,1)] = f_10[i]*c0 + f_11[i]*c1;
                    } /* end of component loop i */
                } /* end of k1 loop */
            } /* end of face 1 */
            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < m_NE0; ++k0) {
                    const register double* f_10 = in.getSampleDataRO(INDEX2(k0+1,0, m_N0));
                    const register double* f_00 = in.getSampleDataRO(INDEX2(k0,0, m_N0));
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = f_00[i]*c1 + f_10[i]*c0;
                        o[INDEX2(i,numComp,1)] = f_00[i]*c0 + f_10[i]*c1;
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of face 2 */
            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0=0; k0 < m_NE0; ++k0) {
                    const register double* f_11 = in.getSampleDataRO(INDEX2(k0+1,m_N1-1, m_N0));
                    const register double* f_01 = in.getSampleDataRO(INDEX2(k0,m_N1-1, m_N0));
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k0);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX2(i,numComp,0)] = f_01[i]*c1 + f_11[i]*c0;
                        o[INDEX2(i,numComp,1)] = f_01[i]*c0 + f_11[i]*c1;
                    } /* end of component loop i */
                } /* end of k0 loop */
            } /* end of face 3 */
            /* GENERATOR SNIP_INTERPOLATE_FACES BOTTOM */
        } // end of parallel section
    }
}

//protected
void Rectangle::nodesToDOF(escript::Data& out, escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    out.requireWrite();

    const index_t left = (m_offset0==0 ? 0 : 1);
    const index_t bottom = (m_offset1==0 ? 0 : 1);
    const index_t right = (m_mpiInfo->rank%m_NX==m_NX-1 ? m_N0 : m_N0-1);
    const index_t top = (m_mpiInfo->rank/m_NX==m_NY-1 ? m_N1 : m_N1-1);
    index_t n=0;
    for (index_t i=bottom; i<top; i++) {
        for (index_t j=left; j<right; j++, n++) {
            const double* src=in.getSampleDataRO(j+i*m_N0);
            copy(src, src+numComp, out.getSampleDataRW(n));
        }
    }
}

//protected
void Rectangle::assemblePDESingle(Paso_SystemMatrix* mat,
        escript::Data& rhs, const escript::Data& A, const escript::Data& B,
        const escript::Data& C, const escript::Data& D,
        const escript::Data& X, const escript::Data& Y,
        const escript::Data& d, const escript::Data& y) const
{
    const double h0 = m_l0/m_gNE0;
    const double h1 = m_l1/m_gNE1;
    /* GENERATOR SNIP_PDE_SINGLE_PRE TOP */
    const double w0 = -0.1555021169820365539*h1/h0;
    const double w1 = 0.041666666666666666667;
    const double w10 = -0.041666666666666666667*h0/h1;
    const double w11 = 0.1555021169820365539*h1/h0;
    const double w12 = 0.1555021169820365539*h0/h1;
    const double w13 = 0.01116454968463011277*h0/h1;
    const double w14 = 0.01116454968463011277*h1/h0;
    const double w15 = 0.041666666666666666667*h1/h0;
    const double w16 = -0.01116454968463011277*h0/h1;
    const double w17 = -0.1555021169820365539*h0/h1;
    const double w18 = -0.33333333333333333333*h1/h0;
    const double w19 = 0.25000000000000000000;
    const double w2 = -0.15550211698203655390;
    const double w20 = -0.25000000000000000000;
    const double w21 = 0.16666666666666666667*h0/h1;
    const double w22 = -0.16666666666666666667*h1/h0;
    const double w23 = -0.16666666666666666667*h0/h1;
    const double w24 = 0.33333333333333333333*h1/h0;
    const double w25 = 0.33333333333333333333*h0/h1;
    const double w26 = 0.16666666666666666667*h1/h0;
    const double w27 = -0.33333333333333333333*h0/h1;
    const double w28 = -0.032861463941450536761*h1;
    const double w29 = -0.032861463941450536761*h0;
    const double w3 = 0.041666666666666666667*h0/h1;
    const double w30 = -0.12264065304058601714*h1;
    const double w31 = -0.0023593469594139828636*h1;
    const double w32 = -0.008805202725216129906*h0;
    const double w33 = -0.008805202725216129906*h1;
    const double w34 = 0.032861463941450536761*h1;
    const double w35 = 0.008805202725216129906*h1;
    const double w36 = 0.008805202725216129906*h0;
    const double w37 = 0.0023593469594139828636*h1;
    const double w38 = 0.12264065304058601714*h1;
    const double w39 = 0.032861463941450536761*h0;
    const double w4 = 0.15550211698203655390;
    const double w40 = -0.12264065304058601714*h0;
    const double w41 = -0.0023593469594139828636*h0;
    const double w42 = 0.0023593469594139828636*h0;
    const double w43 = 0.12264065304058601714*h0;
    const double w44 = -0.16666666666666666667*h1;
    const double w45 = -0.083333333333333333333*h0;
    const double w46 = 0.083333333333333333333*h1;
    const double w47 = 0.16666666666666666667*h1;
    const double w48 = 0.083333333333333333333*h0;
    const double w49 = -0.16666666666666666667*h0;
    const double w5 = -0.041666666666666666667;
    const double w50 = 0.16666666666666666667*h0;
    const double w51 = -0.083333333333333333333*h1;
    const double w52 = 0.025917019497006092316*h0*h1;
    const double w53 = 0.0018607582807716854616*h0*h1;
    const double w54 = 0.0069444444444444444444*h0*h1;
    const double w55 = 0.09672363354357992482*h0*h1;
    const double w56 = 0.00049858867864229740201*h0*h1;
    const double w57 = 0.055555555555555555556*h0*h1;
    const double w58 = 0.027777777777777777778*h0*h1;
    const double w59 = 0.11111111111111111111*h0*h1;
    const double w6 = -0.01116454968463011277*h1/h0;
    const double w60 = -0.19716878364870322056*h1;
    const double w61 = -0.19716878364870322056*h0;
    const double w62 = -0.052831216351296779436*h0;
    const double w63 = -0.052831216351296779436*h1;
    const double w64 = 0.19716878364870322056*h1;
    const double w65 = 0.052831216351296779436*h1;
    const double w66 = 0.19716878364870322056*h0;
    const double w67 = 0.052831216351296779436*h0;
    const double w68 = -0.5*h1;
    const double w69 = -0.5*h0;
    const double w7 = 0.011164549684630112770;
    const double w70 = 0.5*h1;
    const double w71 = 0.5*h0;
    const double w72 = 0.1555021169820365539*h0*h1;
    const double w73 = 0.041666666666666666667*h0*h1;
    const double w74 = 0.01116454968463011277*h0*h1;
    const double w75 = 0.25*h0*h1;
    const double w8 = -0.011164549684630112770;
    const double w9 = -0.041666666666666666667*h1/h0;
    /* GENERATOR SNIP_PDE_SINGLE_PRE BOTTOM */

    rhs.requireWrite();
#pragma omp parallel
    {
        for (index_t k1_0=0; k1_0<2; k1_0++) { // coloring
#pragma omp for
            for (index_t k1=k1_0; k1<m_NE1; k1+=2) {
                for (index_t k0=0; k0<m_NE0; ++k0)  {
                    bool add_EM_S=false;
                    bool add_EM_F=false;
                    vector<double> EM_S(4*4, 0);
                    vector<double> EM_F(4, 0);
                    const index_t e = k0 + m_NE0*k1;
                    /* GENERATOR SNIP_PDE_SINGLE TOP */
                    ///////////////
                    // process A //
                    ///////////////
                    if (!A.isEmpty()) {
                        add_EM_S=true;
                        const double* A_p=const_cast<escript::Data*>(&A)->getSampleDataRO(e);
                        if (A.actsExpanded()) {
                            const register double A_00_0 = A_p[INDEX3(0,0,0,2,2)];
                            const register double A_01_0 = A_p[INDEX3(0,1,0,2,2)];
                            const register double A_10_0 = A_p[INDEX3(1,0,0,2,2)];
                            const register double A_11_0 = A_p[INDEX3(1,1,0,2,2)];
                            const register double A_00_1 = A_p[INDEX3(0,0,1,2,2)];
                            const register double A_01_1 = A_p[INDEX3(0,1,1,2,2)];
                            const register double A_10_1 = A_p[INDEX3(1,0,1,2,2)];
                            const register double A_11_1 = A_p[INDEX3(1,1,1,2,2)];
                            const register double A_00_2 = A_p[INDEX3(0,0,2,2,2)];
                            const register double A_01_2 = A_p[INDEX3(0,1,2,2,2)];
                            const register double A_10_2 = A_p[INDEX3(1,0,2,2,2)];
                            const register double A_11_2 = A_p[INDEX3(1,1,2,2,2)];
                            const register double A_00_3 = A_p[INDEX3(0,0,3,2,2)];
                            const register double A_01_3 = A_p[INDEX3(0,1,3,2,2)];
                            const register double A_10_3 = A_p[INDEX3(1,0,3,2,2)];
                            const register double A_11_3 = A_p[INDEX3(1,1,3,2,2)];
                            const register double tmp4_0 = A_10_1 + A_10_2;
                            const register double tmp12_0 = A_11_0 + A_11_2;
                            const register double tmp2_0 = A_11_0 + A_11_1 + A_11_2 + A_11_3;
                            const register double tmp10_0 = A_01_3 + A_10_3;
                            const register double tmp14_0 = A_01_0 + A_01_3 + A_10_0 + A_10_3;
                            const register double tmp0_0 = A_01_0 + A_01_3;
                            const register double tmp13_0 = A_01_2 + A_10_1;
                            const register double tmp3_0 = A_00_2 + A_00_3;
                            const register double tmp11_0 = A_11_1 + A_11_3;
                            const register double tmp18_0 = A_01_1 + A_10_1;
                            const register double tmp1_0 = A_00_0 + A_00_1;
                            const register double tmp15_0 = A_01_1 + A_10_2;
                            const register double tmp5_0 = A_00_0 + A_00_1 + A_00_2 + A_00_3;
                            const register double tmp16_0 = A_10_0 + A_10_3;
                            const register double tmp6_0 = A_01_3 + A_10_0;
                            const register double tmp17_0 = A_01_1 + A_01_2;
                            const register double tmp9_0 = A_01_0 + A_10_0;
                            const register double tmp7_0 = A_01_0 + A_10_3;
                            const register double tmp8_0 = A_01_1 + A_01_2 + A_10_1 + A_10_2;
                            const register double tmp19_0 = A_01_2 + A_10_2;
                            const register double tmp14_1 = A_10_0*w8;
                            const register double tmp23_1 = tmp3_0*w14;
                            const register double tmp35_1 = A_01_0*w8;
                            const register double tmp54_1 = tmp13_0*w8;
                            const register double tmp20_1 = tmp9_0*w4;
                            const register double tmp25_1 = tmp12_0*w12;
                            const register double tmp2_1 = A_01_1*w4;
                            const register double tmp44_1 = tmp7_0*w7;
                            const register double tmp26_1 = tmp10_0*w4;
                            const register double tmp52_1 = tmp18_0*w8;
                            const register double tmp48_1 = A_10_1*w7;
                            const register double tmp46_1 = A_01_3*w8;
                            const register double tmp50_1 = A_01_0*w2;
                            const register double tmp8_1 = tmp4_0*w5;
                            const register double tmp56_1 = tmp19_0*w8;
                            const register double tmp9_1 = tmp2_0*w10;
                            const register double tmp19_1 = A_10_3*w2;
                            const register double tmp47_1 = A_10_2*w4;
                            const register double tmp16_1 = tmp3_0*w0;
                            const register double tmp18_1 = tmp1_0*w6;
                            const register double tmp31_1 = tmp11_0*w12;
                            const register double tmp55_1 = tmp15_0*w2;
                            const register double tmp39_1 = A_10_2*w7;
                            const register double tmp11_1 = tmp6_0*w7;
                            const register double tmp40_1 = tmp11_0*w17;
                            const register double tmp34_1 = tmp15_0*w8;
                            const register double tmp33_1 = tmp14_0*w5;
                            const register double tmp24_1 = tmp11_0*w13;
                            const register double tmp3_1 = tmp1_0*w0;
                            const register double tmp5_1 = tmp2_0*w3;
                            const register double tmp43_1 = tmp17_0*w5;
                            const register double tmp15_1 = A_01_2*w4;
                            const register double tmp53_1 = tmp19_0*w2;
                            const register double tmp27_1 = tmp3_0*w11;
                            const register double tmp32_1 = tmp13_0*w2;
                            const register double tmp10_1 = tmp5_0*w9;
                            const register double tmp37_1 = A_10_1*w4;
                            const register double tmp38_1 = tmp5_0*w15;
                            const register double tmp17_1 = A_01_1*w7;
                            const register double tmp12_1 = tmp7_0*w4;
                            const register double tmp22_1 = tmp10_0*w7;
                            const register double tmp57_1 = tmp18_0*w2;
                            const register double tmp28_1 = tmp9_0*w7;
                            const register double tmp29_1 = tmp1_0*w14;
                            const register double tmp51_1 = tmp11_0*w16;
                            const register double tmp42_1 = tmp12_0*w16;
                            const register double tmp49_1 = tmp12_0*w17;
                            const register double tmp21_1 = tmp1_0*w11;
                            const register double tmp1_1 = tmp0_0*w1;
                            const register double tmp45_1 = tmp6_0*w4;
                            const register double tmp7_1 = A_10_0*w2;
                            const register double tmp6_1 = tmp3_0*w6;
                            const register double tmp13_1 = tmp8_0*w1;
                            const register double tmp36_1 = tmp16_0*w1;
                            const register double tmp41_1 = A_01_3*w2;
                            const register double tmp30_1 = tmp12_0*w13;
                            const register double tmp4_1 = A_01_2*w7;
                            const register double tmp0_1 = A_10_3*w8;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1;
                            EM_S[INDEX2(1,2,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp9_1;
                            EM_S[INDEX2(3,2,4)]+=tmp14_1 + tmp15_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp1_1 + tmp5_1 + tmp8_1;
                            EM_S[INDEX2(0,0,4)]+=tmp13_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp24_1 + tmp25_1;
                            EM_S[INDEX2(3,3,4)]+=tmp13_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1;
                            EM_S[INDEX2(3,0,4)]+=tmp10_1 + tmp32_1 + tmp33_1 + tmp34_1 + tmp9_1;
                            EM_S[INDEX2(3,1,4)]+=tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp43_1;
                            EM_S[INDEX2(2,1,4)]+=tmp10_1 + tmp13_1 + tmp44_1 + tmp45_1 + tmp9_1;
                            EM_S[INDEX2(0,2,4)]+=tmp36_1 + tmp38_1 + tmp43_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                            EM_S[INDEX2(2,0,4)]+=tmp0_1 + tmp15_1 + tmp17_1 + tmp1_1 + tmp38_1 + tmp49_1 + tmp51_1 + tmp7_1 + tmp8_1;
                            EM_S[INDEX2(1,3,4)]+=tmp14_1 + tmp19_1 + tmp1_1 + tmp2_1 + tmp38_1 + tmp40_1 + tmp42_1 + tmp4_1 + tmp8_1;
                            EM_S[INDEX2(2,3,4)]+=tmp16_1 + tmp18_1 + tmp35_1 + tmp36_1 + tmp41_1 + tmp43_1 + tmp47_1 + tmp48_1 + tmp5_1;
                            EM_S[INDEX2(2,2,4)]+=tmp24_1 + tmp25_1 + tmp27_1 + tmp29_1 + tmp33_1 + tmp52_1 + tmp53_1;
                            EM_S[INDEX2(1,0,4)]+=tmp36_1 + tmp37_1 + tmp39_1 + tmp3_1 + tmp43_1 + tmp46_1 + tmp50_1 + tmp5_1 + tmp6_1;
                            EM_S[INDEX2(0,3,4)]+=tmp10_1 + tmp33_1 + tmp54_1 + tmp55_1 + tmp9_1;
                            EM_S[INDEX2(1,1,4)]+=tmp21_1 + tmp23_1 + tmp30_1 + tmp31_1 + tmp33_1 + tmp56_1 + tmp57_1;
                        } else { /* constant data */
                            const register double A_00 = A_p[INDEX2(0,0,2)];
                            const register double A_01 = A_p[INDEX2(0,1,2)];
                            const register double A_10 = A_p[INDEX2(1,0,2)];
                            const register double A_11 = A_p[INDEX2(1,1,2)];
                            const register double tmp0_0 = A_01 + A_10;
                            const register double tmp0_1 = A_00*w18;
                            const register double tmp10_1 = A_01*w20;
                            const register double tmp12_1 = A_00*w26;
                            const register double tmp4_1 = A_00*w22;
                            const register double tmp8_1 = A_00*w24;
                            const register double tmp13_1 = A_10*w19;
                            const register double tmp9_1 = tmp0_0*w20;
                            const register double tmp3_1 = A_11*w21;
                            const register double tmp11_1 = A_11*w27;
                            const register double tmp1_1 = A_01*w19;
                            const register double tmp6_1 = A_11*w23;
                            const register double tmp7_1 = A_11*w25;
                            const register double tmp2_1 = A_10*w20;
                            const register double tmp5_1 = tmp0_0*w19;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(1,2,4)]+=tmp4_1 + tmp5_1 + tmp6_1;
                            EM_S[INDEX2(3,2,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(0,0,4)]+=tmp5_1 + tmp7_1 + tmp8_1;
                            EM_S[INDEX2(3,3,4)]+=tmp5_1 + tmp7_1 + tmp8_1;
                            EM_S[INDEX2(3,0,4)]+=tmp4_1 + tmp6_1 + tmp9_1;
                            EM_S[INDEX2(3,1,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1;
                            EM_S[INDEX2(2,1,4)]+=tmp4_1 + tmp5_1 + tmp6_1;
                            EM_S[INDEX2(0,2,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1;
                            EM_S[INDEX2(2,0,4)]+=tmp11_1 + tmp12_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(1,3,4)]+=tmp11_1 + tmp12_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(2,3,4)]+=tmp0_1 + tmp10_1 + tmp13_1 + tmp3_1;
                            EM_S[INDEX2(2,2,4)]+=tmp7_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(1,0,4)]+=tmp0_1 + tmp10_1 + tmp13_1 + tmp3_1;
                            EM_S[INDEX2(0,3,4)]+=tmp4_1 + tmp6_1 + tmp9_1;
                            EM_S[INDEX2(1,1,4)]+=tmp7_1 + tmp8_1 + tmp9_1;
                        }
                    }
                    ///////////////
                    // process B //
                    ///////////////
                    if (!B.isEmpty()) {
                        add_EM_S=true;
                        const double* B_p=const_cast<escript::Data*>(&B)->getSampleDataRO(e);
                        if (B.actsExpanded()) {
                            const register double B_0_0 = B_p[INDEX2(0,0,2)];
                            const register double B_1_0 = B_p[INDEX2(1,0,2)];
                            const register double B_0_1 = B_p[INDEX2(0,1,2)];
                            const register double B_1_1 = B_p[INDEX2(1,1,2)];
                            const register double B_0_2 = B_p[INDEX2(0,2,2)];
                            const register double B_1_2 = B_p[INDEX2(1,2,2)];
                            const register double B_0_3 = B_p[INDEX2(0,3,2)];
                            const register double B_1_3 = B_p[INDEX2(1,3,2)];
                            const register double tmp3_0 = B_0_0 + B_0_2;
                            const register double tmp1_0 = B_1_2 + B_1_3;
                            const register double tmp2_0 = B_0_1 + B_0_3;
                            const register double tmp0_0 = B_1_0 + B_1_1;
                            const register double tmp63_1 = B_1_1*w42;
                            const register double tmp79_1 = B_1_1*w40;
                            const register double tmp37_1 = tmp3_0*w35;
                            const register double tmp8_1 = tmp0_0*w32;
                            const register double tmp71_1 = B_0_1*w34;
                            const register double tmp19_1 = B_0_3*w31;
                            const register double tmp15_1 = B_0_3*w34;
                            const register double tmp9_1 = tmp3_0*w34;
                            const register double tmp35_1 = B_1_0*w36;
                            const register double tmp66_1 = B_0_3*w28;
                            const register double tmp28_1 = B_1_0*w42;
                            const register double tmp22_1 = B_1_0*w40;
                            const register double tmp16_1 = B_1_2*w29;
                            const register double tmp6_1 = tmp2_0*w35;
                            const register double tmp55_1 = B_1_3*w40;
                            const register double tmp50_1 = B_1_3*w42;
                            const register double tmp7_1 = tmp1_0*w29;
                            const register double tmp1_1 = tmp1_0*w32;
                            const register double tmp57_1 = B_0_3*w30;
                            const register double tmp18_1 = B_1_1*w32;
                            const register double tmp53_1 = B_1_0*w41;
                            const register double tmp61_1 = B_1_3*w36;
                            const register double tmp27_1 = B_0_3*w38;
                            const register double tmp64_1 = B_0_2*w30;
                            const register double tmp76_1 = B_0_1*w38;
                            const register double tmp39_1 = tmp2_0*w34;
                            const register double tmp62_1 = B_0_1*w31;
                            const register double tmp56_1 = B_0_0*w31;
                            const register double tmp49_1 = B_1_1*w36;
                            const register double tmp2_1 = B_0_2*w31;
                            const register double tmp23_1 = B_0_2*w33;
                            const register double tmp38_1 = B_1_1*w43;
                            const register double tmp74_1 = B_1_2*w41;
                            const register double tmp43_1 = B_1_1*w41;
                            const register double tmp58_1 = B_0_2*w28;
                            const register double tmp67_1 = B_0_0*w33;
                            const register double tmp33_1 = tmp0_0*w39;
                            const register double tmp4_1 = B_0_0*w28;
                            const register double tmp20_1 = B_0_0*w30;
                            const register double tmp13_1 = B_0_2*w38;
                            const register double tmp65_1 = B_1_2*w43;
                            const register double tmp0_1 = tmp0_0*w29;
                            const register double tmp41_1 = tmp3_0*w33;
                            const register double tmp73_1 = B_0_2*w37;
                            const register double tmp69_1 = B_0_0*w38;
                            const register double tmp48_1 = B_1_2*w39;
                            const register double tmp59_1 = B_0_1*w33;
                            const register double tmp17_1 = B_1_3*w41;
                            const register double tmp5_1 = B_0_3*w33;
                            const register double tmp3_1 = B_0_1*w30;
                            const register double tmp21_1 = B_0_1*w28;
                            const register double tmp42_1 = B_1_0*w29;
                            const register double tmp54_1 = B_1_2*w32;
                            const register double tmp60_1 = B_1_0*w39;
                            const register double tmp32_1 = tmp1_0*w36;
                            const register double tmp10_1 = B_0_1*w37;
                            const register double tmp14_1 = B_0_0*w35;
                            const register double tmp29_1 = B_0_1*w35;
                            const register double tmp26_1 = B_1_2*w36;
                            const register double tmp30_1 = B_1_3*w43;
                            const register double tmp70_1 = B_0_2*w35;
                            const register double tmp34_1 = B_1_3*w39;
                            const register double tmp51_1 = B_1_0*w43;
                            const register double tmp31_1 = B_0_2*w34;
                            const register double tmp45_1 = tmp3_0*w28;
                            const register double tmp11_1 = tmp1_0*w39;
                            const register double tmp52_1 = B_1_1*w29;
                            const register double tmp44_1 = B_1_3*w32;
                            const register double tmp25_1 = B_1_1*w39;
                            const register double tmp47_1 = tmp2_0*w33;
                            const register double tmp72_1 = B_1_3*w29;
                            const register double tmp40_1 = tmp2_0*w28;
                            const register double tmp46_1 = B_1_2*w40;
                            const register double tmp36_1 = B_1_2*w42;
                            const register double tmp24_1 = B_0_0*w37;
                            const register double tmp77_1 = B_0_3*w35;
                            const register double tmp68_1 = B_0_3*w37;
                            const register double tmp78_1 = B_0_0*w34;
                            const register double tmp12_1 = tmp0_0*w36;
                            const register double tmp75_1 = B_1_0*w32;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(1,2,4)]+=tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(3,2,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1;
                            EM_S[INDEX2(0,0,4)]+=tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1;
                            EM_S[INDEX2(3,3,4)]+=tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1;
                            EM_S[INDEX2(3,0,4)]+=tmp32_1 + tmp33_1 + tmp6_1 + tmp9_1;
                            EM_S[INDEX2(3,1,4)]+=tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1;
                            EM_S[INDEX2(2,1,4)]+=tmp32_1 + tmp33_1 + tmp40_1 + tmp41_1;
                            EM_S[INDEX2(0,2,4)]+=tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1;
                            EM_S[INDEX2(2,0,4)]+=tmp45_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                            EM_S[INDEX2(1,3,4)]+=tmp37_1 + tmp39_1 + tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1;
                            EM_S[INDEX2(2,3,4)]+=tmp11_1 + tmp12_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1;
                            EM_S[INDEX2(2,2,4)]+=tmp60_1 + tmp61_1 + tmp62_1 + tmp63_1 + tmp64_1 + tmp65_1 + tmp66_1 + tmp67_1;
                            EM_S[INDEX2(1,0,4)]+=tmp0_1 + tmp1_1 + tmp68_1 + tmp69_1 + tmp70_1 + tmp71_1;
                            EM_S[INDEX2(0,3,4)]+=tmp40_1 + tmp41_1 + tmp7_1 + tmp8_1;
                            EM_S[INDEX2(1,1,4)]+=tmp72_1 + tmp73_1 + tmp74_1 + tmp75_1 + tmp76_1 + tmp77_1 + tmp78_1 + tmp79_1;
                        } else { /* constant data */
                            const register double B_0 = B_p[0];
                            const register double B_1 = B_p[1];
                            const register double tmp6_1 = B_1*w50;
                            const register double tmp1_1 = B_1*w45;
                            const register double tmp5_1 = B_1*w49;
                            const register double tmp4_1 = B_1*w48;
                            const register double tmp0_1 = B_0*w44;
                            const register double tmp2_1 = B_0*w46;
                            const register double tmp7_1 = B_0*w51;
                            const register double tmp3_1 = B_0*w47;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1;
                            EM_S[INDEX2(1,2,4)]+=tmp1_1 + tmp2_1;
                            EM_S[INDEX2(3,2,4)]+=tmp3_1 + tmp4_1;
                            EM_S[INDEX2(0,0,4)]+=tmp0_1 + tmp5_1;
                            EM_S[INDEX2(3,3,4)]+=tmp3_1 + tmp6_1;
                            EM_S[INDEX2(3,0,4)]+=tmp2_1 + tmp4_1;
                            EM_S[INDEX2(3,1,4)]+=tmp2_1 + tmp6_1;
                            EM_S[INDEX2(2,1,4)]+=tmp4_1 + tmp7_1;
                            EM_S[INDEX2(0,2,4)]+=tmp5_1 + tmp7_1;
                            EM_S[INDEX2(2,0,4)]+=tmp6_1 + tmp7_1;
                            EM_S[INDEX2(1,3,4)]+=tmp2_1 + tmp5_1;
                            EM_S[INDEX2(2,3,4)]+=tmp0_1 + tmp4_1;
                            EM_S[INDEX2(2,2,4)]+=tmp0_1 + tmp6_1;
                            EM_S[INDEX2(1,0,4)]+=tmp1_1 + tmp3_1;
                            EM_S[INDEX2(0,3,4)]+=tmp1_1 + tmp7_1;
                            EM_S[INDEX2(1,1,4)]+=tmp3_1 + tmp5_1;
                        }
                    }
                    ///////////////
                    // process C //
                    ///////////////
                    if (!C.isEmpty()) {
                        add_EM_S=true;
                        const double* C_p=const_cast<escript::Data*>(&C)->getSampleDataRO(e);
                        if (C.actsExpanded()) {
                            const register double C_0_0 = C_p[INDEX2(0,0,2)];
                            const register double C_1_0 = C_p[INDEX2(1,0,2)];
                            const register double C_0_1 = C_p[INDEX2(0,1,2)];
                            const register double C_1_1 = C_p[INDEX2(1,1,2)];
                            const register double C_0_2 = C_p[INDEX2(0,2,2)];
                            const register double C_1_2 = C_p[INDEX2(1,2,2)];
                            const register double C_0_3 = C_p[INDEX2(0,3,2)];
                            const register double C_1_3 = C_p[INDEX2(1,3,2)];
                            const register double tmp2_0 = C_0_1 + C_0_3;
                            const register double tmp1_0 = C_1_2 + C_1_3;
                            const register double tmp3_0 = C_0_0 + C_0_2;
                            const register double tmp0_0 = C_1_0 + C_1_1;
                            const register double tmp64_1 = C_0_2*w30;
                            const register double tmp14_1 = C_0_2*w28;
                            const register double tmp19_1 = C_0_3*w31;
                            const register double tmp22_1 = C_1_0*w40;
                            const register double tmp37_1 = tmp3_0*w35;
                            const register double tmp29_1 = C_0_1*w35;
                            const register double tmp73_1 = C_0_2*w37;
                            const register double tmp74_1 = C_1_2*w41;
                            const register double tmp52_1 = C_1_3*w39;
                            const register double tmp25_1 = C_1_1*w39;
                            const register double tmp62_1 = C_0_1*w31;
                            const register double tmp79_1 = C_1_1*w40;
                            const register double tmp43_1 = C_1_1*w36;
                            const register double tmp27_1 = C_0_3*w38;
                            const register double tmp28_1 = C_1_0*w42;
                            const register double tmp63_1 = C_1_1*w42;
                            const register double tmp59_1 = C_0_3*w34;
                            const register double tmp72_1 = C_1_3*w29;
                            const register double tmp40_1 = tmp2_0*w35;
                            const register double tmp13_1 = C_0_3*w30;
                            const register double tmp51_1 = C_1_2*w40;
                            const register double tmp54_1 = C_1_2*w42;
                            const register double tmp12_1 = C_0_0*w31;
                            const register double tmp2_1 = tmp1_0*w32;
                            const register double tmp68_1 = C_0_2*w31;
                            const register double tmp75_1 = C_1_0*w32;
                            const register double tmp49_1 = C_1_1*w41;
                            const register double tmp4_1 = C_0_2*w35;
                            const register double tmp66_1 = C_0_3*w28;
                            const register double tmp56_1 = C_0_1*w37;
                            const register double tmp5_1 = C_0_1*w34;
                            const register double tmp38_1 = tmp2_0*w34;
                            const register double tmp76_1 = C_0_1*w38;
                            const register double tmp21_1 = C_0_1*w28;
                            const register double tmp69_1 = C_0_1*w30;
                            const register double tmp53_1 = C_1_0*w36;
                            const register double tmp42_1 = C_1_2*w39;
                            const register double tmp32_1 = tmp1_0*w29;
                            const register double tmp45_1 = C_1_0*w43;
                            const register double tmp33_1 = tmp0_0*w32;
                            const register double tmp35_1 = C_1_0*w41;
                            const register double tmp26_1 = C_1_2*w36;
                            const register double tmp67_1 = C_0_0*w33;
                            const register double tmp31_1 = C_0_2*w34;
                            const register double tmp20_1 = C_0_0*w30;
                            const register double tmp70_1 = C_0_0*w28;
                            const register double tmp8_1 = tmp0_0*w39;
                            const register double tmp30_1 = C_1_3*w43;
                            const register double tmp0_1 = tmp0_0*w29;
                            const register double tmp17_1 = C_1_3*w41;
                            const register double tmp58_1 = C_0_0*w35;
                            const register double tmp9_1 = tmp3_0*w33;
                            const register double tmp61_1 = C_1_3*w36;
                            const register double tmp41_1 = tmp3_0*w34;
                            const register double tmp50_1 = C_1_3*w32;
                            const register double tmp18_1 = C_1_1*w32;
                            const register double tmp6_1 = tmp1_0*w36;
                            const register double tmp3_1 = C_0_0*w38;
                            const register double tmp34_1 = C_1_1*w29;
                            const register double tmp77_1 = C_0_3*w35;
                            const register double tmp65_1 = C_1_2*w43;
                            const register double tmp71_1 = C_0_3*w33;
                            const register double tmp55_1 = C_1_1*w43;
                            const register double tmp46_1 = tmp3_0*w28;
                            const register double tmp24_1 = C_0_0*w37;
                            const register double tmp10_1 = tmp1_0*w39;
                            const register double tmp48_1 = C_1_0*w29;
                            const register double tmp15_1 = C_0_1*w33;
                            const register double tmp36_1 = C_1_2*w32;
                            const register double tmp60_1 = C_1_0*w39;
                            const register double tmp47_1 = tmp2_0*w33;
                            const register double tmp16_1 = C_1_2*w29;
                            const register double tmp1_1 = C_0_3*w37;
                            const register double tmp7_1 = tmp2_0*w28;
                            const register double tmp39_1 = C_1_3*w40;
                            const register double tmp44_1 = C_1_3*w42;
                            const register double tmp57_1 = C_0_2*w38;
                            const register double tmp78_1 = C_0_0*w34;
                            const register double tmp11_1 = tmp0_0*w36;
                            const register double tmp23_1 = C_0_2*w33;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(1,2,4)]+=tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(3,2,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1;
                            EM_S[INDEX2(0,0,4)]+=tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1;
                            EM_S[INDEX2(3,3,4)]+=tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1;
                            EM_S[INDEX2(3,0,4)]+=tmp32_1 + tmp33_1 + tmp7_1 + tmp9_1;
                            EM_S[INDEX2(3,1,4)]+=tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1;
                            EM_S[INDEX2(2,1,4)]+=tmp32_1 + tmp33_1 + tmp40_1 + tmp41_1;
                            EM_S[INDEX2(0,2,4)]+=tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1;
                            EM_S[INDEX2(2,0,4)]+=tmp46_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                            EM_S[INDEX2(1,3,4)]+=tmp37_1 + tmp38_1 + tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1;
                            EM_S[INDEX2(2,3,4)]+=tmp10_1 + tmp11_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1;
                            EM_S[INDEX2(2,2,4)]+=tmp60_1 + tmp61_1 + tmp62_1 + tmp63_1 + tmp64_1 + tmp65_1 + tmp66_1 + tmp67_1;
                            EM_S[INDEX2(1,0,4)]+=tmp0_1 + tmp2_1 + tmp68_1 + tmp69_1 + tmp70_1 + tmp71_1;
                            EM_S[INDEX2(0,3,4)]+=tmp40_1 + tmp41_1 + tmp6_1 + tmp8_1;
                            EM_S[INDEX2(1,1,4)]+=tmp72_1 + tmp73_1 + tmp74_1 + tmp75_1 + tmp76_1 + tmp77_1 + tmp78_1 + tmp79_1;
                        } else { /* constant data */
                            const register double C_0 = C_p[0];
                            const register double C_1 = C_p[1];
                            const register double tmp1_1 = C_1*w45;
                            const register double tmp3_1 = C_0*w51;
                            const register double tmp4_1 = C_0*w44;
                            const register double tmp7_1 = C_0*w46;
                            const register double tmp5_1 = C_1*w49;
                            const register double tmp2_1 = C_1*w48;
                            const register double tmp0_1 = C_0*w47;
                            const register double tmp6_1 = C_1*w50;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1;
                            EM_S[INDEX2(1,2,4)]+=tmp2_1 + tmp3_1;
                            EM_S[INDEX2(3,2,4)]+=tmp2_1 + tmp4_1;
                            EM_S[INDEX2(0,0,4)]+=tmp4_1 + tmp5_1;
                            EM_S[INDEX2(3,3,4)]+=tmp0_1 + tmp6_1;
                            EM_S[INDEX2(3,0,4)]+=tmp1_1 + tmp3_1;
                            EM_S[INDEX2(3,1,4)]+=tmp5_1 + tmp7_1;
                            EM_S[INDEX2(2,1,4)]+=tmp1_1 + tmp7_1;
                            EM_S[INDEX2(0,2,4)]+=tmp3_1 + tmp6_1;
                            EM_S[INDEX2(2,0,4)]+=tmp3_1 + tmp5_1;
                            EM_S[INDEX2(1,3,4)]+=tmp6_1 + tmp7_1;
                            EM_S[INDEX2(2,3,4)]+=tmp0_1 + tmp2_1;
                            EM_S[INDEX2(2,2,4)]+=tmp4_1 + tmp6_1;
                            EM_S[INDEX2(1,0,4)]+=tmp1_1 + tmp4_1;
                            EM_S[INDEX2(0,3,4)]+=tmp2_1 + tmp7_1;
                            EM_S[INDEX2(1,1,4)]+=tmp0_1 + tmp5_1;
                        }
                    }
                    ///////////////
                    // process D //
                    ///////////////
                    if (!D.isEmpty()) {
                        add_EM_S=true;
                        const double* D_p=const_cast<escript::Data*>(&D)->getSampleDataRO(e);
                        if (D.actsExpanded()) {
                            const register double D_0 = D_p[0];
                            const register double D_1 = D_p[1];
                            const register double D_2 = D_p[2];
                            const register double D_3 = D_p[3];
                            const register double tmp4_0 = D_1 + D_3;
                            const register double tmp2_0 = D_0 + D_1 + D_2 + D_3;
                            const register double tmp5_0 = D_0 + D_2;
                            const register double tmp0_0 = D_0 + D_1;
                            const register double tmp6_0 = D_0 + D_3;
                            const register double tmp1_0 = D_2 + D_3;
                            const register double tmp3_0 = D_1 + D_2;
                            const register double tmp16_1 = D_1*w56;
                            const register double tmp14_1 = tmp6_0*w54;
                            const register double tmp8_1 = D_3*w55;
                            const register double tmp2_1 = tmp2_0*w54;
                            const register double tmp12_1 = tmp5_0*w52;
                            const register double tmp4_1 = tmp0_0*w53;
                            const register double tmp3_1 = tmp1_0*w52;
                            const register double tmp13_1 = tmp4_0*w53;
                            const register double tmp10_1 = tmp4_0*w52;
                            const register double tmp1_1 = tmp1_0*w53;
                            const register double tmp7_1 = D_3*w56;
                            const register double tmp0_1 = tmp0_0*w52;
                            const register double tmp11_1 = tmp5_0*w53;
                            const register double tmp9_1 = D_0*w56;
                            const register double tmp5_1 = tmp3_0*w54;
                            const register double tmp18_1 = D_2*w56;
                            const register double tmp17_1 = D_1*w55;
                            const register double tmp6_1 = D_0*w55;
                            const register double tmp15_1 = D_2*w55;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1;
                            EM_S[INDEX2(1,2,4)]+=tmp2_1;
                            EM_S[INDEX2(3,2,4)]+=tmp3_1 + tmp4_1;
                            EM_S[INDEX2(0,0,4)]+=tmp5_1 + tmp6_1 + tmp7_1;
                            EM_S[INDEX2(3,3,4)]+=tmp5_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(3,0,4)]+=tmp2_1;
                            EM_S[INDEX2(3,1,4)]+=tmp10_1 + tmp11_1;
                            EM_S[INDEX2(2,1,4)]+=tmp2_1;
                            EM_S[INDEX2(0,2,4)]+=tmp12_1 + tmp13_1;
                            EM_S[INDEX2(2,0,4)]+=tmp12_1 + tmp13_1;
                            EM_S[INDEX2(1,3,4)]+=tmp10_1 + tmp11_1;
                            EM_S[INDEX2(2,3,4)]+=tmp3_1 + tmp4_1;
                            EM_S[INDEX2(2,2,4)]+=tmp14_1 + tmp15_1 + tmp16_1;
                            EM_S[INDEX2(1,0,4)]+=tmp0_1 + tmp1_1;
                            EM_S[INDEX2(0,3,4)]+=tmp2_1;
                            EM_S[INDEX2(1,1,4)]+=tmp14_1 + tmp17_1 + tmp18_1;
                        } else { /* constant data */
                            const register double D_0 = D_p[0];
                            const register double tmp2_1 = D_0*w59;
                            const register double tmp1_1 = D_0*w58;
                            const register double tmp0_1 = D_0*w57;
                            EM_S[INDEX2(0,1,4)]+=tmp0_1;
                            EM_S[INDEX2(1,2,4)]+=tmp1_1;
                            EM_S[INDEX2(3,2,4)]+=tmp0_1;
                            EM_S[INDEX2(0,0,4)]+=tmp2_1;
                            EM_S[INDEX2(3,3,4)]+=tmp2_1;
                            EM_S[INDEX2(3,0,4)]+=tmp1_1;
                            EM_S[INDEX2(3,1,4)]+=tmp0_1;
                            EM_S[INDEX2(2,1,4)]+=tmp1_1;
                            EM_S[INDEX2(0,2,4)]+=tmp0_1;
                            EM_S[INDEX2(2,0,4)]+=tmp0_1;
                            EM_S[INDEX2(1,3,4)]+=tmp0_1;
                            EM_S[INDEX2(2,3,4)]+=tmp0_1;
                            EM_S[INDEX2(2,2,4)]+=tmp2_1;
                            EM_S[INDEX2(1,0,4)]+=tmp0_1;
                            EM_S[INDEX2(0,3,4)]+=tmp1_1;
                            EM_S[INDEX2(1,1,4)]+=tmp2_1;
                        }
                    }
                    ///////////////
                    // process X //
                    ///////////////
                    if (!X.isEmpty()) {
                        add_EM_F=true;
                        const double* X_p=const_cast<escript::Data*>(&X)->getSampleDataRO(e);
                        if (X.actsExpanded()) {
                            const register double X_0_0 = X_p[INDEX2(0,0,2)];
                            const register double X_1_0 = X_p[INDEX2(1,0,2)];
                            const register double X_0_1 = X_p[INDEX2(0,1,2)];
                            const register double X_1_1 = X_p[INDEX2(1,1,2)];
                            const register double X_0_2 = X_p[INDEX2(0,2,2)];
                            const register double X_1_2 = X_p[INDEX2(1,2,2)];
                            const register double X_0_3 = X_p[INDEX2(0,3,2)];
                            const register double X_1_3 = X_p[INDEX2(1,3,2)];
                            const register double tmp1_0 = X_1_1 + X_1_3;
                            const register double tmp3_0 = X_0_0 + X_0_1;
                            const register double tmp2_0 = X_1_0 + X_1_2;
                            const register double tmp0_0 = X_0_2 + X_0_3;
                            const register double tmp8_1 = tmp2_0*w66;
                            const register double tmp5_1 = tmp3_0*w64;
                            const register double tmp14_1 = tmp0_0*w64;
                            const register double tmp3_1 = tmp3_0*w60;
                            const register double tmp9_1 = tmp3_0*w63;
                            const register double tmp13_1 = tmp3_0*w65;
                            const register double tmp12_1 = tmp1_0*w66;
                            const register double tmp10_1 = tmp0_0*w60;
                            const register double tmp2_1 = tmp2_0*w61;
                            const register double tmp6_1 = tmp2_0*w62;
                            const register double tmp4_1 = tmp0_0*w65;
                            const register double tmp11_1 = tmp1_0*w67;
                            const register double tmp1_1 = tmp1_0*w62;
                            const register double tmp7_1 = tmp1_0*w61;
                            const register double tmp0_1 = tmp0_0*w63;
                            const register double tmp15_1 = tmp2_0*w67;
                            EM_F[0]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                            EM_F[1]+=tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1;
                            EM_F[2]+=tmp10_1 + tmp11_1 + tmp8_1 + tmp9_1;
                            EM_F[3]+=tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1;
                        } else { /* constant data */
                            const register double X_0 = X_p[0];
                            const register double X_1 = X_p[1];
                            const register double tmp3_1 = X_1*w71;
                            const register double tmp2_1 = X_0*w70;
                            const register double tmp1_1 = X_0*w68;
                            const register double tmp0_1 = X_1*w69;
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
                            const register double Y_0 = Y_p[0];
                            const register double Y_1 = Y_p[1];
                            const register double Y_2 = Y_p[2];
                            const register double Y_3 = Y_p[3];
                            const register double tmp0_0 = Y_1 + Y_2;
                            const register double tmp1_0 = Y_0 + Y_3;
                            const register double tmp9_1 = Y_0*w74;
                            const register double tmp4_1 = tmp1_0*w73;
                            const register double tmp5_1 = Y_2*w74;
                            const register double tmp7_1 = Y_1*w74;
                            const register double tmp6_1 = Y_2*w72;
                            const register double tmp2_1 = Y_3*w74;
                            const register double tmp8_1 = Y_3*w72;
                            const register double tmp3_1 = Y_1*w72;
                            const register double tmp0_1 = Y_0*w72;
                            const register double tmp1_1 = tmp0_0*w73;
                            EM_F[0]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_F[1]+=tmp3_1 + tmp4_1 + tmp5_1;
                            EM_F[2]+=tmp4_1 + tmp6_1 + tmp7_1;
                            EM_F[3]+=tmp1_1 + tmp8_1 + tmp9_1;
                        } else { /* constant data */
                            const register double Y_0 = Y_p[0];
                            const register double tmp0_1 = Y_0*w75;
                            EM_F[0]+=tmp0_1;
                            EM_F[1]+=tmp0_1;
                            EM_F[2]+=tmp0_1;
                            EM_F[3]+=tmp0_1;
                        }
                    }
                    /* GENERATOR SNIP_PDE_SINGLE BOTTOM */

                    // add to matrix (if add_EM_S) and RHS (if add_EM_F)
                    const index_t firstNode=m_N0*k1+k0;
                    IndexVector rowIndex;
                    rowIndex.push_back(m_dofMap[firstNode]);
                    rowIndex.push_back(m_dofMap[firstNode+1]);
                    rowIndex.push_back(m_dofMap[firstNode+m_N0]);
                    rowIndex.push_back(m_dofMap[firstNode+m_N0+1]);
                    if (add_EM_F) {
                        //cout << "-- AddtoRHS -- " << endl;
                        double *F_p=rhs.getSampleDataRW(0);
                        for (index_t i=0; i<4; i++) {
                            if (rowIndex[i]<getNumDOF()) {
                                F_p[rowIndex[i]]+=EM_F[i];
                                //cout << "F[" << rowIndex[i] << "]=" << F_p[rowIndex[i]] << endl;
                            }
                        }
                        //cout << "---"<<endl;
                    }
                    if (add_EM_S) {
                        //cout << "-- AddtoSystem -- " << endl;
                        addToSystemMatrix(mat, rowIndex, 1, rowIndex, 1, EM_S);
                    }

                } // end k0 loop
            } // end k1 loop
        } // end of coloring
    } // end of parallel region
}

void Rectangle::addToSystemMatrix(Paso_SystemMatrix* mat, 
       const IndexVector& nodes_Eq, dim_t num_Eq, const IndexVector& nodes_Sol,
       dim_t num_Sol, const vector<double>& array) const
{
    const dim_t numMyCols = mat->pattern->mainPattern->numInput;
    const dim_t numMyRows = mat->pattern->mainPattern->numOutput;

    const index_t* mainBlock_ptr = mat->mainBlock->pattern->ptr;
    const index_t* mainBlock_index = mat->mainBlock->pattern->index;
    double* mainBlock_val = mat->mainBlock->val;
    const index_t* col_coupleBlock_ptr = mat->col_coupleBlock->pattern->ptr;
    const index_t* col_coupleBlock_index = mat->col_coupleBlock->pattern->index;
    double* col_coupleBlock_val = mat->col_coupleBlock->val;
    const index_t* row_coupleBlock_ptr = mat->row_coupleBlock->pattern->ptr;
    const index_t* row_coupleBlock_index = mat->row_coupleBlock->pattern->index;
    double* row_coupleBlock_val = mat->row_coupleBlock->val;

    for (dim_t k_Eq = 0; k_Eq < nodes_Eq.size(); ++k_Eq) {
        // down columns of array
        const dim_t j_Eq = nodes_Eq[k_Eq];
        const dim_t i_row = j_Eq;
//printf("row:%d\n", i_row);
        // only look at the matrix rows stored on this processor
        if (i_row < numMyRows) {
            for (dim_t k_Sol = 0; k_Sol < nodes_Sol.size(); ++k_Sol) {
                const dim_t i_col = nodes_Sol[k_Sol];
                if (i_col < numMyCols) {
                    for (dim_t k = mainBlock_ptr[i_row]; k < mainBlock_ptr[i_row + 1]; ++k) {
                        if (mainBlock_index[k] == i_col) {
                            mainBlock_val[k] += array[INDEX2(k_Eq, k_Sol, nodes_Eq.size())];
                            break;
                        }
                    }
                } else {
                    for (dim_t k = col_coupleBlock_ptr[i_row]; k < col_coupleBlock_ptr[i_row + 1]; ++k) {
//cout << "col:" << i_col-numMyCols << " colIdx:" << col_coupleBlock_index[k] << endl;
                        if (col_coupleBlock_index[k] == i_col - numMyCols) {
                            col_coupleBlock_val[k] += array[INDEX2(k_Eq, k_Sol, nodes_Eq.size())];
                            break;
                        }
                    }
                }
            }
        } else {
            for (dim_t k_Sol = 0; k_Sol < nodes_Sol.size(); ++k_Sol) {
                // across rows of array
                const dim_t i_col = nodes_Sol[k_Sol];
                if (i_col < numMyCols) {
                    for (dim_t k = row_coupleBlock_ptr[i_row - numMyRows];
                         k < row_coupleBlock_ptr[i_row - numMyRows + 1]; ++k)
                    {
//cout << "col:" << i_col << " rowIdx:" << row_coupleBlock_index[k] << endl;
                        if (row_coupleBlock_index[k] == i_col) {
                            row_coupleBlock_val[k] += array[INDEX2(k_Eq, k_Sol, nodes_Eq.size())];
                            break;
                        }
                    }
                }
            }
        }
    }
}

} // end of namespace ripley

