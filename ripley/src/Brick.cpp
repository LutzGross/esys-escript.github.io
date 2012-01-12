
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
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

    if ((n0+1)%m_NX > 0 || (n1+1)%m_NY > 0 || (n2+1)%m_NZ > 0)
        throw RipleyException("Number of elements+1 must be separable into number of ranks in each dimension");

    if ((m_NX > 1 && (n0+1)/m_NX<2) || (m_NY > 1 && (n1+1)/m_NY<2) || (m_NZ > 1 && (n2+1)/m_NZ<2))
        throw RipleyException("Too few elements for the number of ranks");

    // local number of elements (including overlap)
    m_NE0 = m_ownNE0 = (m_NX>1 ? (n0+1)/m_NX : n0);
    if (m_mpiInfo->rank%m_NX>0 && m_mpiInfo->rank%m_NX<m_NX-1)
        m_NE0++;
    else if (m_NX>1 && m_mpiInfo->rank%m_NX==m_NX-1)
        m_ownNE0--;

    m_NE1 = m_ownNE1 = (m_NY>1 ? (n1+1)/m_NY : n1);
    if (m_mpiInfo->rank%(m_NX*m_NY)/m_NX>0 && m_mpiInfo->rank%(m_NX*m_NY)/m_NX<m_NY-1)
        m_NE1++;
    else if (m_NY>1 && m_mpiInfo->rank%(m_NX*m_NY)/m_NX==m_NY-1)
        m_ownNE1--;

    m_NE2 = m_ownNE2 = (m_NZ>1 ? (n2+1)/m_NZ : n2);
    if (m_mpiInfo->rank/(m_NX*m_NY)>0 && m_mpiInfo->rank/(m_NX*m_NY)<m_NZ-1)
        m_NE2++;
    else if (m_NZ>1 && m_mpiInfo->rank/(m_NX*m_NY)==m_NZ-1)
        m_ownNE2--;

    // local number of nodes
    m_N0 = m_NE0+1;
    m_N1 = m_NE1+1;
    m_N2 = m_NE2+1;

    // bottom-left-front node is at (offset0,offset1,offset2) in global mesh
    m_offset0 = (n0+1)/m_NX*(m_mpiInfo->rank%m_NX);
    if (m_offset0 > 0)
        m_offset0--;
    m_offset1 = (n1+1)/m_NY*(m_mpiInfo->rank%(m_NX*m_NY)/m_NX);
    if (m_offset1 > 0)
        m_offset1--;
    m_offset2 = (n2+1)/m_NZ*(m_mpiInfo->rank/(m_NX*m_NY));
    if (m_offset2 > 0)
        m_offset2--;

    populateSampleIds();
    createPattern();
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
    const Brick* o=dynamic_cast<const Brick*>(&other);
    if (o) {
        return (RipleyDomain::operator==(other) &&
                m_gNE0==o->m_gNE0 && m_gNE1==o->m_gNE1 && m_gNE2==o->m_gNE2
                && m_l0==o->m_l0 && m_l1==o->m_l1 && m_l2==o->m_l2
                && m_NX==o->m_NX && m_NY==o->m_NY && m_NZ==o->m_NZ);
    }

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
        << fsType;
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
                const index_t x=id%m_NE0 + 1;
                const index_t y=id%(m_NE0*m_NE1)/m_NE0 + 1;
                const index_t z=id/(m_NE0*m_NE1) + 1;
                return (m_dofMap[x + m_N0*y +m_N0*m_N1*z] < getNumDOF());
            }
        case FaceElements:
        case ReducedFaceElements:
            {
                // check ownership of face element's last node
                const IndexVector faces = getNumFacesPerBoundary();
                dim_t n=0;
                for (size_t i=0; i<faces.size(); i++) {
                    n+=faces[i];
                    if (id<n) {
                        const index_t j=id-n+faces[i];
                        if (i>=4) { // front or back
                            const index_t first=(i==4 ? 0 : m_N0*m_N1*(m_N2-1));
                            return (m_dofMap[first+j%m_NE0+1+(j/m_NE0+1)*m_N0] < getNumDOF());
                        } else if (i>=2) { // bottom or top
                            const index_t first=(i==2 ? 0 : m_N0*(m_N1-1));
                            return (m_dofMap[first+j%m_NE0+1+(j/m_NE0+1)*m_N0*m_N1] < getNumDOF());
                        } else { // left or right
                            const index_t first=(i==0 ? 0 : m_N0-1);
                            return (m_dofMap[first+(j%m_NE1+1)*m_N0+(j/m_NE1+1)*m_N0*m_N1] < getNumDOF());
                        }
                    }
                }
                return false;
            }
        default:
            break;
    }

    stringstream msg;
    msg << "ownSample() not implemented for "
        << functionSpaceTypeAsString(fsType);
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
        out.requireWrite();
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

void Brick::setToSize(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == Elements
            || out.getFunctionSpace().getTypeCode() == ReducedElements) {
        out.requireWrite();
        const dim_t numQuad=out.getNumDataPointsPerSample();
        const double xSize=getFirstCoordAndSpacing(0).second;
        const double ySize=getFirstCoordAndSpacing(1).second;
        const double zSize=getFirstCoordAndSpacing(2).second;
        const double size=min(min(xSize,ySize),zSize);
#pragma omp parallel for
        for (index_t k = 0; k < getNumElements(); ++k) {
            double* o = out.getSampleDataRW(k);
            fill(o, o+numQuad, size);
        }
    } else if (out.getFunctionSpace().getTypeCode() == FaceElements
            || out.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        out.requireWrite();
        const dim_t numQuad=out.getNumDataPointsPerSample();
        const double xSize=getFirstCoordAndSpacing(0).second;
        const double ySize=getFirstCoordAndSpacing(1).second;
        const double zSize=getFirstCoordAndSpacing(2).second;
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
                const double size=min(ySize,zSize);
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE2; ++k2) {
                    for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                        double* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE1));
                        fill(o, o+numQuad, size);
                    }
                }
            }

            if (m_faceOffset[1] > -1) {
                const double size=min(ySize,zSize);
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE2; ++k2) {
                    for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                        double* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE1));
                        fill(o, o+numQuad, size);
                    }
                }
            }

            if (m_faceOffset[2] > -1) {
                const double size=min(xSize,zSize);
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE2; ++k2) {
                    for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE0));
                        fill(o, o+numQuad, size);
                    }
                }
            }

            if (m_faceOffset[3] > -1) {
                const double size=min(xSize,zSize);
#pragma omp for nowait
                for (index_t k2 = 0; k2 < m_NE2; ++k2) {
                    for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE0));
                        fill(o, o+numQuad, size);
                    }
                }
            }

            if (m_faceOffset[4] > -1) {
                const double size=min(xSize,ySize);
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                    for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE0));
                        fill(o, o+numQuad, size);
                    }
                }
            }

            if (m_faceOffset[5] > -1) {
                const double size=min(xSize,ySize);
#pragma omp for nowait
                for (index_t k1 = 0; k1 < m_NE1; ++k1) {
                    for (index_t k0 = 0; k0 < m_NE0; ++k0) {
                        double* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE0));
                        fill(o, o+numQuad, size);
                    }
                }
            }
        } // end of parallel section

    } else {
        stringstream msg;
        msg << "setToSize() not implemented for "
            << functionSpaceTypeAsString(out.getFunctionSpace().getTypeCode());
        throw RipleyException(msg.str());
    }
}

Paso_SystemMatrixPattern* Brick::getPattern(bool reducedRowOrder,
                                            bool reducedColOrder) const
{
    if (reducedRowOrder || reducedColOrder)
        throw RipleyException("getPattern() not implemented for reduced order");

    return m_pattern;
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

IndexVector Brick::getNumSubdivisionsPerDim() const
{
    IndexVector ret;
    ret.push_back(m_NX);
    ret.push_back(m_NY);
    ret.push_back(m_NZ);
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
dim_t Brick::getNumDOF() const
{
    return (m_gNE0+1)/m_NX*(m_gNE1+1)/m_NY*(m_gNE2+1)/m_NZ;
}

//protected
dim_t Brick::getNumFaceElements() const
{
    const IndexVector faces = getNumFacesPerBoundary();
    dim_t n=0;
    for (size_t i=0; i<faces.size(); i++)
        n+=faces[i];
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

//protected
void Brick::assembleGradient(escript::Data& out, escript::Data& in) const
{
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
        out.requireWrite();
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
                    } // end of component loop i
                } // end of k0 loop
            } // end of k1 loop
        } // end of k2 loop
    } else if (out.getFunctionSpace().getTypeCode() == ReducedElements) {
        out.requireWrite();
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
                    } // end of component loop i
                } // end of k0 loop
            } // end of k1 loop
        } // end of k2 loop
    } else if (out.getFunctionSpace().getTypeCode() == FaceElements) {
        out.requireWrite();
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
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 0
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
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 1
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
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 2
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
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 3
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
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of face 4
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
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of face 5
        } // end of parallel section
    } else if (out.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        out.requireWrite();
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
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 0
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
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 1
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
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 2
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
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 3
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
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of face 4
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
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of face 5
        } // end of parallel section
    }
}

//protected
void Brick::assembleIntegrate(vector<double>& integrals, escript::Data& arg) const
{
    const dim_t numComp = arg.getDataPointSize();
    const double h0 = m_l0/m_gNE0;
    const double h1 = m_l1/m_gNE1;
    const double h2 = m_l2/m_gNE2;
    const index_t left = (m_offset0==0 ? 0 : 1);
    const index_t bottom = (m_offset1==0 ? 0 : 1);
    const index_t front = (m_offset2==0 ? 0 : 1);
    if (arg.getFunctionSpace().getTypeCode() == Elements) {
        const double w_0 = h0*h1*h2/8.;
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
#pragma omp for nowait
            for (index_t k2 = front; k2 < front+m_ownNE2; ++k2) {
                for (index_t k1 = bottom; k1 < bottom+m_ownNE1; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE0; ++k0) {
                        const double* f = arg.getSampleDataRO(INDEX3(k0, k1, k2, m_NE0, m_NE1));
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
                        }  // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of k2 loop

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
            for (index_t k2 = front; k2 < front+m_ownNE2; ++k2) {
                for (index_t k1 = bottom; k1 < bottom+m_ownNE1; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE0; ++k0) {
                        const double* f = arg.getSampleDataRO(INDEX3(k0, k1, k2, m_NE0, m_NE1));
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
    } else if (arg.getFunctionSpace().getTypeCode() == FaceElements) {
        const double w_0 = h1*h2/4.;
        const double w_1 = h0*h2/4.;
        const double w_2 = h0*h1/4.;
#pragma omp parallel
        {
            vector<double> int_local(numComp, 0);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE2; ++k2) {
                    for (index_t k1 = bottom; k1 < bottom+m_ownNE1; ++k1) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[0]+INDEX2(k1,k2,m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            const register double f_0 = f[INDEX2(i,0,numComp)];
                            const register double f_1 = f[INDEX2(i,1,numComp)];
                            const register double f_2 = f[INDEX2(i,2,numComp)];
                            const register double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_0;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE2; ++k2) {
                    for (index_t k1 = bottom; k1 < bottom+m_ownNE1; ++k1) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[1]+INDEX2(k1,k2,m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            const register double f_0 = f[INDEX2(i,0,numComp)];
                            const register double f_1 = f[INDEX2(i,1,numComp)];
                            const register double f_2 = f[INDEX2(i,2,numComp)];
                            const register double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_0;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE2; ++k2) {
                    for (index_t k0 = left; k0 < left+m_ownNE0; ++k0) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[2]+INDEX2(k0,k2,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            const register double f_0 = f[INDEX2(i,0,numComp)];
                            const register double f_1 = f[INDEX2(i,1,numComp)];
                            const register double f_2 = f[INDEX2(i,2,numComp)];
                            const register double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_1;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE2; ++k2) {
                    for (index_t k0 = left; k0 < left+m_ownNE0; ++k0) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[3]+INDEX2(k0,k2,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            const register double f_0 = f[INDEX2(i,0,numComp)];
                            const register double f_1 = f[INDEX2(i,1,numComp)];
                            const register double f_2 = f[INDEX2(i,2,numComp)];
                            const register double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_1;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE1; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE0; ++k0) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[4]+INDEX2(k0,k1,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            const register double f_0 = f[INDEX2(i,0,numComp)];
                            const register double f_1 = f[INDEX2(i,1,numComp)];
                            const register double f_2 = f[INDEX2(i,2,numComp)];
                            const register double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_2;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE1; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE0; ++k0) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[5]+INDEX2(k0,k1,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            const register double f_0 = f[INDEX2(i,0,numComp)];
                            const register double f_1 = f[INDEX2(i,1,numComp)];
                            const register double f_2 = f[INDEX2(i,2,numComp)];
                            const register double f_3 = f[INDEX2(i,3,numComp)];
                            int_local[i]+=(f_0+f_1+f_2+f_3)*w_2;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
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
                for (index_t k2 = front; k2 < front+m_ownNE2; ++k2) {
                    for (index_t k1 = bottom; k1 < bottom+m_ownNE1; ++k1) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[0]+INDEX2(k1,k2,m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_0;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE2; ++k2) {
                    for (index_t k1 = bottom; k1 < bottom+m_ownNE1; ++k1) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[1]+INDEX2(k1,k2,m_NE1));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_0;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE2; ++k2) {
                    for (index_t k0 = left; k0 < left+m_ownNE0; ++k0) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[2]+INDEX2(k0,k2,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_1;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k2 = front; k2 < front+m_ownNE2; ++k2) {
                    for (index_t k0 = left; k0 < left+m_ownNE0; ++k0) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[3]+INDEX2(k0,k2,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_1;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE1; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE0; ++k0) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[4]+INDEX2(k0,k1,m_NE0));
                        for (index_t i=0; i < numComp; ++i) {
                            int_local[i]+=f[i]*w_2;
                        }  // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            }

            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k1 = bottom; k1 < bottom+m_ownNE1; ++k1) {
                    for (index_t k0 = left; k0 < left+m_ownNE0; ++k0) {
                        const double* f = arg.getSampleDataRO(m_faceOffset[5]+INDEX2(k0,k1,m_NE0));
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

    }
}

//protected
dim_t Brick::insertNeighbourNodes(IndexVector& index, index_t node) const
{
    const dim_t nDOF0 = (m_gNE0+1)/m_NX;
    const dim_t nDOF1 = (m_gNE1+1)/m_NY;
    const dim_t nDOF2 = (m_gNE2+1)/m_NZ;
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
void Brick::nodesToDOF(escript::Data& out, escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    out.requireWrite();

    const index_t left = (m_offset0==0 ? 0 : 1);
    const index_t bottom = (m_offset1==0 ? 0 : 1);
    const index_t front = (m_offset2==0 ? 0 : 1);
    const dim_t nDOF0 = (m_gNE0+1)/m_NX;
    const dim_t nDOF1 = (m_gNE1+1)/m_NY;
    const dim_t nDOF2 = (m_gNE2+1)/m_NZ;
#pragma omp parallel for
    for (index_t i=0; i<nDOF2; i++) {
        for (index_t j=0; j<nDOF1; j++) {
            for (index_t k=0; k<nDOF0; k++) {
                const index_t n=k+left+(j+bottom)*m_N0+(i+front)*m_N0*m_N1;
                const double* src=in.getSampleDataRO(n);
                copy(src, src+numComp, out.getSampleDataRW(k+j*nDOF0+i*nDOF0*nDOF1));
            }
        }
    }
}

//protected
void Brick::dofToNodes(escript::Data& out, escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    Paso_Coupler* coupler = Paso_Coupler_alloc(m_connector, numComp);
    in.requireWrite();
    Paso_Coupler_startCollect(coupler, in.getSampleDataRW(0));

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
}

//private
void Brick::populateSampleIds()
{
    // identifiers are ordered from left to right, bottom to top, front to back
    // globally

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
#pragma omp for nowait
        // nodes
        for (dim_t i2=0; i2<m_N2; i2++) {
            for (dim_t i1=0; i1<m_N1; i1++) {
                for (dim_t i0=0; i0<m_N0; i0++) {
                    m_nodeId[i0+i1*m_N0+i2*m_N0*m_N1] =
                        (m_offset2+i2)*(m_gNE0+1)*(m_gNE1+1)
                        +(m_offset1+i1)*(m_gNE0+1)
                        +m_offset0+i0;
                }
            }
        }

        // degrees of freedom
#pragma omp for nowait
        for (dim_t k=0; k<numDOF; k++)
            m_dofId[k] = m_nodeDistribution[m_mpiInfo->rank]+k;

        // elements
#pragma omp for nowait
        for (dim_t i2=0; i2<m_NE2; i2++) {
            for (dim_t i1=0; i1<m_NE1; i1++) {
                for (dim_t i0=0; i0<m_NE0; i0++) {
                    m_elementId[i0+i1*m_NE0+i2*m_NE0*m_NE1] =
                        (m_offset2+i2)*m_gNE0*m_gNE1
                        +(m_offset1+i1)*m_gNE0
                        +m_offset0+i0;
                }
            }
        }

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

//private
void Brick::createPattern()
{
    const dim_t nDOF0 = (m_gNE0+1)/m_NX;
    const dim_t nDOF1 = (m_gNE1+1)/m_NY;
    const dim_t nDOF2 = (m_gNE2+1)/m_NZ;
    const index_t left = (m_offset0==0 ? 0 : 1);
    const index_t bottom = (m_offset1==0 ? 0 : 1);
    const index_t front = (m_offset2==0 ? 0 : 1);

    // populate node->DOF mapping with own degrees of freedom.
    // The rest is assigned in the loop further down
    m_dofMap.assign(getNumNodes(), 0);
#pragma omp parallel for
    for (index_t i=front; i<front+nDOF2; i++) {
        for (index_t j=bottom; j<bottom+nDOF1; j++) {
            for (index_t k=left; k<left+nDOF0; k++) {
                m_dofMap[i*m_N0*m_N1+j*m_N0+k]=(i-front)*nDOF0*nDOF1+(j-bottom)*nDOF0+k-left;
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
    const int x=m_mpiInfo->rank%m_NX;
    const int y=m_mpiInfo->rank%(m_NX*m_NY)/m_NX;
    const int z=m_mpiInfo->rank/(m_NX*m_NY);
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
                if (nx>=0 && ny>=0 && nz>=0 && nx<m_NX && ny<m_NY && nz<m_NZ) {
                    neighbour.push_back(nz*m_NX*m_NY+ny*m_NX+nx);
                    if (i0==0 && i1==0) {
                        // sharing front or back plane
                        offsetInShared.push_back(offsetInShared.back()+nDOF0*nDOF1);
                        for (dim_t i=0; i<nDOF1; i++) {
                            const int firstDOF=(i2==-1 ? i*nDOF0
                                    : i*nDOF0 + nDOF0*nDOF1*(nDOF2-1));
                            const int firstNode=(i2==-1 ? left+(i+bottom)*m_N0
                                    : left+(i+bottom)*m_N0+m_N0*m_N1*(m_N2-1));
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
                                    left+(i+front)*m_N0*m_N1
                                    : left+m_N0*((i+1+front)*m_N1-1));
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
                                    (bottom+(i+front)*m_N1)*m_N0
                                    : (bottom+1+(i+front)*m_N1)*m_N0-1);
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
                                m_dofMap[firstNode+j*m_N0]=numDOF+numShared;
                            }
                        }
                    } else if (i0==0) {
                        // sharing an edge in x direction
                        offsetInShared.push_back(offsetInShared.back()+nDOF0);
                        const int firstDOF=(i1+1)/2*nDOF0*(nDOF1-1)
                                           +(i2+1)/2*nDOF0*nDOF1*(nDOF2-1);
                        const int firstNode=(i1+1)/2*m_N0*(m_N1-1)
                                            +(i2+1)/2*m_N0*m_N1*(m_N2-1);
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
                        const int firstNode=(i0+1)/2*(m_N0-1)
                                            +(i2+1)/2*m_N0*m_N1*(m_N2-1);
                        for (dim_t i=0; i<nDOF1; i++, numShared++) {
                            sendShared.push_back(firstDOF+i*nDOF0);
                            recvShared.push_back(numDOF+numShared);
                            if (i>0)
                                colIndices[firstDOF+(i-1)*nDOF0].push_back(numShared);
                            colIndices[firstDOF+i*nDOF0].push_back(numShared);
                            if (i<nDOF1-1)
                                colIndices[firstDOF+(i+1)*nDOF0].push_back(numShared);
                            m_dofMap[firstNode+i*m_N0]=numDOF+numShared;
                        }
                    } else if (i2==0) {
                        // sharing an edge in z direction
                        offsetInShared.push_back(offsetInShared.back()+nDOF2);
                        const int firstDOF=(i0+1)/2*(nDOF0-1)
                                           +(i1+1)/2*nDOF0*(nDOF1-1);
                        const int firstNode=(i0+1)/2*(m_N0-1)
                                            +(i1+1)/2*m_N0*(m_N1-1);
                        for (dim_t i=0; i<nDOF2; i++, numShared++) {
                            sendShared.push_back(firstDOF+i*nDOF0*nDOF1);
                            recvShared.push_back(numDOF+numShared);
                            if (i>0)
                                colIndices[firstDOF+(i-1)*nDOF0*nDOF1].push_back(numShared);
                            colIndices[firstDOF+i*nDOF0*nDOF1].push_back(numShared);
                            if (i<nDOF2-1)
                                colIndices[firstDOF+(i+1)*nDOF0*nDOF1].push_back(numShared);
                            m_dofMap[firstNode+i*m_N0*m_N1]=numDOF+numShared;
                        }
                    } else {
                        // sharing a node
                        const int dof=(i0+1)/2*(nDOF0-1)
                                      +(i1+1)/2*nDOF0*(nDOF1-1)
                                      +(i2+1)/2*nDOF0*nDOF1*(nDOF2-1);
                        const int node=(i0+1)/2*(m_N0-1)
                                       +(i1+1)/2*m_N0*(m_N1-1)
                                       +(i2+1)/2*m_N0*m_N1*(m_N2-1);
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

//protected
void Brick::interpolateNodesOnElements(escript::Data& out, escript::Data& in,
                                       bool reduced) const
{
    const dim_t numComp = in.getDataPointSize();
    if (reduced) {
        out.requireWrite();
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
                    } // end of component loop i
                } // end of k0 loop
            } // end of k1 loop
        } // end of k2 loop
    } else {
        out.requireWrite();
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
                    } // end of component loop i
                } // end of k0 loop
            } // end of k1 loop
        } // end of k2 loop
    }
}

//protected
void Brick::interpolateNodesOnFaces(escript::Data& out, escript::Data& in,
                                    bool reduced) const
{
    const dim_t numComp = in.getDataPointSize();
    if (reduced) {
        out.requireWrite();
        const double c0 = .25;
#pragma omp parallel
        {
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
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 0
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
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 1
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
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 2
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
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 3
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
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of face 4
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
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 0
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
                        } // end of component loop i
                    } // end of k1 loop
                } // end of k2 loop
            } // end of face 1
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
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 2
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
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k2 loop
            } // end of face 3
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
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of face 4
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
                        } // end of component loop i
                    } // end of k0 loop
                } // end of k1 loop
            } // end of face 5
        } // end of parallel section
    }
}

//protected
void Brick::assemblePDESingle(Paso_SystemMatrix* mat, escript::Data& rhs,
        const escript::Data& A, const escript::Data& B,
        const escript::Data& C, const escript::Data& D,
        const escript::Data& X, const escript::Data& Y,
        const escript::Data& d, const escript::Data& y) const
{
    const double h0 = m_l0/m_gNE0;
    const double h1 = m_l1/m_gNE1;
    const double h2 = m_l2/m_gNE2;
    const double w0 = 0.0009303791403858427308*h1*h2/h0;
    const double w1 = 0.0009303791403858427308*h2;
    const double w2 = -0.00024929433932114870101*h1;
    const double w3 = 0.0009303791403858427308*h0*h2/h1;
    const double w4 = -0.00024929433932114870101*h0;
    const double w5 = 0.0009303791403858427308*h1;
    const double w6 = 0.0009303791403858427308*h0;
    const double w7 = -0.00024929433932114870101*h0*h1/h2;
    const double w8 = 0.0034722222222222222222*h2;
    const double w9 = -0.0009303791403858427308*h1;
    const double w10 = 0.012958509748503046158*h0*h2/h1;
    const double w11 = -0.0034722222222222222222*h0;
    const double w12 = 0.0034722222222222222222*h1;
    const double w13 = 0.012958509748503046158*h0;
    const double w14 = -0.0034722222222222222222*h0*h1/h2;
    const double w15 = 0.012958509748503046158*h1*h2/h0;
    const double w16 = -0.0034722222222222222222*h1;
    const double w17 = -0.0009303791403858427308*h0;
    const double w18 = 0.012958509748503046158*h1;
    const double w19 = 0.0034722222222222222222*h0;
    const double w20 = 0.012958509748503046158*h2;
    const double w21 = -0.012958509748503046158*h1;
    const double w22 = -0.012958509748503046158*h0;
    const double w23 = 0.04836181677178996241*h1;
    const double w24 = 0.04836181677178996241*h0;
    const double w25 = -0.04836181677178996241*h0*h1/h2;
    const double w26 = 0.00024929433932114870101*h1;
    const double w27 = 0.00024929433932114870101*h0;
    const double w28 = -0.04836181677178996241*h1;
    const double w29 = -0.04836181677178996241*h0;
    const double w30 = -0.0009303791403858427308*h1*h2/h0;
    const double w31 = -0.0009303791403858427308*h2;
    const double w32 = -0.0009303791403858427308*h0*h2/h1;
    const double w33 = 0.0034722222222222222222*h0*h1/h2;
    const double w34 = -0.0034722222222222222222*h2;
    const double w35 = -0.00024929433932114870101*h2;
    const double w36 = -0.012958509748503046158*h1*h2/h0;
    const double w37 = -0.012958509748503046158*h2;
    const double w38 = -0.012958509748503046158*h0*h2/h1;
    const double w39 = -0.04836181677178996241*h2;
    const double w40 = -0.0034722222222222222222*h0*h2/h1;
    const double w41 = 0.0009303791403858427308*h0*h1/h2;
    const double w42 = 0.04836181677178996241*h2;
    const double w43 = -0.04836181677178996241*h0*h2/h1;
    const double w44 = 0.012958509748503046158*h0*h1/h2;
    const double w45 = -0.00024929433932114870101*h0*h2/h1;
    const double w46 = 0.00024929433932114870101*h2;
    const double w47 = -0.0034722222222222222222*h1*h2/h0;
    const double w48 = -0.00024929433932114870101*h1*h2/h0;
    const double w49 = -0.04836181677178996241*h1*h2/h0;
    const double w50 = 0.0034722222222222222222*h0*h2/h1;
    const double w51 = -0.0009303791403858427308*h0*h1/h2;
    const double w52 = -0.012958509748503046158*h0*h1/h2;
    const double w53 = 0.0034722222222222222222*h1*h2/h0;
    const double w54 = 0.00024929433932114870101*h0*h1/h2;
    const double w55 = 0.04836181677178996241*h0*h2/h1;
    const double w56 = 0.04836181677178996241*h1*h2/h0;
    const double w57 = 0.04836181677178996241*h0*h1/h2;
    const double w58 = 0.00024929433932114870101*h1*h2/h0;
    const double w59 = 0.00024929433932114870101*h0*h2/h1;
    const double w60 = 0.055555555555555555556*h1*h2/h0;
    const double w61 = 0.041666666666666666667*h2;
    const double w62 = -0.083333333333333333333*h1;
    const double w63 = 0.055555555555555555556*h0*h2/h1;
    const double w64 = -0.083333333333333333333*h0;
    const double w65 = 0.083333333333333333333*h1;
    const double w66 = 0.083333333333333333333*h0;
    const double w67 = -0.11111111111111111111*h0*h1/h2;
    const double w68 = -0.055555555555555555556*h1*h2/h0;
    const double w69 = -0.083333333333333333333*h2;
    const double w70 = -0.041666666666666666667*h1;
    const double w71 = -0.055555555555555555556*h0*h2/h1;
    const double w72 = -0.041666666666666666667*h0;
    const double w73 = 0.041666666666666666667*h1;
    const double w74 = 0.041666666666666666667*h0;
    const double w75 = 0.027777777777777777778*h0*h1/h2;
    const double w76 = 0.083333333333333333333*h2;
    const double w77 = -0.11111111111111111111*h0*h2/h1;
    const double w78 = 0.055555555555555555556*h0*h1/h2;
    const double w79 = -0.11111111111111111111*h1*h2/h0;
    const double w80 = -0.027777777777777777778*h1*h2/h0;
    const double w81 = -0.041666666666666666667*h2;
    const double w82 = -0.027777777777777777778*h0*h2/h1;
    const double w83 = -0.027777777777777777778*h0*h1/h2;
    const double w84 = 0.027777777777777777778*h0*h2/h1;
    const double w85 = -0.055555555555555555556*h0*h1/h2;
    const double w86 = 0.11111111111111111111*h1*h2/h0;
    const double w87 = 0.11111111111111111111*h0*h2/h1;
    const double w88 = 0.11111111111111111111*h0*h1/h2;
    const double w89 = 0.027777777777777777778*h1*h2/h0;
    const double w90 = 0.0001966122466178319053*h1*h2;
    const double w91 = 0.0001966122466178319053*h0*h2;
    const double w92 = 0.0001966122466178319053*h0*h1;
    const double w93 = 0.0007337668937680108255*h1*h2;
    const double w94 = 0.0027384553284542113967*h0*h2;
    const double w95 = 0.0027384553284542113967*h0*h1;
    const double w96 = 0.0027384553284542113967*h1*h2;
    const double w97 = 0.0007337668937680108255*h0*h2;
    const double w98 = 0.010220054420048834761*h1*h2;
    const double w99 = 0.010220054420048834761*h0*h2;
    const double w100 = 0.038141762351741127649*h0*h1;
    const double w101 = 0.000052682092703316795705*h0*h1;
    const double w102 = 0.0007337668937680108255*h0*h1;
    const double w103 = 0.010220054420048834761*h0*h1;
    const double w104 = -0.0001966122466178319053*h1*h2;
    const double w105 = -0.0001966122466178319053*h0*h2;
    const double w106 = -0.0007337668937680108255*h1*h2;
    const double w107 = -0.0007337668937680108255*h0*h2;
    const double w108 = -0.0027384553284542113967*h1*h2;
    const double w109 = -0.0027384553284542113967*h0*h2;
    const double w110 = -0.010220054420048834761*h1*h2;
    const double w111 = -0.010220054420048834761*h0*h2;
    const double w112 = -0.0007337668937680108255*h0*h1;
    const double w113 = -0.010220054420048834761*h0*h1;
    const double w114 = -0.038141762351741127649*h0*h2;
    const double w115 = -0.000052682092703316795705*h0*h2;
    const double w116 = -0.0001966122466178319053*h0*h1;
    const double w117 = -0.0027384553284542113967*h0*h1;
    const double w118 = 0.000052682092703316795705*h0*h2;
    const double w119 = 0.038141762351741127649*h0*h2;
    const double w120 = 0.000052682092703316795705*h1*h2;
    const double w121 = 0.038141762351741127649*h1*h2;
    const double w122 = -0.000052682092703316795705*h0*h1;
    const double w123 = -0.038141762351741127649*h0*h1;
    const double w124 = -0.000052682092703316795705*h1*h2;
    const double w125 = -0.038141762351741127649*h1*h2;
    const double w126 = 0.027777777777777777778*h1*h2;
    const double w127 = 0.027777777777777777778*h0*h2;
    const double w128 = 0.055555555555555555556*h0*h1;
    const double w129 = -0.027777777777777777778*h1*h2;
    const double w130 = -0.027777777777777777778*h0*h2;
    const double w131 = 0.013888888888888888889*h0*h1;
    const double w132 = -0.055555555555555555556*h0*h2;
    const double w133 = -0.027777777777777777778*h0*h1;
    const double w134 = 0.055555555555555555556*h0*h2;
    const double w135 = 0.027777777777777777778*h0*h1;
    const double w136 = -0.013888888888888888889*h0*h1;
    const double w137 = 0.055555555555555555556*h1*h2;
    const double w138 = -0.013888888888888888889*h1*h2;
    const double w139 = -0.013888888888888888889*h0*h2;
    const double w140 = -0.055555555555555555556*h0*h1;
    const double w141 = 0.013888888888888888889*h1*h2;
    const double w142 = 0.013888888888888888889*h0*h2;
    const double w143 = -0.055555555555555555556*h1*h2;
    const double w144 = 0.000041549056553524783501*h0*h1*h2;
    const double w145 = 0.0005787037037037037037*h0*h1*h2;
    const double w146 = 0.0080603027952983270684*h0*h1*h2;
    const double w147 = 0.0001550631900643071218*h0*h1*h2;
    const double w148 = 0.002159751624750507693*h0*h1*h2;
    const double w149 = 0.03008145955644280058*h0*h1*h2;
    const double w150 = 0.000011133036149792012204*h0*h1*h2;
    const double w151 = 0.018518518518518518519*h0*h1*h2;
    const double w152 = 0.0092592592592592592592*h0*h1*h2;
    const double w153 = 0.0046296296296296296296*h0*h1*h2;
    const double w154 = 0.037037037037037037037*h0*h1*h2;
    const double w155 = -0.077751058491018276949*h1*h2;
    const double w156 = -0.077751058491018276949*h0*h2;
    const double w157 = -0.077751058491018276949*h0*h1;
    const double w158 = -0.020833333333333333333*h0*h2;
    const double w159 = -0.020833333333333333333*h0*h1;
    const double w160 = -0.020833333333333333333*h1*h2;
    const double w161 = -0.0055822748423150563848*h0*h1;
    const double w162 = -0.0055822748423150563848*h0*h2;
    const double w163 = -0.0055822748423150563848*h1*h2;
    const double w164 = 0.077751058491018276949*h1*h2;
    const double w165 = 0.020833333333333333333*h1*h2;
    const double w166 = 0.0055822748423150563848*h1*h2;
    const double w167 = 0.077751058491018276949*h0*h2;
    const double w168 = 0.020833333333333333333*h0*h2;
    const double w169 = 0.0055822748423150563848*h0*h2;
    const double w170 = 0.077751058491018276949*h0*h1;
    const double w171 = 0.020833333333333333333*h0*h1;
    const double w172 = 0.0055822748423150563848*h0*h1;
    const double w173 = -0.25*h1*h2;
    const double w174 = -0.25*h0*h2;
    const double w175 = -0.25*h0*h1;
    const double w176 = 0.25*h1*h2;
    const double w177 = 0.25*h0*h2;
    const double w178 = 0.25*h0*h1;
    const double w179 = 0.061320326520293008568*h0*h1*h2;
    const double w180 = 0.01643073197072526838*h0*h1*h2;
    const double w181 = 0.004402601362608064953*h0*h1*h2;
    const double w182 = 0.0011796734797069914318*h0*h1*h2;
    const double w183 = 0.125*h0*h1*h2;

    rhs.requireWrite();
#pragma omp parallel
    {
        for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
#pragma omp for
            for (index_t k2=k2_0; k2<m_NE2; k2+=2) {
                for (index_t k1=0; k1<m_NE1; ++k1) {
                    for (index_t k0=0; k0<m_NE0; ++k0)  {
                        bool add_EM_S=false;
                        bool add_EM_F=false;
                        vector<double> EM_S(8*8, 0);
                        vector<double> EM_F(8, 0);
                        const index_t e = k0 + m_NE0*k1 + m_NE0*m_NE1*k2;
                        ///////////////
                        // process A //
                        ///////////////
                        if (!A.isEmpty()) {
                            add_EM_S=true;
                            const double* A_p=const_cast<escript::Data*>(&A)->getSampleDataRO(e);
                            if (A.actsExpanded()) {
                                const register double A_00_0 = A_p[INDEX3(0,0,0,3,3)];
                                const register double A_01_0 = A_p[INDEX3(0,1,0,3,3)];
                                const register double A_02_0 = A_p[INDEX3(0,2,0,3,3)];
                                const register double A_10_0 = A_p[INDEX3(1,0,0,3,3)];
                                const register double A_11_0 = A_p[INDEX3(1,1,0,3,3)];
                                const register double A_12_0 = A_p[INDEX3(1,2,0,3,3)];
                                const register double A_20_0 = A_p[INDEX3(2,0,0,3,3)];
                                const register double A_21_0 = A_p[INDEX3(2,1,0,3,3)];
                                const register double A_22_0 = A_p[INDEX3(2,2,0,3,3)];
                                const register double A_00_1 = A_p[INDEX3(0,0,1,3,3)];
                                const register double A_01_1 = A_p[INDEX3(0,1,1,3,3)];
                                const register double A_02_1 = A_p[INDEX3(0,2,1,3,3)];
                                const register double A_10_1 = A_p[INDEX3(1,0,1,3,3)];
                                const register double A_11_1 = A_p[INDEX3(1,1,1,3,3)];
                                const register double A_12_1 = A_p[INDEX3(1,2,1,3,3)];
                                const register double A_20_1 = A_p[INDEX3(2,0,1,3,3)];
                                const register double A_21_1 = A_p[INDEX3(2,1,1,3,3)];
                                const register double A_22_1 = A_p[INDEX3(2,2,1,3,3)];
                                const register double A_00_2 = A_p[INDEX3(0,0,2,3,3)];
                                const register double A_01_2 = A_p[INDEX3(0,1,2,3,3)];
                                const register double A_02_2 = A_p[INDEX3(0,2,2,3,3)];
                                const register double A_10_2 = A_p[INDEX3(1,0,2,3,3)];
                                const register double A_11_2 = A_p[INDEX3(1,1,2,3,3)];
                                const register double A_12_2 = A_p[INDEX3(1,2,2,3,3)];
                                const register double A_20_2 = A_p[INDEX3(2,0,2,3,3)];
                                const register double A_21_2 = A_p[INDEX3(2,1,2,3,3)];
                                const register double A_22_2 = A_p[INDEX3(2,2,2,3,3)];
                                const register double A_00_3 = A_p[INDEX3(0,0,3,3,3)];
                                const register double A_01_3 = A_p[INDEX3(0,1,3,3,3)];
                                const register double A_02_3 = A_p[INDEX3(0,2,3,3,3)];
                                const register double A_10_3 = A_p[INDEX3(1,0,3,3,3)];
                                const register double A_11_3 = A_p[INDEX3(1,1,3,3,3)];
                                const register double A_12_3 = A_p[INDEX3(1,2,3,3,3)];
                                const register double A_20_3 = A_p[INDEX3(2,0,3,3,3)];
                                const register double A_21_3 = A_p[INDEX3(2,1,3,3,3)];
                                const register double A_22_3 = A_p[INDEX3(2,2,3,3,3)];
                                const register double A_00_4 = A_p[INDEX3(0,0,4,3,3)];
                                const register double A_01_4 = A_p[INDEX3(0,1,4,3,3)];
                                const register double A_02_4 = A_p[INDEX3(0,2,4,3,3)];
                                const register double A_10_4 = A_p[INDEX3(1,0,4,3,3)];
                                const register double A_11_4 = A_p[INDEX3(1,1,4,3,3)];
                                const register double A_12_4 = A_p[INDEX3(1,2,4,3,3)];
                                const register double A_20_4 = A_p[INDEX3(2,0,4,3,3)];
                                const register double A_21_4 = A_p[INDEX3(2,1,4,3,3)];
                                const register double A_22_4 = A_p[INDEX3(2,2,4,3,3)];
                                const register double A_00_5 = A_p[INDEX3(0,0,5,3,3)];
                                const register double A_01_5 = A_p[INDEX3(0,1,5,3,3)];
                                const register double A_02_5 = A_p[INDEX3(0,2,5,3,3)];
                                const register double A_10_5 = A_p[INDEX3(1,0,5,3,3)];
                                const register double A_11_5 = A_p[INDEX3(1,1,5,3,3)];
                                const register double A_12_5 = A_p[INDEX3(1,2,5,3,3)];
                                const register double A_20_5 = A_p[INDEX3(2,0,5,3,3)];
                                const register double A_21_5 = A_p[INDEX3(2,1,5,3,3)];
                                const register double A_22_5 = A_p[INDEX3(2,2,5,3,3)];
                                const register double A_00_6 = A_p[INDEX3(0,0,6,3,3)];
                                const register double A_01_6 = A_p[INDEX3(0,1,6,3,3)];
                                const register double A_02_6 = A_p[INDEX3(0,2,6,3,3)];
                                const register double A_10_6 = A_p[INDEX3(1,0,6,3,3)];
                                const register double A_11_6 = A_p[INDEX3(1,1,6,3,3)];
                                const register double A_12_6 = A_p[INDEX3(1,2,6,3,3)];
                                const register double A_20_6 = A_p[INDEX3(2,0,6,3,3)];
                                const register double A_21_6 = A_p[INDEX3(2,1,6,3,3)];
                                const register double A_22_6 = A_p[INDEX3(2,2,6,3,3)];
                                const register double A_00_7 = A_p[INDEX3(0,0,7,3,3)];
                                const register double A_01_7 = A_p[INDEX3(0,1,7,3,3)];
                                const register double A_02_7 = A_p[INDEX3(0,2,7,3,3)];
                                const register double A_10_7 = A_p[INDEX3(1,0,7,3,3)];
                                const register double A_11_7 = A_p[INDEX3(1,1,7,3,3)];
                                const register double A_12_7 = A_p[INDEX3(1,2,7,3,3)];
                                const register double A_20_7 = A_p[INDEX3(2,0,7,3,3)];
                                const register double A_21_7 = A_p[INDEX3(2,1,7,3,3)];
                                const register double A_22_7 = A_p[INDEX3(2,2,7,3,3)];
                                const register double tmp160_0 = A_12_0 + A_12_6 + A_21_0 + A_21_6;
                                const register double tmp8_0 = A_21_0 + A_21_6;
                                const register double tmp135_0 = A_10_1 + A_10_2 + A_10_5 + A_10_6;
                                const register double tmp67_0 = A_02_2 + A_02_7;
                                const register double tmp211_0 = A_12_6 + A_21_6;
                                const register double tmp180_0 = A_10_2 + A_10_6;
                                const register double tmp37_0 = A_00_0 + A_00_1 + A_00_2 + A_00_3;
                                const register double tmp92_0 = A_11_0 + A_11_1 + A_11_2 + A_11_3 + A_11_4 + A_11_5 + A_11_6 + A_11_7;
                                const register double tmp195_0 = A_02_2 + A_20_2;
                                const register double tmp70_0 = A_01_0 + A_01_7;
                                const register double tmp139_0 = A_02_3 + A_02_4 + A_20_1 + A_20_6;
                                const register double tmp200_0 = A_12_3 + A_12_5 + A_21_3 + A_21_5;
                                const register double tmp60_0 = A_22_0 + A_22_2 + A_22_4 + A_22_6;
                                const register double tmp192_0 = A_01_5 + A_10_5;
                                const register double tmp46_0 = A_21_0 + A_21_7;
                                const register double tmp48_0 = A_10_0 + A_10_7;
                                const register double tmp166_0 = A_11_5 + A_11_7;
                                const register double tmp221_0 = A_02_1 + A_02_6 + A_20_3 + A_20_4;
                                const register double tmp50_0 = A_02_4 + A_02_6 + A_20_4 + A_20_6;
                                const register double tmp217_0 = A_02_3 + A_02_4 + A_20_3 + A_20_4;
                                const register double tmp216_0 = A_01_2 + A_01_5 + A_10_2 + A_10_5;
                                const register double tmp104_0 = A_22_2 + A_22_6;
                                const register double tmp72_0 = A_20_3 + A_20_6;
                                const register double tmp79_0 = A_10_4 + A_10_7;
                                const register double tmp86_0 = A_01_2 + A_01_6 + A_10_1 + A_10_5;
                                const register double tmp214_0 = A_12_0 + A_12_7 + A_21_0 + A_21_7;
                                const register double tmp32_0 = A_02_0 + A_02_2;
                                const register double tmp112_0 = A_01_0 + A_01_4 + A_10_3 + A_10_7;
                                const register double tmp197_0 = A_12_0 + A_21_0;
                                const register double tmp106_0 = A_22_1 + A_22_5;
                                const register double tmp2_0 = A_00_0 + A_00_1 + A_00_4 + A_00_5;
                                const register double tmp115_0 = A_02_5 + A_02_7 + A_20_0 + A_20_2;
                                const register double tmp175_0 = A_01_3 + A_01_7;
                                const register double tmp126_0 = A_01_2 + A_01_5 + A_10_1 + A_10_6;
                                const register double tmp90_0 = A_00_0 + A_00_1 + A_00_2 + A_00_3 + A_00_4 + A_00_5 + A_00_6 + A_00_7;
                                const register double tmp47_0 = A_12_0 + A_12_6;
                                const register double tmp205_0 = A_02_7 + A_20_7;
                                const register double tmp148_0 = A_01_3 + A_01_4;
                                const register double tmp113_0 = A_01_3 + A_01_7 + A_10_0 + A_10_4;
                                const register double tmp43_0 = A_20_4 + A_20_6;
                                const register double tmp161_0 = A_02_1 + A_02_6 + A_20_1 + A_20_6;
                                const register double tmp69_0 = A_12_0 + A_12_1 + A_12_6 + A_12_7 + A_21_0 + A_21_1 + A_21_6 + A_21_7;
                                const register double tmp176_0 = A_01_1 + A_01_2 + A_01_5 + A_01_6;
                                const register double tmp105_0 = A_01_2 + A_01_6 + A_10_2 + A_10_6;
                                const register double tmp22_0 = A_01_5 + A_10_6;
                                const register double tmp91_0 = A_02_4 + A_02_6 + A_20_1 + A_20_3;
                                const register double tmp206_0 = A_12_7 + A_21_7;
                                const register double tmp188_0 = A_02_5 + A_20_5;
                                const register double tmp117_0 = A_21_1 + A_21_6;
                                const register double tmp165_0 = A_01_1 + A_01_6;
                                const register double tmp66_0 = A_00_4 + A_00_5;
                                const register double tmp57_0 = A_02_0 + A_02_2 + A_02_5 + A_02_7 + A_20_0 + A_20_2 + A_20_5 + A_20_7;
                                const register double tmp31_0 = A_21_4 + A_21_5;
                                const register double tmp3_0 = A_11_0 + A_11_2 + A_11_4 + A_11_6;
                                const register double tmp183_0 = A_12_0 + A_12_7;
                                const register double tmp61_0 = A_02_1 + A_02_3 + A_20_1 + A_20_3;
                                const register double tmp54_0 = A_10_5 + A_10_6;
                                const register double tmp18_0 = A_02_3 + A_02_6;
                                const register double tmp119_0 = A_12_2 + A_12_3 + A_12_4 + A_12_5 + A_21_2 + A_21_3 + A_21_4 + A_21_5;
                                const register double tmp29_0 = A_21_2 + A_21_3;
                                const register double tmp17_0 = A_01_3 + A_01_7 + A_10_3 + A_10_7;
                                const register double tmp212_0 = A_02_6 + A_20_6;
                                const register double tmp220_0 = A_02_3 + A_20_6;
                                const register double tmp78_0 = A_20_0 + A_20_7;
                                const register double tmp215_0 = A_01_6 + A_10_6;
                                const register double tmp203_0 = A_01_7 + A_10_7;
                                const register double tmp87_0 = A_12_2 + A_12_3 + A_21_4 + A_21_5;
                                const register double tmp114_0 = A_02_0 + A_02_2 + A_20_5 + A_20_7;
                                const register double tmp0_0 = A_01_0 + A_01_4 + A_10_0 + A_10_4;
                                const register double tmp202_0 = A_01_3 + A_01_4 + A_10_3 + A_10_4;
                                const register double tmp4_0 = A_20_0 + A_20_5;
                                const register double tmp65_0 = A_00_2 + A_00_3;
                                const register double tmp24_0 = A_20_1 + A_20_3;
                                const register double tmp64_0 = A_10_0 + A_10_3;
                                const register double tmp170_0 = A_02_0 + A_02_2 + A_20_0 + A_20_2;
                                const register double tmp11_0 = A_20_1 + A_20_6;
                                const register double tmp82_0 = A_12_4 + A_12_5 + A_21_4 + A_21_5;
                                const register double tmp99_0 = A_01_4 + A_10_7;
                                const register double tmp49_0 = A_12_1 + A_12_7;
                                const register double tmp130_0 = A_12_0 + A_12_1 + A_12_6 + A_12_7;
                                const register double tmp144_0 = A_01_0 + A_10_3;
                                const register double tmp109_0 = A_22_0 + A_22_3 + A_22_4 + A_22_7;
                                const register double tmp185_0 = A_02_0 + A_02_7 + A_20_2 + A_20_5;
                                const register double tmp157_0 = A_01_4 + A_10_4;
                                const register double tmp51_0 = A_22_1 + A_22_3 + A_22_5 + A_22_7;
                                const register double tmp146_0 = A_00_6 + A_00_7;
                                const register double tmp147_0 = A_12_0 + A_12_1 + A_21_0 + A_21_1;
                                const register double tmp150_0 = A_00_2 + A_00_3 + A_00_4 + A_00_5;
                                const register double tmp62_0 = A_21_3 + A_21_5;
                                const register double tmp223_0 = A_12_2 + A_21_4;
                                const register double tmp16_0 = A_02_2 + A_02_5;
                                const register double tmp168_0 = A_11_1 + A_11_3 + A_11_4 + A_11_6;
                                const register double tmp88_0 = A_12_4 + A_12_5 + A_21_2 + A_21_3;
                                const register double tmp142_0 = A_01_7 + A_10_4;
                                const register double tmp34_0 = A_20_0 + A_20_2 + A_20_5 + A_20_7;
                                const register double tmp71_0 = A_00_0 + A_00_1 + A_00_6 + A_00_7;
                                const register double tmp213_0 = A_02_1 + A_20_1;
                                const register double tmp227_0 = A_12_2 + A_12_5 + A_21_3 + A_21_4;
                                const register double tmp228_0 = A_12_1 + A_21_7;
                                const register double tmp140_0 = A_01_2 + A_01_6;
                                const register double tmp74_0 = A_22_0 + A_22_1 + A_22_4 + A_22_5;
                                const register double tmp167_0 = A_11_0 + A_11_2;
                                const register double tmp143_0 = A_01_3 + A_01_4 + A_10_0 + A_10_7;
                                const register double tmp83_0 = A_02_0 + A_02_5;
                                const register double tmp14_0 = A_22_1 + A_22_2 + A_22_5 + A_22_6;
                                const register double tmp5_0 = A_12_1 + A_12_6;
                                const register double tmp94_0 = A_02_1 + A_02_3;
                                const register double tmp193_0 = A_01_1 + A_01_6 + A_10_1 + A_10_6;
                                const register double tmp97_0 = A_02_0 + A_02_2 + A_02_5 + A_02_7;
                                const register double tmp131_0 = A_01_1 + A_01_5;
                                const register double tmp124_0 = A_01_6 + A_10_5;
                                const register double tmp149_0 = A_12_6 + A_12_7 + A_21_6 + A_21_7;
                                const register double tmp187_0 = A_01_2 + A_10_2;
                                const register double tmp93_0 = A_01_1 + A_01_2 + A_10_1 + A_10_2;
                                const register double tmp25_0 = A_01_4 + A_01_7 + A_10_4 + A_10_7;
                                const register double tmp156_0 = A_12_2 + A_12_5 + A_21_2 + A_21_5;
                                const register double tmp20_0 = A_21_2 + A_21_5;
                                const register double tmp55_0 = A_21_2 + A_21_4;
                                const register double tmp208_0 = A_12_1 + A_12_6 + A_21_0 + A_21_7;
                                const register double tmp125_0 = A_12_4 + A_12_5;
                                const register double tmp158_0 = A_01_0 + A_01_7 + A_10_0 + A_10_7;
                                const register double tmp108_0 = A_01_1 + A_01_5 + A_10_1 + A_10_5;
                                const register double tmp199_0 = A_12_2 + A_12_4 + A_21_2 + A_21_4;
                                const register double tmp10_0 = A_02_1 + A_02_4;
                                const register double tmp182_0 = A_02_3 + A_02_6 + A_20_3 + A_20_6;
                                const register double tmp132_0 = A_02_1 + A_20_4;
                                const register double tmp191_0 = A_12_3 + A_12_4 + A_21_3 + A_21_4;
                                const register double tmp35_0 = A_11_0 + A_11_1 + A_11_2 + A_11_3;
                                const register double tmp164_0 = A_10_3 + A_10_4;
                                const register double tmp190_0 = A_12_5 + A_21_5;
                                const register double tmp73_0 = A_02_1 + A_02_6;
                                const register double tmp98_0 = A_01_0 + A_01_7 + A_10_3 + A_10_4;
                                const register double tmp225_0 = A_12_4 + A_21_2;
                                const register double tmp103_0 = A_02_4 + A_02_6;
                                const register double tmp194_0 = A_02_0 + A_02_7 + A_20_0 + A_20_7;
                                const register double tmp207_0 = A_12_0 + A_21_6;
                                const register double tmp102_0 = A_20_5 + A_20_7;
                                const register double tmp1_0 = A_22_3 + A_22_7;
                                const register double tmp172_0 = A_10_1 + A_10_5;
                                const register double tmp222_0 = A_12_5 + A_21_3;
                                const register double tmp201_0 = A_02_2 + A_02_5 + A_20_2 + A_20_5;
                                const register double tmp155_0 = A_12_4 + A_21_4;
                                const register double tmp174_0 = A_02_1 + A_02_4 + A_20_1 + A_20_4;
                                const register double tmp59_0 = A_01_0 + A_01_3;
                                const register double tmp21_0 = A_20_2 + A_20_7;
                                const register double tmp141_0 = A_02_2 + A_02_7 + A_20_2 + A_20_7;
                                const register double tmp210_0 = A_01_1 + A_10_1;
                                const register double tmp145_0 = A_00_0 + A_00_1;
                                const register double tmp121_0 = A_12_0 + A_12_1 + A_21_6 + A_21_7;
                                const register double tmp224_0 = A_12_3 + A_12_4 + A_21_2 + A_21_5;
                                const register double tmp186_0 = A_02_2 + A_20_7;
                                const register double tmp53_0 = A_11_4 + A_11_6;
                                const register double tmp184_0 = A_02_5 + A_20_0;
                                const register double tmp38_0 = A_12_0 + A_12_1;
                                const register double tmp12_0 = A_01_1 + A_01_2 + A_01_5 + A_01_6 + A_10_1 + A_10_2 + A_10_5 + A_10_6;
                                const register double tmp230_0 = A_12_6 + A_21_0;
                                const register double tmp23_0 = A_11_4 + A_11_5 + A_11_6 + A_11_7;
                                const register double tmp81_0 = A_20_1 + A_20_4;
                                const register double tmp134_0 = A_10_3 + A_10_7;
                                const register double tmp129_0 = A_21_0 + A_21_1;
                                const register double tmp137_0 = A_01_0 + A_01_3 + A_01_4 + A_01_7;
                                const register double tmp198_0 = A_01_0 + A_10_0;
                                const register double tmp9_0 = A_21_1 + A_21_7;
                                const register double tmp179_0 = A_01_0 + A_01_4;
                                const register double tmp100_0 = A_20_1 + A_20_3 + A_20_4 + A_20_6;
                                const register double tmp173_0 = A_02_0 + A_20_5;
                                const register double tmp42_0 = A_21_0 + A_21_1 + A_21_6 + A_21_7;
                                const register double tmp226_0 = A_12_3 + A_21_5;
                                const register double tmp6_0 = A_22_0 + A_22_4;
                                const register double tmp218_0 = A_12_1 + A_21_1;
                                const register double tmp28_0 = A_01_2 + A_10_1;
                                const register double tmp133_0 = A_02_6 + A_20_3;
                                const register double tmp13_0 = A_00_2 + A_00_3 + A_00_6 + A_00_7;
                                const register double tmp27_0 = A_12_2 + A_12_3 + A_12_4 + A_12_5;
                                const register double tmp75_0 = A_10_1 + A_10_6;
                                const register double tmp36_0 = A_01_0 + A_01_3 + A_10_0 + A_10_3;
                                const register double tmp138_0 = A_10_0 + A_10_4;
                                const register double tmp189_0 = A_12_2 + A_21_2;
                                const register double tmp181_0 = A_02_7 + A_20_2;
                                const register double tmp85_0 = A_02_1 + A_02_3 + A_20_4 + A_20_6;
                                const register double tmp122_0 = A_01_1 + A_10_2;
                                const register double tmp95_0 = A_01_3 + A_10_0;
                                const register double tmp120_0 = A_12_6 + A_12_7 + A_21_0 + A_21_1;
                                const register double tmp196_0 = A_02_0 + A_20_0;
                                const register double tmp171_0 = A_02_3 + A_02_4;
                                const register double tmp204_0 = A_12_1 + A_12_6 + A_21_1 + A_21_6;
                                const register double tmp45_0 = A_10_1 + A_10_2;
                                const register double tmp101_0 = A_01_5 + A_01_6 + A_10_5 + A_10_6;
                                const register double tmp58_0 = A_11_0 + A_11_2 + A_11_5 + A_11_7;
                                const register double tmp107_0 = A_20_3 + A_20_4;
                                const register double tmp30_0 = A_01_1 + A_01_6 + A_10_2 + A_10_5;
                                const register double tmp63_0 = A_12_2 + A_12_5;
                                const register double tmp127_0 = A_12_2 + A_12_3;
                                const register double tmp177_0 = A_02_2 + A_02_5 + A_20_0 + A_20_7;
                                const register double tmp178_0 = A_10_0 + A_10_3 + A_10_4 + A_10_7;
                                const register double tmp76_0 = A_01_1 + A_01_2;
                                const register double tmp80_0 = A_22_2 + A_22_3 + A_22_6 + A_22_7;
                                const register double tmp41_0 = A_12_6 + A_12_7;
                                const register double tmp89_0 = A_01_0 + A_01_3 + A_01_4 + A_01_7 + A_10_0 + A_10_3 + A_10_4 + A_10_7;
                                const register double tmp116_0 = A_02_1 + A_02_3 + A_02_4 + A_02_6 + A_20_1 + A_20_3 + A_20_4 + A_20_6;
                                const register double tmp33_0 = A_22_0 + A_22_1 + A_22_2 + A_22_3 + A_22_4 + A_22_5 + A_22_6 + A_22_7;
                                const register double tmp169_0 = A_21_3 + A_21_4;
                                const register double tmp96_0 = A_20_0 + A_20_2;
                                const register double tmp111_0 = A_12_3 + A_12_4;
                                const register double tmp118_0 = A_20_2 + A_20_5;
                                const register double tmp19_0 = A_12_3 + A_12_5;
                                const register double tmp68_0 = A_01_5 + A_01_6;
                                const register double tmp7_0 = A_11_1 + A_11_3 + A_11_5 + A_11_7;
                                const register double tmp154_0 = A_12_3 + A_21_3;
                                const register double tmp152_0 = A_02_4 + A_20_4;
                                const register double tmp153_0 = A_02_3 + A_20_3;
                                const register double tmp163_0 = A_02_5 + A_02_7 + A_20_5 + A_20_7;
                                const register double tmp44_0 = A_01_4 + A_01_7;
                                const register double tmp39_0 = A_02_1 + A_02_3 + A_02_4 + A_02_6;
                                const register double tmp123_0 = A_21_2 + A_21_3 + A_21_4 + A_21_5;
                                const register double tmp40_0 = A_02_5 + A_02_7;
                                const register double tmp110_0 = A_02_0 + A_02_7;
                                const register double tmp77_0 = A_12_2 + A_12_3 + A_21_2 + A_21_3;
                                const register double tmp209_0 = A_12_7 + A_21_1;
                                const register double tmp219_0 = A_02_4 + A_20_1;
                                const register double tmp84_0 = A_01_1 + A_01_5 + A_10_2 + A_10_6;
                                const register double tmp162_0 = A_12_1 + A_12_7 + A_21_1 + A_21_7;
                                const register double tmp159_0 = A_01_3 + A_10_3;
                                const register double tmp56_0 = A_11_1 + A_11_3;
                                const register double tmp52_0 = A_01_2 + A_01_5;
                                const register double tmp26_0 = A_00_4 + A_00_5 + A_00_6 + A_00_7;
                                const register double tmp229_0 = A_12_0 + A_12_7 + A_21_1 + A_21_6;
                                const register double tmp151_0 = A_10_2 + A_10_5;
                                const register double tmp136_0 = A_02_0 + A_02_5 + A_20_0 + A_20_5;
                                const register double tmp128_0 = A_21_6 + A_21_7;
                                const register double tmp15_0 = A_12_2 + A_12_4;
                                const register double tmp296_1 = tmp159_0*w42;
                                const register double tmp130_1 = tmp67_0*w5;
                                const register double tmp98_1 = A_01_6*w42;
                                const register double tmp231_1 = tmp125_0*w6;
                                const register double tmp42_1 = tmp34_0*w12;
                                const register double tmp199_1 = A_02_5*w28;
                                const register double tmp113_1 = tmp29_0*w13;
                                const register double tmp330_1 = tmp152_0*w28;
                                const register double tmp90_1 = A_01_1*w46;
                                const register double tmp446_1 = tmp77_0*w22;
                                const register double tmp108_1 = tmp43_0*w5;
                                const register double tmp524_1 = A_12_6*w29;
                                const register double tmp232_1 = tmp126_0*w34;
                                const register double tmp33_1 = tmp25_0*w37;
                                const register double tmp461_1 = tmp180_0*w1;
                                const register double tmp14_1 = tmp8_0*w6;
                                const register double tmp447_1 = tmp205_0*w26;
                                const register double tmp452_1 = tmp198_0*w42;
                                const register double tmp217_1 = tmp81_0*w9;
                                const register double tmp76_1 = tmp59_0*w20;
                                const register double tmp421_1 = tmp134_0*w31;
                                const register double tmp485_1 = tmp51_0*w51;
                                const register double tmp240_1 = tmp131_0*w1;
                                const register double tmp160_1 = tmp91_0*w9;
                                const register double tmp174_1 = A_20_1*w26;
                                const register double tmp273_1 = A_10_1*w46;
                                const register double tmp159_1 = tmp90_0*w47;
                                const register double tmp228_1 = tmp103_0*w5;
                                const register double tmp313_1 = tmp166_0*w45;
                                const register double tmp45_1 = tmp37_0*w30;
                                const register double tmp512_1 = tmp147_0*w13;
                                const register double tmp73_1 = tmp56_0*w43;
                                const register double tmp61_1 = A_01_6*w46;
                                const register double tmp316_1 = tmp167_0*w43;
                                const register double tmp189_1 = tmp112_0*w20;
                                const register double tmp455_1 = tmp215_0*w39;
                                const register double tmp360_1 = A_21_5*w24;
                                const register double tmp258_1 = A_20_7*w2;
                                const register double tmp196_1 = A_20_6*w26;
                                const register double tmp37_1 = tmp29_0*w6;
                                const register double tmp9_1 = A_12_7*w29;
                                const register double tmp80_1 = tmp63_0*w19;
                                const register double tmp312_1 = tmp165_0*w8;
                                const register double tmp264_1 = tmp101_0*w1;
                                const register double tmp124_1 = A_02_3*w26;
                                const register double tmp229_1 = tmp123_0*w11;
                                const register double tmp333_1 = tmp159_0*w46;
                                const register double tmp533_1 = tmp222_0*w4;
                                const register double tmp201_1 = tmp108_0*w37;
                                const register double tmp444_1 = tmp35_0*w10;
                                const register double tmp51_1 = tmp43_0*w18;
                                const register double tmp214_1 = A_21_7*w29;
                                const register double tmp518_1 = tmp86_0*w37;
                                const register double tmp192_1 = tmp115_0*w5;
                                const register double tmp355_1 = A_21_2*w27;
                                const register double tmp156_1 = tmp87_0*w22;
                                const register double tmp516_1 = tmp230_0*w27;
                                const register double tmp366_1 = tmp104_0*w57;
                                const register double tmp271_1 = tmp146_0*w49;
                                const register double tmp437_1 = tmp218_0*w24;
                                const register double tmp436_1 = tmp104_0*w54;
                                const register double tmp167_1 = tmp98_0*w8;
                                const register double tmp136_1 = tmp70_0*w34;
                                const register double tmp406_1 = tmp207_0*w27;
                                const register double tmp193_1 = tmp116_0*w12;
                                const register double tmp486_1 = tmp225_0*w29;
                                const register double tmp469_1 = tmp224_0*w11;
                                const register double tmp287_1 = tmp71_0*w53;
                                const register double tmp430_1 = tmp213_0*w28;
                                const register double tmp462_1 = tmp220_0*w2;
                                const register double tmp294_1 = tmp53_0*w59;
                                const register double tmp218_1 = tmp118_0*w16;
                                const register double tmp116_1 = tmp25_0*w31;
                                const register double tmp495_1 = tmp76_0*w37;
                                const register double tmp501_1 = tmp99_0*w46;
                                const register double tmp0_1 = tmp0_0*w1;
                                const register double tmp99_1 = tmp62_0*w17;
                                const register double tmp429_1 = tmp212_0*w2;
                                const register double tmp249_1 = tmp136_0*w9;
                                const register double tmp504_1 = tmp229_0*w19;
                                const register double tmp197_1 = A_12_2*w27;
                                const register double tmp531_1 = tmp122_0*w35;
                                const register double tmp265_1 = tmp142_0*w46;
                                const register double tmp488_1 = tmp226_0*w4;
                                const register double tmp528_1 = tmp115_0*w18;
                                const register double tmp438_1 = tmp219_0*w2;
                                const register double tmp233_1 = tmp127_0*w13;
                                const register double tmp491_1 = tmp79_0*w1;
                                const register double tmp215_1 = A_21_0*w4;
                                const register double tmp24_1 = tmp18_0*w21;
                                const register double tmp538_1 = tmp209_0*w27;
                                const register double tmp379_1 = tmp167_0*w55;
                                const register double tmp332_1 = tmp154_0*w4;
                                const register double tmp498_1 = tmp68_0*w31;
                                const register double tmp41_1 = tmp33_0*w33;
                                const register double tmp464_1 = tmp179_0*w37;
                                const register double tmp317_1 = tmp168_0*w40;
                                const register double tmp378_1 = tmp106_0*w54;
                                const register double tmp184_1 = tmp109_0*w14;
                                const register double tmp292_1 = tmp14_0*w33;
                                const register double tmp11_1 = tmp5_0*w11;
                                const register double tmp354_1 = A_02_6*w26;
                                const register double tmp84_1 = tmp37_0*w0;
                                const register double tmp422_1 = tmp13_0*w30;
                                const register double tmp132_1 = tmp69_0*w11;
                                const register double tmp251_1 = tmp138_0*w31;
                                const register double tmp18_1 = tmp12_0*w8;
                                const register double tmp88_1 = A_21_1*w4;
                                const register double tmp188_1 = A_12_2*w24;
                                const register double tmp465_1 = tmp175_0*w31;
                                const register double tmp235_1 = tmp128_0*w17;
                                const register double tmp323_1 = A_02_1*w26;
                                const register double tmp31_1 = tmp23_0*w38;
                                const register double tmp397_1 = tmp170_0*w5;
                                const register double tmp175_1 = tmp7_0*w3;
                                const register double tmp148_1 = tmp81_0*w21;
                                const register double tmp238_1 = tmp130_0*w19;
                                const register double tmp59_1 = tmp46_0*w11;
                                const register double tmp432_1 = tmp215_0*w35;
                                const register double tmp398_1 = A_01_2*w46;
                                const register double tmp497_1 = A_10_5*w46;
                                const register double tmp28_1 = tmp21_0*w18;
                                const register double tmp115_1 = tmp23_0*w32;
                                const register double tmp441_1 = tmp23_0*w3;
                                const register double tmp131_1 = tmp68_0*w37;
                                const register double tmp289_1 = tmp155_0*w4;
                                const register double tmp278_1 = tmp80_0*w44;
                                const register double tmp5_1 = A_21_4*w27;
                                const register double tmp254_1 = tmp140_0*w20;
                                const register double tmp183_1 = tmp108_0*w31;
                                const register double tmp279_1 = tmp151_0*w8;
                                const register double tmp298_1 = tmp161_0*w16;
                                const register double tmp505_1 = tmp230_0*w24;
                                const register double tmp246_1 = tmp80_0*w52;
                                const register double tmp100_1 = tmp53_0*w43;
                                const register double tmp440_1 = tmp221_0*w16;
                                const register double tmp481_1 = tmp188_0*w23;
                                const register double tmp480_1 = tmp187_0*w35;
                                const register double tmp384_1 = tmp150_0*w53;
                                const register double tmp142_1 = tmp76_0*w31;
                                const register double tmp372_1 = tmp191_0*w11;
                                const register double tmp307_1 = A_10_7*w35;
                                const register double tmp186_1 = tmp111_0*w19;
                                const register double tmp127_1 = A_20_2*w2;
                                const register double tmp391_1 = tmp167_0*w59;
                                const register double tmp223_1 = tmp113_0*w20;
                                const register double tmp454_1 = tmp197_0*w24;
                                const register double tmp241_1 = tmp74_0*w51;
                                const register double tmp529_1 = tmp114_0*w5;
                                const register double tmp202_1 = tmp104_0*w7;
                                const register double tmp236_1 = tmp96_0*w21;
                                const register double tmp358_1 = tmp183_0*w11;
                                const register double tmp102_1 = tmp51_0*w41;
                                const register double tmp493_1 = A_20_5*w2;
                                const register double tmp468_1 = tmp223_0*w4;
                                const register double tmp435_1 = tmp217_0*w16;
                                const register double tmp110_1 = tmp37_0*w36;
                                const register double tmp479_1 = tmp189_0*w4;
                                const register double tmp120_1 = tmp38_0*w22;
                                const register double tmp16_1 = tmp10_0*w9;
                                const register double tmp407_1 = tmp90_0*w53;
                                const register double tmp442_1 = tmp66_0*w48;
                                const register double tmp60_1 = A_10_4*w35;
                                const register double tmp69_1 = tmp53_0*w45;
                                const register double tmp144_1 = tmp77_0*w17;
                                const register double tmp507_1 = tmp146_0*w48;
                                const register double tmp424_1 = tmp174_0*w18;
                                const register double tmp352_1 = tmp181_0*w23;
                                const register double tmp451_1 = tmp199_0*w13;
                                const register double tmp253_1 = tmp139_0*w16;
                                const register double tmp353_1 = tmp182_0*w18;
                                const register double tmp521_1 = tmp88_0*w22;
                                const register double tmp346_1 = tmp175_0*w37;
                                const register double tmp416_1 = tmp138_0*w37;
                                const register double tmp324_1 = A_10_0*w35;
                                const register double tmp152_1 = tmp84_0*w37;
                                const register double tmp119_1 = tmp32_0*w21;
                                const register double tmp86_1 = A_21_6*w29;
                                const register double tmp290_1 = tmp156_0*w11;
                                const register double tmp382_1 = tmp196_0*w26;
                                const register double tmp91_1 = tmp49_0*w6;
                                const register double tmp499_1 = A_10_2*w42;
                                const register double tmp226_1 = tmp121_0*w13;
                                const register double tmp477_1 = tmp195_0*w26;
                                const register double tmp150_1 = A_02_4*w23;
                                const register double tmp318_1 = tmp15_0*w22;
                                const register double tmp396_1 = tmp206_0*w24;
                                const register double tmp474_1 = A_02_0*w28;
                                const register double tmp245_1 = tmp134_0*w37;
                                const register double tmp3_1 = A_20_4*w26;
                                const register double tmp44_1 = tmp36_0*w31;
                                const register double tmp487_1 = tmp60_0*w52;
                                const register double tmp293_1 = tmp158_0*w8;
                                const register double tmp314_1 = A_01_2*w42;
                                const register double tmp414_1 = tmp80_0*w51;
                                const register double tmp472_1 = A_21_3*w27;
                                const register double tmp321_1 = A_21_2*w24;
                                const register double tmp225_1 = tmp120_0*w6;
                                const register double tmp377_1 = tmp166_0*w59;
                                const register double tmp413_1 = tmp186_0*w26;
                                const register double tmp385_1 = tmp166_0*w55;
                                const register double tmp310_1 = tmp164_0*w34;
                                const register double tmp158_1 = tmp89_0*w34;
                                const register double tmp449_1 = tmp203_0*w46;
                                const register double tmp439_1 = tmp220_0*w28;
                                const register double tmp22_1 = tmp16_0*w16;
                                const register double tmp164_1 = tmp95_0*w46;
                                const register double tmp417_1 = tmp74_0*w52;
                                const register double tmp257_1 = tmp6_0*w25;
                                const register double tmp203_1 = tmp18_0*w9;
                                const register double tmp286_1 = tmp153_0*w28;
                                const register double tmp155_1 = tmp33_0*w14;
                                const register double tmp389_1 = tmp201_0*w12;
                                const register double tmp508_1 = tmp145_0*w49;
                                const register double tmp300_1 = tmp56_0*w55;
                                const register double tmp299_1 = tmp162_0*w22;
                                const register double tmp173_1 = tmp104_0*w25;
                                const register double tmp32_1 = tmp24_0*w5;
                                const register double tmp227_1 = tmp122_0*w39;
                                const register double tmp484_1 = tmp3_0*w38;
                                const register double tmp171_1 = tmp102_0*w21;
                                const register double tmp478_1 = tmp190_0*w29;
                                const register double tmp320_1 = tmp170_0*w18;
                                const register double tmp327_1 = tmp6_0*w57;
                                const register double tmp490_1 = tmp7_0*w32;
                                const register double tmp419_1 = tmp127_0*w6;
                                const register double tmp463_1 = tmp219_0*w28;
                                const register double tmp12_1 = tmp6_0*w7;
                                const register double tmp49_1 = tmp41_0*w22;
                                const register double tmp344_1 = tmp173_0*w26;
                                const register double tmp243_1 = tmp132_0*w2;
                                const register double tmp83_1 = A_10_4*w39;
                                const register double tmp297_1 = tmp160_0*w17;
                                const register double tmp275_1 = tmp148_0*w34;
                                const register double tmp168_1 = tmp99_0*w42;
                                const register double tmp409_1 = tmp3_0*w32;
                                const register double tmp1_1 = tmp1_0*w25;
                                const register double tmp426_1 = tmp210_0*w39;
                                const register double tmp375_1 = tmp109_0*w33;
                                const register double tmp50_1 = tmp42_0*w19;
                                const register double tmp513_1 = A_10_1*w42;
                                const register double tmp97_1 = tmp45_0*w31;
                                const register double tmp403_1 = tmp140_0*w1;
                                const register double tmp71_1 = A_01_1*w42;
                                const register double tmp520_1 = tmp84_0*w31;
                                const register double tmp510_1 = A_10_6*w46;
                                const register double tmp302_1 = A_10_0*w39;
                                const register double tmp364_1 = tmp128_0*w22;
                                const register double tmp515_1 = tmp142_0*w42;
                                const register double tmp283_1 = tmp65_0*w56;
                                const register double tmp222_1 = tmp112_0*w1;
                                const register double tmp428_1 = tmp211_0*w27;
                                const register double tmp371_1 = tmp190_0*w4;
                                const register double tmp423_1 = tmp184_0*w23;
                                const register double tmp276_1 = tmp149_0*w13;
                                const register double tmp65_1 = tmp50_0*w9;
                                const register double tmp305_1 = A_12_0*w29;
                                const register double tmp170_1 = tmp101_0*w20;
                                const register double tmp350_1 = tmp179_0*w31;
                                const register double tmp466_1 = tmp172_0*w20;
                                const register double tmp361_1 = tmp184_0*w26;
                                const register double tmp431_1 = tmp214_0*w19;
                                const register double tmp363_1 = tmp129_0*w17;
                                const register double tmp178_1 = A_02_2*w28;
                                const register double tmp527_1 = tmp120_0*w13;
                                const register double tmp415_1 = tmp182_0*w5;
                                const register double tmp450_1 = tmp200_0*w6;
                                const register double tmp269_1 = A_01_7*w39;
                                const register double tmp285_1 = tmp152_0*w2;
                                const register double tmp272_1 = A_01_0*w35;
                                const register double tmp339_1 = tmp136_0*w21;
                                const register double tmp502_1 = tmp95_0*w42;
                                const register double tmp38_1 = tmp30_0*w34;
                                const register double tmp514_1 = tmp144_0*w46;
                                const register double tmp96_1 = tmp56_0*w45;
                                const register double tmp399_1 = tmp167_0*w45;
                                const register double tmp483_1 = tmp173_0*w23;
                                const register double tmp522_1 = tmp87_0*w17;
                                const register double tmp519_1 = tmp91_0*w21;
                                const register double tmp209_1 = A_12_5*w24;
                                const register double tmp126_1 = tmp65_0*w48;
                                const register double tmp367_1 = tmp187_0*w39;
                                const register double tmp221_1 = tmp67_0*w18;
                                const register double tmp381_1 = tmp146_0*w56;
                                const register double tmp70_1 = tmp54_0*w31;
                                const register double tmp216_1 = tmp117_0*w11;
                                const register double tmp473_1 = A_02_7*w2;
                                const register double tmp149_1 = tmp82_0*w22;
                                const register double tmp357_1 = A_12_6*w4;
                                const register double tmp534_1 = tmp226_0*w29;
                                const register double tmp95_1 = tmp26_0*w15;
                                const register double tmp500_1 = tmp64_0*w20;
                                const register double tmp387_1 = tmp199_0*w6;
                                const register double tmp471_1 = A_20_4*w23;
                                const register double tmp281_1 = tmp74_0*w41;
                                const register double tmp351_1 = tmp180_0*w20;
                                const register double tmp63_1 = tmp48_0*w34;
                                const register double tmp365_1 = tmp186_0*w23;
                                const register double tmp448_1 = tmp206_0*w27;
                                const register double tmp39_1 = tmp31_0*w13;
                                const register double tmp453_1 = tmp196_0*w23;
                                const register double tmp402_1 = tmp163_0*w18;
                                const register double tmp137_1 = tmp71_0*w47;
                                const register double tmp6_1 = A_02_0*w2;
                                const register double tmp34_1 = tmp26_0*w36;
                                const register double tmp383_1 = tmp197_0*w27;
                                const register double tmp166_1 = tmp97_0*w12;
                                const register double tmp114_1 = tmp40_0*w9;
                                const register double tmp306_1 = A_12_7*w4;
                                const register double tmp530_1 = tmp124_0*w39;
                                const register double tmp388_1 = tmp200_0*w13;
                                const register double tmp252_1 = tmp2_0*w30;
                                const register double tmp210_1 = A_02_4*w26;
                                const register double tmp200_1 = tmp21_0*w5;
                                const register double tmp181_1 = tmp3_0*w10;
                                const register double tmp425_1 = tmp106_0*w57;
                                const register double tmp261_1 = A_21_7*w4;
                                const register double tmp64_1 = tmp49_0*w13;
                                const register double tmp506_1 = A_01_0*w39;
                                const register double tmp457_1 = tmp213_0*w2;
                                const register double tmp2_1 = tmp2_0*w0;
                                const register double tmp393_1 = tmp203_0*w42;
                                const register double tmp133_1 = A_01_3*w35;
                                const register double tmp147_1 = tmp80_0*w41;
                                const register double tmp8_1 = tmp4_0*w5;
                                const register double tmp267_1 = tmp144_0*w42;
                                const register double tmp17_1 = tmp11_0*w12;
                                const register double tmp284_1 = tmp58_0*w50;
                                const register double tmp328_1 = tmp66_0*w56;
                                const register double tmp405_1 = tmp60_0*w51;
                                const register double tmp467_1 = tmp222_0*w29;
                                const register double tmp535_1 = tmp225_0*w4;
                                const register double tmp356_1 = A_12_1*w29;
                                const register double tmp274_1 = tmp147_0*w6;
                                const register double tmp476_1 = tmp192_0*w39;
                                const register double tmp206_1 = tmp10_0*w21;
                                const register double tmp334_1 = tmp141_0*w9;
                                const register double tmp482_1 = tmp181_0*w26;
                                const register double tmp212_1 = A_20_7*w28;
                                const register double tmp219_1 = tmp72_0*w21;
                                const register double tmp47_1 = tmp39_0*w16;
                                const register double tmp89_1 = A_10_3*w35;
                                const register double tmp52_1 = tmp44_0*w1;
                                const register double tmp492_1 = A_01_3*w39;
                                const register double tmp81_1 = A_12_3*w24;
                                const register double tmp77_1 = tmp60_0*w41;
                                const register double tmp153_1 = tmp85_0*w21;
                                const register double tmp304_1 = tmp163_0*w5;
                                const register double tmp489_1 = tmp227_0*w11;
                                const register double tmp107_1 = tmp35_0*w38;
                                const register double tmp30_1 = tmp22_0*w39;
                                const register double tmp260_1 = A_21_0*w29;
                                const register double tmp343_1 = tmp172_0*w1;
                                const register double tmp511_1 = tmp149_0*w6;
                                const register double tmp139_1 = tmp73_0*w12;
                                const register double tmp66_1 = tmp51_0*w44;
                                const register double tmp208_1 = tmp4_0*w18;
                                const register double tmp134_1 = tmp23_0*w10;
                                const register double tmp205_1 = tmp105_0*w31;
                                const register double tmp349_1 = tmp178_0*w8;
                                const register double tmp341_1 = tmp53_0*w55;
                                const register double tmp72_1 = tmp55_0*w17;
                                const register double tmp79_1 = tmp62_0*w22;
                                const register double tmp26_1 = tmp20_0*w19;
                                const register double tmp141_1 = tmp75_0*w8;
                                const register double tmp118_1 = tmp41_0*w17;
                                const register double tmp259_1 = A_20_0*w28;
                                const register double tmp458_1 = tmp212_0*w28;
                                const register double tmp68_1 = tmp37_0*w15;
                                const register double tmp154_1 = tmp86_0*w31;
                                const register double tmp335_1 = tmp56_0*w59;
                                const register double tmp359_1 = A_02_1*w23;
                                const register double tmp56_1 = A_21_1*w29;
                                const register double tmp392_1 = tmp145_0*w58;
                                const register double tmp270_1 = tmp145_0*w48;
                                const register double tmp92_1 = tmp47_0*w13;
                                const register double tmp433_1 = tmp216_0*w34;
                                const register double tmp420_1 = tmp125_0*w13;
                                const register double tmp408_1 = tmp51_0*w52;
                                const register double tmp494_1 = A_20_2*w28;
                                const register double tmp362_1 = tmp185_0*w12;
                                const register double tmp411_1 = tmp208_0*w19;
                                const register double tmp336_1 = tmp65_0*w58;
                                const register double tmp475_1 = A_21_4*w24;
                                const register double tmp85_1 = A_12_3*w27;
                                const register double tmp19_1 = tmp13_0*w15;
                                const register double tmp537_1 = tmp132_0*w28;
                                const register double tmp67_1 = tmp52_0*w8;
                                const register double tmp459_1 = tmp210_0*w35;
                                const register double tmp248_1 = tmp135_0*w34;
                                const register double tmp326_1 = A_02_6*w23;
                                const register double tmp23_1 = tmp17_0*w20;
                                const register double tmp35_1 = tmp27_0*w11;
                                const register double tmp62_1 = tmp47_0*w6;
                                const register double tmp180_1 = tmp106_0*w7;
                                const register double tmp277_1 = tmp150_0*w47;
                                const register double tmp373_1 = tmp192_0*w35;
                                const register double tmp337_1 = tmp157_0*w42;
                                const register double tmp106_1 = tmp28_0*w39;
                                const register double tmp369_1 = tmp168_0*w50;
                                const register double tmp434_1 = tmp146_0*w58;
                                const register double tmp331_1 = tmp155_0*w29;
                                const register double tmp503_1 = tmp228_0*w27;
                                const register double tmp93_1 = tmp61_0*w9;
                                const register double tmp25_1 = tmp19_0*w22;
                                const register double tmp146_1 = tmp79_0*w20;
                                const register double tmp280_1 = A_10_6*w42;
                                const register double tmp94_1 = tmp60_0*w44;
                                const register double tmp400_1 = A_01_5*w42;
                                const register double tmp151_1 = tmp83_0*w18;
                                const register double tmp78_1 = tmp61_0*w21;
                                const register double tmp301_1 = tmp6_0*w54;
                                const register double tmp48_1 = tmp40_0*w21;
                                const register double tmp75_1 = tmp58_0*w40;
                                const register double tmp82_1 = tmp59_0*w1;
                                const register double tmp74_1 = tmp57_0*w16;
                                const register double tmp36_1 = tmp28_0*w35;
                                const register double tmp370_1 = tmp189_0*w29;
                                const register double tmp224_1 = tmp119_0*w19;
                                const register double tmp109_1 = tmp36_0*w37;
                                const register double tmp345_1 = tmp174_0*w5;
                                const register double tmp101_1 = tmp44_0*w20;
                                const register double tmp308_1 = A_01_5*w46;
                                const register double tmp295_1 = tmp66_0*w58;
                                const register double tmp117_1 = tmp26_0*w30;
                                const register double tmp125_1 = tmp35_0*w3;
                                const register double tmp309_1 = tmp9_0*w6;
                                const register double tmp412_1 = tmp209_0*w24;
                                const register double tmp46_1 = tmp38_0*w17;
                                const register double tmp7_1 = A_02_7*w28;
                                const register double tmp40_1 = tmp32_0*w9;
                                const register double tmp386_1 = tmp198_0*w46;
                                const register double tmp517_1 = tmp228_0*w24;
                                const register double tmp532_1 = tmp223_0*w29;
                                const register double tmp220_1 = A_02_3*w23;
                                const register double tmp268_1 = tmp93_0*w20;
                                const register double tmp322_1 = A_10_7*w39;
                                const register double tmp311_1 = tmp8_0*w13;
                                const register double tmp123_1 = A_01_4*w39;
                                const register double tmp187_1 = A_20_6*w23;
                                const register double tmp177_1 = A_02_5*w2;
                                const register double tmp58_1 = A_21_6*w4;
                                const register double tmp404_1 = tmp7_0*w38;
                                const register double tmp122_1 = tmp64_0*w1;
                                const register double tmp163_1 = tmp94_0*w5;
                                const register double tmp15_1 = tmp9_0*w13;
                                const register double tmp128_1 = A_20_5*w28;
                                const register double tmp204_1 = tmp2_0*w15;
                                const register double tmp539_1 = tmp207_0*w24;
                                const register double tmp57_1 = tmp45_0*w37;
                                const register double tmp53_1 = A_10_3*w39;
                                const register double tmp157_1 = tmp88_0*w17;
                                const register double tmp169_1 = tmp100_0*w16;
                                const register double tmp162_1 = tmp93_0*w1;
                                const register double tmp325_1 = tmp171_0*w12;
                                const register double tmp179_1 = tmp105_0*w37;
                                const register double tmp207_1 = A_20_1*w23;
                                const register double tmp427_1 = tmp145_0*w56;
                                const register double tmp368_1 = tmp188_0*w26;
                                const register double tmp460_1 = tmp211_0*w24;
                                const register double tmp347_1 = tmp176_0*w34;
                                const register double tmp234_1 = tmp102_0*w9;
                                const register double tmp21_1 = tmp15_0*w17;
                                const register double tmp112_1 = tmp31_0*w6;
                                const register double tmp319_1 = tmp169_0*w19;
                                const register double tmp509_1 = A_01_7*w35;
                                const register double tmp418_1 = tmp2_0*w36;
                                const register double tmp266_1 = tmp143_0*w8;
                                const register double tmp244_1 = tmp133_0*w28;
                                const register double tmp138_1 = tmp72_0*w9;
                                const register double tmp20_1 = tmp14_0*w14;
                                const register double tmp29_1 = A_21_3*w24;
                                const register double tmp190_1 = tmp113_0*w1;
                                const register double tmp338_1 = tmp162_0*w17;
                                const register double tmp87_1 = tmp54_0*w37;
                                const register double tmp374_1 = tmp193_0*w34;
                                const register double tmp195_1 = tmp13_0*w0;
                                const register double tmp194_1 = tmp106_0*w25;
                                const register double tmp111_1 = tmp22_0*w35;
                                const register double tmp213_1 = tmp83_0*w5;
                                const register double tmp291_1 = tmp157_0*w46;
                                const register double tmp329_1 = tmp153_0*w2;
                                const register double tmp256_1 = tmp17_0*w1;
                                const register double tmp342_1 = tmp1_0*w54;
                                const register double tmp376_1 = tmp194_0*w12;
                                const register double tmp239_1 = tmp94_0*w18;
                                const register double tmp340_1 = tmp160_0*w22;
                                const register double tmp262_1 = tmp1_0*w7;
                                const register double tmp54_1 = tmp26_0*w0;
                                const register double tmp470_1 = A_20_3*w26;
                                const register double tmp121_1 = tmp24_0*w18;
                                const register double tmp523_1 = tmp85_0*w9;
                                const register double tmp182_1 = tmp107_0*w12;
                                const register double tmp140_1 = tmp74_0*w44;
                                const register double tmp250_1 = tmp137_0*w8;
                                const register double tmp104_1 = tmp55_0*w22;
                                const register double tmp303_1 = A_21_5*w27;
                                const register double tmp4_1 = tmp3_0*w3;
                                const register double tmp526_1 = tmp121_0*w6;
                                const register double tmp410_1 = tmp131_0*w20;
                                const register double tmp255_1 = tmp141_0*w21;
                                const register double tmp394_1 = tmp204_0*w19;
                                const register double tmp172_1 = tmp103_0*w18;
                                const register double tmp211_1 = A_20_0*w2;
                                const register double tmp263_1 = tmp0_0*w20;
                                const register double tmp230_1 = tmp124_0*w35;
                                const register double tmp496_1 = A_01_4*w35;
                                const register double tmp176_1 = A_12_5*w27;
                                const register double tmp401_1 = tmp166_0*w43;
                                const register double tmp43_1 = tmp35_0*w32;
                                const register double tmp395_1 = tmp205_0*w23;
                                const register double tmp348_1 = tmp177_0*w12;
                                const register double tmp390_1 = tmp202_0*w8;
                                const register double tmp525_1 = A_12_1*w4;
                                const register double tmp237_1 = tmp129_0*w22;
                                const register double tmp105_1 = A_12_4*w24;
                                const register double tmp242_1 = tmp92_0*w50;
                                const register double tmp161_1 = tmp92_0*w40;
                                const register double tmp10_1 = A_12_0*w4;
                                const register double tmp536_1 = tmp133_0*w2;
                                const register double tmp55_1 = A_12_4*w27;
                                const register double tmp445_1 = tmp82_0*w17;
                                const register double tmp198_1 = A_02_2*w2;
                                const register double tmp27_1 = A_20_3*w23;
                                const register double tmp143_1 = A_10_5*w42;
                                const register double tmp165_1 = tmp96_0*w9;
                                const register double tmp145_1 = tmp78_0*w16;
                                const register double tmp282_1 = tmp1_0*w57;
                                const register double tmp13_1 = tmp7_0*w10;
                                const register double tmp129_1 = tmp66_0*w49;
                                const register double tmp443_1 = tmp65_0*w49;
                                const register double tmp135_1 = A_10_2*w46;
                                const register double tmp103_1 = tmp50_0*w21;
                                const register double tmp247_1 = tmp13_0*w36;
                                const register double tmp185_1 = tmp110_0*w16;
                                const register double tmp191_1 = tmp114_0*w18;
                                const register double tmp315_1 = tmp19_0*w17;
                                const register double tmp380_1 = tmp195_0*w23;
                                const register double tmp456_1 = tmp218_0*w27;
                                const register double tmp288_1 = tmp154_0*w29;
                                EM_S[INDEX2(7,3,8)]+=tmp0_1 + tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp1_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                EM_S[INDEX2(4,7,8)]+=tmp30_1 + tmp31_1 + tmp32_1 + tmp33_1 + tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                                EM_S[INDEX2(1,3,8)]+=tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1 + tmp60_1 + tmp61_1 + tmp62_1 + tmp63_1 + tmp64_1 + tmp65_1 + tmp66_1 + tmp67_1 + tmp68_1 + tmp69_1 + tmp70_1 + tmp71_1 + tmp72_1 + tmp73_1 + tmp74_1 + tmp75_1 + tmp76_1 + tmp77_1 + tmp78_1 + tmp79_1 + tmp80_1 + tmp81_1;
                                EM_S[INDEX2(6,4,8)]+=tmp100_1 + tmp101_1 + tmp102_1 + tmp103_1 + tmp104_1 + tmp105_1 + tmp59_1 + tmp63_1 + tmp67_1 + tmp74_1 + tmp75_1 + tmp80_1 + tmp82_1 + tmp83_1 + tmp84_1 + tmp85_1 + tmp86_1 + tmp87_1 + tmp88_1 + tmp89_1 + tmp90_1 + tmp91_1 + tmp92_1 + tmp93_1 + tmp94_1 + tmp95_1 + tmp96_1 + tmp97_1 + tmp98_1 + tmp99_1;
                                EM_S[INDEX2(3,0,8)]+=tmp106_1 + tmp107_1 + tmp108_1 + tmp109_1 + tmp110_1 + tmp111_1 + tmp112_1 + tmp113_1 + tmp114_1 + tmp115_1 + tmp116_1 + tmp117_1 + tmp118_1 + tmp119_1 + tmp120_1 + tmp121_1 + tmp35_1 + tmp38_1 + tmp41_1 + tmp42_1 + tmp47_1 + tmp50_1;
                                EM_S[INDEX2(5,4,8)]+=tmp122_1 + tmp123_1 + tmp124_1 + tmp125_1 + tmp126_1 + tmp127_1 + tmp128_1 + tmp129_1 + tmp130_1 + tmp131_1 + tmp132_1 + tmp133_1 + tmp134_1 + tmp135_1 + tmp136_1 + tmp137_1 + tmp138_1 + tmp139_1 + tmp140_1 + tmp141_1 + tmp142_1 + tmp143_1 + tmp144_1 + tmp145_1 + tmp146_1 + tmp147_1 + tmp148_1 + tmp149_1 + tmp150_1 + tmp151_1;
                                EM_S[INDEX2(0,7,8)]+=tmp132_1 + tmp152_1 + tmp153_1 + tmp154_1 + tmp155_1 + tmp156_1 + tmp157_1 + tmp158_1 + tmp159_1 + tmp160_1 + tmp161_1 + tmp74_1;
                                EM_S[INDEX2(5,6,8)]+=tmp162_1 + tmp163_1 + tmp164_1 + tmp165_1 + tmp166_1 + tmp167_1 + tmp168_1 + tmp169_1 + tmp170_1 + tmp171_1 + tmp172_1 + tmp31_1 + tmp34_1 + tmp35_1 + tmp37_1 + tmp39_1 + tmp41_1 + tmp43_1 + tmp45_1 + tmp46_1 + tmp49_1 + tmp50_1;
                                EM_S[INDEX2(2,6,8)]+=tmp104_1 + tmp158_1 + tmp16_1 + tmp173_1 + tmp174_1 + tmp175_1 + tmp176_1 + tmp177_1 + tmp178_1 + tmp179_1 + tmp180_1 + tmp181_1 + tmp182_1 + tmp183_1 + tmp184_1 + tmp185_1 + tmp186_1 + tmp187_1 + tmp188_1 + tmp19_1 + tmp24_1 + tmp28_1 + tmp2_1 + tmp59_1 + tmp86_1 + tmp88_1 + tmp8_1 + tmp91_1 + tmp92_1 + tmp99_1;
                                EM_S[INDEX2(1,6,8)]+=tmp132_1 + tmp155_1 + tmp156_1 + tmp157_1 + tmp159_1 + tmp161_1 + tmp189_1 + tmp18_1 + tmp190_1 + tmp191_1 + tmp192_1 + tmp193_1;
                                EM_S[INDEX2(5,1,8)]+=tmp13_1 + tmp158_1 + tmp182_1 + tmp184_1 + tmp185_1 + tmp186_1 + tmp194_1 + tmp195_1 + tmp196_1 + tmp197_1 + tmp198_1 + tmp199_1 + tmp200_1 + tmp201_1 + tmp202_1 + tmp203_1 + tmp204_1 + tmp205_1 + tmp206_1 + tmp207_1 + tmp208_1 + tmp209_1 + tmp4_1 + tmp56_1 + tmp58_1 + tmp59_1 + tmp62_1 + tmp64_1 + tmp72_1 + tmp79_1;
                                EM_S[INDEX2(3,7,8)]+=tmp0_1 + tmp12_1 + tmp139_1 + tmp13_1 + tmp18_1 + tmp19_1 + tmp1_1 + tmp20_1 + tmp210_1 + tmp211_1 + tmp212_1 + tmp213_1 + tmp214_1 + tmp215_1 + tmp216_1 + tmp217_1 + tmp218_1 + tmp219_1 + tmp220_1 + tmp221_1 + tmp23_1 + tmp2_1 + tmp4_1 + tmp55_1 + tmp62_1 + tmp64_1 + tmp72_1 + tmp79_1 + tmp80_1 + tmp81_1;
                                EM_S[INDEX2(2,5,8)]+=tmp153_1 + tmp155_1 + tmp159_1 + tmp160_1 + tmp161_1 + tmp18_1 + tmp222_1 + tmp223_1 + tmp224_1 + tmp225_1 + tmp226_1 + tmp74_1;
                                EM_S[INDEX2(0,3,8)]+=tmp107_1 + tmp109_1 + tmp110_1 + tmp115_1 + tmp116_1 + tmp117_1 + tmp166_1 + tmp169_1 + tmp227_1 + tmp228_1 + tmp229_1 + tmp230_1 + tmp231_1 + tmp232_1 + tmp233_1 + tmp234_1 + tmp235_1 + tmp236_1 + tmp237_1 + tmp238_1 + tmp239_1 + tmp41_1;
                                EM_S[INDEX2(7,2,8)]+=tmp112_1 + tmp113_1 + tmp240_1 + tmp241_1 + tmp242_1 + tmp243_1 + tmp244_1 + tmp245_1 + tmp246_1 + tmp247_1 + tmp248_1 + tmp249_1 + tmp250_1 + tmp251_1 + tmp252_1 + tmp253_1 + tmp254_1 + tmp255_1 + tmp35_1 + tmp46_1 + tmp49_1 + tmp50_1;
                                EM_S[INDEX2(4,0,8)]+=tmp104_1 + tmp105_1 + tmp124_1 + tmp130_1 + tmp138_1 + tmp139_1 + tmp148_1 + tmp150_1 + tmp151_1 + tmp175_1 + tmp181_1 + tmp18_1 + tmp195_1 + tmp204_1 + tmp20_1 + tmp216_1 + tmp218_1 + tmp256_1 + tmp257_1 + tmp258_1 + tmp259_1 + tmp260_1 + tmp261_1 + tmp262_1 + tmp263_1 + tmp80_1 + tmp85_1 + tmp91_1 + tmp92_1 + tmp99_1;
                                EM_S[INDEX2(1,2,8)]+=tmp107_1 + tmp108_1 + tmp110_1 + tmp114_1 + tmp115_1 + tmp117_1 + tmp119_1 + tmp121_1 + tmp229_1 + tmp231_1 + tmp233_1 + tmp235_1 + tmp237_1 + tmp238_1 + tmp264_1 + tmp265_1 + tmp266_1 + tmp267_1 + tmp268_1 + tmp41_1 + tmp42_1 + tmp47_1;
                                EM_S[INDEX2(6,7,8)]+=tmp122_1 + tmp125_1 + tmp131_1 + tmp134_1 + tmp142_1 + tmp146_1 + tmp16_1 + tmp174_1 + tmp182_1 + tmp187_1 + tmp224_1 + tmp22_1 + tmp24_1 + tmp269_1 + tmp270_1 + tmp271_1 + tmp272_1 + tmp273_1 + tmp274_1 + tmp275_1 + tmp276_1 + tmp277_1 + tmp278_1 + tmp279_1 + tmp280_1 + tmp281_1 + tmp28_1 + tmp6_1 + tmp7_1 + tmp8_1;
                                EM_S[INDEX2(3,3,8)]+=tmp249_1 + tmp255_1 + tmp264_1 + tmp268_1 + tmp282_1 + tmp283_1 + tmp284_1 + tmp285_1 + tmp286_1 + tmp287_1 + tmp288_1 + tmp289_1 + tmp290_1 + tmp291_1 + tmp292_1 + tmp293_1 + tmp294_1 + tmp295_1 + tmp296_1 + tmp297_1 + tmp298_1 + tmp299_1 + tmp300_1 + tmp301_1;
                                EM_S[INDEX2(2,0,8)]+=tmp102_1 + tmp11_1 + tmp193_1 + tmp302_1 + tmp303_1 + tmp304_1 + tmp305_1 + tmp306_1 + tmp307_1 + tmp308_1 + tmp309_1 + tmp310_1 + tmp311_1 + tmp312_1 + tmp313_1 + tmp314_1 + tmp315_1 + tmp316_1 + tmp317_1 + tmp318_1 + tmp319_1 + tmp320_1 + tmp321_1 + tmp52_1 + tmp54_1 + tmp57_1 + tmp68_1 + tmp70_1 + tmp76_1 + tmp94_1;
                                EM_S[INDEX2(7,6,8)]+=tmp101_1 + tmp125_1 + tmp134_1 + tmp211_1 + tmp212_1 + tmp213_1 + tmp217_1 + tmp218_1 + tmp219_1 + tmp221_1 + tmp224_1 + tmp270_1 + tmp271_1 + tmp274_1 + tmp276_1 + tmp277_1 + tmp278_1 + tmp281_1 + tmp310_1 + tmp322_1 + tmp323_1 + tmp324_1 + tmp325_1 + tmp326_1 + tmp67_1 + tmp82_1 + tmp87_1 + tmp90_1 + tmp97_1 + tmp98_1;
                                EM_S[INDEX2(4,4,8)]+=tmp162_1 + tmp170_1 + tmp284_1 + tmp287_1 + tmp290_1 + tmp292_1 + tmp293_1 + tmp298_1 + tmp327_1 + tmp328_1 + tmp329_1 + tmp330_1 + tmp331_1 + tmp332_1 + tmp333_1 + tmp334_1 + tmp335_1 + tmp336_1 + tmp337_1 + tmp338_1 + tmp339_1 + tmp340_1 + tmp341_1 + tmp342_1;
                                EM_S[INDEX2(6,3,8)]+=tmp112_1 + tmp113_1 + tmp241_1 + tmp242_1 + tmp246_1 + tmp247_1 + tmp252_1 + tmp343_1 + tmp344_1 + tmp345_1 + tmp346_1 + tmp347_1 + tmp348_1 + tmp349_1 + tmp350_1 + tmp351_1 + tmp352_1 + tmp353_1 + tmp35_1 + tmp46_1 + tmp49_1 + tmp50_1;
                                EM_S[INDEX2(1,5,8)]+=tmp127_1 + tmp128_1 + tmp130_1 + tmp138_1 + tmp13_1 + tmp145_1 + tmp148_1 + tmp14_1 + tmp151_1 + tmp158_1 + tmp15_1 + tmp184_1 + tmp194_1 + tmp195_1 + tmp201_1 + tmp202_1 + tmp204_1 + tmp205_1 + tmp21_1 + tmp25_1 + tmp319_1 + tmp325_1 + tmp354_1 + tmp355_1 + tmp356_1 + tmp357_1 + tmp358_1 + tmp359_1 + tmp360_1 + tmp4_1;
                                EM_S[INDEX2(3,6,8)]+=tmp229_1 + tmp231_1 + tmp233_1 + tmp238_1 + tmp240_1 + tmp241_1 + tmp242_1 + tmp245_1 + tmp246_1 + tmp247_1 + tmp248_1 + tmp250_1 + tmp251_1 + tmp252_1 + tmp254_1 + tmp345_1 + tmp353_1 + tmp361_1 + tmp362_1 + tmp363_1 + tmp364_1 + tmp365_1;
                                EM_S[INDEX2(2,2,8)]+=tmp109_1 + tmp116_1 + tmp283_1 + tmp287_1 + tmp295_1 + tmp338_1 + tmp340_1 + tmp345_1 + tmp353_1 + tmp366_1 + tmp367_1 + tmp368_1 + tmp369_1 + tmp370_1 + tmp371_1 + tmp372_1 + tmp373_1 + tmp374_1 + tmp375_1 + tmp376_1 + tmp377_1 + tmp378_1 + tmp379_1 + tmp380_1;
                                EM_S[INDEX2(7,7,8)]+=tmp162_1 + tmp170_1 + tmp282_1 + tmp292_1 + tmp301_1 + tmp345_1 + tmp353_1 + tmp369_1 + tmp381_1 + tmp382_1 + tmp383_1 + tmp384_1 + tmp385_1 + tmp386_1 + tmp387_1 + tmp388_1 + tmp389_1 + tmp390_1 + tmp391_1 + tmp392_1 + tmp393_1 + tmp394_1 + tmp395_1 + tmp396_1;
                                EM_S[INDEX2(5,7,8)]+=tmp101_1 + tmp10_1 + tmp11_1 + tmp14_1 + tmp15_1 + tmp193_1 + tmp21_1 + tmp25_1 + tmp310_1 + tmp312_1 + tmp317_1 + tmp319_1 + tmp322_1 + tmp324_1 + tmp355_1 + tmp360_1 + tmp397_1 + tmp398_1 + tmp399_1 + tmp400_1 + tmp401_1 + tmp402_1 + tmp66_1 + tmp77_1 + tmp82_1 + tmp84_1 + tmp87_1 + tmp95_1 + tmp97_1 + tmp9_1;
                                EM_S[INDEX2(5,3,8)]+=tmp108_1 + tmp121_1 + tmp245_1 + tmp248_1 + tmp250_1 + tmp251_1 + tmp387_1 + tmp388_1 + tmp403_1 + tmp404_1 + tmp405_1 + tmp406_1 + tmp407_1 + tmp408_1 + tmp409_1 + tmp40_1 + tmp410_1 + tmp411_1 + tmp412_1 + tmp42_1 + tmp47_1 + tmp48_1;
                                EM_S[INDEX2(4,1,8)]+=tmp229_1 + tmp235_1 + tmp237_1 + tmp238_1 + tmp242_1 + tmp248_1 + tmp250_1 + tmp362_1 + tmp403_1 + tmp410_1 + tmp413_1 + tmp414_1 + tmp415_1 + tmp416_1 + tmp417_1 + tmp418_1 + tmp419_1 + tmp420_1 + tmp421_1 + tmp422_1 + tmp423_1 + tmp424_1;
                                EM_S[INDEX2(1,1,8)]+=tmp109_1 + tmp116_1 + tmp284_1 + tmp294_1 + tmp300_1 + tmp334_1 + tmp339_1 + tmp375_1 + tmp384_1 + tmp387_1 + tmp388_1 + tmp425_1 + tmp426_1 + tmp427_1 + tmp428_1 + tmp429_1 + tmp430_1 + tmp431_1 + tmp432_1 + tmp433_1 + tmp434_1 + tmp435_1 + tmp436_1 + tmp437_1;
                                EM_S[INDEX2(2,7,8)]+=tmp229_1 + tmp231_1 + tmp233_1 + tmp238_1 + tmp241_1 + tmp242_1 + tmp246_1 + tmp247_1 + tmp249_1 + tmp252_1 + tmp255_1 + tmp343_1 + tmp346_1 + tmp347_1 + tmp349_1 + tmp350_1 + tmp351_1 + tmp363_1 + tmp364_1 + tmp438_1 + tmp439_1 + tmp440_1;
                                EM_S[INDEX2(3,2,8)]+=tmp132_1 + tmp137_1 + tmp16_1 + tmp177_1 + tmp178_1 + tmp17_1 + tmp185_1 + tmp24_1 + tmp278_1 + tmp27_1 + tmp281_1 + tmp28_1 + tmp308_1 + tmp312_1 + tmp314_1 + tmp3_1 + tmp441_1 + tmp442_1 + tmp443_1 + tmp444_1 + tmp445_1 + tmp446_1 + tmp52_1 + tmp53_1 + tmp57_1 + tmp60_1 + tmp63_1 + tmp70_1 + tmp76_1 + tmp8_1;
                                EM_S[INDEX2(0,0,8)]+=tmp264_1 + tmp268_1 + tmp292_1 + tmp327_1 + tmp342_1 + tmp369_1 + tmp377_1 + tmp379_1 + tmp384_1 + tmp389_1 + tmp390_1 + tmp394_1 + tmp415_1 + tmp424_1 + tmp427_1 + tmp434_1 + tmp447_1 + tmp448_1 + tmp449_1 + tmp450_1 + tmp451_1 + tmp452_1 + tmp453_1 + tmp454_1;
                                EM_S[INDEX2(6,6,8)]+=tmp249_1 + tmp255_1 + tmp284_1 + tmp335_1 + tmp33_1 + tmp341_1 + tmp366_1 + tmp375_1 + tmp378_1 + tmp381_1 + tmp384_1 + tmp392_1 + tmp431_1 + tmp433_1 + tmp435_1 + tmp44_1 + tmp450_1 + tmp451_1 + tmp455_1 + tmp456_1 + tmp457_1 + tmp458_1 + tmp459_1 + tmp460_1;
                                EM_S[INDEX2(5,0,8)]+=tmp229_1 + tmp235_1 + tmp237_1 + tmp238_1 + tmp242_1 + tmp334_1 + tmp339_1 + tmp347_1 + tmp349_1 + tmp414_1 + tmp417_1 + tmp418_1 + tmp419_1 + tmp420_1 + tmp422_1 + tmp440_1 + tmp461_1 + tmp462_1 + tmp463_1 + tmp464_1 + tmp465_1 + tmp466_1;
                                EM_S[INDEX2(7,1,8)]+=tmp108_1 + tmp121_1 + tmp297_1 + tmp299_1 + tmp346_1 + tmp347_1 + tmp349_1 + tmp350_1 + tmp404_1 + tmp405_1 + tmp407_1 + tmp408_1 + tmp409_1 + tmp40_1 + tmp42_1 + tmp461_1 + tmp466_1 + tmp467_1 + tmp468_1 + tmp469_1 + tmp47_1 + tmp48_1;
                                EM_S[INDEX2(4,5,8)]+=tmp101_1 + tmp125_1 + tmp126_1 + tmp129_1 + tmp132_1 + tmp134_1 + tmp137_1 + tmp140_1 + tmp144_1 + tmp147_1 + tmp149_1 + tmp17_1 + tmp185_1 + tmp198_1 + tmp199_1 + tmp200_1 + tmp203_1 + tmp206_1 + tmp208_1 + tmp312_1 + tmp398_1 + tmp400_1 + tmp470_1 + tmp471_1 + tmp63_1 + tmp82_1 + tmp83_1 + tmp87_1 + tmp89_1 + tmp97_1;
                                EM_S[INDEX2(0,4,8)]+=tmp11_1 + tmp175_1 + tmp17_1 + tmp181_1 + tmp18_1 + tmp195_1 + tmp200_1 + tmp203_1 + tmp204_1 + tmp206_1 + tmp208_1 + tmp20_1 + tmp22_1 + tmp256_1 + tmp257_1 + tmp262_1 + tmp263_1 + tmp26_1 + tmp305_1 + tmp306_1 + tmp309_1 + tmp311_1 + tmp315_1 + tmp318_1 + tmp470_1 + tmp471_1 + tmp472_1 + tmp473_1 + tmp474_1 + tmp475_1;
                                EM_S[INDEX2(5,5,8)]+=tmp287_1 + tmp297_1 + tmp299_1 + tmp328_1 + tmp336_1 + tmp33_1 + tmp369_1 + tmp372_1 + tmp374_1 + tmp375_1 + tmp376_1 + tmp385_1 + tmp391_1 + tmp415_1 + tmp424_1 + tmp425_1 + tmp436_1 + tmp44_1 + tmp476_1 + tmp477_1 + tmp478_1 + tmp479_1 + tmp480_1 + tmp481_1;
                                EM_S[INDEX2(1,4,8)]+=tmp118_1 + tmp120_1 + tmp242_1 + tmp347_1 + tmp348_1 + tmp349_1 + tmp35_1 + tmp37_1 + tmp39_1 + tmp414_1 + tmp415_1 + tmp417_1 + tmp418_1 + tmp422_1 + tmp424_1 + tmp461_1 + tmp464_1 + tmp465_1 + tmp466_1 + tmp482_1 + tmp483_1 + tmp50_1;
                                EM_S[INDEX2(6,0,8)]+=tmp163_1 + tmp166_1 + tmp169_1 + tmp172_1 + tmp234_1 + tmp236_1 + tmp240_1 + tmp248_1 + tmp250_1 + tmp254_1 + tmp338_1 + tmp340_1 + tmp407_1 + tmp416_1 + tmp421_1 + tmp484_1 + tmp485_1 + tmp486_1 + tmp487_1 + tmp488_1 + tmp489_1 + tmp490_1;
                                EM_S[INDEX2(7,5,8)]+=tmp122_1 + tmp131_1 + tmp135_1 + tmp141_1 + tmp142_1 + tmp143_1 + tmp146_1 + tmp186_1 + tmp193_1 + tmp197_1 + tmp209_1 + tmp214_1 + tmp215_1 + tmp216_1 + tmp269_1 + tmp272_1 + tmp275_1 + tmp317_1 + tmp397_1 + tmp399_1 + tmp401_1 + tmp402_1 + tmp62_1 + tmp64_1 + tmp66_1 + tmp72_1 + tmp77_1 + tmp79_1 + tmp84_1 + tmp95_1;
                                EM_S[INDEX2(2,3,8)]+=tmp132_1 + tmp136_1 + tmp137_1 + tmp139_1 + tmp141_1 + tmp145_1 + tmp210_1 + tmp213_1 + tmp217_1 + tmp219_1 + tmp220_1 + tmp221_1 + tmp278_1 + tmp281_1 + tmp441_1 + tmp442_1 + tmp443_1 + tmp444_1 + tmp445_1 + tmp446_1 + tmp491_1 + tmp492_1 + tmp493_1 + tmp494_1 + tmp495_1 + tmp496_1 + tmp497_1 + tmp498_1 + tmp499_1 + tmp500_1;
                                EM_S[INDEX2(2,1,8)]+=tmp107_1 + tmp110_1 + tmp112_1 + tmp113_1 + tmp115_1 + tmp117_1 + tmp118_1 + tmp120_1 + tmp166_1 + tmp167_1 + tmp169_1 + tmp228_1 + tmp234_1 + tmp236_1 + tmp239_1 + tmp264_1 + tmp268_1 + tmp35_1 + tmp41_1 + tmp501_1 + tmp502_1 + tmp50_1;
                                EM_S[INDEX2(4,2,8)]+=tmp163_1 + tmp166_1 + tmp169_1 + tmp172_1 + tmp234_1 + tmp236_1 + tmp343_1 + tmp347_1 + tmp349_1 + tmp351_1 + tmp407_1 + tmp450_1 + tmp451_1 + tmp464_1 + tmp465_1 + tmp484_1 + tmp485_1 + tmp487_1 + tmp490_1 + tmp503_1 + tmp504_1 + tmp505_1;
                                EM_S[INDEX2(1,0,8)]+=tmp140_1 + tmp147_1 + tmp182_1 + tmp196_1 + tmp200_1 + tmp203_1 + tmp206_1 + tmp207_1 + tmp208_1 + tmp224_1 + tmp22_1 + tmp275_1 + tmp277_1 + tmp279_1 + tmp441_1 + tmp444_1 + tmp473_1 + tmp474_1 + tmp491_1 + tmp495_1 + tmp498_1 + tmp500_1 + tmp506_1 + tmp507_1 + tmp508_1 + tmp509_1 + tmp510_1 + tmp511_1 + tmp512_1 + tmp513_1;
                                EM_S[INDEX2(6,5,8)]+=tmp162_1 + tmp170_1 + tmp229_1 + tmp238_1 + tmp266_1 + tmp31_1 + tmp32_1 + tmp34_1 + tmp363_1 + tmp364_1 + tmp40_1 + tmp419_1 + tmp41_1 + tmp420_1 + tmp42_1 + tmp43_1 + tmp45_1 + tmp47_1 + tmp48_1 + tmp514_1 + tmp515_1 + tmp51_1;
                                EM_S[INDEX2(3,5,8)]+=tmp165_1 + tmp166_1 + tmp169_1 + tmp171_1 + tmp228_1 + tmp239_1 + tmp346_1 + tmp347_1 + tmp349_1 + tmp350_1 + tmp387_1 + tmp388_1 + tmp404_1 + tmp405_1 + tmp407_1 + tmp408_1 + tmp409_1 + tmp461_1 + tmp466_1 + tmp504_1 + tmp516_1 + tmp517_1;
                                EM_S[INDEX2(0,1,8)]+=tmp130_1 + tmp138_1 + tmp140_1 + tmp147_1 + tmp148_1 + tmp151_1 + tmp218_1 + tmp224_1 + tmp258_1 + tmp259_1 + tmp277_1 + tmp302_1 + tmp307_1 + tmp310_1 + tmp325_1 + tmp354_1 + tmp359_1 + tmp441_1 + tmp444_1 + tmp507_1 + tmp508_1 + tmp511_1 + tmp512_1 + tmp52_1 + tmp57_1 + tmp61_1 + tmp67_1 + tmp70_1 + tmp71_1 + tmp76_1;
                                EM_S[INDEX2(7,0,8)]+=tmp132_1 + tmp155_1 + tmp158_1 + tmp159_1 + tmp161_1 + tmp518_1 + tmp519_1 + tmp520_1 + tmp521_1 + tmp522_1 + tmp523_1 + tmp74_1;
                                EM_S[INDEX2(4,6,8)]+=tmp100_1 + tmp102_1 + tmp103_1 + tmp122_1 + tmp123_1 + tmp131_1 + tmp133_1 + tmp136_1 + tmp142_1 + tmp146_1 + tmp26_1 + tmp273_1 + tmp279_1 + tmp280_1 + tmp309_1 + tmp311_1 + tmp315_1 + tmp318_1 + tmp358_1 + tmp472_1 + tmp475_1 + tmp524_1 + tmp525_1 + tmp74_1 + tmp75_1 + tmp84_1 + tmp93_1 + tmp94_1 + tmp95_1 + tmp96_1;
                                EM_S[INDEX2(5,2,8)]+=tmp155_1 + tmp159_1 + tmp161_1 + tmp189_1 + tmp18_1 + tmp190_1 + tmp224_1 + tmp519_1 + tmp523_1 + tmp526_1 + tmp527_1 + tmp74_1;
                                EM_S[INDEX2(6,1,8)]+=tmp132_1 + tmp155_1 + tmp159_1 + tmp161_1 + tmp18_1 + tmp193_1 + tmp222_1 + tmp223_1 + tmp521_1 + tmp522_1 + tmp528_1 + tmp529_1;
                                EM_S[INDEX2(3,1,8)]+=tmp136_1 + tmp14_1 + tmp15_1 + tmp21_1 + tmp25_1 + tmp26_1 + tmp279_1 + tmp29_1 + tmp356_1 + tmp357_1 + tmp358_1 + tmp491_1 + tmp492_1 + tmp495_1 + tmp496_1 + tmp498_1 + tmp500_1 + tmp510_1 + tmp513_1 + tmp54_1 + tmp5_1 + tmp65_1 + tmp66_1 + tmp68_1 + tmp69_1 + tmp73_1 + tmp74_1 + tmp75_1 + tmp77_1 + tmp78_1;
                                EM_S[INDEX2(0,2,8)]+=tmp102_1 + tmp104_1 + tmp141_1 + tmp176_1 + tmp186_1 + tmp188_1 + tmp193_1 + tmp216_1 + tmp260_1 + tmp261_1 + tmp275_1 + tmp304_1 + tmp313_1 + tmp316_1 + tmp317_1 + tmp320_1 + tmp491_1 + tmp495_1 + tmp497_1 + tmp498_1 + tmp499_1 + tmp500_1 + tmp506_1 + tmp509_1 + tmp54_1 + tmp68_1 + tmp91_1 + tmp92_1 + tmp94_1 + tmp99_1;
                                EM_S[INDEX2(7,4,8)]+=tmp163_1 + tmp165_1 + tmp166_1 + tmp169_1 + tmp171_1 + tmp172_1 + tmp229_1 + tmp232_1 + tmp238_1 + tmp31_1 + tmp33_1 + tmp34_1 + tmp363_1 + tmp364_1 + tmp419_1 + tmp41_1 + tmp420_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp530_1 + tmp531_1;
                                EM_S[INDEX2(0,6,8)]+=tmp114_1 + tmp119_1 + tmp32_1 + tmp338_1 + tmp340_1 + tmp343_1 + tmp347_1 + tmp349_1 + tmp351_1 + tmp407_1 + tmp42_1 + tmp464_1 + tmp465_1 + tmp469_1 + tmp47_1 + tmp484_1 + tmp485_1 + tmp487_1 + tmp490_1 + tmp51_1 + tmp532_1 + tmp533_1;
                                EM_S[INDEX2(6,2,8)]+=tmp145_1 + tmp158_1 + tmp173_1 + tmp175_1 + tmp179_1 + tmp180_1 + tmp181_1 + tmp183_1 + tmp184_1 + tmp19_1 + tmp213_1 + tmp217_1 + tmp219_1 + tmp221_1 + tmp2_1 + tmp303_1 + tmp309_1 + tmp311_1 + tmp315_1 + tmp318_1 + tmp319_1 + tmp321_1 + tmp323_1 + tmp325_1 + tmp326_1 + tmp358_1 + tmp493_1 + tmp494_1 + tmp524_1 + tmp525_1;
                                EM_S[INDEX2(4,3,8)]+=tmp152_1 + tmp154_1 + tmp155_1 + tmp158_1 + tmp159_1 + tmp161_1 + tmp193_1 + tmp224_1 + tmp526_1 + tmp527_1 + tmp528_1 + tmp529_1;
                                EM_S[INDEX2(1,7,8)]+=tmp165_1 + tmp166_1 + tmp169_1 + tmp171_1 + tmp228_1 + tmp239_1 + tmp245_1 + tmp248_1 + tmp250_1 + tmp251_1 + tmp297_1 + tmp299_1 + tmp403_1 + tmp404_1 + tmp405_1 + tmp407_1 + tmp408_1 + tmp409_1 + tmp410_1 + tmp489_1 + tmp534_1 + tmp535_1;
                                EM_S[INDEX2(0,5,8)]+=tmp118_1 + tmp120_1 + tmp242_1 + tmp248_1 + tmp250_1 + tmp253_1 + tmp334_1 + tmp339_1 + tmp35_1 + tmp37_1 + tmp39_1 + tmp403_1 + tmp410_1 + tmp414_1 + tmp416_1 + tmp417_1 + tmp418_1 + tmp421_1 + tmp422_1 + tmp50_1 + tmp536_1 + tmp537_1;
                                EM_S[INDEX2(3,4,8)]+=tmp155_1 + tmp158_1 + tmp159_1 + tmp161_1 + tmp191_1 + tmp192_1 + tmp193_1 + tmp224_1 + tmp225_1 + tmp226_1 + tmp518_1 + tmp520_1;
                                EM_S[INDEX2(2,4,8)]+=tmp114_1 + tmp119_1 + tmp240_1 + tmp248_1 + tmp250_1 + tmp254_1 + tmp32_1 + tmp407_1 + tmp411_1 + tmp416_1 + tmp421_1 + tmp42_1 + tmp450_1 + tmp451_1 + tmp47_1 + tmp484_1 + tmp485_1 + tmp487_1 + tmp490_1 + tmp51_1 + tmp538_1 + tmp539_1;
                            } else { // constant data
                                const register double A_00 = A_p[INDEX2(0,0,3)];
                                const register double A_01 = A_p[INDEX2(0,1,3)];
                                const register double A_02 = A_p[INDEX2(0,2,3)];
                                const register double A_10 = A_p[INDEX2(1,0,3)];
                                const register double A_11 = A_p[INDEX2(1,1,3)];
                                const register double A_12 = A_p[INDEX2(1,2,3)];
                                const register double A_20 = A_p[INDEX2(2,0,3)];
                                const register double A_21 = A_p[INDEX2(2,1,3)];
                                const register double A_22 = A_p[INDEX2(2,2,3)];
                                const register double tmp0_0 = A_01 + A_10;
                                const register double tmp1_0 = A_02 + A_20;
                                const register double tmp2_0 = A_12 + A_21;
                                const register double tmp25_1 = A_01*w69;
                                const register double tmp2_1 = tmp0_0*w61;
                                const register double tmp33_1 = A_20*w70;
                                const register double tmp23_1 = A_02*w65;
                                const register double tmp41_1 = A_01*w61;
                                const register double tmp34_1 = A_02*w73;
                                const register double tmp8_1 = A_11*w71;
                                const register double tmp50_1 = A_10*w61;
                                const register double tmp15_1 = A_22*w75;
                                const register double tmp9_1 = A_21*w74;
                                const register double tmp19_1 = A_10*w69;
                                const register double tmp11_1 = A_00*w68;
                                const register double tmp52_1 = tmp2_0*w66;
                                const register double tmp37_1 = tmp2_0*w74;
                                const register double tmp0_1 = A_00*w60;
                                const register double tmp17_1 = A_21*w64;
                                const register double tmp26_1 = A_00*w79;
                                const register double tmp5_1 = A_21*w66;
                                const register double tmp29_1 = A_00*w80;
                                const register double tmp7_1 = A_22*w67;
                                const register double tmp48_1 = A_11*w87;
                                const register double tmp44_1 = A_11*w84;
                                const register double tmp27_1 = tmp2_0*w72;
                                const register double tmp42_1 = A_22*w85;
                                const register double tmp18_1 = A_11*w77;
                                const register double tmp35_1 = tmp0_0*w76;
                                const register double tmp46_1 = A_00*w86;
                                const register double tmp32_1 = A_22*w83;
                                const register double tmp22_1 = A_01*w76;
                                const register double tmp4_1 = A_02*w62;
                                const register double tmp10_1 = A_02*w70;
                                const register double tmp3_1 = A_20*w65;
                                const register double tmp39_1 = A_21*w72;
                                const register double tmp51_1 = tmp1_0*w65;
                                const register double tmp12_1 = A_20*w73;
                                const register double tmp40_1 = A_10*w81;
                                const register double tmp43_1 = tmp1_0*w62;
                                const register double tmp28_1 = A_10*w76;
                                const register double tmp45_1 = tmp2_0*w64;
                                const register double tmp49_1 = A_01*w81;
                                const register double tmp36_1 = tmp1_0*w73;
                                const register double tmp53_1 = A_00*w89;
                                const register double tmp6_1 = A_11*w63;
                                const register double tmp31_1 = A_11*w82;
                                const register double tmp21_1 = A_22*w78;
                                const register double tmp16_1 = tmp1_0*w70;
                                const register double tmp14_1 = A_12*w72;
                                const register double tmp38_1 = A_12*w74;
                                const register double tmp30_1 = tmp0_0*w81;
                                const register double tmp47_1 = A_22*w88;
                                const register double tmp13_1 = tmp0_0*w69;
                                const register double tmp20_1 = A_12*w66;
                                const register double tmp1_1 = A_12*w64;
                                const register double tmp24_1 = A_20*w62;
                                EM_S[INDEX2(7,3,8)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1;
                                EM_S[INDEX2(4,7,8)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp8_1 + tmp9_1;
                                EM_S[INDEX2(1,3,8)]+=tmp0_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1;
                                EM_S[INDEX2(6,4,8)]+=tmp0_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1;
                                EM_S[INDEX2(3,0,8)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp8_1 + tmp9_1;
                                EM_S[INDEX2(5,4,8)]+=tmp21_1 + tmp23_1 + tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp6_1;
                                EM_S[INDEX2(0,7,8)]+=tmp16_1 + tmp27_1 + tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1;
                                EM_S[INDEX2(5,6,8)]+=tmp11_1 + tmp14_1 + tmp15_1 + tmp33_1 + tmp34_1 + tmp35_1 + tmp8_1 + tmp9_1;
                                EM_S[INDEX2(2,6,8)]+=tmp0_1 + tmp17_1 + tmp20_1 + tmp30_1 + tmp3_1 + tmp4_1 + tmp6_1 + tmp7_1;
                                EM_S[INDEX2(1,6,8)]+=tmp27_1 + tmp29_1 + tmp2_1 + tmp31_1 + tmp32_1 + tmp36_1;
                                EM_S[INDEX2(5,1,8)]+=tmp0_1 + tmp17_1 + tmp20_1 + tmp30_1 + tmp3_1 + tmp4_1 + tmp6_1 + tmp7_1;
                                EM_S[INDEX2(3,7,8)]+=tmp0_1 + tmp17_1 + tmp20_1 + tmp23_1 + tmp24_1 + tmp2_1 + tmp6_1 + tmp7_1;
                                EM_S[INDEX2(2,5,8)]+=tmp16_1 + tmp29_1 + tmp2_1 + tmp31_1 + tmp32_1 + tmp37_1;
                                EM_S[INDEX2(0,3,8)]+=tmp11_1 + tmp13_1 + tmp15_1 + tmp33_1 + tmp34_1 + tmp38_1 + tmp39_1 + tmp8_1;
                                EM_S[INDEX2(7,2,8)]+=tmp11_1 + tmp14_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp9_1;
                                EM_S[INDEX2(4,0,8)]+=tmp0_1 + tmp17_1 + tmp20_1 + tmp23_1 + tmp24_1 + tmp2_1 + tmp6_1 + tmp7_1;
                                EM_S[INDEX2(1,2,8)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp15_1 + tmp35_1 + tmp38_1 + tmp39_1 + tmp8_1;
                                EM_S[INDEX2(6,7,8)]+=tmp21_1 + tmp25_1 + tmp26_1 + tmp28_1 + tmp37_1 + tmp3_1 + tmp4_1 + tmp6_1;
                                EM_S[INDEX2(3,3,8)]+=tmp35_1 + tmp43_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp48_1;
                                EM_S[INDEX2(2,0,8)]+=tmp0_1 + tmp18_1 + tmp19_1 + tmp1_1 + tmp21_1 + tmp22_1 + tmp36_1 + tmp5_1;
                                EM_S[INDEX2(7,6,8)]+=tmp19_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp24_1 + tmp26_1 + tmp37_1 + tmp6_1;
                                EM_S[INDEX2(4,4,8)]+=tmp35_1 + tmp43_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp48_1;
                                EM_S[INDEX2(6,3,8)]+=tmp11_1 + tmp14_1 + tmp42_1 + tmp44_1 + tmp49_1 + tmp50_1 + tmp51_1 + tmp9_1;
                                EM_S[INDEX2(1,5,8)]+=tmp0_1 + tmp1_1 + tmp23_1 + tmp24_1 + tmp30_1 + tmp5_1 + tmp6_1 + tmp7_1;
                                EM_S[INDEX2(3,6,8)]+=tmp11_1 + tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp44_1 + tmp51_1;
                                EM_S[INDEX2(2,2,8)]+=tmp13_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp51_1;
                                EM_S[INDEX2(7,7,8)]+=tmp35_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp51_1 + tmp52_1;
                                EM_S[INDEX2(5,7,8)]+=tmp0_1 + tmp18_1 + tmp19_1 + tmp1_1 + tmp21_1 + tmp22_1 + tmp36_1 + tmp5_1;
                                EM_S[INDEX2(5,3,8)]+=tmp10_1 + tmp12_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp52_1 + tmp53_1 + tmp8_1;
                                EM_S[INDEX2(4,1,8)]+=tmp11_1 + tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp44_1 + tmp51_1;
                                EM_S[INDEX2(1,1,8)]+=tmp13_1 + tmp43_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp52_1;
                                EM_S[INDEX2(2,7,8)]+=tmp11_1 + tmp38_1 + tmp39_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp49_1 + tmp50_1;
                                EM_S[INDEX2(3,2,8)]+=tmp19_1 + tmp21_1 + tmp22_1 + tmp26_1 + tmp27_1 + tmp3_1 + tmp4_1 + tmp6_1;
                                EM_S[INDEX2(0,0,8)]+=tmp35_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp51_1 + tmp52_1;
                                EM_S[INDEX2(6,6,8)]+=tmp13_1 + tmp43_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp52_1;
                                EM_S[INDEX2(5,0,8)]+=tmp11_1 + tmp38_1 + tmp39_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp49_1 + tmp50_1;
                                EM_S[INDEX2(7,1,8)]+=tmp10_1 + tmp12_1 + tmp42_1 + tmp45_1 + tmp49_1 + tmp50_1 + tmp53_1 + tmp8_1;
                                EM_S[INDEX2(4,5,8)]+=tmp19_1 + tmp21_1 + tmp22_1 + tmp26_1 + tmp27_1 + tmp3_1 + tmp4_1 + tmp6_1;
                                EM_S[INDEX2(0,4,8)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1;
                                EM_S[INDEX2(5,5,8)]+=tmp13_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp51_1;
                                EM_S[INDEX2(1,4,8)]+=tmp11_1 + tmp14_1 + tmp42_1 + tmp44_1 + tmp49_1 + tmp50_1 + tmp51_1 + tmp9_1;
                                EM_S[INDEX2(6,0,8)]+=tmp33_1 + tmp34_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp45_1 + tmp53_1 + tmp8_1;
                                EM_S[INDEX2(7,5,8)]+=tmp0_1 + tmp17_1 + tmp18_1 + tmp20_1 + tmp21_1 + tmp25_1 + tmp28_1 + tmp36_1;
                                EM_S[INDEX2(2,3,8)]+=tmp21_1 + tmp23_1 + tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp6_1;
                                EM_S[INDEX2(2,1,8)]+=tmp11_1 + tmp14_1 + tmp15_1 + tmp33_1 + tmp34_1 + tmp35_1 + tmp8_1 + tmp9_1;
                                EM_S[INDEX2(4,2,8)]+=tmp33_1 + tmp34_1 + tmp42_1 + tmp49_1 + tmp50_1 + tmp52_1 + tmp53_1 + tmp8_1;
                                EM_S[INDEX2(1,0,8)]+=tmp21_1 + tmp25_1 + tmp26_1 + tmp28_1 + tmp37_1 + tmp3_1 + tmp4_1 + tmp6_1;
                                EM_S[INDEX2(6,5,8)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp15_1 + tmp35_1 + tmp38_1 + tmp39_1 + tmp8_1;
                                EM_S[INDEX2(3,5,8)]+=tmp33_1 + tmp34_1 + tmp42_1 + tmp49_1 + tmp50_1 + tmp52_1 + tmp53_1 + tmp8_1;
                                EM_S[INDEX2(0,1,8)]+=tmp19_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp24_1 + tmp26_1 + tmp37_1 + tmp6_1;
                                EM_S[INDEX2(7,0,8)]+=tmp16_1 + tmp27_1 + tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1;
                                EM_S[INDEX2(4,6,8)]+=tmp0_1 + tmp16_1 + tmp18_1 + tmp1_1 + tmp21_1 + tmp25_1 + tmp28_1 + tmp5_1;
                                EM_S[INDEX2(5,2,8)]+=tmp16_1 + tmp29_1 + tmp2_1 + tmp31_1 + tmp32_1 + tmp37_1;
                                EM_S[INDEX2(6,1,8)]+=tmp27_1 + tmp29_1 + tmp2_1 + tmp31_1 + tmp32_1 + tmp36_1;
                                EM_S[INDEX2(3,1,8)]+=tmp0_1 + tmp16_1 + tmp18_1 + tmp1_1 + tmp21_1 + tmp25_1 + tmp28_1 + tmp5_1;
                                EM_S[INDEX2(0,2,8)]+=tmp0_1 + tmp17_1 + tmp18_1 + tmp20_1 + tmp21_1 + tmp25_1 + tmp28_1 + tmp36_1;
                                EM_S[INDEX2(7,4,8)]+=tmp11_1 + tmp13_1 + tmp15_1 + tmp33_1 + tmp34_1 + tmp38_1 + tmp39_1 + tmp8_1;
                                EM_S[INDEX2(0,6,8)]+=tmp10_1 + tmp12_1 + tmp42_1 + tmp45_1 + tmp49_1 + tmp50_1 + tmp53_1 + tmp8_1;
                                EM_S[INDEX2(6,2,8)]+=tmp0_1 + tmp1_1 + tmp23_1 + tmp24_1 + tmp30_1 + tmp5_1 + tmp6_1 + tmp7_1;
                                EM_S[INDEX2(4,3,8)]+=tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp36_1 + tmp37_1;
                                EM_S[INDEX2(1,7,8)]+=tmp33_1 + tmp34_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp45_1 + tmp53_1 + tmp8_1;
                                EM_S[INDEX2(0,5,8)]+=tmp11_1 + tmp14_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp9_1;
                                EM_S[INDEX2(3,4,8)]+=tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp36_1 + tmp37_1;
                                EM_S[INDEX2(2,4,8)]+=tmp10_1 + tmp12_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp52_1 + tmp53_1 + tmp8_1;
                            }
                        }
                        ///////////////
                        // process B //
                        ///////////////
                        if (!B.isEmpty()) {
                            add_EM_S=true;
                            const double* B_p=const_cast<escript::Data*>(&B)->getSampleDataRO(e);
                            if (B.actsExpanded()) {
                                const register double B_0_0 = B_p[INDEX2(0,0,3)];
                                const register double B_1_0 = B_p[INDEX2(1,0,3)];
                                const register double B_2_0 = B_p[INDEX2(2,0,3)];
                                const register double B_0_1 = B_p[INDEX2(0,1,3)];
                                const register double B_1_1 = B_p[INDEX2(1,1,3)];
                                const register double B_2_1 = B_p[INDEX2(2,1,3)];
                                const register double B_0_2 = B_p[INDEX2(0,2,3)];
                                const register double B_1_2 = B_p[INDEX2(1,2,3)];
                                const register double B_2_2 = B_p[INDEX2(2,2,3)];
                                const register double B_0_3 = B_p[INDEX2(0,3,3)];
                                const register double B_1_3 = B_p[INDEX2(1,3,3)];
                                const register double B_2_3 = B_p[INDEX2(2,3,3)];
                                const register double B_0_4 = B_p[INDEX2(0,4,3)];
                                const register double B_1_4 = B_p[INDEX2(1,4,3)];
                                const register double B_2_4 = B_p[INDEX2(2,4,3)];
                                const register double B_0_5 = B_p[INDEX2(0,5,3)];
                                const register double B_1_5 = B_p[INDEX2(1,5,3)];
                                const register double B_2_5 = B_p[INDEX2(2,5,3)];
                                const register double B_0_6 = B_p[INDEX2(0,6,3)];
                                const register double B_1_6 = B_p[INDEX2(1,6,3)];
                                const register double B_2_6 = B_p[INDEX2(2,6,3)];
                                const register double B_0_7 = B_p[INDEX2(0,7,3)];
                                const register double B_1_7 = B_p[INDEX2(1,7,3)];
                                const register double B_2_7 = B_p[INDEX2(2,7,3)];
                                const register double tmp24_0 = B_2_0 + B_2_2;
                                const register double tmp19_0 = B_2_0 + B_2_1 + B_2_2 + B_2_3;
                                const register double tmp20_0 = B_1_0 + B_1_5;
                                const register double tmp22_0 = B_2_1 + B_2_3;
                                const register double tmp29_0 = B_2_2 + B_2_3;
                                const register double tmp34_0 = B_1_0 + B_1_1 + B_1_4 + B_1_5;
                                const register double tmp15_0 = B_2_4 + B_2_5 + B_2_6 + B_2_7;
                                const register double tmp7_0 = B_0_0 + B_0_4;
                                const register double tmp40_0 = B_1_1 + B_1_4;
                                const register double tmp14_0 = B_1_4 + B_1_5;
                                const register double tmp35_0 = B_1_2 + B_1_3 + B_1_6 + B_1_7;
                                const register double tmp17_0 = B_1_6 + B_1_7;
                                const register double tmp8_0 = B_1_2 + B_1_6;
                                const register double tmp28_0 = B_2_4 + B_2_5;
                                const register double tmp10_0 = B_0_1 + B_0_3;
                                const register double tmp9_0 = B_2_5 + B_2_6;
                                const register double tmp30_0 = B_2_0 + B_2_1;
                                const register double tmp27_0 = B_0_0 + B_0_6;
                                const register double tmp32_0 = B_0_0 + B_0_2 + B_0_4 + B_0_6;
                                const register double tmp16_0 = B_0_0 + B_0_2;
                                const register double tmp2_0 = B_1_3 + B_1_7;
                                const register double tmp3_0 = B_2_1 + B_2_2;
                                const register double tmp33_0 = B_0_1 + B_0_3 + B_0_5 + B_0_7;
                                const register double tmp23_0 = B_2_5 + B_2_7;
                                const register double tmp36_0 = B_2_4 + B_2_7;
                                const register double tmp39_0 = B_0_3 + B_0_5;
                                const register double tmp41_0 = B_1_3 + B_1_6;
                                const register double tmp5_0 = B_1_0 + B_1_4;
                                const register double tmp18_0 = B_1_0 + B_1_1;
                                const register double tmp0_0 = B_0_1 + B_0_5;
                                const register double tmp37_0 = B_2_0 + B_2_3;
                                const register double tmp25_0 = B_2_4 + B_2_6;
                                const register double tmp38_0 = B_0_2 + B_0_4;
                                const register double tmp31_0 = B_2_6 + B_2_7;
                                const register double tmp1_0 = B_0_2 + B_0_6;
                                const register double tmp11_0 = B_0_4 + B_0_6;
                                const register double tmp21_0 = B_1_2 + B_1_7;
                                const register double tmp6_0 = B_1_1 + B_1_5;
                                const register double tmp13_0 = B_1_2 + B_1_3;
                                const register double tmp12_0 = B_0_5 + B_0_7;
                                const register double tmp26_0 = B_0_1 + B_0_7;
                                const register double tmp4_0 = B_0_3 + B_0_7;
                                const register double tmp324_1 = B_0_6*w104;
                                const register double tmp209_1 = B_1_4*w114;
                                const register double tmp255_1 = B_2_3*w103;
                                const register double tmp37_1 = B_1_1*w111;
                                const register double tmp326_1 = tmp38_0*w108;
                                const register double tmp179_1 = tmp21_0*w94;
                                const register double tmp102_1 = tmp4_0*w93;
                                const register double tmp251_1 = B_0_2*w125;
                                const register double tmp321_1 = B_0_4*w90;
                                const register double tmp198_1 = tmp41_0*w97;
                                const register double tmp15_1 = tmp11_0*w108;
                                const register double tmp158_1 = tmp17_0*w99;
                                const register double tmp138_1 = tmp35_0*w94;
                                const register double tmp5_1 = B_2_3*w100;
                                const register double tmp51_1 = tmp25_0*w103;
                                const register double tmp258_1 = B_2_7*w100;
                                const register double tmp221_1 = B_0_4*w125;
                                const register double tmp70_1 = B_0_5*w98;
                                const register double tmp202_1 = tmp38_0*w96;
                                const register double tmp247_1 = tmp41_0*w94;
                                const register double tmp122_1 = tmp34_0*w94;
                                const register double tmp349_1 = B_2_1*w101;
                                const register double tmp167_1 = tmp31_0*w103;
                                const register double tmp408_1 = tmp8_0*w111;
                                const register double tmp20_1 = tmp16_0*w104;
                                const register double tmp2_1 = B_2_7*w103;
                                const register double tmp398_1 = B_1_7*w99;
                                const register double tmp262_1 = B_1_7*w119;
                                const register double tmp97_1 = tmp5_0*w94;
                                const register double tmp157_1 = tmp30_0*w92;
                                const register double tmp67_1 = tmp18_0*w107;
                                const register double tmp144_1 = tmp7_0*w110;
                                const register double tmp264_1 = B_0_0*w120;
                                const register double tmp396_1 = tmp6_0*w97;
                                const register double tmp218_1 = B_1_6*w111;
                                const register double tmp147_1 = tmp8_0*w109;
                                const register double tmp39_1 = tmp12_0*w108;
                                const register double tmp58_1 = tmp15_0*w112;
                                const register double tmp214_1 = B_1_3*w115;
                                const register double tmp69_1 = tmp30_0*w95;
                                const register double tmp54_1 = tmp18_0*w99;
                                const register double tmp261_1 = B_2_0*w101;
                                const register double tmp390_1 = tmp24_0*w103;
                                const register double tmp374_1 = B_2_1*w103;
                                const register double tmp151_1 = tmp2_0*w105;
                                const register double tmp274_1 = tmp22_0*w95;
                                const register double tmp126_1 = tmp16_0*w108;
                                const register double tmp302_1 = B_2_6*w122;
                                const register double tmp14_1 = tmp10_0*w106;
                                const register double tmp190_1 = B_1_2*w99;
                                const register double tmp9_1 = B_2_4*w101;
                                const register double tmp150_1 = B_2_7*w101;
                                const register double tmp66_1 = tmp29_0*w92;
                                const register double tmp334_1 = B_0_7*w124;
                                const register double tmp252_1 = tmp38_0*w93;
                                const register double tmp392_1 = tmp5_0*w99;
                                const register double tmp4_1 = tmp2_0*w99;
                                const register double tmp63_1 = B_0_4*w121;
                                const register double tmp152_1 = tmp5_0*w111;
                                const register double tmp133_1 = tmp4_0*w96;
                                const register double tmp195_1 = tmp40_0*w94;
                                const register double tmp145_1 = B_2_3*w92;
                                const register double tmp35_1 = tmp25_0*w116;
                                const register double tmp316_1 = B_0_3*w98;
                                const register double tmp208_1 = B_2_0*w103;
                                const register double tmp373_1 = B_1_2*w115;
                                const register double tmp62_1 = tmp26_0*w93;
                                const register double tmp85_1 = tmp10_0*w90;
                                const register double tmp187_1 = tmp11_0*w106;
                                const register double tmp342_1 = B_2_5*w92;
                                const register double tmp77_1 = tmp33_0*w108;
                                const register double tmp121_1 = B_2_4*w116;
                                const register double tmp379_1 = B_2_2*w101;
                                const register double tmp391_1 = tmp23_0*w92;
                                const register double tmp451_1 = tmp22_0*w116;
                                const register double tmp73_1 = B_0_2*w90;
                                const register double tmp436_1 = B_1_6*w114;
                                const register double tmp233_1 = tmp29_0*w117;
                                const register double tmp426_1 = tmp23_0*w113;
                                const register double tmp124_1 = tmp11_0*w104;
                                const register double tmp148_1 = B_2_0*w100;
                                const register double tmp314_1 = tmp29_0*w113;
                                const register double tmp163_1 = B_0_0*w124;
                                const register double tmp301_1 = B_1_6*w115;
                                const register double tmp171_1 = B_2_4*w122;
                                const register double tmp7_1 = tmp4_0*w98;
                                const register double tmp11_1 = tmp7_0*w90;
                                const register double tmp308_1 = B_2_2*w116;
                                const register double tmp117_1 = B_2_3*w113;
                                const register double tmp403_1 = B_0_5*w104;
                                const register double tmp178_1 = tmp9_0*w112;
                                const register double tmp128_1 = tmp17_0*w107;
                                const register double tmp89_1 = tmp1_0*w110;
                                const register double tmp434_1 = tmp19_0*w95;
                                const register double tmp270_1 = B_1_7*w114;
                                const register double tmp358_1 = tmp6_0*w99;
                                const register double tmp411_1 = tmp14_0*w107;
                                const register double tmp433_1 = B_0_1*w125;
                                const register double tmp376_1 = B_2_6*w92;
                                const register double tmp42_1 = B_1_6*w99;
                                const register double tmp72_1 = tmp13_0*w105;
                                const register double tmp287_1 = tmp2_0*w111;
                                const register double tmp191_1 = tmp24_0*w113;
                                const register double tmp6_1 = tmp3_0*w95;
                                const register double tmp309_1 = B_1_4*w105;
                                const register double tmp367_1 = B_0_5*w125;
                                const register double tmp234_1 = tmp28_0*w112;
                                const register double tmp339_1 = B_1_2*w111;
                                const register double tmp335_1 = tmp41_0*w107;
                                const register double tmp173_1 = B_2_7*w113;
                                const register double tmp384_1 = tmp28_0*w113;
                                const register double tmp81_1 = tmp19_0*w112;
                                const register double tmp328_1 = B_2_7*w122;
                                const register double tmp413_1 = tmp30_0*w113;
                                const register double tmp111_1 = B_2_6*w101;
                                const register double tmp378_1 = B_0_4*w98;
                                const register double tmp211_1 = tmp21_0*w107;
                                const register double tmp79_1 = tmp35_0*w109;
                                const register double tmp354_1 = tmp1_0*w93;
                                const register double tmp333_1 = B_2_0*w123;
                                const register double tmp169_1 = B_0_3*w121;
                                const register double tmp300_1 = B_0_1*w121;
                                const register double tmp118_1 = B_2_7*w123;
                                const register double tmp372_1 = B_0_5*w121;
                                const register double tmp104_1 = B_2_5*w103;
                                const register double tmp304_1 = B_2_5*w113;
                                const register double tmp397_1 = tmp22_0*w102;
                                const register double tmp61_1 = tmp14_0*w97;
                                const register double tmp290_1 = tmp4_0*w106;
                                const register double tmp269_1 = tmp23_0*w103;
                                const register double tmp361_1 = tmp5_0*w97;
                                const register double tmp189_1 = tmp16_0*w110;
                                const register double tmp116_1 = B_2_0*w122;
                                const register double tmp449_1 = tmp24_0*w117;
                                const register double tmp268_1 = tmp11_0*w96;
                                const register double tmp325_1 = tmp39_0*w106;
                                const register double tmp27_1 = tmp20_0*w107;
                                const register double tmp174_1 = B_2_3*w123;
                                const register double tmp422_1 = tmp14_0*w99;
                                const register double tmp146_1 = tmp6_0*w107;
                                const register double tmp277_1 = B_1_2*w105;
                                const register double tmp443_1 = B_1_6*w118;
                                const register double tmp446_1 = B_1_7*w105;
                                const register double tmp91_1 = tmp8_0*w99;
                                const register double tmp405_1 = B_0_3*w125;
                                const register double tmp57_1 = tmp17_0*w91;
                                const register double tmp455_1 = B_2_6*w103;
                                const register double tmp110_1 = tmp0_0*w98;
                                const register double tmp341_1 = B_1_6*w119;
                                const register double tmp203_1 = B_0_7*w98;
                                const register double tmp371_1 = B_2_7*w116;
                                const register double tmp382_1 = B_0_3*w90;
                                const register double tmp320_1 = tmp28_0*w116;
                                const register double tmp107_1 = tmp2_0*w109;
                                const register double tmp149_1 = tmp4_0*w104;
                                const register double tmp159_1 = tmp29_0*w95;
                                const register double tmp253_1 = B_0_7*w121;
                                const register double tmp90_1 = B_2_1*w122;
                                const register double tmp360_1 = tmp2_0*w94;
                                const register double tmp414_1 = tmp28_0*w117;
                                const register double tmp26_1 = tmp16_0*w96;
                                const register double tmp84_1 = tmp11_0*w98;
                                const register double tmp225_1 = tmp1_0*w108;
                                const register double tmp1_1 = tmp1_0*w96;
                                const register double tmp344_1 = B_2_6*w100;
                                const register double tmp400_1 = B_1_5*w119;
                                const register double tmp345_1 = tmp36_0*w95;
                                const register double tmp407_1 = tmp5_0*w109;
                                const register double tmp228_1 = B_2_2*w122;
                                const register double tmp113_1 = tmp8_0*w105;
                                const register double tmp347_1 = B_0_0*w104;
                                const register double tmp286_1 = tmp23_0*w95;
                                const register double tmp395_1 = tmp8_0*w94;
                                const register double tmp401_1 = B_1_2*w118;
                                const register double tmp380_1 = B_1_7*w111;
                                const register double tmp327_1 = B_0_1*w110;
                                const register double tmp55_1 = tmp19_0*w117;
                                const register double tmp291_1 = tmp7_0*w108;
                                const register double tmp31_1 = tmp10_0*w98;
                                const register double tmp125_1 = tmp12_0*w106;
                                const register double tmp330_1 = tmp40_0*w109;
                                const register double tmp323_1 = tmp17_0*w97;
                                const register double tmp136_1 = tmp31_0*w95;
                                const register double tmp16_1 = tmp12_0*w110;
                                const register double tmp418_1 = B_0_6*w90;
                                const register double tmp428_1 = tmp25_0*w112;
                                const register double tmp385_1 = tmp30_0*w117;
                                const register double tmp351_1 = B_1_1*w118;
                                const register double tmp165_1 = B_0_7*w125;
                                const register double tmp298_1 = tmp29_0*w102;
                                const register double tmp295_1 = tmp34_0*w109;
                                const register double tmp296_1 = tmp28_0*w95;
                                const register double tmp283_1 = tmp25_0*w92;
                                const register double tmp230_1 = B_2_5*w123;
                                const register double tmp350_1 = B_0_1*w124;
                                const register double tmp293_1 = tmp31_0*w92;
                                const register double tmp8_1 = tmp5_0*w91;
                                const register double tmp215_1 = tmp9_0*w95;
                                const register double tmp329_1 = B_1_0*w114;
                                const register double tmp115_1 = tmp6_0*w111;
                                const register double tmp387_1 = tmp29_0*w116;
                                const register double tmp442_1 = B_1_1*w119;
                                const register double tmp281_1 = tmp33_0*w96;
                                const register double tmp415_1 = B_0_1*w98;
                                const register double tmp311_1 = B_1_3*w111;
                                const register double tmp50_1 = tmp23_0*w102;
                                const register double tmp294_1 = tmp35_0*w107;
                                const register double tmp362_1 = tmp27_0*w106;
                                const register double tmp340_1 = B_2_2*w103;
                                const register double tmp282_1 = tmp22_0*w103;
                                const register double tmp254_1 = tmp39_0*w96;
                                const register double tmp186_1 = tmp12_0*w104;
                                const register double tmp106_1 = tmp5_0*w107;
                                const register double tmp170_1 = tmp26_0*w96;
                                const register double tmp427_1 = tmp22_0*w117;
                                const register double tmp452_1 = B_2_1*w92;
                                const register double tmp437_1 = B_1_1*w115;
                                const register double tmp431_1 = B_0_0*w110;
                                const register double tmp201_1 = B_0_6*w121;
                                const register double tmp310_1 = B_0_6*w120;
                                const register double tmp331_1 = B_2_4*w113;
                                const register double tmp212_1 = tmp20_0*w109;
                                const register double tmp47_1 = B_1_4*w119;
                                const register double tmp420_1 = B_0_7*w120;
                                const register double tmp38_1 = tmp16_0*w106;
                                const register double tmp184_1 = tmp20_0*w97;
                                const register double tmp280_1 = tmp32_0*w93;
                                const register double tmp123_1 = tmp35_0*w97;
                                const register double tmp359_1 = tmp8_0*w91;
                                const register double tmp95_1 = tmp6_0*w91;
                                const register double tmp246_1 = tmp36_0*w112;
                                const register double tmp241_1 = B_2_6*w113;
                                const register double tmp3_1 = B_2_0*w92;
                                const register double tmp162_1 = tmp14_0*w94;
                                const register double tmp406_1 = tmp2_0*w107;
                                const register double tmp175_1 = tmp3_0*w117;
                                const register double tmp205_1 = B_0_1*w120;
                                const register double tmp416_1 = tmp29_0*w112;
                                const register double tmp48_1 = tmp21_0*w97;
                                const register double tmp101_1 = tmp32_0*w96;
                                const register double tmp231_1 = B_2_6*w116;
                                const register double tmp366_1 = B_0_2*w124;
                                const register double tmp34_1 = tmp11_0*w90;
                                const register double tmp200_1 = tmp39_0*w93;
                                const register double tmp44_1 = tmp10_0*w104;
                                const register double tmp346_1 = B_0_6*w125;
                                const register double tmp65_1 = tmp28_0*w103;
                                const register double tmp248_1 = B_1_2*w119;
                                const register double tmp288_1 = tmp5_0*w105;
                                const register double tmp227_1 = tmp4_0*w110;
                                const register double tmp343_1 = B_1_4*w99;
                                const register double tmp139_1 = tmp0_0*w90;
                                const register double tmp41_1 = tmp22_0*w92;
                                const register double tmp305_1 = B_2_1*w123;
                                const register double tmp307_1 = B_0_7*w90;
                                const register double tmp217_1 = B_2_3*w101;
                                const register double tmp219_1 = B_0_3*w124;
                                const register double tmp105_1 = B_2_2*w92;
                                const register double tmp423_1 = tmp13_0*w91;
                                const register double tmp381_1 = B_1_0*w105;
                                const register double tmp166_1 = tmp28_0*w102;
                                const register double tmp259_1 = B_0_6*w98;
                                const register double tmp402_1 = B_0_2*w110;
                                const register double tmp263_1 = B_0_1*w90;
                                const register double tmp276_1 = B_1_5*w111;
                                const register double tmp243_1 = B_2_2*w123;
                                const register double tmp131_1 = tmp13_0*w111;
                                const register double tmp265_1 = B_1_0*w118;
                                const register double tmp337_1 = B_1_5*w105;
                                const register double tmp172_1 = B_1_1*w99;
                                const register double tmp161_1 = tmp18_0*w91;
                                const register double tmp88_1 = tmp4_0*w108;
                                const register double tmp236_1 = B_1_5*w118;
                                const register double tmp273_1 = tmp41_0*w109;
                                const register double tmp154_1 = tmp38_0*w106;
                                const register double tmp440_1 = B_1_3*w99;
                                const register double tmp33_1 = B_1_3*w114;
                                const register double tmp141_1 = tmp30_0*w102;
                                const register double tmp393_1 = tmp25_0*w95;
                                const register double tmp289_1 = tmp24_0*w102;
                                const register double tmp180_1 = B_1_3*w119;
                                const register double tmp86_1 = tmp36_0*w117;
                                const register double tmp444_1 = B_1_5*w115;
                                const register double tmp447_1 = B_1_0*w111;
                                const register double tmp315_1 = tmp31_0*w117;
                                const register double tmp120_1 = tmp3_0*w112;
                                const register double tmp355_1 = tmp0_0*w96;
                                const register double tmp223_1 = B_0_5*w110;
                                const register double tmp60_1 = tmp12_0*w90;
                                const register double tmp239_1 = B_2_5*w122;
                                const register double tmp194_1 = tmp22_0*w112;
                                const register double tmp28_1 = tmp21_0*w109;
                                const register double tmp24_1 = tmp12_0*w93;
                                const register double tmp424_1 = tmp17_0*w94;
                                const register double tmp237_1 = B_0_4*w104;
                                const register double tmp322_1 = B_0_5*w120;
                                const register double tmp453_1 = B_2_2*w100;
                                const register double tmp238_1 = B_0_3*w110;
                                const register double tmp370_1 = B_2_4*w123;
                                const register double tmp155_1 = tmp39_0*w108;
                                const register double tmp23_1 = tmp19_0*w102;
                                const register double tmp181_1 = B_2_0*w116;
                                const register double tmp164_1 = tmp13_0*w97;
                                const register double tmp168_1 = tmp27_0*w93;
                                const register double tmp285_1 = tmp6_0*w109;
                                const register double tmp32_1 = tmp24_0*w112;
                                const register double tmp210_1 = B_2_7*w92;
                                const register double tmp399_1 = B_1_0*w91;
                                const register double tmp183_1 = B_0_4*w120;
                                const register double tmp185_1 = B_1_4*w118;
                                const register double tmp40_1 = tmp11_0*w110;
                                const register double tmp59_1 = tmp13_0*w94;
                                const register double tmp279_1 = tmp25_0*w102;
                                const register double tmp242_1 = B_1_7*w91;
                                const register double tmp454_1 = B_2_5*w101;
                                const register double tmp114_1 = tmp1_0*w90;
                                const register double tmp319_1 = tmp18_0*w94;
                                const register double tmp394_1 = tmp2_0*w91;
                                const register double tmp229_1 = B_2_1*w113;
                                const register double tmp93_1 = B_2_6*w123;
                                const register double tmp278_1 = tmp16_0*w90;
                                const register double tmp46_1 = tmp20_0*w94;
                                const register double tmp87_1 = tmp7_0*w106;
                                const register double tmp421_1 = tmp18_0*w111;
                                const register double tmp313_1 = tmp13_0*w99;
                                const register double tmp364_1 = B_0_4*w110;
                                const register double tmp412_1 = tmp13_0*w109;
                                const register double tmp30_1 = tmp23_0*w117;
                                const register double tmp365_1 = B_0_3*w104;
                                const register double tmp425_1 = tmp18_0*w97;
                                const register double tmp134_1 = tmp29_0*w103;
                                const register double tmp143_1 = tmp0_0*w108;
                                const register double tmp213_1 = B_2_4*w100;
                                const register double tmp435_1 = tmp15_0*w102;
                                const register double tmp78_1 = tmp34_0*w107;
                                const register double tmp76_1 = tmp32_0*w106;
                                const register double tmp132_1 = tmp7_0*w93;
                                const register double tmp348_1 = B_1_3*w91;
                                const register double tmp68_1 = tmp17_0*w109;
                                const register double tmp130_1 = tmp14_0*w105;
                                const register double tmp142_1 = tmp1_0*w106;
                                const register double tmp153_1 = B_2_4*w103;
                                const register double tmp432_1 = B_0_6*w124;
                                const register double tmp256_1 = B_2_4*w92;
                                const register double tmp204_1 = B_0_0*w90;
                                const register double tmp82_1 = tmp16_0*w93;
                                const register double tmp222_1 = tmp3_0*w102;
                                const register double tmp272_1 = tmp40_0*w107;
                                const register double tmp19_1 = tmp15_0*w95;
                                const register double tmp226_1 = tmp7_0*w104;
                                const register double tmp25_1 = B_1_4*w115;
                                const register double tmp199_1 = B_1_7*w118;
                                const register double tmp303_1 = B_1_1*w114;
                                const register double tmp318_1 = tmp30_0*w112;
                                const register double tmp100_1 = tmp33_0*w93;
                                const register double tmp439_1 = B_1_4*w111;
                                const register double tmp17_1 = tmp13_0*w107;
                                const register double tmp83_1 = tmp12_0*w96;
                                const register double tmp45_1 = B_1_1*w91;
                                const register double tmp75_1 = tmp14_0*w111;
                                const register double tmp368_1 = B_2_3*w122;
                                const register double tmp216_1 = B_0_2*w104;
                                const register double tmp429_1 = tmp24_0*w116;
                                const register double tmp388_1 = tmp33_0*w106;
                                const register double tmp13_1 = tmp9_0*w102;
                                const register double tmp18_1 = tmp14_0*w109;
                                const register double tmp112_1 = tmp36_0*w102;
                                const register double tmp140_1 = tmp34_0*w97;
                                const register double tmp383_1 = B_0_2*w120;
                                const register double tmp119_1 = tmp9_0*w117;
                                const register double tmp36_1 = B_1_6*w105;
                                const register double tmp266_1 = tmp10_0*w93;
                                const register double tmp336_1 = B_2_3*w116;
                                const register double tmp419_1 = tmp31_0*w116;
                                const register double tmp363_1 = tmp26_0*w108;
                                const register double tmp196_1 = B_1_0*w119;
                                const register double tmp182_1 = B_0_5*w90;
                                const register double tmp0_1 = tmp0_0*w93;
                                const register double tmp297_1 = tmp1_0*w104;
                                const register double tmp12_1 = tmp8_0*w97;
                                const register double tmp244_1 = tmp37_0*w117;
                                const register double tmp332_1 = B_1_7*w115;
                                const register double tmp156_1 = B_0_6*w110;
                                const register double tmp127_1 = tmp10_0*w110;
                                const register double tmp257_1 = B_1_5*w99;
                                const register double tmp64_1 = tmp27_0*w96;
                                const register double tmp389_1 = tmp32_0*w108;
                                const register double tmp10_1 = tmp6_0*w94;
                                const register double tmp317_1 = tmp14_0*w91;
                                const register double tmp275_1 = tmp12_0*w98;
                                const register double tmp99_1 = tmp2_0*w97;
                                const register double tmp192_1 = tmp25_0*w117;
                                const register double tmp224_1 = tmp0_0*w106;
                                const register double tmp410_1 = B_0_0*w121;
                                const register double tmp409_1 = tmp6_0*w105;
                                const register double tmp441_1 = B_1_4*w91;
                                const register double tmp206_1 = tmp26_0*w106;
                                const register double tmp232_1 = tmp31_0*w113;
                                const register double tmp108_1 = B_2_1*w100;
                                const register double tmp21_1 = tmp17_0*w111;
                                const register double tmp52_1 = tmp11_0*w93;
                                const register double tmp103_1 = tmp7_0*w96;
                                const register double tmp197_1 = tmp23_0*w116;
                                const register double tmp438_1 = B_1_3*w105;
                                const register double tmp235_1 = tmp30_0*w116;
                                const register double tmp29_1 = tmp22_0*w113;
                                const register double tmp98_1 = B_2_5*w116;
                                const register double tmp369_1 = B_2_0*w113;
                                const register double tmp74_1 = B_0_3*w120;
                                const register double tmp137_1 = tmp1_0*w98;
                                const register double tmp284_1 = tmp8_0*w107;
                                const register double tmp109_1 = tmp37_0*w95;
                                const register double tmp448_1 = tmp25_0*w113;
                                const register double tmp56_1 = tmp16_0*w98;
                                const register double tmp375_1 = B_1_5*w114;
                                const register double tmp188_1 = tmp10_0*w108;
                                const register double tmp94_1 = tmp0_0*w104;
                                const register double tmp260_1 = B_1_2*w91;
                                const register double tmp49_1 = B_1_3*w118;
                                const register double tmp404_1 = B_0_4*w124;
                                const register double tmp430_1 = B_0_7*w104;
                                const register double tmp299_1 = tmp0_0*w110;
                                const register double tmp417_1 = tmp17_0*w105;
                                const register double tmp445_1 = B_1_2*w114;
                                const register double tmp353_1 = B_0_7*w110;
                                const register double tmp135_1 = tmp28_0*w92;
                                const register double tmp357_1 = tmp4_0*w90;
                                const register double tmp356_1 = tmp7_0*w98;
                                const register double tmp267_1 = B_1_0*w115;
                                const register double tmp220_1 = B_1_1*w105;
                                const register double tmp338_1 = B_0_0*w125;
                                const register double tmp292_1 = tmp30_0*w103;
                                const register double tmp53_1 = tmp10_0*w96;
                                const register double tmp43_1 = tmp24_0*w95;
                                const register double tmp386_1 = tmp31_0*w112;
                                const register double tmp92_1 = B_2_2*w113;
                                const register double tmp245_1 = B_0_5*w124;
                                const register double tmp306_1 = B_0_0*w98;
                                const register double tmp96_1 = tmp37_0*w112;
                                const register double tmp177_1 = B_1_6*w91;
                                const register double tmp176_1 = B_0_2*w98;
                                const register double tmp71_1 = tmp31_0*w102;
                                const register double tmp129_1 = tmp18_0*w109;
                                const register double tmp352_1 = tmp37_0*w102;
                                const register double tmp312_1 = B_0_2*w121;
                                const register double tmp80_1 = tmp15_0*w117;
                                const register double tmp450_1 = tmp23_0*w112;
                                const register double tmp271_1 = tmp24_0*w92;
                                const register double tmp22_1 = tmp18_0*w105;
                                const register double tmp250_1 = tmp40_0*w97;
                                const register double tmp249_1 = B_2_1*w116;
                                const register double tmp207_1 = tmp27_0*w108;
                                const register double tmp193_1 = B_1_5*w91;
                                const register double tmp240_1 = B_1_0*w99;
                                const register double tmp160_1 = B_0_1*w104;
                                const register double tmp377_1 = B_2_5*w100;
                                EM_S[INDEX2(7,3,8)]+=tmp0_1 + tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                EM_S[INDEX2(4,7,8)]+=tmp14_1 + tmp15_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1;
                                EM_S[INDEX2(1,3,8)]+=tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp33_1 + tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1;
                                EM_S[INDEX2(6,4,8)]+=tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                                EM_S[INDEX2(3,0,8)]+=tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1 + tmp60_1 + tmp61_1;
                                EM_S[INDEX2(5,4,8)]+=tmp62_1 + tmp63_1 + tmp64_1 + tmp65_1 + tmp66_1 + tmp67_1 + tmp68_1 + tmp69_1 + tmp70_1 + tmp71_1 + tmp72_1 + tmp73_1 + tmp74_1 + tmp75_1;
                                EM_S[INDEX2(0,7,8)]+=tmp76_1 + tmp77_1 + tmp78_1 + tmp79_1 + tmp80_1 + tmp81_1;
                                EM_S[INDEX2(5,6,8)]+=tmp17_1 + tmp18_1 + tmp19_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp82_1 + tmp83_1 + tmp84_1 + tmp85_1;
                                EM_S[INDEX2(2,6,8)]+=tmp86_1 + tmp87_1 + tmp88_1 + tmp89_1 + tmp90_1 + tmp91_1 + tmp92_1 + tmp93_1 + tmp94_1 + tmp95_1 + tmp96_1 + tmp97_1 + tmp98_1 + tmp99_1;
                                EM_S[INDEX2(1,6,8)]+=tmp100_1 + tmp101_1 + tmp78_1 + tmp79_1 + tmp80_1 + tmp81_1;
                                EM_S[INDEX2(5,1,8)]+=tmp102_1 + tmp103_1 + tmp104_1 + tmp105_1 + tmp106_1 + tmp107_1 + tmp108_1 + tmp109_1 + tmp110_1 + tmp111_1 + tmp112_1 + tmp113_1 + tmp114_1 + tmp115_1;
                                EM_S[INDEX2(3,7,8)]+=tmp0_1 + tmp10_1 + tmp116_1 + tmp117_1 + tmp118_1 + tmp119_1 + tmp11_1 + tmp120_1 + tmp121_1 + tmp12_1 + tmp1_1 + tmp4_1 + tmp7_1 + tmp8_1;
                                EM_S[INDEX2(2,5,8)]+=tmp122_1 + tmp123_1 + tmp76_1 + tmp77_1 + tmp80_1 + tmp81_1;
                                EM_S[INDEX2(0,3,8)]+=tmp124_1 + tmp125_1 + tmp126_1 + tmp127_1 + tmp128_1 + tmp129_1 + tmp130_1 + tmp131_1 + tmp55_1 + tmp58_1;
                                EM_S[INDEX2(7,2,8)]+=tmp132_1 + tmp133_1 + tmp134_1 + tmp135_1 + tmp136_1 + tmp137_1 + tmp138_1 + tmp139_1 + tmp140_1 + tmp141_1;
                                EM_S[INDEX2(4,0,8)]+=tmp13_1 + tmp142_1 + tmp143_1 + tmp144_1 + tmp145_1 + tmp146_1 + tmp147_1 + tmp148_1 + tmp149_1 + tmp150_1 + tmp151_1 + tmp152_1 + tmp153_1 + tmp6_1;
                                EM_S[INDEX2(1,2,8)]+=tmp128_1 + tmp129_1 + tmp130_1 + tmp131_1 + tmp52_1 + tmp53_1 + tmp55_1 + tmp56_1 + tmp58_1 + tmp60_1;
                                EM_S[INDEX2(6,7,8)]+=tmp154_1 + tmp155_1 + tmp156_1 + tmp157_1 + tmp158_1 + tmp159_1 + tmp160_1 + tmp161_1 + tmp162_1 + tmp163_1 + tmp164_1 + tmp165_1 + tmp166_1 + tmp167_1;
                                EM_S[INDEX2(3,3,8)]+=tmp168_1 + tmp169_1 + tmp170_1 + tmp171_1 + tmp172_1 + tmp173_1 + tmp174_1 + tmp175_1 + tmp176_1 + tmp177_1 + tmp178_1 + tmp179_1 + tmp180_1 + tmp181_1 + tmp182_1 + tmp183_1 + tmp184_1 + tmp185_1;
                                EM_S[INDEX2(2,0,8)]+=tmp186_1 + tmp187_1 + tmp188_1 + tmp189_1 + tmp190_1 + tmp191_1 + tmp192_1 + tmp193_1 + tmp194_1 + tmp195_1 + tmp196_1 + tmp197_1 + tmp198_1 + tmp199_1;
                                EM_S[INDEX2(7,6,8)]+=tmp157_1 + tmp158_1 + tmp159_1 + tmp161_1 + tmp162_1 + tmp164_1 + tmp166_1 + tmp167_1 + tmp200_1 + tmp201_1 + tmp202_1 + tmp203_1 + tmp204_1 + tmp205_1;
                                EM_S[INDEX2(4,4,8)]+=tmp206_1 + tmp207_1 + tmp208_1 + tmp209_1 + tmp210_1 + tmp211_1 + tmp212_1 + tmp213_1 + tmp214_1 + tmp215_1 + tmp216_1 + tmp217_1 + tmp218_1 + tmp219_1 + tmp220_1 + tmp221_1 + tmp222_1 + tmp223_1;
                                EM_S[INDEX2(6,3,8)]+=tmp134_1 + tmp135_1 + tmp136_1 + tmp138_1 + tmp140_1 + tmp141_1 + tmp224_1 + tmp225_1 + tmp226_1 + tmp227_1;
                                EM_S[INDEX2(1,5,8)]+=tmp102_1 + tmp103_1 + tmp106_1 + tmp107_1 + tmp110_1 + tmp113_1 + tmp114_1 + tmp115_1 + tmp228_1 + tmp229_1 + tmp230_1 + tmp231_1 + tmp86_1 + tmp96_1;
                                EM_S[INDEX2(3,6,8)]+=tmp132_1 + tmp133_1 + tmp137_1 + tmp138_1 + tmp139_1 + tmp140_1 + tmp232_1 + tmp233_1 + tmp234_1 + tmp235_1;
                                EM_S[INDEX2(2,2,8)]+=tmp206_1 + tmp207_1 + tmp236_1 + tmp237_1 + tmp238_1 + tmp239_1 + tmp240_1 + tmp241_1 + tmp242_1 + tmp243_1 + tmp244_1 + tmp245_1 + tmp246_1 + tmp247_1 + tmp248_1 + tmp249_1 + tmp250_1 + tmp251_1;
                                EM_S[INDEX2(7,7,8)]+=tmp215_1 + tmp222_1 + tmp247_1 + tmp250_1 + tmp252_1 + tmp253_1 + tmp254_1 + tmp255_1 + tmp256_1 + tmp257_1 + tmp258_1 + tmp259_1 + tmp260_1 + tmp261_1 + tmp262_1 + tmp263_1 + tmp264_1 + tmp265_1;
                                EM_S[INDEX2(5,7,8)]+=tmp266_1 + tmp267_1 + tmp268_1 + tmp269_1 + tmp270_1 + tmp271_1 + tmp272_1 + tmp273_1 + tmp274_1 + tmp275_1 + tmp276_1 + tmp277_1 + tmp278_1 + tmp279_1;
                                EM_S[INDEX2(5,3,8)]+=tmp280_1 + tmp281_1 + tmp282_1 + tmp283_1 + tmp284_1 + tmp285_1 + tmp286_1 + tmp287_1 + tmp288_1 + tmp289_1;
                                EM_S[INDEX2(4,1,8)]+=tmp290_1 + tmp291_1 + tmp292_1 + tmp293_1 + tmp294_1 + tmp295_1 + tmp296_1 + tmp297_1 + tmp298_1 + tmp299_1;
                                EM_S[INDEX2(1,1,8)]+=tmp211_1 + tmp212_1 + tmp244_1 + tmp246_1 + tmp252_1 + tmp254_1 + tmp300_1 + tmp301_1 + tmp302_1 + tmp303_1 + tmp304_1 + tmp305_1 + tmp306_1 + tmp307_1 + tmp308_1 + tmp309_1 + tmp310_1 + tmp311_1;
                                EM_S[INDEX2(2,7,8)]+=tmp138_1 + tmp140_1 + tmp224_1 + tmp225_1 + tmp226_1 + tmp227_1 + tmp232_1 + tmp233_1 + tmp234_1 + tmp235_1;
                                EM_S[INDEX2(3,2,8)]+=tmp312_1 + tmp313_1 + tmp314_1 + tmp315_1 + tmp316_1 + tmp317_1 + tmp318_1 + tmp319_1 + tmp320_1 + tmp321_1 + tmp322_1 + tmp323_1 + tmp62_1 + tmp64_1;
                                EM_S[INDEX2(0,0,8)]+=tmp175_1 + tmp178_1 + tmp324_1 + tmp325_1 + tmp326_1 + tmp327_1 + tmp328_1 + tmp329_1 + tmp330_1 + tmp331_1 + tmp332_1 + tmp333_1 + tmp334_1 + tmp335_1 + tmp336_1 + tmp337_1 + tmp338_1 + tmp339_1;
                                EM_S[INDEX2(6,6,8)]+=tmp179_1 + tmp184_1 + tmp325_1 + tmp326_1 + tmp340_1 + tmp341_1 + tmp342_1 + tmp343_1 + tmp344_1 + tmp345_1 + tmp346_1 + tmp347_1 + tmp348_1 + tmp349_1 + tmp350_1 + tmp351_1 + tmp352_1 + tmp353_1;
                                EM_S[INDEX2(5,0,8)]+=tmp292_1 + tmp293_1 + tmp294_1 + tmp295_1 + tmp296_1 + tmp298_1 + tmp354_1 + tmp355_1 + tmp356_1 + tmp357_1;
                                EM_S[INDEX2(7,1,8)]+=tmp280_1 + tmp281_1 + tmp282_1 + tmp283_1 + tmp286_1 + tmp289_1 + tmp358_1 + tmp359_1 + tmp360_1 + tmp361_1;
                                EM_S[INDEX2(4,5,8)]+=tmp362_1 + tmp363_1 + tmp364_1 + tmp365_1 + tmp366_1 + tmp367_1 + tmp65_1 + tmp66_1 + tmp67_1 + tmp68_1 + tmp69_1 + tmp71_1 + tmp72_1 + tmp75_1;
                                EM_S[INDEX2(0,4,8)]+=tmp119_1 + tmp120_1 + tmp142_1 + tmp143_1 + tmp144_1 + tmp146_1 + tmp147_1 + tmp149_1 + tmp151_1 + tmp152_1 + tmp368_1 + tmp369_1 + tmp370_1 + tmp371_1;
                                EM_S[INDEX2(5,5,8)]+=tmp168_1 + tmp170_1 + tmp330_1 + tmp335_1 + tmp345_1 + tmp352_1 + tmp372_1 + tmp373_1 + tmp374_1 + tmp375_1 + tmp376_1 + tmp377_1 + tmp378_1 + tmp379_1 + tmp380_1 + tmp381_1 + tmp382_1 + tmp383_1;
                                EM_S[INDEX2(1,4,8)]+=tmp294_1 + tmp295_1 + tmp354_1 + tmp355_1 + tmp356_1 + tmp357_1 + tmp384_1 + tmp385_1 + tmp386_1 + tmp387_1;
                                EM_S[INDEX2(6,0,8)]+=tmp388_1 + tmp389_1 + tmp390_1 + tmp391_1 + tmp392_1 + tmp393_1 + tmp394_1 + tmp395_1 + tmp396_1 + tmp397_1;
                                EM_S[INDEX2(7,5,8)]+=tmp195_1 + tmp198_1 + tmp266_1 + tmp268_1 + tmp269_1 + tmp271_1 + tmp274_1 + tmp275_1 + tmp278_1 + tmp279_1 + tmp398_1 + tmp399_1 + tmp400_1 + tmp401_1;
                                EM_S[INDEX2(2,3,8)]+=tmp313_1 + tmp314_1 + tmp315_1 + tmp317_1 + tmp318_1 + tmp319_1 + tmp320_1 + tmp323_1 + tmp362_1 + tmp363_1 + tmp402_1 + tmp403_1 + tmp404_1 + tmp405_1;
                                EM_S[INDEX2(2,1,8)]+=tmp124_1 + tmp125_1 + tmp126_1 + tmp127_1 + tmp54_1 + tmp55_1 + tmp57_1 + tmp58_1 + tmp59_1 + tmp61_1;
                                EM_S[INDEX2(4,2,8)]+=tmp388_1 + tmp389_1 + tmp390_1 + tmp391_1 + tmp393_1 + tmp397_1 + tmp406_1 + tmp407_1 + tmp408_1 + tmp409_1;
                                EM_S[INDEX2(1,0,8)]+=tmp200_1 + tmp202_1 + tmp410_1 + tmp411_1 + tmp412_1 + tmp413_1 + tmp414_1 + tmp415_1 + tmp416_1 + tmp417_1 + tmp418_1 + tmp419_1 + tmp420_1 + tmp421_1;
                                EM_S[INDEX2(6,5,8)]+=tmp14_1 + tmp15_1 + tmp16_1 + tmp19_1 + tmp20_1 + tmp23_1 + tmp422_1 + tmp423_1 + tmp424_1 + tmp425_1;
                                EM_S[INDEX2(3,5,8)]+=tmp280_1 + tmp281_1 + tmp358_1 + tmp359_1 + tmp360_1 + tmp361_1 + tmp426_1 + tmp427_1 + tmp428_1 + tmp429_1;
                                EM_S[INDEX2(0,1,8)]+=tmp154_1 + tmp155_1 + tmp411_1 + tmp412_1 + tmp413_1 + tmp414_1 + tmp416_1 + tmp417_1 + tmp419_1 + tmp421_1 + tmp430_1 + tmp431_1 + tmp432_1 + tmp433_1;
                                EM_S[INDEX2(7,0,8)]+=tmp100_1 + tmp101_1 + tmp122_1 + tmp123_1 + tmp434_1 + tmp435_1;
                                EM_S[INDEX2(4,6,8)]+=tmp27_1 + tmp28_1 + tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp436_1 + tmp437_1 + tmp438_1 + tmp439_1 + tmp43_1 + tmp44_1 + tmp50_1 + tmp51_1;
                                EM_S[INDEX2(5,2,8)]+=tmp100_1 + tmp101_1 + tmp434_1 + tmp435_1 + tmp78_1 + tmp79_1;
                                EM_S[INDEX2(6,1,8)]+=tmp122_1 + tmp123_1 + tmp434_1 + tmp435_1 + tmp76_1 + tmp77_1;
                                EM_S[INDEX2(3,1,8)]+=tmp24_1 + tmp26_1 + tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp34_1 + tmp35_1 + tmp440_1 + tmp441_1 + tmp442_1 + tmp443_1 + tmp46_1 + tmp48_1;
                                EM_S[INDEX2(0,2,8)]+=tmp186_1 + tmp187_1 + tmp188_1 + tmp189_1 + tmp191_1 + tmp192_1 + tmp194_1 + tmp197_1 + tmp272_1 + tmp273_1 + tmp444_1 + tmp445_1 + tmp446_1 + tmp447_1;
                                EM_S[INDEX2(7,4,8)]+=tmp19_1 + tmp23_1 + tmp422_1 + tmp423_1 + tmp424_1 + tmp425_1 + tmp82_1 + tmp83_1 + tmp84_1 + tmp85_1;
                                EM_S[INDEX2(0,6,8)]+=tmp388_1 + tmp389_1 + tmp406_1 + tmp407_1 + tmp408_1 + tmp409_1 + tmp448_1 + tmp449_1 + tmp450_1 + tmp451_1;
                                EM_S[INDEX2(6,2,8)]+=tmp109_1 + tmp112_1 + tmp452_1 + tmp453_1 + tmp454_1 + tmp455_1 + tmp87_1 + tmp88_1 + tmp89_1 + tmp91_1 + tmp94_1 + tmp95_1 + tmp97_1 + tmp99_1;
                                EM_S[INDEX2(4,3,8)]+=tmp434_1 + tmp435_1 + tmp76_1 + tmp77_1 + tmp78_1 + tmp79_1;
                                EM_S[INDEX2(1,7,8)]+=tmp280_1 + tmp281_1 + tmp284_1 + tmp285_1 + tmp287_1 + tmp288_1 + tmp426_1 + tmp427_1 + tmp428_1 + tmp429_1;
                                EM_S[INDEX2(0,5,8)]+=tmp290_1 + tmp291_1 + tmp294_1 + tmp295_1 + tmp297_1 + tmp299_1 + tmp384_1 + tmp385_1 + tmp386_1 + tmp387_1;
                                EM_S[INDEX2(3,4,8)]+=tmp100_1 + tmp101_1 + tmp122_1 + tmp123_1 + tmp80_1 + tmp81_1;
                                EM_S[INDEX2(2,4,8)]+=tmp388_1 + tmp389_1 + tmp392_1 + tmp394_1 + tmp395_1 + tmp396_1 + tmp448_1 + tmp449_1 + tmp450_1 + tmp451_1;
                            } else { // constant data
                                const register double B_0 = B_p[0];
                                const register double B_1 = B_p[1];
                                const register double B_2 = B_p[2];
                                const register double tmp7_1 = B_2*w133;
                                const register double tmp17_1 = B_0*w143;
                                const register double tmp1_1 = B_1*w127;
                                const register double tmp8_1 = B_2*w135;
                                const register double tmp15_1 = B_0*w141;
                                const register double tmp3_1 = B_0*w129;
                                const register double tmp13_1 = B_1*w139;
                                const register double tmp2_1 = B_0*w126;
                                const register double tmp12_1 = B_0*w138;
                                const register double tmp5_1 = B_2*w131;
                                const register double tmp10_1 = B_2*w136;
                                const register double tmp4_1 = B_1*w130;
                                const register double tmp6_1 = B_1*w132;
                                const register double tmp9_1 = B_1*w134;
                                const register double tmp14_1 = B_2*w140;
                                const register double tmp0_1 = B_2*w128;
                                const register double tmp16_1 = B_1*w142;
                                const register double tmp11_1 = B_0*w137;
                                EM_S[INDEX2(7,3,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                EM_S[INDEX2(4,7,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                EM_S[INDEX2(1,3,8)]+=tmp2_1 + tmp6_1 + tmp7_1;
                                EM_S[INDEX2(6,4,8)]+=tmp3_1 + tmp8_1 + tmp9_1;
                                EM_S[INDEX2(3,0,8)]+=tmp10_1 + tmp1_1 + tmp2_1;
                                EM_S[INDEX2(5,4,8)]+=tmp11_1 + tmp4_1 + tmp8_1;
                                EM_S[INDEX2(0,7,8)]+=tmp10_1 + tmp12_1 + tmp13_1;
                                EM_S[INDEX2(5,6,8)]+=tmp2_1 + tmp4_1 + tmp5_1;
                                EM_S[INDEX2(2,6,8)]+=tmp14_1 + tmp1_1 + tmp3_1;
                                EM_S[INDEX2(1,6,8)]+=tmp10_1 + tmp13_1 + tmp15_1;
                                EM_S[INDEX2(5,1,8)]+=tmp0_1 + tmp2_1 + tmp4_1;
                                EM_S[INDEX2(3,7,8)]+=tmp14_1 + tmp1_1 + tmp2_1;
                                EM_S[INDEX2(2,5,8)]+=tmp10_1 + tmp12_1 + tmp16_1;
                                EM_S[INDEX2(0,3,8)]+=tmp10_1 + tmp3_1 + tmp4_1;
                                EM_S[INDEX2(7,2,8)]+=tmp16_1 + tmp2_1 + tmp8_1;
                                EM_S[INDEX2(4,0,8)]+=tmp0_1 + tmp3_1 + tmp4_1;
                                EM_S[INDEX2(1,2,8)]+=tmp10_1 + tmp2_1 + tmp4_1;
                                EM_S[INDEX2(6,7,8)]+=tmp17_1 + tmp1_1 + tmp8_1;
                                EM_S[INDEX2(3,3,8)]+=tmp11_1 + tmp14_1 + tmp9_1;
                                EM_S[INDEX2(2,0,8)]+=tmp3_1 + tmp7_1 + tmp9_1;
                                EM_S[INDEX2(7,6,8)]+=tmp11_1 + tmp1_1 + tmp8_1;
                                EM_S[INDEX2(4,4,8)]+=tmp0_1 + tmp17_1 + tmp6_1;
                                EM_S[INDEX2(6,3,8)]+=tmp16_1 + tmp3_1 + tmp8_1;
                                EM_S[INDEX2(1,5,8)]+=tmp14_1 + tmp2_1 + tmp4_1;
                                EM_S[INDEX2(3,6,8)]+=tmp16_1 + tmp2_1 + tmp7_1;
                                EM_S[INDEX2(2,2,8)]+=tmp14_1 + tmp17_1 + tmp9_1;
                                EM_S[INDEX2(7,7,8)]+=tmp0_1 + tmp11_1 + tmp9_1;
                                EM_S[INDEX2(5,7,8)]+=tmp2_1 + tmp6_1 + tmp8_1;
                                EM_S[INDEX2(5,3,8)]+=tmp15_1 + tmp4_1 + tmp8_1;
                                EM_S[INDEX2(4,1,8)]+=tmp13_1 + tmp3_1 + tmp8_1;
                                EM_S[INDEX2(1,1,8)]+=tmp11_1 + tmp14_1 + tmp6_1;
                                EM_S[INDEX2(2,7,8)]+=tmp16_1 + tmp3_1 + tmp7_1;
                                EM_S[INDEX2(3,2,8)]+=tmp11_1 + tmp1_1 + tmp7_1;
                                EM_S[INDEX2(0,0,8)]+=tmp14_1 + tmp17_1 + tmp6_1;
                                EM_S[INDEX2(6,6,8)]+=tmp0_1 + tmp17_1 + tmp9_1;
                                EM_S[INDEX2(5,0,8)]+=tmp13_1 + tmp2_1 + tmp8_1;
                                EM_S[INDEX2(7,1,8)]+=tmp15_1 + tmp1_1 + tmp8_1;
                                EM_S[INDEX2(4,5,8)]+=tmp17_1 + tmp4_1 + tmp8_1;
                                EM_S[INDEX2(0,4,8)]+=tmp14_1 + tmp3_1 + tmp4_1;
                                EM_S[INDEX2(5,5,8)]+=tmp0_1 + tmp11_1 + tmp6_1;
                                EM_S[INDEX2(1,4,8)]+=tmp13_1 + tmp2_1 + tmp7_1;
                                EM_S[INDEX2(6,0,8)]+=tmp12_1 + tmp1_1 + tmp8_1;
                                EM_S[INDEX2(7,5,8)]+=tmp2_1 + tmp8_1 + tmp9_1;
                                EM_S[INDEX2(2,3,8)]+=tmp17_1 + tmp1_1 + tmp7_1;
                                EM_S[INDEX2(2,1,8)]+=tmp10_1 + tmp1_1 + tmp3_1;
                                EM_S[INDEX2(4,2,8)]+=tmp12_1 + tmp4_1 + tmp8_1;
                                EM_S[INDEX2(1,0,8)]+=tmp11_1 + tmp4_1 + tmp7_1;
                                EM_S[INDEX2(6,5,8)]+=tmp1_1 + tmp3_1 + tmp5_1;
                                EM_S[INDEX2(3,5,8)]+=tmp15_1 + tmp1_1 + tmp7_1;
                                EM_S[INDEX2(0,1,8)]+=tmp17_1 + tmp4_1 + tmp7_1;
                                EM_S[INDEX2(7,0,8)]+=tmp15_1 + tmp16_1 + tmp5_1;
                                EM_S[INDEX2(4,6,8)]+=tmp3_1 + tmp6_1 + tmp8_1;
                                EM_S[INDEX2(5,2,8)]+=tmp13_1 + tmp15_1 + tmp5_1;
                                EM_S[INDEX2(6,1,8)]+=tmp12_1 + tmp16_1 + tmp5_1;
                                EM_S[INDEX2(3,1,8)]+=tmp2_1 + tmp7_1 + tmp9_1;
                                EM_S[INDEX2(0,2,8)]+=tmp3_1 + tmp6_1 + tmp7_1;
                                EM_S[INDEX2(7,4,8)]+=tmp1_1 + tmp2_1 + tmp5_1;
                                EM_S[INDEX2(0,6,8)]+=tmp12_1 + tmp4_1 + tmp7_1;
                                EM_S[INDEX2(6,2,8)]+=tmp0_1 + tmp1_1 + tmp3_1;
                                EM_S[INDEX2(4,3,8)]+=tmp12_1 + tmp13_1 + tmp5_1;
                                EM_S[INDEX2(1,7,8)]+=tmp15_1 + tmp4_1 + tmp7_1;
                                EM_S[INDEX2(0,5,8)]+=tmp13_1 + tmp3_1 + tmp7_1;
                                EM_S[INDEX2(3,4,8)]+=tmp10_1 + tmp15_1 + tmp16_1;
                                EM_S[INDEX2(2,4,8)]+=tmp12_1 + tmp1_1 + tmp7_1;
                            }
                        }
                        ///////////////
                        // process C //
                        ///////////////
                        if (!C.isEmpty()) {
                            add_EM_S=true;
                            const double* C_p=const_cast<escript::Data*>(&C)->getSampleDataRO(e);
                            if (C.actsExpanded()) {
                                const register double C_0_0 = C_p[INDEX2(0,0,3)];
                                const register double C_1_0 = C_p[INDEX2(1,0,3)];
                                const register double C_2_0 = C_p[INDEX2(2,0,3)];
                                const register double C_0_1 = C_p[INDEX2(0,1,3)];
                                const register double C_1_1 = C_p[INDEX2(1,1,3)];
                                const register double C_2_1 = C_p[INDEX2(2,1,3)];
                                const register double C_0_2 = C_p[INDEX2(0,2,3)];
                                const register double C_1_2 = C_p[INDEX2(1,2,3)];
                                const register double C_2_2 = C_p[INDEX2(2,2,3)];
                                const register double C_0_3 = C_p[INDEX2(0,3,3)];
                                const register double C_1_3 = C_p[INDEX2(1,3,3)];
                                const register double C_2_3 = C_p[INDEX2(2,3,3)];
                                const register double C_0_4 = C_p[INDEX2(0,4,3)];
                                const register double C_1_4 = C_p[INDEX2(1,4,3)];
                                const register double C_2_4 = C_p[INDEX2(2,4,3)];
                                const register double C_0_5 = C_p[INDEX2(0,5,3)];
                                const register double C_1_5 = C_p[INDEX2(1,5,3)];
                                const register double C_2_5 = C_p[INDEX2(2,5,3)];
                                const register double C_0_6 = C_p[INDEX2(0,6,3)];
                                const register double C_1_6 = C_p[INDEX2(1,6,3)];
                                const register double C_2_6 = C_p[INDEX2(2,6,3)];
                                const register double C_0_7 = C_p[INDEX2(0,7,3)];
                                const register double C_1_7 = C_p[INDEX2(1,7,3)];
                                const register double C_2_7 = C_p[INDEX2(2,7,3)];
                                const register double tmp11_0 = C_0_5 + C_0_7;
                                const register double tmp39_0 = C_0_2 + C_0_4;
                                const register double tmp40_0 = C_1_1 + C_1_4;
                                const register double tmp14_0 = C_0_4 + C_0_6;
                                const register double tmp34_0 = C_1_0 + C_1_1 + C_1_4 + C_1_5;
                                const register double tmp23_0 = C_1_0 + C_1_5;
                                const register double tmp12_0 = C_1_4 + C_1_5;
                                const register double tmp35_0 = C_1_2 + C_1_3 + C_1_6 + C_1_7;
                                const register double tmp16_0 = C_1_6 + C_1_7;
                                const register double tmp19_0 = C_2_0 + C_2_1 + C_2_2 + C_2_3;
                                const register double tmp2_0 = C_1_3 + C_1_7;
                                const register double tmp31_0 = C_2_4 + C_2_5;
                                const register double tmp22_0 = C_2_0 + C_2_2;
                                const register double tmp37_0 = C_2_4 + C_2_7;
                                const register double tmp6_0 = C_2_1 + C_2_2;
                                const register double tmp29_0 = C_2_0 + C_2_1;
                                const register double tmp27_0 = C_0_1 + C_0_7;
                                const register double tmp1_0 = C_0_2 + C_0_6;
                                const register double tmp4_0 = C_0_3 + C_0_7;
                                const register double tmp38_0 = C_0_3 + C_0_5;
                                const register double tmp26_0 = C_0_0 + C_0_6;
                                const register double tmp33_0 = C_0_0 + C_0_2 + C_0_4 + C_0_6;
                                const register double tmp10_0 = C_0_0 + C_0_2;
                                const register double tmp25_0 = C_1_2 + C_1_7;
                                const register double tmp7_0 = C_1_1 + C_1_5;
                                const register double tmp18_0 = C_1_0 + C_1_1;
                                const register double tmp8_0 = C_0_0 + C_0_4;
                                const register double tmp13_0 = C_2_4 + C_2_5 + C_2_6 + C_2_7;
                                const register double tmp28_0 = C_2_2 + C_2_3;
                                const register double tmp17_0 = C_0_1 + C_0_3;
                                const register double tmp36_0 = C_2_0 + C_2_3;
                                const register double tmp9_0 = C_1_2 + C_1_6;
                                const register double tmp20_0 = C_2_1 + C_2_3;
                                const register double tmp41_0 = C_1_3 + C_1_6;
                                const register double tmp5_0 = C_1_0 + C_1_4;
                                const register double tmp3_0 = C_2_5 + C_2_6;
                                const register double tmp15_0 = C_1_2 + C_1_3;
                                const register double tmp24_0 = C_2_4 + C_2_6;
                                const register double tmp32_0 = C_0_1 + C_0_3 + C_0_5 + C_0_7;
                                const register double tmp0_0 = C_0_1 + C_0_5;
                                const register double tmp21_0 = C_2_5 + C_2_7;
                                const register double tmp30_0 = C_2_6 + C_2_7;
                                const register double tmp28_1 = tmp21_0*w117;
                                const register double tmp404_1 = C_0_4*w90;
                                const register double tmp153_1 = tmp5_0*w111;
                                const register double tmp208_1 = C_2_0*w103;
                                const register double tmp201_1 = tmp38_0*w108;
                                const register double tmp102_1 = tmp4_0*w93;
                                const register double tmp277_1 = tmp41_0*w97;
                                const register double tmp39_1 = tmp11_0*w108;
                                const register double tmp443_1 = C_1_1*w111;
                                const register double tmp61_1 = tmp15_0*w111;
                                const register double tmp188_1 = tmp17_0*w108;
                                const register double tmp139_1 = tmp35_0*w94;
                                const register double tmp378_1 = C_0_4*w98;
                                const register double tmp10_1 = tmp7_0*w94;
                                const register double tmp31_1 = tmp22_0*w112;
                                const register double tmp209_1 = C_1_4*w114;
                                const register double tmp369_1 = C_2_0*w100;
                                const register double tmp213_1 = C_2_4*w100;
                                const register double tmp254_1 = tmp38_0*w96;
                                const register double tmp350_1 = C_0_1*w124;
                                const register double tmp183_1 = C_0_4*w120;
                                const register double tmp182_1 = C_0_5*w90;
                                const register double tmp83_1 = tmp14_0*w108;
                                const register double tmp247_1 = tmp41_0*w94;
                                const register double tmp79_1 = tmp34_0*w94;
                                const register double tmp332_1 = C_1_7*w115;
                                const register double tmp75_1 = tmp31_0*w103;
                                const register double tmp303_1 = C_1_1*w114;
                                const register double tmp304_1 = C_2_5*w113;
                                const register double tmp420_1 = C_0_1*w125;
                                const register double tmp271_1 = tmp20_0*w95;
                                const register double tmp444_1 = C_1_2*w99;
                                const register double tmp381_1 = C_1_0*w105;
                                const register double tmp96_1 = tmp5_0*w94;
                                const register double tmp385_1 = tmp30_0*w92;
                                const register double tmp371_1 = C_2_4*w103;
                                const register double tmp66_1 = tmp18_0*w107;
                                const register double tmp300_1 = C_0_1*w121;
                                const register double tmp379_1 = C_2_2*w101;
                                const register double tmp283_1 = tmp21_0*w113;
                                const register double tmp242_1 = C_1_7*w91;
                                const register double tmp234_1 = tmp30_0*w95;
                                const register double tmp51_1 = tmp24_0*w103;
                                const register double tmp151_1 = tmp2_0*w105;
                                const register double tmp46_1 = tmp22_0*w95;
                                const register double tmp120_1 = C_2_4*w101;
                                const register double tmp38_1 = tmp10_0*w106;
                                const register double tmp413_1 = tmp15_0*w109;
                                const register double tmp173_1 = C_2_7*w113;
                                const register double tmp116_1 = C_2_7*w103;
                                const register double tmp221_1 = C_0_4*w125;
                                const register double tmp158_1 = tmp29_0*w92;
                                const register double tmp30_1 = C_1_4*w91;
                                const register double tmp154_1 = tmp38_0*w93;
                                const register double tmp406_1 = tmp5_0*w99;
                                const register double tmp204_1 = C_0_0*w124;
                                const register double tmp452_1 = C_2_1*w122;
                                const register double tmp117_1 = C_2_0*w92;
                                const register double tmp323_1 = C_0_3*w125;
                                const register double tmp225_1 = tmp4_0*w96;
                                const register double tmp274_1 = tmp40_0*w94;
                                const register double tmp71_1 = C_0_2*w124;
                                const register double tmp43_1 = tmp23_0*w107;
                                const register double tmp408_1 = tmp9_0*w94;
                                const register double tmp427_1 = tmp24_0*w92;
                                const register double tmp82_1 = tmp17_0*w106;
                                const register double tmp168_1 = tmp26_0*w93;
                                const register double tmp276_1 = tmp10_0*w90;
                                const register double tmp336_1 = C_2_3*w116;
                                const register double tmp53_1 = tmp11_0*w106;
                                const register double tmp161_1 = C_0_7*w98;
                                const register double tmp382_1 = C_0_3*w90;
                                const register double tmp27_1 = tmp20_0*w113;
                                const register double tmp73_1 = C_0_5*w125;
                                const register double tmp389_1 = tmp33_0*w108;
                                const register double tmp288_1 = tmp22_0*w116;
                                const register double tmp258_1 = C_2_7*w100;
                                const register double tmp295_1 = tmp29_0*w117;
                                const register double tmp64_1 = C_0_4*w110;
                                const register double tmp250_1 = tmp40_0*w97;
                                const register double tmp338_1 = C_0_0*w125;
                                const register double tmp359_1 = tmp7_0*w109;
                                const register double tmp400_1 = C_1_5*w111;
                                const register double tmp414_1 = tmp29_0*w113;
                                const register double tmp214_1 = C_1_3*w115;
                                const register double tmp237_1 = C_0_4*w104;
                                const register double tmp7_1 = tmp4_0*w98;
                                const register double tmp145_1 = C_2_3*w122;
                                const register double tmp129_1 = tmp15_0*w94;
                                const register double tmp45_1 = C_1_1*w115;
                                const register double tmp124_1 = tmp14_0*w93;
                                const register double tmp284_1 = tmp20_0*w117;
                                const register double tmp148_1 = C_2_0*w113;
                                const register double tmp334_1 = C_0_7*w124;
                                const register double tmp257_1 = C_1_5*w99;
                                const register double tmp186_1 = tmp11_0*w104;
                                const register double tmp344_1 = C_2_6*w100;
                                const register double tmp433_1 = C_0_7*w120;
                                const register double tmp169_1 = C_0_3*w121;
                                const register double tmp262_1 = C_1_7*w119;
                                const register double tmp47_1 = tmp17_0*w104;
                                const register double tmp278_1 = C_1_2*w118;
                                const register double tmp217_1 = C_2_3*w101;
                                const register double tmp33_1 = C_1_1*w119;
                                const register double tmp157_1 = tmp30_0*w103;
                                const register double tmp44_1 = tmp25_0*w109;
                                const register double tmp454_1 = C_2_6*w123;
                                const register double tmp205_1 = C_0_7*w125;
                                const register double tmp373_1 = C_1_2*w115;
                                const register double tmp361_1 = tmp2_0*w111;
                                const register double tmp392_1 = tmp24_0*w113;
                                const register double tmp308_1 = C_2_2*w116;
                                const register double tmp215_1 = tmp3_0*w95;
                                const register double tmp417_1 = C_0_6*w124;
                                const register double tmp416_1 = tmp28_0*w112;
                                const register double tmp210_1 = C_2_7*w92;
                                const register double tmp335_1 = tmp41_0*w107;
                                const register double tmp429_1 = tmp22_0*w102;
                                const register double tmp315_1 = tmp28_0*w113;
                                const register double tmp435_1 = tmp19_0*w112;
                                const register double tmp407_1 = tmp2_0*w91;
                                const register double tmp365_1 = C_0_5*w98;
                                const register double tmp136_1 = tmp30_0*w113;
                                const register double tmp218_1 = C_1_6*w111;
                                const register double tmp171_1 = C_2_4*w122;
                                const register double tmp253_1 = C_0_7*w121;
                                const register double tmp123_1 = tmp35_0*w109;
                                const register double tmp149_1 = C_2_4*w123;
                                const register double tmp291_1 = tmp0_0*w96;
                                const register double tmp267_1 = tmp14_0*w96;
                                const register double tmp147_1 = tmp9_0*w109;
                                const register double tmp296_1 = tmp8_0*w98;
                                const register double tmp42_1 = tmp20_0*w92;
                                const register double tmp430_1 = C_0_0*w121;
                                const register double tmp399_1 = C_1_7*w114;
                                const register double tmp355_1 = tmp4_0*w106;
                                const register double tmp56_1 = tmp16_0*w107;
                                const register double tmp289_1 = tmp5_0*w97;
                                const register double tmp125_1 = tmp17_0*w96;
                                const register double tmp397_1 = tmp9_0*w111;
                                const register double tmp333_1 = C_2_0*w123;
                                const register double tmp15_1 = tmp11_0*w96;
                                const register double tmp200_1 = tmp39_0*w106;
                                const register double tmp285_1 = tmp9_0*w91;
                                const register double tmp374_1 = C_2_1*w103;
                                const register double tmp354_1 = tmp1_0*w104;
                                const register double tmp97_1 = tmp2_0*w97;
                                const register double tmp104_1 = C_2_2*w122;
                                const register double tmp339_1 = C_1_2*w111;
                                const register double tmp437_1 = C_1_1*w91;
                                const register double tmp259_1 = C_0_6*w98;
                                const register double tmp110_1 = tmp0_0*w98;
                                const register double tmp419_1 = tmp30_0*w116;
                                const register double tmp236_1 = C_1_5*w118;
                                const register double tmp41_1 = C_1_6*w114;
                                const register double tmp320_1 = C_0_4*w124;
                                const register double tmp74_1 = tmp12_0*w111;
                                const register double tmp63_1 = tmp27_0*w108;
                                const register double tmp106_1 = tmp2_0*w109;
                                const register double tmp150_1 = tmp4_0*w104;
                                const register double tmp195_1 = tmp20_0*w112;
                                const register double tmp270_1 = C_1_7*w99;
                                const register double tmp230_1 = C_2_1*w100;
                                const register double tmp68_1 = tmp29_0*w95;
                                const register double tmp287_1 = tmp2_0*w94;
                                const register double tmp132_1 = tmp28_0*w117;
                                const register double tmp423_1 = tmp12_0*w109;
                                const register double tmp272_1 = tmp11_0*w98;
                                const register double tmp134_1 = tmp1_0*w108;
                                const register double tmp1_1 = tmp1_0*w96;
                                const register double tmp301_1 = C_1_6*w115;
                                const register double tmp241_1 = C_2_6*w113;
                                const register double tmp29_1 = tmp17_0*w98;
                                const register double tmp337_1 = C_1_5*w105;
                                const register double tmp92_1 = tmp36_0*w95;
                                const register double tmp40_1 = tmp14_0*w110;
                                const register double tmp203_1 = C_0_1*w104;
                                const register double tmp349_1 = C_2_1*w101;
                                const register double tmp78_1 = tmp19_0*w95;
                                const register double tmp266_1 = tmp17_0*w93;
                                const register double tmp58_1 = tmp19_0*w117;
                                const register double tmp310_1 = C_0_6*w120;
                                const register double tmp127_1 = tmp10_0*w98;
                                const register double tmp307_1 = C_0_7*w90;
                                const register double tmp330_1 = tmp40_0*w109;
                                const register double tmp386_1 = tmp31_0*w95;
                                const register double tmp351_1 = C_1_1*w118;
                                const register double tmp306_1 = C_0_0*w98;
                                const register double tmp312_1 = tmp30_0*w117;
                                const register double tmp108_1 = C_2_5*w123;
                                const register double tmp140_1 = tmp29_0*w116;
                                const register double tmp235_1 = tmp29_0*w102;
                                const register double tmp293_1 = tmp34_0*w109;
                                const register double tmp446_1 = C_1_0*w119;
                                const register double tmp160_1 = tmp28_0*w95;
                                const register double tmp393_1 = tmp22_0*w117;
                                const register double tmp231_1 = C_2_6*w101;
                                const register double tmp233_1 = tmp31_0*w92;
                                const register double tmp8_1 = tmp5_0*w91;
                                const register double tmp229_1 = C_2_2*w92;
                                const register double tmp128_1 = tmp16_0*w91;
                                const register double tmp187_1 = tmp14_0*w106;
                                const register double tmp86_1 = tmp8_0*w106;
                                const register double tmp49_1 = C_1_3*w105;
                                const register double tmp447_1 = C_1_7*w118;
                                const register double tmp442_1 = C_1_6*w105;
                                const register double tmp77_1 = tmp33_0*w96;
                                const register double tmp358_1 = tmp9_0*w107;
                                const register double tmp292_1 = tmp35_0*w107;
                                const register double tmp81_1 = tmp13_0*w102;
                                const register double tmp206_1 = tmp27_0*w106;
                                const register double tmp366_1 = C_0_2*w90;
                                const register double tmp448_1 = tmp22_0*w103;
                                const register double tmp156_1 = tmp39_0*w96;
                                const register double tmp105_1 = tmp5_0*w107;
                                const register double tmp179_1 = tmp25_0*w94;
                                const register double tmp364_1 = tmp26_0*w96;
                                const register double tmp449_1 = tmp21_0*w92;
                                const register double tmp453_1 = C_2_2*w113;
                                const register double tmp12_1 = tmp8_0*w90;
                                const register double tmp59_1 = tmp13_0*w112;
                                const register double tmp99_1 = C_2_6*w103;
                                const register double tmp327_1 = C_0_1*w110;
                                const register double tmp50_1 = C_1_4*w111;
                                const register double tmp126_1 = tmp18_0*w99;
                                const register double tmp445_1 = C_1_5*w91;
                                const register double tmp152_1 = C_2_7*w116;
                                const register double tmp76_1 = tmp32_0*w93;
                                const register double tmp80_1 = tmp35_0*w97;
                                const register double tmp146_1 = tmp7_0*w107;
                                const register double tmp197_1 = tmp21_0*w116;
                                const register double tmp111_1 = tmp36_0*w112;
                                const register double tmp115_1 = tmp7_0*w111;
                                const register double tmp396_1 = tmp7_0*w105;
                                const register double tmp248_1 = C_1_2*w119;
                                const register double tmp390_1 = tmp2_0*w107;
                                const register double tmp119_1 = tmp6_0*w95;
                                const register double tmp440_1 = C_1_4*w115;
                                const register double tmp6_1 = tmp3_0*w117;
                                const register double tmp346_1 = C_0_6*w125;
                                const register double tmp144_1 = tmp8_0*w110;
                                const register double tmp249_1 = C_2_1*w116;
                                const register double tmp318_1 = tmp29_0*w112;
                                const register double tmp265_1 = C_1_0*w118;
                                const register double tmp387_1 = tmp28_0*w102;
                                const register double tmp130_1 = tmp11_0*w90;
                                const register double tmp348_1 = C_1_3*w91;
                                const register double tmp252_1 = tmp39_0*w93;
                                const register double tmp85_1 = tmp10_0*w104;
                                const register double tmp21_1 = tmp17_0*w90;
                                const register double tmp90_1 = tmp9_0*w99;
                                const register double tmp232_1 = tmp28_0*w103;
                                const register double tmp441_1 = C_1_3*w114;
                                const register double tmp360_1 = tmp5_0*w105;
                                const register double tmp4_1 = C_2_3*w113;
                                const register double tmp313_1 = C_0_2*w110;
                                const register double tmp135_1 = tmp4_0*w110;
                                const register double tmp155_1 = C_0_6*w121;
                                const register double tmp20_1 = tmp16_0*w94;
                                const register double tmp432_1 = C_0_6*w90;
                                const register double tmp227_1 = tmp0_0*w90;
                                const register double tmp269_1 = tmp22_0*w92;
                                const register double tmp261_1 = C_2_0*w101;
                                const register double tmp19_1 = tmp15_0*w91;
                                const register double tmp16_1 = tmp12_0*w99;
                                const register double tmp395_1 = tmp20_0*w116;
                                const register double tmp356_1 = tmp8_0*w108;
                                const register double tmp222_1 = tmp6_0*w102;
                                const register double tmp281_1 = tmp32_0*w96;
                                const register double tmp164_1 = C_0_0*w90;
                                const register double tmp375_1 = C_1_5*w114;
                                const register double tmp238_1 = C_0_3*w110;
                                const register double tmp340_1 = C_2_2*w103;
                                const register double tmp290_1 = tmp1_0*w93;
                                const register double tmp311_1 = C_1_3*w111;
                                const register double tmp162_1 = tmp18_0*w91;
                                const register double tmp107_1 = C_2_1*w113;
                                const register double tmp401_1 = C_1_2*w105;
                                const register double tmp87_1 = tmp4_0*w108;
                                const register double tmp172_1 = C_1_1*w99;
                                const register double tmp402_1 = C_0_2*w121;
                                const register double tmp434_1 = tmp13_0*w117;
                                const register double tmp394_1 = tmp21_0*w112;
                                const register double tmp191_1 = tmp41_0*w109;
                                const register double tmp325_1 = tmp38_0*w106;
                                const register double tmp91_1 = C_2_2*w100;
                                const register double tmp35_1 = tmp14_0*w90;
                                const register double tmp279_1 = tmp24_0*w102;
                                const register double tmp436_1 = C_1_6*w99;
                                const register double tmp244_1 = tmp36_0*w117;
                                const register double tmp264_1 = C_0_0*w120;
                                const register double tmp415_1 = tmp31_0*w117;
                                const register double tmp178_1 = tmp3_0*w112;
                                const register double tmp239_1 = C_2_5*w122;
                                const register double tmp219_1 = C_0_3*w124;
                                const register double tmp370_1 = C_2_7*w101;
                                const register double tmp109_1 = tmp37_0*w117;
                                const register double tmp316_1 = C_0_5*w104;
                                const register double tmp174_1 = C_2_3*w123;
                                const register double tmp89_1 = C_2_1*w92;
                                const register double tmp228_1 = C_2_5*w103;
                                const register double tmp36_1 = tmp25_0*w97;
                                const register double tmp412_1 = tmp12_0*w107;
                                const register double tmp326_1 = tmp39_0*w108;
                                const register double tmp256_1 = C_2_4*w92;
                                const register double tmp23_1 = tmp19_0*w102;
                                const register double tmp362_1 = tmp27_0*w93;
                                const register double tmp286_1 = tmp24_0*w112;
                                const register double tmp114_1 = C_2_6*w116;
                                const register double tmp176_1 = C_0_2*w98;
                                const register double tmp84_1 = tmp11_0*w110;
                                const register double tmp328_1 = C_2_7*w122;
                                const register double tmp113_1 = tmp1_0*w90;
                                const register double tmp55_1 = tmp17_0*w110;
                                const register double tmp319_1 = tmp18_0*w94;
                                const register double tmp95_1 = C_2_5*w101;
                                const register double tmp72_1 = tmp15_0*w105;
                                const register double tmp220_1 = C_1_1*w105;
                                const register double tmp60_1 = tmp12_0*w105;
                                const register double tmp118_1 = C_2_3*w100;
                                const register double tmp32_1 = tmp23_0*w94;
                                const register double tmp166_1 = tmp15_0*w97;
                                const register double tmp159_1 = tmp16_0*w99;
                                const register double tmp421_1 = tmp18_0*w111;
                                const register double tmp367_1 = C_0_3*w120;
                                const register double tmp424_1 = tmp16_0*w111;
                                const register double tmp431_1 = C_0_1*w98;
                                const register double tmp165_1 = C_0_1*w120;
                                const register double tmp309_1 = C_1_4*w105;
                                const register double tmp368_1 = C_2_3*w92;
                                const register double tmp223_1 = C_0_5*w110;
                                const register double tmp384_1 = tmp29_0*w103;
                                const register double tmp143_1 = tmp0_0*w108;
                                const register double tmp240_1 = C_1_0*w99;
                                const register double tmp255_1 = C_2_3*w103;
                                const register double tmp122_1 = tmp34_0*w107;
                                const register double tmp388_1 = tmp32_0*w106;
                                const register double tmp185_1 = C_1_4*w118;
                                const register double tmp88_1 = tmp1_0*w110;
                                const register double tmp342_1 = C_2_5*w92;
                                const register double tmp196_1 = C_1_2*w114;
                                const register double tmp142_1 = tmp1_0*w106;
                                const register double tmp131_1 = tmp12_0*w97;
                                const register double tmp418_1 = tmp16_0*w105;
                                const register double tmp263_1 = C_0_1*w90;
                                const register double tmp17_1 = tmp13_0*w95;
                                const register double tmp422_1 = tmp15_0*w107;
                                const register double tmp121_1 = tmp3_0*w102;
                                const register double tmp190_1 = tmp40_0*w107;
                                const register double tmp70_1 = tmp30_0*w102;
                                const register double tmp363_1 = C_0_4*w121;
                                const register double tmp194_1 = tmp24_0*w117;
                                const register double tmp297_1 = tmp30_0*w112;
                                const register double tmp280_1 = tmp33_0*w93;
                                const register double tmp391_1 = tmp5_0*w109;
                                const register double tmp439_1 = C_1_3*w118;
                                const register double tmp69_1 = C_0_3*w104;
                                const register double tmp48_1 = tmp21_0*w102;
                                const register double tmp243_1 = C_2_2*w123;
                                const register double tmp405_1 = C_0_5*w120;
                                const register double tmp2_1 = C_2_0*w122;
                                const register double tmp372_1 = C_0_5*w121;
                                const register double tmp100_1 = tmp33_0*w106;
                                const register double tmp352_1 = tmp36_0*w102;
                                const register double tmp317_1 = tmp12_0*w91;
                                const register double tmp141_1 = tmp34_0*w97;
                                const register double tmp282_1 = tmp7_0*w99;
                                const register double tmp181_1 = C_2_0*w116;
                                const register double tmp177_1 = C_1_6*w91;
                                const register double tmp14_1 = tmp10_0*w93;
                                const register double tmp321_1 = tmp31_0*w116;
                                const register double tmp98_1 = tmp37_0*w102;
                                const register double tmp9_1 = tmp6_0*w112;
                                const register double tmp322_1 = tmp16_0*w97;
                                const register double tmp103_1 = tmp8_0*w96;
                                const register double tmp0_1 = tmp0_0*w93;
                                const register double tmp3_1 = tmp2_0*w99;
                                const register double tmp216_1 = C_0_2*w104;
                                const register double tmp11_1 = C_2_4*w116;
                                const register double tmp189_1 = tmp10_0*w110;
                                const register double tmp275_1 = C_1_5*w119;
                                const register double tmp376_1 = C_2_6*w92;
                                const register double tmp170_1 = tmp27_0*w96;
                                const register double tmp101_1 = tmp32_0*w108;
                                const register double tmp353_1 = C_0_7*w110;
                                const register double tmp13_1 = tmp9_0*w97;
                                const register double tmp133_1 = tmp0_0*w106;
                                const register double tmp251_1 = C_0_2*w125;
                                const register double tmp62_1 = tmp26_0*w106;
                                const register double tmp137_1 = tmp8_0*w104;
                                const register double tmp294_1 = tmp31_0*w113;
                                const register double tmp260_1 = C_1_2*w91;
                                const register double tmp409_1 = tmp7_0*w97;
                                const register double tmp24_1 = tmp11_0*w93;
                                const register double tmp314_1 = tmp15_0*w99;
                                const register double tmp26_1 = C_1_3*w99;
                                const register double tmp211_1 = tmp25_0*w107;
                                const register double tmp331_1 = C_2_4*w113;
                                const register double tmp199_1 = C_1_0*w111;
                                const register double tmp22_1 = tmp18_0*w97;
                                const register double tmp34_1 = tmp24_0*w116;
                                const register double tmp192_1 = tmp22_0*w113;
                                const register double tmp67_1 = tmp16_0*w109;
                                const register double tmp426_1 = tmp20_0*w103;
                                const register double tmp226_1 = tmp1_0*w98;
                                const register double tmp377_1 = C_2_5*w100;
                                const register double tmp94_1 = tmp7_0*w91;
                                const register double tmp345_1 = tmp37_0*w95;
                                const register double tmp411_1 = C_0_0*w110;
                                const register double tmp438_1 = C_1_4*w119;
                                const register double tmp112_1 = tmp9_0*w105;
                                const register double tmp54_1 = tmp10_0*w108;
                                const register double tmp324_1 = C_0_6*w104;
                                const register double tmp93_1 = tmp0_0*w104;
                                const register double tmp37_1 = C_1_6*w118;
                                const register double tmp329_1 = C_1_0*w114;
                                const register double tmp357_1 = tmp0_0*w110;
                                const register double tmp180_1 = C_1_3*w119;
                                const register double tmp273_1 = C_1_0*w91;
                                const register double tmp65_1 = tmp28_0*w92;
                                const register double tmp299_1 = tmp4_0*w90;
                                const register double tmp184_1 = tmp23_0*w97;
                                const register double tmp245_1 = C_0_5*w124;
                                const register double tmp341_1 = C_1_6*w119;
                                const register double tmp224_1 = tmp8_0*w93;
                                const register double tmp202_1 = C_0_6*w110;
                                const register double tmp302_1 = C_2_6*w122;
                                const register double tmp18_1 = tmp14_0*w98;
                                const register double tmp193_1 = C_1_5*w115;
                                const register double tmp25_1 = tmp10_0*w96;
                                const register double tmp398_1 = C_1_0*w115;
                                const register double tmp198_1 = C_1_7*w105;
                                const register double tmp380_1 = C_1_7*w111;
                                const register double tmp450_1 = tmp24_0*w95;
                                const register double tmp138_1 = tmp31_0*w112;
                                const register double tmp403_1 = C_0_3*w98;
                                const register double tmp212_1 = tmp23_0*w109;
                                const register double tmp246_1 = tmp37_0*w112;
                                const register double tmp167_1 = tmp31_0*w102;
                                const register double tmp57_1 = tmp18_0*w109;
                                const register double tmp5_1 = C_2_7*w123;
                                const register double tmp207_1 = tmp26_0*w108;
                                const register double tmp52_1 = tmp14_0*w104;
                                const register double tmp175_1 = tmp6_0*w117;
                                const register double tmp425_1 = tmp18_0*w105;
                                const register double tmp343_1 = C_1_4*w99;
                                const register double tmp383_1 = C_0_2*w120;
                                const register double tmp163_1 = tmp12_0*w94;
                                const register double tmp455_1 = C_2_5*w116;
                                const register double tmp428_1 = tmp21_0*w95;
                                const register double tmp347_1 = C_0_0*w104;
                                const register double tmp451_1 = tmp20_0*w102;
                                const register double tmp298_1 = tmp28_0*w116;
                                const register double tmp410_1 = C_0_7*w104;
                                const register double tmp268_1 = tmp21_0*w103;
                                const register double tmp305_1 = C_2_1*w123;
                                EM_S[INDEX2(7,3,8)]+=tmp0_1 + tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                EM_S[INDEX2(4,7,8)]+=tmp14_1 + tmp15_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1;
                                EM_S[INDEX2(1,3,8)]+=tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp33_1 + tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1;
                                EM_S[INDEX2(6,4,8)]+=tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                                EM_S[INDEX2(3,0,8)]+=tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1 + tmp60_1 + tmp61_1;
                                EM_S[INDEX2(5,4,8)]+=tmp62_1 + tmp63_1 + tmp64_1 + tmp65_1 + tmp66_1 + tmp67_1 + tmp68_1 + tmp69_1 + tmp70_1 + tmp71_1 + tmp72_1 + tmp73_1 + tmp74_1 + tmp75_1;
                                EM_S[INDEX2(0,7,8)]+=tmp76_1 + tmp77_1 + tmp78_1 + tmp79_1 + tmp80_1 + tmp81_1;
                                EM_S[INDEX2(5,6,8)]+=tmp16_1 + tmp17_1 + tmp19_1 + tmp20_1 + tmp22_1 + tmp23_1 + tmp82_1 + tmp83_1 + tmp84_1 + tmp85_1;
                                EM_S[INDEX2(2,6,8)]+=tmp86_1 + tmp87_1 + tmp88_1 + tmp89_1 + tmp90_1 + tmp91_1 + tmp92_1 + tmp93_1 + tmp94_1 + tmp95_1 + tmp96_1 + tmp97_1 + tmp98_1 + tmp99_1;
                                EM_S[INDEX2(1,6,8)]+=tmp100_1 + tmp101_1 + tmp78_1 + tmp79_1 + tmp80_1 + tmp81_1;
                                EM_S[INDEX2(5,1,8)]+=tmp102_1 + tmp103_1 + tmp104_1 + tmp105_1 + tmp106_1 + tmp107_1 + tmp108_1 + tmp109_1 + tmp110_1 + tmp111_1 + tmp112_1 + tmp113_1 + tmp114_1 + tmp115_1;
                                EM_S[INDEX2(3,7,8)]+=tmp0_1 + tmp10_1 + tmp116_1 + tmp117_1 + tmp118_1 + tmp119_1 + tmp120_1 + tmp121_1 + tmp12_1 + tmp13_1 + tmp1_1 + tmp3_1 + tmp7_1 + tmp8_1;
                                EM_S[INDEX2(2,5,8)]+=tmp122_1 + tmp123_1 + tmp76_1 + tmp77_1 + tmp78_1 + tmp81_1;
                                EM_S[INDEX2(0,3,8)]+=tmp124_1 + tmp125_1 + tmp126_1 + tmp127_1 + tmp128_1 + tmp129_1 + tmp130_1 + tmp131_1 + tmp58_1 + tmp59_1;
                                EM_S[INDEX2(7,2,8)]+=tmp132_1 + tmp133_1 + tmp134_1 + tmp135_1 + tmp136_1 + tmp137_1 + tmp138_1 + tmp139_1 + tmp140_1 + tmp141_1;
                                EM_S[INDEX2(4,0,8)]+=tmp142_1 + tmp143_1 + tmp144_1 + tmp145_1 + tmp146_1 + tmp147_1 + tmp148_1 + tmp149_1 + tmp150_1 + tmp151_1 + tmp152_1 + tmp153_1 + tmp6_1 + tmp9_1;
                                EM_S[INDEX2(1,2,8)]+=tmp126_1 + tmp128_1 + tmp129_1 + tmp131_1 + tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1 + tmp58_1 + tmp59_1;
                                EM_S[INDEX2(6,7,8)]+=tmp154_1 + tmp155_1 + tmp156_1 + tmp157_1 + tmp158_1 + tmp159_1 + tmp160_1 + tmp161_1 + tmp162_1 + tmp163_1 + tmp164_1 + tmp165_1 + tmp166_1 + tmp167_1;
                                EM_S[INDEX2(3,3,8)]+=tmp168_1 + tmp169_1 + tmp170_1 + tmp171_1 + tmp172_1 + tmp173_1 + tmp174_1 + tmp175_1 + tmp176_1 + tmp177_1 + tmp178_1 + tmp179_1 + tmp180_1 + tmp181_1 + tmp182_1 + tmp183_1 + tmp184_1 + tmp185_1;
                                EM_S[INDEX2(2,0,8)]+=tmp186_1 + tmp187_1 + tmp188_1 + tmp189_1 + tmp190_1 + tmp191_1 + tmp192_1 + tmp193_1 + tmp194_1 + tmp195_1 + tmp196_1 + tmp197_1 + tmp198_1 + tmp199_1;
                                EM_S[INDEX2(7,6,8)]+=tmp157_1 + tmp158_1 + tmp159_1 + tmp160_1 + tmp162_1 + tmp163_1 + tmp166_1 + tmp167_1 + tmp200_1 + tmp201_1 + tmp202_1 + tmp203_1 + tmp204_1 + tmp205_1;
                                EM_S[INDEX2(4,4,8)]+=tmp206_1 + tmp207_1 + tmp208_1 + tmp209_1 + tmp210_1 + tmp211_1 + tmp212_1 + tmp213_1 + tmp214_1 + tmp215_1 + tmp216_1 + tmp217_1 + tmp218_1 + tmp219_1 + tmp220_1 + tmp221_1 + tmp222_1 + tmp223_1;
                                EM_S[INDEX2(6,3,8)]+=tmp132_1 + tmp136_1 + tmp138_1 + tmp139_1 + tmp140_1 + tmp141_1 + tmp224_1 + tmp225_1 + tmp226_1 + tmp227_1;
                                EM_S[INDEX2(1,5,8)]+=tmp102_1 + tmp103_1 + tmp105_1 + tmp106_1 + tmp110_1 + tmp112_1 + tmp113_1 + tmp115_1 + tmp228_1 + tmp229_1 + tmp230_1 + tmp231_1 + tmp92_1 + tmp98_1;
                                EM_S[INDEX2(3,6,8)]+=tmp133_1 + tmp134_1 + tmp135_1 + tmp137_1 + tmp139_1 + tmp141_1 + tmp232_1 + tmp233_1 + tmp234_1 + tmp235_1;
                                EM_S[INDEX2(2,2,8)]+=tmp206_1 + tmp207_1 + tmp236_1 + tmp237_1 + tmp238_1 + tmp239_1 + tmp240_1 + tmp241_1 + tmp242_1 + tmp243_1 + tmp244_1 + tmp245_1 + tmp246_1 + tmp247_1 + tmp248_1 + tmp249_1 + tmp250_1 + tmp251_1;
                                EM_S[INDEX2(7,7,8)]+=tmp215_1 + tmp222_1 + tmp247_1 + tmp250_1 + tmp252_1 + tmp253_1 + tmp254_1 + tmp255_1 + tmp256_1 + tmp257_1 + tmp258_1 + tmp259_1 + tmp260_1 + tmp261_1 + tmp262_1 + tmp263_1 + tmp264_1 + tmp265_1;
                                EM_S[INDEX2(5,7,8)]+=tmp266_1 + tmp267_1 + tmp268_1 + tmp269_1 + tmp270_1 + tmp271_1 + tmp272_1 + tmp273_1 + tmp274_1 + tmp275_1 + tmp276_1 + tmp277_1 + tmp278_1 + tmp279_1;
                                EM_S[INDEX2(5,3,8)]+=tmp280_1 + tmp281_1 + tmp282_1 + tmp283_1 + tmp284_1 + tmp285_1 + tmp286_1 + tmp287_1 + tmp288_1 + tmp289_1;
                                EM_S[INDEX2(4,1,8)]+=tmp290_1 + tmp291_1 + tmp292_1 + tmp293_1 + tmp294_1 + tmp295_1 + tmp296_1 + tmp297_1 + tmp298_1 + tmp299_1;
                                EM_S[INDEX2(1,1,8)]+=tmp211_1 + tmp212_1 + tmp244_1 + tmp246_1 + tmp252_1 + tmp254_1 + tmp300_1 + tmp301_1 + tmp302_1 + tmp303_1 + tmp304_1 + tmp305_1 + tmp306_1 + tmp307_1 + tmp308_1 + tmp309_1 + tmp310_1 + tmp311_1;
                                EM_S[INDEX2(2,7,8)]+=tmp139_1 + tmp141_1 + tmp224_1 + tmp225_1 + tmp226_1 + tmp227_1 + tmp232_1 + tmp233_1 + tmp234_1 + tmp235_1;
                                EM_S[INDEX2(3,2,8)]+=tmp312_1 + tmp313_1 + tmp314_1 + tmp315_1 + tmp316_1 + tmp317_1 + tmp318_1 + tmp319_1 + tmp320_1 + tmp321_1 + tmp322_1 + tmp323_1 + tmp62_1 + tmp63_1;
                                EM_S[INDEX2(0,0,8)]+=tmp175_1 + tmp178_1 + tmp324_1 + tmp325_1 + tmp326_1 + tmp327_1 + tmp328_1 + tmp329_1 + tmp330_1 + tmp331_1 + tmp332_1 + tmp333_1 + tmp334_1 + tmp335_1 + tmp336_1 + tmp337_1 + tmp338_1 + tmp339_1;
                                EM_S[INDEX2(6,6,8)]+=tmp179_1 + tmp184_1 + tmp325_1 + tmp326_1 + tmp340_1 + tmp341_1 + tmp342_1 + tmp343_1 + tmp344_1 + tmp345_1 + tmp346_1 + tmp347_1 + tmp348_1 + tmp349_1 + tmp350_1 + tmp351_1 + tmp352_1 + tmp353_1;
                                EM_S[INDEX2(5,0,8)]+=tmp292_1 + tmp293_1 + tmp294_1 + tmp295_1 + tmp297_1 + tmp298_1 + tmp354_1 + tmp355_1 + tmp356_1 + tmp357_1;
                                EM_S[INDEX2(7,1,8)]+=tmp280_1 + tmp281_1 + tmp283_1 + tmp284_1 + tmp286_1 + tmp288_1 + tmp358_1 + tmp359_1 + tmp360_1 + tmp361_1;
                                EM_S[INDEX2(4,5,8)]+=tmp362_1 + tmp363_1 + tmp364_1 + tmp365_1 + tmp366_1 + tmp367_1 + tmp65_1 + tmp66_1 + tmp67_1 + tmp68_1 + tmp70_1 + tmp72_1 + tmp74_1 + tmp75_1;
                                EM_S[INDEX2(0,4,8)]+=tmp119_1 + tmp121_1 + tmp142_1 + tmp143_1 + tmp144_1 + tmp146_1 + tmp147_1 + tmp150_1 + tmp151_1 + tmp153_1 + tmp368_1 + tmp369_1 + tmp370_1 + tmp371_1;
                                EM_S[INDEX2(5,5,8)]+=tmp168_1 + tmp170_1 + tmp330_1 + tmp335_1 + tmp345_1 + tmp352_1 + tmp372_1 + tmp373_1 + tmp374_1 + tmp375_1 + tmp376_1 + tmp377_1 + tmp378_1 + tmp379_1 + tmp380_1 + tmp381_1 + tmp382_1 + tmp383_1;
                                EM_S[INDEX2(1,4,8)]+=tmp292_1 + tmp293_1 + tmp354_1 + tmp355_1 + tmp356_1 + tmp357_1 + tmp384_1 + tmp385_1 + tmp386_1 + tmp387_1;
                                EM_S[INDEX2(6,0,8)]+=tmp388_1 + tmp389_1 + tmp390_1 + tmp391_1 + tmp392_1 + tmp393_1 + tmp394_1 + tmp395_1 + tmp396_1 + tmp397_1;
                                EM_S[INDEX2(7,5,8)]+=tmp190_1 + tmp191_1 + tmp266_1 + tmp267_1 + tmp268_1 + tmp269_1 + tmp271_1 + tmp272_1 + tmp276_1 + tmp279_1 + tmp398_1 + tmp399_1 + tmp400_1 + tmp401_1;
                                EM_S[INDEX2(2,3,8)]+=tmp312_1 + tmp314_1 + tmp315_1 + tmp317_1 + tmp318_1 + tmp319_1 + tmp321_1 + tmp322_1 + tmp362_1 + tmp364_1 + tmp402_1 + tmp403_1 + tmp404_1 + tmp405_1;
                                EM_S[INDEX2(2,1,8)]+=tmp124_1 + tmp125_1 + tmp127_1 + tmp130_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1 + tmp60_1 + tmp61_1;
                                EM_S[INDEX2(4,2,8)]+=tmp388_1 + tmp389_1 + tmp392_1 + tmp393_1 + tmp394_1 + tmp395_1 + tmp406_1 + tmp407_1 + tmp408_1 + tmp409_1;
                                EM_S[INDEX2(1,0,8)]+=tmp200_1 + tmp201_1 + tmp410_1 + tmp411_1 + tmp412_1 + tmp413_1 + tmp414_1 + tmp415_1 + tmp416_1 + tmp417_1 + tmp418_1 + tmp419_1 + tmp420_1 + tmp421_1;
                                EM_S[INDEX2(6,5,8)]+=tmp14_1 + tmp15_1 + tmp17_1 + tmp18_1 + tmp21_1 + tmp23_1 + tmp422_1 + tmp423_1 + tmp424_1 + tmp425_1;
                                EM_S[INDEX2(3,5,8)]+=tmp280_1 + tmp281_1 + tmp358_1 + tmp359_1 + tmp360_1 + tmp361_1 + tmp426_1 + tmp427_1 + tmp428_1 + tmp429_1;
                                EM_S[INDEX2(0,1,8)]+=tmp154_1 + tmp156_1 + tmp412_1 + tmp413_1 + tmp414_1 + tmp415_1 + tmp416_1 + tmp418_1 + tmp419_1 + tmp421_1 + tmp430_1 + tmp431_1 + tmp432_1 + tmp433_1;
                                EM_S[INDEX2(7,0,8)]+=tmp100_1 + tmp101_1 + tmp122_1 + tmp123_1 + tmp434_1 + tmp435_1;
                                EM_S[INDEX2(4,6,8)]+=tmp32_1 + tmp36_1 + tmp38_1 + tmp39_1 + tmp40_1 + tmp42_1 + tmp436_1 + tmp437_1 + tmp438_1 + tmp439_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp51_1;
                                EM_S[INDEX2(5,2,8)]+=tmp100_1 + tmp101_1 + tmp434_1 + tmp435_1 + tmp79_1 + tmp80_1;
                                EM_S[INDEX2(6,1,8)]+=tmp122_1 + tmp123_1 + tmp434_1 + tmp435_1 + tmp76_1 + tmp77_1;
                                EM_S[INDEX2(3,1,8)]+=tmp24_1 + tmp25_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp31_1 + tmp34_1 + tmp35_1 + tmp43_1 + tmp440_1 + tmp441_1 + tmp442_1 + tmp443_1 + tmp44_1;
                                EM_S[INDEX2(0,2,8)]+=tmp186_1 + tmp187_1 + tmp188_1 + tmp189_1 + tmp192_1 + tmp194_1 + tmp195_1 + tmp197_1 + tmp274_1 + tmp277_1 + tmp444_1 + tmp445_1 + tmp446_1 + tmp447_1;
                                EM_S[INDEX2(7,4,8)]+=tmp17_1 + tmp23_1 + tmp422_1 + tmp423_1 + tmp424_1 + tmp425_1 + tmp82_1 + tmp83_1 + tmp84_1 + tmp85_1;
                                EM_S[INDEX2(0,6,8)]+=tmp388_1 + tmp389_1 + tmp406_1 + tmp407_1 + tmp408_1 + tmp409_1 + tmp448_1 + tmp449_1 + tmp450_1 + tmp451_1;
                                EM_S[INDEX2(6,2,8)]+=tmp109_1 + tmp111_1 + tmp452_1 + tmp453_1 + tmp454_1 + tmp455_1 + tmp86_1 + tmp87_1 + tmp88_1 + tmp90_1 + tmp93_1 + tmp94_1 + tmp96_1 + tmp97_1;
                                EM_S[INDEX2(4,3,8)]+=tmp434_1 + tmp435_1 + tmp76_1 + tmp77_1 + tmp79_1 + tmp80_1;
                                EM_S[INDEX2(1,7,8)]+=tmp280_1 + tmp281_1 + tmp282_1 + tmp285_1 + tmp287_1 + tmp289_1 + tmp426_1 + tmp427_1 + tmp428_1 + tmp429_1;
                                EM_S[INDEX2(0,5,8)]+=tmp290_1 + tmp291_1 + tmp292_1 + tmp293_1 + tmp296_1 + tmp299_1 + tmp384_1 + tmp385_1 + tmp386_1 + tmp387_1;
                                EM_S[INDEX2(3,4,8)]+=tmp100_1 + tmp101_1 + tmp122_1 + tmp123_1 + tmp78_1 + tmp81_1;
                                EM_S[INDEX2(2,4,8)]+=tmp388_1 + tmp389_1 + tmp390_1 + tmp391_1 + tmp396_1 + tmp397_1 + tmp448_1 + tmp449_1 + tmp450_1 + tmp451_1;
                            } else { // constant data
                                const register double C_0 = C_p[0];
                                const register double C_1 = C_p[1];
                                const register double C_2 = C_p[2];
                                const register double tmp1_1 = C_1*w127;
                                const register double tmp15_1 = C_0*w138;
                                const register double tmp4_1 = C_2*w133;
                                const register double tmp3_1 = C_2*w131;
                                const register double tmp8_1 = C_1*w132;
                                const register double tmp7_1 = C_0*w129;
                                const register double tmp0_1 = C_2*w140;
                                const register double tmp2_1 = C_0*w126;
                                const register double tmp16_1 = C_1*w139;
                                const register double tmp13_1 = C_0*w141;
                                const register double tmp14_1 = C_2*w128;
                                const register double tmp11_1 = C_0*w143;
                                const register double tmp12_1 = C_1*w142;
                                const register double tmp5_1 = C_1*w134;
                                const register double tmp17_1 = C_0*w137;
                                const register double tmp6_1 = C_2*w135;
                                const register double tmp10_1 = C_2*w136;
                                const register double tmp9_1 = C_1*w130;
                                EM_S[INDEX2(7,3,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                EM_S[INDEX2(4,7,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                EM_S[INDEX2(1,3,8)]+=tmp2_1 + tmp4_1 + tmp5_1;
                                EM_S[INDEX2(6,4,8)]+=tmp6_1 + tmp7_1 + tmp8_1;
                                EM_S[INDEX2(3,0,8)]+=tmp10_1 + tmp7_1 + tmp9_1;
                                EM_S[INDEX2(5,4,8)]+=tmp11_1 + tmp6_1 + tmp9_1;
                                EM_S[INDEX2(0,7,8)]+=tmp12_1 + tmp13_1 + tmp3_1;
                                EM_S[INDEX2(5,6,8)]+=tmp1_1 + tmp3_1 + tmp7_1;
                                EM_S[INDEX2(2,6,8)]+=tmp14_1 + tmp1_1 + tmp7_1;
                                EM_S[INDEX2(1,6,8)]+=tmp12_1 + tmp15_1 + tmp3_1;
                                EM_S[INDEX2(5,1,8)]+=tmp0_1 + tmp2_1 + tmp9_1;
                                EM_S[INDEX2(3,7,8)]+=tmp14_1 + tmp1_1 + tmp2_1;
                                EM_S[INDEX2(2,5,8)]+=tmp13_1 + tmp16_1 + tmp3_1;
                                EM_S[INDEX2(0,3,8)]+=tmp10_1 + tmp1_1 + tmp2_1;
                                EM_S[INDEX2(7,2,8)]+=tmp12_1 + tmp4_1 + tmp7_1;
                                EM_S[INDEX2(4,0,8)]+=tmp0_1 + tmp7_1 + tmp9_1;
                                EM_S[INDEX2(1,2,8)]+=tmp10_1 + tmp1_1 + tmp7_1;
                                EM_S[INDEX2(6,7,8)]+=tmp17_1 + tmp1_1 + tmp6_1;
                                EM_S[INDEX2(3,3,8)]+=tmp0_1 + tmp17_1 + tmp5_1;
                                EM_S[INDEX2(2,0,8)]+=tmp4_1 + tmp7_1 + tmp8_1;
                                EM_S[INDEX2(7,6,8)]+=tmp11_1 + tmp1_1 + tmp6_1;
                                EM_S[INDEX2(4,4,8)]+=tmp11_1 + tmp14_1 + tmp8_1;
                                EM_S[INDEX2(6,3,8)]+=tmp12_1 + tmp2_1 + tmp4_1;
                                EM_S[INDEX2(1,5,8)]+=tmp14_1 + tmp2_1 + tmp9_1;
                                EM_S[INDEX2(3,6,8)]+=tmp12_1 + tmp6_1 + tmp7_1;
                                EM_S[INDEX2(2,2,8)]+=tmp0_1 + tmp11_1 + tmp5_1;
                                EM_S[INDEX2(7,7,8)]+=tmp14_1 + tmp17_1 + tmp5_1;
                                EM_S[INDEX2(5,7,8)]+=tmp2_1 + tmp5_1 + tmp6_1;
                                EM_S[INDEX2(5,3,8)]+=tmp13_1 + tmp1_1 + tmp4_1;
                                EM_S[INDEX2(4,1,8)]+=tmp16_1 + tmp2_1 + tmp4_1;
                                EM_S[INDEX2(1,1,8)]+=tmp0_1 + tmp17_1 + tmp8_1;
                                EM_S[INDEX2(2,7,8)]+=tmp12_1 + tmp2_1 + tmp6_1;
                                EM_S[INDEX2(3,2,8)]+=tmp11_1 + tmp1_1 + tmp4_1;
                                EM_S[INDEX2(0,0,8)]+=tmp0_1 + tmp11_1 + tmp8_1;
                                EM_S[INDEX2(6,6,8)]+=tmp11_1 + tmp14_1 + tmp5_1;
                                EM_S[INDEX2(5,0,8)]+=tmp16_1 + tmp4_1 + tmp7_1;
                                EM_S[INDEX2(7,1,8)]+=tmp13_1 + tmp4_1 + tmp9_1;
                                EM_S[INDEX2(4,5,8)]+=tmp17_1 + tmp6_1 + tmp9_1;
                                EM_S[INDEX2(0,4,8)]+=tmp14_1 + tmp7_1 + tmp9_1;
                                EM_S[INDEX2(5,5,8)]+=tmp14_1 + tmp17_1 + tmp8_1;
                                EM_S[INDEX2(1,4,8)]+=tmp16_1 + tmp6_1 + tmp7_1;
                                EM_S[INDEX2(6,0,8)]+=tmp15_1 + tmp4_1 + tmp9_1;
                                EM_S[INDEX2(7,5,8)]+=tmp2_1 + tmp6_1 + tmp8_1;
                                EM_S[INDEX2(2,3,8)]+=tmp17_1 + tmp1_1 + tmp4_1;
                                EM_S[INDEX2(2,1,8)]+=tmp10_1 + tmp2_1 + tmp9_1;
                                EM_S[INDEX2(4,2,8)]+=tmp15_1 + tmp1_1 + tmp4_1;
                                EM_S[INDEX2(1,0,8)]+=tmp11_1 + tmp4_1 + tmp9_1;
                                EM_S[INDEX2(6,5,8)]+=tmp2_1 + tmp3_1 + tmp9_1;
                                EM_S[INDEX2(3,5,8)]+=tmp13_1 + tmp6_1 + tmp9_1;
                                EM_S[INDEX2(0,1,8)]+=tmp17_1 + tmp4_1 + tmp9_1;
                                EM_S[INDEX2(7,0,8)]+=tmp10_1 + tmp15_1 + tmp16_1;
                                EM_S[INDEX2(4,6,8)]+=tmp5_1 + tmp6_1 + tmp7_1;
                                EM_S[INDEX2(5,2,8)]+=tmp10_1 + tmp12_1 + tmp15_1;
                                EM_S[INDEX2(6,1,8)]+=tmp10_1 + tmp13_1 + tmp16_1;
                                EM_S[INDEX2(3,1,8)]+=tmp2_1 + tmp4_1 + tmp8_1;
                                EM_S[INDEX2(0,2,8)]+=tmp4_1 + tmp5_1 + tmp7_1;
                                EM_S[INDEX2(7,4,8)]+=tmp3_1 + tmp7_1 + tmp9_1;
                                EM_S[INDEX2(0,6,8)]+=tmp15_1 + tmp1_1 + tmp6_1;
                                EM_S[INDEX2(6,2,8)]+=tmp0_1 + tmp1_1 + tmp7_1;
                                EM_S[INDEX2(4,3,8)]+=tmp10_1 + tmp12_1 + tmp13_1;
                                EM_S[INDEX2(1,7,8)]+=tmp13_1 + tmp1_1 + tmp6_1;
                                EM_S[INDEX2(0,5,8)]+=tmp16_1 + tmp2_1 + tmp6_1;
                                EM_S[INDEX2(3,4,8)]+=tmp15_1 + tmp16_1 + tmp3_1;
                                EM_S[INDEX2(2,4,8)]+=tmp15_1 + tmp6_1 + tmp9_1;
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
                                const register double D_4 = D_p[4];
                                const register double D_5 = D_p[5];
                                const register double D_6 = D_p[6];
                                const register double D_7 = D_p[7];
                                const register double tmp30_0 = D_0 + D_2 + D_4 + D_6;
                                const register double tmp22_0 = D_1 + D_3 + D_4 + D_6;
                                const register double tmp6_0 = D_4 + D_6;
                                const register double tmp1_0 = D_0 + D_4;
                                const register double tmp27_0 = D_3 + D_5 + D_6;
                                const register double tmp32_0 = D_2 + D_4 + D_7;
                                const register double tmp8_0 = D_0 + D_1 + D_6 + D_7;
                                const register double tmp24_0 = D_0 + D_2;
                                const register double tmp10_0 = D_4 + D_5;
                                const register double tmp16_0 = D_0 + D_1 + D_4 + D_5;
                                const register double tmp21_0 = D_0 + D_5 + D_6;
                                const register double tmp29_0 = D_1 + D_3 + D_5 + D_7;
                                const register double tmp12_0 = D_0 + D_3 + D_4 + D_7;
                                const register double tmp28_0 = D_1 + D_2 + D_4;
                                const register double tmp11_0 = D_0 + D_1 + D_2 + D_3 + D_4 + D_5 + D_6 + D_7;
                                const register double tmp23_0 = D_5 + D_7;
                                const register double tmp0_0 = D_1 + D_2 + D_5 + D_6;
                                const register double tmp26_0 = D_1 + D_4 + D_7;
                                const register double tmp7_0 = D_1 + D_3;
                                const register double tmp4_0 = D_0 + D_1 + D_2 + D_3;
                                const register double tmp31_0 = D_0 + D_3 + D_5;
                                const register double tmp18_0 = D_0 + D_1;
                                const register double tmp3_0 = D_4 + D_5 + D_6 + D_7;
                                const register double tmp25_0 = D_0 + D_3 + D_6;
                                const register double tmp2_0 = D_3 + D_7;
                                const register double tmp19_0 = D_6 + D_7;
                                const register double tmp20_0 = D_1 + D_2 + D_7;
                                const register double tmp15_0 = D_2 + D_3 + D_6 + D_7;
                                const register double tmp5_0 = D_0 + D_2 + D_5 + D_7;
                                const register double tmp9_0 = D_2 + D_3;
                                const register double tmp17_0 = D_2 + D_3 + D_4 + D_5;
                                const register double tmp13_0 = D_1 + D_5;
                                const register double tmp14_0 = D_2 + D_6;
                                const register double tmp39_1 = tmp25_0*w148;
                                const register double tmp40_1 = D_5*w150;
                                const register double tmp4_1 = tmp4_0*w147;
                                const register double tmp42_1 = D_2*w149;
                                const register double tmp71_1 = tmp30_0*w148;
                                const register double tmp36_1 = D_3*w150;
                                const register double tmp66_1 = D_6*w149;
                                const register double tmp19_1 = tmp14_0*w144;
                                const register double tmp29_1 = D_4*w150;
                                const register double tmp73_1 = tmp19_0*w144;
                                const register double tmp28_1 = tmp20_0*w148;
                                const register double tmp24_1 = tmp1_0*w146;
                                const register double tmp14_1 = tmp10_0*w146;
                                const register double tmp21_1 = tmp15_0*w148;
                                const register double tmp57_1 = tmp10_0*w144;
                                const register double tmp10_1 = tmp4_0*w148;
                                const register double tmp58_1 = tmp9_0*w146;
                                const register double tmp69_1 = tmp25_0*w147;
                                const register double tmp13_1 = tmp9_0*w144;
                                const register double tmp26_1 = tmp18_0*w144;
                                const register double tmp52_1 = tmp15_0*w147;
                                const register double tmp72_1 = tmp29_0*w147;
                                const register double tmp18_1 = tmp14_0*w146;
                                const register double tmp25_1 = tmp17_0*w145;
                                const register double tmp55_1 = tmp32_0*w147;
                                const register double tmp65_1 = tmp31_0*w147;
                                const register double tmp56_1 = D_1*w149;
                                const register double tmp45_1 = tmp28_0*w147;
                                const register double tmp68_1 = D_2*w150;
                                const register double tmp23_1 = tmp2_0*w144;
                                const register double tmp2_1 = tmp2_0*w146;
                                const register double tmp27_1 = tmp19_0*w146;
                                const register double tmp41_1 = tmp26_0*w147;
                                const register double tmp33_1 = tmp23_0*w144;
                                const register double tmp32_1 = tmp22_0*w145;
                                const register double tmp17_1 = tmp13_0*w144;
                                const register double tmp15_1 = tmp11_0*w145;
                                const register double tmp20_1 = tmp13_0*w146;
                                const register double tmp5_1 = tmp5_0*w145;
                                const register double tmp74_1 = tmp18_0*w146;
                                const register double tmp63_1 = tmp32_0*w148;
                                const register double tmp43_1 = tmp27_0*w148;
                                const register double tmp48_1 = tmp23_0*w146;
                                const register double tmp8_1 = tmp7_0*w144;
                                const register double tmp7_1 = tmp7_0*w146;
                                const register double tmp53_1 = tmp31_0*w148;
                                const register double tmp59_1 = tmp28_0*w148;
                                const register double tmp22_1 = tmp16_0*w147;
                                const register double tmp60_1 = D_7*w150;
                                const register double tmp61_1 = tmp27_0*w147;
                                const register double tmp3_1 = tmp3_0*w148;
                                const register double tmp6_1 = tmp6_0*w144;
                                const register double tmp9_1 = tmp6_0*w146;
                                const register double tmp34_1 = tmp24_0*w146;
                                const register double tmp47_1 = tmp24_0*w144;
                                const register double tmp30_1 = tmp21_0*w147;
                                const register double tmp50_1 = tmp30_0*w147;
                                const register double tmp67_1 = tmp26_0*w148;
                                const register double tmp11_1 = tmp3_0*w147;
                                const register double tmp64_1 = D_1*w150;
                                const register double tmp46_1 = D_7*w149;
                                const register double tmp37_1 = tmp20_0*w147;
                                const register double tmp38_1 = D_4*w149;
                                const register double tmp51_1 = tmp16_0*w148;
                                const register double tmp54_1 = D_6*w150;
                                const register double tmp16_1 = tmp12_0*w145;
                                const register double tmp31_1 = D_3*w149;
                                const register double tmp12_1 = tmp8_0*w145;
                                const register double tmp62_1 = D_0*w149;
                                const register double tmp49_1 = tmp29_0*w148;
                                const register double tmp1_1 = tmp1_0*w144;
                                const register double tmp0_1 = tmp0_0*w145;
                                const register double tmp70_1 = D_5*w149;
                                const register double tmp35_1 = tmp21_0*w148;
                                const register double tmp44_1 = D_0*w150;
                                EM_S[INDEX2(0,0,8)]+=tmp59_1 + tmp60_1 + tmp61_1 + tmp62_1;
                                EM_S[INDEX2(1,0,8)]+=tmp25_1 + tmp73_1 + tmp74_1;
                                EM_S[INDEX2(2,0,8)]+=tmp32_1 + tmp33_1 + tmp34_1;
                                EM_S[INDEX2(3,0,8)]+=tmp10_1 + tmp11_1;
                                EM_S[INDEX2(4,0,8)]+=tmp0_1 + tmp23_1 + tmp24_1;
                                EM_S[INDEX2(5,0,8)]+=tmp51_1 + tmp52_1;
                                EM_S[INDEX2(6,0,8)]+=tmp71_1 + tmp72_1;
                                EM_S[INDEX2(7,0,8)]+=tmp15_1;
                                EM_S[INDEX2(0,1,8)]+=tmp25_1 + tmp73_1 + tmp74_1;
                                EM_S[INDEX2(1,1,8)]+=tmp53_1 + tmp54_1 + tmp55_1 + tmp56_1;
                                EM_S[INDEX2(2,1,8)]+=tmp10_1 + tmp11_1;
                                EM_S[INDEX2(3,1,8)]+=tmp5_1 + tmp6_1 + tmp7_1;
                                EM_S[INDEX2(4,1,8)]+=tmp51_1 + tmp52_1;
                                EM_S[INDEX2(5,1,8)]+=tmp16_1 + tmp19_1 + tmp20_1;
                                EM_S[INDEX2(6,1,8)]+=tmp15_1;
                                EM_S[INDEX2(7,1,8)]+=tmp49_1 + tmp50_1;
                                EM_S[INDEX2(0,2,8)]+=tmp32_1 + tmp33_1 + tmp34_1;
                                EM_S[INDEX2(1,2,8)]+=tmp10_1 + tmp11_1;
                                EM_S[INDEX2(2,2,8)]+=tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1;
                                EM_S[INDEX2(3,2,8)]+=tmp12_1 + tmp57_1 + tmp58_1;
                                EM_S[INDEX2(4,2,8)]+=tmp71_1 + tmp72_1;
                                EM_S[INDEX2(5,2,8)]+=tmp15_1;
                                EM_S[INDEX2(6,2,8)]+=tmp16_1 + tmp17_1 + tmp18_1;
                                EM_S[INDEX2(7,2,8)]+=tmp21_1 + tmp22_1;
                                EM_S[INDEX2(0,3,8)]+=tmp10_1 + tmp11_1;
                                EM_S[INDEX2(1,3,8)]+=tmp5_1 + tmp6_1 + tmp7_1;
                                EM_S[INDEX2(2,3,8)]+=tmp12_1 + tmp57_1 + tmp58_1;
                                EM_S[INDEX2(3,3,8)]+=tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1;
                                EM_S[INDEX2(4,3,8)]+=tmp15_1;
                                EM_S[INDEX2(5,3,8)]+=tmp49_1 + tmp50_1;
                                EM_S[INDEX2(6,3,8)]+=tmp21_1 + tmp22_1;
                                EM_S[INDEX2(7,3,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                EM_S[INDEX2(0,4,8)]+=tmp0_1 + tmp23_1 + tmp24_1;
                                EM_S[INDEX2(1,4,8)]+=tmp51_1 + tmp52_1;
                                EM_S[INDEX2(2,4,8)]+=tmp71_1 + tmp72_1;
                                EM_S[INDEX2(3,4,8)]+=tmp15_1;
                                EM_S[INDEX2(4,4,8)]+=tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1;
                                EM_S[INDEX2(5,4,8)]+=tmp12_1 + tmp13_1 + tmp14_1;
                                EM_S[INDEX2(6,4,8)]+=tmp5_1 + tmp8_1 + tmp9_1;
                                EM_S[INDEX2(7,4,8)]+=tmp3_1 + tmp4_1;
                                EM_S[INDEX2(0,5,8)]+=tmp51_1 + tmp52_1;
                                EM_S[INDEX2(1,5,8)]+=tmp16_1 + tmp19_1 + tmp20_1;
                                EM_S[INDEX2(2,5,8)]+=tmp15_1;
                                EM_S[INDEX2(3,5,8)]+=tmp49_1 + tmp50_1;
                                EM_S[INDEX2(4,5,8)]+=tmp12_1 + tmp13_1 + tmp14_1;
                                EM_S[INDEX2(5,5,8)]+=tmp67_1 + tmp68_1 + tmp69_1 + tmp70_1;
                                EM_S[INDEX2(6,5,8)]+=tmp3_1 + tmp4_1;
                                EM_S[INDEX2(7,5,8)]+=tmp32_1 + tmp47_1 + tmp48_1;
                                EM_S[INDEX2(0,6,8)]+=tmp71_1 + tmp72_1;
                                EM_S[INDEX2(1,6,8)]+=tmp15_1;
                                EM_S[INDEX2(2,6,8)]+=tmp16_1 + tmp17_1 + tmp18_1;
                                EM_S[INDEX2(3,6,8)]+=tmp21_1 + tmp22_1;
                                EM_S[INDEX2(4,6,8)]+=tmp5_1 + tmp8_1 + tmp9_1;
                                EM_S[INDEX2(5,6,8)]+=tmp3_1 + tmp4_1;
                                EM_S[INDEX2(6,6,8)]+=tmp63_1 + tmp64_1 + tmp65_1 + tmp66_1;
                                EM_S[INDEX2(7,6,8)]+=tmp25_1 + tmp26_1 + tmp27_1;
                                EM_S[INDEX2(0,7,8)]+=tmp15_1;
                                EM_S[INDEX2(1,7,8)]+=tmp49_1 + tmp50_1;
                                EM_S[INDEX2(2,7,8)]+=tmp21_1 + tmp22_1;
                                EM_S[INDEX2(3,7,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                EM_S[INDEX2(4,7,8)]+=tmp3_1 + tmp4_1;
                                EM_S[INDEX2(5,7,8)]+=tmp32_1 + tmp47_1 + tmp48_1;
                                EM_S[INDEX2(6,7,8)]+=tmp25_1 + tmp26_1 + tmp27_1;
                                EM_S[INDEX2(7,7,8)]+=tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1;
                            } else { // constant data
                                const register double D_0 = D_p[0];
                                const register double tmp3_1 = D_0*w154;
                                const register double tmp0_1 = D_0*w151;
                                const register double tmp2_1 = D_0*w153;
                                const register double tmp1_1 = D_0*w152;
                                EM_S[INDEX2(7,3,8)]+=tmp0_1;
                                EM_S[INDEX2(4,7,8)]+=tmp1_1;
                                EM_S[INDEX2(1,3,8)]+=tmp0_1;
                                EM_S[INDEX2(6,4,8)]+=tmp0_1;
                                EM_S[INDEX2(3,0,8)]+=tmp1_1;
                                EM_S[INDEX2(5,4,8)]+=tmp0_1;
                                EM_S[INDEX2(0,7,8)]+=tmp2_1;
                                EM_S[INDEX2(5,6,8)]+=tmp1_1;
                                EM_S[INDEX2(2,6,8)]+=tmp0_1;
                                EM_S[INDEX2(1,6,8)]+=tmp2_1;
                                EM_S[INDEX2(5,1,8)]+=tmp0_1;
                                EM_S[INDEX2(3,7,8)]+=tmp0_1;
                                EM_S[INDEX2(2,5,8)]+=tmp2_1;
                                EM_S[INDEX2(0,3,8)]+=tmp1_1;
                                EM_S[INDEX2(7,2,8)]+=tmp1_1;
                                EM_S[INDEX2(4,0,8)]+=tmp0_1;
                                EM_S[INDEX2(1,2,8)]+=tmp1_1;
                                EM_S[INDEX2(6,7,8)]+=tmp0_1;
                                EM_S[INDEX2(3,3,8)]+=tmp3_1;
                                EM_S[INDEX2(2,0,8)]+=tmp0_1;
                                EM_S[INDEX2(7,6,8)]+=tmp0_1;
                                EM_S[INDEX2(4,4,8)]+=tmp3_1;
                                EM_S[INDEX2(6,3,8)]+=tmp1_1;
                                EM_S[INDEX2(1,5,8)]+=tmp0_1;
                                EM_S[INDEX2(3,6,8)]+=tmp1_1;
                                EM_S[INDEX2(2,2,8)]+=tmp3_1;
                                EM_S[INDEX2(7,7,8)]+=tmp3_1;
                                EM_S[INDEX2(5,7,8)]+=tmp0_1;
                                EM_S[INDEX2(5,3,8)]+=tmp1_1;
                                EM_S[INDEX2(4,1,8)]+=tmp1_1;
                                EM_S[INDEX2(1,1,8)]+=tmp3_1;
                                EM_S[INDEX2(2,7,8)]+=tmp1_1;
                                EM_S[INDEX2(3,2,8)]+=tmp0_1;
                                EM_S[INDEX2(0,0,8)]+=tmp3_1;
                                EM_S[INDEX2(6,6,8)]+=tmp3_1;
                                EM_S[INDEX2(5,0,8)]+=tmp1_1;
                                EM_S[INDEX2(7,1,8)]+=tmp1_1;
                                EM_S[INDEX2(4,5,8)]+=tmp0_1;
                                EM_S[INDEX2(0,4,8)]+=tmp0_1;
                                EM_S[INDEX2(5,5,8)]+=tmp3_1;
                                EM_S[INDEX2(1,4,8)]+=tmp1_1;
                                EM_S[INDEX2(6,0,8)]+=tmp1_1;
                                EM_S[INDEX2(7,5,8)]+=tmp0_1;
                                EM_S[INDEX2(2,3,8)]+=tmp0_1;
                                EM_S[INDEX2(2,1,8)]+=tmp1_1;
                                EM_S[INDEX2(4,2,8)]+=tmp1_1;
                                EM_S[INDEX2(1,0,8)]+=tmp0_1;
                                EM_S[INDEX2(6,5,8)]+=tmp1_1;
                                EM_S[INDEX2(3,5,8)]+=tmp1_1;
                                EM_S[INDEX2(0,1,8)]+=tmp0_1;
                                EM_S[INDEX2(7,0,8)]+=tmp2_1;
                                EM_S[INDEX2(4,6,8)]+=tmp0_1;
                                EM_S[INDEX2(5,2,8)]+=tmp2_1;
                                EM_S[INDEX2(6,1,8)]+=tmp2_1;
                                EM_S[INDEX2(3,1,8)]+=tmp0_1;
                                EM_S[INDEX2(0,2,8)]+=tmp0_1;
                                EM_S[INDEX2(7,4,8)]+=tmp1_1;
                                EM_S[INDEX2(0,6,8)]+=tmp1_1;
                                EM_S[INDEX2(6,2,8)]+=tmp0_1;
                                EM_S[INDEX2(4,3,8)]+=tmp2_1;
                                EM_S[INDEX2(1,7,8)]+=tmp1_1;
                                EM_S[INDEX2(0,5,8)]+=tmp1_1;
                                EM_S[INDEX2(3,4,8)]+=tmp2_1;
                                EM_S[INDEX2(2,4,8)]+=tmp1_1;
                            }
                        }
                        ///////////////
                        // process X //
                        ///////////////
                        if (!X.isEmpty()) {
                            add_EM_F=true;
                            const double* X_p=const_cast<escript::Data*>(&X)->getSampleDataRO(e);
                            if (X.actsExpanded()) {
                                const register double X_0_0 = X_p[INDEX2(0,0,3)];
                                const register double X_1_0 = X_p[INDEX2(1,0,3)];
                                const register double X_2_0 = X_p[INDEX2(2,0,3)];
                                const register double X_0_1 = X_p[INDEX2(0,1,3)];
                                const register double X_1_1 = X_p[INDEX2(1,1,3)];
                                const register double X_2_1 = X_p[INDEX2(2,1,3)];
                                const register double X_0_2 = X_p[INDEX2(0,2,3)];
                                const register double X_1_2 = X_p[INDEX2(1,2,3)];
                                const register double X_2_2 = X_p[INDEX2(2,2,3)];
                                const register double X_0_3 = X_p[INDEX2(0,3,3)];
                                const register double X_1_3 = X_p[INDEX2(1,3,3)];
                                const register double X_2_3 = X_p[INDEX2(2,3,3)];
                                const register double X_0_4 = X_p[INDEX2(0,4,3)];
                                const register double X_1_4 = X_p[INDEX2(1,4,3)];
                                const register double X_2_4 = X_p[INDEX2(2,4,3)];
                                const register double X_0_5 = X_p[INDEX2(0,5,3)];
                                const register double X_1_5 = X_p[INDEX2(1,5,3)];
                                const register double X_2_5 = X_p[INDEX2(2,5,3)];
                                const register double X_0_6 = X_p[INDEX2(0,6,3)];
                                const register double X_1_6 = X_p[INDEX2(1,6,3)];
                                const register double X_2_6 = X_p[INDEX2(2,6,3)];
                                const register double X_0_7 = X_p[INDEX2(0,7,3)];
                                const register double X_1_7 = X_p[INDEX2(1,7,3)];
                                const register double X_2_7 = X_p[INDEX2(2,7,3)];
                                const register double tmp3_0 = X_1_5 + X_1_7;
                                const register double tmp17_0 = X_0_4 + X_0_5;
                                const register double tmp2_0 = X_0_2 + X_0_3 + X_0_4 + X_0_5;
                                const register double tmp14_0 = X_1_1 + X_1_3;
                                const register double tmp5_0 = X_2_1 + X_2_2 + X_2_5 + X_2_6;
                                const register double tmp1_0 = X_0_0 + X_0_1;
                                const register double tmp8_0 = X_1_0 + X_1_2;
                                const register double tmp12_0 = X_2_1 + X_2_5;
                                const register double tmp16_0 = X_0_2 + X_0_3;
                                const register double tmp9_0 = X_1_4 + X_1_6;
                                const register double tmp7_0 = X_1_1 + X_1_3 + X_1_4 + X_1_6;
                                const register double tmp13_0 = X_1_0 + X_1_2 + X_1_5 + X_1_7;
                                const register double tmp6_0 = X_2_0 + X_2_4;
                                const register double tmp4_0 = X_2_3 + X_2_7;
                                const register double tmp15_0 = X_0_0 + X_0_1 + X_0_6 + X_0_7;
                                const register double tmp11_0 = X_2_0 + X_2_3 + X_2_4 + X_2_7;
                                const register double tmp0_0 = X_0_6 + X_0_7;
                                const register double tmp10_0 = X_2_2 + X_2_6;
                                const register double tmp16_1 = tmp13_0*w158;
                                const register double tmp29_1 = tmp15_0*w165;
                                const register double tmp18_1 = tmp15_0*w160;
                                const register double tmp4_1 = tmp4_0*w161;
                                const register double tmp15_1 = tmp0_0*w166;
                                const register double tmp56_1 = tmp8_0*w169;
                                const register double tmp38_1 = tmp14_0*w162;
                                const register double tmp51_1 = tmp10_0*w170;
                                const register double tmp53_1 = tmp9_0*w167;
                                const register double tmp37_1 = tmp5_0*w171;
                                const register double tmp59_1 = tmp3_0*w167;
                                const register double tmp6_1 = tmp6_0*w157;
                                const register double tmp13_1 = tmp11_0*w159;
                                const register double tmp14_1 = tmp12_0*w157;
                                const register double tmp31_1 = tmp9_0*w169;
                                const register double tmp28_1 = tmp13_0*w168;
                                const register double tmp10_1 = tmp9_0*w162;
                                const register double tmp3_1 = tmp3_0*w162;
                                const register double tmp55_1 = tmp0_0*w164;
                                const register double tmp54_1 = tmp6_0*w172;
                                const register double tmp50_1 = tmp14_0*w169;
                                const register double tmp20_1 = tmp12_0*w161;
                                const register double tmp58_1 = tmp1_0*w166;
                                const register double tmp24_1 = tmp17_0*w163;
                                const register double tmp49_1 = tmp12_0*w172;
                                const register double tmp39_1 = tmp6_0*w170;
                                const register double tmp8_1 = tmp8_0*w156;
                                const register double tmp46_1 = tmp16_0*w166;
                                const register double tmp2_1 = tmp2_0*w160;
                                const register double tmp25_1 = tmp8_0*w167;
                                const register double tmp33_1 = tmp14_0*w167;
                                const register double tmp12_1 = tmp2_0*w165;
                                const register double tmp21_1 = tmp7_0*w168;
                                const register double tmp19_1 = tmp16_0*w155;
                                const register double tmp26_1 = tmp16_0*w164;
                                const register double tmp35_1 = tmp17_0*w155;
                                const register double tmp36_1 = tmp4_0*w172;
                                const register double tmp9_1 = tmp1_0*w164;
                                const register double tmp34_1 = tmp9_0*w156;
                                const register double tmp40_1 = tmp16_0*w163;
                                const register double tmp43_1 = tmp11_0*w171;
                                const register double tmp57_1 = tmp4_0*w170;
                                const register double tmp52_1 = tmp1_0*w163;
                                const register double tmp41_1 = tmp10_0*w172;
                                const register double tmp1_1 = tmp1_0*w155;
                                const register double tmp11_1 = tmp10_0*w161;
                                const register double tmp22_1 = tmp10_0*w157;
                                const register double tmp23_1 = tmp3_0*w169;
                                const register double tmp5_1 = tmp5_0*w159;
                                const register double tmp27_1 = tmp6_0*w161;
                                const register double tmp32_1 = tmp17_0*w166;
                                const register double tmp42_1 = tmp17_0*w164;
                                const register double tmp30_1 = tmp4_0*w157;
                                const register double tmp47_1 = tmp3_0*w156;
                                const register double tmp48_1 = tmp0_0*w155;
                                const register double tmp7_1 = tmp7_0*w158;
                                const register double tmp45_1 = tmp12_0*w170;
                                const register double tmp44_1 = tmp8_0*w162;
                                const register double tmp0_1 = tmp0_0*w163;
                                const register double tmp17_1 = tmp14_0*w156;
                                EM_F[0]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1;
                                EM_F[1]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp16_1 + tmp17_1 + tmp9_1;
                                EM_F[2]+=tmp13_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp24_1 + tmp25_1;
                                EM_F[3]+=tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp33_1 + tmp5_1;
                                EM_F[4]+=tmp16_1 + tmp18_1 + tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1 + tmp40_1;
                                EM_F[5]+=tmp29_1 + tmp41_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp7_1;
                                EM_F[6]+=tmp28_1 + tmp2_1 + tmp43_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1 + tmp52_1 + tmp53_1;
                                EM_F[7]+=tmp12_1 + tmp21_1 + tmp37_1 + tmp54_1 + tmp55_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1;
                            } else { // constant data
                                const register double X_0 = X_p[0];
                                const register double X_1 = X_p[1];
                                const register double X_2 = X_p[2];
                                const register double tmp4_1 = X_1*w177;
                                const register double tmp3_1 = X_0*w176;
                                const register double tmp5_1 = X_2*w178;
                                const register double tmp0_1 = X_0*w173;
                                const register double tmp2_1 = X_2*w175;
                                const register double tmp1_1 = X_1*w174;
                                EM_F[0]+=tmp0_1 + tmp1_1 + tmp2_1;
                                EM_F[1]+=tmp1_1 + tmp2_1 + tmp3_1;
                                EM_F[2]+=tmp0_1 + tmp2_1 + tmp4_1;
                                EM_F[3]+=tmp2_1 + tmp3_1 + tmp4_1;
                                EM_F[4]+=tmp0_1 + tmp1_1 + tmp5_1;
                                EM_F[5]+=tmp1_1 + tmp3_1 + tmp5_1;
                                EM_F[6]+=tmp0_1 + tmp4_1 + tmp5_1;
                                EM_F[7]+=tmp3_1 + tmp4_1 + tmp5_1;
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
                                const register double Y_4 = Y_p[4];
                                const register double Y_5 = Y_p[5];
                                const register double Y_6 = Y_p[6];
                                const register double Y_7 = Y_p[7];
                                const register double tmp3_0 = Y_0 + Y_3 + Y_5;
                                const register double tmp5_0 = Y_0 + Y_3 + Y_6;
                                const register double tmp6_0 = Y_0 + Y_5 + Y_6;
                                const register double tmp2_0 = Y_2 + Y_4 + Y_7;
                                const register double tmp4_0 = Y_1 + Y_4 + Y_7;
                                const register double tmp0_0 = Y_3 + Y_5 + Y_6;
                                const register double tmp1_0 = Y_1 + Y_2 + Y_4;
                                const register double tmp7_0 = Y_1 + Y_2 + Y_7;
                                const register double tmp30_1 = tmp0_0*w180;
                                const register double tmp29_1 = tmp1_0*w181;
                                const register double tmp22_1 = tmp4_0*w180;
                                const register double tmp7_1 = Y_6*w182;
                                const register double tmp31_1 = Y_0*w182;
                                const register double tmp28_1 = Y_7*w179;
                                const register double tmp19_1 = Y_3*w182;
                                const register double tmp12_1 = Y_3*w179;
                                const register double tmp3_1 = Y_7*w182;
                                const register double tmp26_1 = tmp2_0*w180;
                                const register double tmp13_1 = tmp6_0*w181;
                                const register double tmp2_1 = tmp1_0*w180;
                                const register double tmp0_1 = Y_0*w179;
                                const register double tmp17_1 = tmp7_0*w181;
                                const register double tmp24_1 = Y_6*w179;
                                const register double tmp10_1 = tmp5_0*w180;
                                const register double tmp16_1 = Y_4*w179;
                                const register double tmp8_1 = Y_2*w179;
                                const register double tmp18_1 = tmp6_0*w180;
                                const register double tmp5_1 = tmp2_0*w181;
                                const register double tmp6_1 = tmp3_0*w180;
                                const register double tmp27_1 = Y_1*w182;
                                const register double tmp20_1 = Y_5*w179;
                                const register double tmp9_1 = tmp4_0*w181;
                                const register double tmp25_1 = tmp3_0*w181;
                                const register double tmp11_1 = Y_5*w182;
                                const register double tmp4_1 = Y_1*w179;
                                const register double tmp1_1 = tmp0_0*w181;
                                const register double tmp21_1 = tmp5_0*w181;
                                const register double tmp14_1 = tmp7_0*w180;
                                const register double tmp23_1 = Y_2*w182;
                                const register double tmp15_1 = Y_4*w182;
                                EM_F[0]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                                EM_F[1]+=tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1;
                                EM_F[2]+=tmp10_1 + tmp11_1 + tmp8_1 + tmp9_1;
                                EM_F[3]+=tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1;
                                EM_F[4]+=tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1;
                                EM_F[5]+=tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1;
                                EM_F[6]+=tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1;
                                EM_F[7]+=tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1;
                            } else { // constant data
                                const register double Y_0 = Y_p[0];
                                const register double tmp0_1 = Y_0*w183;
                                EM_F[0]+=tmp0_1;
                                EM_F[1]+=tmp0_1;
                                EM_F[2]+=tmp0_1;
                                EM_F[3]+=tmp0_1;
                                EM_F[4]+=tmp0_1;
                                EM_F[5]+=tmp0_1;
                                EM_F[6]+=tmp0_1;
                                EM_F[7]+=tmp0_1;
                            }
                        }

                        // add to matrix (if add_EM_S) and RHS (if add_EM_F)
                        const index_t firstNode=m_N0*m_N1*k2+m_N0*k1+k0;
                        IndexVector rowIndex;
                        rowIndex.push_back(m_dofMap[firstNode]);
                        rowIndex.push_back(m_dofMap[firstNode+1]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0+1]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0*m_N1]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0*m_N1+1]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0*(m_N1+1)]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0*(m_N1+1)+1]);
                        if (add_EM_F) {
                            double *F_p=rhs.getSampleDataRW(0);
                            for (index_t i=0; i<rowIndex.size(); i++) {
                                if (rowIndex[i]<getNumDOF()) {
                                    F_p[rowIndex[i]]+=EM_F[i];
                                }
                            }
                        }
                        if (add_EM_S) {
                            addToSystemMatrix(mat, rowIndex, 1, rowIndex, 1, EM_S);
                        }
                    } // end k0 loop
                } // end k1 loop
            } // end k2 loop
        } // end of colouring
    } // end of parallel region
}

//protected
void Brick::assemblePDESingleReduced(Paso_SystemMatrix* mat,
        escript::Data& rhs, const escript::Data& A, const escript::Data& B,
        const escript::Data& C, const escript::Data& D,
        const escript::Data& X, const escript::Data& Y,
        const escript::Data& d, const escript::Data& y) const
{
    const double h0 = m_l0/m_gNE0;
    const double h1 = m_l1/m_gNE1;
    const double h2 = m_l2/m_gNE2;
    const double w0 = 0.0625*h1*h2/h0;
    const double w1 = 0.0625*h2;
    const double w2 = -0.0625*h1;
    const double w3 = 0.0625*h0*h2/h1;
    const double w4 = -0.0625*h0;
    const double w5 = 0.0625*h1;
    const double w6 = 0.0625*h0;
    const double w7 = -0.0625*h0*h1/h2;
    const double w8 = -0.0625*h1*h2/h0;
    const double w9 = -0.0625*h2;
    const double w10 = -0.0625*h0*h2/h1;
    const double w11 = 0.0625*h0*h1/h2;
    const double w12 = 0.03125*h1*h2;
    const double w13 = 0.03125*h0*h2;
    const double w14 = 0.03125*h0*h1;
    const double w15 = -0.03125*h1*h2;
    const double w16 = -0.03125*h0*h2;
    const double w17 = -0.03125*h0*h1;
    const double w18 = 0.015625*h0*h1*h2;
    const double w19 = -0.25*h1*h2;
    const double w20 = -0.25*h0*h2;
    const double w21 = -0.25*h0*h1;
    const double w22 = 0.25*h1*h2;
    const double w23 = 0.25*h0*h2;
    const double w24 = 0.25*h0*h1;
    const double w25 = 0.125*h0*h1*h2;

    rhs.requireWrite();
#pragma omp parallel
    {
        for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
#pragma omp for
            for (index_t k2=k2_0; k2<m_NE2; k2+=2) {
                for (index_t k1=0; k1<m_NE1; ++k1) {
                    for (index_t k0=0; k0<m_NE0; ++k0)  {
                        bool add_EM_S=false;
                        bool add_EM_F=false;
                        vector<double> EM_S(8*8, 0);
                        vector<double> EM_F(8, 0);
                        const index_t e = k0 + m_NE0*k1 + m_NE0*m_NE1*k2;
                        ///////////////
                        // process A //
                        ///////////////
                        if (!A.isEmpty()) {
                            add_EM_S=true;
                            const double* A_p=const_cast<escript::Data*>(&A)->getSampleDataRO(e);
                            const register double A_00 = A_p[INDEX2(0,0,3)];
                            const register double A_01 = A_p[INDEX2(0,1,3)];
                            const register double A_02 = A_p[INDEX2(0,2,3)];
                            const register double A_10 = A_p[INDEX2(1,0,3)];
                            const register double A_11 = A_p[INDEX2(1,1,3)];
                            const register double A_12 = A_p[INDEX2(1,2,3)];
                            const register double A_20 = A_p[INDEX2(2,0,3)];
                            const register double A_21 = A_p[INDEX2(2,1,3)];
                            const register double A_22 = A_p[INDEX2(2,2,3)];
                            const register double tmp0_0 = A_01 + A_10;
                            const register double tmp1_0 = A_02 + A_20;
                            const register double tmp2_0 = A_12 + A_21;
                            const register double tmp3_1 = A_22*w7;
                            const register double tmp10_1 = A_11*w10;
                            const register double tmp21_1 = A_02*w5;
                            const register double tmp2_1 = A_00*w0;
                            const register double tmp23_1 = tmp2_0*w6;
                            const register double tmp19_1 = A_20*w2;
                            const register double tmp4_1 = A_11*w3;
                            const register double tmp22_1 = tmp1_0*w5;
                            const register double tmp13_1 = A_21*w4;
                            const register double tmp5_1 = A_21*w6;
                            const register double tmp8_1 = A_00*w8;
                            const register double tmp7_1 = A_20*w5;
                            const register double tmp18_1 = tmp2_0*w4;
                            const register double tmp6_1 = A_02*w2;
                            const register double tmp9_1 = A_22*w11;
                            const register double tmp15_1 = tmp1_0*w2;
                            const register double tmp12_1 = A_01*w1;
                            const register double tmp0_1 = tmp0_0*w1;
                            const register double tmp20_1 = A_01*w9;
                            const register double tmp14_1 = A_12*w6;
                            const register double tmp1_1 = A_12*w4;
                            const register double tmp16_1 = A_10*w9;
                            const register double tmp11_1 = tmp0_0*w9;
                            const register double tmp17_1 = A_10*w1;
                            EM_S[INDEX2(7,3,8)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1;
                            EM_S[INDEX2(4,7,8)]+=tmp10_1 + tmp11_1 + tmp1_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(1,3,8)]+=tmp10_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp16_1 + tmp2_1 + tmp9_1;
                            EM_S[INDEX2(6,4,8)]+=tmp10_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp16_1 + tmp2_1 + tmp9_1;
                            EM_S[INDEX2(3,0,8)]+=tmp10_1 + tmp11_1 + tmp1_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(5,4,8)]+=tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp4_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(0,7,8)]+=tmp10_1 + tmp11_1 + tmp15_1 + tmp18_1 + tmp3_1 + tmp8_1;
                            EM_S[INDEX2(5,6,8)]+=tmp0_1 + tmp10_1 + tmp19_1 + tmp1_1 + tmp21_1 + tmp5_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(2,6,8)]+=tmp11_1 + tmp13_1 + tmp14_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp6_1 + tmp7_1;
                            EM_S[INDEX2(1,6,8)]+=tmp0_1 + tmp10_1 + tmp18_1 + tmp22_1 + tmp3_1 + tmp8_1;
                            EM_S[INDEX2(5,1,8)]+=tmp11_1 + tmp13_1 + tmp14_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp6_1 + tmp7_1;
                            EM_S[INDEX2(3,7,8)]+=tmp0_1 + tmp13_1 + tmp14_1 + tmp19_1 + tmp21_1 + tmp2_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(2,5,8)]+=tmp0_1 + tmp10_1 + tmp15_1 + tmp23_1 + tmp3_1 + tmp8_1;
                            EM_S[INDEX2(0,3,8)]+=tmp10_1 + tmp11_1 + tmp13_1 + tmp14_1 + tmp19_1 + tmp21_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(7,2,8)]+=tmp12_1 + tmp15_1 + tmp16_1 + tmp1_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp8_1;
                            EM_S[INDEX2(4,0,8)]+=tmp0_1 + tmp13_1 + tmp14_1 + tmp19_1 + tmp21_1 + tmp2_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(1,2,8)]+=tmp0_1 + tmp10_1 + tmp13_1 + tmp14_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(6,7,8)]+=tmp17_1 + tmp20_1 + tmp23_1 + tmp4_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(3,3,8)]+=tmp0_1 + tmp15_1 + tmp18_1 + tmp2_1 + tmp4_1 + tmp9_1;
                            EM_S[INDEX2(2,0,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp1_1 + tmp22_1 + tmp2_1 + tmp5_1 + tmp9_1;
                            EM_S[INDEX2(7,6,8)]+=tmp12_1 + tmp16_1 + tmp19_1 + tmp21_1 + tmp23_1 + tmp4_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(4,4,8)]+=tmp0_1 + tmp15_1 + tmp18_1 + tmp2_1 + tmp4_1 + tmp9_1;
                            EM_S[INDEX2(6,3,8)]+=tmp17_1 + tmp1_1 + tmp20_1 + tmp22_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp8_1;
                            EM_S[INDEX2(1,5,8)]+=tmp11_1 + tmp19_1 + tmp1_1 + tmp21_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(3,6,8)]+=tmp12_1 + tmp13_1 + tmp14_1 + tmp16_1 + tmp22_1 + tmp3_1 + tmp4_1 + tmp8_1;
                            EM_S[INDEX2(2,2,8)]+=tmp11_1 + tmp18_1 + tmp22_1 + tmp2_1 + tmp4_1 + tmp9_1;
                            EM_S[INDEX2(7,7,8)]+=tmp0_1 + tmp22_1 + tmp23_1 + tmp2_1 + tmp4_1 + tmp9_1;
                            EM_S[INDEX2(5,7,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp1_1 + tmp22_1 + tmp2_1 + tmp5_1 + tmp9_1;
                            EM_S[INDEX2(5,3,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp23_1 + tmp2_1 + tmp3_1 + tmp6_1 + tmp7_1;
                            EM_S[INDEX2(4,1,8)]+=tmp12_1 + tmp13_1 + tmp14_1 + tmp16_1 + tmp22_1 + tmp3_1 + tmp4_1 + tmp8_1;
                            EM_S[INDEX2(1,1,8)]+=tmp11_1 + tmp15_1 + tmp23_1 + tmp2_1 + tmp4_1 + tmp9_1;
                            EM_S[INDEX2(2,7,8)]+=tmp13_1 + tmp14_1 + tmp15_1 + tmp17_1 + tmp20_1 + tmp3_1 + tmp4_1 + tmp8_1;
                            EM_S[INDEX2(3,2,8)]+=tmp12_1 + tmp16_1 + tmp18_1 + tmp4_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(0,0,8)]+=tmp0_1 + tmp22_1 + tmp23_1 + tmp2_1 + tmp4_1 + tmp9_1;
                            EM_S[INDEX2(6,6,8)]+=tmp11_1 + tmp15_1 + tmp23_1 + tmp2_1 + tmp4_1 + tmp9_1;
                            EM_S[INDEX2(5,0,8)]+=tmp13_1 + tmp14_1 + tmp15_1 + tmp17_1 + tmp20_1 + tmp3_1 + tmp4_1 + tmp8_1;
                            EM_S[INDEX2(7,1,8)]+=tmp10_1 + tmp17_1 + tmp18_1 + tmp20_1 + tmp2_1 + tmp3_1 + tmp6_1 + tmp7_1;
                            EM_S[INDEX2(4,5,8)]+=tmp12_1 + tmp16_1 + tmp18_1 + tmp4_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(0,4,8)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1;
                            EM_S[INDEX2(5,5,8)]+=tmp11_1 + tmp18_1 + tmp22_1 + tmp2_1 + tmp4_1 + tmp9_1;
                            EM_S[INDEX2(1,4,8)]+=tmp17_1 + tmp1_1 + tmp20_1 + tmp22_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp8_1;
                            EM_S[INDEX2(6,0,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp18_1 + tmp19_1 + tmp21_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(7,5,8)]+=tmp10_1 + tmp13_1 + tmp14_1 + tmp17_1 + tmp20_1 + tmp22_1 + tmp2_1 + tmp9_1;
                            EM_S[INDEX2(2,3,8)]+=tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp4_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(2,1,8)]+=tmp0_1 + tmp10_1 + tmp19_1 + tmp1_1 + tmp21_1 + tmp5_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(4,2,8)]+=tmp10_1 + tmp17_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp23_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(1,0,8)]+=tmp17_1 + tmp20_1 + tmp23_1 + tmp4_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(6,5,8)]+=tmp0_1 + tmp10_1 + tmp13_1 + tmp14_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(3,5,8)]+=tmp10_1 + tmp17_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp23_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(0,1,8)]+=tmp12_1 + tmp16_1 + tmp19_1 + tmp21_1 + tmp23_1 + tmp4_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(7,0,8)]+=tmp10_1 + tmp11_1 + tmp15_1 + tmp18_1 + tmp3_1 + tmp8_1;
                            EM_S[INDEX2(4,6,8)]+=tmp10_1 + tmp15_1 + tmp17_1 + tmp1_1 + tmp20_1 + tmp2_1 + tmp5_1 + tmp9_1;
                            EM_S[INDEX2(5,2,8)]+=tmp0_1 + tmp10_1 + tmp15_1 + tmp23_1 + tmp3_1 + tmp8_1;
                            EM_S[INDEX2(6,1,8)]+=tmp0_1 + tmp10_1 + tmp18_1 + tmp22_1 + tmp3_1 + tmp8_1;
                            EM_S[INDEX2(3,1,8)]+=tmp10_1 + tmp15_1 + tmp17_1 + tmp1_1 + tmp20_1 + tmp2_1 + tmp5_1 + tmp9_1;
                            EM_S[INDEX2(0,2,8)]+=tmp10_1 + tmp13_1 + tmp14_1 + tmp17_1 + tmp20_1 + tmp22_1 + tmp2_1 + tmp9_1;
                            EM_S[INDEX2(7,4,8)]+=tmp10_1 + tmp11_1 + tmp13_1 + tmp14_1 + tmp19_1 + tmp21_1 + tmp8_1 + tmp9_1;
                            EM_S[INDEX2(0,6,8)]+=tmp10_1 + tmp17_1 + tmp18_1 + tmp20_1 + tmp2_1 + tmp3_1 + tmp6_1 + tmp7_1;
                            EM_S[INDEX2(6,2,8)]+=tmp11_1 + tmp19_1 + tmp1_1 + tmp21_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(4,3,8)]+=tmp10_1 + tmp11_1 + tmp22_1 + tmp23_1 + tmp3_1 + tmp8_1;
                            EM_S[INDEX2(1,7,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp18_1 + tmp19_1 + tmp21_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(0,5,8)]+=tmp12_1 + tmp15_1 + tmp16_1 + tmp1_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp8_1;
                            EM_S[INDEX2(3,4,8)]+=tmp10_1 + tmp11_1 + tmp22_1 + tmp23_1 + tmp3_1 + tmp8_1;
                            EM_S[INDEX2(2,4,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp23_1 + tmp2_1 + tmp3_1 + tmp6_1 + tmp7_1;
                        }
                        ///////////////
                        // process B //
                        ///////////////
                        if (!B.isEmpty()) {
                            add_EM_S=true;
                            const double* B_p=const_cast<escript::Data*>(&B)->getSampleDataRO(e);
                            const register double B_0 = B_p[0];
                            const register double B_1 = B_p[1];
                            const register double B_2 = B_p[2];
                            const register double tmp4_1 = B_0*w15;
                            const register double tmp3_1 = B_1*w16;
                            const register double tmp2_1 = B_0*w12;
                            const register double tmp5_1 = B_2*w17;
                            const register double tmp1_1 = B_2*w14;
                            const register double tmp0_1 = B_1*w13;
                            EM_S[INDEX2(7,3,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(4,7,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(1,3,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                            EM_S[INDEX2(6,4,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                            EM_S[INDEX2(3,0,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                            EM_S[INDEX2(5,4,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(0,7,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(5,6,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(2,6,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(1,6,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                            EM_S[INDEX2(5,1,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(3,7,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                            EM_S[INDEX2(2,5,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(0,3,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(7,2,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(4,0,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(1,2,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                            EM_S[INDEX2(6,7,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                            EM_S[INDEX2(3,3,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                            EM_S[INDEX2(2,0,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(7,6,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(4,4,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(6,3,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                            EM_S[INDEX2(1,5,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                            EM_S[INDEX2(3,6,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                            EM_S[INDEX2(2,2,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(7,7,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(5,7,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(5,3,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(4,1,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(1,1,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                            EM_S[INDEX2(2,7,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(3,2,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                            EM_S[INDEX2(0,0,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(6,6,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                            EM_S[INDEX2(5,0,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(7,1,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(4,5,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(0,4,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(5,5,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(1,4,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                            EM_S[INDEX2(6,0,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                            EM_S[INDEX2(7,5,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(2,3,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(2,1,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(4,2,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(1,0,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                            EM_S[INDEX2(6,5,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                            EM_S[INDEX2(3,5,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                            EM_S[INDEX2(0,1,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(7,0,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(4,6,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(5,2,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(6,1,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                            EM_S[INDEX2(3,1,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                            EM_S[INDEX2(0,2,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(7,4,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(0,6,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(6,2,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                            EM_S[INDEX2(4,3,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(1,7,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                            EM_S[INDEX2(0,5,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(3,4,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                            EM_S[INDEX2(2,4,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                        }
                        ///////////////
                        // process C //
                        ///////////////
                        if (!C.isEmpty()) {
                            add_EM_S=true;
                            const double* C_p=const_cast<escript::Data*>(&C)->getSampleDataRO(e);
                            const register double C_0 = C_p[0];
                            const register double C_1 = C_p[1];
                            const register double C_2 = C_p[2];
                            const register double tmp5_1 = C_0*w15;
                            const register double tmp2_1 = C_0*w12;
                            const register double tmp4_1 = C_1*w16;
                            const register double tmp1_1 = C_2*w17;
                            const register double tmp3_1 = C_2*w14;
                            const register double tmp0_1 = C_1*w13;
                            EM_S[INDEX2(7,3,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(4,7,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(1,3,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(6,4,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(3,0,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(5,4,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(0,7,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(5,6,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                            EM_S[INDEX2(2,6,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                            EM_S[INDEX2(1,6,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                            EM_S[INDEX2(5,1,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                            EM_S[INDEX2(3,7,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(2,5,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(0,3,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(7,2,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                            EM_S[INDEX2(4,0,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(1,2,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                            EM_S[INDEX2(6,7,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(3,3,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(2,0,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(7,6,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                            EM_S[INDEX2(4,4,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(6,3,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(1,5,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(3,6,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                            EM_S[INDEX2(2,2,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                            EM_S[INDEX2(7,7,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(5,7,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(5,3,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(4,1,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                            EM_S[INDEX2(1,1,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                            EM_S[INDEX2(2,7,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(3,2,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                            EM_S[INDEX2(0,0,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(6,6,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                            EM_S[INDEX2(5,0,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(7,1,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                            EM_S[INDEX2(4,5,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(0,4,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(5,5,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(1,4,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(6,0,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(7,5,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(2,3,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(2,1,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                            EM_S[INDEX2(4,2,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                            EM_S[INDEX2(1,0,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(6,5,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(3,5,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(0,1,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                            EM_S[INDEX2(7,0,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(4,6,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                            EM_S[INDEX2(5,2,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                            EM_S[INDEX2(6,1,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                            EM_S[INDEX2(3,1,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                            EM_S[INDEX2(0,2,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                            EM_S[INDEX2(7,4,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(0,6,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                            EM_S[INDEX2(6,2,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                            EM_S[INDEX2(4,3,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_S[INDEX2(1,7,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                            EM_S[INDEX2(0,5,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                            EM_S[INDEX2(3,4,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                            EM_S[INDEX2(2,4,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                        }
                        ///////////////
                        // process D //
                        ///////////////
                        if (!D.isEmpty()) {
                            add_EM_S=true;
                            const double* D_p=const_cast<escript::Data*>(&D)->getSampleDataRO(e);
                            const register double D_0 = D_p[0];
                            const register double tmp0_1 = D_0*w18;
                            EM_S[INDEX2(7,3,8)]+=tmp0_1;
                            EM_S[INDEX2(4,7,8)]+=tmp0_1;
                            EM_S[INDEX2(1,3,8)]+=tmp0_1;
                            EM_S[INDEX2(6,4,8)]+=tmp0_1;
                            EM_S[INDEX2(3,0,8)]+=tmp0_1;
                            EM_S[INDEX2(5,4,8)]+=tmp0_1;
                            EM_S[INDEX2(0,7,8)]+=tmp0_1;
                            EM_S[INDEX2(5,6,8)]+=tmp0_1;
                            EM_S[INDEX2(2,6,8)]+=tmp0_1;
                            EM_S[INDEX2(1,6,8)]+=tmp0_1;
                            EM_S[INDEX2(5,1,8)]+=tmp0_1;
                            EM_S[INDEX2(3,7,8)]+=tmp0_1;
                            EM_S[INDEX2(2,5,8)]+=tmp0_1;
                            EM_S[INDEX2(0,3,8)]+=tmp0_1;
                            EM_S[INDEX2(7,2,8)]+=tmp0_1;
                            EM_S[INDEX2(4,0,8)]+=tmp0_1;
                            EM_S[INDEX2(1,2,8)]+=tmp0_1;
                            EM_S[INDEX2(6,7,8)]+=tmp0_1;
                            EM_S[INDEX2(3,3,8)]+=tmp0_1;
                            EM_S[INDEX2(2,0,8)]+=tmp0_1;
                            EM_S[INDEX2(7,6,8)]+=tmp0_1;
                            EM_S[INDEX2(4,4,8)]+=tmp0_1;
                            EM_S[INDEX2(6,3,8)]+=tmp0_1;
                            EM_S[INDEX2(1,5,8)]+=tmp0_1;
                            EM_S[INDEX2(3,6,8)]+=tmp0_1;
                            EM_S[INDEX2(2,2,8)]+=tmp0_1;
                            EM_S[INDEX2(7,7,8)]+=tmp0_1;
                            EM_S[INDEX2(5,7,8)]+=tmp0_1;
                            EM_S[INDEX2(5,3,8)]+=tmp0_1;
                            EM_S[INDEX2(4,1,8)]+=tmp0_1;
                            EM_S[INDEX2(1,1,8)]+=tmp0_1;
                            EM_S[INDEX2(2,7,8)]+=tmp0_1;
                            EM_S[INDEX2(3,2,8)]+=tmp0_1;
                            EM_S[INDEX2(0,0,8)]+=tmp0_1;
                            EM_S[INDEX2(6,6,8)]+=tmp0_1;
                            EM_S[INDEX2(5,0,8)]+=tmp0_1;
                            EM_S[INDEX2(7,1,8)]+=tmp0_1;
                            EM_S[INDEX2(4,5,8)]+=tmp0_1;
                            EM_S[INDEX2(0,4,8)]+=tmp0_1;
                            EM_S[INDEX2(5,5,8)]+=tmp0_1;
                            EM_S[INDEX2(1,4,8)]+=tmp0_1;
                            EM_S[INDEX2(6,0,8)]+=tmp0_1;
                            EM_S[INDEX2(7,5,8)]+=tmp0_1;
                            EM_S[INDEX2(2,3,8)]+=tmp0_1;
                            EM_S[INDEX2(2,1,8)]+=tmp0_1;
                            EM_S[INDEX2(4,2,8)]+=tmp0_1;
                            EM_S[INDEX2(1,0,8)]+=tmp0_1;
                            EM_S[INDEX2(6,5,8)]+=tmp0_1;
                            EM_S[INDEX2(3,5,8)]+=tmp0_1;
                            EM_S[INDEX2(0,1,8)]+=tmp0_1;
                            EM_S[INDEX2(7,0,8)]+=tmp0_1;
                            EM_S[INDEX2(4,6,8)]+=tmp0_1;
                            EM_S[INDEX2(5,2,8)]+=tmp0_1;
                            EM_S[INDEX2(6,1,8)]+=tmp0_1;
                            EM_S[INDEX2(3,1,8)]+=tmp0_1;
                            EM_S[INDEX2(0,2,8)]+=tmp0_1;
                            EM_S[INDEX2(7,4,8)]+=tmp0_1;
                            EM_S[INDEX2(0,6,8)]+=tmp0_1;
                            EM_S[INDEX2(6,2,8)]+=tmp0_1;
                            EM_S[INDEX2(4,3,8)]+=tmp0_1;
                            EM_S[INDEX2(1,7,8)]+=tmp0_1;
                            EM_S[INDEX2(0,5,8)]+=tmp0_1;
                            EM_S[INDEX2(3,4,8)]+=tmp0_1;
                            EM_S[INDEX2(2,4,8)]+=tmp0_1;
                        }
                        ///////////////
                        // process X //
                        ///////////////
                        if (!X.isEmpty()) {
                            add_EM_F=true;
                            const double* X_p=const_cast<escript::Data*>(&X)->getSampleDataRO(e);
                            const register double X_0 = X_p[0];
                            const register double X_1 = X_p[1];
                            const register double X_2 = X_p[2];
                            const register double tmp0_1 = X_2*w21;
                            const register double tmp1_1 = X_0*w19;
                            const register double tmp2_1 = X_1*w20;
                            const register double tmp3_1 = X_0*w22;
                            const register double tmp4_1 = X_1*w23;
                            const register double tmp5_1 = X_2*w24;
                            EM_F[0]+=tmp0_1 + tmp1_1 + tmp2_1;
                            EM_F[1]+=tmp0_1 + tmp2_1 + tmp3_1;
                            EM_F[2]+=tmp0_1 + tmp1_1 + tmp4_1;
                            EM_F[3]+=tmp0_1 + tmp3_1 + tmp4_1;
                            EM_F[4]+=tmp1_1 + tmp2_1 + tmp5_1;
                            EM_F[5]+=tmp2_1 + tmp3_1 + tmp5_1;
                            EM_F[6]+=tmp1_1 + tmp4_1 + tmp5_1;
                            EM_F[7]+=tmp3_1 + tmp4_1 + tmp5_1;
                        }
                        ///////////////
                        // process Y //
                        ///////////////
                        if (!Y.isEmpty()) {
                            add_EM_F=true;
                            const double* Y_p=const_cast<escript::Data*>(&Y)->getSampleDataRO(e);
                            const register double tmp0_1 = Y_p[0]*w25;
                            EM_F[0]+=tmp0_1;
                            EM_F[1]+=tmp0_1;
                            EM_F[2]+=tmp0_1;
                            EM_F[3]+=tmp0_1;
                            EM_F[4]+=tmp0_1;
                            EM_F[5]+=tmp0_1;
                            EM_F[6]+=tmp0_1;
                            EM_F[7]+=tmp0_1;
                        }

                        // add to matrix (if add_EM_S) and RHS (if add_EM_F)
                        const index_t firstNode=m_N0*m_N1*k2+m_N0*k1+k0;
                        IndexVector rowIndex;
                        rowIndex.push_back(m_dofMap[firstNode]);
                        rowIndex.push_back(m_dofMap[firstNode+1]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0+1]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0*m_N1]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0*m_N1+1]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0*(m_N1+1)]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0*(m_N1+1)+1]);
                        if (add_EM_F) {
                            double *F_p=rhs.getSampleDataRW(0);
                            for (index_t i=0; i<rowIndex.size(); i++) {
                                if (rowIndex[i]<getNumDOF()) {
                                    F_p[rowIndex[i]]+=EM_F[i];
                                }
                            }
                        }
                        if (add_EM_S) {
                            addToSystemMatrix(mat, rowIndex, 1, rowIndex, 1, EM_S);
                        }
                    } // end k0 loop
                } // end k1 loop
            } // end k2 loop
        } // end of colouring
    } // end of parallel region
}

//protected
void Brick::assemblePDESystem(Paso_SystemMatrix* mat, escript::Data& rhs,
        const escript::Data& A, const escript::Data& B,
        const escript::Data& C, const escript::Data& D,
        const escript::Data& X, const escript::Data& Y,
        const escript::Data& d, const escript::Data& y) const
{
    const double h0 = m_l0/m_gNE0;
    const double h1 = m_l1/m_gNE1;
    const double h2 = m_l2/m_gNE2;
    dim_t numEq, numComp;
    if (!mat)
        numEq=numComp=(rhs.isEmpty() ? 1 : rhs.getDataPointSize());
    else {
        numEq=mat->logical_row_block_size;
        numComp=mat->logical_col_block_size;
    }

    const double w0 = 0.0009303791403858427308*h1*h2/h0;
    const double w1 = 0.0009303791403858427308*h2;
    const double w2 = -0.00024929433932114870101*h1;
    const double w3 = 0.0009303791403858427308*h0*h2/h1;
    const double w4 = -0.00024929433932114870101*h0;
    const double w5 = 0.0009303791403858427308*h1;
    const double w6 = 0.0009303791403858427308*h0;
    const double w7 = -0.00024929433932114870101*h0*h1/h2;
    const double w8 = 0.0034722222222222222222*h2;
    const double w9 = -0.0009303791403858427308*h1;
    const double w10 = 0.012958509748503046158*h0*h2/h1;
    const double w11 = -0.0034722222222222222222*h0;
    const double w12 = 0.0034722222222222222222*h1;
    const double w13 = 0.012958509748503046158*h0;
    const double w14 = -0.0034722222222222222222*h0*h1/h2;
    const double w15 = 0.012958509748503046158*h1*h2/h0;
    const double w16 = -0.0034722222222222222222*h1;
    const double w17 = -0.0009303791403858427308*h0;
    const double w18 = 0.012958509748503046158*h1;
    const double w19 = 0.0034722222222222222222*h0;
    const double w20 = 0.012958509748503046158*h2;
    const double w21 = -0.012958509748503046158*h1;
    const double w22 = -0.012958509748503046158*h0;
    const double w23 = 0.04836181677178996241*h1;
    const double w24 = 0.04836181677178996241*h0;
    const double w25 = -0.04836181677178996241*h0*h1/h2;
    const double w26 = 0.00024929433932114870101*h1;
    const double w27 = 0.00024929433932114870101*h0;
    const double w28 = -0.04836181677178996241*h1;
    const double w29 = -0.04836181677178996241*h0;
    const double w30 = -0.0009303791403858427308*h1*h2/h0;
    const double w31 = -0.0009303791403858427308*h2;
    const double w32 = -0.0009303791403858427308*h0*h2/h1;
    const double w33 = 0.0034722222222222222222*h0*h1/h2;
    const double w34 = -0.0034722222222222222222*h2;
    const double w35 = -0.00024929433932114870101*h2;
    const double w36 = -0.012958509748503046158*h1*h2/h0;
    const double w37 = -0.012958509748503046158*h2;
    const double w38 = -0.012958509748503046158*h0*h2/h1;
    const double w39 = -0.04836181677178996241*h2;
    const double w40 = -0.0034722222222222222222*h0*h2/h1;
    const double w41 = 0.0009303791403858427308*h0*h1/h2;
    const double w42 = 0.04836181677178996241*h2;
    const double w43 = -0.04836181677178996241*h0*h2/h1;
    const double w44 = 0.012958509748503046158*h0*h1/h2;
    const double w45 = -0.00024929433932114870101*h0*h2/h1;
    const double w46 = 0.00024929433932114870101*h2;
    const double w47 = -0.0034722222222222222222*h1*h2/h0;
    const double w48 = -0.00024929433932114870101*h1*h2/h0;
    const double w49 = -0.04836181677178996241*h1*h2/h0;
    const double w50 = 0.0034722222222222222222*h0*h2/h1;
    const double w51 = -0.0009303791403858427308*h0*h1/h2;
    const double w52 = -0.012958509748503046158*h0*h1/h2;
    const double w53 = 0.0034722222222222222222*h1*h2/h0;
    const double w54 = 0.00024929433932114870101*h0*h1/h2;
    const double w55 = 0.04836181677178996241*h0*h2/h1;
    const double w56 = 0.04836181677178996241*h1*h2/h0;
    const double w57 = 0.04836181677178996241*h0*h1/h2;
    const double w58 = 0.00024929433932114870101*h1*h2/h0;
    const double w59 = 0.00024929433932114870101*h0*h2/h1;
    const double w60 = 0.055555555555555555556*h1*h2/h0;
    const double w61 = 0.041666666666666666667*h2;
    const double w62 = -0.083333333333333333333*h1;
    const double w63 = 0.055555555555555555556*h0*h2/h1;
    const double w64 = -0.083333333333333333333*h0;
    const double w65 = 0.083333333333333333333*h1;
    const double w66 = 0.083333333333333333333*h0;
    const double w67 = -0.11111111111111111111*h0*h1/h2;
    const double w68 = -0.055555555555555555556*h1*h2/h0;
    const double w69 = -0.083333333333333333333*h2;
    const double w70 = -0.041666666666666666667*h1;
    const double w71 = -0.055555555555555555556*h0*h2/h1;
    const double w72 = -0.041666666666666666667*h0;
    const double w73 = 0.041666666666666666667*h1;
    const double w74 = 0.041666666666666666667*h0;
    const double w75 = 0.027777777777777777778*h0*h1/h2;
    const double w76 = 0.083333333333333333333*h2;
    const double w77 = -0.11111111111111111111*h0*h2/h1;
    const double w78 = 0.055555555555555555556*h0*h1/h2;
    const double w79 = -0.11111111111111111111*h1*h2/h0;
    const double w80 = -0.027777777777777777778*h1*h2/h0;
    const double w81 = -0.041666666666666666667*h2;
    const double w82 = -0.027777777777777777778*h0*h2/h1;
    const double w83 = -0.027777777777777777778*h0*h1/h2;
    const double w84 = 0.027777777777777777778*h0*h2/h1;
    const double w85 = -0.055555555555555555556*h0*h1/h2;
    const double w86 = 0.11111111111111111111*h1*h2/h0;
    const double w87 = 0.11111111111111111111*h0*h2/h1;
    const double w88 = 0.11111111111111111111*h0*h1/h2;
    const double w89 = 0.027777777777777777778*h1*h2/h0;
    const double w90 = 0.0001966122466178319053*h1*h2;
    const double w91 = 0.0001966122466178319053*h0*h2;
    const double w92 = 0.0001966122466178319053*h0*h1;
    const double w93 = 0.0007337668937680108255*h1*h2;
    const double w94 = 0.0027384553284542113967*h0*h2;
    const double w95 = 0.0027384553284542113967*h0*h1;
    const double w96 = 0.0027384553284542113967*h1*h2;
    const double w97 = 0.0007337668937680108255*h0*h2;
    const double w98 = 0.010220054420048834761*h1*h2;
    const double w99 = 0.010220054420048834761*h0*h2;
    const double w100 = 0.038141762351741127649*h0*h1;
    const double w101 = 0.000052682092703316795705*h0*h1;
    const double w102 = 0.0007337668937680108255*h0*h1;
    const double w103 = 0.010220054420048834761*h0*h1;
    const double w104 = -0.0001966122466178319053*h1*h2;
    const double w105 = -0.0001966122466178319053*h0*h2;
    const double w106 = -0.0007337668937680108255*h1*h2;
    const double w107 = -0.0007337668937680108255*h0*h2;
    const double w108 = -0.0027384553284542113967*h1*h2;
    const double w109 = -0.0027384553284542113967*h0*h2;
    const double w110 = -0.010220054420048834761*h1*h2;
    const double w111 = -0.010220054420048834761*h0*h2;
    const double w112 = -0.0007337668937680108255*h0*h1;
    const double w113 = -0.010220054420048834761*h0*h1;
    const double w114 = -0.038141762351741127649*h0*h2;
    const double w115 = -0.000052682092703316795705*h0*h2;
    const double w116 = -0.0001966122466178319053*h0*h1;
    const double w117 = -0.0027384553284542113967*h0*h1;
    const double w118 = 0.000052682092703316795705*h0*h2;
    const double w119 = 0.038141762351741127649*h0*h2;
    const double w120 = 0.000052682092703316795705*h1*h2;
    const double w121 = 0.038141762351741127649*h1*h2;
    const double w122 = -0.000052682092703316795705*h0*h1;
    const double w123 = -0.038141762351741127649*h0*h1;
    const double w124 = -0.000052682092703316795705*h1*h2;
    const double w125 = -0.038141762351741127649*h1*h2;
    const double w126 = 0.027777777777777777778*h1*h2;
    const double w127 = 0.027777777777777777778*h0*h2;
    const double w128 = 0.055555555555555555556*h0*h1;
    const double w129 = -0.027777777777777777778*h1*h2;
    const double w130 = -0.027777777777777777778*h0*h2;
    const double w131 = 0.013888888888888888889*h0*h1;
    const double w132 = -0.055555555555555555556*h0*h2;
    const double w133 = -0.027777777777777777778*h0*h1;
    const double w134 = 0.055555555555555555556*h0*h2;
    const double w135 = 0.027777777777777777778*h0*h1;
    const double w136 = -0.013888888888888888889*h0*h1;
    const double w137 = 0.055555555555555555556*h1*h2;
    const double w138 = -0.013888888888888888889*h1*h2;
    const double w139 = -0.013888888888888888889*h0*h2;
    const double w140 = -0.055555555555555555556*h0*h1;
    const double w141 = 0.013888888888888888889*h1*h2;
    const double w142 = 0.013888888888888888889*h0*h2;
    const double w143 = -0.055555555555555555556*h1*h2;
    const double w144 = 0.000041549056553524783501*h0*h1*h2;
    const double w145 = 0.0005787037037037037037*h0*h1*h2;
    const double w146 = 0.0080603027952983270684*h0*h1*h2;
    const double w147 = 0.0001550631900643071218*h0*h1*h2;
    const double w148 = 0.002159751624750507693*h0*h1*h2;
    const double w149 = 0.03008145955644280058*h0*h1*h2;
    const double w150 = 0.000011133036149792012204*h0*h1*h2;
    const double w151 = 0.018518518518518518519*h0*h1*h2;
    const double w152 = 0.0092592592592592592592*h0*h1*h2;
    const double w153 = 0.0046296296296296296296*h0*h1*h2;
    const double w154 = 0.037037037037037037037*h0*h1*h2;
    const double w155 = -0.077751058491018276949*h1*h2;
    const double w156 = -0.077751058491018276949*h0*h2;
    const double w157 = -0.077751058491018276949*h0*h1;
    const double w158 = -0.020833333333333333333*h0*h2;
    const double w159 = -0.020833333333333333333*h0*h1;
    const double w160 = -0.020833333333333333333*h1*h2;
    const double w161 = -0.0055822748423150563848*h0*h1;
    const double w162 = -0.0055822748423150563848*h0*h2;
    const double w163 = -0.0055822748423150563848*h1*h2;
    const double w164 = 0.077751058491018276949*h1*h2;
    const double w165 = 0.020833333333333333333*h1*h2;
    const double w166 = 0.0055822748423150563848*h1*h2;
    const double w167 = 0.077751058491018276949*h0*h2;
    const double w168 = 0.020833333333333333333*h0*h2;
    const double w169 = 0.0055822748423150563848*h0*h2;
    const double w170 = 0.077751058491018276949*h0*h1;
    const double w171 = 0.020833333333333333333*h0*h1;
    const double w172 = 0.0055822748423150563848*h0*h1;
    const double w173 = -0.25*h1*h2;
    const double w174 = -0.25*h0*h2;
    const double w175 = -0.25*h0*h1;
    const double w176 = 0.25*h1*h2;
    const double w177 = 0.25*h0*h2;
    const double w178 = 0.25*h0*h1;
    const double w179 = 0.061320326520293008568*h0*h1*h2;
    const double w180 = 0.01643073197072526838*h0*h1*h2;
    const double w181 = 0.004402601362608064953*h0*h1*h2;
    const double w182 = 0.0011796734797069914318*h0*h1*h2;
    const double w183 = 0.125*h0*h1*h2;

    rhs.requireWrite();
#pragma omp parallel
    {
        for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
#pragma omp for
            for (index_t k2=k2_0; k2<m_NE2; k2+=2) {
                for (index_t k1=0; k1<m_NE1; ++k1) {
                    for (index_t k0=0; k0<m_NE0; ++k0)  {
                        bool add_EM_S=false;
                        bool add_EM_F=false;
                        vector<double> EM_S(8*8*numEq*numComp, 0);
                        vector<double> EM_F(8*numEq, 0);
                        const index_t e = k0 + m_NE0*k1 + m_NE0*m_NE1*k2;
                        ///////////////
                        // process A //
                        ///////////////
                        if (!A.isEmpty()) {
                            add_EM_S=true;
                            const double* A_p=const_cast<escript::Data*>(&A)->getSampleDataRO(e);
                            if (A.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const register double A_00_0 = A_p[INDEX5(k,0,m,0,0, numEq,3,numComp,3)];
                                        const register double A_01_0 = A_p[INDEX5(k,0,m,1,0, numEq,3,numComp,3)];
                                        const register double A_02_0 = A_p[INDEX5(k,0,m,2,0, numEq,3,numComp,3)];
                                        const register double A_10_0 = A_p[INDEX5(k,1,m,0,0, numEq,3,numComp,3)];
                                        const register double A_11_0 = A_p[INDEX5(k,1,m,1,0, numEq,3,numComp,3)];
                                        const register double A_12_0 = A_p[INDEX5(k,1,m,2,0, numEq,3,numComp,3)];
                                        const register double A_20_0 = A_p[INDEX5(k,2,m,0,0, numEq,3,numComp,3)];
                                        const register double A_21_0 = A_p[INDEX5(k,2,m,1,0, numEq,3,numComp,3)];
                                        const register double A_22_0 = A_p[INDEX5(k,2,m,2,0, numEq,3,numComp,3)];
                                        const register double A_00_1 = A_p[INDEX5(k,0,m,0,1, numEq,3,numComp,3)];
                                        const register double A_01_1 = A_p[INDEX5(k,0,m,1,1, numEq,3,numComp,3)];
                                        const register double A_02_1 = A_p[INDEX5(k,0,m,2,1, numEq,3,numComp,3)];
                                        const register double A_10_1 = A_p[INDEX5(k,1,m,0,1, numEq,3,numComp,3)];
                                        const register double A_11_1 = A_p[INDEX5(k,1,m,1,1, numEq,3,numComp,3)];
                                        const register double A_12_1 = A_p[INDEX5(k,1,m,2,1, numEq,3,numComp,3)];
                                        const register double A_20_1 = A_p[INDEX5(k,2,m,0,1, numEq,3,numComp,3)];
                                        const register double A_21_1 = A_p[INDEX5(k,2,m,1,1, numEq,3,numComp,3)];
                                        const register double A_22_1 = A_p[INDEX5(k,2,m,2,1, numEq,3,numComp,3)];
                                        const register double A_00_2 = A_p[INDEX5(k,0,m,0,2, numEq,3,numComp,3)];
                                        const register double A_01_2 = A_p[INDEX5(k,0,m,1,2, numEq,3,numComp,3)];
                                        const register double A_02_2 = A_p[INDEX5(k,0,m,2,2, numEq,3,numComp,3)];
                                        const register double A_10_2 = A_p[INDEX5(k,1,m,0,2, numEq,3,numComp,3)];
                                        const register double A_11_2 = A_p[INDEX5(k,1,m,1,2, numEq,3,numComp,3)];
                                        const register double A_12_2 = A_p[INDEX5(k,1,m,2,2, numEq,3,numComp,3)];
                                        const register double A_20_2 = A_p[INDEX5(k,2,m,0,2, numEq,3,numComp,3)];
                                        const register double A_21_2 = A_p[INDEX5(k,2,m,1,2, numEq,3,numComp,3)];
                                        const register double A_22_2 = A_p[INDEX5(k,2,m,2,2, numEq,3,numComp,3)];
                                        const register double A_00_3 = A_p[INDEX5(k,0,m,0,3, numEq,3,numComp,3)];
                                        const register double A_01_3 = A_p[INDEX5(k,0,m,1,3, numEq,3,numComp,3)];
                                        const register double A_02_3 = A_p[INDEX5(k,0,m,2,3, numEq,3,numComp,3)];
                                        const register double A_10_3 = A_p[INDEX5(k,1,m,0,3, numEq,3,numComp,3)];
                                        const register double A_11_3 = A_p[INDEX5(k,1,m,1,3, numEq,3,numComp,3)];
                                        const register double A_12_3 = A_p[INDEX5(k,1,m,2,3, numEq,3,numComp,3)];
                                        const register double A_20_3 = A_p[INDEX5(k,2,m,0,3, numEq,3,numComp,3)];
                                        const register double A_21_3 = A_p[INDEX5(k,2,m,1,3, numEq,3,numComp,3)];
                                        const register double A_22_3 = A_p[INDEX5(k,2,m,2,3, numEq,3,numComp,3)];
                                        const register double A_00_4 = A_p[INDEX5(k,0,m,0,4, numEq,3,numComp,3)];
                                        const register double A_01_4 = A_p[INDEX5(k,0,m,1,4, numEq,3,numComp,3)];
                                        const register double A_02_4 = A_p[INDEX5(k,0,m,2,4, numEq,3,numComp,3)];
                                        const register double A_10_4 = A_p[INDEX5(k,1,m,0,4, numEq,3,numComp,3)];
                                        const register double A_11_4 = A_p[INDEX5(k,1,m,1,4, numEq,3,numComp,3)];
                                        const register double A_12_4 = A_p[INDEX5(k,1,m,2,4, numEq,3,numComp,3)];
                                        const register double A_20_4 = A_p[INDEX5(k,2,m,0,4, numEq,3,numComp,3)];
                                        const register double A_21_4 = A_p[INDEX5(k,2,m,1,4, numEq,3,numComp,3)];
                                        const register double A_22_4 = A_p[INDEX5(k,2,m,2,4, numEq,3,numComp,3)];
                                        const register double A_00_5 = A_p[INDEX5(k,0,m,0,5, numEq,3,numComp,3)];
                                        const register double A_01_5 = A_p[INDEX5(k,0,m,1,5, numEq,3,numComp,3)];
                                        const register double A_02_5 = A_p[INDEX5(k,0,m,2,5, numEq,3,numComp,3)];
                                        const register double A_10_5 = A_p[INDEX5(k,1,m,0,5, numEq,3,numComp,3)];
                                        const register double A_11_5 = A_p[INDEX5(k,1,m,1,5, numEq,3,numComp,3)];
                                        const register double A_12_5 = A_p[INDEX5(k,1,m,2,5, numEq,3,numComp,3)];
                                        const register double A_20_5 = A_p[INDEX5(k,2,m,0,5, numEq,3,numComp,3)];
                                        const register double A_21_5 = A_p[INDEX5(k,2,m,1,5, numEq,3,numComp,3)];
                                        const register double A_22_5 = A_p[INDEX5(k,2,m,2,5, numEq,3,numComp,3)];
                                        const register double A_00_6 = A_p[INDEX5(k,0,m,0,6, numEq,3,numComp,3)];
                                        const register double A_01_6 = A_p[INDEX5(k,0,m,1,6, numEq,3,numComp,3)];
                                        const register double A_02_6 = A_p[INDEX5(k,0,m,2,6, numEq,3,numComp,3)];
                                        const register double A_10_6 = A_p[INDEX5(k,1,m,0,6, numEq,3,numComp,3)];
                                        const register double A_11_6 = A_p[INDEX5(k,1,m,1,6, numEq,3,numComp,3)];
                                        const register double A_12_6 = A_p[INDEX5(k,1,m,2,6, numEq,3,numComp,3)];
                                        const register double A_20_6 = A_p[INDEX5(k,2,m,0,6, numEq,3,numComp,3)];
                                        const register double A_21_6 = A_p[INDEX5(k,2,m,1,6, numEq,3,numComp,3)];
                                        const register double A_22_6 = A_p[INDEX5(k,2,m,2,6, numEq,3,numComp,3)];
                                        const register double A_00_7 = A_p[INDEX5(k,0,m,0,7, numEq,3,numComp,3)];
                                        const register double A_01_7 = A_p[INDEX5(k,0,m,1,7, numEq,3,numComp,3)];
                                        const register double A_02_7 = A_p[INDEX5(k,0,m,2,7, numEq,3,numComp,3)];
                                        const register double A_10_7 = A_p[INDEX5(k,1,m,0,7, numEq,3,numComp,3)];
                                        const register double A_11_7 = A_p[INDEX5(k,1,m,1,7, numEq,3,numComp,3)];
                                        const register double A_12_7 = A_p[INDEX5(k,1,m,2,7, numEq,3,numComp,3)];
                                        const register double A_20_7 = A_p[INDEX5(k,2,m,0,7, numEq,3,numComp,3)];
                                        const register double A_21_7 = A_p[INDEX5(k,2,m,1,7, numEq,3,numComp,3)];
                                        const register double A_22_7 = A_p[INDEX5(k,2,m,2,7, numEq,3,numComp,3)];
                                        const register double tmp160_0 = A_12_0 + A_12_6 + A_21_0 + A_21_6;
                                        const register double tmp8_0 = A_21_0 + A_21_6;
                                        const register double tmp135_0 = A_10_1 + A_10_2 + A_10_5 + A_10_6;
                                        const register double tmp67_0 = A_02_2 + A_02_7;
                                        const register double tmp211_0 = A_12_6 + A_21_6;
                                        const register double tmp180_0 = A_10_2 + A_10_6;
                                        const register double tmp37_0 = A_00_0 + A_00_1 + A_00_2 + A_00_3;
                                        const register double tmp92_0 = A_11_0 + A_11_1 + A_11_2 + A_11_3 + A_11_4 + A_11_5 + A_11_6 + A_11_7;
                                        const register double tmp195_0 = A_02_2 + A_20_2;
                                        const register double tmp70_0 = A_01_0 + A_01_7;
                                        const register double tmp139_0 = A_02_3 + A_02_4 + A_20_1 + A_20_6;
                                        const register double tmp200_0 = A_12_3 + A_12_5 + A_21_3 + A_21_5;
                                        const register double tmp60_0 = A_22_0 + A_22_2 + A_22_4 + A_22_6;
                                        const register double tmp192_0 = A_01_5 + A_10_5;
                                        const register double tmp46_0 = A_21_0 + A_21_7;
                                        const register double tmp48_0 = A_10_0 + A_10_7;
                                        const register double tmp166_0 = A_11_5 + A_11_7;
                                        const register double tmp221_0 = A_02_1 + A_02_6 + A_20_3 + A_20_4;
                                        const register double tmp50_0 = A_02_4 + A_02_6 + A_20_4 + A_20_6;
                                        const register double tmp217_0 = A_02_3 + A_02_4 + A_20_3 + A_20_4;
                                        const register double tmp216_0 = A_01_2 + A_01_5 + A_10_2 + A_10_5;
                                        const register double tmp104_0 = A_22_2 + A_22_6;
                                        const register double tmp72_0 = A_20_3 + A_20_6;
                                        const register double tmp79_0 = A_10_4 + A_10_7;
                                        const register double tmp86_0 = A_01_2 + A_01_6 + A_10_1 + A_10_5;
                                        const register double tmp214_0 = A_12_0 + A_12_7 + A_21_0 + A_21_7;
                                        const register double tmp32_0 = A_02_0 + A_02_2;
                                        const register double tmp112_0 = A_01_0 + A_01_4 + A_10_3 + A_10_7;
                                        const register double tmp197_0 = A_12_0 + A_21_0;
                                        const register double tmp106_0 = A_22_1 + A_22_5;
                                        const register double tmp2_0 = A_00_0 + A_00_1 + A_00_4 + A_00_5;
                                        const register double tmp115_0 = A_02_5 + A_02_7 + A_20_0 + A_20_2;
                                        const register double tmp175_0 = A_01_3 + A_01_7;
                                        const register double tmp126_0 = A_01_2 + A_01_5 + A_10_1 + A_10_6;
                                        const register double tmp90_0 = A_00_0 + A_00_1 + A_00_2 + A_00_3 + A_00_4 + A_00_5 + A_00_6 + A_00_7;
                                        const register double tmp47_0 = A_12_0 + A_12_6;
                                        const register double tmp205_0 = A_02_7 + A_20_7;
                                        const register double tmp148_0 = A_01_3 + A_01_4;
                                        const register double tmp113_0 = A_01_3 + A_01_7 + A_10_0 + A_10_4;
                                        const register double tmp43_0 = A_20_4 + A_20_6;
                                        const register double tmp161_0 = A_02_1 + A_02_6 + A_20_1 + A_20_6;
                                        const register double tmp69_0 = A_12_0 + A_12_1 + A_12_6 + A_12_7 + A_21_0 + A_21_1 + A_21_6 + A_21_7;
                                        const register double tmp176_0 = A_01_1 + A_01_2 + A_01_5 + A_01_6;
                                        const register double tmp105_0 = A_01_2 + A_01_6 + A_10_2 + A_10_6;
                                        const register double tmp22_0 = A_01_5 + A_10_6;
                                        const register double tmp91_0 = A_02_4 + A_02_6 + A_20_1 + A_20_3;
                                        const register double tmp206_0 = A_12_7 + A_21_7;
                                        const register double tmp188_0 = A_02_5 + A_20_5;
                                        const register double tmp117_0 = A_21_1 + A_21_6;
                                        const register double tmp165_0 = A_01_1 + A_01_6;
                                        const register double tmp66_0 = A_00_4 + A_00_5;
                                        const register double tmp57_0 = A_02_0 + A_02_2 + A_02_5 + A_02_7 + A_20_0 + A_20_2 + A_20_5 + A_20_7;
                                        const register double tmp31_0 = A_21_4 + A_21_5;
                                        const register double tmp3_0 = A_11_0 + A_11_2 + A_11_4 + A_11_6;
                                        const register double tmp183_0 = A_12_0 + A_12_7;
                                        const register double tmp61_0 = A_02_1 + A_02_3 + A_20_1 + A_20_3;
                                        const register double tmp54_0 = A_10_5 + A_10_6;
                                        const register double tmp18_0 = A_02_3 + A_02_6;
                                        const register double tmp119_0 = A_12_2 + A_12_3 + A_12_4 + A_12_5 + A_21_2 + A_21_3 + A_21_4 + A_21_5;
                                        const register double tmp29_0 = A_21_2 + A_21_3;
                                        const register double tmp17_0 = A_01_3 + A_01_7 + A_10_3 + A_10_7;
                                        const register double tmp212_0 = A_02_6 + A_20_6;
                                        const register double tmp220_0 = A_02_3 + A_20_6;
                                        const register double tmp78_0 = A_20_0 + A_20_7;
                                        const register double tmp215_0 = A_01_6 + A_10_6;
                                        const register double tmp203_0 = A_01_7 + A_10_7;
                                        const register double tmp87_0 = A_12_2 + A_12_3 + A_21_4 + A_21_5;
                                        const register double tmp114_0 = A_02_0 + A_02_2 + A_20_5 + A_20_7;
                                        const register double tmp0_0 = A_01_0 + A_01_4 + A_10_0 + A_10_4;
                                        const register double tmp202_0 = A_01_3 + A_01_4 + A_10_3 + A_10_4;
                                        const register double tmp4_0 = A_20_0 + A_20_5;
                                        const register double tmp65_0 = A_00_2 + A_00_3;
                                        const register double tmp24_0 = A_20_1 + A_20_3;
                                        const register double tmp64_0 = A_10_0 + A_10_3;
                                        const register double tmp170_0 = A_02_0 + A_02_2 + A_20_0 + A_20_2;
                                        const register double tmp11_0 = A_20_1 + A_20_6;
                                        const register double tmp82_0 = A_12_4 + A_12_5 + A_21_4 + A_21_5;
                                        const register double tmp99_0 = A_01_4 + A_10_7;
                                        const register double tmp49_0 = A_12_1 + A_12_7;
                                        const register double tmp130_0 = A_12_0 + A_12_1 + A_12_6 + A_12_7;
                                        const register double tmp144_0 = A_01_0 + A_10_3;
                                        const register double tmp109_0 = A_22_0 + A_22_3 + A_22_4 + A_22_7;
                                        const register double tmp185_0 = A_02_0 + A_02_7 + A_20_2 + A_20_5;
                                        const register double tmp157_0 = A_01_4 + A_10_4;
                                        const register double tmp51_0 = A_22_1 + A_22_3 + A_22_5 + A_22_7;
                                        const register double tmp146_0 = A_00_6 + A_00_7;
                                        const register double tmp147_0 = A_12_0 + A_12_1 + A_21_0 + A_21_1;
                                        const register double tmp150_0 = A_00_2 + A_00_3 + A_00_4 + A_00_5;
                                        const register double tmp62_0 = A_21_3 + A_21_5;
                                        const register double tmp223_0 = A_12_2 + A_21_4;
                                        const register double tmp16_0 = A_02_2 + A_02_5;
                                        const register double tmp168_0 = A_11_1 + A_11_3 + A_11_4 + A_11_6;
                                        const register double tmp88_0 = A_12_4 + A_12_5 + A_21_2 + A_21_3;
                                        const register double tmp142_0 = A_01_7 + A_10_4;
                                        const register double tmp34_0 = A_20_0 + A_20_2 + A_20_5 + A_20_7;
                                        const register double tmp71_0 = A_00_0 + A_00_1 + A_00_6 + A_00_7;
                                        const register double tmp213_0 = A_02_1 + A_20_1;
                                        const register double tmp227_0 = A_12_2 + A_12_5 + A_21_3 + A_21_4;
                                        const register double tmp228_0 = A_12_1 + A_21_7;
                                        const register double tmp140_0 = A_01_2 + A_01_6;
                                        const register double tmp74_0 = A_22_0 + A_22_1 + A_22_4 + A_22_5;
                                        const register double tmp167_0 = A_11_0 + A_11_2;
                                        const register double tmp143_0 = A_01_3 + A_01_4 + A_10_0 + A_10_7;
                                        const register double tmp83_0 = A_02_0 + A_02_5;
                                        const register double tmp14_0 = A_22_1 + A_22_2 + A_22_5 + A_22_6;
                                        const register double tmp5_0 = A_12_1 + A_12_6;
                                        const register double tmp94_0 = A_02_1 + A_02_3;
                                        const register double tmp193_0 = A_01_1 + A_01_6 + A_10_1 + A_10_6;
                                        const register double tmp97_0 = A_02_0 + A_02_2 + A_02_5 + A_02_7;
                                        const register double tmp131_0 = A_01_1 + A_01_5;
                                        const register double tmp124_0 = A_01_6 + A_10_5;
                                        const register double tmp149_0 = A_12_6 + A_12_7 + A_21_6 + A_21_7;
                                        const register double tmp187_0 = A_01_2 + A_10_2;
                                        const register double tmp93_0 = A_01_1 + A_01_2 + A_10_1 + A_10_2;
                                        const register double tmp25_0 = A_01_4 + A_01_7 + A_10_4 + A_10_7;
                                        const register double tmp156_0 = A_12_2 + A_12_5 + A_21_2 + A_21_5;
                                        const register double tmp20_0 = A_21_2 + A_21_5;
                                        const register double tmp55_0 = A_21_2 + A_21_4;
                                        const register double tmp208_0 = A_12_1 + A_12_6 + A_21_0 + A_21_7;
                                        const register double tmp125_0 = A_12_4 + A_12_5;
                                        const register double tmp158_0 = A_01_0 + A_01_7 + A_10_0 + A_10_7;
                                        const register double tmp108_0 = A_01_1 + A_01_5 + A_10_1 + A_10_5;
                                        const register double tmp199_0 = A_12_2 + A_12_4 + A_21_2 + A_21_4;
                                        const register double tmp10_0 = A_02_1 + A_02_4;
                                        const register double tmp182_0 = A_02_3 + A_02_6 + A_20_3 + A_20_6;
                                        const register double tmp132_0 = A_02_1 + A_20_4;
                                        const register double tmp191_0 = A_12_3 + A_12_4 + A_21_3 + A_21_4;
                                        const register double tmp35_0 = A_11_0 + A_11_1 + A_11_2 + A_11_3;
                                        const register double tmp164_0 = A_10_3 + A_10_4;
                                        const register double tmp190_0 = A_12_5 + A_21_5;
                                        const register double tmp73_0 = A_02_1 + A_02_6;
                                        const register double tmp98_0 = A_01_0 + A_01_7 + A_10_3 + A_10_4;
                                        const register double tmp225_0 = A_12_4 + A_21_2;
                                        const register double tmp103_0 = A_02_4 + A_02_6;
                                        const register double tmp194_0 = A_02_0 + A_02_7 + A_20_0 + A_20_7;
                                        const register double tmp207_0 = A_12_0 + A_21_6;
                                        const register double tmp102_0 = A_20_5 + A_20_7;
                                        const register double tmp1_0 = A_22_3 + A_22_7;
                                        const register double tmp172_0 = A_10_1 + A_10_5;
                                        const register double tmp222_0 = A_12_5 + A_21_3;
                                        const register double tmp201_0 = A_02_2 + A_02_5 + A_20_2 + A_20_5;
                                        const register double tmp155_0 = A_12_4 + A_21_4;
                                        const register double tmp174_0 = A_02_1 + A_02_4 + A_20_1 + A_20_4;
                                        const register double tmp59_0 = A_01_0 + A_01_3;
                                        const register double tmp21_0 = A_20_2 + A_20_7;
                                        const register double tmp141_0 = A_02_2 + A_02_7 + A_20_2 + A_20_7;
                                        const register double tmp210_0 = A_01_1 + A_10_1;
                                        const register double tmp145_0 = A_00_0 + A_00_1;
                                        const register double tmp121_0 = A_12_0 + A_12_1 + A_21_6 + A_21_7;
                                        const register double tmp224_0 = A_12_3 + A_12_4 + A_21_2 + A_21_5;
                                        const register double tmp186_0 = A_02_2 + A_20_7;
                                        const register double tmp53_0 = A_11_4 + A_11_6;
                                        const register double tmp184_0 = A_02_5 + A_20_0;
                                        const register double tmp38_0 = A_12_0 + A_12_1;
                                        const register double tmp12_0 = A_01_1 + A_01_2 + A_01_5 + A_01_6 + A_10_1 + A_10_2 + A_10_5 + A_10_6;
                                        const register double tmp230_0 = A_12_6 + A_21_0;
                                        const register double tmp23_0 = A_11_4 + A_11_5 + A_11_6 + A_11_7;
                                        const register double tmp81_0 = A_20_1 + A_20_4;
                                        const register double tmp134_0 = A_10_3 + A_10_7;
                                        const register double tmp129_0 = A_21_0 + A_21_1;
                                        const register double tmp137_0 = A_01_0 + A_01_3 + A_01_4 + A_01_7;
                                        const register double tmp198_0 = A_01_0 + A_10_0;
                                        const register double tmp9_0 = A_21_1 + A_21_7;
                                        const register double tmp179_0 = A_01_0 + A_01_4;
                                        const register double tmp100_0 = A_20_1 + A_20_3 + A_20_4 + A_20_6;
                                        const register double tmp173_0 = A_02_0 + A_20_5;
                                        const register double tmp42_0 = A_21_0 + A_21_1 + A_21_6 + A_21_7;
                                        const register double tmp226_0 = A_12_3 + A_21_5;
                                        const register double tmp6_0 = A_22_0 + A_22_4;
                                        const register double tmp218_0 = A_12_1 + A_21_1;
                                        const register double tmp28_0 = A_01_2 + A_10_1;
                                        const register double tmp133_0 = A_02_6 + A_20_3;
                                        const register double tmp13_0 = A_00_2 + A_00_3 + A_00_6 + A_00_7;
                                        const register double tmp27_0 = A_12_2 + A_12_3 + A_12_4 + A_12_5;
                                        const register double tmp75_0 = A_10_1 + A_10_6;
                                        const register double tmp36_0 = A_01_0 + A_01_3 + A_10_0 + A_10_3;
                                        const register double tmp138_0 = A_10_0 + A_10_4;
                                        const register double tmp189_0 = A_12_2 + A_21_2;
                                        const register double tmp181_0 = A_02_7 + A_20_2;
                                        const register double tmp85_0 = A_02_1 + A_02_3 + A_20_4 + A_20_6;
                                        const register double tmp122_0 = A_01_1 + A_10_2;
                                        const register double tmp95_0 = A_01_3 + A_10_0;
                                        const register double tmp120_0 = A_12_6 + A_12_7 + A_21_0 + A_21_1;
                                        const register double tmp196_0 = A_02_0 + A_20_0;
                                        const register double tmp171_0 = A_02_3 + A_02_4;
                                        const register double tmp204_0 = A_12_1 + A_12_6 + A_21_1 + A_21_6;
                                        const register double tmp45_0 = A_10_1 + A_10_2;
                                        const register double tmp101_0 = A_01_5 + A_01_6 + A_10_5 + A_10_6;
                                        const register double tmp58_0 = A_11_0 + A_11_2 + A_11_5 + A_11_7;
                                        const register double tmp107_0 = A_20_3 + A_20_4;
                                        const register double tmp30_0 = A_01_1 + A_01_6 + A_10_2 + A_10_5;
                                        const register double tmp63_0 = A_12_2 + A_12_5;
                                        const register double tmp127_0 = A_12_2 + A_12_3;
                                        const register double tmp177_0 = A_02_2 + A_02_5 + A_20_0 + A_20_7;
                                        const register double tmp178_0 = A_10_0 + A_10_3 + A_10_4 + A_10_7;
                                        const register double tmp76_0 = A_01_1 + A_01_2;
                                        const register double tmp80_0 = A_22_2 + A_22_3 + A_22_6 + A_22_7;
                                        const register double tmp41_0 = A_12_6 + A_12_7;
                                        const register double tmp89_0 = A_01_0 + A_01_3 + A_01_4 + A_01_7 + A_10_0 + A_10_3 + A_10_4 + A_10_7;
                                        const register double tmp116_0 = A_02_1 + A_02_3 + A_02_4 + A_02_6 + A_20_1 + A_20_3 + A_20_4 + A_20_6;
                                        const register double tmp33_0 = A_22_0 + A_22_1 + A_22_2 + A_22_3 + A_22_4 + A_22_5 + A_22_6 + A_22_7;
                                        const register double tmp169_0 = A_21_3 + A_21_4;
                                        const register double tmp96_0 = A_20_0 + A_20_2;
                                        const register double tmp111_0 = A_12_3 + A_12_4;
                                        const register double tmp118_0 = A_20_2 + A_20_5;
                                        const register double tmp19_0 = A_12_3 + A_12_5;
                                        const register double tmp68_0 = A_01_5 + A_01_6;
                                        const register double tmp7_0 = A_11_1 + A_11_3 + A_11_5 + A_11_7;
                                        const register double tmp154_0 = A_12_3 + A_21_3;
                                        const register double tmp152_0 = A_02_4 + A_20_4;
                                        const register double tmp153_0 = A_02_3 + A_20_3;
                                        const register double tmp163_0 = A_02_5 + A_02_7 + A_20_5 + A_20_7;
                                        const register double tmp44_0 = A_01_4 + A_01_7;
                                        const register double tmp39_0 = A_02_1 + A_02_3 + A_02_4 + A_02_6;
                                        const register double tmp123_0 = A_21_2 + A_21_3 + A_21_4 + A_21_5;
                                        const register double tmp40_0 = A_02_5 + A_02_7;
                                        const register double tmp110_0 = A_02_0 + A_02_7;
                                        const register double tmp77_0 = A_12_2 + A_12_3 + A_21_2 + A_21_3;
                                        const register double tmp209_0 = A_12_7 + A_21_1;
                                        const register double tmp219_0 = A_02_4 + A_20_1;
                                        const register double tmp84_0 = A_01_1 + A_01_5 + A_10_2 + A_10_6;
                                        const register double tmp162_0 = A_12_1 + A_12_7 + A_21_1 + A_21_7;
                                        const register double tmp159_0 = A_01_3 + A_10_3;
                                        const register double tmp56_0 = A_11_1 + A_11_3;
                                        const register double tmp52_0 = A_01_2 + A_01_5;
                                        const register double tmp26_0 = A_00_4 + A_00_5 + A_00_6 + A_00_7;
                                        const register double tmp229_0 = A_12_0 + A_12_7 + A_21_1 + A_21_6;
                                        const register double tmp151_0 = A_10_2 + A_10_5;
                                        const register double tmp136_0 = A_02_0 + A_02_5 + A_20_0 + A_20_5;
                                        const register double tmp128_0 = A_21_6 + A_21_7;
                                        const register double tmp15_0 = A_12_2 + A_12_4;
                                        const register double tmp296_1 = tmp159_0*w42;
                                        const register double tmp130_1 = tmp67_0*w5;
                                        const register double tmp98_1 = A_01_6*w42;
                                        const register double tmp231_1 = tmp125_0*w6;
                                        const register double tmp42_1 = tmp34_0*w12;
                                        const register double tmp199_1 = A_02_5*w28;
                                        const register double tmp113_1 = tmp29_0*w13;
                                        const register double tmp330_1 = tmp152_0*w28;
                                        const register double tmp90_1 = A_01_1*w46;
                                        const register double tmp446_1 = tmp77_0*w22;
                                        const register double tmp108_1 = tmp43_0*w5;
                                        const register double tmp524_1 = A_12_6*w29;
                                        const register double tmp232_1 = tmp126_0*w34;
                                        const register double tmp33_1 = tmp25_0*w37;
                                        const register double tmp461_1 = tmp180_0*w1;
                                        const register double tmp14_1 = tmp8_0*w6;
                                        const register double tmp447_1 = tmp205_0*w26;
                                        const register double tmp452_1 = tmp198_0*w42;
                                        const register double tmp217_1 = tmp81_0*w9;
                                        const register double tmp76_1 = tmp59_0*w20;
                                        const register double tmp421_1 = tmp134_0*w31;
                                        const register double tmp485_1 = tmp51_0*w51;
                                        const register double tmp240_1 = tmp131_0*w1;
                                        const register double tmp160_1 = tmp91_0*w9;
                                        const register double tmp174_1 = A_20_1*w26;
                                        const register double tmp273_1 = A_10_1*w46;
                                        const register double tmp159_1 = tmp90_0*w47;
                                        const register double tmp228_1 = tmp103_0*w5;
                                        const register double tmp313_1 = tmp166_0*w45;
                                        const register double tmp45_1 = tmp37_0*w30;
                                        const register double tmp512_1 = tmp147_0*w13;
                                        const register double tmp73_1 = tmp56_0*w43;
                                        const register double tmp61_1 = A_01_6*w46;
                                        const register double tmp316_1 = tmp167_0*w43;
                                        const register double tmp189_1 = tmp112_0*w20;
                                        const register double tmp455_1 = tmp215_0*w39;
                                        const register double tmp360_1 = A_21_5*w24;
                                        const register double tmp258_1 = A_20_7*w2;
                                        const register double tmp196_1 = A_20_6*w26;
                                        const register double tmp37_1 = tmp29_0*w6;
                                        const register double tmp9_1 = A_12_7*w29;
                                        const register double tmp80_1 = tmp63_0*w19;
                                        const register double tmp312_1 = tmp165_0*w8;
                                        const register double tmp264_1 = tmp101_0*w1;
                                        const register double tmp124_1 = A_02_3*w26;
                                        const register double tmp229_1 = tmp123_0*w11;
                                        const register double tmp333_1 = tmp159_0*w46;
                                        const register double tmp533_1 = tmp222_0*w4;
                                        const register double tmp201_1 = tmp108_0*w37;
                                        const register double tmp444_1 = tmp35_0*w10;
                                        const register double tmp51_1 = tmp43_0*w18;
                                        const register double tmp214_1 = A_21_7*w29;
                                        const register double tmp518_1 = tmp86_0*w37;
                                        const register double tmp192_1 = tmp115_0*w5;
                                        const register double tmp355_1 = A_21_2*w27;
                                        const register double tmp156_1 = tmp87_0*w22;
                                        const register double tmp516_1 = tmp230_0*w27;
                                        const register double tmp366_1 = tmp104_0*w57;
                                        const register double tmp271_1 = tmp146_0*w49;
                                        const register double tmp437_1 = tmp218_0*w24;
                                        const register double tmp436_1 = tmp104_0*w54;
                                        const register double tmp167_1 = tmp98_0*w8;
                                        const register double tmp136_1 = tmp70_0*w34;
                                        const register double tmp406_1 = tmp207_0*w27;
                                        const register double tmp193_1 = tmp116_0*w12;
                                        const register double tmp486_1 = tmp225_0*w29;
                                        const register double tmp469_1 = tmp224_0*w11;
                                        const register double tmp287_1 = tmp71_0*w53;
                                        const register double tmp430_1 = tmp213_0*w28;
                                        const register double tmp462_1 = tmp220_0*w2;
                                        const register double tmp294_1 = tmp53_0*w59;
                                        const register double tmp218_1 = tmp118_0*w16;
                                        const register double tmp116_1 = tmp25_0*w31;
                                        const register double tmp495_1 = tmp76_0*w37;
                                        const register double tmp501_1 = tmp99_0*w46;
                                        const register double tmp0_1 = tmp0_0*w1;
                                        const register double tmp99_1 = tmp62_0*w17;
                                        const register double tmp429_1 = tmp212_0*w2;
                                        const register double tmp249_1 = tmp136_0*w9;
                                        const register double tmp504_1 = tmp229_0*w19;
                                        const register double tmp197_1 = A_12_2*w27;
                                        const register double tmp531_1 = tmp122_0*w35;
                                        const register double tmp265_1 = tmp142_0*w46;
                                        const register double tmp488_1 = tmp226_0*w4;
                                        const register double tmp528_1 = tmp115_0*w18;
                                        const register double tmp438_1 = tmp219_0*w2;
                                        const register double tmp233_1 = tmp127_0*w13;
                                        const register double tmp491_1 = tmp79_0*w1;
                                        const register double tmp215_1 = A_21_0*w4;
                                        const register double tmp24_1 = tmp18_0*w21;
                                        const register double tmp538_1 = tmp209_0*w27;
                                        const register double tmp379_1 = tmp167_0*w55;
                                        const register double tmp332_1 = tmp154_0*w4;
                                        const register double tmp498_1 = tmp68_0*w31;
                                        const register double tmp41_1 = tmp33_0*w33;
                                        const register double tmp464_1 = tmp179_0*w37;
                                        const register double tmp317_1 = tmp168_0*w40;
                                        const register double tmp378_1 = tmp106_0*w54;
                                        const register double tmp184_1 = tmp109_0*w14;
                                        const register double tmp292_1 = tmp14_0*w33;
                                        const register double tmp11_1 = tmp5_0*w11;
                                        const register double tmp354_1 = A_02_6*w26;
                                        const register double tmp84_1 = tmp37_0*w0;
                                        const register double tmp422_1 = tmp13_0*w30;
                                        const register double tmp132_1 = tmp69_0*w11;
                                        const register double tmp251_1 = tmp138_0*w31;
                                        const register double tmp18_1 = tmp12_0*w8;
                                        const register double tmp88_1 = A_21_1*w4;
                                        const register double tmp188_1 = A_12_2*w24;
                                        const register double tmp465_1 = tmp175_0*w31;
                                        const register double tmp235_1 = tmp128_0*w17;
                                        const register double tmp323_1 = A_02_1*w26;
                                        const register double tmp31_1 = tmp23_0*w38;
                                        const register double tmp397_1 = tmp170_0*w5;
                                        const register double tmp175_1 = tmp7_0*w3;
                                        const register double tmp148_1 = tmp81_0*w21;
                                        const register double tmp238_1 = tmp130_0*w19;
                                        const register double tmp59_1 = tmp46_0*w11;
                                        const register double tmp432_1 = tmp215_0*w35;
                                        const register double tmp398_1 = A_01_2*w46;
                                        const register double tmp497_1 = A_10_5*w46;
                                        const register double tmp28_1 = tmp21_0*w18;
                                        const register double tmp115_1 = tmp23_0*w32;
                                        const register double tmp441_1 = tmp23_0*w3;
                                        const register double tmp131_1 = tmp68_0*w37;
                                        const register double tmp289_1 = tmp155_0*w4;
                                        const register double tmp278_1 = tmp80_0*w44;
                                        const register double tmp5_1 = A_21_4*w27;
                                        const register double tmp254_1 = tmp140_0*w20;
                                        const register double tmp183_1 = tmp108_0*w31;
                                        const register double tmp279_1 = tmp151_0*w8;
                                        const register double tmp298_1 = tmp161_0*w16;
                                        const register double tmp505_1 = tmp230_0*w24;
                                        const register double tmp246_1 = tmp80_0*w52;
                                        const register double tmp100_1 = tmp53_0*w43;
                                        const register double tmp440_1 = tmp221_0*w16;
                                        const register double tmp481_1 = tmp188_0*w23;
                                        const register double tmp480_1 = tmp187_0*w35;
                                        const register double tmp384_1 = tmp150_0*w53;
                                        const register double tmp142_1 = tmp76_0*w31;
                                        const register double tmp372_1 = tmp191_0*w11;
                                        const register double tmp307_1 = A_10_7*w35;
                                        const register double tmp186_1 = tmp111_0*w19;
                                        const register double tmp127_1 = A_20_2*w2;
                                        const register double tmp391_1 = tmp167_0*w59;
                                        const register double tmp223_1 = tmp113_0*w20;
                                        const register double tmp454_1 = tmp197_0*w24;
                                        const register double tmp241_1 = tmp74_0*w51;
                                        const register double tmp529_1 = tmp114_0*w5;
                                        const register double tmp202_1 = tmp104_0*w7;
                                        const register double tmp236_1 = tmp96_0*w21;
                                        const register double tmp358_1 = tmp183_0*w11;
                                        const register double tmp102_1 = tmp51_0*w41;
                                        const register double tmp493_1 = A_20_5*w2;
                                        const register double tmp468_1 = tmp223_0*w4;
                                        const register double tmp435_1 = tmp217_0*w16;
                                        const register double tmp110_1 = tmp37_0*w36;
                                        const register double tmp479_1 = tmp189_0*w4;
                                        const register double tmp120_1 = tmp38_0*w22;
                                        const register double tmp16_1 = tmp10_0*w9;
                                        const register double tmp407_1 = tmp90_0*w53;
                                        const register double tmp442_1 = tmp66_0*w48;
                                        const register double tmp60_1 = A_10_4*w35;
                                        const register double tmp69_1 = tmp53_0*w45;
                                        const register double tmp144_1 = tmp77_0*w17;
                                        const register double tmp507_1 = tmp146_0*w48;
                                        const register double tmp424_1 = tmp174_0*w18;
                                        const register double tmp352_1 = tmp181_0*w23;
                                        const register double tmp451_1 = tmp199_0*w13;
                                        const register double tmp253_1 = tmp139_0*w16;
                                        const register double tmp353_1 = tmp182_0*w18;
                                        const register double tmp521_1 = tmp88_0*w22;
                                        const register double tmp346_1 = tmp175_0*w37;
                                        const register double tmp416_1 = tmp138_0*w37;
                                        const register double tmp324_1 = A_10_0*w35;
                                        const register double tmp152_1 = tmp84_0*w37;
                                        const register double tmp119_1 = tmp32_0*w21;
                                        const register double tmp86_1 = A_21_6*w29;
                                        const register double tmp290_1 = tmp156_0*w11;
                                        const register double tmp382_1 = tmp196_0*w26;
                                        const register double tmp91_1 = tmp49_0*w6;
                                        const register double tmp499_1 = A_10_2*w42;
                                        const register double tmp226_1 = tmp121_0*w13;
                                        const register double tmp477_1 = tmp195_0*w26;
                                        const register double tmp150_1 = A_02_4*w23;
                                        const register double tmp318_1 = tmp15_0*w22;
                                        const register double tmp396_1 = tmp206_0*w24;
                                        const register double tmp474_1 = A_02_0*w28;
                                        const register double tmp245_1 = tmp134_0*w37;
                                        const register double tmp3_1 = A_20_4*w26;
                                        const register double tmp44_1 = tmp36_0*w31;
                                        const register double tmp487_1 = tmp60_0*w52;
                                        const register double tmp293_1 = tmp158_0*w8;
                                        const register double tmp314_1 = A_01_2*w42;
                                        const register double tmp414_1 = tmp80_0*w51;
                                        const register double tmp472_1 = A_21_3*w27;
                                        const register double tmp321_1 = A_21_2*w24;
                                        const register double tmp225_1 = tmp120_0*w6;
                                        const register double tmp377_1 = tmp166_0*w59;
                                        const register double tmp413_1 = tmp186_0*w26;
                                        const register double tmp385_1 = tmp166_0*w55;
                                        const register double tmp310_1 = tmp164_0*w34;
                                        const register double tmp158_1 = tmp89_0*w34;
                                        const register double tmp449_1 = tmp203_0*w46;
                                        const register double tmp439_1 = tmp220_0*w28;
                                        const register double tmp22_1 = tmp16_0*w16;
                                        const register double tmp164_1 = tmp95_0*w46;
                                        const register double tmp417_1 = tmp74_0*w52;
                                        const register double tmp257_1 = tmp6_0*w25;
                                        const register double tmp203_1 = tmp18_0*w9;
                                        const register double tmp286_1 = tmp153_0*w28;
                                        const register double tmp155_1 = tmp33_0*w14;
                                        const register double tmp389_1 = tmp201_0*w12;
                                        const register double tmp508_1 = tmp145_0*w49;
                                        const register double tmp300_1 = tmp56_0*w55;
                                        const register double tmp299_1 = tmp162_0*w22;
                                        const register double tmp173_1 = tmp104_0*w25;
                                        const register double tmp32_1 = tmp24_0*w5;
                                        const register double tmp227_1 = tmp122_0*w39;
                                        const register double tmp484_1 = tmp3_0*w38;
                                        const register double tmp171_1 = tmp102_0*w21;
                                        const register double tmp478_1 = tmp190_0*w29;
                                        const register double tmp320_1 = tmp170_0*w18;
                                        const register double tmp327_1 = tmp6_0*w57;
                                        const register double tmp490_1 = tmp7_0*w32;
                                        const register double tmp419_1 = tmp127_0*w6;
                                        const register double tmp463_1 = tmp219_0*w28;
                                        const register double tmp12_1 = tmp6_0*w7;
                                        const register double tmp49_1 = tmp41_0*w22;
                                        const register double tmp344_1 = tmp173_0*w26;
                                        const register double tmp243_1 = tmp132_0*w2;
                                        const register double tmp83_1 = A_10_4*w39;
                                        const register double tmp297_1 = tmp160_0*w17;
                                        const register double tmp275_1 = tmp148_0*w34;
                                        const register double tmp168_1 = tmp99_0*w42;
                                        const register double tmp409_1 = tmp3_0*w32;
                                        const register double tmp1_1 = tmp1_0*w25;
                                        const register double tmp426_1 = tmp210_0*w39;
                                        const register double tmp375_1 = tmp109_0*w33;
                                        const register double tmp50_1 = tmp42_0*w19;
                                        const register double tmp513_1 = A_10_1*w42;
                                        const register double tmp97_1 = tmp45_0*w31;
                                        const register double tmp403_1 = tmp140_0*w1;
                                        const register double tmp71_1 = A_01_1*w42;
                                        const register double tmp520_1 = tmp84_0*w31;
                                        const register double tmp510_1 = A_10_6*w46;
                                        const register double tmp302_1 = A_10_0*w39;
                                        const register double tmp364_1 = tmp128_0*w22;
                                        const register double tmp515_1 = tmp142_0*w42;
                                        const register double tmp283_1 = tmp65_0*w56;
                                        const register double tmp222_1 = tmp112_0*w1;
                                        const register double tmp428_1 = tmp211_0*w27;
                                        const register double tmp371_1 = tmp190_0*w4;
                                        const register double tmp423_1 = tmp184_0*w23;
                                        const register double tmp276_1 = tmp149_0*w13;
                                        const register double tmp65_1 = tmp50_0*w9;
                                        const register double tmp305_1 = A_12_0*w29;
                                        const register double tmp170_1 = tmp101_0*w20;
                                        const register double tmp350_1 = tmp179_0*w31;
                                        const register double tmp466_1 = tmp172_0*w20;
                                        const register double tmp361_1 = tmp184_0*w26;
                                        const register double tmp431_1 = tmp214_0*w19;
                                        const register double tmp363_1 = tmp129_0*w17;
                                        const register double tmp178_1 = A_02_2*w28;
                                        const register double tmp527_1 = tmp120_0*w13;
                                        const register double tmp415_1 = tmp182_0*w5;
                                        const register double tmp450_1 = tmp200_0*w6;
                                        const register double tmp269_1 = A_01_7*w39;
                                        const register double tmp285_1 = tmp152_0*w2;
                                        const register double tmp272_1 = A_01_0*w35;
                                        const register double tmp339_1 = tmp136_0*w21;
                                        const register double tmp502_1 = tmp95_0*w42;
                                        const register double tmp38_1 = tmp30_0*w34;
                                        const register double tmp514_1 = tmp144_0*w46;
                                        const register double tmp96_1 = tmp56_0*w45;
                                        const register double tmp399_1 = tmp167_0*w45;
                                        const register double tmp483_1 = tmp173_0*w23;
                                        const register double tmp522_1 = tmp87_0*w17;
                                        const register double tmp519_1 = tmp91_0*w21;
                                        const register double tmp209_1 = A_12_5*w24;
                                        const register double tmp126_1 = tmp65_0*w48;
                                        const register double tmp367_1 = tmp187_0*w39;
                                        const register double tmp221_1 = tmp67_0*w18;
                                        const register double tmp381_1 = tmp146_0*w56;
                                        const register double tmp70_1 = tmp54_0*w31;
                                        const register double tmp216_1 = tmp117_0*w11;
                                        const register double tmp473_1 = A_02_7*w2;
                                        const register double tmp149_1 = tmp82_0*w22;
                                        const register double tmp357_1 = A_12_6*w4;
                                        const register double tmp534_1 = tmp226_0*w29;
                                        const register double tmp95_1 = tmp26_0*w15;
                                        const register double tmp500_1 = tmp64_0*w20;
                                        const register double tmp387_1 = tmp199_0*w6;
                                        const register double tmp471_1 = A_20_4*w23;
                                        const register double tmp281_1 = tmp74_0*w41;
                                        const register double tmp351_1 = tmp180_0*w20;
                                        const register double tmp63_1 = tmp48_0*w34;
                                        const register double tmp365_1 = tmp186_0*w23;
                                        const register double tmp448_1 = tmp206_0*w27;
                                        const register double tmp39_1 = tmp31_0*w13;
                                        const register double tmp453_1 = tmp196_0*w23;
                                        const register double tmp402_1 = tmp163_0*w18;
                                        const register double tmp137_1 = tmp71_0*w47;
                                        const register double tmp6_1 = A_02_0*w2;
                                        const register double tmp34_1 = tmp26_0*w36;
                                        const register double tmp383_1 = tmp197_0*w27;
                                        const register double tmp166_1 = tmp97_0*w12;
                                        const register double tmp114_1 = tmp40_0*w9;
                                        const register double tmp306_1 = A_12_7*w4;
                                        const register double tmp530_1 = tmp124_0*w39;
                                        const register double tmp388_1 = tmp200_0*w13;
                                        const register double tmp252_1 = tmp2_0*w30;
                                        const register double tmp210_1 = A_02_4*w26;
                                        const register double tmp200_1 = tmp21_0*w5;
                                        const register double tmp181_1 = tmp3_0*w10;
                                        const register double tmp425_1 = tmp106_0*w57;
                                        const register double tmp261_1 = A_21_7*w4;
                                        const register double tmp64_1 = tmp49_0*w13;
                                        const register double tmp506_1 = A_01_0*w39;
                                        const register double tmp457_1 = tmp213_0*w2;
                                        const register double tmp2_1 = tmp2_0*w0;
                                        const register double tmp393_1 = tmp203_0*w42;
                                        const register double tmp133_1 = A_01_3*w35;
                                        const register double tmp147_1 = tmp80_0*w41;
                                        const register double tmp8_1 = tmp4_0*w5;
                                        const register double tmp267_1 = tmp144_0*w42;
                                        const register double tmp17_1 = tmp11_0*w12;
                                        const register double tmp284_1 = tmp58_0*w50;
                                        const register double tmp328_1 = tmp66_0*w56;
                                        const register double tmp405_1 = tmp60_0*w51;
                                        const register double tmp467_1 = tmp222_0*w29;
                                        const register double tmp535_1 = tmp225_0*w4;
                                        const register double tmp356_1 = A_12_1*w29;
                                        const register double tmp274_1 = tmp147_0*w6;
                                        const register double tmp476_1 = tmp192_0*w39;
                                        const register double tmp206_1 = tmp10_0*w21;
                                        const register double tmp334_1 = tmp141_0*w9;
                                        const register double tmp482_1 = tmp181_0*w26;
                                        const register double tmp212_1 = A_20_7*w28;
                                        const register double tmp219_1 = tmp72_0*w21;
                                        const register double tmp47_1 = tmp39_0*w16;
                                        const register double tmp89_1 = A_10_3*w35;
                                        const register double tmp52_1 = tmp44_0*w1;
                                        const register double tmp492_1 = A_01_3*w39;
                                        const register double tmp81_1 = A_12_3*w24;
                                        const register double tmp77_1 = tmp60_0*w41;
                                        const register double tmp153_1 = tmp85_0*w21;
                                        const register double tmp304_1 = tmp163_0*w5;
                                        const register double tmp489_1 = tmp227_0*w11;
                                        const register double tmp107_1 = tmp35_0*w38;
                                        const register double tmp30_1 = tmp22_0*w39;
                                        const register double tmp260_1 = A_21_0*w29;
                                        const register double tmp343_1 = tmp172_0*w1;
                                        const register double tmp511_1 = tmp149_0*w6;
                                        const register double tmp139_1 = tmp73_0*w12;
                                        const register double tmp66_1 = tmp51_0*w44;
                                        const register double tmp208_1 = tmp4_0*w18;
                                        const register double tmp134_1 = tmp23_0*w10;
                                        const register double tmp205_1 = tmp105_0*w31;
                                        const register double tmp349_1 = tmp178_0*w8;
                                        const register double tmp341_1 = tmp53_0*w55;
                                        const register double tmp72_1 = tmp55_0*w17;
                                        const register double tmp79_1 = tmp62_0*w22;
                                        const register double tmp26_1 = tmp20_0*w19;
                                        const register double tmp141_1 = tmp75_0*w8;
                                        const register double tmp118_1 = tmp41_0*w17;
                                        const register double tmp259_1 = A_20_0*w28;
                                        const register double tmp458_1 = tmp212_0*w28;
                                        const register double tmp68_1 = tmp37_0*w15;
                                        const register double tmp154_1 = tmp86_0*w31;
                                        const register double tmp335_1 = tmp56_0*w59;
                                        const register double tmp359_1 = A_02_1*w23;
                                        const register double tmp56_1 = A_21_1*w29;
                                        const register double tmp392_1 = tmp145_0*w58;
                                        const register double tmp270_1 = tmp145_0*w48;
                                        const register double tmp92_1 = tmp47_0*w13;
                                        const register double tmp433_1 = tmp216_0*w34;
                                        const register double tmp420_1 = tmp125_0*w13;
                                        const register double tmp408_1 = tmp51_0*w52;
                                        const register double tmp494_1 = A_20_2*w28;
                                        const register double tmp362_1 = tmp185_0*w12;
                                        const register double tmp411_1 = tmp208_0*w19;
                                        const register double tmp336_1 = tmp65_0*w58;
                                        const register double tmp475_1 = A_21_4*w24;
                                        const register double tmp85_1 = A_12_3*w27;
                                        const register double tmp19_1 = tmp13_0*w15;
                                        const register double tmp537_1 = tmp132_0*w28;
                                        const register double tmp67_1 = tmp52_0*w8;
                                        const register double tmp459_1 = tmp210_0*w35;
                                        const register double tmp248_1 = tmp135_0*w34;
                                        const register double tmp326_1 = A_02_6*w23;
                                        const register double tmp23_1 = tmp17_0*w20;
                                        const register double tmp35_1 = tmp27_0*w11;
                                        const register double tmp62_1 = tmp47_0*w6;
                                        const register double tmp180_1 = tmp106_0*w7;
                                        const register double tmp277_1 = tmp150_0*w47;
                                        const register double tmp373_1 = tmp192_0*w35;
                                        const register double tmp337_1 = tmp157_0*w42;
                                        const register double tmp106_1 = tmp28_0*w39;
                                        const register double tmp369_1 = tmp168_0*w50;
                                        const register double tmp434_1 = tmp146_0*w58;
                                        const register double tmp331_1 = tmp155_0*w29;
                                        const register double tmp503_1 = tmp228_0*w27;
                                        const register double tmp93_1 = tmp61_0*w9;
                                        const register double tmp25_1 = tmp19_0*w22;
                                        const register double tmp146_1 = tmp79_0*w20;
                                        const register double tmp280_1 = A_10_6*w42;
                                        const register double tmp94_1 = tmp60_0*w44;
                                        const register double tmp400_1 = A_01_5*w42;
                                        const register double tmp151_1 = tmp83_0*w18;
                                        const register double tmp78_1 = tmp61_0*w21;
                                        const register double tmp301_1 = tmp6_0*w54;
                                        const register double tmp48_1 = tmp40_0*w21;
                                        const register double tmp75_1 = tmp58_0*w40;
                                        const register double tmp82_1 = tmp59_0*w1;
                                        const register double tmp74_1 = tmp57_0*w16;
                                        const register double tmp36_1 = tmp28_0*w35;
                                        const register double tmp370_1 = tmp189_0*w29;
                                        const register double tmp224_1 = tmp119_0*w19;
                                        const register double tmp109_1 = tmp36_0*w37;
                                        const register double tmp345_1 = tmp174_0*w5;
                                        const register double tmp101_1 = tmp44_0*w20;
                                        const register double tmp308_1 = A_01_5*w46;
                                        const register double tmp295_1 = tmp66_0*w58;
                                        const register double tmp117_1 = tmp26_0*w30;
                                        const register double tmp125_1 = tmp35_0*w3;
                                        const register double tmp309_1 = tmp9_0*w6;
                                        const register double tmp412_1 = tmp209_0*w24;
                                        const register double tmp46_1 = tmp38_0*w17;
                                        const register double tmp7_1 = A_02_7*w28;
                                        const register double tmp40_1 = tmp32_0*w9;
                                        const register double tmp386_1 = tmp198_0*w46;
                                        const register double tmp517_1 = tmp228_0*w24;
                                        const register double tmp532_1 = tmp223_0*w29;
                                        const register double tmp220_1 = A_02_3*w23;
                                        const register double tmp268_1 = tmp93_0*w20;
                                        const register double tmp322_1 = A_10_7*w39;
                                        const register double tmp311_1 = tmp8_0*w13;
                                        const register double tmp123_1 = A_01_4*w39;
                                        const register double tmp187_1 = A_20_6*w23;
                                        const register double tmp177_1 = A_02_5*w2;
                                        const register double tmp58_1 = A_21_6*w4;
                                        const register double tmp404_1 = tmp7_0*w38;
                                        const register double tmp122_1 = tmp64_0*w1;
                                        const register double tmp163_1 = tmp94_0*w5;
                                        const register double tmp15_1 = tmp9_0*w13;
                                        const register double tmp128_1 = A_20_5*w28;
                                        const register double tmp204_1 = tmp2_0*w15;
                                        const register double tmp539_1 = tmp207_0*w24;
                                        const register double tmp57_1 = tmp45_0*w37;
                                        const register double tmp53_1 = A_10_3*w39;
                                        const register double tmp157_1 = tmp88_0*w17;
                                        const register double tmp169_1 = tmp100_0*w16;
                                        const register double tmp162_1 = tmp93_0*w1;
                                        const register double tmp325_1 = tmp171_0*w12;
                                        const register double tmp179_1 = tmp105_0*w37;
                                        const register double tmp207_1 = A_20_1*w23;
                                        const register double tmp427_1 = tmp145_0*w56;
                                        const register double tmp368_1 = tmp188_0*w26;
                                        const register double tmp460_1 = tmp211_0*w24;
                                        const register double tmp347_1 = tmp176_0*w34;
                                        const register double tmp234_1 = tmp102_0*w9;
                                        const register double tmp21_1 = tmp15_0*w17;
                                        const register double tmp112_1 = tmp31_0*w6;
                                        const register double tmp319_1 = tmp169_0*w19;
                                        const register double tmp509_1 = A_01_7*w35;
                                        const register double tmp418_1 = tmp2_0*w36;
                                        const register double tmp266_1 = tmp143_0*w8;
                                        const register double tmp244_1 = tmp133_0*w28;
                                        const register double tmp138_1 = tmp72_0*w9;
                                        const register double tmp20_1 = tmp14_0*w14;
                                        const register double tmp29_1 = A_21_3*w24;
                                        const register double tmp190_1 = tmp113_0*w1;
                                        const register double tmp338_1 = tmp162_0*w17;
                                        const register double tmp87_1 = tmp54_0*w37;
                                        const register double tmp374_1 = tmp193_0*w34;
                                        const register double tmp195_1 = tmp13_0*w0;
                                        const register double tmp194_1 = tmp106_0*w25;
                                        const register double tmp111_1 = tmp22_0*w35;
                                        const register double tmp213_1 = tmp83_0*w5;
                                        const register double tmp291_1 = tmp157_0*w46;
                                        const register double tmp329_1 = tmp153_0*w2;
                                        const register double tmp256_1 = tmp17_0*w1;
                                        const register double tmp342_1 = tmp1_0*w54;
                                        const register double tmp376_1 = tmp194_0*w12;
                                        const register double tmp239_1 = tmp94_0*w18;
                                        const register double tmp340_1 = tmp160_0*w22;
                                        const register double tmp262_1 = tmp1_0*w7;
                                        const register double tmp54_1 = tmp26_0*w0;
                                        const register double tmp470_1 = A_20_3*w26;
                                        const register double tmp121_1 = tmp24_0*w18;
                                        const register double tmp523_1 = tmp85_0*w9;
                                        const register double tmp182_1 = tmp107_0*w12;
                                        const register double tmp140_1 = tmp74_0*w44;
                                        const register double tmp250_1 = tmp137_0*w8;
                                        const register double tmp104_1 = tmp55_0*w22;
                                        const register double tmp303_1 = A_21_5*w27;
                                        const register double tmp4_1 = tmp3_0*w3;
                                        const register double tmp526_1 = tmp121_0*w6;
                                        const register double tmp410_1 = tmp131_0*w20;
                                        const register double tmp255_1 = tmp141_0*w21;
                                        const register double tmp394_1 = tmp204_0*w19;
                                        const register double tmp172_1 = tmp103_0*w18;
                                        const register double tmp211_1 = A_20_0*w2;
                                        const register double tmp263_1 = tmp0_0*w20;
                                        const register double tmp230_1 = tmp124_0*w35;
                                        const register double tmp496_1 = A_01_4*w35;
                                        const register double tmp176_1 = A_12_5*w27;
                                        const register double tmp401_1 = tmp166_0*w43;
                                        const register double tmp43_1 = tmp35_0*w32;
                                        const register double tmp395_1 = tmp205_0*w23;
                                        const register double tmp348_1 = tmp177_0*w12;
                                        const register double tmp390_1 = tmp202_0*w8;
                                        const register double tmp525_1 = A_12_1*w4;
                                        const register double tmp237_1 = tmp129_0*w22;
                                        const register double tmp105_1 = A_12_4*w24;
                                        const register double tmp242_1 = tmp92_0*w50;
                                        const register double tmp161_1 = tmp92_0*w40;
                                        const register double tmp10_1 = A_12_0*w4;
                                        const register double tmp536_1 = tmp133_0*w2;
                                        const register double tmp55_1 = A_12_4*w27;
                                        const register double tmp445_1 = tmp82_0*w17;
                                        const register double tmp198_1 = A_02_2*w2;
                                        const register double tmp27_1 = A_20_3*w23;
                                        const register double tmp143_1 = A_10_5*w42;
                                        const register double tmp165_1 = tmp96_0*w9;
                                        const register double tmp145_1 = tmp78_0*w16;
                                        const register double tmp282_1 = tmp1_0*w57;
                                        const register double tmp13_1 = tmp7_0*w10;
                                        const register double tmp129_1 = tmp66_0*w49;
                                        const register double tmp443_1 = tmp65_0*w49;
                                        const register double tmp135_1 = A_10_2*w46;
                                        const register double tmp103_1 = tmp50_0*w21;
                                        const register double tmp247_1 = tmp13_0*w36;
                                        const register double tmp185_1 = tmp110_0*w16;
                                        const register double tmp191_1 = tmp114_0*w18;
                                        const register double tmp315_1 = tmp19_0*w17;
                                        const register double tmp380_1 = tmp195_0*w23;
                                        const register double tmp456_1 = tmp218_0*w27;
                                        const register double tmp288_1 = tmp154_0*w29;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=tmp0_1 + tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp1_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp30_1 + tmp31_1 + tmp32_1 + tmp33_1 + tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1 + tmp60_1 + tmp61_1 + tmp62_1 + tmp63_1 + tmp64_1 + tmp65_1 + tmp66_1 + tmp67_1 + tmp68_1 + tmp69_1 + tmp70_1 + tmp71_1 + tmp72_1 + tmp73_1 + tmp74_1 + tmp75_1 + tmp76_1 + tmp77_1 + tmp78_1 + tmp79_1 + tmp80_1 + tmp81_1;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=tmp100_1 + tmp101_1 + tmp102_1 + tmp103_1 + tmp104_1 + tmp105_1 + tmp59_1 + tmp63_1 + tmp67_1 + tmp74_1 + tmp75_1 + tmp80_1 + tmp82_1 + tmp83_1 + tmp84_1 + tmp85_1 + tmp86_1 + tmp87_1 + tmp88_1 + tmp89_1 + tmp90_1 + tmp91_1 + tmp92_1 + tmp93_1 + tmp94_1 + tmp95_1 + tmp96_1 + tmp97_1 + tmp98_1 + tmp99_1;
                                        EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp106_1 + tmp107_1 + tmp108_1 + tmp109_1 + tmp110_1 + tmp111_1 + tmp112_1 + tmp113_1 + tmp114_1 + tmp115_1 + tmp116_1 + tmp117_1 + tmp118_1 + tmp119_1 + tmp120_1 + tmp121_1 + tmp35_1 + tmp38_1 + tmp41_1 + tmp42_1 + tmp47_1 + tmp50_1;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=tmp122_1 + tmp123_1 + tmp124_1 + tmp125_1 + tmp126_1 + tmp127_1 + tmp128_1 + tmp129_1 + tmp130_1 + tmp131_1 + tmp132_1 + tmp133_1 + tmp134_1 + tmp135_1 + tmp136_1 + tmp137_1 + tmp138_1 + tmp139_1 + tmp140_1 + tmp141_1 + tmp142_1 + tmp143_1 + tmp144_1 + tmp145_1 + tmp146_1 + tmp147_1 + tmp148_1 + tmp149_1 + tmp150_1 + tmp151_1;
                                        EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp132_1 + tmp152_1 + tmp153_1 + tmp154_1 + tmp155_1 + tmp156_1 + tmp157_1 + tmp158_1 + tmp159_1 + tmp160_1 + tmp161_1 + tmp74_1;
                                        EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp162_1 + tmp163_1 + tmp164_1 + tmp165_1 + tmp166_1 + tmp167_1 + tmp168_1 + tmp169_1 + tmp170_1 + tmp171_1 + tmp172_1 + tmp31_1 + tmp34_1 + tmp35_1 + tmp37_1 + tmp39_1 + tmp41_1 + tmp43_1 + tmp45_1 + tmp46_1 + tmp49_1 + tmp50_1;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=tmp104_1 + tmp158_1 + tmp16_1 + tmp173_1 + tmp174_1 + tmp175_1 + tmp176_1 + tmp177_1 + tmp178_1 + tmp179_1 + tmp180_1 + tmp181_1 + tmp182_1 + tmp183_1 + tmp184_1 + tmp185_1 + tmp186_1 + tmp187_1 + tmp188_1 + tmp19_1 + tmp24_1 + tmp28_1 + tmp2_1 + tmp59_1 + tmp86_1 + tmp88_1 + tmp8_1 + tmp91_1 + tmp92_1 + tmp99_1;
                                        EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp132_1 + tmp155_1 + tmp156_1 + tmp157_1 + tmp159_1 + tmp161_1 + tmp189_1 + tmp18_1 + tmp190_1 + tmp191_1 + tmp192_1 + tmp193_1;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=tmp13_1 + tmp158_1 + tmp182_1 + tmp184_1 + tmp185_1 + tmp186_1 + tmp194_1 + tmp195_1 + tmp196_1 + tmp197_1 + tmp198_1 + tmp199_1 + tmp200_1 + tmp201_1 + tmp202_1 + tmp203_1 + tmp204_1 + tmp205_1 + tmp206_1 + tmp207_1 + tmp208_1 + tmp209_1 + tmp4_1 + tmp56_1 + tmp58_1 + tmp59_1 + tmp62_1 + tmp64_1 + tmp72_1 + tmp79_1;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=tmp0_1 + tmp12_1 + tmp139_1 + tmp13_1 + tmp18_1 + tmp19_1 + tmp1_1 + tmp20_1 + tmp210_1 + tmp211_1 + tmp212_1 + tmp213_1 + tmp214_1 + tmp215_1 + tmp216_1 + tmp217_1 + tmp218_1 + tmp219_1 + tmp220_1 + tmp221_1 + tmp23_1 + tmp2_1 + tmp4_1 + tmp55_1 + tmp62_1 + tmp64_1 + tmp72_1 + tmp79_1 + tmp80_1 + tmp81_1;
                                        EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp153_1 + tmp155_1 + tmp159_1 + tmp160_1 + tmp161_1 + tmp18_1 + tmp222_1 + tmp223_1 + tmp224_1 + tmp225_1 + tmp226_1 + tmp74_1;
                                        EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp107_1 + tmp109_1 + tmp110_1 + tmp115_1 + tmp116_1 + tmp117_1 + tmp166_1 + tmp169_1 + tmp227_1 + tmp228_1 + tmp229_1 + tmp230_1 + tmp231_1 + tmp232_1 + tmp233_1 + tmp234_1 + tmp235_1 + tmp236_1 + tmp237_1 + tmp238_1 + tmp239_1 + tmp41_1;
                                        EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp112_1 + tmp113_1 + tmp240_1 + tmp241_1 + tmp242_1 + tmp243_1 + tmp244_1 + tmp245_1 + tmp246_1 + tmp247_1 + tmp248_1 + tmp249_1 + tmp250_1 + tmp251_1 + tmp252_1 + tmp253_1 + tmp254_1 + tmp255_1 + tmp35_1 + tmp46_1 + tmp49_1 + tmp50_1;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=tmp104_1 + tmp105_1 + tmp124_1 + tmp130_1 + tmp138_1 + tmp139_1 + tmp148_1 + tmp150_1 + tmp151_1 + tmp175_1 + tmp181_1 + tmp18_1 + tmp195_1 + tmp204_1 + tmp20_1 + tmp216_1 + tmp218_1 + tmp256_1 + tmp257_1 + tmp258_1 + tmp259_1 + tmp260_1 + tmp261_1 + tmp262_1 + tmp263_1 + tmp80_1 + tmp85_1 + tmp91_1 + tmp92_1 + tmp99_1;
                                        EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp107_1 + tmp108_1 + tmp110_1 + tmp114_1 + tmp115_1 + tmp117_1 + tmp119_1 + tmp121_1 + tmp229_1 + tmp231_1 + tmp233_1 + tmp235_1 + tmp237_1 + tmp238_1 + tmp264_1 + tmp265_1 + tmp266_1 + tmp267_1 + tmp268_1 + tmp41_1 + tmp42_1 + tmp47_1;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=tmp122_1 + tmp125_1 + tmp131_1 + tmp134_1 + tmp142_1 + tmp146_1 + tmp16_1 + tmp174_1 + tmp182_1 + tmp187_1 + tmp224_1 + tmp22_1 + tmp24_1 + tmp269_1 + tmp270_1 + tmp271_1 + tmp272_1 + tmp273_1 + tmp274_1 + tmp275_1 + tmp276_1 + tmp277_1 + tmp278_1 + tmp279_1 + tmp280_1 + tmp281_1 + tmp28_1 + tmp6_1 + tmp7_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=tmp249_1 + tmp255_1 + tmp264_1 + tmp268_1 + tmp282_1 + tmp283_1 + tmp284_1 + tmp285_1 + tmp286_1 + tmp287_1 + tmp288_1 + tmp289_1 + tmp290_1 + tmp291_1 + tmp292_1 + tmp293_1 + tmp294_1 + tmp295_1 + tmp296_1 + tmp297_1 + tmp298_1 + tmp299_1 + tmp300_1 + tmp301_1;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=tmp102_1 + tmp11_1 + tmp193_1 + tmp302_1 + tmp303_1 + tmp304_1 + tmp305_1 + tmp306_1 + tmp307_1 + tmp308_1 + tmp309_1 + tmp310_1 + tmp311_1 + tmp312_1 + tmp313_1 + tmp314_1 + tmp315_1 + tmp316_1 + tmp317_1 + tmp318_1 + tmp319_1 + tmp320_1 + tmp321_1 + tmp52_1 + tmp54_1 + tmp57_1 + tmp68_1 + tmp70_1 + tmp76_1 + tmp94_1;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=tmp101_1 + tmp125_1 + tmp134_1 + tmp211_1 + tmp212_1 + tmp213_1 + tmp217_1 + tmp218_1 + tmp219_1 + tmp221_1 + tmp224_1 + tmp270_1 + tmp271_1 + tmp274_1 + tmp276_1 + tmp277_1 + tmp278_1 + tmp281_1 + tmp310_1 + tmp322_1 + tmp323_1 + tmp324_1 + tmp325_1 + tmp326_1 + tmp67_1 + tmp82_1 + tmp87_1 + tmp90_1 + tmp97_1 + tmp98_1;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=tmp162_1 + tmp170_1 + tmp284_1 + tmp287_1 + tmp290_1 + tmp292_1 + tmp293_1 + tmp298_1 + tmp327_1 + tmp328_1 + tmp329_1 + tmp330_1 + tmp331_1 + tmp332_1 + tmp333_1 + tmp334_1 + tmp335_1 + tmp336_1 + tmp337_1 + tmp338_1 + tmp339_1 + tmp340_1 + tmp341_1 + tmp342_1;
                                        EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp112_1 + tmp113_1 + tmp241_1 + tmp242_1 + tmp246_1 + tmp247_1 + tmp252_1 + tmp343_1 + tmp344_1 + tmp345_1 + tmp346_1 + tmp347_1 + tmp348_1 + tmp349_1 + tmp350_1 + tmp351_1 + tmp352_1 + tmp353_1 + tmp35_1 + tmp46_1 + tmp49_1 + tmp50_1;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=tmp127_1 + tmp128_1 + tmp130_1 + tmp138_1 + tmp13_1 + tmp145_1 + tmp148_1 + tmp14_1 + tmp151_1 + tmp158_1 + tmp15_1 + tmp184_1 + tmp194_1 + tmp195_1 + tmp201_1 + tmp202_1 + tmp204_1 + tmp205_1 + tmp21_1 + tmp25_1 + tmp319_1 + tmp325_1 + tmp354_1 + tmp355_1 + tmp356_1 + tmp357_1 + tmp358_1 + tmp359_1 + tmp360_1 + tmp4_1;
                                        EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp229_1 + tmp231_1 + tmp233_1 + tmp238_1 + tmp240_1 + tmp241_1 + tmp242_1 + tmp245_1 + tmp246_1 + tmp247_1 + tmp248_1 + tmp250_1 + tmp251_1 + tmp252_1 + tmp254_1 + tmp345_1 + tmp353_1 + tmp361_1 + tmp362_1 + tmp363_1 + tmp364_1 + tmp365_1;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=tmp109_1 + tmp116_1 + tmp283_1 + tmp287_1 + tmp295_1 + tmp338_1 + tmp340_1 + tmp345_1 + tmp353_1 + tmp366_1 + tmp367_1 + tmp368_1 + tmp369_1 + tmp370_1 + tmp371_1 + tmp372_1 + tmp373_1 + tmp374_1 + tmp375_1 + tmp376_1 + tmp377_1 + tmp378_1 + tmp379_1 + tmp380_1;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=tmp162_1 + tmp170_1 + tmp282_1 + tmp292_1 + tmp301_1 + tmp345_1 + tmp353_1 + tmp369_1 + tmp381_1 + tmp382_1 + tmp383_1 + tmp384_1 + tmp385_1 + tmp386_1 + tmp387_1 + tmp388_1 + tmp389_1 + tmp390_1 + tmp391_1 + tmp392_1 + tmp393_1 + tmp394_1 + tmp395_1 + tmp396_1;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=tmp101_1 + tmp10_1 + tmp11_1 + tmp14_1 + tmp15_1 + tmp193_1 + tmp21_1 + tmp25_1 + tmp310_1 + tmp312_1 + tmp317_1 + tmp319_1 + tmp322_1 + tmp324_1 + tmp355_1 + tmp360_1 + tmp397_1 + tmp398_1 + tmp399_1 + tmp400_1 + tmp401_1 + tmp402_1 + tmp66_1 + tmp77_1 + tmp82_1 + tmp84_1 + tmp87_1 + tmp95_1 + tmp97_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp108_1 + tmp121_1 + tmp245_1 + tmp248_1 + tmp250_1 + tmp251_1 + tmp387_1 + tmp388_1 + tmp403_1 + tmp404_1 + tmp405_1 + tmp406_1 + tmp407_1 + tmp408_1 + tmp409_1 + tmp40_1 + tmp410_1 + tmp411_1 + tmp412_1 + tmp42_1 + tmp47_1 + tmp48_1;
                                        EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp229_1 + tmp235_1 + tmp237_1 + tmp238_1 + tmp242_1 + tmp248_1 + tmp250_1 + tmp362_1 + tmp403_1 + tmp410_1 + tmp413_1 + tmp414_1 + tmp415_1 + tmp416_1 + tmp417_1 + tmp418_1 + tmp419_1 + tmp420_1 + tmp421_1 + tmp422_1 + tmp423_1 + tmp424_1;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=tmp109_1 + tmp116_1 + tmp284_1 + tmp294_1 + tmp300_1 + tmp334_1 + tmp339_1 + tmp375_1 + tmp384_1 + tmp387_1 + tmp388_1 + tmp425_1 + tmp426_1 + tmp427_1 + tmp428_1 + tmp429_1 + tmp430_1 + tmp431_1 + tmp432_1 + tmp433_1 + tmp434_1 + tmp435_1 + tmp436_1 + tmp437_1;
                                        EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp229_1 + tmp231_1 + tmp233_1 + tmp238_1 + tmp241_1 + tmp242_1 + tmp246_1 + tmp247_1 + tmp249_1 + tmp252_1 + tmp255_1 + tmp343_1 + tmp346_1 + tmp347_1 + tmp349_1 + tmp350_1 + tmp351_1 + tmp363_1 + tmp364_1 + tmp438_1 + tmp439_1 + tmp440_1;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=tmp132_1 + tmp137_1 + tmp16_1 + tmp177_1 + tmp178_1 + tmp17_1 + tmp185_1 + tmp24_1 + tmp278_1 + tmp27_1 + tmp281_1 + tmp28_1 + tmp308_1 + tmp312_1 + tmp314_1 + tmp3_1 + tmp441_1 + tmp442_1 + tmp443_1 + tmp444_1 + tmp445_1 + tmp446_1 + tmp52_1 + tmp53_1 + tmp57_1 + tmp60_1 + tmp63_1 + tmp70_1 + tmp76_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=tmp264_1 + tmp268_1 + tmp292_1 + tmp327_1 + tmp342_1 + tmp369_1 + tmp377_1 + tmp379_1 + tmp384_1 + tmp389_1 + tmp390_1 + tmp394_1 + tmp415_1 + tmp424_1 + tmp427_1 + tmp434_1 + tmp447_1 + tmp448_1 + tmp449_1 + tmp450_1 + tmp451_1 + tmp452_1 + tmp453_1 + tmp454_1;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=tmp249_1 + tmp255_1 + tmp284_1 + tmp335_1 + tmp33_1 + tmp341_1 + tmp366_1 + tmp375_1 + tmp378_1 + tmp381_1 + tmp384_1 + tmp392_1 + tmp431_1 + tmp433_1 + tmp435_1 + tmp44_1 + tmp450_1 + tmp451_1 + tmp455_1 + tmp456_1 + tmp457_1 + tmp458_1 + tmp459_1 + tmp460_1;
                                        EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp229_1 + tmp235_1 + tmp237_1 + tmp238_1 + tmp242_1 + tmp334_1 + tmp339_1 + tmp347_1 + tmp349_1 + tmp414_1 + tmp417_1 + tmp418_1 + tmp419_1 + tmp420_1 + tmp422_1 + tmp440_1 + tmp461_1 + tmp462_1 + tmp463_1 + tmp464_1 + tmp465_1 + tmp466_1;
                                        EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp108_1 + tmp121_1 + tmp297_1 + tmp299_1 + tmp346_1 + tmp347_1 + tmp349_1 + tmp350_1 + tmp404_1 + tmp405_1 + tmp407_1 + tmp408_1 + tmp409_1 + tmp40_1 + tmp42_1 + tmp461_1 + tmp466_1 + tmp467_1 + tmp468_1 + tmp469_1 + tmp47_1 + tmp48_1;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=tmp101_1 + tmp125_1 + tmp126_1 + tmp129_1 + tmp132_1 + tmp134_1 + tmp137_1 + tmp140_1 + tmp144_1 + tmp147_1 + tmp149_1 + tmp17_1 + tmp185_1 + tmp198_1 + tmp199_1 + tmp200_1 + tmp203_1 + tmp206_1 + tmp208_1 + tmp312_1 + tmp398_1 + tmp400_1 + tmp470_1 + tmp471_1 + tmp63_1 + tmp82_1 + tmp83_1 + tmp87_1 + tmp89_1 + tmp97_1;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=tmp11_1 + tmp175_1 + tmp17_1 + tmp181_1 + tmp18_1 + tmp195_1 + tmp200_1 + tmp203_1 + tmp204_1 + tmp206_1 + tmp208_1 + tmp20_1 + tmp22_1 + tmp256_1 + tmp257_1 + tmp262_1 + tmp263_1 + tmp26_1 + tmp305_1 + tmp306_1 + tmp309_1 + tmp311_1 + tmp315_1 + tmp318_1 + tmp470_1 + tmp471_1 + tmp472_1 + tmp473_1 + tmp474_1 + tmp475_1;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=tmp287_1 + tmp297_1 + tmp299_1 + tmp328_1 + tmp336_1 + tmp33_1 + tmp369_1 + tmp372_1 + tmp374_1 + tmp375_1 + tmp376_1 + tmp385_1 + tmp391_1 + tmp415_1 + tmp424_1 + tmp425_1 + tmp436_1 + tmp44_1 + tmp476_1 + tmp477_1 + tmp478_1 + tmp479_1 + tmp480_1 + tmp481_1;
                                        EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp118_1 + tmp120_1 + tmp242_1 + tmp347_1 + tmp348_1 + tmp349_1 + tmp35_1 + tmp37_1 + tmp39_1 + tmp414_1 + tmp415_1 + tmp417_1 + tmp418_1 + tmp422_1 + tmp424_1 + tmp461_1 + tmp464_1 + tmp465_1 + tmp466_1 + tmp482_1 + tmp483_1 + tmp50_1;
                                        EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp163_1 + tmp166_1 + tmp169_1 + tmp172_1 + tmp234_1 + tmp236_1 + tmp240_1 + tmp248_1 + tmp250_1 + tmp254_1 + tmp338_1 + tmp340_1 + tmp407_1 + tmp416_1 + tmp421_1 + tmp484_1 + tmp485_1 + tmp486_1 + tmp487_1 + tmp488_1 + tmp489_1 + tmp490_1;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=tmp122_1 + tmp131_1 + tmp135_1 + tmp141_1 + tmp142_1 + tmp143_1 + tmp146_1 + tmp186_1 + tmp193_1 + tmp197_1 + tmp209_1 + tmp214_1 + tmp215_1 + tmp216_1 + tmp269_1 + tmp272_1 + tmp275_1 + tmp317_1 + tmp397_1 + tmp399_1 + tmp401_1 + tmp402_1 + tmp62_1 + tmp64_1 + tmp66_1 + tmp72_1 + tmp77_1 + tmp79_1 + tmp84_1 + tmp95_1;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=tmp132_1 + tmp136_1 + tmp137_1 + tmp139_1 + tmp141_1 + tmp145_1 + tmp210_1 + tmp213_1 + tmp217_1 + tmp219_1 + tmp220_1 + tmp221_1 + tmp278_1 + tmp281_1 + tmp441_1 + tmp442_1 + tmp443_1 + tmp444_1 + tmp445_1 + tmp446_1 + tmp491_1 + tmp492_1 + tmp493_1 + tmp494_1 + tmp495_1 + tmp496_1 + tmp497_1 + tmp498_1 + tmp499_1 + tmp500_1;
                                        EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp107_1 + tmp110_1 + tmp112_1 + tmp113_1 + tmp115_1 + tmp117_1 + tmp118_1 + tmp120_1 + tmp166_1 + tmp167_1 + tmp169_1 + tmp228_1 + tmp234_1 + tmp236_1 + tmp239_1 + tmp264_1 + tmp268_1 + tmp35_1 + tmp41_1 + tmp501_1 + tmp502_1 + tmp50_1;
                                        EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp163_1 + tmp166_1 + tmp169_1 + tmp172_1 + tmp234_1 + tmp236_1 + tmp343_1 + tmp347_1 + tmp349_1 + tmp351_1 + tmp407_1 + tmp450_1 + tmp451_1 + tmp464_1 + tmp465_1 + tmp484_1 + tmp485_1 + tmp487_1 + tmp490_1 + tmp503_1 + tmp504_1 + tmp505_1;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=tmp140_1 + tmp147_1 + tmp182_1 + tmp196_1 + tmp200_1 + tmp203_1 + tmp206_1 + tmp207_1 + tmp208_1 + tmp224_1 + tmp22_1 + tmp275_1 + tmp277_1 + tmp279_1 + tmp441_1 + tmp444_1 + tmp473_1 + tmp474_1 + tmp491_1 + tmp495_1 + tmp498_1 + tmp500_1 + tmp506_1 + tmp507_1 + tmp508_1 + tmp509_1 + tmp510_1 + tmp511_1 + tmp512_1 + tmp513_1;
                                        EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp162_1 + tmp170_1 + tmp229_1 + tmp238_1 + tmp266_1 + tmp31_1 + tmp32_1 + tmp34_1 + tmp363_1 + tmp364_1 + tmp40_1 + tmp419_1 + tmp41_1 + tmp420_1 + tmp42_1 + tmp43_1 + tmp45_1 + tmp47_1 + tmp48_1 + tmp514_1 + tmp515_1 + tmp51_1;
                                        EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp165_1 + tmp166_1 + tmp169_1 + tmp171_1 + tmp228_1 + tmp239_1 + tmp346_1 + tmp347_1 + tmp349_1 + tmp350_1 + tmp387_1 + tmp388_1 + tmp404_1 + tmp405_1 + tmp407_1 + tmp408_1 + tmp409_1 + tmp461_1 + tmp466_1 + tmp504_1 + tmp516_1 + tmp517_1;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=tmp130_1 + tmp138_1 + tmp140_1 + tmp147_1 + tmp148_1 + tmp151_1 + tmp218_1 + tmp224_1 + tmp258_1 + tmp259_1 + tmp277_1 + tmp302_1 + tmp307_1 + tmp310_1 + tmp325_1 + tmp354_1 + tmp359_1 + tmp441_1 + tmp444_1 + tmp507_1 + tmp508_1 + tmp511_1 + tmp512_1 + tmp52_1 + tmp57_1 + tmp61_1 + tmp67_1 + tmp70_1 + tmp71_1 + tmp76_1;
                                        EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp132_1 + tmp155_1 + tmp158_1 + tmp159_1 + tmp161_1 + tmp518_1 + tmp519_1 + tmp520_1 + tmp521_1 + tmp522_1 + tmp523_1 + tmp74_1;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=tmp100_1 + tmp102_1 + tmp103_1 + tmp122_1 + tmp123_1 + tmp131_1 + tmp133_1 + tmp136_1 + tmp142_1 + tmp146_1 + tmp26_1 + tmp273_1 + tmp279_1 + tmp280_1 + tmp309_1 + tmp311_1 + tmp315_1 + tmp318_1 + tmp358_1 + tmp472_1 + tmp475_1 + tmp524_1 + tmp525_1 + tmp74_1 + tmp75_1 + tmp84_1 + tmp93_1 + tmp94_1 + tmp95_1 + tmp96_1;
                                        EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp155_1 + tmp159_1 + tmp161_1 + tmp189_1 + tmp18_1 + tmp190_1 + tmp224_1 + tmp519_1 + tmp523_1 + tmp526_1 + tmp527_1 + tmp74_1;
                                        EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp132_1 + tmp155_1 + tmp159_1 + tmp161_1 + tmp18_1 + tmp193_1 + tmp222_1 + tmp223_1 + tmp521_1 + tmp522_1 + tmp528_1 + tmp529_1;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=tmp136_1 + tmp14_1 + tmp15_1 + tmp21_1 + tmp25_1 + tmp26_1 + tmp279_1 + tmp29_1 + tmp356_1 + tmp357_1 + tmp358_1 + tmp491_1 + tmp492_1 + tmp495_1 + tmp496_1 + tmp498_1 + tmp500_1 + tmp510_1 + tmp513_1 + tmp54_1 + tmp5_1 + tmp65_1 + tmp66_1 + tmp68_1 + tmp69_1 + tmp73_1 + tmp74_1 + tmp75_1 + tmp77_1 + tmp78_1;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=tmp102_1 + tmp104_1 + tmp141_1 + tmp176_1 + tmp186_1 + tmp188_1 + tmp193_1 + tmp216_1 + tmp260_1 + tmp261_1 + tmp275_1 + tmp304_1 + tmp313_1 + tmp316_1 + tmp317_1 + tmp320_1 + tmp491_1 + tmp495_1 + tmp497_1 + tmp498_1 + tmp499_1 + tmp500_1 + tmp506_1 + tmp509_1 + tmp54_1 + tmp68_1 + tmp91_1 + tmp92_1 + tmp94_1 + tmp99_1;
                                        EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp163_1 + tmp165_1 + tmp166_1 + tmp169_1 + tmp171_1 + tmp172_1 + tmp229_1 + tmp232_1 + tmp238_1 + tmp31_1 + tmp33_1 + tmp34_1 + tmp363_1 + tmp364_1 + tmp419_1 + tmp41_1 + tmp420_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp530_1 + tmp531_1;
                                        EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp114_1 + tmp119_1 + tmp32_1 + tmp338_1 + tmp340_1 + tmp343_1 + tmp347_1 + tmp349_1 + tmp351_1 + tmp407_1 + tmp42_1 + tmp464_1 + tmp465_1 + tmp469_1 + tmp47_1 + tmp484_1 + tmp485_1 + tmp487_1 + tmp490_1 + tmp51_1 + tmp532_1 + tmp533_1;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=tmp145_1 + tmp158_1 + tmp173_1 + tmp175_1 + tmp179_1 + tmp180_1 + tmp181_1 + tmp183_1 + tmp184_1 + tmp19_1 + tmp213_1 + tmp217_1 + tmp219_1 + tmp221_1 + tmp2_1 + tmp303_1 + tmp309_1 + tmp311_1 + tmp315_1 + tmp318_1 + tmp319_1 + tmp321_1 + tmp323_1 + tmp325_1 + tmp326_1 + tmp358_1 + tmp493_1 + tmp494_1 + tmp524_1 + tmp525_1;
                                        EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp152_1 + tmp154_1 + tmp155_1 + tmp158_1 + tmp159_1 + tmp161_1 + tmp193_1 + tmp224_1 + tmp526_1 + tmp527_1 + tmp528_1 + tmp529_1;
                                        EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp165_1 + tmp166_1 + tmp169_1 + tmp171_1 + tmp228_1 + tmp239_1 + tmp245_1 + tmp248_1 + tmp250_1 + tmp251_1 + tmp297_1 + tmp299_1 + tmp403_1 + tmp404_1 + tmp405_1 + tmp407_1 + tmp408_1 + tmp409_1 + tmp410_1 + tmp489_1 + tmp534_1 + tmp535_1;
                                        EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp118_1 + tmp120_1 + tmp242_1 + tmp248_1 + tmp250_1 + tmp253_1 + tmp334_1 + tmp339_1 + tmp35_1 + tmp37_1 + tmp39_1 + tmp403_1 + tmp410_1 + tmp414_1 + tmp416_1 + tmp417_1 + tmp418_1 + tmp421_1 + tmp422_1 + tmp50_1 + tmp536_1 + tmp537_1;
                                        EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp155_1 + tmp158_1 + tmp159_1 + tmp161_1 + tmp191_1 + tmp192_1 + tmp193_1 + tmp224_1 + tmp225_1 + tmp226_1 + tmp518_1 + tmp520_1;
                                        EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp114_1 + tmp119_1 + tmp240_1 + tmp248_1 + tmp250_1 + tmp254_1 + tmp32_1 + tmp407_1 + tmp411_1 + tmp416_1 + tmp421_1 + tmp42_1 + tmp450_1 + tmp451_1 + tmp47_1 + tmp484_1 + tmp485_1 + tmp487_1 + tmp490_1 + tmp51_1 + tmp538_1 + tmp539_1;
                                    }
                                }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const register double A_00 = A_p[INDEX4(k,0,m,0, numEq,3, numComp)];
                                        const register double A_01 = A_p[INDEX4(k,0,m,1, numEq,3, numComp)];
                                        const register double A_02 = A_p[INDEX4(k,0,m,2, numEq,3, numComp)];
                                        const register double A_10 = A_p[INDEX4(k,1,m,0, numEq,3, numComp)];
                                        const register double A_11 = A_p[INDEX4(k,1,m,1, numEq,3, numComp)];
                                        const register double A_12 = A_p[INDEX4(k,1,m,2, numEq,3, numComp)];
                                        const register double A_20 = A_p[INDEX4(k,2,m,0, numEq,3, numComp)];
                                        const register double A_21 = A_p[INDEX4(k,2,m,1, numEq,3, numComp)];
                                        const register double A_22 = A_p[INDEX4(k,2,m,2, numEq,3, numComp)];
                                        const register double tmp0_0 = A_01 + A_10;
                                        const register double tmp1_0 = A_02 + A_20;
                                        const register double tmp2_0 = A_12 + A_21;
                                        const register double tmp25_1 = A_01*w69;
                                        const register double tmp2_1 = tmp0_0*w61;
                                        const register double tmp33_1 = A_20*w70;
                                        const register double tmp23_1 = A_02*w65;
                                        const register double tmp41_1 = A_01*w61;
                                        const register double tmp34_1 = A_02*w73;
                                        const register double tmp8_1 = A_11*w71;
                                        const register double tmp50_1 = A_10*w61;
                                        const register double tmp15_1 = A_22*w75;
                                        const register double tmp9_1 = A_21*w74;
                                        const register double tmp19_1 = A_10*w69;
                                        const register double tmp11_1 = A_00*w68;
                                        const register double tmp52_1 = tmp2_0*w66;
                                        const register double tmp37_1 = tmp2_0*w74;
                                        const register double tmp0_1 = A_00*w60;
                                        const register double tmp17_1 = A_21*w64;
                                        const register double tmp26_1 = A_00*w79;
                                        const register double tmp5_1 = A_21*w66;
                                        const register double tmp29_1 = A_00*w80;
                                        const register double tmp7_1 = A_22*w67;
                                        const register double tmp48_1 = A_11*w87;
                                        const register double tmp44_1 = A_11*w84;
                                        const register double tmp27_1 = tmp2_0*w72;
                                        const register double tmp42_1 = A_22*w85;
                                        const register double tmp18_1 = A_11*w77;
                                        const register double tmp35_1 = tmp0_0*w76;
                                        const register double tmp46_1 = A_00*w86;
                                        const register double tmp32_1 = A_22*w83;
                                        const register double tmp22_1 = A_01*w76;
                                        const register double tmp4_1 = A_02*w62;
                                        const register double tmp10_1 = A_02*w70;
                                        const register double tmp3_1 = A_20*w65;
                                        const register double tmp39_1 = A_21*w72;
                                        const register double tmp51_1 = tmp1_0*w65;
                                        const register double tmp12_1 = A_20*w73;
                                        const register double tmp40_1 = A_10*w81;
                                        const register double tmp43_1 = tmp1_0*w62;
                                        const register double tmp28_1 = A_10*w76;
                                        const register double tmp45_1 = tmp2_0*w64;
                                        const register double tmp49_1 = A_01*w81;
                                        const register double tmp36_1 = tmp1_0*w73;
                                        const register double tmp53_1 = A_00*w89;
                                        const register double tmp6_1 = A_11*w63;
                                        const register double tmp31_1 = A_11*w82;
                                        const register double tmp21_1 = A_22*w78;
                                        const register double tmp16_1 = tmp1_0*w70;
                                        const register double tmp14_1 = A_12*w72;
                                        const register double tmp38_1 = A_12*w74;
                                        const register double tmp30_1 = tmp0_0*w81;
                                        const register double tmp47_1 = A_22*w88;
                                        const register double tmp13_1 = tmp0_0*w69;
                                        const register double tmp20_1 = A_12*w66;
                                        const register double tmp1_1 = A_12*w64;
                                        const register double tmp24_1 = A_20*w62;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp8_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=tmp0_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=tmp0_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1;
                                        EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp8_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=tmp21_1 + tmp23_1 + tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp16_1 + tmp27_1 + tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1;
                                        EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp11_1 + tmp14_1 + tmp15_1 + tmp33_1 + tmp34_1 + tmp35_1 + tmp8_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=tmp0_1 + tmp17_1 + tmp20_1 + tmp30_1 + tmp3_1 + tmp4_1 + tmp6_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp27_1 + tmp29_1 + tmp2_1 + tmp31_1 + tmp32_1 + tmp36_1;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=tmp0_1 + tmp17_1 + tmp20_1 + tmp30_1 + tmp3_1 + tmp4_1 + tmp6_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=tmp0_1 + tmp17_1 + tmp20_1 + tmp23_1 + tmp24_1 + tmp2_1 + tmp6_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp16_1 + tmp29_1 + tmp2_1 + tmp31_1 + tmp32_1 + tmp37_1;
                                        EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp11_1 + tmp13_1 + tmp15_1 + tmp33_1 + tmp34_1 + tmp38_1 + tmp39_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp11_1 + tmp14_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=tmp0_1 + tmp17_1 + tmp20_1 + tmp23_1 + tmp24_1 + tmp2_1 + tmp6_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp15_1 + tmp35_1 + tmp38_1 + tmp39_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=tmp21_1 + tmp25_1 + tmp26_1 + tmp28_1 + tmp37_1 + tmp3_1 + tmp4_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=tmp35_1 + tmp43_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp48_1;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=tmp0_1 + tmp18_1 + tmp19_1 + tmp1_1 + tmp21_1 + tmp22_1 + tmp36_1 + tmp5_1;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=tmp19_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp24_1 + tmp26_1 + tmp37_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=tmp35_1 + tmp43_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp48_1;
                                        EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp11_1 + tmp14_1 + tmp42_1 + tmp44_1 + tmp49_1 + tmp50_1 + tmp51_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp23_1 + tmp24_1 + tmp30_1 + tmp5_1 + tmp6_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp11_1 + tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp44_1 + tmp51_1;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=tmp13_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp51_1;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=tmp35_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp51_1 + tmp52_1;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=tmp0_1 + tmp18_1 + tmp19_1 + tmp1_1 + tmp21_1 + tmp22_1 + tmp36_1 + tmp5_1;
                                        EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp10_1 + tmp12_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp52_1 + tmp53_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp11_1 + tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp44_1 + tmp51_1;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=tmp13_1 + tmp43_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp52_1;
                                        EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp11_1 + tmp38_1 + tmp39_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp49_1 + tmp50_1;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=tmp19_1 + tmp21_1 + tmp22_1 + tmp26_1 + tmp27_1 + tmp3_1 + tmp4_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=tmp35_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp51_1 + tmp52_1;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=tmp13_1 + tmp43_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp52_1;
                                        EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp11_1 + tmp38_1 + tmp39_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp49_1 + tmp50_1;
                                        EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp10_1 + tmp12_1 + tmp42_1 + tmp45_1 + tmp49_1 + tmp50_1 + tmp53_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=tmp19_1 + tmp21_1 + tmp22_1 + tmp26_1 + tmp27_1 + tmp3_1 + tmp4_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=tmp13_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp51_1;
                                        EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp11_1 + tmp14_1 + tmp42_1 + tmp44_1 + tmp49_1 + tmp50_1 + tmp51_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp33_1 + tmp34_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp45_1 + tmp53_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=tmp0_1 + tmp17_1 + tmp18_1 + tmp20_1 + tmp21_1 + tmp25_1 + tmp28_1 + tmp36_1;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=tmp21_1 + tmp23_1 + tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp11_1 + tmp14_1 + tmp15_1 + tmp33_1 + tmp34_1 + tmp35_1 + tmp8_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp33_1 + tmp34_1 + tmp42_1 + tmp49_1 + tmp50_1 + tmp52_1 + tmp53_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=tmp21_1 + tmp25_1 + tmp26_1 + tmp28_1 + tmp37_1 + tmp3_1 + tmp4_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp15_1 + tmp35_1 + tmp38_1 + tmp39_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp33_1 + tmp34_1 + tmp42_1 + tmp49_1 + tmp50_1 + tmp52_1 + tmp53_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=tmp19_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp24_1 + tmp26_1 + tmp37_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp16_1 + tmp27_1 + tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=tmp0_1 + tmp16_1 + tmp18_1 + tmp1_1 + tmp21_1 + tmp25_1 + tmp28_1 + tmp5_1;
                                        EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp16_1 + tmp29_1 + tmp2_1 + tmp31_1 + tmp32_1 + tmp37_1;
                                        EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp27_1 + tmp29_1 + tmp2_1 + tmp31_1 + tmp32_1 + tmp36_1;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=tmp0_1 + tmp16_1 + tmp18_1 + tmp1_1 + tmp21_1 + tmp25_1 + tmp28_1 + tmp5_1;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=tmp0_1 + tmp17_1 + tmp18_1 + tmp20_1 + tmp21_1 + tmp25_1 + tmp28_1 + tmp36_1;
                                        EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp11_1 + tmp13_1 + tmp15_1 + tmp33_1 + tmp34_1 + tmp38_1 + tmp39_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp10_1 + tmp12_1 + tmp42_1 + tmp45_1 + tmp49_1 + tmp50_1 + tmp53_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp23_1 + tmp24_1 + tmp30_1 + tmp5_1 + tmp6_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp36_1 + tmp37_1;
                                        EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp33_1 + tmp34_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp45_1 + tmp53_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp11_1 + tmp14_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp36_1 + tmp37_1;
                                        EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp10_1 + tmp12_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp52_1 + tmp53_1 + tmp8_1;
                                    }
                                }
                            }
                        }
                        ///////////////
                        // process B //
                        ///////////////
                        if (!B.isEmpty()) {
                            add_EM_S=true;
                            const double* B_p=const_cast<escript::Data*>(&B)->getSampleDataRO(e);
                            if (B.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const register double B_0_0 = B_p[INDEX4(k,0,m,0, numEq,3,numComp)];
                                        const register double B_1_0 = B_p[INDEX4(k,1,m,0, numEq,3,numComp)];
                                        const register double B_2_0 = B_p[INDEX4(k,2,m,0, numEq,3,numComp)];
                                        const register double B_0_1 = B_p[INDEX4(k,0,m,1, numEq,3,numComp)];
                                        const register double B_1_1 = B_p[INDEX4(k,1,m,1, numEq,3,numComp)];
                                        const register double B_2_1 = B_p[INDEX4(k,2,m,1, numEq,3,numComp)];
                                        const register double B_0_2 = B_p[INDEX4(k,0,m,2, numEq,3,numComp)];
                                        const register double B_1_2 = B_p[INDEX4(k,1,m,2, numEq,3,numComp)];
                                        const register double B_2_2 = B_p[INDEX4(k,2,m,2, numEq,3,numComp)];
                                        const register double B_0_3 = B_p[INDEX4(k,0,m,3, numEq,3,numComp)];
                                        const register double B_1_3 = B_p[INDEX4(k,1,m,3, numEq,3,numComp)];
                                        const register double B_2_3 = B_p[INDEX4(k,2,m,3, numEq,3,numComp)];
                                        const register double B_0_4 = B_p[INDEX4(k,0,m,4, numEq,3,numComp)];
                                        const register double B_1_4 = B_p[INDEX4(k,1,m,4, numEq,3,numComp)];
                                        const register double B_2_4 = B_p[INDEX4(k,2,m,4, numEq,3,numComp)];
                                        const register double B_0_5 = B_p[INDEX4(k,0,m,5, numEq,3,numComp)];
                                        const register double B_1_5 = B_p[INDEX4(k,1,m,5, numEq,3,numComp)];
                                        const register double B_2_5 = B_p[INDEX4(k,2,m,5, numEq,3,numComp)];
                                        const register double B_0_6 = B_p[INDEX4(k,0,m,6, numEq,3,numComp)];
                                        const register double B_1_6 = B_p[INDEX4(k,1,m,6, numEq,3,numComp)];
                                        const register double B_2_6 = B_p[INDEX4(k,2,m,6, numEq,3,numComp)];
                                        const register double B_0_7 = B_p[INDEX4(k,0,m,7, numEq,3,numComp)];
                                        const register double B_1_7 = B_p[INDEX4(k,1,m,7, numEq,3,numComp)];
                                        const register double B_2_7 = B_p[INDEX4(k,2,m,7, numEq,3,numComp)];
                                        const register double tmp24_0 = B_2_0 + B_2_2;
                                        const register double tmp19_0 = B_2_0 + B_2_1 + B_2_2 + B_2_3;
                                        const register double tmp20_0 = B_1_0 + B_1_5;
                                        const register double tmp22_0 = B_2_1 + B_2_3;
                                        const register double tmp29_0 = B_2_2 + B_2_3;
                                        const register double tmp34_0 = B_1_0 + B_1_1 + B_1_4 + B_1_5;
                                        const register double tmp15_0 = B_2_4 + B_2_5 + B_2_6 + B_2_7;
                                        const register double tmp7_0 = B_0_0 + B_0_4;
                                        const register double tmp40_0 = B_1_1 + B_1_4;
                                        const register double tmp14_0 = B_1_4 + B_1_5;
                                        const register double tmp35_0 = B_1_2 + B_1_3 + B_1_6 + B_1_7;
                                        const register double tmp17_0 = B_1_6 + B_1_7;
                                        const register double tmp8_0 = B_1_2 + B_1_6;
                                        const register double tmp28_0 = B_2_4 + B_2_5;
                                        const register double tmp10_0 = B_0_1 + B_0_3;
                                        const register double tmp9_0 = B_2_5 + B_2_6;
                                        const register double tmp30_0 = B_2_0 + B_2_1;
                                        const register double tmp27_0 = B_0_0 + B_0_6;
                                        const register double tmp32_0 = B_0_0 + B_0_2 + B_0_4 + B_0_6;
                                        const register double tmp16_0 = B_0_0 + B_0_2;
                                        const register double tmp2_0 = B_1_3 + B_1_7;
                                        const register double tmp3_0 = B_2_1 + B_2_2;
                                        const register double tmp33_0 = B_0_1 + B_0_3 + B_0_5 + B_0_7;
                                        const register double tmp23_0 = B_2_5 + B_2_7;
                                        const register double tmp36_0 = B_2_4 + B_2_7;
                                        const register double tmp39_0 = B_0_3 + B_0_5;
                                        const register double tmp41_0 = B_1_3 + B_1_6;
                                        const register double tmp5_0 = B_1_0 + B_1_4;
                                        const register double tmp18_0 = B_1_0 + B_1_1;
                                        const register double tmp0_0 = B_0_1 + B_0_5;
                                        const register double tmp37_0 = B_2_0 + B_2_3;
                                        const register double tmp25_0 = B_2_4 + B_2_6;
                                        const register double tmp38_0 = B_0_2 + B_0_4;
                                        const register double tmp31_0 = B_2_6 + B_2_7;
                                        const register double tmp1_0 = B_0_2 + B_0_6;
                                        const register double tmp11_0 = B_0_4 + B_0_6;
                                        const register double tmp21_0 = B_1_2 + B_1_7;
                                        const register double tmp6_0 = B_1_1 + B_1_5;
                                        const register double tmp13_0 = B_1_2 + B_1_3;
                                        const register double tmp12_0 = B_0_5 + B_0_7;
                                        const register double tmp26_0 = B_0_1 + B_0_7;
                                        const register double tmp4_0 = B_0_3 + B_0_7;
                                        const register double tmp324_1 = B_0_6*w104;
                                        const register double tmp209_1 = B_1_4*w114;
                                        const register double tmp255_1 = B_2_3*w103;
                                        const register double tmp37_1 = B_1_1*w111;
                                        const register double tmp326_1 = tmp38_0*w108;
                                        const register double tmp179_1 = tmp21_0*w94;
                                        const register double tmp102_1 = tmp4_0*w93;
                                        const register double tmp251_1 = B_0_2*w125;
                                        const register double tmp321_1 = B_0_4*w90;
                                        const register double tmp198_1 = tmp41_0*w97;
                                        const register double tmp15_1 = tmp11_0*w108;
                                        const register double tmp158_1 = tmp17_0*w99;
                                        const register double tmp138_1 = tmp35_0*w94;
                                        const register double tmp5_1 = B_2_3*w100;
                                        const register double tmp51_1 = tmp25_0*w103;
                                        const register double tmp258_1 = B_2_7*w100;
                                        const register double tmp221_1 = B_0_4*w125;
                                        const register double tmp70_1 = B_0_5*w98;
                                        const register double tmp202_1 = tmp38_0*w96;
                                        const register double tmp247_1 = tmp41_0*w94;
                                        const register double tmp122_1 = tmp34_0*w94;
                                        const register double tmp349_1 = B_2_1*w101;
                                        const register double tmp167_1 = tmp31_0*w103;
                                        const register double tmp408_1 = tmp8_0*w111;
                                        const register double tmp20_1 = tmp16_0*w104;
                                        const register double tmp2_1 = B_2_7*w103;
                                        const register double tmp398_1 = B_1_7*w99;
                                        const register double tmp262_1 = B_1_7*w119;
                                        const register double tmp97_1 = tmp5_0*w94;
                                        const register double tmp157_1 = tmp30_0*w92;
                                        const register double tmp67_1 = tmp18_0*w107;
                                        const register double tmp144_1 = tmp7_0*w110;
                                        const register double tmp264_1 = B_0_0*w120;
                                        const register double tmp396_1 = tmp6_0*w97;
                                        const register double tmp218_1 = B_1_6*w111;
                                        const register double tmp147_1 = tmp8_0*w109;
                                        const register double tmp39_1 = tmp12_0*w108;
                                        const register double tmp58_1 = tmp15_0*w112;
                                        const register double tmp214_1 = B_1_3*w115;
                                        const register double tmp69_1 = tmp30_0*w95;
                                        const register double tmp54_1 = tmp18_0*w99;
                                        const register double tmp261_1 = B_2_0*w101;
                                        const register double tmp390_1 = tmp24_0*w103;
                                        const register double tmp374_1 = B_2_1*w103;
                                        const register double tmp151_1 = tmp2_0*w105;
                                        const register double tmp274_1 = tmp22_0*w95;
                                        const register double tmp126_1 = tmp16_0*w108;
                                        const register double tmp302_1 = B_2_6*w122;
                                        const register double tmp14_1 = tmp10_0*w106;
                                        const register double tmp190_1 = B_1_2*w99;
                                        const register double tmp9_1 = B_2_4*w101;
                                        const register double tmp150_1 = B_2_7*w101;
                                        const register double tmp66_1 = tmp29_0*w92;
                                        const register double tmp334_1 = B_0_7*w124;
                                        const register double tmp252_1 = tmp38_0*w93;
                                        const register double tmp392_1 = tmp5_0*w99;
                                        const register double tmp4_1 = tmp2_0*w99;
                                        const register double tmp63_1 = B_0_4*w121;
                                        const register double tmp152_1 = tmp5_0*w111;
                                        const register double tmp133_1 = tmp4_0*w96;
                                        const register double tmp195_1 = tmp40_0*w94;
                                        const register double tmp145_1 = B_2_3*w92;
                                        const register double tmp35_1 = tmp25_0*w116;
                                        const register double tmp316_1 = B_0_3*w98;
                                        const register double tmp208_1 = B_2_0*w103;
                                        const register double tmp373_1 = B_1_2*w115;
                                        const register double tmp62_1 = tmp26_0*w93;
                                        const register double tmp85_1 = tmp10_0*w90;
                                        const register double tmp187_1 = tmp11_0*w106;
                                        const register double tmp342_1 = B_2_5*w92;
                                        const register double tmp77_1 = tmp33_0*w108;
                                        const register double tmp121_1 = B_2_4*w116;
                                        const register double tmp379_1 = B_2_2*w101;
                                        const register double tmp391_1 = tmp23_0*w92;
                                        const register double tmp451_1 = tmp22_0*w116;
                                        const register double tmp73_1 = B_0_2*w90;
                                        const register double tmp436_1 = B_1_6*w114;
                                        const register double tmp233_1 = tmp29_0*w117;
                                        const register double tmp426_1 = tmp23_0*w113;
                                        const register double tmp124_1 = tmp11_0*w104;
                                        const register double tmp148_1 = B_2_0*w100;
                                        const register double tmp314_1 = tmp29_0*w113;
                                        const register double tmp163_1 = B_0_0*w124;
                                        const register double tmp301_1 = B_1_6*w115;
                                        const register double tmp171_1 = B_2_4*w122;
                                        const register double tmp7_1 = tmp4_0*w98;
                                        const register double tmp11_1 = tmp7_0*w90;
                                        const register double tmp308_1 = B_2_2*w116;
                                        const register double tmp117_1 = B_2_3*w113;
                                        const register double tmp403_1 = B_0_5*w104;
                                        const register double tmp178_1 = tmp9_0*w112;
                                        const register double tmp128_1 = tmp17_0*w107;
                                        const register double tmp89_1 = tmp1_0*w110;
                                        const register double tmp434_1 = tmp19_0*w95;
                                        const register double tmp270_1 = B_1_7*w114;
                                        const register double tmp358_1 = tmp6_0*w99;
                                        const register double tmp411_1 = tmp14_0*w107;
                                        const register double tmp433_1 = B_0_1*w125;
                                        const register double tmp376_1 = B_2_6*w92;
                                        const register double tmp42_1 = B_1_6*w99;
                                        const register double tmp72_1 = tmp13_0*w105;
                                        const register double tmp287_1 = tmp2_0*w111;
                                        const register double tmp191_1 = tmp24_0*w113;
                                        const register double tmp6_1 = tmp3_0*w95;
                                        const register double tmp309_1 = B_1_4*w105;
                                        const register double tmp367_1 = B_0_5*w125;
                                        const register double tmp234_1 = tmp28_0*w112;
                                        const register double tmp339_1 = B_1_2*w111;
                                        const register double tmp335_1 = tmp41_0*w107;
                                        const register double tmp173_1 = B_2_7*w113;
                                        const register double tmp384_1 = tmp28_0*w113;
                                        const register double tmp81_1 = tmp19_0*w112;
                                        const register double tmp328_1 = B_2_7*w122;
                                        const register double tmp413_1 = tmp30_0*w113;
                                        const register double tmp111_1 = B_2_6*w101;
                                        const register double tmp378_1 = B_0_4*w98;
                                        const register double tmp211_1 = tmp21_0*w107;
                                        const register double tmp79_1 = tmp35_0*w109;
                                        const register double tmp354_1 = tmp1_0*w93;
                                        const register double tmp333_1 = B_2_0*w123;
                                        const register double tmp169_1 = B_0_3*w121;
                                        const register double tmp300_1 = B_0_1*w121;
                                        const register double tmp118_1 = B_2_7*w123;
                                        const register double tmp372_1 = B_0_5*w121;
                                        const register double tmp104_1 = B_2_5*w103;
                                        const register double tmp304_1 = B_2_5*w113;
                                        const register double tmp397_1 = tmp22_0*w102;
                                        const register double tmp61_1 = tmp14_0*w97;
                                        const register double tmp290_1 = tmp4_0*w106;
                                        const register double tmp269_1 = tmp23_0*w103;
                                        const register double tmp361_1 = tmp5_0*w97;
                                        const register double tmp189_1 = tmp16_0*w110;
                                        const register double tmp116_1 = B_2_0*w122;
                                        const register double tmp449_1 = tmp24_0*w117;
                                        const register double tmp268_1 = tmp11_0*w96;
                                        const register double tmp325_1 = tmp39_0*w106;
                                        const register double tmp27_1 = tmp20_0*w107;
                                        const register double tmp174_1 = B_2_3*w123;
                                        const register double tmp422_1 = tmp14_0*w99;
                                        const register double tmp146_1 = tmp6_0*w107;
                                        const register double tmp277_1 = B_1_2*w105;
                                        const register double tmp443_1 = B_1_6*w118;
                                        const register double tmp446_1 = B_1_7*w105;
                                        const register double tmp91_1 = tmp8_0*w99;
                                        const register double tmp405_1 = B_0_3*w125;
                                        const register double tmp57_1 = tmp17_0*w91;
                                        const register double tmp455_1 = B_2_6*w103;
                                        const register double tmp110_1 = tmp0_0*w98;
                                        const register double tmp341_1 = B_1_6*w119;
                                        const register double tmp203_1 = B_0_7*w98;
                                        const register double tmp371_1 = B_2_7*w116;
                                        const register double tmp382_1 = B_0_3*w90;
                                        const register double tmp320_1 = tmp28_0*w116;
                                        const register double tmp107_1 = tmp2_0*w109;
                                        const register double tmp149_1 = tmp4_0*w104;
                                        const register double tmp159_1 = tmp29_0*w95;
                                        const register double tmp253_1 = B_0_7*w121;
                                        const register double tmp90_1 = B_2_1*w122;
                                        const register double tmp360_1 = tmp2_0*w94;
                                        const register double tmp414_1 = tmp28_0*w117;
                                        const register double tmp26_1 = tmp16_0*w96;
                                        const register double tmp84_1 = tmp11_0*w98;
                                        const register double tmp225_1 = tmp1_0*w108;
                                        const register double tmp1_1 = tmp1_0*w96;
                                        const register double tmp344_1 = B_2_6*w100;
                                        const register double tmp400_1 = B_1_5*w119;
                                        const register double tmp345_1 = tmp36_0*w95;
                                        const register double tmp407_1 = tmp5_0*w109;
                                        const register double tmp228_1 = B_2_2*w122;
                                        const register double tmp113_1 = tmp8_0*w105;
                                        const register double tmp347_1 = B_0_0*w104;
                                        const register double tmp286_1 = tmp23_0*w95;
                                        const register double tmp395_1 = tmp8_0*w94;
                                        const register double tmp401_1 = B_1_2*w118;
                                        const register double tmp380_1 = B_1_7*w111;
                                        const register double tmp327_1 = B_0_1*w110;
                                        const register double tmp55_1 = tmp19_0*w117;
                                        const register double tmp291_1 = tmp7_0*w108;
                                        const register double tmp31_1 = tmp10_0*w98;
                                        const register double tmp125_1 = tmp12_0*w106;
                                        const register double tmp330_1 = tmp40_0*w109;
                                        const register double tmp323_1 = tmp17_0*w97;
                                        const register double tmp136_1 = tmp31_0*w95;
                                        const register double tmp16_1 = tmp12_0*w110;
                                        const register double tmp418_1 = B_0_6*w90;
                                        const register double tmp428_1 = tmp25_0*w112;
                                        const register double tmp385_1 = tmp30_0*w117;
                                        const register double tmp351_1 = B_1_1*w118;
                                        const register double tmp165_1 = B_0_7*w125;
                                        const register double tmp298_1 = tmp29_0*w102;
                                        const register double tmp295_1 = tmp34_0*w109;
                                        const register double tmp296_1 = tmp28_0*w95;
                                        const register double tmp283_1 = tmp25_0*w92;
                                        const register double tmp230_1 = B_2_5*w123;
                                        const register double tmp350_1 = B_0_1*w124;
                                        const register double tmp293_1 = tmp31_0*w92;
                                        const register double tmp8_1 = tmp5_0*w91;
                                        const register double tmp215_1 = tmp9_0*w95;
                                        const register double tmp329_1 = B_1_0*w114;
                                        const register double tmp115_1 = tmp6_0*w111;
                                        const register double tmp387_1 = tmp29_0*w116;
                                        const register double tmp442_1 = B_1_1*w119;
                                        const register double tmp281_1 = tmp33_0*w96;
                                        const register double tmp415_1 = B_0_1*w98;
                                        const register double tmp311_1 = B_1_3*w111;
                                        const register double tmp50_1 = tmp23_0*w102;
                                        const register double tmp294_1 = tmp35_0*w107;
                                        const register double tmp362_1 = tmp27_0*w106;
                                        const register double tmp340_1 = B_2_2*w103;
                                        const register double tmp282_1 = tmp22_0*w103;
                                        const register double tmp254_1 = tmp39_0*w96;
                                        const register double tmp186_1 = tmp12_0*w104;
                                        const register double tmp106_1 = tmp5_0*w107;
                                        const register double tmp170_1 = tmp26_0*w96;
                                        const register double tmp427_1 = tmp22_0*w117;
                                        const register double tmp452_1 = B_2_1*w92;
                                        const register double tmp437_1 = B_1_1*w115;
                                        const register double tmp431_1 = B_0_0*w110;
                                        const register double tmp201_1 = B_0_6*w121;
                                        const register double tmp310_1 = B_0_6*w120;
                                        const register double tmp331_1 = B_2_4*w113;
                                        const register double tmp212_1 = tmp20_0*w109;
                                        const register double tmp47_1 = B_1_4*w119;
                                        const register double tmp420_1 = B_0_7*w120;
                                        const register double tmp38_1 = tmp16_0*w106;
                                        const register double tmp184_1 = tmp20_0*w97;
                                        const register double tmp280_1 = tmp32_0*w93;
                                        const register double tmp123_1 = tmp35_0*w97;
                                        const register double tmp359_1 = tmp8_0*w91;
                                        const register double tmp95_1 = tmp6_0*w91;
                                        const register double tmp246_1 = tmp36_0*w112;
                                        const register double tmp241_1 = B_2_6*w113;
                                        const register double tmp3_1 = B_2_0*w92;
                                        const register double tmp162_1 = tmp14_0*w94;
                                        const register double tmp406_1 = tmp2_0*w107;
                                        const register double tmp175_1 = tmp3_0*w117;
                                        const register double tmp205_1 = B_0_1*w120;
                                        const register double tmp416_1 = tmp29_0*w112;
                                        const register double tmp48_1 = tmp21_0*w97;
                                        const register double tmp101_1 = tmp32_0*w96;
                                        const register double tmp231_1 = B_2_6*w116;
                                        const register double tmp366_1 = B_0_2*w124;
                                        const register double tmp34_1 = tmp11_0*w90;
                                        const register double tmp200_1 = tmp39_0*w93;
                                        const register double tmp44_1 = tmp10_0*w104;
                                        const register double tmp346_1 = B_0_6*w125;
                                        const register double tmp65_1 = tmp28_0*w103;
                                        const register double tmp248_1 = B_1_2*w119;
                                        const register double tmp288_1 = tmp5_0*w105;
                                        const register double tmp227_1 = tmp4_0*w110;
                                        const register double tmp343_1 = B_1_4*w99;
                                        const register double tmp139_1 = tmp0_0*w90;
                                        const register double tmp41_1 = tmp22_0*w92;
                                        const register double tmp305_1 = B_2_1*w123;
                                        const register double tmp307_1 = B_0_7*w90;
                                        const register double tmp217_1 = B_2_3*w101;
                                        const register double tmp219_1 = B_0_3*w124;
                                        const register double tmp105_1 = B_2_2*w92;
                                        const register double tmp423_1 = tmp13_0*w91;
                                        const register double tmp381_1 = B_1_0*w105;
                                        const register double tmp166_1 = tmp28_0*w102;
                                        const register double tmp259_1 = B_0_6*w98;
                                        const register double tmp402_1 = B_0_2*w110;
                                        const register double tmp263_1 = B_0_1*w90;
                                        const register double tmp276_1 = B_1_5*w111;
                                        const register double tmp243_1 = B_2_2*w123;
                                        const register double tmp131_1 = tmp13_0*w111;
                                        const register double tmp265_1 = B_1_0*w118;
                                        const register double tmp337_1 = B_1_5*w105;
                                        const register double tmp172_1 = B_1_1*w99;
                                        const register double tmp161_1 = tmp18_0*w91;
                                        const register double tmp88_1 = tmp4_0*w108;
                                        const register double tmp236_1 = B_1_5*w118;
                                        const register double tmp273_1 = tmp41_0*w109;
                                        const register double tmp154_1 = tmp38_0*w106;
                                        const register double tmp440_1 = B_1_3*w99;
                                        const register double tmp33_1 = B_1_3*w114;
                                        const register double tmp141_1 = tmp30_0*w102;
                                        const register double tmp393_1 = tmp25_0*w95;
                                        const register double tmp289_1 = tmp24_0*w102;
                                        const register double tmp180_1 = B_1_3*w119;
                                        const register double tmp86_1 = tmp36_0*w117;
                                        const register double tmp444_1 = B_1_5*w115;
                                        const register double tmp447_1 = B_1_0*w111;
                                        const register double tmp315_1 = tmp31_0*w117;
                                        const register double tmp120_1 = tmp3_0*w112;
                                        const register double tmp355_1 = tmp0_0*w96;
                                        const register double tmp223_1 = B_0_5*w110;
                                        const register double tmp60_1 = tmp12_0*w90;
                                        const register double tmp239_1 = B_2_5*w122;
                                        const register double tmp194_1 = tmp22_0*w112;
                                        const register double tmp28_1 = tmp21_0*w109;
                                        const register double tmp24_1 = tmp12_0*w93;
                                        const register double tmp424_1 = tmp17_0*w94;
                                        const register double tmp237_1 = B_0_4*w104;
                                        const register double tmp322_1 = B_0_5*w120;
                                        const register double tmp453_1 = B_2_2*w100;
                                        const register double tmp238_1 = B_0_3*w110;
                                        const register double tmp370_1 = B_2_4*w123;
                                        const register double tmp155_1 = tmp39_0*w108;
                                        const register double tmp23_1 = tmp19_0*w102;
                                        const register double tmp181_1 = B_2_0*w116;
                                        const register double tmp164_1 = tmp13_0*w97;
                                        const register double tmp168_1 = tmp27_0*w93;
                                        const register double tmp285_1 = tmp6_0*w109;
                                        const register double tmp32_1 = tmp24_0*w112;
                                        const register double tmp210_1 = B_2_7*w92;
                                        const register double tmp399_1 = B_1_0*w91;
                                        const register double tmp183_1 = B_0_4*w120;
                                        const register double tmp185_1 = B_1_4*w118;
                                        const register double tmp40_1 = tmp11_0*w110;
                                        const register double tmp59_1 = tmp13_0*w94;
                                        const register double tmp279_1 = tmp25_0*w102;
                                        const register double tmp242_1 = B_1_7*w91;
                                        const register double tmp454_1 = B_2_5*w101;
                                        const register double tmp114_1 = tmp1_0*w90;
                                        const register double tmp319_1 = tmp18_0*w94;
                                        const register double tmp394_1 = tmp2_0*w91;
                                        const register double tmp229_1 = B_2_1*w113;
                                        const register double tmp93_1 = B_2_6*w123;
                                        const register double tmp278_1 = tmp16_0*w90;
                                        const register double tmp46_1 = tmp20_0*w94;
                                        const register double tmp87_1 = tmp7_0*w106;
                                        const register double tmp421_1 = tmp18_0*w111;
                                        const register double tmp313_1 = tmp13_0*w99;
                                        const register double tmp364_1 = B_0_4*w110;
                                        const register double tmp412_1 = tmp13_0*w109;
                                        const register double tmp30_1 = tmp23_0*w117;
                                        const register double tmp365_1 = B_0_3*w104;
                                        const register double tmp425_1 = tmp18_0*w97;
                                        const register double tmp134_1 = tmp29_0*w103;
                                        const register double tmp143_1 = tmp0_0*w108;
                                        const register double tmp213_1 = B_2_4*w100;
                                        const register double tmp435_1 = tmp15_0*w102;
                                        const register double tmp78_1 = tmp34_0*w107;
                                        const register double tmp76_1 = tmp32_0*w106;
                                        const register double tmp132_1 = tmp7_0*w93;
                                        const register double tmp348_1 = B_1_3*w91;
                                        const register double tmp68_1 = tmp17_0*w109;
                                        const register double tmp130_1 = tmp14_0*w105;
                                        const register double tmp142_1 = tmp1_0*w106;
                                        const register double tmp153_1 = B_2_4*w103;
                                        const register double tmp432_1 = B_0_6*w124;
                                        const register double tmp256_1 = B_2_4*w92;
                                        const register double tmp204_1 = B_0_0*w90;
                                        const register double tmp82_1 = tmp16_0*w93;
                                        const register double tmp222_1 = tmp3_0*w102;
                                        const register double tmp272_1 = tmp40_0*w107;
                                        const register double tmp19_1 = tmp15_0*w95;
                                        const register double tmp226_1 = tmp7_0*w104;
                                        const register double tmp25_1 = B_1_4*w115;
                                        const register double tmp199_1 = B_1_7*w118;
                                        const register double tmp303_1 = B_1_1*w114;
                                        const register double tmp318_1 = tmp30_0*w112;
                                        const register double tmp100_1 = tmp33_0*w93;
                                        const register double tmp439_1 = B_1_4*w111;
                                        const register double tmp17_1 = tmp13_0*w107;
                                        const register double tmp83_1 = tmp12_0*w96;
                                        const register double tmp45_1 = B_1_1*w91;
                                        const register double tmp75_1 = tmp14_0*w111;
                                        const register double tmp368_1 = B_2_3*w122;
                                        const register double tmp216_1 = B_0_2*w104;
                                        const register double tmp429_1 = tmp24_0*w116;
                                        const register double tmp388_1 = tmp33_0*w106;
                                        const register double tmp13_1 = tmp9_0*w102;
                                        const register double tmp18_1 = tmp14_0*w109;
                                        const register double tmp112_1 = tmp36_0*w102;
                                        const register double tmp140_1 = tmp34_0*w97;
                                        const register double tmp383_1 = B_0_2*w120;
                                        const register double tmp119_1 = tmp9_0*w117;
                                        const register double tmp36_1 = B_1_6*w105;
                                        const register double tmp266_1 = tmp10_0*w93;
                                        const register double tmp336_1 = B_2_3*w116;
                                        const register double tmp419_1 = tmp31_0*w116;
                                        const register double tmp363_1 = tmp26_0*w108;
                                        const register double tmp196_1 = B_1_0*w119;
                                        const register double tmp182_1 = B_0_5*w90;
                                        const register double tmp0_1 = tmp0_0*w93;
                                        const register double tmp297_1 = tmp1_0*w104;
                                        const register double tmp12_1 = tmp8_0*w97;
                                        const register double tmp244_1 = tmp37_0*w117;
                                        const register double tmp332_1 = B_1_7*w115;
                                        const register double tmp156_1 = B_0_6*w110;
                                        const register double tmp127_1 = tmp10_0*w110;
                                        const register double tmp257_1 = B_1_5*w99;
                                        const register double tmp64_1 = tmp27_0*w96;
                                        const register double tmp389_1 = tmp32_0*w108;
                                        const register double tmp10_1 = tmp6_0*w94;
                                        const register double tmp317_1 = tmp14_0*w91;
                                        const register double tmp275_1 = tmp12_0*w98;
                                        const register double tmp99_1 = tmp2_0*w97;
                                        const register double tmp192_1 = tmp25_0*w117;
                                        const register double tmp224_1 = tmp0_0*w106;
                                        const register double tmp410_1 = B_0_0*w121;
                                        const register double tmp409_1 = tmp6_0*w105;
                                        const register double tmp441_1 = B_1_4*w91;
                                        const register double tmp206_1 = tmp26_0*w106;
                                        const register double tmp232_1 = tmp31_0*w113;
                                        const register double tmp108_1 = B_2_1*w100;
                                        const register double tmp21_1 = tmp17_0*w111;
                                        const register double tmp52_1 = tmp11_0*w93;
                                        const register double tmp103_1 = tmp7_0*w96;
                                        const register double tmp197_1 = tmp23_0*w116;
                                        const register double tmp438_1 = B_1_3*w105;
                                        const register double tmp235_1 = tmp30_0*w116;
                                        const register double tmp29_1 = tmp22_0*w113;
                                        const register double tmp98_1 = B_2_5*w116;
                                        const register double tmp369_1 = B_2_0*w113;
                                        const register double tmp74_1 = B_0_3*w120;
                                        const register double tmp137_1 = tmp1_0*w98;
                                        const register double tmp284_1 = tmp8_0*w107;
                                        const register double tmp109_1 = tmp37_0*w95;
                                        const register double tmp448_1 = tmp25_0*w113;
                                        const register double tmp56_1 = tmp16_0*w98;
                                        const register double tmp375_1 = B_1_5*w114;
                                        const register double tmp188_1 = tmp10_0*w108;
                                        const register double tmp94_1 = tmp0_0*w104;
                                        const register double tmp260_1 = B_1_2*w91;
                                        const register double tmp49_1 = B_1_3*w118;
                                        const register double tmp404_1 = B_0_4*w124;
                                        const register double tmp430_1 = B_0_7*w104;
                                        const register double tmp299_1 = tmp0_0*w110;
                                        const register double tmp417_1 = tmp17_0*w105;
                                        const register double tmp445_1 = B_1_2*w114;
                                        const register double tmp353_1 = B_0_7*w110;
                                        const register double tmp135_1 = tmp28_0*w92;
                                        const register double tmp357_1 = tmp4_0*w90;
                                        const register double tmp356_1 = tmp7_0*w98;
                                        const register double tmp267_1 = B_1_0*w115;
                                        const register double tmp220_1 = B_1_1*w105;
                                        const register double tmp338_1 = B_0_0*w125;
                                        const register double tmp292_1 = tmp30_0*w103;
                                        const register double tmp53_1 = tmp10_0*w96;
                                        const register double tmp43_1 = tmp24_0*w95;
                                        const register double tmp386_1 = tmp31_0*w112;
                                        const register double tmp92_1 = B_2_2*w113;
                                        const register double tmp245_1 = B_0_5*w124;
                                        const register double tmp306_1 = B_0_0*w98;
                                        const register double tmp96_1 = tmp37_0*w112;
                                        const register double tmp177_1 = B_1_6*w91;
                                        const register double tmp176_1 = B_0_2*w98;
                                        const register double tmp71_1 = tmp31_0*w102;
                                        const register double tmp129_1 = tmp18_0*w109;
                                        const register double tmp352_1 = tmp37_0*w102;
                                        const register double tmp312_1 = B_0_2*w121;
                                        const register double tmp80_1 = tmp15_0*w117;
                                        const register double tmp450_1 = tmp23_0*w112;
                                        const register double tmp271_1 = tmp24_0*w92;
                                        const register double tmp22_1 = tmp18_0*w105;
                                        const register double tmp250_1 = tmp40_0*w97;
                                        const register double tmp249_1 = B_2_1*w116;
                                        const register double tmp207_1 = tmp27_0*w108;
                                        const register double tmp193_1 = B_1_5*w91;
                                        const register double tmp240_1 = B_1_0*w99;
                                        const register double tmp160_1 = B_0_1*w104;
                                        const register double tmp377_1 = B_2_5*w100;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=tmp0_1 + tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp14_1 + tmp15_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp33_1 + tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                                        EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1 + tmp60_1 + tmp61_1;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=tmp62_1 + tmp63_1 + tmp64_1 + tmp65_1 + tmp66_1 + tmp67_1 + tmp68_1 + tmp69_1 + tmp70_1 + tmp71_1 + tmp72_1 + tmp73_1 + tmp74_1 + tmp75_1;
                                        EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp76_1 + tmp77_1 + tmp78_1 + tmp79_1 + tmp80_1 + tmp81_1;
                                        EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp17_1 + tmp18_1 + tmp19_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp82_1 + tmp83_1 + tmp84_1 + tmp85_1;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=tmp86_1 + tmp87_1 + tmp88_1 + tmp89_1 + tmp90_1 + tmp91_1 + tmp92_1 + tmp93_1 + tmp94_1 + tmp95_1 + tmp96_1 + tmp97_1 + tmp98_1 + tmp99_1;
                                        EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp100_1 + tmp101_1 + tmp78_1 + tmp79_1 + tmp80_1 + tmp81_1;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=tmp102_1 + tmp103_1 + tmp104_1 + tmp105_1 + tmp106_1 + tmp107_1 + tmp108_1 + tmp109_1 + tmp110_1 + tmp111_1 + tmp112_1 + tmp113_1 + tmp114_1 + tmp115_1;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=tmp0_1 + tmp10_1 + tmp116_1 + tmp117_1 + tmp118_1 + tmp119_1 + tmp11_1 + tmp120_1 + tmp121_1 + tmp12_1 + tmp1_1 + tmp4_1 + tmp7_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp122_1 + tmp123_1 + tmp76_1 + tmp77_1 + tmp80_1 + tmp81_1;
                                        EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp124_1 + tmp125_1 + tmp126_1 + tmp127_1 + tmp128_1 + tmp129_1 + tmp130_1 + tmp131_1 + tmp55_1 + tmp58_1;
                                        EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp132_1 + tmp133_1 + tmp134_1 + tmp135_1 + tmp136_1 + tmp137_1 + tmp138_1 + tmp139_1 + tmp140_1 + tmp141_1;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=tmp13_1 + tmp142_1 + tmp143_1 + tmp144_1 + tmp145_1 + tmp146_1 + tmp147_1 + tmp148_1 + tmp149_1 + tmp150_1 + tmp151_1 + tmp152_1 + tmp153_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp128_1 + tmp129_1 + tmp130_1 + tmp131_1 + tmp52_1 + tmp53_1 + tmp55_1 + tmp56_1 + tmp58_1 + tmp60_1;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=tmp154_1 + tmp155_1 + tmp156_1 + tmp157_1 + tmp158_1 + tmp159_1 + tmp160_1 + tmp161_1 + tmp162_1 + tmp163_1 + tmp164_1 + tmp165_1 + tmp166_1 + tmp167_1;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=tmp168_1 + tmp169_1 + tmp170_1 + tmp171_1 + tmp172_1 + tmp173_1 + tmp174_1 + tmp175_1 + tmp176_1 + tmp177_1 + tmp178_1 + tmp179_1 + tmp180_1 + tmp181_1 + tmp182_1 + tmp183_1 + tmp184_1 + tmp185_1;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=tmp186_1 + tmp187_1 + tmp188_1 + tmp189_1 + tmp190_1 + tmp191_1 + tmp192_1 + tmp193_1 + tmp194_1 + tmp195_1 + tmp196_1 + tmp197_1 + tmp198_1 + tmp199_1;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=tmp157_1 + tmp158_1 + tmp159_1 + tmp161_1 + tmp162_1 + tmp164_1 + tmp166_1 + tmp167_1 + tmp200_1 + tmp201_1 + tmp202_1 + tmp203_1 + tmp204_1 + tmp205_1;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=tmp206_1 + tmp207_1 + tmp208_1 + tmp209_1 + tmp210_1 + tmp211_1 + tmp212_1 + tmp213_1 + tmp214_1 + tmp215_1 + tmp216_1 + tmp217_1 + tmp218_1 + tmp219_1 + tmp220_1 + tmp221_1 + tmp222_1 + tmp223_1;
                                        EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp134_1 + tmp135_1 + tmp136_1 + tmp138_1 + tmp140_1 + tmp141_1 + tmp224_1 + tmp225_1 + tmp226_1 + tmp227_1;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=tmp102_1 + tmp103_1 + tmp106_1 + tmp107_1 + tmp110_1 + tmp113_1 + tmp114_1 + tmp115_1 + tmp228_1 + tmp229_1 + tmp230_1 + tmp231_1 + tmp86_1 + tmp96_1;
                                        EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp132_1 + tmp133_1 + tmp137_1 + tmp138_1 + tmp139_1 + tmp140_1 + tmp232_1 + tmp233_1 + tmp234_1 + tmp235_1;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=tmp206_1 + tmp207_1 + tmp236_1 + tmp237_1 + tmp238_1 + tmp239_1 + tmp240_1 + tmp241_1 + tmp242_1 + tmp243_1 + tmp244_1 + tmp245_1 + tmp246_1 + tmp247_1 + tmp248_1 + tmp249_1 + tmp250_1 + tmp251_1;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=tmp215_1 + tmp222_1 + tmp247_1 + tmp250_1 + tmp252_1 + tmp253_1 + tmp254_1 + tmp255_1 + tmp256_1 + tmp257_1 + tmp258_1 + tmp259_1 + tmp260_1 + tmp261_1 + tmp262_1 + tmp263_1 + tmp264_1 + tmp265_1;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=tmp266_1 + tmp267_1 + tmp268_1 + tmp269_1 + tmp270_1 + tmp271_1 + tmp272_1 + tmp273_1 + tmp274_1 + tmp275_1 + tmp276_1 + tmp277_1 + tmp278_1 + tmp279_1;
                                        EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp280_1 + tmp281_1 + tmp282_1 + tmp283_1 + tmp284_1 + tmp285_1 + tmp286_1 + tmp287_1 + tmp288_1 + tmp289_1;
                                        EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp290_1 + tmp291_1 + tmp292_1 + tmp293_1 + tmp294_1 + tmp295_1 + tmp296_1 + tmp297_1 + tmp298_1 + tmp299_1;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=tmp211_1 + tmp212_1 + tmp244_1 + tmp246_1 + tmp252_1 + tmp254_1 + tmp300_1 + tmp301_1 + tmp302_1 + tmp303_1 + tmp304_1 + tmp305_1 + tmp306_1 + tmp307_1 + tmp308_1 + tmp309_1 + tmp310_1 + tmp311_1;
                                        EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp138_1 + tmp140_1 + tmp224_1 + tmp225_1 + tmp226_1 + tmp227_1 + tmp232_1 + tmp233_1 + tmp234_1 + tmp235_1;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=tmp312_1 + tmp313_1 + tmp314_1 + tmp315_1 + tmp316_1 + tmp317_1 + tmp318_1 + tmp319_1 + tmp320_1 + tmp321_1 + tmp322_1 + tmp323_1 + tmp62_1 + tmp64_1;
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=tmp175_1 + tmp178_1 + tmp324_1 + tmp325_1 + tmp326_1 + tmp327_1 + tmp328_1 + tmp329_1 + tmp330_1 + tmp331_1 + tmp332_1 + tmp333_1 + tmp334_1 + tmp335_1 + tmp336_1 + tmp337_1 + tmp338_1 + tmp339_1;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=tmp179_1 + tmp184_1 + tmp325_1 + tmp326_1 + tmp340_1 + tmp341_1 + tmp342_1 + tmp343_1 + tmp344_1 + tmp345_1 + tmp346_1 + tmp347_1 + tmp348_1 + tmp349_1 + tmp350_1 + tmp351_1 + tmp352_1 + tmp353_1;
                                        EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp292_1 + tmp293_1 + tmp294_1 + tmp295_1 + tmp296_1 + tmp298_1 + tmp354_1 + tmp355_1 + tmp356_1 + tmp357_1;
                                        EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp280_1 + tmp281_1 + tmp282_1 + tmp283_1 + tmp286_1 + tmp289_1 + tmp358_1 + tmp359_1 + tmp360_1 + tmp361_1;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=tmp362_1 + tmp363_1 + tmp364_1 + tmp365_1 + tmp366_1 + tmp367_1 + tmp65_1 + tmp66_1 + tmp67_1 + tmp68_1 + tmp69_1 + tmp71_1 + tmp72_1 + tmp75_1;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=tmp119_1 + tmp120_1 + tmp142_1 + tmp143_1 + tmp144_1 + tmp146_1 + tmp147_1 + tmp149_1 + tmp151_1 + tmp152_1 + tmp368_1 + tmp369_1 + tmp370_1 + tmp371_1;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=tmp168_1 + tmp170_1 + tmp330_1 + tmp335_1 + tmp345_1 + tmp352_1 + tmp372_1 + tmp373_1 + tmp374_1 + tmp375_1 + tmp376_1 + tmp377_1 + tmp378_1 + tmp379_1 + tmp380_1 + tmp381_1 + tmp382_1 + tmp383_1;
                                        EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp294_1 + tmp295_1 + tmp354_1 + tmp355_1 + tmp356_1 + tmp357_1 + tmp384_1 + tmp385_1 + tmp386_1 + tmp387_1;
                                        EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp388_1 + tmp389_1 + tmp390_1 + tmp391_1 + tmp392_1 + tmp393_1 + tmp394_1 + tmp395_1 + tmp396_1 + tmp397_1;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=tmp195_1 + tmp198_1 + tmp266_1 + tmp268_1 + tmp269_1 + tmp271_1 + tmp274_1 + tmp275_1 + tmp278_1 + tmp279_1 + tmp398_1 + tmp399_1 + tmp400_1 + tmp401_1;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=tmp313_1 + tmp314_1 + tmp315_1 + tmp317_1 + tmp318_1 + tmp319_1 + tmp320_1 + tmp323_1 + tmp362_1 + tmp363_1 + tmp402_1 + tmp403_1 + tmp404_1 + tmp405_1;
                                        EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp124_1 + tmp125_1 + tmp126_1 + tmp127_1 + tmp54_1 + tmp55_1 + tmp57_1 + tmp58_1 + tmp59_1 + tmp61_1;
                                        EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp388_1 + tmp389_1 + tmp390_1 + tmp391_1 + tmp393_1 + tmp397_1 + tmp406_1 + tmp407_1 + tmp408_1 + tmp409_1;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=tmp200_1 + tmp202_1 + tmp410_1 + tmp411_1 + tmp412_1 + tmp413_1 + tmp414_1 + tmp415_1 + tmp416_1 + tmp417_1 + tmp418_1 + tmp419_1 + tmp420_1 + tmp421_1;
                                        EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp14_1 + tmp15_1 + tmp16_1 + tmp19_1 + tmp20_1 + tmp23_1 + tmp422_1 + tmp423_1 + tmp424_1 + tmp425_1;
                                        EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp280_1 + tmp281_1 + tmp358_1 + tmp359_1 + tmp360_1 + tmp361_1 + tmp426_1 + tmp427_1 + tmp428_1 + tmp429_1;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=tmp154_1 + tmp155_1 + tmp411_1 + tmp412_1 + tmp413_1 + tmp414_1 + tmp416_1 + tmp417_1 + tmp419_1 + tmp421_1 + tmp430_1 + tmp431_1 + tmp432_1 + tmp433_1;
                                        EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp100_1 + tmp101_1 + tmp122_1 + tmp123_1 + tmp434_1 + tmp435_1;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=tmp27_1 + tmp28_1 + tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp436_1 + tmp437_1 + tmp438_1 + tmp439_1 + tmp43_1 + tmp44_1 + tmp50_1 + tmp51_1;
                                        EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp100_1 + tmp101_1 + tmp434_1 + tmp435_1 + tmp78_1 + tmp79_1;
                                        EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp122_1 + tmp123_1 + tmp434_1 + tmp435_1 + tmp76_1 + tmp77_1;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=tmp24_1 + tmp26_1 + tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp34_1 + tmp35_1 + tmp440_1 + tmp441_1 + tmp442_1 + tmp443_1 + tmp46_1 + tmp48_1;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=tmp186_1 + tmp187_1 + tmp188_1 + tmp189_1 + tmp191_1 + tmp192_1 + tmp194_1 + tmp197_1 + tmp272_1 + tmp273_1 + tmp444_1 + tmp445_1 + tmp446_1 + tmp447_1;
                                        EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp19_1 + tmp23_1 + tmp422_1 + tmp423_1 + tmp424_1 + tmp425_1 + tmp82_1 + tmp83_1 + tmp84_1 + tmp85_1;
                                        EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp388_1 + tmp389_1 + tmp406_1 + tmp407_1 + tmp408_1 + tmp409_1 + tmp448_1 + tmp449_1 + tmp450_1 + tmp451_1;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=tmp109_1 + tmp112_1 + tmp452_1 + tmp453_1 + tmp454_1 + tmp455_1 + tmp87_1 + tmp88_1 + tmp89_1 + tmp91_1 + tmp94_1 + tmp95_1 + tmp97_1 + tmp99_1;
                                        EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp434_1 + tmp435_1 + tmp76_1 + tmp77_1 + tmp78_1 + tmp79_1;
                                        EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp280_1 + tmp281_1 + tmp284_1 + tmp285_1 + tmp287_1 + tmp288_1 + tmp426_1 + tmp427_1 + tmp428_1 + tmp429_1;
                                        EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp290_1 + tmp291_1 + tmp294_1 + tmp295_1 + tmp297_1 + tmp299_1 + tmp384_1 + tmp385_1 + tmp386_1 + tmp387_1;
                                        EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp100_1 + tmp101_1 + tmp122_1 + tmp123_1 + tmp80_1 + tmp81_1;
                                        EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp388_1 + tmp389_1 + tmp392_1 + tmp394_1 + tmp395_1 + tmp396_1 + tmp448_1 + tmp449_1 + tmp450_1 + tmp451_1;
                                    }
                                }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const register double B_0 = B_p[INDEX3(k,0,m, numEq, 3)];
                                        const register double B_1 = B_p[INDEX3(k,1,m, numEq, 3)];
                                        const register double B_2 = B_p[INDEX3(k,2,m, numEq, 3)];
                                        const register double tmp7_1 = B_2*w133;
                                        const register double tmp17_1 = B_0*w143;
                                        const register double tmp1_1 = B_1*w127;
                                        const register double tmp8_1 = B_2*w135;
                                        const register double tmp15_1 = B_0*w141;
                                        const register double tmp3_1 = B_0*w129;
                                        const register double tmp13_1 = B_1*w139;
                                        const register double tmp2_1 = B_0*w126;
                                        const register double tmp12_1 = B_0*w138;
                                        const register double tmp5_1 = B_2*w131;
                                        const register double tmp10_1 = B_2*w136;
                                        const register double tmp4_1 = B_1*w130;
                                        const register double tmp6_1 = B_1*w132;
                                        const register double tmp9_1 = B_1*w134;
                                        const register double tmp14_1 = B_2*w140;
                                        const register double tmp0_1 = B_2*w128;
                                        const register double tmp16_1 = B_1*w142;
                                        const register double tmp11_1 = B_0*w137;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                        EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=tmp2_1 + tmp6_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=tmp3_1 + tmp8_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp10_1 + tmp1_1 + tmp2_1;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=tmp11_1 + tmp4_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp10_1 + tmp12_1 + tmp13_1;
                                        EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp2_1 + tmp4_1 + tmp5_1;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=tmp14_1 + tmp1_1 + tmp3_1;
                                        EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp10_1 + tmp13_1 + tmp15_1;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp4_1;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=tmp14_1 + tmp1_1 + tmp2_1;
                                        EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp10_1 + tmp12_1 + tmp16_1;
                                        EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp10_1 + tmp3_1 + tmp4_1;
                                        EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp16_1 + tmp2_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=tmp0_1 + tmp3_1 + tmp4_1;
                                        EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp10_1 + tmp2_1 + tmp4_1;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=tmp17_1 + tmp1_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=tmp11_1 + tmp14_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=tmp3_1 + tmp7_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=tmp11_1 + tmp1_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=tmp0_1 + tmp17_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp16_1 + tmp3_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=tmp14_1 + tmp2_1 + tmp4_1;
                                        EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp16_1 + tmp2_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=tmp14_1 + tmp17_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=tmp0_1 + tmp11_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=tmp2_1 + tmp6_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp15_1 + tmp4_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp13_1 + tmp3_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=tmp11_1 + tmp14_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp16_1 + tmp3_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=tmp11_1 + tmp1_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=tmp14_1 + tmp17_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=tmp0_1 + tmp17_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp13_1 + tmp2_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp15_1 + tmp1_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=tmp17_1 + tmp4_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=tmp14_1 + tmp3_1 + tmp4_1;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=tmp0_1 + tmp11_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp13_1 + tmp2_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp12_1 + tmp1_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=tmp2_1 + tmp8_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=tmp17_1 + tmp1_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp10_1 + tmp1_1 + tmp3_1;
                                        EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp12_1 + tmp4_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=tmp11_1 + tmp4_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp1_1 + tmp3_1 + tmp5_1;
                                        EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp15_1 + tmp1_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=tmp17_1 + tmp4_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp15_1 + tmp16_1 + tmp5_1;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=tmp3_1 + tmp6_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp13_1 + tmp15_1 + tmp5_1;
                                        EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp12_1 + tmp16_1 + tmp5_1;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=tmp2_1 + tmp7_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=tmp3_1 + tmp6_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp5_1;
                                        EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp12_1 + tmp4_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp3_1;
                                        EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp12_1 + tmp13_1 + tmp5_1;
                                        EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp15_1 + tmp4_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp13_1 + tmp3_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp10_1 + tmp15_1 + tmp16_1;
                                        EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp12_1 + tmp1_1 + tmp7_1;
                                    }
                                }
                            }
                        }
                        ///////////////
                        // process C //
                        ///////////////
                        if (!C.isEmpty()) {
                            add_EM_S=true;
                            const double* C_p=const_cast<escript::Data*>(&C)->getSampleDataRO(e);
                            if (C.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const register double C_0_0 = C_p[INDEX4(k,m,0, 0, numEq,numComp,3)];
                                        const register double C_1_0 = C_p[INDEX4(k,m,1, 0, numEq,numComp,3)];
                                        const register double C_2_0 = C_p[INDEX4(k,m,2, 0, numEq,numComp,3)];
                                        const register double C_0_1 = C_p[INDEX4(k,m,0, 1, numEq,numComp,3)];
                                        const register double C_1_1 = C_p[INDEX4(k,m,1, 1, numEq,numComp,3)];
                                        const register double C_2_1 = C_p[INDEX4(k,m,2, 1, numEq,numComp,3)];
                                        const register double C_0_2 = C_p[INDEX4(k,m,0, 2, numEq,numComp,3)];
                                        const register double C_1_2 = C_p[INDEX4(k,m,1, 2, numEq,numComp,3)];
                                        const register double C_2_2 = C_p[INDEX4(k,m,2, 2, numEq,numComp,3)];
                                        const register double C_0_3 = C_p[INDEX4(k,m,0, 3, numEq,numComp,3)];
                                        const register double C_1_3 = C_p[INDEX4(k,m,1, 3, numEq,numComp,3)];
                                        const register double C_2_3 = C_p[INDEX4(k,m,2, 3, numEq,numComp,3)];
                                        const register double C_0_4 = C_p[INDEX4(k,m,0, 4, numEq,numComp,3)];
                                        const register double C_1_4 = C_p[INDEX4(k,m,1, 4, numEq,numComp,3)];
                                        const register double C_2_4 = C_p[INDEX4(k,m,2, 4, numEq,numComp,3)];
                                        const register double C_0_5 = C_p[INDEX4(k,m,0, 5, numEq,numComp,3)];
                                        const register double C_1_5 = C_p[INDEX4(k,m,1, 5, numEq,numComp,3)];
                                        const register double C_2_5 = C_p[INDEX4(k,m,2, 5, numEq,numComp,3)];
                                        const register double C_0_6 = C_p[INDEX4(k,m,0, 6, numEq,numComp,3)];
                                        const register double C_1_6 = C_p[INDEX4(k,m,1, 6, numEq,numComp,3)];
                                        const register double C_2_6 = C_p[INDEX4(k,m,2, 6, numEq,numComp,3)];
                                        const register double C_0_7 = C_p[INDEX4(k,m,0, 7, numEq,numComp,3)];
                                        const register double C_1_7 = C_p[INDEX4(k,m,1, 7, numEq,numComp,3)];
                                        const register double C_2_7 = C_p[INDEX4(k,m,2, 7, numEq,numComp,3)];
                                        const register double tmp11_0 = C_0_5 + C_0_7;
                                        const register double tmp39_0 = C_0_2 + C_0_4;
                                        const register double tmp40_0 = C_1_1 + C_1_4;
                                        const register double tmp14_0 = C_0_4 + C_0_6;
                                        const register double tmp34_0 = C_1_0 + C_1_1 + C_1_4 + C_1_5;
                                        const register double tmp23_0 = C_1_0 + C_1_5;
                                        const register double tmp12_0 = C_1_4 + C_1_5;
                                        const register double tmp35_0 = C_1_2 + C_1_3 + C_1_6 + C_1_7;
                                        const register double tmp16_0 = C_1_6 + C_1_7;
                                        const register double tmp19_0 = C_2_0 + C_2_1 + C_2_2 + C_2_3;
                                        const register double tmp2_0 = C_1_3 + C_1_7;
                                        const register double tmp31_0 = C_2_4 + C_2_5;
                                        const register double tmp22_0 = C_2_0 + C_2_2;
                                        const register double tmp37_0 = C_2_4 + C_2_7;
                                        const register double tmp6_0 = C_2_1 + C_2_2;
                                        const register double tmp29_0 = C_2_0 + C_2_1;
                                        const register double tmp27_0 = C_0_1 + C_0_7;
                                        const register double tmp1_0 = C_0_2 + C_0_6;
                                        const register double tmp4_0 = C_0_3 + C_0_7;
                                        const register double tmp38_0 = C_0_3 + C_0_5;
                                        const register double tmp26_0 = C_0_0 + C_0_6;
                                        const register double tmp33_0 = C_0_0 + C_0_2 + C_0_4 + C_0_6;
                                        const register double tmp10_0 = C_0_0 + C_0_2;
                                        const register double tmp25_0 = C_1_2 + C_1_7;
                                        const register double tmp7_0 = C_1_1 + C_1_5;
                                        const register double tmp18_0 = C_1_0 + C_1_1;
                                        const register double tmp8_0 = C_0_0 + C_0_4;
                                        const register double tmp13_0 = C_2_4 + C_2_5 + C_2_6 + C_2_7;
                                        const register double tmp28_0 = C_2_2 + C_2_3;
                                        const register double tmp17_0 = C_0_1 + C_0_3;
                                        const register double tmp36_0 = C_2_0 + C_2_3;
                                        const register double tmp9_0 = C_1_2 + C_1_6;
                                        const register double tmp20_0 = C_2_1 + C_2_3;
                                        const register double tmp41_0 = C_1_3 + C_1_6;
                                        const register double tmp5_0 = C_1_0 + C_1_4;
                                        const register double tmp3_0 = C_2_5 + C_2_6;
                                        const register double tmp15_0 = C_1_2 + C_1_3;
                                        const register double tmp24_0 = C_2_4 + C_2_6;
                                        const register double tmp32_0 = C_0_1 + C_0_3 + C_0_5 + C_0_7;
                                        const register double tmp0_0 = C_0_1 + C_0_5;
                                        const register double tmp21_0 = C_2_5 + C_2_7;
                                        const register double tmp30_0 = C_2_6 + C_2_7;
                                        const register double tmp28_1 = tmp21_0*w117;
                                        const register double tmp404_1 = C_0_4*w90;
                                        const register double tmp153_1 = tmp5_0*w111;
                                        const register double tmp208_1 = C_2_0*w103;
                                        const register double tmp201_1 = tmp38_0*w108;
                                        const register double tmp102_1 = tmp4_0*w93;
                                        const register double tmp277_1 = tmp41_0*w97;
                                        const register double tmp39_1 = tmp11_0*w108;
                                        const register double tmp443_1 = C_1_1*w111;
                                        const register double tmp61_1 = tmp15_0*w111;
                                        const register double tmp188_1 = tmp17_0*w108;
                                        const register double tmp139_1 = tmp35_0*w94;
                                        const register double tmp378_1 = C_0_4*w98;
                                        const register double tmp10_1 = tmp7_0*w94;
                                        const register double tmp31_1 = tmp22_0*w112;
                                        const register double tmp209_1 = C_1_4*w114;
                                        const register double tmp369_1 = C_2_0*w100;
                                        const register double tmp213_1 = C_2_4*w100;
                                        const register double tmp254_1 = tmp38_0*w96;
                                        const register double tmp350_1 = C_0_1*w124;
                                        const register double tmp183_1 = C_0_4*w120;
                                        const register double tmp182_1 = C_0_5*w90;
                                        const register double tmp83_1 = tmp14_0*w108;
                                        const register double tmp247_1 = tmp41_0*w94;
                                        const register double tmp79_1 = tmp34_0*w94;
                                        const register double tmp332_1 = C_1_7*w115;
                                        const register double tmp75_1 = tmp31_0*w103;
                                        const register double tmp303_1 = C_1_1*w114;
                                        const register double tmp304_1 = C_2_5*w113;
                                        const register double tmp420_1 = C_0_1*w125;
                                        const register double tmp271_1 = tmp20_0*w95;
                                        const register double tmp444_1 = C_1_2*w99;
                                        const register double tmp381_1 = C_1_0*w105;
                                        const register double tmp96_1 = tmp5_0*w94;
                                        const register double tmp385_1 = tmp30_0*w92;
                                        const register double tmp371_1 = C_2_4*w103;
                                        const register double tmp66_1 = tmp18_0*w107;
                                        const register double tmp300_1 = C_0_1*w121;
                                        const register double tmp379_1 = C_2_2*w101;
                                        const register double tmp283_1 = tmp21_0*w113;
                                        const register double tmp242_1 = C_1_7*w91;
                                        const register double tmp234_1 = tmp30_0*w95;
                                        const register double tmp51_1 = tmp24_0*w103;
                                        const register double tmp151_1 = tmp2_0*w105;
                                        const register double tmp46_1 = tmp22_0*w95;
                                        const register double tmp120_1 = C_2_4*w101;
                                        const register double tmp38_1 = tmp10_0*w106;
                                        const register double tmp413_1 = tmp15_0*w109;
                                        const register double tmp173_1 = C_2_7*w113;
                                        const register double tmp116_1 = C_2_7*w103;
                                        const register double tmp221_1 = C_0_4*w125;
                                        const register double tmp158_1 = tmp29_0*w92;
                                        const register double tmp30_1 = C_1_4*w91;
                                        const register double tmp154_1 = tmp38_0*w93;
                                        const register double tmp406_1 = tmp5_0*w99;
                                        const register double tmp204_1 = C_0_0*w124;
                                        const register double tmp452_1 = C_2_1*w122;
                                        const register double tmp117_1 = C_2_0*w92;
                                        const register double tmp323_1 = C_0_3*w125;
                                        const register double tmp225_1 = tmp4_0*w96;
                                        const register double tmp274_1 = tmp40_0*w94;
                                        const register double tmp71_1 = C_0_2*w124;
                                        const register double tmp43_1 = tmp23_0*w107;
                                        const register double tmp408_1 = tmp9_0*w94;
                                        const register double tmp427_1 = tmp24_0*w92;
                                        const register double tmp82_1 = tmp17_0*w106;
                                        const register double tmp168_1 = tmp26_0*w93;
                                        const register double tmp276_1 = tmp10_0*w90;
                                        const register double tmp336_1 = C_2_3*w116;
                                        const register double tmp53_1 = tmp11_0*w106;
                                        const register double tmp161_1 = C_0_7*w98;
                                        const register double tmp382_1 = C_0_3*w90;
                                        const register double tmp27_1 = tmp20_0*w113;
                                        const register double tmp73_1 = C_0_5*w125;
                                        const register double tmp389_1 = tmp33_0*w108;
                                        const register double tmp288_1 = tmp22_0*w116;
                                        const register double tmp258_1 = C_2_7*w100;
                                        const register double tmp295_1 = tmp29_0*w117;
                                        const register double tmp64_1 = C_0_4*w110;
                                        const register double tmp250_1 = tmp40_0*w97;
                                        const register double tmp338_1 = C_0_0*w125;
                                        const register double tmp359_1 = tmp7_0*w109;
                                        const register double tmp400_1 = C_1_5*w111;
                                        const register double tmp414_1 = tmp29_0*w113;
                                        const register double tmp214_1 = C_1_3*w115;
                                        const register double tmp237_1 = C_0_4*w104;
                                        const register double tmp7_1 = tmp4_0*w98;
                                        const register double tmp145_1 = C_2_3*w122;
                                        const register double tmp129_1 = tmp15_0*w94;
                                        const register double tmp45_1 = C_1_1*w115;
                                        const register double tmp124_1 = tmp14_0*w93;
                                        const register double tmp284_1 = tmp20_0*w117;
                                        const register double tmp148_1 = C_2_0*w113;
                                        const register double tmp334_1 = C_0_7*w124;
                                        const register double tmp257_1 = C_1_5*w99;
                                        const register double tmp186_1 = tmp11_0*w104;
                                        const register double tmp344_1 = C_2_6*w100;
                                        const register double tmp433_1 = C_0_7*w120;
                                        const register double tmp169_1 = C_0_3*w121;
                                        const register double tmp262_1 = C_1_7*w119;
                                        const register double tmp47_1 = tmp17_0*w104;
                                        const register double tmp278_1 = C_1_2*w118;
                                        const register double tmp217_1 = C_2_3*w101;
                                        const register double tmp33_1 = C_1_1*w119;
                                        const register double tmp157_1 = tmp30_0*w103;
                                        const register double tmp44_1 = tmp25_0*w109;
                                        const register double tmp454_1 = C_2_6*w123;
                                        const register double tmp205_1 = C_0_7*w125;
                                        const register double tmp373_1 = C_1_2*w115;
                                        const register double tmp361_1 = tmp2_0*w111;
                                        const register double tmp392_1 = tmp24_0*w113;
                                        const register double tmp308_1 = C_2_2*w116;
                                        const register double tmp215_1 = tmp3_0*w95;
                                        const register double tmp417_1 = C_0_6*w124;
                                        const register double tmp416_1 = tmp28_0*w112;
                                        const register double tmp210_1 = C_2_7*w92;
                                        const register double tmp335_1 = tmp41_0*w107;
                                        const register double tmp429_1 = tmp22_0*w102;
                                        const register double tmp315_1 = tmp28_0*w113;
                                        const register double tmp435_1 = tmp19_0*w112;
                                        const register double tmp407_1 = tmp2_0*w91;
                                        const register double tmp365_1 = C_0_5*w98;
                                        const register double tmp136_1 = tmp30_0*w113;
                                        const register double tmp218_1 = C_1_6*w111;
                                        const register double tmp171_1 = C_2_4*w122;
                                        const register double tmp253_1 = C_0_7*w121;
                                        const register double tmp123_1 = tmp35_0*w109;
                                        const register double tmp149_1 = C_2_4*w123;
                                        const register double tmp291_1 = tmp0_0*w96;
                                        const register double tmp267_1 = tmp14_0*w96;
                                        const register double tmp147_1 = tmp9_0*w109;
                                        const register double tmp296_1 = tmp8_0*w98;
                                        const register double tmp42_1 = tmp20_0*w92;
                                        const register double tmp430_1 = C_0_0*w121;
                                        const register double tmp399_1 = C_1_7*w114;
                                        const register double tmp355_1 = tmp4_0*w106;
                                        const register double tmp56_1 = tmp16_0*w107;
                                        const register double tmp289_1 = tmp5_0*w97;
                                        const register double tmp125_1 = tmp17_0*w96;
                                        const register double tmp397_1 = tmp9_0*w111;
                                        const register double tmp333_1 = C_2_0*w123;
                                        const register double tmp15_1 = tmp11_0*w96;
                                        const register double tmp200_1 = tmp39_0*w106;
                                        const register double tmp285_1 = tmp9_0*w91;
                                        const register double tmp374_1 = C_2_1*w103;
                                        const register double tmp354_1 = tmp1_0*w104;
                                        const register double tmp97_1 = tmp2_0*w97;
                                        const register double tmp104_1 = C_2_2*w122;
                                        const register double tmp339_1 = C_1_2*w111;
                                        const register double tmp437_1 = C_1_1*w91;
                                        const register double tmp259_1 = C_0_6*w98;
                                        const register double tmp110_1 = tmp0_0*w98;
                                        const register double tmp419_1 = tmp30_0*w116;
                                        const register double tmp236_1 = C_1_5*w118;
                                        const register double tmp41_1 = C_1_6*w114;
                                        const register double tmp320_1 = C_0_4*w124;
                                        const register double tmp74_1 = tmp12_0*w111;
                                        const register double tmp63_1 = tmp27_0*w108;
                                        const register double tmp106_1 = tmp2_0*w109;
                                        const register double tmp150_1 = tmp4_0*w104;
                                        const register double tmp195_1 = tmp20_0*w112;
                                        const register double tmp270_1 = C_1_7*w99;
                                        const register double tmp230_1 = C_2_1*w100;
                                        const register double tmp68_1 = tmp29_0*w95;
                                        const register double tmp287_1 = tmp2_0*w94;
                                        const register double tmp132_1 = tmp28_0*w117;
                                        const register double tmp423_1 = tmp12_0*w109;
                                        const register double tmp272_1 = tmp11_0*w98;
                                        const register double tmp134_1 = tmp1_0*w108;
                                        const register double tmp1_1 = tmp1_0*w96;
                                        const register double tmp301_1 = C_1_6*w115;
                                        const register double tmp241_1 = C_2_6*w113;
                                        const register double tmp29_1 = tmp17_0*w98;
                                        const register double tmp337_1 = C_1_5*w105;
                                        const register double tmp92_1 = tmp36_0*w95;
                                        const register double tmp40_1 = tmp14_0*w110;
                                        const register double tmp203_1 = C_0_1*w104;
                                        const register double tmp349_1 = C_2_1*w101;
                                        const register double tmp78_1 = tmp19_0*w95;
                                        const register double tmp266_1 = tmp17_0*w93;
                                        const register double tmp58_1 = tmp19_0*w117;
                                        const register double tmp310_1 = C_0_6*w120;
                                        const register double tmp127_1 = tmp10_0*w98;
                                        const register double tmp307_1 = C_0_7*w90;
                                        const register double tmp330_1 = tmp40_0*w109;
                                        const register double tmp386_1 = tmp31_0*w95;
                                        const register double tmp351_1 = C_1_1*w118;
                                        const register double tmp306_1 = C_0_0*w98;
                                        const register double tmp312_1 = tmp30_0*w117;
                                        const register double tmp108_1 = C_2_5*w123;
                                        const register double tmp140_1 = tmp29_0*w116;
                                        const register double tmp235_1 = tmp29_0*w102;
                                        const register double tmp293_1 = tmp34_0*w109;
                                        const register double tmp446_1 = C_1_0*w119;
                                        const register double tmp160_1 = tmp28_0*w95;
                                        const register double tmp393_1 = tmp22_0*w117;
                                        const register double tmp231_1 = C_2_6*w101;
                                        const register double tmp233_1 = tmp31_0*w92;
                                        const register double tmp8_1 = tmp5_0*w91;
                                        const register double tmp229_1 = C_2_2*w92;
                                        const register double tmp128_1 = tmp16_0*w91;
                                        const register double tmp187_1 = tmp14_0*w106;
                                        const register double tmp86_1 = tmp8_0*w106;
                                        const register double tmp49_1 = C_1_3*w105;
                                        const register double tmp447_1 = C_1_7*w118;
                                        const register double tmp442_1 = C_1_6*w105;
                                        const register double tmp77_1 = tmp33_0*w96;
                                        const register double tmp358_1 = tmp9_0*w107;
                                        const register double tmp292_1 = tmp35_0*w107;
                                        const register double tmp81_1 = tmp13_0*w102;
                                        const register double tmp206_1 = tmp27_0*w106;
                                        const register double tmp366_1 = C_0_2*w90;
                                        const register double tmp448_1 = tmp22_0*w103;
                                        const register double tmp156_1 = tmp39_0*w96;
                                        const register double tmp105_1 = tmp5_0*w107;
                                        const register double tmp179_1 = tmp25_0*w94;
                                        const register double tmp364_1 = tmp26_0*w96;
                                        const register double tmp449_1 = tmp21_0*w92;
                                        const register double tmp453_1 = C_2_2*w113;
                                        const register double tmp12_1 = tmp8_0*w90;
                                        const register double tmp59_1 = tmp13_0*w112;
                                        const register double tmp99_1 = C_2_6*w103;
                                        const register double tmp327_1 = C_0_1*w110;
                                        const register double tmp50_1 = C_1_4*w111;
                                        const register double tmp126_1 = tmp18_0*w99;
                                        const register double tmp445_1 = C_1_5*w91;
                                        const register double tmp152_1 = C_2_7*w116;
                                        const register double tmp76_1 = tmp32_0*w93;
                                        const register double tmp80_1 = tmp35_0*w97;
                                        const register double tmp146_1 = tmp7_0*w107;
                                        const register double tmp197_1 = tmp21_0*w116;
                                        const register double tmp111_1 = tmp36_0*w112;
                                        const register double tmp115_1 = tmp7_0*w111;
                                        const register double tmp396_1 = tmp7_0*w105;
                                        const register double tmp248_1 = C_1_2*w119;
                                        const register double tmp390_1 = tmp2_0*w107;
                                        const register double tmp119_1 = tmp6_0*w95;
                                        const register double tmp440_1 = C_1_4*w115;
                                        const register double tmp6_1 = tmp3_0*w117;
                                        const register double tmp346_1 = C_0_6*w125;
                                        const register double tmp144_1 = tmp8_0*w110;
                                        const register double tmp249_1 = C_2_1*w116;
                                        const register double tmp318_1 = tmp29_0*w112;
                                        const register double tmp265_1 = C_1_0*w118;
                                        const register double tmp387_1 = tmp28_0*w102;
                                        const register double tmp130_1 = tmp11_0*w90;
                                        const register double tmp348_1 = C_1_3*w91;
                                        const register double tmp252_1 = tmp39_0*w93;
                                        const register double tmp85_1 = tmp10_0*w104;
                                        const register double tmp21_1 = tmp17_0*w90;
                                        const register double tmp90_1 = tmp9_0*w99;
                                        const register double tmp232_1 = tmp28_0*w103;
                                        const register double tmp441_1 = C_1_3*w114;
                                        const register double tmp360_1 = tmp5_0*w105;
                                        const register double tmp4_1 = C_2_3*w113;
                                        const register double tmp313_1 = C_0_2*w110;
                                        const register double tmp135_1 = tmp4_0*w110;
                                        const register double tmp155_1 = C_0_6*w121;
                                        const register double tmp20_1 = tmp16_0*w94;
                                        const register double tmp432_1 = C_0_6*w90;
                                        const register double tmp227_1 = tmp0_0*w90;
                                        const register double tmp269_1 = tmp22_0*w92;
                                        const register double tmp261_1 = C_2_0*w101;
                                        const register double tmp19_1 = tmp15_0*w91;
                                        const register double tmp16_1 = tmp12_0*w99;
                                        const register double tmp395_1 = tmp20_0*w116;
                                        const register double tmp356_1 = tmp8_0*w108;
                                        const register double tmp222_1 = tmp6_0*w102;
                                        const register double tmp281_1 = tmp32_0*w96;
                                        const register double tmp164_1 = C_0_0*w90;
                                        const register double tmp375_1 = C_1_5*w114;
                                        const register double tmp238_1 = C_0_3*w110;
                                        const register double tmp340_1 = C_2_2*w103;
                                        const register double tmp290_1 = tmp1_0*w93;
                                        const register double tmp311_1 = C_1_3*w111;
                                        const register double tmp162_1 = tmp18_0*w91;
                                        const register double tmp107_1 = C_2_1*w113;
                                        const register double tmp401_1 = C_1_2*w105;
                                        const register double tmp87_1 = tmp4_0*w108;
                                        const register double tmp172_1 = C_1_1*w99;
                                        const register double tmp402_1 = C_0_2*w121;
                                        const register double tmp434_1 = tmp13_0*w117;
                                        const register double tmp394_1 = tmp21_0*w112;
                                        const register double tmp191_1 = tmp41_0*w109;
                                        const register double tmp325_1 = tmp38_0*w106;
                                        const register double tmp91_1 = C_2_2*w100;
                                        const register double tmp35_1 = tmp14_0*w90;
                                        const register double tmp279_1 = tmp24_0*w102;
                                        const register double tmp436_1 = C_1_6*w99;
                                        const register double tmp244_1 = tmp36_0*w117;
                                        const register double tmp264_1 = C_0_0*w120;
                                        const register double tmp415_1 = tmp31_0*w117;
                                        const register double tmp178_1 = tmp3_0*w112;
                                        const register double tmp239_1 = C_2_5*w122;
                                        const register double tmp219_1 = C_0_3*w124;
                                        const register double tmp370_1 = C_2_7*w101;
                                        const register double tmp109_1 = tmp37_0*w117;
                                        const register double tmp316_1 = C_0_5*w104;
                                        const register double tmp174_1 = C_2_3*w123;
                                        const register double tmp89_1 = C_2_1*w92;
                                        const register double tmp228_1 = C_2_5*w103;
                                        const register double tmp36_1 = tmp25_0*w97;
                                        const register double tmp412_1 = tmp12_0*w107;
                                        const register double tmp326_1 = tmp39_0*w108;
                                        const register double tmp256_1 = C_2_4*w92;
                                        const register double tmp23_1 = tmp19_0*w102;
                                        const register double tmp362_1 = tmp27_0*w93;
                                        const register double tmp286_1 = tmp24_0*w112;
                                        const register double tmp114_1 = C_2_6*w116;
                                        const register double tmp176_1 = C_0_2*w98;
                                        const register double tmp84_1 = tmp11_0*w110;
                                        const register double tmp328_1 = C_2_7*w122;
                                        const register double tmp113_1 = tmp1_0*w90;
                                        const register double tmp55_1 = tmp17_0*w110;
                                        const register double tmp319_1 = tmp18_0*w94;
                                        const register double tmp95_1 = C_2_5*w101;
                                        const register double tmp72_1 = tmp15_0*w105;
                                        const register double tmp220_1 = C_1_1*w105;
                                        const register double tmp60_1 = tmp12_0*w105;
                                        const register double tmp118_1 = C_2_3*w100;
                                        const register double tmp32_1 = tmp23_0*w94;
                                        const register double tmp166_1 = tmp15_0*w97;
                                        const register double tmp159_1 = tmp16_0*w99;
                                        const register double tmp421_1 = tmp18_0*w111;
                                        const register double tmp367_1 = C_0_3*w120;
                                        const register double tmp424_1 = tmp16_0*w111;
                                        const register double tmp431_1 = C_0_1*w98;
                                        const register double tmp165_1 = C_0_1*w120;
                                        const register double tmp309_1 = C_1_4*w105;
                                        const register double tmp368_1 = C_2_3*w92;
                                        const register double tmp223_1 = C_0_5*w110;
                                        const register double tmp384_1 = tmp29_0*w103;
                                        const register double tmp143_1 = tmp0_0*w108;
                                        const register double tmp240_1 = C_1_0*w99;
                                        const register double tmp255_1 = C_2_3*w103;
                                        const register double tmp122_1 = tmp34_0*w107;
                                        const register double tmp388_1 = tmp32_0*w106;
                                        const register double tmp185_1 = C_1_4*w118;
                                        const register double tmp88_1 = tmp1_0*w110;
                                        const register double tmp342_1 = C_2_5*w92;
                                        const register double tmp196_1 = C_1_2*w114;
                                        const register double tmp142_1 = tmp1_0*w106;
                                        const register double tmp131_1 = tmp12_0*w97;
                                        const register double tmp418_1 = tmp16_0*w105;
                                        const register double tmp263_1 = C_0_1*w90;
                                        const register double tmp17_1 = tmp13_0*w95;
                                        const register double tmp422_1 = tmp15_0*w107;
                                        const register double tmp121_1 = tmp3_0*w102;
                                        const register double tmp190_1 = tmp40_0*w107;
                                        const register double tmp70_1 = tmp30_0*w102;
                                        const register double tmp363_1 = C_0_4*w121;
                                        const register double tmp194_1 = tmp24_0*w117;
                                        const register double tmp297_1 = tmp30_0*w112;
                                        const register double tmp280_1 = tmp33_0*w93;
                                        const register double tmp391_1 = tmp5_0*w109;
                                        const register double tmp439_1 = C_1_3*w118;
                                        const register double tmp69_1 = C_0_3*w104;
                                        const register double tmp48_1 = tmp21_0*w102;
                                        const register double tmp243_1 = C_2_2*w123;
                                        const register double tmp405_1 = C_0_5*w120;
                                        const register double tmp2_1 = C_2_0*w122;
                                        const register double tmp372_1 = C_0_5*w121;
                                        const register double tmp100_1 = tmp33_0*w106;
                                        const register double tmp352_1 = tmp36_0*w102;
                                        const register double tmp317_1 = tmp12_0*w91;
                                        const register double tmp141_1 = tmp34_0*w97;
                                        const register double tmp282_1 = tmp7_0*w99;
                                        const register double tmp181_1 = C_2_0*w116;
                                        const register double tmp177_1 = C_1_6*w91;
                                        const register double tmp14_1 = tmp10_0*w93;
                                        const register double tmp321_1 = tmp31_0*w116;
                                        const register double tmp98_1 = tmp37_0*w102;
                                        const register double tmp9_1 = tmp6_0*w112;
                                        const register double tmp322_1 = tmp16_0*w97;
                                        const register double tmp103_1 = tmp8_0*w96;
                                        const register double tmp0_1 = tmp0_0*w93;
                                        const register double tmp3_1 = tmp2_0*w99;
                                        const register double tmp216_1 = C_0_2*w104;
                                        const register double tmp11_1 = C_2_4*w116;
                                        const register double tmp189_1 = tmp10_0*w110;
                                        const register double tmp275_1 = C_1_5*w119;
                                        const register double tmp376_1 = C_2_6*w92;
                                        const register double tmp170_1 = tmp27_0*w96;
                                        const register double tmp101_1 = tmp32_0*w108;
                                        const register double tmp353_1 = C_0_7*w110;
                                        const register double tmp13_1 = tmp9_0*w97;
                                        const register double tmp133_1 = tmp0_0*w106;
                                        const register double tmp251_1 = C_0_2*w125;
                                        const register double tmp62_1 = tmp26_0*w106;
                                        const register double tmp137_1 = tmp8_0*w104;
                                        const register double tmp294_1 = tmp31_0*w113;
                                        const register double tmp260_1 = C_1_2*w91;
                                        const register double tmp409_1 = tmp7_0*w97;
                                        const register double tmp24_1 = tmp11_0*w93;
                                        const register double tmp314_1 = tmp15_0*w99;
                                        const register double tmp26_1 = C_1_3*w99;
                                        const register double tmp211_1 = tmp25_0*w107;
                                        const register double tmp331_1 = C_2_4*w113;
                                        const register double tmp199_1 = C_1_0*w111;
                                        const register double tmp22_1 = tmp18_0*w97;
                                        const register double tmp34_1 = tmp24_0*w116;
                                        const register double tmp192_1 = tmp22_0*w113;
                                        const register double tmp67_1 = tmp16_0*w109;
                                        const register double tmp426_1 = tmp20_0*w103;
                                        const register double tmp226_1 = tmp1_0*w98;
                                        const register double tmp377_1 = C_2_5*w100;
                                        const register double tmp94_1 = tmp7_0*w91;
                                        const register double tmp345_1 = tmp37_0*w95;
                                        const register double tmp411_1 = C_0_0*w110;
                                        const register double tmp438_1 = C_1_4*w119;
                                        const register double tmp112_1 = tmp9_0*w105;
                                        const register double tmp54_1 = tmp10_0*w108;
                                        const register double tmp324_1 = C_0_6*w104;
                                        const register double tmp93_1 = tmp0_0*w104;
                                        const register double tmp37_1 = C_1_6*w118;
                                        const register double tmp329_1 = C_1_0*w114;
                                        const register double tmp357_1 = tmp0_0*w110;
                                        const register double tmp180_1 = C_1_3*w119;
                                        const register double tmp273_1 = C_1_0*w91;
                                        const register double tmp65_1 = tmp28_0*w92;
                                        const register double tmp299_1 = tmp4_0*w90;
                                        const register double tmp184_1 = tmp23_0*w97;
                                        const register double tmp245_1 = C_0_5*w124;
                                        const register double tmp341_1 = C_1_6*w119;
                                        const register double tmp224_1 = tmp8_0*w93;
                                        const register double tmp202_1 = C_0_6*w110;
                                        const register double tmp302_1 = C_2_6*w122;
                                        const register double tmp18_1 = tmp14_0*w98;
                                        const register double tmp193_1 = C_1_5*w115;
                                        const register double tmp25_1 = tmp10_0*w96;
                                        const register double tmp398_1 = C_1_0*w115;
                                        const register double tmp198_1 = C_1_7*w105;
                                        const register double tmp380_1 = C_1_7*w111;
                                        const register double tmp450_1 = tmp24_0*w95;
                                        const register double tmp138_1 = tmp31_0*w112;
                                        const register double tmp403_1 = C_0_3*w98;
                                        const register double tmp212_1 = tmp23_0*w109;
                                        const register double tmp246_1 = tmp37_0*w112;
                                        const register double tmp167_1 = tmp31_0*w102;
                                        const register double tmp57_1 = tmp18_0*w109;
                                        const register double tmp5_1 = C_2_7*w123;
                                        const register double tmp207_1 = tmp26_0*w108;
                                        const register double tmp52_1 = tmp14_0*w104;
                                        const register double tmp175_1 = tmp6_0*w117;
                                        const register double tmp425_1 = tmp18_0*w105;
                                        const register double tmp343_1 = C_1_4*w99;
                                        const register double tmp383_1 = C_0_2*w120;
                                        const register double tmp163_1 = tmp12_0*w94;
                                        const register double tmp455_1 = C_2_5*w116;
                                        const register double tmp428_1 = tmp21_0*w95;
                                        const register double tmp347_1 = C_0_0*w104;
                                        const register double tmp451_1 = tmp20_0*w102;
                                        const register double tmp298_1 = tmp28_0*w116;
                                        const register double tmp410_1 = C_0_7*w104;
                                        const register double tmp268_1 = tmp21_0*w103;
                                        const register double tmp305_1 = C_2_1*w123;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=tmp0_1 + tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp14_1 + tmp15_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp33_1 + tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                                        EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1 + tmp60_1 + tmp61_1;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=tmp62_1 + tmp63_1 + tmp64_1 + tmp65_1 + tmp66_1 + tmp67_1 + tmp68_1 + tmp69_1 + tmp70_1 + tmp71_1 + tmp72_1 + tmp73_1 + tmp74_1 + tmp75_1;
                                        EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp76_1 + tmp77_1 + tmp78_1 + tmp79_1 + tmp80_1 + tmp81_1;
                                        EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp16_1 + tmp17_1 + tmp19_1 + tmp20_1 + tmp22_1 + tmp23_1 + tmp82_1 + tmp83_1 + tmp84_1 + tmp85_1;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=tmp86_1 + tmp87_1 + tmp88_1 + tmp89_1 + tmp90_1 + tmp91_1 + tmp92_1 + tmp93_1 + tmp94_1 + tmp95_1 + tmp96_1 + tmp97_1 + tmp98_1 + tmp99_1;
                                        EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp100_1 + tmp101_1 + tmp78_1 + tmp79_1 + tmp80_1 + tmp81_1;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=tmp102_1 + tmp103_1 + tmp104_1 + tmp105_1 + tmp106_1 + tmp107_1 + tmp108_1 + tmp109_1 + tmp110_1 + tmp111_1 + tmp112_1 + tmp113_1 + tmp114_1 + tmp115_1;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=tmp0_1 + tmp10_1 + tmp116_1 + tmp117_1 + tmp118_1 + tmp119_1 + tmp120_1 + tmp121_1 + tmp12_1 + tmp13_1 + tmp1_1 + tmp3_1 + tmp7_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp122_1 + tmp123_1 + tmp76_1 + tmp77_1 + tmp78_1 + tmp81_1;
                                        EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp124_1 + tmp125_1 + tmp126_1 + tmp127_1 + tmp128_1 + tmp129_1 + tmp130_1 + tmp131_1 + tmp58_1 + tmp59_1;
                                        EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp132_1 + tmp133_1 + tmp134_1 + tmp135_1 + tmp136_1 + tmp137_1 + tmp138_1 + tmp139_1 + tmp140_1 + tmp141_1;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=tmp142_1 + tmp143_1 + tmp144_1 + tmp145_1 + tmp146_1 + tmp147_1 + tmp148_1 + tmp149_1 + tmp150_1 + tmp151_1 + tmp152_1 + tmp153_1 + tmp6_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp126_1 + tmp128_1 + tmp129_1 + tmp131_1 + tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1 + tmp58_1 + tmp59_1;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=tmp154_1 + tmp155_1 + tmp156_1 + tmp157_1 + tmp158_1 + tmp159_1 + tmp160_1 + tmp161_1 + tmp162_1 + tmp163_1 + tmp164_1 + tmp165_1 + tmp166_1 + tmp167_1;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=tmp168_1 + tmp169_1 + tmp170_1 + tmp171_1 + tmp172_1 + tmp173_1 + tmp174_1 + tmp175_1 + tmp176_1 + tmp177_1 + tmp178_1 + tmp179_1 + tmp180_1 + tmp181_1 + tmp182_1 + tmp183_1 + tmp184_1 + tmp185_1;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=tmp186_1 + tmp187_1 + tmp188_1 + tmp189_1 + tmp190_1 + tmp191_1 + tmp192_1 + tmp193_1 + tmp194_1 + tmp195_1 + tmp196_1 + tmp197_1 + tmp198_1 + tmp199_1;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=tmp157_1 + tmp158_1 + tmp159_1 + tmp160_1 + tmp162_1 + tmp163_1 + tmp166_1 + tmp167_1 + tmp200_1 + tmp201_1 + tmp202_1 + tmp203_1 + tmp204_1 + tmp205_1;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=tmp206_1 + tmp207_1 + tmp208_1 + tmp209_1 + tmp210_1 + tmp211_1 + tmp212_1 + tmp213_1 + tmp214_1 + tmp215_1 + tmp216_1 + tmp217_1 + tmp218_1 + tmp219_1 + tmp220_1 + tmp221_1 + tmp222_1 + tmp223_1;
                                        EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp132_1 + tmp136_1 + tmp138_1 + tmp139_1 + tmp140_1 + tmp141_1 + tmp224_1 + tmp225_1 + tmp226_1 + tmp227_1;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=tmp102_1 + tmp103_1 + tmp105_1 + tmp106_1 + tmp110_1 + tmp112_1 + tmp113_1 + tmp115_1 + tmp228_1 + tmp229_1 + tmp230_1 + tmp231_1 + tmp92_1 + tmp98_1;
                                        EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp133_1 + tmp134_1 + tmp135_1 + tmp137_1 + tmp139_1 + tmp141_1 + tmp232_1 + tmp233_1 + tmp234_1 + tmp235_1;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=tmp206_1 + tmp207_1 + tmp236_1 + tmp237_1 + tmp238_1 + tmp239_1 + tmp240_1 + tmp241_1 + tmp242_1 + tmp243_1 + tmp244_1 + tmp245_1 + tmp246_1 + tmp247_1 + tmp248_1 + tmp249_1 + tmp250_1 + tmp251_1;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=tmp215_1 + tmp222_1 + tmp247_1 + tmp250_1 + tmp252_1 + tmp253_1 + tmp254_1 + tmp255_1 + tmp256_1 + tmp257_1 + tmp258_1 + tmp259_1 + tmp260_1 + tmp261_1 + tmp262_1 + tmp263_1 + tmp264_1 + tmp265_1;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=tmp266_1 + tmp267_1 + tmp268_1 + tmp269_1 + tmp270_1 + tmp271_1 + tmp272_1 + tmp273_1 + tmp274_1 + tmp275_1 + tmp276_1 + tmp277_1 + tmp278_1 + tmp279_1;
                                        EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp280_1 + tmp281_1 + tmp282_1 + tmp283_1 + tmp284_1 + tmp285_1 + tmp286_1 + tmp287_1 + tmp288_1 + tmp289_1;
                                        EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp290_1 + tmp291_1 + tmp292_1 + tmp293_1 + tmp294_1 + tmp295_1 + tmp296_1 + tmp297_1 + tmp298_1 + tmp299_1;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=tmp211_1 + tmp212_1 + tmp244_1 + tmp246_1 + tmp252_1 + tmp254_1 + tmp300_1 + tmp301_1 + tmp302_1 + tmp303_1 + tmp304_1 + tmp305_1 + tmp306_1 + tmp307_1 + tmp308_1 + tmp309_1 + tmp310_1 + tmp311_1;
                                        EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp139_1 + tmp141_1 + tmp224_1 + tmp225_1 + tmp226_1 + tmp227_1 + tmp232_1 + tmp233_1 + tmp234_1 + tmp235_1;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=tmp312_1 + tmp313_1 + tmp314_1 + tmp315_1 + tmp316_1 + tmp317_1 + tmp318_1 + tmp319_1 + tmp320_1 + tmp321_1 + tmp322_1 + tmp323_1 + tmp62_1 + tmp63_1;
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=tmp175_1 + tmp178_1 + tmp324_1 + tmp325_1 + tmp326_1 + tmp327_1 + tmp328_1 + tmp329_1 + tmp330_1 + tmp331_1 + tmp332_1 + tmp333_1 + tmp334_1 + tmp335_1 + tmp336_1 + tmp337_1 + tmp338_1 + tmp339_1;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=tmp179_1 + tmp184_1 + tmp325_1 + tmp326_1 + tmp340_1 + tmp341_1 + tmp342_1 + tmp343_1 + tmp344_1 + tmp345_1 + tmp346_1 + tmp347_1 + tmp348_1 + tmp349_1 + tmp350_1 + tmp351_1 + tmp352_1 + tmp353_1;
                                        EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp292_1 + tmp293_1 + tmp294_1 + tmp295_1 + tmp297_1 + tmp298_1 + tmp354_1 + tmp355_1 + tmp356_1 + tmp357_1;
                                        EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp280_1 + tmp281_1 + tmp283_1 + tmp284_1 + tmp286_1 + tmp288_1 + tmp358_1 + tmp359_1 + tmp360_1 + tmp361_1;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=tmp362_1 + tmp363_1 + tmp364_1 + tmp365_1 + tmp366_1 + tmp367_1 + tmp65_1 + tmp66_1 + tmp67_1 + tmp68_1 + tmp70_1 + tmp72_1 + tmp74_1 + tmp75_1;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=tmp119_1 + tmp121_1 + tmp142_1 + tmp143_1 + tmp144_1 + tmp146_1 + tmp147_1 + tmp150_1 + tmp151_1 + tmp153_1 + tmp368_1 + tmp369_1 + tmp370_1 + tmp371_1;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=tmp168_1 + tmp170_1 + tmp330_1 + tmp335_1 + tmp345_1 + tmp352_1 + tmp372_1 + tmp373_1 + tmp374_1 + tmp375_1 + tmp376_1 + tmp377_1 + tmp378_1 + tmp379_1 + tmp380_1 + tmp381_1 + tmp382_1 + tmp383_1;
                                        EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp292_1 + tmp293_1 + tmp354_1 + tmp355_1 + tmp356_1 + tmp357_1 + tmp384_1 + tmp385_1 + tmp386_1 + tmp387_1;
                                        EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp388_1 + tmp389_1 + tmp390_1 + tmp391_1 + tmp392_1 + tmp393_1 + tmp394_1 + tmp395_1 + tmp396_1 + tmp397_1;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=tmp190_1 + tmp191_1 + tmp266_1 + tmp267_1 + tmp268_1 + tmp269_1 + tmp271_1 + tmp272_1 + tmp276_1 + tmp279_1 + tmp398_1 + tmp399_1 + tmp400_1 + tmp401_1;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=tmp312_1 + tmp314_1 + tmp315_1 + tmp317_1 + tmp318_1 + tmp319_1 + tmp321_1 + tmp322_1 + tmp362_1 + tmp364_1 + tmp402_1 + tmp403_1 + tmp404_1 + tmp405_1;
                                        EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp124_1 + tmp125_1 + tmp127_1 + tmp130_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1 + tmp60_1 + tmp61_1;
                                        EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp388_1 + tmp389_1 + tmp392_1 + tmp393_1 + tmp394_1 + tmp395_1 + tmp406_1 + tmp407_1 + tmp408_1 + tmp409_1;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=tmp200_1 + tmp201_1 + tmp410_1 + tmp411_1 + tmp412_1 + tmp413_1 + tmp414_1 + tmp415_1 + tmp416_1 + tmp417_1 + tmp418_1 + tmp419_1 + tmp420_1 + tmp421_1;
                                        EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp14_1 + tmp15_1 + tmp17_1 + tmp18_1 + tmp21_1 + tmp23_1 + tmp422_1 + tmp423_1 + tmp424_1 + tmp425_1;
                                        EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp280_1 + tmp281_1 + tmp358_1 + tmp359_1 + tmp360_1 + tmp361_1 + tmp426_1 + tmp427_1 + tmp428_1 + tmp429_1;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=tmp154_1 + tmp156_1 + tmp412_1 + tmp413_1 + tmp414_1 + tmp415_1 + tmp416_1 + tmp418_1 + tmp419_1 + tmp421_1 + tmp430_1 + tmp431_1 + tmp432_1 + tmp433_1;
                                        EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp100_1 + tmp101_1 + tmp122_1 + tmp123_1 + tmp434_1 + tmp435_1;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=tmp32_1 + tmp36_1 + tmp38_1 + tmp39_1 + tmp40_1 + tmp42_1 + tmp436_1 + tmp437_1 + tmp438_1 + tmp439_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp51_1;
                                        EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp100_1 + tmp101_1 + tmp434_1 + tmp435_1 + tmp79_1 + tmp80_1;
                                        EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp122_1 + tmp123_1 + tmp434_1 + tmp435_1 + tmp76_1 + tmp77_1;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=tmp24_1 + tmp25_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp31_1 + tmp34_1 + tmp35_1 + tmp43_1 + tmp440_1 + tmp441_1 + tmp442_1 + tmp443_1 + tmp44_1;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=tmp186_1 + tmp187_1 + tmp188_1 + tmp189_1 + tmp192_1 + tmp194_1 + tmp195_1 + tmp197_1 + tmp274_1 + tmp277_1 + tmp444_1 + tmp445_1 + tmp446_1 + tmp447_1;
                                        EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp17_1 + tmp23_1 + tmp422_1 + tmp423_1 + tmp424_1 + tmp425_1 + tmp82_1 + tmp83_1 + tmp84_1 + tmp85_1;
                                        EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp388_1 + tmp389_1 + tmp406_1 + tmp407_1 + tmp408_1 + tmp409_1 + tmp448_1 + tmp449_1 + tmp450_1 + tmp451_1;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=tmp109_1 + tmp111_1 + tmp452_1 + tmp453_1 + tmp454_1 + tmp455_1 + tmp86_1 + tmp87_1 + tmp88_1 + tmp90_1 + tmp93_1 + tmp94_1 + tmp96_1 + tmp97_1;
                                        EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp434_1 + tmp435_1 + tmp76_1 + tmp77_1 + tmp79_1 + tmp80_1;
                                        EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp280_1 + tmp281_1 + tmp282_1 + tmp285_1 + tmp287_1 + tmp289_1 + tmp426_1 + tmp427_1 + tmp428_1 + tmp429_1;
                                        EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp290_1 + tmp291_1 + tmp292_1 + tmp293_1 + tmp296_1 + tmp299_1 + tmp384_1 + tmp385_1 + tmp386_1 + tmp387_1;
                                        EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp100_1 + tmp101_1 + tmp122_1 + tmp123_1 + tmp78_1 + tmp81_1;
                                        EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp388_1 + tmp389_1 + tmp390_1 + tmp391_1 + tmp396_1 + tmp397_1 + tmp448_1 + tmp449_1 + tmp450_1 + tmp451_1;
                                    }
                                }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const register double C_0 = C_p[INDEX3(k, m, 0, numEq, numComp)];
                                        const register double C_1 = C_p[INDEX3(k, m, 1, numEq, numComp)];
                                        const register double C_2 = C_p[INDEX3(k, m, 2, numEq, numComp)];
                                        const register double tmp1_1 = C_1*w127;
                                        const register double tmp15_1 = C_0*w138;
                                        const register double tmp4_1 = C_2*w133;
                                        const register double tmp3_1 = C_2*w131;
                                        const register double tmp8_1 = C_1*w132;
                                        const register double tmp7_1 = C_0*w129;
                                        const register double tmp0_1 = C_2*w140;
                                        const register double tmp2_1 = C_0*w126;
                                        const register double tmp16_1 = C_1*w139;
                                        const register double tmp13_1 = C_0*w141;
                                        const register double tmp14_1 = C_2*w128;
                                        const register double tmp11_1 = C_0*w143;
                                        const register double tmp12_1 = C_1*w142;
                                        const register double tmp5_1 = C_1*w134;
                                        const register double tmp17_1 = C_0*w137;
                                        const register double tmp6_1 = C_2*w135;
                                        const register double tmp10_1 = C_2*w136;
                                        const register double tmp9_1 = C_1*w130;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                        EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=tmp2_1 + tmp4_1 + tmp5_1;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=tmp6_1 + tmp7_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp10_1 + tmp7_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=tmp11_1 + tmp6_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp12_1 + tmp13_1 + tmp3_1;
                                        EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp1_1 + tmp3_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=tmp14_1 + tmp1_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp12_1 + tmp15_1 + tmp3_1;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=tmp14_1 + tmp1_1 + tmp2_1;
                                        EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp13_1 + tmp16_1 + tmp3_1;
                                        EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp10_1 + tmp1_1 + tmp2_1;
                                        EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp12_1 + tmp4_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=tmp0_1 + tmp7_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp10_1 + tmp1_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=tmp17_1 + tmp1_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=tmp0_1 + tmp17_1 + tmp5_1;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=tmp4_1 + tmp7_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=tmp11_1 + tmp1_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=tmp11_1 + tmp14_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp12_1 + tmp2_1 + tmp4_1;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=tmp14_1 + tmp2_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp12_1 + tmp6_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=tmp0_1 + tmp11_1 + tmp5_1;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=tmp14_1 + tmp17_1 + tmp5_1;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=tmp2_1 + tmp5_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp13_1 + tmp1_1 + tmp4_1;
                                        EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp16_1 + tmp2_1 + tmp4_1;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=tmp0_1 + tmp17_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp12_1 + tmp2_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=tmp11_1 + tmp1_1 + tmp4_1;
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=tmp0_1 + tmp11_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=tmp11_1 + tmp14_1 + tmp5_1;
                                        EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp16_1 + tmp4_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp13_1 + tmp4_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=tmp17_1 + tmp6_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=tmp14_1 + tmp7_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=tmp14_1 + tmp17_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp16_1 + tmp6_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp15_1 + tmp4_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=tmp2_1 + tmp6_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=tmp17_1 + tmp1_1 + tmp4_1;
                                        EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp10_1 + tmp2_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp15_1 + tmp1_1 + tmp4_1;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=tmp11_1 + tmp4_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp2_1 + tmp3_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp13_1 + tmp6_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=tmp17_1 + tmp4_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp10_1 + tmp15_1 + tmp16_1;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=tmp5_1 + tmp6_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp10_1 + tmp12_1 + tmp15_1;
                                        EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp10_1 + tmp13_1 + tmp16_1;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=tmp2_1 + tmp4_1 + tmp8_1;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=tmp4_1 + tmp5_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp3_1 + tmp7_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp15_1 + tmp1_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp10_1 + tmp12_1 + tmp13_1;
                                        EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp13_1 + tmp1_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp16_1 + tmp2_1 + tmp6_1;
                                        EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp15_1 + tmp16_1 + tmp3_1;
                                        EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp15_1 + tmp6_1 + tmp9_1;
                                    }
                                }
                            }
                        }
                        ///////////////
                        // process D //
                        ///////////////
                        if (!D.isEmpty()) {
                            add_EM_S=true;
                            const double* D_p=const_cast<escript::Data*>(&D)->getSampleDataRO(e);
                            if (D.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const register double D_0 = D_p[INDEX3(k, m, 0, numEq, numComp)];
                                        const register double D_1 = D_p[INDEX3(k, m, 1, numEq, numComp)];
                                        const register double D_2 = D_p[INDEX3(k, m, 2, numEq, numComp)];
                                        const register double D_3 = D_p[INDEX3(k, m, 3, numEq, numComp)];
                                        const register double D_4 = D_p[INDEX3(k, m, 4, numEq, numComp)];
                                        const register double D_5 = D_p[INDEX3(k, m, 5, numEq, numComp)];
                                        const register double D_6 = D_p[INDEX3(k, m, 6, numEq, numComp)];
                                        const register double D_7 = D_p[INDEX3(k, m, 7, numEq, numComp)];
                                        const register double tmp30_0 = D_0 + D_2 + D_4 + D_6;
                                        const register double tmp22_0 = D_1 + D_3 + D_4 + D_6;
                                        const register double tmp6_0 = D_4 + D_6;
                                        const register double tmp1_0 = D_0 + D_4;
                                        const register double tmp27_0 = D_3 + D_5 + D_6;
                                        const register double tmp32_0 = D_2 + D_4 + D_7;
                                        const register double tmp8_0 = D_0 + D_1 + D_6 + D_7;
                                        const register double tmp24_0 = D_0 + D_2;
                                        const register double tmp10_0 = D_4 + D_5;
                                        const register double tmp16_0 = D_0 + D_1 + D_4 + D_5;
                                        const register double tmp21_0 = D_0 + D_5 + D_6;
                                        const register double tmp29_0 = D_1 + D_3 + D_5 + D_7;
                                        const register double tmp12_0 = D_0 + D_3 + D_4 + D_7;
                                        const register double tmp28_0 = D_1 + D_2 + D_4;
                                        const register double tmp11_0 = D_0 + D_1 + D_2 + D_3 + D_4 + D_5 + D_6 + D_7;
                                        const register double tmp23_0 = D_5 + D_7;
                                        const register double tmp0_0 = D_1 + D_2 + D_5 + D_6;
                                        const register double tmp26_0 = D_1 + D_4 + D_7;
                                        const register double tmp7_0 = D_1 + D_3;
                                        const register double tmp4_0 = D_0 + D_1 + D_2 + D_3;
                                        const register double tmp31_0 = D_0 + D_3 + D_5;
                                        const register double tmp18_0 = D_0 + D_1;
                                        const register double tmp3_0 = D_4 + D_5 + D_6 + D_7;
                                        const register double tmp25_0 = D_0 + D_3 + D_6;
                                        const register double tmp2_0 = D_3 + D_7;
                                        const register double tmp19_0 = D_6 + D_7;
                                        const register double tmp20_0 = D_1 + D_2 + D_7;
                                        const register double tmp15_0 = D_2 + D_3 + D_6 + D_7;
                                        const register double tmp5_0 = D_0 + D_2 + D_5 + D_7;
                                        const register double tmp9_0 = D_2 + D_3;
                                        const register double tmp17_0 = D_2 + D_3 + D_4 + D_5;
                                        const register double tmp13_0 = D_1 + D_5;
                                        const register double tmp14_0 = D_2 + D_6;
                                        const register double tmp39_1 = tmp25_0*w148;
                                        const register double tmp40_1 = D_5*w150;
                                        const register double tmp4_1 = tmp4_0*w147;
                                        const register double tmp42_1 = D_2*w149;
                                        const register double tmp71_1 = tmp30_0*w148;
                                        const register double tmp36_1 = D_3*w150;
                                        const register double tmp66_1 = D_6*w149;
                                        const register double tmp19_1 = tmp14_0*w144;
                                        const register double tmp29_1 = D_4*w150;
                                        const register double tmp73_1 = tmp19_0*w144;
                                        const register double tmp28_1 = tmp20_0*w148;
                                        const register double tmp24_1 = tmp1_0*w146;
                                        const register double tmp14_1 = tmp10_0*w146;
                                        const register double tmp21_1 = tmp15_0*w148;
                                        const register double tmp57_1 = tmp10_0*w144;
                                        const register double tmp10_1 = tmp4_0*w148;
                                        const register double tmp58_1 = tmp9_0*w146;
                                        const register double tmp69_1 = tmp25_0*w147;
                                        const register double tmp13_1 = tmp9_0*w144;
                                        const register double tmp26_1 = tmp18_0*w144;
                                        const register double tmp52_1 = tmp15_0*w147;
                                        const register double tmp72_1 = tmp29_0*w147;
                                        const register double tmp18_1 = tmp14_0*w146;
                                        const register double tmp25_1 = tmp17_0*w145;
                                        const register double tmp55_1 = tmp32_0*w147;
                                        const register double tmp65_1 = tmp31_0*w147;
                                        const register double tmp56_1 = D_1*w149;
                                        const register double tmp45_1 = tmp28_0*w147;
                                        const register double tmp68_1 = D_2*w150;
                                        const register double tmp23_1 = tmp2_0*w144;
                                        const register double tmp2_1 = tmp2_0*w146;
                                        const register double tmp27_1 = tmp19_0*w146;
                                        const register double tmp41_1 = tmp26_0*w147;
                                        const register double tmp33_1 = tmp23_0*w144;
                                        const register double tmp32_1 = tmp22_0*w145;
                                        const register double tmp17_1 = tmp13_0*w144;
                                        const register double tmp15_1 = tmp11_0*w145;
                                        const register double tmp20_1 = tmp13_0*w146;
                                        const register double tmp5_1 = tmp5_0*w145;
                                        const register double tmp74_1 = tmp18_0*w146;
                                        const register double tmp63_1 = tmp32_0*w148;
                                        const register double tmp43_1 = tmp27_0*w148;
                                        const register double tmp48_1 = tmp23_0*w146;
                                        const register double tmp8_1 = tmp7_0*w144;
                                        const register double tmp7_1 = tmp7_0*w146;
                                        const register double tmp53_1 = tmp31_0*w148;
                                        const register double tmp59_1 = tmp28_0*w148;
                                        const register double tmp22_1 = tmp16_0*w147;
                                        const register double tmp60_1 = D_7*w150;
                                        const register double tmp61_1 = tmp27_0*w147;
                                        const register double tmp3_1 = tmp3_0*w148;
                                        const register double tmp6_1 = tmp6_0*w144;
                                        const register double tmp9_1 = tmp6_0*w146;
                                        const register double tmp34_1 = tmp24_0*w146;
                                        const register double tmp47_1 = tmp24_0*w144;
                                        const register double tmp30_1 = tmp21_0*w147;
                                        const register double tmp50_1 = tmp30_0*w147;
                                        const register double tmp67_1 = tmp26_0*w148;
                                        const register double tmp11_1 = tmp3_0*w147;
                                        const register double tmp64_1 = D_1*w150;
                                        const register double tmp46_1 = D_7*w149;
                                        const register double tmp37_1 = tmp20_0*w147;
                                        const register double tmp38_1 = D_4*w149;
                                        const register double tmp51_1 = tmp16_0*w148;
                                        const register double tmp54_1 = D_6*w150;
                                        const register double tmp16_1 = tmp12_0*w145;
                                        const register double tmp31_1 = D_3*w149;
                                        const register double tmp12_1 = tmp8_0*w145;
                                        const register double tmp62_1 = D_0*w149;
                                        const register double tmp49_1 = tmp29_0*w148;
                                        const register double tmp1_1 = tmp1_0*w144;
                                        const register double tmp0_1 = tmp0_0*w145;
                                        const register double tmp70_1 = D_5*w149;
                                        const register double tmp35_1 = tmp21_0*w148;
                                        const register double tmp44_1 = D_0*w150;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                        EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp3_1 + tmp4_1;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=tmp5_1 + tmp6_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=tmp5_1 + tmp8_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp10_1 + tmp11_1;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=tmp12_1 + tmp13_1 + tmp14_1;
                                        EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp15_1;
                                        EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp3_1 + tmp4_1;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=tmp16_1 + tmp17_1 + tmp18_1;
                                        EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp15_1;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=tmp16_1 + tmp19_1 + tmp20_1;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                        EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp15_1;
                                        EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp10_1 + tmp11_1;
                                        EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp21_1 + tmp22_1;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=tmp0_1 + tmp23_1 + tmp24_1;
                                        EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp10_1 + tmp11_1;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=tmp25_1 + tmp26_1 + tmp27_1;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=tmp32_1 + tmp33_1 + tmp34_1;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=tmp25_1 + tmp26_1 + tmp27_1;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1;
                                        EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp21_1 + tmp22_1;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=tmp16_1 + tmp19_1 + tmp20_1;
                                        EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp21_1 + tmp22_1;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=tmp32_1 + tmp47_1 + tmp48_1;
                                        EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp49_1 + tmp50_1;
                                        EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp51_1 + tmp52_1;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=tmp53_1 + tmp54_1 + tmp55_1 + tmp56_1;
                                        EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp21_1 + tmp22_1;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=tmp12_1 + tmp57_1 + tmp58_1;
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=tmp59_1 + tmp60_1 + tmp61_1 + tmp62_1;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=tmp63_1 + tmp64_1 + tmp65_1 + tmp66_1;
                                        EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp51_1 + tmp52_1;
                                        EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp49_1 + tmp50_1;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=tmp12_1 + tmp13_1 + tmp14_1;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=tmp0_1 + tmp23_1 + tmp24_1;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=tmp67_1 + tmp68_1 + tmp69_1 + tmp70_1;
                                        EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp51_1 + tmp52_1;
                                        EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp71_1 + tmp72_1;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=tmp32_1 + tmp47_1 + tmp48_1;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=tmp12_1 + tmp57_1 + tmp58_1;
                                        EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp10_1 + tmp11_1;
                                        EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp71_1 + tmp72_1;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=tmp25_1 + tmp73_1 + tmp74_1;
                                        EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp3_1 + tmp4_1;
                                        EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp49_1 + tmp50_1;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=tmp25_1 + tmp73_1 + tmp74_1;
                                        EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp15_1;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=tmp5_1 + tmp8_1 + tmp9_1;
                                        EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp15_1;
                                        EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp15_1;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=tmp5_1 + tmp6_1 + tmp7_1;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=tmp32_1 + tmp33_1 + tmp34_1;
                                        EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp3_1 + tmp4_1;
                                        EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp71_1 + tmp72_1;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=tmp16_1 + tmp17_1 + tmp18_1;
                                        EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp15_1;
                                        EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp49_1 + tmp50_1;
                                        EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp51_1 + tmp52_1;
                                        EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp15_1;
                                        EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp71_1 + tmp72_1;
                                    }
                                 }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const register double D_0 = D_p[INDEX2(k, m, numEq)];
                                        const register double tmp3_1 = D_0*w154;
                                        const register double tmp0_1 = D_0*w151;
                                        const register double tmp2_1 = D_0*w153;
                                        const register double tmp1_1 = D_0*w152;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp2_1;
                                        EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp2_1;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp2_1;
                                        EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=tmp3_1;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=tmp3_1;
                                        EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=tmp3_1;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=tmp3_1;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=tmp3_1;
                                        EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=tmp3_1;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=tmp3_1;
                                        EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=tmp3_1;
                                        EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp2_1;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp2_1;
                                        EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp2_1;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=tmp0_1;
                                        EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp2_1;
                                        EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp1_1;
                                        EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp2_1;
                                        EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp1_1;
                                    }
                                }
                            }
                        }
                        ///////////////
                        // process X //
                        ///////////////
                        if (!X.isEmpty()) {
                            add_EM_F=true;
                            const double* X_p=const_cast<escript::Data*>(&X)->getSampleDataRO(e);
                            if (X.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    const register double X_0_0 = X_p[INDEX3(k, 0, 0, numEq, 3)];
                                    const register double X_1_0 = X_p[INDEX3(k, 1, 0, numEq, 3)];
                                    const register double X_2_0 = X_p[INDEX3(k, 2, 0, numEq, 3)];
                                    const register double X_0_1 = X_p[INDEX3(k, 0, 1, numEq, 3)];
                                    const register double X_1_1 = X_p[INDEX3(k, 1, 1, numEq, 3)];
                                    const register double X_2_1 = X_p[INDEX3(k, 2, 1, numEq, 3)];
                                    const register double X_0_2 = X_p[INDEX3(k, 0, 2, numEq, 3)];
                                    const register double X_1_2 = X_p[INDEX3(k, 1, 2, numEq, 3)];
                                    const register double X_2_2 = X_p[INDEX3(k, 2, 2, numEq, 3)];
                                    const register double X_0_3 = X_p[INDEX3(k, 0, 3, numEq, 3)];
                                    const register double X_1_3 = X_p[INDEX3(k, 1, 3, numEq, 3)];
                                    const register double X_2_3 = X_p[INDEX3(k, 2, 3, numEq, 3)];
                                    const register double X_0_4 = X_p[INDEX3(k, 0, 4, numEq, 3)];
                                    const register double X_1_4 = X_p[INDEX3(k, 1, 4, numEq, 3)];
                                    const register double X_2_4 = X_p[INDEX3(k, 2, 4, numEq, 3)];
                                    const register double X_0_5 = X_p[INDEX3(k, 0, 5, numEq, 3)];
                                    const register double X_1_5 = X_p[INDEX3(k, 1, 5, numEq, 3)];
                                    const register double X_2_5 = X_p[INDEX3(k, 2, 5, numEq, 3)];
                                    const register double X_0_6 = X_p[INDEX3(k, 0, 6, numEq, 3)];
                                    const register double X_1_6 = X_p[INDEX3(k, 1, 6, numEq, 3)];
                                    const register double X_2_6 = X_p[INDEX3(k, 2, 6, numEq, 3)];
                                    const register double X_0_7 = X_p[INDEX3(k, 0, 7, numEq, 3)];
                                    const register double X_1_7 = X_p[INDEX3(k, 1, 7, numEq, 3)];
                                    const register double X_2_7 = X_p[INDEX3(k, 2, 7, numEq, 3)];
                                    const register double tmp3_0 = X_1_5 + X_1_7;
                                    const register double tmp17_0 = X_0_4 + X_0_5;
                                    const register double tmp2_0 = X_0_2 + X_0_3 + X_0_4 + X_0_5;
                                    const register double tmp14_0 = X_1_1 + X_1_3;
                                    const register double tmp5_0 = X_2_1 + X_2_2 + X_2_5 + X_2_6;
                                    const register double tmp1_0 = X_0_0 + X_0_1;
                                    const register double tmp8_0 = X_1_0 + X_1_2;
                                    const register double tmp12_0 = X_2_1 + X_2_5;
                                    const register double tmp16_0 = X_0_2 + X_0_3;
                                    const register double tmp9_0 = X_1_4 + X_1_6;
                                    const register double tmp7_0 = X_1_1 + X_1_3 + X_1_4 + X_1_6;
                                    const register double tmp13_0 = X_1_0 + X_1_2 + X_1_5 + X_1_7;
                                    const register double tmp6_0 = X_2_0 + X_2_4;
                                    const register double tmp4_0 = X_2_3 + X_2_7;
                                    const register double tmp15_0 = X_0_0 + X_0_1 + X_0_6 + X_0_7;
                                    const register double tmp11_0 = X_2_0 + X_2_3 + X_2_4 + X_2_7;
                                    const register double tmp0_0 = X_0_6 + X_0_7;
                                    const register double tmp10_0 = X_2_2 + X_2_6;
                                    const register double tmp16_1 = tmp13_0*w158;
                                    const register double tmp29_1 = tmp15_0*w165;
                                    const register double tmp18_1 = tmp15_0*w160;
                                    const register double tmp4_1 = tmp4_0*w161;
                                    const register double tmp15_1 = tmp0_0*w166;
                                    const register double tmp56_1 = tmp8_0*w169;
                                    const register double tmp38_1 = tmp14_0*w162;
                                    const register double tmp51_1 = tmp10_0*w170;
                                    const register double tmp53_1 = tmp9_0*w167;
                                    const register double tmp37_1 = tmp5_0*w171;
                                    const register double tmp59_1 = tmp3_0*w167;
                                    const register double tmp6_1 = tmp6_0*w157;
                                    const register double tmp13_1 = tmp11_0*w159;
                                    const register double tmp14_1 = tmp12_0*w157;
                                    const register double tmp31_1 = tmp9_0*w169;
                                    const register double tmp28_1 = tmp13_0*w168;
                                    const register double tmp10_1 = tmp9_0*w162;
                                    const register double tmp3_1 = tmp3_0*w162;
                                    const register double tmp55_1 = tmp0_0*w164;
                                    const register double tmp54_1 = tmp6_0*w172;
                                    const register double tmp50_1 = tmp14_0*w169;
                                    const register double tmp20_1 = tmp12_0*w161;
                                    const register double tmp58_1 = tmp1_0*w166;
                                    const register double tmp24_1 = tmp17_0*w163;
                                    const register double tmp49_1 = tmp12_0*w172;
                                    const register double tmp39_1 = tmp6_0*w170;
                                    const register double tmp8_1 = tmp8_0*w156;
                                    const register double tmp46_1 = tmp16_0*w166;
                                    const register double tmp2_1 = tmp2_0*w160;
                                    const register double tmp25_1 = tmp8_0*w167;
                                    const register double tmp33_1 = tmp14_0*w167;
                                    const register double tmp12_1 = tmp2_0*w165;
                                    const register double tmp21_1 = tmp7_0*w168;
                                    const register double tmp19_1 = tmp16_0*w155;
                                    const register double tmp26_1 = tmp16_0*w164;
                                    const register double tmp35_1 = tmp17_0*w155;
                                    const register double tmp36_1 = tmp4_0*w172;
                                    const register double tmp9_1 = tmp1_0*w164;
                                    const register double tmp34_1 = tmp9_0*w156;
                                    const register double tmp40_1 = tmp16_0*w163;
                                    const register double tmp43_1 = tmp11_0*w171;
                                    const register double tmp57_1 = tmp4_0*w170;
                                    const register double tmp52_1 = tmp1_0*w163;
                                    const register double tmp41_1 = tmp10_0*w172;
                                    const register double tmp1_1 = tmp1_0*w155;
                                    const register double tmp11_1 = tmp10_0*w161;
                                    const register double tmp22_1 = tmp10_0*w157;
                                    const register double tmp23_1 = tmp3_0*w169;
                                    const register double tmp5_1 = tmp5_0*w159;
                                    const register double tmp27_1 = tmp6_0*w161;
                                    const register double tmp32_1 = tmp17_0*w166;
                                    const register double tmp42_1 = tmp17_0*w164;
                                    const register double tmp30_1 = tmp4_0*w157;
                                    const register double tmp47_1 = tmp3_0*w156;
                                    const register double tmp48_1 = tmp0_0*w155;
                                    const register double tmp7_1 = tmp7_0*w158;
                                    const register double tmp45_1 = tmp12_0*w170;
                                    const register double tmp44_1 = tmp8_0*w162;
                                    const register double tmp0_1 = tmp0_0*w163;
                                    const register double tmp17_1 = tmp14_0*w156;
                                    EM_F[INDEX2(k,0,numEq)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1;
                                    EM_F[INDEX2(k,1,numEq)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp16_1 + tmp17_1 + tmp9_1;
                                    EM_F[INDEX2(k,2,numEq)]+=tmp13_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp24_1 + tmp25_1;
                                    EM_F[INDEX2(k,3,numEq)]+=tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp33_1 + tmp5_1;
                                    EM_F[INDEX2(k,4,numEq)]+=tmp16_1 + tmp18_1 + tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1 + tmp40_1;
                                    EM_F[INDEX2(k,5,numEq)]+=tmp29_1 + tmp41_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp7_1;
                                    EM_F[INDEX2(k,6,numEq)]+=tmp28_1 + tmp2_1 + tmp43_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1 + tmp52_1 + tmp53_1;
                                    EM_F[INDEX2(k,7,numEq)]+=tmp12_1 + tmp21_1 + tmp37_1 + tmp54_1 + tmp55_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1;
                                }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    const register double X_0 = X_p[INDEX2(k, 0, numEq)];
                                    const register double X_1 = X_p[INDEX2(k, 1, numEq)];
                                    const register double X_2 = X_p[INDEX2(k, 2, numEq)];
                                    const register double tmp4_1 = X_1*w177;
                                    const register double tmp3_1 = X_0*w176;
                                    const register double tmp5_1 = X_2*w178;
                                    const register double tmp0_1 = X_0*w173;
                                    const register double tmp2_1 = X_2*w175;
                                    const register double tmp1_1 = X_1*w174;
                                    EM_F[INDEX2(k,0,numEq)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                    EM_F[INDEX2(k,1,numEq)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                    EM_F[INDEX2(k,2,numEq)]+=tmp0_1 + tmp2_1 + tmp4_1;
                                    EM_F[INDEX2(k,3,numEq)]+=tmp2_1 + tmp3_1 + tmp4_1;
                                    EM_F[INDEX2(k,4,numEq)]+=tmp0_1 + tmp1_1 + tmp5_1;
                                    EM_F[INDEX2(k,5,numEq)]+=tmp1_1 + tmp3_1 + tmp5_1;
                                    EM_F[INDEX2(k,6,numEq)]+=tmp0_1 + tmp4_1 + tmp5_1;
                                    EM_F[INDEX2(k,7,numEq)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                }
                            }
                        }
                        ///////////////
                        // process Y //
                        ///////////////
                        if (!Y.isEmpty()) {
                            add_EM_F=true;
                            const double* Y_p=const_cast<escript::Data*>(&Y)->getSampleDataRO(e);
                            if (Y.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    const register double Y_0 = Y_p[INDEX2(k, 0, numEq)];
                                    const register double Y_1 = Y_p[INDEX2(k, 1, numEq)];
                                    const register double Y_2 = Y_p[INDEX2(k, 2, numEq)];
                                    const register double Y_3 = Y_p[INDEX2(k, 3, numEq)];
                                    const register double Y_4 = Y_p[INDEX2(k, 4, numEq)];
                                    const register double Y_5 = Y_p[INDEX2(k, 5, numEq)];
                                    const register double Y_6 = Y_p[INDEX2(k, 6, numEq)];
                                    const register double Y_7 = Y_p[INDEX2(k, 7, numEq)];
                                    const register double tmp3_0 = Y_0 + Y_3 + Y_5;
                                    const register double tmp5_0 = Y_0 + Y_3 + Y_6;
                                    const register double tmp6_0 = Y_0 + Y_5 + Y_6;
                                    const register double tmp2_0 = Y_2 + Y_4 + Y_7;
                                    const register double tmp4_0 = Y_1 + Y_4 + Y_7;
                                    const register double tmp0_0 = Y_3 + Y_5 + Y_6;
                                    const register double tmp1_0 = Y_1 + Y_2 + Y_4;
                                    const register double tmp7_0 = Y_1 + Y_2 + Y_7;
                                    const register double tmp30_1 = tmp0_0*w180;
                                    const register double tmp29_1 = tmp1_0*w181;
                                    const register double tmp22_1 = tmp4_0*w180;
                                    const register double tmp7_1 = Y_6*w182;
                                    const register double tmp31_1 = Y_0*w182;
                                    const register double tmp28_1 = Y_7*w179;
                                    const register double tmp19_1 = Y_3*w182;
                                    const register double tmp12_1 = Y_3*w179;
                                    const register double tmp3_1 = Y_7*w182;
                                    const register double tmp26_1 = tmp2_0*w180;
                                    const register double tmp13_1 = tmp6_0*w181;
                                    const register double tmp2_1 = tmp1_0*w180;
                                    const register double tmp0_1 = Y_0*w179;
                                    const register double tmp17_1 = tmp7_0*w181;
                                    const register double tmp24_1 = Y_6*w179;
                                    const register double tmp10_1 = tmp5_0*w180;
                                    const register double tmp16_1 = Y_4*w179;
                                    const register double tmp8_1 = Y_2*w179;
                                    const register double tmp18_1 = tmp6_0*w180;
                                    const register double tmp5_1 = tmp2_0*w181;
                                    const register double tmp6_1 = tmp3_0*w180;
                                    const register double tmp27_1 = Y_1*w182;
                                    const register double tmp20_1 = Y_5*w179;
                                    const register double tmp9_1 = tmp4_0*w181;
                                    const register double tmp25_1 = tmp3_0*w181;
                                    const register double tmp11_1 = Y_5*w182;
                                    const register double tmp4_1 = Y_1*w179;
                                    const register double tmp1_1 = tmp0_0*w181;
                                    const register double tmp21_1 = tmp5_0*w181;
                                    const register double tmp14_1 = tmp7_0*w180;
                                    const register double tmp23_1 = Y_2*w182;
                                    const register double tmp15_1 = Y_4*w182;
                                    EM_F[INDEX2(k,0,numEq)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                                    EM_F[INDEX2(k,1,numEq)]+=tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1;
                                    EM_F[INDEX2(k,2,numEq)]+=tmp10_1 + tmp11_1 + tmp8_1 + tmp9_1;
                                    EM_F[INDEX2(k,3,numEq)]+=tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1;
                                    EM_F[INDEX2(k,4,numEq)]+=tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1;
                                    EM_F[INDEX2(k,5,numEq)]+=tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1;
                                    EM_F[INDEX2(k,6,numEq)]+=tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1;
                                    EM_F[INDEX2(k,7,numEq)]+=tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1;
                                }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    const register double Y_0 = Y_p[k];
                                    const register double tmp0_1 = Y_0*w183;
                                    EM_F[INDEX2(k,0,numEq)]+=tmp0_1;
                                    EM_F[INDEX2(k,1,numEq)]+=tmp0_1;
                                    EM_F[INDEX2(k,2,numEq)]+=tmp0_1;
                                    EM_F[INDEX2(k,3,numEq)]+=tmp0_1;
                                    EM_F[INDEX2(k,4,numEq)]+=tmp0_1;
                                    EM_F[INDEX2(k,5,numEq)]+=tmp0_1;
                                    EM_F[INDEX2(k,6,numEq)]+=tmp0_1;
                                    EM_F[INDEX2(k,7,numEq)]+=tmp0_1;
                                }
                            }
                        }

                        // add to matrix (if add_EM_S) and RHS (if add_EM_F)
                        const index_t firstNode=m_N0*m_N1*k2+m_N0*k1+k0;
                        IndexVector rowIndex;
                        rowIndex.push_back(m_dofMap[firstNode]);
                        rowIndex.push_back(m_dofMap[firstNode+1]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0+1]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0*m_N1]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0*m_N1+1]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0*(m_N1+1)]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0*(m_N1+1)+1]);
                        if (add_EM_F) {
                            double *F_p=rhs.getSampleDataRW(0);
                            for (index_t i=0; i<rowIndex.size(); i++) {
                                if (rowIndex[i]<getNumDOF()) {
                                    for (index_t eq=0; eq<numEq; eq++) {
                                        F_p[INDEX2(eq,rowIndex[i],numEq)]+=EM_F[INDEX2(eq,i,numEq)];
                                    }
                                }
                            }
                        }
                        if (add_EM_S) {
                            addToSystemMatrix(mat, rowIndex, numEq, rowIndex, numComp, EM_S);
                        }
                    } // end k0 loop
                } // end k1 loop
            } // end k2 loop
        } // end of colouring
    } // end of parallel region
}

//protected
void Brick::assemblePDESystemReduced(Paso_SystemMatrix* mat,
        escript::Data& rhs, const escript::Data& A, const escript::Data& B,
        const escript::Data& C, const escript::Data& D,
        const escript::Data& X, const escript::Data& Y,
        const escript::Data& d, const escript::Data& y) const
{
    const double h0 = m_l0/m_gNE0;
    const double h1 = m_l1/m_gNE1;
    const double h2 = m_l2/m_gNE2;
    dim_t numEq, numComp;
    if (!mat)
        numEq=numComp=(rhs.isEmpty() ? 1 : rhs.getDataPointSize());
    else {
        numEq=mat->logical_row_block_size;
        numComp=mat->logical_col_block_size;
    }

    const double w0 = .0625*h1*h2/h0;
    const double w1 = .0625*h2;
    const double w2 = -.0625*h1;
    const double w3 = .0625*h0*h2/h1;
    const double w4 = -.0625*h0;
    const double w5 = .0625*h1;
    const double w6 = .0625*h0;
    const double w7 = -.0625*h0*h1/h2;
    const double w8 = -.0625*h1*h2/h0;
    const double w9 = -.0625*h2;
    const double w10 = -.0625*h0*h2/h1;
    const double w11 = .0625*h0*h1/h2;
    const double w12 = .03125*h1*h2;
    const double w13 = .03125*h0*h2;
    const double w14 = .03125*h0*h1;
    const double w15 = -.03125*h1*h2;
    const double w16 = -.03125*h0*h2;
    const double w17 = -.03125*h0*h1;
    const double w18 = .015625*h0*h1*h2;
    const double w19 = -.25*h1*h2;
    const double w20 = -.25*h0*h2;
    const double w21 = -.25*h0*h1;
    const double w22 = .25*h1*h2;
    const double w23 = .25*h0*h2;
    const double w24 = .25*h0*h1;
    const double w25 = .125*h0*h1*h2;

    rhs.requireWrite();
#pragma omp parallel
    {
        for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
#pragma omp for
            for (index_t k2=k2_0; k2<m_NE2; k2+=2) {
                for (index_t k1=0; k1<m_NE1; ++k1) {
                    for (index_t k0=0; k0<m_NE0; ++k0)  {
                        bool add_EM_S=false;
                        bool add_EM_F=false;
                        vector<double> EM_S(8*8*numEq*numComp, 0);
                        vector<double> EM_F(8*numEq, 0);
                        const index_t e = k0 + m_NE0*k1 + m_NE0*m_NE1*k2;
                        ///////////////
                        // process A //
                        ///////////////
                        if (!A.isEmpty()) {
                            add_EM_S=true;
                            const double* A_p=const_cast<escript::Data*>(&A)->getSampleDataRO(e);
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const register double A_00 = A_p[INDEX4(k,0,m,0, numEq,3, numComp)];
                                    const register double A_01 = A_p[INDEX4(k,0,m,1, numEq,3, numComp)];
                                    const register double A_02 = A_p[INDEX4(k,0,m,2, numEq,3, numComp)];
                                    const register double A_10 = A_p[INDEX4(k,1,m,0, numEq,3, numComp)];
                                    const register double A_11 = A_p[INDEX4(k,1,m,1, numEq,3, numComp)];
                                    const register double A_12 = A_p[INDEX4(k,1,m,2, numEq,3, numComp)];
                                    const register double A_20 = A_p[INDEX4(k,2,m,0, numEq,3, numComp)];
                                    const register double A_21 = A_p[INDEX4(k,2,m,1, numEq,3, numComp)];
                                    const register double A_22 = A_p[INDEX4(k,2,m,2, numEq,3, numComp)];
                                    const register double tmp0_0 = A_01 + A_10;
                                    const register double tmp1_0 = A_02 + A_20;
                                    const register double tmp2_0 = A_12 + A_21;
                                    const register double tmp3_1 = A_22*w7;
                                    const register double tmp10_1 = A_11*w10;
                                    const register double tmp21_1 = A_02*w5;
                                    const register double tmp2_1 = A_00*w0;
                                    const register double tmp23_1 = tmp2_0*w6;
                                    const register double tmp19_1 = A_20*w2;
                                    const register double tmp4_1 = A_11*w3;
                                    const register double tmp22_1 = tmp1_0*w5;
                                    const register double tmp13_1 = A_21*w4;
                                    const register double tmp5_1 = A_21*w6;
                                    const register double tmp8_1 = A_00*w8;
                                    const register double tmp7_1 = A_20*w5;
                                    const register double tmp18_1 = tmp2_0*w4;
                                    const register double tmp6_1 = A_02*w2;
                                    const register double tmp9_1 = A_22*w11;
                                    const register double tmp15_1 = tmp1_0*w2;
                                    const register double tmp12_1 = A_01*w1;
                                    const register double tmp0_1 = tmp0_0*w1;
                                    const register double tmp20_1 = A_01*w9;
                                    const register double tmp14_1 = A_12*w6;
                                    const register double tmp1_1 = A_12*w4;
                                    const register double tmp16_1 = A_10*w9;
                                    const register double tmp11_1 = tmp0_0*w9;
                                    const register double tmp17_1 = A_10*w1;
                                    EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp10_1 + tmp11_1 + tmp1_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=tmp10_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp16_1 + tmp2_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=tmp10_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp16_1 + tmp2_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp10_1 + tmp11_1 + tmp1_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp4_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp10_1 + tmp11_1 + tmp15_1 + tmp18_1 + tmp3_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp0_1 + tmp10_1 + tmp19_1 + tmp1_1 + tmp21_1 + tmp5_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=tmp11_1 + tmp13_1 + tmp14_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp6_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp0_1 + tmp10_1 + tmp18_1 + tmp22_1 + tmp3_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=tmp11_1 + tmp13_1 + tmp14_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp6_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=tmp0_1 + tmp13_1 + tmp14_1 + tmp19_1 + tmp21_1 + tmp2_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp0_1 + tmp10_1 + tmp15_1 + tmp23_1 + tmp3_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp10_1 + tmp11_1 + tmp13_1 + tmp14_1 + tmp19_1 + tmp21_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp12_1 + tmp15_1 + tmp16_1 + tmp1_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=tmp0_1 + tmp13_1 + tmp14_1 + tmp19_1 + tmp21_1 + tmp2_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp0_1 + tmp10_1 + tmp13_1 + tmp14_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=tmp17_1 + tmp20_1 + tmp23_1 + tmp4_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=tmp0_1 + tmp15_1 + tmp18_1 + tmp2_1 + tmp4_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp1_1 + tmp22_1 + tmp2_1 + tmp5_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=tmp12_1 + tmp16_1 + tmp19_1 + tmp21_1 + tmp23_1 + tmp4_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=tmp0_1 + tmp15_1 + tmp18_1 + tmp2_1 + tmp4_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp17_1 + tmp1_1 + tmp20_1 + tmp22_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=tmp11_1 + tmp19_1 + tmp1_1 + tmp21_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp12_1 + tmp13_1 + tmp14_1 + tmp16_1 + tmp22_1 + tmp3_1 + tmp4_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=tmp11_1 + tmp18_1 + tmp22_1 + tmp2_1 + tmp4_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=tmp0_1 + tmp22_1 + tmp23_1 + tmp2_1 + tmp4_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp1_1 + tmp22_1 + tmp2_1 + tmp5_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp23_1 + tmp2_1 + tmp3_1 + tmp6_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp12_1 + tmp13_1 + tmp14_1 + tmp16_1 + tmp22_1 + tmp3_1 + tmp4_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=tmp11_1 + tmp15_1 + tmp23_1 + tmp2_1 + tmp4_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp13_1 + tmp14_1 + tmp15_1 + tmp17_1 + tmp20_1 + tmp3_1 + tmp4_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=tmp12_1 + tmp16_1 + tmp18_1 + tmp4_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=tmp0_1 + tmp22_1 + tmp23_1 + tmp2_1 + tmp4_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=tmp11_1 + tmp15_1 + tmp23_1 + tmp2_1 + tmp4_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp13_1 + tmp14_1 + tmp15_1 + tmp17_1 + tmp20_1 + tmp3_1 + tmp4_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp10_1 + tmp17_1 + tmp18_1 + tmp20_1 + tmp2_1 + tmp3_1 + tmp6_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=tmp12_1 + tmp16_1 + tmp18_1 + tmp4_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=tmp11_1 + tmp18_1 + tmp22_1 + tmp2_1 + tmp4_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp17_1 + tmp1_1 + tmp20_1 + tmp22_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp18_1 + tmp19_1 + tmp21_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=tmp10_1 + tmp13_1 + tmp14_1 + tmp17_1 + tmp20_1 + tmp22_1 + tmp2_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp4_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp0_1 + tmp10_1 + tmp19_1 + tmp1_1 + tmp21_1 + tmp5_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp10_1 + tmp17_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp23_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=tmp17_1 + tmp20_1 + tmp23_1 + tmp4_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp0_1 + tmp10_1 + tmp13_1 + tmp14_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp10_1 + tmp17_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp23_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=tmp12_1 + tmp16_1 + tmp19_1 + tmp21_1 + tmp23_1 + tmp4_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp10_1 + tmp11_1 + tmp15_1 + tmp18_1 + tmp3_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=tmp10_1 + tmp15_1 + tmp17_1 + tmp1_1 + tmp20_1 + tmp2_1 + tmp5_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp0_1 + tmp10_1 + tmp15_1 + tmp23_1 + tmp3_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp0_1 + tmp10_1 + tmp18_1 + tmp22_1 + tmp3_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=tmp10_1 + tmp15_1 + tmp17_1 + tmp1_1 + tmp20_1 + tmp2_1 + tmp5_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=tmp10_1 + tmp13_1 + tmp14_1 + tmp17_1 + tmp20_1 + tmp22_1 + tmp2_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp10_1 + tmp11_1 + tmp13_1 + tmp14_1 + tmp19_1 + tmp21_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp10_1 + tmp17_1 + tmp18_1 + tmp20_1 + tmp2_1 + tmp3_1 + tmp6_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=tmp11_1 + tmp19_1 + tmp1_1 + tmp21_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp10_1 + tmp11_1 + tmp22_1 + tmp23_1 + tmp3_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp18_1 + tmp19_1 + tmp21_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp12_1 + tmp15_1 + tmp16_1 + tmp1_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp10_1 + tmp11_1 + tmp22_1 + tmp23_1 + tmp3_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp23_1 + tmp2_1 + tmp3_1 + tmp6_1 + tmp7_1;
                                }
                            }
                        }
                        ///////////////
                        // process B //
                        ///////////////
                        if (!B.isEmpty()) {
                            add_EM_S=true;
                            const double* B_p=const_cast<escript::Data*>(&B)->getSampleDataRO(e);
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const register double B_0 = B_p[INDEX3(k,0,m, numEq, 3)];
                                    const register double B_1 = B_p[INDEX3(k,1,m, numEq, 3)];
                                    const register double B_2 = B_p[INDEX3(k,2,m, numEq, 3)];
                                    const register double tmp4_1 = B_0*w15;
                                    const register double tmp3_1 = B_1*w16;
                                    const register double tmp2_1 = B_0*w12;
                                    const register double tmp5_1 = B_2*w17;
                                    const register double tmp1_1 = B_2*w14;
                                    const register double tmp0_1 = B_1*w13;
                                    EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                                }
                            }
                        }
                        ///////////////
                        // process C //
                        ///////////////
                        if (!C.isEmpty()) {
                            add_EM_S=true;
                            const double* C_p=const_cast<escript::Data*>(&C)->getSampleDataRO(e);
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const register double C_0 = C_p[INDEX3(k, m, 0, numEq, numComp)];
                                    const register double C_1 = C_p[INDEX3(k, m, 1, numEq, numComp)];
                                    const register double C_2 = C_p[INDEX3(k, m, 2, numEq, numComp)];
                                    const register double tmp5_1 = C_0*w15;
                                    const register double tmp2_1 = C_0*w12;
                                    const register double tmp4_1 = C_1*w16;
                                    const register double tmp1_1 = C_2*w17;
                                    const register double tmp3_1 = C_2*w14;
                                    const register double tmp0_1 = C_1*w13;
                                    EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                }
                            }
                        }
                        ///////////////
                        // process D //
                        ///////////////
                        if (!D.isEmpty()) {
                            add_EM_S=true;
                            const double* D_p=const_cast<escript::Data*>(&D)->getSampleDataRO(e);
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const register double D_0 = D_p[INDEX2(k, m, numEq)];
                                    const register double tmp0_1 = D_0*w18;
                                    EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp0_1;
                                }
                            }
                        }
                        ///////////////
                        // process X //
                        ///////////////
                        if (!X.isEmpty()) {
                            add_EM_F=true;
                            const double* X_p=const_cast<escript::Data*>(&X)->getSampleDataRO(e);
                            for (index_t k=0; k<numEq; k++) {
                                const register double X_0 = X_p[INDEX2(k, 0, numEq)];
                                const register double X_1 = X_p[INDEX2(k, 1, numEq)];
                                const register double X_2 = X_p[INDEX2(k, 2, numEq)];
                                const register double tmp1_1 = X_0*w19;
                                const register double tmp2_1 = X_1*w20;
                                const register double tmp3_1 = X_0*w22;
                                const register double tmp4_1 = X_1*w23;
                                const register double tmp5_1 = X_2*w24;
                                const register double tmp0_1 = X_2*w21;
                                EM_F[INDEX2(k,0,numEq)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                EM_F[INDEX2(k,1,numEq)]+=tmp0_1 + tmp2_1 + tmp3_1;
                                EM_F[INDEX2(k,2,numEq)]+=tmp0_1 + tmp1_1 + tmp4_1;
                                EM_F[INDEX2(k,3,numEq)]+=tmp0_1 + tmp3_1 + tmp4_1;
                                EM_F[INDEX2(k,4,numEq)]+=tmp1_1 + tmp2_1 + tmp5_1;
                                EM_F[INDEX2(k,5,numEq)]+=tmp2_1 + tmp3_1 + tmp5_1;
                                EM_F[INDEX2(k,6,numEq)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                EM_F[INDEX2(k,7,numEq)]+=tmp3_1 + tmp4_1 + tmp5_1;
                            }
                        }
                        ///////////////
                        // process Y //
                        ///////////////
                        if (!Y.isEmpty()) {
                            add_EM_F=true;
                            const double* Y_p=const_cast<escript::Data*>(&Y)->getSampleDataRO(e);
                            for (index_t k=0; k<numEq; k++) {
                                const register double Y_0 = Y_p[k];
                                const register double tmp0_1 = Y_0*w25;
                                EM_F[INDEX2(k,0,numEq)]+=tmp0_1;
                                EM_F[INDEX2(k,1,numEq)]+=tmp0_1;
                                EM_F[INDEX2(k,2,numEq)]+=tmp0_1;
                                EM_F[INDEX2(k,3,numEq)]+=tmp0_1;
                                EM_F[INDEX2(k,4,numEq)]+=tmp0_1;
                                EM_F[INDEX2(k,5,numEq)]+=tmp0_1;
                                EM_F[INDEX2(k,6,numEq)]+=tmp0_1;
                                EM_F[INDEX2(k,7,numEq)]+=tmp0_1;
                            }
                        }

                        // add to matrix (if add_EM_S) and RHS (if add_EM_F)
                        const index_t firstNode=m_N0*m_N1*k2+m_N0*k1+k0;
                        IndexVector rowIndex;
                        rowIndex.push_back(m_dofMap[firstNode]);
                        rowIndex.push_back(m_dofMap[firstNode+1]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0+1]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0*m_N1]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0*m_N1+1]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0*(m_N1+1)]);
                        rowIndex.push_back(m_dofMap[firstNode+m_N0*(m_N1+1)+1]);
                        if (add_EM_F) {
                            double *F_p=rhs.getSampleDataRW(0);
                            for (index_t i=0; i<rowIndex.size(); i++) {
                                if (rowIndex[i]<getNumDOF()) {
                                    for (index_t eq=0; eq<numEq; eq++) {
                                        F_p[INDEX2(eq,rowIndex[i],numEq)]+=EM_F[INDEX2(eq,i,numEq)];
                                    }
                                }
                            }
                        }
                        if (add_EM_S) {
                            addToSystemMatrix(mat, rowIndex, numEq, rowIndex, numComp, EM_S);
                        }
                    } // end k0 loop
                } // end k1 loop
            } // end k2 loop
        } // end of colouring
    } // end of parallel region
}

} // end of namespace ripley

