
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
#pragma omp parallel
    {
#pragma omp for
        for (dim_t i0 = 0; i0 < m_N0; i0++) {
            coords[0][i0]=(m_l0*(i0+m_offset0))/m_gNE0;
        }
#pragma omp for
        for (dim_t i1 = 0; i1 < m_N1; i1++) {
            coords[1][i1]=(m_l1*(i1+m_offset1))/m_gNE1;
        }
#pragma omp for
        for (dim_t i2 = 0; i2 < m_N2; i2++) {
            coords[2][i2]=(m_l2*(i2+m_offset2))/m_gNE2;
        }
    }
    int dims[] = { m_N0, m_N1, m_N2 };
    DBPutQuadmesh(dbfile, "mesh", NULL, coords, dims, 3, DB_DOUBLE,
            DB_COLLINEAR, NULL);

    DBPutQuadvar1(dbfile, "nodeId", "mesh", (void*)&m_nodeId[0], dims, 3, NULL, 0,
            DB_INT, DB_NODECENT, NULL);

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
    if (fsType == Nodes)
        return &m_nodeId[0];

    throw RipleyException("borrowSampleReferenceIDs() only implemented for Nodes");
}

bool Brick::ownSample(int fsCode, index_t id) const
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

void Brick::interpolateOnDomain(escript::Data& target,
                                const escript::Data& in) const
{
    const Brick& inDomain=dynamic_cast<const Brick&>(*(in.getFunctionSpace().getDomain()));
    const Brick& targetDomain=dynamic_cast<const Brick&>(*(target.getFunctionSpace().getDomain()));
    if (inDomain != *this)
        throw RipleyException("Illegal domain of interpolant");
    if (targetDomain != *this)
        throw RipleyException("Illegal domain of interpolation target");

    throw RipleyException("interpolateOnDomain() not implemented");
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
        for (index_t i=0; i < getNumNodes(); i++) {
            cout << "  " << setw(5) << m_nodeId[i]
                << "  " << (m_l0*(i%m_N0+m_offset0))/m_gNE0
                << "  " << (m_l1*(i%(m_N0*m_N1)/m_N0+m_offset1))/m_gNE1
                << "  " << (m_l2*(i/(m_N0*m_N1)+m_offset2))/m_gNE2 << endl;
        }
    }
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

    arg.requireWrite();
#pragma omp parallel for
    for (dim_t i2 = 0; i2 < m_N2; i2++) {
        for (dim_t i1 = 0; i1 < m_N1; i1++) {
            for (dim_t i0 = 0; i0 < m_N0; i0++) {
                double* point = arg.getSampleDataRW(i0+m_N0*i1+m_N0*m_N1*i2);
                point[0] = (m_l0*(i0+m_offset0))/m_gNE0;
                point[1] = (m_l1*(i1+m_offset1))/m_gNE1;
                point[2] = (m_l2*(i2+m_offset2))/m_gNE2;
            }
        }
    }
}

//private
void Brick::populateSampleIds()
{
    const index_t firstId = getNumNodes()*m_mpiInfo->rank;
    const index_t diff0 = m_N0*(m_N1*m_N2-1)+1;
    const index_t diff1 = m_N0*m_N1*(m_N2*m_NX-1)+m_N0;
    const index_t diff2 = m_N0*m_N1*m_N2*m_NX*m_NY-m_N0*m_N1*(m_N2-1);
    m_nodeId.resize(getNumNodes());
#pragma omp parallel for
    for (dim_t k=0; k<getNumNodes(); k++) {
        index_t id = firstId+k;
        if (m_offset0 > 0 && k%m_N0==0)
            id -= diff0;
        if (m_offset1 > 0 && k%(m_N0*m_N1)<m_N0)
            id -= diff1;
        if (m_offset2 > 0 && k/(m_N0*m_N1)==0)
            id -= diff2;
        m_nodeId[k]=id;
    }
}

} // end of namespace ripley

