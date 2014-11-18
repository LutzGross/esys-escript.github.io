
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

#include <algorithm>
#include <limits>

#include <speckley/Rectangle.h>
#include <esysUtils/esysFileWriter.h>
#include <esysUtils/index.h>
#include <esysUtils/Esys_MPI.h>
#include <speckley/DefaultAssembler2D.h>
#include <boost/scoped_array.hpp>
#include <escript/FunctionSpaceFactory.h>
#include "esysUtils/EsysRandom.h"

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

using esysUtils::FileWriter;

namespace speckley {

Rectangle::Rectangle(int order, dim_t n0, dim_t n1, double x0, double y0, double x1,
                     double y1, int d0, int d1,
                     const std::vector<double>& points,
                     const std::vector<int>& tags,
                     const simap_t& tagnamestonums,
		    escript::SubWorld_ptr w
		    ) :
    SpeckleyDomain(2, order, w)
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
        std::cout << "Warning: Automatic domain subdivision (d0=" << d0 << ", d1="
            << d1 << "). This may not be optimal!" << std::endl;
    }

    double l0 = x1-x0;
    double l1 = y1-y0;
    m_dx[0] = l0/n0;
    m_dx[1] = l1/n1;

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

    for (simap_t::const_iterator i = tagnamestonums.begin();
            i != tagnamestonums.end(); i++) {
        setTagMap(i->first, i->second);
    }
    addPoints(points, tags);
}

Rectangle::~Rectangle()
{
}

std::string Rectangle::getDescription() const
{
    return "speckley::Rectangle";
}

bool Rectangle::operator==(const AbstractDomain& other) const
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

const dim_t* Rectangle::borrowSampleReferenceIDs(int fsType) const
{
    switch (fsType) {
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom: // FIXME: reduced
//            return &m_dofId[0];
        case Nodes:
        case ReducedNodes: // FIXME: reduced
            return &m_nodeId[0];
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
    throw SpeckleyException(msg.str());
}

bool Rectangle::ownSample(int fsType, index_t id) const
{
    throw SpeckleyException("ownSample not implemented");
}

void Rectangle::setToNormal(escript::Data& out) const
{
    throw SpeckleyException("setToNormal not implemented");
}

void Rectangle::setToSize(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == Elements) {
        out.requireWrite();
        const dim_t numQuad = out.getNumDataPointsPerSample();
        const dim_t numElements = getNumElements();
        const double size=sqrt(m_dx[0]*m_dx[0]+m_dx[1]*m_dx[1]);
#pragma omp parallel for
        for (index_t k = 0; k < numElements; ++k) {
            double* o = out.getSampleDataRW(k);
            std::fill(o, o+numQuad, size);
        }
    } else if (out.getFunctionSpace().getTypeCode() == FaceElements) {
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
                    std::fill(o, o+numQuad, m_dx[1]);
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k1 = 0; k1 < NE1; ++k1) {
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k1);
                    std::fill(o, o+numQuad, m_dx[1]);
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < NE0; ++k0) {
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k0);
                    std::fill(o, o+numQuad, m_dx[0]);
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k0 = 0; k0 < NE0; ++k0) {
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k0);
                    std::fill(o, o+numQuad, m_dx[0]);
                }
            }
        } // end of parallel section

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
    escriptDataC x = arg.getDataC();
    int numDim = m_numDim;
    if (!isDataPointShapeEqual(&x, 1, &numDim))
        throw SpeckleyException("setToX: Invalid Data object shape");
    if (!numSamplesEqual(&x, 1, getNumNodes()))
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
void Rectangle::assembleIntegrate(std::vector<double>& integrals,
                                  const escript::Data& arg) const
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

//protected
dim_t Rectangle::insertNeighbourNodes(IndexVector& index, index_t node) const
{
    std::cerr << "insertNeighbourNodes not updated\n";
    const dim_t nDOF0 = (m_gNE[0]+1)/m_NX[0];
    const dim_t nDOF1 = (m_gNE[1]+1)/m_NX[1];
    const int x=node%nDOF0;
    const int y=node/nDOF0;
    dim_t num=0;
    // loop through potential neighbours and add to index if positions are
    // within bounds
    for (int i1=-1; i1<2; i1++) {
        for (int i0=-1; i0<2; i0++) {
            // skip node itself
            if (i0==0 && i1==0)
                continue;
            // location of neighbour node
            const int nx=x+i0;
            const int ny=y+i1;
            if (nx>=0 && ny>=0 && nx<nDOF0 && ny<nDOF1) {
                index.push_back(ny*nDOF0+nx);
                num++;
            }
        }
    }

    return num;
}

/* This is a wrapper for filtered (and non-filtered) randoms
 * For detailed doco see randomFillWorker
*/ 
escript::Data Rectangle::randomFill(const escript::DataTypes::ShapeType& shape,
           const escript::FunctionSpace& fs,
           long seed, const boost::python::tuple& filter) const {
    int numvals=escript::DataTypes::noValues(shape);
    int per_element = (m_order+1)*(m_order+1)*numvals;
    if (len(filter)>0) {
        throw SpeckleyException("Speckley does not support filters.");
    }

    double* src=new double[m_NE[0]*m_NE[1]*per_element*numvals];
    esysUtils::randomFillArray(seed, src, m_NE[0]*m_NE[1]*per_element);
    escript::Data res(0, shape, escript::function(*this), true);
    int current = 0;
    for (int ei = 0; ei < m_NE[1]; ++ei) {
        for (int ej = 0; ej < m_NE[0]; ++ej) {
            double *e = res.getSampleDataRW(INDEX2(ej,ei,m_NE[0]));
            memcpy(e, &src[current], sizeof(double)*per_element);
            current += per_element;
        }
    }
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

    const dim_t NFE = getNumFaceElements();
    m_faceId.resize(NFE);


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
}

//private
void Rectangle::addToMatrixAndRHS(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<double>& EM_S, const std::vector<double>& EM_F, bool addS,
         bool addF, index_t firstNode, dim_t nEq, dim_t nComp) const
{
    throw SpeckleyException("Rectangle::addToMatrixAndRHS, adding to matrix not supported");
}

//protected
void Rectangle::interpolateNodesOnElements(escript::Data& out,
                                           const escript::Data& in,
                                           bool reduced) const
{
    if (reduced) {
        throw SpeckleyException("Speckley domains do not support reduced function spaces");
    }
    const dim_t numComp = in.getDataPointSize();
    const dim_t NE0 = m_NE[0];
    const dim_t NE1 = m_NE[1];
    const int quads = m_order + 1;
    const int max_x = m_NN[0];
    out.requireWrite();
#pragma omp parallel for
    for (dim_t ey = 0; ey < NE1; ey++) {
        for (dim_t ex = 0; ex < NE0; ex++) {
            double *e_out = out.getSampleDataRW(ex + ey*NE0);
            int start = ex*m_order + ey*max_x*m_order;
            int quad = 0;
            for (int qy = 0; qy < quads; qy++) {
                for (int qx = 0; qx < quads; qx++, quad++) {
                    const double *n_in = in.getSampleDataRO(start + max_x*qy + qx);
                    memcpy(e_out+quad*numComp, n_in, sizeof(double) * numComp);
                }
            }
        }
    }
}

#ifdef ESYS_MPI
//protected
void Rectangle::balanceNeighbours(escript::Data& data, bool average) const {
    if (m_NX[0] * m_NX[1] == 1) {
        return;
    }
    const int numComp = data.getDataPointSize();
    const int rx = m_mpiInfo->rank % m_NX[0];
    const int ry = m_mpiInfo->rank / m_NX[0];
    //include bordering ranks in summation
    if (m_NX[1] != 1)
        shareVertical(data, rx, ry);
    if (m_NX[0] != 1)
        shareSides(data, rx, ry);
    if (m_NX[0] != 1 && m_NX[1] != 1) {
        shareCorners(data, rx, ry);
        if (!average)
            return;
        //averaging out corners
        // bottom left
        if (rx && ry) {
            double *values = data.getSampleDataRW(0);
            for (int comp = 0; comp < numComp; comp++) {
                values[comp] /= 2;
            }
        }
        // bottom right
        if (rx < (m_NX[0] - 1) && ry) {
            double *values = data.getSampleDataRW(m_NN[0]-1);
            for (int comp = 0; comp < numComp; comp++) {
                values[comp] /= 2;
            }
        }
        // top left
        if (rx && ry < (m_NX[0] - 1)) {
            double *values = data.getSampleDataRW((m_NN[1]-1)*m_NN[0]);
            for (int comp = 0; comp < numComp; comp++) {
                values[comp] /= 2;
            }
        }
        // top right
        if (rx < (m_NX[0] - 1) && ry < (m_NX[0] - 1)) {
            double *values = data.getSampleDataRW(m_NN[1]*m_NN[0] - 1);
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
            double *values = data.getSampleDataRW(qy*m_NN[0]);
            for (int comp = 0; comp < numComp; comp++) {
                values[comp] /= 2;
            }
        }
    }
    //right
    if (rx < m_NX[0] - 1) {
#pragma omp parallel for
        for (dim_t qy = 0; qy < m_NN[1]; qy++) {
            double *values = data.getSampleDataRW(qy*m_NN[0] + m_NN[0] - 1);
            for (int comp = 0; comp < numComp; comp++) {
                values[comp] /= 2;
            }
        }
    }
    //bottom
    if (ry) {
#pragma omp parallel for
        for (dim_t qx = 0; qx < m_NN[0]; qx++) {
            double *values = data.getSampleDataRW(qx);
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
            double *values = data.getSampleDataRW(start + qx);
            for (int comp = 0; comp < numComp; comp++) {
                values[comp] /= 2;
            }
        }
    }
}
#endif //#ifdef ESYS_MPI

//protected
void Rectangle::interpolateElementsOnNodes(escript::Data& out,
        const escript::Data& in, bool reduced) const {
    const dim_t numComp = in.getDataPointSize();
    const dim_t NE0 = m_NE[0];
    const dim_t NE1 = m_NE[1];
    const int quads = m_order + 1;
    const dim_t max_x = (m_order*NE0) + 1;
    const dim_t max_y = (m_order*NE1) + 1;
    out.requireWrite();
    if (reduced) {
        throw SpeckleyException("Speckley domains do not support reduced function spaces");
    }
    // the summation portion
    for (dim_t colouring = 0; colouring < 2; colouring++) {
#pragma omp parallel for
        for (dim_t ey = colouring; ey < NE1; ey += 2) {
            for (dim_t ex = 0; ex < NE0; ex++) {
                int start = ex*m_order + ey*max_x*m_order;
                const double *e_in = in.getSampleDataRO(ex + ey*NE0);
                for (int qy = 0; qy < quads; qy++) {
                    for (int qx = 0; qx < quads; qx++) {
                        double *n_out = out.getSampleDataRW(start + max_x*qy + qx);
                        for (int comp = 0; comp < numComp; comp++) {
                            n_out[comp] += e_in[INDEX3(comp, qx, qy, numComp, quads)];
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
            double *n_out = out.getSampleDataRW(qx + qy*max_x);
            for (int comp = 0; comp < numComp; comp++) {
                n_out[comp] /= 2;
            }
        }
    }

    // for every internal edge in y
#pragma omp parallel for
    for (dim_t qy = m_order; qy < max_y - m_order; qy += m_order) {
        for (dim_t qx = 0; qx < max_x; qx ++) {
            double *n_out = out.getSampleDataRW(qx + qy*max_x);
            for (int comp = 0; comp < numComp; comp++) {
                n_out[comp] /= 2;
            }
        }
    }
}

#ifdef ESYS_MPI
//private
void Rectangle::shareCorners(escript::Data& out, int rx, int ry) const
{
    //setup
    const int tag = 0;
    MPI_Status status;
    const int numComp = out.getDataPointSize();
    const int count = 4 * numComp;
    std::vector<double> outbuf(count, 0);
    std::vector<double> inbuf(count, 0);
    const int rank = m_mpiInfo->rank;
    //precalc bounds so we can loop nicely, can probably be cleaned up
    const bool conds[4] = {rx && ry, rx < (m_NX[0] - 1) && ry,
            rx && ry < (m_NX[1] - 1),rx < (m_NX[0] - 1) && ry < (m_NX[1] - 1)};
    const int ranks[4] = {rank-m_NX[0]-1,rank-m_NX[0]+1,
                        rank+m_NX[0]-1,rank+m_NX[0]+1};
    //fill everything, regardless of whether we're sharing that corner or not
    for (int y = 0; y < 2; y++) {
        for (int x = 0; x < 2; x++) {
            const double *data = out.getSampleDataRO(x*(m_NN[0]-1) + y*(m_NN[1]-1)*m_NN[0]);
            std::copy(data, data + numComp, &outbuf[(x + 2*y)*numComp]);
        }
    }
    
    //share
    for (int i = 0; i < 4; i++) {
        if (conds[i]) {
            if (rx % 2) {
                MPI_Send(&outbuf[i], numComp, MPI_DOUBLE, ranks[i], tag,
                        m_mpiInfo->comm);
                MPI_Recv(&inbuf[i], count, MPI_DOUBLE, ranks[i], tag,
                        m_mpiInfo->comm, &status);
            } else {
                MPI_Recv(&inbuf[i], count, MPI_DOUBLE, ranks[i], tag,
                        m_mpiInfo->comm, &status);
                MPI_Send(&outbuf[i], numComp, MPI_DOUBLE, ranks[i], tag,
                        m_mpiInfo->comm);
            }
        }
    }

    //unpack
    for (int y = 0; y < 2; y++) {
        for (int x = 0; x < 2; x++) {
            double *data = out.getSampleDataRW(x*(m_NN[0]-1) + y*(m_NN[1]-1)*m_NN[0]);
            for (int comp = 0; comp < numComp; comp++) {
                data[comp] += inbuf[(x + 2*y)*numComp + comp];
            }
        }
    }
}

//private
void Rectangle::shareVertical(escript::Data& out, int rx, int ry) const
{
    const int tag = 0;
    MPI_Status status;
    const int numComp = out.getDataPointSize();
    const dim_t count = m_NN[0]*numComp;
    const int up_neighbour = m_mpiInfo->rank + m_NX[0];
    const int down_neighbour = m_mpiInfo->rank - m_NX[0];
    //allocate some space to recieve
    std::vector<double> recv(count);
    //get our sources
    double *top = out.getSampleDataRW((m_NN[1]-1) * m_NN[0]);
    double *bottom = out.getSampleDataRW(0);
    
    if (ry % 2) {
        //send to up, read from up
        if (ry < m_NX[1] - 1) {
            MPI_Send(&top[0], count, MPI_DOUBLE, up_neighbour, tag,
                    m_mpiInfo->comm);
            MPI_Recv(&recv[0], count, MPI_DOUBLE, up_neighbour, tag,
                    m_mpiInfo->comm, &status);
            //unpack up
            for (dim_t i = 0; i < count; i++) {
                top[i] += recv[i];
            }
        }
        //send down, read down
        MPI_Send(&bottom[0], count, MPI_DOUBLE, down_neighbour, tag,
                    m_mpiInfo->comm);
        MPI_Recv(&recv[0], count, MPI_DOUBLE, down_neighbour, tag,
                    m_mpiInfo->comm, &status);
        //unpack bottom
        for (dim_t i = 0; i < count; i++) {
            bottom[i] += recv[i];
        }
    } else {
        //read down, send down
        if (ry) {
            MPI_Recv(&recv[0], count, MPI_DOUBLE, down_neighbour, tag,
                    m_mpiInfo->comm, &status);
            MPI_Send(&bottom[0], count, MPI_DOUBLE, down_neighbour, tag,
                    m_mpiInfo->comm);
            //unpack bottom
            for (dim_t i = 0; i < count; i++) {
                bottom[i] += recv[i];
            }
        }

        //read up, send up
        if (ry < m_NX[1] - 1) {
            MPI_Recv(&recv[0], count, MPI_DOUBLE, up_neighbour, tag,
                    m_mpiInfo->comm, &status);
            MPI_Send(&top[0], count, MPI_DOUBLE, up_neighbour, tag,
                    m_mpiInfo->comm);
            //unpack up
            for (dim_t i = 0; i < count; i++) {
                top[i] += recv[i];
            }
        }
    }
}

//private
void Rectangle::shareSides(escript::Data& out, int rx, int ry) const
{
    const int tag = 0;
    MPI_Status status;
    const int numComp = out.getDataPointSize();
    const int count = m_NN[1]*numComp;
    const int left_neighbour = m_mpiInfo->rank - 1;
    const int right_neighbour = m_mpiInfo->rank + 1;
    //allocate some space
    std::vector<double> left(count,170);
    std::vector<double> right(count,17000);
    std::vector<double> recv(count,1700);
    //fill those values
    for (dim_t n = 0; n < m_NN[1]; n++) {
        index_t index = n*m_NN[0];
        const double *leftData = out.getSampleDataRO(index);
        std::copy(leftData, leftData + numComp, &left[n*numComp]);
        const double *rightData = out.getSampleDataRO(index+m_NN[0]-1);
        std::copy(rightData, rightData + numComp, &right[n*numComp]);
    }
    
    if (rx % 2) {
        //send to right, read from right
        if (rx < m_NX[0] - 1) {
            MPI_Send(&right[0], count, MPI_DOUBLE, right_neighbour, tag,
                    m_mpiInfo->comm);
            MPI_Recv(&recv[0], count, MPI_DOUBLE, right_neighbour, tag,
                    m_mpiInfo->comm, &status);
            //unpack to right
            for (dim_t i = 0; i < m_NN[1]; i++) {
                double *data = out.getSampleDataRW((i + 1) * m_NN[0] - 1);
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += recv[i*numComp+comp];
                }
            }
        }
        //send left, read left
        MPI_Send(&left[0], count, MPI_DOUBLE, left_neighbour, tag,
                    m_mpiInfo->comm);
        MPI_Recv(&recv[0], count, MPI_DOUBLE, left_neighbour, tag,
                    m_mpiInfo->comm, &status);
        //unpack to left
        for (dim_t i = 0; i < m_NN[1]; i++) {
            double *data = out.getSampleDataRW(i*m_NN[0]);
            for (int comp = 0; comp < numComp; comp++) {
                data[comp] += recv[i*numComp+comp];
            }
        }
    } else {
        //read left, send left
        if (rx) {
            MPI_Recv(&recv[0], count, MPI_DOUBLE, left_neighbour, tag,
                    m_mpiInfo->comm, &status);
            MPI_Send(&left[0], count, MPI_DOUBLE, left_neighbour, tag,
                    m_mpiInfo->comm);
            //unpack to left
            for (dim_t i = 0; i < m_NN[1]; i++) {
                double *data = out.getSampleDataRW(i*m_NN[0]);
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += recv[i*numComp+comp];
                }
            }
        }

        //read right, send right
        if (rx < m_NX[0] - 1) {
            MPI_Recv(&recv[0], count, MPI_DOUBLE, right_neighbour, tag,
                    m_mpiInfo->comm, &status);
            MPI_Send(&right[0], count, MPI_DOUBLE, right_neighbour, tag,
                    m_mpiInfo->comm);
            //unpack to right
            for (dim_t i = 0; i < m_NN[1]; i++) {
                double *data = out.getSampleDataRW((i + 1) * m_NN[0] - 1);
                for (int comp = 0; comp < numComp; comp++) {
                    data[comp] += recv[i*numComp+comp];
                }
            }
        }
    }
}
#endif //#ifdef ESYS_MPI

//protected
void Rectangle::interpolateNodesOnFaces(escript::Data& out,
                                        const escript::Data& in,
                                        bool reduced) const
{
    throw SpeckleyException("Rectangle::interpolateNodesOnFaces not implemented");
}



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

    // distance in elements
    dim_t ex = (dim_t) floor(x / m_dx[0] + 0.01*m_dx[0]);
    dim_t ey = (dim_t) floor(y / m_dx[1] + 0.01*m_dx[1]);
    dim_t start = ex*m_order + ey*m_order*m_NN[0];
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
                closest = start + dx*m_order + dy*m_NN[0]*m_order;
				std::cerr << "Rectangle::findNode not updated for MPI\n";
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
        const DataMap& options) const {
    if (type.compare("DefaultAssembler") == 0) {
        return Assembler_ptr(new DefaultAssembler2D(shared_from_this(), m_dx, m_NX, m_NE, m_NN));
    }
    throw SpeckleyException("Speckley::Rectangle does not support the"
            " requested assembler");
}

} // end of namespace speckley
