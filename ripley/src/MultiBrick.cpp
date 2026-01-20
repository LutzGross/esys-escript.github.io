
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

#include <ripley/MultiBrick.h>
#include <ripley/blocktools.h>
#include <ripley/domainhelpers.h>

#include <escript/DataFactory.h>
#include <escript/FunctionSpaceFactory.h>
#include <escript/index.h>

#define FIRST_QUAD 0.21132486540518711775
#define SECOND_QUAD 0.78867513459481288225

#include <iomanip>
#include <iostream>
#include <limits>

using std::vector;
using std::string;
using std::min;
using std::ios;
using std::fill;

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

MultiBrick::MultiBrick(dim_t n0, dim_t n1, dim_t n2, double x0, double y0, double z0,
             double x1, double y1, double z1, int d0, int d1, int d2,
             const vector<double>& points, const vector<int>& tags,
             const TagMap& tagnamestonums,
             unsigned int subdivisions
         ) :
    Brick(n0, n1, n2, x0, y0, z0, x1, y1, z1, d0, d1, d2, points, tags, tagnamestonums),
    m_subdivisions(subdivisions)
{
    if (m_mpiInfo->size != 1)
        throw RipleyException("Multiresolution Brick domains don't currently support multiple processes");

    if (subdivisions == 0 || (subdivisions & (subdivisions - 1)) != 0)
        throw RipleyException("Element subdivisions must be a power of two");

    if (d0 == 0 || d1 == 0)
        throw RipleyException("Domain subdivisions must be positive");

    dim_t oldNN[3] = {0};

    for (int i = 0; i < 3; i++) {
        m_NE[i] *= subdivisions;
        oldNN[i] = m_NN[i];
        m_NN[i] = m_NE[i] + 1;
        m_gNE[i] *= subdivisions;
        m_ownNE[i] *= subdivisions;
        m_dx[i] /= subdivisions;
        m_faceCount[i] *= subdivisions;
        m_faceCount[2+i] *= subdivisions;
    }

    // bottom-left node is at (offset0,offset1) in global mesh
    m_offset[0] = m_gNE[0]*subdivisions/d0*(m_mpiInfo->rank%d0);
    m_offset[1] = m_gNE[1]*subdivisions/d1*(m_mpiInfo->rank/d0);
    m_offset[2] = m_gNE[2]*subdivisions/d2*(m_mpiInfo->rank/(d0*d1));
    populateSampleIds();
    
    const dim_t nDirac = m_diracPoints.size();
#pragma omp parallel for
    for (int i = 0; i < nDirac; i++) {
        const dim_t node = m_diracPoints[i].node;
        const dim_t x = node % oldNN[0];
        const dim_t y = (node % (oldNN[0]*oldNN[1])) / oldNN[0];
        const dim_t z = node / (oldNN[0]*oldNN[1]);
        m_diracPoints[i].node = INDEX3(x*subdivisions, y*subdivisions, z*subdivisions,
                m_NN[0], m_NN[1]);
        m_diracPointNodeIDs[i] = m_diracPoints[i].node;
    }
}

MultiBrick::MultiBrick(escript::JMPI jmpi, dim_t n0, dim_t n1, dim_t n2, double x0, double y0, double z0,
             double x1, double y1, double z1, int d0, int d1, int d2,
             const vector<double>& points, const vector<int>& tags,
             const TagMap& tagnamestonums,
             unsigned int subdivisions
         ) :
    Brick(jmpi, n0, n1, n2, x0, y0, z0, x1, y1, z1, d0, d1, d2, points, tags, tagnamestonums),
    m_subdivisions(subdivisions)
{
    if (m_mpiInfo->size != 1)
        throw RipleyException("Multiresolution Brick domains don't currently support multiple processes");

    if (subdivisions == 0 || (subdivisions & (subdivisions - 1)) != 0)
        throw RipleyException("Element subdivisions must be a power of two");

    if (d0 == 0 || d1 == 0)
        throw RipleyException("Domain subdivisions must be positive");

    dim_t oldNN[3] = {0};

    for (int i = 0; i < 3; i++) {
        m_NE[i] *= subdivisions;
        oldNN[i] = m_NN[i];
        m_NN[i] = m_NE[i] + 1;
        m_gNE[i] *= subdivisions;
        m_ownNE[i] *= subdivisions;
        m_dx[i] /= subdivisions;
        m_faceCount[i] *= subdivisions;
        m_faceCount[2+i] *= subdivisions;
    }

    // bottom-left node is at (offset0,offset1) in global mesh
    m_offset[0] = m_gNE[0]*subdivisions/d0*(m_mpiInfo->rank%d0);
    m_offset[1] = m_gNE[1]*subdivisions/d1*(m_mpiInfo->rank/d0);
    m_offset[2] = m_gNE[2]*subdivisions/d2*(m_mpiInfo->rank/(d0*d1));
    populateSampleIds();

    const dim_t nDirac = m_diracPoints.size();
#pragma omp parallel for
    for (int i = 0; i < nDirac; i++) {
        const dim_t node = m_diracPoints[i].node;
        const dim_t x = node % oldNN[0];
        const dim_t y = (node % (oldNN[0]*oldNN[1])) / oldNN[0];
        const dim_t z = node / (oldNN[0]*oldNN[1]);
        m_diracPoints[i].node = INDEX3(x*subdivisions, y*subdivisions, z*subdivisions,
                m_NN[0], m_NN[1]);
        m_diracPointNodeIDs[i] = m_diracPoints[i].node;
    }
}

MultiBrick::~MultiBrick()
{
}


void MultiBrick::validateInterpolationAcross(int fsType_source,
        const escript::AbstractDomain& domain, int fsType_target) const
{
    const MultiBrick *other = dynamic_cast<const MultiBrick *>(&domain);
    if (other == NULL)
        throw RipleyException("Invalid interpolation: domains must both be instances of MultiBrick");

    const double *len = other->getLength();
    const int *subdivs = other->getNumSubdivisionsPerDim();
    const dim_t *elements = other->getNumElementsPerDim();
    const unsigned int level = other->getNumSubdivisionsPerElement();
    const unsigned int factor = m_subdivisions > level ? m_subdivisions/level : level/m_subdivisions;
    if ((factor & (factor - 1)) != 0) //factor == 2**x
        throw RipleyException("Invalid interpolation: elemental subdivisions of each domain must be powers of two");

    if (other->getMPIComm() != m_mpiInfo->comm)
        throw RipleyException("Invalid interpolation: Domains are on different communicators");

    for (int i = 0; i < m_numDim; i++) {
        if (m_length[i] != len[i])
            throw RipleyException("Invalid interpolation: domain length mismatch");
        if (m_NX[i] != subdivs[i])
            throw RipleyException("Invalid interpolation: domain process subdivision mismatch");
        if (m_subdivisions > level) {
            if (m_ownNE[i]/elements[i] != factor)
                throw RipleyException("Invalid interpolation: element factor mismatch");
        } else {
            if (elements[i]/m_ownNE[i] != factor)
                throw RipleyException("Invalid interpolation: element factor mismatch");
        }
    }
}

void MultiBrick::interpolateNodesToNodesFiner(const escript::Data& source,
        escript::Data& target, const MultiBrick& other) const
{
    if (source.isComplex()!=target.isComplex())
    {
        throw RipleyException("Programmer Error: in and out parameters do not have the same complexity.");        
    }
    if (source.isComplex())
    {
        interpolateNodesToNodesFinerWorker(source, target, other, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateNodesToNodesFinerWorker(source, target, other, escript::DataTypes::real_t(0));      
    }
}

template <typename S>
void MultiBrick::interpolateNodesToNodesFinerWorker(const escript::Data& source,
        escript::Data& target, const MultiBrick& other, S sentinel) const
{
    const int scaling = other.getNumSubdivisionsPerElement()/m_subdivisions;
    const dim_t NN0 = m_NN[0], NN1 = m_NN[1], NN2 = m_NN[2], *otherNN = other.getNumNodesPerDim();
    const dim_t numComp = source.getDataPointSize();
    target.requireWrite();
#pragma omp parallel for
    for (dim_t nz = 0; nz < NN2 - 1; nz++) { //source nodes
        for (dim_t ny = 0; ny < NN1 - 1; ny++) {
            for (dim_t nx = 0; nx < NN0 - 1; nx++) {
                const S *x0y0z0 = source.getSampleDataRO(INDEX3(nx,   ny,   nz, NN0, NN1), sentinel);
                const S *x0y1z0 = source.getSampleDataRO(INDEX3(nx,   ny+1, nz, NN0, NN1), sentinel);
                const S *x1y0z0 = source.getSampleDataRO(INDEX3(nx+1, ny,   nz, NN0, NN1), sentinel);
                const S *x1y1z0 = source.getSampleDataRO(INDEX3(nx+1, ny+1, nz, NN0, NN1), sentinel);
                const S *x0y0z1 = source.getSampleDataRO(INDEX3(nx,   ny,   nz+1, NN0, NN1), sentinel);
                const S *x0y1z1 = source.getSampleDataRO(INDEX3(nx,   ny+1, nz+1, NN0, NN1), sentinel);
                const S *x1y0z1 = source.getSampleDataRO(INDEX3(nx+1, ny,   nz+1, NN0, NN1), sentinel);
                const S *x1y1z1 = source.getSampleDataRO(INDEX3(nx+1, ny+1, nz+1, NN0, NN1), sentinel);
                const S origin[3] = {getLocalCoordinate(nx, 0),
                                          getLocalCoordinate(ny, 1),
                                          getLocalCoordinate(nz, 2)};
                for (int sz = 0; sz < scaling + 1; sz++) { //target nodes                
                    const S z = (other.getLocalCoordinate(nz*scaling+sz, 2) - origin[2]) / m_dx[2];
		    const S mz = static_cast<S>(1)-z;
                    for (int sy = 0; sy < scaling + 1; sy++) {
                        const S y = (other.getLocalCoordinate(ny*scaling+sy, 1) - origin[1]) / m_dx[1];
			const S my = static_cast<S>(1)-y;
                        for (int sx = 0; sx < scaling + 1; sx++) {
                            const S x = (other.getLocalCoordinate(nx*scaling+sx, 0) - origin[0]) / m_dx[0];
			    const S mx = static_cast<S>(1)-x;
                            S *out = target.getSampleDataRW(INDEX3(nx*scaling+sx, ny*scaling+sy, nz*scaling+sz, otherNN[0], otherNN[1]), sentinel);
			    
                            for (int comp = 0; comp < numComp; comp++) {
                                out[comp] = x0y0z0[comp]*(mx)*(my)*(mz) + x1y0z0[comp]*x*(my)*(mz) + x0y1z0[comp]*(mx)*y*(mz) + x1y1z0[comp]*x*y*(mz)
                                          + x0y0z1[comp]*(mx)*(my)*z     + x1y0z1[comp]*x*(my)*z     + x0y1z1[comp]*(mx)*y*z     + x1y1z1[comp]*x*y*z;
                            }
                        }
                    }
                }
            }
        }
    }
}

void MultiBrick::interpolateReducedToElementsFiner(const escript::Data& source,
        escript::Data& target, const MultiBrick& other) const
{
    if (source.isComplex()!=target.isComplex())
    {
        throw RipleyException("Programmer Error: in and out parameters do not have the same complexity.");        
    }
    if (source.isComplex())
    {
        interpolateReducedToElementsFinerWorker(source, target, other, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateReducedToElementsFinerWorker(source, target, other, escript::DataTypes::real_t(0));      
    }
}        

template <typename S>
void MultiBrick::interpolateReducedToElementsFinerWorker(const escript::Data& source,
        escript::Data& target, const MultiBrick& other, S sentinel) const
{
    const int scaling = other.getNumSubdivisionsPerElement()/m_subdivisions;
    const dim_t numComp = source.getDataPointSize();
    target.requireWrite();
    //for each of ours
#pragma omp parallel for
    for (dim_t ez = 0; ez < m_NE[2]; ez++) {
        for (dim_t ey = 0; ey < m_NE[1]; ey++) {
            for (dim_t ex = 0; ex < m_NE[0]; ex++) {
                const S *in = source.getSampleDataRO(INDEX3(ex, ey, ez, m_NE[0], m_NE[1]), sentinel);
                //for each subelement
                for (dim_t sz = 0; sz < scaling; sz++) {
                    const dim_t tz = ez*scaling + sz;
                    for (dim_t sy = 0; sy < scaling; sy++) {
                        const dim_t ty = ey*scaling + sy;
                        for (dim_t sx = 0; sx < scaling; sx++) {
                            const dim_t tx = ex*scaling + sx;
                            S *out = target.getSampleDataRW(INDEX3(tx, ty, tz, m_NE[0]*scaling, m_NE[1]*scaling), sentinel);
                            for (dim_t comp = 0; comp < numComp; comp++) {
                                const S quadvalue = in[comp];
                                for (int i = 0; i < 8; i++) {
                                    out[comp + i*numComp] = quadvalue;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void MultiBrick::interpolateReducedToReducedFiner(const escript::Data& source,
        escript::Data& target, const MultiBrick& other) const
{
    if (source.isComplex()!=target.isComplex())
    {
        throw RipleyException("Programmer Error: in and out parameters do not have the same complexity.");        
    }
    if (source.isComplex())
    {
        interpolateReducedToReducedFinerWorker(source, target, other, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateReducedToReducedFinerWorker(source, target, other, escript::DataTypes::real_t(0));      
    }
}   



template <typename S>
void MultiBrick::interpolateReducedToReducedFinerWorker(const escript::Data& source,
        escript::Data& target, const MultiBrick& other, S sentinel) const
{
    const int scaling = other.getNumSubdivisionsPerElement()/m_subdivisions;
    const dim_t numComp = source.getDataPointSize();
    target.requireWrite();
    //for each of ours
#pragma omp parallel for
    for (dim_t ey = 0; ey < m_NE[1]; ey++) {
        for (dim_t ex = 0; ex < m_NE[0]; ex++) {
            const S *in = source.getSampleDataRO(ex + ey*m_NE[0], sentinel);
            //for each subelement
            for (dim_t sy = 0; sy < scaling; sy++) {
                const dim_t ty = ey*scaling + sy;
                for (dim_t sx = 0; sx < scaling; sx++) {
                    const dim_t tx = ex*scaling + sx;
                    S *out = target.getSampleDataRW(tx + ty*m_NE[0]*scaling, sentinel);
                    for (dim_t comp = 0; comp < numComp; comp++) {
                        out[comp] = in[comp];
                    }
                }
            }
        }
    }
}

void MultiBrick::interpolateNodesToElementsFiner(const escript::Data& source,
        escript::Data& target, const MultiBrick& other) const
{
    if (source.isComplex()!=target.isComplex())
    {
        throw RipleyException("Programmer Error: in and out parameters do not have the same complexity.");        
    }
    if (source.isComplex())
    {
        interpolateNodesToElementsFinerWorker(source, target, other, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateNodesToElementsFinerWorker(source, target, other, escript::DataTypes::real_t(0));      
    }
}   

template <typename S>
void MultiBrick::interpolateNodesToElementsFinerWorker(const escript::Data& source,
        escript::Data& target, const MultiBrick& other, S sentinel) const
{
    const int scaling = other.getNumSubdivisionsPerElement()/m_subdivisions;
    const dim_t NE0 = m_NE[0], NE1 = m_NE[1], NE2 = m_NE[2], *theirNE = other.getNumElementsPerDim();
    const dim_t numComp = source.getDataPointSize();
    target.requireWrite();
#pragma omp parallel for
    for (dim_t ez = 0; ez < NE2; ez++) { //source nodes
        for (dim_t ey = 0; ey < NE1; ey++) {
            for (dim_t ex = 0; ex < NE0; ex++) {
                const S *points[8] = {
                        source.getSampleDataRO(INDEX3(ex,   ey,   ez, NE0+1, NE1+1), sentinel),
                        source.getSampleDataRO(INDEX3(ex+1, ey,   ez, NE0+1, NE1+1), sentinel),
                        source.getSampleDataRO(INDEX3(ex,   ey+1, ez, NE0+1, NE1+1), sentinel),
                        source.getSampleDataRO(INDEX3(ex+1, ey+1, ez, NE0+1, NE1+1), sentinel),
                        source.getSampleDataRO(INDEX3(ex,   ey,   ez+1, NE0+1, NE1+1), sentinel),
                        source.getSampleDataRO(INDEX3(ex+1, ey,   ez+1, NE0+1, NE1+1), sentinel),
                        source.getSampleDataRO(INDEX3(ex,   ey+1, ez+1, NE0+1, NE1+1), sentinel),
                        source.getSampleDataRO(INDEX3(ex+1, ey+1, ez+1, NE0+1, NE1+1), sentinel),
                    };
                const S origin[3] = {getLocalCoordinate(ex, 0),
                                          getLocalCoordinate(ey, 1),
                                          getLocalCoordinate(ez, 2)
                                          };
                for (int sz = 0; sz < scaling; sz++) { //target elements
                    for (int sy = 0; sy < scaling; sy++) {
                        for (int sx = 0; sx < scaling; sx++) {
                            const S x1 = (other.getLocalCoordinate(ex*scaling+sx, 0) - origin[0]) / m_dx[0] + FIRST_QUAD/scaling;
                            const S x2 = x1 + (SECOND_QUAD - FIRST_QUAD)/scaling;
                            const S y1 = (other.getLocalCoordinate(ey*scaling+sy, 1) - origin[1]) / m_dx[1] + FIRST_QUAD/scaling;
                            const S y2 = y1 + (SECOND_QUAD - FIRST_QUAD)/scaling;
                            const S z1 = (other.getLocalCoordinate(ez*scaling+sz, 2) - origin[2]) / m_dx[2] + FIRST_QUAD/scaling;
                            const S z2 = z1 + (SECOND_QUAD - FIRST_QUAD)/scaling;
			    const S mx1 = static_cast<S>(1)-x1;
			    const S mx2 = static_cast<S>(1)-x2;
			    const S my1 = static_cast<S>(1)-y1;
			    const S my2 = static_cast<S>(1)-y2;
			    const S mz1 = static_cast<S>(1)-z1;
			    const S mz2 = static_cast<S>(1)-z2;
			    
                            S *out = target.getSampleDataRW(INDEX3(ex*scaling+sx, ey*scaling+sy, ez*scaling+sz, theirNE[0], theirNE[1]), sentinel);
                            for (int comp = 0; comp < numComp; comp++) {
                                out[INDEX4(comp, 0, 0, 0, numComp, 2, 2)] = 
                                              points[0][comp]*(mx1)*(my1)*(mz1)
                                            + points[1][comp]*x1*(my1)*(mz1) 
                                            + points[2][comp]*(mx1)*y1*(mz1)
                                            + points[3][comp]*x1*y1*(mz1)
                                            + points[4][comp]*(mx1)*(my1)*z1
                                            + points[5][comp]*x1*(my1)*z1 
                                            + points[6][comp]*(mx1)*y1*z1
                                            + points[7][comp]*x1*y1*z1;
                                out[INDEX4(comp, 1, 0, 0, numComp, 2, 2)] = 
                                              points[0][comp]*(mx2)*(my1)*(mz1)
                                            + points[1][comp]*x2*(my1)*(mz1) 
                                            + points[2][comp]*(mx2)*y1*(mz1)
                                            + points[3][comp]*x2*y1*(mz1)
                                            + points[4][comp]*(mx2)*(my1)*z1
                                            + points[5][comp]*x2*(my1)*z1 
                                            + points[6][comp]*(mx2)*y1*z1
                                            + points[7][comp]*x2*y1*z1;
                                out[INDEX4(comp, 0, 1, 0, numComp, 2, 2)] = 
                                              points[0][comp]*(mx1)*(my2)*(mz1)
                                            + points[1][comp]*x1*(my2)*(mz1) 
                                            + points[2][comp]*(mx1)*y2*(mz1)
                                            + points[3][comp]*x1*y2*(mz1)
                                            + points[4][comp]*(mx1)*(my2)*z1
                                            + points[5][comp]*x1*(my2)*z1 
                                            + points[6][comp]*(mx1)*y2*z1
                                            + points[7][comp]*x1*y2*z1;
                                out[INDEX4(comp, 1, 1, 0, numComp, 2, 2)] = 
                                              points[0][comp]*(mx2)*(my2)*(mz1)
                                            + points[1][comp]*x2*(my2)*(mz1) 
                                            + points[2][comp]*(mx2)*y2*(mz1)
                                            + points[3][comp]*x2*y2*(mz1)
                                            + points[4][comp]*(mx2)*(my2)*z1
                                            + points[5][comp]*x2*(my2)*z1 
                                            + points[6][comp]*(mx2)*y2*z1
                                            + points[7][comp]*x2*y2*z1;
                                out[INDEX4(comp, 0, 0, 1, numComp, 2, 2)] = 
                                              points[0][comp]*(mx1)*(my1)*(mz2)
                                            + points[1][comp]*x1*(my1)*(mz2) 
                                            + points[2][comp]*(mx1)*y1*(mz2)
                                            + points[3][comp]*x1*y1*(mz2)
                                            + points[4][comp]*(mx1)*(my1)*z2
                                            + points[5][comp]*x1*(my1)*z2 
                                            + points[6][comp]*(mx1)*y1*z2
                                            + points[7][comp]*x1*y1*z2;
                                out[INDEX4(comp, 1, 0, 1, numComp, 2, 2)] = 
                                              points[0][comp]*(mx2)*(my1)*(mz2)
                                            + points[1][comp]*x2*(my1)*(mz2) 
                                            + points[2][comp]*(mx2)*y1*(mz2)
                                            + points[3][comp]*x2*y1*(mz2)
                                            + points[4][comp]*(mx2)*(my1)*z2
                                            + points[5][comp]*x2*(my1)*z2 
                                            + points[6][comp]*(mx2)*y1*z2
                                            + points[7][comp]*x2*y1*z2;
                                out[INDEX4(comp, 0, 1, 1, numComp, 2, 2)] = 
                                              points[0][comp]*(mx1)*(my2)*(mz2)
                                            + points[1][comp]*x1*(my2)*(mz2) 
                                            + points[2][comp]*(mx1)*y2*(mz2)
                                            + points[3][comp]*x1*y2*(mz2)
                                            + points[4][comp]*(mx1)*(my2)*z2
                                            + points[5][comp]*x1*(my2)*z2 
                                            + points[6][comp]*(mx1)*y2*z2
                                            + points[7][comp]*x1*y2*z2;
                                out[INDEX4(comp, 1, 1, 1, numComp, 2, 2)] = 
                                              points[0][comp]*(mx2)*(my2)*(mz2)
                                            + points[1][comp]*x2*(my2)*(mz2) 
                                            + points[2][comp]*(mx2)*y2*(mz2)
                                            + points[3][comp]*x2*y2*(mz2)
                                            + points[4][comp]*(mx2)*(my2)*z2
                                            + points[5][comp]*x2*(my2)*z2 
                                            + points[6][comp]*(mx2)*y2*z2
                                            + points[7][comp]*x2*y2*z2;
                            }
                        }
                    }
                }
            }
        }
    }
}

void MultiBrick::interpolateElementsToElementsCoarser(const escript::Data& source,
        escript::Data& target, const MultiBrick& other) const
{
    if (source.isComplex()!=target.isComplex())
    {
        throw RipleyException("Programmer Error: in and out parameters do not have the same complexity.");        
    }
    if (source.isComplex())
    {
        interpolateElementsToElementsCoarserWorker(source, target, other, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateElementsToElementsCoarserWorker(source, target, other, escript::DataTypes::real_t(0));      
    }
}   

template <typename S>
void MultiBrick::interpolateElementsToElementsCoarserWorker(const escript::Data& source,
        escript::Data& target, const MultiBrick& other, S sentinel) const
{
    const int scaling = m_subdivisions/other.getNumSubdivisionsPerElement();
    const double scaling_volume = (1./scaling)*(1./scaling)*(1./scaling);
    const dim_t *theirNE = other.getNumElementsPerDim();
    const dim_t numComp = source.getDataPointSize();

    vector<double> points(scaling*2, 0);
    vector<double> first_lagrange(scaling*2, 1);
    vector<double> second_lagrange(scaling*2, 1);
    
    for (int i = 0; i < scaling*2; i+=2) {
        points[i] = (i/2 + FIRST_QUAD)/scaling;
        points[i+1] = (i/2 + SECOND_QUAD)/scaling;
    }
    
    for (int i = 0; i < scaling*2; i++) {
        first_lagrange[i] = (points[i] - SECOND_QUAD) / (FIRST_QUAD - SECOND_QUAD);
        second_lagrange[i] = (points[i] - FIRST_QUAD) / (SECOND_QUAD - FIRST_QUAD);
    }
    target.requireWrite();
    //for each of theirs
#pragma omp parallel for
    for (dim_t tz = 0; tz < theirNE[2]; tz++) {
        for (dim_t ty = 0; ty < theirNE[1]; ty++) {
            for (dim_t tx = 0; tx < theirNE[0]; tx++) {
                S *out = target.getSampleDataRW(INDEX3(tx, ty, tz, theirNE[0], theirNE[1]), sentinel);
                //for each subelement
                for (dim_t sz = 0; sz < scaling; sz++) {
                    const dim_t ez = tz*scaling + sz;
                    for (dim_t sy = 0; sy < scaling; sy++) {
                        const dim_t ey = ty*scaling + sy;
                        for (dim_t sx = 0; sx < scaling; sx++) {
                            const dim_t ex = tx*scaling + sx;
                            const S *in = source.getSampleDataRO(INDEX3(ex, ey, ez, m_NE[0], m_NE[1]), sentinel);
                            for (int quad = 0; quad < 8; quad++) {
                                int lx = sx*2 + quad%2;
                                int ly = sy*2 + (quad%4)/2;
                                int lz = sz*2 + quad/4;
                                for (dim_t comp = 0; comp < numComp; comp++) {
                                    const S quadvalue = scaling_volume * in[comp + quad*numComp];
                                    out[INDEX4(comp, 0, 0, 0, numComp, 2, 2)] += quadvalue * first_lagrange[lx] * first_lagrange[ly] * first_lagrange[lz];
                                    out[INDEX4(comp, 1, 0, 0, numComp, 2, 2)] += quadvalue * second_lagrange[lx] * first_lagrange[ly] * first_lagrange[lz];
                                    out[INDEX4(comp, 0, 1, 0, numComp, 2, 2)] += quadvalue * first_lagrange[lx] * second_lagrange[ly] * first_lagrange[lz];
                                    out[INDEX4(comp, 1, 1, 0, numComp, 2, 2)] += quadvalue * second_lagrange[lx] * second_lagrange[ly] * first_lagrange[lz];
                                    out[INDEX4(comp, 0, 0, 1, numComp, 2, 2)] += quadvalue * first_lagrange[lx] * first_lagrange[ly] * second_lagrange[lz];
                                    out[INDEX4(comp, 1, 0, 1, numComp, 2, 2)] += quadvalue * second_lagrange[lx] * first_lagrange[ly] * second_lagrange[lz];
                                    out[INDEX4(comp, 0, 1, 1, numComp, 2, 2)] += quadvalue * first_lagrange[lx] * second_lagrange[ly] * second_lagrange[lz];
                                    out[INDEX4(comp, 1, 1, 1, numComp, 2, 2)] += quadvalue * second_lagrange[lx] * second_lagrange[ly] * second_lagrange[lz];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void MultiBrick::interpolateElementsToElementsFiner(const escript::Data& source,
        escript::Data& target, const MultiBrick& other) const
{
    if (source.isComplex()!=target.isComplex())
    {
        throw RipleyException("Programmer Error: in and out parameters do not have the same complexity.");        
    }
    if (source.isComplex())
    {
        interpolateElementsToElementsFinerWorker(source, target, other, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateElementsToElementsFinerWorker(source, target, other, escript::DataTypes::real_t(0));      
    }
}  

template <typename S>
void MultiBrick::interpolateElementsToElementsFinerWorker(const escript::Data& source,
        escript::Data& target, const MultiBrick& other, S sentinel) const
{
    const int scaling = other.getNumSubdivisionsPerElement()/m_subdivisions;
    const dim_t numComp = source.getDataPointSize();

    vector<double> points(scaling*2, 0);
    vector<double> lagranges(scaling*4, 1);

    for (int i = 0; i < scaling*2; i+=2) {
        points[i] = (i/2 + FIRST_QUAD)/scaling;
        points[i+1] = (i/2 + SECOND_QUAD)/scaling;
    }
    for (int i = 0; i < scaling*2; i++) {
        lagranges[i] = (points[i] - SECOND_QUAD) / (FIRST_QUAD - SECOND_QUAD);
        lagranges[i + 2*scaling] = (points[i] - FIRST_QUAD) / (SECOND_QUAD - FIRST_QUAD);
    }
    target.requireWrite();
    //for each of ours
#pragma omp parallel for
    for (dim_t ez = 0; ez < m_NE[2]; ez++) {
        for (dim_t ey = 0; ey < m_NE[1]; ey++) {
            for (dim_t ex = 0; ex < m_NE[0]; ex++) {
                const S *in = source.getSampleDataRO(INDEX3(ex, ey, ez, m_NE[0], m_NE[1]), sentinel);
                //for each subelement
                for (dim_t sz = 0; sz < scaling; sz++) {
                    const dim_t tz = ez*scaling + sz;
                    for (dim_t sy = 0; sy < scaling; sy++) {
                        const dim_t ty = ey*scaling + sy;
                        for (dim_t sx = 0; sx < scaling; sx++) {
                            const dim_t tx = ex*scaling + sx;
                            S *out = target.getSampleDataRW(INDEX3(tx, ty, tz, m_NE[0]*scaling, m_NE[1]*scaling), sentinel);
                            for (int quad = 0; quad < 8; quad++) {
                                const int lx = scaling*2*(quad%2) + sx*2;
                                const int ly = scaling*2*((quad%4)/2) + sy*2;
                                const int lz = scaling*2*(quad/4) + sz*2;
                                for (dim_t comp = 0; comp < numComp; comp++) {
                                    const S quadvalue = in[comp + quad*numComp];
                                    out[INDEX4(comp, 0, 0, 0, numComp, 2, 2)] += quadvalue * lagranges[lx] * lagranges[ly] * lagranges[lz];
                                    out[INDEX4(comp, 0, 1, 0, numComp, 2, 2)] += quadvalue * lagranges[lx] * lagranges[ly+1] * lagranges[lz];
                                    out[INDEX4(comp, 1, 0, 0, numComp, 2, 2)] += quadvalue * lagranges[lx+1] * lagranges[ly] * lagranges[lz];
                                    out[INDEX4(comp, 1, 1, 0, numComp, 2, 2)] += quadvalue * lagranges[lx+1] * lagranges[ly+1] * lagranges[lz];
                                    out[INDEX4(comp, 0, 0, 1, numComp, 2, 2)] += quadvalue * lagranges[lx] * lagranges[ly] * lagranges[lz+1];
                                    out[INDEX4(comp, 0, 1, 1, numComp, 2, 2)] += quadvalue * lagranges[lx] * lagranges[ly+1] * lagranges[lz+1];
                                    out[INDEX4(comp, 1, 0, 1, numComp, 2, 2)] += quadvalue * lagranges[lx+1] * lagranges[ly] * lagranges[lz+1];
                                    out[INDEX4(comp, 1, 1, 1, numComp, 2, 2)] += quadvalue * lagranges[lx+1] * lagranges[ly+1] * lagranges[lz+1];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void MultiBrick::interpolateAcross(escript::Data& target,
                                     const escript::Data& source) const
{
    const MultiBrick *other =
                dynamic_cast<const MultiBrick *>(target.getDomain().get());
    if (other == NULL)
        throw RipleyException("Invalid interpolation: Domains must both be instances of MultiBrick");
    //shouldn't ever happen, but I want to know if it does
    if (other == this)
        throw RipleyException("interpolateAcross: this domain is the target");
        
    validateInterpolationAcross(source.getFunctionSpace().getTypeCode(),
            *(target.getDomain().get()), target.getFunctionSpace().getTypeCode());
    int fsSource = source.getFunctionSpace().getTypeCode();
    int fsTarget = target.getFunctionSpace().getTypeCode();

    std::stringstream msg;
    msg << "Invalid interpolation: interpolation not implemented for function space "
        << functionSpaceTypeAsString(fsSource)
        << " -> "
        << functionSpaceTypeAsString(fsTarget);
    if (other->getNumSubdivisionsPerElement() > getNumSubdivisionsPerElement()) {
        switch (fsSource) {
            case Nodes:
                switch (fsTarget) {
                    case Nodes:
                    case ReducedNodes:
                    case DegreesOfFreedom:
                    case ReducedDegreesOfFreedom:
                        interpolateNodesToNodesFiner(source, target, *other);
                        return;
                    case Elements:
                        interpolateNodesToElementsFiner(source, target, *other);
                        return;
                }
                break;
            case Elements:
                switch (fsTarget) {
                    case Elements:
                        interpolateElementsToElementsFiner(source, target, *other);
                        return;
                }
                break;
            case ReducedElements:
                switch (fsTarget) {
                    case Elements:
                        interpolateReducedToElementsFiner(source, target, *other);
                        return;
                }
                break;
        }
        msg << " when target is a finer mesh";
    } else {
        switch (fsSource) {
            case Nodes:
                switch (fsTarget) {
                    case Elements:
                        escript::Data elements=escript::Vector(0., escript::function(*this), true);
                        interpolateNodesOnElements(elements, source, false);
                        interpolateElementsToElementsCoarser(elements, target, *other);
                        return;
                }
                break;
            case Elements:
                switch (fsTarget) {
                    case Elements:
                        interpolateElementsToElementsCoarser(source, target, *other);
                        return;
                }
                break;
        }
        msg << " when target is a coarser mesh";
    }
    throw RipleyException(msg.str());
}

string MultiBrick::getDescription() const
{
    return "ripley::MultiBrick";
}

bool MultiBrick::operator==(const AbstractDomain& other) const
{
    const MultiBrick* o=dynamic_cast<const MultiBrick*>(&other);
    if (o) {
        return (RipleyDomain::operator==(other) &&
                m_gNE[0]==o->m_gNE[0] && m_gNE[1]==o->m_gNE[1] && m_gNE[2]==o->m_gNE[2]
                && m_origin[0]==o->m_origin[0] && m_origin[1]==o->m_origin[1] && m_origin[2]==o->m_origin[2]
                && m_length[0]==o->m_length[0] && m_length[1]==o->m_length[1] && m_length[2]==o->m_length[2]
                && m_NX[0]==o->m_NX[0] && m_NX[1]==o->m_NX[1] && m_NX[2]==o->m_NX[2]
                && m_subdivisions == o->m_subdivisions);
    }

    return false;
}

void MultiBrick::readNcGrid(escript::Data& out, string filename, string varname,
            const ReaderParameters& params) const
{
    if (m_subdivisions != 1)
        throw RipleyException("Non-parent MultiBricks cannot read datafiles");
    Brick::readNcGrid(out, filename, varname, params);
}

void MultiBrick::readBinaryGrid(escript::Data& out, string filename,
                               const ReaderParameters& params) const
{
    if (m_subdivisions != 1)
        throw RipleyException("Non-parent MultiBricks cannot read datafiles");
    Brick::readBinaryGrid(out, filename, params);
}

void MultiBrick::readBinaryGridFromZipped(escript::Data& out, string filename,
                               const ReaderParameters& params) const
{
    if (m_subdivisions != 1)
        throw RipleyException("Non-parent MultiBricks cannot read datafiles");
    Brick::readBinaryGridFromZipped(out, filename, params);
}

void MultiBrick::writeBinaryGrid(const escript::Data& in, string filename,
                                int byteOrder, int dataType) const
{
    if (m_subdivisions != 1)
        throw RipleyException("Non-parent MultiBricks cannot read datafiles");
    Brick::writeBinaryGrid(in, filename, byteOrder, dataType);
}

void MultiBrick::dump(const string& fileName) const
{
    if (m_subdivisions != 1)
        throw RipleyException("Non-parent MultiBricks dump not implemented");
    Brick::dump(fileName);
}

const dim_t* MultiBrick::borrowSampleReferenceIDs(int fsType) const
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
    throw RipleyException(msg.str());
}

bool MultiBrick::ownSample(int fsType, index_t id) const
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
    throw RipleyException(msg.str());
}

void MultiBrick::setToNormal(escript::Data& out) const
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
        throw RipleyException(msg.str());
    }
}

void MultiBrick::setToSize(escript::Data& out) const
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
        throw RipleyException(msg.str());
    }
}

void MultiBrick::Print_Mesh_Info(const bool full) const
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
IndexVector MultiBrick::getDiagonalIndices(bool upperOnly) const
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

//private
void MultiBrick::populateSampleIds()
{
    // degrees of freedom are numbered from left to right, bottom to top, front
    // to back in each rank, continuing on the next rank (ranks also go
    // left-right, bottom-top, front-back).
    // This means rank 0 has id 0...n0-1, rank 1 has id n0...n1-1 etc. which
    // helps when writing out data rank after rank.

    // build node distribution vector first.
    // rank i owns m_nodeDistribution[i+1]-nodeDistribution[i] nodes.
    // Unlike regular ripley domains this is NOT constant for all ranks so
    // we do an Allgather (we could have also computed per rank but it's a bit
    // involved)
    m_nodeDistribution.assign(m_mpiInfo->size+1, 0);
    dim_t numDOF=getNumDOF();
    if (m_mpiInfo->size > 1) {
#if ESYS_MPI
        MPI_Allgather(&numDOF, 1, MPI_DIM_T, &m_nodeDistribution[0], 1,
                      MPI_DIM_T, m_mpiInfo->comm);

        // accumulate
        dim_t accu = 0;
        for (int rank=0; rank<m_mpiInfo->size; rank++) {
            const dim_t n = m_nodeDistribution[rank];
            m_nodeDistribution[rank] = accu;
            accu += n;
        }
        ESYS_ASSERT(accu == getNumDataPointsGlobal(),
                "something went wrong computing the DOF distribution!");

        m_nodeDistribution[m_mpiInfo->size] = accu;
#endif
    } else {
        m_nodeDistribution[m_mpiInfo->size] = numDOF;
    }

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

//private
vector<IndexVector> MultiBrick::getConnections(bool includeShared) const
{
    // returns a vector v of size numDOF where v[i] is a vector with indices
    // of DOFs connected to i (up to 27 in 3D)
    const dim_t nDOF0 = (m_gNE[0]+1)/m_NX[0];
    const dim_t nDOF1 = (m_gNE[1]+1)/m_NX[1];
    const dim_t nDOF2 = (m_gNE[2]+1)/m_NX[2];
    const dim_t M = nDOF0*nDOF1*nDOF2;
    vector<IndexVector> indices(M);

#pragma omp parallel for
    for (index_t i=0; i < M; i++) {
        const index_t x = i % nDOF0;
        const index_t y = i % (nDOF0*nDOF1)/nDOF0;
        const index_t z = i / (nDOF0*nDOF1);
        // loop through potential neighbours and add to index if positions are
        // within bounds
        for (int i2=z-1; i2<z+2; i2++) {
            for (int i1=y-1; i1<y+2; i1++) {
                for (int i0=x-1; i0<x+2; i0++) {
                    if (i0>=0 && i1>=0 && i2>=0
                            && i0<nDOF0 && i1<nDOF1 && i2<nDOF2) {
                        indices[i].push_back(i2*nDOF0*nDOF1 + i1*nDOF0 + i0);
                    }
                }
            }
        }
    }
    return indices;
}

RankVector MultiBrick::getOwnerVector(int fsType) const
{
    if (m_subdivisions != 1)
        throw RipleyException("Multiresolution domains only support ownership for the coarsest level");
    return Brick::getOwnerVector(fsType);
}

dim_t MultiBrick::findNode(const double *coords) const
{
    const dim_t NOT_MINE = -1;
    //is the found element even owned by this rank
    // (inside owned or shared elements but will map to an owned element)
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
                         "node, implementation problem in MultiBrick::findNode()");
    }
    return closest;
}

} // end of namespace ripley

