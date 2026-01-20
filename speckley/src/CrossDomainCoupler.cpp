
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <speckley/CrossDomainCoupler.h>
#include <speckley/lagrange_functions.h>

#include <escript/index.h>

#define MINE 1
#define SHARED 0
#define THEIRS -1

namespace speckley {

//this should probably come from ripley dynamically
const double ripleyLocations[2] = {.21132486540518711775, .78867513459481288225};

bool probeInterpolationAcross(int fsType_source,
        const escript::AbstractDomain& domain, int fsType_target, int dim)
{
    try {
        const ripley::RipleyDomain& other =
                dynamic_cast<const ripley::RipleyDomain &>(domain);
        if (other.getDim() != dim)
            return false;
    } catch (const std::bad_cast& e) {
        return false;
    }
    return (fsType_source == Elements && fsType_target == ripley::Elements);
}

RipleyCoupler::RipleyCoupler(const speckley::SpeckleyDomain *speck,
        const double s_dx[2], int rank) : speck(speck)
#ifdef ESYS_MPI
, rank(rank)
#endif
{
    const int *splits = speck->getNumSubdivisionsPerDim();
    const dim_t *elements = speck->getNumElementsPerDim();
    const dim_t *edges = speck->getNumFacesPerBoundary();//left,right,bottom,top
    for (int i = 0; i < speck->getDim(); i++) {
        this->s_dx[i] = s_dx[i];
        s_NX[i] = splits[i];
        s_NE[i] = elements[i];
        speckley_origin[i] = speck->getLocalCoordinate(0,i);
        hasLower[i] = (edges[2*i] == 0);
        hasUpper[i] = (edges[2*i + 1] == 0);
    }
    if (speck->getDim() == 2) {
        hasLower[2] = false;
        hasUpper[2] = false;
        s_NX[2] = 1;
    }
    order = speck->getOrder();
    numQuads = order + 1;
#ifdef ESYS_MPI
    comm = speck->getMPIComm();
#endif
}

bool RipleyCoupler::validInterpolation(escript::Data& target,
            const escript::Data& source,
            const SpeckleyDomain *speck, const double *s_dx,
            const ripley::RipleyDomain *other) const
{
    if (source.getDomain().get() != speck)
        throw SpeckleyException(
            "ripleyCoupler: interpolation from unsupported domain");

    if (speck->getDim() != other->getDim())
        throw SpeckleyException(
            "ripleyCoupler: domains must have the same dimensions");

    //validate functionspaces
    const int tFS = target.getFunctionSpace().getTypeCode();
    const int sFS = source.getFunctionSpace().getTypeCode();
    if (sFS != Elements)
        throw SpeckleyException(
                "ripleyCoupler: source data must be in Function functionspace");
    if (tFS != ripley::Elements)
        throw SpeckleyException(
                "ripleyCoupler: target data must be in Function functionspace");

    //ensure same division setup
    const int *r_NX = other->getNumSubdivisionsPerDim();
    for (int i = 0; i < speck->getDim(); i++) {
        if (r_NX[i] != s_NX[i]) {
            throw SpeckleyException(
                    "ripleyCoupler: domain subdivisions don't match");
        }
    }

    //ensure same data shape
    if (target.getDataPointSize() != source.getDataPointSize())
        throw SpeckleyException(
                "ripleyCoupler: data point size mismatch");

    const double *r_len = other->getLength();
    const double *s_len = speck->getLength();
    for (int i = 0; i < speck->getDim(); i++) {
        if (r_len[i] != s_len[i]) {
            throw SpeckleyException("ripleyCoupler: domain length mismatch");
        }
    }

#ifdef ESYS_MPI
    int res;
    if (MPI_Comm_compare(speck->getMPIComm(), other->getMPIComm(), &res)
            != MPI_SUCCESS)
        throw SpeckleyException(
                "ripleyCoupler: domains have bad communicators");
    if (res != MPI_IDENT)
        throw SpeckleyException(
                "ripleyCoupler: domain communicators are not identical");
#endif
    return true;
}

void RipleyCoupler::generateLocations(struct Ripley& r, double *positions[3]) const
{
    switch(order) {
        case 2:
            for (int i = 0; i < speck->getDim(); i++) {
                double *axis = positions[i];
                double first = ripleyLocations[0]*r.dx[i];
                double second = ripleyLocations[1]*r.dx[i];
        #pragma omp parallel for
                for (dim_t e = 0; e < r.NE[i]; e++) {
                    const double x = r.domain->getLocalCoordinate(e,i);
                    calculateOrder2(i, x + first, axis + 2*e*numQuads);
                    calculateOrder2(i, x + second, axis + (2*e+1)*numQuads);
                }
            }
            break;
        case 3:
            for (int i = 0; i < speck->getDim(); i++) {
                double *axis = positions[i];
                double first = ripleyLocations[0]*r.dx[i];
                double second = ripleyLocations[1]*r.dx[i];
        #pragma omp parallel for
                for (dim_t e = 0; e < r.NE[i]; e++) {
                    const double x = r.domain->getLocalCoordinate(e,i);
                    calculateOrder3(i, x + first, axis + 2*e*numQuads);
                    calculateOrder3(i, x + second, axis + (2*e+1)*numQuads);
                }
            }
            break;
        case 4:
            for (int i = 0; i < speck->getDim(); i++) {
                double *axis = positions[i];
                double first = ripleyLocations[0]*r.dx[i];
                double second = ripleyLocations[1]*r.dx[i];
        #pragma omp parallel for
                for (dim_t e = 0; e < r.NE[i]; e++) {
                    const double x = r.domain->getLocalCoordinate(e,i);
                    calculateOrder4(i, x + first, axis + 2*e*numQuads);
                    calculateOrder4(i, x + second, axis + (2*e+1)*numQuads);
                }
            }
            break;
        case 5:
            for (int i = 0; i < speck->getDim(); i++) {
                double *axis = positions[i];
                double first = ripleyLocations[0]*r.dx[i];
                double second = ripleyLocations[1]*r.dx[i];
        #pragma omp parallel for
                for (dim_t e = 0; e < r.NE[i]; e++) {
                    const double x = r.domain->getLocalCoordinate(e,i);
                    calculateOrder5(i, x + first, axis + 2*e*numQuads);
                    calculateOrder5(i, x + second, axis + (2*e+1)*numQuads);
                }
            }
            break;
        case 6:
            for (int i = 0; i < speck->getDim(); i++) {
                double *axis = positions[i];
                double first = ripleyLocations[0]*r.dx[i];
                double second = ripleyLocations[1]*r.dx[i];
        #pragma omp parallel for
                for (dim_t e = 0; e < r.NE[i]; e++) {
                    const double x = r.domain->getLocalCoordinate(e,i);
                    calculateOrder6(i, x + first, axis + 2*e*numQuads);
                    calculateOrder6(i, x + second, axis + (2*e+1)*numQuads);
                }
            }
            break;
        case 7:
            for (int i = 0; i < speck->getDim(); i++) {
                double *axis = positions[i];
                double first = ripleyLocations[0]*r.dx[i];
                double second = ripleyLocations[1]*r.dx[i];
        #pragma omp parallel for
                for (dim_t e = 0; e < r.NE[i]; e++) {
                    const double x = r.domain->getLocalCoordinate(e,i);
                    calculateOrder7(i, x + first, axis + 2*e*numQuads);
                    calculateOrder7(i, x + second, axis + (2*e+1)*numQuads);
                }
            }
            break;
        case 8:
            for (int i = 0; i < speck->getDim(); i++) {
                double *axis = positions[i];
                double first = ripleyLocations[0]*r.dx[i];
                double second = ripleyLocations[1]*r.dx[i];
        #pragma omp parallel for
                for (dim_t e = 0; e < r.NE[i]; e++) {
                    const double x = r.domain->getLocalCoordinate(e,i);
                    calculateOrder8(i, x + first, axis + 2*e*numQuads);
                    calculateOrder8(i, x + second, axis + (2*e+1)*numQuads);
                }
            }
            break;
        case 9:
            for (int i = 0; i < speck->getDim(); i++) {
                double *axis = positions[i];
                double first = ripleyLocations[0]*r.dx[i];
                double second = ripleyLocations[1]*r.dx[i];
        #pragma omp parallel for
                for (dim_t e = 0; e < r.NE[i]; e++) {
                    const double x = r.domain->getLocalCoordinate(e,i);
                    calculateOrder9(i, x + first, axis + 2*e*numQuads);
                    calculateOrder9(i, x + second, axis + (2*e+1)*numQuads);
                }
            }
            break;
        case 10:
            for (int i = 0; i < speck->getDim(); i++) {
                double *axis = positions[i];
                double first = ripleyLocations[0]*r.dx[i];
                double second = ripleyLocations[1]*r.dx[i];
        #pragma omp parallel for
                for (dim_t e = 0; e < r.NE[i]; e++) {
                    const double x = r.domain->getLocalCoordinate(e,i);
                    calculateOrder10(i, x + first, axis + 2*e*numQuads);
                    calculateOrder10(i, x + second, axis + (2*e+1)*numQuads);
                }
            }
            break;
        default:
            throw SpeckleyException(
                    "RipleyCoupler:: unexpected order of domain");
    }
}

void RipleyCoupler::calculateOrder2(int dim, double loc, double *results) const
{
    loc -= speckley_origin[dim];
    dim_t source_e = loc / s_dx[dim];
    double local = ((loc - source_e * s_dx[dim]) / s_dx[dim]) * 2 - 1;
    results[0] = lagrange_degree2_0(local);
    results[1] = lagrange_degree2_1(local);
    results[2] = lagrange_degree2_2(local);
}

void RipleyCoupler::calculateOrder3(int dim, double loc, double *results) const
{
    loc -= speckley_origin[dim];
    dim_t source_e = loc / s_dx[dim];
    double local = ((loc - source_e * s_dx[dim]) / s_dx[dim]) * 2 - 1;
    results[0] = lagrange_degree3_0(local);
    results[1] = lagrange_degree3_1(local);
    results[2] = lagrange_degree3_2(local);
    results[3] = lagrange_degree3_3(local);
}

void RipleyCoupler::calculateOrder4(int dim, double loc, double *results) const
{
    loc -= speckley_origin[dim];
    dim_t source_e = loc / s_dx[dim];
    double local = ((loc - source_e * s_dx[dim]) / s_dx[dim]) * 2 - 1;
    results[0] = lagrange_degree4_0(local);
    results[1] = lagrange_degree4_1(local);
    results[2] = lagrange_degree4_2(local);
    results[3] = lagrange_degree4_3(local);
    results[4] = lagrange_degree4_4(local);
}

void RipleyCoupler::calculateOrder5(int dim, double loc, double *results) const
{
    loc -= speckley_origin[dim];
    dim_t source_e = loc / s_dx[dim];
    double local = ((loc - source_e * s_dx[dim]) / s_dx[dim]) * 2 - 1;
    results[0] = lagrange_degree5_0(local);
    results[1] = lagrange_degree5_1(local);
    results[2] = lagrange_degree5_2(local);
    results[3] = lagrange_degree5_3(local);
    results[4] = lagrange_degree5_4(local);
    results[5] = lagrange_degree5_5(local);
}

void RipleyCoupler::calculateOrder6(int dim, double loc, double *results) const
{
    loc -= speckley_origin[dim];
    dim_t source_e = loc / s_dx[dim];
    double local = ((loc - source_e * s_dx[dim]) / s_dx[dim]) * 2 - 1;
    results[0] = lagrange_degree6_0(local);
    results[1] = lagrange_degree6_1(local);
    results[2] = lagrange_degree6_2(local);
    results[3] = lagrange_degree6_3(local);
    results[4] = lagrange_degree6_4(local);
    results[5] = lagrange_degree6_5(local);
    results[6] = lagrange_degree6_6(local);
}

void RipleyCoupler::calculateOrder7(int dim, double loc, double *results) const
{
    loc -= speckley_origin[dim];
    dim_t source_e = loc / s_dx[dim];
    double local = ((loc - source_e * s_dx[dim]) / s_dx[dim]) * 2 - 1;
    results[0] = lagrange_degree7_0(local);
    results[1] = lagrange_degree7_1(local);
    results[2] = lagrange_degree7_2(local);
    results[3] = lagrange_degree7_3(local);
    results[4] = lagrange_degree7_4(local);
    results[5] = lagrange_degree7_5(local);
    results[6] = lagrange_degree7_6(local);
    results[7] = lagrange_degree7_7(local);
}

void RipleyCoupler::calculateOrder8(int dim, double loc, double *results) const
{
    loc -= speckley_origin[dim];
    dim_t source_e = loc / s_dx[dim];
    double local = ((loc - source_e * s_dx[dim]) / s_dx[dim]) * 2 - 1;
    results[0] = lagrange_degree8_0(local);
    results[1] = lagrange_degree8_1(local);
    results[2] = lagrange_degree8_2(local);
    results[3] = lagrange_degree8_3(local);
    results[4] = lagrange_degree8_4(local);
    results[5] = lagrange_degree8_5(local);
    results[6] = lagrange_degree8_6(local);
    results[7] = lagrange_degree8_7(local);
    results[8] = lagrange_degree8_8(local);
}

void RipleyCoupler::calculateOrder9(int dim, double loc, double *results) const
{
    loc -= speckley_origin[dim];
    dim_t source_e = loc / s_dx[dim];
    double local = ((loc - source_e * s_dx[dim]) / s_dx[dim]) * 2 - 1;
    results[0] = lagrange_degree9_0(local);
    results[1] = lagrange_degree9_1(local);
    results[2] = lagrange_degree9_2(local);
    results[3] = lagrange_degree9_3(local);
    results[4] = lagrange_degree9_4(local);
    results[5] = lagrange_degree9_5(local);
    results[6] = lagrange_degree9_6(local);
    results[7] = lagrange_degree9_7(local);
    results[8] = lagrange_degree9_8(local);
    results[9] = lagrange_degree9_9(local);
}

void RipleyCoupler::calculateOrder10(int dim, double loc, double *results) const
{
    loc -= speckley_origin[dim];
    dim_t source_e = loc / s_dx[dim];
    double local = ((loc - source_e * s_dx[dim]) / s_dx[dim]) * 2 - 1;
    results[0] = lagrange_degree10_0(local);
    results[1] = lagrange_degree10_1(local);
    results[2] = lagrange_degree10_2(local);
    results[3] = lagrange_degree10_3(local);
    results[4] = lagrange_degree10_4(local);
    results[5] = lagrange_degree10_5(local);
    results[6] = lagrange_degree10_6(local);
    results[7] = lagrange_degree10_7(local);
    results[8] = lagrange_degree10_8(local);
    results[9] = lagrange_degree10_9(local);
    results[10] = lagrange_degree10_10(local);
}

void RipleyCoupler::getEdgeSpacing(struct Ripley r, int *lower, int *upper) const
{
    for (int i = 0; i < speck->getDim(); i++) {
        const double first = ripleyLocations[0]*r.dx[i];
        const double second = ripleyLocations[1]*r.dx[i];
        double start = r.domain->getLocalCoordinate(0, i) - speckley_origin[i];

        lower[i] = SHARED;

        if (start + first > 0)
            lower[i] = MINE;
        else if (start + second < 0)
            lower[i] = THEIRS;

        start = r.domain->getLocalCoordinate(r.NE[i]-1, i) - speckley_origin[i];
        upper[i] = 0;
        if ((start + first)/s_dx[i] >= s_NE[i])
            upper[i] = THEIRS;
        else if ((start + second)/s_dx[i] < s_NE[i])
            upper[i] = MINE;
    }
}


inline void RipleyCoupler::calculate(struct Ripley& r,
        dim_t ex, dim_t ey, dim_t ez, int oqx, int oqy, int oqz, double *out,
        const double *factor_x, const double *factor_y, const double *factor_z,
        const escript::Data& source) const
{
    const double x = r.domain->getLocalCoordinate(ex, 0) - speckley_origin[0]
                + ripleyLocations[oqx]*r.dx[0];
    const double y = r.domain->getLocalCoordinate(ey, 1) - speckley_origin[1]
                + ripleyLocations[oqy]*r.dx[1];
    const dim_t source_ex = x / s_dx[0];
    const dim_t source_ey = y / s_dx[1];

    if (speck->getDim() == 2) {
        const double *sdata = source.getSampleDataRO(
                                    INDEX2(source_ex, source_ey, s_NE[0]));
        for (int sqy = 0; sqy < numQuads; sqy++) {
            for (int sqx = 0; sqx < numQuads; sqx++) {
                for (int comp = 0; comp < numComp; comp++) {
                    out[INDEX3(comp,oqx,oqy,numComp,2)]
                          += sdata[INDEX3(comp,sqx,sqy,numComp,numQuads)]
                             * factor_x[sqx] * factor_y[sqy];
                }
            }
        }
    } else { //dim == 3
        const double z = r.domain->getLocalCoordinate(ez, 2) - speckley_origin[2]
                + ripleyLocations[oqz]*r.dx[2];
        const dim_t source_ez = z / s_dx[2];
        const double *sdata = source.getSampleDataRO(
                    INDEX3(source_ex, source_ey, source_ez, s_NE[0], s_NE[1]));
        for (int sqz = 0; sqz < numQuads; sqz++) {
            for (int sqy = 0; sqy < numQuads; sqy++) {
                for (int sqx = 0; sqx < numQuads; sqx++) {
                    for (int comp = 0; comp < numComp; comp++) {
                        out[INDEX4(comp,oqx,oqy,oqz,numComp,2,2)]
                              += sdata[INDEX4(comp,sqx,sqy,sqz,numComp,numQuads,numQuads)]
                                * factor_x[sqx] * factor_y[sqy] * factor_z[sqz];
                    }
                }
            }
        }
    }
}

void RipleyCoupler::interpolate(escript::Data& target,
        const escript::Data& source) const
{
    //can only interpolate to ripley right now, so check it
    const ripley::RipleyDomain *other =
                            dynamic_cast<const ripley::RipleyDomain *>
                                                (target.getDomain().get());
    if (other == NULL)
        throw SpeckleyException("interpolation to unsupported domain");
    validInterpolation(target, source, speck, s_dx, other); //throws if bad

    //gather details from the ripley side
    const dim_t *r_NE = other->getNumElementsPerDim();
    const double *r_dx = other->getElementLength();
    const dim_t *edges = other->getNumFacesPerBoundary();//left,right,bottom,top
    struct Ripley r = {other, {0,0,0},{0,0,0},{0,0,0}, {0,0,0}};
    for (int i = 0; i < speck->getDim(); i++) {
        r.NE[i] = r_NE[i];
        r.dx[i] = r_dx[i];
        r.mins[i] =           (edges[2*i] == 0 ? 1 : 0);
        r.maxs[i] = r_NE[i] - (edges[2*i + 1] == 0 ? 1 : 0);
    }

    std::vector<double> factors_x_vec(2*r.NE[0]*numQuads);
    std::vector<double> factors_y_vec(2*r.NE[1]*numQuads);
    std::vector<double> factors_z_vec(speck->getDim() == 3 ? 2*r.NE[2]*numQuads : numQuads, 1);
    double *factors[3] = {&factors_x_vec[0], &factors_y_vec[0], &factors_z_vec[0]};
    generateLocations(r, factors);
    double *factors_x = factors[0];
    double *factors_y = factors[1];
    double *factors_z = factors[2];

    numComp = source.getDataPointSize();
    target.requireWrite();

    //positional help, 0 = split, -1 = both on neighbour, 1 = both local
    int upper[3] = {1,1,1};
    int lower[3] = {1,1,1};
    getEdgeSpacing(r, lower, upper);
    for (int i = 0; i < speck->getDim(); i++) {
        if (hasUpper[i] && upper[i] == MINE) {
            r.maxs[i] += 1;
        }
        if (hasLower[i] && lower[i] == MINE) {
            r.mins[i] = 0;
        }
    }
    if (speck->getDim() == 2) {
        r.mins[2] = 0;
        r.maxs[2] = 1;
    }
    //end of setup

    const dim_t xmin = (lower[0] == SHARED) ? 1 : (lower[0] == MINE) ? 0 : 2;
    const dim_t xmax = r.NE[0]*2 - ((upper[0] == SHARED) ? 1 : (upper[0] == MINE) ? 0 : 2);
    const dim_t ymin = (lower[1] == SHARED) ? 1 : (lower[1] == MINE) ? 0 : 2;
    const dim_t ymax = r.NE[1]*2 - ((upper[1] == SHARED) ? 1 : (upper[1] == MINE) ? 0 : 2);
    const dim_t zmin = (speck->getDim() == 1) ? 1 : (lower[2] == SHARED) ? 1 : (lower[2] == MINE) ? 0 : 2;
    const dim_t zmax = r.NE[2]*2 - ((speck->getDim() == 2) ? -1 : ((upper[2] == SHARED) ? 1 : (upper[2] == MINE) ? 0 : 2));

    // inefficient cache use, but the alternative isn't simple
    for (dim_t z = zmin; z < zmax; z++) {
        const dim_t ez = z/2;
        const int oqz = z%2;
        double *precalc_z = factors_z + z * numQuads;
#pragma omp parallel for
        for (dim_t y = ymin; y < ymax; y++) {
            const dim_t ey = y/2;
            const dim_t oqy = y%2;
            double *precalc_y = factors_y + y * numQuads;
            for (dim_t x = xmin; x < xmax; x++) {
                const dim_t ex = x/2;
                const dim_t oqx = x%2;
                double *precalc_x = factors_x + x * numQuads;
                double *out = target.getSampleDataRW(INDEX3(ex,ey,ez,r.NE[0],r.NE[1]));
                calculate(r, ex, ey, ez, oqx, oqy, oqz, out, precalc_x, precalc_y, precalc_z, source);
            }
        }
    }
    //return early if we don't even need to do this next work
    if (s_NX[0] * s_NX[1] * s_NX[2] == 1) {
        return;
    }

    if (speck->getDim() == 2) {
        if (hasLower[0] || hasUpper[0])
            shareRectangleXEdges(r, hasLower[0], hasUpper[0], lower[0], upper[0], target);
        if (hasLower[1] || hasUpper[1])
            shareRectangleYEdges(r, hasLower[1], hasUpper[1], lower[1], upper[1], target);
    } else {
        if (hasLower[0] || hasUpper[0])
            shareBrickXFaces(r, hasLower[0], hasUpper[0], lower[0], upper[0], target);
        if (hasLower[1] || hasUpper[1])
            shareBrickYFaces(r, hasLower[1], hasUpper[1], lower[1], upper[1], target);
        if (hasLower[2] || hasUpper[2])
            shareBrickZFaces(r, hasLower[2], hasUpper[2], lower[2], upper[2], target);
    }
}

void RipleyCoupler::shareBrickXFaces(struct Ripley& r, int hasLower,
        int hasUpper, int lower, int upper, escript::Data& target) const
{
#ifdef ESYS_MPI
    const dim_t leftCount = r.NE[2]*r.NE[1]*numComp*4 * (1 + lower*lower);
    const dim_t rightCount = r.NE[2]*r.NE[1]*numComp*4 * (1 + upper*upper);
    std::vector<double> left(leftCount, 0);
    std::vector<double> right(rightCount, 0);
    std::vector<double> rrecv(rightCount, 0);
    std::vector<double> lrecv(leftCount, 0);
    const size_t pointsize = numComp * sizeof(double);

    //fill outbound arrays
    if (lower == SHARED) { //single
        const dim_t z_size = r.NE[1]*numComp*4;
#pragma omp parallel for
        for (dim_t ez = 0; ez < r.NE[2]; ez++) {
            const dim_t z_offset = ez * z_size;
            for (dim_t ey = 0; ey < r.NE[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX3(0, ey, ez, r.NE[0], r.NE[1]));
                memcpy(&left[z_offset + ey*4*numComp], out + numComp, pointsize);
                memcpy(&left[z_offset + (ey*4 + 1)*numComp], out + 3*numComp, pointsize);
                memcpy(&left[z_offset + (ey*4 + 2)*numComp], out + 5*numComp, pointsize);
                memcpy(&left[z_offset + (ey*4 + 3)*numComp], out + 7*numComp, pointsize);
            }
        }
    } else if (hasLower && lower == MINE) { //double
        const dim_t z_size = r.NE[1]*numComp*8;
#pragma omp parallel for
        for (dim_t ez = 0; ez < r.NE[2]; ez++) {
            const dim_t z_offset = ez * z_size;
            for (dim_t ey = 0; ey < r.NE[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX3(0, ey, ez, r.NE[0], r.NE[1]));
                memcpy(&left[z_offset + ey*8*numComp], out, pointsize*8);
            }
        }
    }
    if (upper == SHARED) { //single
        const dim_t z_size = r.NE[1]*numComp*4;
#pragma omp parallel for
        for (dim_t ez = 0; ez < r.NE[2]; ez++) {
            const dim_t z_offset = ez * z_size;
            for (dim_t ey = 0; ey < r.NE[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX3(r.NE[0]-1, ey, ez, r.NE[0], r.NE[1]));
                memcpy(&right[z_offset + ey*4*numComp], out, pointsize);
                memcpy(&right[z_offset + (ey*4 + 1)*numComp], out + 2*numComp, pointsize);
                memcpy(&right[z_offset + (ey*4 + 2)*numComp], out + 4*numComp, pointsize);
                memcpy(&right[z_offset + (ey*4 + 3)*numComp], out + 6*numComp, pointsize);
            }
        }
    } else if (hasUpper && upper == MINE) { //double
        const dim_t z_size = r.NE[1]*numComp*8;
#pragma omp parallel for
        for (dim_t ez = 0; ez < r.NE[2]; ez++) {
            const dim_t z_offset = ez * z_size;
            for (dim_t ey = 0; ey < r.NE[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX3(r.NE[0]-1, ey, ez, r.NE[0], r.NE[1]));
                memcpy(&right[z_offset + ey*8*numComp], out, pointsize*8);
            }
        }
    }

    shareWithNeighbours((rank % s_NX[0]) % 2, hasLower, hasUpper, &left[0],
            &right[0], &lrecv[0], &rrecv[0], leftCount, rightCount, 1);

    //unpacking
    if (lower == SHARED) {
        const dim_t z_size = r.NE[1]*numComp*4;
#pragma omp parallel for
        for (dim_t ez = 0; ez < r.NE[2]; ez++) {
            const dim_t z_offset = ez * z_size;
            for (dim_t ey = 0; ey < r.NE[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX3(0, ey, ez, r.NE[0], r.NE[1]));
                memcpy(out, &lrecv[z_offset + ey*4*numComp], pointsize);
                memcpy(out + 2*numComp, &lrecv[z_offset + (ey*4 + 1)*numComp], pointsize);
                memcpy(out + 4*numComp, &lrecv[z_offset + (ey*4 + 2)*numComp], pointsize);
                memcpy(out + 6*numComp, &lrecv[z_offset + (ey*4 + 3)*numComp], pointsize);
            }
        }
    } else if (lower == THEIRS) {
        const dim_t z_size = r.NE[1]*numComp*8;
#pragma omp parallel for
        for (dim_t ez = 0; ez < r.NE[2]; ez++) {
            const dim_t z_offset = ez * z_size;
            for (dim_t ey = 0; ey < r.NE[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX3(0, ey, ez, r.NE[0], r.NE[1]));
                memcpy(out, &lrecv[z_offset + ey*8*numComp], pointsize*8);
            }
        }
    }
    if (upper == SHARED) {
        const dim_t z_size = r.NE[1]*numComp*4;
#pragma omp parallel for
        for (dim_t ez = 0; ez < r.NE[2]; ez++) {
            const dim_t z_offset = ez * z_size;
            for (dim_t ey = 0; ey < r.NE[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX3(r.NE[0]-1, ey, ez, r.NE[0], r.NE[1]));
                memcpy(out + numComp, &rrecv[z_offset + ey*4*numComp], pointsize);
                memcpy(out + 3*numComp, &rrecv[z_offset + (ey*4 + 1)*numComp], pointsize);
                memcpy(out + 5*numComp, &rrecv[z_offset + (ey*4 + 2)*numComp], pointsize);
                memcpy(out + 7*numComp, &rrecv[z_offset + (ey*4 + 3)*numComp], pointsize);
            }
        }
    } else if (upper == THEIRS) {
        const dim_t z_size = r.NE[1]*numComp*8;
#pragma omp parallel for
        for (dim_t ez = 0; ez < r.NE[2]; ez++) {
            const dim_t z_offset = ez * z_size;
            for (dim_t ey = 0; ey < r.NE[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX3(r.NE[0]-1, ey, ez, r.NE[0], r.NE[1]));
                memcpy(out, &rrecv[z_offset + ey*8*numComp], pointsize*8);
            }
        }
    }
#endif
}

void RipleyCoupler::shareBrickYFaces(struct Ripley& r, int hasLower,
        int hasUpper, int lower, int upper, escript::Data& target) const
{

#ifdef ESYS_MPI
    const size_t pointsize = numComp * sizeof(double);
    //fill outbound arrays
    const dim_t bottomCount = r.NE[2]*r.NE[0]*numComp*4 * (1 + lower*lower);
    const dim_t topCount = r.NE[2]*r.NE[0]*numComp*4 * (1 + upper*upper);
    std::vector<double> bottom(bottomCount, 0);
    std::vector<double> top(topCount, 0);
    std::vector<double> brecv(bottomCount, 0);
    std::vector<double> trecv(topCount, 0);

    if (lower == SHARED) { //single
#pragma omp parallel for
        for (dim_t ez = 0; ez < r.NE[2]; ez++) {
            for (dim_t ex = 0; ex < r.NE[0]; ex++) {
                double *out = target.getSampleDataRW(INDEX3(ex, 0, ez, r.NE[0], r.NE[1]));
                memcpy(&bottom[ez*4*numComp*r.NE[0] + ex*4*numComp], out + 2*numComp, pointsize*2);
                memcpy(&bottom[ez*4*numComp*r.NE[0] + ex*4*numComp + 2*numComp], out + 6*numComp, pointsize*2);
            }
        }
    } else if (hasLower && lower == MINE) { //double
        for (dim_t ez = 0; ez < r.NE[2]; ez++) {
            double *out = target.getSampleDataRW(INDEX3(0, 0, ez, r.NE[0], r.NE[1]));
            memcpy(&bottom[ez*8*numComp*r.NE[0]], out, pointsize*8*r.NE[0]);
        }
    }
    if (upper == SHARED) { //single
#pragma omp parallel for
        for (dim_t ez = 0; ez < r.NE[2]; ez++) {
            for (dim_t ex = 0; ex < r.NE[0]; ex++) {
                double *out = target.getSampleDataRW(INDEX3(ex, r.NE[1]-1, ez, r.NE[0], r.NE[1]));
                memcpy(&top[ez*4*numComp*r.NE[0] + ex*4*numComp], out, pointsize*4);
                memcpy(&top[ez*4*numComp*r.NE[0] + ex*4*numComp + 2*numComp], out + 4*numComp, pointsize*2);
            }
        }
    } else if (hasUpper && upper == MINE) { //double
        for (dim_t ez = 0; ez < r.NE[2]; ez++) {
            double *out = target.getSampleDataRW(INDEX3(0, r.NE[1]-1, ez, r.NE[0], r.NE[1]));
            memcpy(&top[ez*8*numComp*r.NE[0]], out, pointsize*8*r.NE[0]);
        }
    }

    shareWithNeighbours((rank / s_NX[0]) % 2, hasLower, hasUpper, &bottom[0],
            &top[0], &brecv[0], &trecv[0], bottomCount, topCount, s_NX[0]);

    //unpacking
    if (lower == SHARED) {
#pragma omp parallel for
        for (dim_t ez = 0; ez < r.NE[2]; ez++) {
            for (dim_t ex = 0; ex < r.NE[0]; ex++) {
                double *out = target.getSampleDataRW(INDEX3(ex, 0, ez, r.NE[0], r.NE[1]));
                memcpy(out, &brecv[INDEX4(0,0,ex,ez,numComp,4,r.NE[0])], pointsize*2);
                memcpy(out + 4*numComp, &brecv[INDEX4(0,2,ex,ez,numComp,4,r.NE[0])], pointsize*2);
            }
        }
    } else if (lower == THEIRS) {
        for (dim_t ez = 0; ez < r.NE[2]; ez++) {
            double *out = target.getSampleDataRW(INDEX3(0, 0, ez, r.NE[0], r.NE[1]));
            memcpy(out, &brecv[ez*8*numComp*r.NE[0]], pointsize*8*r.NE[0]);
        }
    }
    if (upper == SHARED) {
#pragma omp parallel for
        for (dim_t ez = 0; ez < r.NE[2]; ez++) {
            for (dim_t ex = 0; ex < r.NE[0]; ex++) {
                double *out = target.getSampleDataRW(INDEX3(ex, r.NE[1]-1, ez, r.NE[0], r.NE[1]));

                memcpy(out + 2*numComp, &trecv[INDEX4(0,0,ex,ez,numComp,4,r.NE[0])], pointsize*2);
                memcpy(out + 6*numComp, &trecv[INDEX4(0,2,ex,ez,numComp,4,r.NE[0])], pointsize*2);
            }
        }
    } else if (upper == THEIRS) {
        for (dim_t ez = 0; ez < r.NE[2]; ez++) {
            double *out = target.getSampleDataRW(INDEX3(0, r.NE[1]-1, ez, r.NE[0], r.NE[1]));
            memcpy(out, &trecv[ez*8*numComp*r.NE[0]], pointsize*8*r.NE[0]);
        }
    }
#endif
}


void RipleyCoupler::shareBrickZFaces(struct Ripley& r, int hasLower,
        int hasUpper, int lower, int upper, escript::Data& target) const
{
#ifdef ESYS_MPI
    const size_t pointsize = numComp * sizeof(double);
    //fill outbound arrays
    const dim_t bottomCount = r.NE[1]*r.NE[0]*numComp*4 * (1 + lower*lower);
    const dim_t topCount = r.NE[1]*r.NE[0]*numComp*4 * (1 + upper*upper);
    std::vector<double> bottom(bottomCount, 0);
    std::vector<double> top(topCount, 0);
    std::vector<double> brecv(bottomCount, 0);
    std::vector<double> trecv(topCount, 0);

    if (lower == SHARED) { //single
#pragma omp parallel for
        for (dim_t ey = 0; ey < r.NE[1]; ey++) {
            for (dim_t ex = 0; ex < r.NE[0]; ex++) {
                double *out = target.getSampleDataRW(INDEX3(ex, ey, 0, r.NE[0], r.NE[1]));
                memcpy(&bottom[INDEX4(0,0,ex,ey,numComp,4,r.NE[0])], out + 4*numComp, pointsize*4);
            }
        }
    } else if (hasLower && lower == MINE) { //double
        double *out = target.getSampleDataRW(0);
        memcpy(&bottom[0], out, pointsize*8*r.NE[0]*r.NE[1]);
    }
    if (upper == SHARED) { //single
#pragma omp parallel for
        for (dim_t ey = 0; ey < r.NE[1]; ey++) {
            for (dim_t ex = 0; ex < r.NE[0]; ex++) {
                double *out = target.getSampleDataRW(INDEX3(ex,ey,r.NE[2]-1,r.NE[0],r.NE[1]));
                memcpy(&top[INDEX4(0,0,ex,ey,numComp,4,r.NE[0])], out, pointsize*4);
            }
        }
    } else if (hasUpper && upper == MINE) { //double
        double *out = target.getSampleDataRW((r.NE[2]-1)*r.NE[0]*r.NE[1]);
        memcpy(&top[0], out, pointsize*8*r.NE[0]*r.NE[1]);
    }

    shareWithNeighbours((rank / s_NX[0]*s_NX[1]) % 2, hasLower, hasUpper, &bottom[0],
            &top[0], &brecv[0], &trecv[0], bottomCount, topCount, s_NX[0]*s_NX[1]);

    //unpacking
    if (lower == SHARED) {
#pragma omp parallel for
        for (dim_t ey = 0; ey < r.NE[1]; ey++) {
            for (dim_t ex = 0; ex < r.NE[0]; ex++) {
                double *out = target.getSampleDataRW(ex + ey*r.NE[0]);
                memcpy(out, &brecv[INDEX4(0,0,ex,ey,numComp,4,r.NE[0])], pointsize*4);
            }
        }
    } else if (lower == THEIRS) {
        double *out = target.getSampleDataRW(0);
        memcpy(out, &brecv[0], pointsize*8*r.NE[0]*r.NE[1]);
    }
    if (upper == SHARED) {
#pragma omp parallel for
        for (dim_t ey = 0; ey < r.NE[1]; ey++) {
            for (dim_t ex = 0; ex < r.NE[0]; ex++) {
                double *out = target.getSampleDataRW(INDEX3(ex,ey,r.NE[2]-1,r.NE[0],r.NE[1]));
                memcpy(out + 4*numComp, &trecv[INDEX4(0,0,ex,ey,numComp,4,r.NE[0])], pointsize*4);
            }
        }
    } else if (upper == THEIRS) {
        double *out = target.getSampleDataRW((r.NE[2]-1)*r.NE[0]*r.NE[1]);
        memcpy(out, &trecv[0], pointsize*8*r.NE[0]*r.NE[1]);
    }
#endif
}

void RipleyCoupler::shareRectangleXEdges(struct Ripley& r,
        int hasLower, int hasUpper, int lower, int upper,
        escript::Data& target) const
{
#ifdef ESYS_MPI
    const dim_t leftCount = r.NE[1]*numComp*2 * (1 + lower*lower);
    const dim_t rightCount = r.NE[1]*numComp*2 * (1 + upper*upper);
    std::vector<double> left(leftCount, 0);
    std::vector<double> right(rightCount, 0);
    std::vector<double> rrecv(rightCount, 0);
    std::vector<double> lrecv(leftCount, 0);
    const size_t pointsize = numComp * sizeof(double);

    //fill outbound arrays
    if (lower == SHARED) { //single
#pragma omp parallel for
        for (dim_t ey = 0; ey < r.NE[1]; ey++) {
            double *out = target.getSampleDataRW(INDEX2(0, ey, r.NE[0]));
            memcpy(&left[ey*2*numComp], out + numComp, pointsize);
            memcpy(&left[(ey*2 + 1)*numComp], out + 3*numComp, pointsize);
        }
    } else if (hasLower && lower == MINE) { //double
#pragma omp parallel for
        for (dim_t ey = 0; ey < r.NE[1]; ey++) {
            double *out = target.getSampleDataRW(INDEX2(0, ey, r.NE[0]));
            memcpy(&left[ey*4*numComp], out, pointsize*4);
        }
    }
    if (upper == SHARED) { //single
#pragma omp parallel for
        for (dim_t ey = 0; ey < r.NE[1]; ey++) {
            double *out = target.getSampleDataRW(INDEX2(r.NE[0]-1, ey, r.NE[0]));
            memcpy(&right[ey*2*numComp], out, pointsize);
            memcpy(&right[(ey*2 + 1)*numComp], out + 2*numComp, pointsize);
        }
    } else if (hasUpper && upper == MINE) { //double
#pragma omp parallel for
        for (dim_t ey = 0; ey < r.NE[1]; ey++) {
            double *out = target.getSampleDataRW(INDEX2(r.NE[0]-1, ey, r.NE[0]));
            memcpy(&right[ey*4*numComp], out, pointsize*4);
        }
    }

    shareWithNeighbours((rank % s_NX[0]) % 2, hasLower, hasUpper, &left[0],
            &right[0], &lrecv[0], &rrecv[0], leftCount, rightCount, 1);

    //unpacking
    if (lower == SHARED) {
#pragma omp parallel for
        for (dim_t ey = 0; ey < r.NE[1]; ey++) {
            double *out = target.getSampleDataRW(INDEX2(0, ey, r.NE[0]));
            memcpy(out, &lrecv[ey*2*numComp], pointsize);
            memcpy(out + 2*numComp, &lrecv[(ey*2 + 1)*numComp], pointsize);
        }
    } else if (lower == THEIRS) {
#pragma omp parallel for
        for (dim_t ey = 0; ey < r.NE[1]; ey++) {
            double *out = target.getSampleDataRW(INDEX2(0, ey, r.NE[0]));
            memcpy(out, &lrecv[ey*4*numComp], pointsize*4);
        }
    }
    if (upper == SHARED) {
#pragma omp parallel for
        for (dim_t ey = 0; ey < r.NE[1]; ey++) {
            double *out = target.getSampleDataRW(INDEX2(r.NE[0]-1, ey, r.NE[0]));
            memcpy(out + numComp, &rrecv[ey*2*numComp], pointsize);
            memcpy(out + 3*numComp, &rrecv[(ey*2 + 1)*numComp], pointsize);
        }
    } else if (upper == THEIRS) {
#pragma omp parallel for
        for (dim_t ey = 0; ey < r.NE[1]; ey++) {
            double *out = target.getSampleDataRW(INDEX2(r.NE[0]-1, ey, r.NE[0]));
            memcpy(out, &rrecv[ey*4*numComp], pointsize*4);
        }
    }
#endif
}

void RipleyCoupler::shareRectangleYEdges(struct Ripley& r, int hasLower,
        int hasUpper, int lower, int upper, escript::Data& target) const
{
#ifdef ESYS_MPI
    const size_t pointsize = numComp * sizeof(double);
    //fill outbound arrays
    const dim_t bottomCount = r.NE[0]*numComp*2 * (1 + lower*lower);
    const dim_t topCount = r.NE[0]*numComp*2 * (1 + upper*upper);
    std::vector<double> bottom(bottomCount, 0);
    std::vector<double> top(topCount, 0);
    std::vector<double> brecv(bottomCount, 0);
    std::vector<double> trecv(topCount, 0);

    if (lower == SHARED) { //single
#pragma omp parallel for
        for (dim_t ex = 0; ex < r.NE[0]; ex++) {
            double *out = target.getSampleDataRW(ex);
            memcpy(&bottom[ex*2*numComp], out + 2*numComp, pointsize*2);
        }
    } else if (hasLower && lower == MINE) { //double
        double *out = target.getSampleDataRW(0);
        memcpy(&bottom[0], out, pointsize*4*r.NE[0]);
    }
    if (upper == SHARED) { //single
#pragma omp parallel for
        for (dim_t ex = 0; ex < r.NE[0]; ex++) {
            double *out = target.getSampleDataRW(INDEX2(ex,r.NE[1]-1,r.NE[0]));
            memcpy(&top[ex*2*numComp], out, pointsize*2);
        }
    } else if (hasUpper && upper == MINE) { //double
        double *out = target.getSampleDataRW((r.NE[1]-1)*r.NE[0]);
        memcpy(&top[0], out, pointsize*4*r.NE[0]);
    }
    //now share
    shareWithNeighbours((rank / s_NX[0]) % 2, hasLower, hasUpper, &bottom[0],
            &top[0], &brecv[0], &trecv[0], bottomCount, topCount, s_NX[0]);
    //unpacking
    if (lower == SHARED) {
#pragma omp parallel for
        for (dim_t ex = 0; ex < r.NE[0]; ex++) {
            double *out = target.getSampleDataRW(ex);
            memcpy(out, &brecv[ex*2*numComp], pointsize*2);
        }
    } else if (lower == THEIRS) {
        double *out = target.getSampleDataRW(0);
        memcpy(out, &brecv[0], pointsize*4*r.NE[0]);
    }
    if (upper == SHARED) {
#pragma omp parallel for
        for (dim_t ex = 0; ex < r.NE[0]; ex++) {
            double *out = target.getSampleDataRW(INDEX2(ex,r.NE[1]-1,r.NE[0]));
            memcpy(out + 2*numComp, &trecv[ex*2*numComp], pointsize*2);
        }
    } else if (upper == THEIRS) {
        double *out = target.getSampleDataRW((r.NE[1]-1)*r.NE[0]);
        memcpy(out, &trecv[0], pointsize*4*r.NE[0]);
    }
#endif
}

void RipleyCoupler::shareWithNeighbours(bool lowerFirst, int hasLower,
        int hasUpper, double *bottom, double *top, double *brecv, double *trecv,
        int bSize, int tSize, int distance) const
{
#ifdef ESYS_MPI
    const int above = rank + distance;
    const int below = rank - distance;
    MPI_Status status;

    if (lowerFirst) {
        //down first
        if (hasLower) {
            MPI_Sendrecv(bottom, bSize, MPI_DOUBLE, below, below,
                    brecv, bSize, MPI_DOUBLE, below, rank,
                    comm, &status);

        }
        if (hasUpper) {
            MPI_Sendrecv(top, tSize, MPI_DOUBLE, above, above,
                    trecv, tSize, MPI_DOUBLE, above, rank,
                    comm, &status);
        }
    } else {
        //up first
        if (hasUpper) {
            MPI_Sendrecv(top, tSize, MPI_DOUBLE, above, above,
                    trecv, tSize, MPI_DOUBLE, above, rank,
                    comm, &status);
        }
        if (hasLower) {
            MPI_Sendrecv(bottom, bSize, MPI_DOUBLE, below, below,
                    brecv, bSize, MPI_DOUBLE, below, rank,
                    comm, &status);
        }
    }
#endif
}

} // end of namespace speckley

