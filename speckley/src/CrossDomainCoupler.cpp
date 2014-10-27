
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
#include <speckley/CrossDomainCoupler.h>

#include <speckley/lagrange_functions.h>
#include <esysUtils/index.h>
#include <esysUtils/Esys_MPI.h>

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

RectangleCoupler::RectangleCoupler(const speckley::Rectangle *speck, 
        const double s_dx[2], int rank
#ifdef ESYS_MPI
        , MPI_Comm comm) : speck(speck), rank(rank), comm(comm)
#else
        ) : speck(speck), rank(rank)
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
    order = speck->getOrder();
    numQuads = order + 1;
}

bool RectangleCoupler::validInterpolation(escript::Data& target,
            const escript::Data& source,
            const SpeckleyDomain *speck, const double *s_dx,
            const ripley::RipleyDomain *other) const
{
    if (source.getDomain().get() != speck)
        throw SpeckleyException(
            "ripleyCoupler: interpolation from unsupported domain");

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

    return true;
}

inline void RectangleCoupler::calculate(struct Ripley& r, dim_t ex, dim_t ey,
        int oqx, int oqy, double *out, double *factor_x, double *factor_y,
        const escript::Data& source) const
{
    double x = r.domain->getLocalCoordinate(ex, 0) - speckley_origin[0]
                + ripleyLocations[oqx]*r.dx[0];
    double y = r.domain->getLocalCoordinate(ey, 1) - speckley_origin[1]
                + ripleyLocations[oqy]*r.dx[1];
    dim_t source_ex = x / s_dx[0];
    dim_t source_ey = y / s_dx[1];
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
}

void RectangleCoupler::generateLocations(struct Ripley& r, double **positions) const
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
            throw SpeckleyException("RectangleCoupler:: unexpected order of domain, can only be implementation error");
    }
}

void RectangleCoupler::calculateOrder2(int dim, double loc, double *results) const
{
    loc -= speckley_origin[dim];
    dim_t source_e = loc / s_dx[dim];
    double local = ((loc - source_e * s_dx[dim]) / s_dx[dim]) * 2 - 1;
    results[0] = lagrange_degree2_0(local);
    results[1] = lagrange_degree2_1(local);
    results[2] = lagrange_degree2_2(local);
}

void RectangleCoupler::calculateOrder3(int dim, double loc, double *results) const
{
    loc -= speckley_origin[dim];
    dim_t source_e = loc / s_dx[dim];
    double local = ((loc - source_e * s_dx[dim]) / s_dx[dim]) * 2 - 1;
    results[0] = lagrange_degree3_0(local);
    results[1] = lagrange_degree3_1(local);
    results[2] = lagrange_degree3_2(local);
    results[3] = lagrange_degree3_3(local);
}

void RectangleCoupler::calculateOrder4(int dim, double loc, double *results) const
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

void RectangleCoupler::calculateOrder5(int dim, double loc, double *results) const
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

void RectangleCoupler::calculateOrder6(int dim, double loc, double *results) const
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

void RectangleCoupler::calculateOrder7(int dim, double loc, double *results) const
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

void RectangleCoupler::calculateOrder8(int dim, double loc, double *results) const
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

void RectangleCoupler::calculateOrder9(int dim, double loc, double *results) const
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

void RectangleCoupler::calculateOrder10(int dim, double loc, double *results) const
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

void RectangleCoupler::getEdgeSpacing(struct Ripley r, int *lower, int *upper) const
{
    for (int i = 0; i < speck->getDim(); i++) {
        const double first = ripleyLocations[0]*r.dx[i];
        const double second = ripleyLocations[1]*r.dx[i];
        double start = r.domain->getLocalCoordinate(0, i) - speckley_origin[i];
        lower[i] = 0;
        if (start + first > 0)
            lower[i] = 1; //both in this
        else if (start + second < 0)
            lower[i] = -1; //both in neighbour
        start = r.domain->getLocalCoordinate(r.NE[i]-1, i) - speckley_origin[i];
        upper[i] = 0;
        if ((start + first)/s_dx[i] >= s_NE[i])
            upper[i] = -1; //both in neighbour
        else if ((start + second)/s_dx[i] < s_NE[i])
            upper[i] = 1; //both in this
    }
}

void RectangleCoupler::interpolate(escript::Data& target,
        const escript::Data& source) const
{
    //can only interpolate to ripley right now, so check it
    const ripley::Rectangle *other = dynamic_cast<const ripley::Rectangle *>
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
    double *factors_x = &factors_x_vec[0];
    std::vector<double> factors_y_vec(2*r.NE[1]*numQuads);
    double *factors_y = &factors_y_vec[0];
    double *positions[2] = {factors_x, factors_y};
    generateLocations(r, positions);
    
    numComp = source.getDataPointSize();
    target.requireWrite();

    //positional help, 0 = split, -1 = both on neighbour, 1 = both local
    int upper[3];
    int lower[3];
    getEdgeSpacing(r, lower, upper);
    for (int i = 0; i < speck->getDim(); i++) {
        if (hasUpper[i] && upper[i] == 1) {
            r.maxs[i] += 1;
        }
        if (hasLower[i] && lower[i] == 1) {
            r.mins[i] = 0;
        }
    }
    
    //elements that are entirely internal to this rank
#pragma omp parallel for
    for (dim_t ey = r.mins[1]; ey < r.maxs[1]; ey++) {
        for (dim_t ex = r.mins[0]; ex < r.maxs[0]; ex++) {
            double *out = target.getSampleDataRW(INDEX2(ex,ey,r.NE[0]));
            for (int oqy = 0; oqy < 2; oqy++) {
                double *precalc_y = factors_y + (2*ey + oqy) * numQuads;
                for (int oqx = 0; oqx < 2; oqx++) {
                    double *precalc_x = factors_x + (2*ex + oqx) * numQuads;
                    calculate(r, ex, ey, oqx, oqy, out, precalc_x, precalc_y, source);
                }
            }
        }
    }
#ifdef ESYS_MPI
    if (s_NX[0] * s_NX[1] == 1) {
        return;
    }

    const size_t pointsize = numComp * sizeof(double);
    MPI_Status status;
    // communicate across X splits
    if (hasLower[0] || hasUpper[0]) {
        const dim_t leftCount = r.NE[1]*numComp*2 * (1 + lower[0]*lower[0]);
        const dim_t rightCount = r.NE[1]*numComp*2 * (1 + upper[0]*upper[0]);
        std::vector<double> left(leftCount, 0);
        std::vector<double> right(rightCount, 0);
        std::vector<double> rrecv(rightCount, 0);
        std::vector<double> lrecv(leftCount, 0);
        if (lower[0] == 0) { //left, only if an element is split in two
            const dim_t ex = 0;
            const int oqx = 1;
#pragma omp parallel for
            for (dim_t ey = r.mins[1]; ey < r.maxs[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX2(ex,ey,r.NE[0]));
                for (int oqy = 0; oqy < 2; oqy++) {
                    double *precalc_y = factors_y + (2*ey + oqy) * numQuads;
                    double *precalc_x = factors_x + (2*ex + oqx) * numQuads;
                    calculate(r, ex, ey, oqx, oqy, out, precalc_x, precalc_y, source);
                }
            }
            //fill the corners while we're here, if required
            if (lower[1] == 0) { //bottom left element, top right quad
                const dim_t ey = 0;
                double *out = target.getSampleDataRW(0);
                calculate(r, ex, ey, oqx, 1, out, factors_x + numQuads, factors_y + numQuads, source);
            }
            if (upper[1] == 0) { //top left element, bottom right quad
                const dim_t ey = r.NE[1] - 1;
                double *out = target.getSampleDataRW(INDEX2(ex,ey,r.NE[0]));
                calculate(r, ex, ey, oqx, 0, out, factors_x + (2*ex+1)*numQuads, factors_y + 2*ey*numQuads, source);
                
            }
        }
        if (upper[0] == 0) { //right, only if an element is split in two
            const dim_t ex = r.NE[0] - 1;
            const int oqx = 0;
#pragma omp parallel for
            for (dim_t ey = r.mins[1]; ey < r.maxs[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX2(ex,ey,r.NE[0]));
                for (int oqy = 0; oqy < 2; oqy++) {     //for each remote quad
                    double *precalc_y = factors_y + (2*ey + oqy) * numQuads;
                    calculate(r, ex, ey, oqx, oqy, out, factors_x + 2*ex*numQuads, precalc_y, source);
                }
            }
            //fill the corners while we're here, if required
            if (lower[1] == 0) { //bottom right element, top left quad
                const dim_t ey = 0;
                double *out = target.getSampleDataRW(ex);
                calculate(r, ex, ey, oqx, 1, out, factors_x + 2*ex*numQuads, factors_y + numQuads, source);
            }
            if (upper[1] == 0) { //top right element, bottom left quad
                const dim_t ey = r.maxs[1];
                double *out = target.getSampleDataRW(INDEX2(ex,ey,r.NE[0]));
                calculate(r, ex, ey, oqx, 0, out, factors_x + 2*ex*numQuads, factors_y + 2*ey*numQuads, source);
            }
        }
        //fill outbound arrays
        if (lower[0] == 0) { //single 
#pragma omp parallel for
            for (dim_t ey = 0; ey < r.NE[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX2(0, ey, r.NE[0]));
                memcpy(&left[ey*2*numComp], out + numComp, pointsize);
                memcpy(&left[(ey*2 + 1)*numComp], out + 3*numComp, pointsize);
            }
        } else if (lower[0] == 1) { //double
#pragma omp parallel for
            for (dim_t ey = 0; ey < r.NE[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX2(0, ey, r.NE[0]));
                memcpy(&left[ey*4*numComp], out, pointsize*4);
            }
        }
        if (upper[0] == 0) { //single
#pragma omp parallel for
            for (dim_t ey = 0; ey < r.NE[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX2(r.NE[0]-1, ey, r.NE[0]));
                memcpy(&right[ey*2*numComp], out, pointsize);
                memcpy(&right[(ey*2 + 1)*numComp], out + 2*numComp, pointsize);
            }
        } else if (upper[0] == 1) { //double
#pragma omp parallel for
            for (dim_t ey = 0; ey < r.NE[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX2(r.NE[0]-1, ey, r.NE[0]));
                memcpy(&right[ey*4*numComp], out, pointsize*4);
            }
        }
        //now share
        if ((rank % s_NX[0]) % 2) {
            //left first
            if (hasLower[0]) {
                MPI_Sendrecv(&left[0], leftCount, MPI_DOUBLE, rank - 1, 0,
                        &lrecv[0], leftCount, MPI_DOUBLE, rank - 1, 0,
                        comm, &status);

            }
            if (hasUpper[0]) {
                MPI_Sendrecv(&right[0], rightCount, MPI_DOUBLE, rank + 1, 0,
                        &rrecv[0], rightCount, MPI_DOUBLE, rank+1, 0,
                        comm, &status);
            }
        } else {
            //right first
            if (hasUpper[0]) {
                MPI_Sendrecv(&right[0], rightCount, MPI_DOUBLE, rank + 1, 0,
                        &rrecv[0], rightCount, MPI_DOUBLE, rank+1, 0,
                        comm, &status);
            }
            if (hasLower[0]) {
                MPI_Sendrecv(&left[0], leftCount, MPI_DOUBLE, rank - 1, 0,
                        &lrecv[0], leftCount, MPI_DOUBLE, rank - 1, 0,
                        comm, &status);
            }
        }
        //unpacking
        if (lower[0] == 0) {
#pragma omp parallel for
            for (dim_t ey = 0; ey < r.NE[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX2(0, ey, r.NE[0]));
                memcpy(out, &lrecv[ey*2*numComp], pointsize);
                memcpy(out + 2*numComp, &lrecv[(ey*2 + 1)*numComp], pointsize);
            }
        } else if (lower[0] == -1) {
#pragma omp parallel for
            for (dim_t ey = 0; ey < r.NE[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX2(0, ey, r.NE[0]));
                memcpy(out, &lrecv[ey*4*numComp], pointsize*4);
            }
        }
        if (upper[0] == 0) {
#pragma omp parallel for
            for (dim_t ey = 0; ey < r.NE[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX2(r.NE[0]-1, ey, r.NE[0]));
                memcpy(out + numComp, &rrecv[ey*2*numComp], pointsize);
                memcpy(out + 3*numComp, &rrecv[(ey*2 + 1)*numComp], pointsize);
            }
        } else if (upper[0] == -1) {
#pragma omp parallel for
            for (dim_t ey = 0; ey < r.NE[1]; ey++) {
                double *out = target.getSampleDataRW(INDEX2(r.NE[0]-1, ey, r.NE[0]));
                memcpy(out, &rrecv[ey*4*numComp], pointsize*4);
            }
        }
    }

    // communicate across Y splits
    if (hasLower[1] || hasUpper[1]) {
        const dim_t bottomCount = r.NE[0]*numComp*2 * (1 + lower[1]*lower[1]);
        const dim_t topCount = r.NE[0]*numComp*2 * (1 + upper[1]*upper[1]);
        std::vector<double> bottom(bottomCount, 0);
        std::vector<double> top(topCount, 0);
        std::vector<double> brecv(bottomCount, 0);
        std::vector<double> trecv(topCount, 0);
        if (lower[1] == 0) { //bottom, only if an element is split in two
            const dim_t ey = 0;
            const int oqy = 1;
            double *precalc_y = factors_y + (2*ey + oqy) * numQuads;
#pragma omp parallel for
            for (dim_t ex = r.mins[0]; ex < r.maxs[0]; ex++) {
                double *out = target.getSampleDataRW(INDEX2(ex,ey,r.NE[0]));
                for (int oqx = 0; oqx < 2; oqx++) {
                    double *precalc_x = factors_x + (2*ex + oqx) * numQuads;
                    calculate(r, ex, ey, oqx, oqy, out, precalc_x, precalc_y, source);
                }
            }
        }
        if (upper[1] == 0) { //top, only if an element is split in two
            const dim_t ey = r.NE[1] - 1;
            const int oqy = 0;
            double *precalc_y = factors_y + (2*ey + oqy) * numQuads;
#pragma omp parallel for
            for (dim_t ex = r.mins[0]; ex < r.maxs[0]; ex++) {
                double *out = target.getSampleDataRW(INDEX2(ex,ey,r.NE[0]));
                for (int oqx = 0; oqx < 2; oqx++) {     //for each remote quad
                    double *precalc_x = factors_x + (2*ex + oqx) * numQuads;
                    calculate(r, ex, ey, oqx, oqy, out, precalc_x, precalc_y, source);
                }
            }
        }
        //fill outbound arrays
        if (lower[1] == 0) { //single 
#pragma omp parallel for
            for (dim_t ex = 0; ex < r.NE[0]; ex++) {
                double *out = target.getSampleDataRW(ex);
                memcpy(&bottom[ex*2*numComp], out + 2*numComp, pointsize*2);
            }
        } else if (lower[1] == 1) { //double
            double *out = target.getSampleDataRW(0);
            memcpy(&bottom[0], out, pointsize*4*r.NE[0]);
        }
        if (upper[1] == 0) { //single
#pragma omp parallel for
            for (dim_t ex = 0; ex < r.NE[0]; ex++) {
                double *out = target.getSampleDataRW(INDEX2(ex,r.NE[1]-1,r.NE[0]));
                memcpy(&top[ex*2*numComp], out, pointsize*2);
            }
        } else if (upper[1] == 1) { //double
            double *out = target.getSampleDataRW((r.NE[1]-1)*r.NE[0]);
            memcpy(&top[0], out, pointsize*4*r.NE[0]);
        }
        //now share
        if ((rank / s_NX[0]) % 2) {
            //down first
            if (hasLower[1]) {
                MPI_Sendrecv(&bottom[0], bottomCount, MPI_DOUBLE, rank - s_NX[0], 0,
                        &brecv[0], bottomCount, MPI_DOUBLE, rank - s_NX[0], 0,
                        comm, &status);

            }
            if (hasUpper[1]) {
                MPI_Sendrecv(&top[0], topCount, MPI_DOUBLE, rank + s_NX[0], 0,
                        &trecv[0], topCount, MPI_DOUBLE, rank + s_NX[0], 0,
                        comm, &status);
            }
        } else {
            //up first
            if (hasUpper[1]) {
                MPI_Sendrecv(&top[0], topCount, MPI_DOUBLE, rank + s_NX[0], 0,
                        &trecv[0], topCount, MPI_DOUBLE, rank+s_NX[0], 0,
                        comm, &status);
            }
            if (hasLower[1]) {
                MPI_Sendrecv(&bottom[0], bottomCount, MPI_DOUBLE, rank - s_NX[0], 0,
                        &brecv[0], bottomCount, MPI_DOUBLE, rank - s_NX[0], 0,
                        comm, &status);
            }
        }
        //unpacking
        if (lower[1] == 0) {
#pragma omp parallel for
            for (dim_t ex = 0; ex < r.NE[0]; ex++) {
                double *out = target.getSampleDataRW(ex);
                memcpy(out, &brecv[ex*2*numComp], pointsize*2);
            }
        } else if (lower[1] == -1) {
            double *out = target.getSampleDataRW(0);
            memcpy(out, &brecv[0], pointsize*4*r.NE[0]);
        }
        if (upper[1] == 0) {
#pragma omp parallel for
            for (dim_t ex = 0; ex < r.NE[0]; ex++) {
                double *out = target.getSampleDataRW(INDEX2(ex,r.NE[1]-1,r.NE[0]));
                memcpy(out + 2*numComp, &trecv[ex*2*numComp], pointsize*2);
            }
        } else if (upper[1] == -1) {
            double *out = target.getSampleDataRW((r.NE[1]-1)*r.NE[0]);
            memcpy(out, &trecv[0], pointsize*4*r.NE[0]);
        }
    }
#endif //ESYS_MPI
}



//these should be unrolled
double (*interpolationFuncs[9][11]) (double) = {
    {lagrange_degree2_0, lagrange_degree2_1, lagrange_degree2_2},
    {lagrange_degree3_0, lagrange_degree3_1, lagrange_degree3_2, lagrange_degree3_3},
    {lagrange_degree4_0, lagrange_degree4_1, lagrange_degree4_2, lagrange_degree4_3, lagrange_degree4_4},
    {lagrange_degree5_0, lagrange_degree5_1, lagrange_degree5_2, lagrange_degree5_3, lagrange_degree5_4, lagrange_degree5_5},
    {lagrange_degree6_0, lagrange_degree6_1, lagrange_degree6_2, lagrange_degree6_3, lagrange_degree6_4, lagrange_degree6_5, lagrange_degree6_6},
    {lagrange_degree7_0, lagrange_degree7_1, lagrange_degree7_2, lagrange_degree7_3, lagrange_degree7_4, lagrange_degree7_5, lagrange_degree7_6, lagrange_degree7_7},
    {lagrange_degree8_0, lagrange_degree8_1, lagrange_degree8_2, lagrange_degree8_3, lagrange_degree8_4, lagrange_degree8_5, lagrange_degree8_6, lagrange_degree8_7, lagrange_degree8_8},
    {lagrange_degree9_0, lagrange_degree9_1, lagrange_degree9_2, lagrange_degree9_3, lagrange_degree9_4, lagrange_degree9_5, lagrange_degree9_6, lagrange_degree9_7, lagrange_degree9_8, lagrange_degree9_9},
    {lagrange_degree10_0, lagrange_degree10_1, lagrange_degree10_2, lagrange_degree10_3, lagrange_degree10_4, lagrange_degree10_5, lagrange_degree10_6, lagrange_degree10_7, lagrange_degree10_8, lagrange_degree10_9, lagrange_degree10_10}
    };
   

void interpolateAcross3D(escript::Data& target, const escript::Data& source,
        const Brick *speck, const double s_dx[3], int rank,
        MPI_Comm comm)
{
    const int *s_NX = speck->getNumSubdivisionsPerDim();
    if (s_NX[0] * s_NX[1] * s_NX[2] > 1)
        throw SpeckleyException("speckley::Brick doesn't support multiple ranks");
    //can only interpolate to ripley right now, so check it
    const ripley::Brick& other = dynamic_cast<const ripley::Brick &>
                                                (*(target.getDomain().get()));
    // this throw is defensive in case the validation function at some point
    // returns false
//    if (!validInterpolation(target, source, speck, s_dx, other))
//        throw SpeckleyException("ripleyCoupler: unable to interpolate");
    const dim_t *r_NE = other.getNumElementsPerDim();
    const dim_t *s_NE = speck->getNumElementsPerDim();
    const double r_dx[3] = {other.getLength()[0] / r_NE[0],
                            other.getLength()[1] / r_NE[1],
                            other.getLength()[2] / r_NE[2]};
    const int numComp = source.getDataPointSize();
    const int order = speck->getOrder();
    const int numQuads = order + 1;

    target.requireWrite();

#pragma omp parallel for
    for (dim_t ez = 0; ez < r_NE[2]; ez++) {
        for (dim_t ey = 0; ey < r_NE[1]; ey++) {
            for (dim_t ex = 0; ex < r_NE[0]; ex++) {
                double *out = target.getSampleDataRW(INDEX3(ex,ey,ez,r_NE[0],r_NE[1]));
                //initial position
                for (int oqz = 0; oqz < 2; oqz++) { //for each remote quad
                    for (int oqy = 0; oqy < 2; oqy++) {
                        for (int oqx = 0; oqx < 2; oqx++) {
                            double x = other.getLocalCoordinate(ex, 0) + ripleyLocations[oqx]*r_dx[0];
                            double y = other.getLocalCoordinate(ey, 1) + ripleyLocations[oqy]*r_dx[1];
                            double z = other.getLocalCoordinate(ez, 2) + ripleyLocations[oqz]*r_dx[2];
                            //which source element does it live in
                            dim_t source_ex = x / s_dx[0];
                            dim_t source_ey = y / s_dx[1];
                            dim_t source_ez = z / s_dx[2];
                            //now modify coords to the element reference
                            x = ((x - source_ex * s_dx[0]) / s_dx[0]) * 2 - 1;
                            y = ((y - source_ey * s_dx[1]) / s_dx[1]) * 2 - 1;
                            z = ((z - source_ez * s_dx[2]) / s_dx[2]) * 2 - 1;
                            //MPI TODO: check element belongs to this rank
                            const double *sdata = source.getSampleDataRO(INDEX3(source_ex, source_ey, source_ez,
                                                             s_NE[0], s_NE[1]));
                            //and do the actual interpolation
                            for (int sqz = 0; sqz < numQuads; sqz++) {
                                for (int sqy = 0; sqy < numQuads; sqy++) {
                                    for (int sqx = 0; sqx < numQuads; sqx++) {
                                        for (int comp = 0; comp < numComp; comp++) {
                                            out[INDEX4(comp,oqx,oqy,oqz,numComp,2,2)]
                                                  += sdata[INDEX4(comp,sqx,sqy,sqz,numComp,numQuads,numQuads)]
                                                     * interpolationFuncs[order-2][sqx](x)
                                                     * interpolationFuncs[order-2][sqy](y)
                                                     * interpolationFuncs[order-2][sqz](z);
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
}

} // end of namespace speckley

