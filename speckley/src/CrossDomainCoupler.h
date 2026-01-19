
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __SPECKLEY_CROSSDOMAINCOUPLER_H__
#define __SPECKLEY_CROSSDOMAINCOUPLER_H__

#include <speckley/Brick.h>
#include <speckley/Rectangle.h>

#include <ripley/Brick.h>
#include <ripley/Rectangle.h>

namespace speckley {

class RipleyCoupler {
public:
    RipleyCoupler(const SpeckleyDomain *speck, const double s_dx[2], int rank);


    void interpolate(escript::Data& target, const escript::Data& source) const;
private:
    // a struct type to hold all the relevant info on the target domain
    struct Ripley {
        const ripley::RipleyDomain *domain;
        double dx[3];
        dim_t NE[3];
        dim_t mins[3];
        dim_t maxs[3];
    };
    void calculateOrder2(int dim, double loc, double *results) const;
    void calculateOrder3(int dim, double loc, double *results) const;
    void calculateOrder4(int dim, double loc, double *results) const;
    void calculateOrder5(int dim, double loc, double *results) const;
    void calculateOrder6(int dim, double loc, double *results) const;
    void calculateOrder7(int dim, double loc, double *results) const;
    void calculateOrder8(int dim, double loc, double *results) const;
    void calculateOrder9(int dim, double loc, double *results) const;
    void calculateOrder10(int dim, double loc, double *results) const;

    void generateLocations(struct Ripley& r, double **positions) const;

    bool validInterpolation(escript::Data& target, const escript::Data& source,
            const SpeckleyDomain *speck, const double *s_dx,
            const ripley::RipleyDomain *other) const;
    void calculate(struct Ripley& r, dim_t ex, dim_t ey, dim_t ez,
            int oqx, int oqy, int oqz, double *out, const double *factor_x,
            const double *factor_y, const double *factor_z,
            const escript::Data& source) const;

    void shareWithNeighbours(bool lowerFirst, int hasLower, int hasUpper,
            double *bottom, double *top, double *brecv, double *trecv,
            int bSize, int tSize, int distance) const;

    void getEdgeSpacing(struct Ripley r, int *lower, int *upper) const;

    void shareBrickXFaces(struct Ripley& r, int hasLower,
            int hasUpper, int lower, int upper, escript::Data& target) const;
    void shareBrickYFaces(struct Ripley& r, int hasLower,
            int hasUpper, int lower, int upper, escript::Data& target) const;
    void shareBrickZFaces(struct Ripley& r, int hasLower,
            int hasUpper, int lower, int upper, escript::Data& target) const;

    void shareRectangleXEdges(struct Ripley& r, int hasLower,
            int hasUpper, int lower, int upper, escript::Data& target) const;
    void shareRectangleYEdges(struct Ripley& r, int hasLower,
            int hasUpper, int lower, int upper, escript::Data& target) const;
    //speckley info
    const SpeckleyDomain *speck;
    dim_t s_NE[3];
    double s_dx[3];
    int s_NX[3];
    double speckley_origin[3];
    int order;
    int numQuads;

    //coupling info
    bool hasLower[3];
    bool hasUpper[3];

    //per interpolation
    mutable int numComp;

#ifdef ESYS_MPI
    int rank;
    MPI_Comm comm;
#endif

};

/**
   \brief
   interpolates data given on source onto target where source and target
   are given on different domains
*/
void interpolateAcross3D(escript::Data& target, const escript::Data& source,
                    const Brick *speck, const double s_dx[3], int rank,
                    MPI_Comm comm);

bool probeInterpolationAcross(int fsType_source,
        const escript::AbstractDomain& domain, int fsType_target, int dim);

} // end of namespace speckley

#endif // __SPECKLEY_CROSSDOMAINCOUPLER_H__

