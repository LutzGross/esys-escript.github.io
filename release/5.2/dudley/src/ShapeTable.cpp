/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include "ShapeTable.h"

namespace dudley {

// Joel Fenwick - derived from info in Finley's Quadrature and shape files

// This method is not thread-safe unless the initial call has completed.
// Evaluates the shape functions at nodes (This is the S value from the finley
// ShapeFunctions). The dim argument is the dimension of the element not the
// dimension of the embedding space. The reduced arg is whether the elements
// are reduced or not
bool getQuadShape(dim_t dim, bool reduced, const double **shapearr)
{
#define _dudley_s_alpha 0.58541019662496852
#define _dudley_s_beta  0.1381966011250105

    // {Line, TRI, TET} X {single_quad_point, more} X max number of quadpoints
    static const double _dudley_V[3 * 2][12] = {
        {0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // Line single
        {(1. - .577350269189626) / 2., (1. + .577350269189626) / 2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},     // Line 2 points
        {1 / 3., 1 / 3., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // Tri single
        {0.5, 0, 0, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0, 0},   // Tri 3 points
        {0.25, 0.25, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0},  // Tet single
        {_dudley_s_beta, _dudley_s_beta, _dudley_s_beta,
         _dudley_s_alpha, _dudley_s_beta, _dudley_s_beta,
         _dudley_s_beta, _dudley_s_alpha, _dudley_s_beta,
         _dudley_s_beta, _dudley_s_beta, _dudley_s_alpha} // Tet 4 points
    };

#undef _dudley_s_alpha
#undef _dudley_s_beta

    static double **arr = NULL;

    if (!arr) {
        int i;
        arr = new double*[8];   // point occupies two slots to make things simpler
        arr[0] = new double[1];
        arr[0][0] = 1.;         // point
        arr[1] = arr[0];
        arr[2] = new double[4]; // Line Single
        arr[3] = new double[4]; // Line 2

        for (i = 0; i < 2; ++i) {
            arr[2][2 * i] = 1 - _dudley_V[0][i];
            arr[2][2 * i + 1] = _dudley_V[0][i];
            arr[3][2 * i] = 1 - _dudley_V[1][i];
            arr[3][2 * i + 1] = _dudley_V[1][i];
        }

        arr[4] = new double[3]; // Tri single
        arr[4][0] = 1. - _dudley_V[2][0] - _dudley_V[2][1];
        arr[4][1] = _dudley_V[2][0];
        arr[4][2] = _dudley_V[2][1];

        arr[5] = new double[9]; // Tri 3
        for (i = 0; i < 3; ++i) {
            arr[5][3 * i] = 1 - _dudley_V[3][2 * i] - _dudley_V[3][2 * i + 1];
            arr[5][3 * i + 1] = _dudley_V[3][2 * i];
            arr[5][3 * i + 2] = _dudley_V[3][2 * i + 1];
        }
        arr[6] = new double[4]; // Tet single
        arr[6][0] = 1 - _dudley_V[4][0] - _dudley_V[4][1] - _dudley_V[4][2];
        arr[6][1] = _dudley_V[4][0];
        arr[6][2] = _dudley_V[4][1];
        arr[6][3] = _dudley_V[4][2];

        arr[7] = new double[16]; // Tet 4
        for (i = 0; i < 4; ++i) {
            const double x = _dudley_V[5][3 * i];
            const double y = _dudley_V[5][3 * i + 1];
            const double z = _dudley_V[5][3 * i + 2];
            arr[7][4 * i] = 1 - x - y - z;
            arr[7][4 * i + 1] = x;
            arr[7][4 * i + 2] = y;
            arr[7][4 * i + 3] = z;
        }
    } // end if

    if (dim > -1 && dim < 4) {
        *shapearr = arr[(!reduced) ? (2 * dim + 1) : (2 * dim)];
        return 1;
    }
    *shapearr = NULL;
    return 0;
}

const char *getElementName(ElementTypeId id)
{
    switch (id) {
        case Dudley_Point1:
            return "Point1";
        case Dudley_Line2:
            return "Line2";
        case Dudley_Tri3:
            return "Tri3";
        case Dudley_Tet4:
            return "Tet4";
        case Dudley_Line2Face:
            return "Line2Face";
        case Dudley_Tri3Face:
            return "Tri3Face";
        case Dudley_Tet4Face:
            return "Tet4Face";
        default:
            return "noElement";
    }
}

} // namespace dudley

