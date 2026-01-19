
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/


/****************************************************************************/

/* Paso: Pattern: Pattern_mis

   Searches for a maximal independent set MIS in the matrix pattern.
   Vertices in the maximal independent set are marked in mis_marker.
   Nodes to be considered are marked by -1 on the input in mis_marker.

*/
/****************************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005                      */
/* Author: Lutz Gross, l.gross@uq.edu.au                              */

/****************************************************************************/

#include "Pattern.h"
#include "PasoException.h"
#include "PasoUtil.h"

namespace paso {

// used to generate pseudo random numbers
static double Pattern_mis_seed=.4142135623730951;

#define IS_AVAILABLE -1
#define IS_IN_MIS_NOW -2
#define IS_IN_MIS -3
#define IS_CONNECTED_TO_MIS -4

void Pattern::mis(index_t* mis_marker) const
{
    const index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    if (numOutput != numInput) {
        throw PasoException("Pattern::mis: pattern must be square.");
    }

    const dim_t n = numOutput;
    double* value = new double[n];

    // is there any vertex available?
    while (util::isAny(n, mis_marker, IS_AVAILABLE)) {
        /* step 1: assign a random number in [0,1] to each vertex
         * step 2: if the vertex is available, check if its value is smaller
         *         than all values of its neighbours. If the answer is yes,
         *         the vertex is put into the independent set and all its
         *         neighbours are removed from the graph by setting
         *         mis_marker to false
         */
        /* assign random number in [0,1] to each vertex */
#pragma omp parallel for schedule(static)
        for (dim_t i=0; i < n; ++i) {
            if (mis_marker[i]==IS_AVAILABLE) {
                    value[i]=fmod(Pattern_mis_seed*(i+1),1.);
            } else {
                    value[i]=2.;
            }
        }
        // update the seed
        // Pattern_mis_seed=fmod(sqrt(Pattern_mis_seed*(n+1)),1.);

        // detect independent vertices as those vertices that have a value
        // less than all values of its neighbours
#pragma omp parallel for schedule(static)
        for (dim_t i=0; i < n; ++i) {
            if (mis_marker[i]==IS_AVAILABLE) {
                index_t flag=IS_IN_MIS_NOW;
                for (index_t iptr=ptr[i]-index_offset; iptr<ptr[i+1]-index_offset; ++iptr) {
                    const index_t naib=index[iptr]-index_offset;
                    if (naib != i && value[naib] <= value[i]) {
                        flag=IS_AVAILABLE;
                        break;
                    }
                }
                mis_marker[i]=flag;
            }
        }
        // detect independent vertices as those vertices that have a value
        // less than all values of its neighbours
#pragma omp parallel for schedule(static)
        for (dim_t i=0; i < n; i++) {
            if (mis_marker[i]==IS_IN_MIS_NOW) {
                for (index_t iptr=ptr[i]-index_offset; iptr<ptr[i+1]-index_offset; ++iptr) {
                    const index_t naib=index[iptr]-index_offset;
                    if (naib != i)
                        mis_marker[naib]=IS_CONNECTED_TO_MIS;
                }
                mis_marker[i]=IS_IN_MIS;
            }
        }
    }
    // swap to true/false in mis_marker
#pragma omp parallel for schedule(static)
    for (dim_t i=0; i < n; i++)
        mis_marker[i] = (mis_marker[i]==IS_IN_MIS);
    delete[] value;
}
#undef IS_AVAILABLE
#undef IS_IN_MIS_NOW
#undef IS_IN_MIS
#undef IS_CONNECTED_TO_MIS

} // namespace paso

