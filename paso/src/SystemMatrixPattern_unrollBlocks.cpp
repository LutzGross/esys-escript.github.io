
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


/****************************************************************************/

/* Paso: SystemMatrixPattern_unrollBlocks */

/****************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "SystemMatrixPattern.h"
#include "Paso.h"
#include "esysUtils/error.h"

namespace paso {

SystemMatrixPattern* SystemMatrixPattern_unrollBlocks(
        SystemMatrixPattern* pattern, int type, dim_t output_block_size,
        dim_t input_block_size)
{
    SystemMatrixPattern* out = NULL;
    Pattern *new_mainPattern=NULL, *new_col_couplePattern=NULL, *new_row_couplePattern=NULL;
    Distribution_ptr new_output_distribution, new_input_distribution;
    Connector *new_col_connector=NULL, *new_row_connector=NULL;

    if ( (output_block_size == 1) && (input_block_size == 1) &&
            ((pattern->type & MATRIX_FORMAT_OFFSET1) == (type & MATRIX_FORMAT_OFFSET1)) ) {
        out = SystemMatrixPattern_getReference(pattern);
    } else {
        new_mainPattern = Pattern_unrollBlocks(pattern->mainPattern, type,
                output_block_size, input_block_size);
        new_col_couplePattern = Pattern_unrollBlocks(
                pattern->col_couplePattern, type, output_block_size,
                input_block_size);
        new_row_couplePattern = Pattern_unrollBlocks(
                pattern->row_couplePattern, type, output_block_size,
                input_block_size);
        if (output_block_size > 1) {
            new_output_distribution.reset(new Distribution(
                    pattern->output_distribution->mpi_info,
                    pattern->output_distribution->first_component,
                    output_block_size, 0));
            new_row_connector = Connector_unroll(pattern->row_connector,
                                                      output_block_size);
        } else {
            new_output_distribution = pattern->output_distribution;
            new_row_connector = Connector_getReference(pattern->row_connector);
        }
        if (input_block_size > 1) {
            new_input_distribution.reset(new Distribution(
                    pattern->input_distribution->mpi_info,
                    pattern->input_distribution->first_component,
                    input_block_size, 0));
            new_col_connector = Connector_unroll(pattern->col_connector,
                    input_block_size);
        } else {
            new_input_distribution = pattern->input_distribution;
            new_col_connector = Connector_getReference(pattern->col_connector);
        }

        if (Esys_noError()) {
            out = new SystemMatrixPattern(type, new_output_distribution,
                                          new_input_distribution,
                                          new_mainPattern,
                                          new_col_couplePattern,
                                          new_row_couplePattern,
                                          new_col_connector,
                                          new_row_connector);
        }
        Pattern_free(new_mainPattern);
        Pattern_free(new_col_couplePattern);
        Pattern_free(new_row_couplePattern);
        Connector_free(new_row_connector);
        Connector_free(new_col_connector);
    }

    if (Esys_noError()) {
        return out;
    } else {
        SystemMatrixPattern_free(out);
        return NULL;
    }
}

} // namespace paso

