
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

/* Paso: SystemMatrixPattern::unrollBlocks */

/****************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "SystemMatrixPattern.h"
#include "esysUtils/error.h"

namespace paso {

SystemMatrixPattern_ptr SystemMatrixPattern::unrollBlocks(
                        int newType, dim_t output_block_size,
                        dim_t input_block_size)
{
    SystemMatrixPattern_ptr out;
    Pattern *new_mainPattern=NULL, *new_col_couplePattern=NULL, *new_row_couplePattern=NULL;
    Distribution_ptr new_output_distribution, new_input_distribution;
    Connector_ptr new_col_connector, new_row_connector;

    if ( (output_block_size == 1) && (input_block_size == 1) &&
            ((type & MATRIX_FORMAT_OFFSET1) == (newType & MATRIX_FORMAT_OFFSET1)) ) {
        out = shared_from_this();
    } else {
        new_mainPattern = Pattern_unrollBlocks(mainPattern, newType,
                output_block_size, input_block_size);
        new_col_couplePattern = Pattern_unrollBlocks(
                col_couplePattern, newType, output_block_size,
                input_block_size);
        new_row_couplePattern = Pattern_unrollBlocks(
                row_couplePattern, newType, output_block_size,
                input_block_size);
        if (output_block_size > 1) {
            new_output_distribution.reset(new Distribution(
                    output_distribution->mpi_info,
                    output_distribution->first_component,
                    output_block_size, 0));
            new_row_connector = row_connector->unroll(output_block_size);
        } else {
            new_output_distribution = output_distribution;
            new_row_connector = row_connector;
        }
        if (input_block_size > 1) {
            new_input_distribution.reset(new Distribution(
                    input_distribution->mpi_info,
                    input_distribution->first_component,
                    input_block_size, 0));
            new_col_connector = col_connector->unroll(input_block_size);
        } else {
            new_input_distribution = input_distribution;
            new_col_connector = col_connector;
        }

        if (Esys_noError()) {
            out.reset(new SystemMatrixPattern(newType, new_output_distribution,
                                              new_input_distribution,
                                              new_mainPattern,
                                              new_col_couplePattern,
                                              new_row_couplePattern,
                                              new_col_connector,
                                              new_row_connector));
        }
        Pattern_free(new_mainPattern);
        Pattern_free(new_col_couplePattern);
        Pattern_free(new_row_couplePattern);
    }

    if (!Esys_noError()) {
        return SystemMatrixPattern_ptr();
    }
    return out;
}

} // namespace paso

