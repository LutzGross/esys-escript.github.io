
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

/* Paso: SystemMatrixPattern::unrollBlocks */

/****************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "SystemMatrixPattern.h"

namespace paso {

SystemMatrixPattern_ptr SystemMatrixPattern::unrollBlocks(
                        int newType, dim_t output_block_size,
                        dim_t input_block_size)
{
    SystemMatrixPattern_ptr out;
    escript::Distribution_ptr new_output_distribution, new_input_distribution;
    Connector_ptr new_col_connector, new_row_connector;

    if ( (output_block_size == 1) && (input_block_size == 1) &&
            ((type & MATRIX_FORMAT_OFFSET1) == (newType & MATRIX_FORMAT_OFFSET1)) ) {
        out = shared_from_this();
    } else {
        Pattern_ptr new_mainPattern(mainPattern->unrollBlocks(newType,
                output_block_size, input_block_size));
        Pattern_ptr new_col_couplePattern(col_couplePattern->unrollBlocks(
                newType, output_block_size, input_block_size));
        Pattern_ptr new_row_couplePattern(row_couplePattern->unrollBlocks(
                newType, output_block_size, input_block_size));
        if (output_block_size > 1) {
            new_output_distribution.reset(new escript::Distribution(
                    output_distribution->mpi_info,
                    output_distribution->first_component,
                    output_block_size, 0));
            new_row_connector = row_connector->unroll(output_block_size);
        } else {
            new_output_distribution = output_distribution;
            new_row_connector = row_connector;
        }
        if (input_block_size > 1) {
            new_input_distribution.reset(new escript::Distribution(
                    input_distribution->mpi_info,
                    input_distribution->first_component,
                    input_block_size, 0));
            new_col_connector = col_connector->unroll(input_block_size);
        } else {
            new_input_distribution = input_distribution;
            new_col_connector = col_connector;
        }

        out.reset(new SystemMatrixPattern(newType, new_output_distribution,
                                          new_input_distribution,
                                          new_mainPattern,
                                          new_col_couplePattern,
                                          new_row_couplePattern,
                                          new_col_connector,
                                          new_row_connector));
    }

    return out;
}

} // namespace paso

