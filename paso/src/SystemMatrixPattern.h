
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


/****************************************************************************/

/*   Paso: system matrix pattern                              */

/****************************************************************************/

/*   Copyrights by ACcESS Australia 2004,2005 */
/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_SYSTEMMATRIXPATTERN_H__
#define __PASO_SYSTEMMATRIXPATTERN_H__

#include "Coupler.h"
#include "Pattern.h"

#include <escript/Distribution.h>

namespace paso {

struct SystemMatrixPattern;
typedef boost::shared_ptr<SystemMatrixPattern> SystemMatrixPattern_ptr;
typedef boost::shared_ptr<const SystemMatrixPattern> const_SystemMatrixPattern_ptr;

struct PASO_DLL_API SystemMatrixPattern : boost::enable_shared_from_this<SystemMatrixPattern>
{
    // constructor
    SystemMatrixPattern(int type, escript::Distribution_ptr output_distribution,
        escript::Distribution_ptr input_distribution, Pattern_ptr mainPattern,
        Pattern_ptr col_couplePattern, Pattern_ptr row_couplePattern,
        Connector_ptr col_connector, Connector_ptr row_connector);

    ~SystemMatrixPattern() {}

    inline index_t getNumOutput() const {
        return mainPattern->numOutput;
    }

    SystemMatrixPattern_ptr unrollBlocks(int type, dim_t output_block_size,
                                         dim_t input_block_size);

    int type;
    escript::JMPI mpi_info;
    Pattern_ptr mainPattern;
    Pattern_ptr col_couplePattern;
    Pattern_ptr row_couplePattern;
    Connector_ptr col_connector;
    Connector_ptr row_connector;
    escript::Distribution_ptr output_distribution;
    escript::Distribution_ptr input_distribution;
};


} // namespace paso

#endif // __PASO_SYSTEMMATRIXPATTERN_H__

