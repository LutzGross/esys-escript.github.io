
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

/*   Paso: system matrix pattern                              */

/****************************************************************************/

/*   Copyrights by ACcESS Australia 2004,2005 */
/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_SYSTEMMATRIXPATTERN_H__
#define __PASO_SYSTEMMATRIXPATTERN_H__

#include "Distribution.h"
#include "Pattern.h"
#include "Coupler.h"

namespace paso {

struct SystemMatrixPattern;
typedef boost::shared_ptr<SystemMatrixPattern> SystemMatrixPattern_ptr;
typedef boost::shared_ptr<const SystemMatrixPattern> const_SystemMatrixPattern_ptr;

PASO_DLL_API
struct SystemMatrixPattern : boost::enable_shared_from_this<SystemMatrixPattern>
{
    // constructor
    SystemMatrixPattern(int type, Distribution_ptr output_distribution,
        Distribution_ptr input_distribution, Pattern* mainPattern,
        Pattern* col_couplePattern, Pattern* row_couplePattern,
        Connector_ptr col_connector, Connector_ptr row_connector);

    ~SystemMatrixPattern()
    {
        Pattern_free(mainPattern);
        Pattern_free(row_couplePattern);
        Pattern_free(col_couplePattern);
        Esys_MPIInfo_free(mpi_info);
    }

    inline index_t getNumOutput() const {
        return mainPattern->numOutput;
    }

    SystemMatrixPattern_ptr unrollBlocks(int type, dim_t output_block_size,
                                         dim_t input_block_size);

    int type;
    Esys_MPIInfo* mpi_info;
    Pattern* mainPattern;
    Pattern* col_couplePattern;
    Pattern* row_couplePattern;
    Connector_ptr col_connector;
    Connector_ptr row_connector;
    Distribution_ptr output_distribution; 
    Distribution_ptr input_distribution; 
};


} // namespace paso

#endif // __PASO_SYSTEMMATRIXPATTERN_H__

