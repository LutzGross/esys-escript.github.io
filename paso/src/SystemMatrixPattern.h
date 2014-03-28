
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

PASO_DLL_API
struct SystemMatrixPattern
{
    // constructor
    SystemMatrixPattern(int type, Distribution_ptr output_distribution,
        Distribution_ptr input_distribution, Pattern* mainPattern,
        Pattern* col_couplePattern, Pattern* row_couplePattern,
        Connector* col_connector, Connector* row_connector);

    int type;

    Esys_MPIInfo* mpi_info;

    Pattern* mainPattern;
    Pattern* col_couplePattern;
    Pattern* row_couplePattern;
    Connector* col_connector;
    Connector* row_connector;
    Distribution_ptr output_distribution; 
    Distribution_ptr input_distribution; 

    dim_t reference_counter;
};


/*  interfaces: */

PASO_DLL_API
SystemMatrixPattern* SystemMatrixPattern_getReference(SystemMatrixPattern*);

PASO_DLL_API
void SystemMatrixPattern_free(SystemMatrixPattern*);

SystemMatrixPattern* SystemMatrixPattern_unrollBlocks(
                                           SystemMatrixPattern* pattern,
                                           int type, dim_t output_block_size,
                                           dim_t input_block_size);

index_t SystemMatrixPattern_getNumOutput(const SystemMatrixPattern*);

} // namespace paso

#endif // __PASO_SYSTEMMATRIXPATTERN_H__

