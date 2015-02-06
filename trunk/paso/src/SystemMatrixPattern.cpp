
/*****************************************************************************
*
* Copyright (c) 2003-2015 by University of Queensland
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

/* Paso: SystemMatrixPattern */

/****************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "SystemMatrixPattern.h"

namespace paso {

SystemMatrixPattern::SystemMatrixPattern(int patType, Distribution_ptr outDist,
        Distribution_ptr inDist, Pattern_ptr mainPat, Pattern_ptr colPat,
        Pattern_ptr rowPat, Connector_ptr colConn, Connector_ptr rowConn) :
    type(patType),
    mainPattern(mainPat),
    col_couplePattern(colPat),
    row_couplePattern(rowPat),
    col_connector(colConn),
    row_connector(rowConn),
    output_distribution(outDist),
    input_distribution(inDist)
{
    Esys_resetError();

    if (outDist->mpi_info != inDist->mpi_info) {
        Esys_setError(SYSTEM_ERROR, "SystemMatrixPattern: output distribution and input distribution MPI communicators don't match.");
    }
    if (outDist->mpi_info != colConn->mpi_info) {
        Esys_setError(SYSTEM_ERROR, "SystemMatrixPattern: output distribution and col connector MPI communicators don't match.");
    }
    if (outDist->mpi_info != rowConn->mpi_info ) {
        Esys_setError(SYSTEM_ERROR, "SystemMatrixPattern: output distribution and row connector MPI communicators don't match.");
    }
    if (mainPat->type != patType)  {
        Esys_setError(VALUE_ERROR, "SystemMatrixPattern: type of mainPattern does not match expected type.");
    }
    if (colPat->type != patType)  {
        Esys_setError(VALUE_ERROR, "SystemMatrixPattern: type of col couplePattern does not match expected type.");
    }
    if (rowPat->type != patType)  {
        Esys_setError(VALUE_ERROR, "SystemMatrixPattern: type of row couplePattern does not match expected type.");
    }
    if (colPat->numOutput != mainPat->numOutput) {
        Esys_setError(VALUE_ERROR, "SystemMatrixPattern: number of outputs for couple and main pattern don't match.");
    }
    if (mainPat->numOutput != outDist->getMyNumComponents()) {
        Esys_setError(VALUE_ERROR, "SystemMatrixPattern: number of outputs and given distribution don't match.");
    }
    if (mainPat->numInput != inDist->getMyNumComponents()) {
        Esys_setError(VALUE_ERROR, "SystemMatrixPattern: number of input for main pattern and number of send components in connector don't match.");
    }
    if (colPat->numInput != colConn->recv->numSharedComponents) {
        Esys_setError(VALUE_ERROR, "SystemMatrixPattern: number of inputs for column couple pattern and number of received components in connector don't match.");
    }
    if (rowPat->numOutput != rowConn->recv->numSharedComponents) {
        Esys_setError(VALUE_ERROR, "SystemMatrixPattern: number of inputs for row couple pattern and number of received components in connector don't match.");
    }
    mpi_info = outDist->mpi_info;

}

} // namespace paso


