
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

/* Paso: SystemMatrixPattern */

/****************************************************************************/
 
/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "SystemMatrixPattern.h"

namespace paso {

// constructor for a SystemMatrixPattern
SystemMatrixPattern::SystemMatrixPattern(int patType,
        Paso_Distribution* outDist, Paso_Distribution* inDist,
        Paso_Pattern* mainPat, Paso_Pattern* colPat, Paso_Pattern* rowPat,
        Paso_Connector* colConn, Paso_Connector* rowConn) 
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
    if (mainPat->numOutput != Paso_Distribution_getMyNumComponents(outDist)) {
        Esys_setError(VALUE_ERROR, "SystemMatrixPattern: number of outputs and given distribution don't match.");
    }
    if (mainPat->numInput != Paso_Distribution_getMyNumComponents(inDist)) {
        Esys_setError(VALUE_ERROR, "SystemMatrixPattern: number of input for main pattern and number of send components in connector don't match.");
    }
    if (colPat->numInput != colConn->recv->numSharedComponents) {
        Esys_setError(VALUE_ERROR, "SystemMatrixPattern: number of inputs for column couple pattern and number of received components in connector don't match.");
    }
    if (rowPat->numOutput != rowConn->recv->numSharedComponents) {
        Esys_setError(VALUE_ERROR, "SystemMatrixPattern: number of inputs for row couple pattern and number of received components in connector don't match.");
    }

    type=patType;
    reference_counter=1;
    mainPattern=Paso_Pattern_getReference(mainPat);
    row_couplePattern=Paso_Pattern_getReference(rowPat);
    col_couplePattern=Paso_Pattern_getReference(colPat);
    row_connector=Paso_Connector_getReference(rowConn);
    col_connector=Paso_Connector_getReference(colConn);
    output_distribution=Paso_Distribution_getReference(outDist);
    input_distribution=Paso_Distribution_getReference(inDist);
    mpi_info = Esys_MPIInfo_getReference(outDist->mpi_info);
#ifdef Paso_TRACE
    printf("SystemMatrixPattern: system matrix pattern has been allocated.\n");
#endif
}

// returns a reference to in
SystemMatrixPattern* SystemMatrixPattern_getReference(SystemMatrixPattern* in)
{
    if (in != NULL) {
        ++(in->reference_counter);
    }
    return in;
}
  
// deallocates a SystemMatrixPattern
void SystemMatrixPattern_free(SystemMatrixPattern* in)
{
    if (in != NULL) {
        in->reference_counter--;
        if (in->reference_counter <= 0) {
            Paso_Pattern_free(in->mainPattern);
            Paso_Pattern_free(in->row_couplePattern);
            Paso_Pattern_free(in->col_couplePattern);
            Paso_Connector_free(in->row_connector);
            Paso_Connector_free(in->col_connector);
            Paso_Distribution_free(in->output_distribution);
            Paso_Distribution_free(in->input_distribution);
            Esys_MPIInfo_free(in->mpi_info);
            delete in;
#ifdef Paso_TRACE
            printf("SystemMatrixPattern_free: system matrix pattern has been deallocated.\n");
#endif
        }
    }
}

dim_t SystemMatrixPattern_getNumOutput(const SystemMatrixPattern* in)
{
    if (in != NULL) {
        return in->mainPattern->numOutput;
    }
    return 0;
}

} // namespace paso

