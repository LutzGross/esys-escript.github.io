
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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
    std::stringstream ss;

    if (outDist->mpi_info != inDist->mpi_info) {
        ss << "SystemMatrixPattern: output distribution and input distribution MPI communicators don't match.";
    } else if (outDist->mpi_info != colConn->mpi_info) {
        ss << "SystemMatrixPattern: output distribution and col connector MPI communicators don't match.";
    } else if (outDist->mpi_info != rowConn->mpi_info) {
        ss << "SystemMatrixPattern: output distribution and row connector MPI communicators don't match.";
    } else if (mainPat->type != patType)  {
        ss << "SystemMatrixPattern: type of mainPattern (" << mainPat->type
           << ") does not match expected type (" << patType << ")";
    } else if (colPat->type != patType)  {
        ss << "SystemMatrixPattern: type of col couplePattern (" << colPat->type
           << ") does not match expected type (" << patType << ")";
    } else if (rowPat->type != patType)  {
        ss << "SystemMatrixPattern: type of row couplePattern (" << rowPat->type
           << ") does not match expected type (" << patType << ")";
    } else if (colPat->numOutput != mainPat->numOutput) {
        ss << "SystemMatrixPattern: number of outputs for couple and main "
              "pattern don't match: " << colPat->numOutput << " != "
           << mainPat->numOutput;
    } else if (mainPat->numOutput != outDist->getMyNumComponents()) {
        ss << "SystemMatrixPattern: number of outputs and given distribution "
              "don't match: " << mainPat->numOutput << " != "
           << outDist->getMyNumComponents();
    } else if (mainPat->numInput != inDist->getMyNumComponents()) {
        ss << "SystemMatrixPattern: number of input for main pattern and "
              "number of send components in connector don't match: "
           << mainPat->numInput << " != " << inDist->getMyNumComponents();
    } else if (colPat->numInput != colConn->recv->numSharedComponents) {
        ss << "SystemMatrixPattern: number of inputs for column couple pattern"
              " and number of received components in connector don't match: "
           << colPat->numInput << " != " << colConn->recv->numSharedComponents;
    } else if (rowPat->numOutput != rowConn->recv->numSharedComponents) {
        ss << "SystemMatrixPattern: number of inputs for row couple pattern "
              "and number of received components in connector don't match: "
           << rowPat->numOutput << " != " << rowConn->recv->numSharedComponents;
    }
    const std::string msg(ss.str());
    int error = msg.length(); // proxy for error condition
    int gerror = error;
    escript::checkResult(error, gerror, outDist->mpi_info);
    if (gerror > 0) {
        char* gmsg;
        escript::shipString(msg.c_str(), &gmsg, outDist->mpi_info->comm);
        std::string realmsg(error? msg : gmsg);
        //throw PasoException(realmsg);
    }
    std::cerr << msg << std::endl;

    mpi_info = outDist->mpi_info;
}

} // namespace paso


