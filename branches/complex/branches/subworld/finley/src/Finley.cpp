
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

#include "Finley.h"
#include "esysUtils/error.h"
#include "finley/CppAdapter/FinleyAdapterException.h" // temporary

namespace finley {

/// returns a time mark
double timer()
{
    return Esys_timer();
}

/// checks if the pointer ptr has a target. If not an error is raised and
/// TRUE is returned.
bool checkPtr(void* arg)
{
    return Esys_checkPtr(arg);
}

/// resets the error to NO_ERROR
void resetError()
{
    Esys_resetError();
}

/// sets an error
void setError(ErrorCodeType err, const char* msg)
{
    Esys_setError(err,msg);
}

/// checks if there is no error
bool noError()
{
    return Esys_noError();
}

/// returns the error code
ErrorCodeType getErrorType()
{
    return Esys_getErrorType();
}

/// returns the error message
char* getErrorMessage(void)
{
    return Esys_getErrorMessage();
}

void checkFinleyError() 
{
    if (!noError()) {
        // reset the error code to no error otherwise the next call to
        // this function may resurrect a previous error
        resetError();
        throw FinleyAdapterException(getErrorMessage());
    }
}

/* checks that there is no error across all processes in a communicator */
/* NOTE : does not guarantee consistency of error string on each process */
bool MPI_noError(esysUtils::JMPI& mpi_info)
{
    return esysUtils::Esys_MPIInfo_noError(mpi_info);
}

} // namespace finley

