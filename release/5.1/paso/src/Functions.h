
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


#ifndef __PASO_FUNCTIONS_H__
#define __PASO_FUNCTIONS_H__

#include "Paso.h"
#include "performance.h"
#include "SystemMatrix.h"

namespace paso {

struct Function
{
    Function(const escript::JMPI& mpi_info);
    virtual ~Function();

    /// sets value=F(arg)
    virtual SolverResult call(double* value, const double* arg, Performance* pp) = 0;

    /// numerical calculation of the directional derivative J0w of F at x0 in
    /// the direction w. f0 is the value of F at x0. setoff is workspace
    SolverResult derivative(double* J0w, const double* w, const double* f0,
                            const double* x0, double* setoff, Performance* pp);

    /// returns the length of the vectors used by this function
    virtual dim_t getLen() = 0;

    const escript::JMPI mpi_info;
};

struct LinearSystem : public Function
{
    LinearSystem(SystemMatrix_ptr A, double* b, Options* options);
    virtual ~LinearSystem();

    virtual SolverResult call(double* value, const double* arg, Performance* pp);

    virtual dim_t getLen() { return n; }

    SystemMatrix_ptr mat;
    double* tmp;
    double* b;
    dim_t n;
};

} // namespace paso

#endif // __PASO_FUNCTIONS_H__

