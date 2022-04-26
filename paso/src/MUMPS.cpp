
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
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


/****************************************************************************/

/* Paso: interface to the MUMPS library */

/****************************************************************************/

#include "Paso.h"
#include "MUMPS.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>

namespace paso {

std::ostream& operator<<(std::ostream& os, const cplx_t& c)
{
    os << "(" << c.real() << ", " << c.imag() << ")";
    return os;
}

} // namespace paso
