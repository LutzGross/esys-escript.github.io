
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


/****************************************************************************

  Assembles the system of numEqu PDEs into the stiffness matrix S and right
  hand side F

      d_dirac_{k,m} u_m and y_dirac_k

  u has p.numEqu components in a 3D domain.
  The shape functions for test and solution must be identical and
  row_NS == row_NN.

  Shape of the coefficients:

      d_dirac = p.numEqu x p.numEqu
      y_dirac = p.numEqu


*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

#include <escript/index.h>

namespace dudley {

template<typename Scalar>
void Assemble_PDE_Points(const AssembleParameters& p,
                         const escript::Data& d_dirac,
                         const escript::Data& y_dirac)
{
    Scalar* F_p = NULL;
    const Scalar zero = static_cast<Scalar>(0);
    if (!p.F.isEmpty()) {
        p.F.requireWrite();
        F_p = p.F.getSampleDataRW(0, zero);
    }

#pragma omp parallel
    {
        std::vector<index_t> rowIndex(1);
        std::vector<Scalar> values(p.numEqu*p.numEqu);

        for (index_t color = p.elements->minColor; color <= p.elements->maxColor; color++) {
            // loop over all elements
#pragma omp for
            for (index_t e = 0; e < p.elements->numElements; e++) {
                if (p.elements->Color[e] == color) {
                    rowIndex[0] = p.DOF[p.elements->Nodes[INDEX2(0,e,p.NN)]];
                    if (!y_dirac.isEmpty()) {
                        const Scalar* y_dirac_p = y_dirac.getSampleDataRO(e, zero);
                        util::addScatter(1, &rowIndex[0], p.numEqu,
                                         y_dirac_p, F_p, p.DOF_UpperBound);
                    }
                   
                    if (!d_dirac.isEmpty()) {
                        const Scalar* EM_S = d_dirac.getSampleDataRO(e, zero);
                        values.assign(EM_S, EM_S+p.numEqu*p.numEqu);
                        Assemble_addToSystemMatrix(p.S, rowIndex, p.numEqu,
                                                   values);
                    }
                } // end color check
            } // end element loop
        } // end color loop
    } // end parallel region
}

// instantiate our two supported versions
template void Assemble_PDE_Points<escript::DataTypes::real_t>(
                            const AssembleParameters& p,
                            const escript::Data& d, const escript::Data& y);
template void Assemble_PDE_Points<escript::DataTypes::cplx_t>(
                            const AssembleParameters& p,
                            const escript::Data& d, const escript::Data& y);

} // namespace dudley

