
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


/****************************************************************************

  Assembles the system of numEq PDEs into the stiffness matrix S and right
  hand side F

      d_dirac_{k,m} u_m and y_dirac_k

  u has p.numComp components in a 3D domain.
  The shape functions for test and solution must be identical and
  row_NS == row_NN.

  Shape of the coefficients:

      d_dirac = p.numEqu x p.numComp
      y_dirac = p.numEqu

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

namespace finley {

void Assemble_PDE_Points(const AssembleParameters& p,
                         const escript::Data& d_dirac,
                         const escript::Data& y_dirac)
{

    double *F_p=NULL;
    if(!p.F.isEmpty()) {
        p.F.requireWrite();
        F_p=p.F.getSampleDataRW(0);
    }
#pragma omp parallel
    {
        for (int color=p.elements->minColor; color<=p.elements->maxColor; color++) {
            // loop over all elements
#pragma omp for
            for (int e=0; e<p.elements->numElements; e++) {
                if (p.elements->Color[e]==color) {
                    int row_index=p.row_DOF[p.elements->Nodes[INDEX2(0,e,p.NN)]];
                    if (!y_dirac.isEmpty()) {
                        const double *y_dirac_p=y_dirac.getSampleDataRO(e);
                        util::addScatter(1, &row_index, p.numEqu,
                                         y_dirac_p, F_p, p.row_DOF_UpperBound);
                    }
                    if (!d_dirac.isEmpty()) {
                        const double *d_dirac_p=d_dirac.getSampleDataRO(e);
                        Assemble_addToSystemMatrix(p.S, 1, &row_index,
                                p.numEqu, 1, &row_index, p.numComp, d_dirac_p);
                    }
                } // end color check
            } // end element loop
        } // end color loop
    } // end parallel section
}

} // namespace finley

