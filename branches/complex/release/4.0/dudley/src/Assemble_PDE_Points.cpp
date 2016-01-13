
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


/************************************************************************************/

/*    assembles the system of numEq PDEs into the stiffness matrix S right hand side F  */
/*    the shape functions for test and solution must be identical */


/*      d_dirac_{k,m} u_m yand _dirac_k */

/*    u has p.numComp components in a 3D domain. The shape functions for test and solution must be identical  */
/*    and row_NS == row_NN                                                                                  */

/*    Shape of the coefficients: */

/*      d_dirac = p.numEqu x p.numComp  */
/*      y_dirac = p.numEqu   */


/************************************************************************************/


#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif


/************************************************************************************/

void  Dudley_Assemble_PDE_Points(Dudley_Assemble_Parameters p,
                                 Dudley_ElementFile* elements,
                                 paso::SystemMatrix_ptr Mat, escriptDataC* F,
                                 escriptDataC* d_dirac, escriptDataC* y_dirac) {

    index_t color, e, row_index;
    __const double  *d_dirac_p, *y_dirac_p;
    
    double *F_p=(requireWrite(F), getSampleDataRW(F,0));	/* use comma, to get around the mixed code and declarations thing */

    #pragma omp parallel private(color, d_dirac_p, y_dirac_p)
    {
          for (color=elements->minColor;color<=elements->maxColor;color++) {
             /*  open loop over all elements: */
             #pragma omp for private(e) schedule(static)
             for(e=0;e<elements->numElements;e++){
                if (elements->Color[e]==color) {
                   
		   d_dirac_p=getSampleDataRO(d_dirac, e);
                   y_dirac_p=getSampleDataRO(y_dirac, e);
		   
                   row_index=p.row_DOF[elements->Nodes[INDEX2(0,e,p.NN)]];
		   
		   if (NULL!=y_dirac_p)  Dudley_Util_AddScatter(1,
                                                        &row_index,
                                                        p.numEqu,
                                                        y_dirac_p,
                                                        F_p, 
                                                        p.row_DOF_UpperBound);
		   
                   if (NULL!=d_dirac_p) Dudley_Assemble_addToSystemMatrix(Mat,
                                                                   1,
                                                                   &row_index,
                                                                   p.numEqu,
                                                                   1,
                                                                   &row_index,
                                                                   p.numComp,
                                                                   d_dirac_p);
                } /* end color check */
             } /* end element loop */
         } /* end color loop */
   } /* end parallel region */
}
