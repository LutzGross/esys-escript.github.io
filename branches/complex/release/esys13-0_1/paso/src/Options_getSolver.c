/* $Id$ */


/*
********************************************************************************
*               Copyright   2006 by ACcESS MNRF                                *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

/**************************************************************/

/* Paso: returns the package to be used                   */

/**************************************************************/

/* Copyrights by ACcESS Australia 2004,2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "Options.h"

/**************************************************************/

index_t Paso_Options_getSolver(index_t solver,index_t package, bool_t symmetry) {
  index_t out=PASO_DEFAULT;
  /* PASO */
  if (package==PASO_PASO) {
     switch (solver) {
        case PASO_BICGSTAB:
            out=PASO_BICGSTAB;
            break;
        case PASO_PCG:
            out=PASO_PCG;
            break;
        case PASO_PRES20:
            out=PASO_PRES20;
            break;
        case PASO_GMRES:
            out=PASO_GMRES;
            break;
        default:
            if (symmetry) {
               out=PASO_PCG;
            } else {
               out=PASO_BICGSTAB;
            }
            break;
     }
  /* SCSL */
  } else if (package==PASO_SCSL) {
    switch (solver) {
      case PASO_PCG:
        out=PASO_PCG;
        break;
      case PASO_CR:
        out=PASO_CR;
        break;
      case PASO_CGS:
        out=PASO_CGS;
        break;
      case PASO_BICGSTAB:
        out=PASO_BICGSTAB;
        break;
      case PASO_ITERATIVE:
        if (symmetry) {
          out=PASO_PCG;
        } else {
          out=PASO_BICGSTAB;
        }
        break;
      case PASO_CHOLEVSKY:
        out=PASO_CHOLEVSKY;
        break;
      case PASO_DIRECT:
        out=PASO_DIRECT;
        break;
      default:
        if (symmetry) {
          out=PASO_CHOLEVSKY;
        } else {
          out=PASO_DIRECT;
        }
        break;
     }
  /* MKL */
  } else if (package==PASO_MKL) {
    switch (solver) {
      case PASO_CHOLEVSKY:
        out=PASO_CHOLEVSKY;
        break;
      case PASO_DIRECT:
        out=PASO_DIRECT;
        break;
      default:
        if (symmetry) {
          out=PASO_CHOLEVSKY;
        } else {
          out=PASO_DIRECT;
        }
        break;
    }
  /* UMFPACK */
  } else if (package==PASO_UMFPACK) {
      out=PASO_DIRECT;
  } else {
      Paso_setError(VALUE_ERROR,"Unidentified package.");
  }
  return out;
}

/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:38  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.1  2005/09/05 06:29:46  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
