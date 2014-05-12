
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

/*    Paso finite element solver library                                    */

/****************************************************************************/

/*  Copyrights by ACcESS Australia, 2003,2004,2005 */
/*  Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_H__
#define __PASO_H__

#include <cfloat>
#include <esysUtils/error.h>
#include <esysUtils/Esys_MPI.h>
#include <esysUtils/index.h>
#include <esysUtils/maths.h>

#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>

#define PASO_DLL_API
#ifdef _WIN32
#   ifndef PASO_STATIC_LIB
#      undef PASO_DLL_API
#      ifdef PASO_EXPORTS
#         define PASO_DLL_API __declspec(dllexport)
#      else
#         define PASO_DLL_API __declspec(dllimport)
#      endif
#   endif
#endif

#define MATRIX_FORMAT_DEFAULT 1
#define MATRIX_FORMAT_CSC 2
#define MATRIX_FORMAT_BLK1 4
#define MATRIX_FORMAT_OFFSET1 8
#define MATRIX_FORMAT_TRILINOS_CRS 16
#define MATRIX_FORMAT_DIAGONAL_BLOCK 32

#define PASO_ONE (double)(1.0)
#define PASO_ZERO (double)(0.0)

#endif // __PASO_H__

