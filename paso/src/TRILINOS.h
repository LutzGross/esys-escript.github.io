/* $Id$ */

/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

/* Interface to Sandia TRILINOS sparse solver */

/* Author: k.steube@uq.edu.au */


#ifndef INC_PASO_TRILINOS
#define INC_PASO_TRILINOS

void Initialize_TrilinosData(void* p, Paso_SystemMatrixPattern *pattern);

void Paso_TRILINOS(Paso_SystemMatrix* A,
                          double* out,
                          double* in,
                          Paso_Options* options,
                          Paso_Performance* pp);

void Paso_TRILINOS_free(double* in);

#endif /* ifndef INC_PASO_TRILINOS */

