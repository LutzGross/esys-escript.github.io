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

/* TODO
	New TrilinosData object
	Create Map
	Begin insert global values
	Insert values
	End submit entries
	Fill complete
	Clean up
	Solve system given a RHS
	Get result
*/


#ifdef TRILINOS

#include "Paso.h"
#include "performance.h"
#include "escript/system_dep.h"

#include "Epetra_ConfigDefs.h"
#ifdef PASO_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_IntSerialDenseVector.h"



/* This C++ class is accessed from a C program via the extern "C" functions that are included below */



/**
  \brief
  Defines the class that stores TRILINOS map and epetra matrix

  Description:
  This file is only necessary because we have to use C++ objects from C
*/

class ESCRIPT_DLL_API TrilinosData {

    Paso_MPIInfo *MPIInfo;
    Epetra_Map *epetra_map;
    Epetra_CrsMatrix *epetra_crs_matrix;

  public:

  TrilinosData(Paso_MPIInfo *info, Paso_SystemMatrixPattern *pattern) {
    MPIInfo = info;
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    Epetra_Map temp_map(row_distribution->numDOF * row_blocksize, row_distribution->numLocal * row_blocksize, 0, Comm);
    epetra_map = &temp_map;
    int *numNz = new int[row_distribution->numDOF]; // probably need numLocal
    for (i=0; i<row_distribution->numLocal; i++) {
      for (ii=0; i<row_blocksize; i++) {
	numNz[ii+row_blocksize*i] = (pattern->ptr[i+1] - pattern->ptr[i]) * col_blocksize
      }
    }
    Epetra_CrsMatrix temp(Copy, *epetra_map, numNz);
    epetra_crs_matrix = &temp;
  }

  void SumIntoMyValues(int row, int num, double *value, int *col) {
    printf("ksteube TrilinosData::SumIntoMyValues num=%d value=%f col=%d\n", num, *value, *col);
    epetra_crs_matrix->SumIntoMyValues(row, num, value, col);
  }

  protected:

  private:

};



/* Below here are the methods we use to access the class above from a C program */



extern "C"
void Initialize_TrilinosData(TrilinosData *p, Paso_SystemMatrixPattern *pattern) {
    printf("ksteube in Initialize_TrilinosData\n");
    p = new TrilinosData(pattern->MPIInfo, pattern);
    printf("ksteube Paso_SystemMatrixPattern len=%d\n", pattern->len);
}

extern "C"
void Trilinos_SumIntoMyValues(TrilinosData *p, int row, int col, double *value) {
    printf("ksteube in Trilinos_SumIntoMyValues row=%d col=%d value=%lf\n", row, col, value);
    p->SumIntoMyValues(row, 1, &value, &col);
}

extern "C"
void Paso_TRILINOS(Paso_SystemMatrix* A,
                          double* out,
                          double* in,
                          Paso_Options* options,
                          Paso_Performance* pp) {
#ifdef TRILINOS
#else
    Paso_setError(SYSTEM_ERROR,"Paso_TRILINOS: TRILINOS is not available.");
#endif
}

extern "C"
void Paso_TRILINOS_free(double* in) {
#ifdef TRILINOS
#endif
}

#endif /* ifdef TRILINOS */

