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



extern "C" {
#include "TRILINOS.h"
}

#ifdef TRILINOS
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
#endif



/* This C++ class is accessed from a C program via the extern "C" functions that are included below */



/**
  \brief
  Defines the class that stores TRILINOS map and epetra matrix

  Description:
  This file is only necessary because we have to use C++ objects from C
*/

class ESCRIPT_DLL_API TrilinosData {

    Epetra_Map *epetra_map;
    Epetra_CrsMatrix *epetra_crs_matrix;

  public:

  TrilinosData(Paso_SystemMatrixPattern *pattern,dim_t row_block_size,dim_t col_block_size) {
    #ifdef TRILINOS
       Paso_Distribution *row_distribution=pattern->output_distribution;
       Paso_MPIInfo *mpi_info=row_distribution->mpi_info;
       Epetra_MpiComm Comm(mpi_info->comm);
       Epetra_Map temp_map(row_distribution->numComponents,row_distribution->myNumComponents, 0,Comm);
       epetra_map = &temp_map;
       int *numNz = new int[row_distribution->myNumComponents];
       for (int i=0; i<row_distribution->myNumComponents; i++) numNz[i] = pattern->ptr[i+1] - pattern->ptr[i];
       Epetra_CrsMatrix temp(Copy, *epetra_map, numNz, true);
       epetra_crs_matrix = &temp;
    #else
       Paso_setError(SYSTEM_ERROR,"Paso_TRILINOS: TRILINOS is not available.");
    #endif
  }

  ~TrilinosData() {
      delete &epetra_map;
      delete &epetra_crs_matrix;
  }

  void SumIntoMyValues(int row, int num, double *value, int *col) {
    printf("ksteube TrilinosData::SumIntoMyValues num=%d value=%f col=%d\n", num, *value, *col);
    #ifdef TRILINOS
       epetra_crs_matrix->SumIntoMyValues(row, num, value, col);
    #else
       Paso_setError(SYSTEM_ERROR,"Paso_TRILINOS: TRILINOS is not available.");
    #endif
  }

  protected:

  private:

};


/* Below here are the methods we use to access the class above from a C program */



extern "C"
void Paso_TRILINOS_alloc(void* trilinos_data, Paso_SystemMatrixPattern *pattern,dim_t row_block_size,dim_t col_block_size) {
    printf("ksteube in Initialize_TrilinosData\n");
    trilinos_data = new TrilinosData(pattern,row_block_size,col_block_size);
    printf("ksteube Paso_SystemMatrixPattern len=%d\n", pattern->myLen);
}

extern "C"
void Trilinos_SumIntoMyValues(void *p, int row, int col, double *value) {
    printf("ksteube in Trilinos_SumIntoMyValues row=%d col=%d value=%lf\n", row, col, value);
    ((TrilinosData *)p)->SumIntoMyValues(row, 1, value, &col);
}


extern "C"
void Paso_TRILINOS(Paso_SystemMatrix* A,
                   double* out,
                   double* in,
                   Paso_Options* options,
                   Paso_Performance* pp) {
#ifdef TRILINOS
printf("no solver called.\n");
#else
    Paso_setError(SYSTEM_ERROR,"Paso_TRILINOS: TRILINOS is not available.");
#endif
}

extern "C"
void Paso_TRILINOS_free(void* in) {
   ((TrilinosData *)in)->~TrilinosData();
}
