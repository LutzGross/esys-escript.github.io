/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:50 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  **************  MPI_File_read_ordered_begin.c  *********************/
/****************************************************************************/
/* Author : Lisa Alano July 22 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_File_read_ordered_begin(MPI_File fh, void *buf, int count, 
                             MPI_Datatype datatype)
{
  _MPI_COVERAGE();
  return PMPI_File_read_ordered_begin (fh, buf, count, datatype); 
}

