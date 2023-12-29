/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:51 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************  MPI_Reduce_scatter.c    ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Reduce_scatter ( void *sendbuf, void *recvbuf, int *recvcnts, 
                       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm )
{
  _MPI_COVERAGE();
  return PMPI_Reduce_scatter (sendbuf, recvbuf, recvcnts, datatype, op, comm);
}

