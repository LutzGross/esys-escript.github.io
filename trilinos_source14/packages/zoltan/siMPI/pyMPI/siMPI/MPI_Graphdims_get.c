/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:50 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  *******************  MPI_Graphdims_get.c    ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Graphdims_get ( MPI_Comm comm, int *nnodes, int *nedges )
{
  _MPI_COVERAGE();
  return PMPI_Graphdims_get (comm, nnodes, nedges);
}

