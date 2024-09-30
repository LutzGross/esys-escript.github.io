/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:50 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  *********************  MPI_Graph_create.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Graph_create ( MPI_Comm comm_old, int nnodes, int *index, int *edges, 
                      int reorder, MPI_Comm *comm_graph )
{
  _MPI_COVERAGE();
  return PMPI_Graph_create (comm_old, nnodes, index, edges, reorder, comm_graph);
}

