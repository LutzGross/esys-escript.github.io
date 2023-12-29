/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2006/04/25 19:03:07 $
 *    Revision: 1.1 $
 ****************************************************************************/
/* FILE  ******************      MPI_Bcast.c        ************************/
/****************************************************************************/
/* Author : Lisa Alano July 1 2002                                          */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/*==========================================================================*/
int MPI_Bcast ( void* buffer, int count, MPI_Datatype datatype, int root, 
  MPI_Comm comm )
{
  _MPI_COVERAGE();
  return PMPI_Bcast(buffer, count, datatype, root, comm);
}

