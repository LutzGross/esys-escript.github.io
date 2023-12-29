/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:52 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ********************    PMPI_Group_size.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int PMPI_Group_size ( MPI_Group group, int *size )
{
  if (group == MPI_GROUP_EMPTY)
    *size = 0;
  else
    *size = _MPI_SIZE;
  return MPI_SUCCESS;
}

