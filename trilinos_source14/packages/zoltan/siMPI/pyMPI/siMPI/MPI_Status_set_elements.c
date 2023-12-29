/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:51 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  **************  MPI_Status_set_elements.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Status_set_elements( MPI_Status *status, MPI_Datatype datatype,
                           int count )
{
  _MPI_COVERAGE();
  return PMPI_Status_set_elements (status, datatype, count);
}

