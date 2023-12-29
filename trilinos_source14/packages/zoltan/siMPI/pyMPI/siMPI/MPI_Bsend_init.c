/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:49 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************     MPI_Bsend_init.c     ************************/
/****************************************************************************/
/* Author : Lisa Alano July 18 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int MPI_Bsend_init( void *buf, int count, MPI_Datatype datatype, int dest,
                   int tag, MPI_Comm comm, MPI_Request *request )
{
  _MPI_COVERAGE();
  return PMPI_Bsend_init (buf, count, datatype, dest, tag, comm, request); 
}

