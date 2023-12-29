/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:49 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************    MPI_Alltoall.c        ************************/
/****************************************************************************/
/* Author : Lisa Alano July 17 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int MPI_Alltoall( void *sendbuf, int sendcount, MPI_Datatype sendtype, 
                  void *recvbuf, int recvcnt, MPI_Datatype recvtype, 
                 MPI_Comm comm )
{
  _MPI_COVERAGE();
  return PMPI_Alltoall (sendbuf, sendcount, sendtype, recvbuf, recvcnt, recvtype, comm);
}

