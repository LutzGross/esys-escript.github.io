/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:50 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************   MPI_File_iread_at.c    ************************/
/****************************************************************************/
/* Author : Lisa Alano July 22 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int MPI_File_iread_at(MPI_File fh, MPI_Offset offset, void *buf,
                      int count, MPI_Datatype datatype, 
                      MPIO_Request *request)
{
  _MPI_COVERAGE();
  return PMPI_File_iread_at(fh, offset, buf, count, datatype, request); 
}

