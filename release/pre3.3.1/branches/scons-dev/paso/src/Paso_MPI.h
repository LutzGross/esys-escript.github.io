
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#ifndef INC_PASO_MPI
#define INC_PASO_MPI

#include "Common.h"
#include "Paso.h"

#ifdef PASO_MPI
   #include "mpi_C.h"
#else
   typedef int MPI_Comm;
   typedef int MPI_Request;
   #define MPI_INT 6
   #define MPI_DOUBLE 11
   #define MPI_COMM_WORLD 91
#endif

typedef int Paso_MPI_rank;

#define PASO_MPI_TODO 	{ fprintf( stdout, "\nTODO : %s:%d\n", __FILE__, __LINE__);	MPI_Finalize(); exit(1); }

/* Datatypes */
struct Paso_MPIInfo {
  dim_t reference_counter;
  int size;
  Paso_MPI_rank rank;
  MPI_Comm comm;
  int msg_tag_counter;
};

typedef struct Paso_MPIInfo Paso_MPIInfo;

/* Function prototypes */
Paso_MPIInfo* Paso_MPIInfo_alloc( MPI_Comm comm );
void          Paso_MPIInfo_free( Paso_MPIInfo* );
Paso_MPIInfo *Paso_MPIInfo_getReference( Paso_MPIInfo* in );
int           Paso_MPIInfo_initialized( void );
index_t Paso_MPIInfo_mod(index_t n, index_t k);
dim_t Paso_MPIInfo_setDistribution(Paso_MPIInfo* in ,index_t min_id,index_t max_id,index_t* distribution);
void Paso_MPIInfo_Split( Paso_MPIInfo *mpi_info, dim_t n, dim_t* local_N,index_t* offset); 
bool_t Paso_MPIInfo_noError( Paso_MPIInfo *mpi_info);
char *Paso_MPI_appendRankToFileName(const char *, int, int);

#endif /* INC_PASO_MPI */
