#ifndef INC_MPI
#define INC_MPI

#include "mpi_C.h"
#include "Common.h"

#define PASO_MPI_TODO 	{ fprintf( stdout, "\nTODO : %s:%d\n", __FILE__, __LINE__);	MPI_Finalize(); exit(1); }

/* Datatypes */
struct Paso_MPIInfo{
  dim_t reference_counter;
  int size;
  int rank;
  MPI_Comm comm;
};

typedef struct Paso_MPIInfo Paso_MPIInfo;

/* Function prototypes */
Paso_MPIInfo* Paso_MPIInfo_alloc( MPI_Comm comm );
void          Paso_MPIInfo_dealloc( Paso_MPIInfo* );
Paso_MPIInfo *Paso_MPIInfo_getReference( Paso_MPIInfo* in );
int           Paso_MPI_initialized( void );

#endif
