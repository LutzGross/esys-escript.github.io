#ifndef INC_MPI
#define INC_MPI

#include "Common.h"

#ifdef PASO_MPI

#include "mpi_C.h"

#define PASO_MPI_TODO 	{ fprintf( stdout, "\nTODO : %s:%d\n", __FILE__, __LINE__);	MPI_Finalize(); exit(1); }
#define PASO_INFO_ERRORMSG( err, msg ) { char _msg__[256]; sprintf( _msg__, "%s : %s:%d\n", msg, __FILE__, __LINE__ ); Paso_setError( err, _msg__ );  }

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

#endif // INC_PASO_MPI

#endif // PASO_MPI
