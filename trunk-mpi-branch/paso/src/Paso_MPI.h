#ifndef INC_PASO_MPI
#define INC_PASO_MPI

#include "Common.h"
#include "Paso.h"

#ifdef PASO_MPI
   #include "mpi_C.h"
#else
   typedef int MPI_Comm;
   #define MPI_INT 6
   #define MPI_DOUBLE 11
   #define MPI_COMM_WORLD 91
#endif

#define PASO_MPI_TODO 	{ fprintf( stdout, "\nTODO : %s:%d\n", __FILE__, __LINE__);	MPI_Finalize(); exit(1); }

/* Datatypes */
struct Paso_MPIInfo {
  dim_t reference_counter;
  int size;
  int rank;
  MPI_Comm comm;
  int msg_tag_counter;
};

typedef struct Paso_MPIInfo Paso_MPIInfo;

/* Function prototypes */
Paso_MPIInfo* Paso_MPIInfo_alloc( MPI_Comm comm );
void          Paso_MPIInfo_dealloc( Paso_MPIInfo* );
Paso_MPIInfo *Paso_MPIInfo_getReference( Paso_MPIInfo* in );
int           Paso_MPI_initialized( void );
bool_t Paso_MPI_noError( Paso_MPIInfo *mpi_info);
index_t PASO_MPI_mod(index_t k, index_t n);

#endif /* INC_PASO_MPI */
