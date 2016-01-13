
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   Paso: coupler                                            */ 

/**************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#ifndef INC_PASO_COUPLER
#define INC_PASO_COUPLER

#include "SharedComponents.h"
/**************************************************************/

typedef struct Paso_Connector {

  Paso_SharedComponents* send;
  Paso_SharedComponents* recv;
  dim_t reference_counter;
  Esys_MPIInfo *mpi_info;

} Paso_Connector;
typedef struct Paso_Coupler {

  dim_t block_size;

  Paso_Connector* connector;

  double *data; /* unmanaged pointer to data be send */
  double *send_buffer;
  double *recv_buffer;
  #ifdef ESYS_MPI
    MPI_Request* mpi_requests;
    MPI_Status* mpi_stati;
 #else
    void* mpi_requests;
    void* mpi_stati;
 #endif
  
  dim_t reference_counter;
  Esys_MPIInfo *mpi_info;

} Paso_Coupler;



PASO_DLL_API
Paso_Connector* Paso_Connector_alloc(Paso_SharedComponents * send, Paso_SharedComponents* recv);

PASO_DLL_API
Paso_Connector* Paso_Connector_getReference(Paso_Connector*);

PASO_DLL_API
Paso_Connector* Paso_Connector_unroll(Paso_Connector* in, index_t block_size);

PASO_DLL_API
Paso_Connector* Paso_Connector_copy(Paso_Connector* in);

PASO_DLL_API
void Paso_Connector_free(Paso_Connector*);


PASO_DLL_API
Paso_Coupler* Paso_Coupler_alloc(Paso_Connector*, dim_t blockSize);

PASO_DLL_API
Paso_Coupler* Paso_Coupler_getReference(Paso_Coupler*);

PASO_DLL_API
void Paso_Coupler_startCollect(Paso_Coupler* self,const double* in);

PASO_DLL_API
double* Paso_Coupler_finishCollect(Paso_Coupler* self);

PASO_DLL_API
void Paso_Coupler_free(Paso_Coupler* in);

#define Paso_Coupler_borrowLocalData(_in_) (_in_)->data
#define Paso_Coupler_borrowRemoteData(_in_) (_in_)->recv_buffer

PASO_DLL_API
void Paso_Coupler_copyAll(const Paso_Coupler* src, Paso_Coupler* target);

PASO_DLL_API
dim_t Paso_Coupler_getLocalLength(const Paso_Coupler* in);
#endif 
