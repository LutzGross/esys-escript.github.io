
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/************************************************************************************/

/*   Paso: coupler                                            */ 

/************************************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/************************************************************************************/

#ifndef INC_PASO_COUPLER
#define INC_PASO_COUPLER

#include "SharedComponents.h"
/************************************************************************************/

typedef struct Paso_Connector {

  Paso_SharedComponents* send;
  Paso_SharedComponents* recv;
  dim_t reference_counter;
  Esys_MPIInfo *mpi_info;

} Paso_Connector;
typedef struct Paso_Coupler {

  dim_t block_size;

  Paso_Connector* connector;

  bool in_use;
  double *data; /* unmanaged pointer to data to be sent */
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
#define Paso_Coupler_getNumSharedComponents(_in_)  ((_in_)->connector->send->numSharedComponents)
#define Paso_Coupler_getNumOverlapComponents(_in_) ((_in_)->connector->recv->numSharedComponents)
#define Paso_Coupler_getNumSharedValues(_in_) ( Paso_Coupler_getNumSharedComponents(_in_) * (_in_)->block_size )
#define Paso_Coupler_getNumOverlapValues(_in_) ( Paso_Coupler_getNumOverlapComponents(_in_) * (_in_)->block_size )
#define Paso_Coupler_getLocalLength(_in_) ( (_in_)->connector->send->local_length )

PASO_DLL_API
void Paso_Coupler_copyAll(const Paso_Coupler* src, Paso_Coupler* target);


PASO_DLL_API
void Paso_Coupler_fillOverlap(const dim_t n, double* x, Paso_Coupler *coupler);

PASO_DLL_API
void Paso_Coupler_max(const dim_t n, double* x, Paso_Coupler *coupler);

#endif

