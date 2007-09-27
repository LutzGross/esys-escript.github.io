
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

/**************************************************************/

/*   Paso: coupler                                            */ 

/**************************************************************/

/*   Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_COUPLER
#define INC_PASO_COUPLER

#include "SharedComponents.h"
/**************************************************************/

typedef struct Paso_Coupler {

  dim_t block_size;

  Paso_SharedComponents* send;
  double *send_buffer;
  Paso_SharedComponents* recv;
  double *recv_buffer;
  #ifdef PASO_MPI
    MPI_Request* mpi_requests;
    MPI_Status* mpi_stati;
 #else
    void* mpi_requests;
    void* mpi_stati;
 #endif
  
  dim_t reference_counter;
  Paso_MPIInfo *mpi_info;

} Paso_Coupler;


Paso_Coupler* Paso_Coupler_alloc(Paso_SharedComponents * send,
                                 Paso_SharedComponents* recv);
Paso_Coupler* Paso_Coupler_getReference(Paso_Coupler*);
void Paso_Coupler_free(Paso_Coupler*);
void Paso_Coupler_startCollect(Paso_Coupler* self, double* in);
double* Paso_Coupler_finishCollect(Paso_Coupler* self);
void Paso_Coupler_freeBuffer(Paso_Coupler* coupler);
void Paso_Coupler_allocBuffer(Paso_Coupler* coupler,dim_t blockSize);
Paso_Coupler* Paso_Coupler_unroll(Paso_Coupler* in, index_t block_size);


#endif 
