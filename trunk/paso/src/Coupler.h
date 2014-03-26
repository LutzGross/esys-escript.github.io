
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


/****************************************************************************/

/*   Paso: coupler                                            */ 

/****************************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_COUPLER_H__
#define __PASO_COUPLER_H__

#include "SharedComponents.h"

namespace paso {

struct Connector
{
    Paso_SharedComponents* send;
    Paso_SharedComponents* recv;
    dim_t reference_counter;
    Esys_MPIInfo* mpi_info;
};

struct Coupler
{
    Connector* connector;

    dim_t block_size;
    bool in_use;

    // unmanaged pointer to data to be sent
    double *data;
    double *send_buffer;
    double *recv_buffer;
    MPI_Request* mpi_requests;
    MPI_Status* mpi_stati;
  
    dim_t reference_counter;
    Esys_MPIInfo *mpi_info;
};


PASO_DLL_API
Connector* Connector_alloc(Paso_SharedComponents* send, Paso_SharedComponents* recv);

PASO_DLL_API
Connector* Connector_getReference(Connector*);

PASO_DLL_API
Connector* Connector_unroll(Connector* in, index_t block_size);

PASO_DLL_API
Connector* Connector_copy(Connector* in);

PASO_DLL_API
void Connector_free(Connector*);


PASO_DLL_API
Coupler* Coupler_alloc(Connector*, dim_t blockSize);

PASO_DLL_API
Coupler* Coupler_getReference(Coupler*);

PASO_DLL_API
void Coupler_startCollect(Coupler*, const double* in);

PASO_DLL_API
double* Coupler_finishCollect(Coupler*);

PASO_DLL_API
void Coupler_free(Coupler*);

PASO_DLL_API
void Coupler_copyAll(const Coupler* src, Coupler* target);

PASO_DLL_API
void Coupler_fillOverlap(dim_t n, double* x, Coupler* coupler);

PASO_DLL_API
void Coupler_max(dim_t n, double* x, Coupler* coupler);

inline const double* Coupler_borrowLocalData(const Coupler* in)
{
    return in->data;
}

inline const double* Coupler_borrowRemoteData(const Coupler* in)
{
    return in->recv_buffer;
}

inline dim_t Coupler_getNumSharedComponents(const Coupler* in)
{
    return in->connector->send->numSharedComponents;
}

inline dim_t Coupler_getNumOverlapComponents(const Coupler* in)
{
    return in->connector->recv->numSharedComponents;
}

inline dim_t Coupler_getNumSharedValues(const Coupler* in)
{
    return Coupler_getNumSharedComponents(in) * in->block_size;
}

inline dim_t Coupler_getNumOverlapValues(const Coupler* in)
{
    return Coupler_getNumOverlapComponents(in) * in->block_size;
}

inline dim_t Coupler_getLocalLength(const Coupler* in)
{
    return in->connector->send->local_length;
}

} // namespace paso

#endif // __PASO_COUPLER_H__

