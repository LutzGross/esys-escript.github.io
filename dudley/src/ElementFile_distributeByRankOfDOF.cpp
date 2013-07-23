
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/

/************************************************************************************/

/*   Dudley: ElementFile: this will redistribute the Elements including overlap by */

/************************************************************************************/

#include "ElementFile.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/************************************************************************************/

void Dudley_ElementFile_distributeByRankOfDOF(Dudley_ElementFile * self, Esys_MPI_rank * mpiRankOfDOF, index_t * Id)
{
    size_t size_size;
    Esys_MPI_rank myRank, p, *Owner_buffer = NULL, loc_proc_mask_max;
    dim_t e, j, i, size, *send_count = NULL, *recv_count = NULL, *newOwner = NULL, *loc_proc_mask =
	NULL, *loc_send_count = NULL, newNumElements, numElementsInBuffer, numNodes, NN;
    index_t *send_offset = NULL, *recv_offset = NULL, *Id_buffer = NULL, *Tag_buffer = NULL, *Nodes_buffer = NULL, k;
    bool_t *proc_mask = NULL;
#ifdef ESYS_MPI
    dim_t numRequests = 0;
    MPI_Request *mpi_requests = NULL;
    MPI_Status *mpi_stati = NULL;
#endif
    if (self == NULL)
	return;
    myRank = self->MPIInfo->rank;
    size = self->MPIInfo->size;
    size_size = size * sizeof(dim_t);
    numNodes = self->numNodes;
    NN = self->numNodes;
    if (size > 1)
    {
#ifdef ESYS_MPI
	mpi_requests = new  MPI_Request[8 * size];
	mpi_stati = new  MPI_Status[8 * size];
	Dudley_checkPtr(mpi_requests);
	Dudley_checkPtr(mpi_stati);
#endif

	/* count the number elements that have to be send to each processor (send_count) 
	   and define a new element owner as the processor with the largest number of DOFs and the smallest id */
	send_count = new  dim_t[size];
	recv_count = new  dim_t[size];
	newOwner = new  Esys_MPI_rank[self->numElements];
	if (!(Dudley_checkPtr(send_count) || Dudley_checkPtr(recv_count) || Dudley_checkPtr(newOwner)))
	{
	    memset(send_count, 0, size_size);
#pragma omp parallel private(p,loc_proc_mask,loc_send_count)
	    {
		loc_proc_mask = new  dim_t[size];
		loc_send_count = new  dim_t[size];
		memset(loc_send_count, 0, size_size);
#pragma omp for private(e,j,loc_proc_mask_max) schedule(static)
		for (e = 0; e < self->numElements; e++)
		{
		    if (self->Owner[e] == myRank)
		    {
			newOwner[e] = myRank;
			memset(loc_proc_mask, 0, size_size);
			for (j = 0; j < numNodes; j++)
			{
			    p = mpiRankOfDOF[self->Nodes[INDEX2(j, e, NN)]];
			    loc_proc_mask[p]++;
			}
			loc_proc_mask_max = 0;
			for (p = 0; p < size; ++p)
			{
			    if (loc_proc_mask[p] > 0)
				loc_send_count[p]++;
			    if (loc_proc_mask[p] > loc_proc_mask_max)
			    {
				newOwner[e] = p;
				loc_proc_mask_max = loc_proc_mask[p];
			    }
			}
		    }
		    else
		    {
			newOwner[e] = -1;
		    }
		}
#pragma omp critical
		{
		    for (p = 0; p < size; ++p)
			send_count[p] += loc_send_count[p];
		}
		delete[] loc_proc_mask;
		delete[] loc_send_count;
	    }
#ifdef ESYS_MPI
	    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, self->MPIInfo->comm);
#else
	    for (p = 0; p < size; ++p)
		recv_count[p] = send_count[p];
#endif
	    /* get the new number of elements for this processor */
	    newNumElements = 0;
	    for (p = 0; p < size; ++p)
		newNumElements += recv_count[p];

	    /* get the new number of elements for this processor */
	    numElementsInBuffer = 0;
	    for (p = 0; p < size; ++p)
		numElementsInBuffer += send_count[p];
	    /* allocate buffers */
	    Id_buffer = new  index_t[numElementsInBuffer];
	    Tag_buffer = new  index_t[numElementsInBuffer];
	    Owner_buffer = new  Esys_MPI_rank[numElementsInBuffer];
	    Nodes_buffer = new  index_t[numElementsInBuffer * NN];
	    send_offset = new  index_t[size];
	    recv_offset = new  index_t[size];
	    proc_mask = new  bool_t[size];
	    if (!(Dudley_checkPtr(Id_buffer) || Dudley_checkPtr(Tag_buffer) || Dudley_checkPtr(Owner_buffer) ||
		  Dudley_checkPtr(Nodes_buffer) || Dudley_checkPtr(send_offset) || Dudley_checkPtr(recv_offset) ||
		  Dudley_checkPtr(proc_mask)))
	    {

		/* calculate the offsets for the processor buffers */
		recv_offset[0] = 0;
		for (p = 0; p < size - 1; ++p)
		    recv_offset[p + 1] = recv_offset[p] + recv_count[p];
		send_offset[0] = 0;
		for (p = 0; p < size - 1; ++p)
		    send_offset[p + 1] = send_offset[p] + send_count[p];

		memset(send_count, 0, size_size);
		/* copy element into buffers. proc_mask makes sure that an element is copied once only for each processor */
		for (e = 0; e < self->numElements; e++)
		{
		    if (self->Owner[e] == myRank)
		    {
			memset(proc_mask, TRUE, size_size);
			for (j = 0; j < numNodes; j++)
			{
			    p = mpiRankOfDOF[self->Nodes[INDEX2(j, e, NN)]];
			    if (proc_mask[p])
			    {
				k = send_offset[p] + send_count[p];
				Id_buffer[k] = self->Id[e];
				Tag_buffer[k] = self->Tag[e];
				Owner_buffer[k] = newOwner[e];
				for (i = 0; i < numNodes; i++)
				    Nodes_buffer[INDEX2(i, k, NN)] = Id[self->Nodes[INDEX2(i, e, NN)]];
				send_count[p]++;
				proc_mask[p] = FALSE;
			    }
			}
		    }
		}
		/* allocate new tables */
		Dudley_ElementFile_allocTable(self, newNumElements);

		/* start to receive new elements */
		for (p = 0; p < size; ++p)
		{
		    if (recv_count[p] > 0)
		    {
#ifdef ESYS_MPI
			MPI_Irecv(&(self->Id[recv_offset[p]]), recv_count[p],
				  MPI_INT, p, self->MPIInfo->msg_tag_counter + myRank,
				  self->MPIInfo->comm, &mpi_requests[numRequests]);
			numRequests++;
			MPI_Irecv(&(self->Tag[recv_offset[p]]), recv_count[p],
				  MPI_INT, p, self->MPIInfo->msg_tag_counter + size + myRank,
				  self->MPIInfo->comm, &mpi_requests[numRequests]);
			numRequests++;
			MPI_Irecv(&(self->Owner[recv_offset[p]]), recv_count[p],
				  MPI_INT, p, self->MPIInfo->msg_tag_counter + 2 * size + myRank,
				  self->MPIInfo->comm, &mpi_requests[numRequests]);
			numRequests++;
			MPI_Irecv(&(self->Nodes[recv_offset[p] * NN]), recv_count[p] * NN,
				  MPI_INT, p, self->MPIInfo->msg_tag_counter + 3 * size + myRank,
				  self->MPIInfo->comm, &mpi_requests[numRequests]);
			numRequests++;
#endif
		    }
		}
		/* now the buffers can be send away */
		for (p = 0; p < size; ++p)
		{
		    if (send_count[p] > 0)
		    {
#ifdef ESYS_MPI
			MPI_Issend(&(Id_buffer[send_offset[p]]), send_count[p],
				   MPI_INT, p, self->MPIInfo->msg_tag_counter + p,
				   self->MPIInfo->comm, &mpi_requests[numRequests]);
			numRequests++;
			MPI_Issend(&(Tag_buffer[send_offset[p]]), send_count[p],
				   MPI_INT, p, self->MPIInfo->msg_tag_counter + size + p,
				   self->MPIInfo->comm, &mpi_requests[numRequests]);
			numRequests++;
			MPI_Issend(&(Owner_buffer[send_offset[p]]), send_count[p],
				   MPI_INT, p, self->MPIInfo->msg_tag_counter + 2 * size + p,
				   self->MPIInfo->comm, &mpi_requests[numRequests]);
			numRequests++;
			MPI_Issend(&(Nodes_buffer[send_offset[p] * NN]), send_count[p] * NN,
				   MPI_INT, p, self->MPIInfo->msg_tag_counter + 3 * size + p,
				   self->MPIInfo->comm, &mpi_requests[numRequests]);
			numRequests++;
#endif

		    }
		}
		ESYS_MPI_INC_COUNTER(*(self->MPIInfo), 4 * size);
		/* wait for the requests to be finalized */
#ifdef ESYS_MPI
		MPI_Waitall(numRequests, mpi_requests, mpi_stati);
#endif
	    }
	    /* clear buffer */
	    delete[] Id_buffer;
	    delete[] Tag_buffer;
	    delete[] Owner_buffer;
	    delete[] Nodes_buffer;
	    delete[] send_offset;
	    delete[] recv_offset;
	    delete[] proc_mask;
	}
#ifdef ESYS_MPI
	delete[] mpi_requests;
	delete[] mpi_stati;
#endif
	delete[] send_count;
	delete[] recv_count;
	delete[] newOwner;
    }
    else
    {
#pragma omp for private(e,i) schedule(static)
	for (e = 0; e < self->numElements; e++)
	{
	    self->Owner[e] = myRank;
	    for (i = 0; i < numNodes; i++)
		self->Nodes[INDEX2(i, e, NN)] = Id[self->Nodes[INDEX2(i, e, NN)]];
	}
    }
    return;
}
