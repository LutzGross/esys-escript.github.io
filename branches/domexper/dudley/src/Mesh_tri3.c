
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

/*   Dudley: generates triangular meshes by splitting rectangles */

/*   Generates a numElements[0] x numElements[1] x 2 mesh with first order elements (Tri3) in the rectangle */
/*   [0,Length[0]] x [0,Length[1]]. order is the desired accuracy of the integration scheme. */

/**************************************************************/

#include "TriangularMesh.h"

Dudley_Mesh *Dudley_TriangularMesh_Tri3(dim_t * numElements,
					double *Length, index_t order, index_t reduced_order, bool_t optimize)
{
#define N_PER_E 1
#define DIM 2
    dim_t N0, N1, NE0, NE1, i0, i1, Nstride0 = 0, Nstride1 = 0, local_NE0, local_NE1, local_N0 = 0, local_N1 = 0;
    index_t offset0 = 0, offset1 = 0, e_offset0 = 0, e_offset1 = 0;
    dim_t totalNECount, faceNECount, NDOF0 = 0, NDOF1 = 0, NFaceElements;
    index_t myRank;
    Dudley_Mesh *out;
    Esys_MPIInfo *mpi_info = NULL;
    char name[50];
    const int LEFTTAG = 1;	/* boundary x1=0 */
    const int RIGHTTAG = 2;	/* boundary x1=1 */
    const int BOTTOMTAG = 10;	/* boundary x2=0 */
    const int TOPTAG = 20;	/* boundary x2=1 */

#ifdef Dudley_TRACE
    double time0 = Dudley_timer();
#endif

    /* get MPI information */
    mpi_info = Esys_MPIInfo_alloc(MPI_COMM_WORLD);
    if (!Dudley_noError())
    {
	return NULL;
    }
    myRank = mpi_info->rank;

    /* set up the global dimensions of the mesh */

    NE0 = MAX(1, numElements[0]);
    NE1 = MAX(1, numElements[1]);
    N0 = N_PER_E * NE0 + 1;
    N1 = N_PER_E * NE1 + 1;

    /* This code was originally copied from Finley's Rec4 constructor.
       NE? refers to the number of rectangular elements in each direction.
       The number of nodes produced is the same but the number of non-face elements
       will double.
     */

    /*  allocate mesh: */
    sprintf(name, "Triangular %d x %d (x 2) mesh", N0, N1);
    out = Dudley_Mesh_alloc(name, DIM, mpi_info);
    if (!Dudley_noError())
    {
	Esys_MPIInfo_free(mpi_info);
	return NULL;
    }
    if (Dudley_noError())
    {

	Dudley_Mesh_setPoints(out, Dudley_ElementFile_alloc(Point1, mpi_info));
	Dudley_Mesh_setFaceElements(out, Dudley_ElementFile_alloc(Line2, mpi_info));
	Dudley_Mesh_setElements(out, Dudley_ElementFile_alloc(Tri3, mpi_info));
	Nstride0 = 1;
	Nstride1 = N0;
	if (N1 == MAX(N0, N1))
	{
	    local_NE0 = NE0;
	    e_offset0 = 0;
	    Esys_MPIInfo_Split(mpi_info, NE1, &local_NE1, &e_offset1);
	}
	else
	{
	    Esys_MPIInfo_Split(mpi_info, NE0, &local_NE0, &e_offset0);
	    local_NE1 = NE1;
	    e_offset1 = 0;
	}
	offset0 = e_offset0 * N_PER_E;
	offset1 = e_offset1 * N_PER_E;
	local_N0 = local_NE0 > 0 ? local_NE0 * N_PER_E + 1 : 0;
	local_N1 = local_NE1 > 0 ? local_NE1 * N_PER_E + 1 : 0;

	/* get the number of surface elements */

	NFaceElements = 0;
	if (local_NE0 > 0)
	{
	    NDOF0 = N0;
	    if (e_offset0 == 0)
		NFaceElements += local_NE1;
	    if (local_NE0 + e_offset0 == NE0)
		NFaceElements += local_NE1;
	}
	else
	{
	    NDOF0 = N0 - 1;
	}
	if (local_NE1 > 0)
	{
	    NDOF1 = N1;
	    if (e_offset1 == 0)
		NFaceElements += local_NE0;
	    if (local_NE1 + e_offset1 == NE1)
		NFaceElements += local_NE0;
	}
	else
	{
	    NDOF1 = N1 - 1;
	}

	/*  allocate tables: */

	Dudley_NodeFile_allocTable(out->Nodes, local_N0 * local_N1);

	/* This code was oringinally copied from Finley's rec4 generator 
	   We double these numbers because each "rectangle" will be split into
	   two triangles. So the number of nodes is the same but the 
	   number of elements will double */
	Dudley_ElementFile_allocTable(out->Elements, local_NE0 * local_NE1 * 2);
	Dudley_ElementFile_allocTable(out->FaceElements, NFaceElements);

    }
    if (Dudley_noError())
    {
	/* create nodes */
#pragma omp parallel for private(i0,i1)
	for (i1 = 0; i1 < local_N1; i1++)
	{
	    for (i0 = 0; i0 < local_N0; i0++)
	    {
		dim_t k = i0 + local_N0 * i1;
		dim_t global_i0 = i0 + offset0;
		dim_t global_i1 = i1 + offset1;
		out->Nodes->Coordinates[INDEX2(0, k, DIM)] = DBLE(global_i0) / DBLE(N0 - 1) * Length[0];
		out->Nodes->Coordinates[INDEX2(1, k, DIM)] = DBLE(global_i1) / DBLE(N1 - 1) * Length[1];
		out->Nodes->Id[k] = Nstride0 * global_i0 + Nstride1 * global_i1;
		out->Nodes->Tag[k] = 0;
		out->Nodes->globalDegreesOfFreedom[k] = Nstride0 * (global_i0 % NDOF0) + Nstride1 * (global_i1 % NDOF1);
	    }
	}
	/*   set the elements: */
	dim_t NN = out->Elements->numNodes;
	index_t global_adjustment = (offset0 + offset1) % 2;
#pragma omp parallel for private(i0,i1)
	for (i1 = 0; i1 < local_NE1; i1++)
	{
	    for (i0 = 0; i0 < local_NE0; i0++)
	    {
		/* we will split this "rectangle" into two triangles */
		dim_t k = 2 * (i0 + local_NE0 * i1);
		index_t node0 = Nstride0 * N_PER_E * (i0 + e_offset0) + Nstride1 * N_PER_E * (i1 + e_offset1);

		out->Elements->Id[k] = 2 * ((i0 + e_offset0) + NE0 * (i1 + e_offset1));
		out->Elements->Tag[k] = 0;
		out->Elements->Owner[k] = myRank;
		out->Elements->Id[k + 1] = out->Elements->Id[k] + 1;
		out->Elements->Tag[k + 1] = 0;
		out->Elements->Owner[k + 1] = myRank;

		/* a,b,c,d gives the nodes in the rectangle in clockwise order */
		index_t a = node0, b = node0 + Nstride0, c = node0 + Nstride1 + Nstride0, d = node0 + Nstride1;
		/* For a little bit of variety  */
		if ((global_adjustment + node0) % 2)
		{
		    out->Elements->Nodes[INDEX2(0, k, NN)] = a;
		    out->Elements->Nodes[INDEX2(1, k, NN)] = b;
		    out->Elements->Nodes[INDEX2(2, k, NN)] = d;
		    out->Elements->Nodes[INDEX2(0, k + 1, NN)] = b;
		    out->Elements->Nodes[INDEX2(1, k + 1, NN)] = c;
		    out->Elements->Nodes[INDEX2(2, k + 1, NN)] = d;
		}
		else
		{
		    out->Elements->Nodes[INDEX2(0, k, NN)] = a;
		    out->Elements->Nodes[INDEX2(1, k, NN)] = b;
		    out->Elements->Nodes[INDEX2(2, k, NN)] = c;
		    out->Elements->Nodes[INDEX2(0, k + 1, NN)] = a;
		    out->Elements->Nodes[INDEX2(1, k + 1, NN)] = c;
		    out->Elements->Nodes[INDEX2(2, k + 1, NN)] = d;
		}
	    }
	}
	/* face elements */
	NN = out->FaceElements->numNodes;
	totalNECount = 2 * NE0 * NE1;	/* because we have split the rectangles */
	faceNECount = 0;
	if (local_NE0 > 0)
	{
	    /* **  elements on boundary 001 (x1=0): */

	    if (e_offset0 == 0)
	    {
#pragma omp parallel for private(i1)
		for (i1 = 0; i1 < local_NE1; i1++)
		{

		    dim_t k = i1 + faceNECount;
		    index_t node0 = Nstride1 * N_PER_E * (i1 + e_offset1);

		    out->FaceElements->Id[k] = i1 + e_offset1 + totalNECount;
		    out->FaceElements->Tag[k] = LEFTTAG;
		    out->FaceElements->Owner[k] = myRank;
		    out->FaceElements->Nodes[INDEX2(0, k, NN)] = node0 + Nstride1;
		    out->FaceElements->Nodes[INDEX2(1, k, NN)] = node0;
		}
		faceNECount += local_NE1;
	    }
	    totalNECount += NE1;
	    /* **  elements on boundary 002 (x1=1): */
	    if (local_NE0 + e_offset0 == NE0)
	    {
#pragma omp parallel for private(i1)
		for (i1 = 0; i1 < local_NE1; i1++)
		{
		    dim_t k = i1 + faceNECount;
		    index_t node0 = Nstride0 * N_PER_E * (NE0 - 1) + Nstride1 * N_PER_E * (i1 + e_offset1);

		    out->FaceElements->Id[k] = (i1 + e_offset1) + totalNECount;
		    out->FaceElements->Tag[k] = RIGHTTAG;
		    out->FaceElements->Owner[k] = myRank;
		    out->FaceElements->Nodes[INDEX2(0, k, NN)] = node0 + Nstride0;
		    out->FaceElements->Nodes[INDEX2(1, k, NN)] = node0 + Nstride1 + Nstride0;
		}
		faceNECount += local_NE1;
	    }
	    totalNECount += NE1;
	}
	if (local_NE1 > 0)
	{
	    /* **  elements on boundary 010 (x2=0): */
	    if (e_offset1 == 0)
	    {
#pragma omp parallel for private(i0)
		for (i0 = 0; i0 < local_NE0; i0++)
		{
		    dim_t k = i0 + faceNECount;
		    index_t node0 = Nstride0 * N_PER_E * (i0 + e_offset0);

		    out->FaceElements->Id[k] = e_offset0 + i0 + totalNECount;
		    out->FaceElements->Tag[k] = BOTTOMTAG;
		    out->FaceElements->Owner[k] = myRank;

		    out->FaceElements->Nodes[INDEX2(0, k, NN)] = node0;
		    out->FaceElements->Nodes[INDEX2(1, k, NN)] = node0 + Nstride0;
		}
		faceNECount += local_NE0;
	    }
	    totalNECount += NE0;
	    /* **  elements on boundary 020 (x2=1): */
	    if (local_NE1 + e_offset1 == NE1)
	    {
#pragma omp parallel for private(i0)
		for (i0 = 0; i0 < local_NE0; i0++)
		{
		    dim_t k = i0 + faceNECount;
		    index_t node0 = Nstride0 * N_PER_E * (i0 + e_offset0) + Nstride1 * N_PER_E * (NE1 - 1);

		    out->FaceElements->Id[k] = i0 + e_offset0 + totalNECount;
		    out->FaceElements->Tag[k] = TOPTAG;
		    out->FaceElements->Owner[k] = myRank;

		    out->FaceElements->Nodes[INDEX2(0, k, NN)] = node0 + Nstride1 + Nstride0;
		    out->FaceElements->Nodes[INDEX2(1, k, NN)] = node0 + Nstride1;
/*printf("E=%d: %d=%d %d=%d\n",k,INDEX2(0,k,NN),out->FaceElements->Nodes[INDEX2(0,k,NN)], 
INDEX2(1,k,NN),out->FaceElements->Nodes[INDEX2(1,k,NN)]); */
		}
		faceNECount += local_NE0;
	    }
	    totalNECount += NE0;
	}
    }
    if (Dudley_noError())
    {
	/* add tag names */
	Dudley_Mesh_addTagMap(out, "top", TOPTAG);
	Dudley_Mesh_addTagMap(out, "bottom", BOTTOMTAG);
	Dudley_Mesh_addTagMap(out, "left", LEFTTAG);
	Dudley_Mesh_addTagMap(out, "right", RIGHTTAG);
    }
    /* prepare mesh for further calculatuions: */
    if (Dudley_noError())
    {
	Dudley_Mesh_resolveNodeIds(out);
    }
    if (Dudley_noError())
    {
	Dudley_Mesh_prepare(out, optimize);
    }

    /* free up memory */
    Esys_MPIInfo_free(mpi_info);

    return out;
}
