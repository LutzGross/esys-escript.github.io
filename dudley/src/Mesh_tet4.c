
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

/*   Dudley: generates rectangular meshes */

/*   Generates a numElements[0] x numElements[1] x numElements[2] mesh with first order elements (Hex8) in the brick */
/*   [0,Length[0]] x [0,Length[1]] x [0,Length[2]]. order is the desired accuracy of the */
/*   integration scheme. */

/**************************************************************/

#include "TriangularMesh.h"


/* Be careful reading this function. The X? and NStride? are 1,2,3 but the loop vars are 0,1,2 */
Dudley_Mesh* Dudley_TriangularMesh_Tet4(dim_t* numElements,
                                          double* Length,
					  index_t order,
                                          index_t reduced_order, 
                                          bool_t optimize)
{
  #define N_PER_E 1
  #define DIM 3
  dim_t N0,N1,N2,NE0,NE1,NE2,i0,i1,i2,k,Nstride0=0, Nstride1=0,Nstride2=0, local_NE0, local_NE1, local_NE2, local_N0=0, local_N1=0, local_N2=0;
  dim_t totalNECount,faceNECount,NDOF0=0,NDOF1=0,NDOF2=0,NFaceElements=0, NN;
  index_t node0, myRank, e_offset2, e_offset1, e_offset0=0, offset1=0, offset2=0, offset0=0, global_i0, global_i1, global_i2;
  Dudley_ReferenceElementSet *refPoints=NULL, *refFaceElements=NULL, *refElements=NULL;
  Dudley_Mesh* out;
  Paso_MPIInfo *mpi_info = NULL;
  char name[50];
  #ifdef Dudley_TRACE
  double time0=Dudley_timer();
  #endif

  const int LEFTTAG=1;		/* boundary x1=0 */
  const int RIGHTTAG=2;		/* boundary x1=1 */
  const int BOTTOMTAG=100;	/* boundary x3=1 */
  const int TOPTAG=200;		/* boundary x3=0 */
  const int FRONTTAG=10;	/* boundary x2=0 */
  const int BACKTAG=20;		/* boundary x2=1 */




  /* get MPI information */
  mpi_info = Paso_MPIInfo_alloc( MPI_COMM_WORLD );
  if (! Dudley_noError()) {
        return NULL;
  }
  myRank=mpi_info->rank;

  /* set up the global dimensions of the mesh */

  NE0=MAX(1,numElements[0]);
  NE1=MAX(1,numElements[1]);
  NE2=MAX(1,numElements[2]);
  N0=N_PER_E*NE0+1;
  N1=N_PER_E*NE1+1;
  N2=N_PER_E*NE2+1;

  /*  allocate mesh: */  
  sprintf(name,"Triangular %d x %d x %d (x 5) mesh",N0,N1,N2);
  out=Dudley_Mesh_alloc(name,DIM, mpi_info);
  if (! Dudley_noError()) { 
      Paso_MPIInfo_free( mpi_info );
      return NULL;
  }
  refElements= Dudley_ReferenceElementSet_alloc(Tet4,order,reduced_order);
  refFaceElements=Dudley_ReferenceElementSet_alloc(Tri3, order, reduced_order);
  refPoints=Dudley_ReferenceElementSet_alloc(Point1, order, reduced_order);
  
 
  if ( Dudley_noError()) {
  
	  Dudley_Mesh_setPoints(out,Dudley_ElementFile_alloc(refPoints, mpi_info));
	  Dudley_Mesh_setFaceElements(out,Dudley_ElementFile_alloc(refFaceElements, mpi_info));
	  Dudley_Mesh_setElements(out,Dudley_ElementFile_alloc(refElements, mpi_info));

	  /* work out the largest dimension */
	  if (N2==MAX3(N0,N1,N2)) {
		  Nstride0=1;
		  Nstride1=N0;
		  Nstride2=N0*N1;
		  local_NE0=NE0;
		  e_offset0=0;
		  local_NE1=NE1;
		  e_offset1=0;
		  Paso_MPIInfo_Split(mpi_info,NE2,&local_NE2,&e_offset2);
	  } else if (N1==MAX3(N0,N1,N2)) {
		  Nstride0=N2;
		  Nstride1=N0*N2;
		  Nstride2=1;
		  local_NE0=NE0;
 	    e_offset0=0;
		Paso_MPIInfo_Split(mpi_info,NE1,&local_NE1,&e_offset1);
		local_NE2=NE2;
		e_offset2=0;
	  } else {
		  Nstride0=N1*N2;
		  Nstride1=1;
		  Nstride2=N1;
		  Paso_MPIInfo_Split(mpi_info,NE0,&local_NE0,&e_offset0);
		  local_NE1=NE1;
		  e_offset1=0;
		  local_NE2=NE2;
		  e_offset2=0;
	  }
	  offset0=e_offset0*N_PER_E;
	  offset1=e_offset1*N_PER_E;
	  offset2=e_offset2*N_PER_E;
	  local_N0=local_NE0>0 ? local_NE0*N_PER_E+1 : 0;
	  local_N1=local_NE1>0 ? local_NE1*N_PER_E+1 : 0;
	  local_N2=local_NE2>0 ? local_NE2*N_PER_E+1 : 0;

	  /* get the number of surface elements */

	  NFaceElements=0;
	  if (local_NE2>0) {
		  NDOF2=N2;
		  if (offset2==0) NFaceElements+=2*local_NE1*local_NE0;		// each face is split
		  if (local_NE2+e_offset2 == NE2) NFaceElements+=2*local_NE1*local_NE0;
	  } else {
		  NDOF2=N2-1;
	  }
 
	  if (local_NE0>0) {
		  NDOF0=N0;
		  if (e_offset0 == 0) NFaceElements+=2*local_NE1*local_NE2;
		  if (local_NE0+e_offset0 == NE0) NFaceElements+=2*local_NE1*local_NE2;
	  } else {
		  NDOF0=N0-1;
	  }
	 	 
	  if (local_NE1>0) {
		NDOF1=N1;
		if (e_offset1 == 0) NFaceElements+=2*local_NE0*local_NE2;
		if (local_NE1+e_offset1 == NE1) NFaceElements+=2*local_NE0*local_NE2;
	  } else {
		NDOF1=N1-1;
	  }
  }
  
  
  /*  allocate tables: */
  if (Dudley_noError()) {

	  Dudley_NodeFile_allocTable(out->Nodes,local_N0*local_N1*local_N2);
	/* we split the rectangular prism this code used to produce into 5 tetrahedrons */
	  Dudley_ElementFile_allocTable(out->Elements,local_NE0*local_NE1*local_NE2*5);
	/* each border face will be split in half */
	  Dudley_ElementFile_allocTable(out->FaceElements,NFaceElements);
  }
  
  if (Dudley_noError()) {
     /* create nodes */
   
     #pragma omp parallel for private(i0,i1,i2,k,global_i0,global_i1,global_i2)
     for (i2=0;i2<local_N2;i2++) {
       for (i1=0;i1<local_N1;i1++) {
         for (i0=0;i0<local_N0;i0++) {
           k=i0+local_N0*i1+local_N0*local_N1*i2;
           global_i0=i0+offset0;
           global_i1=i1+offset1;
           global_i2=i2+offset2;
           out->Nodes->Coordinates[INDEX2(0,k,DIM)]=DBLE(global_i0)/DBLE(N0-1)*Length[0];
           out->Nodes->Coordinates[INDEX2(1,k,DIM)]=DBLE(global_i1)/DBLE(N1-1)*Length[1];
           out->Nodes->Coordinates[INDEX2(2,k,DIM)]=DBLE(global_i2)/DBLE(N2-1)*Length[2];
           out->Nodes->Id[k]=Nstride0*global_i0+Nstride1*global_i1+Nstride2*global_i2;
           out->Nodes->Tag[k]=0;
           out->Nodes->globalDegreesOfFreedom[k]=Nstride0*(global_i0%NDOF0) 
                                               +Nstride1*(global_i1%NDOF1) 
                                               +Nstride2*(global_i2%NDOF2);
// 	   printf("N=%d: %f,%f,%f\n", k, out->Nodes->Coordinates[INDEX2(0,k,DIM)],
// 		out->Nodes->Coordinates[INDEX2(1,k,DIM)],
// 		out->Nodes->Coordinates[INDEX2(2,k,DIM)]);
         }
       }
     }
// printf("Now listing elements\n");
     /*   set the elements: */

     int global_adjustment=0;	// If we are not the only rank we may need to shift our pattern to match neighbours

     NN=out->Elements->numNodes;
     #pragma omp parallel for private(i0,i1,i2,k,node0) 
     for (i2=0;i2<local_NE2;i2++) {
       for (i1=0;i1<local_NE1;i1++) {
         for (i0=0;i0<local_NE0;i0++) {
           
           k=5*(i0+local_NE0*i1+local_NE0*local_NE1*i2);
           node0=Nstride0*N_PER_E*(i0+e_offset0)+Nstride1*N_PER_E*(i1+e_offset1)+Nstride2*N_PER_E*(i2+e_offset2);

	   index_t res=5*((i0+e_offset0)+NE0*(i1+e_offset1)+NE0*NE1*(i2+e_offset2));
	   for (int j=0;j<5;++j)
	   {
               out->Elements->Id[k+j]=res+j;
               out->Elements->Tag[k+j]=0;
               out->Elements->Owner[k+j]=myRank;
	   }

	   index_t v[8];
/*	   in non-rotated orientation the points are numbered as follows:
	   The bottom face (anticlockwise= 0,1,3,2), top face (anticlockwise 4,5,7,6)*/
	   if ((global_adjustment+i0+i1+i2)%2==0)
	   {
// printf("Type0 %d, %d, %d\n",i0,i1,i2);
	   v[0]=node0;
	   v[1]=node0+Nstride0;
	   v[2]=node0+Nstride1;
	   v[3]=node0+Nstride1+Nstride0;
	   v[4]=node0+Nstride2;
	   v[5]=node0+Nstride0+Nstride2;
	   v[6]=node0+Nstride1+Nstride2;
	   v[7]=node0+Nstride2+Nstride1+Nstride0;
	   }
	   else
	   {
// printf("Type1 %d, %d, %d\n",i0,i1,i2);
	   // this form is rotated around the 0,2,4,6 face clockwise 90 degrees

	   v[0]=node0+Nstride1;			// node 0 ends up in position 2
	   v[2]=node0+Nstride1+Nstride2;	// node 2 ends up in position 6
	   v[6]=node0+Nstride2;			// node 6 ends up in position 4
	   v[4]=node0;				// node 4 ends up in position 0

	   v[1]=node0+Nstride0+Nstride1;	// node 1 -> pos 3
	   v[3]=node0+Nstride2+Nstride1+Nstride0;	// node 3-> pos 7
	   v[7]=node0+Nstride0+Nstride2;	// node 7 -> pos 5
	   v[5]=node0+Nstride0;			// node 5 -> pos 1
	   }

// for (int z=0;z<8;++z)
// {
//     printf("z[%d]=%d\n", z, v[z]);
// 
// }


// 	   index_t a=node0, b=node0+Nstride0, c=node0+Nstride1+Nstride0, d=node0+Nstride1;
// 	   index_t e=node0+Nstride2, f=node0+Nstride2+Nstride0, g=node0+Nstride2+Nstride1+Nstride0,
// 		 h=node0+Nstride2+Nstride1;
// 	   
// 
// 
// 	   a=0, b=1, c=3, d=2
// 	   e=4, f=5, g=7, h=6

	   // elements nodes are numbered: centre, x, y, z

	   out->Elements->Nodes[INDEX2(0,k,NN)] =v[4];
	   out->Elements->Nodes[INDEX2(1,k,NN)] =v[5];
	   out->Elements->Nodes[INDEX2(2,k,NN)] =v[6];
	   out->Elements->Nodes[INDEX2(3,k,NN)] =v[0];

	   out->Elements->Nodes[INDEX2(0,k+1,NN)] =v[7];
	   out->Elements->Nodes[INDEX2(1,k+1,NN)] =v[6];
	   out->Elements->Nodes[INDEX2(2,k+1,NN)] =v[5];
	   out->Elements->Nodes[INDEX2(3,k+1,NN)] =v[3];

	   out->Elements->Nodes[INDEX2(0,k+2,NN)] =v[2];
	   out->Elements->Nodes[INDEX2(1,k+2,NN)] =v[3];
	   out->Elements->Nodes[INDEX2(2,k+2,NN)] =v[0];
	   out->Elements->Nodes[INDEX2(3,k+2,NN)] =v[6];


	   out->Elements->Nodes[INDEX2(0,k+3,NN)] =v[1];
	   out->Elements->Nodes[INDEX2(1,k+3,NN)] =v[0];
	   out->Elements->Nodes[INDEX2(2,k+3,NN)] =v[3];
	   out->Elements->Nodes[INDEX2(3,k+3,NN)] =v[5];

	   // I can't work out where the center is for this one
	   out->Elements->Nodes[INDEX2(0,k+4,NN)] =v[5];
	   out->Elements->Nodes[INDEX2(1,k+4,NN)] =v[0];
	   out->Elements->Nodes[INDEX2(2,k+4,NN)] =v[6];
	   out->Elements->Nodes[INDEX2(3,k+4,NN)] =v[3];


// for (int z=0;z<5;++z)
// {
// printf("E %d:",z);
//    for (int q=0;q<4;++q)
//    {
// 	index_t id=out->Elements->Nodes[INDEX2(q,z,NN)];
// 	printf("   %d = %f, %f, %f\n", id, out->Nodes->Coordinates[INDEX2(0,id,DIM)],
// 		out->Nodes->Coordinates[INDEX2(1,id,DIM)], out->Nodes->Coordinates[INDEX2(2,id,DIM)]);
//    }
// printf("\n");}


/*
for (int j=0;j<4;++j)
{

  printf("Elt %d",j);
  for (int m=0;m<4;++m)
  {
     printf(" %d",out->Elements->Nodes[INDEX2(m,k+j,NN)]);
  }
  printf("\n");
}*/



         }

       }
     } /* end for */
     /* face elements */
// printf("Starting face elements\n");
     NN=out->FaceElements->numNodes;
     totalNECount=5*NE0*NE1*NE2;
     faceNECount=0;
//printf("Bottom\n");
     /*   these are the quadrilateral elements on boundary 1 (x3=0): */
     if (local_NE2>0) 
     {
       /* **  elements on boundary 100 (x3=0): */
       if (e_offset2==0)
       {
          #pragma omp parallel for private(i0,i1,k,node0) 
          for (i1=0;i1<local_NE1;i1++)
	  {
            for (i0=0;i0<local_NE0;i0++) 
	    {
              k=2*(i0+local_NE0*i1)+faceNECount;
              node0=Nstride0*N_PER_E*(i0+e_offset0)+Nstride1*N_PER_E*(i1+e_offset1);
	      index_t res=2*((i0+e_offset0)+NE0*(i1+e_offset1))+totalNECount;
              out->FaceElements->Id[k]=res;
              out->FaceElements->Tag[k]=BOTTOMTAG;
              out->FaceElements->Owner[k]=myRank;
              out->FaceElements->Id[k+1]=res+1;
              out->FaceElements->Tag[k+1]=BOTTOMTAG;
              out->FaceElements->Owner[k+1]=myRank;

	      index_t n0=node0;
	      index_t n1=node0+Nstride0;
	      index_t n2=node0+Nstride1;
	      index_t n3=node0+Nstride0+Nstride1;

	      if ((global_adjustment+i0+i1)%2==0)
	      {
	        out->FaceElements->Nodes[INDEX2(0,k,NN)]=n0;
	        out->FaceElements->Nodes[INDEX2(1,k,NN)]=n3;
	        out->FaceElements->Nodes[INDEX2(2,k,NN)]=n1;

	        out->FaceElements->Nodes[INDEX2(0,k+1,NN)]=n0;
	        out->FaceElements->Nodes[INDEX2(1,k+1,NN)]=n2;
	        out->FaceElements->Nodes[INDEX2(2,k+1,NN)]=n3;

	      }
	      else
	      {
	        out->FaceElements->Nodes[INDEX2(0,k,NN)]=n0;
	        out->FaceElements->Nodes[INDEX2(1,k,NN)]=n2;
	        out->FaceElements->Nodes[INDEX2(2,k,NN)]=n1;

	        out->FaceElements->Nodes[INDEX2(0,k+1,NN)]=n1;
	        out->FaceElements->Nodes[INDEX2(1,k+1,NN)]=n2;
	        out->FaceElements->Nodes[INDEX2(2,k+1,NN)]=n3;



	      }
// printf("%d: %d, %d, %d\n", k,out->FaceElements->Nodes[INDEX2(0,k,NN)],  out->FaceElements->Nodes[INDEX2(1,k,NN)],
// out->FaceElements->Nodes[INDEX2(2,k,NN)]);
// printf("%d: %d, %d, %d\n", k+1,out->FaceElements->Nodes[INDEX2(0,k+1,NN)],  out->FaceElements->Nodes[INDEX2(1,k+1,NN)],
// out->FaceElements->Nodes[INDEX2(2,k+1,NN)]);

            }
          }
// printf("\n");
          faceNECount+=2*local_NE1*local_NE0;
       }
       totalNECount+=2*NE1*NE0;

//printf("Top\n");
       /* **  elements on boundary 200 (x3=1) */
       if (local_NE2+e_offset2 == NE2) {
          #pragma omp parallel for private(i0,i1,k,node0) 
          for (i1=0;i1<local_NE1;i1++) {
            for (i0=0;i0<local_NE0;i0++) {
      
              k=2*(i0+local_NE0*i1)+faceNECount;
              node0=Nstride0*N_PER_E*(i0+e_offset0)+Nstride1*N_PER_E*(i1+e_offset1)+Nstride2*N_PER_E*(NE2-1);

	      index_t res=2*((i0+e_offset0)+NE0*(i1+e_offset1))+totalNECount;
              out->FaceElements->Id[k]=res;
              out->FaceElements->Tag[k]=TOPTAG;
              out->FaceElements->Owner[k]=myRank;
              out->FaceElements->Id[k+1]=res+1;
              out->FaceElements->Tag[k+1]=TOPTAG;
              out->FaceElements->Owner[k+1]=myRank;

	      index_t n4=node0+Nstride2;
	      index_t n5=node0+Nstride0+Nstride2;
	      index_t n6=node0+Nstride1+Nstride2;
	      index_t n7=node0+Nstride1+Nstride0+Nstride2;

	      if ((global_adjustment+i0+i1+local_NE2-1)%2==0)
	      {
	        out->FaceElements->Nodes[INDEX2(0,k,NN)]=n4;
	        out->FaceElements->Nodes[INDEX2(1,k,NN)]=n5;
	        out->FaceElements->Nodes[INDEX2(2,k,NN)]=n6;

	        out->FaceElements->Nodes[INDEX2(0,k+1,NN)]=n5;
	        out->FaceElements->Nodes[INDEX2(1,k+1,NN)]=n7;
	        out->FaceElements->Nodes[INDEX2(2,k+1,NN)]=n6;
	      }
	      else
	      {

	        out->FaceElements->Nodes[INDEX2(0,k,NN)]=n4;
	        out->FaceElements->Nodes[INDEX2(1,k,NN)]=n5;
	        out->FaceElements->Nodes[INDEX2(2,k,NN)]=n7;

	        out->FaceElements->Nodes[INDEX2(0,k+1,NN)]=n4;
	        out->FaceElements->Nodes[INDEX2(1,k+1,NN)]=n7;
	        out->FaceElements->Nodes[INDEX2(2,k+1,NN)]=n6;
	      }

// printf("%d: %d, %d, %d\n", k,out->FaceElements->Nodes[INDEX2(0,k,NN)],  out->FaceElements->Nodes[INDEX2(1,k,NN)],
// out->FaceElements->Nodes[INDEX2(2,k,NN)]);
// printf("%d: %d, %d, %d\n", k+1,out->FaceElements->Nodes[INDEX2(0,k+1,NN)],  out->FaceElements->Nodes[INDEX2(1,k+1,NN)],
// out->FaceElements->Nodes[INDEX2(2,k+1,NN)]);

            }
          }
/*printf("\n");*/
          faceNECount+=2*local_NE1*local_NE0;
       }
       totalNECount+=2*NE1*NE0;
     }
//printf("Left\n");
     if (local_NE0>0) {
        /* **  elements on boundary 001 (x1=0): */
     
        if (e_offset0 == 0) {
           #pragma omp parallel for private(i1,i2,k,node0) 
           for (i2=0;i2<local_NE2;i2++) {
             for (i1=0;i1<local_NE1;i1++) {
      
               k=2*(i1+local_NE1*i2)+faceNECount;

//printf("%d, %d ",k,k+1);
               node0=Nstride1*N_PER_E*(i1+e_offset1)+Nstride2*N_PER_E*(i2+e_offset2);
	       index_t res=2*((i1+e_offset1)+NE1*(i2+e_offset2))+totalNECount;
               out->FaceElements->Id[k]=res;
               out->FaceElements->Tag[k]=LEFTTAG;
               out->FaceElements->Owner[k]=myRank;
               out->FaceElements->Id[k+1]=res+1;
               out->FaceElements->Tag[k+1]=LEFTTAG;
               out->FaceElements->Owner[k+1]=myRank;
      
	      index_t n0=node0;
	      index_t n2=node0+Nstride1;
	      index_t n4=node0+Nstride2;
	      index_t n6=node0+Nstride1+Nstride2;
	      
	       if ((global_adjustment+0+i1+i2)%2==0)
	   {
	        out->FaceElements->Nodes[INDEX2(0,k,NN)]=n0;
	        out->FaceElements->Nodes[INDEX2(1,k,NN)]=n4;
	        out->FaceElements->Nodes[INDEX2(2,k,NN)]=n6;

	        out->FaceElements->Nodes[INDEX2(0,k+1,NN)]=n0;
	        out->FaceElements->Nodes[INDEX2(1,k+1,NN)]=n6;
	        out->FaceElements->Nodes[INDEX2(2,k+1,NN)]=n2;
	   }
	   else
	   {
	   // this form is rotated around the 0,2,4,6 face clockwise 90 degrees

	        out->FaceElements->Nodes[INDEX2(0,k,NN)]=n0;
	        out->FaceElements->Nodes[INDEX2(1,k,NN)]=n4;
	        out->FaceElements->Nodes[INDEX2(2,k,NN)]=n2;

	        out->FaceElements->Nodes[INDEX2(0,k+1,NN)]=n4;
	        out->FaceElements->Nodes[INDEX2(1,k+1,NN)]=n6;
	        out->FaceElements->Nodes[INDEX2(2,k+1,NN)]=n2;
	   }
// printf("%d: %d, %d, %d\n", k,out->FaceElements->Nodes[INDEX2(0,k,NN)],  out->FaceElements->Nodes[INDEX2(1,k,NN)],
// out->FaceElements->Nodes[INDEX2(2,k,NN)]);
// printf("%d: %d, %d, %d\n", k+1,out->FaceElements->Nodes[INDEX2(0,k+1,NN)],  out->FaceElements->Nodes[INDEX2(1,k+1,NN)],
// out->FaceElements->Nodes[INDEX2(2,k+1,NN)]);
             }
           }
           faceNECount+=2*local_NE1*local_NE2;
        }
/*printf("\n");*/
//printf("Right\n");
        totalNECount+=2*NE1*NE2;
        /* **  elements on boundary 002 (x1=1): */
        if (local_NE0+e_offset0 == NE0) {
           #pragma omp parallel for private(i1,i2,k,node0) 
           for (i2=0;i2<local_NE2;i2++) {
             for (i1=0;i1<local_NE1;i1++) {
               k=2*(i1+local_NE1*i2)+faceNECount;

//printf("%d, %d ",k,k+1);
               node0=Nstride0*N_PER_E*(NE0-1)+Nstride1*N_PER_E*(i1+e_offset1)+Nstride2*N_PER_E*(i2+e_offset2);
	       index_t res=2*((i1+e_offset1)+NE1*(i2+e_offset2))+totalNECount;
               out->FaceElements->Id[k]=res;
               out->FaceElements->Tag[k]=RIGHTTAG;
               out->FaceElements->Owner[k]=myRank;
               out->FaceElements->Id[k+1]=res+1;
               out->FaceElements->Tag[k+1]=RIGHTTAG;
               out->FaceElements->Owner[k+1]=myRank;

	      index_t n1=node0+Nstride0;
	      index_t n3=node0+Nstride0+Nstride1;
	      index_t n5=node0+Nstride0+Nstride2;
	      index_t n7=node0+Nstride0+Nstride1+Nstride2;
	       if ((global_adjustment+local_NE0-1+i1+i2)%2==0)
	   {
	        out->FaceElements->Nodes[INDEX2(0,k,NN)]=n1;
	        out->FaceElements->Nodes[INDEX2(1,k,NN)]=n3;
	        out->FaceElements->Nodes[INDEX2(2,k,NN)]=n5;

	        out->FaceElements->Nodes[INDEX2(0,k+1,NN)]=n3;
	        out->FaceElements->Nodes[INDEX2(1,k+1,NN)]=n7;
	        out->FaceElements->Nodes[INDEX2(2,k+1,NN)]=n5;
	   }
	   else
	   {
	   // this form is rotated around the 0,2,4,6 face clockwise 90 degrees

	        out->FaceElements->Nodes[INDEX2(0,k,NN)]=n1;
	        out->FaceElements->Nodes[INDEX2(1,k,NN)]=n7;
	        out->FaceElements->Nodes[INDEX2(2,k,NN)]=n5;

	        out->FaceElements->Nodes[INDEX2(0,k+1,NN)]=n1;
	        out->FaceElements->Nodes[INDEX2(1,k+1,NN)]=n3;
	        out->FaceElements->Nodes[INDEX2(2,k+1,NN)]=n7;
	   }
// printf("%d: %d, %d, %d\n", k,out->FaceElements->Nodes[INDEX2(0,k,NN)],  out->FaceElements->Nodes[INDEX2(1,k,NN)],
// out->FaceElements->Nodes[INDEX2(2,k,NN)]);
// printf("%d: %d, %d, %d\n", k+1,out->FaceElements->Nodes[INDEX2(0,k+1,NN)],  out->FaceElements->Nodes[INDEX2(1,k+1,NN)],
// out->FaceElements->Nodes[INDEX2(2,k+1,NN)]);

             }
           }
           faceNECount+=2*local_NE1*local_NE2;
         }
         totalNECount+=2*NE1*NE2;
     }
//printf("\n");
//printf("Front\n");
     if (local_NE1>0) {
        /* **  elements on boundary 010 (x2=0): */
        if (e_offset1 == 0) {
           #pragma omp parallel for private(i0,i2,k,node0) 
           for (i2=0;i2<local_NE2;i2++) {
             for (i0=0;i0<local_NE0;i0++) {
               k=2*(i0+local_NE0*i2)+faceNECount;
               node0=Nstride0*N_PER_E*(i0+e_offset0)+Nstride2*N_PER_E*(i2+e_offset2);
//printf("%d, %d ",k,k+1);         
	       index_t res=2*((i2+e_offset2)+NE2*(e_offset0+i0))+totalNECount;
               out->FaceElements->Id[k]=res;
               out->FaceElements->Tag[k]=FRONTTAG;
               out->FaceElements->Owner[k]=myRank;
               out->FaceElements->Id[k+1]=res+1;
               out->FaceElements->Tag[k+1]=FRONTTAG;
               out->FaceElements->Owner[k+1]=myRank;

	      index_t n0=node0;
	      index_t n1=node0+Nstride0;
	      index_t n4=node0+Nstride2;
	      index_t n5=node0+Nstride0+Nstride2;

	       if ((global_adjustment+i0+0+i2)%2==0)
	   {
	        out->FaceElements->Nodes[INDEX2(0,k,NN)]=n0;
	        out->FaceElements->Nodes[INDEX2(1,k,NN)]=n1;
	        out->FaceElements->Nodes[INDEX2(2,k,NN)]=n5;

	        out->FaceElements->Nodes[INDEX2(0,k+1,NN)]=n0;
	        out->FaceElements->Nodes[INDEX2(1,k+1,NN)]=n5;
	        out->FaceElements->Nodes[INDEX2(2,k+1,NN)]=n4;


	   }
	   else
	   {
	   // this form is rotated around the 0,2,4,6 face clockwise 90 degrees

	        out->FaceElements->Nodes[INDEX2(0,k,NN)]=n0;
	        out->FaceElements->Nodes[INDEX2(1,k,NN)]=n1;
	        out->FaceElements->Nodes[INDEX2(2,k,NN)]=n4;

	        out->FaceElements->Nodes[INDEX2(0,k+1,NN)]=n1;
	        out->FaceElements->Nodes[INDEX2(1,k+1,NN)]=n5;
	        out->FaceElements->Nodes[INDEX2(2,k+1,NN)]=n4;

	   }
/*printf("%d: %d, %d, %d\n", k,out->FaceElements->Nodes[INDEX2(0,k,NN)],  out->FaceElements->Nodes[INDEX2(1,k,NN)],
out->FaceElements->Nodes[INDEX2(2,k,NN)]);
printf("%d: %d, %d, %d\n", k+1,out->FaceElements->Nodes[INDEX2(0,k+1,NN)],  out->FaceElements->Nodes[INDEX2(1,k+1,NN)],
out->FaceElements->Nodes[INDEX2(2,k+1,NN)]);*/
             }
           }
           faceNECount+=2*local_NE0*local_NE2;
        }
// printf("\n");
        totalNECount+=2*NE0*NE2;
//printf("Back\n");
        /* **  elements on boundary 020 (x2=1): */
        if (local_NE1+e_offset1 == NE1) {
           #pragma omp parallel for private(i0,i2,k,node0) 
           for (i2=0;i2<local_NE2;i2++) {
             for (i0=0;i0<local_NE0;i0++) {
               k=2*(i0+local_NE0*i2)+faceNECount;
               node0=Nstride0*N_PER_E*(i0+e_offset0)+Nstride1*N_PER_E*(NE1-1)+Nstride2*N_PER_E*(i2+e_offset2);
// printf("%d, %d ",k,k+1);   
	       index_t res=2*((i2+e_offset2)+NE2*(i0+e_offset0))+totalNECount;
               out->FaceElements->Id[k]=res;
               out->FaceElements->Tag[k]=BACKTAG;
               out->FaceElements->Owner[k]=myRank;
               out->FaceElements->Id[k+1]=res+1;
               out->FaceElements->Tag[k+1]=BACKTAG;
               out->FaceElements->Owner[k+1]=myRank;

	      index_t n2=node0+Nstride1;
	      index_t n6=node0+Nstride1+Nstride2;
	      index_t n7=node0+Nstride0+Nstride1+Nstride2;
	      index_t n3=node0+Nstride0+Nstride1;

	       if ((global_adjustment+i0+local_NE1-1+i2)%2==0)
	   {
	        out->FaceElements->Nodes[INDEX2(0,k,NN)]=n2;
	        out->FaceElements->Nodes[INDEX2(1,k,NN)]=n6;
	        out->FaceElements->Nodes[INDEX2(2,k,NN)]=n3;

	        out->FaceElements->Nodes[INDEX2(0,k+1,NN)]=n6;
	        out->FaceElements->Nodes[INDEX2(1,k+1,NN)]=n7;
	        out->FaceElements->Nodes[INDEX2(2,k+1,NN)]=n3;


	   }
	   else
	   {
	   // this form is rotated around the 0,2,4,6 face clockwise 90 degrees
	        out->FaceElements->Nodes[INDEX2(0,k,NN)]=n2;
	        out->FaceElements->Nodes[INDEX2(1,k,NN)]=n6;
	        out->FaceElements->Nodes[INDEX2(2,k,NN)]=n7;

	        out->FaceElements->Nodes[INDEX2(0,k+1,NN)]=n2;
	        out->FaceElements->Nodes[INDEX2(1,k+1,NN)]=n7;
	        out->FaceElements->Nodes[INDEX2(2,k+1,NN)]=n3;

	   }
/*printf("%d: %d, %d, %d\n", k,out->FaceElements->Nodes[INDEX2(0,k,NN)],  out->FaceElements->Nodes[INDEX2(1,k,NN)],
out->FaceElements->Nodes[INDEX2(2,k,NN)]);
printf("%d: %d, %d, %d\n", k+1,out->FaceElements->Nodes[INDEX2(0,k+1,NN)],  out->FaceElements->Nodes[INDEX2(1,k+1,NN)],
out->FaceElements->Nodes[INDEX2(2,k+1,NN)]);*/
             }
           }
           faceNECount+=2*local_NE0*local_NE2;
        }
        totalNECount+=2*NE0*NE2;
     }
  }
// printf("\n");
  if (Dudley_noError()) {
     /* add tag names */
     Dudley_Mesh_addTagMap(out,"top", TOPTAG);
     Dudley_Mesh_addTagMap(out,"bottom", BOTTOMTAG);
     Dudley_Mesh_addTagMap(out,"left", LEFTTAG);
     Dudley_Mesh_addTagMap(out,"right", RIGHTTAG);
     Dudley_Mesh_addTagMap(out,"front", FRONTTAG);
     Dudley_Mesh_addTagMap(out,"back", BACKTAG);
  }
  /* prepare mesh for further calculatuions:*/
  if (Dudley_noError()) {
         Dudley_Mesh_resolveNodeIds(out);
  }
  if (Dudley_noError()) {
         Dudley_Mesh_prepare(out, optimize);
  }

  if (!Dudley_noError()) {
      Dudley_Mesh_free(out);
  }
    /* free up memory */
  Dudley_ReferenceElementSet_dealloc(refPoints);
  Dudley_ReferenceElementSet_dealloc(refFaceElements);
  Dudley_ReferenceElementSet_dealloc(refElements);
  Paso_MPIInfo_free( mpi_info );  

  return out;
}


