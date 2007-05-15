/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/


/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif



/**************************************************************/
/*                                                            */
/*  Jacobean 1D                                               */
/*                                                            */
void Assemble_jacobeans_1D(double* coordinates, dim_t numQuad,double* QuadWeights,
                           dim_t numShape, dim_t numElements, dim_t numNodes, index_t* nodes,
                           double* DSDv, dim_t numTest, double* DTDv, 
                           double* dTdX, double* volume, index_t* element_id) {
     #define DIM 1
     #define LOCDIM 1
     register int e,q,s;
     char error_msg[LenErrorMsg_MAX];
     #pragma omp parallel 
     {
       register double D,invD, X0_loc;
       #pragma omp for private(e,q,s,D,invD,X0_loc) schedule(static) 
       for(e=0;e<numElements;e++){
           for (q=0;q<numQuad;q++) {
              D=0;
              for (s=0;s<numShape; s++) {
                 X0_loc=coordinates[INDEX2(0,nodes[INDEX2(s,e,numNodes)],DIM)];
                 D+=X0_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
              }
              if (D==0.) {
                  sprintf(error_msg,"Assemble_jacobeans_1D: element %d (id %d) has length zero.",e,element_id[e]);
                  Finley_setError(ZERO_DIVISION_ERROR,error_msg);
              } else {
                 invD=1./D;
                 for (s=0;s<numTest; s++) dTdX[INDEX4(s,0,q,e,numTest,DIM,numQuad)]=DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*invD;
              }
	      volume[INDEX2(q,e,numQuad)]=ABS(D)*QuadWeights[q];
           }
       }

     }
     #undef DIM 
     #undef LOCDIM 

}
/**************************************************************/
/*                                                            */
/*  Jacobean 2D with area element                             */
/*                                                            */
void Assemble_jacobeans_2D(double* coordinates, dim_t numQuad,double* QuadWeights,
                           dim_t numShape, dim_t numElements, dim_t numNodes, index_t* nodes,
                           double* DSDv, dim_t numTest, double* DTDv, 
                           double* dTdX, double* volume, index_t* element_id) {
     #define DIM 2
     #define LOCDIM 2
     register int e,q,s;
     char error_msg[LenErrorMsg_MAX];
     #pragma omp parallel 
     {
       register double dXdv00,dXdv10,dXdv01,dXdv11,
                       dvdX00,dvdX10,dvdX01,dvdX11, D,invD,
                       X0_loc, X1_loc;
       #pragma omp for private(e,q,s,dXdv00,dXdv10,dXdv01,dXdv11,dvdX00,dvdX10,dvdX01,dvdX11, D,invD,X0_loc, X1_loc) schedule(static) 
       for(e=0;e<numElements;e++){
           for (q=0;q<numQuad;q++) {
              dXdv00=0;
              dXdv10=0;
              dXdv01=0;
              dXdv11=0;
              for (s=0;s<numShape; s++) {
                 X0_loc=coordinates[INDEX2(0,nodes[INDEX2(s,e,numNodes)],DIM)];
                 X1_loc=coordinates[INDEX2(1,nodes[INDEX2(s,e,numNodes)],DIM)];
                 dXdv00+=X0_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv10+=X1_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv01+=X0_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv11+=X1_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
              }
              D  =  dXdv00*dXdv11 - dXdv01*dXdv10;
              if (D==0.) {
                  sprintf(error_msg,"Assemble_jacobeans_2D: element %d (id %d) has area zero.",e,element_id[e]);
                  Finley_setError(ZERO_DIVISION_ERROR,error_msg);
              } else {
                 invD=1./D;
                 dvdX00= dXdv11*invD;
                 dvdX10=-dXdv10*invD;
                 dvdX01=-dXdv01*invD;
                 dvdX11= dXdv00*invD;

                 for (s=0;s<numTest; s++) {
                   dTdX[INDEX4(s,0,q,e,numTest,DIM,numQuad)]=DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10;
                   dTdX[INDEX4(s,1,q,e,numTest,DIM,numQuad)]=DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11;
                 }
              }
              volume[INDEX2(q,e,numQuad)]=ABS(D)*QuadWeights[q];
           }
       }

     }
     #undef DIM 
     #undef LOCDIM 
}
/**************************************************************/
/*                                                            */
/*  Jacobean 1D manifold in 2D and 1D elements                */
/*                                                            */
void Assemble_jacobeans_2D_M1D_E1D(double* coordinates, dim_t numQuad,double* QuadWeights,
                                   dim_t numShape, dim_t numElements, dim_t numNodes, index_t* nodes,
                                   double* DSDv, dim_t numTest, double* DTDv, 
                                   double* dTdX, double* volume, index_t* element_id) {
     #define DIM 2
     #define LOCDIM 1
     register int e,q,s;
     char error_msg[LenErrorMsg_MAX];
     #pragma omp parallel 
     {
       register double dXdv00,dXdv10,dvdX00,dvdX01,D,invD,
                       X0_loc, X1_loc;
       #pragma omp for private(e,q,s,dXdv00,dXdv10,dvdX00,dvdX01,D,invD,X0_loc, X1_loc) schedule(static) 
       for(e=0;e<numElements;e++){
           for (q=0;q<numQuad;q++) {
              dXdv00=0;
              dXdv10=0;
              for (s=0;s<numShape; s++) {
                 X0_loc=coordinates[INDEX2(0,nodes[INDEX2(s,e,numNodes)],DIM)];
                 X1_loc=coordinates[INDEX2(1,nodes[INDEX2(s,e,numNodes)],DIM)];
                 dXdv00+=X0_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv10+=X1_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
              }
              D=dXdv00*dXdv00+dXdv10*dXdv10;
              if (D==0.) {
                  sprintf(error_msg,"Assemble_jacobeans_2D_M1D_E1D: element %d (id %d) has length zero.",e,element_id[e]);
                  Finley_setError(ZERO_DIVISION_ERROR,error_msg);
              } else {
                 invD=1./D;
                 dvdX00=dXdv00*invD;
                 dvdX01=dXdv10*invD;
                 for (s=0;s<numTest; s++) {
                   dTdX[INDEX4(s,0,q,e,numTest,DIM,numQuad)]=DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00;
                   dTdX[INDEX4(s,1,q,e,numTest,DIM,numQuad)]=DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01;
                 }
	         volume[INDEX2(q,e,numQuad)]=sqrt(D)*QuadWeights[q];
              }
           }
       }

     }
     #undef DIM 
     #undef LOCDIM 
}
/**************************************************************/
/*                                                            */
/*  Jacobean 1D manifold in 2D and 1D elements woth contact   */
/*                                                            */
void Assemble_jacobeans_2D_M1D_E1D_C(double* coordinates, dim_t numQuad,double* QuadWeights,
                                   dim_t numShape, dim_t numElements, dim_t numNodes, index_t* nodes,
                                   double* DSDv, dim_t numTest, double* DTDv, 
                                   double* dTdX, double* volume, index_t* element_id) {
     #define DIM 2
     #define LOCDIM 1
     register int e,q,s;
     char error_msg[LenErrorMsg_MAX];
     #pragma omp parallel 
     {
       register double dXdv00_0,dXdv10_0,dvdX00_0,dvdX01_0,D_0,invD_0,
                       dXdv00_1,dXdv10_1,dvdX00_1,dvdX01_1,D_1,invD_1,
                       X0_loc_0, X1_loc_0, X0_loc_1, X1_loc_1;
       #pragma omp for private(e,q,s,dXdv00_0,dXdv10_0,dvdX00_0,dvdX01_0,D_0,invD_0,dXdv00_1,dXdv10_1,dvdX00_1,dvdX01_1,D_1,invD_1,X0_loc_0, X1_loc_0, X0_loc_1, X1_loc_1) schedule(static) 
       for(e=0;e<numElements;e++){
           for (q=0;q<numQuad;q++) {
              dXdv00_0=0;
              dXdv10_0=0;
              dXdv00_1=0;
              dXdv10_1=0;
              for (s=0;s<numShape; s++) {
                 X0_loc_0=coordinates[INDEX2(0,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                 X1_loc_0=coordinates[INDEX2(1,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                 X0_loc_1=coordinates[INDEX2(0,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                 X1_loc_1=coordinates[INDEX2(1,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                 dXdv00_0+=X0_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv10_0+=X1_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv00_1+=X0_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv10_1+=X1_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
              }
              D_0=dXdv00_0*dXdv00_0+dXdv10_0*dXdv10_0;
              D_1=dXdv00_1*dXdv00_1+dXdv10_1*dXdv10_1;
              if (D_0 == 0.  || D_1 == 0.) {
                  sprintf(error_msg,"Assemble_jacobeans_2D_M1D_E1D: element %d (id %d) has length zero.",e,element_id[e]);
                  Finley_setError(ZERO_DIVISION_ERROR,error_msg);
              } else {
                 invD_0=1./D_0;
                 dvdX00_0=dXdv00_0*invD_0;
                 dvdX01_0=dXdv10_0*invD_0;
                 dvdX00_1=dXdv00_1*invD_1;
                 dvdX01_1=dXdv10_1*invD_1;
                 for (s=0;s<numTest; s++) {
                   dTdX[INDEX4(        s,0,q,e,2*numTest,DIM,numQuad)]=DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00_0;
                   dTdX[INDEX4(        s,1,q,e,2*numTest,DIM,numQuad)]=DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01_0;
                   dTdX[INDEX4(numTest+s,0,q,e,2*numTest,DIM,numQuad)]=DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00_1;
                   dTdX[INDEX4(numTest+s,1,q,e,2*numTest,DIM,numQuad)]=DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01_1;
                 }
	         volume[INDEX2(q,e,numQuad)]=(sqrt(D_0)+sqrt(D_1))/2.*QuadWeights[q];
              }
           }
       }

     }
     #undef DIM 
     #undef LOCDIM 
}
/**************************************************************/
/*                                                            */
/*  Jacobean 1D manifold in 2D and 2D elements                */
/*                                                            */
void Assemble_jacobeans_2D_M1D_E2D(double* coordinates, dim_t numQuad,double* QuadWeights,
                                   dim_t numShape, dim_t numElements, dim_t numNodes, index_t* nodes,
                                   double* DSDv, dim_t numTest,double* DTDv,
                                   double* dTdX, double* volume, index_t* element_id) {
     #define DIM 2
     #define LOCDIM 2
     register int e,q,s;
     char error_msg[LenErrorMsg_MAX];
     #pragma omp parallel 
     {
       register double dXdv00,dXdv10,dXdv01,dXdv11,
                       dvdX00,dvdX10,dvdX01,dvdX11, D,invD,
                       X0_loc, X1_loc;
       #pragma omp for private(e,q,s,dXdv00,dXdv10,dXdv01,dXdv11,dvdX00,dvdX10,dvdX01,dvdX11, D,invD,X0_loc, X1_loc) schedule(static) 
       for(e=0;e<numElements;e++){
           for (q=0;q<numQuad;q++) {
              dXdv00=0;
              dXdv10=0;
              dXdv01=0;
              dXdv11=0;
              for (s=0;s<numShape; s++) {
                 X0_loc=coordinates[INDEX2(0,nodes[INDEX2(s,e,numNodes)],DIM)];
                 X1_loc=coordinates[INDEX2(1,nodes[INDEX2(s,e,numNodes)],DIM)];
                 dXdv00+=X0_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv10+=X1_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv01+=X0_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv11+=X1_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
              }
              D  =  dXdv00*dXdv11 - dXdv01*dXdv10;
              if (D==0.) {
                  sprintf(error_msg,"Assemble_jacobeans_2D_E2D: element %d (id %d) has area zero.",e,element_id[e]);
                  Finley_setError(ZERO_DIVISION_ERROR,error_msg);
              } else {
                 invD=1./D;
                 dvdX00= dXdv11*invD;
                 dvdX10=-dXdv10*invD;
                 dvdX01=-dXdv01*invD;
                 dvdX11= dXdv00*invD;

                 for (s=0;s<numTest; s++) {
                   dTdX[INDEX4(s,0,q,e,numTest,DIM,numQuad)]=DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10;
                   dTdX[INDEX4(s,1,q,e,numTest,DIM,numQuad)]=DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11;
                 }
              }
	      volume[INDEX2(q,e,numQuad)]=sqrt(dXdv00*dXdv00+dXdv10*dXdv10)*QuadWeights[q];
           }
       }

     }
     #undef DIM 
     #undef LOCDIM 
}
/**************************************************************/
/*                                                            */
/*  Jacobean 1D manifold in 2D and 2D elements with contact   */
/*                                                            */
void Assemble_jacobeans_2D_M1D_E2D_C(double* coordinates, dim_t numQuad,double* QuadWeights,
                                     dim_t numShape, dim_t numElements, dim_t numNodes, index_t* nodes,
                                     double* DSDv, dim_t numTest,double* DTDv,
                                     double* dTdX, double* volume, index_t* element_id) {
     #define DIM 2
     #define LOCDIM 2
     register int e,q,s;
     char error_msg[LenErrorMsg_MAX];
     #pragma omp parallel 
     {
       register double dXdv00_0,dXdv10_0,dXdv01_0,dXdv11_0,dvdX00_0,dvdX10_0,dvdX01_0,dvdX11_0, D_0,invD_0,
                       dXdv00_1,dXdv10_1,dXdv01_1,dXdv11_1,dvdX00_1,dvdX10_1,dvdX01_1,dvdX11_1, D_1,invD_1,
                       X0_loc_0, X1_loc_0, X0_loc_1, X1_loc_1;
       #pragma omp for private(e,q,s,dXdv00_0,dXdv10_0,dXdv01_0,dXdv11_0,dvdX00_0,dvdX10_0,dvdX01_0,dvdX11_0, D_0,invD_0,dXdv00_1,dXdv10_1,dXdv01_1,dXdv11_1,dvdX00_1,dvdX10_1,dvdX01_1,dvdX11_1, D_1,invD_1,X0_loc_0, X1_loc_0, X0_loc_1, X1_loc_1) schedule(static) 
       for(e=0;e<numElements;e++){
           for (q=0;q<numQuad;q++) {
              dXdv00_0=0;
              dXdv10_0=0;
              dXdv01_0=0;
              dXdv11_0=0;
              dXdv00_1=0;
              dXdv10_1=0;
              dXdv01_1=0;
              dXdv11_1=0;
              for (s=0;s<numShape; s++) {
                 X0_loc_0=coordinates[INDEX2(0,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                 X1_loc_0=coordinates[INDEX2(1,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                 X0_loc_1=coordinates[INDEX2(0,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                 X1_loc_1=coordinates[INDEX2(1,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                 dXdv00_0+=X0_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv10_0+=X1_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv01_0+=X0_loc_0*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv11_0+=X1_loc_0*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv00_1+=X0_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv10_1+=X1_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv01_1+=X0_loc_1*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv11_1+=X1_loc_1*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
              }
              D_0  =  dXdv00_0*dXdv11_0 - dXdv01_0*dXdv10_0;
              D_1  =  dXdv00_1*dXdv11_1 - dXdv01_1*dXdv10_1;
              if ( (D_0 ==0.) || (D_1 ==0.) ) {
                  sprintf(error_msg,"Assemble_jacobeans_2D_E2D_C: element %d (id %d) has area zero.",e,element_id[e]);
                  Finley_setError(ZERO_DIVISION_ERROR,error_msg);
              } else {
                 invD_0=1./D_0;
                 dvdX00_0= dXdv11_0*invD_0;
                 dvdX10_0=-dXdv10_0*invD_0;
                 dvdX01_0=-dXdv01_0*invD_0;
                 dvdX11_0= dXdv00_0*invD_0;
                 invD_1=1./D_1;
                 dvdX00_1= dXdv11_1*invD_1;
                 dvdX10_1=-dXdv10_1*invD_1;
                 dvdX01_1=-dXdv01_1*invD_1;
                 dvdX11_1= dXdv00_1*invD_1;
                 for (s=0;s<numTest; s++) {
                   dTdX[INDEX4(        s,0,q,e,2*numTest,DIM,numQuad)]=
                                                 DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00_0+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10_0;
                   dTdX[INDEX4(        s,1,q,e,2*numTest,DIM,numQuad)]=
                                                 DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01_0+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11_0;
                   dTdX[INDEX4(numTest+s,0,q,e,2*numTest,DIM,numQuad)]=
                                                 DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00_1+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10_1;
                   dTdX[INDEX4(numTest+s,1,q,e,2*numTest,DIM,numQuad)]=
                                                 DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01_1+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11_1;
                 }
              }
	      volume[INDEX2(q,e,numQuad)]=(sqrt(dXdv00_0*dXdv00_0+dXdv10_0*dXdv10_0)+sqrt(dXdv00_1*dXdv00_1+dXdv10_1*dXdv10_1))/2.*QuadWeights[q];
           }
       }

     }
     #undef DIM 
     #undef LOCDIM 
}
/**************************************************************/
/*                                                            */
/*  Jacobean 3D                                               */
/*                                                            */
void Assemble_jacobeans_3D(double* coordinates, dim_t numQuad,double* QuadWeights,
                           dim_t numShape, dim_t numElements, dim_t numNodes, index_t* nodes,
                           double* DSDv, dim_t numTest, double* DTDv, 
                           double* dTdX, double* volume, index_t* element_id) {
     #define DIM 3
     #define LOCDIM 3
     int e,q,s;
     char error_msg[LenErrorMsg_MAX];
     #pragma omp parallel 
     {
       register double dXdv00,dXdv10,dXdv20,dXdv01,dXdv11,dXdv21,dXdv02,dXdv12,dXdv22,
                       dvdX00,dvdX10,dvdX20,dvdX01,dvdX11,dvdX21,dvdX02,dvdX12,dvdX22, D,invD,
                       X0_loc,X1_loc,X2_loc;

      #pragma omp for private(e,q,s,dXdv00,dXdv10,dXdv20,dXdv01,dXdv11,dXdv21,dXdv02,dXdv12,dXdv22,dvdX00,dvdX10,dvdX20,dvdX01,dvdX11,dvdX21,dvdX02,dvdX12,dvdX22,D,invD,X0_loc,X1_loc,X2_loc) schedule(static)
       for(e=0;e<numElements;e++){
           for (q=0;q<numQuad;q++) {
              dXdv00=0;
              dXdv10=0;
              dXdv20=0;
              dXdv01=0;
              dXdv11=0;
              dXdv21=0;
              dXdv02=0;
              dXdv12=0;
              dXdv22=0;
              for (s=0;s<numShape; s++) {
                 X0_loc=coordinates[INDEX2(0,nodes[INDEX2(s,e,numNodes)],DIM)];
                 X1_loc=coordinates[INDEX2(1,nodes[INDEX2(s,e,numNodes)],DIM)];
                 X2_loc=coordinates[INDEX2(2,nodes[INDEX2(s,e,numNodes)],DIM)];
                 dXdv00+=X0_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv10+=X1_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv20+=X2_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv01+=X0_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv11+=X1_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv21+=X2_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv02+=X0_loc*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                 dXdv12+=X1_loc*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                 dXdv22+=X2_loc*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
              }
              D  =  dXdv00*(dXdv11*dXdv22-dXdv12*dXdv21)+ dXdv01*(dXdv20*dXdv12-dXdv10*dXdv22)+dXdv02*(dXdv10*dXdv21-dXdv20*dXdv11);
              if (D==0.) {
                  sprintf(error_msg,"Assemble_jacobeans_3D: element %d (id %d) has volume zero.",e,element_id[e]);
                  Finley_setError(ZERO_DIVISION_ERROR,error_msg);
              } else {
                 invD=1./D;
                 dvdX00=(dXdv11*dXdv22-dXdv12*dXdv21)*invD;
                 dvdX10=(dXdv20*dXdv12-dXdv10*dXdv22)*invD;
                 dvdX20=(dXdv10*dXdv21-dXdv20*dXdv11)*invD;
                 dvdX01=(dXdv02*dXdv21-dXdv01*dXdv22)*invD;
                 dvdX11=(dXdv00*dXdv22-dXdv20*dXdv02)*invD;
                 dvdX21=(dXdv01*dXdv20-dXdv00*dXdv21)*invD;
                 dvdX02=(dXdv01*dXdv12-dXdv02*dXdv11)*invD;
                 dvdX12=(dXdv02*dXdv10-dXdv00*dXdv12)*invD;
                 dvdX22=(dXdv00*dXdv11-dXdv01*dXdv10)*invD;

                 for (s=0;s<numTest; s++) {
                   dTdX[INDEX4(s,0,q,e,numTest,DIM,numQuad)]=DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10+DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX20;
                   dTdX[INDEX4(s,1,q,e,numTest,DIM,numQuad)]=DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11+DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX21;
                   dTdX[INDEX4(s,2,q,e,numTest,DIM,numQuad)]=DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX02+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX12+DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX22;
                 }
	         volume[INDEX2(q,e,numQuad)]=ABS(D)*QuadWeights[q];
              }
           }
       }

     }
     #undef DIM 
     #undef LOCDIM 
}
/**************************************************************/
/*                                                            */
/*  Jacobean 2D manifold in 3D with 3D elements               */
/*                                                            */
void Assemble_jacobeans_3D_M2D_E3D(double* coordinates, dim_t numQuad,double* QuadWeights,
                                   dim_t numShape, dim_t numElements, dim_t numNodes, index_t* nodes,
                                   double* DSDv, dim_t numTest,double* DTDv,
                                   double* dTdX, double* volume, index_t* element_id) {
     #define DIM 3
     #define LOCDIM 3
     register int e,q,s;
     char error_msg[LenErrorMsg_MAX];
     #pragma omp parallel 
     {
       register double dXdv00,dXdv10,dXdv20,dXdv01,dXdv11,dXdv21,dXdv02,dXdv12,dXdv22, m0, m1, m2,
                       dvdX00,dvdX10,dvdX20,dvdX01,dvdX11,dvdX21,dvdX02,dvdX12,dvdX22, D,invD,
                       X0_loc, X1_loc, X2_loc;
       #pragma omp for private(e,q,s,dXdv00,dXdv10,dXdv20,dXdv01,dXdv11,dXdv21,dXdv02,dXdv12,dXdv22, m0, m1, m2,dvdX00,dvdX10,dvdX20,dvdX01,dvdX11,dvdX21,dvdX02,dvdX12,dvdX22, D,invD,X0_loc, X1_loc, X2_loc) schedule(static) 
       for(e=0;e<numElements;e++){
           for (q=0;q<numQuad;q++) {
              dXdv00=0;
              dXdv10=0;
              dXdv20=0;
              dXdv01=0;
              dXdv11=0;
              dXdv21=0;
              dXdv02=0;
              dXdv12=0;
              dXdv22=0;
              for (s=0;s<numShape; s++) {
                 X0_loc=coordinates[INDEX2(0,nodes[INDEX2(s,e,numNodes)],DIM)];
                 X1_loc=coordinates[INDEX2(1,nodes[INDEX2(s,e,numNodes)],DIM)];
                 X2_loc=coordinates[INDEX2(2,nodes[INDEX2(s,e,numNodes)],DIM)];
                 dXdv00+=X0_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv10+=X1_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv20+=X2_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv01+=X0_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv11+=X1_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv21+=X2_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv02+=X0_loc*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                 dXdv12+=X1_loc*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                 dXdv22+=X2_loc*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
              }
              D  =  dXdv00*(dXdv11*dXdv22-dXdv12*dXdv21)+ dXdv01*(dXdv20*dXdv12-dXdv10*dXdv22)+dXdv02*(dXdv10*dXdv21-dXdv20*dXdv11);
              if (D==0.) {
                  sprintf(error_msg,"Assemble_jacobeans_3D: element %d (id %d) has volume zero.",e,element_id[e]);
                  Finley_setError(ZERO_DIVISION_ERROR,error_msg);
              } else {
                 invD=1./D;
                 dvdX00=(dXdv11*dXdv22-dXdv12*dXdv21)*invD;
                 dvdX10=(dXdv20*dXdv12-dXdv10*dXdv22)*invD;
                 dvdX20=(dXdv10*dXdv21-dXdv20*dXdv11)*invD;
                 dvdX01=(dXdv02*dXdv21-dXdv01*dXdv22)*invD;
                 dvdX11=(dXdv00*dXdv22-dXdv20*dXdv02)*invD;
                 dvdX21=(dXdv01*dXdv20-dXdv00*dXdv21)*invD;
                 dvdX02=(dXdv01*dXdv12-dXdv02*dXdv11)*invD;
                 dvdX12=(dXdv02*dXdv10-dXdv00*dXdv12)*invD;
                 dvdX22=(dXdv00*dXdv11-dXdv01*dXdv10)*invD;

                 for (s=0;s<numTest; s++) {
                   dTdX[INDEX4(s,0,q,e,numTest,DIM,numQuad)]=
                       DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10+DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX20;
                   dTdX[INDEX4(s,1,q,e,numTest,DIM,numQuad)]=
                       DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11+DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX21;
                   dTdX[INDEX4(s,2,q,e,numTest,DIM,numQuad)]=
                       DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX02+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX12+DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX22;
                 }
              }
              m0=dXdv10*dXdv21-dXdv20*dXdv11;
              m1=dXdv20*dXdv01-dXdv00*dXdv21;
              m2=dXdv00*dXdv11-dXdv10*dXdv01;
	      volume[INDEX2(q,e,numQuad)]=sqrt(m0*m0+m1*m1+m2*m2)*QuadWeights[q];
           }
       }

     }
     #undef DIM 
     #undef LOCDIM 
}
/**************************************************************/
/*                                                            */
/*  Jacobean 2D manifold in 3D with 3D elements on contact    */
/*                                                            */
void Assemble_jacobeans_3D_M2D_E3D_C(double* coordinates, dim_t numQuad,double* QuadWeights,
                                     dim_t numShape, dim_t numElements, dim_t numNodes, index_t* nodes,
                                     double* DSDv, dim_t numTest,double* DTDv,
                                     double* dTdX, double* volume, index_t* element_id) {
     #define DIM 3
     #define LOCDIM 3
     register int e,q,s;
     char error_msg[LenErrorMsg_MAX];
     #pragma omp parallel 
     {
       register double dXdv00_0,dXdv10_0,dXdv20_0,dXdv01_0,dXdv11_0,dXdv21_0,dXdv02_0,dXdv12_0,dXdv22_0, m0_0, m1_0, m2_0,
                       dvdX00_0,dvdX10_0,dvdX20_0,dvdX01_0,dvdX11_0,dvdX21_0,dvdX02_0,dvdX12_0,dvdX22_0, D_0,invD_0,
                       dXdv00_1,dXdv10_1,dXdv20_1,dXdv01_1,dXdv11_1,dXdv21_1,dXdv02_1,dXdv12_1,dXdv22_1, m0_1, m1_1, m2_1,
                       dvdX00_1,dvdX10_1,dvdX20_1,dvdX01_1,dvdX11_1,dvdX21_1,dvdX02_1,dvdX12_1,dvdX22_1, D_1,invD_1,
                       X0_loc_0, X1_loc_0, X2_loc_0, X0_loc_1, X1_loc_1, X2_loc_1;
       #pragma omp for private(e,q,s,dXdv00_0,dXdv10_0,dXdv20_0,dXdv01_0,dXdv11_0,dXdv21_0,dXdv02_0,dXdv12_0,dXdv22_0, m0_0, m1_0, m2_0,dvdX00_0,dvdX10_0,dvdX20_0,dvdX01_0,dvdX11_0,dvdX21_0,dvdX02_0,dvdX12_0,dvdX22_0, D_0,invD_0,dXdv00_1,dXdv10_1,dXdv20_1,dXdv01_1,dXdv11_1,dXdv21_1,dXdv02_1,dXdv12_1,dXdv22_1, m0_1, m1_1, m2_1,dvdX00_1,dvdX10_1,dvdX20_1,dvdX01_1,dvdX11_1,dvdX21_1,dvdX02_1,dvdX12_1,dvdX22_1, D_1,invD_1,X0_loc_0, X1_loc_0, X2_loc_0, X0_loc_1, X1_loc_1, X2_loc_1) schedule(static) 
       for(e=0;e<numElements;e++){
           for (q=0;q<numQuad;q++) {
              dXdv00_0=0;
              dXdv10_0=0;
              dXdv20_0=0;
              dXdv01_0=0;
              dXdv11_0=0;
              dXdv21_0=0;
              dXdv02_0=0;
              dXdv12_0=0;
              dXdv22_0=0;
              dXdv00_1=0;
              dXdv10_1=0;
              dXdv20_1=0;
              dXdv01_1=0;
              dXdv11_1=0;
              dXdv21_1=0;
              dXdv02_1=0;
              dXdv12_1=0;
              dXdv22_1=0;
              for (s=0;s<numShape; s++) {
                 X0_loc_0=coordinates[INDEX2(0,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                 X1_loc_0=coordinates[INDEX2(1,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                 X2_loc_0=coordinates[INDEX2(2,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                 X0_loc_1=coordinates[INDEX2(0,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                 X1_loc_1=coordinates[INDEX2(1,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                 X2_loc_1=coordinates[INDEX2(2,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                 dXdv00_0+=X0_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv10_0+=X1_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv20_0+=X2_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv01_0+=X0_loc_0*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv11_0+=X1_loc_0*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv21_0+=X2_loc_0*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv02_0+=X0_loc_0*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                 dXdv12_0+=X1_loc_0*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                 dXdv22_0+=X2_loc_0*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                 dXdv00_1+=X0_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv10_1+=X1_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv20_1+=X2_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv01_1+=X0_loc_1*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv11_1+=X1_loc_1*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv21_1+=X2_loc_1*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv02_1+=X0_loc_1*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                 dXdv12_1+=X1_loc_1*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                 dXdv22_1+=X2_loc_1*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
              }

              D_0=dXdv00_0*(dXdv11_0*dXdv22_0-dXdv12_0*dXdv21_0)+dXdv01_0*(dXdv20_0*dXdv12_0-dXdv10_0*dXdv22_0)+dXdv02_0*(dXdv10_0*dXdv21_0-dXdv20_0*dXdv11_0);
              D_1=dXdv00_1*(dXdv11_1*dXdv22_1-dXdv12_1*dXdv21_1)+dXdv01_1*(dXdv20_1*dXdv12_1-dXdv10_1*dXdv22_1)+dXdv02_1*(dXdv10_1*dXdv21_1-dXdv20_1*dXdv11_1);
              if ( (D_0==0.) || (D_1 == 0.)) {
                  sprintf(error_msg,"Assemble_jacobeans_3D_C: element %d (id %d) has volume zero.",e,element_id[e]);
                  Finley_setError(ZERO_DIVISION_ERROR,error_msg);
              } else {
                 invD_0=1./D_0;
                 dvdX00_0=(dXdv11_0*dXdv22_0-dXdv12_0*dXdv21_0)*invD_0;
                 dvdX10_0=(dXdv20_0*dXdv12_0-dXdv10_0*dXdv22_0)*invD_0;
                 dvdX20_0=(dXdv10_0*dXdv21_0-dXdv20_0*dXdv11_0)*invD_0;
                 dvdX01_0=(dXdv02_0*dXdv21_0-dXdv01_0*dXdv22_0)*invD_0;
                 dvdX11_0=(dXdv00_0*dXdv22_0-dXdv20_0*dXdv02_0)*invD_0;
                 dvdX21_0=(dXdv01_0*dXdv20_0-dXdv00_0*dXdv21_0)*invD_0;
                 dvdX02_0=(dXdv01_0*dXdv12_0-dXdv02_0*dXdv11_0)*invD_0;
                 dvdX12_0=(dXdv02_0*dXdv10_0-dXdv00_0*dXdv12_0)*invD_0;
                 dvdX22_0=(dXdv00_0*dXdv11_0-dXdv01_0*dXdv10_0)*invD_0;
                 invD_1=1./D_1;
                 dvdX00_1=(dXdv11_1*dXdv22_1-dXdv12_1*dXdv21_1)*invD_1;
                 dvdX10_1=(dXdv20_1*dXdv12_1-dXdv10_1*dXdv22_1)*invD_1;
                 dvdX20_1=(dXdv10_1*dXdv21_1-dXdv20_1*dXdv11_1)*invD_1;
                 dvdX01_1=(dXdv02_1*dXdv21_1-dXdv01_1*dXdv22_1)*invD_1;
                 dvdX11_1=(dXdv00_1*dXdv22_1-dXdv20_1*dXdv02_1)*invD_1;
                 dvdX21_1=(dXdv01_1*dXdv20_1-dXdv00_1*dXdv21_1)*invD_1;
                 dvdX02_1=(dXdv01_1*dXdv12_1-dXdv02_1*dXdv11_1)*invD_1;
                 dvdX12_1=(dXdv02_1*dXdv10_1-dXdv00_1*dXdv12_1)*invD_1;
                 dvdX22_1=(dXdv00_1*dXdv11_1-dXdv01_1*dXdv10_1)*invD_1;

                 for (s=0;s<numTest; s++) {
                   dTdX[INDEX4(        s,0,q,e,2*numTest,DIM,numQuad)]=
                      DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00_0+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10_0+DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX20_0;
                   dTdX[INDEX4(        s,1,q,e,2*numTest,DIM,numQuad)]=
                      DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01_0+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11_0+DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX21_0;
                   dTdX[INDEX4(        s,2,q,e,2*numTest,DIM,numQuad)]=
                      DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX02_0+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX12_0+DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX22_0;
                   dTdX[INDEX4(numTest+s,0,q,e,2*numTest,DIM,numQuad)]=
                      DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00_1+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10_1+DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX20_1;
                   dTdX[INDEX4(numTest+s,1,q,e,2*numTest,DIM,numQuad)]=
                      DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01_1+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11_1+DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX21_1;
                   dTdX[INDEX4(numTest+s,2,q,e,2*numTest,DIM,numQuad)]=
                      DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX02_1+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX12_1+DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX22_1;
                 }
              }
              m0_0=dXdv10_0*dXdv21_0-dXdv20_0*dXdv11_0;
              m1_0=dXdv20_0*dXdv01_0-dXdv00_0*dXdv21_0;
              m2_0=dXdv00_0*dXdv11_0-dXdv10_0*dXdv01_0;
              m0_1=dXdv10_1*dXdv21_1-dXdv20_1*dXdv11_1;
              m1_1=dXdv20_1*dXdv01_1-dXdv00_1*dXdv21_1;
              m2_1=dXdv00_1*dXdv11_1-dXdv10_1*dXdv01_1;
	      volume[INDEX2(q,e,numQuad)]=(sqrt(m0_0*m0_0+m1_0*m1_0+m2_0*m2_0)+sqrt(m0_1*m0_1+m1_1*m1_1+m2_1*m2_1))/2.*QuadWeights[q];
           }
       }

     }
     #undef DIM 
     #undef LOCDIM 
}
/**************************************************************/
/*                                                            */
/*  Jacobean 2D manifold in 3D with 2D elements               */
/*                                                            */
void Assemble_jacobeans_3D_M2D_E2D(double* coordinates, dim_t numQuad,double* QuadWeights,
                                   dim_t numShape, dim_t numElements, dim_t numNodes, index_t* nodes,
                                   double* DSDv, dim_t numTest,double* DTDv,
                                   double* dTdX, double* volume, index_t* element_id) {
     #define DIM 3
     #define LOCDIM 2
     register int e,q,s;
     char error_msg[LenErrorMsg_MAX];
     #pragma omp parallel 
     {
       register double dXdv00,dXdv10,dXdv20,dXdv01,dXdv11,dXdv21,m00,m01,m11,
                       dvdX00,dvdX01,dvdX02,dvdX10,dvdX11,dvdX12,D,invD,
                       X0_loc, X1_loc, X2_loc;
       #pragma omp for private(e,q,s,dXdv00,dXdv10,dXdv20,dXdv01,dXdv11,dXdv21,m00,m01,m11,dvdX00,dvdX01,dvdX02,dvdX10,dvdX11,dvdX12,D,invD, X0_loc, X1_loc, X2_loc) schedule(static) 
       for(e=0;e<numElements;e++){
           for (q=0;q<numQuad;q++) {
              dXdv00=0;
              dXdv10=0;
              dXdv20=0;
              dXdv01=0;
              dXdv11=0;
              dXdv21=0;
              for (s=0;s<numShape; s++) {
                 X0_loc=coordinates[INDEX2(0,nodes[INDEX2(s,e,numNodes)],DIM)];
                 X1_loc=coordinates[INDEX2(1,nodes[INDEX2(s,e,numNodes)],DIM)];
                 X2_loc=coordinates[INDEX2(2,nodes[INDEX2(s,e,numNodes)],DIM)];
                 dXdv00+=X0_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv10+=X1_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv20+=X2_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv01+=X0_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv11+=X1_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv21+=X2_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
              }
              m00=dXdv00*dXdv00+dXdv10*dXdv10+dXdv20*dXdv20;
              m01=dXdv00*dXdv01+dXdv10*dXdv11+dXdv20*dXdv21;
              m11=dXdv01*dXdv01+dXdv11*dXdv11+dXdv21*dXdv21;
              D=m00*m11-m01*m01;
              if (D==0.) {
                  sprintf(error_msg,"Assemble_jacobeans_3D_M2D: element %d (id %d) has area zero.",e,element_id[e]);
                  Finley_setError(ZERO_DIVISION_ERROR,error_msg);
              } else {
                 invD=1./D;
                 dvdX00=( m00*dXdv00-m01*dXdv01)*invD;
                 dvdX01=( m00*dXdv10-m01*dXdv11)*invD;
                 dvdX02=( m00*dXdv20-m01*dXdv21)*invD;
                 dvdX10=(-m01*dXdv00+m11*dXdv01)*invD;
                 dvdX11=(-m01*dXdv10+m11*dXdv11)*invD;
                 dvdX12=(-m01*dXdv20+m11*dXdv21)*invD;
                 for (s=0;s<numTest; s++) {
                   dTdX[INDEX4(s,0,q,e,numTest,DIM,numQuad)]=DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10;
                   dTdX[INDEX4(s,1,q,e,numTest,DIM,numQuad)]=DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11;
                   dTdX[INDEX4(s,2,q,e,numTest,DIM,numQuad)]=DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX02+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX12;
                 }
	         volume[INDEX2(q,e,numQuad)]=sqrt(D)*QuadWeights[q];
              }
           }
       }

     }
     #undef DIM 
     #undef LOCDIM 
}
/**************************************************************/
/*                                                            */
/*  Jacobean 2D manifold in 3D with 2D elements  with contact */
/*                                                            */
void Assemble_jacobeans_3D_M2D_E2D_C(double* coordinates, dim_t numQuad,double* QuadWeights,
                                     dim_t numShape, dim_t numElements, dim_t numNodes, index_t* nodes,
                                     double* DSDv, dim_t numTest,double* DTDv,
                                     double* dTdX, double* volume, index_t* element_id) {
     #define DIM 3
     #define LOCDIM 2
     register int e,q,s;
     char error_msg[LenErrorMsg_MAX];
     #pragma omp parallel 
     {
       register double dXdv00_0,dXdv10_0,dXdv20_0,dXdv01_0,dXdv11_0,dXdv21_0,m00_0,m01_0,m11_0,
                       dvdX00_0,dvdX01_0,dvdX02_0,dvdX10_0,dvdX11_0,dvdX12_0,D_0,invD_0,
                       dXdv00_1,dXdv10_1,dXdv20_1,dXdv01_1,dXdv11_1,dXdv21_1,m00_1,m01_1,m11_1,
                       dvdX00_1,dvdX01_1,dvdX02_1,dvdX10_1,dvdX11_1,dvdX12_1,D_1,invD_1,
                       X0_loc_0, X1_loc_0, X2_loc_0, X0_loc_1, X1_loc_1, X2_loc_1;
       #pragma omp for private(e,q,s,dXdv00_0,dXdv10_0,dXdv20_0,dXdv01_0,dXdv11_0,dXdv21_0,m00_0,m01_0,m11_0,dvdX00_0,dvdX01_0,dvdX02_0,dvdX10_0,dvdX11_0,dvdX12_0,D_0,invD_0,dXdv00_1,dXdv10_1,dXdv20_1,dXdv01_1,dXdv11_1,dXdv21_1,m00_1,m01_1,m11_1,dvdX00_1,dvdX01_1,dvdX02_1,dvdX10_1,dvdX11_1,dvdX12_1,D_1,invD_1,X0_loc_0, X1_loc_0, X2_loc_0, X0_loc_1, X1_loc_1, X2_loc_1) schedule(static) 
       for(e=0;e<numElements;e++){
           for (q=0;q<numQuad;q++) {
              dXdv00_0=0;
              dXdv10_0=0;
              dXdv20_0=0;
              dXdv01_0=0;
              dXdv11_0=0;
              dXdv21_0=0;
              dXdv00_1=0;
              dXdv10_1=0;
              dXdv20_1=0;
              dXdv01_1=0;
              dXdv11_1=0;
              dXdv21_1=0;
              for (s=0;s<numShape; s++) {
                 X0_loc_0=coordinates[INDEX2(0,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                 X1_loc_0=coordinates[INDEX2(1,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                 X2_loc_0=coordinates[INDEX2(2,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                 X0_loc_1=coordinates[INDEX2(0,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                 X1_loc_1=coordinates[INDEX2(1,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                 X2_loc_1=coordinates[INDEX2(2,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                 dXdv00_0+=X0_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv10_0+=X1_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv20_0+=X2_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv01_0+=X0_loc_0*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv11_0+=X1_loc_0*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv21_0+=X2_loc_0*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv00_1+=X0_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv10_1+=X1_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv20_1+=X2_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                 dXdv01_1+=X0_loc_1*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv11_1+=X1_loc_1*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                 dXdv21_1+=X2_loc_1*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
              }
              m00_0=dXdv00_0*dXdv00_0+dXdv10_0*dXdv10_0+dXdv20_0*dXdv20_0;
              m01_0=dXdv00_0*dXdv01_0+dXdv10_0*dXdv11_0+dXdv20_0*dXdv21_0;
              m11_0=dXdv01_0*dXdv01_0+dXdv11_0*dXdv11_0+dXdv21_0*dXdv21_0;
              D_0=m00_0*m11_0-m01_0*m01_0;
              m00_1=dXdv00_1*dXdv00_1+dXdv10_1*dXdv10_1+dXdv20_1*dXdv20_1;
              m01_1=dXdv00_1*dXdv01_1+dXdv10_1*dXdv11_1+dXdv20_1*dXdv21_1;
              m11_1=dXdv01_1*dXdv01_1+dXdv11_1*dXdv11_1+dXdv21_1*dXdv21_1;
              D_1=m00_1*m11_1-m01_1*m01_1;
              if ( (D_0==0.) || (D_1 == 0.) ) {
                  sprintf(error_msg,"Assemble_jacobeans_3D_M2D: element %d (id %d) has area zero.",e,element_id[e]);
                  Finley_setError(ZERO_DIVISION_ERROR,error_msg);
              } else {
                 invD_0=1./D_0;
                 dvdX00_0=( m00_0*dXdv00_0-m01_0*dXdv01_0)*invD_0;
                 dvdX01_0=( m00_0*dXdv10_0-m01_0*dXdv11_0)*invD_0;
                 dvdX02_0=( m00_0*dXdv20_0-m01_0*dXdv21_0)*invD_0;
                 dvdX10_0=(-m01_0*dXdv00_0+m11_0*dXdv01_0)*invD_0;
                 dvdX11_0=(-m01_0*dXdv10_0+m11_0*dXdv11_0)*invD_0;
                 dvdX12_0=(-m01_0*dXdv20_0+m11_0*dXdv21_0)*invD_0;
                 invD_1=1./D_1;
                 dvdX00_1=( m00_1*dXdv00_1-m01_1*dXdv01_1)*invD_1;
                 dvdX01_1=( m00_1*dXdv10_1-m01_1*dXdv11_1)*invD_1;
                 dvdX02_1=( m00_1*dXdv20_1-m01_1*dXdv21_1)*invD_1;
                 dvdX10_1=(-m01_1*dXdv00_1+m11_1*dXdv01_1)*invD_1;
                 dvdX11_1=(-m01_1*dXdv10_1+m11_1*dXdv11_1)*invD_1;
                 dvdX12_1=(-m01_1*dXdv20_1+m11_1*dXdv21_1)*invD_1;
                 for (s=0;s<numTest; s++) {
                   dTdX[INDEX4(        s,0,q,e,2*numTest,DIM,numQuad)]=
                                 DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00_0+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10_0;
                   dTdX[INDEX4(        s,1,q,e,2*numTest,DIM,numQuad)]=
                                 DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01_0+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11_0;
                   dTdX[INDEX4(        s,2,q,e,2*numTest,DIM,numQuad)]=
                                 DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX02_0+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX12_0;
                   dTdX[INDEX4(numTest+s,0,q,e,2*numTest,DIM,numQuad)]=
                                 DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00_1+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10_1;
                   dTdX[INDEX4(numTest+s,1,q,e,2*numTest,DIM,numQuad)]=
                                 DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01_1+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11_1;
                   dTdX[INDEX4(numTest+s,2,q,e,2*numTest,DIM,numQuad)]=
                                 DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX02_1+DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX12_1;
                 }
	         volume[INDEX2(q,e,numQuad)]=(sqrt(D_0)+sqrt(D_1))/2.*QuadWeights[q];
              }
           }
       }

     }
     #undef DIM 
     #undef LOCDIM 
}
/*
 * $Log:$
 *
 */
