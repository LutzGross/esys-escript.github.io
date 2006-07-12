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
                           dim_t numShape, dim_t numElements,index_t* nodes,
                           double* DSDv, dim_t numTest,double* DTDv,
                           double* dTdX, double* volume, index_t* element_id) {
     #define DIM 1
     #define LOCALDIM 1
     register int e,q,s;
     char error_msg[LenErrorMsg_MAX];
     #pragma omp parallel 
     {
       register double D,invD;
       double X0[numShape];
       #pragma omp for private(e,q,s) schedule(static) 
       for(e=0;e<numElements;e++){
           for (s=0;s<numShape; s++) X0[s]=coordinates[INDEX2(0,nodes[INDEX2(s,e,numShape)],DIM)];
           for (q=0;q<numQuad;q++) {
              D=0;
              for (s=0;s<numShape; s++) D+=X0[s]*DSDv[INDEX3(s,0,q,numShape,LOCALDIM)];
              if (D==0.) {
                  sprintf(error_msg,"Assemble_jacobeans_1D: element %d (id %d) has length zero.",e,element_id[e]);
                  Finley_setError(ZERO_DIVISION_ERROR,error_msg);
              } else {
                 invD=1./D;
                 for (s=0;s<numTest; s++) dTdX[INDEX4(s,0,q,e,numTest,LOCALDIM,numQuad)]=DTDv[INDEX3(s,0,q,numShape,LOCALDIM)]*invD;
		 volume[INDEX2(q,e,numQuad)]=ABS(D)*QuadWeights[q];
              }
           }
       }

     }
     #undef DIM 
     #undef LOCALDIM 

}
/**************************************************************/
/*                                                            */
/*  Jacobean 2D                                               */
/*                                                            */
void Assemble_jacobeans_2D(double* coordinates, dim_t numQuad,double* QuadWeights,
                           dim_t numShape, dim_t numElements,index_t* nodes,
                           double* DSDv, dim_t numTest,double* DTDv,
                           double* dTdX, double* volume, index_t* element_id) {
     #define DIM 2
     #define LOCALDIM 2
     register int e,q,s;
     char error_msg[LenErrorMsg_MAX];
     #pragma omp parallel 
     {
       register double dXdv00,dXdv10,dXdv01,dXdv11,
                       dvdX00,dvdX10,dvdX01,dvdX11, D,invD;
       double X0[numShape], X1[numShape];
       #pragma omp for private(e,q,s) schedule(static) 
       for(e=0;e<numElements;e++){
           for (s=0;s<numShape; s++) {
             X0[s]=coordinates[INDEX2(0,nodes[INDEX2(s,e,numShape)],DIM)];
             X1[s]=coordinates[INDEX2(1,nodes[INDEX2(s,e,numShape)],DIM)];
           }
           for (q=0;q<numQuad;q++) {
              dXdv00=0;
              dXdv10=0;
              dXdv01=0;
              dXdv11=0;
              for (s=0;s<numShape; s++) {
                 dXdv00+=X0[s]*DSDv[INDEX3(s,0,q,numShape,LOCALDIM)];
                 dXdv10+=X1[s]*DSDv[INDEX3(s,0,q,numShape,LOCALDIM)];
                 dXdv01+=X0[s]*DSDv[INDEX3(s,1,q,numShape,LOCALDIM)];
                 dXdv11+=X1[s]*DSDv[INDEX3(s,1,q,numShape,LOCALDIM)];
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
                   dTdX[INDEX4(s,0,q,e,numTest,LOCALDIM,numQuad)]=DTDv[INDEX3(s,0,q,numShape,LOCALDIM)]*dvdX00+DTDv[INDEX3(s,1,q,numShape,LOCALDIM)]*dvdX10;
                   dTdX[INDEX4(s,1,q,e,numTest,LOCALDIM,numQuad)]=DTDv[INDEX3(s,0,q,numShape,LOCALDIM)]*dvdX01+DTDv[INDEX3(s,1,q,numShape,LOCALDIM)]*dvdX11;
                 }
		 volume[INDEX2(q,e,numQuad)]=ABS(D)*QuadWeights[q];
              }
           }
       }

     }
     #undef DIM 
     #undef LOCALDIM 
}
/**************************************************************/
/*                                                            */
/*  Jacobean 3D                                               */
/*                                                            */
void Assemble_jacobeans_3D(double* coordinates, dim_t numQuad,double* QuadWeights,
                           dim_t numShape, dim_t numElements,index_t* nodes,
                           double* DSDv, dim_t numTest,double* DTDv,
                           double* dTdX, double* volume, index_t* element_id) {
     #define DIM 3
     #define LOCALDIM 3
     register int e,q,s;
     char error_msg[LenErrorMsg_MAX];
     #pragma omp parallel 
     {
       register double dXdv00,dXdv10,dXdv20,dXdv01,dXdv11,dXdv21,dXdv02,dXdv12,dXdv22,
                       dvdX00,dvdX10,dvdX20,dvdX01,dvdX11,dvdX21,dvdX02,dvdX12,dvdX22, D,invD;
       double X0[numShape], X1[numShape], X2[numShape];
       #pragma omp for private(e,q,s) schedule(static) 
       for(e=0;e<numElements;e++){
           for (s=0;s<numShape; s++) {
             X0[s]=coordinates[INDEX2(0,nodes[INDEX2(s,e,numShape)],DIM)];
             X1[s]=coordinates[INDEX2(1,nodes[INDEX2(s,e,numShape)],DIM)];
             X2[s]=coordinates[INDEX2(2,nodes[INDEX2(s,e,numShape)],DIM)];
           }
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
                 dXdv00+=X0[s]*DSDv[INDEX3(s,0,q,numShape,LOCALDIM)];
                 dXdv10+=X1[s]*DSDv[INDEX3(s,0,q,numShape,LOCALDIM)];
                 dXdv20+=X2[s]*DSDv[INDEX3(s,0,q,numShape,LOCALDIM)];
                 dXdv01+=X0[s]*DSDv[INDEX3(s,1,q,numShape,LOCALDIM)];
                 dXdv11+=X1[s]*DSDv[INDEX3(s,1,q,numShape,LOCALDIM)];
                 dXdv21+=X2[s]*DSDv[INDEX3(s,1,q,numShape,LOCALDIM)];
                 dXdv02+=X0[s]*DSDv[INDEX3(s,2,q,numShape,LOCALDIM)];
                 dXdv12+=X1[s]*DSDv[INDEX3(s,2,q,numShape,LOCALDIM)];
                 dXdv22+=X2[s]*DSDv[INDEX3(s,2,q,numShape,LOCALDIM)];
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
                   dTdX[INDEX4(s,0,q,e,numTest,LOCALDIM,numQuad)]=DTDv[INDEX3(s,0,q,numShape,LOCALDIM)]*dvdX00+DTDv[INDEX3(s,1,q,numShape,LOCALDIM)]*dvdX10+DTDv[INDEX3(s,2,q,numShape,LOCALDIM)]*dvdX20;
                   dTdX[INDEX4(s,1,q,e,numTest,LOCALDIM,numQuad)]=DTDv[INDEX3(s,0,q,numShape,LOCALDIM)]*dvdX01+DTDv[INDEX3(s,1,q,numShape,LOCALDIM)]*dvdX11+DTDv[INDEX3(s,2,q,numShape,LOCALDIM)]*dvdX21;
                   dTdX[INDEX4(s,2,q,e,numTest,LOCALDIM,numQuad)]=DTDv[INDEX3(s,0,q,numShape,LOCALDIM)]*dvdX02+DTDv[INDEX3(s,1,q,numShape,LOCALDIM)]*dvdX12+DTDv[INDEX3(s,2,q,numShape,LOCALDIM)]*dvdX22;
                 }
	         volume[INDEX2(q,e,numQuad)]=ABS(D)*QuadWeights[q];
              }
           }
       }

     }
     #undef DIM 
     #undef LOCALDIM 
}
/**************************************************************/
/*                                                            */
/*  Jacobean 2D manifold in 3D                                */
/*                                                            */
void Assemble_jacobeans_3D_M2D(double* coordinates, dim_t numQuad,double* QuadWeights,
                               dim_t numShape, dim_t numElements,index_t* nodes,
                               double* DSDv, dim_t numTest,double* DTDv,
                               double* dTdX, double* volume, index_t* element_id) {
     #define DIM 3
     #define LOCALDIM 2
     register int e,q,s;
     char error_msg[LenErrorMsg_MAX];
     #pragma omp parallel 
     {
       register double dXdv00,dXdv10,dXdv20,dXdv01,dXdv11,dXdv21,m00,m01,m11,
                       dvdX00,dvdX01,dvdX02,dvdX10,dvdX11,dvdX12,D,invD;
       double X0[numShape], X1[numShape], X2[numShape];
       #pragma omp for private(e,q,s) schedule(static) 
       for(e=0;e<numElements;e++){
           for (s=0;s<numShape; s++) {
             X0[s]=coordinates[INDEX2(0,nodes[INDEX2(s,e,numShape)],DIM)];
             X1[s]=coordinates[INDEX2(1,nodes[INDEX2(s,e,numShape)],DIM)];
             X2[s]=coordinates[INDEX2(2,nodes[INDEX2(s,e,numShape)],DIM)];
           }
           for (q=0;q<numQuad;q++) {
              dXdv00=0;
              dXdv10=0;
              dXdv20=0;
              dXdv01=0;
              dXdv11=0;
              dXdv21=0;
              for (s=0;s<numShape; s++) {
                 dXdv00+=X0[s]*DSDv[INDEX3(s,0,q,numShape,LOCALDIM)];
                 dXdv10+=X1[s]*DSDv[INDEX3(s,0,q,numShape,LOCALDIM)];
                 dXdv20+=X2[s]*DSDv[INDEX3(s,0,q,numShape,LOCALDIM)];
                 dXdv01+=X0[s]*DSDv[INDEX3(s,1,q,numShape,LOCALDIM)];
                 dXdv11+=X1[s]*DSDv[INDEX3(s,1,q,numShape,LOCALDIM)];
                 dXdv21+=X2[s]*DSDv[INDEX3(s,1,q,numShape,LOCALDIM)];
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
                   dTdX[INDEX4(s,0,q,e,numTest,LOCALDIM,numQuad)]=DTDv[INDEX3(s,0,q,numShape,LOCALDIM)]*dvdX00+DTDv[INDEX3(s,1,q,numShape,LOCALDIM)]*dvdX10;
                   dTdX[INDEX4(s,1,q,e,numTest,LOCALDIM,numQuad)]=DTDv[INDEX3(s,0,q,numShape,LOCALDIM)]*dvdX01+DTDv[INDEX3(s,1,q,numShape,LOCALDIM)]*dvdX11;
                   dTdX[INDEX4(s,2,q,e,numTest,LOCALDIM,numQuad)]=DTDv[INDEX3(s,0,q,numShape,LOCALDIM)]*dvdX02+DTDv[INDEX3(s,1,q,numShape,LOCALDIM)]*dvdX12;
                 }
	         volume[INDEX2(q,e,numQuad)]=sqrt(D)*QuadWeights[q];
              }
           }
       }

     }
     #undef DIM 
     #undef LOCALDIM 
}
/**************************************************************/
/*                                                            */
/*  Jacobean 1D manifold in 3D                                */
/*                                                            */
void Assemble_jacobeans_3D_M1D(double* coordinates, dim_t numQuad,double* QuadWeights,
                               dim_t numShape, dim_t numElements,index_t* nodes,
                               double* DSDv, dim_t numTest,double* DTDv,
                               double* dTdX, double* volume, index_t* element_id) {
     #define DIM 3
     #define LOCALDIM 1
     register int e,q,s;
     char error_msg[LenErrorMsg_MAX];
     #pragma omp parallel 
     {
       register double dXdv00,dXdv10,dXdv20,dvdX00,dvdX01,dvdX02,D,invD;
       double X0[numShape], X1[numShape], X2[numShape];
       #pragma omp for private(e,q,s) schedule(static) 
       for(e=0;e<numElements;e++){
           for (s=0;s<numShape; s++) {
             X0[s]=coordinates[INDEX2(0,nodes[INDEX2(s,e,numShape)],DIM)];
             X1[s]=coordinates[INDEX2(1,nodes[INDEX2(s,e,numShape)],DIM)];
             X2[s]=coordinates[INDEX2(2,nodes[INDEX2(s,e,numShape)],DIM)];
           }
           for (q=0;q<numQuad;q++) {
              dXdv00=0;
              dXdv10=0;
              dXdv20=0;
              for (s=0;s<numShape; s++) {
                 dXdv00+=X0[s]*DSDv[INDEX3(s,0,q,numShape,LOCALDIM)];
                 dXdv10+=X1[s]*DSDv[INDEX3(s,0,q,numShape,LOCALDIM)];
                 dXdv20+=X2[s]*DSDv[INDEX3(s,0,q,numShape,LOCALDIM)];
              }
              D=dXdv00*dXdv00+dXdv10*dXdv10+dXdv20*dXdv20;
              if (D==0.) {
                  sprintf(error_msg,"Assemble_jacobeans_3D_M1D: element %d (id %d) has length zero.",e,element_id[e]);
                  Finley_setError(ZERO_DIVISION_ERROR,error_msg);
              } else {
                 invD=1./D;
                 dvdX00=dXdv00*invD;
                 dvdX01=dXdv10*invD;
                 dvdX02=dXdv20*invD;
                 for (s=0;s<numTest; s++) {
                   dTdX[INDEX4(s,0,q,e,numTest,LOCALDIM,numQuad)]=DTDv[INDEX3(s,0,q,numShape,LOCALDIM)]*dvdX00;
                   dTdX[INDEX4(s,1,q,e,numTest,LOCALDIM,numQuad)]=DTDv[INDEX3(s,0,q,numShape,LOCALDIM)]*dvdX01;
                   dTdX[INDEX4(s,2,q,e,numTest,LOCALDIM,numQuad)]=DTDv[INDEX3(s,0,q,numShape,LOCALDIM)]*dvdX02;
                 }
	         volume[INDEX2(q,e,numQuad)]=sqrt(D)*QuadWeights[q];
              }
           }
       }

     }
     #undef DIM 
     #undef LOCALDIM 
}
/**************************************************************/
/*                                                            */
/*  Jacobean 1D manifold in 2D                                */
/*                                                            */
void Assemble_jacobeans_2D_M1D(double* coordinates, dim_t numQuad,double* QuadWeights,
                               dim_t numShape, dim_t numElements,index_t* nodes,
                               double* DSDv, dim_t numTest,double* DTDv,
                               double* dTdX, double* volume, index_t* element_id) {
     #define DIM 2
     #define LOCALDIM 1
     register int e,q,s;
     char error_msg[LenErrorMsg_MAX];
     #pragma omp parallel 
     {
       register double dXdv00,dXdv10,dvdX00,dvdX01,D,invD;
       double X0[numShape], X1[numShape];
       #pragma omp for private(e,q,s) schedule(static) 
       for(e=0;e<numElements;e++){
           for (s=0;s<numShape; s++) {
             X0[s]=coordinates[INDEX2(0,nodes[INDEX2(s,e,numShape)],DIM)];
             X1[s]=coordinates[INDEX2(1,nodes[INDEX2(s,e,numShape)],DIM)];
           }
           for (q=0;q<numQuad;q++) {
              dXdv00=0;
              dXdv10=0;
              for (s=0;s<numShape; s++) {
                 dXdv00+=X0[s]*DSDv[INDEX3(s,0,q,numShape,LOCALDIM)];
                 dXdv10+=X1[s]*DSDv[INDEX3(s,0,q,numShape,LOCALDIM)];
              }
              D=dXdv00*dXdv00+dXdv10*dXdv10;
              if (D==0.) {
                  sprintf(error_msg,"Assemble_jacobeans_2D_M1D: element %d (id %d) has length zero.",e,element_id[e]);
                  Finley_setError(ZERO_DIVISION_ERROR,error_msg);
              } else {
                 invD=1./D;
                 dvdX00=dXdv00*invD;
                 dvdX01=dXdv10*invD;
                 for (s=0;s<numTest; s++) {
                   dTdX[INDEX4(s,0,q,e,numTest,LOCALDIM,numQuad)]=DTDv[INDEX3(s,0,q,numShape,LOCALDIM)]*dvdX00;
                   dTdX[INDEX4(s,1,q,e,numTest,LOCALDIM,numQuad)]=DTDv[INDEX3(s,0,q,numShape,LOCALDIM)]*dvdX01;
                 }
	         volume[INDEX2(q,e,numQuad)]=sqrt(D)*QuadWeights[q];
              }
           }
       }

     }
     #undef DIM 
     #undef LOCALDIM 
}
/*
 * $Log$
 *
 */
