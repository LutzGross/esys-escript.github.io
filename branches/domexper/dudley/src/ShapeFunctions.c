
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

/*   Dudley: Shape functions */

/**************************************************************/

#include "ShapeFunctions.h"



Dudley_ShapeFunctionInfo Dudley_ShapeFunction_InfoList[]={
	{Point1Shape, "Point1", 0,  1, 1, 1,	Dudley_Shape_Point1 } ,
        {Line2Shape,  "Line2",  1,  2, 1, 2,	Dudley_Shape_Line2  } ,
	{Line3Shape,  "Line3",  1,  3, 2, 2,	Dudley_Shape_Line3  },
	{Line4Shape,  "Line4",  1,  4, 3, 2,	Dudley_Shape_Line4  },
	{Tri3Shape,   "Tri3",   2,  3, 1, 3,	Dudley_Shape_Tri3   },
	{Tri6Shape,   "Tri6",   2,  6, 2, 3,	Dudley_Shape_Tri6   },
	{Tri9Shape,   "Tri9",   2,  9, 3, 3,	Dudley_Shape_Tri9   },
	{Tri10Shape,  "Tri10",  2, 10, 3, 3,	Dudley_Shape_Tri10, },
	{Tet4Shape,   "Tet4",   3,  4, 1, 4,	Dudley_Shape_Tet4,  },
	{Tet10Shape,  "Tet10",  3, 10, 2, 4,	Dudley_Shape_Tet10, },
	{Tet16Shape,  "Tet16",  3, 16, 3, 4,	Dudley_Shape_Tet16, },
	{NoShape, "NoType", 0,  1, 1, 1,	Dudley_Shape_Point1  }
};


/******************************************************************************************************************************
   
    creates an evaluation of the ShapeFunction on the given quadrature scheme. 
    if the spatial dimension of the scheme and the shape functions don't match
    
    if QuadNodes==Null or QuadWeights==Null the shape functions method is used to generate a quadrature scheme with numQuasNodes
    nodes. otherwise its assumed that a quadraure scheme is given on these array and copy is created within the structure.

*/
Dudley_ShapeFunction* Dudley_ShapeFunction_alloc(Dudley_ShapeFunctionTypeId id,int numQuadDim, int numQuadNodes, double *QuadNodes, double *QuadWeights) {
	Dudley_ShapeFunction *out=NULL;
	int numDim, numShapes, i, q;
  
	numDim=Dudley_ShapeFunction_InfoList[id].numDim;
        numShapes=Dudley_ShapeFunction_InfoList[id].numShapes;
  
    if (numQuadDim>numDim) {    
	    Dudley_setError(VALUE_ERROR,"Dudley_ShapeFunction_alloc: spatial dimension of quadrature scheme is bigger then spatial dimension of shape function.");
	    return NULL;
    }
	
	/*  allocate the Dudley_ShapeFunction to be returned: */
  
	out=MEMALLOC(1,Dudley_ShapeFunction);
	if (Dudley_checkPtr(out)) return NULL;

  
	out->Type=Dudley_ShapeFunction_getInfo(id);
	out->numQuadNodes=numQuadNodes;
	out->QuadNodes=NULL;
	out->QuadWeights=NULL;
	out->S=NULL;
	out->dSdv=NULL;
	out->reference_counter=0;
  
	/*  allocate memory: */
  
	out->QuadNodes=MEMALLOC(numQuadNodes*numDim,double);
	out->QuadWeights=MEMALLOC(numQuadNodes,double);
	out->S=MEMALLOC(numShapes*numQuadNodes,double);
	out->dSdv=MEMALLOC(numShapes*numDim*numQuadNodes,double);
	if ( Dudley_checkPtr(out->QuadNodes) || Dudley_checkPtr(out->QuadWeights) || Dudley_checkPtr(out->S) || Dudley_checkPtr(out->dSdv) ) {
         Dudley_ShapeFunction_dealloc(out);
         return NULL;
	}
  
	/*  set the quadrature nodes (missing values are filled with 0): */

    for (q=0;q<numQuadNodes;q++) {
       for (i=0;i<numQuadDim;i++) 	 out->QuadNodes[INDEX2(i,q,numDim)]=QuadNodes[INDEX2(i,q,numQuadDim)];
       for (i=numQuadDim;i<numDim;i++) out->QuadNodes[INDEX2(i,q,numDim)]=0;
       out->QuadWeights[q]=QuadWeights[q];
    }
	
	/*  eval shape functions on quadrature node: */
  
	out->Type->getValues(numQuadNodes,out->QuadNodes,out->S,out->dSdv);

	if (! Dudley_noError()) {
         Dudley_ShapeFunction_dealloc(out);
         return NULL;
	} 
  
	/*  all done: */
	out->reference_counter=1;
	return out;
}

Dudley_ShapeFunction* Dudley_ShapeFunction_reference(Dudley_ShapeFunction* in) {
     if (in!=NULL) ++(in->reference_counter);
     return in;
}
/**************************************************************/

void Dudley_ShapeFunction_dealloc(Dudley_ShapeFunction* in) {
  if (in!=NULL) { 
	  in->reference_counter--;
	  if (in->reference_counter<1) {
		  MEMFREE(in->QuadNodes);
		  MEMFREE(in->QuadWeights);
		  MEMFREE(in->S);
		  MEMFREE(in->dSdv);
		  MEMFREE(in);
	  }
  }
}

/**************************************************************/

Dudley_ShapeFunctionTypeId Dudley_ShapeFunction_getTypeId(char* element_type) 
{
    int ptr=0;
    Dudley_ShapeFunctionTypeId out=NoShape;
    while (Dudley_ShapeFunction_InfoList[ptr].TypeId!=NoShape && out==NoShape) {
       if (strcmp(element_type,Dudley_ShapeFunction_InfoList[ptr].Name)==0) out=Dudley_ShapeFunction_InfoList[ptr].TypeId;
       ptr++;
    }
    return out;
}

Dudley_ShapeFunctionInfo* Dudley_ShapeFunction_getInfo(Dudley_ShapeFunctionTypeId id)
{
    int ptr=0;
    Dudley_ShapeFunctionInfo* out=NULL;
    while (Dudley_ShapeFunction_InfoList[ptr].TypeId!=NoShape && out==NULL) {
       if (Dudley_ShapeFunction_InfoList[ptr].TypeId==id) out=&(Dudley_ShapeFunction_InfoList[ptr]);
       ptr++;
    }
    if (out==NULL) {
        Dudley_setError(VALUE_ERROR,"Dudley_ShapeFunctionInfo_getInfo: canot find requested shape function");
    }
    return out;
}

/**************************************************************/

#define V(_K_,_I_) v[INDEX2((_K_)-1,(_I_),DIM)]
#define S(_J_,_I_) s[S_INDEX((_J_)-1,(_I_),NUMSHAPES)]
#define DSDV(_J_,_K_,_I_) dsdv[DSDV_INDEX((_J_)-1,(_K_)-1,(_I_),NUMSHAPES,DIM)]

/**************************************************************/

void Dudley_Shape_Point1(int NumV,double* v,double* s,double* dsdv) {
  #define NUMSHAPES 1
  #define DIM 0
  int i;
  #pragma ivdep
  for (i=0;i<NumV;i++) {
    S(1,i)=1.;
  }
  #undef NUMSHAPES
  #undef DIM
}

/**************************************************************/

void Dudley_Shape_Line2(int NumV,double* v,double* s,double* dsdv) {
  #define NUMSHAPES 2
  #define DIM 1
  register double x;
  int i;
  #pragma ivdep
  for (i=0;i<NumV;i++) {
    x=V(1,i);
    S(1,i)=1.-x;
    S(2,i)=   x;
    DSDV(1,1,i)=-1.;
    DSDV(2,1,i)= 1.;
  }
  #undef NUMSHAPES
  #undef DIM
}

/**************************************************************/

void Dudley_Shape_Line3(int NumV,double* v,double* s,double* dsdv) {
  #define NUMSHAPES 3
  #define DIM 1
  register double x;
  int i;
  #pragma ivdep
  for (i=0;i<NumV;i++) {
    x=V(1,i);
    S(1,i)=(2.*x -1. )*(x -1.);
    S(2,i)=(2.*x -1.)*x;
    S(3,i)=  4.*x*(1. -x );
    DSDV(1,1,i)= 4.*x -3.;
    DSDV(2,1,i)= 4.*x -1.;
    DSDV(3,1,i)=-8.*x+4.;
 }
  #undef NUMSHAPES
  #undef DIM
}

/**************************************************************/

void Dudley_Shape_Line4(int NumV,double* v,double* s,double* dsdv) {
  #define NUMSHAPES 4
  #define DIM 1
  register double x;
  int i;
  #pragma ivdep
  for (i=0;i<NumV;i++) {
    x=V(1,i);
    S(1,i)=(10.)+(-5.5)*x+(9.)*x*x+(-4.5)*x*x*x ;
    S(2,i)=(10.)*x+(-4.5)*x*x+(4.5)*x*x*x  ;
    S(3,i)=(9.)*x+(-22.5)*x*x+(13.5)*x*x*x ;
    S(4,i)=(-4.5)*x+(18.)*x*x+(-13.5)*x*x*x;
    DSDV(1,1,i)=(-5.5)+(18.)*x+(-13.5)*x*x;
    DSDV(2,1,i)=(10.)+(-9.)*x+(13.5)*x*x;
    DSDV(3,1,i)=(9.)+(-45.)*x+(0.405e2)*x*x;
    DSDV(4,1,i)=(-4.5)+(36.)*x+(-0.405e2)*x*x;
  }
  #undef NUMSHAPES
  #undef DIM
}

/**************************************************************/

void Dudley_Shape_Tri3(int NumV,double* v,double* s,double* dsdv) {
  #define NUMSHAPES 3
  #define DIM 2
  register double x,y;
  int i;
  #pragma ivdep
  for (i=0;i<NumV;i++) {
    x=V(1,i);
    y=V(2,i);
    S(1,i)=1.-x-y;
    S(2,i)=   x;
    S(3,i)=   y;
    DSDV(1,1,i)=-1.;
    DSDV(1,2,i)=-1.;
    DSDV(2,1,i)= 1.;
    DSDV(2,2,i)= 0.;
    DSDV(3,1,i)= 0.;
    DSDV(3,2,i)= 1.;
  }
  #undef NUMSHAPES
  #undef DIM
}

/**************************************************************/

void Dudley_Shape_Tri6(int NumV,double* v,double* s,double* dsdv) {
  #define NUMSHAPES 6
  #define DIM 2
  register double x,y;
  int i;
  #pragma ivdep
  for (i=0;i<NumV;i++) {
    x=V(1,i);
    y=V(2,i);
    S(1,i)=  (1. -x -y)*(1. -2.*x -2.* y);
    S(2,i)=  x*(2.* x -1.);
    S(3,i)=  y*(2.* y -1.);
    S(4,i)=  (1. -x -y)*4.* x;
    S(5,i)=  4.*x*y;
    S(6,i)=  (1. -x -y)*4.* y;
    DSDV(1,1,i)= -3.+4.*x+4.*y;
    DSDV(1,2,i)= -3.+4.*x+4.*y;
    DSDV(2,1,i)= -1.+4.*x;
    DSDV(2,2,i)=     0.;
    DSDV(3,1,i)=     0.;
    DSDV(3,2,i)= -1.           +4.*y;
    DSDV(4,1,i)=     4. -8.*x -4.*y;
    DSDV(4,2,i)=        -4.*x;
    DSDV(5,1,i)=                     4.*y;
    DSDV(5,2,i)=          4.*x;
    DSDV(6,1,i)=                   -4.*y;
    DSDV(6,2,i)=     4. -4.*x -8.*y;
  }
  #undef NUMSHAPES
  #undef DIM
}

/**************************************************************/

void Dudley_Shape_Tri9(int NumV,double* v,double* s,double* dsdv) {
  #define NUMSHAPES 9
  #define DIM 2
  register double x,y;
  int i;
  #pragma ivdep
  for (i=0;i<NumV;i++) {
    x=V(1,i);
    y=V(2,i);
    S(1,i)=(10.)+(-5.5)*x+(-5.5)*y+(9.)*x*x+(-4.5)*x*x*x+(9.)*y*y+(-4.5)*y*y*y+(4.5)*x*y*y+(4.5)*x*x*y;
    S(2,i)=(10.)*x+(-4.5)*x*x+(4.5)*x*x*x;
    S(3,i)=(10.)*y+(-4.5)*y*y+(4.5)*y*y*y;
    S(4,i)=(9.)*x+(-22.5)*x*x+(13.5)*x*x*x+(-9.)*x*y*y+(4.5)*x*x*y;
    S(5,i)=(-4.5)*x+(18.)*x*x+(-13.5)*x*x*x+(4.5)*x*y*y+(-9.)*x*x*y;
    S(6,i)=(-4.5)*x*y*y+(9.)*x*x*y;
    S(7,i)=(9.)*x*y*y+(-4.5)*x*x*y;
    S(8,i)=(-4.5)*y+(18.)*y*y+(-13.5)*y*y*y+(-9.)*x*y*y+(4.5)*x*x*y;
    S(9,i)=(9.)*y+(-22.5)*y*y+(13.5)*y*y*y+(4.5)*x*y*y+(-9.)*x*x*y;
    DSDV(1, 1,i)=(-5.5)+(18.)*x+(-13.5)*x*x+(4.5)*y*y+(9.)*x*y;
    DSDV(2, 1,i)=(10.)+(-9.)*x+(13.5)*x*x;
    DSDV(3, 1,i)=           0.;
    DSDV(4, 1,i)=(9.)+(-45.)*x+(0.405e2)*x*x+(-9.)*y*y+(9.)*x*y;
    DSDV(5, 1,i)=(-4.5)+(36.)*x+(-0.405e2)*x*x+(4.5)*y*y+(-18.)*x*y;
    DSDV(6, 1,i)=(-4.5)*y*y+(18.)*x*y;
    DSDV(7, 1,i)=(9.)*y*y+(-9.)*x*y;
    DSDV(8, 1,i)=(-9.)*y*y+(9.)*x*y;
    DSDV(9, 1,i)=(4.5)*y*y+(-18.)*x*y;
    DSDV(1, 2,i)=(-5.5)+(18.)*y+(-13.5)*y*y+(9.)*x*y+(4.5)*x*x;
    DSDV(2, 2,i)=           0.;
    DSDV(3, 2,i)=(10.)+(-9.)*y+(13.5)*y*y;
    DSDV(4, 2,i)=(-18.)*x*y+(4.5)*x*x;
    DSDV(5, 2,i)=(9.)*x*y+(-9.)*x*x;
    DSDV(6, 2,i)=(-9.)*x*y+(9.)*x*x;
    DSDV(7, 2,i)=(18.)*x*y+(-4.5)*x*x;
    DSDV(8, 2,i)=(-4.5)+(36.)*y+(-0.405e2)*y*y+(-18.)*x*y+(4.5)*x*x;
    DSDV(9, 2,i)=(9.)+(-45.)*y+(0.405e2)*y*y+(9.)*x*y+(-9.)*x*x;
  }
  #undef NUMSHAPES
  #undef DIM
}

/**************************************************************/

void Dudley_Shape_Tri10(int NumV,double* v,double* s,double* dsdv) {
  #define NUMSHAPES 10
  #define DIM 2
  register double x,y;
  int i;
  #pragma ivdep
  for (i=0;i<NumV;i++) {
    x=V(1,i);
    y=V(2,i);
    S(1,i)=(10.)+(-5.5)*x+(-5.5)*y+(9.)*x*x+(-4.5)*x*x*x+(9.)*y*y+(-4.5)*y*y*y+(-13.5)*x*y*y+(-13.5)*x*x*y+(18.)*x*y;
    S(2,i)=(10.)*x+(-4.5)*x*x+(4.5)*x*x*x;
    S(3,i)=(10.)*y+(-4.5)*y*y+(4.5)*y*y*y;
    S(4,i)=(9.)*x+(-22.5)*x*x+(13.5)*x*x*x+(13.5)*x*y*y+(0.27e2)*x*x*y+(-22.5)*x*y;
    S(5,i)=(-4.5)*x+(18.)*x*x+(-13.5)*x*x*x+(-13.5)*x*x*y+(4.5)*x*y;
    S(6,i)=(13.5)*x*x*y+(-4.5)*x*y;
    S(7,i)=(13.5)*x*y*y+(-4.5)*x*y;
    S(8,i)=(-4.5)*y+(18.)*y*y+(-13.5)*y*y*y+(-13.5)*x*y*y+(4.5)*x*y;
    S(9,i)=(9.)*y+(-22.5)*y*y+(13.5)*y*y*y+(0.27e2)*x*y*y+(13.5)*x*x*y+(-22.5)*x*y;
    S(10,i)=(-0.27e2)*x*y*y+(-0.27e2)*x*x*y+(0.27e2)*x*y;
    DSDV(1, 1,i)=(-5.5)+(18.)*x+(-13.5)*x*x+(-13.5)*y*y+(-0.27e2)*x*y+(18.)*y;
    DSDV(2, 1,i)=(10.)+(-9.)*x+(13.5)*x*x;
    DSDV(3, 1,i)=            0.;
    DSDV(4, 1,i)=(9.)+(-45.)*x+(0.405e2)*x*x+(13.5)*y*y+(0.54e2)*x*y+(-22.5)*y;
    DSDV(5, 1,i)=(-4.5)+(36.)*x+(-0.405e2)*x*x+(-0.27e2)*x*y+(4.5)*y;
    DSDV(6, 1,i)=(0.27e2)*x*y+(-4.5)*y;
    DSDV(7, 1,i)=(13.5)*y*y+(-4.5)*y;
    DSDV(8, 1,i)=(-13.5)*y*y+(4.5)*y;
    DSDV(9, 1,i)=(0.27e2)*y*y+(0.27e2)*x*y+(-22.5)*y;
    DSDV(10, 1,i)=(-0.27e2)*y*y+(-0.54e2)*x*y+(0.27e2)*y;
    DSDV(1, 2,i)=(-5.5)+(18.)*y+(-13.5)*y*y+(-0.27e2)*x*y+(-13.5)*x*x+(18.)*x;
    DSDV(2, 2,i)=0.;
    DSDV(3, 2,i)=(10.)+(-9.)*y+(13.5)*y*y;
    DSDV(4, 2,i)=(0.27e2)*x*y+(0.27e2)*x*x+(-22.5)*x;
    DSDV(5, 2,i)=(-13.5)*x*x+(4.5)*x;
    DSDV(6, 2,i)=(13.5)*x*x+(-4.5)*x;
    DSDV(7, 2,i)=(0.27e2)*x*y+(-4.5)*x;
    DSDV(8, 2,i)=(-4.5)+(36.)*y+(-0.405e2)*y*y+(-0.27e2)*x*y+(4.5)*x;
    DSDV(9, 2,i)=(9.)+(-45.)*y+(0.405e2)*y*y+(0.54e2)*x*y+(13.5)*x*x+(-22.5)*x;
    DSDV(10, 2,i)=(-0.54e2)*x*y+(-0.27e2)*x*x+(0.27e2)*x;
  }
  #undef NUMSHAPES
  #undef DIM
}


/**************************************************************/

void Dudley_Shape_Tet4(int NumV,double* v,double* s,double* dsdv) {
  #define NUMSHAPES 4
  #define DIM 3
  register double x,y,z;
  int i;
  #pragma ivdep
  for (i=0;i<NumV;i++) {
    x=V(1,i);
    y=V(2,i);
    z=V(3,i);
    S(1,i)=1.-x-y-z;
    S(2,i)=x;
    S(3,i)=y;
    S(4,i)=z;
    DSDV(1,1,i)=-1.;
    DSDV(1,2,i)=-1.;
    DSDV(1,3,i)=-1.;
    DSDV(2,1,i)= 1.;
    DSDV(2,2,i)= 0.;
    DSDV(2,3,i)= 0.;
    DSDV(3,1,i)= 0.;
    DSDV(3,2,i)= 1.;
    DSDV(3,3,i)= 0.;
    DSDV(4,1,i)= 0.;
    DSDV(4,2,i)= 0.;
    DSDV(4,3,i)= 1.;
  }
  #undef NUMSHAPES
  #undef DIM
}

/**************************************************************/

void Dudley_Shape_Tet10(int NumV,double* v,double* s,double* dsdv) {
  #define NUMSHAPES 10
  #define DIM 3
  register double x,y,z;
  int i;
  #pragma ivdep
  for (i=0;i<NumV;i++) {
    x=V(1,i);
    y=V(2,i);
    z=V(3,i);
    S(1,i) = (1.-x-y-z)*(1.-2.*x-2.*y-2.*z);
    S(2,i) = x*(2.*x-1.);
    S(3,i) = y*(2.*y-1.);
    S(4,i) = z*(2.*z-1.);
    S(5,i) = (1.-x-y-z)*4.*x;
    S(6,i) = 4.*x*y;
    S(7,i) = (1.-x-y-z)*4.*y;
    S(8,i) = (1.-x-y-z)*4.*z;
    S(9,i) = 4.*x*z;
    S(10,i)= 4.*y*z;

    DSDV(1,1,i)= -3.+4.*x+4.*y+4.*z;
    DSDV(1,2,i)= -3.+4.*x+4.*y+4.*z;
    DSDV(1,3,i)= -3.+4.*x+4.*y+4.*z;


    DSDV(2,1,i)= -1.+4.*x;
    DSDV(2,2,i)=   0.;
    DSDV(2,3,i)=   0.;

    DSDV(3,1,i)=   0.;
    DSDV(3,2,i)=  -1.           +4.*y;
    DSDV(3,3,i)=   0.;

    DSDV(4,1,i)=   0.;
    DSDV(4,2,i)=   0.;
    DSDV(4,3,i)=  -1.                        +4.*z;

    DSDV(5,1,i)=   4. -8.*x -4.*y -4.*z;
    DSDV(5,2,i)=      -4.*x;
    DSDV(5,3,i)=      -4.*x;

    DSDV(6,1,i)=             4.*y;
    DSDV(6,2,i)=       4.*x;
    DSDV(6,3,i)=   0.;

    DSDV(7,1,i)=            -4.*y;
    DSDV(7,2,i)=   4. -4.*x -8.*y -4.*z;
    DSDV(7,3,i)=            -4.*y;

    DSDV(8,1,i)=                                -4.*z;
    DSDV(8,2,i)=                                -4.*z;
    DSDV(8,3,i)=   4. -4.*x -4.*y -8.*z;

    DSDV(9,1,i)=                                 4.*z;
    DSDV(9,2,i)=   0.;
    DSDV(9,3,i)=         4.*x;

    DSDV(10,1,i)=  0.;
    DSDV(10,2,i)=                                  4.*z;
    DSDV(10,3,i)=                      4.*y;
  }
  #undef NUMSHAPES
  #undef DIM
}

/**************************************************************/

void Dudley_Shape_Tet16(int NumV,double* v,double* s,double* dsdv) {
  #define NUMSHAPES 16
  #define DIM 3
  register double x,y,z;
  int i;
  #pragma ivdep
  for (i=0;i<NumV;i++) {
    x=V(1,i);
    y=V(2,i);
    z=V(3,i);
    S(1,i)=(10.)+(-5.5)*x+(-5.5)*y+(-5.5)*z+(9.)*x*x+(-4.5)*x*x*x+(4.5)*x*x*y+(4.5)*x*y*y+(-4.5)*y*y*y+(9.)*y*y+(9.)*z*z+(4.5)*x*x*z+(4.5)*y*y*z+(-4.5)*z*z*z+(4.5)*x*z*z+(4.5)*y*z*z;
    S(2,i)=(1.e0)*x+(-4.5)*x*x+(4.5)*x*x*x;
    S(3,i)=(1.e0)*y+(4.5)*y*y*y+(-4.5)*y*y;
    S(4,i)=(1.e0)*z+(-4.5)*z*z+(4.5)*z*z*z;
    S(5,i)=(9.)*x+(-22.5)*x*x+(13.5)*x*x*x+(4.5)*x*x*y+(-9.)*x*y*y+(4.5)*x*x*z+(-9.)*x*z*z;
    S(6,i)=(-4.5)*x+(18.)*x*x+(-13.5)*x*x*x+(-9.)*x*x*y+(4.5)*x*y*y+(-9.)*x*x*z+(4.5)*x*z*z;
    S(7,i)=(9.)*x*x*y+(-4.5)*x*y*y;
    S(8,i)=(-4.5)*x*x*y+(9.)*x*y*y;
    S(9,i)=(-4.5)*y+(4.5)*x*x*y+(-9.)*x*y*y+(-13.5)*y*y*y+(18.)*y*y+(-9.)*y*y*z+(4.5)*y*z*z;
    S(10,i)=(9.)*y+(-9.)*x*x*y+(4.5)*x*y*y+(13.5)*y*y*y+(-22.5)*y*y+(4.5)*y*y*z+(-9.)*y*z*z;
    S(11,i)=(9.)*z+(-22.5)*z*z+(-9.)*x*x*z+(-9.)*y*y*z+(13.5)*z*z*z+(4.5)*x*z*z+(4.5)*y*z*z;
    S(12,i)=(9.)*x*x*z+(-4.5)*x*z*z;
    S(13,i)=(9.)*y*y*z+(-4.5)*y*z*z;
    S(14,i)=(-4.5)*z+(18.)*z*z+(4.5)*x*x*z+(4.5)*y*y*z+(-13.5)*z*z*z+(-9.)*x*z*z+(-9.)*y*z*z;
    S(15,i)=(-4.5)*x*x*z+(9.)*x*z*z;
    S(16,i)=(-4.5)*y*y*z+(9.)*y*z*z;
    DSDV(1, 1,i)=(-5.5)+(18.)*x+(-13.5)*x*x+(9.)*x*y+(4.5)*y*y+(9.)*x*z+(4.5)*z*z;
    DSDV(2, 1,i)=(1.e0)+(-9.)*x+(13.5)*x*x;
    DSDV(3, 1,i)=            0.;
    DSDV(4, 1,i)=            0.;
    DSDV(5, 1,i)=(9.)+(-45.)*x+(0.405e2)*x*x+(9.)*x*y+(-9.)*y*y+(9.)*x*z+(-9.)*z*z;
    DSDV(6, 1,i)=(-4.5)+(36.)*x+(-0.405e2)*x*x+(-18.)*x*y+(4.5)*y*y+(-18.)*x*z+(4.5)*z*z;
    DSDV(7, 1,i)=(18.)*x*y+(-4.5)*y*y;
    DSDV(8, 1,i)=(-9.)*x*y+(9.)*y*y;
    DSDV(9, 1,i)=(9.)*x*y+(-9.)*y*y;
    DSDV(10, 1,i)=(-18.)*x*y+(4.5)*y*y;
    DSDV(11, 1,i)=(-18.)*x*z+(4.5)*z*z;
    DSDV(12, 1,i)=(18.)*x*z+(-4.5)*z*z;
    DSDV(13, 1,i)=0.;
    DSDV(14, 1,i)=(9.)*x*z+(-9.)*z*z;
    DSDV(15, 1,i)=(-9.)*x*z+(9.)*z*z;
    DSDV(16, 1,i)=0.;
    DSDV(1, 2,i)=(-5.5)+(4.5)*x*x+(9.)*x*y+(-13.5)*y*y+(18.)*y+(9.)*y*z+(4.5)*z*z;
    DSDV(2, 2,i)=0.;
    DSDV(3, 2,i)=(1.e0)+(13.5)*y*y+(-9.)*y;
    DSDV(4, 2,i)=0.;
    DSDV(5, 2,i)=(4.5)*x*x+(-18.)*x*y;
    DSDV(6, 2,i)=(-9.)*x*x+(9.)*x*y;
    DSDV(7, 2,i)=(9.)*x*x+(-9.)*x*y;
    DSDV(8, 2,i)=(-4.5)*x*x+(18.)*x*y;
    DSDV(9, 2,i)=(-4.5)+(4.5)*x*x+(-18.)*x*y+(-0.405e2)*y*y+(36.)*y+(-18.)*y*z+(4.5)*z*z;
    DSDV(10, 2,i)=(9.)+(-9.)*x*x+(9.)*x*y+(0.405e2)*y*y+(-45.)*y+(9.)*y*z+(-9.)*z*z;
    DSDV(11, 2,i)=(-18.)*y*z+(4.5)*z*z;
    DSDV(12, 2,i)=0.;
    DSDV(13, 2,i)=(18.)*y*z+(-4.5)*z*z;
    DSDV(14, 2,i)=(9.)*y*z+(-9.)*z*z;
    DSDV(15, 2,i)=0.;
    DSDV(16, 2,i)=(-9.)*y*z+(9.)*z*z;
    DSDV(1, 3,i)=(-5.5)+(18.)*z+(4.5)*x*x+(4.5)*y*y+(-13.5)*z*z+(.9e1)*x*z+(9.)*y*z;
    DSDV(2, 3,i)=           0.;
    DSDV(3, 3,i)=           0.;
    DSDV(4, 3,i)=(1.e0)+(-9.)*z+(13.5)*z*z;
    DSDV(5, 3,i)=(4.5)*x*x+(-18.)*x*z;
    DSDV(6, 3,i)=(-9.)*x*x+(9.)*x*z;
    DSDV(7, 3,i)=           0.;
    DSDV(8, 3,i)=           0.;
    DSDV(9, 3,i)=(-9.)*y*y+(9.)*y*z;
    DSDV(10, 3,i)=(4.5)*y*y+(-18.)*y*z;
    DSDV(11, 3,i)=(9.)+(-45.)*z+(-9.)*x*x+(-9.)*y*y+(0.405e2)*z*z+(.9e1)*x*z+(9.)*y*z;
    DSDV(12, 3,i)=(9.)*x*x+(-9.)*x*z;
    DSDV(13, 3,i)=(9.)*y*y+(-9.)*y*z;
    DSDV(14, 3,i)=(-4.5)+(36.)*z+(4.5)*x*x+(4.5)*y*y+(-0.405e2)*z*z+(-18.)*x*z+(-18.)*y*z;
    DSDV(15, 3,i)=(-4.5)*x*x+(18.)*x*z;
    DSDV(16, 3,i)=(-4.5)*y*y+(18.)*y*z;
  }
  #undef NUMSHAPES
  #undef DIM
}


#undef V
#undef S
#undef DSDV


