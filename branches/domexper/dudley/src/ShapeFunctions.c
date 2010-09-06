
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
	{Tri3Shape,   "Tri3",   2,  3, 1, 3,	Dudley_Shape_Tri3   },
	{Tet4Shape,   "Tet4",   3,  4, 1, 4,	Dudley_Shape_Tet4,  },
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

#undef V
#undef S
#undef DSDV


