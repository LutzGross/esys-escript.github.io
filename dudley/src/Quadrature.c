
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

/*   Dudley: quadrature schemes */

/**************************************************************/

#include "Quadrature.h"


#define QUADNODES(_K_,_I_) quadNodes[INDEX2(_K_,_I_,DIM)]
#define QUADWEIGHTS(_I_) quadWeights[_I_]

/**************************************************************/

Dudley_QuadInfo Dudley_QuadInfoList[]={
	{PointQuad, "Point", 0,  1, 	Dudley_Quad_getNodesPoint,		Dudley_Quad_getNumNodesPoint} ,
        {LineQuad,  "Line",  1,  2,	Dudley_Quad_getNodesLine,		Dudley_Quad_getNumNodesLine} ,
	{TriQuad,   "Tri",   2,  3,   	Dudley_Quad_getNodesTri,   	       Dudley_Quad_getNumNodesTri},
	{TetQuad,   "Tet",   3,  4,   	Dudley_Quad_getNodesTet,        Dudley_Quad_getNumNodesTet},
	{NoQuad, "NoType", 0,  1,	Dudley_Quad_getNodesPoint,		Dudley_Quad_getNumNodesPoint}
};

Dudley_QuadInfo* Dudley_QuadInfo_getInfo(Dudley_QuadTypeId id) 
{
    int ptr=0;
    Dudley_QuadInfo* out=NULL;
    while (Dudley_QuadInfoList[ptr].TypeId!=NoQuad && out==NULL) {
       if (Dudley_QuadInfoList[ptr].TypeId==id) out=&(Dudley_QuadInfoList[ptr]);
       ptr++;
    }
    if (out==NULL) {
        Dudley_setError(VALUE_ERROR,"Dudley_QuadInfo_getInfo: canot find requested quadrature scheme.");
    }
    return out;
}

/**************************************************************/

/*   get a quadrature scheme with numQuadNodes quadrature nodes for the tri  */
/*   as a queezed scheme on a quad [0,1]^2 */

void Dudley_Quad_getNodesTri(int numQuadNodes,double* quadNodes,double* quadWeights) {
//  int i;
//  double Q1,Q2,a,b,c,d,e,f,g,u,v,w;
  #define DIM 2
  
  /*  the easy cases: */
  if (numQuadNodes==1) {
    QUADNODES(0,0)=1./3.;
    QUADNODES(1,0)=1./3.;
    QUADWEIGHTS(0)=1./2.;
  } else if (numQuadNodes==3){
    QUADNODES(0,0)=1./2.;
    QUADNODES(1,0)=0.;
    QUADWEIGHTS(0)=1./6.;
    QUADNODES(0,1)=0.;
    QUADNODES(1,1)=1./2.;
    QUADWEIGHTS(1)=1./6.;
    QUADNODES(0,2)=1./2.;
    QUADNODES(1,2)=1./2.;
    QUADWEIGHTS(2)=1./6.;
  } 
  #undef DIM


}

/**************************************************************/

/*   get a quadrature scheme with numQuadNodes quadrature nodes for the tet */
/*   as a queezed scheme on a hex [0,1]^3 */

void Dudley_Quad_getNodesTet(int numQuadNodes,double* quadNodes,double* quadWeights) {
  double alpha=0.58541019662496852;
  double beta =0.1381966011250105;
  #define DIM 3
  
  /*  the easy cases: */
  if (numQuadNodes==1) {
    QUADNODES(0,0)=0.25;
    QUADNODES(1,0)=0.25;
    QUADNODES(2,0)=0.25;
    QUADWEIGHTS(0)=1./6.;
  } else if (numQuadNodes==4){
    QUADNODES(0,0)=beta;
    QUADNODES(1,0)=beta;
    QUADNODES(2,0)=beta;
    QUADWEIGHTS(0)=1./24.;
    QUADNODES(0,1)=alpha;
    QUADNODES(1,1)=beta;
    QUADNODES(2,1)=beta;
    QUADWEIGHTS(1)=1./24.;
    QUADNODES(0,2)=beta;
    QUADNODES(1,2)=alpha;
    QUADNODES(2,2)=beta;
    QUADWEIGHTS(2)=1./24.;
    QUADNODES(0,3)=beta;
    QUADNODES(1,3)=beta;
    QUADNODES(2,3)=alpha;
    QUADWEIGHTS(3)=1./24.;
  }
  #undef DIM
}

/**************************************************************/

/*   get a quadrature scheme with     DSDV(2,2,i)= 0.;
    DSDV(2,3,i)= 0.;
    DSDV(3,1,i)= 0.;numQuadNodes quadrature nodes for a point. As there */
/*   in no quadrature scheme for a point any value for numQuadNodes other than 0 throws */
/*   an error. */

void Dudley_Quad_getNodesPoint(int numQuadNodes,double* quadNodes,double* quadWeights) {
  if (numQuadNodes==0) {
        return;
  } else {
       Dudley_setError(VALUE_ERROR,"Dudley_Quad_getNodesPoint: Illegal number of quadrature nodes.");
  }
}

/**************************************************************/

/*   get a quadrature scheme with numQuadNodes quadrature nodes on the line [0,1]: */
/*   The nodes and weights are set from a table. */

void Dudley_Quad_getNodesLine(int numQuadNodes,double* quadNodes,double* quadWeights) {
  switch(numQuadNodes) {
      case 1:
        quadNodes[0]=0.5;
        quadWeights[0]=1.;
        break;

      case 2:
        quadNodes[0]=(1.-.577350269189626)/2.;
        quadNodes[1]=(1.+.577350269189626)/2.;
        quadWeights[0]=.5;
        quadWeights[1]=.5;
        break;
      default:
        Dudley_setError(VALUE_ERROR,"Dudley_Quad_getNodesLine: Negative intergration order.");
        break;
  }
}


/**************************************************************/

/*    The following functions Dudley_Quad_getNumNodes* return the nmber of quadrature points needed to */
/*    achieve a certain accuracy. Notice that for Tet and Tri the the order is increased */
/*    to consider the accuracy reduction through the construction process.  */


int Dudley_Quad_getNumNodesPoint(int order) {
    return 0;
}

int Dudley_Quad_getNumNodesLine(int order) {
  char error_msg[LenErrorMsg_MAX];
  if (order <0 ) {
    Dudley_setError(VALUE_ERROR,"Dudley_Quad_getNumNodesPoint: Negative intergration order.");
    return -1;
  } else { 
    if (order > 2*MAX_numQuadNodesLine-1) {
      sprintf(error_msg,"Dudley_Quad_getNumNodesPoint: requested integration order %d on line is too large (>%d).",
                                                           order,2*MAX_numQuadNodesLine-1);
      Dudley_setError(VALUE_ERROR,error_msg);
      return -1;
    } else { 
       Dudley_resetError();
       return order/2+1;
    }
  }
}

int Dudley_Quad_getNumNodesTri(int order) {
  if (order<=1) {
      return 1;
  } else if (order<=2){
      return 3;
  } else {
      return -1;
  }
}

int Dudley_Quad_getNumNodesTet(int order) {
  if (order<=1) {
      return 1;
  } else if (order<=2){
      return 4;
  } else {
      return -1;
  }
}

