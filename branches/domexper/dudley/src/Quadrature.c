
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


/*
else if (numQuadNodes==4){
    QUADNODES(0,0)=1./3.;
    QUADNODES(1,0)=1./3.;
    QUADWEIGHTS(0)=-27./96.;
    QUADNODES(0,1)=0.2;
    QUADNODES(1,1)=0.2;
    QUADWEIGHTS(1)=25./96.;
    QUADNODES(0,2)=0.6;
    QUADNODES(1,2)=0.2;
    QUADWEIGHTS(2)=25./96.;
    QUADNODES(0,3)=0.2;
    QUADNODES(1,3)=0.6;
    QUADWEIGHTS(3)=25./96.;
  } else if (numQuadNodes==6){
      QUADWEIGHTS(0) = 0.109951743655322/2.;
      QUADWEIGHTS(1) = 0.109951743655322/2.;
      QUADWEIGHTS(2) = 0.109951743655322/2.;
      QUADWEIGHTS(3) = 0.223381589678011/2.;
      QUADWEIGHTS(4) = 0.223381589678011/2.;
      QUADWEIGHTS(5) = 0.223381589678011/2.;

      QUADNODES(0,0) = 0.816847572980459;
      QUADNODES(0,1) = 0.091576213509771;
      QUADNODES(0,2) = 0.091576213509771;
      QUADNODES(0,3) = 0.108103018168070;
      QUADNODES(0,4) = 0.445948490915965;
      QUADNODES(0,5) = 0.445948490915965;

      QUADNODES(1,0) = 0.091576213509771;
      QUADNODES(1,1) = 0.816847572980459;
      QUADNODES(1,2) = 0.091576213509771;
      QUADNODES(1,3) = 0.445948490915965;
      QUADNODES(1,4) = 0.108103018168070;
      QUADNODES(1,5) = 0.445948490915965;

  } else if (numQuadNodes==7){
      QUADNODES(0,0) = 0.33333333333333333;
      QUADNODES(0,1) = 0.7974269853530872;
      QUADNODES(0,2) = 0.10128650732345633;
      QUADNODES(0,3) = 0.10128650732345633;
      QUADNODES(0,4) = 0.059715871789769809;
      QUADNODES(0,5) = 0.47014206410511505;
      QUADNODES(0,6) = 0.47014206410511505;

      QUADNODES(1,0) = 0.33333333333333333;
      QUADNODES(1,1) = 0.10128650732345633;
      QUADNODES(1,2) = 0.7974269853530872;
      QUADNODES(1,3) = 0.10128650732345633;
      QUADNODES(1,4) = 0.47014206410511505;
      QUADNODES(1,5) = 0.059715871789769809;
      QUADNODES(1,6) = 0.47014206410511505;

      QUADWEIGHTS(0) = 0.225/2.;
      QUADWEIGHTS(1) = 0.12593918054482717/2.;
      QUADWEIGHTS(2) = 0.12593918054482717/2.;
      QUADWEIGHTS(3) = 0.12593918054482717/2.;
      QUADWEIGHTS(4) = 0.13239415278850616/2.;
      QUADWEIGHTS(5) = 0.13239415278850616/2.;
      QUADWEIGHTS(6) = 0.13239415278850616/2.;

   } else if (numQuadNodes==12){
       a = 0.873821971016996;
       b = 0.063089014491502;
       c = 0.501426509658179;
       d = 0.249286745170910;
       e = 0.636502499121399;
       f = 0.310352451033785;
       g = 0.053145049844816;

       u = 0.050844906370207/2.;
       v = 0.116786275726379/2.;
       w = 0.082851075618374/2.;

          QUADNODES(0,0) = a;
          QUADNODES(0,1) =  b;
          QUADNODES(0,2) =  b;
          QUADNODES(0,3) =  c;
          QUADNODES(0,4) =  d;
          QUADNODES(0,5) =  d;
          QUADNODES(0,6) =  e;
          QUADNODES(0,7) =  e;
          QUADNODES(0,8) =  f;
          QUADNODES(0,9) =  f;
          QUADNODES(0,10) =  g;
          QUADNODES(0,11) =  g;

          QUADNODES(1,0) = b;
          QUADNODES(1,1) =  a;
          QUADNODES(1,2) =  b;
          QUADNODES(1,3) =  d;
          QUADNODES(1,4) =  c;
          QUADNODES(1,5) =  d;
          QUADNODES(1,6) =  f;
          QUADNODES(1,7) =  g;
          QUADNODES(1,8) =  e;
          QUADNODES(1,9) =  g;
          QUADNODES(1,10) =  e;
          QUADNODES(1,11) =  f;

          QUADWEIGHTS(0)= u;
          QUADWEIGHTS(1)= u;
          QUADWEIGHTS(2)= u;
          QUADWEIGHTS(3)= v;
          QUADWEIGHTS(4)= v;
          QUADWEIGHTS(5)= v;
          QUADWEIGHTS(6)= w;
          QUADWEIGHTS(7)= w;
          QUADWEIGHTS(8)= w;
          QUADWEIGHTS(9)= w;
          QUADWEIGHTS(10)= w;
          QUADWEIGHTS(11)= w;

  } else if (numQuadNodes==13){
      QUADWEIGHTS(0) =-0.149570044467670/2.;
      QUADWEIGHTS(1) = 0.175615257433204/2.;
      QUADWEIGHTS(2) = 0.175615257433204/2.;
      QUADWEIGHTS(3) = 0.175615257433204/2.;
      QUADWEIGHTS(4) = 0.053347235608839/2.;
      QUADWEIGHTS(5) = 0.053347235608839/2.;
      QUADWEIGHTS(6) = 0.053347235608839/2.;
      QUADWEIGHTS(7) = 0.077113760890257/2.;
      QUADWEIGHTS(8) = 0.077113760890257/2.;
      QUADWEIGHTS(9) = 0.077113760890257/2.;
      QUADWEIGHTS(10) = 0.077113760890257/2.;
      QUADWEIGHTS(11) = 0.077113760890257/2.;
      QUADWEIGHTS(12) = 0.077113760890257/2.;

      QUADNODES(0,0) = 0.3333333333333333;
      QUADNODES(0,1) = 0.479308067841923;
      QUADNODES(0,2) = 0.260345966079038;
      QUADNODES(0,3) = 0.260345966079038;
      QUADNODES(0,4) = 0.869739794195568;
      QUADNODES(0,5) = 0.065130102902216;
      QUADNODES(0,6) = 0.065130102902216;
      QUADNODES(0,7) = 0.638444188569809;
      QUADNODES(0,8) = 0.638444188569809;
      QUADNODES(0,9) = 0.048690315425316;
      QUADNODES(0,10) = 0.048690315425316;
      QUADNODES(0,11) = 0.312865496004875;
      QUADNODES(0,12) = 0.312865496004875;

      QUADNODES(1,0) = 0.3333333333333333;
      QUADNODES(1,1) = 0.260345966079038;
      QUADNODES(1,2) = 0.479308067841923;
      QUADNODES(1,3) = 0.260345966079038;
      QUADNODES(1,4) = 0.065130102902216;
      QUADNODES(1,5) = 0.869739794195568;
      QUADNODES(1,6) = 0.065130102902216;
      QUADNODES(1,7) = 0.048690315425316;
      QUADNODES(1,8) = 0.312865496004875;
      QUADNODES(1,9) = 0.638444188569809;
      QUADNODES(1,10) = 0.312865496004875;
      QUADNODES(1,11) = 0.638444188569809;
      QUADNODES(1,12) = 0.048690315425316;

  } else if (numQuadNodes==16){
      QUADWEIGHTS(0) = 0.07215780;
      QUADWEIGHTS(1) = 0.04754582;
      QUADWEIGHTS(2) = 0.04754582;
      QUADWEIGHTS(3) = 0.04754582;
      QUADWEIGHTS(4) = 0.01622925;
      QUADWEIGHTS(5) = 0.01622925;
      QUADWEIGHTS(6) = 0.01622925;
      QUADWEIGHTS(7) = 0.05160869;
      QUADWEIGHTS(8) = 0.05160869;
      QUADWEIGHTS(9) = 0.05160869;
      QUADWEIGHTS(10) = 0.01361516;
      QUADWEIGHTS(11) = 0.01361516;
      QUADWEIGHTS(12) = 0.01361516;
      QUADWEIGHTS(13) = 0.01361516;
      QUADWEIGHTS(14) = 0.01361516;
      QUADWEIGHTS(15) = 0.01361516;
 
      QUADNODES(0,0) = 0.3333333;
      QUADNODES(0,1) = 0.08141482;
      QUADNODES(0,2) = 0.4592926;
      QUADNODES(0,3) = 0.4592926;
      QUADNODES(0,4) = 0.8989055;
      QUADNODES(0,5) = 0.05054723;
      QUADNODES(0,6) = 0.05054723;
      QUADNODES(0,7) = 0.6588614;
      QUADNODES(0,8) = 0.1705693;
      QUADNODES(0,9) = 0.1705693;
      QUADNODES(0,10) = 0.008394777;
      QUADNODES(0,11) = 0.008394777;
      QUADNODES(0,12) = 0.7284924;
      QUADNODES(0,13) = 0.7284924;
      QUADNODES(0,14) = 0.2631128;
      QUADNODES(0,15) = 0.2631128;
 
      QUADNODES(1,0) = 0.3333333;
      QUADNODES(1,1) = 0.4592926;
      QUADNODES(1,2) = 0.08141482;
      QUADNODES(1,3) = 0.4592926;
      QUADNODES(1,4) = 0.05054723;
      QUADNODES(1,5) = 0.8989055;
      QUADNODES(1,6) = 0.05054723;
      QUADNODES(1,7) = 0.1705693;
      QUADNODES(1,8) = 0.6588614;
      QUADNODES(1,9) = 0.1705693;
      QUADNODES(1,10) = 0.7284924;
      QUADNODES(1,11) = 0.2631128;
      QUADNODES(1,12) = 0.008394777;
      QUADNODES(1,13) = 0.2631128;
      QUADNODES(1,14) = 0.008394777;
      QUADNODES(1,15) = 0.7284924;

  } else if (numQuadNodes==19){
      QUADWEIGHTS(0) = 0.04856790;
      QUADWEIGHTS(1) = 0.01566735;
      QUADWEIGHTS(2) = 0.01566735;
      QUADWEIGHTS(3) = 0.01566735;
      QUADWEIGHTS(4) = 0.03891377;
      QUADWEIGHTS(5) = 0.03891377;
      QUADWEIGHTS(6) = 0.03891377;
      QUADWEIGHTS(7) = 0.03982387;
      QUADWEIGHTS(8) = 0.03982387;
      QUADWEIGHTS(9) = 0.03982387;
      QUADWEIGHTS(10) = 0.01278884;
      QUADWEIGHTS(11) = 0.01278884;
      QUADWEIGHTS(12) = 0.01278884;
      QUADWEIGHTS(13) = 0.02164177;
      QUADWEIGHTS(14) = 0.02164177;
      QUADWEIGHTS(15) = 0.02164177;
      QUADWEIGHTS(16) = 0.02164177;
      QUADWEIGHTS(17) = 0.02164177;
      QUADWEIGHTS(18) = 0.02164177;
 
      QUADNODES(0,0) = 0.3333333;
      QUADNODES(0,1) = 0.02063496;
      QUADNODES(0,2) = 0.4896825;
      QUADNODES(0,3) = 0.4896825;
      QUADNODES(0,4) = 0.1258208;
      QUADNODES(0,5) = 0.4370896;
      QUADNODES(0,6) = 0.4370896;
      QUADNODES(0,7) = 0.6235929;
      QUADNODES(0,8) = 0.1882035;
      QUADNODES(0,9) = 0.1882035;
      QUADNODES(0,10) = 0.9105410;
      QUADNODES(0,11) = 0.04472951;
      QUADNODES(0,12) = 0.04472951;
      QUADNODES(0,13) = 0.03683841;
      QUADNODES(0,14) = 0.03683841;
      QUADNODES(0,15) = 0.7411986;
      QUADNODES(0,16) = 0.7411986;
      QUADNODES(0,17) = 0.2219630;
      QUADNODES(0,18) = 0.2219630;
 
      QUADNODES(1,0) = 0.3333333;
      QUADNODES(1,1) = 0.4896825;
      QUADNODES(1,2) = 0.02063496;
      QUADNODES(1,3) = 0.4896825;
      QUADNODES(1,4) = 0.4370896;
      QUADNODES(1,5) = 0.1258208;
      QUADNODES(1,6) = 0.4370896;
      QUADNODES(1,7) = 0.1882035;
      QUADNODES(1,8) = 0.6235929;
      QUADNODES(1,9) = 0.1882035;
      QUADNODES(1,10) = 0.04472951;
      QUADNODES(1,11) = 0.9105410;
      QUADNODES(1,12) = 0.04472951;
      QUADNODES(1,13) = 0.7411986;
      QUADNODES(1,14) = 0.2219630;
      QUADNODES(1,15) = 0.03683841;
      QUADNODES(1,16) = 0.2219630;
      QUADNODES(1,17) = 0.03683841;
      QUADNODES(1,18) = 0.7411986;
  } else {
*/
//    /*  get scheme on [0.1]^2 */
/*    Dudley_Quad_getNodesRec(numQuadNodes,quadNodes,quadWeights);
    if (! Dudley_noError()) return;
    
    
    for (i=0;i<numQuadNodes;i++) {
        Q1=QUADNODES(0,i);
        Q2=QUADNODES(1,i);
        QUADWEIGHTS(i)=QUADWEIGHTS(i)*(1.-(1./2.)*(Q1+Q2));
        QUADNODES(0,i)=Q1*(1.-(1./2.)*Q2);
        QUADNODES(1,i)=Q2*(1.-(1./2.)*Q1);
    }
  }*/
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

/*   get a quadrature scheme with numQuadNodes quadrature nodes for a point. As there */
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

