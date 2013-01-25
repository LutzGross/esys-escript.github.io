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

/*   Finley:  */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Quadrature.h"

#define QUADNODES(_K_,_I_) quadNodes[INDEX2(_K_,_I_,DIM)]
#define QUADWEIGHTS(_I_) quadWeights[_I_]

/**************************************************************/

/*   get a quadrature scheme with numQuadNodes quadrature nodes for the tri  */
/*   as a queezed scheme on a quad [0,1]^2 */

void Finley_Quad_getNodesTri(int numQuadNodes,double* quadNodes,double* quadWeights) {
  int i;
  double Q1,Q2;
  #define DIM 2
  
  /*  the easy case: */
  
  if (numQuadNodes==1) {
    QUADNODES(0,0)=1./3.;
    QUADNODES(1,0)=1./3.;
    QUADWEIGHTS(0)= .5;
  } else {
    
    /*  get scheme on [0.1]^2 */
    
    Finley_Quad_getNodesRec(numQuadNodes,quadNodes,quadWeights);
    if (! Finley_noError()) return;
    
    /*  squeeze it: */
    
    for (i=0;i<numQuadNodes;i++) {
        Q1=QUADNODES(0,i);
        Q2=QUADNODES(1,i);
        QUADWEIGHTS(i)=QUADWEIGHTS(i)*(1.-(1./2.)*(Q1+Q2));
        QUADNODES(0,i)=Q1*(1.-(1./2.)*Q2);
        QUADNODES(1,i)=Q2*(1.-(1./2.)*Q1);
    }
  }
  #undef DIM
}

/**************************************************************/

/*   get a quadrature scheme with numQuadNodes quadrature nodes for the tet */
/*   as a queezed scheme on a hex [0,1]^3 */

void Finley_Quad_getNodesTet(int numQuadNodes,double* quadNodes,double* quadWeights) {
  int i;
  double Q1,Q2,Q3,JA11,JA12,JA13,JA21,JA22,JA23,JA31,JA32,JA33,DET;
  #define DIM 3
  
  /*  the easy case: */
  
  if (numQuadNodes==1) {
    QUADNODES(0,0)= .25;
    QUADNODES(1,0)= .25;
    QUADNODES(2,0)= .25;
    QUADWEIGHTS(0)=1./6.;
  } else {
    
    /*  get scheme on [0.1]^3 */
    
    Finley_Quad_getNodesHex(numQuadNodes,quadNodes,quadWeights);
    if (! Finley_noError()) return;
    
    /*  squeeze it: */
    
    for (i=0;i<numQuadNodes;i++) {
      Q1=QUADNODES(0,i);
      Q2=QUADNODES(1,i);
      Q3=QUADNODES(2,i);

      JA11= (1./3.)*Q2*Q3-(1./2.)*(Q2+Q3) +1.;
      JA12= (1./3.)*Q1*Q3-(1./2.)*Q1;
      JA13= (1./3.)*Q1*Q2-(1./2.)*Q1;
      JA21= (1./3.)*Q2*Q3-(1./2.)*Q2;
      JA22= (1./3.)*Q1*Q3-(1./2.)*(Q1+Q3) +1.;
      JA23= (1./3.)*Q1*Q2-(1./2.)*Q2;
      JA31= (1./3.)*Q2*Q3-(1./2.)*Q3;
      JA32= (1./3.)*Q1*Q3-(1./2.)*Q3;
      JA33= (1./3.)*Q1*Q2-(1./2.)*(Q1+Q2) +1.;
      DET=JA11*JA22*JA33+JA12*JA23*JA31+JA13*JA21*JA32-JA13*JA22*JA31-JA11*JA23*JA32-JA12*JA21*JA33;
      quadWeights[i]=quadWeights[i]*ABS(DET);
      QUADNODES(0,i)=Q1*((1./3.)*Q2*Q3-(1./2.)*(Q2+Q3)+1.);
      QUADNODES(1,i)=Q2*((1./3.)*Q1*Q3-(1./2.)*(Q1+Q3)+1.);
      QUADNODES(2,i)=Q3*((1./3.)*Q1*Q2-(1./2.)*(Q1+Q2)+1.);
    }
  }
  #undef DIM
}

/**************************************************************/

/*   get a quadrature scheme with numQuadNodes quadrature nodes for the quad [0.1]^2 */
/*   as a X-product of a 1D scheme. */

void Finley_Quad_getNodesRec(int numQuadNodes,double* quadNodes,double* quadWeights) {
  char error_msg[LenErrorMsg_MAX];
  int numQuadNodes1d,i,j,l;
  double *quadNodes1d=NULL,*quadWeights1d=NULL;
  bool_t set=FALSE;
  #define DIM 2
  
  quadNodes1d=TMPMEMALLOC(numQuadNodes,double);
  quadWeights1d=TMPMEMALLOC(numQuadNodes,double);
  if (! ( Finley_checkPtr(quadNodes1d) || Finley_checkPtr(quadWeights1d) ) ) {
     /*  find numQuadNodes1d with numQuadNodes1d**2==numQuadNodes: */
     
     for (numQuadNodes1d=1;numQuadNodes1d<=MAX_numQuadNodesLine;numQuadNodes1d++) {
       if (numQuadNodes1d*numQuadNodes1d==numQuadNodes) {
      
         /*  get 1D scheme: */
         
         Finley_Quad_getNodesLine(numQuadNodes1d,quadNodes1d,quadWeights1d);
      
         /*  make 2D scheme: */
      
         if (Finley_noError()) {
           l=0;
           for (i=0;i<numQuadNodes1d;i++) {
             for (j=0;j<numQuadNodes1d;j++) {
               QUADNODES(0,l)=quadNodes1d[i];
               QUADNODES(1,l)=quadNodes1d[j];
               QUADWEIGHTS(l)=quadWeights1d[i]*quadWeights1d[j];
               l++;
             }
           }
           set=TRUE;
           break;
         }
       }
     }
     if (!set) {
         sprintf(error_msg,"Finley_Quad_getNodesRec: Illegal number of quadrature nodes %d on hexahedron.",numQuadNodes);
         Finley_setError(VALUE_ERROR,error_msg);
     }
     TMPMEMFREE(quadNodes1d);
     TMPMEMFREE(quadWeights1d);
   }
   #undef DIM
}

/**************************************************************/

/*   get a quadrature scheme with numQuadNodes quadrature nodes for the hex [0.1]^3 */
/*   as a X-product of a 1D scheme. */

void Finley_Quad_getNodesHex(int numQuadNodes,double* quadNodes,double* quadWeights) {
  char error_msg[LenErrorMsg_MAX];
  int numQuadNodes1d,i,j,k,l;
  double *quadNodes1d=NULL,*quadWeights1d=NULL;
  bool_t set;
  #define DIM 3
  
  /*  find numQuadNodes1d with numQuadNodes1d**3==numQuadNodes: */
  
  quadNodes1d=TMPMEMALLOC(numQuadNodes,double);
  quadWeights1d=TMPMEMALLOC(numQuadNodes,double);
  if (! ( Finley_checkPtr(quadNodes1d) || Finley_checkPtr(quadWeights1d) ) ) {
     for (numQuadNodes1d=1;numQuadNodes1d<=MAX_numQuadNodesLine;numQuadNodes1d++) {
       if (numQuadNodes1d*numQuadNodes1d*numQuadNodes1d==numQuadNodes) {
      
         /*  get 1D scheme: */
      
         Finley_Quad_getNodesLine(numQuadNodes1d,quadNodes1d,quadWeights1d);
      
         /*  make 3D scheme: */
      
         if (Finley_noError()) {
           l=0;
           for (i=0;i<numQuadNodes1d;i++) {
             for (j=0;j<numQuadNodes1d;j++) {
               for (k=0;k<numQuadNodes1d;k++) {
                 QUADNODES(0,l)=quadNodes1d[i];
                 QUADNODES(1,l)=quadNodes1d[j];
                 QUADNODES(2,l)=quadNodes1d[k];
                 QUADWEIGHTS(l)=quadWeights1d[i]*quadWeights1d[j]*quadWeights1d[k];
                 l++;
               }
             }
           }
           set=TRUE;
           break;
         }
       }
     }
     if (!set) {
          sprintf(error_msg,"Finley_Quad_getNodesHex: Illegal number of quadrature nodes %d on hexahedron.",numQuadNodes);
          Finley_setError(VALUE_ERROR,error_msg);
     }
     TMPMEMFREE(quadNodes1d);
     TMPMEMFREE(quadWeights1d);
  }
  #undef DIM
}

/**************************************************************/

/*   get a quadrature scheme with numQuadNodes quadrature nodes for a point. As there */
/*   in no quadrature scheme for a point any value for numQuadNodes other than 0 throws */
/*   an error. */

void Finley_Quad_getNodesPoint(int numQuadNodes,double* quadNodes,double* quadWeights) {
  if (numQuadNodes!=0) {
       Finley_setError(VALUE_ERROR,"Finley_Quad_getNodesPoint: There is no quadrature scheme on points.");
  }
}

/**************************************************************/

/*   get a quadrature scheme with numQuadNodes quadrature nodes on the line [0,1]: */
/*   The nodes and weights are set from a table. */

void Finley_Quad_getNodesLine(int numQuadNodes,double* quadNodes,double* quadWeights) {
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

      case 3:
        quadNodes[0]=(1.-.774596669241483)/2.;
        quadNodes[1]=.5;
        quadNodes[2]=(1.+.774596669241483)/2.;
        quadWeights[0]=5./18.;
        quadWeights[1]=4./ 9.;
        quadWeights[2]=5./18.;
        break;

      case 4:
        quadNodes[0]=(1.-.861136311594053)/2.;
        quadNodes[1]=(1.-.339981043584856)/2.;
        quadNodes[2]=(1.+.339981043584856)/2.;
        quadNodes[3]=(1.+.861136311594053)/2.;
        quadWeights[0]=.347854845137454/2.;
        quadWeights[1]=.652145154862546/2.;
        quadWeights[2]=.652145154862546/2.;
        quadWeights[3]=.347854845137454/2.;
        break;

      case 5:
        quadNodes[0]=(1.-.906179845938664)/2.;
        quadNodes[1]=(1.-.538469310105683)/2.;
        quadNodes[2]= .5;
        quadNodes[3]=(1.+.538469310105683)/2.;
        quadNodes[4]=(1.+.906179845938664)/2.;
        quadWeights[0]=.236926885056189/2.;
        quadWeights[1]=.478628670499366/2.;
        quadWeights[2]=.568888888888889/2.;
        quadWeights[3]=.478628670499366/2.;
        quadWeights[4]=.236926885056189/2.;
        break;

      case 6:
        quadNodes[0]=(1.-.932469514203152)/2.;
        quadNodes[1]=(1.-.661209386466265)/2.;
        quadNodes[2]=(1.-.238619186083197)/2.;
        quadNodes[3]=(1.+.238619186083197)/2.;
        quadNodes[4]=(1.+.661209386466265)/2.;
        quadNodes[5]=(1.+.932469514203152)/2.;
        quadWeights[0]=.171324492379170/2.;
        quadWeights[1]=.360761573048139/2.;
        quadWeights[2]=.467913934572691/2.;
        quadWeights[3]=.467913934572691/2.;
        quadWeights[4]=.360761573048139/2.;
        quadWeights[5]=.171324492379170/2.;
        break;

      case 7:
        quadNodes[0]=(1.-.949107912342759)/2.;
        quadNodes[1]=(1.-.741531185599394)/2.;
        quadNodes[2]=(1.-.405845151377397)/2.;
        quadNodes[3]=0.5;
        quadNodes[4]=(1.+.405845151377397)/2.;
        quadNodes[5]=(1.+.741531185599394)/2.;
        quadNodes[6]=(1.+.949107912342759)/2.;
        quadWeights[0]= .129484966168870/2.;
        quadWeights[1]= .279705391489277/2.;
        quadWeights[2]= .381830050505119/2.;
        quadWeights[3]= .417959183673469/2.;
        quadWeights[4]= .381830050505119/2.;
        quadWeights[5]= .279705391489277/2.;
        quadWeights[6]= .129484966168870/2.;
        break;

      case 8:
        quadNodes[0]=(1.-.960289856497536)/2.;
        quadNodes[1]=(1.-.796666477413627)/2.;
        quadNodes[2]=(1.-.525532409916329)/2.;
        quadNodes[3]=(1.-.183434642495650)/2.;
        quadNodes[4]=(1.+.183434642495650)/2.;
        quadNodes[5]=(1.+.525532409916329)/2.;
        quadNodes[6]=(1.+.796666477413627)/2.;
        quadNodes[7]=(1.+.960289856497536)/2.;
        quadWeights[0]= .101228536290376/2.;
        quadWeights[1]= .222381034453374/2.;
        quadWeights[2]= .313706645877887/2.;
        quadWeights[3]= .362683783378362/2.;
        quadWeights[4]= .362683783378362/2.;
        quadWeights[5]= .313706645877887/2.;
        quadWeights[6]= .222381034453374/2.;
        quadWeights[7]= .101228536290376/2.;
        break;

      case 9:
        quadNodes[0]=(1.-.968160239507626)/2.;
        quadNodes[1]=(1.-.836031107326636)/2.;
        quadNodes[2]=(1.-.613371432700590)/2.;
        quadNodes[3]=(1.-.324253423403809)/2.;
        quadNodes[4]= .5;
        quadNodes[5]=(1.+.324253423403809)/2.;
        quadNodes[6]=(1.+.613371432700590)/2.;
        quadNodes[7]=(1.+.836031107326636)/2.;
        quadNodes[8]=(1.+.968160239507626)/2.;
        quadWeights[0]= .081274388361574/2.;
        quadWeights[1]= .180648160694857/2.;
        quadWeights[2]= .260610696402935/2.;
        quadWeights[3]= .312347077040003/2.;
        quadWeights[4]= .330239355001260/2.;
        quadWeights[5]= .312347077040003/2.;
        quadWeights[6]= .260610696402935/2.;
        quadWeights[7]= .180648160694857/2.;
        quadWeights[8]= .081274388361574/2.;
        break;

      case 10:
        quadNodes[0]=(1.-.973906528517172)/2.;
        quadNodes[1]=(1.-.865063366688985)/2.;
        quadNodes[2]=(1.-.679409568299024)/2.;
        quadNodes[3]=(1.-.433395394129247)/2.;
        quadNodes[4]=(1.-.148874338981631)/2.;
        quadNodes[5]=(1.+.148874338981631)/2.;
        quadNodes[6]=(1.+.433395394129247)/2.;
        quadNodes[7]=(1.+.679409568299024)/2.;
        quadNodes[8]=(1.+.865063366688985)/2.;
        quadNodes[9]=(1.+.973906528517172)/2.;
        quadWeights[0]= .066671344308688/2.;
        quadWeights[1]= .149451349150581/2.;
        quadWeights[2]= .219086362515982/2.;
        quadWeights[3]= .269266719309996/2.;
        quadWeights[4]= .295524224714753/2.;
        quadWeights[5]= .295524224714753/2.;
        quadWeights[6]= .269266719309996/2.;
        quadWeights[7]= .219086362515982/2.;
        quadWeights[8]= .149451349150581/2.;
        quadWeights[9]= .066671344308688/2.;
        break;

      default:
        Finley_setError(VALUE_ERROR,"__FILE__: Negative intergration order.");
        break;
  }
}
/**************************************************************/
/*                                                            */
/*  the following function are used define the meshes on the surface in the xy-plane */

/* triangle surface on a tetrahedron */
void Finley_Quad_getNodesTriOnFace(int numQuadNodes,double* quadNodes,double* quadWeights) {
       Finley_Quad_makeNodesOnFace(3,numQuadNodes,quadNodes,quadWeights,Finley_Quad_getNodesTri);
}
/* rectangular surface on a hexahedron */
void Finley_Quad_getNodesRecOnFace(int numQuadNodes,double* quadNodes,double* quadWeights) {
       Finley_Quad_makeNodesOnFace(3,numQuadNodes,quadNodes,quadWeights,Finley_Quad_getNodesRec);
}
/* line surface on a triangle or rectangle */
void Finley_Quad_getNodesLineOnFace(int numQuadNodes,double* quadNodes,double* quadWeights) {
       Finley_Quad_makeNodesOnFace(2,numQuadNodes,quadNodes,quadWeights,Finley_Quad_getNodesLine);
}
/* point surface on a line */
void Finley_Quad_getNodesPointOnFace(int numQuadNodes,double* quadNodes,double* quadWeights) {
     Finley_Quad_makeNodesOnFace(1,numQuadNodes,quadNodes,quadWeights,Finley_Quad_getNodesPoint);
}

void Finley_Quad_makeNodesOnFace(int dim, int numQuadNodes,double* quadNodes,double* quadWeights, Finley_Quad_getNodes getFaceNodes) {
    int q,i;
    double *quadNodesOnFace=NULL;
    #define DIM dim
    quadNodesOnFace=TMPMEMALLOC(numQuadNodes*(dim-1),double);

    if (! Finley_checkPtr(quadNodesOnFace) ) {
       getFaceNodes(numQuadNodes,quadNodesOnFace,quadWeights);
       if (Finley_noError()) {
          for (q=0;q<numQuadNodes;q++) {
             for (i=0;i<dim-1;i++) QUADNODES(i,q)=quadNodesOnFace[INDEX2(i,q,dim-1)];
             QUADNODES(dim-1,q)=0;
          }
       }
       TMPMEMFREE(quadNodesOnFace);
    }
    #undef DIM
}

/**************************************************************/

/*    The following functions Finley_Quad_getNumNodes* return the nmber of quadrature points needed to */
/*    achieve a certain accuracy. Notice that for Tet and Tri the the order is increased */
/*    to consider the accuracy reduction through the construction process.  */


int Finley_Quad_getNumNodesPoint(int order) {
  if (order <0 ) {
    Finley_setError(VALUE_ERROR,"Finley_Quad_getNumNodesPoint: Negative intergration order.");
    return -1;
  } else { 
    return 0;
  }
}

int Finley_Quad_getNumNodesLine(int order) {
  char error_msg[LenErrorMsg_MAX];
  if (order <0 ) {
    Finley_setError(VALUE_ERROR,"Finley_Quad_getNumNodesPoint: Negative intergration order.");
    return -1;
  } else { 
    if (order > 2*MAX_numQuadNodesLine-1) {
      sprintf(error_msg,"Finley_Quad_getNumNodesPoint: requested integration order %d on line is too large (>%d).",
                                                           order,2*MAX_numQuadNodesLine-1);
      Finley_setError(VALUE_ERROR,error_msg);
      return -1;
    } else { 
       Finley_resetError();
       return order/2+1;
    }
  }
}

int Finley_Quad_getNumNodesTri(int order) {
  int numQuadNodesLine;
  if (order<=1) {
      return 1;
  } else {
      numQuadNodesLine=Finley_Quad_getNumNodesLine(order+1);
      if (Finley_noError()) {
         return numQuadNodesLine*numQuadNodesLine;
      } else {
         return -1;
      }
  }
}

int Finley_Quad_getNumNodesRec(int order) {
  int numQuadNodesLine;
  numQuadNodesLine=Finley_Quad_getNumNodesLine(order);
  if (Finley_noError()) {
      return numQuadNodesLine*numQuadNodesLine;
  } else {
      return -1;
  }
}

int Finley_Quad_getNumNodesTet(int order) {
  int numQuadNodesLine;
  if (order<=1) {
      return 1;
  } else {
     numQuadNodesLine=Finley_Quad_getNumNodesLine(order+2);
     if (Finley_noError()) {
         return numQuadNodesLine*numQuadNodesLine*numQuadNodesLine;
     } else {
         return -1;
     }
  }
}

int Finley_Quad_getNumNodesHex(int order) {
  int numQuadNodesLine;
  numQuadNodesLine=Finley_Quad_getNumNodesLine(order);
  if (Finley_noError()) {
      return numQuadNodesLine*numQuadNodesLine*numQuadNodesLine;
  } else {
      return -1;
  }
}
/* 
* $Log$
* Revision 1.2  2005/09/15 03:44:23  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.1.1.1.6.1  2005/09/07 06:26:20  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.2  2004/08/03 04:49:06  gross
* bug in Quadrature.c fixed
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/
