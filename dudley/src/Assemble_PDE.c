
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

/*    assembles the system of numEq PDEs into the stiffness matrix S and right hand side F */

/*     -div(A*grad u)-div(B*u)+C*grad u + D*u= -div X + Y */

/*      -(A_{k,i,m,j} u_m,j)_i-(B_{k,i,m} u_m)_i+C_{k,m,j} u_m,j-D_{k,m} u_m = -(X_{k,i})_i + Y_k */

/*    u has numComp components. */

/*    Shape of the coefficients: */

/*      A = numEqu x numDim x numComp x numDim */
/*      B = numDim x numEqu x numComp  */
/*      C = numEqu x numDim x numComp  */
/*      D = numEqu x numComp  */
/*      X = numEqu x numDim   */
/*      Y = numEqu */

/*    The coefficients A,B,C,D,X and Y have to be defined on the integartion points or not present (=NULL). */

/*    S and F have to be initialized before the routine is called. S or F can be NULL. In this case the left or */
/*    the right hand side of the PDE is not processed.  */

/*    The routine does not consider any boundary conditions. */

/**************************************************************/

#include "Assemble.h"
#include "Util.h"
#include "esysUtils/blocktimer.h"
#ifdef _OPENMP
#include <omp.h>
#endif


/**************************************************************/

void Dudley_Assemble_PDE(Dudley_NodeFile* nodes,Dudley_ElementFile* elements,Paso_SystemMatrix* S, escriptDataC* F,
			 escriptDataC* A, escriptDataC* B, escriptDataC* C, escriptDataC* D, escriptDataC* X, escriptDataC* Y ) {

  bool_t reducedIntegrationOrder=FALSE;
  char error_msg[LenErrorMsg_MAX];
  Assemble_Parameters p;
  double time0;
  dim_t dimensions[ESCRIPT_MAX_DATA_RANK];
  type_t funcspace;
  double blocktimer_start = blocktimer_time();

  Dudley_resetError();

  {
#ifdef Paso_MPI
  int iam, numCPUs;
  iam = elements->elementDistribution->MPIInfo->rank;
  numCPUs = elements->elementDistribution->MPIInfo->size;
#endif
  }

  if (nodes==NULL || elements==NULL) return;
  if (S==NULL && isEmpty(F)) return;

  if (isEmpty(F) && !isEmpty(X) && !isEmpty(F)) {
        Dudley_setError(TYPE_ERROR,"Dudley_Assemble_PDE: right hand side coefficients are non-zero bat no right hand side vector given.");
  }

  if (S==NULL && !isEmpty(A) && !isEmpty(B) && !isEmpty(C) && !isEmpty(D)) {
        Dudley_setError(TYPE_ERROR,"Dudley_Assemble_PDE: coefficients are non-zero but no matrix is given.");
  }

  /*  get the functionspace for this assemblage call */
  funcspace=UNKNOWN;
  updateFunctionSpaceType(funcspace,A);
  updateFunctionSpaceType(funcspace,B);
  updateFunctionSpaceType(funcspace,C);
  updateFunctionSpaceType(funcspace,D);
  updateFunctionSpaceType(funcspace,X);
  updateFunctionSpaceType(funcspace,Y);
  if (funcspace==UNKNOWN) return; /* all  data are empty */

  /* check if all function spaces are the same */
  if (! functionSpaceTypeEqual(funcspace,A) ) {
        Dudley_setError(TYPE_ERROR,"Dudley_Assemble_PDE: unexpected function space type for coefficient A");
  }
  if (! functionSpaceTypeEqual(funcspace,B) ) {
        Dudley_setError(TYPE_ERROR,"Dudley_Assemble_PDE: unexpected function space type for coefficient B");
  }
  if (! functionSpaceTypeEqual(funcspace,C) ) {
        Dudley_setError(TYPE_ERROR,"Dudley_Assemble_PDE: unexpected function space type for coefficient C");
  }
  if (! functionSpaceTypeEqual(funcspace,D) ) {
        Dudley_setError(TYPE_ERROR,"Dudley_Assemble_PDE: unexpected function space type for coefficient D");
  }
  if (! functionSpaceTypeEqual(funcspace,X) ) {
        Dudley_setError(TYPE_ERROR,"Dudley_Assemble_PDE: unexpected function space type for coefficient X");
  }
  if (! functionSpaceTypeEqual(funcspace,Y) ) {
        Dudley_setError(TYPE_ERROR,"Dudley_Assemble_PDE: unexpected function space type for coefficient Y");
  }
  if (! Dudley_noError()) return;

  /* check if all function spaces are the same */
  if (funcspace==DUDLEY_ELEMENTS) {
       reducedIntegrationOrder=FALSE;
  } else if (funcspace==DUDLEY_FACE_ELEMENTS)  {
       reducedIntegrationOrder=FALSE;
  } else if (funcspace==DUDLEY_REDUCED_ELEMENTS) {
       reducedIntegrationOrder=TRUE;
  } else if (funcspace==DUDLEY_REDUCED_FACE_ELEMENTS)  {
       reducedIntegrationOrder=TRUE;
  } else {
       Dudley_setError(TYPE_ERROR,"Dudley_Assemble_PDE: assemblage failed because of illegal function space.");
  }
  if (! Dudley_noError()) return;

  /* set all parameters in p*/
  Assemble_getAssembleParameters(nodes,elements,S,F, reducedIntegrationOrder, &p);
  if (! Dudley_noError()) return;

  /* check if all function spaces are the same */

  if (! numSamplesEqual(A,p.numQuadTotal,elements->numElements) ) {
        sprintf(error_msg,"Dudley_Assemble_PDE: sample points of coefficient A don't match (%d,%d)",p.numQuadTotal,elements->numElements);
        Dudley_setError(TYPE_ERROR,error_msg);
  }

  if (! numSamplesEqual(B,p.numQuadTotal,elements->numElements) ) {
        sprintf(error_msg,"Dudley_Assemble_PDE: sample points of coefficient B don't match (%d,%d)",p.numQuadTotal,elements->numElements);
        Dudley_setError(TYPE_ERROR,error_msg);
  }

  if (! numSamplesEqual(C,p.numQuadTotal,elements->numElements) ) {
        sprintf(error_msg,"Dudley_Assemble_PDE: sample points of coefficient C don't match (%d,%d)",p.numQuadTotal,elements->numElements);
        Dudley_setError(TYPE_ERROR,error_msg);
  }

  if (! numSamplesEqual(D,p.numQuadTotal,elements->numElements) ) {
        sprintf(error_msg,"Dudley_Assemble_PDE: sample points of coefficient D don't match (%d,%d)",p.numQuadTotal,elements->numElements);
        Dudley_setError(TYPE_ERROR,error_msg);
  }

  if (! numSamplesEqual(X,p.numQuadTotal,elements->numElements) ) {
        sprintf(error_msg,"Dudley_Assemble_PDE: sample points of coefficient X don't match (%d,%d)",p.numQuadTotal,elements->numElements);
        Dudley_setError(TYPE_ERROR,error_msg);
  }

  if (! numSamplesEqual(Y,p.numQuadTotal,elements->numElements) ) {
        sprintf(error_msg,"Dudley_Assemble_PDE: sample points of coefficient Y don't match (%d,%d)",p.numQuadTotal,elements->numElements);
        Dudley_setError(TYPE_ERROR,error_msg);
  }

  /*  check the dimensions: */

  if (p.numEqu==1 && p.numComp==1) {
    if (!isEmpty(A)) {
      dimensions[0]=p.numDim;
      dimensions[1]=p.numDim;
      if (!isDataPointShapeEqual(A,2,dimensions)) {
          sprintf(error_msg,"Dudley_Assemble_PDE: coefficient A: illegal shape, expected shape (%d,%d)",dimensions[0],dimensions[1]);
          Dudley_setError(TYPE_ERROR,error_msg);
      }
    }
    if (!isEmpty(B)) {
       dimensions[0]=p.numDim;
       if (!isDataPointShapeEqual(B,1,dimensions)) {
          sprintf(error_msg,"Dudley_Assemble_PDE: coefficient B: illegal shape (%d,)",dimensions[0]);
          Dudley_setError(TYPE_ERROR,error_msg);
       }
    }
    if (!isEmpty(C)) {
       dimensions[0]=p.numDim;
       if (!isDataPointShapeEqual(C,1,dimensions)) {
          sprintf(error_msg,"Dudley_Assemble_PDE: coefficient C, expected shape (%d,)",dimensions[0]);
          Dudley_setError(TYPE_ERROR,error_msg);
       }
    }
    if (!isEmpty(D)) {
       if (!isDataPointShapeEqual(D,0,dimensions)) {
          Dudley_setError(TYPE_ERROR,"Dudley_Assemble_PDE: coefficient D, rank 0 expected.");
       }
    }
    if (!isEmpty(X)) {
       dimensions[0]=p.numDim;
       if (!isDataPointShapeEqual(X,1,dimensions)) {
          sprintf(error_msg,"Dudley_Assemble_PDE: coefficient X, expected shape (%d,",dimensions[0]);
          Dudley_setError(TYPE_ERROR,error_msg);
       }
    }
    if (!isEmpty(Y)) {
       if (!isDataPointShapeEqual(Y,0,dimensions)) {
          Dudley_setError(TYPE_ERROR,"Dudley_Assemble_PDE: coefficient Y, rank 0 expected.");
       }
    }
  } else {
    if (!isEmpty(A)) {
      dimensions[0]=p.numEqu;
      dimensions[1]=p.numDim;
      dimensions[2]=p.numComp;
      dimensions[3]=p.numDim;
      if (!isDataPointShapeEqual(A,4,dimensions)) {
          sprintf(error_msg,"Dudley_Assemble_PDE: coefficient A, expected shape (%d,%d,%d,%d)",dimensions[0],dimensions[1],dimensions[2],dimensions[3]);
          Dudley_setError(TYPE_ERROR,error_msg);
      }
    }
    if (!isEmpty(B)) {
      dimensions[0]=p.numEqu;
      dimensions[1]=p.numDim;
      dimensions[2]=p.numComp;
      if (!isDataPointShapeEqual(B,3,dimensions)) {
          sprintf(error_msg,"Dudley_Assemble_PDE: coefficient B, expected shape (%d,%d,%d)",dimensions[0],dimensions[1],dimensions[2]);
          Dudley_setError(TYPE_ERROR,error_msg);
      }
    }
    if (!isEmpty(C)) {
      dimensions[0]=p.numEqu;
      dimensions[1]=p.numComp;
      dimensions[2]=p.numDim;
      if (!isDataPointShapeEqual(C,3,dimensions)) {
          sprintf(error_msg,"Dudley_Assemble_PDE: coefficient C, expected shape (%d,%d,%d)",dimensions[0],dimensions[1],dimensions[2]);
          Dudley_setError(TYPE_ERROR,error_msg);
      }
    }
    if (!isEmpty(D)) {
      dimensions[0]=p.numEqu;
      dimensions[1]=p.numComp;
      if (!isDataPointShapeEqual(D,2,dimensions)) {
          sprintf(error_msg,"Dudley_Assemble_PDE: coefficient D, expected shape (%d,%d)",dimensions[0],dimensions[1]);
          Dudley_setError(TYPE_ERROR,error_msg);
      }
    }
    if (!isEmpty(X)) {
      dimensions[0]=p.numEqu;
      dimensions[1]=p.numDim;
      if (!isDataPointShapeEqual(X,2,dimensions)) {
          sprintf(error_msg,"Dudley_Assemble_PDE: coefficient X, expected shape (%d,%d)",dimensions[0],dimensions[1]);
          Dudley_setError(TYPE_ERROR,error_msg);
      }
    }
    if (!isEmpty(Y)) {
      dimensions[0]=p.numEqu;
      if (!isDataPointShapeEqual(Y,1,dimensions)) {
          sprintf(error_msg,"Dudley_Assemble_PDE: coefficient Y, expected shape (%d,)",dimensions[0]);
          Dudley_setError(TYPE_ERROR,error_msg);
      }
    }
  }
  if (Dudley_noError()) {
     time0=Dudley_timer();
     if (p.numEqu == p. numComp) {
        if (p.numEqu > 1) {
          /* system of PDESs */
          if (p.numDim==3) {
               Dudley_Assemble_PDE_System2_3D(p,elements,S,F,A,B,C,D,X,Y);
          } else if (p.numDim==2) {
               Dudley_Assemble_PDE_System2_2D(p,elements,S,F,A,B,C,D,X,Y);
          } else {
            Dudley_setError(VALUE_ERROR,"Dudley_Assemble_PDE supports spatial dimensions 2 and 3 only.");
          }
        } else {
          /* single PDES */
          if (p.numDim==3) {
               Dudley_Assemble_PDE_Single2_3D(p,elements,S,F,A,B,C,D,X,Y);
          } else if (p.numDim==2) {
               Dudley_Assemble_PDE_Single2_2D(p,elements,S,F,A,B,C,D,X,Y);
          } else {
            Dudley_setError(VALUE_ERROR,"Dudley_Assemble_PDE supports spatial dimensions 2 and 3 only.");
          }
        }
     } else {
          Dudley_setError(VALUE_ERROR,"Dudley_Assemble_PDE requires number of equations == number of solutions  .");
     }
  }
  blocktimer_increment("Dudley_Assemble_PDE()", blocktimer_start);
}
