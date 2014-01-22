
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/****************************************************************************

  Assembles the system of numEqu PDEs into the stiffness matrix S and right
  hand side F:

      -div(A*grad u)-div(B*u)+C*grad u + D*u= -div X + Y

      -(A_{k,i,m,j} u_m,j)_i-(B_{k,i,m} u_m)_i+C_{k,m,j} u_m,j-D_{k,m} u_m = -(X_{k,i})_i + Y_k

  u has numComp components.
  Shape of the coefficients:

      A = numEqu x numDim x numComp x numDim
      B = numDim x numEqu x numComp
      C = numEqu x numDim x numComp
      D = numEqu x numComp
      X = numEqu x numDim
      Y = numEqu

  The coefficients A,B,C,D,X and Y have to be defined on the integration
  points or not present (i.e. empty).

  S and F have to be initialized before the routine is called. S or F can
  be NULL. In this case the left or the right hand side of the PDE is not
  processed.

  The routine does not consider any boundary conditions.

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"
#include "esysUtils/blocktimer.h"

#include <sstream>

namespace finley {

inline void setNumSamplesError(const char* c, int n0, int n1)
{
    std::stringstream ss;
    ss << "Assemble_PDE: number of sample points of coefficient " << c
        << " don't match (" << n0 << "," << n1 << ").";
    std::string errorMsg(ss.str());
    setError(TYPE_ERROR, errorMsg.c_str());
}

inline void setShapeError(const char* c, int num, const int *dims)
{
    std::stringstream ss;
    ss << "Assemble_PDE: shape of coefficient " << c
        << " does not match (" << dims[0] << ",";
    if (num > 1) {
       ss << dims[1];
       if (num > 2) {
           ss << "," << dims[2];
           if (num > 3) {
               ss << "," << dims[3];
           }
       }
    }
    ss << ").";
    std::string errorMsg(ss.str());
    setError(TYPE_ERROR, errorMsg.c_str());
}

void Assemble_PDE(const NodeFile* nodes, const ElementFile* elements,
                  Paso_SystemMatrix* S, escript::Data& F,
                  const escript::Data& A, const escript::Data& B,
                  const escript::Data& C, const escript::Data& D,
                  const escript::Data& X, const escript::Data& Y)
{
    resetError();
    if (!nodes || !elements || (S==NULL && F.isEmpty()))
        return;

    if (F.isEmpty() && (!X.isEmpty() || !Y.isEmpty())) {
        setError(TYPE_ERROR, "Assemble_PDE: right hand side coefficients are non-zero but no right hand side vector given.");
        return;
    }

    if (S==NULL && !A.isEmpty() && !B.isEmpty() && !C.isEmpty() && !D.isEmpty()) {
        setError(TYPE_ERROR, "Assemble_PDE: coefficients are non-zero but no matrix is given.");
        return;
    }

    // get the function space for this assemblage call
    int funcspace = UNKNOWN;
    if (!A.isEmpty()) funcspace=A.getFunctionSpace().getTypeCode();
    if (!B.isEmpty()) funcspace=B.getFunctionSpace().getTypeCode();
    if (!C.isEmpty()) funcspace=C.getFunctionSpace().getTypeCode();
    if (!D.isEmpty()) funcspace=D.getFunctionSpace().getTypeCode();
    if (!X.isEmpty()) funcspace=X.getFunctionSpace().getTypeCode();
    if (!Y.isEmpty()) funcspace=Y.getFunctionSpace().getTypeCode();
    if (funcspace==UNKNOWN)
        return; // all data are empty

    // check if all function spaces are the same
    if (!A.isEmpty() && A.getFunctionSpace().getTypeCode()!=funcspace) {
        setError(TYPE_ERROR, "Assemble_PDE: unexpected function space type for coefficient A");
    } else if (!B.isEmpty() && B.getFunctionSpace().getTypeCode()!=funcspace) {
        setError(TYPE_ERROR, "Assemble_PDE: unexpected function space type for coefficient B");
    } else if (!C.isEmpty() && C.getFunctionSpace().getTypeCode()!=funcspace) {
        setError(TYPE_ERROR, "Assemble_PDE: unexpected function space type for coefficient C");
    } else if (!D.isEmpty() && D.getFunctionSpace().getTypeCode()!=funcspace) {
        setError(TYPE_ERROR, "Assemble_PDE: unexpected function space type for coefficient D");
    } else if (!X.isEmpty() && X.getFunctionSpace().getTypeCode()!=funcspace) {
        setError(TYPE_ERROR, "Assemble_PDE: unexpected function space type for coefficient X");
    } else if (!Y.isEmpty() && Y.getFunctionSpace().getTypeCode()!=funcspace) {
        setError(TYPE_ERROR, "Assemble_PDE: unexpected function space type for coefficient Y");
    }
    if (!noError())
        return;

    bool reducedIntegrationOrder;
    if (funcspace==FINLEY_ELEMENTS) {
       reducedIntegrationOrder=false;
    } else if (funcspace==FINLEY_FACE_ELEMENTS)  {
       reducedIntegrationOrder=false;
    } else if (funcspace==FINLEY_CONTACT_ELEMENTS_1)  {
       reducedIntegrationOrder=false;
    } else if (funcspace==FINLEY_CONTACT_ELEMENTS_2)  {
       reducedIntegrationOrder=false;
    } else if (funcspace==FINLEY_REDUCED_ELEMENTS) {
       reducedIntegrationOrder=true;
    } else if (funcspace==FINLEY_REDUCED_FACE_ELEMENTS)  {
       reducedIntegrationOrder=true;
    } else if (funcspace==FINLEY_REDUCED_CONTACT_ELEMENTS_1)  {
       reducedIntegrationOrder=true;
    } else if (funcspace==FINLEY_REDUCED_CONTACT_ELEMENTS_2)  {
       reducedIntegrationOrder=true;
    } else if (funcspace==FINLEY_POINTS)  {
       reducedIntegrationOrder=false;
    } else {
       setError(TYPE_ERROR, "Assemble_PDE: assemblage failed because of illegal function space.");
       return;
    }

    // get assemblage parameters
    AssembleParameters p(nodes, elements, S, F, reducedIntegrationOrder);
    if (!noError())
        return;

    // check if sample numbers are the same
    if (!A.numSamplesEqual(p.numQuadTotal, elements->numElements)) {
        setNumSamplesError("A", p.numQuadTotal, elements->numElements);
    } else if (!B.numSamplesEqual(p.numQuadTotal, elements->numElements)) {
        setNumSamplesError("B", p.numQuadTotal, elements->numElements);
    } else if (!C.numSamplesEqual(p.numQuadTotal, elements->numElements)) {
        setNumSamplesError("C", p.numQuadTotal, elements->numElements);
    } else if (!D.numSamplesEqual(p.numQuadTotal, elements->numElements)) {
        setNumSamplesError("D", p.numQuadTotal, elements->numElements);
    } else if (!X.numSamplesEqual(p.numQuadTotal, elements->numElements)) {
        setNumSamplesError("X", p.numQuadTotal, elements->numElements);
    } else if (!Y.numSamplesEqual(p.numQuadTotal, elements->numElements)) {
        setNumSamplesError("Y", p.numQuadTotal, elements->numElements);
    }
    if (!noError())
        return;

    // check the dimensions:
    if (p.numEqu != p. numComp) {
        setError(VALUE_ERROR, "Assemble_PDE requires number of equations == number of solutions.");
    } else if (p.numEqu==1) {
        const int dimensions[2] = { p.numDim, p.numDim };
        if (!A.isDataPointShapeEqual(2, dimensions)) {
            setShapeError("A", 2, dimensions);
        } else if (!B.isDataPointShapeEqual(1, dimensions)) {
            setShapeError("B", 1, dimensions);
        } else if (!C.isDataPointShapeEqual(1, dimensions)) {
            setShapeError("C", 1, dimensions);
        } else if (!D.isDataPointShapeEqual(0, dimensions)) {
            setError(TYPE_ERROR, "Assemble_PDE: coefficient D must be rank 0.");
        } else if (!X.isDataPointShapeEqual(1, dimensions)) {
            setShapeError("X", 1, dimensions);
        } else if (!Y.isDataPointShapeEqual(0, dimensions)) {
            setError(TYPE_ERROR, "Assemble_PDE: coefficient Y must be rank 0.");
        }
    } else {
        const int dimAB[4] = { p.numEqu, p.numDim, p.numComp, p.numDim };
        const int dimCD[3] = { p.numEqu, p.numComp, p.numDim };
        if (!A.isDataPointShapeEqual(4, dimAB)) {
            setShapeError("A", 4, dimAB);
        } else if (!B.isDataPointShapeEqual(3, dimAB)) {
            setShapeError("B", 3, dimAB);
        } else if (!C.isDataPointShapeEqual(3, dimCD)) {
            setShapeError("C", 3, dimCD);
        } else if (!D.isDataPointShapeEqual(2, dimCD)) {
            setShapeError("D", 2, dimCD);
        } else if (!X.isDataPointShapeEqual(2, dimAB)) {
            setShapeError("X", 2, dimAB);
        } else if (!Y.isDataPointShapeEqual(1, dimAB)) {
            setShapeError("Y", 1, dimAB);
        }
    }
    if (!noError())
        return;

    double blocktimer_start = blocktimer_time();

    if (p.numSides == 1) {
        if (funcspace==FINLEY_POINTS) {
            if (!A.isEmpty() || !B.isEmpty() || !C.isEmpty() || !X.isEmpty()) {
                setError(TYPE_ERROR, "Assemble_PDE: Point elements require A, B, C and X to be empty.");
            } else {
                Assemble_PDE_Points(p, D, Y);
            }
        } else if (p.numEqu > 1) { // system of PDEs
            if (p.numDim==3) {
                Assemble_PDE_System_3D(p, A, B, C, D, X, Y);
            } else if (p.numDim==2) {
                Assemble_PDE_System_2D(p, A, B, C, D, X, Y);
            } else if (p.numDim==1) {
                Assemble_PDE_System_1D(p, A, B, C, D, X, Y);
            } else {
                setError(VALUE_ERROR, "Assemble_PDE supports spatial dimensions 1,2,3 only.");
            }
        } else { // single PDE
            if (p.numDim==3) {
                Assemble_PDE_Single_3D(p, A, B, C, D, X, Y);
            } else if (p.numDim==2) {
                Assemble_PDE_Single_2D(p, A, B, C, D, X, Y);
            } else if (p.numDim==1) {
                Assemble_PDE_Single_1D(p, A, B, C, D, X, Y);
            } else {
                setError(VALUE_ERROR, "Assemble_PDE supports spatial dimensions 1,2,3 only.");
            }
        }
    } else if (p.numSides == 2) {
        if (!A.isEmpty() || !B.isEmpty() || !C.isEmpty() || !X.isEmpty()) {
            setError(TYPE_ERROR, "Assemble_PDE: Contact elements require A, B, C and X to be empty.");
        } else if (p.numEqu > 1) { // system of PDEs
            Assemble_PDE_System_C(p, D, Y);
        } else { // single PDE
            Assemble_PDE_Single_C(p, D, Y);
        }
    } else {
        setError(TYPE_ERROR,"Assemble_PDE supports numShape=NumNodes or 2*numShape=NumNodes only.");
    }

    blocktimer_increment("Assemble_PDE()", blocktimer_start);
}

} // namespace finley

