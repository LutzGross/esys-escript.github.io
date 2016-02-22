
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

/************************************************************************************/

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

/*    The coefficients A,B,C,D,X and Y have to be defined on the integration points or not present (=NULL). */

/*    S and F have to be initialised before the routine is called. S or F can be NULL. In this case the left or */
/*    the right hand side of the PDE is not processed.  */

/*    The routine does not consider any boundary conditions. */

/************************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "Assemble.h"
#include "Util.h"

inline void setNumSamplesError(const char* c, int n0, int n1)
{
    std::stringstream ss;
    ss << "Assemble_PDE: number of sample points of coefficient " << c
        << " don't match (" << n0 << "," << n1 << ").";
    std::string errorMsg(ss.str());
    Dudley_setError(TYPE_ERROR, errorMsg.c_str());
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
    Dudley_setError(TYPE_ERROR, errorMsg.c_str());
}

void Dudley_Assemble_PDE(Dudley_NodeFile* nodes, Dudley_ElementFile* elements,
                         escript::ASM_ptr S, escript::Data* F,
                         const escript::Data* A, const escript::Data* B, const escript::Data* C,
                         const escript::Data* D, const escript::Data* X, const escript::Data* Y)
{
    bool reducedIntegrationOrder = false;
    Dudley_Assemble_Parameters p;
    int funcspace;

    Dudley_resetError();

    if (nodes == NULL || elements == NULL)
        return;
    if (S == NULL && isEmpty(F))
        return;
    if (isEmpty(F) && ( !isEmpty(X) || !isEmpty(Y) ) ) 
    {
        Dudley_setError(TYPE_ERROR,
                        "Dudley_Assemble_PDE: right hand side coefficients are non-zero but no right hand side vector given.");
    }

    if (S == NULL && !isEmpty(A) && !isEmpty(B) && !isEmpty(C) && !isEmpty(D))
    {
        Dudley_setError(TYPE_ERROR, "Dudley_Assemble_PDE: coefficients are non-zero but no matrix is given.");
    }

    /*  get the functionspace for this assemblage call */
    funcspace = -1;
    updateFunctionSpaceType(funcspace, A);
    updateFunctionSpaceType(funcspace, B);
    updateFunctionSpaceType(funcspace, C);
    updateFunctionSpaceType(funcspace, D);
    updateFunctionSpaceType(funcspace, X);
    updateFunctionSpaceType(funcspace, Y);
    if (funcspace == -1)
        return;                 /* all  data are empty */

    /* check if all function spaces are the same */
    if (!functionSpaceTypeEqual(funcspace, A))
    {
        Dudley_setError(TYPE_ERROR, "Dudley_Assemble_PDE: unexpected function space type for coefficient A");
    }
    if (!functionSpaceTypeEqual(funcspace, B))
    {
        Dudley_setError(TYPE_ERROR, "Dudley_Assemble_PDE: unexpected function space type for coefficient B");
    }
    if (!functionSpaceTypeEqual(funcspace, C))
    {
        Dudley_setError(TYPE_ERROR, "Dudley_Assemble_PDE: unexpected function space type for coefficient C");
    }
    if (!functionSpaceTypeEqual(funcspace, D))
    {
        Dudley_setError(TYPE_ERROR, "Dudley_Assemble_PDE: unexpected function space type for coefficient D");
    }
    if (!functionSpaceTypeEqual(funcspace, X))
    {
        Dudley_setError(TYPE_ERROR, "Dudley_Assemble_PDE: unexpected function space type for coefficient X");
    }
    if (!functionSpaceTypeEqual(funcspace, Y))
    {
        Dudley_setError(TYPE_ERROR, "Dudley_Assemble_PDE: unexpected function space type for coefficient Y");
    }
    if (!Dudley_noError())
        return;

    /* check if all function spaces are the same */
    if (funcspace == DUDLEY_ELEMENTS)
    {
        reducedIntegrationOrder = false;
    }
    else if (funcspace == DUDLEY_FACE_ELEMENTS)
    {
        reducedIntegrationOrder = false;
    }
    else if (funcspace == DUDLEY_REDUCED_ELEMENTS)
    {
        reducedIntegrationOrder = true;
    }
    else if (funcspace == DUDLEY_REDUCED_FACE_ELEMENTS)
    {
        reducedIntegrationOrder = true;
    }
    else if (funcspace == DUDLEY_POINTS)
    {
        reducedIntegrationOrder = true;
    }
    else
    {
        Dudley_setError(TYPE_ERROR, "Dudley_Assemble_PDE: assemblage failed because of illegal function space.");
    }
    if (!Dudley_noError())
        return;

    /* set all parameters in p */
    Dudley_Assemble_getAssembleParameters(nodes, elements, S, F, reducedIntegrationOrder, &p);
    if (!Dudley_noError())
        return;

    // check if all function spaces are the same

    if (!numSamplesEqual(A, p.numQuad, elements->numElements)) {
        setNumSamplesError("A", p.numQuad, elements->numElements);
    } else if (!numSamplesEqual(B, p.numQuad, elements->numElements)) {
        setNumSamplesError("B", p.numQuad, elements->numElements);
    } else if (!numSamplesEqual(C, p.numQuad, elements->numElements)) {
        setNumSamplesError("C", p.numQuad, elements->numElements);
    } else if (!numSamplesEqual(D, p.numQuad, elements->numElements)) {
        setNumSamplesError("D", p.numQuad, elements->numElements);
    } else if (!numSamplesEqual(X, p.numQuad, elements->numElements)) {
        setNumSamplesError("X", p.numQuad, elements->numElements);
    } else if (!numSamplesEqual(Y, p.numQuad, elements->numElements)) {
        setNumSamplesError("Y", p.numQuad, elements->numElements);
    }

    // check the dimensions

    if (p.numEqu == 1 && p.numComp == 1)
    {
        const int dimensions[2] = { p.numDim, p.numDim };
        if (!isEmpty(A)) {
            if (!isDataPointShapeEqual(A, 2, dimensions)) {
                setShapeError("A", 2, dimensions);
            }
        }
        if (!isEmpty(B)) {
            if (!isDataPointShapeEqual(B, 1, dimensions)) {
                setShapeError("B", 1, dimensions);
            }
        }
        if (!isEmpty(C)) {
            if (!isDataPointShapeEqual(C, 1, dimensions)) {
                setShapeError("C", 1, dimensions);
            }
        }
        if (!isEmpty(D)) {
            if (!isDataPointShapeEqual(D, 0, dimensions)) {
                Dudley_setError(TYPE_ERROR, "Assemble_PDE: coefficient D must be rank 0.");
            }
        }
        if (!isEmpty(X)) {
            if (!isDataPointShapeEqual(X, 1, dimensions)) {
                setShapeError("X", 1, dimensions);
            }
        }
        if (!isEmpty(Y)) {
            if (!isDataPointShapeEqual(Y, 0, dimensions)) {
                Dudley_setError(TYPE_ERROR, "Assemble_PDE: coefficient Y must be rank 0.");
            }
        }
    } else {
        const int dimAB[4] = { p.numEqu, p.numDim, p.numComp, p.numDim };
        const int dimCD[3] = { p.numEqu, p.numComp, p.numDim };
        if (!isEmpty(A)) {
            if (!isDataPointShapeEqual(A, 4, dimAB)) {
                setShapeError("A", 4, dimAB);
            }
        }
        if (!isEmpty(B)) {
            if (!isDataPointShapeEqual(B, 3, dimAB)) {
                setShapeError("B", 3, dimAB);
            }
        }
        if (!isEmpty(C)) {
            if (!isDataPointShapeEqual(C, 3, dimCD)) {
                setShapeError("C", 3, dimCD);
            }
        }
        if (!isEmpty(D)) {
            if (!isDataPointShapeEqual(D, 2, dimCD)) {
                setShapeError("D", 2, dimCD);
            }
        }
        if (!isEmpty(X)) {
            if (!isDataPointShapeEqual(X, 2, dimAB)) {
                setShapeError("X", 2, dimAB);
            }
        }
        if (!isEmpty(Y)) {
            if (!isDataPointShapeEqual(Y, 1, dimAB)) {
                setShapeError("Y", 1, dimAB);
            }
        }
    }
    if (Dudley_noError())
    {
        if (funcspace==DUDLEY_POINTS) {
            if ( !isEmpty(A) || !isEmpty(B) || !isEmpty(C) || !isEmpty(X) ) {
                Dudley_setError(TYPE_ERROR,"Finley_Assemble_PDE: Point elements require A, B, C and X to be empty.");
            } else {
                Dudley_Assemble_PDE_Points(p, elements,S,F, D, Y);
            }
        }
        else
        {
            if (p.numEqu == p.numComp)
            {
                if (p.numEqu > 1)
                {
                    // system of PDEs
                    if (p.numDim == 3)
                    {
                        Dudley_Assemble_PDE_System2_3D(p, elements, S, F, A, B, C, D, X, Y);
                    }
                    else if (p.numDim == 2)
                    {
                        Dudley_Assemble_PDE_System2_2D(p, elements, S, F, A, B, C, D, X, Y);
                    }
                    else
                    {
                        Dudley_setError(VALUE_ERROR, "Dudley_Assemble_PDE supports spatial dimensions 2 and 3 only.");
                    }
                }
                else
                {
                    // single PDE
                    if (p.numDim == 3)
                    {
                        Dudley_Assemble_PDE_Single2_3D(p, elements, S, F, A, B, C, D, X, Y);
                    }
                    else if (p.numDim == 2)
                    {
                        Dudley_Assemble_PDE_Single2_2D(p, elements, S, F, A, B, C, D, X, Y);
                    }
                    else
                    {
                        Dudley_setError(VALUE_ERROR, "Dudley_Assemble_PDE supports spatial dimensions 2 and 3 only.");
                    }
                }
            }
            else
            {
                Dudley_setError(VALUE_ERROR, "Dudley_Assemble_PDE requires number of equations == number of solutions  .");
            }
        }
    }
}

