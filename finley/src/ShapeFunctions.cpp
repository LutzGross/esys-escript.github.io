
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


/****************************************************************************

  Finley: Shape functions

*****************************************************************************/

#include "ShapeFunctions.h"

#include <escript/index.h>

#include <cstring>

namespace finley {

const ShapeFunctionInfo ShapeFunction_InfoList[] = {
    { Point1Shape, "Point1", 0,  1, 1, 1, Shape_Point1 },
    { Line2Shape,  "Line2",  1,  2, 1, 2, Shape_Line2  },
    { Line3Shape,  "Line3",  1,  3, 2, 2, Shape_Line3  },
    { Line4Shape,  "Line4",  1,  4, 3, 2, Shape_Line4  },
    { Tri3Shape,   "Tri3",   2,  3, 1, 3, Shape_Tri3   },
    { Tri6Shape,   "Tri6",   2,  6, 2, 3, Shape_Tri6   },
    { Tri9Shape,   "Tri9",   2,  9, 3, 3, Shape_Tri9   },
    { Tri10Shape,  "Tri10",  2, 10, 3, 3, Shape_Tri10, },
    { Rec4Shape,   "Rec4",   2,  4, 1, 4, Shape_Rec4,  },
    { Rec8Shape,   "Rec8",   2,  8, 2, 4, Shape_Rec8,  },
    { Rec9Shape,   "Rec9",   2,  9, 2, 4, Shape_Rec9,  },
    { Rec12Shape,  "Rec12",  2, 12, 3, 4, Shape_Rec12, },
    { Rec16Shape,  "Rec16",  2, 16, 3, 4, Shape_Rec16, },
    { Tet4Shape,   "Tet4",   3,  4, 1, 4, Shape_Tet4,  },
    { Tet10Shape,  "Tet10",  3, 10, 2, 4, Shape_Tet10, },
    { Tet16Shape,  "Tet16",  3, 16, 3, 4, Shape_Tet16, },
    { Hex8Shape,   "Hex8",   3,  8, 1, 8, Shape_Hex8,  },
    { Hex20Shape,  "Hex20",  3, 20, 2, 8, Shape_Hex20, },
    { Hex27Shape,  "Hex27",  3, 27, 2, 8, Shape_Hex27, },
    { Hex32Shape,  "Hex32",  3, 32, 3, 8, Shape_Hex32, },
    { NoShape,     "NoType", 0,  1, 1, 1, Shape_Point1 }
};


/// Creates an evaluation of the ShapeFunction on the given quadrature scheme.
/// If QuadNodes==Null or QuadWeights==Null the shape functions method is used
/// to generate a quadrature scheme with numQuadNodes nodes. Otherwise it is
/// assumed that a quadrature scheme is given on this array and a copy is
/// created within the structure.
ShapeFunction::ShapeFunction(ShapeFunctionTypeId id, int numQDim,
                             int numQNodes, const std::vector<double>& qNodes,
                             const std::vector<double>& qWeights)
{
    const int numDim=ShapeFunction_InfoList[id].numDim;
    const int numShapes=ShapeFunction_InfoList[id].numShapes;

    if (numQDim>numDim) {
        throw escript::ValueError("ShapeFunction: number of spatial dimensions of quadrature scheme is larger than the spatial dimensionality of shape function.");
    }

    Type=getInfo(id);
    numQuadNodes=numQNodes;

    // allocate memory
    QuadNodes.assign(numQuadNodes*numDim, 0);
    QuadWeights=qWeights;
    S.assign(numShapes*numQuadNodes, 0);
    dSdv.assign(numShapes*numDim*numQuadNodes, 0);

    // set the quadrature nodes (missing values are filled with 0)
    for (int q=0; q<numQuadNodes; q++) {
       for (int i=0; i<numQDim; i++)
           QuadNodes[INDEX2(i,q,numDim)]=qNodes[INDEX2(i,q,numQDim)];
    }

    // evaluate shape functions on quadrature nodes
    Type->getValues(numQuadNodes, QuadNodes, S, dSdv);
}

ShapeFunctionTypeId ShapeFunction::getTypeId(const char* element_type)
{
    int idx=0;
    ShapeFunctionTypeId out=NoShape;
    while (ShapeFunction_InfoList[idx].TypeId!=NoShape && out==NoShape) {
        if (!strcmp(element_type, ShapeFunction_InfoList[idx].Name))
            out=ShapeFunction_InfoList[idx].TypeId;
        idx++;
    }
    return out;
}

const ShapeFunctionInfo* ShapeFunction::getInfo(ShapeFunctionTypeId id)
{
    int idx=0;
    const ShapeFunctionInfo* out=NULL;
    while (ShapeFunction_InfoList[idx].TypeId!=NoShape && out==NULL) {
       if (ShapeFunction_InfoList[idx].TypeId==id)
           out=&ShapeFunction_InfoList[idx];
       idx++;
    }
    if (out==NULL) {
        throw escript::ValueError("ShapeFunction::getInfo: cannot find requested shape function");
    }
    return out;
}


#define V(_K_,_I_) v[INDEX2((_K_)-1,(_I_),DIM)]
#define S(_J_,_I_) s[S_INDEX((_J_)-1,(_I_),NUMSHAPES)]
#define DSDV(_J_,_K_,_I_) dsdv[DSDV_INDEX((_J_)-1,(_K_)-1,(_I_),NUMSHAPES,DIM)]

/****************************************************************************/
void Shape_Point1(int NumV, std::vector<double>& v, std::vector<double>& s,
                  std::vector<double>& dsdv)
{
#define NUMSHAPES 1
#define DIM 0
  for (int i=0; i<NumV; i++)
        S(1,i)=1.;
#undef NUMSHAPES
#undef DIM
}

/****************************************************************************/
void Shape_Line2(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 2
#define DIM 1
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
        S(1,i)=1.-x;
        S(2,i)=   x;
        DSDV(1,1,i)=-1.;
        DSDV(2,1,i)= 1.;
    }
#undef NUMSHAPES
#undef DIM
}

/****************************************************************************/
void Shape_Line3(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 3
#define DIM 1
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
        S(1,i)=(2.*x -1.)*(x-1.);
        S(2,i)=(2.*x -1.)*x;
        S(3,i)= 4.*x*(1.-x);
        DSDV(1,1,i)= 4.*x-3.;
        DSDV(2,1,i)= 4.*x-1.;
        DSDV(3,1,i)=-8.*x+4.;
    }
#undef NUMSHAPES
#undef DIM
}

/****************************************************************************/
void Shape_Line4(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 4
#define DIM 1
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
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

/****************************************************************************/
void Shape_Tri3(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 3
#define DIM 2
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
        const double y=V(2,i);
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

/****************************************************************************/
void Shape_Tri6(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 6
#define DIM 2
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
        const double y=V(2,i);
        S(1,i)=  (1. -x -y)*(1. -2.*x -2.* y);
        S(2,i)=  x*(2.* x -1.);
        S(3,i)=  y*(2.* y -1.);
        S(4,i)=  (1. -x -y)*4.* x;
        S(5,i)=  4.*x*y;
        S(6,i)=  (1. -x -y)*4.* y;
        DSDV(1,1,i)= -3.+4.*x+4.*y;
        DSDV(1,2,i)= -3.+4.*x+4.*y;
        DSDV(2,1,i)= -1.+4.*x;
        DSDV(2,2,i)=  0.;
        DSDV(3,1,i)=  0.;
        DSDV(3,2,i)= -1.          +4.*y;
        DSDV(4,1,i)=     4. -8.*x -4.*y;
        DSDV(4,2,i)=        -4.*x;
        DSDV(5,1,i)=               4.*y;
        DSDV(5,2,i)=         4.*x;
        DSDV(6,1,i)=              -4.*y;
        DSDV(6,2,i)=     4. -4.*x -8.*y;
    }
#undef NUMSHAPES
#undef DIM
}

/****************************************************************************/
void Shape_Tri9(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 9
#define DIM 2
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
        const double y=V(2,i);
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

/****************************************************************************/
void Shape_Tri10(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 10
#define DIM 2
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
        const double y=V(2,i);
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

/****************************************************************************/
void Shape_Rec4(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 4
#define DIM 2
  #pragma ivdep
  for (int i=0; i<NumV; i++) {
    const double x=V(1,i);
    const double y=V(2,i);
    S(1,i)=(1.-x)*(1.-y);
    S(2,i)= x*(1.-y);
    S(3,i)= x*y;
    S(4,i)= (1.-x)*y;
    DSDV(1,1,i)=y-1.;
    DSDV(1,2,i)=x-1.;
    DSDV(2,1,i)= 1.-y;
    DSDV(2,2,i)=-x;
    DSDV(3,1,i)= y;
    DSDV(3,2,i)= x;
    DSDV(4,1,i)=-y;
    DSDV(4,2,i)= 1.-x;
  }
#undef NUMSHAPES
#undef DIM
}

/****************************************************************************/
void Shape_Rec8(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 8
#define DIM 2
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
        const double y=V(2,i);
        S(1,i)= 1.-3.*(x+y)+2.*x*x*(1.-y)+2.*y*y*(1.-x)+5.*x*y;
        S(2,i)= x*(-1.-y+2.*x+2.*y*y-2.*x*y);
        S(3,i)= x*y*(-3.+2.*(x+y));
        S(4,i)= y*(-1.-x+2.*y+2.*x*x-2.*x*y);
        S(5,i)=4.*x*(1.-x-y+x* y);
        S(6,i)= 4.*x*y*(1.-y);
        S(7,i)= 4.*x*y*(1.-x);
        S(8,i)=4.*y*(1.-x-y+x* y);
        DSDV(1,1,i)=-3.+4.*x*(1.-y)+y*(5.-2.*y);
        DSDV(1,2,i)=-3.+4.*y*(1.-x)+x*(5.-2.*x);
        DSDV(2,1,i)=-1.+4.*x*(1.-y)+y*(-1.+2.*y);
        DSDV(2,2,i)= x*(-1.-2.*x+4.*y);
        DSDV(3,1,i)= y*(-3.+4.*x+2.*y);
        DSDV(3,2,i)= x*(-3.+4.*y+2.*x);
        DSDV(4,1,i)= y*(-1.-2.*y+4.*x);
        DSDV(4,2,i)=-1.+4.*y*(1.-x)+x*(-1.+2.*x);
        DSDV(5,1,i)= 4.*(1.-y)+8.*x*(y -1.);
        DSDV(5,2,i)= 4.*x*(x -1.);
        DSDV(6,1,i)= 4.*y*(1.-y);
        DSDV(6,2,i)= 4.*x*(1.-2.*y);
        DSDV(7,1,i)= 4.*y*(1.-2.*x);
        DSDV(7,2,i)= 4.*x*(1.-x);
        DSDV(8,1,i)= 4.*y*(y-1.);
        DSDV(8,2,i)= 4.*(1.-x)+8.*y*(x-1.);
    }
#undef NUMSHAPES
#undef DIM
}

/****************************************************************************/
void Shape_Rec9(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 9
#define DIM 2
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
        const double y=V(2,i);
        S(1,i)= +1. - 3.*x + 2.*x*x - 3.*y + 9.*x*y - 6.*x*x*y + 2.*y*y - 6.*x*y*y + 4.*x*x*y*y;
        S(2,i)= -1.*x + 2.*x*x + 3.*x*y - 6.*x*x*y - 2.*x*y*y + 4.*x*x*y*y;
        S(3,i)= 1.*x*y - 2.*x*x*y - 2.*x*y*y + 4.*x*x*y*y;
        S(4,i)= -1.*y + 3.*x*y - 2.*x*x*y + 2.*y*y - 6.*x*y*y + 4.*x*x*y*y;
        S(5,i)= 4.*x - 4.*x*x - 12.*x*y + 12.*x*x*y + 8.*x*y*y - 8.*x*x*y*y;
        S(6,i)= -4.*x*y + 8.*x*x*y + 4.*x*y*y - 8.*x*x*y*y;
        S(7,i)= -4.*x*y + 4.*x*x*y + 8.*x*y*y - 8.*x*x*y*y;
        S(8,i)= 4.*y - 12.*x*y + 8.*x*x*y - 4.*y*y + 12.*x*y*y - 8.*x*x*y*y;
        S(9,i)= 16.*x*y - 16.*x*x*y - 16.*x*y*y + 16.*x*x*y*y;
        DSDV(1,1,i)= -3. + 4.*x + 9.*y - 12.*x*y - 6.*y*y + 8.*x*y*y;
        DSDV(1,2,i)= -3. + 9.*x - 6.*x*x + 4.*y - 12.*x*y + 8.*x*x*y;
        DSDV(2,1,i)= -1. + 4.*x + 3.*y - 12.*x*y - 2.*y*y + 8.*x*y*y;
        DSDV(2,2,i)= 3.*x - 6.*x*x - 4.*x*y + 8.*x*x*y;
        DSDV(3,1,i)= 1.*y - 4.*x*y - 2.*y*y + 8.*x*y*y;
        DSDV(3,2,i)= 1.*x - 2.*x*x - 4.*x*y + 8.*x*x*y;
        DSDV(4,1,i)= 3.*y - 4.*x*y - 6.*y*y + 8.*x*y*y;
        DSDV(4,2,i)= -1. + 3.*x - 2.*x*x + 4.*y - 12.*x*y + 8.*x*x*y;
        DSDV(5,1,i)= 4. - 8.*x - 12.*y + 24.*x*y + 8.*y*y - 16.*x*y*y;
        DSDV(5,2,i)= -12.*x + 12.*x*x + 16.*x*y - 16.*x*x*y;
        DSDV(6,1,i)= -4.*y + 16.*x*y + 4.*y*y - 16.*x*y*y;
        DSDV(6,2,i)= -4.*x + 8.*x*x + 8.*x*y - 16.*x*x*y;
        DSDV(7,1,i)= -4.*y + 8.*x*y + 8.*y*y - 16.*x*y*y;
        DSDV(7,2,i)= -4.*x + 4.*x*x + 16.*x*y - 16.*x*x*y;
        DSDV(8,1,i)= -12.*y + 16.*x*y + 12.*y*y - 16.*x*y*y;
        DSDV(8,2,i)= 4. - 12.*x + 8.*x*x - 8.*y + 24.*x*y - 16.*x*x*y;
        DSDV(9,1,i)= 16.*y - 32.*x*y - 16.*y*y + 32.*x*y*y;
        DSDV(9,2,i)= 16.*x - 16.*x*x - 32.*x*y + 32.*x*x*y;
    }
#undef NUMSHAPES
#undef DIM
}

/****************************************************************************/

void Shape_Rec12(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 12
#define DIM 2
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
        const double y=V(2,i);
        S(1,i)=10.-5.5*x+(10.)*x*y+(-5.5)*y+(9.)*x*x+(-4.5)*x*x*x+(-9.)*x*y*y+(4.5)*x*y*y*y+(4.5)*x*x*x*y+(-9.)*x*x*y+(-4.5)*y*y*y+(9.)*y*y;
        S(2,i)=10.*x+(-5.5)*x*y+(-4.5)*x*x+(4.5)*x*x*x+(9.)*x*y*y+(-4.5)*x*y*y*y+(-4.5)*x*x*x*y+(4.5)*x*x*y;
        S(3,i)=10.*x*y+(-4.5)*x*y*y+(4.5)*x*y*y*y+(4.5)*x*x*x*y+(-4.5)*x*x*y;
        S(4,i)=-5.5*x*y+10.*y+4.5*x*y*y-4.5*x*y*y*y-4.5*x*x*x*y+9.*x*x*y+4.5*y*y*y-4.5*y*y;
        S(5,i)=9.*x-9.*x*y-22.5*x*x+13.5*x*x*x-13.5*x*x*x*y+22.5*x*x*y;
        S(6,i)=-4.5*x+4.5*x*y+18.*x*x-13.5*x*x*x+13.5*x*x*x*y-18.*x*x*y;
        S(7,i)=9.*x*y+(-22.5)*x*y*y+(13.5)*x*y*y*y;
        S(8,i)=-4.5*x*y+(18.)*x*y*y+(-13.5)*x*y*y*y;
        S(9,i)=-4.5*x*y+(-13.5)*x*x*x*y+(18.)*x*x*y;
        S(10,i)=9.*x*y+(13.5)*x*x*x*y+(-22.5)*x*x*y;
        S(11,i)=4.5*x*y+(-4.5)*y+(-18.)*x*y*y+(13.5)*x*y*y*y+(-13.5)*y*y*y+(18.)*y*y;
        S(12,i)=-9.*x*y+9.*y+22.5*x*y*y-13.5*x*y*y*y+13.5*y*y*y-22.5*y*y;
        DSDV(1,1,i)=-5.5+10.*y+18.*x-13.5*x*x-9.*y*y+4.5*y*y*y+13.5*x*x*y-18.*x*y;
        DSDV(2,1,i)=10.-5.5*y-9.*x+13.5*x*x+9.*y*y-4.5*y*y*y-13.5*x*x*y+9.*x*y;
        DSDV(3,1,i)=10.*y+(-4.5)*y*y+(4.5)*y*y*y+(13.5)*x*x*y+(-9.)*x*y;
        DSDV(4,1,i)=-5.5*y+(4.5)*y*y+(-4.5)*y*y*y+(-13.5)*x*x*y+(18.)*x*y;
        DSDV(5,1,i)=9.-9.*y-45.*x+(0.405e2)*x*x+(-0.405e2)*x*x*y+(45.)*x*y;
        DSDV(6,1,i)=-4.5+4.5*y+36.*x+(-0.405e2)*x*x+(0.405e2)*x*x*y-36.*x*y;
        DSDV(7,1,i)=9.*y+(-22.5)*y*y+(13.5)*y*y*y;
        DSDV(8,1,i)=-4.5*y+(18.)*y*y+(-13.5)*y*y*y;
        DSDV(9,1,i)=-4.5*y+(-0.405e2)*x*x*y+(36.)*x*y;
        DSDV(10,1,i)=9.*y+(0.405e2)*x*x*y+(-45.)*x*y;
        DSDV(11,1,i)=4.5*y+(-18.)*y*y+(13.5)*y*y*y;
        DSDV(12,1,i)=-9.*y+(22.5)*y*y+(-13.5)*y*y*y;
        DSDV(1,2,i)=10*x-5.5-18.*x*y+13.5*x*y*y+4.5*x*x*x-9.*x*x-13.5*y*y+18.*y;
        DSDV(2,2,i)=-5.5*x+(18.)*x*y+(-13.5)*x*y*y+(-4.5)*x*x*x+(4.5)*x*x;
        DSDV(3,2,i)=10.*x+(-9.)*x*y+(13.5)*x*y*y+(4.5)*x*x*x+(-4.5)*x*x;
        DSDV(4,2,i)=-5.5*x+10.+9.*x*y-13.5*x*y*y-4.5*x*x*x+9.*x*x+13.5*y*y-9.*y;
        DSDV(5,2,i)=-9.*x+(-13.5)*x*x*x+(22.5)*x*x;
        DSDV(6,2,i)=4.5*x+(13.5)*x*x*x+(-18.)*x*x;
        DSDV(7,2,i)=9.*x+(-45.)*x*y+(0.405e2)*x*y*y;
        DSDV(8,2,i)=-4.5*x+(36.)*x*y+(-0.405e2)*x*y*y;
        DSDV(9,2,i)=-4.5*x+(-13.5)*x*x*x+(18.)*x*x;
        DSDV(10,2,i)=9.*x+(13.5)*x*x*x+(-22.5)*x*x;
        DSDV(11,2,i)=4.5*x-4.5-36.*x*y+(0.405e2)*x*y*y+(-0.405e2)*y*y+36.*y;
        DSDV(12,2,i)=-9.*x+9.+45.*x*y+(-0.405e2)*x*y*y+(0.405e2)*y*y+(-45.)*y;
    }
#undef NUMSHAPES
#undef DIM
}

/****************************************************************************/
void Shape_Rec16(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 16
#define DIM 2
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
        const double y=V(2,i);
        S(1,i)=(10.)+(-5.5)*x+(0.3025e2)*x*y+(-5.5)*y+(9.)*x*x+(-4.5)*x*x*x+(-0.495e2)*x*y*y+(0.2475e2)*x*y*y*y+(0.2475e2)*x*x*x*y+(-0.495e2)*x*x*y+(-4.5)*y*y*y+(9.)*y*y+(0.81e2)*x*x*y*y+(-0.405e2)*x*x*x*y*y+(0.2025e2)*x*x*x*y*y*y+(-0.405e2)*x*x*y*y*y;
        S(2,i)=(10.)*x+(-5.5)*x*y+(-4.5)*x*x+(4.5)*x*x*x+(9.)*x*y*y+(-4.5)*x*y*y*y+(-0.2475e2)*x*x*x*y+(0.2475e2)*x*x*y+(-0.405e2)*x*x*y*y+(0.405e2)*x*x*x*y*y+(-0.2025e2)*x*x*x*y*y*y+(0.2025e2)*x*x*y*y*y;
        S(3,i)=(10.)*x*y+(-4.5)*x*y*y+(4.5)*x*y*y*y+(4.5)*x*x*x*y+(-4.5)*x*x*y+(0.2025e2)*x*x*y*y+(-0.2025e2)*x*x*x*y*y+(0.2025e2)*x*x*x*y*y*y+(-0.2025e2)*x*x*y*y*y;
        S(4,i)=(-5.5)*x*y+(10.)*y+(0.2475e2)*x*y*y+(-0.2475e2)*x*y*y*y+(-4.5)*x*x*x*y+(9.)*x*x*y+(4.5)*y*y*y+(-4.5)*y*y+(-0.405e2)*x*x*y*y+(0.2025e2)*x*x*x*y*y+(-0.2025e2)*x*x*x*y*y*y+(0.405e2)*x*x*y*y*y;
        S(5,i)=(9.)*x+(-0.495e2)*x*y+(-22.5)*x*x+(13.5)*x*x*x+(0.81e2)*x*y*y+(-0.405e2)*x*y*y*y+(-0.7425e2)*x*x*x*y+(0.12375e3)*x*x*y+(-0.2025e3)*x*x*y*y+(0.1215e3)*x*x*x*y*y+(-0.6075e2)*x*x*x*y*y*y+(0.10125e3)*x*x*y*y*y;
        S(6,i)=(-4.5)*x+(0.2475e2)*x*y+(18.)*x*x+(-13.5)*x*x*x+(-0.405e2)*x*y*y+(0.2025e2)*x*y*y*y+(0.7425e2)*x*x*x*y+(-0.99e2)*x*x*y+(0.162e3)*x*x*y*y+(-0.1215e3)*x*x*x*y*y+(0.6075e2)*x*x*x*y*y*y+(-0.81e2)*x*x*y*y*y;
        S(7,i)=(9.)*x*y+(-22.5)*x*y*y+(13.5)*x*y*y*y+(0.405e2)*x*x*x*y+(-0.405e2)*x*x*y+(0.10125e3)*x*x*y*y+(-0.10125e3)*x*x*x*y*y+(0.6075e2)*x*x*x*y*y*y+(-0.6075e2)*x*x*y*y*y;
        S(8,i)=(-4.5)*x*y+(18.)*x*y*y+(-13.5)*x*y*y*y+(-0.2025e2)*x*x*x*y+(0.2025e2)*x*x*y+(-0.81e2)*x*x*y*y+(0.81e2)*x*x*x*y*y+(-0.6075e2)*x*x*x*y*y*y+(0.6075e2)*x*x*y*y*y;
        S(9,i)=(-4.5)*x*y+(0.2025e2)*x*y*y+(-0.2025e2)*x*y*y*y+(-13.5)*x*x*x*y+(18.)*x*x*y+(-0.81e2)*x*x*y*y+(0.6075e2)*x*x*x*y*y+(-0.6075e2)*x*x*x*y*y*y+(0.81e2)*x*x*y*y*y;
        S(10,i)=(9.)*x*y+(-0.405e2)*x*y*y+(0.405e2)*x*y*y*y+(13.5)*x*x*x*y+(-22.5)*x*x*y+(0.10125e3)*x*x*y*y+(-0.6075e2)*x*x*x*y*y+(0.6075e2)*x*x*x*y*y*y+(-0.10125e3)*x*x*y*y*y;
        S(11,i)=(0.2475e2)*x*y+(-4.5)*y+(-0.99e2)*x*y*y+(0.7425e2)*x*y*y*y+(0.2025e2)*x*x*x*y+(-0.405e2)*x*x*y+(-13.5)*y*y*y+(18.)*y*y+(0.162e3)*x*x*y*y+(-0.81e2)*x*x*x*y*y+(0.6075e2)*x*x*x*y*y*y+(-0.1215e3)*x*x*y*y*y;
        S(12,i)=(-0.495e2)*x*y+(9.)*y+(0.12375e3)*x*y*y+(-0.7425e2)*x*y*y*y+(-0.405e2)*x*x*x*y+(0.81e2)*x*x*y+(13.5)*y*y*y+(-22.5)*y*y+(-0.2025e3)*x*x*y*y+(0.10125e3)*x*x*x*y*y+(-0.6075e2)*x*x*x*y*y*y+(0.1215e3)*x*x*y*y*y;
        S(13,i)=(0.81e2)*x*y+(-0.2025e3)*x*y*y+(0.1215e3)*x*y*y*y+(0.1215e3)*x*x*x*y+(-0.2025e3)*x*x*y+(0.50625e3)*x*x*y*y+(-0.30375e3)*x*x*x*y*y+(0.18225e3)*x*x*x*y*y*y+(-0.30375e3)*x*x*y*y*y;
        S(14,i)=(-0.405e2)*x*y+(0.10125e3)*x*y*y+(-0.6075e2)*x*y*y*y+(-0.1215e3)*x*x*x*y+(0.162e3)*x*x*y+(-0.405e3)*x*x*y*y+(0.30375e3)*x*x*x*y*y+(-0.18225e3)*x*x*x*y*y*y+(0.243e3)*x*x*y*y*y;
        S(15,i)=(0.2025e2)*x*y+(-0.81e2)*x*y*y+(0.6075e2)*x*y*y*y+(0.6075e2)*x*x*x*y+(-0.81e2)*x*x*y+(0.324e3)*x*x*y*y+(-0.243e3)*x*x*x*y*y+(0.18225e3)*x*x*x*y*y*y+(-0.243e3)*x*x*y*y*y;
        S(16,i)=(-0.405e2)*x*y+(0.162e3)*x*y*y+(-0.1215e3)*x*y*y*y+(-0.6075e2)*x*x*x*y+(0.10125e3)*x*x*y+(-0.405e3)*x*x*y*y+(0.243e3)*x*x*x*y*y+(-0.18225e3)*x*x*x*y*y*y+(0.30375e3)*x*x*y*y*y;
        DSDV(1, 1,i)=(-5.5)+(0.3025e2)*y+(18.)*x+(-13.5)*x*x+(-0.495e2)*y*y+(0.2475e2)*y*y*y+(0.7425e2)*x*x*y+(-0.99e2)*x*y+(0.162e3)*x*y*y+(-0.1215e3)*x*x*y*y+(0.6075e2)*x*x*y*y*y+(-0.81e2)*x*y*y*y;
        DSDV(2, 1,i)=(10.)+(-5.5)*y+(-9.)*x+(13.5)*x*x+(9.)*y*y+(-4.5)*y*y*y+(-0.7425e2)*x*x*y+(0.495e2)*x*y+(-0.81e2)*x*y*y+(0.1215e3)*x*x*y*y+(-0.6075e2)*x*x*y*y*y+(0.405e2)*x*y*y*y;
        DSDV(3, 1,i)=(10.)*y+(-4.5)*y*y+(4.5)*y*y*y+(13.5)*x*x*y+(-9.)*x*y+(0.405e2)*x*y*y+(-0.6075e2)*x*x*y*y+(0.6075e2)*x*x*y*y*y+(-0.405e2)*x*y*y*y;
        DSDV(4, 1,i)=(-5.5)*y+(0.2475e2)*y*y+(-0.2475e2)*y*y*y+(-13.5)*x*x*y+(18.)*x*y+(-0.81e2)*x*y*y+(0.6075e2)*x*x*y*y+(-0.6075e2)*x*x*y*y*y+(0.81e2)*x*y*y*y;
        DSDV(5, 1,i)=(9.)+(-0.495e2)*y+(-45.)*x+(0.405e2)*x*x+(0.81e2)*y*y+(-0.405e2)*y*y*y+(-0.22275e3)*x*x*y+(0.2475e3)*x*y+(-0.405e3)*x*y*y+(0.3645e3)*x*x*y*y+(-0.18225e3)*x*x*y*y*y+(0.2025e3)*x*y*y*y;
        DSDV(6, 1,i)=(-4.5)+(0.2475e2)*y+(36.)*x+(-0.405e2)*x*x+(-0.405e2)*y*y+(0.2025e2)*y*y*y+(0.22275e3)*x*x*y+(-0.198e3)*x*y+(0.324e3)*x*y*y+(-0.3645e3)*x*x*y*y+(0.18225e3)*x*x*y*y*y+(-0.162e3)*x*y*y*y;
        DSDV(7, 1,i)=(9.)*y+(-22.5)*y*y+(13.5)*y*y*y+(0.1215e3)*x*x*y+(-0.81e2)*x*y+(0.2025e3)*x*y*y+(-0.30375e3)*x*x*y*y+(0.18225e3)*x*x*y*y*y+(-0.1215e3)*x*y*y*y;
        DSDV(8, 1,i)=(-4.5)*y+(18.)*y*y+(-13.5)*y*y*y+(-0.6075e2)*x*x*y+(0.405e2)*x*y+(-0.162e3)*x*y*y+(0.243e3)*x*x*y*y+(-0.18225e3)*x*x*y*y*y+(0.1215e3)*x*y*y*y;
        DSDV(9, 1,i)=(-4.5)*y+(0.2025e2)*y*y+(-0.2025e2)*y*y*y+(-0.405e2)*x*x*y+(36.)*x*y+(-0.162e3)*x*y*y+(0.18225e3)*x*x*y*y+(-0.18225e3)*x*x*y*y*y+(0.162e3)*x*y*y*y;
        DSDV(10, 1,i)=(9.)*y+(-0.405e2)*y*y+(0.405e2)*y*y*y+(0.405e2)*x*x*y+(-45.)*x*y+(0.2025e3)*x*y*y+(-0.18225e3)*x*x*y*y+(0.18225e3)*x*x*y*y*y+(-0.2025e3)*x*y*y*y;
        DSDV(11, 1,i)=(0.2475e2)*y+(-0.99e2)*y*y+(0.7425e2)*y*y*y+(0.6075e2)*x*x*y+(-0.81e2)*x*y+(0.324e3)*x*y*y+(-0.243e3)*x*x*y*y+(0.18225e3)*x*x*y*y*y+(-0.243e3)*x*y*y*y;
        DSDV(12, 1,i)=(-0.495e2)*y+(0.12375e3)*y*y+(-0.7425e2)*y*y*y+(-0.1215e3)*x*x*y+(0.162e3)*x*y+(-0.405e3)*x*y*y+(0.30375e3)*x*x*y*y+(-0.18225e3)*x*x*y*y*y+(0.243e3)*x*y*y*y;
        DSDV(13, 1,i)=(0.81e2)*y+(-0.2025e3)*y*y+(0.1215e3)*y*y*y+(0.3645e3)*x*x*y+(-0.405e3)*x*y+(0.10125e4)*x*y*y+(-0.91125e3)*x*x*y*y+(0.54675e3)*x*x*y*y*y+(-0.6075e3)*x*y*y*y;
        DSDV(14, 1,i)=(-0.405e2)*y+(0.10125e3)*y*y+(-0.6075e2)*y*y*y+(-0.3645e3)*x*x*y+(0.324e3)*x*y+(-0.81e3)*x*y*y+(0.91125e3)*x*x*y*y+(-0.54675e3)*x*x*y*y*y+(0.486e3)*x*y*y*y;
        DSDV(15, 1,i)=(0.2025e2)*y+(-0.81e2)*y*y+(0.6075e2)*y*y*y+(0.18225e3)*x*x*y+(-0.162e3)*x*y+(0.648e3)*x*y*y+(-0.729e3)*x*x*y*y+(0.54675e3)*x*x*y*y*y+(-0.486e3)*x*y*y*y;
        DSDV(16, 1,i)=(-0.405e2)*y+(0.162e3)*y*y+(-0.1215e3)*y*y*y+(-0.18225e3)*x*x*y+(0.2025e3)*x*y+(-0.81e3)*x*y*y+(0.729e3)*x*x*y*y+(-0.54675e3)*x*x*y*y*y+(0.6075e3)*x*y*y*y;
        DSDV(1, 2,i)=(0.3025e2)*x+(-5.5)+(-0.99e2)*x*y+(0.7425e2)*x*y*y+(0.2475e2)*x*x*x+(-0.495e2)*x*x+(-13.5)*y*y+(18.)*y+(0.162e3)*x*x*y+(-0.81e2)*x*x*x*y+(0.6075e2)*x*x*x*y*y+(-0.1215e3)*x*x*y*y;
        DSDV(2, 2,i)=(-5.5)*x+(18.)*x*y+(-13.5)*x*y*y+(-0.2475e2)*x*x*x+(0.2475e2)*x*x+(-0.81e2)*x*x*y+(0.81e2)*x*x*x*y+(-0.6075e2)*x*x*x*y*y+(0.6075e2)*x*x*y*y;
        DSDV(3, 2,i)=(10.)*x+(-9.)*x*y+(13.5)*x*y*y+(4.5)*x*x*x+(-4.5)*x*x+(0.405e2)*x*x*y+(-0.405e2)*x*x*x*y+(0.6075e2)*x*x*x*y*y+(-0.6075e2)*x*x*y*y;
        DSDV(4, 2,i)=(-5.5)*x+(10.)+(0.495e2)*x*y+(-0.7425e2)*x*y*y+(-4.5)*x*x*x+(9.)*x*x+(13.5)*y*y+(-9.)*y+(-0.81e2)*x*x*y+(0.405e2)*x*x*x*y+(-0.6075e2)*x*x*x*y*y+(0.1215e3)*x*x*y*y;
        DSDV(5, 2,i)=(-0.495e2)*x+(0.162e3)*x*y+(-0.1215e3)*x*y*y+(-0.7425e2)*x*x*x+(0.12375e3)*x*x+(-0.405e3)*x*x*y+(0.243e3)*x*x*x*y+(-0.18225e3)*x*x*x*y*y+(0.30375e3)*x*x*y*y;
        DSDV(6, 2,i)=(0.2475e2)*x+(-0.81e2)*x*y+(0.6075e2)*x*y*y+(0.7425e2)*x*x*x+(-0.99e2)*x*x+(0.324e3)*x*x*y+(-0.243e3)*x*x*x*y+(0.18225e3)*x*x*x*y*y+(-0.243e3)*x*x*y*y;
        DSDV(7, 2,i)=(9.)*x+(-45.)*x*y+(0.405e2)*x*y*y+(0.405e2)*x*x*x+(-.405e2)*x*x+(0.2025e3)*x*x*y+(-0.2025e3)*x*x*x*y+(0.18225e3)*x*x*x*y*y+(-0.18225e3)*x*x*y*y;
        DSDV(8, 2,i)=(-4.5)*x+(36.)*x*y+(-0.405e2)*x*y*y+(-0.2025e2)*x*x*x+(0.2025e2)*x*x+(-0.162e3)*x*x*y+(0.162e3)*x*x*x*y+(-0.18225e3)*x*x*x*y*y+(0.18225e3)*x*x*y*y;
        DSDV(9, 2,i)=(-4.5)*x+(0.405e2)*x*y+(-0.6075e2)*x*y*y+(-13.5)*x*x*x+(18.)*x*x+(-0.162e3)*x*x*y+(0.1215e3)*x*x*x*y+(-0.18225e3)*x*x*x*y*y+(0.243e3)*x*x*y*y;
        DSDV(10, 2,i)=(9.)*x+(-0.81e2)*x*y+(0.1215e3)*x*y*y+(13.5)*x*x*x+(-22.5)*x*x+(0.2025e3)*x*x*y+(-0.1215e3)*x*x*x*y+(0.18225e3)*x*x*x*y*y+(-0.30375e3)*x*x*y*y;
        DSDV(11, 2,i)=(0.2475e2)*x+(-4.5)+(-0.198e3)*x*y+(0.22275e3)*x*y*y+(0.2025e2)*x*x*x+(-0.405e2)*x*x+(-0.405e2)*y*y+(36.)*y+(0.324e3)*x*x*y+(-0.162e3)*x*x*x*y+(0.18225e3)*x*x*x*y*y+(-0.3645e3)*x*x*y*y;
        DSDV(12, 2,i)=(-0.495e2)*x+(9.)+(0.2475e3)*x*y+(-0.22275e3)*x*y*y+(-0.405e2)*x*x*x+(0.81e2)*x*x+(0.405e2)*y*y+(-45.)*y+(-0.405e3)*x*x*y+(0.2025e3)*x*x*x*y+(-0.18225e3)*x*x*x*y*y+(0.3645e3)*x*x*y*y;
        DSDV(13, 2,i)=(0.81e2)*x+(-0.405e3)*x*y+(0.3645e3)*x*y*y+(0.1215e3)*x*x*x+(-0.2025e3)*x*x+(0.10125e4)*x*x*y+(-0.6075e3)*x*x*x*y+(0.54675e3)*x*x*x*y*y+(-0.91125e3)*x*x*y*y;
        DSDV(14, 2,i)=(-0.405e2)*x+(0.2025e3)*x*y+(-0.18225e3)*x*y*y+(-0.1215e3)*x*x*x+(0.162e3)*x*x+(-0.81e3)*x*x*y+(0.6075e3)*x*x*x*y+(-0.54675e3)*x*x*x*y*y+(0.729e3)*x*x*y*y;
        DSDV(15, 2,i)=(0.2025e2)*x+(-0.162e3)*x*y+(0.18225e3)*x*y*y+(0.6075e2)*x*x*x+(-0.81e2)*x*x+(0.648e3)*x*x*y+(-0.486e3)*x*x*x*y+(0.54675e3)*x*x*x*y*y+(-0.729e3)*x*x*y*y;
        DSDV(16, 2,i)=(-0.405e2)*x+(0.324e3)*x*y+(-0.3645e3)*x*y*y+(-0.6075e2)*x*x*x+(0.10125e3)*x*x+(-0.81e3)*x*x*y+(0.486e3)*x*x*x*y+(-0.54675e3)*x*x*x*y*y+(0.91125e3)*x*x*y*y;
    }
#undef NUMSHAPES
#undef DIM
}

/****************************************************************************/
void Shape_Tet4(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 4
#define DIM 3
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
        const double y=V(2,i);
        const double z=V(3,i);
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

/****************************************************************************/
void Shape_Tet10(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 10
#define DIM 3
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
        const double y=V(2,i);
        const double z=V(3,i);
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
        DSDV(2,2,i)= 0.;
        DSDV(2,3,i)= 0.;

        DSDV(3,1,i)= 0.;
        DSDV(3,2,i)= -1.      +4.*y;
        DSDV(3,3,i)= 0.;

        DSDV(4,1,i)= 0.;
        DSDV(4,2,i)= 0.;
        DSDV(4,3,i)= -1.            +4.*z;

        DSDV(5,1,i)= 4. -8.*x -4.*y -4.*z;
        DSDV(5,2,i)=    -4.*x;
        DSDV(5,3,i)=    -4.*x;

        DSDV(6,1,i)=           4.*y;
        DSDV(6,2,i)=     4.*x;
        DSDV(6,3,i)= 0.;

        DSDV(7,1,i)=          -4.*y;
        DSDV(7,2,i)= 4. -4.*x -8.*y -4.*z;
        DSDV(7,3,i)=          -4.*y;

        DSDV(8,1,i)=                -4.*z;
        DSDV(8,2,i)=                -4.*z;
        DSDV(8,3,i)= 4. -4.*x -4.*y -8.*z;

        DSDV(9,1,i)=                 4.*z;
        DSDV(9,2,i)= 0.;
        DSDV(9,3,i)=     4.*x;

        DSDV(10,1,i)=0.;
        DSDV(10,2,i)=                4.*z;
        DSDV(10,3,i)=          4.*y;
    }
#undef NUMSHAPES
#undef DIM
}

/****************************************************************************/
void Shape_Tet16(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 16
#define DIM 3
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
        const double y=V(2,i);
        const double z=V(3,i);
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

/****************************************************************************/
void Shape_Hex8(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 8
#define DIM 3
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
        const double y=V(2,i);
        const double z=V(3,i);
        S(1,i)=(1.-x)*(1.-y)*(1.-z);
        S(2,i)= x*(1.-z)*(1.-y);
        S(3,i)= x*(1.-z)*y;
        S(4,i)= (1.-z)*(1.-x)*y;
        S(5,i)= (1.-x)*z*(1.-y);
        S(6,i)= x*z*(1.-y);
        S(7,i)= x*y*z;
        S(8,i)= y*z*(1.-x);
        DSDV(1,1,i)= (1.-z)*(y-1.);
        DSDV(1,2,i)= (1.-x)*(z-1.);
        DSDV(1,3,i)= (1.-x)*(y-1.);
        DSDV(2,1,i)= (1.-z)*(1.-y);
        DSDV(2,2,i)= (z-1.)*x;
        DSDV(2,3,i)= (y-1.)*x;
        DSDV(3,1,i)= (1.-z)*y;
        DSDV(3,2,i)= (1.-z)*x;
        DSDV(3,3,i)=-y*x;
        DSDV(4,1,i)= y*(z-1.);
        DSDV(4,2,i)= (1.-z)*(1.-x);
        DSDV(4,3,i)= y*(x-1.);
        DSDV(5,1,i)= z*(y-1.);
        DSDV(5,2,i)= z*(x-1.);
        DSDV(5,3,i)= (x-1.)*(y-1.);
        DSDV(6,1,i)= z*(1.-y);
        DSDV(6,2,i)= -x*z;
        DSDV(6,3,i)= (1.-y)*x;
        DSDV(7,1,i)= y*z;
        DSDV(7,2,i)= x*z;
        DSDV(7,3,i)= x*y;
        DSDV(8,1,i)=-y*z;
        DSDV(8,2,i)= z*(1.-x);
        DSDV(8,3,i)= y*(1.-x);
    }
#undef NUMSHAPES
#undef DIM
}

/****************************************************************************/
void Shape_Hex20(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 20
#define DIM 3
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
        const double y=V(2,i);
        const double z=V(3,i);
        S(1,i)=1.+(-3.)*x+(-3.)*y+(-3.)*z+(5.)*x*y+(5.)*x*z+(5.)*y*z+(2.)*x*x+(2.)*y*y+(2.)*z*z+(-2.)*x*x*y+(-2.)*x*x*z+(-2.)*x*y*y+(-2.)*y*y*z+(-2.)*x*z*z+(-2.)*y*z*z+(-7.)*x*y*z+(2.)*x*x*y*z+(2.)*x*y*y*z+(2.)*x*y*z*z;
        S(2,i)=(-1.)*x+(-1.)*x*y+(-1.)*x*z+(2.)*x*x+(-2.)*x*x*y+(-2.)*x*x*z+(2.)*x*y*y+(2.)*x*z*z+(3.)*x*y*z+(2.)*x*x*y*z+(-2.)*x*y*y*z+(-2.)*x*y*z*z;
        S(3,i)=(-3.)*x*y+(2.)*x*x*y+(2.)*x*y*y+1.*x*y*z+(-2.)*x*x*y*z+(-2.)*x*y*y*z+(2.)*x*y*z*z;
        S(4,i)=(-1.)*y+(-1.)*x*y+(-1.)*y*z+(2.)*y*y+(2.)*x*x*y+(-2.)*x*y*y+(-2.)*y*y*z+(2.)*y*z*z+(3.)*x*y*z+(-2.)*x*x*y*z+(2.)*x*y*y*z+(-2.)*x*y*z*z;
        S(5,i)=(-1.)*z+(-1.)*x*z+(-1.)*y*z+(2.)*z*z+(2.)*x*x*z+(2.)*y*y*z+(-2.)*x*z*z+(-2.)*y*z*z+(3.)*x*y*z+(-2.)*x*x*y*z+(-2.)*x*y*y*z+(2.)*x*y*z*z;
        S(6,i)=(-3.)*x*z+(2.)*x*x*z+(2.)*x*z*z+1.*x*y*z+(-2.)*x*x*y*z+(2.)*x*y*y*z+(-2.)*x*y*z*z;
        S(7,i)=(-5.)*x*y*z+(2.)*x*x*y*z+(2.)*x*y*y*z+(2.)*x*y*z*z;
        S(8,i)=(-3.)*y*z+(2.)*y*y*z+(2.)*y*z*z+1.*x*y*z+(2.)*x*x*y*z+(-2.)*x*y*y*z+(-2.)*x*y*z*z;
        S(9,i)=(4.)*x+(-4.)*x*y+(-4.)*x*z+(-4.)*x*x+(4.)*x*x*y+(4.)*x*x*z+(4.)*x*y*z+(-4.)*x*x*y*z;
        S(10,i)=(4.)*x*y+(-4.)*x*y*y+(-4.)*x*y*z+(4.)*x*y*y*z;
        S(11,i)=(4.)*x*y+(-4.)*x*x*y+(-4.)*x*y*z+(4.)*x*x*y*z;
        S(12,i)=(4.)*y+(-4.)*x*y+(-4.)*y*z+(-4.)*y*y+(4.)*x*y*y+(4.)*y*y*z+(4.)*x*y*z+(-4.)*x*y*y*z;
        S(13,i)=(4.)*z+(-4.)*x*z+(-4.)*y*z+(-4.)*z*z+(4.)*x*z*z+(4.)*y*z*z+(4.)*x*y*z+(-4.)*x*y*z*z;
        S(14,i)=(4.)*x*z+(-4.)*x*z*z+(-4.)*x*y*z+(4.)*x*y*z*z;
        S(15,i)=(4.)*x*y*z+(-4.)*x*y*z*z;
        S(16,i)=(4.)*y*z+(-4.)*y*z*z+(-4.)*x*y*z+(4.)*x*y*z*z;
        S(17,i)=(4.)*x*z+(-4.)*x*x*z+(-4.)*x*y*z+(4.)*x*x*y*z;
        S(18,i)=(4.)*x*y*z+(-4.)*x*y*y*z;
        S(19,i)=(4.)*x*y*z+(-4.)*x*x*y*z;
        S(20,i)=(4.)*y*z+(-4.)*y*y*z+(-4.)*x*y*z+(4.)*x*y*y*z;
        DSDV(1,1,i)=(-3.)+(5.)*y+(5.)*z+(4.)*x+(-4.)*x*y+(-4.)*x*z+(-2.)*y*y+(-2.)*z*z+(-7.)*y*z+(4.)*x*y*z+(2.)*y*y*z+(2.)*y*z*z;
        DSDV(2,1,i)=(-1.)+(-1.)*y+(-1.)*z+(4.)*x+(-4.)*x*y+(-4.)*x*z+(2.)*y*y+(2.)*z*z+(3.)*y*z+(4.)*x*y*z+(-2.)*y*y*z+(-2.)*y*z*z;
        DSDV(3,1,i)=(-3.)*y+(4.)*x*y+(2.)*y*y+1.*y*z+(-4.)*x*y*z+(-2.)*y*y*z+(2.)*y*z*z;
        DSDV(4,1,i)=(-1.)*y+(4.)*x*y+(-2.)*y*y+(3.)*y*z+(-4.)*x*y*z+(2.)*y*y*z+(-2.)*y*z*z;
        DSDV(5,1,i)=(-1.)*z+(4.)*x*z+(-2.)*z*z+(3.)*y*z+(-4.)*x*y*z+(-2.)*y*y*z+(2.)*y*z*z;
        DSDV(6,1,i)=(-3.)*z+(4.)*x*z+(2.)*z*z+1.*y*z+(-4.)*x*y*z+(2.)*y*y*z+(-2.)*y*z*z;
        DSDV(7,1,i)=(-5.)*y*z+(4.)*x*y*z+(2.)*y*y*z+(2.)*y*z*z;
        DSDV(8,1,i)=1.*y*z+(4.)*x*y*z+(-2.)*y*y*z+(-2.)*y*z*z;
        DSDV(9,1,i)=(4.)+(-4.)*y+(-4.)*z+(-8.)*x+(8.)*x*y+(8.)*x*z+(4.)*y*z+(-8.)*x*y*z;
        DSDV(10,1,i)=(4.)*y+(-4.)*y*y+(-4.)*y*z+(4.)*y*y*z;
        DSDV(11,1,i)=(4.)*y+(-8.)*x*y+(-4.)*y*z+(8.)*x*y*z;
        DSDV(12,1,i)=(-4.)*y+(4.)*y*y+(4.)*y*z+(-4.)*y*y*z;
        DSDV(13,1,i)=(-4.)*z+(4.)*z*z+(4.)*y*z+(-4.)*y*z*z;
        DSDV(14,1,i)=(4.)*z+(-4.)*z*z+(-4.)*y*z+(4.)*y*z*z;
        DSDV(15,1,i)=(4.)*y*z+(-4.)*y*z*z;
        DSDV(16,1,i)=(-4.)*y*z+(4.)*y*z*z;
        DSDV(17,1,i)=(4.)*z+(-8.)*x*z+(-4.)*y*z+(8.)*x*y*z;
        DSDV(18,1,i)=(4.)*y*z+(-4.)*y*y*z;
        DSDV(19,1,i)=(4.)*y*z+(-8.)*x*y*z;
        DSDV(20,1,i)=(-4.)*y*z+(4.)*y*y*z;
        DSDV(1,2,i)=(-3.)+(5.)*x+(5.)*z+(4.)*y+(-2.)*x*x+(-4.)*x*y+(-4.)*y*z+(-2.)*z*z+(-7.)*x*z+(2.)*x*x*z+(4.)*x*y*z+(2.)*x*z*z;
        DSDV(2,2,i)=(-1.)*x+(-2.)*x*x+(4.)*x*y+(3.)*x*z+(2.)*x*x*z+(-4.)*x*y*z+(-2.)*x*z*z;
        DSDV(3,2,i)=(-3.)*x+(2.)*x*x+(4.)*x*y+1.*x*z+(-2.)*x*x*z+(-4.)*x*y*z+(2.)*x*z*z;
        DSDV(4,2,i)=(-1.)+(-1.)*x+(-1.)*z+(4.)*y+(2.)*x*x+(-4.)*x*y+(-4.)*y*z+(2.)*z*z+(3.)*x*z+(-2.)*x*x*z+(4.)*x*y*z+(-2.)*x*z*z;
        DSDV(5,2,i)=(-1.)*z+(4.)*y*z+(-2.)*z*z+(3.)*x*z+(-2.)*x*x*z+(-4.)*x*y*z+(2.)*x*z*z;
        DSDV(6,2,i)=1.*x*z+(-2.)*x*x*z+(4.)*x*y*z+(-2.)*x*z*z;
        DSDV(7,2,i)=(-5.)*x*z+(2.)*x*x*z+(4.)*x*y*z+(2.)*x*z*z;
        DSDV(8,2,i)=(-3.)*z+(4.)*y*z+(2.)*z*z+1.*x*z+(2.)*x*x*z+(-4.)*x*y*z+(-2.)*x*z*z;
        DSDV(9,2,i)=(-4.)*x+(4.)*x*x+(4.)*x*z+(-4.)*x*x*z;
        DSDV(10,2,i)=(4.)*x+(-8.)*x*y+(-4.)*x*z+(8.)*x*y*z;
        DSDV(11,2,i)=(4.)*x+(-4.)*x*x+(-4.)*x*z+(4.)*x*x*z;
        DSDV(12,2,i)=(4.)+(-4.)*x+(-4.)*z+(-8.)*y+(8.)*x*y+(8.)*y*z+(4.)*x*z+(-8.)*x*y*z;
        DSDV(13,2,i)=(-4.)*z+(4.)*z*z+(4.)*x*z+(-4.)*x*z*z;
        DSDV(14,2,i)=(-4.)*x*z+(4.)*x*z*z;
        DSDV(15,2,i)=(4.)*x*z+(-4.)*x*z*z;
        DSDV(16,2,i)=(4.)*z+(-4.)*z*z+(-4.)*x*z+(4.)*x*z*z;
        DSDV(17,2,i)=(-4.)*x*z+(4.)*x*x*z;
        DSDV(18,2,i)=(4.)*x*z+(-8.)*x*y*z;
        DSDV(19,2,i)=(4.)*x*z+(-4.)*x*x*z;
        DSDV(20,2,i)=(4.)*z+(-8.)*y*z+(-4.)*x*z+(8.)*x*y*z;
        DSDV(1,3,i)=(-3.)+(5.)*x+(5.)*y+(4.)*z+(-2.)*x*x+(-2.)*y*y+(-4.)*x*z+(-4.)*y*z+(-7.)*x*y+(2.)*x*x*y+(2.)*x*y*y+(4.)*x*y*z;
        DSDV(2,3,i)=(-1.)*x+(-2.)*x*x+(4.)*x*z+(3.)*x*y+(2.)*x*x*y+(-2.)*x*y*y+(-4.)*x*y*z;
        DSDV(3,3,i)=1.*x*y+(-2.)*x*x*y+(-2.)*x*y*y+(4.)*x*y*z;
        DSDV(4,3,i)=(-1.)*y+(-2.)*y*y+(4.)*y*z+(3.)*x*y+(-2.)*x*x*y+(2.)*x*y*y+(-4.)*x*y*z;
        DSDV(5,3,i)=(-1.)+(-1.)*x+(-1.)*y+(4.)*z+(2.)*x*x+(2.)*y*y+(-4.)*x*z+(-4.)*y*z+(3.)*x*y+(-2.)*x*x*y+(-2.)*x*y*y+(4.)*x*y*z;
        DSDV(6,3,i)=(-3.)*x+(2.)*x*x+(4.)*x*z+1.*x*y+(-2.)*x*x*y+(2.)*x*y*y+(-4.)*x*y*z;
        DSDV(7,3,i)=(-5.)*x*y+(2.)*x*x*y+(2.)*x*y*y+(4.)*x*y*z;
        DSDV(8,3,i)=(-3.)*y+(2.)*y*y+(4.)*y*z+1.*x*y+(2.)*x*x*y+(-2.)*x*y*y+(-4.)*x*y*z;
        DSDV(9,3,i)=(-4.)*x+(4.)*x*x+(4.)*x*y+(-4.)*x*x*y;
        DSDV(10,3,i)=(-4.)*x*y+(4.)*x*y*y;
        DSDV(11,3,i)=(-4.)*x*y+(4.)*x*x*y;
        DSDV(12,3,i)=(-4.)*y+(4.)*y*y+(4.)*x*y+(-4.)*x*y*y;
        DSDV(13,3,i)=(4.)+(-4.)*x+(-4.)*y+(-8.)*z+(8.)*x*z+(8.)*y*z+(4.)*x*y+(-8.)*x*y*z;
        DSDV(14,3,i)=(4.)*x+(-8.)*x*z+(-4.)*x*y+(8.)*x*y*z;
        DSDV(15,3,i)=(4.)*x*y+(-8.)*x*y*z;
        DSDV(16,3,i)=(4.)*y+(-8.)*y*z+(-4.)*x*y+(8.)*x*y*z;
        DSDV(17,3,i)=(4.)*x+(-4.)*x*x+(-4.)*x*y+(4.)*x*x*y;
        DSDV(18,3,i)=(4.)*x*y+(-4.)*x*y*y;
        DSDV(19,3,i)=(4.)*x*y+(-4.)*x*x*y;
        DSDV(20,3,i)=(4.)*y+(-4.)*y*y+(-4.)*x*y+(4.)*x*y*y;
    }
#undef NUMSHAPES
#undef DIM
}

/****************************************************************************/
void Shape_Hex27(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 27
#define DIM 3
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
        const double y=V(2,i);
        const double z=V(3,i);
        S(1,i)= + 1.0 - 3.0*x + 2.0*x*x - 3.0*y + 9.0*x*y - 6.0*x*x*y + 2.0*y*y - 6.0*x*y*y + 4.0*x*x*y*y - 3.0*z + 9.0*x*z - 6.0*x*x*z + 9.0*y*z - 27.0*x*y*z + 18.0*x*x*y*z - 6.0*y*y*z + 18.0*x*y*y*z - 12.0*x*x*y*y*z + 2.0*z*z - 6.0*x*z*z + 4.0*x*x*z*z - 6.0*y*z*z + 18.0*x*y*z*z - 12.0*x*x*y*z*z + 4.0*y*y*z*z - 12.0*x*y*y*z*z + 8.0*x*x*y*y*z*z;
        S(2,i)= - 1.0*x + 2.0*x*x + 3.0*x*y - 6.0*x*x*y - 2.0*x*y*y + 4.0*x*x*y*y + 3.0*x*z - 6.0*x*x*z - 9.0*x*y*z + 18.0*x*x*y*z + 6.0*x*y*y*z - 12.0*x*x*y*y*z - 2.0*x*z*z + 4.0*x*x*z*z + 6.0*x*y*z*z - 12.0*x*x*y*z*z - 4.0*x*y*y*z*z + 8.0*x*x*y*y*z*z;
        S(3,i)= + 1.0*x*y - 2.0*x*x*y - 2.0*x*y*y + 4.0*x*x*y*y - 3.0*x*y*z + 6.0*x*x*y*z + 6.0*x*y*y*z - 12.0*x*x*y*y*z + 2.0*x*y*z*z - 4.0*x*x*y*z*z - 4.0*x*y*y*z*z + 8.0*x*x*y*y*z*z;
        S(4,i)= - 1.0*y + 3.0*x*y - 2.0*x*x*y + 2.0*y*y - 6.0*x*y*y + 4.0*x*x*y*y + 3.0*y*z - 9.0*x*y*z + 6.0*x*x*y*z - 6.0*y*y*z + 18.0*x*y*y*z - 12.0*x*x*y*y*z - 2.0*y*z*z + 6.0*x*y*z*z - 4.0*x*x*y*z*z + 4.0*y*y*z*z - 12.0*x*y*y*z*z + 8.0*x*x*y*y*z*z;
        S(5,i)= - 1.0*z + 3.0*x*z - 2.0*x*x*z + 3.0*y*z - 9.0*x*y*z + 6.0*x*x*y*z - 2.0*y*y*z + 6.0*x*y*y*z - 4.0*x*x*y*y*z + 2.0*z*z - 6.0*x*z*z + 4.0*x*x*z*z - 6.0*y*z*z + 18.0*x*y*z*z - 12.0*x*x*y*z*z + 4.0*y*y*z*z - 12.0*x*y*y*z*z + 8.0*x*x*y*y*z*z;
        S(6,i)= + 1.0*x*z - 2.0*x*x*z - 3.0*x*y*z + 6.0*x*x*y*z + 2.0*x*y*y*z - 4.0*x*x*y*y*z - 2.0*x*z*z + 4.0*x*x*z*z + 6.0*x*y*z*z - 12.0*x*x*y*z*z - 4.0*x*y*y*z*z + 8.0*x*x*y*y*z*z;
        S(7,i)= - 1.0*x*y*z + 2.0*x*x*y*z + 2.0*x*y*y*z - 4.0*x*x*y*y*z + 2.0*x*y*z*z - 4.0*x*x*y*z*z - 4.0*x*y*y*z*z + 8.0*x*x*y*y*z*z;
        S(8,i)= + 1.0*y*z - 3.0*x*y*z + 2.0*x*x*y*z - 2.0*y*y*z + 6.0*x*y*y*z - 4.0*x*x*y*y*z - 2.0*y*z*z + 6.0*x*y*z*z - 4.0*x*x*y*z*z + 4.0*y*y*z*z - 12.0*x*y*y*z*z + 8.0*x*x*y*y*z*z;
        S(9,i)= + 4.0*x - 4.0*x*x - 12.0*x*y + 12.0*x*x*y + 8.0*x*y*y - 8.0*x*x*y*y - 12.0*x*z + 12.0*x*x*z + 36.0*x*y*z - 36.0*x*x*y*z - 24.0*x*y*y*z + 24.0*x*x*y*y*z + 8.0*x*z*z - 8.0*x*x*z*z - 24.0*x*y*z*z + 24.0*x*x*y*z*z + 16.0*x*y*y*z*z - 16.0*x*x*y*y*z*z;
        S(10,i)= - 4.0*x*y + 8.0*x*x*y + 4.0*x*y*y - 8.0*x*x*y*y + 12.0*x*y*z - 24.0*x*x*y*z - 12.0*x*y*y*z + 24.0*x*x*y*y*z - 8.0*x*y*z*z + 16.0*x*x*y*z*z + 8.0*x*y*y*z*z - 16.0*x*x*y*y*z*z;
        S(11,i)= - 4.0*x*y + 4.0*x*x*y + 8.0*x*y*y - 8.0*x*x*y*y + 12.0*x*y*z - 12.0*x*x*y*z - 24.0*x*y*y*z + 24.0*x*x*y*y*z - 8.0*x*y*z*z + 8.0*x*x*y*z*z + 16.0*x*y*y*z*z - 16.0*x*x*y*y*z*z;
        S(12,i)= + 4.0*y - 12.0*x*y + 8.0*x*x*y - 4.0*y*y + 12.0*x*y*y - 8.0*x*x*y*y - 12.0*y*z + 36.0*x*y*z - 24.0*x*x*y*z + 12.0*y*y*z - 36.0*x*y*y*z + 24.0*x*x*y*y*z + 8.0*y*z*z - 24.0*x*y*z*z + 16.0*x*x*y*z*z - 8.0*y*y*z*z + 24.0*x*y*y*z*z - 16.0*x*x*y*y*z*z;
        S(13,i)= + 4.0*z - 12.0*x*z + 8.0*x*x*z - 12.0*y*z + 36.0*x*y*z - 24.0*x*x*y*z + 8.0*y*y*z - 24.0*x*y*y*z + 16.0*x*x*y*y*z - 4.0*z*z + 12.0*x*z*z - 8.0*x*x*z*z + 12.0*y*z*z - 36.0*x*y*z*z + 24.0*x*x*y*z*z - 8.0*y*y*z*z + 24.0*x*y*y*z*z - 16.0*x*x*y*y*z*z;
        S(14,i)= - 4.0*x*z + 8.0*x*x*z + 12.0*x*y*z - 24.0*x*x*y*z - 8.0*x*y*y*z + 16.0*x*x*y*y*z + 4.0*x*z*z - 8.0*x*x*z*z - 12.0*x*y*z*z + 24.0*x*x*y*z*z + 8.0*x*y*y*z*z - 16.0*x*x*y*y*z*z;
        S(15,i)= + 4.0*x*y*z - 8.0*x*x*y*z - 8.0*x*y*y*z + 16.0*x*x*y*y*z - 4.0*x*y*z*z + 8.0*x*x*y*z*z + 8.0*x*y*y*z*z - 16.0*x*x*y*y*z*z;
        S(16,i)= - 4.0*y*z + 12.0*x*y*z - 8.0*x*x*y*z + 8.0*y*y*z - 24.0*x*y*y*z + 16.0*x*x*y*y*z + 4.0*y*z*z - 12.0*x*y*z*z + 8.0*x*x*y*z*z - 8.0*y*y*z*z + 24.0*x*y*y*z*z - 16.0*x*x*y*y*z*z;
        S(17,i)= - 4.0*x*z + 4.0*x*x*z + 12.0*x*y*z - 12.0*x*x*y*z - 8.0*x*y*y*z + 8.0*x*x*y*y*z + 8.0*x*z*z - 8.0*x*x*z*z - 24.0*x*y*z*z + 24.0*x*x*y*z*z + 16.0*x*y*y*z*z - 16.0*x*x*y*y*z*z;
        S(18,i)= + 4.0*x*y*z - 8.0*x*x*y*z - 4.0*x*y*y*z + 8.0*x*x*y*y*z - 8.0*x*y*z*z + 16.0*x*x*y*z*z + 8.0*x*y*y*z*z - 16.0*x*x*y*y*z*z;
        S(19,i)= + 4.0*x*y*z - 4.0*x*x*y*z - 8.0*x*y*y*z + 8.0*x*x*y*y*z - 8.0*x*y*z*z + 8.0*x*x*y*z*z + 16.0*x*y*y*z*z - 16.0*x*x*y*y*z*z;
        S(20,i)= - 4.0*y*z + 12.0*x*y*z - 8.0*x*x*y*z + 4.0*y*y*z - 12.0*x*y*y*z + 8.0*x*x*y*y*z + 8.0*y*z*z - 24.0*x*y*z*z + 16.0*x*x*y*z*z - 8.0*y*y*z*z + 24.0*x*y*y*z*z - 16.0*x*x*y*y*z*z;
        S(21,i)= + 16.0*x*y - 16.0*x*x*y - 16.0*x*y*y + 16.0*x*x*y*y - 48.0*x*y*z + 48.0*x*x*y*z + 48.0*x*y*y*z - 48.0*x*x*y*y*z + 32.0*x*y*z*z - 32.0*x*x*y*z*z - 32.0*x*y*y*z*z + 32.0*x*x*y*y*z*z;
        S(22,i)= + 16.0*x*z - 16.0*x*x*z - 48.0*x*y*z + 48.0*x*x*y*z + 32.0*x*y*y*z - 32.0*x*x*y*y*z - 16.0*x*z*z + 16.0*x*x*z*z + 48.0*x*y*z*z - 48.0*x*x*y*z*z - 32.0*x*y*y*z*z + 32.0*x*x*y*y*z*z;
        S(23,i)= - 16.0*x*y*z + 32.0*x*x*y*z + 16.0*x*y*y*z - 32.0*x*x*y*y*z + 16.0*x*y*z*z - 32.0*x*x*y*z*z - 16.0*x*y*y*z*z + 32.0*x*x*y*y*z*z;
        S(24,i)= - 16.0*x*y*z + 16.0*x*x*y*z + 32.0*x*y*y*z - 32.0*x*x*y*y*z + 16.0*x*y*z*z - 16.0*x*x*y*z*z - 32.0*x*y*y*z*z + 32.0*x*x*y*y*z*z;
        S(25,i)= + 16.0*y*z - 48.0*x*y*z + 32.0*x*x*y*z - 16.0*y*y*z + 48.0*x*y*y*z - 32.0*x*x*y*y*z - 16.0*y*z*z + 48.0*x*y*z*z - 32.0*x*x*y*z*z + 16.0*y*y*z*z - 48.0*x*y*y*z*z + 32.0*x*x*y*y*z*z;
        S(26,i)= - 16.0*x*y*z + 16.0*x*x*y*z + 16.0*x*y*y*z - 16.0*x*x*y*y*z + 32.0*x*y*z*z - 32.0*x*x*y*z*z - 32.0*x*y*y*z*z + 32.0*x*x*y*y*z*z;
        S(27,i)= + 64.0*x*y*z - 64.0*x*x*y*z - 64.0*x*y*y*z + 64.0*x*x*y*y*z - 64.0*x*y*z*z + 64.0*x*x*y*z*z + 64.0*x*y*y*z*z - 64.0*x*x*y*y*z*z;
        DSDV(1,1,i)= - 3.0 + 4.0*x + 9.0*y - 12.0*x*y - 6.0*y*y + 8.0*x*y*y + 9.0*z - 12.0*x*z - 27.0*y*z + 36.0*x*y*z + 18.0*y*y*z - 24.0*x*y*y*z - 6.0*z*z + 8.0*x*z*z + 18.0*y*z*z - 24.0*x*y*z*z - 12.0*y*y*z*z + 16.0*x*y*y*z*z;
        DSDV(1,2,i)= - 3.0 + 9.0*x - 6.0*x*x + 4.0*y - 12.0*x*y + 8.0*x*x*y + 9.0*z - 27.0*x*z + 18.0*x*x*z - 12.0*y*z + 36.0*x*y*z - 24.0*x*x*y*z - 6.0*z*z + 18.0*x*z*z - 12.0*x*x*z*z + 8.0*y*z*z - 24.0*x*y*z*z + 16.0*x*x*y*z*z;
        DSDV(1,3,i)= - 3.0 + 9.0*x - 6.0*x*x + 9.0*y - 27.0*x*y + 18.0*x*x*y - 6.0*y*y + 18.0*x*y*y - 12.0*x*x*y*y + 4.0*z - 12.0*x*z + 8.0*x*x*z - 12.0*y*z + 36.0*x*y*z - 24.0*x*x*y*z + 8.0*y*y*z - 24.0*x*y*y*z + 16.0*x*x*y*y*z;
        DSDV(2,1,i)= - 1.0 + 4.0*x + 3.0*y - 12.0*x*y - 2.0*y*y + 8.0*x*y*y + 3.0*z - 12.0*x*z - 9.0*y*z + 36.0*x*y*z + 6.0*y*y*z - 24.0*x*y*y*z - 2.0*z*z + 8.0*x*z*z + 6.0*y*z*z - 24.0*x*y*z*z - 4.0*y*y*z*z + 16.0*x*y*y*z*z;
        DSDV(2,2,i)= + 3.0*x - 6.0*x*x - 4.0*x*y + 8.0*x*x*y - 9.0*x*z + 18.0*x*x*z + 12.0*x*y*z - 24.0*x*x*y*z + 6.0*x*z*z - 12.0*x*x*z*z - 8.0*x*y*z*z + 16.0*x*x*y*z*z;
        DSDV(2,3,i)= + 3.0*x - 6.0*x*x - 9.0*x*y + 18.0*x*x*y + 6.0*x*y*y - 12.0*x*x*y*y - 4.0*x*z + 8.0*x*x*z + 12.0*x*y*z - 24.0*x*x*y*z - 8.0*x*y*y*z + 16.0*x*x*y*y*z;
        DSDV(3,1,i)= + 1.0*y - 4.0*x*y - 2.0*y*y + 8.0*x*y*y - 3.0*y*z + 12.0*x*y*z + 6.0*y*y*z - 24.0*x*y*y*z + 2.0*y*z*z - 8.0*x*y*z*z - 4.0*y*y*z*z + 16.0*x*y*y*z*z;
        DSDV(3,2,i)= + 1.0*x - 2.0*x*x - 4.0*x*y + 8.0*x*x*y - 3.0*x*z + 6.0*x*x*z + 12.0*x*y*z - 24.0*x*x*y*z + 2.0*x*z*z - 4.0*x*x*z*z - 8.0*x*y*z*z + 16.0*x*x*y*z*z;
        DSDV(3,3,i)= - 3.0*x*y + 6.0*x*x*y + 6.0*x*y*y - 12.0*x*x*y*y + 4.0*x*y*z - 8.0*x*x*y*z - 8.0*x*y*y*z + 16.0*x*x*y*y*z;
        DSDV(4,1,i)= + 3.0*y - 4.0*x*y - 6.0*y*y + 8.0*x*y*y - 9.0*y*z + 12.0*x*y*z + 18.0*y*y*z - 24.0*x*y*y*z + 6.0*y*z*z - 8.0*x*y*z*z - 12.0*y*y*z*z + 16.0*x*y*y*z*z;
        DSDV(4,2,i)= - 1.0 + 3.0*x - 2.0*x*x + 4.0*y - 12.0*x*y + 8.0*x*x*y + 3.0*z - 9.0*x*z + 6.0*x*x*z - 12.0*y*z + 36.0*x*y*z - 24.0*x*x*y*z - 2.0*z*z + 6.0*x*z*z - 4.0*x*x*z*z + 8.0*y*z*z - 24.0*x*y*z*z + 16.0*x*x*y*z*z;
        DSDV(4,3,i)= + 3.0*y - 9.0*x*y + 6.0*x*x*y - 6.0*y*y + 18.0*x*y*y - 12.0*x*x*y*y - 4.0*y*z + 12.0*x*y*z - 8.0*x*x*y*z + 8.0*y*y*z - 24.0*x*y*y*z + 16.0*x*x*y*y*z;
        DSDV(5,1,i)= + 3.0*z - 4.0*x*z - 9.0*y*z + 12.0*x*y*z + 6.0*y*y*z - 8.0*x*y*y*z - 6.0*z*z + 8.0*x*z*z + 18.0*y*z*z - 24.0*x*y*z*z - 12.0*y*y*z*z + 16.0*x*y*y*z*z;
        DSDV(5,2,i)= + 3.0*z - 9.0*x*z + 6.0*x*x*z - 4.0*y*z + 12.0*x*y*z - 8.0*x*x*y*z - 6.0*z*z + 18.0*x*z*z - 12.0*x*x*z*z + 8.0*y*z*z - 24.0*x*y*z*z + 16.0*x*x*y*z*z;
        DSDV(5,3,i)= - 1.0 + 3.0*x - 2.0*x*x + 3.0*y - 9.0*x*y + 6.0*x*x*y - 2.0*y*y + 6.0*x*y*y - 4.0*x*x*y*y + 4.0*z - 12.0*x*z + 8.0*x*x*z - 12.0*y*z + 36.0*x*y*z - 24.0*x*x*y*z + 8.0*y*y*z - 24.0*x*y*y*z + 16.0*x*x*y*y*z;
        DSDV(6,1,i)= + 1.0*z - 4.0*x*z - 3.0*y*z + 12.0*x*y*z + 2.0*y*y*z - 8.0*x*y*y*z - 2.0*z*z + 8.0*x*z*z + 6.0*y*z*z - 24.0*x*y*z*z - 4.0*y*y*z*z + 16.0*x*y*y*z*z;
        DSDV(6,2,i)= - 3.0*x*z + 6.0*x*x*z + 4.0*x*y*z - 8.0*x*x*y*z + 6.0*x*z*z - 12.0*x*x*z*z - 8.0*x*y*z*z + 16.0*x*x*y*z*z;
        DSDV(6,3,i)= + 1.0*x - 2.0*x*x - 3.0*x*y + 6.0*x*x*y + 2.0*x*y*y - 4.0*x*x*y*y - 4.0*x*z + 8.0*x*x*z + 12.0*x*y*z - 24.0*x*x*y*z - 8.0*x*y*y*z + 16.0*x*x*y*y*z;
        DSDV(7,1,i)= - 1.0*y*z + 4.0*x*y*z + 2.0*y*y*z - 8.0*x*y*y*z + 2.0*y*z*z - 8.0*x*y*z*z - 4.0*y*y*z*z + 16.0*x*y*y*z*z;
        DSDV(7,2,i)= - 1.0*x*z + 2.0*x*x*z + 4.0*x*y*z - 8.0*x*x*y*z + 2.0*x*z*z - 4.0*x*x*z*z - 8.0*x*y*z*z + 16.0*x*x*y*z*z;
        DSDV(7,3,i)= - 1.0*x*y + 2.0*x*x*y + 2.0*x*y*y - 4.0*x*x*y*y + 4.0*x*y*z - 8.0*x*x*y*z - 8.0*x*y*y*z + 16.0*x*x*y*y*z;
        DSDV(8,1,i)= - 3.0*y*z + 4.0*x*y*z + 6.0*y*y*z - 8.0*x*y*y*z + 6.0*y*z*z - 8.0*x*y*z*z - 12.0*y*y*z*z + 16.0*x*y*y*z*z;
        DSDV(8,2,i)= + 1.0*z - 3.0*x*z + 2.0*x*x*z - 4.0*y*z + 12.0*x*y*z - 8.0*x*x*y*z - 2.0*z*z + 6.0*x*z*z - 4.0*x*x*z*z + 8.0*y*z*z - 24.0*x*y*z*z + 16.0*x*x*y*z*z;
        DSDV(8,3,i)= + 1.0*y - 3.0*x*y + 2.0*x*x*y - 2.0*y*y + 6.0*x*y*y - 4.0*x*x*y*y - 4.0*y*z + 12.0*x*y*z - 8.0*x*x*y*z + 8.0*y*y*z - 24.0*x*y*y*z + 16.0*x*x*y*y*z;
        DSDV(9,1,i)= + 4.0 - 8.0*x - 12.0*y + 24.0*x*y + 8.0*y*y - 16.0*x*y*y - 12.0*z + 24.0*x*z + 36.0*y*z - 72.0*x*y*z - 24.0*y*y*z + 48.0*x*y*y*z + 8.0*z*z - 16.0*x*z*z - 24.0*y*z*z + 48.0*x*y*z*z + 16.0*y*y*z*z - 32.0*x*y*y*z*z;
        DSDV(9,2,i)= - 12.0*x + 12.0*x*x + 16.0*x*y - 16.0*x*x*y + 36.0*x*z - 36.0*x*x*z - 48.0*x*y*z + 48.0*x*x*y*z - 24.0*x*z*z + 24.0*x*x*z*z + 32.0*x*y*z*z - 32.0*x*x*y*z*z;
        DSDV(9,3,i)= - 12.0*x + 12.0*x*x + 36.0*x*y - 36.0*x*x*y - 24.0*x*y*y + 24.0*x*x*y*y + 16.0*x*z - 16.0*x*x*z - 48.0*x*y*z + 48.0*x*x*y*z + 32.0*x*y*y*z - 32.0*x*x*y*y*z;
        DSDV(10,1,i)= - 4.0*y + 16.0*x*y + 4.0*y*y - 16.0*x*y*y + 12.0*y*z - 48.0*x*y*z - 12.0*y*y*z + 48.0*x*y*y*z - 8.0*y*z*z + 32.0*x*y*z*z + 8.0*y*y*z*z - 32.0*x*y*y*z*z;
        DSDV(10,2,i)= - 4.0*x + 8.0*x*x + 8.0*x*y - 16.0*x*x*y + 12.0*x*z - 24.0*x*x*z - 24.0*x*y*z + 48.0*x*x*y*z - 8.0*x*z*z + 16.0*x*x*z*z + 16.0*x*y*z*z - 32.0*x*x*y*z*z;
        DSDV(10,3,i)= + 12.0*x*y - 24.0*x*x*y - 12.0*x*y*y + 24.0*x*x*y*y - 16.0*x*y*z + 32.0*x*x*y*z + 16.0*x*y*y*z - 32.0*x*x*y*y*z;
        DSDV(11,1,i)= - 4.0*y + 8.0*x*y + 8.0*y*y - 16.0*x*y*y + 12.0*y*z - 24.0*x*y*z - 24.0*y*y*z + 48.0*x*y*y*z - 8.0*y*z*z + 16.0*x*y*z*z + 16.0*y*y*z*z - 32.0*x*y*y*z*z;
        DSDV(11,2,i)= - 4.0*x + 4.0*x*x + 16.0*x*y - 16.0*x*x*y + 12.0*x*z - 12.0*x*x*z - 48.0*x*y*z + 48.0*x*x*y*z - 8.0*x*z*z + 8.0*x*x*z*z + 32.0*x*y*z*z - 32.0*x*x*y*z*z;
        DSDV(11,3,i)= + 12.0*x*y - 12.0*x*x*y - 24.0*x*y*y + 24.0*x*x*y*y - 16.0*x*y*z + 16.0*x*x*y*z + 32.0*x*y*y*z - 32.0*x*x*y*y*z;
        DSDV(12,1,i)= - 12.0*y + 16.0*x*y + 12.0*y*y - 16.0*x*y*y + 36.0*y*z - 48.0*x*y*z - 36.0*y*y*z + 48.0*x*y*y*z - 24.0*y*z*z + 32.0*x*y*z*z + 24.0*y*y*z*z - 32.0*x*y*y*z*z;
        DSDV(12,2,i)= + 4.0 - 12.0*x + 8.0*x*x - 8.0*y + 24.0*x*y - 16.0*x*x*y - 12.0*z + 36.0*x*z - 24.0*x*x*z + 24.0*y*z - 72.0*x*y*z + 48.0*x*x*y*z + 8.0*z*z - 24.0*x*z*z + 16.0*x*x*z*z - 16.0*y*z*z + 48.0*x*y*z*z - 32.0*x*x*y*z*z;
        DSDV(12,3,i)= - 12.0*y + 36.0*x*y - 24.0*x*x*y + 12.0*y*y - 36.0*x*y*y + 24.0*x*x*y*y + 16.0*y*z - 48.0*x*y*z + 32.0*x*x*y*z - 16.0*y*y*z + 48.0*x*y*y*z - 32.0*x*x*y*y*z;
        DSDV(13,1,i)= - 12.0*z + 16.0*x*z + 36.0*y*z - 48.0*x*y*z - 24.0*y*y*z + 32.0*x*y*y*z + 12.0*z*z - 16.0*x*z*z - 36.0*y*z*z + 48.0*x*y*z*z + 24.0*y*y*z*z - 32.0*x*y*y*z*z;
        DSDV(13,2,i)= - 12.0*z + 36.0*x*z - 24.0*x*x*z + 16.0*y*z - 48.0*x*y*z + 32.0*x*x*y*z + 12.0*z*z - 36.0*x*z*z + 24.0*x*x*z*z - 16.0*y*z*z + 48.0*x*y*z*z - 32.0*x*x*y*z*z;
        DSDV(13,3,i)= + 4.0 - 12.0*x + 8.0*x*x - 12.0*y + 36.0*x*y - 24.0*x*x*y + 8.0*y*y - 24.0*x*y*y + 16.0*x*x*y*y - 8.0*z + 24.0*x*z - 16.0*x*x*z + 24.0*y*z - 72.0*x*y*z + 48.0*x*x*y*z - 16.0*y*y*z + 48.0*x*y*y*z - 32.0*x*x*y*y*z;
        DSDV(14,1,i)= - 4.0*z + 16.0*x*z + 12.0*y*z - 48.0*x*y*z - 8.0*y*y*z + 32.0*x*y*y*z + 4.0*z*z - 16.0*x*z*z - 12.0*y*z*z + 48.0*x*y*z*z + 8.0*y*y*z*z - 32.0*x*y*y*z*z;
        DSDV(14,2,i)= + 12.0*x*z - 24.0*x*x*z - 16.0*x*y*z + 32.0*x*x*y*z - 12.0*x*z*z + 24.0*x*x*z*z + 16.0*x*y*z*z - 32.0*x*x*y*z*z;
        DSDV(14,3,i)= - 4.0*x + 8.0*x*x + 12.0*x*y - 24.0*x*x*y - 8.0*x*y*y + 16.0*x*x*y*y + 8.0*x*z - 16.0*x*x*z - 24.0*x*y*z + 48.0*x*x*y*z + 16.0*x*y*y*z - 32.0*x*x*y*y*z;
        DSDV(15,1,i)= + 4.0*y*z - 16.0*x*y*z - 8.0*y*y*z + 32.0*x*y*y*z - 4.0*y*z*z + 16.0*x*y*z*z + 8.0*y*y*z*z - 32.0*x*y*y*z*z;
        DSDV(15,2,i)= + 4.0*x*z - 8.0*x*x*z - 16.0*x*y*z + 32.0*x*x*y*z - 4.0*x*z*z + 8.0*x*x*z*z + 16.0*x*y*z*z - 32.0*x*x*y*z*z;
        DSDV(15,3,i)= + 4.0*x*y - 8.0*x*x*y - 8.0*x*y*y + 16.0*x*x*y*y - 8.0*x*y*z + 16.0*x*x*y*z + 16.0*x*y*y*z - 32.0*x*x*y*y*z;
        DSDV(16,1,i)= + 12.0*y*z - 16.0*x*y*z - 24.0*y*y*z + 32.0*x*y*y*z - 12.0*y*z*z + 16.0*x*y*z*z + 24.0*y*y*z*z - 32.0*x*y*y*z*z;
        DSDV(16,2,i)= - 4.0*z + 12.0*x*z - 8.0*x*x*z + 16.0*y*z - 48.0*x*y*z + 32.0*x*x*y*z + 4.0*z*z - 12.0*x*z*z + 8.0*x*x*z*z - 16.0*y*z*z + 48.0*x*y*z*z - 32.0*x*x*y*z*z;
        DSDV(16,3,i)= - 4.0*y + 12.0*x*y - 8.0*x*x*y + 8.0*y*y - 24.0*x*y*y + 16.0*x*x*y*y + 8.0*y*z - 24.0*x*y*z + 16.0*x*x*y*z - 16.0*y*y*z + 48.0*x*y*y*z - 32.0*x*x*y*y*z;
        DSDV(17,1,i)= - 4.0*z + 8.0*x*z + 12.0*y*z - 24.0*x*y*z - 8.0*y*y*z + 16.0*x*y*y*z + 8.0*z*z - 16.0*x*z*z - 24.0*y*z*z + 48.0*x*y*z*z + 16.0*y*y*z*z - 32.0*x*y*y*z*z;
        DSDV(17,2,i)= + 12.0*x*z - 12.0*x*x*z - 16.0*x*y*z + 16.0*x*x*y*z - 24.0*x*z*z + 24.0*x*x*z*z + 32.0*x*y*z*z - 32.0*x*x*y*z*z;
        DSDV(17,3,i)= - 4.0*x + 4.0*x*x + 12.0*x*y - 12.0*x*x*y - 8.0*x*y*y + 8.0*x*x*y*y + 16.0*x*z - 16.0*x*x*z - 48.0*x*y*z + 48.0*x*x*y*z + 32.0*x*y*y*z - 32.0*x*x*y*y*z;
        DSDV(18,1,i)= + 4.0*y*z - 16.0*x*y*z - 4.0*y*y*z + 16.0*x*y*y*z - 8.0*y*z*z + 32.0*x*y*z*z + 8.0*y*y*z*z - 32.0*x*y*y*z*z;
        DSDV(18,2,i)= + 4.0*x*z - 8.0*x*x*z - 8.0*x*y*z + 16.0*x*x*y*z - 8.0*x*z*z + 16.0*x*x*z*z + 16.0*x*y*z*z - 32.0*x*x*y*z*z;
        DSDV(18,3,i)= + 4.0*x*y - 8.0*x*x*y - 4.0*x*y*y + 8.0*x*x*y*y - 16.0*x*y*z + 32.0*x*x*y*z + 16.0*x*y*y*z - 32.0*x*x*y*y*z;
        DSDV(19,1,i)= + 4.0*y*z - 8.0*x*y*z - 8.0*y*y*z + 16.0*x*y*y*z - 8.0*y*z*z + 16.0*x*y*z*z + 16.0*y*y*z*z - 32.0*x*y*y*z*z;
        DSDV(19,2,i)= + 4.0*x*z - 4.0*x*x*z - 16.0*x*y*z + 16.0*x*x*y*z - 8.0*x*z*z + 8.0*x*x*z*z + 32.0*x*y*z*z - 32.0*x*x*y*z*z;
        DSDV(19,3,i)= + 4.0*x*y - 4.0*x*x*y - 8.0*x*y*y + 8.0*x*x*y*y - 16.0*x*y*z + 16.0*x*x*y*z + 32.0*x*y*y*z - 32.0*x*x*y*y*z;
        DSDV(20,1,i)= + 12.0*y*z - 16.0*x*y*z - 12.0*y*y*z + 16.0*x*y*y*z - 24.0*y*z*z + 32.0*x*y*z*z + 24.0*y*y*z*z - 32.0*x*y*y*z*z;
        DSDV(20,2,i)= - 4.0*z + 12.0*x*z - 8.0*x*x*z + 8.0*y*z - 24.0*x*y*z + 16.0*x*x*y*z + 8.0*z*z - 24.0*x*z*z + 16.0*x*x*z*z - 16.0*y*z*z + 48.0*x*y*z*z - 32.0*x*x*y*z*z;
        DSDV(20,3,i)= - 4.0*y + 12.0*x*y - 8.0*x*x*y + 4.0*y*y - 12.0*x*y*y + 8.0*x*x*y*y + 16.0*y*z - 48.0*x*y*z + 32.0*x*x*y*z - 16.0*y*y*z + 48.0*x*y*y*z - 32.0*x*x*y*y*z;
        DSDV(21,1,i)= + 16.0*y - 32.0*x*y - 16.0*y*y + 32.0*x*y*y - 48.0*y*z + 96.0*x*y*z + 48.0*y*y*z - 96.0*x*y*y*z + 32.0*y*z*z - 64.0*x*y*z*z - 32.0*y*y*z*z + 64.0*x*y*y*z*z;
        DSDV(21,2,i)= + 16.0*x - 16.0*x*x - 32.0*x*y + 32.0*x*x*y - 48.0*x*z + 48.0*x*x*z + 96.0*x*y*z - 96.0*x*x*y*z + 32.0*x*z*z - 32.0*x*x*z*z - 64.0*x*y*z*z + 64.0*x*x*y*z*z;
        DSDV(21,3,i)= - 48.0*x*y + 48.0*x*x*y + 48.0*x*y*y - 48.0*x*x*y*y + 64.0*x*y*z - 64.0*x*x*y*z - 64.0*x*y*y*z + 64.0*x*x*y*y*z;
        DSDV(22,1,i)= + 16.0*z - 32.0*x*z - 48.0*y*z + 96.0*x*y*z + 32.0*y*y*z - 64.0*x*y*y*z - 16.0*z*z + 32.0*x*z*z + 48.0*y*z*z - 96.0*x*y*z*z - 32.0*y*y*z*z + 64.0*x*y*y*z*z;
        DSDV(22,2,i)= - 48.0*x*z + 48.0*x*x*z + 64.0*x*y*z - 64.0*x*x*y*z + 48.0*x*z*z - 48.0*x*x*z*z - 64.0*x*y*z*z + 64.0*x*x*y*z*z;
        DSDV(22,3,i)= + 16.0*x - 16.0*x*x - 48.0*x*y + 48.0*x*x*y + 32.0*x*y*y - 32.0*x*x*y*y - 32.0*x*z + 32.0*x*x*z + 96.0*x*y*z - 96.0*x*x*y*z - 64.0*x*y*y*z + 64.0*x*x*y*y*z;
        DSDV(23,1,i)= - 16.0*y*z + 64.0*x*y*z + 16.0*y*y*z - 64.0*x*y*y*z + 16.0*y*z*z - 64.0*x*y*z*z - 16.0*y*y*z*z + 64.0*x*y*y*z*z;
        DSDV(23,2,i)= - 16.0*x*z + 32.0*x*x*z + 32.0*x*y*z - 64.0*x*x*y*z + 16.0*x*z*z - 32.0*x*x*z*z - 32.0*x*y*z*z + 64.0*x*x*y*z*z;
        DSDV(23,3,i)= - 16.0*x*y + 32.0*x*x*y + 16.0*x*y*y - 32.0*x*x*y*y + 32.0*x*y*z - 64.0*x*x*y*z - 32.0*x*y*y*z + 64.0*x*x*y*y*z;
        DSDV(24,1,i)= - 16.0*y*z + 32.0*x*y*z + 32.0*y*y*z - 64.0*x*y*y*z + 16.0*y*z*z - 32.0*x*y*z*z - 32.0*y*y*z*z + 64.0*x*y*y*z*z;
        DSDV(24,2,i)= - 16.0*x*z + 16.0*x*x*z + 64.0*x*y*z - 64.0*x*x*y*z + 16.0*x*z*z - 16.0*x*x*z*z - 64.0*x*y*z*z + 64.0*x*x*y*z*z;
        DSDV(24,3,i)= - 16.0*x*y + 16.0*x*x*y + 32.0*x*y*y - 32.0*x*x*y*y + 32.0*x*y*z - 32.0*x*x*y*z - 64.0*x*y*y*z + 64.0*x*x*y*y*z;
        DSDV(25,1,i)= - 48.0*y*z + 64.0*x*y*z + 48.0*y*y*z - 64.0*x*y*y*z + 48.0*y*z*z - 64.0*x*y*z*z - 48.0*y*y*z*z + 64.0*x*y*y*z*z;
        DSDV(25,2,i)= + 16.0*z - 48.0*x*z + 32.0*x*x*z - 32.0*y*z + 96.0*x*y*z - 64.0*x*x*y*z - 16.0*z*z + 48.0*x*z*z - 32.0*x*x*z*z + 32.0*y*z*z - 96.0*x*y*z*z + 64.0*x*x*y*z*z;
        DSDV(25,3,i)= + 16.0*y - 48.0*x*y + 32.0*x*x*y - 16.0*y*y + 48.0*x*y*y - 32.0*x*x*y*y - 32.0*y*z + 96.0*x*y*z - 64.0*x*x*y*z + 32.0*y*y*z - 96.0*x*y*y*z + 64.0*x*x*y*y*z;
        DSDV(26,1,i)= - 16.0*y*z + 32.0*x*y*z + 16.0*y*y*z - 32.0*x*y*y*z + 32.0*y*z*z - 64.0*x*y*z*z - 32.0*y*y*z*z + 64.0*x*y*y*z*z;
        DSDV(26,2,i)= - 16.0*x*z + 16.0*x*x*z + 32.0*x*y*z - 32.0*x*x*y*z + 32.0*x*z*z - 32.0*x*x*z*z - 64.0*x*y*z*z + 64.0*x*x*y*z*z;
        DSDV(26,3,i)= - 16.0*x*y + 16.0*x*x*y + 16.0*x*y*y - 16.0*x*x*y*y + 64.0*x*y*z - 64.0*x*x*y*z - 64.0*x*y*y*z + 64.0*x*x*y*y*z;
        DSDV(27,1,i)= + 64.0*y*z - 128.0*x*y*z - 64.0*y*y*z + 128.0*x*y*y*z - 64.0*y*z*z + 128.0*x*y*z*z + 64.0*y*y*z*z - 128.0*x*y*y*z*z;
        DSDV(27,2,i)= + 64.0*x*z - 64.0*x*x*z - 128.0*x*y*z + 128.0*x*x*y*z - 64.0*x*z*z + 64.0*x*x*z*z + 128.0*x*y*z*z - 128.0*x*x*y*z*z;
        DSDV(27,3,i)= + 64.0*x*y - 64.0*x*x*y - 64.0*x*y*y + 64.0*x*x*y*y - 128.0*x*y*z + 128.0*x*x*y*z + 128.0*x*y*y*z - 128.0*x*x*y*y*z;
    }
#undef NUMSHAPES
#undef DIM
}

/****************************************************************************/
void Shape_Hex32(int NumV, std::vector<double>& v, std::vector<double>& s, std::vector<double>& dsdv)
{
#define NUMSHAPES 32
#define DIM 3
    #pragma ivdep
    for (int i=0; i<NumV; i++) {
        const double x=V(1,i);
        const double y=V(2,i);
        const double z=V(3,i);
        S(1,i)=(10.)+(-5.5)*x+(10.)*x*y+(-5.5)*y+(-5.5)*z+(10.)*x*z+(-0.145e2)*x*y*z+(10.)*y*z+(9.)*x*x+(-4.5)*x*x*x+(-9.)*x*y*y+(4.5)*x*y*y*y+(4.5)*x*x*x*y+(-9.)*x*x*y+(-4.5)*y*y*y+(9.)*y*y+(9.)*z*z+(-9.)*x*z*z+(9.)*x*y*z*z+(-9.)*y*z*z+(-4.5)*z*z*z+(4.5)*x*z*z*z+(-4.5)*x*y*z*z*z+(4.5)*y*z*z*z+(-9.)*x*x*z+(4.5)*x*x*x*z+(9.)*x*y*y*z+(-4.5)*x*y*y*y*z+(-4.5)*x*x*x*y*z+(9.)*x*x*y*z+(4.5)*y*y*y*z+(-9.)*y*y*z;
        S(2,i)=(10.)*x+(-5.5)*x*y+(-5.5)*x*z+(10.)*x*y*z+(-.45e1)*x*x+(4.5)*x*x*x+(9.)*x*y*y+(-4.5)*x*y*y*y+(-4.5)*x*x*x*y+(4.5)*x*x*y+(9.)*x*z*z+(-9.)*x*y*z*z+(-4.5)*x*z*z*z+(4.5)*x*y*z*z*z+(4.5)*x*x*z+(-4.5)*x*x*x*z+(-9.)*x*y*y*z+(4.5)*x*y*y*y*z+(4.5)*x*x*x*y*z+(-4.5)*x*x*y*z;
        S(3,i)=(10.)*x*y+(-5.5)*x*y*z+(-4.5)*x*y*y+(4.5)*x*y*y*y+(4.5)*x*x*x*y+(-4.5)*x*x*y+(9.)*x*y*z*z+(-4.5)*x*y*z*z*z+(4.5)*x*y*y*z+(-4.5)*x*y*y*y*z+(-4.5)*x*x*x*y*z+(4.5)*x*x*y*z;
        S(4,i)=(-5.5)*x*y+(10.)*y+(10.)*x*y*z+(-5.5)*y*z+(.45e1)*x*y*y+(-4.5)*x*y*y*y+(-4.5)*x*x*x*y+(9.)*x*x*y+(4.5)*y*y*y+(-4.5)*y*y+(-9.)*x*y*z*z+(9.)*y*z*z+(4.5)*x*y*z*z*z+(-4.5)*y*z*z*z+(-4.5)*x*y*y*z+(4.5)*x*y*y*y*z+(4.5)*x*x*x*y*z+(-9.)*x*x*y*z+(-4.5)*y*y*y*z+(4.5)*y*y*z;
        S(5,i)=(10.)*z+(-5.5)*x*z+(10.)*x*y*z+(-5.5)*y*z+(-.45e1)*z*z+(4.5)*x*z*z+(-4.5)*x*y*z*z+(4.5)*y*z*z+(4.5)*z*z*z+(-4.5)*x*z*z*z+(4.5)*x*y*z*z*z+(-4.5)*y*z*z*z+(9.)*x*x*z+(-4.5)*x*x*x*z+(-9.)*x*y*y*z+(4.5)*x*y*y*y*z+(4.5)*x*x*x*y*z+(-9.)*x*x*y*z+(-4.5)*y*y*y*z+(9.)*y*y*z;
        S(6,i)=(10.)*x*z+(-5.5)*x*y*z+(-4.5)*x*z*z+(4.5)*x*y*z*z+(4.5)*x*z*z*z+(-4.5)*x*y*z*z*z+(-4.5)*x*x*z+(4.5)*x*x*x*z+(9.)*x*y*y*z+(-.45+01)*x*y*y*y*z+(-4.5)*x*x*x*y*z+(4.5)*x*x*y*z;
        S(7,i)=(10.)*x*y*z+(-4.5)*x*y*z*z+(4.5)*x*y*z*z*z+(-4.5)*x*y*y*z+(4.5)*x*y*y*y*z+(4.5)*x*x*x*y*z+(-4.5)*x*x*y*z;
        S(8,i)=(-5.5)*x*y*z+(10.)*y*z+(4.5)*x*y*z*z+(-4.5)*y*z*z+(-4.5)*x*y*z*z*z+(4.5)*y*z*z*z+(4.5)*x*y*y*z+(-4.5)*x*y*y*y*z+(-4.5)*x*x*x*y*z+(9.)*x*x*y*z+(4.5)*y*y*y*z+(-4.5)*y*y*z;
        S(9,i)=(9.)*x+(-9.)*x*y+(-9.)*x*z+(9.)*x*y*z+(-.225e2)*x*x+(13.5)*x*x*x+(-13.5)*x*x*x*y+(22.5)*x*x*y+(22.5)*x*x*z+(-13.5)*x*x*x*z+(13.5)*x*x*x*y*z+(-22.5)*x*x*y*z;
        S(10,i)=(-4.5)*x+(4.5)*x*y+(4.5)*x*z+(-4.5)*x*y*z+(18.)*x*x+(-13.5)*x*x*x+(13.5)*x*x*x*y+(-18.)*x*x*y+(-18.)*x*x*z+(13.5)*x*x*x*z+(-13.5)*x*x*x*y*z+(18.)*x*x*y*z;
        S(11,i)=(9.)*x*y+(-9.)*x*y*z+(-22.5)*x*y*y+(13.5)*x*y*y*y+(22.5)*x*y*y*z+(-13.5)*x*y*y*y*z;
        S(12,i)=(-4.5)*x*y+(4.5)*x*y*z+(18.)*x*y*y+(-13.5)*x*y*y*y+(-18.)*x*y*y*z+(13.5)*x*y*y*y*z;
        S(13,i)=(-4.5)*x*y+(4.5)*x*y*z+(-13.5)*x*x*x*y+(18.)*x*x*y+(13.5)*x*x*x*y*z+(-18.)*x*x*y*z;
        S(14,i)=(9.)*x*y+(-9.)*x*y*z+(13.5)*x*x*x*y+(-22.5)*x*x*y+(-13.5)*x*x*x*y*z+(22.5)*x*x*y*z;
        S(15,i)=(4.5)*x*y+(-4.5)*y+(-4.5)*x*y*z+(4.5)*y*z+(-18.)*x*y*y+(13.5)*x*y*y*y+(-13.5)*y*y*y+(18.)*y*y+(18.)*x*y*y*z+(-13.5)*x*y*y*y*z+(.135e2)*y*y*y*z+(-18.)*y*y*z;
        S(16,i)=(-9.)*x*y+(9.)*y+(9.)*x*y*z+(-9.)*y*z+(.225e2)*x*y*y+(-13.5)*x*y*y*y+(13.5)*y*y*y+(-22.5)*y*y+(-22.5)*x*y*y*z+(13.5)*x*y*y*y*z+(-13.5)*y*y*y*z+(22.5)*y*y*z;
        S(17,i)=(9.)*z+(-9.)*x*z+(9.)*x*y*z+(-9.)*y*z+(-.225e2)*z*z+(22.5)*x*z*z+(-22.5)*x*y*z*z+(22.5)*y*z*z+(13.5)*z*z*z+(-13.5)*x*z*z*z+(13.5)*x*y*z*z*z+(-13.5)*y*z*z*z;
        S(18,i)=(9.)*x*z+(-9.)*x*y*z+(-22.5)*x*z*z+(22.5)*x*y*z*z+(13.5)*x*z*z*z+(-13.5)*x*y*z*z*z;
        S(19,i)=(9.)*x*y*z+(-22.5)*x*y*z*z+(13.5)*x*y*z*z*z;
        S(20,i)=(-9.)*x*y*z+(9.)*y*z+(22.5)*x*y*z*z+(-22.5)*y*z*z+(-13.5)*x*y*z*z*z+(13.5)*y*z*z*z;
        S(21,i)=(-4.5)*z+(4.5)*x*z+(-4.5)*x*y*z+(4.5)*y*z+(18.)*z*z+(-18.)*x*z*z+(18.)*x*y*z*z+(-18.)*y*z*z+(-13.5)*z*z*z+(13.5)*x*z*z*z+(-13.5)*x*y*z*z*z+(13.5)*y*z*z*z;
        S(22,i)=(-4.5)*x*z+(4.5)*x*y*z+(18.)*x*z*z+(-18.)*x*y*z*z+(-13.5)*x*z*z*z+(13.5)*x*y*z*z*z;
        S(23,i)=(-4.5)*x*y*z+(18.)*x*y*z*z+(-13.5)*x*y*z*z*z;
        S(24,i)=(4.5)*x*y*z+(-4.5)*y*z+(-18.)*x*y*z*z+(18.)*y*z*z+(13.5)*x*y*z*z*z+(-13.5)*y*z*z*z;
        S(25,i)=(9.)*x*z+(-9.)*x*y*z+(-22.5)*x*x*z+(13.5)*x*x*x*z+(-13.5)*x*x*x*y*z+(22.5)*x*x*y*z;
        S(26,i)=(-4.5)*x*z+(4.5)*x*y*z+(18.)*x*x*z+(-13.5)*x*x*x*z+(13.5)*x*x*x*y*z+(-18.)*x*x*y*z;
        S(27,i)=(9.)*x*y*z+(-22.5)*x*y*y*z+(13.5)*x*y*y*y*z;
        S(28,i)=(-4.5)*x*y*z+(18.)*x*y*y*z+(-13.5)*x*y*y*y*z;
        S(29,i)=(-4.5)*x*y*z+(-13.5)*x*x*x*y*z+(18.)*x*x*y*z;
        S(30,i)=(9.)*x*y*z+(13.5)*x*x*x*y*z+(-22.5)*x*x*y*z;
        S(31,i)=(4.5)*x*y*z+(-4.5)*y*z+(-18.)*x*y*y*z+(13.5)*x*y*y*y*z+(-13.5)*y*y*y*z+(18.)*y*y*z;
        S(32,i)=(-9.)*x*y*z+(9.)*y*z+(22.5)*x*y*y*z+(-13.5)*x*y*y*y*z+(13.5)*y*y*y*z+(-22.5)*y*y*z;
        DSDV(1, 1,i)=(-5.5)+(10.)*y+(10.)*z+(-0.145e2)*y*z+(18.)*x+(-13.5)*x*x+(-9.)*y*y+(4.5)*y*y*y+(13.5)*x*x*y+(-18.)*x*y+(-9.)*z*z+(9.)*y*z*z+(4.5)*z*z*z+(-4.5)*y*z*z*z+(-18.)*x*z+(13.5)*x*x*z+(9.)*y*y*z+(-4.5)*y*y*y*z+(-13.5)*x*x*y*z+(.18e2)*x*y*z;
        DSDV(2, 1,i)=(10.)+(-5.5)*y+(-5.5)*z+(10.)*y*z+(-9.)*x+(13.5)*x*x+(9.)*y*y+(-4.5)*y*y*y+(-13.5)*x*x*y+(9.)*x*y+(9.)*z*z+(-9.)*y*z*z+(-4.5)*z*z*z+(.45e1)*y*z*z*z+(9.)*x*z+(-13.5)*x*x*z+(-9.)*y*y*z+(4.5)*y*y*y*z+(13.5)*x*x*y*z+(-.9e1)*x*y*z;
        DSDV(3, 1,i)=(10.)*y+(-5.5)*y*z+(-4.5)*y*y+(4.5)*y*y*y+(13.5)*x*x*y+(-9.)*x*y+(9.)*y*z*z+(-4.5)*y*z*z*z+(4.5)*y*y*z+(-4.5)*y*y*y*z+(-13.5)*x*x*y*z+(.9e1)*x*y*z;
        DSDV(4, 1,i)=(-5.5)*y+(10.)*y*z+(4.5)*y*y+(-4.5)*y*y*y+(-13.5)*x*x*y+(18.)*x*y+(-9.)*y*z*z+(4.5)*y*z*z*z+(-4.5)*y*y*z+(4.5)*y*y*y*z+(13.5)*x*x*y*z+(-18.)*x*y*z;
        DSDV(5, 1,i)=(-5.5)*z+(10.)*y*z+(4.5)*z*z+(-4.5)*y*z*z+(-.45e1)*z*z*z+(4.5)*y*z*z*z+(18.)*x*z+(-13.5)*x*x*z+(-9.)*y*y*z+(4.5)*y*y*y*z+(13.5)*x*x*y*z+(-.18e2)*x*y*z;
        DSDV(6, 1,i)=(10.)*z+(-5.5)*y*z+(-4.5)*z*z+(4.5)*y*z*z+(.45e1)*z*z*z+(-4.5)*y*z*z*z+(-9.)*x*z+(13.5)*x*x*z+(9.)*y*y*z+(-4.5)*y*y*y*z+(-13.5)*x*x*y*z+(9.)*x*y*z;
        DSDV(7, 1,i)=(10.)*y*z+(-4.5)*y*z*z+(4.5)*y*z*z*z+(-4.5)*y*y*z+(4.5)*y*y*y*z+(13.5)*x*x*y*z+(-9.)*x*y*z;
        DSDV(8, 1,i)=(-5.5)*y*z+(4.5)*y*z*z+(-4.5)*y*z*z*z+(4.5)*y*y*z+(-4.5)*y*y*y*z+(-13.5)*x*x*y*z+(18.)*x*y*z;
        DSDV(9, 1,i)=(9.)+(-9.)*y+(-9.)*z+(9.)*y*z+(-45.)*x+(0.405e2)*x*x+(-0.405e2)*x*x*y+(45.)*x*y+(45.)*x*z+(-0.405e2)*x*x*z+(0.405e2)*x*x*y*z+(-45.)*x*y*z;
        DSDV(10, 1,i)=(-4.5)+(4.5)*y+(4.5)*z+(-4.5)*y*z+(36.)*x+(-.405e2)*x*x+(0.405e2)*x*x*y+(-36.)*x*y+(-36.)*x*z+(0.405e2)*x*x*z+(-0.405e2)*x*x*y*z+(36.)*x*y*z;
        DSDV(11, 1,i)=(9.)*y+(-9.)*y*z+(-22.5)*y*y+(13.5)*y*y*y+(22.5)*y*y*z+(-13.5)*y*y*y*z;
        DSDV(12, 1,i)=(-4.5)*y+(4.5)*y*z+(18.)*y*y+(-13.5)*y*y*y+(-18.)*y*y*z+(13.5)*y*y*y*z;
        DSDV(13, 1,i)=(-4.5)*y+(4.5)*y*z+(-0.405e2)*x*x*y+(36.)*x*y+(0.405e2)*x*x*y*z+(-36.)*x*y*z;
        DSDV(14, 1,i)=(9.)*y+(-9.)*y*z+(0.405e2)*x*x*y+(-45.)*x*y+(-0.405e2)*x*x*y*z+(45.)*x*y*z;
        DSDV(15, 1,i)=(4.5)*y+(-4.5)*y*z+(-18.)*y*y+(13.5)*y*y*y+(18.)*y*y*z+(-13.5)*y*y*y*z;
        DSDV(16, 1,i)=(-9.)*y+(9.)*y*z+(22.5)*y*y+(-13.5)*y*y*y+(-22.5)*y*y*z+(13.5)*y*y*y*z;
        DSDV(17, 1,i)=(-9.)*z+(9.)*y*z+(22.5)*z*z+(-22.5)*y*z*z+(-13.5)*z*z*z+(13.5)*y*z*z*z;
        DSDV(18, 1,i)=(9.)*z+(-9.)*y*z+(-22.5)*z*z+(22.5)*y*z*z+(.135e2)*z*z*z+(-13.5)*y*z*z*z;
        DSDV(19, 1,i)=(9.)*y*z+(-22.5)*y*z*z+(13.5)*y*z*z*z;
        DSDV(20, 1,i)=(-9.)*y*z+(22.5)*y*z*z+(-13.5)*y*z*z*z;
        DSDV(21, 1,i)=(4.5)*z+(-4.5)*y*z+(-18.)*z*z+(18.)*y*z*z+(.135e2)*z*z*z+(-13.5)*y*z*z*z;
        DSDV(22, 1,i)=(-4.5)*z+(4.5)*y*z+(18.)*z*z+(-18.)*y*z*z+(-13.5)*z*z*z+(13.5)*y*z*z*z;
        DSDV(23, 1,i)=(-4.5)*y*z+(18.)*y*z*z+(-13.5)*y*z*z*z;
        DSDV(24, 1,i)=(4.5)*y*z+(-18.)*y*z*z+(13.5)*y*z*z*z;
        DSDV(25, 1,i)=(9.)*z+(-9.)*y*z+(-45.)*x*z+(0.405e2)*x*x*z+(-0.405e2)*x*x*y*z+(45.)*x*y*z;
        DSDV(26, 1,i)=(-4.5)*z+(4.5)*y*z+(36.)*x*z+(-0.405e2)*x*x*z+(0.405e2)*x*x*y*z+(-36.)*x*y*z;
        DSDV(27, 1,i)=(9.)*y*z+(-22.5)*y*y*z+(13.5)*y*y*y*z;
        DSDV(28, 1,i)=(-4.5)*y*z+(18.)*y*y*z+(-13.5)*y*y*y*z;
        DSDV(29, 1,i)=(-4.5)*y*z+(-0.405e2)*x*x*y*z+(36.)*x*y*z;
        DSDV(30, 1,i)=(9.)*y*z+(0.405e2)*x*x*y*z+(-45.)*x*y*z;
        DSDV(31, 1,i)=(4.5)*y*z+(-18.)*y*y*z+(13.5)*y*y*y*z;
        DSDV(32, 1,i)=(-9.)*y*z+(22.5)*y*y*z+(-13.5)*y*y*y*z;
        DSDV(1, 2,i)=(10.)*x+(-5.5)+(-0.145e2)*x*z+(10.)*z+(-18.)*x*y+(13.5)*x*y*y+(4.5)*x*x*x+(-9.)*x*x+(-13.5)*y*y+(18.)*y+(9.)*x*z*z+(-9.)*z*z+(-4.5)*x*z*z*z+(4.5)*z*z*z+(18.)*x*y*z+(-13.5)*x*y*y*z+(-4.5)*x*x*x*z+(9.)*x*x*z+(.135e2)*y*y*z+(-18.)*y*z;
        DSDV(2, 2,i)=(-5.5)*x+(10.)*x*z+(18.)*x*y+(-13.5)*x*y*y+(-4.5)*x*x*x+(4.5)*x*x+(-9.)*x*z*z+(4.5)*x*z*z*z+(-18.)*x*y*z+(13.5)*x*y*y*z+(4.5)*x*x*x*z+(-4.5)*x*x*z;
        DSDV(3, 2,i)=(10.)*x+(-5.5)*x*z+(-9.)*x*y+(13.5)*x*y*y+(4.5)*x*x*x+(-4.5)*x*x+(9.)*x*z*z+(-4.5)*x*z*z*z+(9.)*x*y*z+(-13.5)*x*y*y*z+(-4.5)*x*x*x*z+(4.5)*x*x*z;
        DSDV(4, 2,i)=(-5.5)*x+(10.)+(10.)*x*z+(-5.5)*z+(9.)*x*y+(-13.5)*x*y*y+(-4.5)*x*x*x+(9.)*x*x+(13.5)*y*y+(-9.)*y+(-9.)*x*z*z+(9.)*z*z+(4.5)*x*z*z*z+(-4.5)*z*z*z+(-9.)*x*y*z+(13.5)*x*y*y*z+(4.5)*x*x*x*z+(-9.)*x*x*z+(-.135e2)*y*y*z+(9.)*y*z;
        DSDV(5, 2,i)=(10.)*x*z+(-5.5)*z+(-4.5)*x*z*z+(4.5)*z*z+(.45e1)*x*z*z*z+(-4.5)*z*z*z+(-18.)*x*y*z+(13.5)*x*y*y*z+(4.5)*x*x*x*z+(-9.)*x*x*z+(-.135e2)*y*y*z+(18.)*y*z;
        DSDV(6, 2,i)=(-5.5)*x*z+(4.5)*x*z*z+(-4.5)*x*z*z*z+(18.)*x*y*z+(-13.5)*x*y*y*z+(-4.5)*x*x*x*z+(4.5)*x*x*z;
        DSDV(7, 2,i)=(10.)*x*z+(-4.5)*x*z*z+(4.5)*x*z*z*z+(-9.)*x*y*z+(13.5)*x*y*y*z+(4.5)*x*x*x*z+(-4.5)*x*x*z;
        DSDV(8, 2,i)=(-5.5)*x*z+(10.)*z+(4.5)*x*z*z+(-4.5)*z*z+(-.45e1)*x*z*z*z+(4.5)*z*z*z+(9.)*x*y*z+(-13.5)*x*y*y*z+(-4.5)*x*x*x*z+(9.)*x*x*z+(.135e2)*y*y*z+(-9.)*y*z;
        DSDV(9, 2,i)=(-9.)*x+(9.)*x*z+(-13.5)*x*x*x+(22.5)*x*x+(13.5)*x*x*x*z+(-22.5)*x*x*z;
        DSDV(10, 2,i)=(4.5)*x+(-4.5)*x*z+(13.5)*x*x*x+(-18.)*x*x+(-13.5)*x*x*x*z+(18.)*x*x*z;
        DSDV(11, 2,i)=(9.)*x+(-9.)*x*z+(-45.)*x*y+(0.405e2)*x*y*y+(45.)*x*y*z+(-0.405e2)*x*y*y*z;
        DSDV(12, 2,i)=(-4.5)*x+(4.5)*x*z+(36.)*x*y+(-0.405e2)*x*y*y+(-36.)*x*y*z+(0.405e2)*x*y*y*z;
        DSDV(13, 2,i)=(-4.5)*x+(4.5)*x*z+(-13.5)*x*x*x+(18.)*x*x+(13.5)*x*x*x*z+(-18.)*x*x*z;
        DSDV(14, 2,i)=(9.)*x+(-9.)*x*z+(13.5)*x*x*x+(-22.5)*x*x+(-13.5)*x*x*x*z+(22.5)*x*x*z;
        DSDV(15, 2,i)=(4.5)*x+(-4.5)+(-4.5)*x*z+(4.5)*z+(-36.)*x*y+(0.405e2)*x*y*y+(-0.405e2)*y*y+(36.)*y+(36.)*x*y*z+(-0.405e2)*x*y*y*z+(0.405e2)*y*y*z+(-36.)*y*z;
        DSDV(16, 2,i)=(-9.)*x+(9.)+(9.)*x*z+(-9.)*z+(45.)*x*y+(-0.405e2)*x*y*y+(0.405e2)*y*y+(-45.)*y+(-45.)*x*y*z+(0.405e2)*x*y*y*z+(-0.405e2)*y*y*z+(45.)*y*z;
        DSDV(17, 2,i)=(9.)*x*z+(-9.)*z+(-22.5)*x*z*z+(22.5)*z*z+(13.5)*x*z*z*z+(-13.5)*z*z*z;
        DSDV(18, 2,i)=(-9.)*x*z+(22.5)*x*z*z+(-13.5)*x*z*z*z;
        DSDV(19, 2,i)=(9.)*x*z+(-22.5)*x*z*z+(13.5)*x*z*z*z;
        DSDV(20, 2,i)=(-9.)*x*z+(9.)*z+(22.5)*x*z*z+(-22.5)*z*z+(-13.5)*x*z*z*z+(13.5)*z*z*z;
        DSDV(21, 2,i)=(-4.5)*x*z+(4.5)*z+(18.)*x*z*z+(-18.)*z*z+(-13.5)*x*z*z*z+(13.5)*z*z*z;
        DSDV(22, 2,i)=(4.5)*x*z+(-18.)*x*z*z+(13.5)*x*z*z*z;
        DSDV(23, 2,i)=(-4.5)*x*z+(18.)*x*z*z+(-13.5)*x*z*z*z;
        DSDV(24, 2,i)=(4.5)*x*z+(-4.5)*z+(-18.)*x*z*z+(18.)*z*z+(13.5)*x*z*z*z+(-13.5)*z*z*z;
        DSDV(25, 2,i)=(-9.)*x*z+(-13.5)*x*x*x*z+(22.5)*x*x*z;
        DSDV(26, 2,i)=(4.5)*x*z+(13.5)*x*x*x*z+(-18.)*x*x*z;
        DSDV(27, 2,i)=(9.)*x*z+(-45.)*x*y*z+(0.405e2)*x*y*y*z;
        DSDV(28, 2,i)=(-4.5)*x*z+(36.)*x*y*z+(-0.405e2)*x*y*y*z;
        DSDV(29, 2,i)=(-4.5)*x*z+(-13.5)*x*x*x*z+(18.)*x*x*z;
        DSDV(30, 2,i)=(9.)*x*z+(13.5)*x*x*x*z+(-22.5)*x*x*z;
        DSDV(31, 2,i)=(4.5)*x*z+(-4.5)*z+(-36.)*x*y*z+(0.405e2)*x*y*y*z+(-0.405e2)*y*y*z+(36.)*y*z;
        DSDV(32, 2,i)=(-9.)*x*z+(9.)*z+(45.)*x*y*z+(-0.405e2)*x*y*y*z+(0.405e2)*y*y*z+(-45.)*y*z;
        DSDV(1, 3,i)=(-5.5)+(10.)*x+(-0.145e2)*x*y+(10.)*y+(18.)*z+(-18.)*x*z+(18.)*x*y*z+(-18.)*y*z+(-13.5)*z*z+(13.5)*x*z*z+(-13.5)*x*y*z*z+(13.5)*y*z*z+(-9.)*x*x+(4.5)*x*x*x+(9.)*x*y*y+(-4.5)*x*y*y*y+(-4.5)*x*x*x*y+(9.)*x*x*y+(4.5)*y*y*y+(-9.)*y*y;
        DSDV(2, 3,i)=(-5.5)*x+(10.)*x*y+(18.)*x*z+(-18.)*x*y*z+(-13.5)*x*z*z+(13.5)*x*y*z*z+(4.5)*x*x+(-4.5)*x*x*x+(-9.)*x*y*y+(4.5)*x*y*y*y+(4.5)*x*x*x*y+(-4.5)*x*x*y;
        DSDV(3, 3,i)=(-5.5)*x*y+(18.)*x*y*z+(-13.5)*x*y*z*z+(4.5)*x*y*y+(-4.5)*x*y*y*y+(-4.5)*x*x*x*y+(4.5)*x*x*y;
        DSDV(4, 3,i)=(10.)*x*y+(-5.5)*y+(-18.)*x*y*z+(18.)*y*z+(13.5)*x*y*z*z+(-13.5)*y*z*z+(-4.5)*x*y*y+(4.5)*x*y*y*y+(4.5)*x*x*x*y+(-9.)*x*x*y+(-4.5)*y*y*y+(4.5)*y*y;
        DSDV(5, 3,i)=(10.)+(-5.5)*x+(10.)*x*y+(-5.5)*y+(-9.)*z+(9.)*x*z+(-9.)*x*y*z+(9.)*y*z+(13.5)*z*z+(-13.5)*x*z*z+(13.5)*x*y*z*z+(-13.5)*y*z*z+(9.)*x*x+(-4.5)*x*x*x+(-9.)*x*y*y+(4.5)*x*y*y*y+(4.5)*x*x*x*y+(-9.)*x*x*y+(-4.5)*y*y*y+(9.)*y*y;
        DSDV(6, 3,i)=(10.)*x+(-5.5)*x*y+(-9.)*x*z+(9.)*x*y*z+(13.5)*x*z*z+(-13.5)*x*y*z*z+(-4.5)*x*x+(4.5)*x*x*x+(9.)*x*y*y+(-4.5)*x*y*y*y+(-4.5)*x*x*x*y+(4.5)*x*x*y;
        DSDV(7, 3,i)=(10.)*x*y+(-9.)*x*y*z+(13.5)*x*y*z*z+(-4.5)*x*y*y+(4.5)*x*y*y*y+(4.5)*x*x*x*y+(-4.5)*x*x*y;
        DSDV(8, 3,i)=(-5.5)*x*y+(10.)*y+(9.)*x*y*z+(-9.)*y*z+(-13.5)*x*y*z*z+(13.5)*y*z*z+(4.5)*x*y*y+(-4.5)*x*y*y*y+(-4.5)*x*x*x*y+(9.)*x*x*y+(4.5)*y*y*y+(-4.5)*y*y;
        DSDV(9, 3,i)=(-9.)*x+(9.)*x*y+(22.5)*x*x+(-13.5)*x*x*x+(13.5)*x*x*x*y+(-22.5)*x*x*y;
        DSDV(10, 3,i)=(4.5)*x+(-4.5)*x*y+(-18.)*x*x+(13.5)*x*x*x+(-13.5)*x*x*x*y+(18.)*x*x*y;
        DSDV(11, 3,i)=(-9.)*x*y+(22.5)*x*y*y+(-13.5)*x*y*y*y;
        DSDV(12, 3,i)=(4.5)*x*y+(-18.)*x*y*y+(13.5)*x*y*y*y;
        DSDV(13, 3,i)=(4.5)*x*y+(13.5)*x*x*x*y+(-18.)*x*x*y;
        DSDV(14, 3,i)=(-9.)*x*y+(-13.5)*x*x*x*y+(22.5)*x*x*y;
        DSDV(15, 3,i)=(-4.5)*x*y+(4.5)*y+(18.)*x*y*y+(-13.5)*x*y*y*y+(13.5)*y*y*y+(-18.)*y*y;
        DSDV(16, 3,i)=(9.)*x*y+(-9.)*y+(-22.5)*x*y*y+(13.5)*x*y*y*y+(-13.5)*y*y*y+(22.5)*y*y;
        DSDV(17, 3,i)=(9.)+(-9.)*x+(9.)*x*y+(-9.)*y+(-45.)*z+(45.)*x*z+(-45.)*x*y*z+(45.)*y*z+(0.405e2)*z*z+(-0.405e2)*x*z*z+(0.405e2)*x*y*z*z+(-0.405e2)*y*z*z;
        DSDV(18, 3,i)=(9.)*x+(-9.)*x*y+(-45.)*x*z+(45.)*x*y*z+(0.405e2)*x*z*z+(-0.405e2)*x*y*z*z;
        DSDV(19, 3,i)=(9.)*x*y+(-45.)*x*y*z+(0.405e2)*x*y*z*z;
        DSDV(20, 3,i)=(-9.)*x*y+(9.)*y+(45.)*x*y*z+(-45.)*y*z+(-0.405e2)*x*y*z*z+(0.405e2)*y*z*z;
        DSDV(21, 3,i)=(-4.5)+(4.5)*x+(-4.5)*x*y+(4.5)*y+(36.)*z+(-36.)*x*z+(36.)*x*y*z+(-36.)*y*z+(-0.405e2)*z*z+(0.405e2)*x*z*z+(-0.405e2)*x*y*z*z+(0.405e2)*y*z*z;
        DSDV(22, 3,i)=(-4.5)*x+(4.5)*x*y+(36.)*x*z+(-36.)*x*y*z+(-0.405e2)*x*z*z+(0.405e2)*x*y*z*z;
        DSDV(23, 3,i)=(-4.5)*x*y+(36.)*x*y*z+(-0.405e2)*x*y*z*z;
        DSDV(24, 3,i)=(4.5)*x*y+(-4.5)*y+(-36.)*x*y*z+(36.)*y*z+(0.405e2)*x*y*z*z+(-0.405e2)*y*z*z;
        DSDV(25, 3,i)=(9.)*x+(-9.)*x*y+(-22.5)*x*x+(13.5)*x*x*x+(-13.5)*x*x*x*y+(22.5)*x*x*y;
        DSDV(26, 3,i)=(-4.5)*x+(4.5)*x*y+(18.)*x*x+(-13.5)*x*x*x+(13.5)*x*x*x*y+(-18.)*x*x*y;
        DSDV(27, 3,i)=(9.)*x*y+(-22.5)*x*y*y+(13.5)*x*y*y*y;
        DSDV(28, 3,i)=(-4.5)*x*y+(18.)*x*y*y+(-13.5)*x*y*y*y;
        DSDV(29, 3,i)=(-4.5)*x*y+(-13.5)*x*x*x*y+(18.)*x*x*y;
        DSDV(30, 3,i)=(9.)*x*y+(13.5)*x*x*x*y+(-22.5)*x*x*y;
        DSDV(31, 3,i)=(4.5)*x*y+(-4.5)*y+(-18.)*x*y*y+(13.5)*x*y*y*y+(-13.5)*y*y*y+(18.)*y*y;
        DSDV(32, 3,i)=(-9.)*x*y+(9.)*y+(22.5)*x*y*y+(-13.5)*x*y*y*y+(13.5)*y*y*y+(-22.5)*y*y;
    }
#undef NUMSHAPES
#undef DIM
}

#undef V
#undef S
#undef DSDV

} // namespace finley

