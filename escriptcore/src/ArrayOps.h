
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#ifndef __ESCRIPT_LOCALOPS_H__
#define __ESCRIPT_LOCALOPS_H__

#include "DataTypes.h"
#include "DataException.h"
#include <iostream>
#include <cmath>
#include <complex>
#ifdef _WIN32
# include <algorithm>
#endif // _WIN32

#ifdef ESYS_USE_BOOST_ACOS
#include <boost/math/complex/acos.hpp>	// std::acos for complex on OSX (elcapitan) is wrong
#endif

#ifndef M_PI
#   define M_PI           3.14159265358979323846  /* pi */
#endif


/**
\file LocalOps.h 
\brief Describes binary operations performed on double*.

For operations on DataAbstract see BinaryOp.h.
For operations on DataVector see DataMaths.h.
*/

#include "ES_optype.h"

namespace escript {



bool always_real(escript::ES_optype operation);


/**
   \brief
   Return the maximum value of the two given values.
*/
struct FMax 
{
  inline DataTypes::real_t operator()(DataTypes::real_t x, DataTypes::real_t y) const
  {
    return std::max(x,y);
  }
  typedef DataTypes::real_t first_argument_type;
  typedef DataTypes::real_t second_argument_type;
  typedef DataTypes::real_t result_type;
};

/**
   \brief
   Return the minimum value of the two given values.
*/
struct FMin
{
  inline DataTypes::real_t operator()(DataTypes::real_t x, DataTypes::real_t y) const
  {
    return std::min(x,y);
  }
  typedef DataTypes::real_t first_argument_type;
  typedef DataTypes::real_t second_argument_type;
  typedef DataTypes::real_t result_type;  
};

/**
   \brief
   Return the absolute maximum value of the two given values.
*/
template<typename T>
struct AbsMax 
{
  inline DataTypes::real_t operator()(T x, T y) const
  {
    return std::max(std::abs(x),std::abs(y));
  }
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef DataTypes::real_t result_type;
};


inline
DataTypes::real_t
fsign(DataTypes::real_t x)
{
  if (x == 0) {
    return 0;
  } else {
    return x/fabs(x);
  }
}

/**
\brief acts as a wrapper to isnan.
\warning if compiler does not support FP_NAN this function will always return false.
*/
inline
bool nancheck(DataTypes::real_t d)
{
                // Q: so why not just test d!=d?
                // A: Coz it doesn't always work [I've checked].
                // One theory is that the optimizer skips the test.
    return std::isnan(d);       // isNan should be a function in C++ land
}

// I'm not sure if there is agreement about what a complex nan would be
// so, in this case I've just checked both components
inline
bool nancheck(DataTypes::cplx_t d)
{
                // Q: so why not just test d!=d?
                // A: Coz it doesn't always work [I've checked].
                // One theory is that the optimizer skips the test.
    return std::isnan( real(d) ) || std::isnan( imag(d) );       // isNan should be a function in C++ land
}


/**
\brief returns a NaN.
\warning Should probably only used where you know you can test for NaNs
*/
inline
DataTypes::real_t makeNaN()
{
#ifdef nan
    return nan("");
#else
    return sqrt(-1.);
#endif

}


/**
   \brief
   solves a 1x1 eigenvalue A*V=ev*V problem

   \param A00 Input - A_00
   \param ev0 Output - eigenvalue
*/
inline
void eigenvalues1(const DataTypes::real_t A00,DataTypes::real_t* ev0) {
   *ev0=A00;
}

inline
void eigenvalues1(const DataTypes::cplx_t A00,DataTypes::cplx_t* ev0) {
   *ev0=A00;
}


/**
   \brief
   solves a 2x2 eigenvalue A*V=ev*V problem for symmetric A

   \param A00 Input - A_00
   \param A01 Input - A_01
   \param A11 Input - A_11
   \param ev0 Output - smallest eigenvalue
   \param ev1 Output - largest eigenvalue
*/
template <class T>
inline
void eigenvalues2(const T A00,const T A01,const T A11,
                 T* ev0, T* ev1) {
      const T trA=(A00+A11)/2.;
      const T A_00=A00-trA;
      const T A_11=A11-trA;
      const T s=sqrt(A01*A01-A_00*A_11);
      *ev0=trA-s;
      *ev1=trA+s;
}
/**
   \brief
   solves a 3x3 eigenvalue A*V=ev*V problem for symmetric A

   \param A00 Input - A_00
   \param A01 Input - A_01
   \param A02 Input - A_02
   \param A11 Input - A_11
   \param A12 Input - A_12
   \param A22 Input - A_22
   \param ev0 Output - smallest eigenvalue
   \param ev1 Output - eigenvalue
   \param ev2 Output - largest eigenvalue
*/
inline
void eigenvalues3(const DataTypes::real_t A00, const DataTypes::real_t A01, const DataTypes::real_t A02,
                                   const DataTypes::real_t A11, const DataTypes::real_t A12,
                                                     const DataTypes::real_t A22,
                 DataTypes::real_t* ev0, DataTypes::real_t* ev1,DataTypes::real_t* ev2) {

      const DataTypes::real_t trA=(A00+A11+A22)/3.;
      const DataTypes::real_t A_00=A00-trA;
      const DataTypes::real_t A_11=A11-trA;
      const DataTypes::real_t A_22=A22-trA;
      const DataTypes::real_t A01_2=A01*A01;
      const DataTypes::real_t A02_2=A02*A02;
      const DataTypes::real_t A12_2=A12*A12;
      const DataTypes::real_t p=A02_2+A12_2+A01_2+(A_00*A_00+A_11*A_11+A_22*A_22)/2.;
      if (p<=0.) {
         *ev2=trA;
         *ev1=trA;
         *ev0=trA;

      } else {
         const DataTypes::real_t q=(A02_2*A_11+A12_2*A_00+A01_2*A_22)-(A_00*A_11*A_22+2*A01*A12*A02);
         const DataTypes::real_t sq_p=sqrt(p/3.);
         DataTypes::real_t z=-q/(2*pow(sq_p,3));
         if (z<-1.) {
            z=-1.;
         } else if (z>1.) {
            z=1.;
         }
         const DataTypes::real_t alpha_3=acos(z)/3.;
         *ev2=trA+2.*sq_p*cos(alpha_3);
         *ev1=trA-2.*sq_p*cos(alpha_3+M_PI/3.);
         *ev0=trA-2.*sq_p*cos(alpha_3-M_PI/3.);
      }
}
/**
   \brief
   solves a 1x1 eigenvalue A*V=ev*V problem for symmetric A

   \param A00 Input - A_00
   \param ev0 Output - eigenvalue
   \param V00 Output - eigenvector
   \param tol Input - tolerance to identify to eigenvalues
*/
inline
void  eigenvalues_and_eigenvectors1(const DataTypes::real_t A00,DataTypes::real_t* ev0,DataTypes::real_t* V00,const DataTypes::real_t tol)
{
      eigenvalues1(A00,ev0);
      *V00=1.;
      return;
}
/**
   \brief
   returns a non-zero vector in the kernel of [[A00,A01],[A01,A11]] assuming that the kernel dimension is at least 1.

   \param A00 Input - matrix component
   \param A10 Input - matrix component
   \param A01 Input - matrix component
   \param A11 Input - matrix component
   \param V0 Output - vector component
   \param V1 Output - vector component
*/
inline
void  vectorInKernel2(const DataTypes::real_t A00,const DataTypes::real_t A10,const DataTypes::real_t A01,const DataTypes::real_t A11,
                      DataTypes::real_t* V0, DataTypes::real_t*V1)
{
      DataTypes::real_t absA00=fabs(A00);
      DataTypes::real_t absA10=fabs(A10);
      DataTypes::real_t absA01=fabs(A01);
      DataTypes::real_t absA11=fabs(A11);
      DataTypes::real_t m=absA11>absA10 ? absA11 : absA10;
      if (absA00>m || absA01>m) {
         *V0=-A01;
         *V1=A00;
      } else {
         if (m<=0) {
           *V0=1.;
           *V1=0.;
         } else {
           *V0=A11;
           *V1=-A10;
         }
     }
}
/**
   \brief
   returns a non-zero vector in the kernel of [[A00,A01,A02],[A10,A11,A12],[A20,A21,A22]]
   assuming that the kernel dimension is at least 1 and A00 is non zero.

   \param A00 Input - matrix component
   \param A10 Input - matrix component
   \param A20 Input - matrix component
   \param A01 Input - matrix component
   \param A11 Input - matrix component
   \param A21 Input - matrix component
   \param A02 Input - matrix component
   \param A12 Input - matrix component
   \param A22 Input - matrix component
   \param V0 Output - vector component
   \param V1 Output - vector component
   \param V2 Output - vector component
*/
inline
void  vectorInKernel3__nonZeroA00(const DataTypes::real_t A00,const DataTypes::real_t A10,const DataTypes::real_t A20,
                                const DataTypes::real_t A01,const DataTypes::real_t A11,const DataTypes::real_t A21,
                                const DataTypes::real_t A02,const DataTypes::real_t A12,const DataTypes::real_t A22,
                                DataTypes::real_t* V0,DataTypes::real_t* V1,DataTypes::real_t* V2)
{
    DataTypes::real_t TEMP0,TEMP1;
    const DataTypes::real_t I00=1./A00;
    const DataTypes::real_t IA10=I00*A10;
    const DataTypes::real_t IA20=I00*A20;
    vectorInKernel2(A11-IA10*A01,A12-IA10*A02,
                    A21-IA20*A01,A22-IA20*A02,&TEMP0,&TEMP1);
    *V0=-(A10*TEMP0+A20*TEMP1);
    *V1=A00*TEMP0;
    *V2=A00*TEMP1;
}

/**
   \brief
   solves a 2x2 eigenvalue A*V=ev*V problem for symmetric A. Eigenvectors are
   ordered by increasing value and eigen vectors are normalizeVector3d such that
   length is zero and first non-zero component is positive.

   \param A00 Input - A_00
   \param A01 Input - A_01
   \param A11 Input - A_11
   \param ev0 Output - smallest eigenvalue
   \param ev1 Output - eigenvalue
   \param V00 Output - eigenvector componenent coresponding to ev0
   \param V10 Output - eigenvector componenent coresponding to ev0
   \param V01 Output - eigenvector componenent coresponding to ev1
   \param V11 Output - eigenvector componenent coresponding to ev1
   \param tol Input - tolerance to identify to eigenvalues
*/
inline
void  eigenvalues_and_eigenvectors2(const DataTypes::real_t A00,const DataTypes::real_t A01,const DataTypes::real_t A11,
                                    DataTypes::real_t* ev0, DataTypes::real_t* ev1,
                                    DataTypes::real_t* V00, DataTypes::real_t* V10, DataTypes::real_t* V01, DataTypes::real_t* V11,
                                    const DataTypes::real_t tol)
{
     DataTypes::real_t TEMP0,TEMP1;
     eigenvalues2(A00,A01,A11,ev0,ev1);
     const DataTypes::real_t absev0=fabs(*ev0);
     const DataTypes::real_t absev1=fabs(*ev1);
     DataTypes::real_t max_ev=absev0>absev1 ? absev0 : absev1;
     if (fabs((*ev0)-(*ev1))<tol*max_ev) {
        *V00=1.;
        *V10=0.;
        *V01=0.;
        *V11=1.;
     } else {
        vectorInKernel2(A00-(*ev0),A01,A01,A11-(*ev0),&TEMP0,&TEMP1);
        const DataTypes::real_t scale=1./sqrt(TEMP0*TEMP0+TEMP1*TEMP1);
        if (TEMP0<0.) {
            *V00=-TEMP0*scale;
            *V10=-TEMP1*scale;
            if (TEMP1<0.) {
               *V01=  *V10;
               *V11=-(*V00);
            } else {
               *V01=-(*V10);
               *V11= (*V00);
            }
        } else if (TEMP0>0.) {
            *V00=TEMP0*scale;
            *V10=TEMP1*scale;
            if (TEMP1<0.) {
               *V01=-(*V10);
               *V11= (*V00);
            } else {
               *V01= (*V10);
               *V11=-(*V00);
            }
        } else {
           *V00=0.;
           *V10=1;
           *V11=0.;
           *V01=1.;
       }
   }
}
/**
   \brief
   nomalizes a 3-d vector such that length is one and first non-zero component is positive.

   \param V0 - vector componenent
   \param V1 - vector componenent
   \param V2 - vector componenent
*/
inline
void  normalizeVector3(DataTypes::real_t* V0,DataTypes::real_t* V1,DataTypes::real_t* V2)
{
    DataTypes::real_t s;
    if (*V0>0) {
        s=1./sqrt((*V0)*(*V0)+(*V1)*(*V1)+(*V2)*(*V2));
        *V0*=s;
        *V1*=s;
        *V2*=s;
    } else if (*V0<0)  {
        s=-1./sqrt((*V0)*(*V0)+(*V1)*(*V1)+(*V2)*(*V2));
        *V0*=s;
        *V1*=s;
        *V2*=s;
    } else {
        if (*V1>0) {
            s=1./sqrt((*V1)*(*V1)+(*V2)*(*V2));
            *V1*=s;
            *V2*=s;
        } else if (*V1<0)  {
            s=-1./sqrt((*V1)*(*V1)+(*V2)*(*V2));
            *V1*=s;
            *V2*=s;
        } else {
            *V2=1.;
        }
    }
}
/**
   \brief
   solves a 2x2 eigenvalue A*V=ev*V problem for symmetric A. Eigenvectors are
   ordered by increasing value and eigen vectors are normalizeVector3d such that
   length is zero and first non-zero component is positive.

   \param A00 Input - A_00
   \param A01 Input - A_01
   \param A02 Input - A_02
   \param A11 Input - A_11
   \param A12 Input - A_12
   \param A22 Input - A_22
   \param ev0 Output - smallest eigenvalue
   \param ev1 Output - eigenvalue
   \param ev2 Output -
   \param V00 Output - eigenvector componenent coresponding to ev0
   \param V10 Output - eigenvector componenent coresponding to ev0
   \param V20 Output -
   \param V01 Output - eigenvector componenent coresponding to ev1
   \param V11 Output - eigenvector componenent coresponding to ev1
   \param V21 Output -
   \param V02 Output -
   \param V12 Output -
   \param V22 Output -
   \param tol Input - tolerance to identify to eigenvalues
*/
inline
void  eigenvalues_and_eigenvectors3(const DataTypes::real_t A00, const DataTypes::real_t A01, const DataTypes::real_t A02,
                                    const DataTypes::real_t A11, const DataTypes::real_t A12, const DataTypes::real_t A22,
                                    DataTypes::real_t* ev0, DataTypes::real_t* ev1, DataTypes::real_t* ev2,
                                    DataTypes::real_t* V00, DataTypes::real_t* V10, DataTypes::real_t* V20,
                                    DataTypes::real_t* V01, DataTypes::real_t* V11, DataTypes::real_t* V21,
                                    DataTypes::real_t* V02, DataTypes::real_t* V12, DataTypes::real_t* V22,
                                    const DataTypes::real_t tol)
{
      const DataTypes::real_t absA01=fabs(A01);
      const DataTypes::real_t absA02=fabs(A02);
      const DataTypes::real_t m=absA01>absA02 ? absA01 : absA02;
      if (m<=0) {
        DataTypes::real_t TEMP_V00,TEMP_V10,TEMP_V01,TEMP_V11,TEMP_EV0,TEMP_EV1;
        eigenvalues_and_eigenvectors2(A11,A12,A22,
                                      &TEMP_EV0,&TEMP_EV1,
                                      &TEMP_V00,&TEMP_V10,&TEMP_V01,&TEMP_V11,tol);
        if (A00<=TEMP_EV0) {
            *V00=1.;
            *V10=0.;
            *V20=0.;
            *V01=0.;
            *V11=TEMP_V00;
            *V21=TEMP_V10;
            *V02=0.;
            *V12=TEMP_V01;
            *V22=TEMP_V11;
            *ev0=A00;
            *ev1=TEMP_EV0;
            *ev2=TEMP_EV1;
        } else if (A00>TEMP_EV1) {
            *V02=1.;
            *V12=0.;
            *V22=0.;
            *V00=0.;
            *V10=TEMP_V00;
            *V20=TEMP_V10;
            *V01=0.;
            *V11=TEMP_V01;
            *V21=TEMP_V11;
            *ev0=TEMP_EV0;
            *ev1=TEMP_EV1;
            *ev2=A00;
        } else {
            *V01=1.;
            *V11=0.;
            *V21=0.;
            *V00=0.;
            *V10=TEMP_V00;
            *V20=TEMP_V10;
            *V02=0.;
            *V12=TEMP_V01;
            *V22=TEMP_V11;
            *ev0=TEMP_EV0;
            *ev1=A00;
            *ev2=TEMP_EV1;
        }
      } else {
         eigenvalues3(A00,A01,A02,A11,A12,A22,ev0,ev1,ev2);
         const DataTypes::real_t absev0=fabs(*ev0);
         const DataTypes::real_t absev1=fabs(*ev1);
         const DataTypes::real_t absev2=fabs(*ev2);
         DataTypes::real_t max_ev=absev0>absev1 ? absev0 : absev1;
         max_ev=max_ev>absev2 ? max_ev : absev2;
         const DataTypes::real_t d_01=fabs((*ev0)-(*ev1));
         const DataTypes::real_t d_12=fabs((*ev1)-(*ev2));
         const DataTypes::real_t max_d=d_01>d_12 ? d_01 : d_12;
         if (max_d<=tol*max_ev) {
             *V00=1.;
             *V10=0;
             *V20=0;
             *V01=0;
             *V11=1.;
             *V21=0;
             *V02=0;
             *V12=0;
             *V22=1.;
         } else {
            const DataTypes::real_t S00=A00-(*ev0);
            const DataTypes::real_t absS00=fabs(S00);
            if (absS00>m) {
                vectorInKernel3__nonZeroA00(S00,A01,A02,A01,A11-(*ev0),A12,A02,A12,A22-(*ev0),V00,V10,V20);
            } else if (absA02<m) {
                vectorInKernel3__nonZeroA00(A01,A11-(*ev0),A12,S00,A01,A02,A02,A12,A22-(*ev0),V00,V10,V20);
            } else {
                vectorInKernel3__nonZeroA00(A02,A12,A22-(*ev0),S00,A01,A02,A01,A11-(*ev0),A12,V00,V10,V20);
            }
            normalizeVector3(V00,V10,V20);;
            const DataTypes::real_t T00=A00-(*ev2);
            const DataTypes::real_t absT00=fabs(T00);
            if (absT00>m) {
                 vectorInKernel3__nonZeroA00(T00,A01,A02,A01,A11-(*ev2),A12,A02,A12,A22-(*ev2),V02,V12,V22);
            } else if (absA02<m) {
                 vectorInKernel3__nonZeroA00(A01,A11-(*ev2),A12,T00,A01,A02,A02,A12,A22-(*ev2),V02,V12,V22);
            } else {
                 vectorInKernel3__nonZeroA00(A02,A12,A22-(*ev2),T00,A01,A02,A01,A11-(*ev2),A12,V02,V12,V22);
            }
            const DataTypes::real_t dot=(*V02)*(*V00)+(*V12)*(*V10)+(*V22)*(*V20);
            *V02-=dot*(*V00);
            *V12-=dot*(*V10);
            *V22-=dot*(*V20);
            normalizeVector3(V02,V12,V22);
            *V01=(*V10)*(*V22)-(*V12)*(*V20);
            *V11=(*V20)*(*V02)-(*V00)*(*V22);
            *V21=(*V00)*(*V12)-(*V02)*(*V10);
            normalizeVector3(V01,V11,V21);
         }
   }
}

// General tensor product: arg_2(SL x SR) = arg_0(SL x SM) * arg_1(SM x SR)
// SM is the product of the last axis_offset entries in arg_0.getShape().
template <class LEFT, class RIGHT, class RES>
inline
void matrix_matrix_product(const int SL, const int SM, const int SR, const LEFT* A, const RIGHT* B, RES* C, int transpose)
{
  if (transpose == 0) {
    for (int i=0; i<SL; i++) {
      for (int j=0; j<SR; j++) {
        RES sum = 0.0;
        for (int l=0; l<SM; l++) {
          sum += A[i+SL*l] * B[l+SM*j];
        }
        C[i+SL*j] = sum;
      }
    }
  }
  else if (transpose == 1) {
    for (int i=0; i<SL; i++) {
      for (int j=0; j<SR; j++) {
        RES sum = 0.0;
        for (int l=0; l<SM; l++) {
          sum += A[i*SM+l] * B[l+SM*j];
        }
        C[i+SL*j] = sum;
      }
    }
  }
  else if (transpose == 2) {
    for (int i=0; i<SL; i++) {
      for (int j=0; j<SR; j++) {
        RES sum = 0.0;
        for (int l=0; l<SM; l++) {
          sum += A[i+SL*l] * B[l*SR+j];
        }
        C[i+SL*j] = sum;
      }
    }
  }
}

inline
DataTypes::real_t calc_erf(DataTypes::real_t x)
{
    return ::erf(x);
}

inline
DataTypes::cplx_t calc_erf(DataTypes::cplx_t x)
{
    return makeNaN();
}

inline DataTypes::real_t calc_sign(DataTypes::real_t x)
{
    return escript::fsign(x);
}

inline DataTypes::cplx_t calc_sign(DataTypes::cplx_t x)
{
    return makeNaN();
}

inline 
DataTypes::real_t calc_acos(DataTypes::real_t x)
{
    return acos(x);
}

inline 
DataTypes::cplx_t calc_acos(DataTypes::cplx_t x)
{
#ifdef ESYS_USE_BOOST_ACOS
              return boost::math::acos(x);
#else
              return acos(x);
#endif  
}


inline escript::DataTypes::real_t fabs(const escript::DataTypes::cplx_t c)
{
    return abs(c);
}



inline DataTypes::real_t calc_gtzero(const DataTypes::real_t& x) {return x>0;}
inline DataTypes::cplx_t calc_gtzero(const DataTypes::cplx_t& x) {return makeNaN();}


inline DataTypes::real_t calc_gezero(const DataTypes::real_t& x) {return x>=0;}
inline DataTypes::cplx_t calc_gezero(const DataTypes::cplx_t& x) {return makeNaN();}


inline DataTypes::real_t calc_ltzero(const DataTypes::real_t& x) {return x<0;}
inline DataTypes::cplx_t calc_ltzero(const DataTypes::cplx_t& x) {return makeNaN();}

inline DataTypes::real_t calc_lezero(const DataTypes::real_t& x) {return x<=0;}
inline DataTypes::cplx_t calc_lezero(const DataTypes::cplx_t& x) {return makeNaN();}

template <typename T>
inline DataTypes::real_t abs_f(T i)
{
    return fabs(i);
}

template <>
inline DataTypes::real_t abs_f(DataTypes::cplx_t i)
{
    return abs(i);
}




// deals with unary operations which return real, regardless of
// their input type
template <class T>
inline void tensor_unary_array_operation_real(const size_t size,
                             const T *arg1,
                             DataTypes::real_t * argRes,
                             escript::ES_optype operation,
                             DataTypes::real_t tol=0)
{
   switch (operation)
   {
     case REAL: 
          for (int i = 0; i < size; ++i) {
              argRes[i] = std::real(arg1[i]);
          }
          break;          
     case IMAG: 
          for (int i = 0; i < size; ++i) {
              argRes[i] = std::imag(arg1[i]);
          }
          break;  
    case EZ:   
          for (size_t i = 0; i < size; ++i) {
              argRes[i] = (fabs(arg1[i])<=tol);
          }
          break;
    case NEZ: 
          for (size_t i = 0; i < size; ++i) {
              argRes[i] = (fabs(arg1[i])>tol);
          }
          break;
    case ABS: 
          for (size_t i = 0; i < size; ++i) {
              argRes[i] = abs_f(arg1[i]);
          }
          break;  
    case PHS:
          for (size_t i = 0; i < size; ++i) {
              argRes[i] = std::arg(arg1[i]);
          }
          break;
     default:
	  std::ostringstream oss;
          oss << "Unsupported unary operation=";
	  oss << opToString(operation);
	  oss << '/';
	  oss << operation;
	  oss << " (Was expecting an operation with real results)";
          throw DataException(oss.str());      
   }  
}



template <typename T1, typename T2>
inline T1 conjugate(const T2 i)
{
    return conj(i);
}

// This should never actually be called
template <>
inline DataTypes::real_t conjugate(const DataTypes::real_t r)
{
    return r;
}

inline void tensor_unary_promote(const size_t size,
                             const DataTypes::real_t *arg1,
                             DataTypes::cplx_t * argRes)
{
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = arg1[i];
      }
}

// No openmp because it's called by Lazy
// In most cases, T1 and T2 will be the same
// but not ruling out putting Re() and Im()
// through this
template <class T1, typename T2>
inline void tensor_unary_array_operation(const size_t size,
                             const T1 *arg1,
                             T2 * argRes,
                             escript::ES_optype operation,
                             DataTypes::real_t tol=0)
{
  switch (operation)
  {
    case NEG:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = -arg1[i];
          }
          break;
    case SIN: 
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = sin(arg1[i]);
          }
          break;
    case COS:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = cos(arg1[i]);
          }
          break;
    case TAN:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = tan(arg1[i]);
          }
          break;
    case ASIN: 
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = asin(arg1[i]);
          }
          break;
    case ACOS:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i]=calc_acos(arg1[i]);
          }
          break;
    case ATAN:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = atan(arg1[i]);
          }
          break;
    case ABS:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = std::abs(arg1[i]);
          }
          break;      
    case SINH:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = sinh(arg1[i]);
          }
          break;
    case COSH:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = cosh(arg1[i]);
          }
          break;
    case TANH:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = tanh(arg1[i]);
          }
          break;
    case ERF: 
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = calc_erf(arg1[i]);
          }
          break;
    case ASINH:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = asinh(arg1[i]);
          }
          break;
    case ACOSH:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = acosh(arg1[i]);
          }
          break;
    case ATANH:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = atanh(arg1[i]);
          }
          break;
    case LOG10:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = log10(arg1[i]);
          }
          break;
    case LOG:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = log(arg1[i]);
          }
          break;      
    case SIGN:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = calc_sign(arg1[i]);
          }
          break;      
    case EXP:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = exp(arg1[i]);
          }
          break;      
    case SQRT:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = sqrt(arg1[i]);
          }
          break;      
    case GZ:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = calc_gtzero(arg1[i]);
          }
          break;      
    case GEZ:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = calc_gezero(arg1[i]);
          }
          break;            
    case LZ:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = calc_ltzero(arg1[i]);
          }
          break;            
    case LEZ:
	  for (size_t i = 0; i < size; ++i) {
              argRes[i] = calc_lezero(arg1[i]);
          }
          break;            
    case CONJ: 
          for (size_t i = 0; i < size; ++i) {
              argRes[i] = conjugate<T2,T1>(arg1[i]);
          }
          break; 
    case RECIP: 
          for (size_t i = 0; i < size; ++i) {
              argRes[i] = 1.0/arg1[i];
          }
          break; 
    case EZ:
          for (size_t i = 0; i < size; ++i) {
              argRes[i] = fabs(arg1[i])<=tol;
          }	  
	  break;
    case NEZ:
          for (size_t i = 0; i < size; ++i) {
              argRes[i] = fabs(arg1[i])>tol;
          }	  
	  break;      
      
    default:
      std::ostringstream oss;
      oss << "Unsupported unary operation=";
      oss << opToString(operation);
      oss << '/';
      oss << operation;
      throw DataException(oss.str());
  }
  return;
}

bool supports_cplx(escript::ES_optype operation);


} // end of namespace

#endif // __ESCRIPT_LOCALOPS_H__

