// $Id$
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

#if !defined escript_LocalOps_H
#define escript_LocalOps_H
#ifdef __INTEL_COMPILER
#include <mathimf.h>
# define M_PI           3.14159265358979323846  /* pi */
#else
#include <math.h>
#endif

namespace escript {


/**
   \brief
   solves a 1x1 eigenvalue A*V=ev*V problem

   \param A00 Input - A_00 
   \param ev0 Output - eigenvalue
*/
inline
void eigenvalues1(const double A00,double* ev0) {

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
inline
void eigenvalues2(const double A00,const double A01,const double A11,
                 double* ev0, double* ev1) {
      const register double trA=(A00+A11)/2.;
      const register double A_00=A00-trA;
      const register double A_11=A11-trA;
      const register double s=sqrt(A01*A01-A_00*A_11);
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
void eigenvalues3(const double A00, const double A01, const double A02,
                                   const double A11, const double A12,
                                                     const double A22,
                 double* ev0, double* ev1,double* ev2) {

      const register double trA=(A00+A11+A22)/3.;
      const register double A_00=A00-trA;
      const register double A_11=A11-trA;
      const register double A_22=A22-trA;
      const register double A01_2=A01*A01;
      const register double A02_2=A02*A02;
      const register double A12_2=A12*A12;
      const register double p=A02_2+A12_2+A01_2+(A_00*A_00+A_11*A_11+A_22*A_22)/2.;
      if (p<=0.) {
         *ev2=trA;
         *ev1=trA;
         *ev0=trA;

      } else {
         const register double q=(A02_2*A_11+A12_2*A_00+A01_2*A_22)-(A_00*A_11*A_22+2*A01*A12*A02);
         const register double sq_p=sqrt(p/3.);
         register double z=-q/(2*pow(sq_p,3));
         if (z<-1.) {
            z=-1.;
         } else if (z>1.) {
            z=1.;
         }
         const register double alpha_3=acos(z)/3.;
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
void  eigenvalues_and_eigenvectors1(const double A00,double* ev0,double* V00,const double tol)
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
void  vectorInKernel2(const double A00,const double A10,const double A01,const double A11,
                      double* V0, double*V1)
{
      register double absA00=fabs(A00);
      register double absA10=fabs(A10);
      register double absA01=fabs(A01);
      register double absA11=fabs(A11);
      register double m=absA11>absA10 ? absA11 : absA10;
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
void  vectorInKernel3__nonZeroA00(const double A00,const double A10,const double A20,
                                const double A01,const double A11,const double A21,
                                const double A02,const double A12,const double A22,
                                double* V0,double* V1,double* V2)
{
    double TEMP0,TEMP1;
    register const double I00=1./A00;
    register const double IA10=I00*A10;
    register const double IA20=I00*A20;
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
void  eigenvalues_and_eigenvectors2(const double A00,const double A01,const double A11,
                                    double* ev0, double* ev1,
                                    double* V00, double* V10, double* V01, double* V11,
                                    const double tol) 
{
     double TEMP0,TEMP1;
     eigenvalues2(A00,A01,A11,ev0,ev1);
     const register double absev0=fabs(*ev0);
     const register double absev1=fabs(*ev1);
     register double max_ev=absev0>absev1 ? absev0 : absev1;
     if (fabs((*ev0)-(*ev1))<tol*max_ev) {
        *V00=1.;
        *V10=0.;
        *V01=0.;
        *V11=1.;
     } else {
        vectorInKernel2(A00-(*ev0),A01,A01,A11-(*ev0),&TEMP0,&TEMP1);
        const register double scale=1./sqrt(TEMP0*TEMP0+TEMP1*TEMP1);
        if (TEMP0<0.) {
            *V00=-TEMP0*scale;
            *V10=-TEMP1*scale;
            if (TEMP1<0.) {
               *V01=  *V10; 
               *V11=-(*V00);
            } else {
               *V01=-(*V10);
               *V11= (*V10);
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
void  normalizeVector3(double* V0,double* V1,double* V2)
{
    register double s;
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
void  eigenvalues_and_eigenvectors3(const double A00, const double A01, const double A02,
                                    const double A11, const double A12, const double A22,
                                    double* ev0, double* ev1, double* ev2,
                                    double* V00, double* V10, double* V20, 
                                    double* V01, double* V11, double* V21, 
                                    double* V02, double* V12, double* V22, 
                                    const double tol)
{
      register const double absA01=fabs(A01);
      register const double absA02=fabs(A02);
      register const double m=absA01>absA02 ? absA01 : absA02;
      if (m<=0) {
        double TEMP_V00,TEMP_V10,TEMP_V01,TEMP_V11,TEMP_EV0,TEMP_EV1;
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
         const register double absev0=fabs(*ev0);
         const register double absev1=fabs(*ev1);
         const register double absev2=fabs(*ev2);
         register double max_ev=absev0>absev1 ? absev0 : absev1;
         max_ev=max_ev>absev2 ? max_ev : absev2;
         const register double d_01=fabs((*ev0)-(*ev1));
         const register double d_12=fabs((*ev1)-(*ev2));
         const register double max_d=d_01>d_12 ? d_01 : d_12;
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
            const register double S00=A00-(*ev0);
            const register double absS00=fabs(S00);
            if (fabs(S00)>m) {
                vectorInKernel3__nonZeroA00(S00,A01,A02,A01,A11-(*ev0),A12,A02,A12,A22-(*ev0),V00,V10,V20);
            } else if (absA02<m) {
                vectorInKernel3__nonZeroA00(A01,A11-(*ev0),A12,S00,A01,A02,A02,A12,A22-(*ev0),V00,V10,V20);
            } else {
                vectorInKernel3__nonZeroA00(A02,A12,A22-(*ev0),S00,A01,A02,A01,A11-(*ev0),A12,V00,V10,V20);
            }
            normalizeVector3(V00,V10,V20);;
            const register double T00=A00-(*ev2);
            const register double absT00=fabs(T00);
            if (fabs(T00)>m) {
                 vectorInKernel3__nonZeroA00(T00,A01,A02,A01,A11-(*ev2),A12,A02,A12,A22-(*ev2),V02,V12,V22);
            } else if (absA02<m) {
                 vectorInKernel3__nonZeroA00(A01,A11-(*ev2),A12,T00,A01,A02,A02,A12,A22-(*ev2),V02,V12,V22);
            } else {
                 vectorInKernel3__nonZeroA00(A02,A12,A22-(*ev2),T00,A01,A02,A01,A11-(*ev2),A12,V02,V12,V22);
            }
            const register double dot=(*V02)*(*V00)+(*V12)*(*V10)+(*V22)*(*V20);
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
} // end of namespace
#endif
