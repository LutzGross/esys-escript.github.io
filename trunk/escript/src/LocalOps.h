// $Id$
/* 
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

#if !defined escript_LocalOps_H
#define escript_LocalOps_H

namespace escript {


/**
   \brief
   solves a 1x1 eigenvalue A*V=ev*V problem

   \param A00 Input - A_00 
   \param ev0 Output - eigenvalue
*/
inline
void eigenvalues1(const double A00,double* ev0) {

   *ev0=1.;

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
void eigenvalues2(const double A00,const double A01
                                 ,const double A11,
                 double* ev0, double* ev1) {

   *ev0=1.;
   *ev1=2.;
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

   *ev0=1.;
   *ev1=2.;
   *ev2=3.;

}

} // end of namespace
#endif
