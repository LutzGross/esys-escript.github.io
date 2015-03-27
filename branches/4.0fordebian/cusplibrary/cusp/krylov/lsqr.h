/*
 *  Copyright 2011 The Regents of the University of California
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#pragma once

#include <cusp/detail/config.h>

namespace cusp
{
   namespace krylov
   {

      /*! \addtogroup iterative_solvers Iterative Solvers
       *  \addtogroup krylov_methods Krylov Methods
       *  \ingroup iterative_solvers
       *  \{
       */     
     /*!
       This class is used to specify more detailed parameters to the LSQR solver.
       
       The following quantities are used in discussing the structure:
       \f[   \mathrm{Abar} =
\left[ {\begin{array}{cc}
 A  \\
 \mathrm{damp} \cdot I\\
 \end{array} } \right]
     \f]
       
     relpr  =  the relative precision of floating-point arithmetic
     on the machine being used.  On most machines,
     relpr is about 1.0e-7 and 1.0d-16 in single and double
     precision respectively.
     */
     template <typename Real>
       struct lsqr_parameters
       {
         /*!
           \param damp The damping parameter for problem 3 above.
           (damp should be 0.0 for problems 1 and 2.)
           If the system A*x = b is incompatible, values
           of damp in the range 0 to sqrt(relpr)*norm(A)
           will probably have a negligible effect.
           Larger values of damp will tend to decrease
           the norm of x and reduce the number of 
           iterations required by LSQR.    
           The work per iteration and the storage needed
           by LSQR are the same for all values of damp.
           \param atol An estimate of the relative error in the data
           defining the matrix A.  For example,
           if A is accurate to about 6 digits, set
           atol = 1.0e-6.  If 0 machine precision will be used.
           \param  btol An estimate of the relative error in the data
           defining the rhs vector b.  For example,
           if b is accurate to about 6 digits, set
           btol = 1.0e-6. If 0 machine precision will be used.
           If both atol and btol are 0 the relative tolerance from the monitor
           parameter to lsqr will be used instead.
           \param conlim An upper limit on cond(Abar), the apparent
           condition number of the matrix Abar.
           Iterations will be terminated if a computed
           estimate of cond(Abar) exceeds conlim.
           This is intended to prevent certain small or
           zero singular values of A or Abar from
           coming into effect and causing unwanted growth
           in the computed solution.
           conlim and damp may be used separately or
           together to regularize ill-conditioned systems.         
           Normally, conlim should be in the range
           1000 to 1/relpr.
           Suggested value:
           conlim = 1/(100*relpr)  for compatible systems,
           conlim = 1/(10*sqrt(relpr)) for least squares.
           If set to 0 then 1/(machine precision) will be used.
         */
         lsqr_parameters(Real _damp = 0, Real _atol = 0, Real _btol = 0, Real _conlim = 0):
         damp(_damp),    
           atol(_atol),  
           btol(_btol),  
           conlim(_conlim)
           {
           }
         Real damp;
         Real atol;
         Real btol;
         Real conlim;
       };

     /*! This class holds several results from of an LSQR run 
      
       The following quantities are used in discussing the structure:
       \f[   \mathrm{Abar} =
       \left[ {\begin{array}{cc}
       A  \\
       \mathrm{damp} \cdot I \\
       \end{array} } \right]
       \hspace{1cm}
       \mathrm{rbar} =  \left[ {\begin{array}{cc}
       b  \\
       0 \\
       \end{array} } \right] - \mathrm{Abar} \cdot x
       \f]

     relpr  =  the relative precision of floating-point arithmetic
     on the machine being used.  On most machines,
     relpr is about 1.0e-7 and 1.0d-16 in single and double
     precision respectively.
     */
     template <typename Real>
       struct lsqr_results
       {
         /*! This enumeration describes the reason why the lsqr solver stopped.
          */
         enum StopCondition{
           /*! x = 0  is the exact solution.
             No iterations were performed. */
           TrivialSolution=0, 
           /*! The equations A*x = b are probably
             compatible.  Norm(A*x - b) is sufficiently
             small, given the values of atol and btol. */
           ExactSolution, 
           /*! damp is zero.  The system A*x = b is probably
             not compatible.  A least-squares solution has
             been obtained that is sufficiently accurate,
             given the value of atol. */
           LeastSquaresSolution, 
           /*! damp is nonzero.  A damped 
             least-squares solution has been
             obtained that is sufficiently
             accurate, given the value of atol.*/
           DampedLeastSquaresSolution, 
           /*! An estimate of cond(Abar) has exceeded
             the condition number limit.  The system A*x = b appears to 
             be ill-conditioned. */
           AbarConditionNumberExceeded, 
           /*! The iteration limit was reached */
           IterationLimitReached 
         };
       lsqr_results(int _istop, Real _anorm, Real _acond, Real _rnorm, 
                   Real _arnorm, Real _xnorm):
           anorm(_anorm),
           acond(_acond),
           rnorm(_rnorm),
           arnorm(_arnorm),
           xnorm(_xnorm)
           {
             istop = StopCondition(_istop);
           }
         /*!
           The reason why the lsqr solver stopped.
           \sa StopCondition
         */
         StopCondition istop;
         /*! An estimate of the Frobenius norm of  Abar.
           This is the square-root of the sum of squares
           of the elements of Abar.
           If damp is small and if the columns of A
           have all been scaled to have length 1.0,
           anorm should increase to roughly sqrt(n).
           A radically different value for anorm may
           indicate an error in subroutine aprod (there
           may be an inconsistency between modes 1 and 2).
         */
         Real anorm;
         /*! An estimate of cond(Abar), the condition
           number of Abar.  A very high value of acond
           may again indicate a problem in the spmv.
         */
         Real acond;
         /*! An estimate of the final value of
           norm( Abar(transpose)*rbar ), the norm of
           the residual for the usual normal equations.
           This should be small in all cases.  (arnorm
           will often be smaller than the true value
           computed from the output vector x.)
         */
         Real rnorm;
         /*! An estimate of the final value of
           norm( Abar(transpose)*rbar ), the norm of
           the residual for the usual normal equations.
           This should be small in all cases.  (arnorm
           will often be smaller than the true value
           computed from the output vector x.)
         */
         Real arnorm;
         /*! An estimate of the norm of the final
           solution vector x.
         */
         Real xnorm;     
       };


      /*! \p lsqr : LSQR solver
       *
       * Solves:
       *    Ax = b
       *    in the case of square A
       *  or
       *    minimize ||Ax - b||^2 (least squares)
       *    for rectangular A when d == 0
       *  or
       *    minimize ||Ax - b||^2 + d^2 ||x||^2 (dampened least squares)
       *    for rectangular A when d > 0
       *    
       * using the default convergence criteria.
       */
     template <class LinearOperator, class Vector, typename RealType>
       lsqr_results<RealType> lsqr(LinearOperator& A,
                                  Vector& x,
                                  Vector& b,
                                  const lsqr_parameters<RealType> & p);
     
      /*! \p lsqr : LSQR solver
       *
       * Solves:
       *    Ax = b
       *    for square A when d == 0
       *  or
       *    minimize ||Ax - b||^2 (least squares)
       *    for rectangular A when d == 0
       *  or
       *    minimize ||Ax - b||^2 + d^2 ||x||^2 (dampened least squares)
       *    when d > 0
       *
       * \param A matrix of the linear system 
       * \param At transpose of the matrix of the linear system 
       * \param x approximate solution of the linear system
       * \param b right-hand side of the linear system
       * \param d the damping factor
       * \param monitor montiors iteration and determines stopping conditions
       *
       * \tparam LinearOperator is a matrix or subclass of \p linear_operator
       * \tparam Vector vector
       * \tparam Monitor is a monitor such as \p default_monitor or \p verbose_monitor
       *
       *  The following code snippet demonstrates how to use \p lsqr to 
       *  solve a 10x10 Poisson problem.
       *
       *  \code
       *  #include <cusp/csr_matrix.h>
       *  #include <cusp/monitor.h>
       *  #include <cusp/krylov/lsqr.h>
       *  #include <cusp/gallery/poisson.h>
       *  
       *  int main(void)
       *  {
       *      // create an empty sparse matrix structure (CSR format)
       *      cusp::csr_matrix<int, float, cusp::device_memory> A;
       *
       *      // initialize matrix
       *      cusp::gallery::poisson5pt(A, 10, 10);
       *
       *      // allocate storage for solution (x) and right hand side (b)
       *      cusp::array1d<float, cusp::device_memory> x(A.num_rows, 0);
       *      cusp::array1d<float, cusp::device_memory> b(A.num_rows, 1);
       *
       *      // set stopping criteria:
       *      //  iteration_limit    = 100
       *      //  relative_tolerance = 1e-6
       *      cusp::verbose_monitor<float> monitor(b, 100, 1e-6);
       *
       *      // solve the linear system A x = b
       *      // because A is hermitian we can use 
       *      // it for its own conjugate transpose
       *      cusp::krylov::lsqr(A, A, x, b, cusp::krylov::lsqr_paramters<float>(), monitor);
       *
       *      return 0;
       *  }
       *  \endcode

       *  \see \p default_monitor
       *  \see \p verbose_monitor
       *
       */
      template <class LinearOperator,
                class Vector,
                typename RealType,
                class Monitor>
       lsqr_results<RealType> lsqr(LinearOperator& A,
                                   Vector& x,
                                   Vector& b,
                                   const lsqr_parameters<RealType> & p,
                                   Monitor& monitor);

      /*! \}
      */

   } // end namespace krylov
} // end namespace cusp

#include <cusp/krylov/detail/lsqr.inl>

