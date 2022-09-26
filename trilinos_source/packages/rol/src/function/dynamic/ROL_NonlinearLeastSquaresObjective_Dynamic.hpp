// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_NONLINEARLEASTSQUARESOBJECTIVE_DYNAMIC_H
#define ROL_NONLINEARLEASTSQUARESOBJECTIVE_DYNAMIC_H

#include "ROL_Objective.hpp"
#include "ROL_DynamicConstraint.hpp"
#include "ROL_Types.hpp"

/** @ingroup func_group
    \class ROL::NonlinearLeastSquaresObjective_Dynamic
    \brief Provides the interface to evaluate nonlinear least squares objective
           functions.

    ROL's nonlinear least squares objective function interface constructs the
    the nonlinear least squares objective function associated with the equality
    constraint \f$c_n(u_{n-1},u_n,z)=0\f$.  That is, given \f$z\f$ and \f$u_{n-1}\f$,
    \f[
       J(u) = \langle \mathfrak{R} c_n(u_{n-1},u,z),c_n(u_{n-1},u,z) \rangle_{\mathcal{C}^*,\mathcal{C}}
    \f]
    where \f$c_n:\mathcal{U}\times\mathcal{U}\times\mathcal{Z}\to\mathcal{C}\f$ and
    \f$\mathfrak{R}\in\mathcal{L}(\mathcal{C},\mathcal{C}^*)\f$ denotes the
    Riesz map from \f$\mathcal{C}\f$ into \f$\mathcal{C}^*\f$.

    ---
*/


namespace ROL {

template <class Real>
class DynamicConstraint;

template <class Real>
class NonlinearLeastSquaresObjective_Dynamic : public Objective<Real> {
private:
  const Ptr<DynamicConstraint<Real>> con_;
  const Ptr<const Vector<Real>> uo_;
  const Ptr<const Vector<Real>> z_;
  const Ptr<const TimeStamp<Real>> ts_;
  const bool GaussNewtonHessian_;

  Ptr<Vector<Real> > c1_, c2_, cdual_, udual_;

public:
  /** \brief Constructor. 

      This function constructs a nonlinear least squares objective function. 
      @param[in]          con   is the nonlinear equation to be solved.
      @param[in]          vec   is a constraint space vector used for cloning.
      @param[in]          GHN   is a flag dictating whether or not to use the Gauss-Newton Hessian.
  */
  NonlinearLeastSquaresObjective_Dynamic(const Ptr<DynamicConstraint<Real>> &con,
                                         const Vector<Real> &c,
                                         const Ptr<const Vector<Real>> &uo,
                                         const Ptr<const Vector<Real>> &z,
                                         const Ptr<const TimeStamp<Real>> &ts,
                                         const bool GNH = false)
    : con_(con), uo_(uo), z_(z), ts_(ts), GaussNewtonHessian_(GNH) {
    c1_ = c.clone();
    c2_ = c.clone();
    cdual_ = c.dual().clone();
    udual_ = uo->dual().clone();
  }

  void update( const Vector<Real> &u, bool flag = true, int iter = -1 ) {
    //con_->update_un(u,*ts_);
    con_->update(*uo_,u,*z_,*ts_);
    con_->value(*c1_,*uo_,u,*z_,*ts_);
    cdual_->set(c1_->dual());
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    Real half(0.5);
    return half*(c1_->dot(*cdual_));
  }

  void gradient( Vector<Real> &g, const Vector<Real> &u, Real &tol ) {
    con_->applyAdjointJacobian_un(g,*cdual_,*uo_,u,*z_,*ts_);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, Real &tol ) {
    con_->applyJacobian_un(*c2_,v,*uo_,u,*z_,*ts_);
    con_->applyAdjointJacobian_un(hv,c2_->dual(),*uo_,u,*z_,*ts_);
    if ( !GaussNewtonHessian_ ) {
      con_->applyAdjointHessian_un_un(*udual_,*cdual_,v,*uo_,u,*z_,*ts_);
      hv.plus(*udual_);
    }
  }

  void precond( Vector<Real> &pv, const Vector<Real> &v, const Vector<Real> &u, Real &tol ) {
    con_->applyInverseAdjointJacobian_un(*cdual_,v,*uo_,u,*z_,*ts_);
    con_->applyInverseJacobian_un(pv,cdual_->dual(),*uo_,u,*z_,*ts_);
  }

// Definitions for parametrized (stochastic) equality constraints
//public:
//  void setParameter(const std::vector<Real> &param) {
//    Objective<Real>::setParameter(param);
//    con_->setParameter(param);
//  }
};

} // namespace ROL

#endif
