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

#ifndef PDEOPT_MODEL_H
#define PDEOPT_MODEL_H

#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"

template <class Real>
class Model : public ROL::Constraint<Real> {
private:
  const std::vector<ROL::Ptr<ROL::Objective<Real>>> objs_;
  ROL::Ptr<ROL::Vector<Real>> xdual_;

public:
  Model(const std::vector<ROL::Ptr<ROL::Objective<Real>>> &objs)
    : objs_(objs) {}

  void setParameter(const std::vector<Real> &param) {
    ROL::Constraint<Real>::setParameter(param);
    for (const auto obj : objs_) obj->setParameter(param);
  }

  void update(const ROL::Vector<Real> &x, bool flag = true, int iter = -1) {
    for (const auto obj : objs_) obj->update(x,flag,iter);
  }

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &x, Real &tol ) {
    ROL::Ptr<std::vector<Real>> cdata = static_cast<ROL::StdVector<Real>&>(c).getVector();
    int nobs = objs_.size();
    for (int i = 0; i < nobs; ++i) {
      (*cdata)[i] = objs_[i]->value(x,tol);
    }
  }

  void applyJacobian(ROL::Vector<Real> &jv,
                     const ROL::Vector<Real> &v,
                     const ROL::Vector<Real> &x,
                     Real &tol) {
    if (xdual_==ROL::nullPtr) xdual_ = x.dual().clone();
    ROL::Ptr<std::vector<Real>> jdata = static_cast<ROL::StdVector<Real>&>(jv).getVector();
    int nobs = objs_.size();
    for (int i = 0; i < nobs; ++i) {
      xdual_->zero();
      objs_[i]->gradient(*xdual_,x,tol);
      (*jdata)[i] = xdual_->dot(v.dual());
    }
  }

  void applyAdjointJacobian(ROL::Vector<Real> &ajv,
                            const ROL::Vector<Real> &v,
                            const ROL::Vector<Real> &x,
                            Real &tol) {
    if (xdual_==ROL::nullPtr) xdual_ = x.dual().clone();
    ROL::Ptr<const std::vector<Real>> vdata = static_cast<const ROL::StdVector<Real>&>(v).getVector();
    ajv.zero();
    int nobs = objs_.size();
    for (int i = 0; i < nobs; ++i) {
      xdual_->zero();
      objs_[i]->gradient(*xdual_,x,tol);
      ajv.axpy((*vdata)[i],*xdual_);
    }
  }

}; // class Model

#endif
