//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
//@HEADER

#ifndef NOX_OBSERVER_VECTOR_HPP
#define NOX_OBSERVER_VECTOR_HPP

#include "NOX_Common.H"
#include "Teuchos_RCP.hpp"
#include "NOX_Observer.hpp"
#include <vector>

namespace NOX {

  /** \brief Concrete implementation of NOX::Observer that stores a vector of Observers.

      The intent of this object to to aggregate a set of Observer objects.
  */
  class ObserverVector : public NOX::Observer {
  public:
    // Derived methods
    void runPreIterate(const NOX::Solver::Generic& solver);
    void runPostIterate(const NOX::Solver::Generic& solver);
    void runPreSolve(const NOX::Solver::Generic& solver);
    void runPostSolve(const NOX::Solver::Generic& solver);
    void runPreSolutionUpdate(const NOX::Abstract::Vector& update, const NOX::Solver::Generic& solver);
    void runPostSolutionUpdate(const NOX::Solver::Generic& solver);
    void runPreLineSearch(const NOX::Solver::Generic& solver);
    void runPostLineSearch(const NOX::Solver::Generic& solver);

    //! Add observer to end of vector.
    void pushBack(const Teuchos::RCP<NOX::Observer>& observer);
    //! Remove observer from end of vector.
    void popBack();
    //! Clear the vector of observers.
    void clear();

  private:
    //! std::vector of observer objects
    std::vector<Teuchos::RCP<NOX::Observer>> vec_;
  };
}

#endif
