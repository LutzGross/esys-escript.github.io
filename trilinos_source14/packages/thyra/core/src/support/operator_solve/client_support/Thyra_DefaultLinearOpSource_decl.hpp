// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFUALT_LINEAR_OP_SOURCE_DECL_HPP
#define THYRA_DEFUALT_LINEAR_OP_SOURCE_DECL_HPP

#include "Thyra_LinearOpSourceBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"

namespace Thyra {

/** \brief Default implementation of a <tt>LinearOpSourceBase</tt> that just
 * accepts and gives up a single linear operator object.
 */
template<class Scalar>
class DefaultLinearOpSource : virtual public LinearOpSourceBase<Scalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized.
   */
  DefaultLinearOpSource();

  /** \brief Construct with a non-const linear operator.
   */
  DefaultLinearOpSource(
    const Teuchos::RCP<LinearOpBase<Scalar> >    &op
    );

  /** \brief Construct with a const linear operator.
   */
  DefaultLinearOpSource(
    const Teuchos::RCP<const LinearOpBase<Scalar> >    &op
    );

  /** \brief Initialize with a non-const linear operator.
   */
  void initialize(
    const Teuchos::RCP<LinearOpBase<Scalar> >    &op
    );

  /** \brief Initialize with a const linear operator.
   */
  void initialize(
    const Teuchos::RCP<const LinearOpBase<Scalar> >    &op
    );

  /** \brief Uninitialize.
   *
   * Note: If the client wants to access the underlying linear operator, then
   * it had better grab them with the below access functions before calling
   * this function.
   */
  void uninitialize();

  //@}

  // ToDo: Override functions from Describable!

  /** @name Overridden from LinearOpSourceBase */
  //@{
  /** \brief . */
  bool isOpConst() const;
  /** \brief . */
  Teuchos::RCP<LinearOpBase<Scalar> > getNonconstOp();
  /** \brief . */
  Teuchos::RCP<const LinearOpBase<Scalar> > getOp() const;
  //@}
  
private:
  
  Teuchos::ConstNonconstObjectContainer<LinearOpBase<Scalar> >  op_;
  
};

// //////////////////////////////
// Related functions

/** \brief Create a <tt>DefaultLinearOpSource</tt> object out of a
 * <tt>LinearOpBase</tt> object.
 *
 * \relates DefaultLinearOpSource
 */
template <class Scalar>
Teuchos::RCP<const DefaultLinearOpSource<Scalar> >
defaultLinearOpSource(
  const Teuchos::RCP<const LinearOpBase<Scalar> >    &op
  )
{
  return Teuchos::rcp(new DefaultLinearOpSource<Scalar>(op));
}

} // namespace Thyra

#endif // THYRA_DEFUALT_LINEAR_OP_SOURCE_DECL_HPP
