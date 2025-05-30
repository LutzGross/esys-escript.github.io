// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
#ifndef TPETRAEXAMPLES_FEM_ASSEMBLY_UTILITY_HPP
#define TPETRAEXAMPLES_FEM_ASSEMBLY_UTILITY_HPP

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include <iostream>

namespace TpetraExamples {

template<typename V>
Kokkos::View<typename V::data_type, typename V::array_layout, typename V::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
makeUnmanaged(const V& v)
{
  return v;
}

// Return a pointer (RCP is like std::shared_ptr) to an output
// stream.  It prints on Process 0 of the given MPI communicator,
// but ignores all output on other MPI processes.
Teuchos::RCP<Teuchos::FancyOStream>
getOutputStream (const Teuchos::Comm<int>& comm)
{
  using Teuchos::getFancyOStream;

  const int myRank = comm.getRank ();
  if (0 == myRank) {
    // Process 0 of the given communicator prints to std::cout.
    return getFancyOStream (Teuchos::rcpFromRef (std::cout));
  }
  else {
    // A "black hole output stream" ignores all output directed to it.
    return getFancyOStream (Teuchos::rcp (new Teuchos::oblackholestream ()));
  }
}

} // namespace TpetraExamples

#endif  // TPETRAEXAMPLES_FEM_ASSEMBLY_UTILITY_HPP

