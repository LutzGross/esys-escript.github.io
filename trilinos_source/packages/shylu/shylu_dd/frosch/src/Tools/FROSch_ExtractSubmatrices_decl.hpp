//@HEADER
// ************************************************************************
//
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
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
// Questions? Contact Alexander Heinlein (alexander.heinlein@uni-koeln.de)
//
// ************************************************************************
//@HEADER

#ifndef _FROSCH_EXTRACTSUBMATRICES_DECL_HPP
#define _FROSCH_EXTRACTSUBMATRICES_DECL_HPP

#include <ShyLU_DDFROSch_config.h>

#include <FROSch_Timers.h>

#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>
#include <Xpetra_ImportFactory_fwd.hpp>

#include <KokkosKernels_Utils.hpp>

namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    RCP<const Matrix<SC,LO,GO,NO> > ExtractLocalSubdomainMatrix(RCP<const Matrix<SC,LO,GO,NO> > globalMatrix,
                                                                RCP<const Map<LO,GO,NO> > map);

    // ----------------------------------------------------------- //
    // split ExtractLocalSubdomainMatrix into symbolic / compute
    template <class SC,class LO,class GO,class NO>
    void ExtractLocalSubdomainMatrix_Symbolic(RCP<Matrix<SC,LO,GO,NO> > subdomainMatrix,        // input  : globalMatrix, re-distributed with map
                                              RCP<Matrix<SC,LO,GO,NO> > localSubdomainMatrix);  // output : local submatrix

    template <class SC,class LO,class GO,class NO>
    void ExtractLocalSubdomainMatrix_Compute(RCP<const Matrix<SC,LO,GO,NO> > globalMatrix,
                                             RCP<      Matrix<SC,LO,GO,NO> > subdomainMatrix,
                                             RCP<      Matrix<SC,LO,GO,NO> > repeatedMatrix);
    // ----------------------------------------------------------- //

    template <class SC,class LO,class GO,class NO>
    RCP<const Matrix<SC,LO,GO,NO> > ExtractLocalSubdomainMatrix(RCP<const Matrix<SC,LO,GO,NO> > globalMatrix,
                                                                RCP<const Map<LO,GO,NO> > map,
                                                                SC value);

    template <class SC,class LO,class GO,class NO>
    int UpdateLocalSubdomainMatrix(RCP<Matrix<SC,LO,GO,NO> > globalMatrix,
                                   RCP<Map<LO,GO,NO> > &map,
                                   RCP<Matrix<SC,LO,GO,NO> > &localSubdomainMatrix);

    template <class SC,class LO,class GO,class NO>
    int BuildSubmatrices(RCP<const Matrix<SC,LO,GO,NO> > k,
                         ArrayView<GO> indI,
                         RCP<const Matrix<SC,LO,GO,NO> > &kII,
                         RCP<const Matrix<SC,LO,GO,NO> > &kIJ,
                         RCP<const Matrix<SC,LO,GO,NO> > &kJI,
                         RCP<const Matrix<SC,LO,GO,NO> > &kJJ);

    template <class SC,class LO,class GO,class NO>
    int BuildSubmatrix(RCP<const Matrix<SC,LO,GO,NO> > k,
                       ArrayView<GO> indI,
                       RCP<const Matrix<SC,LO,GO,NO> > &kII);

    template <class LO,class GO,class NO>
    int BuildSubgraph(RCP<const CrsGraph<LO,GO,NO> > k,
                      ArrayView<GO> indI,
                      RCP<const CrsGraph<LO,GO,NO> > &kII);
}

#endif
