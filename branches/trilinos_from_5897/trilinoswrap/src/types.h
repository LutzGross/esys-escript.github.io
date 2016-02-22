
/*****************************************************************************
*
* Copyright (c) 2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#ifndef __TRILINOSWRAP_TYPES_H__
#define __TRILINOSWRAP_TYPES_H__

#include <escript/DataTypes.h>

#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Amesos2_Solver_decl.hpp>
#include <BelosSolverManager.hpp>

namespace esys_trilinos {

/// Scalar type
typedef double  ST;
/// Global Ordinal type
typedef escript::DataTypes::index_t GO;
/// Local Ordinal type
typedef escript::DataTypes::index_t LO;
/// Kokkos Node type
#ifdef _OPENMP
typedef Kokkos::Compat::KokkosOpenMPWrapperNode NT;
#elif USE_CUDA
typedef Kokkos::Compat::KokkosCudaWrapperNode NT;
#else
typedef Kokkos::Compat::KokkosSerialWrapperNode NT;
#endif

typedef Tpetra::CrsGraph<LO,GO,NT>                 GraphType;
typedef Tpetra::CrsMatrix<ST,LO,GO,NT>             MatrixType;
typedef Tpetra::MultiVector<ST,LO,GO,NT>           VectorType;
typedef Tpetra::Operator<ST,LO,GO,NT>              OpType;
typedef Belos::LinearProblem<ST,VectorType,OpType> ProblemType;
typedef Belos::SolverManager<ST,VectorType,OpType> SolverType;
typedef Amesos2::Solver<MatrixType,VectorType>     DirectSolverType;
typedef MatrixType::map_type                       MapType;
typedef Tpetra::Import<LO,GO,NT>                   ImportType;

typedef Teuchos::RCP<MapType> TrilinosMap_ptr;
typedef Teuchos::RCP<const MapType> const_TrilinosMap_ptr;

typedef Teuchos::RCP<GraphType> TrilinosGraph_ptr;
typedef Teuchos::RCP<const GraphType> const_TrilinosGraph_ptr;

inline
Teuchos::RCP<const Teuchos::Comm<int> > TeuchosCommFromEsysComm(MPI_Comm comm)
{
    return Teuchos::rcp(new Teuchos::MpiComm<int>(comm));
}

} // namespace esys_trilinos

#endif // __TRILINOSWRAP_TYPES_H__
