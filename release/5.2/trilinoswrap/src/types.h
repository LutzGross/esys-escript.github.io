
/*****************************************************************************
*
* Copyright (c) 2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#ifndef __ESYS_TRILINOSWRAP_TYPES_H__
#define __ESYS_TRILINOSWRAP_TYPES_H__

#include <escript/DataTypes.h>

#ifndef ESYS_MPI
#include <escript/EsysMPI.h> // for MPI_Comm typedef
#include <Teuchos_DefaultComm.hpp>
#endif
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_RowMatrix.hpp>

#include <Tpetra_Experimental_BlockVector.hpp>

namespace esys_trilinos {

/// Scalar types
typedef escript::DataTypes::real_t  real_t;
typedef escript::DataTypes::cplx_t  cplx_t;

/// Global Ordinal type
typedef escript::DataTypes::index_t GO;
/// Local Ordinal type
typedef escript::DataTypes::index_t LO;
/// Kokkos Node type
#ifdef _OPENMP
typedef Kokkos::Compat::KokkosOpenMPWrapperNode NT;
#elif ESYS_HAVE_CUDA
typedef Kokkos::Compat::KokkosCudaWrapperNode   NT;
#else
typedef Kokkos::Compat::KokkosSerialWrapperNode NT;
#endif

typedef Tpetra::CrsGraph<LO,GO,NT>    GraphType;
typedef Tpetra::Import<LO,GO,NT>      ImportType;
typedef Teuchos::RCP<GraphType>       TrilinosGraph_ptr;
typedef Teuchos::RCP<const GraphType> const_TrilinosGraph_ptr;
typedef GraphType::map_type           MapType;
typedef Teuchos::RCP<MapType>         TrilinosMap_ptr;
typedef Teuchos::RCP<const MapType>   const_TrilinosMap_ptr;

template<typename ST> using MatrixType = Tpetra::RowMatrix<ST,LO,GO,NT>;
template<typename ST> using VectorType = Tpetra::MultiVector<ST,LO,GO,NT>;
template<typename ST> using OpType     = Tpetra::Operator<ST,LO,GO,NT>;

typedef VectorType<real_t> RealVector;
typedef OpType<real_t>     RealOperator;

typedef VectorType<cplx_t> ComplexVector;
typedef OpType<cplx_t>     ComplexOperator;

// experimental block types
template<typename ST> using BlockVectorType =
                               Tpetra::Experimental::BlockVector<ST,LO,GO,NT>;
typedef BlockVectorType<real_t> RealBlockVector;
typedef BlockVectorType<cplx_t> ComplexBlockVector;


/// converts an MPI communicator to a Teuchos communicator
inline
Teuchos::RCP<const Teuchos::Comm<int> > TeuchosCommFromEsysComm(MPI_Comm comm)
{
#ifdef ESYS_MPI
    return Teuchos::rcp(new Teuchos::MpiComm<int>(comm));
#else
    return Teuchos::DefaultComm<int>::getComm();
#endif
}

} // namespace esys_trilinos

#endif // __ESYS_TRILINOSWRAP_TYPES_H__

