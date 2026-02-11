
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __ESYS_TRILINOSWRAP_TYPES_H__
#define __ESYS_TRILINOSWRAP_TYPES_H__

#include <escript/DataTypes.h>
#include <escript/EsysMPI.h> // for MPI_Comm typedef (always needed)

// Get Trilinos MPI configuration (defines HAVE_MPI if Trilinos has MPI)
#include <Teuchos_config.h>

// Only include MpiComm if both escript AND Trilinos have MPI
#if defined(ESYS_MPI) && defined(HAVE_MPI)
#include <Teuchos_DefaultMpiComm.hpp>
#define TRILINOS_HAS_MPI 1
#endif
#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_RowMatrix.hpp>
// Use standard BlockVector (experimental version deprecated in Trilinos 16.2+)
#include <Tpetra_BlockVector.hpp>

namespace esys_trilinos {

/// Scalar types
typedef escript::DataTypes::real_t  real_t;
typedef escript::DataTypes::cplx_t  cplx_t;

/// Global Ordinal type
#ifdef MANUALLY_SET_GO
    #ifdef SET_GO_INT
        typedef int GO;
    #elif SET_GO_LONG
        typedef long GO;
    #elif SET_GO_LONG_LONG
        typedef long long int GO;
    #elif SET_GO_COMPLEX_DOUBLE
        typedef std::complex<double> GO;
    #elif SET_GO_REALT
        typedef real_t GO;
    #elif SET_GO_CPLXT
        typedef cplx_t GO;
    #endif
#else
typedef escript::DataTypes::index_t GO;
#endif
/// Local Ordinal type
#ifdef MANUALLY_SET_LO
    #ifdef SET_LO_INT
        typedef int LO;
    #elif SET_LO_LONG
        typedef long LO;
    #elif SET_LO_LONG_LONG
        typedef long long int LO;
    #elif SET_LO_COMPLEX_DOUBLE
        typedef std::complex<double> LO;
    #elif SET_LO_REALT
        typedef real_t LO;
    #elif SET_LO_CPLXT
        typedef cplx_t LO;
    #endif
#else
typedef escript::DataTypes::index_t LO;
#endif
/// Kokkos Node type
#ifdef _OPENMP
    typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode NT;
#elif ESYS_HAVE_CUDA
    typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode   NT;
#else
    typedef Tpetra::KokkosCompat::KokkosSerialWrapperNode NT;
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

// Block vector type (Tpetra::BlockVector is standard in Trilinos 16.2+)
template<typename ST> using BlockVectorType =
                               Tpetra::BlockVector<ST,LO,GO,NT>;
typedef BlockVectorType<real_t> RealBlockVector;
typedef BlockVectorType<cplx_t> ComplexBlockVector;


/// converts an MPI communicator to a Teuchos communicator
inline
Teuchos::RCP<const Teuchos::Comm<int> > TeuchosCommFromEsysComm(MPI_Comm comm)
{
#ifdef TRILINOS_HAS_MPI
    // Both escript and Trilinos have MPI - wrap the communicator
    return Teuchos::rcp(new Teuchos::MpiComm<int>(comm));
#else
    // Either escript or Trilinos (or both) lack MPI - use default serial comm
    (void)comm; // suppress unused parameter warning
    return Teuchos::DefaultComm<int>::getComm();
#endif
}

} // namespace esys_trilinos

#endif // __ESYS_TRILINOSWRAP_TYPES_H__

