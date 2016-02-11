
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

#ifndef __TRILINOSMATRIXADAPTER_H__
#define __TRILINOSMATRIXADAPTER_H__

#include <escript/AbstractSystemMatrix.h>
#include <escript/FunctionSpace.h>

#include <Amesos2_Solver_decl.hpp>
#include <BelosSolverManager.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>

namespace escript {
    class SolverBuddy;
}

namespace esys_trilinos {

/// Scalar type
typedef double  ST;
/// Global Ordinal type
typedef index_t GO;
/// Local Ordinal type
typedef index_t LO;
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

class TrilinosMatrixAdapter : public escript::AbstractSystemMatrix
{
public:
    /**
       \brief
       Creates a new Trilinos CRS matrix adapter using a compatible
       fill-complete Trilinos matrix graph.
    */
    TrilinosMatrixAdapter(esysUtils::JMPI mpiInfo, int blocksize,
                          const escript::FunctionSpace& fs,
                          const_TrilinosGraph_ptr graph);

    virtual ~TrilinosMatrixAdapter() {}

    virtual void nullifyRowsAndCols(escript::Data& row_q,
                                    escript::Data& col_q,
                                    double mdv);  
  
    virtual void saveMM(const std::string& filename) const;
    virtual void saveHB(const std::string& filename) const;
    virtual void resetValues();

    /// notifies the matrix that changes are about to happen.
    void resumeFill() { mat->resumeFill(); }

    /// notifies the matrix that a set of changes has occured.
    void fillComplete(bool localOnly = true);

    void add(const std::vector<LO>& rowIndex, const std::vector<double>& array);

    inline int getBlockSize() const { return getRowBlockSize(); }

private:
    Teuchos::RCP<OpType> createPreconditioner(
                                        const escript::SolverBuddy& sb) const;

    Teuchos::RCP<SolverType> createSolver(const escript::SolverBuddy& sb) const;

    Teuchos::RCP<DirectSolverType> createDirectSolver(
                                      const escript::SolverBuddy& sb,
                                      Teuchos::RCP<VectorType> X,
                                      Teuchos::RCP<const VectorType> B) const;

    virtual void setToSolution(escript::Data& out, escript::Data& in,
                               boost::python::object& options) const;

    virtual void ypAx(escript::Data& y, escript::Data& x) const;

    esysUtils::JMPI m_mpiInfo;
    Teuchos::RCP<MatrixType> mat;
    Teuchos::RCP<const ImportType> importer;
};

} // namespace esys_trilinos

#endif // __TRILINOSMATRIXADAPTER_H__

