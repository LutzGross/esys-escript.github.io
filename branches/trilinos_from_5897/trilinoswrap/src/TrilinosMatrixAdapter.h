
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

#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>

namespace escript {
class SolverBuddy;
}

namespace esys_trilinos {

/// Scalar Type
typedef double  ST;
/// Global Ordinal Type
typedef index_t GO;
/// Local Ordinal Type
typedef index_t LO;

typedef std::vector<GO> IndexVector;
typedef Tpetra::CrsGraph<LO,GO> GraphType;
typedef Tpetra::CrsMatrix<ST,LO,GO> MatrixType;
typedef MatrixType::map_type MapType;
typedef Tpetra::MultiVector<ST,LO,GO> VectorType;
typedef Tpetra::Import<LO,GO> ImportType;

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

    void add(const IndexVector& rowIndex, const std::vector<double>& array);

    inline int getBlockSize() const { return getRowBlockSize(); }

private:
    virtual void setToSolution(escript::Data& out, escript::Data& in,
                               boost::python::object& options) const;

    virtual void ypAx(escript::Data& y, escript::Data& x) const;

    esysUtils::JMPI m_mpiInfo;
    Teuchos::RCP<MatrixType> mat;
    Teuchos::RCP<const ImportType> importer;
};

} // namespace esys_trilinos

#endif // __TRILINOSMATRIXADAPTER_H__

