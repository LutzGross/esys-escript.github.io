
# for TpetraTestingUtilities.hpp
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

ADD_SUBDIRECTORIES(
  BasicPerfTest
  Behavior
  Blas
  Block
  BugTests
  Comm
  Core
  CrsGraph
  CrsGraphTransposer
  CrsMatrix
  Directory
  Distributor
  FECrsGraph
  FECrsMatrix
  FEMultiVector
  HashTable
  ImportExport
  ImportExport2
  inout
  KokkosIntegration
  Map
  MatrixMatrix
  Merge
  MPIAdvance
  MultiVector
  Node
  Operator
  PerformanceCGSolve
  Sort
  Spaces
  Utils
  RowMatrixTransposer
  Tsqr
  )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Epetra)
IF (${PACKAGE_NAME}_ENABLE_Epetra)
  ADD_SUBDIRECTORY(EpetraRowMatrix)
ENDIF()
