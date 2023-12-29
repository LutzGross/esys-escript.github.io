SET(LIB_REQUIRED_DEP_PACKAGES Teuchos Xpetra KokkosCore KokkosContainers KokkosKernels)
SET(LIB_OPTIONAL_DEP_PACKAGES Amesos Amesos2 AvatarT Belos Epetra EpetraExt Teko
                              Ifpack Ifpack2 Intrepid2 ML Tpetra
                              Zoltan Zoltan2Core Stratimikos Thyra ThyraTpetraAdapters
                              Isorropia)
SET(TEST_REQUIRED_DEP_PACKAGES Galeri)
SET(TEST_OPTIONAL_DEP_PACKAGES AztecOO Pamgen)
SET(LIB_REQUIRED_DEP_TPLS BLAS LAPACK)
SET(LIB_OPTIONAL_DEP_TPLS Boost MATLAB AmgX CGAL ViennaCL MKL Avatar CUSPARSE MAGMASparse mlpack)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS HYPRE PETSC)
