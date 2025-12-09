#ifndef TPETRACORE_CONFIG_H
#define TPETRACORE_CONFIG_H
/* CMake uses this file to generate TpetraCore_config.h automatically */

/* Define whether Tpetra enables deprecated code at compile time. */
#define TPETRA_ENABLE_DEPRECATED_CODE

/* BlockCrs has default layout right and it can be compile time changed with the following cmake option */
/* #undef TPETRA_ENABLE_BLOCKCRS_LITTLEBLOCK_LAYOUTLEFT */

/* Define if want to build with epetra enabled */
#define HAVE_TPETRACORE_EPETRA
#ifdef HAVE_TPETRACORE_EPETRA
#  define HAVE_TPETRA_EPETRA
#endif // HAVE_TPETRACORE_EPETRA

/* Define if want to build teuchos-debug */
/* #undef HAVE_TPETRA_DEBUG */

/* #undef HAVE_TPETRA_ENABLE_SS_TESTING */

/* Define if we have MPI */
/* #undef HAVE_TPETRACORE_MPI */
#ifdef HAVE_TPETRACORE_MPI
#  define HAVE_TPETRA_MPI
#endif // HAVE_TPETRACORE_MPI

/* Determine if the CUDA TPL is enabled */
/* #undef HAVE_TPETRACORE_CUDA */

/* Determine if we have the quadmath TPL */
/* #undef HAVE_TPETRACORE_QUADMATH */

/* Determine if we have the mpi_advance TPL */
/* #undef HAVE_TPETRACORE_MPI_ADVANCE */

/* Determine if we have QD */
/* #undef HAVE_TPETRACORE_QD */
#ifdef HAVE_TPETRACORE_QD
#  define HAVE_TPETRA_QD
#endif // HAVE_TPETRACORE_QD

/* Define if want to build tpetra-throw_warnings */
/* #undef HAVE_TPETRA_THROW_WARNINGS */

/* Define if want to build tpetra-print_warnings */
/* #undef HAVE_TPETRA_PRINT_WARNINGS */

/* Define if want to build tpetra-print_efficiency_warnings */
/* #undef HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS */

/* Define if want to build tpetra-throw_abuse_warnings */
/* #undef HAVE_TPETRA_THROW_ABUSE_WARNINGS */

/* Define if want to build tpetra-print_abuse_warnings */
/* #undef HAVE_TPETRA_PRINT_ABUSE_WARNINGS */

/* Define when building Tpetra with the TpetraTSQR subpackage enabled.
   The HAVE_TPETRA_TSQR macro tells you whether it is safe to use TSQR
   in Tpetra (and downstream packages).  TSQR is enabled by default in
   Tpetra if KokkosTSQR is enabled, but users may turn off TSQR by
   setting TpetraCore_ENABLE_TSQR to OFF. */
#define HAVE_TPETRACORE_TPETRATSQR
#ifdef HAVE_TPETRACORE_TPETRATSQR
#  define HAVE_TPETRA_KOKKOSTSQR
#  define HAVE_TPETRA_TPETRATSQR
#endif // HAVE_TPETRACORE_TPETRATSQR

/* Define when building TpetraCore with TSQR enabled.  TSQR is enabled
   by default in TpetraCore if the TpetraTSQR subpackage is enabled.
   Users may turn off TSQR even if TpetraTSQR is enabled, by setting
   TpetraCore_ENABLE_TSQR to OFF. */
#define HAVE_TPETRA_TSQR

/* Define when enabling the Murmur hash function in Tpetra */
/* #undef TPETRA_USE_MURMUR_HASH */

/* Tpetra timers */
/* #undef HAVE_TPETRA_MMM_TIMINGS */
/* #undef HAVE_TPETRA_DISTRIBUTOR_TIMINGS */

/* Define when enabling Kokkos in Tpetra */
#define HAVE_TPETRACORE_KOKKOS
#ifdef HAVE_TPETRACORE_KOKKOS
// For backwards compatibility
#  define HAVE_TPETRA_KOKKOS
#endif // HAVE_TPETRACORE_KOKKOS

/* Define when enabling KokkosCompat in TpetraCore */
#define HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
#ifdef HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
// For backwards compatibility
#  define HAVE_TPETRACORE_KOKKOSCOMPAT
#  define HAVE_TPETRA_TEUCHOSKOKKOSCOMPAT
#  define HAVE_TPETRA_KOKKOSCOMPAT
#endif // HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT

/* Define when enabling KokkosComm in TpetraCore */
#define HAVE_TPETRACORE_TEUCHOSKOKKOSCOMM
#ifdef HAVE_TPETRACORE_TEUCHOSKOKKOSCOMM
// For backwards compatibility
#  define HAVE_TPETRACORE_KOKKOSMPICOMM
#  define HAVE_TPETRA_KOKKOSMPICOMM
#endif // HAVE_TPETRACORE_TEUCHOSKOKKOSCOMM

/* Define when enabling TeuchosNumerics in TpetraCore */
#define HAVE_TPETRACORE_TEUCHOSNUMERICS

/* Define when enabling Kokkos::View DistObject in Tpetra */
/* #undef TPETRA_ENABLE_KOKKOS_DISTOBJECT */

/* Define if user requested explicit instantiation of classes into libtpetra */
#define HAVE_TPETRA_EXPLICIT_INSTANTIATION

/* Define when Tpetra assumes that the MPI implementation is GPU aware */
/* #undef TPETRA_ASSUME_GPU_AWARE_MPI */

/* List of enabled (LocalOrdinal, GlobalOrdinal) pairs */
#define HAVE_TPETRA_INST_INT_INT
/* #undef HAVE_TPETRA_INST_INT_LONG */
/* #undef HAVE_TPETRA_INST_INT_LONG_LONG */
/* #undef HAVE_TPETRA_INST_INT_UNSIGNED */
/* #undef HAVE_TPETRA_INST_INT_UNSIGNED_LONG */

/* List of enabled Scalar types */
/* #undef HAVE_TPETRA_INST_FLOAT */
#define HAVE_TPETRA_INST_DOUBLE
/* #undef HAVE_TPETRA_INST_LONG_DOUBLE */
/* #undef HAVE_TPETRA_INST_COMPLEX_FLOAT */
#define HAVE_TPETRA_INST_COMPLEX_DOUBLE
/* #undef HAVE_TPETRA_INST_DD_REAL */
/* #undef HAVE_TPETRA_INST_QD_REAL */
/* #undef HAVE_TPETRA_INST_FLOAT128 */

/* List of enabled Node types */
/* #undef HAVE_TPETRA_INST_SERIALCLASSIC */
#define HAVE_TPETRA_INST_SERIAL
/* #undef HAVE_TPETRA_INST_PTHREAD */
#define HAVE_TPETRA_INST_OPENMP
/* #undef HAVE_TPETRA_INST_CUDA */
/* #undef HAVE_TPETRA_INST_HIP */
/* #undef HAVE_TPETRA_INST_SYCL */

#define HAVE_TPETRA_INT_INT
/* #undef HAVE_TPETRA_INT_LONG */
/* #undef HAVE_TPETRA_INT_LONG_LONG */
/* #undef HAVE_TPETRA_INT_UNSIGNED */
/* #undef HAVE_TPETRA_INT_UNSIGNED_LONG */

/* #undef HAVE_TPETRA_FLOAT */
#define HAVE_TPETRA_DOUBLE
/* #undef HAVE_TPETRA_LONG_DOUBLE */
/* #undef HAVE_TPETRA_COMPLEX_FLOAT */
#define HAVE_TPETRA_COMPLEX_DOUBLE
/* #undef HAVE_TPETRA_DD_REAL */
/* #undef HAVE_TPETRA_QD_REAL */

/* #undef HAVE_TPETRA_SERIALCLASSIC */
#define HAVE_TPETRA_SERIAL
/* #undef HAVE_TPETRA_PTHREAD */
#define HAVE_TPETRA_OPENMP
/* #undef HAVE_TPETRA_CUDA */
/* #undef HAVE_TPETRA_HIP */
/* #undef HAVE_TPETRA_SYCL */

/* #undef HAVE_TPETRA_REDUCED_ETI */

/* #undef HAVE_TPETRA_THREADED_MKL */

/* #undef HAVE_TPETRA_BCRS_DO_POINT_IMPORT */

/* Define the Fortran name mangling to be used for the BLAS */
#ifndef TPETRACORE_F77_BLAS_MANGLE
#  define TPETRACORE_F77_BLAS_MANGLE(name,NAME) name ## _
#endif


/* #undef HAVE_TPETRA_SHARED_ALLOCS */

/* #undef HAVE_TPETRA_EXPERIMENTAL */

#ifndef TPETRA_DEPRECATED
#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define TPETRA_DEPRECATED  __attribute__((__deprecated__))
#  else
#    define TPETRA_DEPRECATED
#  endif
#endif

#ifndef TPETRA_DEPRECATED_MSG
#  if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5))
#    define TPETRA_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__ (#MSG) ))
#  elif (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define TPETRA_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__))
#  else
#    define TPETRA_DEPRECATED_MSG(MSG)
#  endif
#endif


#define TPETRA_DEPRECATED_DECLARATIONS

#endif // TPETRACORE_CONFIG_H
