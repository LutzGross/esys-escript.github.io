#ifndef MUELU_CONFIG_HPP
#define MUELU_CONFIG_HPP

#ifndef F77_BLAS_MANGLE
 #define F77_BLAS_MANGLE(name,NAME) name ## _
#endif

/* Options */

/* #undef HAVE_MPI */

/* #undef HAVE_MUELU_DEBUG */

#define HAVE_MUELU_EXPLICIT_INSTANTIATION

/* #undef HAVE_MUELU_GOOGLE_PERFTOOLS */

/* #undef HAVE_MUELU_EXPERIMENTAL */

/* #undef HAVE_MUELU_ADDITIVE_VARIANT */

/* #undef HAVE_MUELU_BOOST_FOR_REAL */

#define HAVE_MUELU_DEPRECATED_CODE

/* #undef HAVE_MUELU_DEPRECATED_TESTS */

/* Optional Dependencies */

#define HAVE_MUELU_AMESOS

#define HAVE_MUELU_AMESOS2

#define HAVE_MUELU_AZTECOO

#define HAVE_MUELU_BELOS

#define HAVE_MUELU_EPETRA

#define HAVE_MUELU_EPETRAEXT

#define HAVE_MUELU_GALERI

#define HAVE_MUELU_IFPACK

#define HAVE_MUELU_IFPACK2

/* #undef HAVE_MUELU_INTREPID2 */

/* #undef HAVE_MUELU_ISORROPIA */

/* #undef HAVE_MUELU_KOKKOSCLASSIC */

#define HAVE_MUELU_ML

/* #undef HAVE_MUELU_PAMGEN */

/* #undef HAVE_MUELU_TEKO */

/* #undef HAVE_MUELU_THYRA */

/* #undef HAVE_MUELU_THYRATPETRAADAPTERS */

#define HAVE_MUELU_TPETRA

/* Whether Tpetra is enabled with LocalOrdinal = int and GlobalOrdinal = int */
#define HAVE_MUELU_TPETRA_INST_INT_INT

/* #undef HAVE_MUELU_STRATIMIKOS */

/* #undef HAVE_MUELU_ZOLTAN */

/* #undef HAVE_MUELU_ZOLTAN2CORE */
#ifdef HAVE_MUELU_ZOLTAN2CORE
#  define HAVE_MUELU_ZOLTAN2
#endif


/* Optional TPL */

#define HAVE_MUELU_BOOST

/* #undef HAVE_MUELU_MATLAB */

/* #undef HAVE_MUELU_AMGX */

/* #undef HAVE_MUELU_VIENNACL */

/* #undef HAVE_MUELU_MKL */

/* Special logic for avatar */
/* #undef HAVE_MUELU_AVATART */
/* #undef HAVE_MUELU_AVATAR */
#if defined(HAVE_MUELU_AVATART)
 #define HAVE_MUELU_AVATAR
#endif


/* #undef HAVE_MUELU_MLPACK */

/* #undef HAVE_MUELU_CUSPARSE */

/* #undef HAVE_MUELU_MAGMASPARSE */


/* #undef HAVE_MUELU_HYPRE */

/* #undef HAVE_MUELU_PETSC */

/* Flags for active Tpetra nodes if Tpetra is enabled */
#define HAVE_MUELU_SERIAL
#define HAVE_MUELU_OPENMP
/* #undef HAVE_MUELU_CUDA */
/* #undef HAVE_MUELU_HIP */
/* #undef HAVE_MUELU_SYCL */

/* #undef HAVE_MUELU_DEFAULT_GO_LONGLONG */
/* #undef HAVE_MUELU_DEFAULT_GO_LONG */

/*
 If deprecated warnings are on, and the compiler supports them, then
 define MUELU_DEPRECATED to emit deprecated warnings.  Otherwise,
 give it an empty definition.
*/
#ifndef MUELU_DEPRECATED
#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define MUELU_DEPRECATED  __attribute__((__deprecated__))
#  else
#    define MUELU_DEPRECATED
#  endif
#endif

#ifndef MUELU_DEPRECATED_MSG
#  if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5))
#    define MUELU_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__ (#MSG) ))
#  elif (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define MUELU_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__))
#  else
#    define MUELU_DEPRECATED_MSG(MSG)
#  endif
#endif


#endif /* MUELU_CONFIG_HPP */
