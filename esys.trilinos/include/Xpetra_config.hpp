#ifndef XPETRA_CONFIG_HPP
#define XPETRA_CONFIG_HPP

/* #undef HAVE_XPETRA_DEBUG */

#define HAVE_XPETRA_TPETRA

/* Define whether Xpetra enables deprecated code at compile time. */
#define XPETRA_ENABLE_DEPRECATED_CODE

#define HAVE_XPETRA_EPETRA

#define HAVE_XPETRA_EPETRAEXT

/* #undef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES */

#define XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES

/* #undef HAVE_XPETRA_EXPERIMENTAL */

#define HAVE_XPETRA_EXPLICIT_INSTANTIATION

/* #undef HAVE_XPETRA_KOKKOSCLASSIC */

#define HAVE_XPETRA_KOKKOS

#define HAVE_XPETRA_KOKKOSKERNELS

/* #undef HAVE_XPETRA_PROFILING */

/* #undef HAVE_XPETRA_THYRA */

#ifndef XPETRA_DEPRECATED
#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define XPETRA_DEPRECATED  __attribute__((__deprecated__))
#  else
#    define XPETRA_DEPRECATED
#  endif
#endif

#ifndef XPETRA_DEPRECATED_MSG
#  if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5))
#    define XPETRA_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__ (#MSG) ))
#  elif (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define XPETRA_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__))
#  else
#    define XPETRA_DEPRECATED_MSG(MSG)
#  endif
#endif


#endif // XPETRA_CONFIG_HPP
