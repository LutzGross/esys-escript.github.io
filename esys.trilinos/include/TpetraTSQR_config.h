#ifndef TPETRATSQR_CONFIG_H
#define TPETRATSQR_CONFIG_H

/* Define if building TSQR with std::complex<T> support */
#define HAVE_TPETRATSQR_COMPLEX
#ifdef HAVE_TPETRATSQR_COMPLEX
   /* For backwards compatibility */
#  define HAVE_KOKKOSTSQR_COMPLEX HAVE_TPETRATSQR_COMPLEX
#endif

/* Define if TSQR supports the CUBLAS TPL */
/* #undef HAVE_TPETRATSQR_CUBLAS */

/* Define if TSQR supports the CUSOLVER TPL */
/* #undef HAVE_TPETRATSQR_CUSOLVER */

#endif // TPETRATSQR_CONFIG_H
