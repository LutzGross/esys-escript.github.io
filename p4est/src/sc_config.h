#ifndef _CONFIG_SC_CONFIG_H
#define _CONFIG_SC_CONFIG_H 1
 
/* config/sc_config.h. Generated automatically at end of configure. */
/* config/sc_config_autotools.h.  Generated from sc_config_autotools.h.in by configure.  */
/* config/sc_config_autotools.h.in.  Generated from configure.ac by autoheader.  */

/* DEPRECATED (use SC_ENABLE_DEBUG instead) */
/* #undef DEBUG */

/* enable debug mode (assertions and extra checks) */
/* #undef ENABLE_DEBUG */

/* Undefine if: disable tests that use file i/o functions */
#ifndef SC_ENABLE_FILE_CHECKS
#define SC_ENABLE_FILE_CHECKS 1
#endif

/* Undefine if: while the default alignment is sizeof (void *), this switch
   will choose the standard system malloc. For custom alignment use
   --enable-memalign=<bytes> */
#ifndef SC_ENABLE_MEMALIGN
#define SC_ENABLE_MEMALIGN 1
#endif

/* Define to 1 if we are using MPI */
#ifndef SC_ENABLE_MPI
#define SC_ENABLE_MPI 1
#endif

/* Define to 1 if we can use MPI_COMM_TYPE_SHARED */
#ifndef SC_ENABLE_MPICOMMSHARED
#define SC_ENABLE_MPICOMMSHARED 1
#endif

/* Define to 1 if we are using MPI I/O */
#ifndef SC_ENABLE_MPIIO
#define SC_ENABLE_MPIIO 1
#endif

/* Define to 1 if we can use MPI split nodes and shared memory */
#ifndef SC_ENABLE_MPISHARED
#define SC_ENABLE_MPISHARED 1
#endif

/* Define to 1 if we are using MPI_Init_thread */
#ifndef SC_ENABLE_MPITHREAD
#define SC_ENABLE_MPITHREAD 1
#endif

/* Define to 1 if we can use MPI_Win_allocate_shared */
#ifndef SC_ENABLE_MPIWINSHARED
#define SC_ENABLE_MPIWINSHARED 1
#endif

/* enable POSIX threads: Using --enable-pthread without arguments does not
   specify any CFLAGS; to supply CFLAGS use --enable-pthread=<PTHREAD_CFLAGS>.
   We check first for linking without any libraries and then with -lpthread;
   to avoid the latter, specify LIBS=<PTHREAD_LIBS> on configure line */
/* #undef ENABLE_PTHREAD */

/* Undefine if: disable non-thread-safe internal debug counters */
#ifndef SC_ENABLE_USE_COUNTERS
#define SC_ENABLE_USE_COUNTERS 1
#endif

/* Undefine if: resize arrays with malloc/copy/free (HISTORIC) */
#ifndef SC_ENABLE_USE_REALLOC
#define SC_ENABLE_USE_REALLOC 1
#endif

/* Development with V4L2 devices works */
#ifndef SC_ENABLE_V4L2
#define SC_ENABLE_V4L2 1
#endif

/* Enable valgrind in executing tests */
/* #undef ENABLE_VALGRIND */

/* DEPRECATED (use SC_ENABLE_FILE_CHECKS instead) */
#ifndef SC_FILE_CHECKS
#define SC_FILE_CHECKS 1
#endif

/* Define to 1 if we have MPI_Aint_diff */
#ifndef SC_HAVE_AINT_DIFF
#define SC_HAVE_AINT_DIFF 1
#endif

/* Define to 1 if you have the `aligned_alloc' function. */
#ifndef SC_HAVE_ALIGNED_ALLOC
#define SC_HAVE_ALIGNED_ALLOC 1
#endif

/* Define to 1 if you have the `backtrace' function. */
#ifndef SC_HAVE_BACKTRACE
#define SC_HAVE_BACKTRACE 1
#endif

/* Define to 1 if you have the `backtrace_symbols' function. */
#ifndef SC_HAVE_BACKTRACE_SYMBOLS
#define SC_HAVE_BACKTRACE_SYMBOLS 1
#endif

/* Define to 1 if you have the `basename' function. */
#ifndef SC_HAVE_BASENAME
#define SC_HAVE_BASENAME 1
#endif

/* Define to 1 if qsort_r conforms to BSD standard */
/* #undef HAVE_BSD_QSORT_R */

/* Define to 1 if you have the `dirname' function. */
#ifndef SC_HAVE_DIRNAME
#define SC_HAVE_DIRNAME 1
#endif

/* Define to 1 if you have the <dlfcn.h> header file. */
#ifndef SC_HAVE_DLFCN_H
#define SC_HAVE_DLFCN_H 1
#endif

/* Define to 1 if you have the <execinfo.h> header file. */
#ifndef SC_HAVE_EXECINFO_H
#define SC_HAVE_EXECINFO_H 1
#endif

/* Define to 1 if you have the <fcntl.h> header file. */
#ifndef SC_HAVE_FCNTL_H
#define SC_HAVE_FCNTL_H 1
#endif

/* Define to 1 if you have the `fsync' function. */
#ifndef SC_HAVE_FSYNC
#define SC_HAVE_FSYNC 1
#endif

/* Define to 1 if you have the `gettimeofday' function. */
#ifndef SC_HAVE_GETTIMEOFDAY
#define SC_HAVE_GETTIMEOFDAY 1
#endif

/* Define to 1 if qsort_r conforms to GNU standard */
#ifndef SC_HAVE_GNU_QSORT_R
#define SC_HAVE_GNU_QSORT_R 1
#endif

/* Define to 1 if you have the <inttypes.h> header file. */
#ifndef SC_HAVE_INTTYPES_H
#define SC_HAVE_INTTYPES_H 1
#endif

/* Define to 1 if json_integer and json_real link */
/* #undef HAVE_JSON */

/* Define to 1 if you have the <libgen.h> header file. */
#ifndef SC_HAVE_LIBGEN_H
#define SC_HAVE_LIBGEN_H 1
#endif

/* Define to 1 if you have the <linux/version.h> header file. */
#ifndef SC_HAVE_LINUX_VERSION_H
#define SC_HAVE_LINUX_VERSION_H 1
#endif

/* Define to 1 if you have the <linux/videodev2.h> header file. */
#ifndef SC_HAVE_LINUX_VIDEODEV2_H
#define SC_HAVE_LINUX_VIDEODEV2_H 1
#endif

/* Have we found function pthread_create. */
/* #undef HAVE_LPTHREAD */

/* Define to 1 if sqrt links successfully */
#ifndef SC_HAVE_MATH
#define SC_HAVE_MATH 1
#endif

/* Define to 1 if you have the <memory.h> header file. */
#ifndef SC_HAVE_MEMORY_H
#define SC_HAVE_MEMORY_H 1
#endif

/* Define to 1 if we have MPI_INT8_T */
#ifndef SC_HAVE_MPI_INT8_T
#define SC_HAVE_MPI_INT8_T 1
#endif

/* Define to 1 if we have MPI_SIGNED_CHAR */
#ifndef SC_HAVE_MPI_SIGNED_CHAR
#define SC_HAVE_MPI_SIGNED_CHAR 1
#endif

/* Define to 1 if we have MPI_UNSIGNED_LONG_LONG */
#ifndef SC_HAVE_MPI_UNSIGNED_LONG_LONG
#define SC_HAVE_MPI_UNSIGNED_LONG_LONG 1
#endif

/* Define to 1 if you have the `posix_memalign' function. */
#ifndef SC_HAVE_POSIX_MEMALIGN
#define SC_HAVE_POSIX_MEMALIGN 1
#endif

/* Define to 1 if you have the `qsort_r' function. */
#ifndef SC_HAVE_QSORT_R
#define SC_HAVE_QSORT_R 1
#endif

/* Define to 1 if you have the <signal.h> header file. */
#ifndef SC_HAVE_SIGNAL_H
#define SC_HAVE_SIGNAL_H 1
#endif

/* Define to 1 if you have the <stdint.h> header file. */
#ifndef SC_HAVE_STDINT_H
#define SC_HAVE_STDINT_H 1
#endif

/* Define to 1 if you have the <stdlib.h> header file. */
#ifndef SC_HAVE_STDLIB_H
#define SC_HAVE_STDLIB_H 1
#endif

/* Define to 1 if you have the <strings.h> header file. */
#ifndef SC_HAVE_STRINGS_H
#define SC_HAVE_STRINGS_H 1
#endif

/* Define to 1 if you have the <string.h> header file. */
#ifndef SC_HAVE_STRING_H
#define SC_HAVE_STRING_H 1
#endif

/* Define to 1 if you have the `strtok_r' function. */
#ifndef SC_HAVE_STRTOK_R
#define SC_HAVE_STRTOK_R 1
#endif

/* Define to 1 if you have the `strtol' function. */
#ifndef SC_HAVE_STRTOL
#define SC_HAVE_STRTOL 1
#endif

/* Define to 1 if you have the `strtoll' function. */
#ifndef SC_HAVE_STRTOLL
#define SC_HAVE_STRTOLL 1
#endif

/* Define to 1 if you have the <sys/ioctl.h> header file. */
#ifndef SC_HAVE_SYS_IOCTL_H
#define SC_HAVE_SYS_IOCTL_H 1
#endif

/* Define to 1 if you have the <sys/select.h> header file. */
#ifndef SC_HAVE_SYS_SELECT_H
#define SC_HAVE_SYS_SELECT_H 1
#endif

/* Define to 1 if you have the <sys/stat.h> header file. */
#ifndef SC_HAVE_SYS_STAT_H
#define SC_HAVE_SYS_STAT_H 1
#endif

/* Define to 1 if you have the <sys/time.h> header file. */
#ifndef SC_HAVE_SYS_TIME_H
#define SC_HAVE_SYS_TIME_H 1
#endif

/* Define to 1 if you have the <sys/types.h> header file. */
#ifndef SC_HAVE_SYS_TYPES_H
#define SC_HAVE_SYS_TYPES_H 1
#endif

/* Define to 1 if you have the <time.h> header file. */
#ifndef SC_HAVE_TIME_H
#define SC_HAVE_TIME_H 1
#endif

/* Define to 1 if you have the <unistd.h> header file. */
#ifndef SC_HAVE_UNISTD_H
#define SC_HAVE_UNISTD_H 1
#endif

/* Define to 1 if zlib's adler32_combine links */
#ifndef SC_HAVE_ZLIB
#define SC_HAVE_ZLIB 1
#endif

/* minimal log priority */
/* #undef LOG_PRIORITY */

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#ifndef SC_LT_OBJDIR
#define SC_LT_OBJDIR ".libs/"
#endif

/* DEPRECATED (use SC_ENABLE_MEMALIGN instead) */
#ifndef SC_MEMALIGN
#define SC_MEMALIGN 1
#endif

/* desired alignment of allocations in bytes */
#ifndef SC_MEMALIGN_BYTES
#define SC_MEMALIGN_BYTES (8)
#endif

/* DEPRECATED (use SC_ENABLE_MPI instead) */
#ifndef SC_MPI
#define SC_MPI 1
#endif

/* DEPRECATED (use SC_ENABLE_MPIIO instead) */
#ifndef SC_MPIIO
#define SC_MPIIO 1
#endif

/* Define to 1 if your C compiler doesn't accept -c and -o together. */
/* #undef NO_MINUS_C_MINUS_O */

/* Name of package */
#ifndef SC_PACKAGE
#define SC_PACKAGE "libsc"
#endif

/* Define to the address where bug reports for this package should be sent. */
#ifndef SC_PACKAGE_BUGREPORT
#define SC_PACKAGE_BUGREPORT "p4est@ins.uni-bonn.de"
#endif

/* Define to the full name of this package. */
#ifndef SC_PACKAGE_NAME
#define SC_PACKAGE_NAME "libsc"
#endif

/* Define to the full name and version of this package. */
#ifndef SC_PACKAGE_STRING
#define SC_PACKAGE_STRING "libsc 2.8.7"
#endif

/* Define to the one symbol short name of this package. */
#ifndef SC_PACKAGE_TARNAME
#define SC_PACKAGE_TARNAME "libsc"
#endif

/* Define to the home page for this package. */
#ifndef SC_PACKAGE_URL
#define SC_PACKAGE_URL ""
#endif

/* Define to the version of this package. */
#ifndef SC_PACKAGE_VERSION
#define SC_PACKAGE_VERSION "2.8.7"
#endif

/* DEPRECATED (use SC_WITH_PAPI instead) */
/* #undef PAPI */

/* Use builtin getopt */
/* #undef PROVIDE_GETOPT */

/* DEPRECATED (use SC_ENABLE_PTHREAD instead) */
/* #undef PTHREAD */

/* Define to 1 if you have the ANSI C header files. */
#ifndef SC_STDC_HEADERS
#define SC_STDC_HEADERS 1
#endif

/* DEPRECATED (use SC_ENABLE_USE_COUNTERS instead) */
#ifndef SC_USE_COUNTERS
#define SC_USE_COUNTERS 1
#endif

/* DEPRECATED (use SC_ENABLE_USE_REALLOC instead) */
#ifndef SC_USE_REALLOC
#define SC_USE_REALLOC 1
#endif

/* Define to 1 if using autoconf build */
#ifndef SC_USING_AUTOCONF
#define SC_USING_AUTOCONF 1
#endif

/* Version number of package */
#ifndef SC_VERSION
#define SC_VERSION "2.8.7"
#endif

/* Package major version */
#ifndef SC_VERSION_MAJOR
#define SC_VERSION_MAJOR 2
#endif

/* Package minor version */
#ifndef SC_VERSION_MINOR
#define SC_VERSION_MINOR 8
#endif

/* Package point version */
#ifndef SC_VERSION_POINT
#define SC_VERSION_POINT 7
#endif

/* enable Flop counting with papi */
/* #undef WITH_PAPI */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to the equivalent of the C99 'restrict' keyword, or to
   nothing if this is not supported.  Do not define if restrict is
   supported directly.  */
#ifndef _sc_restrict
#define _sc_restrict __restrict
#endif
/* Work around a bug in Sun C++: it does not support _Restrict or
   __restrict__, even though the corresponding Sun C compiler ends up with
   "#define restrict _Restrict" or "#define restrict __restrict__" in the
   previous line.  Perhaps some future version of Sun C++ will work with
   restrict; if so, hopefully it defines __RESTRICT like Sun C does.  */
#if defined __SUNPRO_CC && !defined __RESTRICT
# define _Restrict
# define __restrict__
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef ssize_t */
 
/* once: _CONFIG_SC_CONFIG_H */
#endif
