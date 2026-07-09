#ifndef _CONFIG_P_EST_CONFIG_H
#define _CONFIG_P_EST_CONFIG_H 1
 
/* config/p4est_config.h. Generated automatically at end of configure. */
/* config/p4est_autotools_config.h.  Generated from p4est_autotools_config.h.in by configure.  */
/* config/p4est_autotools_config.h.in.  Generated from configure.ac by autoheader.  */

/* DEPRECATED (use P4EST_ENABLE_BUILD_2D instead) */
#ifndef P4EST_BUILD_2D
#define P4EST_BUILD_2D 1
#endif

/* DEPRECATED (use P4EST_ENABLE_BUILD_3D instead) */
#ifndef P4EST_BUILD_3D
#define P4EST_BUILD_3D 1
#endif

/* DEPRECATED (use P4EST_ENABLE_BUILD_P6EST instead) */
#ifndef P4EST_BUILD_P6EST
#define P4EST_BUILD_P6EST 1
#endif

/* DEPRECATED (use P4EST_ENABLE_DEBUG instead) */
/* #undef DEBUG */

/* Undefine if: disable the 2D library */
#ifndef P4EST_ENABLE_BUILD_2D
#define P4EST_ENABLE_BUILD_2D 1
#endif

/* Undefine if: disable the 3D library */
#ifndef P4EST_ENABLE_BUILD_3D
#define P4EST_ENABLE_BUILD_3D 1
#endif

/* Undefine if: disable hybrid 2D+1D p6est library */
#ifndef P4EST_ENABLE_BUILD_P6EST
#define P4EST_ENABLE_BUILD_P6EST 1
#endif

/* enable debug mode (assertions and extra checks) */
/* #undef ENABLE_DEBUG */

/* Undefine if: disable tests that use file i/o functions */
#ifndef P4EST_ENABLE_FILE_CHECKS
#define P4EST_ENABLE_FILE_CHECKS 1
#endif

/* use depreacted data file format */
/* #undef ENABLE_FILE_DEPRECATED */

/* Undefine if: while the default alignment is sizeof (void *), this switch
   will choose the standard system malloc. For custom alignment use
   --enable-memalign=<bytes> */
#ifndef P4EST_ENABLE_MEMALIGN
#define P4EST_ENABLE_MEMALIGN 1
#endif

/* Define to 1 if we are using MPI */
#ifndef P4EST_ENABLE_MPI
#define P4EST_ENABLE_MPI 1
#endif

/* Define to 1 if we can use MPI_COMM_TYPE_SHARED */
#ifndef P4EST_ENABLE_MPICOMMSHARED
#define P4EST_ENABLE_MPICOMMSHARED 1
#endif

/* Define to 1 if we are using MPI I/O */
#ifndef P4EST_ENABLE_MPIIO
#define P4EST_ENABLE_MPIIO 1
#endif

/* Define to 1 if we can use MPI split nodes and shared memory */
#ifndef P4EST_ENABLE_MPISHARED
#define P4EST_ENABLE_MPISHARED 1
#endif

/* Define to 1 if we are using MPI_Init_thread */
#ifndef P4EST_ENABLE_MPITHREAD
#define P4EST_ENABLE_MPITHREAD 1
#endif

/* Define to 1 if we can use MPI_Win_allocate_shared */
#ifndef P4EST_ENABLE_MPIWINSHARED
#define P4EST_ENABLE_MPIWINSHARED 1
#endif

/* enable POSIX threads: Using --enable-pthread without arguments does not
   specify any CFLAGS; to supply CFLAGS use --enable-pthread=<PTHREAD_CFLAGS>.
   We check first for linking without any libraries and then with -lpthread;
   to avoid the latter, specify LIBS=<PTHREAD_LIBS> on configure line */
/* #undef ENABLE_PTHREAD */

/* Development with V4L2 devices works */
/* #undef ENABLE_V4L2 */

/* Enable valgrind in executing tests */
/* #undef ENABLE_VALGRIND */

/* Undefine if: write vtk ascii file data */
#ifndef P4EST_ENABLE_VTK_BINARY
#define P4EST_ENABLE_VTK_BINARY 1
#endif

/* Undefine if: disable zlib compression for vtk binary data */
#ifndef P4EST_ENABLE_VTK_COMPRESSION
#define P4EST_ENABLE_VTK_COMPRESSION 1
#endif

/* use doubles for vtk file data */
/* #undef ENABLE_VTK_DOUBLES */

/* DEPRECATED (use P4EST_ENABLE_FILE_CHECKS instead) */
#ifndef P4EST_FILE_CHECKS
#define P4EST_FILE_CHECKS 1
#endif

/* DEPRECATED (use P4EST_ENABLE_FILE_DEPRECATED instead) */
/* #undef FILE_DEPRECATED */

/* Define to 1 if we have MPI_Aint_diff */
#ifndef P4EST_HAVE_AINT_DIFF
#define P4EST_HAVE_AINT_DIFF 1
#endif

/* Define to 1 if you have the `aligned_alloc' function. */
#ifndef P4EST_HAVE_ALIGNED_ALLOC
#define P4EST_HAVE_ALIGNED_ALLOC 1
#endif

/* Define to 1 if you have the <arpa/inet.h> header file. */
#ifndef P4EST_HAVE_ARPA_INET_H
#define P4EST_HAVE_ARPA_INET_H 1
#endif

/* Define to 1 if qsort_r conforms to BSD standard */
/* #undef HAVE_BSD_QSORT_R */

/* Define to 1 if you have the <dlfcn.h> header file. */
#ifndef P4EST_HAVE_DLFCN_H
#define P4EST_HAVE_DLFCN_H 1
#endif

/* Define to 1 if qsort_r conforms to GNU standard */
#ifndef P4EST_HAVE_GNU_QSORT_R
#define P4EST_HAVE_GNU_QSORT_R 1
#endif

/* Define to 1 if you have the <inttypes.h> header file. */
#ifndef P4EST_HAVE_INTTYPES_H
#define P4EST_HAVE_INTTYPES_H 1
#endif

/* Define to 1 if json_integer and json_real link */
/* #undef HAVE_JSON */

/* Have we found function pthread_create. */
/* #undef HAVE_LPTHREAD */

/* Define to 1 if sqrt links successfully */
#ifndef P4EST_HAVE_MATH
#define P4EST_HAVE_MATH 1
#endif

/* Define to 1 if you have the <memory.h> header file. */
#ifndef P4EST_HAVE_MEMORY_H
#define P4EST_HAVE_MEMORY_H 1
#endif

/* Define to 1 if we have MPI_INT8_T */
#ifndef P4EST_HAVE_MPI_INT8_T
#define P4EST_HAVE_MPI_INT8_T 1
#endif

/* Define to 1 if we have MPI_SIGNED_CHAR */
#ifndef P4EST_HAVE_MPI_SIGNED_CHAR
#define P4EST_HAVE_MPI_SIGNED_CHAR 1
#endif

/* Define to 1 if we have MPI_UNSIGNED_LONG_LONG */
#ifndef P4EST_HAVE_MPI_UNSIGNED_LONG_LONG
#define P4EST_HAVE_MPI_UNSIGNED_LONG_LONG 1
#endif

/* Define to 1 if you have the <netinet/in.h> header file. */
#ifndef P4EST_HAVE_NETINET_IN_H
#define P4EST_HAVE_NETINET_IN_H 1
#endif

/* Define to 1 if you have the `posix_memalign' function. */
#ifndef P4EST_HAVE_POSIX_MEMALIGN
#define P4EST_HAVE_POSIX_MEMALIGN 1
#endif

/* Define to 1 if you have the <stdint.h> header file. */
#ifndef P4EST_HAVE_STDINT_H
#define P4EST_HAVE_STDINT_H 1
#endif

/* Define to 1 if you have the <stdlib.h> header file. */
#ifndef P4EST_HAVE_STDLIB_H
#define P4EST_HAVE_STDLIB_H 1
#endif

/* Define to 1 if you have the <strings.h> header file. */
#ifndef P4EST_HAVE_STRINGS_H
#define P4EST_HAVE_STRINGS_H 1
#endif

/* Define to 1 if you have the <string.h> header file. */
#ifndef P4EST_HAVE_STRING_H
#define P4EST_HAVE_STRING_H 1
#endif

/* Define to 1 if you have the <sys/stat.h> header file. */
#ifndef P4EST_HAVE_SYS_STAT_H
#define P4EST_HAVE_SYS_STAT_H 1
#endif

/* Define to 1 if you have the <sys/types.h> header file. */
#ifndef P4EST_HAVE_SYS_TYPES_H
#define P4EST_HAVE_SYS_TYPES_H 1
#endif

/* Define to 1 if you have the <unistd.h> header file. */
#ifndef P4EST_HAVE_UNISTD_H
#define P4EST_HAVE_UNISTD_H 1
#endif

/* Define to 1 if zlib's adler32_combine links */
#ifndef P4EST_HAVE_ZLIB
#define P4EST_HAVE_ZLIB 1
#endif

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#ifndef P4EST_LT_OBJDIR
#define P4EST_LT_OBJDIR ".libs/"
#endif

/* DEPRECATED (use P4EST_ENABLE_MEMALIGN instead) */
#ifndef P4EST_MEMALIGN
#define P4EST_MEMALIGN 1
#endif

/* desired alignment of allocations in bytes */
#ifndef P4EST_MEMALIGN_BYTES
#define P4EST_MEMALIGN_BYTES (8)
#endif

/* DEPRECATED (use P4EST_WITH_METIS instead) */
/* #undef METIS */

/* DEPRECATED (use P4EST_ENABLE_MPI instead) */
#ifndef P4EST_MPI
#define P4EST_MPI 1
#endif

/* DEPRECATED (use P4EST_ENABLE_MPIIO instead) */
#ifndef P4EST_MPIIO
#define P4EST_MPIIO 1
#endif

/* Define to 1 if your C compiler doesn't accept -c and -o together. */
/* #undef NO_MINUS_C_MINUS_O */

/* Name of package */
#ifndef P4EST_PACKAGE
#define P4EST_PACKAGE "p4est"
#endif

/* Define to the address where bug reports for this package should be sent. */
#ifndef P4EST_PACKAGE_BUGREPORT
#define P4EST_PACKAGE_BUGREPORT "p4est@ins.uni-bonn.de"
#endif

/* Define to the full name of this package. */
#ifndef P4EST_PACKAGE_NAME
#define P4EST_PACKAGE_NAME "p4est"
#endif

/* Define to the full name and version of this package. */
#ifndef P4EST_PACKAGE_STRING
#define P4EST_PACKAGE_STRING "p4est 2.8.7"
#endif

/* Define to the one symbol short name of this package. */
#ifndef P4EST_PACKAGE_TARNAME
#define P4EST_PACKAGE_TARNAME "p4est"
#endif

/* Define to the home page for this package. */
#ifndef P4EST_PACKAGE_URL
#define P4EST_PACKAGE_URL ""
#endif

/* Define to the version of this package. */
#ifndef P4EST_PACKAGE_VERSION
#define P4EST_PACKAGE_VERSION "2.8.7"
#endif

/* DEPRECATED (use P4EST_WITH_PETSC instead) */
/* #undef PETSC */

/* Use builtin getopt */
/* #undef PROVIDE_GETOPT */

/* DEPRECATED (use P4EST_ENABLE_PTHREAD instead) */
/* #undef PTHREAD */

/* DEPRECATED (use P4EST_WITH_SC instead) */
/* #undef SC */

/* Define to 1 if you have the ANSI C header files. */
#ifndef P4EST_STDC_HEADERS
#define P4EST_STDC_HEADERS 1
#endif

/* Define to 1 if using autoconf build */
#ifndef P4EST_USING_AUTOCONF
#define P4EST_USING_AUTOCONF 1
#endif

/* Version number of package */
#ifndef P4EST_VERSION
#define P4EST_VERSION "2.8.7"
#endif

/* Package major version */
#ifndef P4EST_VERSION_MAJOR
#define P4EST_VERSION_MAJOR 2
#endif

/* Package minor version */
#ifndef P4EST_VERSION_MINOR
#define P4EST_VERSION_MINOR 8
#endif

/* Package point version */
#ifndef P4EST_VERSION_POINT
#define P4EST_VERSION_POINT 7
#endif

/* DEPRECATED (use P4EST_ENABLE_VTK_BINARY instead) */
#ifndef P4EST_VTK_BINARY
#define P4EST_VTK_BINARY 1
#endif

/* DEPRECATED (use P4EST_ENABLE_VTK_COMPRESSION instead) */
#ifndef P4EST_VTK_COMPRESSION
#define P4EST_VTK_COMPRESSION 1
#endif

/* DEPRECATED (use P4EST_ENABLE_VTK_DOUBLES instead) */
/* #undef VTK_DOUBLES */

/* enable metis-dependent code */
/* #undef WITH_METIS */

/* enable PETSc-dependent code */
/* #undef WITH_PETSC */

/* specifiy how to depend on package sc (optional). If the option value is
   literal no or the option is not present, use the source subdirectory. If
   the option value is the literal external, expect all necessary environment
   variables to be set to compile and link against sc and to run the examples.
   Otherwise, the option value must be an absolute path to the toplevel
   directory of a make installed sc. */
/* #undef WITH_SC */
 
/* once: _CONFIG_P_EST_CONFIG_H */
#endif
