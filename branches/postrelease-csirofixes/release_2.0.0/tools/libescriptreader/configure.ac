#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([libescriptreader], [1.0])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/esdcreate.cpp])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AM_PROG_LIBTOOL

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stdint.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_C_VOLATILE

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([pow sqrt])

CXXFLAGS="-Wall -O2"

AC_ARG_ENABLE([esdcreate],
              AC_HELP_STRING([--disable-esdcreate],
                             [do not build the esdcreate tool]),
              [
               if test "x$enableval" = "xno"; then
                   wants_esdcreate=false
               else
                   wants_esdcreate=true
               fi
              ],
              [wants_esdcreate=true]
             )

AM_CONDITIONAL([BUILD_ESDCREATE], [test "$wants_esdcreate" = true])

AC_ARG_ENABLE([escript2silo],
              AC_HELP_STRING([--disable-escript2silo],
                             [do not build the escript to Silo converter]),
              [
               if test "x$enableval" = "xno"; then
                   wants_escript2silo=false
               else
                   wants_escript2silo=true
               fi
              ],
              [wants_escript2silo=true]
             )

AC_ARG_WITH([escript-prefix],
            [AS_HELP_STRING([--with-escript-prefix=PATH],
                            [Prefix of escript installation containing headers])],
            [ESCRIPTDIR="$withval"])

AC_MSG_CHECKING([for escript include files])
if test "x$ESCRIPTDIR" = "x"; then
    ESCRIPTDIR=`../..`
fi

if ! test -f "$ESCRIPTDIR/include/finley/ReferenceElements.h"; then
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([
    The escript headers could not be found. Please provide the
    location of your escript installation using --with-escript-prefix
                  ])
else
    ESCRIPT_CFLAGS=" -I$ESCRIPTDIR/include"
    AC_MSG_RESULT([yes])
fi

AC_SUBST([ESCRIPT_CFLAGS])


# Check for NetCDF library and headers

AC_ARG_WITH([netcdf-prefix],
            [AS_HELP_STRING([--with-netcdf-prefix=PATH],
                            [Prefix of NetCDF headers and libraries])],
            [NETCDFPREFIX="$withval"])

AC_MSG_CHECKING([for NetCDF files])
if ! test "x$NETCDFPREFIX" = "x"; then
    CPPFLAGS_SAVED="$CPPFLAGS"
    NETCDF_CFLAGS="-I$NETCDFPREFIX/include"
    NETCDF_LIBS="-L$NETCDFPREFIX/lib -lnetcdf_c++ -lnetcdf"
    CPPFLAGS="$CPPFLAGS $NETCDF_CFLAGS"
    export CPPFLAGS
else
    NETCDF_CFLAGS=""
    NETCDF_LIBS="-lnetcdf_c++ -lnetcdf"
fi
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE(
               [@%:@include <netcdf.hh>],
               [NcError ncerr(NcError::verbose_fatal);],
               [AC_MSG_RESULT([yes])],
               [AC_MSG_ERROR([
    The NetCDF headers could not be found. Please provide the
    location of your NetCDF installation using --with-netcdf-prefix
                  ])
               ]
)
AC_LANG_RESTORE
if ! test "x$NETCDFPREFIX" = "x"; then
    CPPFLAGS="$CPPFLAGS_SAVED"
fi

AC_SUBST([NETCDF_CFLAGS])
AC_SUBST([NETCDF_LIBS])


# Check for Silo library and headers

has_silo_lib=false
AC_ARG_WITH([silo-prefix],
            [AS_HELP_STRING([--with-silo-prefix=PATH],
                            [Prefix of Silo headers and libraries])],
            [SILOPREFIX="$withval"])

AC_MSG_CHECKING([for Silo files])
if test "x$SILOPREFIX" = "x"; then
    if test -f "/usr/include/silo.h"; then
        SILOPREFIX="/usr"
    elif test -f "/usr/local/include/silo.h"; then
        SILOPREFIX="/usr/local"
    fi
fi

if test -f "$SILOPREFIX/include/silo.h"; then
    SILO_CFLAGS=" -I$SILOPREFIX/include"

    if test -f "$SILOPREFIX/lib/libsiloh5.a"; then
        SILO_LIBS="-L$SILOPREFIX/lib -lsiloh5"
        has_silo_lib=true
    elif test -f "$SILOPREFIX/lib/libsilo.a"; then
        SILO_LIBS="-L$SILOPREFIX/lib -lsilo"
        has_silo_lib=true
    fi
fi
AC_SUBST([SILO_CFLAGS])
AC_SUBST([SILO_LIBS])

if test "$has_silo_lib" = "true"; then
    AC_MSG_RESULT([yes])
    AC_DEFINE([HAVE_SILO], [1])
    AC_CHECK_LIB([mpi], [MPI_Wait])
    AC_CHECK_LIB([z], [deflate])
else
    AC_MSG_RESULT([no])
    AC_MSG_WARN([
    The Silo library could not be found. Your escriptreader library will not
    be able to save data in the Silo format and the escript2silo tool will
    not be built. If you need this functionality please download the Silo
    library from http://silo.llnl.gov/ and specify the install location using
    --with-silo-prefix.
               ])
    wants_escript2silo=false
fi

AM_CONDITIONAL([BUILD_ESCRIPT2SILO], [test "$wants_escript2silo" = true])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/escriptreader/Makefile
                 ])
AC_OUTPUT
