
set CMAKE_C_COMPILER=clang
set CMAKE_CXX_COMPILER=clang++

rem :: Set the right msvc version according to Python version
rem if "%PY_VER%"=="2.7" (
rem     set MSVC_VER=9.0
rem     set LIB_VER=90
rem ) else if "%PY_VER%"=="3.4" (
rem     set MSVC_VER=10.0
rem     set LIB_VER=100
rem ) else (
rem     set MSVC_VER=14.1
rem     set LIB_VER=141
rem )

set INCLUDE_PATH="%PREFIX%/include"
set LIBRARY_PATH="%PREFIX%/lib"
set LDFLAGS=-L%LIBRARY_PATH%
set CMAKE_INSTALL_PREFIX=%PREFIX%

cd %SRC_DIR%/escript
scons -j%CPU_COUNT% ^
    options_file="%SRC_DIR%/escript/scons/templates/stretch_options.py" ^
    prefix="%PREFIX%" ^
    cxx_extra="-fPIC" ^
    ld_extra="-L%PREFIX%/lib -lblas -llapack -lumfpack -lnetcdf_c++4" ^
    boost_prefix="%PREFIX%" ^
    boost_libs='boost_python27' ^
    pythoncmd=`which python` ^
    pythonlibpath="%PREFIX%/lib" ^
    pythonincpath="%PREFIX%/include/python%PY_VER%" ^
    pythonlibname="python2.7" ^
    trilinos=0 ^
    trilinos_prefix="%PREFIX%" ^
    umfpack=1 ^
    umfpack_prefix=["%PREFIX%/include","%PREFIX%/lib"] ^
    lapack=1 ^
    lapack_prefix=["%PREFIX%/include/atlas","%PREFIX%/lib"] ^
    lapack_libs=['lapack'] ^
    werror=0

