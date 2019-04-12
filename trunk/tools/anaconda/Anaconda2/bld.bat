
:: Set the right msvc version according to Python version
if "%PY_VER%"=="2.7" (
    set MSVC_VER=9.0
    set LIB_VER=90
) else if "%PY_VER%"=="3.4" (
    set MSVC_VER=10.0
    set LIB_VER=100
) else (
    set MSVC_VER=14.1
    set LIB_VER=141
)

set PREFIX=pwd
set INCLUDE_PATH="%PREFIX%/include"
set LIBRARY_PATH="%PREFIX%/lib"
set LDFLAGS=-L%LIBRARY_PATH%
set CMAKE_INSTALL_PREFIX=%PREFIX%

mkdir %SRC_DIR%/lapack_build
cd %SRC_DIR%/lapack_build
cmake -DCMAKE_INSTALL_PREFIX=%PREFIX% ^
    -DCMAKE_MAKE_PROGRAM=%PREFIX%/bin/make ^
    -DCMAKE_INSTALL_LIBDIR=%LIBRARY_PATH% ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DBUILD_SHARED_LIBS=OFF ^
    -DBUILD_TESTING=OFF ^
    -DCMAKE_CC_FLAGS='-fPIC' ^
    -DCMAKE_CXX_FLAGS='-fPIC' ^
    -DLAPACKE=ON ^
    -DCBLAS=ON ^
    %SRC_DIR%/lapack_source
make -j"%CPU_COUNT%" install

rem mkdir %SRC_DIR%/netcdf_build
rem cd %SRC_DIR%/netcdf_build
rem cmake ^
rem      -DCMAKE_INSTALL_PREFIX=%PREFIX% ^
rem      -DCMAKE_SHARED_LINKER_FLAGS=%LDFLAGS% ^
rem      -DCMAKE_CC_FLAGS='-fPIC' ^
rem      -DCMAKE_CXX_FLAGS='-fPIC' ^
rem      %SRC_DIR%/netcdf_source 
rem make -j"%CPU_COUNT%" install

rem cd %SRC_DIR%/netcdf_cxx_source
rem ./configure --prefix=%PREFIX%
rem make -j"%CPU_COUNT%" install

rem mkdir ${SRC_DIR}/netcdf_cxx_build
rem cd ${SRC_DIR}/netcdf_cxx_build
rem cmake ^
rem      -DCMAKE_INSTALL_PREFIX=%PREFIX% ^
rem      -DCMAKE_SHARED_LINKER_FLAGS=%LDFLAGS% ^
rem      -DCMAKE_CC_FLAGS='-fPIC' ^
rem      -DCMAKE_CXX_FLAGS='-fPIC' ^
rem      ${SRC_DIR}/netcdf_source 
rem make -j%CPU_COUNT% install

rem cd ${SRC_DIR}/netcdf_cxx_source
rem ./configure --prefix=%PREFIX%
rem make -j%CPU_COUNT% install

cd ${SRC_DIR}/boost_source
./bootstrap.sh ^
    --prefix="%PREFIX%" ^
    --with-python="${PYTHON}" ^
    --with-python-root="%PREFIX% : %PREFIX%/include/python${PY_VER}m %PREFIX%/include/python${PY_VER}" ^
    --with-icu="%PREFIX%" 
./b2 -q -d+2 -a ^
    variant=release ^
    address-model="${ARCH}" ^
    architecture=x86 ^
    debug-symbols=off ^
    threading=multi ^
    runtime-link=static ^
    link=static ^
    toolset=gcc ^
    cxxflags='-fPIC -w' ^
    python="${PY_VER}" ^
    include="%PREFIX%/include" ^
    linkflags="-L%LIBRARY_PATH%" ^
    --with-python ^
    --with-iostreams ^
    --with-random ^
    --layout=system ^
    -j%CPU_COUNT% ^
    install

cd ${SRC_DIR}/trilinos_source
git reset --hard c87c9f7224b36c845aa7ef19a636f9c99c9c02d4

mkdir ${SRC_DIR}/trilinos_build
cd ${SRC_DIR}/trilinos_build
cmake -DCMAKE_INSTALL_PREFIX=%PREFIX% ^
    -DCMAKE_MAKE_PROGRAM=%PREFIX%/bin/make ^
    -DCMAKE_VERBOSE_MAKEFILE=FALSE ^
    -DTrilinos_ENABLE_CXX11=ON ^
    -DTrilinos_ENABLE_Fortran=OFF ^
    -DBUILD_SHARED_LIBS=OFF ^
    -DOpenMP_C_FLAGS='-w -fopenmp' ^
    -DOpenMP_CXX_FLAGS='-w -fopenmp' ^
    -DCMAKE_CC_FLAGS='-fPIC' ^
    -DCMAKE_CXX_FLAGS='-fPIC' ^
    -DTPL_ENABLE_BLAS=ON ^
    -DBLAS_LIBRARY_NAMES='blas' ^
    -DTPL_BLAS_INCLUDE_DIRS=%PREFIX%/include ^
    -DTPL_BLAS_LIBRARY_DIRS=%LIBRARY_PATH% ^
    -DTPL_BLAS_LIBRARY_NAMES='blas' ^
    -DTPL_ENABLE_LAPACK=ON ^
    -DTPL_LAPACK_INCLUDE_DIRS=%PREFIX%/include ^
    -DTPL_LAPACK_LIBRARY_DIRS=%LIBRARY_PATH% ^
    -DTPL_LAPACK_LIBRARY_NAMES='lapack' ^
    -DTPL_ENABLE_Boost=ON ^
    -DTPL_Boost_INCLUDE_DIRS=%PREFIX%/include ^
    -DTPL_Boost_LIBRARY_DIRS=%LIBRARY_PATH% ^
    -DTPL_ENABLE_Cholmod=OFF ^
    -DTPL_Cholmod_INCLUDE_DIRS=%PREFIX%/include ^
    -DTPL_Cholmod_LIBRARY_DIRS=%LIBRARY_PATH% ^
    -DTPL_Cholmod_LIBRARIES='libcholmod.so;libamd.so;libcolamd.so' ^
    -DTPL_ENABLE_METIS=ON ^
    -DTPL_METIS_INCLUDE_DIRS=%PREFIX%/include ^
    -DTPL_METIS_LIBRARY_DIRS=%LIBRARY_PATH% ^
    -DMETIS_LIBRARY_NAMES='libmetis.a' ^
    -DTPL_ENABLE_UMFPACK=OFF ^
    -DTPL_UMFPACK_INCLUDE_DIRS=%PREFIX%/include ^
    -DTPL_UMFPACK_LIBRARY_DIRS=%LIBRARY_PATH% ^
    -DTrilinos_ENABLE_Amesos2=ON ^
    -DTrilinos_ENABLE_Belos=ON ^
    -DTrilinos_ENABLE_Ifpack2=ON ^
    -DTrilinos_ENABLE_Kokkos=ON ^
    -DTrilinos_ENABLE_MueLu=ON ^
    -DTrilinos_ENABLE_Tpetra=ON ^
    -DTrilinos_ENABLE_Teuchos=ON ^
    -DTrilinos_ENABLE_COMPLEX=ON ^
    -DTrilinos_ENABLE_OpenMP=ON ^
    -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF ^
    -DTrilinos_ENABLE_EXPLICIT_INSTANTIATION=ON ^
    ${SRC_DIR}/trilinos_source
make -j%CPU_COUNT% install

cd ${SRC_DIR}/escript
scons -j%CPU_COUNT% ^
    options_file="${SRC_DIR}/escript/scons/templates/sid_options.py" ^
    prefix="%PREFIX%" ^
    cxx_extra="-fPIC" ^
    boost_prefix="%PREFIX%" ^
    boost_libs='boost_python27' ^
    trilinos=1 ^
    trilinos_prefix="%PREFIX%" ^
    pythonlibpath="%LIBRARY_PATH%" ^
    pythonincpath="%PREFIX%/include/python${PY_VER}" ^
    pythonlibname="python2.7" ^
    umfpack=0 ^
    umfpack_prefix=["%PREFIX%/include","%LIBRARY_PATH%"] ^
    werror=0

