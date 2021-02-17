:: build nectcdf-cxx4 first
::
set NC_BUILD_TYPE=Release
mkdir %SRC_DIR%\netcdf-cxx4\build
cd %SRC_DIR%\netcdf-cxx4\build
cmake -G "Visual Studio 15 2017 Win64" ^
    -DBUILD_SHARED_LIBS=OFF ^
    -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" ^
    -DCMAKE_LIBRARY_PATH="%LIBRARY_LIB%" ^
    -DCMAKE_PREFIX_PATH="%LIBRARY_PREFIX%" ^
    -DNETCDF_LIB_NAME="netcdf" ^
    -DHDF5_LIB_NAME="hdf5" ^
    %SRC_DIR%\netcdf-cxx4
if errorlevel 1 exit \b 1
cmake --build . --config %NC_BUILD_TYPE%
if errorlevel 1 exit \b 1
cmake --build . --config %NC_BUILD_TYPE% --target install
if errorlevel 1 exit \b 1

:: trilinos
mkdir %SRC_DIR%\trilinos_build
cd %SRC_DIR%\trilinos_build
cmake -G "Visual Studio 15 2017 Win64" ^
    -DBUILD_SHARED_LIBS=OFF ^
    -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" ^
    -DCMAKE_LIBRARY_PATH="%LIBRARY_LIB%" ^
    -DCMAKE_PREFIX_PATH="%LIBRARY_PREFIX%" ^
    -DNETCDF_LIB_NAME="netcdf" ^
    -DHDF5_LIB_NAME="hdf5" ^
    %SRC_DIR%\netcdf-cxx4
cmake -DCMAKE_BUILD_TYPE=RELEASE ^
      -DCMAKE_CXX_FLAGS=' -fPIC ' ^
      -DCMAKE_C_FLAGS=' -fPIC ' ^
      -DTrilinos_ENABLE_CXX11=ON ^
      -DTrilinos_ENABLE_Fortran=OFF ^
      -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" ^
      -DBUILD_SHARED_LIBS=ON ^
      -DTPL_ENABLE_LAPACK=ON ^
      -DTPL_ENABLE_MPI=OFF ^
      -DTPL_ENABLE_UMFPACK=ON ^
      -DTrilinos_ENABLE_Amesos=ON ^
      -DTrilinos_ENABLE_Amesos2=ON ^
      -DTrilinos_ENABLE_AztecOO=ON ^
      -DTrilinos_ENABLE_Belos=ON ^
      -DTrilinos_ENABLE_Ifpack=ON ^
      -DTrilinos_ENABLE_Ifpack2=ON ^
      -DTrilinos_ENABLE_Kokkos=ON ^
      -DTrilinos_ENABLE_Komplex=ON ^
      -DTrilinos_ENABLE_ML=ON ^
      -DTrilinos_ENABLE_MueLu=ON ^
      -DTrilinos_ENABLE_Teuchos=ON ^
      -DTrilinos_ENABLE_Tpetra=ON ^
      -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=ON ^
      -DKOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION=ON ^
      -DTpetra_INST_COMPLEX_DOUBLE=ON ^
      -DTeuchos_ENABLE_COMPLEX=ON ^
      -DTpetraKernels_ENABLE_Experimental=ON ^
      -DTrilinos_ENABLE_OpenMP=ON ^
      -DTrilinos_ENABLE_EXPLICIT_INSTANTIATION=ON ^
      -DKOKKOS_ENABLE_COMPILER_WARNINGS=ON ^
      -DAmesos2_ENABLE_Basker=ON ^
      -DTpetra_INST_SERIAL:BOOL=ON ^
      -DTpetra_INST_INT_INT=ON ^
%SRC_DIR%/trilinos_source
make -j"${CPU_COUNT}" install

make -j%CPU_COUNT% install

:: now build silo
::
set SILO_BUILD_TYPE=Release
cd %SRC_DIR%\silo\SiloWindows\MSVC2012
msbuild Silo.sln ^
    -toolsVersion:15.0 ^
    -property:Configuration=%SILO_BUILD_TYPE% ^
    -property:PlatformToolset=v141;WindowsTargetPlatformVersion=%WindowsSDKVer% ^
    -property:HDF5_INC_DIR_X64=%LIBRARY_INC% ^
    -property:HDF5_LIB_DIR_X64=%LIBRARY_LIB%
if errorlevel 1 exit \b 1
echo Install configuration: "%SILO_BUILD_TYPE%"
copy /y %SRC_DIR%\silo\SiloWindows\include\silo*.h "%LIBRARY_INC%"
if errorlevel 1 exit \b 1
copy /y %SRC_DIR%\silo\src\silo\silo_exports.h "%LIBRARY_INC%"
if errorlevel 1 exit \b 1
copy /y %SRC_DIR%\silo\SiloWindows\MSVC2012\x64\%SILO_BUILD_TYPE%\silohdf5.lib "%LIBRARY_LIB%"
if errorlevel 1 exit \b 1
copy /y %SRC_DIR%\silo\SiloWindows\MSVC2012\x64\%SILO_BUILD_TYPE%\silohdf5.dll "%LIBRARY_BIN%"
if errorlevel 1 exit \b 1

:: workaround for bug in win-64/boost-1.73.0-py38_11
if not exist "%LIBRARY_INC%\boost\python.hpp" (
    if exist "%LIBRARY_INC%\boost\python\python.hpp" (
        copy "%LIBRARY_INC%\boost\python\python.hpp" "%LIBRARY_INC%\boost"
    )
)
:: now build escript
::
cd %SRC_DIR%\escript
call scons -j%CPU_COUNT% ^
    options_file="scons\templates\windows_msvc141_options.py" ^
    compressed_files=0 ^
    silo=1 ^
    prefix="%PREFIX%" ^
    build_dir="%BUILD_PREFIX%\escript_build" ^
    build_full
if errorlevel 1 exit \b 1

:: type config.log
copy /y %SRC_DIR%\escript\LICENSE %SRC_DIR%\LICENSE
copy /y %PREFIX%\esys %SP_DIR%\esys
copy /y %BUILD_PREFIX%\escript_build\scripts\release_sanity.py %TEMP%\release_sanity.py