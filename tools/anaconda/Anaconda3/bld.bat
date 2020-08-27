::@echo off

:: get boost library name
for /r "%LIBRARY_LIB%" %%f in (boost_python%PY_VER:.=%*.*) do (set "ES_BOOST_LIBS=%%~nf")
set "VCPKG_PREFIX=C:\Users\%USERNAME%\vcpkg\packages"
set "ES_LD_EXTRA=/LIBPATH:%VCPKG_PREFIX%\hdf5_x64-windows\lib"
set "ES_LD_EXTRA=%ES_LD_EXTRA% /LIBPATH:%VCPKG_PREFIX%\curl_x64-windows\lib"
set "ES_LD_EXTRA=%ES_LD_EXTRA% /LIBPATH:%VCPKG_PREFIX%\zlib_x64-windows\lib"
set "ES_LD_EXTRA=%ES_LD_EXTRA% /LIBPATH:%VCPKG_PREFIX%\szip_x64-windows\lib"
set "ES_LD_EXTRA=%ES_LD_EXTRA% /LIBPATH:%LIBRARY_PREFIX%\mingw-w64\lib"
set "ES_LD_EXTRA=%ES_LD_EXTRA% /LIBPATH:%LIBRARY_PREFIX%\mingw-w64\bin"
set "ES_LD_EXTRA=%ES_LD_EXTRA% libmumps_common.a libdmumps.dll.a libzmumps.dll.a"
set "PATH=%PATH%;%LIBRARY_PREFIX%\mingw-w64\bin;%VCPKG_PREFIX%\cppunit_x64-windows\bin"
::set "PATH=%PATH%;%LIBRARY_PREFIX%\mingw-w64\bin;%VCPKG_PREFIX%\hdf5_x64-windows\bin;%VCPKG_PREFIX%\curl_x64-windows\bin;%VCPKG_PREFIX%\zlib_x64-windows\bin;%VCPKG_PREFIX%\szip_x64-windows\bin;%VCPKG_PREFIX%\cppunit_x64-windows\bin"

:: netcdf libhdf5 libcurl zlib szip
cd %SRC_DIR%/escript
::call scons -j"%CPU_COUNT%" ^
::    verbose=1 ^
::    prefix="%PREFIX%" ^
::    build_dir="%BUILD_PREFIX%\escript_build" ^
::    cc_flags="/EHsc /MD /DBOOST_ALL_NO_LIB /wd4068" ^
::    ld_extra="%ES_LD_EXTRA%" ^
::    pythoncmd="%PYTHON%" ^
::    boost_libs="%ES_BOOST_LIBS%" ^
::    boost_prefix="%LIBRARY_PREFIX%" ^
::    netcdf=4 ^
::    netcdf_libs="netcdf-cxx4" ^
::    netcdf_prefix="%VCPKG_PREFIX%\netcdf-cxx4_x64-windows" ^
::    cppunit_libs="cppunit_dll" ^
::    cppunit_prefix="%VCPKG_PREFIX%\cppunit_x64-windows" ^
::    mumps=1 ^
::    mumps_libs="" ^
::    mumps_prefix="%LIBRARY_PREFIX%\mingw-w64" ^
::    openmp=1 ^
::    omp_flags="/openmp" ^
::    build_full
call scons -j%CPU_COUNT% ^
    options_file="scons\templates\windows_msvc141_options.py" ^
    prefix="%PREFIX%" ^
    build_dir="%BUILD_PREFIX%\escript_build" ^
    build_full
if errorlevel 1 exit 1

type config.log
copy /y %SRC_DIR%\escript\LICENSE %SRC_DIR%\LICENSE
copy /y %PREFIX%\esys %SP_DIR%\esys
copy /y %BUILD_PREFIX%\escript_build\scripts\release_sanity.py %TEMP%\release_sanity.py
