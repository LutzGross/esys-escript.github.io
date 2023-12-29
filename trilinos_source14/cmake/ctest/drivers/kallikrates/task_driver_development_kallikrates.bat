rem Driver script for dashboards on kallikrates

rem Call this script and pass in the base directory for the dashboard
rem save state before changing anything
setlocal

rem Setting the path to have git on it.
set PATH=%PATH%;C:\Program Files (x86)\Git\cmd;C:\Qt\4.6.3\bin

rem GIT needs to know how to use ssh
set GIT_SSH=C:\Users\bmpersc\Documents\plink_wrap.bat

rem Set location of CTEST_EXE, and GIT_EXE
set GIT_EXE=C:\Program Files (x86)\Git\cmd\git
set CTEST_EXE=C:\Program Files (x86)\CMake 2.8\bin\ctest.exe
set TRILINOS_REPOSITORY_LOCATION=software.sandia.gov:/space/git/nightly/Trilinos

rem Set the base directory which is one above where Trilinos will be 
rem checked out.

set BASEDIR=%~1

rem setup the environment for command line cl to run
call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\Common7\Tools\vsvars32.bat"

rem change into the basedir
cd "%BASEDIR%"

rem checkout the basics from Trilinos needed to run the dashboard including
rem this script.
if exist Trilinos goto update else goto checkout

:update
  echo "Doing update of an existing directory"
  cd Trilinos
  call "%GIT_EXE%" pull
  cd ..
  goto endif

:checkout
  echo "Cloning the repository because none exists yet."
  call "%GIT_EXE%" clone %TRILINOS_REPOSITORY_LOCATION%

:endif

rem Now run ctest on each of the ctest build scripts for this machine

call "%CTEST_EXE%" -S "%BASEDIR%\Trilinos\cmake\ctest\drivers\kallikrates\ctest_windows_nightly_serial_development.cmake" -VV >"%BASEDIR%\ctest_msvc_nightly_serial_development_kallikrates.out" 2>&1

rem Have to set the path so that the tests can find the dlls during runtime. We need a better solution than this long term. Something like -rpath for gnu.
set PATH=%PATH%;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\epetra\src;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\teuchos\src;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\anasazi\src;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\epetra\test\src;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\teuchos\test\FancyOutputting;%BASEDIR%\SERIAL_OPT_DEV_SHARED\BUILD\packages\teuchos\test\ParameterList

call "%CTEST_EXE%" -S "%BASEDIR%\Trilinos\cmake\ctest\drivers\kallikrates\ctest_windows_nightly_serial_development_shared.cmake" -VV >"%BASEDIR%\ctest_msvc_nightly_serial_development_shared_kallikrates.out" 2>&1


endlocal
