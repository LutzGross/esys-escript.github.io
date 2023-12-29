rem Call this script and pass in the base directory for the dashboard
rem save state before changing anything
setlocal
set foo=%1%

rem cygwin CVS needs to know to use ssh
set CVS_RSH=C:/cygwin/bin/ssh.exe
rem Driver script for dashboards on Corrin
rem Set location of CTEST_EXE, and CVS_EXE
set CVS_EXE=C:\cygwin\bin\cvs.exe
set CTEST_EXE="c:\hoffman\My Builds\CMake-build26-rel\bin\ctest.exe"

rem Set the base directory which is one above where Trilinos will be 
rem checked out.

set BASEDIR=%1%

rem setup the environment for command line cl to run
call "C:\Program Files\Microsoft Visual Studio 9.0\Common7\Tools\vsvars32.bat"

rem change into the basedir
cd %BASEDIR%
rem checkout the basics from Trilinos needed to run the dashboard including
rem this script.
%CVS_EXE% -q -d :ext:software.sandia.gov:/space/CVS co Trilinos/cmake Trilinos/CTestConfig.cmake

rem Now run ctest on each of the ctest build scripts for this machine

%CTEST_EXE% -S %BASEDIR%\Trilinos\cmake\ctest\drivers\corrin\ctest_windows_nightly_serial_release.cmake -VV > %BASEDIR%\ctest_msvc_nightly_serial_optimized_corrin.out
endlocal
