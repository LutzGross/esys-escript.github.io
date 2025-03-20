# *esys-escript* installation guide (under review)
This is for Version 6.

(sorry this document is still a mess)

*esys-escript* is primarily developed on Linux systems and macOS systems. 
The current version has not been ported to Windows. 
For support and questions please visit the [Issues site](https://github.com/LutzGross/esys-escript.github.io/issues) for help.

Dependencies are

- c++ compiler
- python3 
- scons
- cmake
- python3-numpy
- boost-python
- boost-random
- boost-iostreams
- [mpi4py](https://mpi4py.readthedocs.io) or openMPI (with MPI)
- parmetis (with MPI)
- hdf5-serial (optional)
- silo (optional,  recommended when using VisIt)
- boost-numpy (optional, recommended when using Jupyter)
- python3-sympy (optional, **currently deactivated**)
- cppunit (optional, for testing)

The source code is shipped with [Trilinos](https://trilinos.github.io/) and [p4est](https://www.p4est.org/).

Recommended tools are

- [gmsh](https://gmsh.info/) for mesh generation,
- [VisIt](https://visit-dav.github.io) for visualization (via Silo or VTK),
- [paraview](https://www.paraview.org/) for visualization (via VTK),
- [mayavi](https://docs.enthought.com/mayavi/mayavi/) for visualization (via VTK),
- [python3-matplotlib](https://matplotlib.org/) for plotting.

They are not required to be installed at compile time. Please visit the particular package for installation instructions.

For parallelization *esys-escript* supports

- [OpenMP](https://www.openmp.org/) for shared memory parallelization (multi-core): needs to be activated at compilation.
- [MPI](https://www.openmp.org/) for distributed memory (compute clusters): requires compilation with mpi4py
- hybrid mode with both OpenMP and MPI parallelization.


## Installation from Source

### General setup 

The first step is to download the newest version of *esys-escript* source from 
[gitHub](https://github.com/LutzGross/esys-escript.github.io) using  

    mkdir esys6
    git clone --single-branch --branch master https://github.com/LutzGross/esys-escript.github.io.git esys6
    cd esys6

Alternatively, you can download a tagged version (zipped) from [tagged versions](https://github.com/LutzGross/esys-escript.github.io/tags)

After you have installed the required dependencies (see below) start *scons* to run the installation from the source directory:     

    scons -j4 options_file=scons/templates/<OS>_options.py

where `<OS>` is your operating system/kernel name. Typically, the templates are assuming installation with MPI. 
By default, libraries and *Python* modules are installed into `./lib/esys` and `./esys`, respectively.
*Trilinos* libraries are installed into `esys.trilions/lib`. So you need to add 
    
    export PYTHONPATH=$ESYSESCRIPT
    export LD_LIBRARY_PATH=$ESYSESCRIPT/lib/esys:$ESYSESCRIPT/esys.trilinos/lib

to your environment where `ESYSESCRIPT` is the root to the source (installation) directory.

If you wish to test your build, you can use the following:

    scons -j4 py_tests options_file=scons/templates/<OS>_options.py

If the options file of your operating system/kernel is not available yet, or you want to modify the options
(for instance switching off MPI) you can create your own option file  `scons/<HOST>_options.py` typically by copying and modifying one of the available
template options files in [scons/templates](scons/templates). `<HOST>` refers to the hostname which can be obtained by

    HOST=`uname -n`
    echo $HOST

Installation can then be run by

    scons -j4    



### Debian
These instructions have been tested on Debian 6.

As a preliminary step, you need to install the dependencies that *esys-escript* requires from the distribution repository:

    sudo apt-get install python3-dev python3-numpy 
    sudo apt-get install g++ scons cmake 
    sudo apt-get install libboost-numpy-dev libboost-python-dev libboost-random-dev libboost-iostreams-dev
    sudo apt-get install libhdf5-serial-dev libsilo-dev libsuitesparse-dev
    #sudo apt-get install python3-sympy

If you want to use MPI you additionally need to install 
    
    sudo apt-get install python3-mpi4py libparmetis-dev
    

### Ubuntu
These instructions have been tested on Ubuntu 24.

You should instead install the following packages:

    sudo apt-get install python3-dev python3-numpy
    sudo apt-get install g++ scons cmake 
    sudo apt-get install libboost-numpy-dev libboost-python-dev libboost-random-dev libboost-iostreams-dev
    sudo apt-get install libhdf5-serial-dev libsilo-dev libsuitesparse-dev
    #sudo apt-get install python3-sympy

If you want to use MPI you additionally need to install 
    
    sudo apt-get install python3-mpi4py libparmetis-dev


--------------
# The following is WORK IN PROGRESS
Chapter~\ref{chap:source} covers installing from source.
Appendix~\ref{app:cxxfeatures} lists some c++ features which your compiler must support in order to compile *esys-escript*.
This version of *esys-escript* has the option of using \texttt{Trilinos} in addition to our regular solvers.
Appendix~\ref{app:trilinos} covers features of \texttt{Trilinos} which *esys-escript* needs.

## Installing from Flatpak}\label{chap:flatpak}

To install *esys-escript* on any linux distribution using flatpak\footnote{For most linux distributions this can be installed from the repository. Otherwise, flatpak is available at \url{https://flathub.org/home}}, type
\begin{shellCode}
flatpak install flathub au.edu.uq.esys.*esys-escript*
\end{shellCode}

This wil download and install *esys-escript* on your machine. The *esys-escript* build installed utilises both the Trilinos and PASO solver libraries, with openMP but without OPENMPI. 

After flatpak has finished downloading and installing *esys-escript*, you can launch an *esys-escript* window from the menu or run *esys-escript* in a terminal with the command:
\begin{shellCode}
flatpak run au.edu.uq.esys.*esys-escript* [other arguments to pass to *esys-escript*]
\end{shellCode}

Finally, to uninstall *esys-escript* from your machine, type
\begin{shellCode}
flatpak uninstall *esys-escript*
\end{shellCode}

#--------------
## Installing inside Anaconda}\label{chap:conda}

There are precompiled binaries of *esys-escript* available for Anaconda python (for Linux), available in the conda-forge channel on Anaconda cloud. They can be installed using the command
\begin{shellCode}
conda install *esys-escript*-c conda-forge
\end{shellCode}
#--------------

## Installing from Docker}\label{chap:docker}

To install an *esys-escript* Docker container on your machine, first install Docker and then type:
\begin{shellCode}
docker pull esys*esys-escript*/esys-*esys-escript*
\end{shellCode}

\vspace{1cm}\
Once installed, you can launch an *esys-escript* session using the command:
\begin{shellCode}
docker run -ti esys*esys-escript*/*esys-escript*run-*esys-escript*
\end{shellCode}

If you would also like to mount a folder, you can use the command:
\begin{shellCode}
docker run -ti -v $(pwd):/app/ esys*esys-escript*/*esys-escript*run-*esys-escript*
\end{shellCode}

Additionally, if you wish to run an *esys-escript* named test_program.py, you can use the command:
\begin{shellCode}
docker run -ti -v $(pwd):/app/ esys*esys-escript*/*esys-escript*\
					run-*esys-escript* [path to test_program.py]
\end{shellCode}

## Installing from Source

### Parallelization
It is likely that the computer you run *esys-escript* on, will have more than one processor core.
*esys-escript* can make use of multiple cores in order to solve problems more quickly if it is told to do so,
but this functionality must be enabled at compile time.
Section~\ref{sec:needpar} gives some rough guidelines to help you determine what you need.

There are two technologies which *esys-escript* can employ here.

 - OpenMP -- the more efficient of the two [thread level parallelism].
 - MPI -- Uses multiple processes (less efficient), needs less help from
   the compiler (not supported on Windows).
\end{itemize}

*esys-escript* is primarily tested on recent versions of the GNU and Intel suites
(``g++'' / ``icpc''), and Microsoft Visual C++ (MSVC).  However, it also passes
our tests when compiled using ``clang++''.  *esys-escript* now requires compiler
support for some features of the C++11 standard.  See
Appendix~\ref{app:cxxfeatures} for a list.


Our current test compilers include:

 - g++ 10.2
 - clang++ 11.0
 - intel icpc v17
 - MSVC 2017 or 2019

Note that:
 - OpenMP will not function correctly for g++ $\leq$ 4.2.1 (and is not currently supported by clang).
 - icpc v11 has a subtle bug involving OpenMP and C++ exception handling, so this combination should not be used.

\subsection{What parallel technology do I need?}\label{sec:needpar} If you are
using any version of Linux released in the past few years, then your system
compiler will support \openmp with no extra work (and give better performance);
so you should use it. MSVC 2017 and 2019 also have \openmp support on Windows
(\openmp 2.0). You will not need MPI unless your computer is some form of
cluster. MPI is not recommended on Windows as it will interfer with Jupyter.

If you are using BSD or MacOSX and you are just experimenting with *esys-escript*, then performance is
probably not a major issue for you at the moment, so you don't need to use either \openmp or MPI.
This also applies if you write and polish your scripts on your computer and then send them to a cluster to execute.
If in the future you find *esys-escript* useful and your scripts take significant time to run, then you may want to recompile
*esys-escript* with more options.



Note that even if your version of *esys-escript* has support for \openmp or MPI, you will still need to tell the system to
use it when you run your scripts.
If you are using the \texttt{run-*esys-escript*} launcher, then this is controlled by
the \texttt{-t}, \texttt{-p}, and \texttt{-n} flags.
If not, then consult the documentation for your MPI libraries (or the compiler documentation in the case of OpenMP
\footnote{It may be enough to set the \texttt{OMP\_NUM\_THREADS} environment variable.}).

If you are using MacOSX, then see the next section, if not, then skip to Section~\ref{sec:build}.

\section{MacOS}
This release of *esys-escript* has only been tested on OSX 10.13.
For this section we assume you are using either \texttt{homebrew} or \texttt{MacPorts} as a package
manager\footnote{Note that package managers will make changes to your computer based on programs configured by other people from
various places around the internet. It is important to satisfy yourself as to the security of those systems.}.
You can of course install prerequisite software in other ways.
For example, we have had \emph{some} success changing the default
compilers used by those systems. However, this is more complicated, and we do not provide a guide here.

\noindent Both of those systems require the XCode command line tools to be installed\footnote{As of OSX10.9, the
command \texttt{xcode-select --install} will allow you to download and install the commandline tools.}.

\section{Building}\label{sec:build}

*esys-escript* is built using \textit{SCons}. To simplify the installation process, we have prepared \textit{SCons} \textit{_options.py} files for a number of common systems\footnote{These are correct a time of writing but later versions of those systems may require tweaks.
Also, these systems represent a cross section of possible platforms rather than meaning those systems get particular support.}.
The options files are located in the \textit{scons/templates} directory. We suggest that the file most relevant to your OS
be copied from the templates directory to the scons directory and renamed to the form XXXX_options.py where XXXX
should be replaced with your computer's (host-)name.
If your particular system is not in the list below, or if you want a more customised
build,
see Section~\ref{sec:othersrc} for instructions.
\begin{itemize}
 - Debian - \ref{sec:debsrc}
 - Ubuntu - \ref{sec:ubsrc}
 - Mint - \ref{sec:mintsrc}
 - Arch Linux - \ref{sec:archsrc}
 - OpenSuse - \ref{sec:susesrc}
 - Centos - \ref{sec:centossrc}
 - Fedora - \ref{sec:fedorasrc}
 - MacOS (macports) - \ref{sec:macportsrc}
 - MacOS (homebrew) - \ref{sec:homebrewsrc}
 - FreeBSD - \ref{sec:freebsdsrc}
 - Windows - \ref{sec:windowssrc}
\end{itemize}

Once these are done proceed to Section~\ref{sec:cleanup} for cleanup steps.



\subsection{Mint}\label{sec:mintsrc}
These instructions were prepared on Mint 20.3. \newline

\noindent As a preliminary step, you should install the dependencies that *esys-escript* requires from the repository.
\begin{shellCode}
sudo apt-get install python3-dev python3-numpy
sudo apt-get install python3-sympy python3-matplotlib python3-scipy
sudo apt-get install libhdf5-serial-dev
sudo apt-get install libboost-random-dev libboost-python-dev libboost-iostreams-dev
sudo apt-get install scons lsb-release libsuitesparse-dev
\end{shellCode}

\noindent Then navigate to the source directory and execute the following 
\begin{shellCode}
scons -j4 options_file=scons/templates/mint_options.py
\end{shellCode}

\subsection{Arch}\label{sec:archsrc}
These instructions were prepared on Arch Linux. \newline

First, install the dependencies that *esys-escript* uses:
\begin{shellCode}
pacman -Sy --noconfirm python python-numpy python-scipy
pacman -Sy --noconfirm community/hdf5-serial
pacman -Sy --noconfirm extra/boost extra/boost-libs suitesparse
\end{shellCode}

Now you can compile *esys-escript* using the command
\begin{shellCode}
scons options_file=scons/templates/arch_py3_options.py -j4 build_full
\end{shellCode}

\subsection{OpenSuse}\label{sec:susesrc}
These instructions were prepared using OpenSUSE Leap 15.2. \newline

\noindent As a preliminary step, you should install the dependencies that *esys-escript* requires from the repository.
\noindent  If you intend to use Python 2.7, then you should install the following packages
\begin{shellCode}
sudo zypper in python-devel python2-numpy
sudo zypper in python2-scipy python2-sympy python2-matplotlib
sudo zypper in gcc gcc-c++ scons hdf5-devel
sudo zypper in libboost_python-py2_7-1_66_0-devel libboost_numpy-py2_7-1_66_0-devel
sudo zypper in libboost_iostreams1_66_0-devel suitesparse-devel
\end{shellCode}

\noindent If you intend to use Python 3.0, then you should instead install the following packages
\begin{shellCode}
sudo zypper in python3-devel python3-numpy
sudo zypper in python3-scipy python3-sympy python3-matplotlib
sudo zypper in gcc gcc-c++ scons hdf5-devel
sudo zypper in libboost_python-py3-1_66_0-devel libboost_numpy-py3-1_66_0-devel
sudo zypper in libboost_random1_66_0-devel libboost_iostreams1_66_0-devel
sudo zypper in suitesparse-devel
\end{shellCode}

\noindent Now to build *esys-escript* itself.
\noindent In the *esys-escript* source directory execute the following (substitute \textit{opensuse_py2} or \textit{opensuse_py3} as appropriate for XXXX):
\begin{shellCode}
scons -j4 options_file=scons/templates/XXXX_options.py
\end{shellCode}

\noindent If you wish to test your build, you can use the following:
\begin{shellCode}
scons -j4 py_tests options_file=scons/templates/XXXX_options.py
\end{shellCode}

\noindent Now go to Section~\ref{sec:cleanup} for cleanup.

\subsection{CentOS}\label{sec:centossrc}
It is possible to install *esys-escript* on both CentOS releases $7$ and $8$. We include separate instructions for each of these CentOS releases in this section.
\subsubsection{CentOS release $7$}
The core of *esys-escript* works, however some functionality is not available because the default packages for some dependencies in CentOS are too old.
At present, it is not possible to compile *esys-escript* using Python 3.0+ on CentOS $7$ as Python 3.0+ versions of many of the dependencies are not currently available in any of the CentOS repositories, but this may change in the future. 
In this section we only outline how to install a version of *esys-escript* that uses Python 2.7.

\noindent First, add the \texttt{EPEL} repository.
\begin{shellCode}
yum install epel-release.noarch
\end{shellCode}

\noindent Install packages:
\begin{shellCode}
yum install hdf5-devel
yum install python-devel numpy scipy sympy python2-scons
yum install python-matplotlib gcc gcc-c++ boost-devel
yum install boost-python suitesparse-devel
\end{shellCode}

\noindent Now to build *esys-escript* itself.
In the *esys-escript* source directory:
\begin{shellCode}
scons -j4 options_file=scons/templates/centos7_0_options.py
\end{shellCode}

\noindent Now go to Section~\ref{sec:cleanup} for cleanup.

\subsubsection{CentOS release $8$}
The core of *esys-escript* works in CentOS $8$, however some functionality is not available because the default packages for some dependencies in CentOS are too old. This install is for Python 3.

First, add the EPEL, PowerTools and Okay repositories:
\begin{shellCode}
yum update
yum install epel-release.noarch dnf-plugins-core
yum config-manager --set-enabled PowerTools
rpm -ivh http://repo.okay.com.mx/centos/8/x86_64/release/okay-release-1-3.el8.noarch.rpm
yum update
\end{shellCode}

Now, install the packages:
\begin{shellCode}
yum install python3-devel python3-numpy python3-scipy
yum install boost-devel boost-python3 boost-python3-devel
yum install gcc gcc-c++ scons
yum install suitesparse suitesparse-devel 
\end{shellCode}

Finally, you can compile *esys-escript* with the command
\begin{shellCode}
scons -j4 options_file=scons/templates/centos8_0_options.py
\end{shellCode}

\subsection{Fedora}\label{sec:fedorasrc}
These instructions were prepared using Fedora $31$ Workstation.

\noindent To build the a version of *esys-escript* that uses Python 2.7, install the following packages:
\begin{shellCode}
yum install gcc-c++ scons suitesparse-devel
yum install python2-devel boost-python2-devel
yum install python2-scipy
yum install hdf5-devel
\end{shellCode}

\noindent To build the a version of *esys-escript* that uses Python 3.0+, install the following packages:
\begin{shellCode}
yum install gcc-c++ scons suitesparse-devel
yum install python3-devel boost-python3-devel
yum install python3-scipy python3-matplotlib
yum install boost-python3 boost-numpy3 boost-iostreams boost-random
yum install hdf5-devel
\end{shellCode}

\noindent In the source directory execute the following (substitute \textit{fedora_py2} or \textit{fedora_py3} for XXXX):
\begin{shellCode}
scons -j4 options_file=scons/templates/XXXX_options.py
\end{shellCode}

\noindent Now go to Section~\ref{sec:cleanup} for cleanup.

\subsection{MacOS 10.10 and later (macports)}\label{sec:macportsrc}

The following will install the capabilities needed for the \texttt{macports_10.10_options.py} file (later versions can use the same options file).

\begin{shellCode}
sudo port install scons
sudo port select --set python python27
sudo port install boost
sudo port install py27-numpy
sudo port install py27-sympy
sudo port select --set py-sympy py27-sympy
sudo port install py27-scipy
sudo port install hdf5
sudo port install silo
\end{shellCode}

\begin{shellCode}
scons -j4 options_file=scons/templates/macports_10.10options.py
\end{shellCode}


\subsection{MacOS 10.13 and later (homebrew)}\label{sec:homebrewsrc}

The following will install the capabilities needed for the \texttt{homebrew_10.13_options.py} file.

\begin{shellCode}
brew install scons
brew install boost-python
brew install hdf5
\end{shellCode}

There do not appear to be formulae for \texttt{sympy}  so if you wish to use those features, then
you will need to install them separately.


\begin{shellCode}
scons -j4 options_file=scons/templates/homebrew_10.13_options.py
\end{shellCode}


\subsection{FreeBSD}\label{sec:freebsdsrc}

At time of writing, \texttt{numpy} does not install correctly on FreeBSD.
Since \texttt{numpy} is a critical dependency for *esys-escript*, we have been unable to test on FreeBSD.

\begin{comment}
\subsubsection{Release 10.0}

Install the following packages:
\begin{itemize}
 - python
 - scons
 - boost-python-libs
 - bash
 - netcdf
 - silo
 - py27-scipy
 - py27-matplotlib
 - py27-sympy
\end{itemize}

\noindent Next choose (or create) your options file.
For the setup as above the *esys-escript* source comes with a prepared file in
\texttt{scons/templates/freebsd10.0_options.py}.
Finally to build *esys-escript* issue the following in the *esys-escript* source directory
(replace the options file as required):
\begin{shellCode}
scons -j4 options_file=scons/templates/freebsd10.0_options.py
\end{shellCode}

\emph{Note:} Some packages installed above are built with gcc 4.7. Somewhere
in the toolchain a system-installed gcc library is pulled in which is
incompatible with the one from version 4.7 and would prevent *esys-escript* from
executing successfully. As explained in the FreeBSD
documentation\footnote{see \url{http://www.freebsd.org/doc/en/articles/custom-gcc/article.html}}
this can be fixed by adding a line to \texttt{/etc/libmap.conf}:
\begin{shellCode}
libgcc_s.so.1 gcc47/libgcc_s.so.1
\end{shellCode}

\end{comment}
\subsection{Windows}\label{sec:windowssrc}

\noindent These instructions were prepared for Microsoft Windows 10.

\noindent Start by installing *esys-escript* dependencies.

\begin{itemize}
\item Microsoft Visual Studio
\begin{enumerate}
\item Download the Microsoft Visual Studio Community 2017 (or VS 2019 if
preferred) installer from
\begin{itemize}
\item[] \url{https://visualstudio.microsoft.com}.
\end{itemize}
\item Launch the Visual Studio installer, selecting Individual components:
\begin{itemize}
\item VS 2017: \textbf{VC++ 2017 latest v141 tools} \\
VS 2019: \textbf{MSVC v142 - VS 2019 C++ build tools}
\item \textbf{Windows 10 SDK}
\item \textbf{MSBuild}
\item \textbf{Visual C++ tools for CMake}
\end{itemize}
\end{enumerate}
\item Anaconda
\begin{enumerate}
\item Download the Python 3.7 64-Bit Graphical Installer for Windows from
\begin{itemize}
\item[] \url{https://anaconda.org/}.
\end{itemize}
\item Launch the Anaconda installer, selecting installation type: \textbf{Just
Me} and destination folder: \newline \verb!C:\Users\%USERNAME%\Anaconda3!.
\end{enumerate}
\end{itemize}

\noindent Next, open Windows Command Prompt (\file{cmd.exe}) and set-up the
*esys-escript* dependencies.

\begin{itemize}
\item Conda environment
\begin{enumerate}
\item Create and activate a new environment
\begin{shellCode}
C:\Users\%USERNAME%\Anaconda3\Scripts\activate
conda create --name *esys-escript* python=3.7
conda deactivate
C:\Users\%USERNAME%\Anaconda3\Scripts\activate *esys-escript*
\end{shellCode}
\item Install required conda modules
\begin{shellCode}
conda install numpy==1.15.4 matplotlib==2.2.2 sympy==1.1.1
    boost git scons scipy m2-patch mumps gmsh
    -c defaults -c conda-forge
\end{shellCode}
\end{enumerate}
\item Vcpkg
\begin{enumerate}
\item Build vcpkg package manager
\begin{shellCode}
cd C:\Users\%USERNAME%
git clone https://github.com/Microsoft/vcpkg.git
cd vcpkg
bootstrap-vcpkg
\end{shellCode}
\item Install the CppUnit vcpkg package
\begin{shellCode}
vcpkg install cppunit:x64-windows
\end{shellCode}
\end{enumerate}
\end{itemize}

\noindent Once the dependencies are installed and set-up, you can download and
build *esys-escript* from source.

\begin{enumerate}
\item Activate the conda environment (if not active).
\begin{shellCode}
C:\Users\%USERNAME%\Anaconda3\Scripts\activate *esys-escript*
\end{shellCode}
\item Set-up the Command Prompt for the 64-bit MSVC command line build environment.
\begin{shellCode}
"C:\Program Files (x86)\Microsoft Visual Studio\2017\
    Community\VC\Auxiliary\Build\vcvars64.bat"
\end{shellCode}
\item Add CppUnit to the Windows System Path.
\begin{shellCode}
set PATH=%PATH%;C:\Users\%USERNAME%\vcpkg\packages\cppunit_x64-windows\bin
\end{shellCode}
\item Download the *esys-escript* source code tarball from
\begin{itemize}
\item[] \url{https://launchpad.net/*esys-escript*-finley}
\end{itemize}
Extract the tarball to \verb!C:\Users\%USERNAME%*esys-escript*!
\item Build and install the netCDF-4 C++ library before starting the *esys-escript*
build.  Download the netCDF-4 C++ v4.3.1 source code tarball from
\begin{itemize}
\item[] \url{https://github.com/Unidata/netcdf-cxx4/archive/v4.3.1.tar.gz}
\end{itemize}
Extract the tarball to \verb!C:\Users\%USERNAME%*esys-escript*!
\item Apply the provided patch.
\begin{shellCode}
cd C:\Users\%USERNAME%*esys-escript*\netcdf-cxx4-4.3.1
patch -p1 < ..\src\tools\anaconda\Anaconda3\netcdf-cxx4.patch
\end{shellCode}
\item Configure, build, and install netcdf-cxx4.
\begin{shellCode}
mkdir build
cd build
cmake -G "Visual Studio 15 2017 Win64" -DBUILD_SHARED_LIBS=OFF
    -DCMAKE_INSTALL_PREFIX="%CONDA_PREFIX%\Library"
    -DCMAKE_LIBRARY_PATH="%CONDA_PREFIX%\Library\lib"
    -DCMAKE_PREFIX_PATH="%CONDA_PREFIX%\Library"
    -DNETCDF_LIB_NAME="netcdf" -DHDF5_LIB_NAME="hdf5" ..
cmake --build . --config Release
cmake --build . --config Release --target install
\end{shellCode}
\item Kick-off the *esys-escript* build when the netCDF-4 C++ install is complete.
\begin{shellCode}
cd C:\Users\%USERNAME%*esys-escript*\src
scons -j4 build_full options_file=scons/templates/windows_msvc141_options.py
\end{shellCode}
\item Once the build completes successfully, you can validate *esys-escript* using
the provided test script.
\begin{shellCode}
python utest.py C:\Users\%USERNAME%*esys-escript*\src\build -t8
\end{shellCode}
\end{enumerate}

\subsection{Other Systems / Custom Builds}\label{sec:othersrc}

*esys-escript* has support for a number of optional packages.
Some, like \texttt{netcdf} need to be enabled at compile time, while others, such as \texttt{sympy} and the projection packages
used in \downunder are checked at run time.
For the second type, you can install them at any time (ensuring that python can find them) and they should work.
For the first type, you need to modify the options file and recompile with scons.
The rest of this section deals with this.

To avoid having to specify the options file each time you run scons, copy an existing \texttt{_options.py} file from the
\texttt{scons/} or \texttt{scons/templates/} directories. Put the file in the \texttt{scons} directory and name
it \textit{yourmachinename}\texttt{_options.py}.\footnote{If the name
has - or other non-alpha characters, they must be replaced with underscores in the filename}.
For example: on a machine named toybox, the file would be \texttt{scons/toybox_options.py}.

Individual lines can be enabled/disabled, by removing or adding \# (the python comment character) to the beginning of the line.
For example, to enable OpenMP, change the line
\begin{verbatim}
#openmp = True
\end{verbatim}
to
\begin{verbatim}
openmp = True
\end{verbatim}

If you are using libraries which are not installed in the standard places (or have different names) you will need to
change the relevant lines.
A common need for this would be using a more recent version of the boost::python library.
You can also change the compiler or the options passed to it by modifying the relevant lines.

\subsubsection*{MPI}
If you wish to enable or disable MPI, or if you wish to use a different implementation of MPI, you can use the \texttt{mpi}
configuration variable.
You will also need to ensure that the \texttt{mpi_prefix} and \texttt{mpi_libs} variables are uncommented and set correctly.
To disable MPI use, \verb|mpi = 'none'|.


\subsubsection{Testing}
As indicated earlier, you can test your build using \texttt{scons py_tests}.
Note however, that some features like \texttt{netCDF} are optional for using *esys-escript*, the tests will report a failure if
they are missing.

\section{Cleaning up}
\label{sec:cleanup}

Once the build (and optional testing) is complete, you can remove everything except:
\begin{itemize}
 - bin
 - esys
 - lib
 - doc
 - CREDITS
 - LICENSE
 - README
\end{itemize}
The last three aren't strictly required for operation.
The \texttt{doc} directory is not required either but does contain examples of *esys-escript* scripts.

You can run *esys-escript* using \texttt{\textit{path_to_*esys-escript*_files}/bin/run-*esys-escript*}.\\
Where \texttt{\textit{path_to_*esys-escript*_files}} is replaced with the real path.

\begin{optionalstep}
You can add the *esys-escript* \texttt{bin} directory to your \texttt{PATH} variable.
The launcher will then take care of the rest of the environment.
\end{optionalstep}

\section{Optional Extras}

Some other packages which might be useful include:
\begin{itemize}
 - Lapack and UMFPACK --- direct solvers (install the relevant libraries and enable them in the options file).
 - support for the Silo file format (install the relevant libraries and enable them in the options file).
 - VisIt --- visualisation package. Can be used independently but our \texttt{weipa} library can make a Visit
plug-in to allow direct visualisation of *esys-escript* simulations.
 - gmsh --- meshing software used by our \texttt{pycad} library.
 - Mayavi2 --- another visualisation tool.
\end{itemize}


\subsection{Trilinos}
*esys-escript* now has some support for Trilinos\footnote{\url{https://trilinos.org/}}
solvers and preconditioners.
The most significant limitation is that the current Trilinos release does not
support block matrices so *esys-escript* can only use Trilinos solvers for single
PDEs (i.e. no PDE systems).

If your distribution does not provide Trilinos packages you can build a working
version from source. (See Appendix~\ref{app:trilinos})


\section{Testing *esys-escript*}\label{chap:utest}

*esys-escript* has extensive testing that can be used to confirm that the program is working correctly. To run the unit testing, compile *esys-escript* with the flag \texttt{build_full}. This will build *esys-escript* normally and then create a shell script named \texttt{utest.sh}. Once this file has been created, you can run unit testing using the command
\begin{shellCode}
sh utest.sh path_to_build_folder '-tT -nN -pP'
\end{shellCode}
where \texttt{T}, \texttt{N} and \texttt{P} represent the number of threads, nodes and processes to run the testing on. Some of these terms can be omitted. For example, to run the testing in serial, you would run
\begin{shellCode}
sh utest.sh path_to_build_folder '-t1'
\end{shellCode}

Note that a careless selection of these parameters may cause the testing program to skip many of the tests. For example, if you compile *esys-escript* with OpenMP enabled but then instruct the testing program to run on a single thread, many of the OpenMP tests will not be run.


\esysappendix
## Required compiler features
\label{app:cxxfeatures}

Building *esys-escript* from source requires that your c++ compiler supports at least the following features:
\begin{itemize}
 - \texttt{std::complex<>}
 - Variables declared with type \texttt{auto}
 - Variables declared with type \texttt{decltype(T)}
 - \texttt{extern template class} to prevent instantiation of templates. 
 - \texttt{template class \textit{classname$<$type$>$};} to force instantiation of templates
 - \texttt{isnan()} is defined in the \texttt{std::} namespace
\end{itemize}
The above is not exhaustive and only lists language features which are more recent that our previous baseline of c++99 (or which
we have recently begun to rely on).
Note that we test on up to date versions of \texttt{g++, icpc \& clang++} so they should be fine.

Note that in future we may use c++14 features as well.

## Trilinos}
\label{app:trilinos}

In order to solve PDEs with complex coefficients, *esys-escript* needs to be compiled with \texttt{Trilinos} support.
This requires that your version of Trilinos has certain features enabled.
Since some precompiled distributions of \texttt{Trilinos} are not built with these features, you may 
need to compile \texttt{Trilinos} yourself as well.

While we can't provide support for building \texttt{Trilinos}, we provide here two configuration files which seem to work for 
Debian 10 ``buster'. One of these cmake script builds \texttt{Trilinos} with MPI support and one builds \texttt{Trilinos} without MPI support.

\section{Debian ``buster'' configuration}


\subsection{Required packages}

The following packages should be installed to attempt this build:
\begin{itemize}
\item[] cmake
\item[] g++
\item[] libsuitesparse-dev
\item[] libmumps-dev
\item[] libboost-dev
\item[] libscotchparmetis-dev
\item[] libmetis-dev
\item[] libcppunit-dev
\end{itemize}
\end{document}

