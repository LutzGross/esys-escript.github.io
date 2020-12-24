# Escript pip module

The `esys-escript` pip module includes Python source files, as listed in
`py_src.lst`, and dependent libraries, as collected from an `esys-escript`
conda-forge installation for the required system. The Python packaging build
process expects to download the required package of libraries from
[GitHub](https://github.com/esys-escript/esys-escript.github.io/releases)

## Collecting the Python source and dependent libraries

The Python sources files listed in `py_src.lst`, along with the escript
README.md and LICENSE files, need to be collected and placed in the same folder
as the `setup.py` script.

Collecting the required libraries from the conda-forge installation involves
manually following the dependency trail of the required libraries.

On Linux systems use:

```
$ ldd libescript.so
        linux-vdso.so.1
        libgomp.so.1
        libpython3.8.so.1.0
        libboost_python38.so.1.72.0
        libboost_numpy38.so.1.72.0
        libstdc++.so.6
...
```

On Windows systems use (`dumpbin` is installed with MSVC):

```
>dumpbin /dependents escript.dll
...
File Type: DLL

  Image has the following dependencies:

    netcdf.dll
    python38.dll
    boost_numpy38.dll
    boost_python38.dll
    MSVCP140.dll
...
```

The sample output shows Boost Python and NumPy dependencies that should be
checked for their sub-dependencies. The process should continue until all
dependencies are identified. Any libraries that can be safely assummed to be
available can be ignored (e.g. `MSVCP140.dll` on Windows systems).

The escript `buildvars` file from the conda-forge installation should also be
included with the libraries. Libraries are packaged in the `lib` folder in a
`zip` file on Windows systems or a `tgz` file on Linux systems.

Python-specific escript libraries, such as `esys\escriptcore\escriptcpp.pyd`
are packaged under the `esys` folder.

## Build process

The pip module build script `setup.py` uses the Python setuptools API. The
highlevel steps in `setup.py` are:

- Check the platform is supported. Currently supported platforms are Windows 10
and Debian 10.

- Download and extract the dependent libraries from
[GitHub](https://github.com/esys-escript/esys-escript.github.io/releases). If
the required `zip`/`tgz` file is already in the working directory it will be
used rather than downloading, which is useful when you are working on a new
release.

- Convert the `buildvars` file to a template that can be populated by
`esys.escript.__init__.py` on the end-user's system.

- Configure the dependent libraries and `buildvars` template to be installed in
the end-user's `site-packages` folder, under `esys_escript_lib`. The escript
libraries must be in the system library path, so `esys.escript.__init__.py` on
the end-user's system will issue a warning if not set and restart Python with
the library path set.

- Set up `run-escript` and `runmodel` script entry points, runable from the
end-user's command line path.

The pip module build script dependencies need to be installed before launching
the build:

```
pip install requests
```

The module build can then be launched as follows:

```
python setup.py bdist_wheel
```

The generated wheel is specific to the current platform, so repeat the process
on each required target platform.

## Installation from wheel file

The pip build generates a wheel file in the `dist` folder. You can install the
pip module from the local wheel file, for example:

```
pip install dist\esys_escript-5.6.2-cp38-cp38-win_amd64.whl
```

## Uploading wheels to Test PyPI

To upload the generated wheel files to Test PyPI, you need to get a Test PyPI
API token and setup your `$HOME/.pypirc` as follows:

```
[testpypi]
  username = __token__
  password = <your token here>
```

See [PyPI Help](https://pypi.org/help) for further details.

Twine needs to be installed next:

```
pip install twine
```

You can then upload your wheels to Test PyPI using Twine:

```
python -m twine upload --repository testpypi dist/*
```

## Installation from Test PyPI

Avoid installing the `esys-escript` Python dependencies from Test PyPI by
installing them manually as the first step:

```
pip install netCDF4 pyproj scipy sympy
```

Then install `esys-escript` from Test PyPI as follows:

```
pip install --index-url https://test.pypi.org/simple/ --no-deps esys-escript
```
