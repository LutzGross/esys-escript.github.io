
##############################################################################
#
# Copyright (c) 2013-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

import os, platform, re, setuptools, setuptools.command.build_ext, setuptools.command.install, sys

esys_version = '5.6.2'
python_version = '%d%d' % sys.version_info[:2]
python_version_site = '%d.%d' % sys.version_info[:2]
system, machine = platform.system().lower(), platform.machine().lower()
plat, dist_name, dist_ver = (None, None, None)
url = 'https://github.com/esys-escript/esys-escript.github.io/releases/download'

if 'linux' in system and '64' in machine:
    plat = 'Linux64'
    py_lib_ext = 'so'
    url_lib_ext = 'tgz'
    try:
        with open('/etc/os-release', 'r') as f:
            os_rel = {k:v.strip('"') for k, v in (l.strip().split('=') for l in f)}
        dist_name, dist_ver = (os_rel['NAME'], int(os_rel['VERSION_ID']))
    except:
        pass
    if 'debian' in dist_name.lower():
        dist_name = 'debian'
    if not (dist_name == 'debian' and dist_ver == 10):
        raise Exception('Linux distro not supported: {} {} - please try Debian 10'.format(dist_name, dist_ver))
elif 'windows' in system and '64' in machine:
    plat = 'Windows64'
    py_lib_ext = 'pyd'
    url_lib_ext = 'zip'
    dist_name = 'windows'
    dist_ver = sys.getwindowsversion()[0]
    if not (dist_ver == 10):
        raise Exception('Windows version not supported: {} {} - please try Windows 10'.format(dist_name, dist_ver))
# elif 'darwin' in system:
#     plat = 'MacOSX'

esys_bin = '-'.join([plat, 'py'+python_version, dist_name, str(dist_ver), 'lib'])

# check we have the required binaries
if plat is None:
    raise Exception('Platform not supported: py{}-{}-{}'.format(python_version, system, machine))

esys_lib = 'esys_escript_lib'

class esys_build_ext(setuptools.command.build_ext.build_ext):
    ''' add platform specific python libs (.so/.pyd) to python module '''

    def _download(self):
        import requests, tarfile, zipfile
        fname = esys_bin + '.' + url_lib_ext
        _url = '/'.join([url, esys_version, fname])
        if not os.path.exists(fname):
            print('downloading {}...'.format(_url))
            libs = requests.get(_url, allow_redirects=True, verify=False)
            with open(fname, 'wb') as f:
                f.write(libs.content)
        print('extracting {}...'.format(fname))
        tar = tarfile.open(fname) if url_lib_ext == 'tgz' else zipfile.ZipFile(fname, 'r')
        tar.extractall()

    def build_extension(self, ext):
        self._download()
        for f in self.get_source_files():
            src = os.path.join(ext.name, f)
            dst = os.path.join(self.build_lib, f)
            self.copy_file(src, dst)

class esys_install(setuptools.command.install.install):
    ''' install platform specific non-python libs in a separate site-packages folder '''

    def _updateBuildvars(self, path):
        prefix_rep = re.compile('^(prefix=).*$')
        python_rep = re.compile('^(python=).*$')
        inc_path_rep = re.compile('^(.*_inc_path=).*$')
        lib_path_rep = re.compile('^(.*_lib_path=).*$')
        fn_1 = path
        fn_2 = path + '.in'

        with open(fn_1) as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            line = line.strip()
            m = prefix_rep.match(line)
            if m:
                lines[i] = m.group(1) + '%(prefix)s\n'
            m = python_rep.match(line)
            if m:
                lines[i] = m.group(1) + '%(python)s\n'
            m = inc_path_rep.match(line)
            if m:
                lines[i] = m.group(1) + '%(inc_prefix)s\n'
            m = lib_path_rep.match(line)
            if m:
                lines[i] = m.group(1) + '%(lib_prefix)s\n'

        with open(fn_2, 'w') as f:
            f.writelines(lines)
        
    def run(self):
        if plat.startswith('Windows'):
            site_path = 'lib/site-packages'
        else:
            site_path = 'lib/python'+python_version_site+'/site-packages'
        lib_path = os.path.join(esys_bin, 'lib')
        self._updateBuildvars(os.path.join(lib_path, 'buildvars'))
        data_files = [
            (
                os.path.join(site_path, esys_lib),
                [os.path.join(dp, f) for dp, dn, fns in os.walk(lib_path) for f in fns if f != 'buildvars']
            )
        ]
        self.distribution.include(data_files=data_files)
        setuptools.command.install.install.run(self)

install_requires = ['netCDF4', 'pyproj', 'scipy']
if python_version == '37':
    if plat.startswith('Windows'):
        install_requires.append('numpy==1.15.4')
    else:
        install_requires.append('numpy')
    install_requires.append('sympy==1.1.1')
else:
    install_requires.append('numpy')
    install_requires.append('sympy')

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name='esys-escript',
    version=esys_version,
    description='numerical modelling library using the finite element method',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://esys-escript.github.io',
    python_requires='>=3.7',
    packages=setuptools.find_packages(),
    install_requires=install_requires,
    entry_points = {
        'console_scripts': [
            'run-escript=esys.escript.run_escript:main',
            'runmodel=esys.escriptcore.runmodel:main'
        ]
    },
    ext_modules = [
        setuptools.Extension(esys_bin, [
            'esys/dudley/dudleycpp.'+py_lib_ext,
            'esys/escriptcore/escriptcpp.'+py_lib_ext,
            'esys/finley/finleycpp.'+py_lib_ext,
            'esys/ripley/ripleycpp.'+py_lib_ext,
            'esys/speckley/speckleycpp.'+py_lib_ext,
            'esys/weipa/weipacpp.'+py_lib_ext,
        ]),
    ],
    cmdclass = {
        'build_ext': esys_build_ext,
        'install': esys_install
    }
)

