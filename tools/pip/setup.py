
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
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
url = 'https://github.com/esys-escript/esys-escript.github.io/releases/download'

if 'linux' in system and '64' in machine:
    plat = 'Linux64'
    py_lib_ext = 'so'
    url_lib_ext = 'tgz'
elif 'windows' in system and '64' in machine:
    plat = 'Windows64'
    py_lib_ext = 'pyd'
    url_lib_ext = 'zip'
# elif 'darwin' in system:
#     plat = 'MacOSX'

# check we have the required binaries
if plat is None:
    raise Exception('Platform not supported: py{}-{}-{}'.format(python_version, system, machine))

esys_bin = '-'.join([plat, 'py'+python_version, 'lib'])
esys_lib = 'esys_escript_lib'

class esys_install(setuptools.command.install.install):
    ''' install platform specific libraries '''

    def _download(self):
        import requests, tarfile, zipfile
        fname = esys_bin + '.' + url_lib_ext
        _url = '/'.join([url, esys_version, fname])
        if not os.path.exists(fname):
            print('downloading {}...'.format(_url))
            resp = requests.get(_url, allow_redirects=True, verify=False)
            if resp.status_code != requests.codes.ok:
                if resp.status_code == requests.codes.not_found:
                    print('\nescript platform-specific libs not found: ' + fname)
                    print('please see available libs at:\n' + url + '\n')
                resp.raise_for_status()
            with open(fname, 'wb') as f:
                f.write(resp.content)
        print('extracting {}...'.format(fname))
        tar = tarfile.open(fname) if url_lib_ext == 'tgz' else zipfile.ZipFile(fname, 'r')
        tar.extractall()

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
        self._download()
        if plat.startswith('Windows'):
            site_path = 'lib/site-packages'
        else:
            site_path = 'lib/python'+python_version_site+'/site-packages'

        # add platform specific python libs (.so/.pyd) to python module
        lib_path = os.path.join(esys_bin, 'esys')
        esys_bin_strip_len = len(esys_bin) + 1
        data_files = [
            (os.path.join(site_path, dp[esys_bin_strip_len:]), [os.path.join(dp, f)])
                for dp, dn, fns in os.walk(lib_path) for f in fns
        ]

        # add platform specific non-python libs to a separate site-packages folder
        lib_path = os.path.join(esys_bin, 'lib')
        self._updateBuildvars(os.path.join(lib_path, 'buildvars'))
        data_files.append(
            (
                os.path.join(site_path, esys_lib),
                [os.path.join(dp, f) for dp, dn, fns in os.walk(lib_path) for f in fns if f != 'buildvars']
            )
        )
        self.distribution.include(data_files=data_files)
        setuptools.command.install.install.run(self)

install_requires = ['matplotlib', 'netCDF4', 'scipy']
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
    cmdclass = { 'install': esys_install },
    setup_requires = ['requests']
)

