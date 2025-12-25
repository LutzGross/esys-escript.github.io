import os, re, subprocess

rep1 = re.compile('^\s*(?P<lib>[^ ]+) => (?P<loc>[^ ]+) [(](?P<addr>[^)]+)[)]$')
rep2 = re.compile('^\s*(?P<lib>[^ ]+) [(](?P<addr>[^)]+)[)]$')
rep3 = re.compile('^\s*(?P<lib>[^ ]+) => not found$')

def which(cmd):
    hits = []
    paths = os.environ['PATH'].split(os.pathsep)
    for path in paths:
        hit = os.path.join(path, cmd)
        if os.path.isfile(hit):
            hits.append(hit)
    return hits

def getLibDependents(lib):
    if os.name == 'nt':
        lib_pth = which(lib)[0]
        txt = subprocess.check_output(['dumpbin', '/dependents', lib_pth]).decode()
        lines = txt.splitlines()
        libs = []
        have_deps = False
        for line in lines:
            if line.endswith('dependencies:'):
                have_deps = True
                continue
            if have_deps:
                dep = line.strip()
                if dep:
                    if dep.startswith('api-ms-win-'):
                        # ignore these
                        continue
                    try:
                        libs.append(which(dep)[0])
                    except:
                        print(txt)
                        raise Exception('lib not found: ' + dep)
                elif len(libs):
                    # blank line with libs found means end of list
                    break
    else:
        txt = subprocess.check_output(['ldd', lib]).decode()
        lines = txt.splitlines()
        libs = []
        for line in lines:
            m = rep1.match(line)
            if m:
                libs.append(os.path.abspath(m.group('loc')))
                continue
            m = rep2.match(line)
            if m:
                # virtual kernel lib
                continue
            m = rep3.match(line)
            if m:
                print('lib not found: ' + m.group('lib'))
                continue
            raise Exception('line not recognised:\n' + line)
    return libs

def searchLibs(lib):
    ''' recursive search for lib dependents '''
    if os.name == 'nt':
        lib_pth = which(lib)[0]
    else:
        lib_pth = os.path.abspath(lib)
    if lib_pth not in libDict:
        libDict[lib_pth] = getLibDependents(lib_pth)
        for l in libDict[lib_pth]:
            searchLibs(l)

esys_linux_libs = [
    'libescript.so',
    'libescriptreader.so',
    'libfinley.so',
    'libpaso.so',
    'libripley.so',
    'libspeckley.so',
    'libweipa.so'
]
esys_win_libs = [
    'escript.dll',
    'escriptreader.dll',
    'finley.dll',
    'paso.dll',
    'ripley.dll',
    'speckley.dll',
    'weipa.dll',
    'libdmumps.dll', # mumps dlls linked at run time
    'libzmumps.dll'
]

libDict = {}

if os.name == 'nt':
    esys_libs = esys_win_libs
else:
    esys_libs = esys_linux_libs

if not os.path.exists(esys_libs[0]):
    raise Exception('lib not found: ' + esys_libs[0] +
        '\nexpect esys libs to be in the current directory')

print('searching for esys lib dependents...\n')
for i in esys_libs:
    searchLibs(i)
print('\n'.join(sorted(list(libDict.keys()))))
print('\ndependent libs found: {}'.format(len(libDict)))
