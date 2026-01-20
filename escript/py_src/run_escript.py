#!/usr/bin/env python
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

import argparse, esys.escript, os, re, sys, subprocess

def which(cmd):
    hits = []
    paths = env['PATH'].split(os.pathsep)
    for path in paths:
        hit = os.path.join(path, cmd)
        if os.path.isfile(hit):
            hits.append(hit)
    return hits

def die(msg):
    print('Error: {msg}'.format(msg=msg), file = sys.stderr)
    exit(1)

def main():
    # Escript wrapper for python
    # Sets LD_LIBRARY_PATH and PYTHONPATH and then runs either python or the MPI
    # launcher.

    # Extra paths can be configured about a page further down
    # Search for EXTRA_PATH = ''

    def vlog(msg):
        if ESCRIPT_VERBOSE != 0:
            print(msg)

    def printExports():
        cmd, env1, env2 = ('set', '%', '%') if IS_WINDOWS else ('export', '$', '')
        if IS_WINDOWS:
            if EXTRA_PATH or EXTRA_LD_LIBRARY_PATH:
                print(cmd+' PATH='+os.pathsep.join(filter(None, [EXTRA_PATH, EXTRA_LD_LIBRARY_PATH, env1+'PATH'+env2])))
        else:
            if EXTRA_LD_LIBRARY_PATH:
                print(cmd+' LD_LIBRARY_PATH='+os.pathsep.join(filter(None, [EXTRA_LD_LIBRARY_PATH,
                    env1+'LD_LIBRARY_PATH'+env2])))
            if EXTRA_PATH:
                print(cmd+' PATH='+os.pathsep.join(filter(None, [EXTRA_PATH, env1+'PATH'+env2])))
        if EXTRA_PYTHONPATH:
            print(cmd+' PYTHONPATH='+os.pathsep.join(filter(None, [EXTRA_PYTHONPATH, env1+'PYTHONPATH'+env2])))
        if IS_DARWIN:
            if EXTRA_DYLD_LIBRARY_PATH:
                print('export DYLD_LIBRARY_PATH = '+os.pathsep.join(filter(None, [EXTRA_DYLD_LIBRARY_PATH,
                    EXTRA_LD_LIBRARY_PATH, '$DYLD_LIBRARY_PATH'])))
        exit(0)

    def get_buildvar(bldvar):
        with open(buildinfo_file, 'r') as f:
            m = re.search(bldvar+'[ ]*=(?P<bldvar>.*)$', f.read(), re.MULTILINE)
            if m:
                return m.group('bldvar').strip()

    def fixEnvVars(cmd):
        # workaround for pip module using pre/post/launch vars generated on wrong platform
        if IS_WINDOWS:
            res = '[$]{([^}]*)}'
            env1, env2 = ('%', '%')
        else:
            res = '%([^%]*)%'
            env1, env2 = ('${', '}')
        rep = re.compile(res)
        return rep.sub(env1+r'\1'+env2, cmd)

    # set to 1 if this is part of a packaged build (.deb) and files will be
    # installed in standard locations rather than everything in a single directory
    stdlocation = '0'

    # Now we find the location of this script
    # Note that this location should be absolute but does not need to be unique
    scriptdir = ''
    curdir = os.getcwd()
    buildinfo_file = None
    pip_lib_prefix = None

    # Environment vars which control operations:
    # ESCRIPT_NUM_NODES, ESCRIPT_NUM_PROCS, ESCRIPT_NUM_THREADS, ESCRIPT_HOSTFILE, ESCRIPT_CREATESTDFILES

    IS_WINDOWS = (os.name == 'nt')
    env = os.environ
    if IS_WINDOWS:
        tempdir, username = env['TEMP'], env['USERNAME']
        IS_DARWIN = False
    else:
        tempdir, username = '/tmp', env['USER']
        IS_DARWIN = (os.uname()[0] == 'Darwin')
    hostfile = os.path.join(tempdir, 'escript.'+username+'.'+str(os.getpid()))

    #Begin finding ESCRIPT_ROOT
    if stdlocation != '0':
        #Package building scripts will replace this line
        ESCRIPT_ROOT = '/usr/lib/python-escript'
    elif esys.escript.__buildinfo_file__ is not None:
        # pip install will set __buildinfo_file__
        buildinfo_file = esys.escript.__buildinfo_file__
        ESCRIPT_ROOT = get_buildvar('prefix')
        pip_lib_prefix = os.path.join(ESCRIPT_ROOT, 'esys_escript_lib')
    else:
        # We don't know the escript root so we need to work it out from the invocation
        # Need to match if the name contains /
        if sys.argv[0].find(os.sep):
            # We are not using the PATH to find the script
            scriptdir = os.path.dirname(os.path.abspath(sys.argv[0]))
        else:
            # name does not contain / therefore we are using
            scriptdirs = which(sys.argv[0])
            if len(scriptdirs) < 1:
                die('Unable to determine script directory!')
            scriptdir = os.path.dirname(scriptdirs[0])
        ESCRIPT_ROOT = os.path.dirname(scriptdir)
        ESCRIPT_PARENT = os.path.dirname(ESCRIPT_ROOT)

    ##### End finding ESCRIPT_ROOT ########

    # if possible please express paths relative to $ESCRIPT_ROOT unless
    # they are in an unrelated location

    EXTRA_DYLD_LIBRARY_PATH = ''
    if pip_lib_prefix is None:
        EXTRA_PYTHONPATH = ESCRIPT_ROOT
        EXTRA_PATH = os.path.join(ESCRIPT_ROOT, 'bin')
        if IS_WINDOWS:
            EXTRA_LD_LIBRARY_PATH = ''
        else:
            EXTRA_LD_LIBRARY_PATH = os.path.join(ESCRIPT_ROOT, 'lib')
    else:
        EXTRA_PYTHONPATH = ''
        if IS_WINDOWS:
            EXTRA_PATH = pip_lib_prefix
            EXTRA_LD_LIBRARY_PATH = ''
        else:
            EXTRA_PATH = ''
            EXTRA_LD_LIBRARY_PATH = pip_lib_prefix

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', metavar='nn', help='number of nodes to use')
    parser.add_argument('-p', metavar='np', help='number of MPI processes to spawn per node')
    parser.add_argument('-t', metavar='nt', help='number of OpenMP threads to use')
    parser.add_argument('-f', metavar='file', help='name of MPI hostfile')
    parser.add_argument('-c', action='store_true', help='print compile information for escript and exit')
    parser.add_argument('-V', action='store_true', help='print escript version and exit')
    parser.add_argument('-i', action='store_true', help='interactive mode')
    parser.add_argument('-b', action='store_true', help='do not invoke python (run non-python programs)')
    parser.add_argument('-e', action='store_true', help='print export statements for environment and exit')
    parser.add_argument('-o', action='store_true', help='redirect output from MPI to files')
    parser.add_argument('-v', action='store_true', help='print diagnostics')
    parser.add_argument('-x', action='store_true', help='run in new xterm instance')
    parser.add_argument('-m', metavar='tool', help='run with valgrind {tool=m[emcheck]/c[allgrind]/[cac]h[egrind]}')
    parser.add_argument('script.py [args...]', nargs=argparse.REMAINDER, help='Your python script and optional command-line arguments')
    args = parser.parse_args()

    if buildinfo_file is None:
        if IS_WINDOWS:
            buildinfo_dir = 'bin'
        else:
            buildinfo_dir = 'lib'
        buildinfo_file = os.path.join(ESCRIPT_ROOT, buildinfo_dir, 'buildvars')
    if not os.access(buildinfo_file, os.R_OK):
        if args.e:
            printExports()
        die('Unable to read escript build information!')

    PYTHON_MPI_NULL = os.path.join(ESCRIPT_ROOT, 'lib', 'esys', 'pythonMPI')
    PYTHON_MPI_REDIRECT = os.path.join(ESCRIPT_ROOT, 'lib', 'esys', 'pythonMPIredirect')
    PYTHON_CMD = get_buildvar('python')

    #==============================================================================

    # Parse the command-line options

    DO_BINARY = args.b
    DO_VALGRIND = args.m
    ESCRIPT_NUM_NODES = 1 if args.n is None else int(args.n)
    ESCRIPT_NUM_PROCS = 1 if args.p is None else int(args.p)
    ESCRIPT_NUM_THREADS = args.t if args.t is None else int(args.t)
    ESCRIPT_HOSTFILE = args.f
    DO_INTERACTIVE = args.i
    ESCRIPT_CREATESTDFILES = args.o
    ESCRIPT_VERBOSE = args.v
    DO_XTERM = args.x

    if args.c:
        with open(buildinfo_file, 'r') as f:
            print(f.read())
        exit(0)
    if args.V:
        print('escript-development(build '+get_buildvar('svn_revision')+')')
        exit(0)
    if args.e:
        printExports()
        exit(0)

    PY_SCRIPT = ' '.join(args.__dict__.get('script.py [args...]', ''))

    #==============================================
    #
    #   Read MPI_FLAVOUR and WITH_OPENMP from the buildvars
    #
    MPI_FLAVOUR = get_buildvar('mpi')
    WITH_OPENMP = get_buildvar('openmp')

    vlog('MPI flavour is '+MPI_FLAVOUR+'.')
    if WITH_OPENMP == '1':
        vlog('OpenMP enabled.')
    else:
        vlog('OpenMP disabled.')

    #
    #   Add VisIt paths if required
    #
    WITH_VISIT = get_buildvar('visit')
    if WITH_VISIT == '1':
        visit_bins = which('visit')
        if len(visit_bins) > 0:
            visit_bin = visit_bins[0]
            m = re.search('LIBPATH[ ]*=(?P<libpath>.*)$', subprocess.Popen(visit_bin+' -env', shell=True,
                stdout=subprocess.PIPE).stdout.read().decode(sys.stdout.encoding), re.MULTILINE)
            VISIT_PY_PATH = m.group('libpath').strip()
            EXTRA_PYTHONPATH = os.pathsep.join(filter(None, [EXTRA_PYTHONPATH, VISIT_PY_PATH]))
            EXTRA_LD_LIBRARY_PATH = os.pathsep.join(filter(None, [EXTRA_LD_LIBRARY_PATH, VISIT_PY_PATH]))
        else:
            vlog('Warning: VisIt module enabled but VisIt not in path!')

    #
    #  extend path variables
    #
    export_env = env
    if IS_WINDOWS:
        export_env['PATH'] = os.pathsep.join(filter(None, [EXTRA_PATH, EXTRA_LD_LIBRARY_PATH, export_env.get('PATH', '')]))
    else:
        export_env['PATH'] = os.pathsep.join(filter(None, [EXTRA_PATH, export_env.get('PATH', '')]))
        export_env['LD_LIBRARY_PATH'] = os.pathsep.join(filter(None, [EXTRA_LD_LIBRARY_PATH,
            export_env.get('LD_LIBRARY_PATH', '')]))
    export_env['PYTHONPATH'] = os.pathsep.join(filter(None, [EXTRA_PYTHONPATH, export_env.get('PYTHONPATH', '')]))
    if IS_DARWIN:
        export_env['DYLD_LIBRARY_PATH'] = os.pathsep.join(filter(None, [EXTRA_DYLD_LIBRARY_PATH, EXTRA_LD_LIBRARY_PATH,
            export_env.get('DYLD_LIBRARY_PATH', '')]))
    vlog('PATH='+export_env['PATH']+'\nLD_LIBRARY_PATH='+export_env.get('LD_LIBRARY_PATH','None')+
        '\nPYTHONPATH='+export_env['PYTHONPATH'])
    if IS_DARWIN:
        vlog('DYLD_LIBRARY_PATH='+export_env['DYLD_LIBRARY_PATH'])

    #==============================================
    #
    #  Ensure the variables have sensible values
    #
    if MPI_FLAVOUR == 'none':
        if ESCRIPT_NUM_NODES is not None:
            if ESCRIPT_NUM_NODES > 1:
                print('Warning: MPI disabled but number of nodes set. Option ignored.')
        if ESCRIPT_NUM_PROCS is not None:
            if ESCRIPT_NUM_PROCS > 1:
                print('Warning: MPI disabled but number of processors per node set. Option ignored.')
        if ESCRIPT_HOSTFILE is not None:
            print('Warning: MPI disabled but host file is given. Option ignored.')
        ESCRIPT_NUM_NODES = 1
        ESCRIPT_NUM_PROCS = 1
    else:
        # use the PBS_NODEFILE if not otherwise specified
        PBS_NODEFILE = env.get('PBS_NODEFILE')
        if (PBS_NODEFILE is not None) and (ESCRIPT_HOSTFILE is None):
            ESCRIPT_HOSTFILE = PBS_NODEFILE
        if ESCRIPT_HOSTFILE is not None:
            if os.path.isfile(ESCRIPT_HOSTFILE):
                with open(ESCRIPT_HOSTFILE, 'r') as f:
                    hosts = [s for s in sorted(set(f.read().split('\n'))) if len(s) > 0] # sort unique no blanks
                with open(hostfile, 'w') as f:
                    f.write('\n'.join(hosts))
                HOSTLIST = ','.join(hosts)
                NUM_HOSTS = len(hosts)
                if ESCRIPT_NUM_NODES is not None:
                    if NUM_HOSTS < ESCRIPT_NUM_NODES:
                        die('Number of requested nodes must not exceed the number of entries selected in the host file '+
                            ESCRIPT_HOSTFILE+'.  You asked for '+str(ESCRIPT_NUM_NODES)+' from '+str(NUM_HOSTS)+'.')
                else:
                    ESCRIPT_NUM_NODES = NUM_HOSTS
            else:
                die('Cannot find hostfile '+ESCRIPT_HOSTFILE+'!')
        else:
            with open(hostfile, 'a') as f:
                f.write('\nlocalhost')
            HOSTLIST = 'localhost'

        vlog('ESCRIPT_NUM_NODES = '+str(ESCRIPT_NUM_NODES)+'\nESCRIPT_NUM_PROCS = '+str(ESCRIPT_NUM_PROCS))

    if WITH_OPENMP == '1':
        if ESCRIPT_NUM_THREADS is None:
            ESCRIPT_NUM_THREADS = int(env.get('OMP_NUM_THREADS', '1'))
        vlog('ESCRIPT_NUM_THREADS = '+str(ESCRIPT_NUM_THREADS))
    else:
        if ESCRIPT_NUM_THREADS is not None:
            if ESCRIPT_NUM_THREADS != 1:
                print('Warning: OpenMP is disabled but number of threads requested is '+
                    str(ESCRIPT_NUM_THREADS)+' != 1. Option ignored.')
        ESCRIPT_NUM_THREADS = 1

    #
    # Now we compute total number of Processes
    #
    try:
        TOTPROC = ESCRIPT_NUM_NODES * ESCRIPT_NUM_PROCS
    except:
        #Some compute error
        #This could happen if the args were not a number
        die('Expression of total number of processors = '+str(ESCRIPT_NUM_NODES)+' * '+
            str(ESCRIPT_NUM_PROCS)+' is not numerical!')

    # set up thread binding if unset -- disabled by default because it interfers
    # with MPI binding
    #if env.get('OMP_PROC_BIND') == '':
    #    #Force OpenMP binding for Intel (and GCC, though GCC is on by default)
    #    export_env['OMP_PROC_BIND'] = 'true'
    #if os.environ.get('KMP_AFFINITY') == '':
    #    #Set the style of binding (overrides OMP_PROC_BIND in many cases)
    #    export_env['KMP_AFFINITY'] = 'verbose,compact'

    #
    # Test to ensure people aren't trying to combine interactive and multi-process
    #
    if (DO_INTERACTIVE or (len(sys.argv) == 0)) and (TOTPROC > 1 ):
        die('Interactive mode cannot be used with more than one process!')

    if TOTPROC > 1:
        if ESCRIPT_CREATESTDFILES == 'y':
            PYTHON_MPI = PYTHON_MPI_REDIRECT
        else:
            PYTHON_MPI = PYTHON_MPI_NULL
    else:
        PYTHON_MPI = PYTHON_MPI_NULL
    #==============================================================================
    # Must have at least one command-line arg: the python script
    if len(sys.argv) == 0:
        if DO_BINARY:
            die('No program to run was specified!')
        else:
            DO_INTERACTIVE = True

    #==============================================================================

    if DO_XTERM:
        EXEC_CMD = 'xterm -e'
    else:
        EXEC_CMD = ''

    if DO_VALGRIND is not None:
        VALGRIND_BINS = which('valgrind')
        if len(VALGRIND_BINS) > 0:
            LOGDIR = os.path.join(ESCRIPT_ROOT, 'valgrind_logs')
            if not os.isdir(LOGDIR):
                os.mkdir(LOGDIR)
            if DO_VALGRIND.startwith('c'):
                # run callgrind
                LOGFILE = os.path.join(LOGDIR, 'callgrind.'+str(os.getpid())+'.xml')
                VALGRIND = 'valgrind --tool=callgrind --callgrind-out-file='+LOGFILE
                EXEC_CMD = (EXEC_CMD+' '+VALGRIND).strip()
            elif DO_VALGRIND.startwith('h'):
                # run cachegrind
                LOGFILE = os.path.join(LOGDIR, 'cachegrind.'+str(os.getpid())+'.xml')
                VALGRIND = 'valgrind --tool=cachegrind --cachegrind-out-file='+LOGFILE
                EXEC_CMD = (EXEC_CMD+' '+VALGRIND).strip()
            else:
                # run memcheck by default
                with open(LOGDIR, 'r') as f:
                    new_n = int(re.findall('^memcheck[^.]*[.](.*)$',f,re.MULTILINE)[-1])+1
                LOGFILE = os.path.join(LOGDIR, 'memcheck.'+new_n+'.'+str(os.getpid())+'.xml')
                VALGRIND = 'valgrind --tool=memcheck --xml=yes --show-reachable=yes --error-limit=no --suppressions=' \
                        +os.path.join(ESCRIPT_ROOT, 'scripts', 'escript.supp --leak-check=full --xml-file='+LOGFILE)
                EXEC_CMD = (EXEC_CMD+' '+VALGRIND).strip()
        else:
            die('Execution with valgrind requested but valgrind not in path!')

    if DO_BINARY:
        # TODO:
        EXEC_CMD = (EXEC_CMD+' '+PY_SCRIPT).strip()
    else:
        # Check to see if the python version we were compiled with matches the
        # one of PYTHON_CMD.
        compfull = get_buildvar('python_version')
        m = re.search('^(?P<major>[^.])[.](?P<minor>[^.])[.].*$', compfull)
        compmajor = m.group('major')
        compversion = compmajor+'.'+m.group('minor')
        if PYTHON_CMD == 'python': # if people have customised the command they
            if compmajor == '3':   # might not want us changing it
                PYTHON_CMD = 'python3'
        cmd = PYTHON_CMD+' -c "import sys;print(str(sys.version_info[0])+\'.\'+str(sys.version_info[1]))"' 
        intversion = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read().decode(
            sys.stdout.encoding).strip()
        if compversion != intversion:
            die("Python versions do not match. Escript was compiled for '"+
                compversion+"'.\nCurrent version of Python appears to be '"+intversion+"'.")
        if MPI_FLAVOUR == 'none':
            if DO_INTERACTIVE:
                EXEC_CMD = (EXEC_CMD+' '+PYTHON_CMD+' -i '+PY_SCRIPT).strip()
            else:
                EXEC_CMD = (EXEC_CMD+' '+PYTHON_CMD+' '+PY_SCRIPT).strip()
        else:
            if DO_INTERACTIVE:
                EXEC_CMD = (EXEC_CMD+' '+PYTHON_MPI+' -i '+PY_SCRIPT).strip()
            else:
                EXEC_CMD = (EXEC_CMD+' '+PYTHON_MPI+' '+PY_SCRIPT).strip()
    export_env['EXEC_CMD'] = EXEC_CMD
    vlog("Command to be executed is '"+EXEC_CMD+"'")

    #==============================================================================
    #
    #   now we start to spawn things:
    #
    if WITH_OPENMP == '1':
        export_env['OMP_NUM_THREADS'] = str(ESCRIPT_NUM_THREADS)

    if MPI_FLAVOUR == 'OPENMPI':
        if len(which('rsh'))+len(which('ssh')) == 0:
            AGENTOVERRIDE = '--gmca plm_rsh_agent /bin/false'

    prelaunch = fixEnvVars('')
    if len(prelaunch) > 0:
        vlog("Pre-launch command: '"+prelaunch+"'")
        subprocess.call(prelaunch, shell=True, env=export_env)

    launch = fixEnvVars('mpirun ${AGENTOVERRIDE} --gmca mpi_warn_on_fork 0 ${EE} --host ${HOSTLIST} --map-by node:pe=${ESCRIPT_NUM_THREADS} -bind-to core -np ${TOTPROC} ${EXEC_CMD}')
    if len(launch) > 0:
        vlog("Launch command: '"+launch+"'")
        exit_code = subprocess.call(launch, shell=True, env=export_env)

    postlaunch = fixEnvVars('')
    if len(postlaunch) > 0:
        vlog("Post-launch command: '"+postlaunch+"'")
        subprocess.call(postlaunch, shell=True, env=export_env)

    if DO_VALGRIND is not None:
        print('Valgrind log file written to '+LOGFILE)

    if os.path.isfile(hostfile):
        os.remove(hostfile)

    exit(exit_code)

if __name__=='__main__':
    main()
