

# Extensions to Scons

import py_compile
import sys
import os
import glob
import fnmatch
import types

from SCons.Script.SConscript import SConsEnvironment

###############################################################################
def matchingFiles(env,directory,includes='*',excludes=None) :

    # Pre-process for more convenient arguments

    if isinstance(includes,str) :
        includes = env.Split(includes)

    if isinstance(excludes,str) :
        excludes = env.Split(excludes)

    def fn_filter(node):
        fn = os.path.basename(str(node))
        match = False
        for include in includes:
            if fnmatch.fnmatchcase( fn, include ):
                match = True
                break

        if match and not excludes is None:
            for exclude in excludes:
                if fnmatch.fnmatchcase( fn, exclude ):
                    match = False
                    break

        return match

    def filter_nodes(where):
        contents = glob.glob( os.path.join(str(where),"*") )
        children = [ x for x in contents if fn_filter(x) ]
        nodes = []
        for f in children:
            nodes.append(gen_node(f))
        return nodes

    def gen_node(n):
        """Checks first to see if the node is a file or a dir, then
        creates the appropriate node. [code seems redundant, if the node
        is a node, then shouldn't it just be left as is?
        """
        if type(n) in (type(''), type(u'')):
            path = n
        else:
            path = n.abspath

        if os.path.isdir(path):
            return env.Dir(n)

        return env.File(n)

    here = env.Dir(directory)
    nodes = filter_nodes(here)

    node_srcs = [n.srcnode() for n in nodes]

    src = here.srcnode()
    if src is not here:
        for s in filter_nodes(src):
            if s not in node_srcs:
                # Probably need to check if this node is a directory
                nodes.append(
                    gen_node(os.path.join(str(directory),
                                          os.path.basename(str(s)))))

    return nodes
#==============================================================================
SConsEnvironment.matchingFiles = matchingFiles
###############################################################################



###############################################################################
def InstallPyModule(env,target_dir,source,shlib=None,excludes=None) :

    # put the .so/.dll over in  <target_dir>
    if shlib :
        shtarg = env.Install(target_dir, shlib)

    # Gather the python sources
    python_src = env.matchingFiles(source, "*.py",excludes=excludes)

    # Here is a hack to deal with the (possibly not yet generated)
    # __init__.py. The problem is that the python_src list is built before
    # __init__.py is updated from __init__.in. The AlwaysBuild call
    # ensures that this does not cause a problem, except on the first
    # build after a clean checkout, in which case there is no old
    # __init__.py in src, and hence it does not make it into python_src!

    init_input = env.matchingFiles(source, "__init__.in")
    if init_input :
        if python_src :
            names = [x.name for x in python_src]
            if "__init__.py" not in names :
                python_src.append(env.File(os.path.join(str(source),
                                                        "__init__.py")))
        else:
            python_src = [env.File(os.path.join(str(source),"__init__.py"))]
            
    # decide if we're doing py or pyc distribn and install.

    if  env['distrib_py_src'] :
        pytarg = env.Install(target_dir, python_src)
    else:
        pyc = env.PyCompile(python_src)
        pytarg = env.Install(target_dir, pyc)

    if shlib :
        targ_ret = env.Flatten([pytarg] + [shtarg])
    else:
        targ_ret = pytarg

    return targ_ret
#==============================================================================
SConsEnvironment.InstallPyModule = InstallPyModule
###############################################################################

###############################################################################
def installDirectory(env,target,source,includes="*", excludes=None,
                     recursive=False):
    
    
    if os.path.isfile(str(target)) :
        raise UserError("target must be a directory")

    if os.path.isfile(str(source)) :
        raise UserError("source must be a directory")

    source_files = env.matchingFiles(source,includes,excludes)

    ret_targ = []

    for f in source_files :
        if f.isfile() :
            targ = env.Install(target,f)
            ret_targ.append(targ)

        if f.isdir() and recursive :
            x = os.path.basename(str(f))
            t = env.Dir(os.path.join(str(target),x))
            targ = env.installDirectory(t,f,includes,excludes,recursive)
            ret_targ += targ
    
    return ret_targ
#==============================================================================
SConsEnvironment.installDirectory = installDirectory
###############################################################################

###############################################################################
# Code to build .pyc from .py
def build_py(target, source, env):
    py_compile.compile(str(source[0]), str(target[0]))
    return 0

def doSubstitution(target,source,env) :
    import product_info as PI
    data = source[0].get_contents()
    data = data.replace('$ProductName$',PI.PRODUCT_NAME)
    data = data.replace('$LowerProductName$',PI.product_name)
    data = data.replace('$ProductVersion$',PI.product_version)
    data = data.replace('$VersionString$',PI.pkg_version_string)
    data = data.replace('$SVNRevision$',PI.svn_revision)
    open(str(target[0]),'w').write(data)
    return 0

# Code to run unit_test executables
def runUnitTest(target, source, env):
    app = str(source[0].abspath)

    olddir = os.getcwd()
    newdir = os.path.dirname(str(source[0]))
    os.chdir(newdir)

    if env.Execute(app) != 0:
        os.chdir(olddir)
        return 1

    os.chdir(olddir)
    open(str(target[0]),'w').write("PASSED\n")
    return 0

def runPyUnitTest(target, source, env): 
    app = env['python_path'] + ' "' + str(source[0].abspath) + '"'

    olddir = os.getcwd()
    newdir = os.path.dirname(str(source[0]))
    os.chdir(newdir)

    if env.Execute(app)  != 0:
        os.chdir(olddir)
        return 1

    os.chdir(olddir)
    open(str(target[0]),'w').write("PASSED\n")
    return 0
  

def addBuilders(env) :
    py_builder = env.Builder(action = build_py,
                             suffix = '.pyc',
                             src_suffix = '.py',
                             single_source=True)

    env.Append(BUILDERS = {'PyCompile' : py_builder});

    substituter = env.Builder(action = doSubstitution,
                              suffix = '',
                              src_suffix = '.in',
                              single_source=True )

    env.Append(BUILDERS = {'VariableSubstitution' : substituter});

    runUnitTest_builder = env.Builder(action = runUnitTest,
                                      suffix = '.passed',
                                      src_suffix=env.get('PROGSUFFIX',''),
                                      single_source=True)

    env.Append(BUILDERS = {'RunUnitTest' : runUnitTest_builder});

    runPyUnitTest_builder = env.Builder(action = runPyUnitTest,
                                        suffix = '.passed',
                                        src_suffix='.py',
                                        single_source=True)

    env.Append(BUILDERS = {'RunPyUnitTest' : runPyUnitTest_builder});
    return
#==============================================================================
SConsEnvironment.addBuilders = addBuilders
###############################################################################

###############################################################################
def epydocAction(target, source, env):

    doc_dir = os.path.dirname(str(target[0].abspath))

    cmd = [
        [env['epydoc_path'], "-qqqq", "-o", doc_dir, str(source[0].tpath)]
          ]
    print 'executing epydoc...'
    return env.Execute(cmd,"Build epydoc documentation")


def epydocDepend(env, target, source,
                 src_pattern="*.py",
                 file_names=['index.html','epydoc.css'],
                 subdirs=['private','public']) :
    """
    \brief add a dependency between the directory containing the doco and that
           containing the source
    \param target - the directory containing index.html for the doco,
                    and the private and public sub-directories.
    \param source - the python module source directory
    """
    the_subdirs = [os.path.join(target.abspath,sub) for sub in subdirs]
    the_targets = [os.path.join(target.abspath,file) for file in file_names]


    ret_target = target
    
    # if absolutely anything goes wrong, turn this on.
    force_build = False

    dst_time = 0
    src_time = 1

    try:

        # have a shot at digging out all the source file times.
        src_time = max(
            [os.path.getmtime(str(x))
             for x in env.matchingFiles(source,src_pattern) +
                      [source.abspath]]
            )

        # now try to find the target files and their mod time.
        a = [os.path.getmtime(os.path.join(target.abspath,str(x)))
             for x in file_names ] 

        for directory in the_subdirs :
            # include the mod time of the directory
            a += [os.path.getmtime(str(directory))]

            # now go for the mod times of all files below the subdirs.
            if os.path.isdir(str(directory)) :
                a += [os.path.getmtime(str(x))
                      for x in env.matchingFiles(directory,"*")]

            else:
                # if it is not a directory, and we expected it to be
                # do something nasty.
                # we're in a try anyway.
                force_build = True
                os.unlink(directory)

        dst_time = max(a)
            
    except:
        # Force an unlink and re-build.
        force_build = True
        
    if src_time > dst_time or force_build :
        for x in the_targets :
            try:
                os.unlink(x)
            except OSError:
                pass

        ret_target = env.Command(the_targets,source,epydocAction)
    

    env.Clean(target, the_subdirs + file_names)

    return ret_target
#==============================================================================
SConsEnvironment.epydocDepend = epydocDepend
###############################################################################


###############################################################################
def genSConscriptCalls(env,dir_list):

    for d in dir_list :
        print 'calling SConscript in "./%s"' %(d)
        env.SConscript(dirs = [d],
                       build_dir='build/$PLATFORM/' + str(d),
                       duplicate=0)


    return
#==============================================================================
SConsEnvironment.genSConscriptCalls = genSConscriptCalls
###############################################################################


###############################################################################
def print_all_nodes(env,dirnode, level=0):
    """Print all the scons nodes that are children of this node, recursively."""
    if type(dirnode)==type(''):
        dirnode=env.Dir(dirnode)
    dt = type(env.Dir('.'))
    for f in dirnode.all_children():
        if type(f) == dt:
            print "%s%s: .............."%(level * ' ', str(f))
            env.print_all_nodes(f, level+2)
        print "%s%s"%(level * ' ', str(f))

#==============================================================================
SConsEnvironment.print_all_nodes = print_all_nodes
###############################################################################
  




















###############################################################################
def Glob(env,dir='.',includes="*",excludes=None, scan_dir=True):

    """Similar to glob.glob, except globs SCons nodes, and thus sees
    generated files and files from build directories.  Basically, it sees
    anything SCons knows about.  A key subtlety is that since this function
    operates on generated nodes as well as source nodes on the filesystem,
    it needs to be called after builders that generate files you want to
    include.

    It will return both Dir entries and File entries
    """

    # Pre-process for more convenient arguments

    if isinstance(includes,str) :
        includes = env.Split(includes)

    if isinstance(excludes,str) :
        excludes = env.Split(excludes)

    def fn_filter(node):
        fn = os.path.basename(str(node))
        match = 0
        for include in includes:
            if fnmatch.fnmatchcase( fn, include ):
                match = 1
                break

        if match == 1 and not excludes is None:
            for exclude in excludes:
                if fnmatch.fnmatchcase( fn, exclude ):
                    match = 0
                    break

        return match

    def filter_nodes(where):
        contents = where.all_children(scan=scan_dir)
        children = [ x for x in contents if fn_filter(x) ]
        nodes = []
        for f in children:
            nodes.append(gen_node(f))
        return nodes

    def gen_node(n):
        """Checks first to see if the node is a file or a dir, then
        creates the appropriate node. [code seems redundant, if the node
        is a node, then shouldn't it just be left as is?
        """
        if type(n) in (type(''), type(u'')):
            path = n
        else:
            path = n.abspath
        if os.path.isdir(path):
            return env.Dir(n)
        else:
            return env.File(n)

    here = env.Dir(dir)
    nodes = filter_nodes(here)

    node_srcs = [n.srcnode() for n in nodes]

    src = here.srcnode()
    if src is not here:
        for s in filter_nodes(src):
            if s not in node_srcs:
                # Probably need to check if this node is a directory
                nodes.append(
                    gen_node(os.path.join(dir,os.path.basename(str(s)))))

    return nodes
