# Extensions to Scons
def build_py(target, source, env):
   # Code to build .pyc from .py
   import py_compile, sys;
   py_compile.compile(str(source[0]), str(target[0]))
   return None
   
def runUnitTest(target, source, env): 
   import os
   app = str(source[0].abspath)
   if not os.system(app):
      open(str(target[0]),'w').write("PASSED\n")
   return None

def runPyUnitTest(target, source, env): 
   import os
   app = 'python '+str(source[0].abspath)
   if not os.system(app):
      open(str(target[0]),'w').write("PASSED\n")
   return None


#import fnmatch
#import os

#def Glob( includes = Split( '*' ), excludes = None, dir = '.'):
#   """Similar to glob.glob, except globs SCons nodes, and thus sees
#   generated files and files from build directories.  Basically, it sees
#   anything SCons knows about.  A key subtlety is that since this function
#   operates on generated nodes as well as source nodes on the filesystem,
#   it needs to be called after builders that generate files you want to
#   include.
#
#   It will return both Dir entries and File entries
#   """
#   def fn_filter(node):
#      fn = os.path.basename(str(node))
#      match = 0
#      for include in includes:
#         if fnmatch.fnmatchcase( fn, include ):
#            match = 1
#            break
#
#      if match == 1 and not excludes is None:
#         for exclude in excludes:
#            if fnmatch.fnmatchcase( fn, exclude ):
#               match = 0
#               break
#
#      return match
#
#   def filter_nodes(where):
#       children = filter(fn_filter, where.all_children())
#       nodes = []
#       for f in children:
#           nodes.append(gen_node(f))
#       return nodes
#
#   def gen_node(n):
#       """Checks first to see if the node is a file or a dir, then
#       creates the appropriate node. [code seems redundant, if the node
#       is a node, then shouldn't it just be left as is?
#       """
#       if type(n) in (type(''), type(u'')):
#           path = n
#       else:
#           path = n.abspath
#       if os.path.isdir(path):
#           return Dir(n)
#       else:
#           return File(n)
#
#   here = Dir(dir)
#   nodes = filter_nodes(here)
#
#   node_srcs = [n.srcnode() for n in nodes]
#
#   src = here.srcnode()
#   if src is not here:
#       for s in filter_nodes(src):
#           if s not in node_srcs:
#               # Probably need to check if this node is a directory
#               nodes.append(gen_node(os.path.join(dir,os.path.basename(str(s)))))
#
#   return nodes