
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


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"


#Extracts Debian build options

import subprocess
import sys

def getdebbuildflags():
  usedflags={'CFLAGS':None, 'CPPFLAGS':'cpp_flags', 'CXXFLAGS':'cxx_extra', 'LDFLAGS':'ld_extra'}
  ignoreflags=['FFLAGS','FCFLAGS', 'GCJFLAGS','OBJCFLAGS','OBJCXXFLAGS']
  mycflags=None
  mycxxflags=None
  try:
    deps=subprocess.check_output("dpkg-buildflags")
  except OSError:
    return []
  res=[]
  if not isinstance(deps, str):
      deps=deps.decode()
  deps=deps.split("\n")
  for line in deps:
    ind=line.find("=")
    if ind==-1:
        continue
    key=line[:ind]
    val=line[ind+1:]
    print(key, val)
    if key in usedflags:
      if key=="CFLAGS":
        mycflags=val
      if key=="CXXFLAGS":
        mycxxflags=val
      if mycflags is not None and mycxxflags is not None and mycflags!=mycxxflags:
        raise RuntimeError("We do not current support different different dpkg-buildflags for C vs C++")
      if usedflags[key] is None:
        continue
      res.append([usedflags[key],val])
  return res    
