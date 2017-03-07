
##############################################################################
#
# Copyright (c)2015-2017 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c)2015-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"


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
  deps=deps.split("\n")
  for line in deps:
    ind=line.find("=")
    if ind==-1:
        continue
    key=line[:ind]
    val=line[ind+1:]
    if key in ignoreflags:
        continue
    if key not in usedflags:
        raise RuntimeError("Unknown key ("+key+") in dpkg-buildflags")
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
