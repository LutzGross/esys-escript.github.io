#!/usr/bin/python

##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the MIT (ExPat) license
# https://opensource.org/licenses/MIT
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

# locates the source of .pyc in the esys directory and copies to the specified dest directory


import os, shutil, sys

if len(sys.argv)!=2:
   print("Please specify source directory", file=sys.stderr)
   exit(1)

print("TESTING")

for dirn, subdir, files in os.walk("esys"):
  if dirn.find('__pycache__')!=-1:
    continue
  first=True
  for n in files:
    if n.endswith(".pyc"):
      n=n[:-1]
      if first:
        first=False
        print("os.makedirs("+dirn+")")
        try:
          os.makedirs(dirn)
        except OSError:
          pass
      lst=dirn.split("/")
      if len(lst)==1:
        continue
      source="/".join([sys.argv[1]]+[lst[1]]+["py_src"]+lst[2:]+[n])
      dest=dirn+"/"+n
      shutil.copyfile(source,dest)
      print("Copy: "+source+"  "+dest)
