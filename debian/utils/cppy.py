#!/usr/bin/python

# locates the source of .pyc in the esys directory and copies to the specified dest directory

from __future__ import print_function, division

import os, shutil, sys

if len(sys.argv)!=2:
   print("Please specify source directory", file=sys.stderr)
   exit(1)

print("TESTING")

for dirn, subdir, files in os.walk("esys"):
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
