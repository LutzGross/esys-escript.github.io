
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################
import py_compile
from sys import * 

if len(argv)!=3:
    print("%s source dest"%argv[0], file=stderr)
    exit(2)
try:
    py_compile.compile(argv[1], argv[2], argv[1], True)
except Exception as e:
   print(e.args)
   exit(1)

