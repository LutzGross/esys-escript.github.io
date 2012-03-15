from __future__ import print_function
import py_compile
from sys import * 

if len(argv)!=3:
    print("%s source dest"%argv[0], file=stderr)
    exit(2)
try:
    py_compile.compile(argv[1], argv[2], argv[1], True)
except Exception as e:
   print(e.args)

