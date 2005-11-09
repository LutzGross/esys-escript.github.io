import sys
import os
import unittest

from esys.escript import *
from esys.escript.linearPDEs import *

from esys import bruce

mydomain=bruce.Rectangle(5,5)

u1=Scalar(1.0,Function(mydomain))
u1.setTaggedValue(1,2.0)

u2=Scalar(1.0,Function(mydomain))
u2.setTaggedValue(1,3.0)


print "u1:"
print u1

print "u2:"
print u2

v=u1+u2
print "v=u1+u2:"
print v

print "Lsup(v):"
print Lsup(v)
