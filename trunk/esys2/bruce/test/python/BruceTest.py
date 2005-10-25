# $Id$

import os
import sys

from esys.escript import *
from esys.bruce import *

"""

Some simple tests of Bruce.

"""

b = Bruce()

assert (b.getDescription()=="Bruce")

assert (b.isValidFunctionSpaceType(0))
assert (b.isValidFunctionSpaceType(1))
assert ( not (b.isValidFunctionSpaceType(2)))

assert (b.getContinuousFunctionCode()==0)
assert (b.getFunctionCode()==1)

assert (b.getDim()==0)

brick = Brick(11,11,11,10,10,10)

assert (brick.getDim()==3)

brick.getX()
brick.getSize()

rectangle = Rectangle(11,11,10,10)

assert (rectangle.getDim()==2)

rectangle.getX()
rectangle.getSize()

sys.exit(0)
