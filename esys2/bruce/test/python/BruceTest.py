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

brick = Brick()

assert (brick.getDim()==3)

brick.getX()

rectangle = Rectangle()

assert (rectangle.getDim()==2)

rectangle.getX()

sys.exit(0)
