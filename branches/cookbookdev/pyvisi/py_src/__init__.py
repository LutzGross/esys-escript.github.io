
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__deprecated__="""The pyvisi module has been deprecated and is no longer supported.
Please see the escript user guide and tutorials for suggestions
 on visualisation."""
"""
:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"

import warnings
warnings.warn("pyvisi is deprecated. Please see the escript user guide and tutorials for visualisation suggestions.", DeprecationWarning, stacklevel=2)

from camera import *
from carpet import *
from contour import *
from datacollector import *
from ellipsoid import *
from image import *
from light import *
from map import *
from position import *
from scene import *
from streamline import *
from text import *
from velocity import *
from imagereader import *
from logo import *
from legend import *
from movie import *
from rectangle import *
from rotation import *


