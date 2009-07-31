# This is the cookbook library for generic scripts.
# Imports come from the folder cblib and contain script definitions.

import sys
import os

#set the subdirectory path
sys.path.insert(0,'cblib')
#import all examples to library
from wavesolver2d import *
from wavesolver2df import *
from phones import *
#exit subdirectory
del sys.path[0]