
##############################################################################
#
# Copyright (c) 2003-2016 by The University of Queensland
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

# This is a template configuration file for escript on Debian/GNU Linux.
# Refer to README_FIRST for usage instructions.

from scons.templates.sid_options import *

pythoncmd = 'python3'

import subprocess
p=subprocess.Popen([pythoncmd, '-c', 'import sysconfig\nprint(sysconfig.get_config_var("LDLIBRARY"))'], stdout=subprocess.PIPE)
pythonlibname = p.stdout.readline().encode().strip()
p.wait()

boost_libs = [p3name[3:-3]]
