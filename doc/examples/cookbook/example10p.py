##############################################################################
#
# Copyright (c) 2009-2018 by The University of Queensland
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
from __future__ import division, print_function

__copyright__="""Copyright (c) 2009-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""

############################################################FILE HEADER
# example10p.py
# Profiles of gravitational potential with increasing domain size.

#######################################################EXTERNAL MODULES
# To solve the problem it is necessary to import the modules we require.
import pylab as pl
import os
from esys.escript import * # This imports everything from the escript library

################################################ESTABLISHING PARAMETERS
#the folder to put our outputs in, leave blank "" for script path 
save_path= os.path.join("data","example10")
#ensure the dir exists
mkDir(save_path)

#load the data files.
dat1=pl.loadtxt(os.path.join(save_path,'example10b_2000.asc'))
dat2=pl.loadtxt(os.path.join(save_path,'example10b_0500.asc'))
dat3=pl.loadtxt(os.path.join(save_path,'example10b_4000.asc'))
dat4=pl.loadtxt(os.path.join(save_path,'example10b_1000.asc'))

#plot the data.
pl.plot(dat2,label='500^2')
pl.plot(dat4,label='1000^2')
pl.plot(dat1,label='2000^2')
pl.plot(dat3,label='4000^2')
pl.legend()
pl.title('Small vs Large Gravity Model')
pl.xlabel('X')
pl.ylabel('Potential')
pl.savefig(os.path.join(save_path,"ex10p_boundeff.pdf"),dpi=150)


