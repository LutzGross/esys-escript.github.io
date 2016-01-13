
##############################################################################
#
# Copyright (c) 2003-2015 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from .sid_options import *

usepython3=True
pythoncmd='python3'

import subprocess
import sysconfig

#pythonlibname='python3.4m'
#pythonlibname=sysconfig.get_config_var("LDLIBRARY")

p=subprocess.Popen([pythoncmd, '-c', 'import sysconfig\nprint(sysconfig.get_config_var("LDLIBRARY"))'], stdout=subprocess.PIPE)
pythonlibname=p.stdout.readline().encode().strip()
p.wait()

#pythonincpath='/usr/include/python3.4'
#pythonincpath=sysconfig.get_config_var("INCLUDEPY")

p=subprocess.Popen([pythoncmd, '-c', 'import sysconfig\nprint(sysconfig.get_config_var("INCLUDEPY"))'], stdout=subprocess.PIPE)
pythonincpath=p.stdout.readline().encode().strip()
p.wait()


import os

p=subprocess.Popen(["ld","--verbose"], stdout=subprocess.PIPE)
out,err=p.communicate()
spath=[x[13:-3] for x in out.split() if 'SEARCH_DIR' in x]
p2name=''
p3name=''
for name in spath:
  try:
    l=os.listdir(name)
    p2res=[x for x in l if x.startswith('libboost_python-py2') and x.endswith('.so')]
    p3res=[x for x in l if x.startswith('libboost_python-py3') and x.endswith('.so')]
    if len(p2name)==0 and len(p2res)>0:
      p2name=p2res[0]
    if len(p3name)==0 and len(p3res)>0:
      p3name=p3res[0]
  except OSError:
    pass

# boost-python library/libraries to link against
boost_libs = [p3name[3:-3]]
