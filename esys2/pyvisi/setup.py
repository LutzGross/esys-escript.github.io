#!/usr/bin/env python

# $Id$

from distutils.core import setup

a=setup(name="pyvisi",
      version="0.1a",
      description="The Python Visualisation Interface",
      author="Paul Cochrane",
      author_email="cochrane@esscc.uq.edu.au",
      url="http://pyvisi.sourceforge.net",
      packages=['pyvisi','pyvisi.lib'],
      scripts=['bin/pyvisi'],
)


