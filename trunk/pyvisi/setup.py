#!/usr/bin/env python

# $Id: setup.py,v 1.6 2005/04/29 00:19:11 paultcochrane Exp $

from distutils.core import setup

a=setup(name="pyvisi",
      version="0.1-pre-alpha-4",
      description="The Python Visualisation Interface",
      author="Paul Cochrane",
      author_email="cochrane@esscc.uq.edu.au",
      url="http://pyvisi.sourceforge.net",
      packages=['pyvisi',
      'pyvisi.renderers',
      'pyvisi.renderers.gnuplot',
      'pyvisi.renderers.vtk',
      ],
)


