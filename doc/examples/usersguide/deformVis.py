##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
##############################################################################

"""
VisIt Visualization Script for Deformation Data

This script demonstrates how to use VisIt's Python interface to create
visualizations of deformation data (stress and displacement) from VTK files.

Usage:
    visit -cli -nowin -s deformVis.py

The -nowin option prevents the visualization window from being shown, which
is not required since the purpose of the script is to save an image file.

Prerequisites:
    - VisIt must be installed
    - A deform.vtu file must exist in the current directory
"""

__copyright__ = """Copyright (c) 2003-2026 by the esys.escript Group
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""

# Open the VTU database
OpenDatabase("deform.vtu")

# Add pseudocolor plot for stress
AddPlot("Pseudocolor", "stress")
DrawPlots()

# Configure pseudocolor plot
p = PseudocolorAttributes()
p.centering = p.Nodal
SetPlotOptions(p)

# Add three-slice operator
AddOperator("ThreeSlice")
DrawPlots()

# Configure three-slice operator
t = ThreeSliceAttributes()
t.x = 0.3
t.y = 0.3
SetOperatorOptions(t)

# Add vector plot for displacement
AddPlot("Vector", "disp")
DrawPlots()

# Configure vector plot
v = VectorAttributes()
v.useStride = 1
SetPlotOptions(v)

# Configure save window attributes
s = SaveWindowAttributes()
# Change settings as required, for example:
# s.fileName = "deform_visualization"
# s.format = s.PNG
# s.width = 1024
# s.height = 768
SetSaveWindowAttributes(s)

# Save the window and exit
SaveWindow()
exit()
