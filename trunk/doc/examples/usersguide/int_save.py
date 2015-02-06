from __future__ import division
##############################################################################
#
# Copyright (c) 2003-2015 by University of Queensland
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

__copyright__="""Copyright (c) 2003-2015 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"



# This example demonstrates both interpolation and saving data in CSV format

from esys.escript import saveDataCSV, sup, interpolateTable
import numpy
try:
    from esys.finley import Rectangle
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False

if not HAVE_FINLEY:
    print("Finley module not available")
else:

    n=4		#Change this to whatever you like
    r=Rectangle(n,n)
    x=r.getX()
    x0=x[0]
    x1=x[1]    #we'll use this later
    toobig=100	#An exception will be thrown if interpolation produces a value larger than this

    #First one dimensional interpolation

    #In this example we will interpolate a sine curve
    #The values we take from the domain will range from 0 to 1 (inclusive)

    sine_table=[0, 0.70710678118654746, 1, 0.70710678118654746, 0, -0.70710678118654746, -1, -0.70710678118654746, 0]

    numslices=len(sine_table)-1

    minval=0
    maxval=1

    step=sup(maxval-minval)/numslices	#The width of the gap between entries in the table

    result=interpolateTable(sine_table, x0, minval, step, toobig)

    #Now we save the input and output for comparison

    saveDataCSV("1d.csv", inp=x0, out=result)

    #Now 2D interpolation

    #This time the sine curve will be at full height along the x (ie x0) axis.
    #Its amplitude will decrease to a flat line along x1=1.1

    #Interpolate works with numpy arrays so we'll use them
    st=numpy.array(sine_table)

    table=[st, 0.5*st, 0*st ]   #Note that this table is 2D

    #The y dimension should be the outer the dimension of the table
    #Note that the order of tuples for the 2nd and 3rd param is (x,y)
    result2=interpolateTable(table, x, (minval,0), (0.55, step), toobig)
    saveDataCSV("2d.csv",inp0=x0, inp2=x1, out=result2)
