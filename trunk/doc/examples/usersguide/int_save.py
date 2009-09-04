
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



# This example demonstrates both interpolation and saving data in CSV format

from esys.escript import saveDataCSV, sup
import numpy
from esys.finley import Rectangle

n=4		#Change this to whatever you like
r=Rectangle(n,n)
x=r.getX()
x0=x[0]
x1=x[1]    #we'll use this later
toobig=100	#An exception will be thrown if interpolation produces a value larger than this

#First one dimensional interpolation

#In this example we will interpolate a sine curve
#The values we take from the domain will range from 0 to 1 (inclusive)
#Because we will actually see the maximum value in the input we need to add
#an extra entry to the table (so we have 1.125 sine cycles)

sine_table=[0, 0.70710678118654746, 1, 0.70710678118654746, 0, -0.70710678118654746, -1, -0.70710678118654746, 0, 0.70710678118654746]

numslices=len(sine_table)-1

minval=0
maxval=1.125	# The extra 0.25 is to account for the extra entry in the table

step=sup(maxval-minval)/numslices	#The width of the gap between entries in the table

result=x0.interpolateTable(sine_table, minval, step, toobig)

#Now we save the input and output for comparison

saveDataCSV("1d.csv", inp=x0, out=result)

#Now 2D interpolation

#This time the sine curve will be at full height along the x (ie x0) axis.
#Its amplitude will decrease to a flat line along x1=1.1
#Since 1.1 is larger than any value in the input we don't need to add an extra row
#to our table

#Interpolate works with numpy arrays so we'll use them
#st=numpy.array(sine_table)

#table=[st, 0.5*st, 0*st ]   #Note that this table is 2D

##note that we call the interpolate table method on the object
##which corresponds to the outer dimension of the table
#result2=x1.interpolateTable(table, 0, 0.55, x0, minval, step, 500)
#saveDataCSV("2d.csv",inp0=x0, inp2=x1, out=result2)