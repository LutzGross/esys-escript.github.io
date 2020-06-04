##############################################################################
#
# Copyright (c) 2015-2018 by The University of Queensland
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

from __future__ import print_function, division
"""
Test script to run test model COMMEMI-4
"""

try:
  import matplotlib
  # The following line is here to allow automated testing. Remove or comment if
  # you would like to display the final plot in a window instead.
  matplotlib.use('agg')
  
  from matplotlib import pyplot
  HAVE_MPL=True
except:
  HAVE_MPL=False

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

import numpy
import datetime
import esys.downunder.magtel2d as mt2d
import esys.escript            as escript
import esys.escript.pdetools   as pdetools
try:
    import esys.finley         as finley
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False
    
try: 
    from scipy.interpolate import InterpolatedUnivariateSpline
    HAVE_SCIPY = True
except:
    HAVE_SCIPY = False

HAVE_GMSH = escript.hasFeature("gmsh")
HAVE_DIRECT = escript.hasFeature("PASO_DIRECT") or escript.hasFeature('trilinos')

# Matplotlib uses outdated code -- ignore the warnings until an update is available:
import warnings
warnings.filterwarnings("ignore")  #, category=DeprecationWarning
import logging
# this is mainly to avoid warning messages
logging.basicConfig(format='%(name)s: %(message)s', level=logging.INFO)


# ==========================================================
# ==========================================================



def setupMesh(mode, coord, elem_sizes):         
    #---------------------------------------------------------------------------
    # DESCRIPTION:
    # -----------
    # This is a utility function which setups the COMMEMI-4 mesh.
    # 
    #
    # ARGUMENTS:                                                              
    # ----------
    # mode       :: TE or TM mode.
    # coord      :: dictionary with coordinate tuples.
    # elem_sizes :: mesh element sizes, large, normal, small. 
    #
    # RETURNS:
    # --------
    # <Nothing> A mesh file is written to the output folder.
    # 
    #
    # AUTHOR:
    # -------
    # Ralf Schaa, 
    # University of Queensland
    #
    #---------------------------------------------------------------------------



    #---------------------------------------------------------------------------
    # Imports.
    #---------------------------------------------------------------------------
        
    import esys.pycad              as pycad     # @UnresolvedImport   
    import esys.finley             as finley    # @UnresolvedImport
    import esys.escript            as escript   # @UnresolvedImport
    import esys.weipa              as weipa     # @UnresolvedImport    
    # <Note>: "@UnresolvedImport" ignores any warnings in Eclipse/PyDev (PyDev has trouble with external libraries).



    model = "COMMEMI-4"

    print("Preparing the mesh " + model + " ...")
    print("")
    
    # Warn about magnetotelluric TM mode:
    if mode.lower() == 'tm':
        print("TM mode not yet supported")
        return


        
    # Path to write the mesh:
    outpath = "../out/commemi4"
    
    
        
     
    # --------------------------------------------------------------------------
    # Initialisations.
    # --------------------------------------------------------------------------

    # Get coordinates from dictionary as list of tuples  
    a0 = coord["air"]   
    l1 = coord["lyr1"]  
    s1 = coord["slab"]  
    b1 = coord["basin"] 
    l2 = coord["lyr2"]  
    l3 = coord["lyr3"]  
    
    # Mesh length from top-boundary.
    x_extent = abs(a0[3][0]-a0[0][0])
    
    

        
    # --------------------------------------------------------------------------
    # Point definitions.
    # --------------------------------------------------------------------------
    
    #<Note>: define all points spanning the mesh, anomalies and layers; 
    #        note also shared domain points must be defined only once.
 
 
    # Mesh top boundary.    
    air = []
    air.append( pycad.Point( *a0[0] ) )    # 0: left  , top    (@ boundary)
    air.append( pycad.Point( *a0[3] ) )    # 3: right , top    (@ boundary)
    
    
    # First-layer.
    ly1 = []
    ly1.append( pycad.Point( *l1[0] ) )    # 0: left  , top    (@ air/earth interface)                       
    ly1.append( pycad.Point( *l1[1] ) )    # 1: left  , bottom (@ boundary)                       
    ly1.append( pycad.Point( *l1[2] ) )    # 2: right , bottom (@ slab/basin)   
    ly1.append( pycad.Point( *l1[3] ) )    # 3: right , bottom (@ boundary)     
    ly1.append( pycad.Point( *l1[4] ) )    # 4: right , top    (@ air/earth interface)                 

   
    # Slab.
    sl1 = []
    sl1.append( ly1[1]                )    # 0: left  , top    (@ boundary)                       
    sl1.append( pycad.Point( *s1[1] ) )    # 1: left  , bottom (@ boundary)                       
    sl1.append( pycad.Point( *s1[2] ) )    # 2: right , bottom (@ slab/basin)                         
    sl1.append( ly1[2]                )    # 3: right , top    (@ slab/basin)                       
    
    
    # Basin.
    bs1 = []
    bs1.append( ly1[2]                )    # 0: left  , top    (@ slab/basin)
    bs1.append( sl1[2]                )    # 1: left  , centre (@ slab/basin) 
    bs1.append( pycad.Point( *b1[2] ) )    # 2: left  , bottom (@ lyr1/basin)                       
    bs1.append( pycad.Point( *b1[3] ) )    # 3: centre, bottom (@ lyr1/basin)                       
    bs1.append( pycad.Point( *b1[4] ) )    # 4: edge  , bottom (@ lyr1/basin)                       
    bs1.append( pycad.Point( *b1[5] ) )    # 5: right , bottom (@ boundary)
    bs1.append( ly1[3]                )    # 6: right , top 
    
    
    # Second-Layer.
    ly2 = []
    ly2.append( sl1[1]                )    # 0: left  , top    (@ lyr2/slab)
    ly2.append( pycad.Point( *l2[1] ) )    # 1: left  , bottom (@ boundary) 
    ly2.append( pycad.Point( *l2[2] ) )    # 2: right , bottom (@ boundary)                       
    ly2.append( bs1[5]                )    # 3: right , top    (@ basin/boundary)                       
    ly2.append( bs1[4]                )    # 4: edge  , top    (@ lyr2/basin)                      
    ly2.append( bs1[3]                )    # 5: centre, top    (@ lyr2/basin)
    ly2.append( bs1[2]                )    # 6: left  , top    (@ lyr2/basin)
    ly2.append( sl1[2]                )    # 7: left  , centre (@ slab/basin) 
    
    
    # Basement layer.       
    ly3 = []    
    ly3.append( ly2[1]                )    # 0: left  , top    (@ boundary)
    ly3.append( pycad.Point( *l3[1] ) )    # 1: left  , bottom (@ boundary) 
    ly3.append( pycad.Point( *l3[2] ) )    # 2: right , bottom (@ boundary) 
    ly3.append( ly2[2]                )    # 3: right , top    (@ boundary)
    #___________________________________________________________________________




    # --------------------------------------------------------------------------
    # Line definitions.
    # --------------------------------------------------------------------------

    #<Note>: connects the points to define lines counterclockwise;    
    #        shared lines are re-used to ensure that all domains  
    #        are recognised as parts of the same mesh. 
        
    # Air.
    ln0 = []
    ln0.append( pycad.Line(air[0], ly1[0]) ) # 0 left-top     to left-bottom.
    ln0.append( pycad.Line(ly1[0], ly1[4]) ) # 1 left-bottom  to right-bottom (air-earth interface).
    ln0.append( pycad.Line(ly1[4], air[1]) ) # 2 right-bottom to right-top.
    ln0.append( pycad.Line(air[1], air[0]) ) # 3 right-top    to left-top.
        
    # Top Layer.
    ln1 = []
    ln1.append( pycad.Line(ly1[0], ly1[1]) ) # 0 left-top         to left-bottom.   
    ln1.append( pycad.Line(ly1[1], ly1[2]) ) # 1 left-bottom      to start-slab/basin.  
    ln1.append( pycad.Line(ly1[2], ly1[3]) ) # 2 start-slab/basin to basin-boundary 
    ln1.append( pycad.Line(ly1[3], ly1[4]) ) # 3 basin-boundary   to right-top.     
    ln1.append( -ln0[1]                    ) # 4 right-top        to left-top.

 
    # Slab.
    ln2 = []
    ln2.append( pycad.Line(sl1[0], sl1[1]) ) # 0 left-top     to left-bottom.   
    ln2.append( pycad.Line(sl1[1], sl1[2]) ) # 1 left-bottom  to right-bottom.         
    ln2.append( pycad.Line(sl1[2], sl1[3]) ) # 2 right-bottom to right-top.            
    ln2.append( -ln1[1]                    ) # 3 right-top    to left-top


    # Basin.
    ln3 = []
    ln3.append( -ln2[2]                    ) # 0 left-top         to left-centre.         
    ln3.append( pycad.Line(bs1[1], bs1[2]) ) # 1 left-centre      to left-bottom.         
    ln3.append( pycad.Line(bs1[2], bs1[3]) ) # 2 left-bottom      to mid-bottom.          
    ln3.append( pycad.Line(bs1[3], bs1[4]) ) # 3 mid-bottom       to right-mid-top.       
    ln3.append( pycad.Line(bs1[4], bs1[5]) ) # 4 right-mid-top    to right-bottom.        
    ln3.append( pycad.Line(bs1[5], bs1[6]) ) # 5 right-bottom     to right-top.           
    ln3.append( -ln1[2]                    ) # 6 right-top        to right-slab/basin.    
    
    
    # Layer below.
    ln4 = []
    ln4.append( pycad.Line(ly2[0], ly2[1]) ) # 0 left-top      to left-bottom.        
    ln4.append( pycad.Line(ly2[1], ly2[2]) ) # 1 left-bottom   to right-bottom.        
    ln4.append( pycad.Line(ly2[2], ly2[3]) ) # 2 right-bottom  to right-top.            
    ln4.append( -ln3[4]                    ) # 3 right-top     to right-mid-top.       
    ln4.append( -ln3[3]                    ) # 4 right-mid-top to mid-bottom.          
    ln4.append( -ln3[2]                    ) # 5 mid-bottom    to left-bottom.         
    ln4.append( -ln3[1]                    ) # 6 left-bottom   to left-centre.         
    ln4.append( -ln2[1]                    ) # 7 left-centre   to left-top.            
        
    # Basement layer.
    ln5 = []
    ln5.append( pycad.Line(ly3[0], ly3[1]) ) # 0 left-top     to left-bottom.
    ln5.append( pycad.Line(ly3[1], ly3[2]) ) # 1 left-bottom  to right-bottom.
    ln5.append( pycad.Line(ly3[2], ly3[3]) ) # 2 right-bottom to right-top.
    ln5.append( -ln4[1]                    ) # 3 right-top    to left-top.
    #___________________________________________________________________________




    # --------------------------------------------------------------------------
    # Domain definitions.
    # --------------------------------------------------------------------------
    
       
    # First define all borders.       
    borders = []   
    borders.append( pycad.CurveLoop(*ln0) )   
    borders.append( pycad.CurveLoop(*ln1) )   
    borders.append( pycad.CurveLoop(*ln2) )   
    borders.append( pycad.CurveLoop(*ln3) )    
    borders.append( pycad.CurveLoop(*ln4) )    
    borders.append( pycad.CurveLoop(*ln5) )    

    # And next the domains.
    domains = []
    for i in range( len(borders) ):        
        domains.append( pycad.PlaneSurface(borders[i]) ) 
    #___________________________________________________________________________




    # --------------------------------------------------------------------------
    # Set element sizes in domains.
    # --------------------------------------------------------------------------
    
    # Horizontal extents of segments along slab and basin:
    x_extents = []
    x_extents.append( l1[2][0] - l1[0][0] ) # 0
    x_extents.append( l1[3][0] - l1[2][0] ) # 1

    # Number of elements in the air-domain, first-layer as well as slab- and basin-domain.
    domains[0].setElementDistribution(     x_extent / elem_sizes["large"]   )
    domains[1].setElementDistribution(     x_extent / (elem_sizes["small"]) )
    domains[2].setElementDistribution( 0.4*x_extent / (elem_sizes["small"]) )
    domains[3].setElementDistribution( 0.5*x_extent / (elem_sizes["small"]) )
    #<Note> slab and basin multiplied by approximate ratio of their x_extent.
    #___________________________________________________________________________




    #---------------------------------------------------------------------------
    # Now define the gmsh 'design' object. 
    #---------------------------------------------------------------------------

    design2D = pycad.gmsh.Design(dim=2, element_size=elem_sizes['large'], keep_files=False)
    
    # Also specify the domains for tagging with property values later on:
    design2D.addItems(   
    pycad.PropertySet( "air"   , domains[0]) ,   
    pycad.PropertySet( "lyr1"  , domains[1]) ,   
    pycad.PropertySet( "slab"  , domains[2]) ,   
    pycad.PropertySet( "basin" , domains[3]) ,
    pycad.PropertySet( "lyr2"  , domains[4]) ,
    pycad.PropertySet( "lyr3"  , domains[5]) )   
    
    # Now define the unstructured finley-mesh..
    model2D = finley.MakeDomain(design2D)  
    #___________________________________________________________________________


    return model2D

def generateCommemi4Mesh():
    #---------------------------------------------------------------------------
    # DESCRIPTION:
    # ------------
    # Script for preparing the COMMEMI-2 2D model.
    #
    # The COMMEMI-4 2D model consist of a 3-layered halfspace,
    # hosting an anomalous horizontal slab and a basin-structure 
    # in the first layer. 
    #
    # References:
    # -----------
    # See Franke A., p.89, 2003 (MSc. Thesis).
    # 
    # Antje Franke, "Zweidimensionale Finite-Elemente-Modellierung 
    # niederfrequenter elektromagnetischer Felder in der Fernzone", 
    # Diplomarbeit (MSc.), 2003, Technische Universtitaet Freiberg.
    #
    # --------------------------------------------------------------------------


    #---------------------------------------------------------------------------
    # Geometric mesh parameters.
    # --------------------------------------------------------------------------

    # Horizontal extent and start point of mesh.
    a_extent = 50000   # 50km - Vertical extent of air-layer in (m).
    z_extent = 50000   # 50km - Vertical extent of subsurface in (m).
    x_extent = 60000   # 60km - Horizontal extent of model in (m).

    # Start point of mesh.
    x_start  = 0 #-x_extent/2.0

    # Mesh elements sizes.
    elem_sizes = { 
                'large' : 4.00 * x_extent/100.0, # 
                'normal': 2.00 * x_extent/100.0, # 
                'small' : 0.25 * x_extent/100.0  # 
                }
   #____________________________________________________________________________





    #---------------------------------------------------------------------------
    # Coordinate definitions.
    # --------------------------------------------------------------------------

    # X-coordinates of all domain corners (in order of appearance, left to right).
    x0 = x_start                          # left         (@ boundary)
    x1 = x_start + 24000                  # centre       (@ slab/basin)
    x2 = x_start + 24000 + 8000           # edge-bottom  (@ slab/lyr1)
    x3 = x_start + 24000 + 8000 + 3000    # edge-top     (@ slab/lyr1)
    x4 = x_start + x_extent               # right        (@ boundary) 

    # Y-coordinates of all domain corners (in order of appearance, top to bottom).
    y0 = a_extent                         # top          
    y1 = 0                                # centre       (@ air/earth)
    y2 =-500                              # lyr1-bottom  (@ boundary-left) 
    y3 =-1000                             # basin-bottom (@ boundary-right) 
    y4 =-2000                             # slab-bottom  (@ boundary-left) 
    y5 =-4000                             # basin-bottom (@ centre)  
    y6 =-25000                            # lyr1-bottom 
    y7 =-z_extent                         # bottom

    # Save in dictionary as a list of tuples for each domain, from left-top corner, counterclockwise.
    coord = {                                 
            'air'  : ([x0, y0, 0],    # 0: left  , top
                        [x0, y1, 0],    # 1: left  , bottom (@ air/earth)
                        [x4, y1, 0],    # 2: right , bottom (@ air/earth)
                        [x4, y0, 0]),   # 3: right , top
                                        
            'lyr1' : ([x0, y1, 0],    # 0: left  , top    
                        [x0, y2, 0],    # 1: left  , bottom 
                        [x1, y2, 0],    # 2: right , bottom (@ slab/basin)
                        [x4, y2, 0],    # 3: right , bottom (@ boundary)
                        [x4, y1, 0]),   # 4: right , top 
                                            
            'slab' : ([x0, y2, 0],    # 0: left  , top    
                        [x0, y4, 0],    # 1: left  , bottom 
                        [x1, y4, 0],    # 2: right , bottom (@ slab/basin)
                        [x1, y2, 0]),   # 3: right , top    (@ slab/basin)
                                    
            'basin': ([x1, y2, 0],    # 0: left  , top    (@ slab/basin)
                        [x1, y4, 0],    # 1: left  , centre (@ slab/basin) 
                        [x1, y5, 0],    # 2: left  , bottom (@ lyr1/basin) 
                        [x2, y5, 0],    # 3: centre, bottom (@ lyr1/basin)        
                        [x3, y3, 0],    # 4: edge  , bottom (@ lyr1/basin)
                        [x4, y3, 0],    # 5: right , bottom (@ boundary)
                        [x4, y2, 0]),   # 6: right , top
                                    
            'lyr2' : ([x0, y4, 0],    # 0: left  , top    
                        [x0, y6, 0],    # 1: left  , bottom 
                        [x4, y6, 0],    # 2: right , bottom 
                        [x4, y3, 0],    # 3: right , top    (@ basin/boundary)
                        [x3, y3, 0],    # 4: edge  , top    (@ lyr2/basin)
                        [x2, y5, 0],    # 5: centre, top    (@ lyr2/basin)
                        [x1, y5, 0],    # 6: left  , top    (@ lyr2/basin)
                        [x1, y4, 0]),   # 7: left  , centre (@ slab/basin)
                                    
            'lyr3' : ([x0, y6, 0],    # 0: left  , top    
                        [x0, y7, 0],    # 1: left  , bottom 
                        [x4, y7, 0],    # 2: right , bottom 
                        [x4, y6, 0]),   # 3: right , top                   
            }
    #___________________________________________________________________________


    #---------------------------------------------------------------------------
    # Setup the COMMEMI-4 mesh.
    #---------------------------------------------------------------------------

    # This creates the mesh and saves it to the output folder.
    return setupMesh("TE", coord, elem_sizes)
    #___________________________________________________________________________



# ==========================================================
# ==========================================================

class Test_COMMEMI4(unittest.TestCase):
    @unittest.skipUnless(HAVE_FINLEY, "Test requires finley to be available")
    @unittest.skipUnless(HAVE_GMSH, "Test requires gmsh to be available")
    @unittest.skipUnless(HAVE_DIRECT, "Missing direct solver")
    @unittest.skipUnless(HAVE_SCIPY, "Test requires scipy to be available")
    def test_comm4(self):
        # ---
        # Initialisations
        # ---

        # Get timing:
        startTime = datetime.datetime.now()

        # Mode (TE includes air-layer, whereas TM does not):
        mode = 'TE'

        # Read the mesh file and define the 'finley' domain: 
        #mesh_file = "mesh/commemi-4/commemi4_tm.msh"
        #domain = finley.ReadGmsh(mesh_file, numDim=2)
        domain = generateCommemi4Mesh()

        #mesh_file = "mesh/commemi-4/commemi4_tm.fly"
        #domain = finley.ReadMesh(mesh_file)

        # Sounding frequencies (in Hz):
        freq_def = {"high":1.0e+0,"low":1.0e-0,"step":1}
        # Frequencies will be mapped on a log-scale from
        # 'high' to 'low' with 'step' points per decade.
        # (also only one frequency must be passed via dict)

        # Step sizes for sampling along vertical and horizontal axis (in m):
        xstep=100
        zstep=100



        # ---
        # Resistivity model
        # ---

        # Resistivity values assigned to tagged regions (in Ohm.m):
        rho  = [
                1.0e+14, # 0: air     1.0e-30
                25.0   , # 1: lyr1    0.04
                10.0   , # 2: slab    0.1
                2.5    , # 3: basin   0.4
                1000.0 , # 4: lyr2    0.001
                5.0      # 5: lyr3    0.2
            ]

        # Tags must match those in the file:
        tags = ["air", "lyr1", "slab", "basin", "lyr2", "lyr3"]

        # Optional user defined map of resistivity:
        def f4(x,z,r): return escript.sqrt(escript.sqrt(x*x+z*z))/r
        maps = [None, None, None, None, f4, None]



        # ---
        # Layer definitions for 1D response at boundaries.
        # ---

        # List with resistivity values for left and right boundary.
        rho_1d_left  = [ rho[0], rho[1], rho[2], rho[4], rho[5] ]
        rho_1d_rght  = [ rho[0], rho[1], rho[3], rho[4], rho[5] ]

        # Associated interfaces for 1D response left and right (must match the mesh file).
        ifc_1d_left = [ 50000, 0, -500, -2000, -25000, -50000]
        ifc_1d_rght = [ 50000, 0, -500, -1000, -25000, -50000]

        # Save in dictionary with layer interfaces and resistivities left and right:
        ifc_1d = {"left":ifc_1d_left , "right":ifc_1d_rght}
        rho_1d = {"left":rho_1d_left , "right":rho_1d_rght}



        # ---
        # Adjust parameters here for TM mode
        # ---

        # Simply delete first element from lists:
        if mode.upper() == 'TM':
            tags.pop(0)
            rho.pop(0)
            rho_1d['left'].pop(0)
            rho_1d['right'].pop(0)
            ifc_1d['left'].pop(0)
            ifc_1d['right'].pop(0)
            if maps is not None:
                maps.pop(0)



        # ---
        # Run MT_2D
        # ---

        # Class options:
        mt2d.MT_2D._solver = "DIRECT" #"ITERATIVE" #"CHOLEVSKY" #"CGLS " #"BICGSTAB" #"DIRECT" "ITERATIVE"
        mt2d.MT_2D._debug   = False

        # Instantiate an MT_2D object with required & optional parameters:
        obj_mt2d = mt2d.MT_2D(domain, mode, freq_def, tags, rho, rho_1d, ifc_1d,
                xstep=xstep ,zstep=zstep, maps=None, plot=False)

        # Solve for fields, apparent resistivity and phase:
        mt2d_fields, arho_2d, aphi_2d = obj_mt2d.pdeSolve()
        
        
        #import random

        #mt2d_fields[0]['real']+=random.random()
        #mt2d_fields[0]['imag']+=50*random.random()

        #print(arho_2d[0][0])
        #for i in range(len(aphi_2d[0])):
            #aphi_2d[0][i]+=(50*random.random())

        #for i in range(len(arho_2d[0])):
            #arho_2d[0][i]-=7*(random.random())    
        

        # ---
        # User defined plots
        # ---

        # Setup abscissas/Ordinates for escript data:
        x  = numpy.array( obj_mt2d.loc.getX() )[:,0]
        y0 = numpy.array( obj_mt2d.loc.getValue(arho_2d[0]) )
        y1 = numpy.array( obj_mt2d.loc.getValue(aphi_2d[0]) )

        # Sort these arrays and delete any duplicates to prevent the call to 
        # InterpolatedUnivariateSpline below from returning an error
        indices = numpy.argsort(x,kind='quicksort')
        x = x[indices]
        y0 = y0[indices]
        y1 = y1[indices]
        indices = numpy.unique(x,return_index=True)
        x = x[indices[1]]
        y0 = y0[indices[1]]
        y1 = y1[indices[1]]

        # Zhdanov et al, 1997, -- Model 2D-1 Table B.33. Model2D-4 (T=1.0, z=0), see 
        # "Methods for modelling electromagnetic fields. Results from COMMEMI -- the
        # international project on the comparison of modelling results for electromag-
        # netic induction", Journal of Applied Geophysics, 133-271
        rte = [12.70, 12.00, 8.80, 6.84, 6.67, 6.25] # TE rho_a (3 Canada)   
        rtm = [11.40, 11.50, 9.03, 6.78, 6.80, 5.71] # TM rho_a (3 Canada) 
        if mode.lower() == 'te':
            ra = rte
        else:
            ra = rtm  
        # Associated stations shifted to match escript coordinates:
        xs = numpy.array( [-10, -7, -6, -5, 2, 5] )*1000 + x.max()/2.0

        # Setup interpolation to get values at specified stations (for comparison):
        fi = InterpolatedUnivariateSpline(x, y0)
        # Save esscript values at comparison points in text file:
        # uncomment to investigate
        #numpy.savetxt("mesh/commemi-4/commemi4_"+mode.lower()+".dat", numpy.column_stack((xs,fi(xs))), fmt='%g')

        # X plot-limits:
        x0lim = [2000,38000]
        #y1lim = [0,120]
        #y2lim = [40,85]
        
        # Plot labels:
        title = '    escript COMMEMI-4 MT-2D ' + '(' + mode.upper() + ')' + ' freq: ' + str(obj_mt2d.frequencies[0]) + ' Hz'
        ylbl0 = r'resistivity $(\Omega m)$'
        ylbl1 = r'phase $(\circ)$'
        xlbl1 = 'X (m)'
        if HAVE_MPL:
            # Setup the plot window with app. res. on top and phase on bottom:
            f, ax = pyplot.subplots(2, figsize=(3.33,3.33), dpi=1200,
                    facecolor='w', edgecolor='k', sharex=True) # Mind shared axis
            f.subplots_adjust(hspace=0.1, top=0.95, left=0.135, bottom=0.125, right=0.975)  
            f.suptitle(title, y=0.99,fontsize=8) # 
                
            # Top: apparent resistivity and points from Weaver for comparison:
            ax[0].plot(x, y0, color='red',  label = 'escript')
            ax[0].plot(xs,ra, linestyle='', markersize=3, marker='o',
                    color='blue', label = 'Weaver') 
            ax[0].grid(b=True, which='both', color='grey',linestyle=':')
            ax[0].set_ylabel( ylbl0)
            ax[0].yaxis.set_label_coords(-0.082, 0.5)
            # Plot limits:
            ax[0].set_xlim(x0lim)      
            #ax[0].set_ylim(y1lim)    

            # Bottom: phase on linear plot
            ax[1].plot(x,y1, color='blue')
            ax[1].grid(b=True, which='both', color='grey',linestyle=':')
            ax[1].set_xlabel( xlbl1 )
            ax[1].set_ylabel( ylbl1 )
            # Plot limits:
            ax[1].set_xlim(x0lim)      
            #ax[1].set_ylim(y2lim)     

            # ask matplotlib for the plotted objects and their labels
            lna, la = ax[0].get_legend_handles_labels()
            ax[0].legend(lna, la, bbox_to_anchor=(0.02, 0.325), loc=2,
                    borderaxespad=0.,prop={'size':8}, frameon=False)

            pyplot.ticklabel_format(style='sci', axis='x', scilimits=(0,0),
                    useMathText=True)
            ax[0].xaxis.major.formatter._useMathText = True
            pyplot.rc('font', **{'size': 8,'family':'sans-serif'})
        #uncomment to allow visual inspection
        #f.savefig("mesh/commemi-4/commemi4_"+mode.lower()+".png", dpi=1200)

        # Now let's see if the points match
        # First, we need to find correspondance between xs and x
        indices=[]
        for i in range(len(xs)):
            mindiff=40000
            mindex=0
            for j in range(len(x)):
                if abs(xs[i]-x[j]) < mindiff: 
                    mindiff=abs(xs[i]-x[j])
                    mindex=j
            indices.append(mindex)
            
        # The following are very simple checks based on the visual shape of the correct result        
        maxdiff=0
        for i in range(len(indices)):
            if abs(y0[indices[i]]-ra[i])>maxdiff:
                maxdiff=abs(y0[indices[i]]-ra[i])
            
        if maxdiff>5:           #Threshold is pretty arbitrary
            raise RuntimeError("Mismatch with reference data")

        c=0
        for y in y1:
            if y<46:
                c+=1

        if not (62 < escript.Lsup(y1) < 64):
            raise RuntimeError("Peak of bottom plot is off.")
            
        if not (0.61 < c/len(y1) < 0.65):
            print(c/len(y1))
            raise RuntimeError("Bottom plot has too many high points")

        #
        print("Runtime:", datetime.datetime.now()-startTime)
        print("Done!")


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
