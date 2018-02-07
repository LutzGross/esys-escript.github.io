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


"""
Test script to run test model COMMEMI-1
"""

from __future__ import print_function, division

import matplotlib
# The following line is here to allow automated testing. Remove or comment if
# you would like to display the final plot in a window instead.
matplotlib.use('agg')

import datetime
import numpy
import esys.downunder.magtel2d as mt2d
import esys.escript            as escript
import esys.escript.pdetools   as pdetools

try:
    import esys.finley         as finley
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False

HAVE_DIRECT = escript.hasFeature("PASO_DIRECT") or escript.hasFeature('trilinos')

#-------------------------------------------------------------
# The following functions create the mesh used by this example
#-------------------------------------------------------------


def makeLayerCake(x_start,x_extent,z_layers):
    # --------------------------------------------------------------------------
    # DESCRIPTION:
    # -----------
    # This is a utility function which sets up a 2D model with N layers.
    # 
    # ARGUMENTS:                                                              
    # ----------
    # x_start             :: start coordinate of mesh.
    # x_extent            :: horizontal extent of mesh.
    # z_layers            :: list with interface coordinates.
    #
    # RETURNS:
    # --------
    # borders             :: borders of layers. 
    # air_earth_interface :: line at the air/earth interface.
    #
    # AUTHOR:
    # -------
    # Ralf Schaa, 
    # University of Queensland
    #
    #
    # HISTORY:
    # --------
    #
    # --------------------------------------------------------------------------

    import esys.pycad   as pycad     # @UnresolvedImport
    import esys.weipa   as weipa     # @UnresolvedImport    
    import esys.finley  as finley    # @UnresolvedImport
    import esys.escript as escript   # @UnresolvedImport

         
    # --------------------------------------------------------------------------
    # Point definitions.
    # --------------------------------------------------------------------------
         
    # Loop through all layers and define the vertices at all interfaces.
    scale = 1.0
    points = []
    for i in range(0,len(z_layers)):
            # Adjust scale at corners of air/earth interface:
            if z_layers[i] == 0:
                scale = 0.15
            else:
                scale = 1.0
            points.append( pycad.Point(x_start           , z_layers[i], 0.0, local_scale = scale) ) # Left-Corner.     
            points.append( pycad.Point(x_start + x_extent, z_layers[i], 0.0, local_scale = scale) ) # Right-Corner. 


    # --------------------------------------------------------------------------
    # Line definitions.
    # --------------------------------------------------------------------------

    # Now connect the points to define the horizontal lines for all interfaces:
    hlines = []
    for i in range(0,len(points),2):
        if i <= len(points)-1:
            hlines.append( pycad.Line(points[i],points[i+1]) )     
    
    # Now connect the points to define the vertical lines for all interfaces:
    vlines_left = []
    for i in range(0,len(points),2):
        if i <= len(points)-3:
            vlines_left.append( pycad.Line(points[i],points[i+2]) )     

    vlines_right = []
    for i in range(0,len(points),2):
        if i <= len(points)-4:
            vlines_right.append( pycad.Line(points[i+1],points[i+3]) )     



    # --------------------------------------------------------------------------
    # Curveloop and Area definitions.
    # --------------------------------------------------------------------------

    # Join line segments for each layer.          
    borders = []
    for i in range(0,len(z_layers)-1):
        border = [ hlines[i],vlines_right[i],-hlines[i+1],-vlines_left[i] ]
        borders.append( pycad.CurveLoop( border) )       



    # --------------------------------------------------------------------------
    # Return values.
    # --------------------------------------------------------------------------

    # Explicitly specify the air-earth-boundary:
    air_earth_interface = hlines[1]
    
    return borders, air_earth_interface
                                      
#_______________________________________________________________________________




def setupMesh(mode, x_start, x_extent, a_extent, z_layers, anomaly_coord, elem_sizes):
    # --------------------------------------------------------------------------
    # DESCRIPTION:
    # -----------
    # This is a utility function which sets up the COMMEMI-1 mesh.
    # 
    #
    # ARGUMENTS:                                                              
    # ----------
    # mode           :: TE or TM mode.
    # x_start        :: horizontal start-point mesh.
    # x_extent       :: horizontal extent of mesh.
    # a_extent       :: vertical extent of air-layer.
    # z_layers       :: list with coordinates of top-interfaces in Z-direction, incl. basement.
    # anomaly_coord  :: dictionary with coordinate tuples of anomalies, counterclockwise.
    # elem_sizes     :: mesh element sizes, large, normal, small. 
    #
    # RETURNS:
    # --------
    # <Nothing> A mesh file is written to the output folder.
    # 
    #
    # AUTHOR:
    # -------
    # Ralf Schaa, 
    # The University of Queensland
    #
    #
    # HISTORY:
    # --------
    #
    # --------------------------------------------------------------------------



    # --------------------------------------------------------------------------
    # Imports.
    # --------------------------------------------------------------------------
    
    # System imports.
    import math
    
    # Escript modules.
    import esys.pycad              as pycad     # @UnresolvedImport   
    import esys.finley             as finley    # @UnresolvedImport
    import esys.escript            as escript   # @UnresolvedImport
    import esys.weipa              as weipa     # @UnresolvedImport    
    # <Note>: "@UnresolvedImport" ignores any warnings in Eclipse/PyDev (PyDev has trouble with external libraries).
    
    # Warn about magnetotelluric TM mode:
    if mode.lower() == 'tm':
        print("TM mode not yet supported")
        return None
        
    # --------------------------------------------------------------------------
    # Anomaly border.
    # --------------------------------------------------------------------------
     
    #<Note>: define the anomaly which must be 'cut out' in the main mesh. 
    
    
    # Prepare list to store the anomaly borders:
    border_anomaly = []
                
    # Cycle anomaly dictionary and define the border for each.
    for anomaly in anomaly_coord:
        
        # Extract the coordinates for current key:
        coord = anomaly_coord[anomaly]
            
        # Points defining the anomaly from left-top.
        points0 = []
        for i in range( 0, len(coord) ):            
            points0.append(pycad.Point(coord[i][0], coord[i][1], 0.0))
    
        # Define the line segments connecting the points.
        lines0 = []
        for i in range( 0, len(points0)-1 ):
            lines0.append(pycad.Line(points0[i],points0[i+1]))
        # Connect the last segment from end to start:    
        lines0.append(pycad.Line(points0[-1], points0[0])) 
        
        # And define the border of the anomalous area.
        border_anomaly.append( pycad.CurveLoop(*lines0) ) 
    #___________________________________________________________________________




    # --------------------------------------------------------------------------
    # Get the borders for each layer (air & host).
    # --------------------------------------------------------------------------

    # Borders around layers and the air/earth interface.
    borders, air_earth_interface = makeLayerCake(x_start,x_extent,z_layers)





    # --------------------------------------------------------------------------
    # Specification of number of elements in domains.
    # --------------------------------------------------------------------------
        
    #<Note>: specifying the number of mesh elements is somewhat heuristic 
    #        and is dependent on the mesh size and the anomaly sizes. 

    coord = anomaly_coord["anomaly_1"]
     
    # First get the max-length of the anomaly to specify the number of elements.
    length = max(( abs(coord[2][0]-coord[0][0]) ),  # X-length
                 ( abs(coord[2][1]-coord[0][1]) ))  # Y-length                 
    
    # Specify number of elements in air, anomaly and on air/earth interface:        
    nr_elements_air       = 1 * x_extent / elem_sizes["large"]
    nr_elements_anomaly   = 2 * length   / elem_sizes["small"]
    nr_elements_interface = 4 * x_extent / elem_sizes["small"]
    #___________________________________________________________________________
    
    
    
     
    # --------------------------------------------------------------------------
    # Domain definitions.
    # --------------------------------------------------------------------------

    # Define the air & layer areas; note the 'holes' specifiers.
    domain_air     = pycad.PlaneSurface( borders[0] )   
    domain_host    = pycad.PlaneSurface( borders[1] , holes = [ border_anomaly[0] ] )    
    domain_anomaly = pycad.PlaneSurface( border_anomaly[0] )    
        
    # Specify the element sizes in the domains and along the interface.
    #<Note>: Sizes must be assigned in the order as they appear below:    
    domain_air.setElementDistribution( nr_elements_air )         
    domain_anomaly.setElementDistribution( nr_elements_anomaly ) 
    air_earth_interface.setElementDistribution( nr_elements_interface )

    # Ready to define the mesh-design..
    design2D = pycad.gmsh.Design(dim=2, element_size=elem_sizes["normal"] , keep_files=False)
    # ..and also specify the domains for tagging with property values later on:
    design2D.addItems( pycad.PropertySet("domain_air"    , domain_air),
                       pycad.PropertySet("domain_host"   , domain_host),
                       pycad.PropertySet("domain_anomaly", domain_anomaly) ) 
    
    # Now define the unstructured finley-mesh..
    model2D = finley.MakeDomain(design2D)
    #___________________________________________________________________________


    return model2D    
    #___________________________________________________________________________
    
def generateCommemi1Mesh():
    # --------------------------------------------------------------------------
    # Geometric mesh parameters.
    # --------------------------------------------------------------------------

    # Mesh extents.
    a_extent = 20000    # 20km - Vertical extent of air-layer in (m).
    z_extent = 20000    # 20km - Vertical extent of subsurface in (m).
    x_extent = 40000    # 40km - Horizontal extent of mesh in (m).

    # Start point of mesh.
    x_start = 0 #-x_extent/2.0

    # Define interface locations in z-direction: top, air/earth, basement. 
    z_layers    = [   a_extent, 0, -z_extent]

    # Mesh elements sizes.
    elem_sizes = { 
                'large' : 10.00 * x_extent/100.0, # 5.00% of x_extent.
                'normal': 05.00 * x_extent/100.0, # 2.50% of x_extent.
                'small' : 00.50 * x_extent/100.0  # 0.25% of x_extent.
                }
    #___________________________________________________________________________




    # --------------------------------------------------------------------------
    # Geometric anomaly parameters.
    # --------------------------------------------------------------------------

    # Extents of the rectangular 2D anomaly.
    x_anomaly = 1000    # 1km - Horizontal extent of anomaly in (m).
    z_anomaly = 2000    # 2km - Vertical extent of anomaly in (m).

    # Coordinates of the rectangular 2D anomaly.
    ya1 = -250                                    # Top
    ya2 = -z_anomaly + ya1                        # Bottom
    xa1 = x_start + x_extent/2.0 - x_anomaly/2.0  # Left
    xa2 = x_start + x_extent/2.0 + x_anomaly/2.0  # Right

    # Save in dictionary as a list of tuples from left-top corner, counterclockwise.
    anomaly_coord = { 
                    'anomaly_1': ([xa1,ya1],[xa1,ya2],[xa2,ya2],[xa2,ya1]) 
                    }
    #___________________________________________________________________________




    # --------------------------------------------------------------------------
    # Setup the COMMEMI-1 mesh.
    # --------------------------------------------------------------------------

    # This creates the mesh and saves it to the output folder.
    return setupMesh("TE", x_start, x_extent, a_extent, z_layers,  anomaly_coord, elem_sizes)
    #___________________________________________________________________________



#-------------------------------------------------------------
# End of mesh set up functions
#-------------------------------------------------------------


if HAVE_FINLEY:
    # ---
    # Initialisations
    # ---

    # Get timing:
    startTime = datetime.datetime.now()

    # Mode (TE includes air-layer, whereas TM does not):
    mode = 'TE'

    # Read the mesh file and define the 'finley' domain:
    #mesh_file = "data/commemi1_te.fly"
    #domain = finley.ReadMesh(mesh_file, numDim=2)
    if escript.hasFeature('gmsh'):
        domain = generateCommemi1Mesh()

    # Sounding frequencies (in Hz):
    freq_def = {"high":1.0e+1,"low":1.0e+1,"step":1}
    # Frequencies will be mapped on a log-scale from
    # 'high' to 'low' with 'step' points per decade.
    # (also only one frequency must be passed via dict)

    # Step sizes for sampling along vertical and horizontal axis (in m):
    xstep=400
    zstep=200

    # ---
    # Resistivity model
    # ---

    # Resistivity values assigned to tagged regions (in Ohm.m):
    rho  = [
            1.0e+14, # 0: air
            100.0  , # 1: host
              0.5    # 2: anomaly
           ]

    # Tags must match those in the file:
    tags = ["domain_air", "domain_host", "domain_anomaly"]


    # ---
    # Layer definitions for 1D response at boundaries.
    # ---

    # List with resistivity values for left and right boundary.
    rho_1d_left  = [ rho[0], rho[1] ]
    rho_1d_rght  = [ rho[0], rho[1] ]

    # Associated interfaces for 1D response left and right (must match the mesh file).
    ifc_1d_left = [ 20000, 0, -20000]
    ifc_1d_rght = [ 20000, 0, -20000]

    # Save in dictionary with layer interfaces and resistivities left and right:
    ifc_1d = {"left":ifc_1d_left , "right":ifc_1d_rght}
    rho_1d = {"left":rho_1d_left , "right":rho_1d_rght}

    # ---
    # Run MT_2D
    # ---

    # Class options:
    mt2d.MT_2D._solver = "DIRECT"
    mt2d.MT_2D._debug   = False

    if mt2d.MT_2D._solver == "DIRECT" and not escript.hasFeature('paso'):
        print("Trilinos direct solvers cannot currently handle PDE systems. Please compile with Paso.")
    elif mt2d.MT_2D._solver == "DIRECT" and not HAVE_DIRECT:
        if escript.getMPISizeWorld() > 1:
            print("Direct solvers and multiple MPI processes are not currently supported.")
        else:
            print("escript was not built with support for direct solvers, aborting.")
    elif not escript.hasFeature('gmsh'):
        print("This example requires gmsh, aborting.")
    else:
        # Instantiate an MT_2D object with required & optional parameters:
        obj_mt2d = mt2d.MT_2D(domain, mode, freq_def, tags, rho, rho_1d, ifc_1d,
                xstep=xstep ,zstep=zstep, maps=None, plot=True)

        # Solve for fields, apparent resistivity and phase:
        mt2d_fields, arho_2d, aphi_2d = obj_mt2d.pdeSolve()

        #
        print(datetime.datetime.now()-startTime)
        print("Done!")

else: # no finley
    print("Finley module not available.")


