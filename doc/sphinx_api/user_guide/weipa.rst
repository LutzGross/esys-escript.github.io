================
The weipa Module
================

The **esys.weipa** Module and Data Visualization
==================================================

.. _chap:weipa:

The *weipa* C++ library and accompanying *Python* module allow exporting **esys.escript** “Data“ objects and their domain in a format suitable for visualization. Besides creating output files, *weipa* can also interface with the *VisIt* visualization software. This allows accessing the latest simulation data while the simulation is still running without the need to save any files.

The “EscriptDataset“ class
--------------------------

.. container:: classdesc

   EscriptDataset holds an *escript* dataset including a domain and data variables for a single time step and offers methods to export the data in various formats. It is preferable to create a dataset object using the “createDataset“ function from **esys.weipa** (see Section :ref:‘sec:weipafuncs‘) rather than using the (non-exposed) *Python* constructor for the class.

The following methods are available:

.. container:: methoddesc

   setDomaindomain sets the “Domain“ for this dataset. Note that the domain can only be set once and all “Data“ objects added to this dataset must be defined on the same domain.

.. container:: methoddesc

   addDatadata, name adds the “Data“ object “data“ to this dataset which will be exported by the given “name“. Some export formats support data units which can be set through the “units“ parameter, e.g. “"km/h"\`‘. Before calling this method a domain must be set with “setDomain“ and all “Data“ objects added must be defined on the same domain. There is no restriction, however, on the “FunctionSpace“ used.

.. container:: methoddesc

   setCycleAndTimecycle, time sets the cycle and time values for this dataset. The cycle is an integer value which usually corresponds with the loop counter of the simulation script. That is, every time a new data file is created this counter is incremented. The value of “time“ on the other hand is a floating point number that encodes some form of simulation time. Both, cycle and time may be read by analysis tools and shown alongside other metadata to the user.

.. container:: methoddesc

   setMeshLabelsx, y sets the labels of the X, Y, and Z axis. By default, visualization tools display default strings such as "X-Axis" or "X" along the axes. Some export formats allow overriding these with more specific strings such as "Width", "Horizontal Distance", etc.

.. container:: methoddesc

   setMeshUnitsx, y sets the units to be displayed along the X, Y, and Z axis in visualization tools (if supported). Not all export formats will use these values.

.. container:: methoddesc

   setMetadataSchemaString adds custom metadata and/or XML schema strings to VTK files. The content of “schema“ is added to the top-level *VTKFile* element so care must be taken to keep the resulting file valid. As an example, “schema“ may contain the string “xmlns:gml="http://www.opengis.net/gml"\`‘. The content of “metadata“ will be written enclosed in “<MetaData>“ tags. Thus, a valid example would be “<dataSource>something</dataSource>“. Note that these values are ignored by other exporters.

.. container:: methoddesc

   saveSilofilename saves the dataset in the *SILO* file format to a file named “filename“. The file extension “.silo“ will be automatically added if not present.

.. container:: methoddesc

   saveVTKfilename saves the dataset in the *VTK* file format to a file named “filename“. The file extension “.vtu“ will be automatically added if not present. Certain combinations of function spaces cannot be written to a single *VTK* file due to format restrictions. In these cases this method will save separate files where each file contains compatible data. The function space name is appended to the filename to distinguish them.

Functions
---------

.. _sec:weipafuncs:

.. container:: funcdesc

   createDatasetdomain, \**data creates an “EscriptDataset“ object, sets its domain, populates it with the given “Data“ objects and returns it. Note that it is not possible to set units for the data variables added with this function. If this is required, it is recommended to call this function with a domain only and use the “addData“ method subsequently.

.. container:: funcdesc

   saveSilofilename , \**data convenience function that creates a dataset with the given domain and “Data“ objects and saves it to a file in the *SILO* file format. If “domain“ is “None“ the domain will be determined by the “Data“ objects. See the “setDomain“, “addData“, and “saveSilo“ methods of the “EscriptDataset“ class for details.

.. container:: funcdesc

   saveVoxetfilename, \**data saves “Data“ objects defined on a *ripley* grid in the file format suitable for import into *GOCAD*. A dataset consists of a header file (extension ``.vo``) and one property file (with no file extension) for each “Data“ object.

.. container:: funcdesc

   visitInitializesimFile initializes the *VisIt* simulation interface which is responsible for the communication with a *VisIt* client. This function will create a file by the name given via “simFile“ (extension “.sim2“) which can be loaded by a compatible *VisIt* client in order to connect to the simulation. The optional “comment“ string is forwarded to the client. Note that this function only succeeds if *escript* was compiled with support for *VisIt* and the appropriate libraries are found in the runtime environment. Clients wanting to connect can only do so if the version number matches the version number used to compile **esys.weipa**. Calling this function does not make any data available yet, see the “visitPublishData“ function.

.. container:: funcdesc

   visitPublishDatadataset publishes an “EscriptDataset“ object through the *VisIt* simulation interface, checks for client requests and handles any outstanding ones. Before publishing any data, the “visitInitialize“ function must be called to set up the interface. Since this function not only publishes new data but polls for incoming connections and handles requests, it should be called as often as practical (even with the same dataset) to avoid timeout errors from clients. On the other hand it should be noted that the same process(es) deal with visualization requests that run your simulation. So a request for an expensive task by a *VisIt* client will pause the simulation code while it is being processed.

Visualizing *escript* Data
--------------------------

This section gives a very brief overview on how data exported through **esys.weipa** can be visualized. While there are many visualization packages available that are compatible with *VTK* and *SILO* files produced by *escript*, this discussion will refer to *VisIt* , an actively maintained open source package optimally suited to visualize and analyze large datasets both interactively and through *Python* scripts. You can find a number of manuals, a wiki page and links to mailing lists on the *VisIt* website. It is assumed that you have a working *VisIt* installation that can be started by entering “visit“ on the command line.

The examples that follow will use the output produced by the Elastic Deformation example from section `[sec:ElasticChap] <#sec:ElasticChap>`__ (“heatedblock.py“ in the example directory) which produces the file “deform.vtu“. This *VTK* file contains a 3D scalar variable called “stress“ and a vector variable called “disp“, among others.

Using the *VisIt* GUI
~~~~~~~~~~~~~~~~~~~~~~~

Start the VisIt graphical user interface and open the file “deform.vtu“ via the ’File’ menu. Alternatively, you can directly open the file on startup by issuing

::

   visit -o deform.vtu

You should see the *VisIt* GUI on the left hand side and an empty visualization window on the right. Click on ’Add’ under Plots in the GUI to bring up a menu of plot types, then click on ’Pseudocolor’ and select ’stress’. This will add a plot to the list which maps values of the ’stress’ variable to colors. Note, that the plot will not be generated until you click on the ’Draw’ button in the GUI. You should now see a coloured box in the visualization window which you can rotate around and inspect from different angles using your mouse. The example uses a coarse mesh of 10 by 10 by 10 elements which are clearly visible in this plot.

We can improve the visual effect by enabling interpolation between the elements. To do so, bring up the plot attributes by double-clicking the ’Pseudocolor - stress’ plot entry in the GUI. Next, select ’Nodal’ under ’Centering’, click on ’Apply’ and dismiss the dialog. Notice how the colours now smoothly blend into each other and the element boundaries are no longer visible.

Now we will add arrows to visualize the displacement vectors. Click on ’Add’ and under ’Vector’ select ’disp’. Once again click on ’Draw’ to execute the new plot. By default only few vectors are shown but since the mesh is very coarse we can tell *VisIt* to draw all available vectors. Bring up the Vector plot attributes (double-click on the plot as before) and under ’Vector amount’ select ’Stride’, leaving the parameter as 1. Click on ’Apply’ and dismiss the dialog.

As a final step we would like to see inside the plot. One possibility to do so is slicing. However, we want to keep all vectors while slicing only the Pseudocolor plot. In *VisIt* slicing is one of the Operators that may be added to plots and by default, Operators are added to *all* plots. To change this behaviour, uncheck the ’Apply operators to all plots’ box which is located underneath the plot list in the GUI. Then select the Pseudocolor plot, bring up the Operators menu by clicking on ’Operators’ and select ’ThreeSlice’ from the ’Slicing’ submenu. Again, click on ’Draw’ to update the plots and notice how the box has now been sliced. We can move the slices to more suitable positions by editing the operator attributes. Click on the little triangle to the left of the Pseudocolor plot to reveal the list of elements that have been applied to it. Next, double-click the ’ThreeSlice’ element to bring up the attribute window. Change the values to :math:`X=0.3` and :math:`Y=0.3`, leaving :math:`Z=0`. Apply the changes and dismiss the dialog to see the result.

You can now create an image of the plots as shown in the window. First, adjust the save options to your needs in the ’Set Save options’ dialog which is accessible from the ’File’ menu. Then select ’Save Window’ and you should find an image file with the name and location as entered in the options dialog.

Using the *VisIt* CLI (command line interface)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We will now perform exactly the same steps as in the last section but using the *Python* interface of *VisIt* instead of the GUI. Start up the CLI by issuing

::

   visit -cli

You should now see an empty visualization window but unlike in the previous section there will be no graphical user interface but a *Python* command line instead. Enter the following commands, one by one, noticing the changes in the visualization window after every block of commands:

.. code:: python

   OpenDatabase("deform.vtu")
     AddPlot("Pseudocolor","stress")
     DrawPlots()

     p=PseudocolorAttributes()
     p.centering=p.Nodal
     SetPlotOptions(p)

     AddOperator("ThreeSlice")
     DrawPlots()

     t=ThreeSliceAttributes()
     t.x=0.3
     t.y=0.3
     SetOperatorOptions(t)

     AddPlot("Vector", "disp")
     DrawPlots()

     v=VectorAttributes()
     v.useStride=1
     SetPlotOptions(v)

     s=SaveWindowAttributes()
     #change settings as required
     SaveWindow()
     exit()

All but the last call to “DrawPlots()“ is not required and was only put there for demonstrating the effects of the commands. You can save these commands to a file, e.g. “deformVis.py“ and let *VisIt* process them non-interactively like so:

::

   visit -cli -nowin -s deformVis.py

The “-nowin“ option prevents the visualization window from being shown which is not required since the purpose of the script is to save an image file.

Obviously, we have barely touched on the powerful features of *VisIt* and this section was only meant to give you a minimal introduction. The *VisIt* website has a reference manual for the *Python* interface that explains how to perform other operations programmatically, such as changing the view.
