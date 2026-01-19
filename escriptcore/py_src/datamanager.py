##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

"""
An escript data import and export manager.

:var __author__: name of authors
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point to documentation
"""


__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"
__author__="Lutz Gross, Cihan Altinay"

import pickle
import os
import shutil
from . import util
from . import escriptcpp as esc

class DataManager(object):
    """
    Escript data import/export manager.

    Example::

        dm=DataManager(formats=[DataManager.RESTART,DataManager.VTK])
        if dm.hasData():
            dom = dm.getDomain()
            time = dm.getValue("time")
            dt = dm.getValue("dt")
            T = dm.getValue("T")
            u = dm.getValue("u")
        else:
            T = ...
            u = ...
        dm.addData(time=time,dt=dt,T=T,u=u) # add data and variables
        dm.setTime(time)                    # set the simulation timestamp
        dm.export()                         # write out data
    """

    RESTART, SILO, VISIT, VTK = list(range(4))

    def __init__(self, formats=[RESTART], work_dir=".", restart_prefix="restart", do_restart=True):
        """
        Initialises the data manager. If do_restart is True and a restart
        directory is found the contained data is loaded (hasData() returns True)
        otherwise restart directories are removed (hasData() returns False).
        Values are only written to disk when export() is called.

        :param formats: A list of export file formats to use. Allowed values
                        are RESTART, SILO, VISIT, VTK.
        :param work_dir: top-level directory where files are exported to
        :param restart_prefix: prefix for restart directories. Will be used to
                               load restart files (if do_restart is True) and
                               store new restart files (if RESTART is used)
        :param do_restart: whether to attempt to load restart files
        """
        self._metadata=""
        self._md_schema=""
        self._data={}
        self._domain=None
        self._meshlabels=["","",""]
        self._meshunits=["","",""]
        self._stamp={}
        self._time=0.
        self._restartdir=None
        self._N=-1
        self._checkpointfreq=1
        self._myrank=esc.getMPIRankWorld()
        self._exportformats=set(formats)
        self._restartprefix=restart_prefix
        self._workdir=work_dir
        util.mkDir(self._workdir)
        if self.VISIT in self._exportformats:
            simFile=os.path.join(self._workdir, "escriptsim.sim2")
            if not self.__initVisit(simFile, "Escript simulation"):
                print("Warning: Could not initialize VisIt interface")
                self._exportformats.remove(self.VISIT)
        if self.RESTART in self._exportformats:
            # find all restart directories
            restart_folders = []
            for f in os.listdir(self._workdir):
                if f.startswith(self._restartprefix):
                    if os.path.isdir(f):
                        restart_folders.append(f)
            # remove unneeded restart directories
            if len(restart_folders)>0:
                restart_folders.sort()
                if do_restart:
                    self._restartdir=restart_folders[-1]
                    print(("Restart from "+os.path.join(self._workdir, self._restartdir)))
                    for f in restart_folders[:-1]:
                        self.__removeDirectory(f)
                    self.__loadState()
                else:
                    for f in restart_folders:
                        self.__removeDirectory(f)

    def addData(self, **data):
        """
        Adds 'escript.Data' objects and other data to be exported to this
        manager.

        :note: This method does not make copies of Data objects so
               any modifications will be carried over until export() is called.
        """
        # if this is the first addition after a restart, clear data first
        if self._restartdir != None:
            self.__clearData()

        for name,var in sorted(data.items(), key=lambda x: x[0]):
            if hasattr(var, "getDomain"):
                if self._domain is None:
                    self._domain=var.getDomain()
                elif self._domain != var.getDomain():
                    raise ValueError("addData: Data must be on the same domain!")
                self._data[name]=var
            else:
                self._stamp[name]=var

    def setDomain(self, domain):
        """
        Sets the domain without adding data.

        :param domain: the escript domain to set
        :type domain: `Domain`
        :raise ValueError: if a different domain has already been set
        """
        if self._domain is None:
            self._domain = domain
        elif self._domain != domain:
            raise ValueError("setDomain: Domain already set!")


    def hasData(self):
        """
        Returns True if the manager holds data for restart.

        :return: ``True`` if restart data is available, ``False`` otherwise
        :rtype: ``bool``
        """
        return self._restartdir != None

    def getDomain(self):
        """
        Returns the domain as recovered from restart files.

        :return: the domain loaded from restart files
        :rtype: `Domain`
        :raise ValueError: if no restart data is available
        """
        if not self.hasData():
            raise ValueError("No restart data available")
        return self._domain 

    def getValue(self, value_name):
        """
        Returns an 'escript.Data' object or other value that has been loaded
        from restart files.

        :param value_name: name of the value to retrieve
        :type value_name: ``str``
        :return: the requested data object or value
        :rtype: `Data` or other type depending on what was stored
        :raise ValueError: if no restart data is available
        """
        if not self.hasData():
            raise ValueError("No restart data available")

        if value_name in self._stamp:
            return self._stamp[value_name]

        ff=self.__getDumpFilename(value_name, self._restartdir)
        var = esc.load(ff, self._domain)
        #print("Value %s recovered from %s."%(value_name, ff))
        return var 

    def getCycle(self):
        """
        Returns the export cycle (=number of times export() has been called).

        :return: the current export cycle number
        :rtype: ``int``
        """
        return self._N

    def setCheckpointFrequency(self, freq):
        """
        Sets the number of calls to export() before new restart files are
        generated.

        :param freq: checkpoint frequency (1 = every export, 2 = every other, etc.)
        :type freq: ``int``
        """
        self._checkpointfreq=freq

    def setTime(self, time):
        """
        Sets the simulation timestamp.

        :param time: the current simulation time
        :type time: ``float``
        """
        self._time = time

    def setMeshLabels(self, x, y, z=""):
        """
        Sets labels for the mesh axes. These are currently only used by the
        Silo exporter.

        :param x: label for the x-axis
        :type x: ``str``
        :param y: label for the y-axis
        :type y: ``str``
        :param z: label for the z-axis (optional for 2D)
        :type z: ``str``
        """
        self._meshlabels=[x,y,z]

    def setMeshUnits(self, x, y, z=""):
        """
        Sets units for the mesh axes. These are currently only used by the
        Silo exporter.

        :param x: unit for the x-axis
        :type x: ``str``
        :param y: unit for the y-axis
        :type y: ``str``
        :param z: unit for the z-axis (optional for 2D)
        :type z: ``str``
        """
        self._meshunits=[x,y,z]

    def setMetadataSchemaString(self, schema, metadata=""):
        """
        Sets metadata namespaces and the corresponding metadata.
        Only used for the VTK file format at the moment.

        :param schema: A dictionary that maps namespace prefixes to namespace
                       names, e.g. {'gml':'http://www.opengis.net/gml'}
        :param metadata: The actual metadata string which will be enclosed in
                         '<MetaData>' tags.
        """
        self._metadata=metadata
        ss=""
        for i,p in sorted(list(schema.items()), key=lambda x: x[0]):
            ss="%s xmlns:%s=\"%s\""%(ss, i, p)
        self._md_schema=ss.strip()

    def export(self):
        """
        Executes the actual data export. Depending on the formats parameter
        used in the constructor all data added by addData() is written to disk
        (RESTART,SILO,VTK) or made available through the VisIt simulation
        interface (VISIT).
        """

        if self._domain is None:
            print("Warning: DataManager.export() called but no domain set!")
            return

        self._N += 1
        ds = None
        nameprefix=os.path.join(self._workdir, "dataset.%04d"%(self._N))

        for f in self._exportformats:
            if f == self.SILO:
                if ds is None:
                    ds=self.__createDataset()
                ds.saveSilo(nameprefix)
            elif f == self.VTK:
                if ds is None:
                    ds=self.__createDataset()
                ds.saveVTK(nameprefix)
            elif f == self.VISIT:
                from esys.weipa.weipacpp import visitPublishData
                if ds is None:
                    ds=self.__createDataset()
                visitPublishData(ds)
            elif f == self.RESTART:
                # only write checkpoint files with the requested frequency
                if self._N % self._checkpointfreq==0:
                    self.__saveState()
            else:
                raise ValueError("export: Unknown export format "+str(f))

        self.__clearData()

    def __createDataset(self):
        from esys.weipa.weipacpp import EscriptDataset
        from esys.weipa import createDataset

        ds = createDataset(self._domain, **self._data)
        ds.setCycleAndTime(self._N, self._time)
        ds.setMetadataSchemaString(self._md_schema, self._metadata)
        ds.setMeshLabels(self._meshlabels[0], self._meshlabels[1], self._meshlabels[2])
        ds.setMeshUnits(self._meshunits[0], self._meshunits[1], self._meshunits[2])
        return ds
    
    def __clearData(self):
        #print("Clearing all data")
        self._restartdir = None
        self._domain = None
        self._stamp = {}
        self._data = {}

    def __getStampFilename(self, dir_name):
        return os.path.join(self._workdir, dir_name, "stamp.%d"%self._myrank)

    def __getDumpFilename(self, data_name, dir_name):
        return os.path.join(self._workdir, dir_name, "%s.h5"%data_name)

    def __initVisit(self, simFile, comment=""):
        """
        Initialises the VisIt interface if available.

        :param simFile: Name of the sim file to be generated which can be
                        loaded into a VisIt client
        :param comment: A short description of this simulation
        """
        from esys.weipa.weipacpp import visitInitialize
        return visitInitialize(simFile, comment)

    def __loadState(self):
        stamp_file=self.__getStampFilename(self._restartdir)
        try:
            self._stamp = pickle.load(open(stamp_file, "rb"))
            self._N = int(self._restartdir[len(self._restartprefix)+1:])
        except:
            raise IOError("Could not load stamp file "+stamp_file)
        # load domain
        ff=self.__getDumpFilename("_domain",self._restartdir)
        modname=self._stamp['__domainmodule']
        clsname=self._stamp['__domainclass']
        try:
            domclass=__import__(modname, fromlist=[clsname])
            self._domain = domclass.LoadMesh(ff)
        except:
            raise ImportError("Unable to load %s using %s.%s!"%(ff, modname, clsname))

    def __saveState(self):
        restartdir = "%s_%04d"%(self._restartprefix, self._N)
        util.mkDir(os.path.join(self._workdir, restartdir))
        stamp_file=self.__getStampFilename(restartdir)
        self._stamp['__domainmodule']=self._domain.__module__
        self._stamp['__domainclass']=type(self._domain).__name__
        pickle.dump(self._stamp, open(stamp_file, "wb"))
        ff=self.__getDumpFilename("_domain", restartdir)
        self._domain.dump(ff)
        for name, var in sorted(list(self._data.items()), key=lambda x: x[0]):
            ff=self.__getDumpFilename(name, restartdir)
            var.dump(ff)
        print(("Restart files saved in "+os.path.join(self._workdir, restartdir)))
        # keep only one restart directory
        old_restartdir = "%s_%04d"%(self._restartprefix, self._N-self._checkpointfreq)
        self.__removeDirectory(os.path.join(self._workdir, old_restartdir))

    def __removeDirectory(self, path):
        if self._myrank==0 and os.path.isdir(path):
            shutil.rmtree(path, True)
            #print("Removed restart directory %s."%path)
        esc.MPIBarrierWorld()

