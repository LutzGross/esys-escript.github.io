########################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"
__author__="Lutz Gross, Cihan Altinay"

"""
an escript data import and export manager (still under development)

:var __author__: name of authors
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point to documentation
"""

import cPickle
import os
import shutil
import util
from esys.escript import getMPIRankWorld, MPIBarrierWorld, load

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

    RESTART, SILO, VISIT, VTK = range(4)

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
        self._stamp={}
        self._time=0.
        self._restartdir=None
        self._N=-1
        self._myrank=getMPIRankWorld()
        self._exportformats=set(formats)
        self._restartprefix=restart_prefix
        self._workdir=work_dir
        util.mkDir(self._workdir)
        if self.VISIT in self._exportformats:
            simFile=os.path.join(self._workdir, "escriptsim.sim2")
            if not self.__initVisit(simFile, "Escript simulation"):
                print "Warning: Could not initialize VisIt interface"
                self._exportformats.remove(self.VISIT)
        # find all restart directories
        restart_folders = []
        for f in os.listdir(self._workdir):
            if f.startswith(self._restartprefix):
                restart_folders.append(f)
        # remove unneeded restart directories
        if len(restart_folders)>0:
            restart_folders.sort()
            if do_restart:
                self._restartdir=restart_folders[-1]
                print "Restart from "+os.path.join(self._workdir, self._restartdir)
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

        for name,var in data.items():
            if hasattr(var, "getDomain"):
                if self._domain==None:
                    self._domain=var.getDomain()
                elif self._domain != var.getDomain():
                    raise ValueError, "addData: Data must be on the same domain!"
                self._data[name]=var
            else:
                self._stamp[name]=var

    def hasData(self):
        """
        Returns True if the manager holds data for restart
        """
        return self._restartdir != None

    def getDomain(self):
        """
        Returns the domain as recovered from restart files.
        """
        if not self.hasData():
            raise ValueError, "No restart data available"
        return self._domain 

    def getValue(self, value_name):
        """
        Returns an 'escript.Data' object or other value that has been loaded
        from restart files.
        """
        if not self.hasData():
            raise ValueError, "No restart data available"

        if value_name in self._stamp:
            return self._stamp[value_name]

        ff=self.__getDumpFilename(value_name, self._restartdir)
        var = load(ff, self._domain)
        #print "Value %s recovered from %s."%(value_name, ff)
        return var 

    def getCycle(self):
        """
        Returns the export cycle (=number of times export() has been called)
        """
        return self._N

    def setTime(self, time):
        """
        Sets the simulation timestamp.
        """
        self._time = time

    def setMetadataSchemaString(self, schema, metadata=""):
        """
        Sets metadata namespaces and the corresponding metadata.
        Only used for the VTK file format at the moment.

        :param schema: A dictionary that maps namespace prefixes to namespace
                       names, e.g. {'gml':'http://www.opengis.net/gml'}
        :param metadata: The actual metadata string which will be enclosed in
                         '<MetaData>' tags.
        """
        self._metadata="<MetaData>"+metadata+"</MetaData>"
        ss=""
        for i,p in schema.items():
            ss="%s xmlns:%s=\"%s\""%(ss, i, p)
        self._md_schema=ss.strip()

    def export(self):
        """
        Executes the actual data export. Depending on the formats parameter
        used in the constructor all data added by addData() is written to disk
        (RESTART,SILO,VTK) or made available through the VisIt simulation
        interface (VISIT).
        """

        if self._domain == None:
            return

        idata = self.__interpolateData()

        self._N += 1
        nameprefix=os.path.join(self._workdir, "dataset%04d"%(self._N))

        for f in self._exportformats:
            if f == self.SILO:
                from esys.weipa.weipacpp import _saveSilo
                filename=nameprefix+".silo"
                _saveSilo(filename, self._N, self._time, self._domain, idata)
            elif f == self.VTK:
                from esys.weipa.weipacpp import _saveVTK
                filename=nameprefix+".vtu"
                _saveVTK(filename, self._N, self._time, self._domain,
                        idata, self._metadata, self._md_schema)
            elif f == self.VISIT:
                from esys.weipa.weipacpp import _visitPublishData
                _visitPublishData(self._N, self._time, self._domain, idata)
            elif f == self.RESTART:
                self.__saveState()
            else:
                raise ValueError, "export: Unknown export format "+f

        self.__clearData()

    def __clearData(self):
        #print "Clearing all data"
        self._restartdir = None
        self._domain = None
        self._stamp = {}
        self._data = {}

    def __getStampFilename(self, dir_name):
        return os.path.join(self._workdir, dir_name, "stamp.%d"%self._myrank)

    def __getDumpFilename(self, data_name, dir_name):
        return os.path.join(self._workdir, dir_name, "%s.nc"%data_name)

    def __initVisit(self, simFile, comment=""):
        """
        Initialises the VisIt interface if available.

        :param simFile: Name of the sim file to be generated which can be
                        loaded into a VisIt client
        :param comment: A short description of this simulation
        """
        from esys.weipa.weipacpp import _visitInitialize
        return _visitInitialize(simFile, comment)

    def __interpolateData(self):
        # (Reduced)Solution is not directly supported so interpolate to
        # different function space
        from esys.escript import Solution, ReducedSolution
        from esys.escript import ContinuousFunction, ReducedContinuousFunction
        from esys.escript.util import interpolate
        new_data={}
        for n,d in self._data.items():
            if not d.isEmpty():
                fs=d.getFunctionSpace()
                domain=fs.getDomain()
                if fs == Solution(domain):
                    new_data[n]=interpolate(d, ContinuousFunction(domain))
                elif fs == ReducedSolution(domain):
                    new_data[n]=interpolate(d, ReducedContinuousFunction(domain))
                else:
                    new_data[n]=d
        return new_data

    def __loadState(self):
        stamp_file=self.__getStampFilename(self._restartdir)
        try:
            self._stamp = cPickle.load(open(stamp_file, "rb"))
            self._N = int(self._restartdir[len(self._restartprefix)+1:])
        except:
            raise IOError, "Could not load stamp file "+stamp_file
        # load domain
        path=os.path.join(self._workdir, self._restartdir)
        for f in os.listdir(path):
            if f.startswith("__mesh_finley"): 
                ff=self.__getDumpFilename("__mesh_finley",self._restartdir)
                import esys.finley
                self._domain = esys.finley.LoadMesh(ff)
                #print "Domain recovered from file %s."%ff
                break

    def __saveState(self):
        restartdir = "%s_%04d"%(self._restartprefix, self._N)
        util.mkDir(os.path.join(self._workdir, restartdir))
        stamp_file=self.__getStampFilename(restartdir)
        cPickle.dump(self._stamp, open(stamp_file, "wb"))
        ff=self.__getDumpFilename("__mesh_finley", restartdir)
        self._domain.dump(ff)
        for name, var in self._data.items():
            ff=self.__getDumpFilename(name, restartdir)
            var.dump(ff)
        print "Restart files saved in "+os.path.join(self._workdir, self._restartdir)
        # keep only one restart directory
        old_restartdir = "%s_%04d"%(self._restartprefix, self._N-1)
        self.__removeDirectory(os.path.join(self._workdir, old_restartdir))

    def __removeDirectory(self, path):
        if self._myrank==0 and os.path.isdir(path):
            shutil.rmtree(path, True)
            #print "Removed restart directory %s."%path
        MPIBarrierWorld()

