# -*- coding: utf-8 -*-
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

"""
a simple restart manager (still under testing)

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

import util
import cPickle
import os
from esys.escript import getMPIRankWorld, MPIBarrierWorld, load

class RestartManager(object):
   """
   this is a simple restart manager

   usage
    
   rm=RestartManager(prefix="myrestart")
   if rm.hasData():
        dom=rm.getDomain()
        t=rm.getStamp("t")
        dt=rm.getStamp("stemp")
        T=rm.getValue("T")
        u=rm.getValues("d")
   else:
        T=...
        u=...
     
   rm.setStamp(t=t,step=dt) # this creates a new restart directory
   rm.checkIn(T=T,d=u)      # adds values to the restart data base
   """
   def __init__(self, new=False, prefix="restart", work_dir=".", myrank=getMPIRankWorld()):
       """
       initializes the restart manager. If new is True a new restart directory is created
       otherwise the last available directory is used for restart.
       
       :param new: fresh start. Otherwise the last available restart file will be used.
       :type new: `bool`
       :param prefix: directory name prefix
       :type prefix: `str`
       :param work_dir: directory where restart files are created
       :type prefix: `str`
       """
       self.PREFIX=prefix
       self.MYRANK=myrank
       self.WORK_DIR=work_dir
       self.NEW_RESTART_DIR=None
       util.mkDir(self.WORK_DIR)
       
       restart_files=[]
       
       for f in os.listdir(self.WORK_DIR):
           if f.startswith(self.PREFIX): restart_files.append(f)
       restart_files.sort()
           
       if len(restart_files)==0 or new:
          for f in restart_files: self.__removeRestartDirectory(f)
          self.RESTART_DIR=None
          self.N=-1
          self.__STAMP=None
          
       else:
          self.RESTART_DIR=restart_files[-1]
          for f in restart_files[:-1]: self.__removeRestartDirectory(f)
          print "Restart from file %s."%self.RESTART_DIR
          self.N=int(self.RESTART_DIR[len(self.PREFIX)+1:])
          self.__STAMP=cPickle.load(open(self.__getStampFileName(self.RESTART_DIR),"rb"))
       self.DOMAIN = None
       
   def __removeRestartDirectory(self, dir_name):
       if self.MYRANK==0 and os.path.isdir(dir_name):
           for root, dirs, files in os.walk(dir_name, topdown=False):
               for name in files: os.remove(os.path.join(root,name))
               for name in dirs: os.remove(os.path.join(root,name))
           os.rmdir(dir_name)
           print "Restart files %s have been removed."%dir_name
       MPIBarrierWorld()

   def __getStampFileName(self, dir_name):
         return os.path.join(self.WORK_DIR,dir_name,"stamp.%d"%self.MYRANK)
   def __getValueFileName(self, value_name, dir_name):
         return os.path.join(self.WORK_DIR,dir_name,"%s.nc"%value_name)

   def hasData(self):
       """
       returns True if manager holds data for restart
       """
       return not self.RESTART_DIR == None

   def setStamp(self, **values):
       """
       sets a time stap defined by a set of pickable values (e.g. float, int, str)
       A new restart directory is created and the directory before the last is removed.
       """
       # complete last saver:
       if not self.NEW_RESTART_DIR == None: 
          if self.RESTART_DIR!= None: self.__removeRestartDirectory(self.RESTART_DIR)
          self.RESTART_DIR=self.NEW_RESTART_DIR
          self.DOMAIN=self.NEW_DOMAIN
          self.N+=1
       # create a new restart directory
       self.NEW_RESTART_DIR="%s_%d"%(self.PREFIX,self.N+1)
       self.NEW_DOMAIN = None
       util.mkDir(self.NEW_RESTART_DIR)
       cPickle.dump(values,open(self.__getStampFileName(self.NEW_RESTART_DIR),"wb"))

   def getStamp(self,value_name):
       """
       returns the value of stamp paramater with name value_name
       
       :param value_name: requested value
       :type value_name: `str`
       :return: requested value from latest restart file.
       """
       if self.__STAMP == None:
          raise ValueError,"no restart data available."
       if not self.__STAMP.has_key(value_name):
           raise ValueError,"no stamp parameter %s."%value_name
       return self.__STAMP[value_name]

   def checkIn(self,**values):
       """
       writes 'escript.Data' objects to restart dierctory.
       """
       if self.NEW_RESTART_DIR == None:
           raise ValueError,"No time stamp set."
       for i,v in values.items():
         self.checkinDomain(v.getDomain())
         ff=self.__getValueFileName(i, self.NEW_RESTART_DIR)
         v.dump(ff)
         print "%s dumped to file %s."%(i,ff)

   def getValue(self,value_name):
       """
       returns the value of value_name in the latest restart file
       
       :param value_name: requested value
       :type value_name: `str`
       :return: requested value from latest restart file.
       """
       if self.RESTART_DIR == None:
          raise ValueError,"no restart data available."
       domain=self.getDomain()
       if domain == None:
	  raise ValueError,"no domain found for restart."
       ff=self.__getValueFileName(value_name, self.RESTART_DIR)
       v=load(ff, domain)
       print "Value %s recovered from %s."%(value_name, ff)
       return v 
   def checkinDomain(self, domain):
       """
       writes the domain to current restart directory
       """
       if self.NEW_RESTART_DIR == None:
           raise ValueError,"No time stamp set."
       if self.NEW_DOMAIN == None:
           ff=self.__getValueFileName("__mesh_finley",self.NEW_RESTART_DIR) # finley extension added to mark domain type
           domain.dump(ff)
           print "Domain dumped to file %s."%(ff,)
           self.NEW_DOMAIN=domain

   def getDomain(self):
       """
       recovers the domain from the current restart directory
       """
       if self.RESTART_DIR == None:
          raise ValueError,"no restart data available."
       if self.DOMAIN == None:
	  for f in os.listdir(os.path.join(self.WORK_DIR, self.RESTART_DIR)):
               if f.startswith("__mesh_finley"): 
                 ff=self.__getValueFileName("__mesh_finley",self.RESTART_DIR)
                 import esys.finley
                 self.DOMAIN=esys.finley.LoadMesh(ff)
                 print "Domain recovered from file %s."%ff
       return self.DOMAIN
