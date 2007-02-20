"""
set up of mining activities in modelframe
                                                                                                                                                                                                     
@var __author__: name of author
@var __licence__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Lutz Gross, l.gross@uq.edu.au"
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision$"
__date__="$Date$"

from esys.escript.modelframe import Model
from mines import parse

class MiningHistory(Model):
    """
    manages the history of mines. It mandels the time steping according to the available 
    data and a dictionary of mass changes per year for all the mines. 

    @ivar history: mine site history file.
    @type history: C{DataSource}
    @ivar mass_changes: current mass change per year.
    @type mass_changes: C{DataSource}

    """
    def __init__(self,**kwargs):
        """
        """
        super(MiningHistory,self).__init__(**kwargs)
        self.declareParameter(t=0.,
                              history=None,
                              mass_changes=None)

    def doInitialization(self):
        """
        initialize time integration
        """
        self.__minesite=parse(open(self.history.getLocalFileName(),'r'))
        self.mass_changes=self.__minesite.getMassChanges(self.t)

    def doStepPreprocessing(self, dt):
        self.mass_changes=self.__minesite.getMassChanges(self.t)
        print self.t,":",self.mass_changes

    def getSafeTimeStepSize(self, dt):
        print "new time marker:",self.t,":",self.__minesite.getNextTimeMarker(self.t)
        return self.__minesite.getNextTimeMarker(self.t)-self.t



