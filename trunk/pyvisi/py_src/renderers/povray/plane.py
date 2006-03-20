# Copyright (C) 2004-2005 Paul Cochrane
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

# $Id: plane.py,v 1.5 2005/02/24 04:34:07 paultcochrane Exp $

## @file plane.py

"""
The classes associated with Planes
"""

# generic imports
from common import debugMsg

# module specific imports
from item import Item

__revision__ = '$Revision: 1.5 $'

class Plane(Item):
    """
    Generic class for Plane objects
    """

    def __init__(self, scene):
        """
        Initialisation of the Plane object
        """
        debugMsg("Called Plane.__init__()")
        Item.__init__()

        self.renderer = scene.renderer

    def mapImageToPlane(self, image):
        """
        Maps an Image object onto a Plane object
        """
        debugMsg("Called Plane.mapImageToPlane()")

        if image is None:
            raise ValueError, "You must specify an image object"

        return

    def render(self):
        """
        Perform Plane object specific (pre)rendering tasks
        """
        debugMsg("Called Plane.mapImageToPlane()")

        return
    
# vim: expandtab shiftwidth=4:
