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

# $Id: item.py,v 1.12 2005/02/24 01:18:23 paultcochrane Exp $

## @file item.py

"""
This is the file of the base class for items within a scene
"""

# generic imports
from pyvisi.renderers.vtk.common import debugMsg
from pyvisi.item import Item as BaseItem

__revision__ = '$Revision: 1.12 $'

class Item(BaseItem):
    """
    This is the base class for items within a scene
    """

    def __init__(self):
        """
        Initialisation of Item class
        """
        debugMsg("Called Item.__init__()")
        BaseItem.__init__(self)

# vim: expandtab shiftwidth=4:
