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

# $Id: item.py,v 1.10 2005/03/07 07:45:48 paultcochrane Exp $

## @file item.py

"""
This is the file of the base class for items within a scene
"""

from pyvisi.common import debugMsg, overrideWarning

__revision__ = '$Revision: 1.10 $'

class Item(object):
    """
    This is the base class for items within a scene
    """

    def __init__(self):
        """
        Initialisation
        """
        object.__init__(self)
        debugMsg("Called Item.__init__()")

    def render(self):
        """
        Render the object
        """
        debugMsg("Called Item.render()")

        # print a warning if get to here
        overrideWarning("Item.render")

        return

# vim: expandtab shiftwidth=4:
