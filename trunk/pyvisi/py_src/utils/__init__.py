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

# $Id: __init__.py,v 1.1 2005/06/24 00:55:10 paultcochrane Exp $

## @file __init__.py

"""
Initialisation of the pyvisi utilities
"""

from pyvisi.common import _pyvisiVersion, _pyvisiRevision
print "This is PyVisi version %s-%s" % (_pyvisiVersion, _pyvisiRevision)

__author__ = 'Paul Cochrane'
__version__ = _pyvisiVersion
__revision__ = _pyvisiRevision

# vim: expandtab shiftwidth=4:
