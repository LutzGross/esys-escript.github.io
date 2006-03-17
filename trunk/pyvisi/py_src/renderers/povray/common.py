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

# $Id: common.py,v 1.5 2005/06/24 00:48:17 paultcochrane Exp $

## @file common.py

"""
Variables common to all classes and functions
"""

from pyvisi.common import _debug
_rendererName = "POVRAY"
_rendererVersion = '0.1'
_rendererRevision = 'pre-alpha-1'

__revision__ = '$Revision: 1.5 $'

def debugMsg(message):
    """
    Convenience function for debugging messages.

    This function will print out a debugging message if the debug variable
    is set.

    @param message: the message to output if the debug flag is set
    @type message: string
    """
    if _debug:
        print "\t%s: %s" % (_rendererName, message)

def unsupportedError():
    """
    Print an error message when a method is called that is defined in pyvisi
    but is not supported at the renderer module level.
    """
    errorString = "Sorry, but %s doesn't support this method." % _rendererName
    raise NotImplementedError, errorString

# vim: expandtab shiftwidth=4:
