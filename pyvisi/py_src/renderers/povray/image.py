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

# $Id: image.py,v 1.5 2005/02/24 04:34:07 paultcochrane Exp $

## @file image.py

"""
Class and functions associated with a pyvisi Image objects
"""

# generic imports
from pyvisi.renderers.povray.common import debugMsg
from pyvisi.common import fileCheck

# module specific imports
from pyvisi.renderers.povray.item import Item

__revision__ = '$Revision: 1.5 $'

class Image(Item):
    """
    Image class.  Generic class to handle image data.
    """
    def __init__(self, format, scene):
        """
        Initialises the Image class object
        
        @param format: The image format
        @type format: string

        @param scene: The Scene object to add to
        @type scene: Scene object
        """
        debugMsg("Called Image.__init__()")
        Item.__init__()

        #if format == "jpeg":
            #if _debug: print "\t%s: Using jpeg image format" % rendererName
            #return JpegImage(scene)
        #else:
            #print "Unknown image format %s" % format
            #return None

        if format is None:
            raise ValueError, "You must specify an image format"

        if scene is None:
            raise ValueError, "You must specify a scene object"

    def load(self, fname):
        """
        Loads image data from file.

        @param fname: The filename from which to load image data
        @type fname: string
        """
        debugMsg("Called Image.load()")

        fileCheck(fname)

        return

class JpegImage(Image):
    """
    Subclass of Image class to explicitly handle jpeg images
    """
    def __init__(self, scene):
        """
        Initialises the JpegImage class object

        @param scene: The Scene object to add to
        @type scene: Scene object
        """
        debugMsg("Called JpegImage.__init__()")
        Image.__init__()
        self.renderer = scene.renderer

    def load(self, fname):
        """
        Loads jpeg image data from file.

        @param fname: The filename from which to load jpeg image data
        @type fname: string
        """
        debugMsg("Called JpegImage.load()")

        fileCheck(fname)
        
        return

    def render(self):
        """
        Does JpegImage object specific (pre)rendering stuff
        """
        debugMsg("Called JpegImage.render()")

        return


# vim: expandtab shiftwidth=4:
