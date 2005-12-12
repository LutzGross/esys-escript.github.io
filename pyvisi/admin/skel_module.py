#!/usr/bin/env python

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

# $Id: skel_module.py,v 1.8 2005/02/24 04:23:12 paultcochrane Exp $

import string,sys,os

if (len(sys.argv) != 2):
    print "Usage: python skel_module.py <moduleName>"
    sys.exit(1)

baseClassNames = [ 
                    'item',
                    'scene',
                    'renderer',
                    'axes',
                    'camera',
                    'image',
                    'plane',
                    'plot',
                    'text',
                    ]

imageClassNames = [ 
                    'JpegImage',
                    'PbmImage',
                    'PdfImage',
                    'PngImage',
                    'PnmImage',
                    'PsImage',
                    'TiffImage',
                    ]

plotClassNames = [
                    'ArrowPlot',
                    'ContourPlot',
                    'LinePlot',
                    'SurfacePlot',
                    ]

# this is the name of the module, and the name of the directory to create
moduleName = sys.argv[1]

def createDirs():
    # create the file structure (module directory plus tests subdirectory)
    os.mkdir(moduleName)
    testDirName = "%s/tests" % moduleName
    os.mkdir(testDirName)

# the copyright string to put at the top of the file
copyrightStr = """# Copyright (C) 2004-2005 Paul Cochrane
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

# $"""

# doing this to get around cvs keyword expansion issues
copyrightStr += "Id"
copyrightStr += "$\n\n"

# this gets vim to "do the right thing" wrt spaces and tabs etc
vimStr = "# vim: expandtab shiftwidth=4:"

def writeInitFile():
    """
    Generate the __init__.py file
    """

    initText = "## @file __init__.py"
    initText += """

\"\"\"
Initialisation of the %s renderer module
\"\"\"
    """ % moduleName

    initText += """
from pyvisi.renderers.%s.common \\
        import _rendererName, _rendererVersion, _rendererRevision
print \"This is the \\\"%%s\\\" renderer module version %%s-%%s\" %% \\
        (_rendererName, _rendererVersion, _rendererRevision)

__author__ = 'Insert author name here'
__version__ = _rendererVersion
__revision__ = _rendererRevision

""" % moduleName

    for baseClass in baseClassNames:
        if baseClass == 'plot':
            initText += "from pyvisi.renderers.%s.%s import %s, \\\n" % \
                    (moduleName, baseClass.lower(), baseClass.capitalize())
            initText += "        "
            for i in range(len(plotClassNames)-1):
                initText += "%s, " % plotClassNames[i]
            initText += "%s\n" % plotClassNames[-1]
        elif baseClass == 'image':
            initText += "from pyvisi.renderers.%s.%s import %s, \\\n" % \
                    (moduleName, baseClass.lower(), baseClass.capitalize())
            initText += "        "
            for i in range(len(imageClassNames)-1):
                initText += "%s, " % imageClassNames[i]
            initText += "%s\n" % imageClassNames[-1]
        else:
            initText += "from pyvisi.renderers.%s.%s import %s\n" % \
                    (moduleName, baseClass.lower(), baseClass.capitalize())
        
    fname = moduleName + '/__init__.py'
    f = open(fname, 'w')
    f.write(copyrightStr)
    f.write(initText + "\n")
    f.write(vimStr + "\n")
    f.close()

def writeCommonFile():
    """
    Generate the common.py file
    """

    commonText = "## @file common.py"
    commonText += """

\"\"\"
Variables common to all classes and functions
\"\"\"

_debug = 1
_rendererName = \'%s\'
_rendererVersion = \'Enter renderer version number\'
_rendererRevision = \'Enter renderer revision text\'

__revision__ = \'$
"""
    # this bit here is to get around cvs keyword expansion issues
    commonText += "Revision"
    commonText += """$\'

    """ % moduleName.upper()

    commonText += """
def debugMsg(message):
    \"\"\"
    Convenience function for debugging messages.

    This function will print out a debugging message if the debug variable
    is set.

    @param message: the message to output if the debug flag is set
    @type message: string
    \"\"\"
    if _debug:
        print \"\\t%s: %s\" % (_rendererName, message)

def unsupportedError():
    \"\"\"
    Print an error message when a method is called that is defined in pyvisi
    but is not supported at the renderer module level.
    \"\"\"
    errorString = \"Sorry, but %s doesn't support this method.\" % _rendererName
    raise NotImplementedError, errorString

def getRevision():
    \"\"\"
    Get the revision string/number
    \"\"\"
    return _rendererRevision
    """

    fname = moduleName + '/common.py'
    f = open(fname, 'w')
    f.write(copyrightStr)
    f.write(commonText + "\n")
    f.write(vimStr + "\n")
    f.close()

def getClassFileHeader(className, baseClassName):
    # epydoc documentation strings for info about the file itself
    fileDoxStr = "## @file %s.py\n" % className.lower()

    # the header of the file
    fileHeaderStr = """
\"\"\"
Brief introduction to what the file contains/does
\"\"\"

from pyvisi.renderers.%s.common import debugMsg, overrideWarning

from pyvisi.renderers.%s.%s import %s

__revision__ = '$
"""
    # this is here to get around cvs keyword expansion problems
    fileHeaderStr += "Revision"
    fileHeaderStr += """$'

""" % (moduleName, moduleName, baseClassName.lower(), baseClassName)
    return copyrightStr + fileDoxStr + fileHeaderStr

def getBaseClassFileHeader(className, baseClassName):
    # epydoc documentation strings for info about the file itself
    fileDoxStr = "## @file %s.py\n" % className.lower()

    # the header of the file
    fileHeaderStr = """
\"\"\"
Brief introduction to what the file contains/does
\"\"\"

from pyvisi.renderers.%s.common import debugMsg, overrideWarning, getRevision

from pyvisi.%s import %s as %s
""" % (moduleName, className.lower(), className.capitalize(), baseClassName)

    if className.lower() == 'scene':
        fileHeaderStr += """
from pyvisi.renderers.%s.renderer import Renderer
""" % moduleName

    fileHeaderStr += """

__revision__ = getRevision()

""" 
    return copyrightStr + fileDoxStr + fileHeaderStr
   
def getClassFileBody(className, baseClassName):
    """
    Gets the text of the class body skeleton

    @param className: the name of the class body skeleton
    @type className: string

    @param baseClassName: the name of the class' base class
    @type baseClassName: string
    """

    # the class body skeleton, with epydoc strings added, and a dummy function
    classStr = """
class %s(%s):
""" % (className, baseClassName)
    classStr += """
    \"\"\"
    Brief introduction to what the class does
    \"\"\"

    def __init__(self, arg):
        \"\"\"
        Brief description of the init function

        @param arg: a description of the argument
        @type arg: the type of the argument
        \"\"\"
        debugMsg(\"Called %s.__init__()\")
        %s.__init__(self)  # initialisation of base class
    
    def myfunc(self, myarg):
        \"\"\"
        Brief description of what the function does
    
        Replace the text given here with an actual description of what
        the function does.  Also change the name of the function and
        the name of the argument.
    
        @param myarg: Description of what the parameter means/does
        @type myarg: the type of the argument
        \"\"\"
        return

    """ % (baseClassName, className)
    return classStr

def writeClassFile(classFileText, className):
    """
    Write the text of the generated class skeleton out to file
    """

    fname = moduleName + '/' + className + '.py'
    f = open(fname, 'w')
    f.write(classFileText + "\n")
    f.write(vimStr + "\n")
    f.close()

# now biff out the stuff for the test file
def getTestFileHeader(className):
    """
    Get the text of the generated header for the test file
    """

    fileDoxStr = '## @file test_' + className.lower() + '.py'

    headerTextStr = """

import unittest
import sys,os,string
here = os.getcwd() + '/../../'
sys.path.append(here)
# this should import all of the pyvisi stuff needed
from pyvisi import *   
# this should import the renderer specific stuff
from pyvisi.renderers.%s import * 
    """ % moduleName
    return copyrightStr + fileDoxStr + headerTextStr

def getTestFileBody(className):
    """
    Get the text of the generated body for the test class
    """

    bodyTextStr = """
\"\"\"
Class and functions for testing the %s class
\"\"\"

class Test%s(unittest.TestCase):
    \"\"\"
    The main test class
    \"\"\"

    def testFunction(self):
        \"\"\"
        A test function
        \"\"\"
        self.assertEqual()
        self.assertRaises()
        self.assert_()

if __name__ == '__main__':
    unittest.main()

    """ % (className, className)
    return bodyTextStr

def writeTestFile(testFileText, className):
    """
    Write the generated test skeleton to file
    """

    fname = moduleName + '/tests/' + 'test_' + className + '.py'
    f = open(fname, 'w')
    f.write(testFileText + "\n")
    f.write(vimStr + "\n")
    f.close()


# right, now do stuff
createDirs()
writeInitFile()
writeCommonFile()
for name in baseClassNames:
    if name == 'item' or name == 'scene' or name == 'renderer':
        # write the class file
        baseName = "Base" + name.capitalize()
        classFileText = getBaseClassFileHeader(name.capitalize(), baseName)
        classFileText += getClassFileBody(name.capitalize(), baseName)
        writeClassFile(classFileText, name)
        # write the test file
        testFileText = getTestFileHeader(name.capitalize())
        testFileText += getTestFileBody(name.capitalize())
        writeTestFile(testFileText, name)
    elif name == 'axes':
        # write the class file
        baseName = 'Plot'
        classFileText = getClassFileHeader(name.capitalize(), baseName)
        classFileText += getClassFileBody(name.capitalize(), baseName)
        writeClassFile(classFileText, name)
        # write the test file
        testFileText = getTestFileHeader(name.capitalize())
        testFileText += getTestFileBody(name.capitalize())
        writeTestFile(testFileText, name)
    else:
        # write the class file
        baseName = 'Item'
        classFileText = getClassFileHeader(name.capitalize(), baseName)
        classFileText += getClassFileBody(name.capitalize(), baseName)
        testFileText = getTestFileHeader(name.capitalize())
        testFileText += getTestFileBody(name.capitalize())
        if name == 'plot':
            for subClass in plotClassNames:
                classFileText += getClassFileBody(subClass, name.capitalize())
                testFileText += getTestFileBody(subClass)
        elif name == 'image':
            for subClass in imageClassNames:
                classFileText += getClassFileBody(subClass, name.capitalize())
                testFileText += getTestFileBody(subClass)
        writeClassFile(classFileText, name)
        # write the test file
        writeTestFile(testFileText, name)

# vim: expandtab shiftwidth=4:
