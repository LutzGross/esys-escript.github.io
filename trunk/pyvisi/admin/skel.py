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

# $Id: skel.py,v 1.3 2005/02/08 05:54:17 paultcochrane Exp $

import string,sys

if (len(sys.argv) != 2):
    print "Usage: python skel.py <className>"
    sys.exit(1)

# this is the filename (minus the .py) and the name of the class
classname = sys.argv[1]  

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

# $
""" % classname
# this is here to get around cvs keyword expansion issues
copyrightStr += "Id"
copyrightStr += "$"

# doxygen documentation strings for info about the file itself
fileDoxStr = "## @file %s.py\n" % classname
fileDoxStr += """
\"\"\"
Brief introduction to what the file contains/does
\"\"\"

from pyvisi.common import debugMsg, overrideWarning

"""

# the class body skeleton, with doxygen strings added, and a dummy function
classStr = "class %s():" % classname.capitalize()
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
        BaseClass.__init__(self)
        debugMsg(\"Called %s.__init__()\")
    
    def myfunc(myarg):
        \"\"\"
        Brief description of what the function does
    
        Replace the text given here with an actual description of what
        the function does.  Also change the name of the function and
        the name of the argument.
    
        @param myarg: Description of what the parameter means/does
        @type myarg: the type of the argument
        \"\"\"
        return
        """ % classname.capitalize()

# this gets vim to "do the right thing" wrt spaces and tabs etc
vimStr = "# vim: expandtab shiftwidth=4:"

fname = classname + '.py'
f = open(fname, 'w')
f.write(copyrightStr + "\n")
f.write(fileDoxStr + "\n")
f.write(classStr + "\n")
f.write(vimStr + "\n")
f.close()

# now biff out the stuff for the test file

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

# \$Id\$

""" % classname

fileDoxStr = '## @file test_' + classname + '.py'

bodyTextStr = """import unittest
import sys,os,string
here = os.getcwd() + '/../../'
sys.path.append(here)
from pyvisi import *   # this should import all of the pyvisi stuff needed

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

""" % (classname.capitalize(), classname.capitalize())

fname = 'test_' + classname + '.py'
f = open(fname, 'w')
f.write(copyrightStr + "\n")
f.write(fileDoxStr + "\n")
f.write(bodyTextStr + "\n")
f.write(vimStr + "\n")
f.close()

# vim: expandtab shiftwidth=4:

