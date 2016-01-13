#!/usr/bin/python
# $Id$

"""
commandline utility to take an xml file, parse it, and run a simulation. 
invoke this by doing ./runmodel.py <filename.xml>

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Elspeth Thorne, e.thorne@uq.edu.au"
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"


from esys.escript import modelframe
import optparse

parser = optparse.OptionParser(usage="\n%prog [options]\n%prog files...")
parser.add_option('-f', '--file', dest='filename', 
        help='the FILE', metavar='FILE')
parser.add_option('-n', '--old-name', action="store",
        help='the old filename, used in aup',
        dest='old_name', default='')
def main():
    (options, args) = parser.parse_args()
    if options.filename and args:
        parser.usage("Please only specifiy 1 file if using --file=")
    if options.filename:
        files = [(file(options.filename), options.filename)]
    elif args:
        files = [(file(arg), arg) for arg in args]
    else:
        parser.usage

    for f, filename in files:
        
        simstring = f.read()
        sim = modelframe.parse(simstring)
        print sim
        sim.run()

if __name__=='__main__':
    main()
    
