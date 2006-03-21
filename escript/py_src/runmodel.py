#!/usr/bin/python
# $Id$

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__licence__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licences/osl-3.0.php"""

from esys.escript import modelframe

#commandline utility to take an xml file, parse it, and run a simulation. 
# invoke this by doing ./runmodel.py <filename.xml>

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
    
