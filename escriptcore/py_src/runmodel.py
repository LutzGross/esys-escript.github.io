
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
commandline utility to take an xml file, parse it, and run a simulation.
invoke this by doing ./runmodel.py <filename.xml>

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Elspeth Thorne, e.thorne@uq.edu.au"


from . import modelframe
import optparse

parser = optparse.OptionParser(usage="%prog [options] <ESySXML files>")
parser.add_option('-f', '--file', dest='filename',
        help='the input ESySXML file', metavar='FILE')
parser.add_option('-d', '--debug', dest='dbg', action="store_true",
        help='switch debug on', default=False)
parser.add_option('-n', '--new', action="store",
        help='output ESySXML file',
        dest='new_file_name', default='')
def main():
    (options, args) = parser.parse_args()
    if options.filename:
       filenames=list(options.filename) + args
    else:
       filenames=args
    if len(filenames)<1:
        parser.error("no input file.")

    files = [(file(arg), arg) for arg in filenames]
    for f, filename in files:
        xml = modelframe.ESySXMLParser(f.read(), debug=options.dbg)
        sims = xml.parse()
        for s in sims:
          if isinstance(s, modelframe.Simulation):
               if options.new_file_name: s.writeXML(file(options.new_file_name,'w'))
               s.run()

if __name__=='__main__':
    main()

