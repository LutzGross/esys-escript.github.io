
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.uq.edu.au/esscc/escript-finley"

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


from esys.escript import modelframe
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
    
