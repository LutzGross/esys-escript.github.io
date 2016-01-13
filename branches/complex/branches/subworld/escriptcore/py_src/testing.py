
##############################################################################
#
# Copyright (c) 2014 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

"""
A small set of functions to automate test discovery and running and also 
allow specific tests to be run without modifying the test scripts themselves.

When imported into a test file, any external script can import the test file
module and run tests as they wish.

Both specific and general forms of the functions return the result state,
allowing access to list all skipped tests and other data kept by the 
`TestTextResults` returned.

Examples:

  Running all test classes
    ```
    import your_test_module_or_file as tests
    tests.run_tests("your_test_module_or_file")
    ```

  Running specific test classes
    ```
    import your_test_module_or_file as tests
    tests.run_tests("your_test_module_or_file", [tests.Test_classA,
            tests.Test_classB])
    ```
    
  Running a specific test within a class
    ```
    import your_test_module_or_file as tests
    tests.run_single_test(tests.Test_classname("test_functionname"))
    ```

  Printing the list of skipped tests
    ```
    results = tests.run_tests("name")
    for skipped_test in results.skipped:
        print(skipped_test)
    ```

"""

__copyright__="""Copyright (c) 2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import esys.escriptcore.utestselect as unittest
import sys

def __add_tests(suite, modulename):
    import inspect
    clsmembers = inspect.getmembers(sys.modules[modulename], inspect.isclass)
    for name, cls in clsmembers:
        if modulename == cls.__module__ and name.startswith("Test") \
                and issubclass(cls, unittest.TestCase):
            suite.addTest(unittest.makeSuite(cls))

def run_single_test(test, exit_on_failure=False):
    suite = unittest.TestSuite()
    suite.addTest(test)
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if exit_on_failure and not s.wasSuccessful():
        sys.exit(1)
    return s

def run_tests(modulename, classes = [], exit_on_failure = False):
    suite = unittest.TestSuite()
    if len(classes) == 0:
        __add_tests(suite, modulename)
    else:
        for test_class in classes:
            suite.addTest(unittest.makeSuite(test_class))
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if exit_on_failure and not s.wasSuccessful():
        sys.exit(1)
    return s

