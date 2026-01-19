
##############################################################################
#
# Copyright (c) 2014-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file in the repository root for contributors and development history
# https://github.com/LutzGross/esys-escript.github.io/blob/master/CREDITS
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

Running all test classes::

    import your_test_module_or_file as tests
    tests.run_tests("your_test_module_or_file")

Running specific test classes::

    import your_test_module_or_file as tests
    tests.run_tests("your_test_module_or_file", [tests.Test_classA, tests.Test_classB])

Running a specific test within a class::

    import your_test_module_or_file as tests
    tests.run_single_test(tests.Test_classname("test_functionname"))

Printing the list of skipped tests::

    results = tests.run_tests("name")
    for skipped_test in results.skipped:
        print(skipped_test)

"""


__copyright__="""Copyright (c) 2014-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

import esys.escriptcore.utestselect as unittest
from .util import getMPIRankWorld
from .escriptcpp import MPIBarrierWorld
import os
import sys

def __add_tests(suite, modulename):
    """
    Discovers and adds all test classes from a module to a test suite.

    Test classes are identified by having a name starting with "Test" and
    being a subclass of ``unittest.TestCase``.

    :param suite: the test suite to add tests to
    :type suite: ``unittest.TestSuite``
    :param modulename: name of the module to discover tests from
    :type modulename: ``str``
    """
    import inspect
    clsmembers = inspect.getmembers(sys.modules[modulename], inspect.isclass)
    for name, cls in clsmembers:
        if modulename == cls.__module__ and name.startswith("Test") \
                and issubclass(cls, unittest.TestCase):
            suite.addTest(unittest.TestLoader().loadTestsFromTestCase(cls))

def run_single_test(test, exit_on_failure=False):
    """
    Runs a single test case and returns the result.

    :param test: the test instance to run
    :type test: ``unittest.TestCase``
    :param exit_on_failure: if ``True``, exits with code 1 on test failure
    :type exit_on_failure: ``bool``
    :return: the test result object containing pass/fail information
    :rtype: ``unittest.TestResult``
    """
    suite = unittest.TestSuite()
    suite.addTest(test)
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if exit_on_failure and not s.wasSuccessful():
        sys.exit(1)
    return s

def run_tests(modulename, classes = [], exit_on_failure = False):
    """
    Runs test classes from a module and returns the result.

    If no specific classes are provided, all test classes (those starting
    with "Test" and subclassing ``unittest.TestCase``) are discovered and
    run automatically.

    In MPI environments, only rank 0 produces output; other ranks run
    silently.

    :param modulename: name of the module containing the tests
    :type modulename: ``str``
    :param classes: specific test classes to run; if empty, all tests
                    in the module are discovered and run
    :type classes: ``list`` of ``unittest.TestCase`` subclasses
    :param exit_on_failure: if ``True``, exits with code 1 on test failure
    :type exit_on_failure: ``bool``
    :return: the test result object containing pass/fail information
    :rtype: ``unittest.TestResult``
    """
    rank = getMPIRankWorld()
    stream=sys.stderr #default
    verb=2
    if rank > 0:
        stream=open(os.devnull,'w')
        verb=0

    suite = unittest.TestSuite()
    if len(classes) == 0:
        __add_tests(suite, modulename)
    else:
        for test_class in classes:
            suite.addTest(unittest.TestLoader().loadTestsFromTestCase(test_class))
    s=unittest.TextTestRunner(stream=stream,verbosity=verb).run(suite)
    if exit_on_failure and not s.wasSuccessful():
        sys.stderr.flush()
        MPIBarrierWorld()
        sys.exit(1)
    return s

