
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################
import sys

fail_format = """======================================================================
FAIL: {0} {1}
----------------------------------------------------------------------
{2}
----------------------------------------------------------------------
"""

def rearrange(string):
    parts = string.split()
    parts = [parts[1][1:-1], parts[0]]
    return ".".join(parts)

def rearrange_to_default(string):
    parts = string.split()
    print(parts)
    return ".".join(parts)

def run_tests(modules, exit_on_failure=False):
    skiplist = []
    faillist = []
    for module in modules:
        if module[-3:] == ".py":
            module = module[:-3]
        m = __import__(module)
        res = m.run_tests(module, exit_on_failure=exit_on_failure)
        skiplist.extend(["%s : %s\n"%(rearrange(str(i[0])),i[1]) for i in res.skipped])
        faillist.extend([fail_format.format(str(i[0]).split()[0],str(i[0]).split()[1], i[1]) for i in res.failures+res.errors])
    return skiplist, faillist

if __name__ == "__main__":
    modules = sys.argv[1:]
    if len(modules) == 0:
        print("%s missing argument, provide module to run tests on"%sys.argv[0])
        sys.exit(1)
    skipfile = None
    skipappendfile = None
    failfile = None
    failappendfile = None
    exit = False
    n = 0
    while n < len(modules):
        m = modules[n]
        if m.startswith("-skipfile="):
            modules.pop(n)
            skipfile = m.split("=")[1]
            continue
        if m.startswith("-skipappendfile="):
            modules.pop(n)
            skipappendfile = m.split("=")[1]
            continue
        if m.startswith("-failfile="):
            modules.pop(n)
            failfile = m.split("=")[1]
            continue
        if m.startswith("-failappendfile="):
            modules.pop(n)
            failappendfile = m.split("=")[1]
            continue
        if m == "-exit":
            modules.pop(n)
            exit = True
            continue
        n += 1

    skipped, failed = run_tests(modules, exit_on_failure=exit)
    if skipfile:
        open(skipfile, "w").writelines("".join(skipped))
    elif skipappendfile:
        open(skipappendfile, "a").writelines("".join(skipped))
    else:
        print("".join(skipped))

    if failfile:
        open(failfile, "w").writelines("".join(failed))
    elif failappendfile:
        open(failappendfile, "a").writelines("".join(failed))
    else:
        print("".join(failed))
