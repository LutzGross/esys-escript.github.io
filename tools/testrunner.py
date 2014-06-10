from __future__ import print_function
import sys

def rearrange(string):
    parts = string.split()
    parts = [parts[1][1:-1], parts[0]]
    return ".".join(parts)

def run_tests(modules):
    skiplist = []
    for module in modules:
        if module[-3:] == ".py":
            module = module[:-3]
        m = __import__(module)
        res = m.run_tests(module)
        skiplist.extend(["%s : %s"%(rearrange(str(i[0])),i[1]) for i in res.skipped])
    return skiplist

if __name__ == "__main__":
    modules = sys.argv[1:]
    if len(modules) == 0:
        print("%s missing argument, provide module to run tests on"%sys.argv[0])
        sys.exit(1)
    outputfile = None
    appendfile = None
    for n, m in enumerate(modules):
        if m.startswith("-outputfile="):
            modules.pop(n)
            outputfile = m.split("=")[1]
            break
        if m.startswith("-appendfile="):
            modules.pop(n)
            appendfile = m.split("=")[1]
            break
    skipped = run_tests(modules)
    if outputfile:
        open(outputfile, "w").writelines("\n".join(skipped))
    elif appendfile:
        open(appendfile, "a").writelines("\n".join(skipped))
    else:
        print("\n".join(skipped))
