
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
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


import lamesource
import sys

def buildTempAndSummation(dim, ids, temps, summations, forced_substitutions = []):
    declarations = {}
    sumStatements = {}
    for k in range(dim):
        for m in range(dim):
            declarations[(k,m)] = []
            sumStatements[(k,m)] = []
            zeroes = [] #tracks tmpvars that will always be zero
            nonzeroes = [] #tracks tmpvars that will possibly be non-zero
            use_counts = {} #tmp var usage counts, used later for culling
#temp declarations
            #hold on to the building expression for every identifier
            var_expressions = {}
            for line in temps:
                newline = line.format(k,m)
                insert = False
                for index in ids.iterkeys():
                    if index in newline:
                        insert = True
                        break
                if not insert:
                    value = newline.split("=")[1]
                    for tmp in nonzeroes:
                        if tmp in value:
                            x = use_counts.get(tmp, 0)
                            use_counts[tmp] = x+1
                            insert = True
                i = newline.index("tmp")
                identifier = newline[i:].split()[0]
                expression = newline.split("= ")[1][:-1]
                if insert:
                    nonzeroes.append(identifier)
                    var_expressions[identifier] = expression
                else:
                    zeroes.append(identifier)
#summations
            with_nonzero = [] #only those with some non-zero temp var
            for line in summations:
                newline = replaceZeroes(line.format(k,m),zeroes)
                for nz in nonzeroes:
                    if nz in newline:
                        with_nonzero.append(newline)
                        components = newline[:-1].split("+=")[1].lstrip().replace(" + ", "|").replace(" - ", "|")
                        if components[0] == "-":
                            components = components[1:]
                        components = components.split("|")
                        for var in components:
                            x = use_counts.get(var, 0)
                            use_counts[var] = x+1
                        break
                for fs in forced_substitutions:
                    if fs in nonzeroes:
                        use_counts[fs] = 1
                for z in zeroes:
                    use_counts[z] = 0
            #only interested in keeping variables declarations with 2 or more uses
            for var in nonzeroes:
                #0 we don't care and 1 we substitute the expression
                if use_counts.get(var, 0) > 1:
                    declarations[(k,m)].append("    const double %s = %s;"%(var, var_expressions[var]))
            #remove zeroes and replace single use tmpvars with their expression
            for line in with_nonzero:
                for key in use_counts.iterkeys():
                    if use_counts[key] == 0:
                        line = line.replace("%s;"%key, ";")
                        line = line.replace("%s "%key, " ")
                    elif use_counts[key] < 2:
                        s = var_expressions[key]
                        line = line.replace("%s;"%key, s+";")
                        line = line.replace("%s "%key, s+" ")
                #print only if there's a right-hand-side of the expression
                if "+=;" not in line and "+=;" not in line:
                    sumStatements[(k,m)].append("    "+line)
    return declarations, sumStatements

def replaceZeroes(line, zeroes):
    for zero in zeroes:
        if zero+"|" in line:
            line = line.replace(" + %s|"%zero, "")
            line = line.replace("%s| + "%zero, "0 + ") #there's a better solution
            line = line.replace(" - %s|"%zero, "")     #it's just to stop
            line = line.replace("%s| - "%zero, "0 - ") # a - b -> b instead of -b
            line = line.replace("=%s|;"%zero,"=;")
    if line:
        return line.replace("|","").replace(" - 0", "").replace(" + 0", "").replace("+= 0 +", "+=").replace("+= 0 - ", "+=-")
    return line

def print2DAExpanded():
    dim = 2
    quads = 2**dim
    ids = {}
    for i in range(dim):
        for j in range(dim):
            ids["{0}{0}{1}{1}".format(i,j)] = None
            ids["{0}{1}{1}{0}".format(i,j)] = None
            ids["{0}{1}{0}{1}".format(i,j)] = None

    for name in sorted(ids.iterkeys()):
        print("double A_{0}[{1}] =".format(name,quads), "{0};")
    #   ijji += mu
    #   ijij += mu
    #   iijj += lambda
    print("if (!mu.isEmpty()) {\n    const double *mu_p = mu.getSampleDataRO(e);")
    completed = {}
    for i in range(dim):
        for j in range(dim):
            for q in range(quads):
                if i == j:
                    print("    A_{0}{0}{0}{0}[{1}] += 2*mu_p[{1}];".format(i,q))
                else:
                    print("    A_{0}{1}{1}{0}[{2}] += mu_p[{2}];".format(i,j,q))
                    print("    A_{0}{1}{0}{1}[{2}] += mu_p[{2}];".format(i,j,q))
    print("}\nif (!lambda.isEmpty()) {\n    const double *lambda_p = lambda.getSampleDataRO(e);")
    for i in range(dim):
        for j in range(dim):
            for q in range(quads):
                print("    A_{0}{0}{1}{1}[{2}] += lambda_p[{2}];".format(i,j,q))
    print("}")

    decl, sums = buildTempAndSummation(dim, ids, lamesource.expanded2Dtemps, lamesource.expanded2Dsummations)
    for k in range(dim):
        for m in range(dim):
            print("{")
            print("\n".join(decl[(k,m)]))
            print("\n".join(sums[(k,m)]))
            print("}")

def print2DAReduced():
    dim = 2
    quads = 2**dim
    ids = {}
    for i in range(dim):
        for j in range(dim):
            ids["{0}{0}{1}{1}".format(i,j)] = None
            ids["{0}{1}{1}{0}".format(i,j)] = None
            ids["{0}{1}{0}{1}".format(i,j)] = None

    for name in sorted(ids.iterkeys()):
        print("double A_{0} =".format(name,quads), "0;")
    #   ijji += mu
    #   ijij += mu
    #   iijj += lambda

    print("if (!mu.isEmpty()) {\n    const double *mu_p = mu.getSampleDataRO(e);")

    for i in range(dim):
        for j in range(dim):
            if i == j:
                print("    A_{0}{0}{0}{0} += 2*mu_p[0];".format(i))
            else:
                print("    A_{0}{1}{1}{0} += mu_p[0];".format(i,j))
                print("    A_{0}{1}{0}{1} += mu_p[0];".format(i,j))
    print("}")
    print("if (!lambda.isEmpty()) {\n    const double *lambda_p = lambda.getSampleDataRO(e);")
    for i in range(dim):
        for j in range(dim):
            print("    A_{0}{0}{1}{1} += lambda_p[0];".format(i,j))
    print("}")

    lines = [
        "const double tmp_0 = 6*w1*(A_{0}0{1}1 - A_{0}1{1}0);",
        "const double tmp_1 = 6*w1*(A_{0}0{1}1 + A_{0}1{1}0);",
        "const double tmp_2 = 6*w1*(-A_{0}0{1}1 - A_{0}1{1}0);",
        "const double tmp_3 = 6*w1*(-A_{0}0{1}1 + A_{0}1{1}0);"
    ]

    zeroes = []
    nonzeroes = []
    for k in range(dim):
        for m in range(dim):
            if k == m:
                for q in range(quads):
                    zeroes.append("tmp{0}{0}_{1}".format(m, q))
                continue
            for line in lines:
                newline = line.format(k,m)
                insert = False
                for index in ids.iterkeys():
                    if index in newline:
                        insert = True
                        break
                i = newline.index("tmp")
                if insert:
                    print(newline)
                    nonzeroes.append(newline[i:].split()[0].rstrip())
                else:
                    zeroes.append(newline[i:].split()[0].rstrip())

    for k in range(dim):
        for m in range(dim):
            for line in lamesource.reduced2Dsummations:
                newline = replaceZeroes(line.format(k,m),zeroes)
                if not newline:
                    continue
                if k != m:
                    st = "A_{0}0{1}0".format(k,m)
                    if st in newline:
                        i = newline.index(st)
                        newline = newline[:i-3] + newline[i+12:]
                    st = "A_{0}1{1}1".format(k,m)
                    if st in newline:
                        i = newline.index(st)
                        newline = newline[:i-2] + newline[i+11:]
                print(newline)

def print3DAExpanded():
    dim = 3
    quads = 2**dim
    ids = {}
    for i in range(dim):
        for j in range(dim):
            ids["{0}{0}{1}{1}".format(i,j)] = None
            ids["{0}{1}{1}{0}".format(i,j)] = None
            ids["{0}{1}{0}{1}".format(i,j)] = None

    for name in sorted(ids.iterkeys()):
        print("double A_{0}[{1}] =".format(name,quads), "{0};")
    #   ijji += mu
    #   ijij += mu
    #   iijj += lambda
    print("if (!mu.isEmpty()) {\n    const double *mu_p = mu.getSampleDataRO(e);")
    for i in range(dim):
        for j in range(dim):
            for q in range(quads):
                if i == j:
                    print("    A_{0}{0}{0}{0}[{1}] += 2*mu_p[{1}];".format(i,q))
                else:
                    print("    A_{0}{1}{1}{0}[{2}] += mu_p[{2}];".format(i,j,q))
                    print("    A_{0}{1}{0}{1}[{2}] += mu_p[{2}];".format(i,j,q))
    print("}")
    print("if (!lambda.isEmpty()) {\n    const double *lambda_p = lambda.getSampleDataRO(e);")
    for i in range(dim):
        for j in range(dim):
            for q in range(quads):
                print("    A_{0}{0}{1}{1}[{2}] += lambda_p[{2}];".format(i,j,q))
    print("}")

    decl, sums = buildTempAndSummation(dim, ids, lamesource.expanded3Dtemps, lamesource.expanded3Dsummations)
    for k in range(dim):
        for m in range(dim):
            print("{")
            print("\n".join(decl[(k,m)]))
            print("\n".join(sums[(k,m)]))
            print("}")

def print3DAReduced():
    dim = 3
    ids = {}
    for i in range(dim):
        for j in range(dim):
            ids["{0}{0}{1}{1}".format(i,j)] = None
            ids["{0}{1}{1}{0}".format(i,j)] = None
            ids["{0}{1}{0}{1}".format(i,j)] = None
            
    for i in sorted(ids.iterkeys()):
        print("double Aw%s = 0;"%i)
        
    print("if (!mu.isEmpty()) {\n    const double *mu_p = mu.getSampleDataRO(e);")
    for i in range(dim):
        for j in range(dim):
            if i == j:
                print("    Aw{0}{0}{0}{0} += 2*mu_p[0];".format(i))
            else:
                print("    Aw{0}{1}{1}{0} += mu_p[0];".format(i,j))
                print("    Aw{0}{1}{0}{1} += mu_p[0];".format(i,j))
    
    print("}\nif (!lambda.isEmpty()) {\n    const double *lambda_p = lambda.getSampleDataRO(e);")
    for i in range(dim):
        for j in range(dim):
            print("    Aw{0}{0}{1}{1} += lambda_p[0];".format(i,j))
    print("}")
    
    for k in range(dim):
        for m in range(dim):
            for line in ["Aw{0}0{1}0 *= 8*w27;","Aw{0}0{1}1 *= 12*w8;","Aw{0}0{1}2 *= 12*w11;",
                    "Aw{0}1{1}0 *= 12*w8;","Aw{0}1{1}1 *= 8*w22;","Aw{0}1{1}2 *= 12*w10;",
                    "Aw{0}2{1}0 *= 12*w11;","Aw{0}2{1}1 *= 12*w10;","Aw{0}2{1}2 *= 8*w13;"]:
                newline = line.format(k,m)
                found = False
                for ident in ids.iterkeys():
                    if ident in newline:
                        found = True
                        break
                if found:
                    print(newline)
    
    decs, sums = buildTempAndSummation(dim, ids, lamesource.reduced3Dtemps,
            lamesource.reduced3Dsummations,
            forced_substitutions = ["tmp12","tmp13","tmp14","tmp21","tmp22","tmp23"])
    for k in range(dim):
        for m in range(dim):
            print("{")
            print("\n".join(decs[(k,m)]))
            print("\n".join(sums[(k,m)]))
            print("}")
    


def printAReduced(dim):
    if dim == 2:
        return print2DAReduced()
    elif dim == 3:
        return print3DAReduced()
    else:
        raise

def printAExpanded(dim):
    if dim == 2:
        return print2DAExpanded()
    elif dim == 3:
        return print3DAExpanded()
    else:
        raise

if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1] not in ["R", "E"] or sys.argv[2] not in ["2","3"]:
        print("Usage: {0} [R]educed/[E]xpanded dimensions\nE.g. {0} R 3")
        exit(1)
    dim = int(sys.argv[2])
    if sys.argv[1] == "R":
        printAReduced(dim)
    else:
        printAExpanded(dim)

