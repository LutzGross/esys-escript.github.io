#/usr/bin/python
# $Id:$

#
#      COPYRIGHT ACcESS 2004 -  All Rights Reserved
#
#   This software is the property of ACcESS.  No part of this code
#   may be copied in any form or by any means without the expressed written
#   consent of ACcESS.  Copying, use or modification of this software
#   by any unauthorised person is illegal unless that
#   person has a software license agreement with ACcESS.
#

"""
some benchmarks for tetsing the finley solver. The idea is to develop a set of standart benchmarks

  * Laplace2Dorder1_?k
  * Laplace3Dorder2_?k

where ? is approximatively the number of unknowns in 1000.

@var __author__: name of author
@var __licence__: licence agreement
var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"
__licence__="contact: esys@access.uq.edu.au"
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision:$"
__date__="$Date:$"

from esys.finley.finleybench import *
from esys.escript.benchmark import BenchmarkSuite,Benchmark

thlist=[1,2,4,8,16,32]
# thlist=[1,2,4,8,16,32,64,128]
# thlist=[1,2,4,8,16,32,64,128]
show=True
ff=FinleyFilter()

opt1=FinleyOptions(solver_method=LinearPDE.PCG,preconditioner=LinearPDE.JACOBI,verbose=show)
opt2=FinleyOptions(solver_method=LinearPDE.PCG,preconditioner=LinearPDE.ILU0,verbose=show)

bm_L2Do1=Benchmark(name="Laplace 2D (order 1)")
bm_L2Do1.addProblem(Laplace2DOrder1_30k())
bm_L2Do1.addProblem(Laplace2DOrder1_60k())
bm_L2Do1.addProblem(Laplace2DOrder1_120k())
bm_L2Do1.addProblem(Laplace2DOrder1_240k())
bm_L2Do1.addProblem(Laplace2DOrder1_480k())
bm_L2Do1.addProblem(Laplace2DOrder1_960k())
bm_L2Do1.addProblem(Laplace2DOrder1_1920k())
bm_L2Do1.addProblem(Laplace2DOrder1_3840k())
bm_L2Do1.addProblem(Laplace2DOrder1_7680k())
bm_L2Do1.addProblem(Laplace2DOrder1_15360k())
bm_L2Do1.addOptions(opt1)
bm_L2Do1.addOptions(opt2)

bm_L2Do2=Benchmark("Laplace 2D (order 2)")
bm_L2Do2.addProblem(Laplace2DOrder2_30k())
bm_L2Do2.addProblem(Laplace2DOrder2_60k())
bm_L2Do2.addProblem(Laplace2DOrder2_120k())
bm_L2Do2.addProblem(Laplace2DOrder2_240k())
bm_L2Do2.addProblem(Laplace2DOrder2_480k())
bm_L2Do2.addProblem(Laplace2DOrder2_960k())
bm_L2Do2.addProblem(Laplace2DOrder2_1920k())
bm_L2Do2.addProblem(Laplace2DOrder2_3840k())
bm_L2Do2.addProblem(Laplace2DOrder2_7680k())
bm_L2Do2.addProblem(Laplace2DOrder2_15360k())
bm_L2Do2.addOptions(opt1)
bm_L2Do2.addOptions(opt2)

bm_L3Do1=Benchmark("Laplace 3D (order 1)")
bm_L3Do1.addProblem(Laplace3DOrder1_30k())
bm_L3Do1.addProblem(Laplace3DOrder1_60k())
bm_L3Do1.addProblem(Laplace3DOrder1_120k())
bm_L3Do1.addProblem(Laplace3DOrder1_240k())
bm_L3Do1.addProblem(Laplace3DOrder1_480k())
bm_L3Do1.addProblem(Laplace3DOrder1_960k())
bm_L3Do1.addProblem(Laplace3DOrder1_1920k())
bm_L3Do1.addProblem(Laplace3DOrder1_3840k())
bm_L3Do1.addProblem(Laplace3DOrder1_7680k())
bm_L3Do1.addProblem(Laplace3DOrder1_15360k())
bm_L3Do1.addOptions(opt1)
bm_L3Do1.addOptions(opt2)

bm_L3Do2=Benchmark("Laplace 3D (order 2)")
bm_L3Do2.addProblem(Laplace3DOrder2_30k())
bm_L3Do2.addProblem(Laplace3DOrder2_60k())
bm_L3Do2.addProblem(Laplace3DOrder2_120k())
bm_L3Do2.addProblem(Laplace3DOrder2_240k())
bm_L3Do2.addProblem(Laplace3DOrder2_480k())
bm_L3Do2.addProblem(Laplace3DOrder2_960k())
bm_L3Do2.addProblem(Laplace3DOrder2_1920k())
bm_L3Do2.addProblem(Laplace3DOrder2_3840k())
bm_L3Do2.addProblem(Laplace3DOrder2_7680k())
bm_L3Do2.addProblem(Laplace3DOrder2_15360k())
bm_L3Do2.addOptions(opt1)
bm_L3Do2.addOptions(opt2)



bms=BenchmarkSuite("Paso/Finley")
bms.addBenchmark(bm_L2Do1)
bms.addBenchmark(bm_L2Do2)
bms.addBenchmark(bm_L3Do1)
bms.addBenchmark(bm_L3Do2)
bms.run(scale=thlist)
print bms.getHTML(filter=ff)
