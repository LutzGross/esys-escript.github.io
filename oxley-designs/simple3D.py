from .tools import *

world_comm = MPI.COMM_WORLD
domain1 = Brick(n0=20, n1=20, n2 = 10, l0=1000, l1=1000, l2=800, origin0=-500, origi1=-500, origin2=-800,  comm=world_comm, framework=None)
# alternatively and then switch to Brick or Rectangle
domain1 = Block(n=(20, 20, 10), l0=(1000., 1000., 1000.), origin=(-500, -500, -800.),  comm=world_comm, framework=None)

refiner = Refiner(levels_max = 5)
refner.add(Sphere(center=(0.,0, -400), radius=200, tagname="Anomaly", resolution=5))
refiner.add(PlaneInterface(origin=(0.,0, 400), normal=(0,0,1), tagname="Deep", resolution=10))
domain2=refiner(domain1)

# now we can do things like this:
# coordinates on all nodes:
x_cf = ContinousFuction(domain2).getX()
# identical to x_cf
x_rcf = ReducedContinousFuction(domain2).getX()
# on 8 quadrature points
x_f = Function(domain2).getX()
# on element center
x_rf = ReducedFunction(domain2).getX()
# on 4 quadrature point of faces
x_fb = FunctionOnBoundary(domain2).getX()
# on element center
x_rfb = ReducedFunctionOnBoundary(domain2).getX()
# coordinates of the DOFs = *nodes that are not hanging nodes*.
# labeling depends on framework!!!!
x_s = Solution(domain2).getX()
# equivalent to x_s
x_rs = ReducedSolution(domain2).getX()

# interpolation: for instance
interpolate(x_cf, ContinuousFunction(domain1))
# etc

# then we want to do things like:
mypde=SingleLinearPDE(domain2)
kappa=Scalar(1., Function(domain2))
kappa.setTaggedValue("Anomaly",10)
kappa.setTaggedValue("Deep",10)

input=Scalar(0., FunctionOnBoundary(domain2))
input.setTaggedValue("bottom",10)

mypde.setValue(A=kappa * kronecker(mydomain),y=input, q=whereZero(x_s[2]))
u=mypde.getSolution()
# etc.

