import matplotlib
matplotlib.use('agg')

from esys.escript import *
import esys.escript.symbolic as esym
from esys.finley import Rectangle
import numpy
import os
try:
    import scipy.interpolate
    HAVE_SCIPY=True
except:
    HAVE_SCIPY=False
#set up domain and symbols
mydomain = Rectangle(l0=1.,l1=1.,n0=10, n1=10)
u = Symbol('u',(2,), dim=2)
q = Symbol('q', (2,2))
sigma = Symbol('sigma',(2,2))
theta = Symbol('theta')
# q is a rotation matrix represented by a Symbol. Values can be substituted for
# theta.
q[0,0]=cos(theta)
q[0,1]=-sin(theta)
q[1,0]=sin(theta)
q[1,1]=cos(theta)
# Theta gets substituted by pi/4 and masked to lie between .3 and .7 in the
# vertical direction. Using this masking means that when q is used it will apply
# only to the specified area of the domain.
x = Function(mydomain).getX()
q=q.subs(theta,(esym.Symconsts.pi/4)*whereNonNegative(x[1]-.30)*whereNegative(x[1]-.70))
# epsilon is defined in terms of u and has the rotation applied.
epsilon0 = symmetric(grad(u))
epsilon = matrixmult(matrixmult(q,epsilon0),q.transpose(1))
# For the purposes of demonstration, an arbitrary c with isotropic constraints
# is chosen here. In order to act as an isotropic material c is chosen such that
# c00 = c11 = c01+c1+2*c55
c00 = 10
c01 = 8; c11 = 10
c05 = 0; c15 = 0; c55 = 1
# sigma is defined in terms of epsilon
sigma[0,0] = c00*epsilon[0,0]+c01*epsilon[1,1]+c05*2*epsilon[1,0]
sigma[1,1] = c01*epsilon[0,0]+c11*epsilon[1,1]+c15*2*epsilon[1,0]
sigma[0,1] = c05*epsilon[0,0]+c15*epsilon[1,1]+c55*2*epsilon[1,0]
sigma[1,0] = sigma[0,1]
sigma0=matrixmult(matrixmult(q.transpose(1),sigma),q)
# set up boundary conditions
x=mydomain.getX()
gammaD=whereZero(x[1])*[1,1]
yconstraint = FunctionOnBoundary(mydomain).getX()[1]
# The nonlinear PDE is set up, the values are substituted in and the solution is
# calculated, y represents an external shearing force acting on the domain.
# In this case a force of magnitude 50 is acting in the x[0] direction.
p = NonlinearPDE(mydomain, u, debug=NonlinearPDE.DEBUG0)
p.setValue(X=sigma0,q=gammaD,y=[-50,0]*whereZero(yconstraint-1),r=[1,1])
v = p.getSolution(u=[0,0])

# Create vector plot of displacement field v
if HAVE_SCIPY:
    import matplotlib.pyplot as plt

    # Get coordinates and velocity components
    x_data = mydomain.getX()[0].toListOfTuples()
    y_data = mydomain.getX()[1].toListOfTuples()
    v0_data = interpolate(v[0], mydomain.getX().getFunctionSpace()).toListOfTuples()
    v1_data = interpolate(v[1], mydomain.getX().getFunctionSpace()).toListOfTuples()

    # Create regular grid for interpolation
    x_grid = numpy.linspace(0., 1., 20)
    y_grid = numpy.linspace(0., 1., 20)
    X_grid, Y_grid = numpy.meshgrid(x_grid, y_grid)

    # Interpolate vector components to grid
    v0_grid = scipy.interpolate.griddata((x_data, y_data), v0_data, (X_grid, Y_grid), method='linear')
    v1_grid = scipy.interpolate.griddata((x_data, y_data), v1_data, (X_grid, Y_grid), method='linear')

    # Create quiver plot
    plt.figure(figsize=(8, 8))
    plt.quiver(X_grid, Y_grid, v0_grid, v1_grid, scale=500.0, scale_units='xy')
    plt.xlabel('x[0]')
    plt.ylabel('x[1]')
    plt.title('Displacement field v')
    plt.axis('equal')

    # Save in the same directory as this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, 'anisotropy_vector.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Vector plot saved to {output_path}")
else:
    print("scipy not available, skipping vector plot")


