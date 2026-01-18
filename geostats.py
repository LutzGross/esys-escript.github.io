import numpy as np
from mpi4py import MPI

from escriptcore.py_src.util import kronecker

L = 100.
nx, ny = 100,100
dx = 5.
VARIO_SILL = 1.
VARIO_RANGE = 300.
NUM_SAMPLES = 10.


#comm = MPI.COMM_WORLD
#myrank = comm.Get_rank()
myrank =0
if myrank == 0:
    xline=np.arange(0, L, dx)
    yline=np.arange(0, L, dx)
    xg, yg=np.meshgrid(xline, yline)
    xx=xg.flatten()
    yy=yg.flatten()

    h=np.zeros((len(xx), len(xx)))
    for i in range(len(xx)):
        h[i,:]=np.sqrt((xx[i]-xx)**2+(yy[i]-yy)**2)
    mean=np.zeros((len(xx),))
    cov=VARIO_SILL*np.exp(-h/(VARIO_RANGE/3))


    Q=np.random.multivariate_normal(mean, cov)
    Qg=Q.reshape(len(xline), len(yline))

    import matplotlib.pyplot as plt
    plt.figure()
    plt.imshow(Qg, extent =(min(xline), max(xline), min(yline), max(yline)) )
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title(f"A Realization for sill ={VARIO_SILL} and range ={VARIO_RANGE}")
    plt.colorbar()
    plt.show()


dom=Rectangle(nx, ny, l0=L, l1=L, comm=sub_comm)

pde = LinearPDE(domain)
pde.setSymmetryOn()
pde.setValue(A= k * kronecker(2))
x = Solution(domain).getX()

pde.setValue(A = Y=source)
pde.setValue(q=whereZero(x[0]) + whereZero(x[0]-1.0) +
               whereZero(x[1]) + whereZero(x[1]-1.0))


result2=interpolateTable(table2, x2, (xmin, ymin), (xstep, ystep), toobig)
