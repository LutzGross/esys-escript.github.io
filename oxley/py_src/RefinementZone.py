
from esys.oxley import *

def RefinementZone(dim):
	if dim == 2:
		return RefinementZone2D
	elif dim == 3:
		return RefinementZone3D
	else:
		print("Invalid dimension")
