
##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
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
this a very simple test for inversion of Dc resistivity Data
"""

#For this hacked version, you will need to create a (non-version controlled)
#subdir called pklMsh with the mesh in it


__copyright__="""Copyright (c) 2003-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import os
import sys
from esys.escript import *
from esys.downunder import *
import esys.ripley as ripley
from esys.escript import unitsSI as U
from esys.weipa import saveSilo
import numpy                     as np
import esys.escript.pdetools     as pdetools
import pickle

import logging
logger=logging.getLogger('inv')
#logger.setLevel(logging.INFO)
logger.setLevel(logging.DEBUG)

"""
Inversion script for 3d dc resistivity data
It is highly recommended to use the same .msh file that was used for generating the
sythetic if you are using synthetic measurements.
"""

def getApparentResistivity(delPhi, n, a, current):
    resistivity=(-delPhi/current)*(n*(n+1))*pi*a
    return resistivity

current=0.5 #500 milli amps 
#load data to be inverted, the pkl file is to be generated using the forward
#modeling code adjust path to suit


pklFileName = "dcResResistiveLayer20-2-24-10.pkl"
pklFilePtr = open(pklFileName, "r")
pklFile = pickle.load(pklFilePtr)
#sourceInfo descibes the current electode setup, it containts pairs of tags as a list of tupples.
sourceInfo = pklFile["srcInfo"]
#VPSQ is a list of potential differences across different measurement pairs s and source pairs p and q
vpqs = pklFile["delPotArr"]
sampleInfo = pklFile["samples"]
#run name contains the name of the msh file to be used
runName = pklFile["runName"]
interval_a= pklFile["interval_a"]
numElectrodes=pklFile["numElectrodes"]
tags=[]
points=[]
electrodeDict = pklFile["electrodeDict"]
for i in range(numElectrodes):
    tags.append("e%d"%i)
    points.append(electrodeDict["e%d"%i][:-1]) #first three elements of electrode dict provide point info
    
    
sw=SplitWorld(getMPISizeWorld())    # Depending on how big the meshes are, we may want to be more efficient about number of worlds
buildDomains(sw,ReadGmsh, runName+".msh", 3 , diracTags=tags, diracPoints=points)


addVariable(sw,"regularization", makeLocalOnly)
addVariable(sw,"mappings", makeLocalOnly)
addVariable(sw,"fwdmodels", makeLocalOnly)

# model count is how many models this Job should create, modelnumstart is where in the global numbering should
# these jobs be?
# This imposes a constraint on the arity of functions accepted for this call
#Note that self, is needed because this will be executed in the context of a job and we want to be able to access it
def makeParamBlob(self, **kwargs):
  dom=self.domain
  DIM=dom.getDim()
  z=ContinuousFunction(dom).getX()[1]
  try:
    modelcount=kwargs["rangestart"]		#Where does our part of the range start
    modelnumstart=kwargs["rangelen"]	#How many things in our part of the range
    nLS=kwargs["numLevelSets"]
  except KeyError as e:
    raise ValueError("The following expected value was not supplied to the function: "+e.message)
  
  #====================== THIS IS WHERE THE INVERSION STARTS:

  # create regularization with two level set functions:
  z=Solution(dom).getX()[DIM-1]
  x=Solution(dom).getX()[0]
  # ask luts about this
  reg_mask = whereZero(z-inf(z)) # fix the top resistivity only.
  #====================== SETUP INITIAL VALUES
  primaryRes=100.0
  Sigma0 = Scalar(0,Function(dom))
  SigmaPrimary = Scalar(0,Function(dom))
  Sigma0.setTaggedValue("volume-1",1/primaryRes)
  Sigma0.setTaggedValue("volume-2",1/primaryRes)
  Sigma0.setTaggedValue("volume-3",1/primaryRes)
  Sigma0.setTaggedValue("volume-4",1/primaryRes)
  SigmaPrimary.setTaggedValue("volume-1",1/primaryRes)
  SigmaPrimary.setTaggedValue("volume-2",1/primaryRes)
  SigmaPrimary.setTaggedValue("volume-3",1/primaryRes)
  SigmaPrimary.setTaggedValue("volume-4",1/primaryRes)
  Sigma0.expand()
  SigmaPrimary.expand()
#===================== SETUP SAMPLES
  samples=[]
  for i in range(len(sampleInfo)):
      samples.append([electrodeDict[sampleInfo[i][0]],electrodeDict[sampleInfo[i][1]]])

  #===================== SETUP PDE FOR PRIMARY POTENTIAL
  pde=LinearPDE(dom, numEquations=1)
  kro=kronecker(dom)
  tol=1e-8
  DIM=dom.getDim()
  x = dom.getX()
  q=whereZero(x[DIM-1]-inf(x[DIM-1]))
  for i in xrange(DIM-1):
      xi=x[i]
      q+=whereZero(xi-inf(xi))+whereZero(xi-sup(xi))
  A=pde.createCoefficient('A')
  A=SigmaPrimary * kro
  dcm=mappings.DcResMapping(2*Sigma0)
  

  # iterate over source info
  n=1
  j=0
  
  # Need to make sure that the numLevelSets here is the same as the one given as a parameter
  # to the overall InversionCostFunction but don't want to pull values out of the subworld  
  regularization=Regularization(dom, numLevelSets=nLS,
				w1=np.ones(DIM), # consider gradient terms
				location_of_set_m=reg_mask) 
  
  
  # Now create a number of forward models to be processed on this subworld
  forwardmodels=[]  # this will hold the list of forward models   
  for k in range(0, modelcount):
      i=modelnumstart+k
      nCount = (numElectrodes-3) - (2*(n-1))
      loc = pdetools.Locator(ContinuousFunction(dom), samples[i])
      y_dirac=Scalar(0,DiracDeltaFunctions(dom))
      y_dirac.setTaggedValue(sourceInfo[i][0],current)
      y_dirac.setTaggedValue(sourceInfo[i][1],-current)
      pde.getSolverOptions().setTolerance(tol)
      pde.setSymmetryOn()
      pde.setValue(A=A,y_dirac=y_dirac,q=q)
      phiPrimary=pde.getSolution()
      if (vpqs[n-1][j]!=0):
	  DCRES=DcRes(dom, loc, (vpqs[n-1][j],), [sampleInfo[i]], phiPrimary, SigmaPrimary, w=(1./vpqs[n-1][j])**2)
      else:
	  DCRES=DcRes(dom, loc, (vpqs[n-1][j],), [sampleInfo[i]], phiPrimary, SigmaPrimary, w=(1.))
      j+=1
      if j == nCount:
	  n+=1
	  j=0
      forwardmodels.append((DCRES, 0))
  dcm=(mappings.DcResMapping(2*Sigma0),)	
      
  dcm, forwardmodels=InversionCostFunction.processMapsAndModels(dcm, forwardmodels, nLS) 
  
  print "dcm is",dcm
  

  trafo = regularization.getCoordinateTransformation()
  for m in forwardmodels:
      print "--------------------------------"
      print type(m[0])
      print "--------------------------------"    
      if not m[0].getCoordinateTransformation() == trafo:
	  raise ValueError("Coordinate transformation for regularization and model don't match.") 
	
	
  muVal=5000
  mu_c=None
     #cf.setTradeOffFactorsRegularization(mu=muVal, mu_c=None)	
  regularization.setTradeOffFactorsForVariation(muVal)
  regularization.setTradeOffFactorsForCrossGradient(mu_c)  
  
  
    # Need to export forward models as well
  self.exportValue("fwdmodels", forwardmodels)
  # Need to export regularization				
  self.exportValue("regularization",regularization)
  
  # Need to export dcm				

  self.exportValue("mappings", dcm)
  


# finally we can set up the cost function:
# Need to pass in subworld and world init function here
# Also need to modify so we don't take in regularization or mappings .... or forward models 
#cf=InversionCostFunction(regularization, (dcm,), tuple(forwardmodels) , splitw=sw, worldsinit_fn=makeParamBlob)
cf=InversionCostFunction(None, None, None , splitw=sw, worldsinit_fn=makeParamBlob, numLevelSets=1, 
       numModels=len(sourceInfo), numMappings=1)

#============== set regularization
#muVal=5000
#cf.setTradeOffFactorsRegularization(mu=muVal, mu_c=None)

# this code can be used to save silo at every step 
#def minimizerCallBack(k, x, Jx, g_Jxx):
 # saveSilo("/scratch/jplessis/mtinversion/silos/%d-%d.silo"%(siloCount,k),est=x,sigma=cf.getProperties(x))

#==================================== setup and run solver
# and then start running the show:
solver=MinimizerLBFGS()
solver.setCostFunction(cf)
# solver.setCallback(minimizerCallBack)
solver.setTolerance(1e-6)
solver.setMaxIterations(180)

try:
  solver.run(Scalar(0, Solution(dom)))
except (MinimizerMaxIterReached ,RuntimeError) as exceptionVal:
    print "Caught exception saving silo and rethrowing."
    print "silo save to,","/scratch/uqjdupl1/dcres/silos/sigmaResistiveLayerfailed%d-%d.silo"%(numElectrodes, muVal)
    m=solver.getResult()
    sigma = cf.getProperties(m)
    saveSilo("/scratch/uqjdupl1/dcres/silos/sigmaResistiveLayerfailed%d-%d-%d.silo"%(interval_a, numElectrodes ,muVal), sigma=sigma, con=1./sigma)
    raise exceptionVal

m=solver.getResult()
sigma = cf.getProperties(m)
saveSilo("/scratch/uqjdupl1/dcres/silos/sigmaResistiveLayer%d-%d-%d.silo"%(interval_a, numElectrodes ,muVal), sigma=sigma, con=1./sigma)
print ("output silo in:","silos/sigmaResistiveLayer%d-%d-%d.silo"%(interval_a, numElectrodes ,muVal))
resIn=[]
resFinal=[]
n=1
j=0
for i in forwardmodels:
    nCount = (numElectrodes-3) - (2*(n-1))
    phi, loc_phi=i[0].getArguments(sigma)
    delPot=loc_phi[1]-loc_phi[0]
    delPhiIn = (i[0].delphiIn)[0]
    resIn.append(primaryRes + getApparentResistivity(delPhiIn, n, interval_a, current))
    resFinal.append(primaryRes + getApparentResistivity(delPot, n, interval_a, current))
    j+=1
    if j == nCount:
        n+=1
        j=0

pklF=open(pklFileName[:-4]+"-%d-inverted.pkl"%muVal,"w")
print ("output in:",pklFileName[:-4]+"-%d-inverted.pkl"%muVal)
pklDict={"resIn":resIn,"resFinal":resFinal,"regVal":muVal}
pickle.dump(pklDict,pklF)
print "m=",m
print "sigma=", sigma
