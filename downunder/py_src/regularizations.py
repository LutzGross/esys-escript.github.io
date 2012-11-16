
##############################################################################
#
# Copyright (c) 2003-2012 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development since 2012 by School of Earth Sciences
#
##############################################################################

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['Regularization']

from costfunctions import CostFunction


import numpy as np
from esys.escript.linearPDEs import LinearPDE, IllegalCoefficientValue
from esys.escript import Data, grad, inner, integrate, kronecker, boundingBoxEdgeLengths, interpolate, vol
from esys.escript.pdetools import ArithmeticTuple

class Regularization(CostFunction):
    """
    The regularization term for the level set function `m` within the cost function J for an inversion:
    
    *J(m)=1/2 * sum_k imtegrate( mu_0[k]*s0[k] * m_k**2 + mu_1[k]*s1[k,i] * m_{k,i}**2) + sum_l<k mu_c[l,k] sc[l,l]* |curl(m_k) x curl(m_l)|^2*
    
    where s0[k], s1[k,i] and  sc[k,l] are non-negative scaling factors and mu_0[k], mu_1[k], mu_c[l,k] are weighting factors
    which may be alter during the inversion. The scaling factors are normalized such that their integrals over the
    domain are constant:
   
    *integrate(s0[k])=1*
    *integrate(inner(s1[k,:],L[:])=1*
    *integrate(inner(sc[l,k]*L**4)=1*
    
    """
    def __init__(self, domain, numLevelSets=1,  
                       s0=None, s1=None, sc=None,
                       location_of_set_m=Data(), 
                       useDiagonalHessianApproximation=True, tol=1e-8):
        """
        initialization
        
        :param domain: domain 
        :type domain: `Domain`
        :param numLevelSets: number of level sets
        :type numLevelSets: ``int``
        :param s0: scaling factor for the m**2 term. If not set zero is assumed.
        :type s0: ``Scalar`` if `numLevelSets`==1 or `Data` object of shape ('numLevelSets`,) if numLevelSets > 1
        :param s1: scaling factor for the grad(m_i) terms. If not set zero is assumed.
        :type s1: ``Vector`` if `numLevelSets`==1 or `Data` object of shape (`numLevelSets`,DIM) if numLevelSets > 1.
        :param sc: scaling factor for the cross correlation terms. If not set zero is assumed. Used for the case if numLevelSets > 1 only.
                   values `sc[l,k]``  in the lower triangle (l<k) are used only.
        :type sc: `Data` object of shape (`numLevelSets`,`numLevelSets`)     
        :param location_of_set_m: marks location of zero values of the level set function `m` by a positive entry.
        :type location_of_set_m: ``Scalar`` if `numLevelSets`==1 or `Data` object of shape ('numLevelSets`,) if numLevelSets > 1. 
        :param useDiagonalHessianApproximation: if True cross correllation terms between level set components are ignored when calculating 
                                                approximations of the inverse of the Hessian Operator. This can speep-up the calculation of 
                                                the inverse but may lead to an increase of the number of iteration steps in the inversion.
        :type useDiagonalHessianApproximation: ``bool``
        :param tol: toleramce when solving the PDE for the inverse of the Hessian Operator
        :type tol: positive ``float``
        
        
        """
        if numLevelSets>1:
	      raise ValueError("Currently only numLevelSets<=1 is supportered.")
	if s0 == None and s1==None:
	      raise ValueError("Values for s0 or s1 must be given.")
	
        self.__domain=domain
        DIM=self.__domain.getDim()
        self.__L2=np.asarray(boundingBoxEdgeLengths(domain))**2
        self.__L4=np.sum(self.__L2)**2
        self.__numLevelSets=numLevelSets
        
        if self.__numLevelSets > 1:
            self.__useDiagonalHessianApproximation=useDiagonalHessianApproximation
        else:
	    self.__useDiagonalHessianApproximation=True 
	self._update_Hessian=True
	
	self.__pde=LinearPDE(self.__domain, numEquations=self.__numLevelSets)
        self.__pde.getSolverOptions().setTolerance(tol)
        self.__pde.setSymmetryOn()
        self.__pde_is_set=False
        try:
	  self.__pde.setValue(q=location_of_set_m)
	except IllegalCoefficientValue: 
	  raise ValueError("Unable to set location of fixed level set function.")
	
	self.__total_num_weights=2*numLevelSets+((numLevelSets-1)*numLevelSets)/2
	self.__weight_index=[]  # this is a mapping from the relevant mu-coefficients to the set of all mu-coefficients
	                   # we count s0, then s1, then sc (k<l).
	# THE S0 weighting factor
	n=0
	VV=vol(domain)
        if not s0 is None:
	    s0 = interpolate(s0,self.__pde.getFunctionSpaceForCoefficient('D'))	    
	    s=s0.getShape()
	    if numLevelSets == 1 :
	         if s == () :
		     V=integrate(s0)
		     if V > 0:
		       self.__weight_index.append(n)
		       s0*=VV/V
	             else:
		       s0=None
		 else:
		     raise ValueError("Unexpected shape %s for weight s0."%s)
            else:
	       	 if s == (numLevelSets,):
		     for k in xrange(numLevelSets):
		        V=integrate(s0[k])
		        if V > 0: 
		            self.__weight_index.append(n+k)
		            s0[k]*=VV/V
		 else:
		     raise ValueError("Unexpected shape %s for weight s0."%s)   
        self.__s0=s0
        n+=numLevelSets
        
        # The S1 weighting factor
        if not s1 is None:
	    s1 = interpolate(s1,self.__pde.getFunctionSpaceForCoefficient('A'))
	    s=s1.getShape()
	    
	    if numLevelSets == 1 :
	         if s == (DIM,) :
		       V=integrate(inner(s1, 1/self.__L2))
		       if V > 0:
			  self.__weight_index.append(n)
		          s1*=VV/V
		       print "REG SCALE = ",s1
		 else:
		     raise ValueError("Unexpected shape %s for weight s1."%s)
            else:
	       	 if s == (numLevelSets,DIM):
		     for k in xrange(numLevelSets):
		        for i in xrange(DIM):
			    ww=s1[k,:]
			    V=integrate(inner(ww,1/self.__L2))
		            if V > 0: 
		               self.__weight_index.append(n+i)
		               s1[k,:]=ww*(VV/V)
		 else:
		     raise ValueError("Unexpected shape %s for weight s1."%s)      
		   
        self.__s1=s1   
        n+=numLevelSets

        # The SC weighting factor
        if not sc is None:
	    if numLevelSets == 1 :
              sc=None
            else:
              sc = interpolate(sc,self.__pde.getFunctionSpaceForCoefficient('A'))
	      s=sc.getShape()
	      if s == (numLevelSets,numLevelSets):
		     for k in xrange(numLevelSets):
		        sc[k,k]=0.
  		        for l in xrange(k):
			    ww=sc[l,k]
			    V=integrate(ww)
		            if V > 0: 
		               self.__weight_index.append(n+k*numLevelSets+l)
		               sc[l,k]=ww*VV/V*self.__L4
		               sc[k,l]=0
              else:
		     raise ValueError("Unexpected shape %s for weight s0."%s)      
		   
        self.__sc=sc
	self.setWeights()
	    
    def getDomain(self):
        """
        return the domain of the regularization term
        :rtype: ``Domain``
        """
        return self.__domain
        
    def getNumLevelSets(self):
        """
        return the number of level set functions
        :rtype: ``int``
        """
        return self.__numLevelSets

    def getPDE(self):
        """
        returns the linear PDE to be solved for the Hessian Operator inverse
        :rtype: ``LinearPDE``
        """
        return self.__pde
     
    def getDualProduct(self, m, r):
        """
        returns the dual product of a gradient represented by X=r[1] and Y=r[0] with a level set function m:
             
             *Y_i*m_i + X_ij*m_{i,j}*
        
        :type m: ``esys.escript.Data``
        :type r: ``ArithmeticTuple``
        :rtype: float
        """
        A=0
        if not r[0].isEmpty(): A+=integrate(inner(r[0], m))
        if not r[1].isEmpty(): A+=integrate(inner(r[1], grad(m)))
	return A

    def getValue(self, m, grad_m):
        """
        return the value of the costfunction J
        
        
        
        :rtype: ``float``
        """
        mu=self.getWeights( uncompress=True)
        DIM=self.getDomain().getDim()
        numLS=self.getNumLevelSets()
                
        A=0
        n=0
        
        if self.__s0 is not None:
            A+=inner(integrate(m**2*self.__s0), mu[:numLS])
        n+=numLS
        
        if self.__s1 is not None:
	    if numLS == 1:
	        A+=integrate(inner(grad_m**2, self.__s1))*mu[n]
	    else:
	        for k in xrange(numLS):
		    A+=mu[n+k]*integrate(inner(grad_m[k,:]**2,self.__s1[k,:]))
        n+=numLS    
        
        if self.__sc is not None:
	      for k in xrange(numLS):
		   gk=grad_m[k,:]
		   len_gk=length(gk)
  		   for l in xrange(k):
		       gl=grad_m[l,:]
		       A+= mu[n+k*numLS+l] * integrate( self.__sc[l,k] * ( len_gk * length(gl) )**2 - inner(gk, gl)**2 ) 
        return A/2
        
    def getGradient(self, m,  grad_m):
        """
        returns the gradient of the costfunction J with respect to m. The function returns Y_k=dPsi/dm_k and
        X_kj=dPsi/dm_kj
        """
        
        mu=self.getWeights( uncompress=True)
        DIM=self.getDomain().getDim()
        numLS=self.getNumLevelSets()
        
        n=0
        
        if self.__s0 is not None:
	    Y = m * self.__s0 * mu[:numLS]
	else:
	    Y = Data()
        n+=numLS

        if self.__s1 is not None:
	    if numLS == 1:
	        X=grad_m* self.__s1*mu[n]
	    else:
	        X=self.getPDE().createCoefficient("X")
	        for k in xrange(numLS):
		    X[k,:]=mu[n+k]*grad_m[k,:]*self.__s1[k,:]
	else:
	    X = Data()
        n+=numLS 
        if self.__sc is not None:
           raise NotImplementedError
        return ArithmeticTuple(Y, X)
        
    def getArguments(self, m):
        """
        """
        return ( grad(m),)
         
    def getInverseHessianApproximation(self, m, r, grad_m):
        """
        """
        if self._new_mu or self._update_Hessian:
	    
	    self._new_mu=False
	    self._update_Hessian=False
	    
	    mu=self.getWeights( uncompress=True)
            DIM=self.getDomain().getDim()
            numLS=self.getNumLevelSets()
            n=0
            if self.__s0 is not None:
	        if numLS == 1:
	             D=self.__s0 * mu[n]
                else:	
                     D=self.getPDE().createCoefficient("D")
	             for k in xrange(numLS): D[k,k]=self.__s0[k] * mu[n+k]
	        self.getPDE().setValue(D=D)
	    n+=numLS

	    A=self.getPDE().createCoefficient("A")
            if self.__s1 is not None:
	       if numLS == 1:
		   for i in xrange(DIM): A[i,i]=self.__s1[i] * mu[n]
	       else:
	           for k in xrange(numLS): 
	                for i in xrange(DIM): A[k,i,k,i]=self.__s1[k,i] * mu[n+k]
            n+=numLS 
            if self.__sc is not None:
                raise NotImplementedError
            print A
            self.getPDE().setValue(A=A)
        self.getPDE().resetRightHandSideCoefficients()
        self.getPDE().setValue(X=r[1])
        print "X only: ",self.getPDE().getSolution()
        self.getPDE().resetRightHandSideCoefficients()
        self.getPDE().setValue(Y=r[0])
        print "Y only: ",self.getPDE().getSolution()
        
        self.getPDE().resetRightHandSideCoefficients()
        self.getPDE().setValue(X=r[1], Y=r[0])
        return self.getPDE().getSolution()
        
    def updateHessian(self):
        """
        notify the class to recalculate the Hessian operator
        """
	if not self.__useDiagonalHessianApproximation:
	    self._update_Hessian=True

   # ================ we should factor these out =============================================================
    def getNumRelevantTerms(self):
        """
        returns the number of terms in the costfunction of regularization with non-zero scaling factors
        :rtype: ``int``
        """
        return len(self.__weight_index)
    
    def getNumTerms(self):
        """
        returns the number of terms in the costfunction of regularization with non-zero scaling factors
        :rtype: ``int``
        """
        return len(self.__weight_index)
        
    def setWeights(self, mu=None):
        """
        sets the values for the weighting terms in the costsfunction. Note that only values
        correspond to entries with non-negative scaling factors are used.
        
        :param mu: weights 
        :type mu: ``list`` of ``float`` or ``np,array``
        """
        if mu == None:
	    mu = np.ones(self.getNumRelevantTerms())
	mu=np.asarray(mu)
	
	
	if len(mu) == self.getNumRelevantTerms():
	     if not mu.shape == (self.getNumRelevantTerms(),):
		raise ValueError("%s values are expected."%self.getNumRelevantTerms())
	     self.__mu=mu
	else:
	     if not mu.shape == (self.__total_num_weights,):
		raise ValueError("%s values are expected."%self.__total_num_weights)
	     self.__mu = np.zeros(self.getNumRelevantTerms())
	     for i in xrange(len(self.__weight_index)): self.__mu[i]=mu[self.__weight_index[i]]
	self._new_mu=True

    def setWeightsForS0(self, mu=None):
         """
         set the weights for s0-terms
         """
         numLS=self.getNumLevelSets()
         my_mu=self.getWeights(uncompress=True)
         if mu is None:
             my_mu[:numLS]=1
         else:
	     my_mu[:numLS]=mu
         self.setWeights(my_mu)
	   
    def setWeightsForS1(self, mu=None):
         """
         set the weights for s1-terms
         """
         numLS=self.getNumLevelSets()
         my_mu=self.getWeights(uncompress=True)
         if mu is None:
              my_mu[numLS:2*numLS]=1
         else:
              my_mu[numLS:2*numLS]=mu
         self.setWeights(my_mu)
      
    def setWeightsForSc(self, mu): 
         """
         set the weights for s1-terms
         """
         numLS=self.getNumLevelSets()
         my_mu=self.getWeights(uncompress=True)
         if mu is None:
              my_mu[2*numLS:]=1
         else:
              my_mu[2*numLS:]=mu
         self.setWeights(my_mu)

    
    def getWeights(self, uncompress=False):
        """
        Returns the weights for the terms in the costsfunction. 
        The first ``numLevelSets`` values used for the 
        regularization terms and the remaining values for the cross correlation terms.
 
       :type mu: ``list`` of ``float``
        """
        if uncompress:
	     mu = np.zeros(self.__total_num_weights) 
	     for i in xrange(len(self.__weight_index)): mu[self.__weight_index[i]] = self.__mu[i]
	     return mu
	else:
           return self.__mu
    
    def getWeightIndex(self):
        """
        returns an iondex to the contributions of terms with non-zero scaling factor.
        """
        return self.__weight_index
	 
