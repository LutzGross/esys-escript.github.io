
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
from esys.escript import ReducedFunction, outer, Data, Scalar, grad, inner, integrate, interpolate, kronecker, boundingBoxEdgeLengths, vol, sqrt, length
from esys.escript.linearPDEs import LinearPDE, IllegalCoefficientValue
from esys.escript.pdetools import ArithmeticTuple

class Regularization(CostFunction):
    """
    The regularization term for the level set function ``m`` within the cost
    function J for an inversion:

    *J(m)=1/2 * sum_k imtegrate( mu[k] * ( w0[k] * m_k**2 * w1[k,i] * m_{k,i}**2) + sum_l<k mu_c[l,k] wc[l,k] * | curl(m_k) x curl(m_l) |^2*

    where w0[k], w1[k,i] and  wc[k,l] are non-negative weighting factors and
    mu[k] and mu_c[l,k] are trade-off factors which may be altered
    during the inversion. The weighting factors are normalized such that their
    integrals over the domain are constant:

    *integrate(w0[k] + inner(w1[k,:],1/L[:]**2))=scale[k]* volume(domain)*
    *integrate(wc[l,k]*1/L**4)=scale_c[k]* volume(domain) *

    """
    def __init__(self, domain, numLevelSets=1,
                       w0=None, w1=None, wc=None,
                       location_of_set_m=Data(),
                       useDiagonalHessianApproximation=False, tol=1e-8,
                       scale=None, scale_c=None):
        """
        initialization.

        :param domain: domain
        :type domain: `Domain`
        :param numLevelSets: number of level sets
        :type numLevelSets: ``int``
        :param w0: weighting factor for the m**2 term. If not set zero is assumed.
        :type w0: ``Scalar`` if ``numLevelSets`` == 1 or `Data` object of shape
                  (``numLevelSets`` ,) if ``numLevelSets`` > 1
        :param w1: weighting factor for the grad(m_i) terms. If not set zero is assumed
        :type w1: ``Vector`` if ``numLevelSets`` == 1 or `Data` object of shape
                  (``numLevelSets`` , DIM) if ``numLevelSets`` > 1
        :param wc: weighting factor for the cross gradient terms. If not set
                   zero is assumed. Used for the case if ``numLevelSets`` > 1
                   only. Only values ``wc[l,k]`` in the lower triangle (l<k)
                   are used.
        :type wc: `Data` object of shape (``numLevelSets`` , ``numLevelSets``)
        :param location_of_set_m: marks location of zero values of the level
                                  set function ``m`` by a positive entry.
        :type location_of_set_m: ``Scalar`` if ``numLevelSets`` == 1 or `Data`
                object of shape (``numLevelSets`` ,) if ``numLevelSets`` > 1
        :param useDiagonalHessianApproximation: if True cross gradient terms
                    between level set components are ignored when calculating
                    approximations of the inverse of the Hessian Operator.
                    This can speed-up the calculation of the inverse but may
                    lead to an increase of the number of iteration steps in the
                    inversion.
        :type useDiagonalHessianApproximation: ``bool``
        :param tol: tolerance when solving the PDE for the inverse of the
                    Hessian Operator
        :type tol: positive ``float``
        
        :param scale: weighting factor for level set function variation terms. If not set one is used. 
        :type scale: ``Scalar`` if ``numLevelSets`` == 1 or `Data` object of shape
                  (``numLevelSets`` ,) if ``numLevelSets`` > 1
        :param scale_c: scale for the cross gradient terms. If not set
                   one is assumed. Used for the case if ``numLevelSets`` > 1
                   only. Only values ``scale_c[l,k]`` in the lower triangle (l<k)
                   are used.
        :type scale_c: `Data` object of shape (``numLevelSets`` , ``numLevelSets``)
        

        """
        if w0 == None and w1==None:
              raise ValueError("Values for w0 or for w1 must be given.")
	if wc == None and  numLevelSets>1:
	      raise ValueError("Values for wc must be given.")

        self.__domain=domain
        DIM=self.__domain.getDim()
        self.__numLevelSets=numLevelSets
        
        self.__pde=LinearPDE(self.__domain, numEquations=self.__numLevelSets)
        self.__pde.getSolverOptions().setTolerance(tol)
        self.__pde.setSymmetryOn()
        try:
            self.__pde.setValue(q=location_of_set_m)
        except IllegalCoefficientValue:
            raise ValueError("Unable to set location of fixed level set function.")
	  
        # =========== check the shape of the scales: =================================
        if scale is None:
	    if numLevelSets == 1 :
	       scale = 1.
	    else:
	       scale = np.ones((numLevelSets,))
	else:
	    scale=np.asarray(scale)
	    if numLevelSets == 1 :
	        if scale.shape == ():
		   if not scale > 0 :
		      raise ValueError("Value for scale must be positive.")
		else:
		   raise ValueError("Unexpected shape %s for scale."%scale.shape)
            else:
	         if scale.shape is (numLevelSets,):
		     if not min(scale) > 0:
		        raise ValueError("All value for scale must be positive.")
		 else:
		   raise ValueError("Unexpected shape %s for scale."%scale.shape)
	
        if scale_c is None or numLevelSets < 2:
	    scale_c = np.ones((numLevelSets,numLevelSets))
	else:
	    scale_c=np.asarray(scale_c)
	    if scale_c.shape == (numLevelSets,numLevelSets):
	        if not all( [ [ scale_c[l,k] > 0. for l in xrange(k) ] for k in xrange(1,numLevelSets) ]):
		        raise ValueError("All values in the lower triangle of scale_c must be positive.")
            else:
		 raise ValueError("Unexpected shape %s for scale."%scale_c.shape)
	# ===== check the shape of the weights: ============================================
        if w0 is not None:
	      w0 = interpolate(w0,self.__pde.getFunctionSpaceForCoefficient('D'))
	      s0=w0.getShape()
	      if numLevelSets == 1 :
		   if  not s0 == () :
		      raise ValueError("Unexpected shape %s for weight w0."%s0)
              else:
		   if not s0 == (numLevelSets,):
		      raise ValueError("Unexpected shape %s for weight w0."%s0) 
	if not w1 is None:
	      w1 = interpolate(w1,self.__pde.getFunctionSpaceForCoefficient('A'))
	      s1=w1.getShape()
	      if numLevelSets is 1 :
		   if not s1 == (DIM,) :
		      raise ValueError("Unexpected shape %s for weight w1."%s1)
              else:
		   if not s1 == (numLevelSets,DIM):
		      raise ValueError("Unexpected shape %s for weight w1."%s1) 
        if numLevelSets == 1 :
             wc=None
        else:
             wc = interpolate(wc,self.__pde.getFunctionSpaceForCoefficient('A'))
             sc=wc.getShape()
             if not sc == (numLevelSets, numLevelSets):
                raise ValueError("Unexpected shape %s for weight wc."%(sc,))
        # ============= now we rescale weights: ======================================
        L2s=np.asarray(boundingBoxEdgeLengths(domain))**2
        L4=1/np.sum(1/L2s)**2
        if numLevelSets == 1 : 
            A=0
            if w0 is not None:
	        A = integrate(w0)
	    if w1 is not None:
	        A += integrate(inner(w1, 1/L2s))
	    if A > 0:
	        f = scale/A
	        if w0 is not None:
	             w0*=f
	        if w1 is not None:
		     w1*=f
            else:
	       raise ValueError("Non-positive weighting factor detected.") 
        else:

	     for k in xrange(numLevelSets):
	         A=0
                 if w0 is not None:
	             A = integrate(w0[k])
	         if w1 is not None:
	              A += integrate(inner(w1[k,:], 1/L2s))
	         if A > 0:
	              f = scale[k]/A
	              if w0 is not None:
	                 w0[k]*=f
	              if w1 is not None:
		         w1[k,:]*=f
                 else:
	           raise ValueError("Non-positive weighting factor for level set component %d detected."%k) 
		 
	         # and now the cross-gradient:
	         if wc is not None:
	             for l in xrange(k):   
	                A = integrate(wc[l,k])/L4
  	                if A > 0:
 	                   f = scale_c[l,k]/A
 	                   wc[l,k]*=f
#                        else:
#	                   raise ValueError("Non-positive weighting factor for cross-gradient level set components %d and %d detected."%(l,k)) 
	            
        self.__w0=w0
        self.__w1=w1
        self.__wc=wc

        self.__pde_is_set=False        
        if self.__numLevelSets > 1:
            self.__useDiagonalHessianApproximation=useDiagonalHessianApproximation
        else:
            self.__useDiagonalHessianApproximation=True
        self._update_Hessian=True

        self.__num_tradeoff_factors=numLevelSets+((numLevelSets-1)*numLevelSets)/2
        self.setTradeOffFactors()
        self.__vol_d=vol(self.__domain)
        
    def getDomain(self):
        """
        returns the domain of the regularization term

        :rtype: ``Domain``
        """
        return self.__domain

    def getNumLevelSets(self):
        """
        returns the number of level set functions

        :rtype: ``int``
        """
        return self.__numLevelSets

    def getPDE(self):
        """
        returns the linear PDE to be solved for the Hessian Operator inverse

        :rtype: `LinearPDE`
        """
        return self.__pde
    
    def getDualProduct(self, m, r):
        """
        returns the dual product of a gradient represented by X=r[1] and Y=r[0]
        with a level set function m:

             *Y_i*m_i + X_ij*m_{i,j}*

        :type m: `Data`
        :type r: `ArithmeticTuple`
        :rtype: ``float``
        """
        A=0
        if not r[0].isEmpty(): A+=integrate(inner(r[0], m))
        if not r[1].isEmpty(): A+=integrate(inner(r[1], grad(m)))
        return A
    def getNumTradeOffFactors(self):
        """
        returns the number of trade-off factors being used.

        :rtype: ``int``
        """
        return self.__num_tradeoff_factors

    def setTradeOffFactors(self, mu=None):
        """
        sets the trade-off factors for the level-set variation and the cross-gradient
        
        :param mu: new values for the trade-off factors where values mu[:numLevelSets] are the 
                   trade-off factors for the level-set variation and the remaining values for 
                   the cross-gradient part with mu_c[l,k]=mu[numLevelSets+l+((k-1)*k)/2] (l<k).
                   If no values for mu is given ones are used. Values must be positive.
        :type mu: ``list`` of ``float`` or ```numpy.array```
        """
        numLS=self.getNumLevelSets()
        numTF=self.getNumTradeOffFactors()
        if mu is None:
	   mu = np.ones((numTF,))
	else:
	   mu = np.asarray(mu)

	if mu.shape == (numTF,):
	    self.setTradeOffFactorsForVariation(mu[:numLS])
	    mu_c2=np.zeros((numLS,numLS))
	    for k in xrange(numLS):
	       for l in xrange(k):
		   mu_c2[l,k] = mu[numLS+l+((k-1)*k)/2] 
	    self.setTradeOffFactorsForCrossGradient(mu_c2)
	elif mu.shape == () and numLS ==1:
	    self.setTradeOffFactorsForVariation(mu)
	else:
	   raise ValueError("Unexpected shape %s for mu."%(mu.shape,)) 
	   
    def setTradeOffFactorsForVariation(self, mu=None):
         """
         sets the trade-off factors for the level-set variation part
         
         :param mu:  new values for the trade-off factors. Values must be positive.
         :type mu: `float``, ``list`` of ``float`` or ```numpy.array```
         """
         numLS=self.getNumLevelSets()
         if mu is None:
	    if numLS == 1:
	       mu = 1.
	    else:
	       mu = np.ones((numLS,))

	 mu=np.asarray(mu)
	 if numLS == 1:
	   if mu.shape == (1,): mu=mu[0] 
	   if mu.shape == ():
	      if mu > 0:
		 self.__mu= mu
		 self._new_mu=True
	      else:
		 raise ValueError("Value for trade-off factor must be positive.") 
	   else:
	      raise ValueError("Unexpected shape %s for mu."%mu.shape)
	 else:
	   if mu.shape == (numLS,):
	       if min(mu) > 0:
		   self.__mu= mu
		   self._new_mu=True
	       else:
		   raise ValueError("All value for mu must be positive.")
           else:
	       raise ValueError("Unexpected shape %s for trade-off factor."%mu.shape) 

    def setTradeOffFactorsForCrossGradient(self, mu_c=None):
        """
        sets the trade-off factors for the cross--gradient terms
         
        :param mu_c:  new values for the trade-off factors for the cross--gradient terms. Values must be positive.
                      if now value is given ones are used. Onky value mu_c[l,k] for l<k are used.
        :type mu: `float``, ``list`` of ``float`` or ```numpy.array```
         
        """
        numLS=self.getNumLevelSets()
        if mu_c is None or numLS < 2:
	    self.__mu_c = np.ones((numLS,numLS))
	if isinstance(mu_c, float) or isinstance(mu_c, int):
	    self.__mu_c = np.zeros((numLS,numLS))
	    self.__mu_c[:,:]=mu_c 
	else:
	    mu_c=np.asarray(mu_c)
	    if mu_c.shape == (numLS,numLS):
	        if not all( [ [ mu_c[l,k] > 0. for l in xrange(k) ] for k in xrange(1,numLS) ]):
		     raise ValueError("All trade-off factors in the lower triangle of mu_c must be positive.")
		else:
		     self.__mu_c =  mu_c
		     self._new_mu=True
            else:
		 raise ValueError("Unexpected shape %s for mu."%mu_c.shape)
    
    def getArguments(self, m):
        """
        """
        return ( grad(m),)
		 

    def getValue(self, m, grad_m):
        """
        returns the value of the cost function J with respect to m.

        :rtype: ``float``
        """
        mu=self.__mu
        mu_c=self.__mu_c
        DIM=self.getDomain().getDim()
        numLS=self.getNumLevelSets()

        A=0
        if self.__w0 is not None:
           A+=inner(integrate(m**2 * self.__w0), mu)

        if self.__w1 is not None:
            if numLS == 1:
                A+=integrate(inner(grad_m**2, self.__w1))*mu
            else:
                for k in xrange(numLS):
                    A+=mu[k]*integrate(inner(grad_m[k,:]**2,self.__w1[k,:]))
                    
        if numLS > 1:
            for k in xrange(numLS):
                gk=grad_m[k,:]
                len_gk=length(gk)
                for l in xrange(k):
                    gl=grad_m[l,:]
                    A+= mu_c[l,k] * integrate( self.__wc[l,k] * ( len_gk * length(gl) )**2 - inner(gk, gl)**2 )
        return A/2

    def getGradient(self, m,  grad_m):
        """
        returns the gradient of the cost function J with respect to m.
        The function returns Y_k=dPsi/dm_k and X_kj=dPsi/dm_kj
        """

        mu=self.__mu
        mu_c=self.__mu_c
        DIM=self.getDomain().getDim()
        numLS=self.getNumLevelSets()

        print "WARNING: WRONG FUNCTION SPACE"
        
        grad_m=grad(m, ReducedFunction(m.getDomain()))
        if self.__w0 is not None:
            Y = m * self.__w0 * mu
        else:
	    if numLS == 1:
	      Y = Scalar(0,  grad_m.getFunctionSpace())
	    else:
              Y = Data(0, (numLS,) , grad_m.getFunctionSpace())

        if self.__w1 is not None:
	    X=grad_m*self.__w1
            if numLS == 1:
                X=grad_m* self.__w1*mu
            else:
                for k in xrange(numLS):
                    X[k,:]*=mu[k]
        else:
            X = Data(0, grad_m.getShape(), grad_m.getFunctionSpace())
        
        # cross gradient terms:
        if numLS > 1:	  
	  for  k in xrange(numLS):
       	     grad_m_k=grad_m[k,:]
  	     l2_grad_m_k = length(grad_m_k)**2
	     for  l in xrange(k):
	       grad_m_l=grad_m[l,:]
	       l2_grad_m_l = length(grad_m_l)**2
	       grad_m_lk = inner(grad_m_l, grad_m_k)
	       f=  mu_c[l,k]* self.__wc[l,k]
	       X[l,:] += f * ( l2_grad_m_l *  grad_m_l - grad_m_lk * grad_m_k )
	       X[k,:] += f * ( l2_grad_m_k *  grad_m_k - grad_m_lk * grad_m_l )
        
        return ArithmeticTuple(Y, X)



    def getInverseHessianApproximation(self, m, r, grad_m):
        """
        """
        if self._new_mu or self._update_Hessian:

            self._new_mu=False
            self._update_Hessian=False
            mu=self.__mu
            mu_c=self.__mu_c

            DIM=self.getDomain().getDim()
            numLS=self.getNumLevelSets()
            if self.__w0 is not None:
                if numLS == 1:
                     D=self.__w0 * mu
                else:
                     D=self.getPDE().createCoefficient("D")
                     for k in xrange(numLS): D[k,k]=self.__w0[k] * mu[k]
                self.getPDE().setValue(D=D)

            A=self.getPDE().createCoefficient("A")
            if self.__w1 is not None:
               if numLS == 1:
                   for i in xrange(DIM): A[i,i]=self.__w1[i] * mu
               else:
                   for k in xrange(numLS):
                        for i in xrange(DIM): A[k,i,k,i]=self.__w1[k,i] * mu[k]
                        
            if numLS > 1:
	       for  k in xrange(numLS):
       	         grad_m_k=grad_m[k,:]
  	         l2_grad_m_k = length(grad_m_k)**2
  	         o_kk=outer(grad_m_k, grad_m_k)
	         for  l in xrange(k):
	            grad_m_l=grad_m[l,:]
	            l2_grad_m_l = length(grad_m_l)**2
	            i_lk = inner(grad_m_l, grad_m_k)
	            o_lk = outer(grad_m_l, grad_m_k)
	            o_kl = outer(grad_m_k, grad_m_l)
	            o_ll=outer(grad_m_l, grad_m_l)
	            f=  mu_c[l,k]* self.__wc[l,k]
	            
	            A[l,:,l,:] += f * ( l2_grad_m_k * kronecker(DIM) - o_kk )
	            A[l,:,k,:] += f * ( 2 * o_lk -   o_kl - i_lk * kronecker(DIM) )
	            A[k,:,l,:] += f * ( 2 * o_kl -   o_lk - i_lk * kronecker(DIM) )
	            A[k,:,k,:] += f * ( l2_grad_m_l * kronecker(DIM) - o_ll )
            self.getPDE().setValue(A=A)
        #self.getPDE().resetRightHandSideCoefficients()
        #self.getPDE().setValue(X=r[1])
        #print "X only: ",self.getPDE().getSolution()
        #self.getPDE().resetRightHandSideCoefficients()
        #self.getPDE().setValue(Y=r[0])
        #print "Y only: ",self.getPDE().getSolution()

        self.getPDE().resetRightHandSideCoefficients()
        self.getPDE().setValue(X=r[1], Y=r[0])
        return self.getPDE().getSolution()

    def updateHessian(self):
        """
        notify the class to recalculate the Hessian operator
        """
        if not self.__useDiagonalHessianApproximation:
            self._update_Hessian=True
            
    def getNorm(self, m):
        """
        returns the norm of ``m``

        :param m: level set function
        :type m: `Data`
        :rtype: ``float``
        """
        return sqrt(integrate(length(m)**2)/self.__vol_d)