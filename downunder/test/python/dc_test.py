from esys.downunder import *
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
import numpy as np
from esys.escript import *
from esys.weipa import saveSilo
from esys.escript.linearPDEs import LinearSinglePDE, LinearPDE


from esys.escript.pdetools import Locator



mpisize = getMPISizeWorld()
class TestDCResistivity(unittest.TestCase):

    def test_PDE2D(self):
         
        dx_tests=0.1 
        electrodes=[(0.5-2*dx_tests,1.), (0.5-dx_tests,1.), (0.5+dx_tests,1.), (0.5+2*dx_tests,1.)]
        from esys.finley import Rectangle
        domain=Rectangle(20,20, d1=mpisize,  diracPoints=electrodes, diracTags=["sl0", "sl1", "sr0", "sr1"] )
        loc=Locator(domain,electrodes[2:])

        # this creates some reference Data:
        x=domain.getX()
        q=whereZero(x[0]-inf(x[0]))+whereZero(x[0]-sup(x[0]))+whereZero(x[1]-inf(x[1]))
        ppde=LinearPDE(domain, numEquations=1)
        s=Scalar(0.,DiracDeltaFunctions(domain))
        s.setTaggedValue("sl0" ,1.)
        s.setTaggedValue("sl1",-1.)
        ppde.setValue(A=kronecker(2), q=q, y_dirac=s)
        pp=ppde.getSolution()
        uu=loc(pp)
        
        # arguments for DcRes
        sigmaPrimary=5
        current = 10        
        delphi_in = [ (uu[1]-uu[0]) * current ]
        sourceInfo = [ "sl0",  "sl1" ]
        sampleTags = [ "-sr0", "sr1" ]


        
        
        acw=DcRes(domain, loc, delphi_in, sourceInfo, current, sampleTags,sigmaPrimary, w=1., coordinates=None, tol=1e-8,saveMemory=True,b=None)
        
        # this should be the primary potential:
        p_primary = pp * current/sigmaPrimary
        print p_primary 
        print acw.u_primary
        
        self.assertTrue(Lsup(p_primary-acw.u_primary) < 1.e-6 * Lsup(p_primary) )    
#        :param domain: the domain of the model
#        :type: escript domain
#        :param locator: should be setup to contain the measurement pairs
#        :type: escript locator
#        :param: delphi_in: this is v_pq, the potential difference for the current source  and a set of measurement pairs. a list of measured potential differences is expected. Note this should be the secondary potential only
#        :type delphi_in: tuple
#        :param sourceInfo: descibes the current electode setup. a pair of tags should be provided for the current source setup. the first tag will be set to current and the second tag to -current
#        :type sourceInfo: tuple
#        :param sigmaPrimary: the conductivity to be used for the Primary Solution.
#        :type sigmaPrimary: ``Data`` of shape (1,)
        print "uu =", uu
        print acw.u_primary
        print pp
        
        P0=10. # matches current        
        args0=acw.getArguments(P0)
        print args0
        p=args0[0]
        u=args0[1]
        print u
        print p
        
        #self.assertTrue(Lsup() < 1.e-8)
        #self.assertTrue(Lsup(u[1]) < 1.e-8)
        #self.assertTrue(Lsup(u[2]-2.5*domain.getX()[2]) < 1.e-8)
        
        dd=acw.getDefect(P0, *args0)
        
        self.assertTrue( dd >= 0.)
        self.assertTrue( dd <= 1e-7 * 2.5 )
    def rtest_Differential2D(self):
    
        lam=2.
        mu=1.
        
        INC=0.01 
        from esys.ripley import Brick
        domain=Brick(20,20,20*mpisize-1 , d2=mpisize)
        
        xb=FunctionOnBoundary(domain).getX()
        m=whereZero(xb[2]-1)
        w=m*[0,0,1]
        d=m*2.5
        acw=DCResistivity(domain, w,d, lam, mu )
        
        
        x=Function(domain).getX()
        P0=x[0]*x[1]
        args0=acw.getArguments(P0)
        d0=acw.getDefect(P0, *args0)
        grad_d=acw.getGradient(P0, *args0)

        
        dP=exp(-(length(x-[0.5,0.5,0.5])/0.06)**2)
        P1=P0+INC*dP
        args1=acw.getArguments(P1)
        d1=acw.getDefect(P1, *args1)
        ref=abs((d1-d0)/INC)
        self.assertTrue(abs((d1-d0)/INC-integrate(grad_d* dP)) < ref * 1.e-5) 

        dP=exp(-(length(x-[0.3,0.3,0.5])/0.06)**2)
        P2=P0-INC*dP
        args2=acw.getArguments(P2)
        d2=acw.getDefect(P2, *args2)
        ref=abs((d2-d0)/INC)
        self.assertTrue(abs((d2-d0)/INC+integrate(grad_d* dP)) < ref * 1.e-5) 

################################
if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
    
    
