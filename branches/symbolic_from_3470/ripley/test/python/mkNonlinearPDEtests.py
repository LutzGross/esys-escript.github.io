print "import unittest"
print "from esys.escript import *"
print "DEBUG=False"


print "from esys.escript.linearPDEs import IllegalCoefficient, IllegalCoefficient"
print "class Test_NonlinearPDE(unittest.TestCase):"
for nsol in [1, 2, 5 ] :
   for c in [ "X", "X_reduced", "Y", "Y_reduced", "y", "y_reduced", "y_contact" ,"y_contact_reduced" ]:
        for j in xrange(nsol):
	   if c.startswith("X"):
	          r=xrange(nsol)
	   else:
	          r=xrange(j,j+1)
	   for j2 in r:
	      if  c.startswith("X"):
                  print "   def test_NonLinearPDE_Unknown%s_%s_Component_%s_%s(self):"%(nsol, c, j, j2)
              else:
		  print "   def test_NonLinearPDE_Unknown%s_%s_Component_%s(self):"%(nsol, c, j)
              print "       DIM = self.domain.getDim()"

              if nsol == 1 :
                 print "       u=Symbol('v', (), dim=DIM)"
              else:
		 print "       u=Symbol('v', (%s,), dim=DIM)"%nsol
              print "       pde=NonlinearPDE(self.domain, u, debug=DEBUG)"
              print "       self.assertEqual(pde.getSolutionSymbol(), u, \"solution symbol is wrong.\")"
              print "       g=grad(u)"
              print "       f=pde.createCoefficient(\"%s\")"%c
              if nsol == 1 :
		 if c.startswith("X"):
		     print "       self.assertEqual(f.getShape(), (DIM,), \"wrong shape of coefficient %s.\")"%c
		     print "       f[:]=g*u"
		 else:
		     print "       self.assertEqual(f.getShape(), () ,\"wrong shape of coefficient %s.\")"%c
		     if c.startswith("Y"):
		         print "       f=length(g)**2*u"
		     else:
		         print "       f=u"
              else:
		  if c.startswith("X"):
		     print "       self.assertEqual(f.getShape(), (%s,DIM), \"wrong shape of coefficient %s.\")"%(nsol,c)
		     print "       f[%s,:]=g[%s,:]*length(u)**2"%(j,j2)
		  else:
		     print "       self.assertEqual(f.getShape(), (%s,), \"wrong shape of coefficient %s.\")"%(nsol,c)
		     if c.startswith("Y"):
		        print "       f[%s]=length(g)**2*length(u)**2"%j
		     else:
    		        print "       f[%s]=length(u)**2"%j

	      print "       pde.setValue(%s=f)"%c
	      if not c == "X" :
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"X\")"
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"A\")"
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"B\")"
              if not c == "Y" :	         
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"Y\")"
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"D\")"
       	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"C\")"
              if not c == "y" :	         
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"y\")"	         
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"d\")"
              if not c == "y_contact" :	         
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"y_contact\")"
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"d_contact\")"
	      if not c == "X_reduced" :
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"X_reduced\")"
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"A_reduced\")"
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"B_reduced\")"
              if not c == "Y_reduced" :	         
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"Y_reduced\")"
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"D_reduced\")"
       	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"C_reduced\")"
              if not c == "y_reduced" :	         
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"y_reduced\")"	         
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"d_reduced\")"
              if not c == "y_contact_reduced" :	         
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"y_contact_reduced\")"
	         print "       self.assertRaises(IllegalCoefficient, pde.getCoefficient, \"d_contact_reduced\")"
	         
	      if c.endswith("reduced"):
		   postfix="_reduced"
	      else:
		   postfix=""

		 
	      if c.startswith("y"):
		    print "       Y=pde.getCoefficient(\"%s\")"%c
		    print "       self.assertEqual(Y, f, \"wrong coefficient %s.\")"%c
		    print "       D=pde.getCoefficient(\"d%s\")"%c[1:]
		    if nsol>1:
                        print "       D_ref=Symbol('D_ref', (%s, %s), dim=DIM)"%(nsol, nsol)
                        for k0 in xrange(nsol):
			    for k1 in xrange(nsol):
			      if k0 == j : 
			         print "       D_ref[%s,%s]=2*u[%s]"%(k0,k1, k1)
			      else:
			         print "       D_ref[%s,%s]=0"%(k0,k1)
                    else:
                        print "       D_ref=Symbol('D_ref', (), dim=DIM)*0+1"
                        
                    print "       self.assertEqual(D_ref.simplify(), D.simplify(), \"wrong coefficient D%s.\")"%c[1:]
	      elif c.startswith("Y"): 
		    print "       Y=pde.getCoefficient(\"%s\")"%c
		    print "       self.assertEqual(Y, f, \"wrong coefficient %s.\")"%c
                    print "       D=pde.getCoefficient(\"D%s\")"%c[1:]
                    print "       C=pde.getCoefficient(\"C%s\")"%c[1:]
                    if nsol>1:
                        print "       D_ref=Symbol('D_ref', (%s,%s), dim=DIM)"%(nsol,nsol)
                        print "       C_ref=Symbol('C_ref', (%s,%s,DIM), dim=DIM)"%(nsol,nsol)
                        for k0 in xrange(nsol):
			    for k1 in xrange(nsol):
			      if k0 == j : 
			         print "       D_ref[%s,%s]=2*length(g)**2*u[%s]"%(k0,k1,k1)
			         print "       C_ref[%s,%s,:]=2*length(u)**2*grad(u)[%s,:]"%(k0,k1, k1)
			      else:
			         print "       D_ref[%s,%s]=0"%(k0,k1)
			         print "       C_ref[%s,%s,:]=0"%(k0,k1)
                    else:
                        print "       D_ref=length(g)**2"
                        print "       C_ref=2*g*u"
                    print "       self.assertEqual(D_ref.simplify(), D.simplify(), \"wrong coefficient D%s.\")"%c[1:]
                    print "       self.assertEqual(C_ref.simplify(), C.simplify(), \"wrong coefficient C%s.\")"%c[1:]
	      else: 
		    print "       X=pde.getCoefficient(\"%s\")"%c
		    print "       self.assertEqual(X, f, \"wrong coefficient %s.\")"%c
                    print "       A=pde.getCoefficient(\"A%s\")"%c[1:]
                    print "       B=pde.getCoefficient(\"B%s\")"%c[1:]
                    if nsol>1:
                        print "       A_ref=Symbol('A_ref', (%s,DIM,%s,DIM), dim=DIM)"%(nsol,nsol)
                        print "       B_ref=Symbol('B_ref', (%s,DIM,%s), dim=DIM)"%(nsol,nsol)
                        for k0 in xrange(nsol):
			    for k1 in xrange(nsol):
			      if k0 == j : 
			         if k1 == j2:
			             print "       A_ref[%s,:,%s,:] = length(u)**2 *kronecker(DIM)"%(k0,k1)
			         else:
				     print "       A_ref[%s,:,%s,:] = 0"%(k0,k1)
			         print "       B_ref[%s,:,%s]=2*grad(u)[%s,:]*u[%s]"%(k0,k1, j2, k1)
			      else:
			         print "       A_ref[%s,:,%s,:]=0"%(k0,k1)
			         print "       B_ref[%s,:, %s]=0"%(k0,k1)
                    else:
                        print "       A_ref=u*kronecker(DIM)"
                        print "       B_ref=g"
                    print "       self.assertEqual(B_ref.simplify(), B.simplify(), \"wrong coefficient B%s.\")"%c[1:]    
                    print "       self.assertEqual(A_ref.simplify(), A.simplify(), \"wrong coefficient A%s.\")"%c[1:]

              #print "       u0=Data(0., u.getShape(), ContinuousFunction(dom))
              #print "       ui=pde.getSolution(u=u0, **params)