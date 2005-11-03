# $Id$

"""
program generates parts of the util.py and the test_util.py script
"""
u_prog ="import esys.escript\n"
u_prog+="import numarray\n"
u_prog+="import math\n"
t_prog=""

func= [ ["atan", [-100.,100.], "math.atan(%a%)"],
        ["log", [1.e-8,10.], "math.log(%a%)"]
      ]
import random 
import numarray
import math
finc=1.e-6

def makeArray(shape,rng):
   l=rng[1]-rng[0]
   out=numarray.zeros(shape,numarray.Float64)
   if len(shape)==0:
       out=random.random()
   elif len(shape)==1:
       for i0 in range(shape[0]):
                   out[i0]=l*random.random()+rng[0]
   elif len(shape)==2:
       for i0 in range(shape[0]):
          for i1 in range(shape[1]):
                   out[i0,i1]=l*random.random()+rng[0]
   elif len(shape)==3:
       for i0 in range(shape[0]):
          for i1 in range(shape[1]):
             for i2 in range(shape[2]):
                   out[i0,i1,i2]=l*random.random()+rng[0]
   elif len(shape)==4:
       for i0 in range(shape[0]):
          for i1 in range(shape[1]):
             for i2 in range(shape[2]):
                for i3 in range(shape[3]):
                   out[i0,i1,i2,i3]=l*random.random()+rng[0]
   else:
       raise SystemError,"rank is restricted to 4"
   return out         

def makeResult2(val,case):
   if isinstance(val,float):
       out=val
   elif len(val.shape)==0:
       out=val[0]
   elif len(val.shape)==1:
       out=val[0]
       for i0 in range(val.shape[0]):
           if case == "Lsup":
               out=max(out,abs(val[i0]))
           elif case == "sup":
               out=max(out,val[i0])
           else:
               out=min(out,val[i0])
   elif len(val.shape)==2:
       out=val[0,0]
       for i0 in range(val.shape[0]):
          for i1 in range(val.shape[1]):
               if case == "Lsup":
                   out=max(out,abs(val[i0,i1]))
               elif case == "sup":
                   out=max(out,val[i0,i1])
               else:
                   out=min(out,val[i0,i1])
   elif len(val.shape)==3:
       out=val[0,0,0]
       for i0 in range(val.shape[0]):
          for i1 in range(val.shape[1]):
             for i2 in range(val.shape[2]):
               if case == "Lsup":
                   out=max(out,abs(val[i0,i1,i2]))
               elif case == "sup":
                   out=max(out,val[i0,i1,i2])
               else:
                   out=min(out,val[i0,i1,i2])
   elif len(val.shape)==4:
       out=val[0,0,0,0]
       for i0 in range(val.shape[0]):
          for i1 in range(val.shape[1]):
             for i2 in range(val.shape[2]):
                for i3 in range(val.shape[3]):
                   if case == "Lsup":
                       out=max(out,abs(val[i0,i1,i2,i3]))
                   elif case == "sup":
                       out=max(out,val[i0,i1,i2,i3])
                   else:
                       out=min(out,val[i0,i1,i2,i3])
   else:
       raise SystemError,"rank is restricted to 4"
   return out         

def makeResult(val,f):
   if isinstance(val,float):
       out=eval(f[2].replace("%a%","val"))
   elif len(val.shape)==0:
       out=eval(f[2].replace("%a%","val"))
   elif len(val.shape)==1:
       out=numarray.zeros(val.shape,numarray.Float64)
       for i0 in range(val.shape[0]):
                   out[i0]=eval(f[2].replace("%a%","val[i0]"))
   elif len(val.shape)==2:
       out=numarray.zeros(val.shape,numarray.Float64)
       for i0 in range(val.shape[0]):
          for i1 in range(val.shape[1]):
                   out[i0,i1]=eval(f[2].replace("%a%","val[i0,i1]"))
   elif len(val.shape)==3:
       out=numarray.zeros(val.shape,numarray.Float64)
       for i0 in range(val.shape[0]):
          for i1 in range(val.shape[1]):
             for i2 in range(val.shape[2]):
                   out[i0,i1,i2]=eval(f[2].replace("%a%","val[i0,i1,i2]"))
   elif len(val.shape)==4:
       out=numarray.zeros(val.shape,numarray.Float64)
       for i0 in range(val.shape[0]):
          for i1 in range(val.shape[1]):
             for i2 in range(val.shape[2]):
                for i3 in range(val.shape[3]):
                   out[i0,i1,i2,i3]=eval(f[2].replace("%a%","val[i0,i1,i2,i3]"))
   else:
       raise SystemError,"rank is restricted to 4"
   return out         
    	
for case in ["Lsup", "sup", "inf"]:
   for args in ["float","array","constData","taggedData","expandedData"]: 
     for sh in [ (),(2,), (4,5), (6,2,2),(3,2,3,4)]:
       if not args=="float" or len(sh)==0:
         t_prog+="   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
         tname="def test_%s_%s_rank%s"%(case,args,len(sh))
         t_prog+="   %s(self):\n"%tname
         if args in ["float","array" ]:
           a=makeArray(sh,[-1,1])
           r=makeResult2(a,case)
           if len(sh)==0:
              t_prog+="      arg=%s\n"%a
           else:
              t_prog+="      arg=numarray.array(%s)\n"%a.tolist()
           t_prog+="      ref=%s\n"%r
           t_prog+="      res=%s(arg)\n"%case
           t_prog+="      self.failUnless(isinstance(res,float),\"wrong type of result.\")\n"
           t_prog+="      self.failUnless(abs(res-ref)<=self.tol*abs(ref),\"wrong result\")\n"
         elif args== "constData":
           a=makeArray(sh,[-1,1])
           r=makeResult2(a,case)
           if len(sh)==0:
              t_prog+="      arg=Data(%s,self.functionspace)\n"%(a)
           else:
              t_prog+="      arg=Data(numarray.array(%s),self.functionspace)\n"%(a.tolist())
           t_prog+="      ref=%s\n"%r
           t_prog+="      res=%s(arg)\n"%case
           t_prog+="      self.failUnless(isinstance(res,float),\"wrong type of result.\")\n"
           t_prog+="      self.failUnless(abs(res-ref)<=self.tol*abs(ref),\"wrong result\")\n"
         elif args in [ "taggedData","expandedData"]:
           a=makeArray(sh,[-1,1])
           r=makeResult2(a,case)
           a1=makeArray(sh,[-1,1])
           r1=makeResult2(a1,case)
           if case in ["Lsup","sup"]:
               r=max(r,r1)
           else:
               r=min(r,r1)
           if len(sh)==0:
              if args=="expandedData":
                 t_prog+="      arg=Data(%s,self.functionspace,True)\n"%(a)
              else:
                 t_prog+="      arg=Data(%s,self.functionspace)\n"%(a)
              t_prog+="      arg.setTaggedValue(1,%s)\n"%a
           else:
              if args=="expandedData":
                 t_prog+="      arg=Data(numarray.array(%s),self.functionspace,True)\n"%(a.tolist())
              else:
                 t_prog+="      arg=Data(numarray.array(%s),self.functionspace)\n"%(a.tolist())
              t_prog+="      arg.setTaggedValue(1,%s)\n"%a1.tolist()
           t_prog+="      res=%s(arg)\n"%case
           t_prog+="      ref=%s\n"%r
           t_prog+="      self.failUnless(isinstance(res,float),\"wrong type of result.\")\n"
           t_prog+="      self.failUnless(abs(res-ref)<=self.tol*abs(ref),\"wrong result\")\n"

print t_prog

1/0


for f in func:
   u_prog+="def %s(arg):\n"%f[0]
   u_prog+="   \"\"\"\n"
   u_prog+="   returns %s of argument arg\n\n"%f[0]
   u_prog+="   @param arg: argument\n"
   u_prog+="   @type arg: C{float}, L{escript.Data}, L{escript.Symbol}, L{numarray.NumArray}.\n"
   u_prog+="   @rtype:C{float}, L{escript.Data}, L{escript.Symbol}, L{numarray.NumArray} depending on the type of arg.\n"
   u_prog+="   @raises TypeError: if the type of the argument is not expected.\n"
   u_prog+="   \"\"\"\n"
   u_prog+="   if isinstance(arg,numarray.NumArray):\n"
   u_prog+="       return numarray.%s(arg)\n"%f[0]
   u_prog+="   elif isinstance(arg,escript.Data):\n"
   u_prog+="       return arg._%s(arg)\n"%f[0]
   u_prog+="   elif isinstance(arg,float):\n"
   u_prog+="       return math.%s(arg)\n"%f[0]
   u_prog+="   elif isinstance(arg,int):\n"
   u_prog+="       return math.%s(float(arg))\n"%f[0]
   u_prog+="   elif isinstance(arg,Symbol):\n"
   u_prog+="       return Symbol_%s(arg)\n"%f[0]
   u_prog+="   else:\n"
   u_prog+="       raise TypeError,\"%s: Unknown argument type.\"\n"%f[0]
   u_prog+="class Symbol_%s(Symbol):\n"%f[0]
   u_prog+="   \"\"\"\n"
   u_prog+="   Symbol of the result of the %s function\n"%f[0]
   u_prog+="   \"\"\"\n"
   u_prog+="   def __init__(self,arg):\n"
   u_prog+="      \"\"\"\n"
   u_prog+="      initialization %s function with argument arg\n"%f[0]
   u_prog+="      @param arg: argument of function\n"
   u_prog+="      @type arg: L{Symbol}.\n"
   u_prog+="      \"\"\"\n"
   u_prog+="      Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])\n"
   u_prog+="   def __str__(self):\n"
   u_prog+="      \"\"\"\n"
   u_prog+="      string representation of the object\n"
   u_prog+="      @rtype: C{str}\n"
   u_prog+="      \"\"\"\n"
   u_prog+="      return \"%s(%%s)\"%%str(self.getArgument(0))\n"%f[0]
   u_prog+="   def substitute(self,argvals):\n"
   u_prog+="      \"\"\"\n"
   u_prog+="      assigns new values to symbols in the definition of the symbol\n      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.\n\n"
   u_prog+="      @param argvals: new values assigned to symbols\n"
   u_prog+="      @type argvals: C{dict} with keywords of type L{Symbol}.\n"
   u_prog+="      @return: result of the substitution process. Operations are executed as much as possible.\""
   u_prog+="      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution\n"
   u_prog+="      raise: TypeError: if a value for a L{Symbol} cannot be substituted.\n"
   u_prog+="      \"\"\"\n"
   u_prog+="      if argval.has_key(self): \n"
   u_prog+="         if self.isAppropriateValue(arg):\n"
   u_prog+="            return argval[self] \n"
   u_prog+="         else:\n"
   u_prog+="            raise TypeError,\"%s: new value is not appropriate.\"\n"
   u_prog+="      else:\n"
   u_prog+="         arg=self.getSubstitutedArguments(argvals)[0]\n"
   u_prog+="         return %s(arg)\n"%f[0]
   u_prog+="   def diff(self,arg):\n"
   u_prog+="      \"\"\"\n"
   u_prog+="      returns the derivative of the symbol with respect to L{Symbol} arg\n"

   u_prog+="      @param arg: the derivative is calculated with respect to arg\n"
   u_prog+="      @type arg: typically L{escript.Symbol} but can also be C{float}, L{escript.Data}, L{numarray.NumArray} depending the involved functions and data.\n"
   u_prog+="      @raise NotImplementedError: derivative not available.\n"
   u_prog+="      \"\"\"\n"
   u_prog+="      if arg==self: \n"
   u_prog+="         return identity(self.getShape())\n"
   u_prog+="      else:\n"
   u_prog+="        darg=self.getDifferentiatedArguments(arg)[0]\n"
   u_prog+="        raise NotImplementedError,\"%s: derivative is not available.\"\n"%f[0]
   u_prog+="        # dself=?\n"
   u_prog+="        return dself.matchShape(darg)*darg\n"


   for args in ["float","array","constData","taggedData","expandedData","symbol"]: 
     for sh in [ (),(2,), (4,5), (6,2,2),(3,2,3,4)]:
       if not args=="float" or len(sh)==0:
         t_prog+="   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
         tname="def test_%s_%s_rank%s"%(f[0],args,len(sh))
         t_prog+="   %s(self):\n"%tname
         if args in ["float","array" ]:
           a=makeArray(sh,f[1])
           r=makeResult(a,f)
           if len(sh)==0:
              t_prog+="      arg=%s\n"%a
              t_prog+="      ref=%s\n"%r
           else:
              t_prog+="      arg=numarray.array(%s)\n"%a.tolist()
              t_prog+="      ref=numarray.array(%s)\n"%r.tolist()
           t_prog+="      res=%s(arg)\n"%f[0]
           if args=="float":
              t_prog+="      self.failUnless(isinstance(res,float),\"wrong type of result.\")\n"
           else:
              t_prog+="      self.failUnlessEqual(res.shape,%s,\"wrong shape of result.\")\n"%str(sh)
           t_prog+="      self.failUnless(Lsup(res-ref)<=self.tol*Lsup(ref),\"wrong result\")\n"
         elif args== "constData":
           a=makeArray(sh,f[1])
           r=makeResult(a,f)
           if len(sh)==0:
              t_prog+="      arg=Data(%s,self.functionspace)\n"%(a)
              t_prog+="      ref=%s\n"%r
           else:
              t_prog+="      arg=Data(numarray.array(%s),self.functionspace)\n"%(a.tolist())
              t_prog+="      ref=numarray.array(%s)\n"%r.tolist()
           t_prog+="      res=%s(arg)\n"%f[0]
           t_prog+="      self.failUnlessEqual(res.getShape(),%s,\"wrong shape of result.\")\n"%str(sh)
           t_prog+="      self.failUnless(Lsup(res-ref)<=self.tol*Lsup(ref),\"wrong result\")\n"
         elif args in [ "taggedData","expandedData"]:
           a=makeArray(sh,f[1])
           r=makeResult(a,f)
           a1=makeArray(sh,f[1])
           r1=makeResult(a1,f)
           if len(sh)==0:
              if args=="expandedData":
                 t_prog+="      arg=Data(%s,self.functionspace,True)\n"%(a)
                 t_prog+="      ref=Data(%s,self.functionspace,True)\n"%(r)
              else:
                 t_prog+="      arg=Data(%s,self.functionspace)\n"%(a)
                 t_prog+="      ref=Data(%s,self.functionspace)\n"%(r)
              t_prog+="      arg.setTaggedValue(1,%s)\n"%a
              t_prog+="      ref.setTaggedValue(1,%s)\n"%r1
           else:
              if args=="expandedData":
                 t_prog+="      arg=Data(numarray.array(%s),self.functionspace,True)\n"%(a.tolist())
                 t_prog+="      ref=Data(numarray.array(%s),self.functionspace,True)\n"%(r.tolist())
              else:
                 t_prog+="      arg=Data(numarray.array(%s),self.functionspace)\n"%(a.tolist())
                 t_prog+="      ref=Data(numarray.array(%s),self.functionspace)\n"%(r.tolist())
              t_prog+="      arg.setTaggedValue(1,%s)\n"%a1.tolist()
              t_prog+="      ref.setTaggedValue(1,%s)\n"%r1.tolist()
           t_prog+="      res=%s(arg)\n"%f[0]
           t_prog+="      self.failUnlessEqual(res.getShape(),%s,\"wrong shape of result.\")\n"%str(sh)
           t_prog+="      self.failUnless(Lsup(res-ref)<=self.tol*Lsup(ref),\"wrong result\")\n"
         elif args=="symbol":
           t_prog+="      arg=Symbol(shape=%s)\n"%str(sh)
           t_prog+="      v=%s(arg)\n"%f[0]
           t_prog+="      self.failUnlessRaises(ValueError,v.substitute,Symbol(shape=(1,1)),\"illegal shape of substitute not identified.\")\n"
           a=makeArray(sh,f[1])
           r=makeResult(a,f)
           if len(sh)==0:
              t_prog+="      res=v.substitute({arg : %s})\n"%a
              t_prog+="      ref=%s\n"%r
              t_prog+="      self.failUnless(isinstance(res,float),\"wrong type of result.\")\n"
           else:
              t_prog+="      res=v.substitute({arg : numarray.array(%s)})\n"%a.tolist()
              t_prog+="      ref=numarray.array(%s)\n"%r.tolist()
              t_prog+="      self.failUnlessEqual(res.getShape(),%s,\"wrong shape of substitution result.\")\n"%str(sh)
           t_prog+="      self.failUnless(Lsup(res-ref)<=self.tol*Lsup(ref),\"wrong result\")\n"
 
           if len(sh)==0:
               t_prog+="      # test derivative with respect to itself:\n"
               t_prog+="      dvdv=v.diff(v)\n"
               t_prog+="      self.failUnlessEqual(dvdv,1.,\"derivative with respect to self is not 1.\")\n"
           elif len(sh)==1:
               t_prog+="      # test derivative with respect to itself:\n"
               t_prog+="      dvdv=v.diff(v)\n"
               t_prog+="      self.failUnlessEqual(dvdv.shape,%s,\"shape of derivative with respect is wrong\")\n"%str(sh+sh)
               for i0_l in range(sh[0]):
                 for i0_r in range(sh[0]):
                     if i0_l == i0_r:
                         v=1.
                     else:
                         v=0.
                     t_prog+="      self.failUnlessEqual(dvdv[%s,%s],%s,\"derivative with respect to self: [%s,%s] is not %s\")\n"%(i0_l,i0_r,v,i0_l,i0_r,v)
           elif len(sh)==2:
               t_prog+="      # test derivative with respect to itself:\n"
               t_prog+="      dvdv=v.diff(v)\n"
               t_prog+="      self.failUnlessEqual(dvdv.shape,%s,\"shape of derivative with respect is wrong\")\n"%str(sh+sh)
               for i0_l in range(sh[0]):
                 for i0_r in range(sh[0]):
                    for i1_l in range(sh[1]):
                       for i1_r in range(sh[1]):
                         if i0_l == i0_r and i1_l == i1_r:
                            v=1.
                         else:
                            v=0.
                         t_prog+="      self.failUnlessEqual(dvdv[%s,%s,%s,%s],%s,\"derivative with respect to self: [%s,%s,%s,%s] is not %s\")\n"%(i0_l,i1_l,i0_r,i1_r,v,i0_l,i1_l,i0_r,i1_r,v)
 
           for sh_in in [ (),(2,), (4,5), (6,2,2),(3,2,3,4)]:
             if len(sh_in)+len(sh)<=4:

               t_prog+="      # test derivative with shape %s as argument\n"%str(sh_in)
               trafo=makeArray(sh+sh_in,[0,1.])
               a_in=makeArray(sh_in,f[1])
               t_prog+="      arg_in=Symbol(shape=%s)\n"%str(sh_in)
               t_prog+="      arg2=Symbol(shape=%s)\n"%str(sh)

               if len(sh)==0:
                  t_prog+="      arg2="
                  if len(sh_in)==0:
                     t_prog+="%s*arg_in\n"%trafo
                  elif len(sh_in)==1:
                     for i0 in range(sh_in[0]):
                        if i0>0: t_prog+="+"
                        t_prog+="%s*arg_in[%s]"%(trafo[i0],i0)
                     t_prog+="\n"
                  elif len(sh_in)==2:
                     for i0 in range(sh_in[0]):
                        for i1 in range(sh_in[1]):
                           if i0+i1>0: t_prog+="+"
                           t_prog+="%s*arg_in[%s,%s]"%(trafo[i0,i1],i0,i1)

                  elif len(sh_in)==3:
                     for i0 in range(sh_in[0]):
                        for i1 in range(sh_in[1]):
                           for i2 in range(sh_in[2]):
                                 if i0+i1+i2>0: t_prog+="+"
                                 t_prog+="%s*arg_in[%s,%s,%s]"%(trafo[i0,i1,i2],i0,i1,i2)
                  elif len(sh_in)==4:
                     for i0 in range(sh_in[0]):
                        for i1 in range(sh_in[1]):
                           for i2 in range(sh_in[2]):
                              for i3 in range(sh_in[3]):
                                 if i0+i1+i2+i3>0: t_prog+="+"
                                 t_prog+="%s*arg_in[%s,%s,%s,%s]"%(trafo[i0,i1,i2,i3],i0,i1,i2,i3)
                     t_prog+="\n"
               elif len(sh)==1:
                  for j0 in range(sh[0]):
                     t_prog+="      arg2[%s]="%j0
                     if len(sh_in)==0:
                        t_prog+="%s*arg_in"%trafo[j0]
                     elif len(sh_in)==1:
                        for i0 in range(sh_in[0]):
                           if i0>0: t_prog+="+"
                           t_prog+="%s*arg_in[%s]"%(trafo[j0,i0],i0)
                     elif len(sh_in)==2:
                        for i0 in range(sh_in[0]):
                           for i1 in range(sh_in[1]):
                              if i0+i1>0: t_prog+="+"
                              t_prog+="%s*arg_in[%s,%s]"%(trafo[j0,i0,i1],i0,i1)
                     elif len(sh_in)==3:
                        for i0 in range(sh_in[0]):
                           for i1 in range(sh_in[1]):
                              for i2 in range(sh_in[2]):
                                    if i0+i1+i2>0: t_prog+="+"
                                    t_prog+="%s*arg_in[%s,%s,%s]"%(trafo[j0,i0,i1,i2],i0,i1,i2)
                     t_prog+="\n"
               elif len(sh)==2:
                  for j0 in range(sh[0]):
                    for j1 in range(sh[1]):
                        t_prog+="      arg2[%s,%s]="%(j0,j1)
                        if len(sh_in)==0:
                           t_prog+="%s*arg_in"%trafo[j0,j1]
                        elif len(sh_in)==1:
                           for i0 in range(sh_in[0]):
                              if i0>0: t_prog+="+"
                              t_prog+="%s*arg_in[%s]"%(trafo[j0,j1,i0],i0)
                        elif len(sh_in)==2:
                           for i0 in range(sh_in[0]):
                              for i1 in range(sh_in[1]):
                                  if i0+i1>0: t_prog+="+"
                                  t_prog+="%s*arg_in[%s,%s]"%(trafo[j0,j1,i0,i1],i0,i1)
                        t_prog+="\n"
               elif len(sh)==3:
                  for j0 in range(sh[0]):
                    for j1 in range(sh[1]):
                       for j2 in range(sh[2]):
                           t_prog+="      arg2[%s,%s,%s]="%(j0,j1,j2)
                           if len(sh_in)==0:
                              t_prog+="%s*arg_in"%trafo[j0,j1,j2]
                           elif len(sh_in)==1:
                              for i0 in range(sh_in[0]):
                                 if i0>0: t_prog+="+"
                                 t_prog+="%s*arg_in[%s]"%(trafo[j0,j1,j2,i0],i0)
                           t_prog+="\n"
               elif len(sh)==4:
                  for j0 in range(sh[0]):
                    for j1 in range(sh[1]):
                       for j2 in range(sh[2]):
                         for j3 in range(sh[3]):
                           t_prog+="      arg2[%s,%s,%s,%s]="%(j0,j1,j2,j3)
                           if len(sh_in)==0:
                              t_prog+="%s*arg_in"%trafo[j0,j1,j2,j3]
                           t_prog+="\n"
               t_prog+="      dvdin=v.substitute({arg : arg2}).diff(arg_in)\n"
               if len(sh_in)==0:
                  t_prog+="      res_in=dvdin.substitute({arg_in : %s})\n"%a_in
               else:
                  t_prog+="      res_in=dvdin.substitute({arg : numarray.array(%s)})\n"%a_in.tolist()

               if len(sh)==0:
                  if len(sh_in)==0:
                     ref_diff=(makeResult(trafo*a_in+finc,f)-makeResult(trafo*a_in,f))/finc
                     t_prog+="      self.failUnlessAlmostEqual(dvdin,%s,self.places,\"%s-derivative: wrong derivative\")\n"%(ref_diff,str(sh_in))
                  elif len(sh_in)==1:
                     s=0
                     for k0 in range(sh_in[0]):
                         s+=trafo[k0]*a_in[k0]
                     for i0 in range(sh_in[0]):
                        ref_diff=(makeResult(s+trafo[i0]*finc,f)-makeResult(s,f))/finc
                        t_prog+="      self.failUnlessAlmostEqual(dvdin[%s],%s,self.places,\"%s-derivative: wrong derivative with respect of %s\")\n"%(i0,ref_diff,str(sh_in),str(i0))
                  elif len(sh_in)==2:
                     s=0
                     for k0 in range(sh_in[0]):
                       for k1 in range(sh_in[1]):
                           s+=trafo[k0,k1]*a_in[k0,k1]
                     for i0 in range(sh_in[0]):
                        for i1 in range(sh_in[1]):
                           ref_diff=(makeResult(s+trafo[i0,i1]*finc,f)-makeResult(s,f))/finc
                           t_prog+="      self.failUnlessAlmostEqual(dvdin[%s,%s],%s,self.places,\"%s-derivative: wrong derivative with respect of %s\")\n"%(i0,i1,ref_diff,str(sh_in),str((i0,i1)))

                  elif len(sh_in)==3:
                     s=0
                     for k0 in range(sh_in[0]):
                       for k1 in range(sh_in[1]):
                          for k2 in range(sh_in[2]):
                             s+=trafo[k0,k1,k2]*a_in[k0,k1,k2]
                     for i0 in range(sh_in[0]):
                        for i1 in range(sh_in[1]):
                           for i2 in range(sh_in[2]):
                                 ref_diff=(makeResult(s+trafo[i0,i1,i2]*finc,f)-makeResult(s,f))/finc
                                 t_prog+="      self.failUnlessAlmostEqual(dvdin[%s,%s,%s],%s,self.places,\"%s-derivative: wrong derivative with respect of %s\")\n"%(i0,i1,i2,ref_diff,str(sh_in),str((i0,i1,i2)))
                  elif len(sh_in)==4:
                     s=0
                     for k0 in range(sh_in[0]):
                       for k1 in range(sh_in[1]):
                          for k2 in range(sh_in[2]):
                             for k3 in range(sh_in[3]):
                                s+=trafo[k0,k1,k2,k3]*a_in[k0,k1,k2,k3]
                     for i0 in range(sh_in[0]):
                        for i1 in range(sh_in[1]):
                           for i2 in range(sh_in[2]):
                              for i3 in range(sh_in[3]):
                                 ref_diff=(makeResult(s+trafo[i0,i1,i2,i3]*finc,f)-makeResult(s,f))/finc
                                 t_prog+="      self.failUnlessAlmostEqual(dvdin[%s,%s,%s,%s],%s,self.places,\"%s-derivative: wrong derivative with respect of %s\")\n"%(i0,i1,i2,i3,ref_diff,str(sh_in),str((i0,i1,i2,i3)))
               elif len(sh)==1:
                  for j0 in range(sh[0]):
                     if len(sh_in)==0:
                        ref_diff=(makeResult(trafo[j0]*a_in+finc,f)-makeResult(trafo[j0]*a_in,f))/finc
                        t_prog+="      self.failUnlessAlmostEqual(dvdin[%s],%s,self.places,\"%s-derivative: wrong derivative of %s\")\n"%(j0,ref_diff,str(sh_in),j0)
                     elif len(sh_in)==1:
                        s=0
                        for k0 in range(sh_in[0]):
                            s+=trafo[j0,k0]*a_in[k0]
                        for i0 in range(sh_in[0]):
                           ref_diff=(makeResult(s+trafo[j0,i0]*finc,f)-makeResult(s,f))/finc
                           t_prog+="      self.failUnlessAlmostEqual(dvdin[%s,%s],%s,self.places,\"%s-derivative: wrong derivative of %s with respect of %s\")\n"%(j0,i0,ref_diff,str(sh_in),str(j0),str(i0))
                     elif len(sh_in)==2:
                        s=0
                        for k0 in range(sh_in[0]):
                           for k1 in range(sh_in[1]):
                               s+=trafo[j0,k0,k1]*a_in[k0,k1]
                        for i0 in range(sh_in[0]):
                           for i1 in range(sh_in[1]):
                              ref_diff=(makeResult(s+trafo[j0,i0,i1]*finc,f)-makeResult(s,f))/finc
                              t_prog+="      self.failUnlessAlmostEqual(dvdin[%s,%s,%s],%s,self.places,\"%s-derivative: wrong derivative of %s with respect of %s\")\n"%(j0,i0,i1,ref_diff,str(sh_in),str(j0),str((i0,i1)))

                     elif len(sh_in)==3:

                        s=0
                        for k0 in range(sh_in[0]):
                          for k1 in range(sh_in[1]):
                             for k2 in range(sh_in[2]):
                                s+=trafo[j0,k0,k1,k2]*a_in[k0,k1,k2]

                        for i0 in range(sh_in[0]):
                           for i1 in range(sh_in[1]):
                              for i2 in range(sh_in[2]):
                                    ref_diff=(makeResult(s+trafo[j0,i0,i1,i2]*finc,f)-makeResult(s,f))/finc
                                    t_prog+="      self.failUnlessAlmostEqual(dvdin[%s,%s,%s,%s],%s,self.places,\"%s-derivative: wrong derivative of %s with respect of %s\")\n"%(j0,i0,i1,i2,ref_diff,str(sh_in),str(j0),str((i0,i1,i2)))
               elif len(sh)==2:
                  for j0 in range(sh[0]):
                    for j1 in range(sh[1]):
                        if len(sh_in)==0:
                           ref_diff=(makeResult(trafo[j0,j1]*a_in+finc,f)-makeResult(trafo[j0,j1]*a_in,f))/finc
                           t_prog+="      self.failUnlessAlmostEqual(dvdin[%s,%s],%s,self.places,\"%s-derivative: wrong derivative of %s\")\n"%(j0,j1,ref_diff,str(sh_in),str((j0,j1)))
                        elif len(sh_in)==1:
                           s=0
                           for k0 in range(sh_in[0]):
                               s+=trafo[j0,j1,k0]*a_in[k0]
                           for i0 in range(sh_in[0]):
                              ref_diff=(makeResult(s+trafo[j0,j1,i0]*finc,f)-makeResult(s,f))/finc
                              t_prog+="      self.failUnlessAlmostEqual(dvdin[%s,%s,%s],%s,self.places,\"%s-derivative: wrong derivative of %s with respect of %s\")\n"%(j0,j1,i0,ref_diff,str(sh_in),str((j0,j1)),str(i0))

                        elif len(sh_in)==2:
                           s=0
                           for k0 in range(sh_in[0]):
                              for k1 in range(sh_in[1]):
                                  s+=trafo[j0,j1,k0,k1]*a_in[k0,k1]
                           for i0 in range(sh_in[0]):
                              for i1 in range(sh_in[1]):
                                  ref_diff=(makeResult(s+trafo[j0,j1,i0,i1]*finc,f)-makeResult(s,f))/finc
                                  t_prog+="      self.failUnlessAlmostEqual(dvdin[%s,%s,%s,%s],%s,self.places,\"%s-derivative: wrong derivative of %s with respect of %s\")\n"%(j0,j1,i0,i1,ref_diff,str(sh_in),str((j0,j1)),str((i0,i1)))
               elif len(sh)==3:
                  for j0 in range(sh[0]):
                    for j1 in range(sh[1]):
                       for j2 in range(sh[2]):
                           if len(sh_in)==0:
                              ref_diff=(makeResult(trafo[j0,j1,j2]*a_in+finc,f)-makeResult(trafo[j0,j1,j2]*a_in,f))/finc
                              t_prog+="      self.failUnlessAlmostEqual(dvdin[%s,%s,%s],%s,self.places,\"%s-derivative: wrong derivative of %s\")\n"%(j0,j1,j2,ref_diff,str(sh_in),str((j0,j1,j2)))
                           elif len(sh_in)==1:
                              s=0
                              for k0 in range(sh_in[0]):
                                  s+=trafo[j0,j1,j2,k0]*a_in[k0]
                              for i0 in range(sh_in[0]):
                                 ref_diff=(makeResult(s+trafo[j0,j1,j2,i0]*finc,f)-makeResult(s,f))/finc
                                 t_prog+="      self.failUnlessAlmostEqual(dvdin[%s,%s,%s,%s],%s,self.places,\"%s-derivative: wrong derivative of %s with respect of %s\")\n"%(j0,j1,j2,i0,ref_diff,str(sh_in),str((j0,j1,j2)),i0)
               elif len(sh)==4:
                  for j0 in range(sh[0]):
                    for j1 in range(sh[1]):
                       for j2 in range(sh[2]):
                         for j3 in range(sh[3]):
                           if len(sh_in)==0:
                              ref_diff=(makeResult(trafo[j0,j1,j2,j3]*a_in+finc,f)-makeResult(trafo[j0,j1,j2,j3]*a_in,f))/finc
                              t_prog+="      self.failUnlessAlmostEqual(dvdin[%s,%s,%s,%s],%s,self.places,\"%s-derivative: wrong derivative of %s\")\n"%(j0,j1,j2,j3,ref_diff,str(sh_in),str((j0,j1,j2,j3)))
               
# print u_prog
print t_prog
