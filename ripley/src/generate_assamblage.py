
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################



"""
this script generates the assemblage routine for the ripley rectangular grid solver


"""
from sympy import *

DIGITS=20

#
#  quadrature form in [0,1]
#

q0=1/sqrt(RealNumber(3))
qD1=[ (1-q0)/2,  (1+q0)/2 ]
wD1=[RealNumber(1)/2, RealNumber(1)/2 ]


q0=1/sqrt(RealNumber(3))
qD1=[ (1-q0)/2,  (1+q0)/2 ]
wD1=[RealNumber(1)/2, RealNumber(1)/2 ]

qD1_r=[ RealNumber(1)/2 ]
wD1_r= [RealNumber(1)]

print "1D quadrature nodes =",qD1
print "1D quadrature weights =",wD1
print "1D reduced quadrature nodes =",qD1_r
print "1D reduced quadrature weights =",wD1_r

def generate(DIM):
   # variables
   h=[Symbol('h%s'%i) for i in range(DIM) ]
   x=[Symbol('x%s'%i) for i in range(DIM) ]
   #  shape functions: labeling differently from finley!!!!
   S=[ 1 ]
   GridOffset=[[]]
   GridOffsetText=[ "" ]
   idx_full=[]
   Q = [ [] ]
   W = [ 1 ]
   Q_r = [ [] ]
   W_r = [ 1 ]
   for i in range(DIM):
       idx_full.append(i)
       S = [ s * (1-x[i]/h[i])  for s in S] + [ s * x[i]/h[i]  for s in S]
       GridOffset = [ s + [0]  for s in GridOffset] + [ s + [1]  for s in GridOffset]
       GridOffsetText = [ s + "0"  for s in GridOffsetText] + [ s + "1"  for s in GridOffsetText]

       # generate quadrature points in element
       Qnew, Wnew =[], []
       for j in range(len(qD1)): 
              Qnew += [ q + [qD1[j]*h[i],] for q in Q ]
              Wnew += [ w * wD1[j]*h[i] for w in W ]
       Q, W=Qnew, Wnew

       Qnew, Wnew =[], []
       for j in range(len(qD1_r)): 
               Qnew += [ q + [qD1_r[j]*h[i],] for q in Q_r ]
               Wnew += [ w * wD1_r[j]*h[i] for w in W_r ]
       Q_r, W_r=Qnew, Wnew

   # now the same thing for the faces:
   Q_faces=[]
   W_faces=[]
   Q_r_faces=[]
   W_r_faces=[]
   idx_faces=[]
   for face in range(DIM):
       Q_left, W_left = [ [] ], [ 1 ]
       Q_left_r, W_left_r = [ [] ], [ 1 ]
       Q_right, W_right = [ [] ], [ 1 ]
       Q_right_r, W_right_r = [ [] ], [ 1 ]
       idx_left,idx_right=[],[]
       for i in range(DIM):  
          # generate quadrature points in element
          Q_left_new, W_left_new =[], []
          Q_right_new, W_right_new =[], []
          
          if face == i :
             idx_left.append(-1)
             idx_right.append(-2)
             Q_left_new += [ q + [0.,] for q in Q_left ]
             W_left_new += [ w * 1. for w in W_left ]
             Q_right_new += [ q + [ h[i], ] for q in Q_right ]
             W_right_new += [ w * 1. for w in W_right ]
          else:
             idx_left.append(i)
             idx_right.append(i)
             for j in range(len(qD1)): 
                 Q_left_new += [ q + [qD1[j]*h[i],] for q in Q_left ]
                 W_left_new += [ w * wD1[j]*h[i] for w in W_left ]
                 Q_right_new += [ q + [qD1[j]*h[i],] for q in Q_right ]
                 W_right_new += [ w * wD1[j]*h[i] for w in W_right ]
          Q_left, W_left=Q_left_new, W_left_new
          Q_right, W_right=Q_right_new, W_right_new
       Q_faces=Q_faces+[Q_left, Q_right]
       W_faces=W_faces+[W_left, W_right]
       idx_faces=idx_faces + [ idx_left, idx_right ]

       Q_left, W_left = [ [] ], [ 1 ]
       Q_left_r, W_left_r = [ [] ], [ 1 ]
       Q_right, W_right = [ [] ], [ 1 ]
       Q_right_r, W_right_r = [ [] ], [ 1 ]
       for i in range(DIM):
          # generate quadrature points in element
          Q_left_new, W_left_new =[], []
          Q_right_new, W_right_new =[], []
          if face == i :
             Q_left_new += [ q + [0.,] for q in Q_left ]
             W_left_new += [ w * 1. for w in W_left ]
             Q_right_new += [ q + [ h[i], ] for q in Q_right ]
             W_right_new += [ w * 1. for w in W_right ]
          else:
             for j in range(len(qD1_r)): 
                 Q_left_new += [ q + [qD1_r[j]*h[i],] for q in Q_left ]
                 W_left_new += [ w * wD1_r[j]*h[i] for w in W_left ]
                 Q_right_new += [ q + [qD1_r[j]*h[i],] for q in Q_right ]
                 W_right_new += [ w * wD1_r[j]*h[i] for w in W_right ]
          Q_left, W_left=Q_left_new, W_left_new
          Q_right, W_right=Q_right_new, W_right_new
       Q_r_faces=Q_r_faces+[Q_left, Q_right]
       W_r_faces=W_r_faces+[W_left, W_right]
       
   # the interpolation function
   f_interpolation = 0
   for i in range(len(GridOffsetText)):
     f_interpolation = f_interpolation + S[i] * Symbol("f_%s"%GridOffsetText[i])

   # interpolation to Quad points
   CODE="\nif (out_data_type==RIPLEY_ELEMENTS) {\n"
   CODE+=createCode(f_interpolation, x, Q, loopindex=idx_full)
   CODE+="} else if (out_data_type==RIPLEY_REDUCED_ELEMENTS) {\n"
   CODE+=createCode(f_interpolation, x, Q_r, loopindex=idx_full)
   CODE+="} else if (out_data_type==RIPLEY_BOUNDARY_ELEMENTS) {\n"
   for i in range(len(Q_r_faces)):
        CODE+="if (face_offset(%s)>-1) {\n"%i
        CODE+=createCode(f_interpolation, x, Q_faces[i], loopindex=idx_faces[i],  gridoffset="face_offset(%s)"%i)
        CODE+="\n} /* end of face %s */\n"%i
   CODE+="} else if (out_data_type==RIPLEY_REDUCED_BOUNDARY_ELEMENTS) {\n"
   for i in range(len(Q_r_faces)):
        CODE+="if (face_offset(%s)>-1) {\n"%i
        CODE+=createCode(f_interpolation, x, Q_r_faces[i], loopindex=idx_faces[i],  gridoffset="face_offset(%s)"%i)
        CODE+="\n} /* end of face %s */\n"%i
   CODE+="\n} /* end of out_data_type branching */\n"
   insertCode("Assemble_Interpolation_%sD.c"%DIM, { "SNIP" : CODE})
   
      # gradient to Quad points
   CODE="\nif (out_data_type==RIPLEY_ELEMENTS) {\n"
   CODE+=createGradientCode(f_interpolation, x, Q, loopindex=idx_full)
   CODE+="} else if (out_data_type==RIPLEY_REDUCED_ELEMENTS) {\n"
   CODE+=createGradientCode(f_interpolation, x, Q_r, loopindex=idx_full)
   CODE+="} else if (out_data_type==RIPLEY_BOUNDARY_ELEMENTS) {\n"
   for i in range(len(Q_r_faces)):
        CODE+="if (face_offset(%s)>-1) {\n"%i
        CODE+=createGradientCode(f_interpolation, x, Q_faces[i], loopindex=idx_faces[i],  gridoffset="face_offset(%s)"%i)
        CODE+="\n} /* end of face %s */\n"%i
   CODE+="} else if (out_data_type==RIPLEY_REDUCED_BOUNDARY_ELEMENTS) {\n"
   for i in range(len(Q_r_faces)):
        CODE+="if (face_offset(%s)>-1) {\n"%i
        CODE+=createGradientCode(f_interpolation, x, Q_r_faces[i], loopindex=idx_faces[i],  gridoffset="face_offset(%s)"%i)
        CODE+="\n} /* end of face %s */\n"%i
   CODE+="\n} /* end of out_data_type branching */\n"
   insertCode("Assemble_Gradient_%sD.c"%DIM, { "SNIP" : CODE})

   #generate PDE assemblage
   CODE,PRECODE = makePDE(S, x, Q, W, DIM=DIM, system=True)
   insertCode("Assemble_PDE_System_%sD.c"%DIM, { "SNIP" : CODE , "SNIP_PRE" : PRECODE } )
   CODE,PRECODE = makePDE(S, x, Q_r, W_r, DIM=DIM, system=True)
   insertCode("Assemble_PDE_System_%sD_reduced.c"%DIM, { "SNIP" : CODE , "SNIP_PRE" : PRECODE })
   CODE,PRECODE = makePDE(S, x, Q, W, DIM=DIM, system=False)
   insertCode("Assemble_PDE_Single_%sD.c"%DIM, { "SNIP" : CODE , "SNIP_PRE" : PRECODE })
   CODE,PRECODE = makePDE(S, x, Q_r, W_r, DIM=DIM, system=False)
   insertCode("Assemble_PDE_Single_%sD_reduced.c"%DIM, { "SNIP" : CODE , "SNIP_PRE" : PRECODE })


    

def insertCode(fn, replacement):
    TOP_LINE_FMT="/* GENERATOR %s TOP"
    BOTTOM_LINE_FMT="/* GENERATOR %s BOTTOM"
    TAG="   "

    # read code and create back up 
    bkp_stream=open(fn+".bkp",'w')
    old_code=""
    for l in open(fn, 'r').readlines():
       s=l.strip()
       if len(s)>0:
          bkp_stream.write(s+"\n")
          old_code+=s+"\n"
    bkp_stream.close()
    
    for k in replacement.keys():
        TOP_LINE=TOP_LINE_FMT%k
        BOTTOM_LINE=BOTTOM_LINE_FMT%k
        new_code=""
        ignore=False
        for l in old_code.split("\n"):           
           s=l.strip()
           if len(s)>0:
              if s.startswith(TOP_LINE): 
                 ignore=True
                 new_code+=s+"\n"+replacement[k]
              elif s.startswith(BOTTOM_LINE): 
                 ignore=False
                 new_code+=s+"\n"
              elif not ignore:
                 new_code+=s+"\n"
        old_code=new_code
    #insert tabs:
    N=0
    new_code=""
    for l in old_code.split("\n"):
       s=l.strip()
       if s.find("/*")>-1 and s.find("*/")>-1 :  
            s2=(s[:s.find("/*")] + s[s.find("*/")+2:]).strip()
       else:
            s2=s
       if len(s)>0:
          if s2.startswith("}"): 
             N=max(N-1, 0)
          new_code+=TAG*N + s + "\n"
          if s2.endswith("{"): 
             N+=1
    open(fn,'w').write(new_code)
    print "file %s generated."%fn

   

def optimizeEvaluate(F, x, Q):
  F2=[]
  for f in F:
     f0=[]
     for q in range(len(Q)):
        tmp=f
        for i in range(len(x)): tmp=tmp.subs(x[i], Q[q][i])
        f0.append(tmp)
     F2.append(f0)
 
  # collect arguments:
  args=[]
  for f in F2:
   for fq in f:
     for s in fq.atoms():
       if isinstance(s, Symbol):
         if collect(fq, s, evaluate=False).has_key(s):
            if s.name.startswith("f_"):
              if args.count(s) == 0 and not (collect(fq, s, evaluate=False)[s]).is_zero: args.append(s)
  c0={}
  cc=0
  F_new=[]
  for f in F2:
   fq_new=[]
   for fq in f:
       s1={}
       for a in args:
         s=collect(fq, a, evaluate=False)
         if s.has_key(a):
             k0=None
             for k, v in c0.items():
                if s[a] == v: k0=k
             if k0==None:
                 k0=Symbol("tmp0_%s"%cc)
                 c0[k0]=s[a]
                 cc+=1
             if not s1.has_key(k0): s1[k0]=0
             s1[k0] = s1[k0] + a
       tmp=0
       for k,v in s1.items():
          tmp = tmp + k * v
       fq_new.append(tmp)
   F_new.append(fq_new)
      
  return F_new, args, c0

def createGradientCode(F, x, Q, gridoffset="", loopindex=(0,1,2)):
   DIM=len(loopindex)
   dF=[ expand(diff(F, x[i])) for i in range(DIM)]
   dF, args, consts = optimizeEvaluate(dF, x, Q)
   k=""
   N=""
   M1=""
   for i in range(len(Q[0])):
     if len(k)==0:
         k+="k%s"%i
         N+="N%s"%i
         M1+="M%s"%i
     else:
         k+=",k%s"%i
         if i < DIM-1: N+=",N%s"%i
         if i < DIM-1: M1+=",M%s"%i
   
   k3=""
   for i in range(DIM-1,-1,-1): 
          if loopindex[i] > -1: k3+=",k%i"%loopindex[i]
   TXT=""
   for a,v in  consts.items():
      TXT+="const double %s = %s;\n"%(ccode(a), ccode(v.evalf(n=DIGITS)))
   TXT+="#pragma omp parallel for private(i%s)"%k3

   for i in range(DIM-1,-1,-1): 
          if loopindex[i] > -1: TXT+="\nfor (k%i =0; k%i < N%i; ++k%i) {"%(loopindex[i],loopindex[i],loopindex[i],loopindex[i])
   TXT+="\nfor (i =0; i < NCOMP; ++i) {\n"
   for s in args:
            k1=""
            for i in range(DIM):
              if s.name[2+i]=="0":
                  if loopindex[i]>-1:
                     k1+=",k%s"%i
                  elif  loopindex[i] == -1:
                     k1+=",0"
                  elif  loopindex[i] == -2:
                     k1+=",M%s-2"%i
              else:
                  if loopindex[i]>-1:
                     k1+=",k%s+1"%i
                  elif  loopindex[i] == -1:
                     k1+=",1"
                  elif  loopindex[i] == -2:
                     k1+=",M%s-1"%i
            TXT+="const double %s = in[INDEX2(i,INDEX%s(%s, %s),NCOMP)];\n"%(s.name,DIM,k1[1:],M1)
   #  interpolation to quadrature points
   for q in range(len(Q)):
         IDX2="INDEX%s(%s,%s)"%(DIM,k,N)
         if len(gridoffset) > 0: IDX2=gridoffset+"+"+IDX2
         for i in range(DIM):
             TXT+="out[INDEX4(i,%s,%s,%s,NCOMP,%s,%s)] = %s;\n"%(i,q,IDX2,DIM,len(Q),ccode(dF[i][q]))
   TXT+="} /* close component loop i */\n"
   for i in range(DIM): 
           if loopindex[i]>-1 : TXT+="} /* close k%i loop */\n"%loopindex[i]
   return TXT

def createCode(F, x, Q, gridoffset="", loopindex=(0,1,2)):
   DIM=len(loopindex)
   F, args, consts = optimizeEvaluate([F, ], x, Q)
   F=F[0]
   k=""
   N=""
   M1=""
   for i in range(len(Q[0])):
     if len(k)==0:
         k+="k%s"%i
         N+="N%s"%i
         M1+="M%s"%i
     else:
         k+=",k%s"%i
         if i < DIM-1: N+=",N%s"%i
         if i < DIM-1: M1+=",M%s"%i
   
   k3=""
   for i in range(DIM-1,-1,-1): 
          if loopindex[i] > -1: k3+=",k%i"%loopindex[i]
   TXT=""
   for a,v in  consts.items():
      TXT+="const double %s = %s;\n"%(ccode(a), ccode(v.evalf(n=DIGITS)))
   TXT+="#pragma omp parallel for private(i%s)"%k3
   for i in range(DIM-1,-1,-1): 
          if loopindex[i] > -1: TXT+="\nfor (k%i =0; k%i < N%i; ++k%i) {"%(loopindex[i],loopindex[i],loopindex[i],loopindex[i])
   TXT+="\nfor (i =0; i < NCOMP; ++i) {\n"
   for s in args:
            k1=""
            for i in range(DIM):
              if s.name[2+i]=="0":
                  if loopindex[i]>-1:
                     k1+=",k%s"%i
                  elif  loopindex[i] == -1:
                     k1+=",0"
                  elif  loopindex[i] == -2:
                     k1+=",M%s-2"%i
              else:
                  if loopindex[i]>-1:
                     k1+=",k%s+1"%i
                  elif  loopindex[i] == -1:
                     k1+=",1"
                  elif  loopindex[i] == -2:
                     k1+=",M%s-1"%i
            TXT+="const double %s = in[INDEX2(i,INDEX%s(%s, %s),NCOMP)];\n"%(s.name,DIM,k1[1:],M1)
   #  interpolation to quadrature points
   for q in range(len(Q)):
      IDX2="INDEX%s(%s,%s)"%(DIM,k,N)
      if len(gridoffset) > 0: IDX2=gridoffset+"+"+IDX2
      TXT+="out[INDEX3(i,%s,%s,NCOMP,%s)] = %s;\n"%(q,IDX2,len(Q),ccode(F[q]))
   TXT+="} /* close component loop i */\n"
   for i in range(DIM): 
           if loopindex[i]>-1 : TXT+="} /* close k%i loop */\n"%loopindex[i]
   return TXT
   
def subPoint(g,x, p):
     out=g
     for i in range(len(x)): out=out.subs(x[i], p[i])
     # return out.evalf()
     return out

def generatePDECode(DATA_A, EM,GLOBAL_TMP, system=False):         
        LOCAL_TMP={}
        LOCAL2_TMP={}
        OUT=""
        for p in EM.keys():
            # first we collect the constant coefficients in the expression (outside the parallel region):
            tmp={}
            for a in DATA_A:
                c=  collect(EM[p], a, evaluate=False)
                if c.has_key(a):
                   if not GLOBAL_TMP.has_key(c[a]):
                           GLOBAL_TMP[c[a]]=Symbol("w%s"%(len(GLOBAL_TMP),))
                   if tmp.has_key(GLOBAL_TMP[c[a]]):
                       tmp[GLOBAL_TMP[c[a]]]+=a
                   else:
                        tmp[GLOBAL_TMP[c[a]]]=a
            # now we find temporary values which are combinations of the input data:
            em_new=[]
            for t in tmp.keys():
                if tmp[t] in DATA_A:
                   tt = t * tmp[t]
                else:
                   if not LOCAL_TMP.has_key(tmp[t]):
                      LOCAL_TMP[tmp[t]]=Symbol("tmp%s_0"%len(LOCAL_TMP))
                   tt= t * LOCAL_TMP[tmp[t]]
                if not LOCAL2_TMP.has_key(tt):
                      LOCAL2_TMP[tt]=Symbol("tmp%s_1"%len(LOCAL2_TMP))
                em_new.append(LOCAL2_TMP[tt])
            EM[p]=em_new
                         
        for p in EM.keys():
            sum_em=0
            for t in EM[p]: sum_em=sum_em + t
            if isinstance(p, int) :
              if system:
                 OUT+="  EM_F[INDEX2(k,%s,p.numEqu)]+=%s;\n"%(p,ccode(sum_em))
              else:
                 OUT+="  EM_F[%s]+=%s;\n"%(p,ccode(sum_em))
            else:
              if system:
                 OUT+="  EM_S[INDEX4(k,m,%s,%s,p.numEqu,p.numComp,%s)]+=%s;\n"%(p[0],p[1],max([ qq[0] for qq in EM.keys()])+1,ccode(sum_em))
              else:
                 OUT+="  EM_S[INDEX2(%s,%s,%s)]+=%s;\n"%(p[0],p[1],max([ qq[0] for qq in EM.keys()])+1,ccode(sum_em)) 

        OUT2=""
        for p in LOCAL_TMP:
             OUT2+="  const double %s = %s;\n"%(LOCAL_TMP[p],ccode(p))
        for p in LOCAL2_TMP:
             OUT2+="  const double %s = %s;\n"%(LOCAL2_TMP[p],ccode(p))
        return OUT2+OUT

def makePDE(S, x, Q, W, DIM=2, system=False):
   GLOBAL_TMP={}
   GLOBAL_N=0
   PRECODE=""
   CODE="""
        /**************************************************************/
        /*   process A: */
        /**************************************************************/   
        if (NULL!=A_p) {
            add_EM_S=TRUE; 
            """
   if len(Q) > 1:
        CODE+="if (extendedA) {\n"
        if system: CODE+= """for (k=0;k<p.numEqu;k++) {
                            for (m=0;m<p.numComp;m++) {  
                         """
        DATA_A=[]
        CODE2=""
        EM = {}
        for i in range(len(S)):
            for j in range(len(S)):
                EM[(i,j)] = 0
                
        for q in range(len(Q)):
            for di in range(DIM):
               for dj in range(DIM):
                    A_name="A_%d%d_%d"%(di,dj,q)
                  A=Symbol(A_name)
                  DATA_A.append(A)
                  if system:
                      CODE2+="const double %s = A_p[INDEX5(k,%s,m,%s,%s, p.numEqu,%s,p.numComp,%s)];\n"%(A_name, di, dj, q, DIM, DIM)
                  else:
                      CODE2+="const double %s = A_p[INDEX3(%s,%s,%s,%s,%s)];\n"%(A_name, di, dj, q, DIM, DIM)
                   for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[i],x[di]) * diff(S[j],x[dj])).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_A, EM, GLOBAL_TMP,system)
   
        if system: CODE+="}\n }\n" 
        CODE+="} else { /* constant data */\n"
   if system: 
       CODE+= """for (k=0;k<p.numEqu;k++) {
                            for (m=0;m<p.numComp;m++) {
            """ 
   DATA_A=[]
   CODE2=""
   EM = {}
   for i in range(len(S)):
       for j in range(len(S)):
          EM[(i,j)] = 0 
   for di in range(DIM):
       for dj in range(DIM):
               A_name="A_%d%d"%(di,dj)
             A=Symbol(A_name)
             DATA_A.append(A)
             if system:
                      CODE2+="const double %s = A_p[INDEX4(k,%s,m,%s p.numEqu,%s, p.numComp)];\n"%(A_name, di, dj, DIM)
             else:
                      CODE2+="const double %s = A_p[INDEX2(%s,%s,%s)];\n"%(A_name, di, dj, DIM)
             for q in range(len(Q)):
                   for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[i],x[di]) * diff(S[j],x[dj])).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_A, EM, GLOBAL_TMP,system)    
   if system:  CODE+="}\n }\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n" 
   #BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
   CODE+="""
        /**************************************************************/
        /*   process B: */
        /**************************************************************/   
        if (NULL!=B_p) {
            add_EM_S=TRUE; 
            """
   if len(Q) > 1:
        CODE+="if (extendedB) {\n"
        if system: CODE+= """for (k=0;k<p.numEqu;k++) {
                            for (m=0;m<p.numComp;m++) {  
                         """
        DATA_B=[]
        CODE2=""
        EM = {}
        for i in range(len(S)):
            for j in range(len(S)):
                EM[(i,j)] = 0
                
        for q in range(len(Q)):
            for di in range(DIM):
                    A_name="B_%d_%d"%(di,q)
                  A=Symbol(A_name)
                  DATA_B.append(A)
                  if system:
                      CODE2+="const double %s = B_p[INDEX4(k,%s,m,%s, p.numEqu,%s,p.numComp)];\n"%(A_name, di,  q, DIM)
                  else:
                      CODE2+="const double %s = B_p[INDEX2(%s,%s,%s)];\n"%(A_name, di, q, DIM)
                   for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[i],x[di]) * S[j]).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_B, EM, GLOBAL_TMP,system)
        
        if system: CODE+="}\n }\n" 
        CODE+="} else { /* constant data */\n"
   if system: 
       CODE+= """for (k=0;k<p.numEqu;k++) {
                            for (m=0;m<p.numComp;m++) {
            """ 
   DATA_B=[]
   CODE2=""
   EM = {}
   for i in range(len(S)):
       for j in range(len(S)):
          EM[(i,j)] = 0 
   for di in range(DIM):
               A_name="B_%d"%(di)
             A=Symbol(A_name)
             DATA_B.append(A)
             if system:
                      CODE2+="const double %s = B_p[INDEX3(k,%s,m, p.numEqu,%s)];\n"%(A_name, di,  DIM)
             else:
                      CODE2+="const double %s = B_p[%s];\n"%(A_name, di)
             for q in range(len(Q)):
                   for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[i],x[di]) * S[j] ).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_B, EM, GLOBAL_TMP,system)    
   if system:  CODE+="}\n }\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n"  
   #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   CODE+="""
        /**************************************************************/
        /*   process C: */
        /**************************************************************/   
        if (NULL!=C_p) {
            add_EM_S=TRUE; 
            """
   if len(Q) > 1:
        CODE+="if (extendedC) {\n"
        if system: CODE+= """for (k=0;k<p.numEqu;k++) {
                            for (m=0;m<p.numComp;m++) {  
                         """
        DATA_C=[]
        CODE2=""
        EM = {}
        for i in range(len(S)):
            for j in range(len(S)):
                EM[(i,j)] = 0
                
        for q in range(len(Q)):
            for dj in range(DIM):
                    A_name="C_%d_%d"%(dj,q)
                  A=Symbol(A_name)
                  DATA_C.append(A)
                  if system:
                      CODE2+="const double %s = C_p[INDEX4(k,m,%s, %s, p.numEqu,p.numComp,%s)];\n"%(A_name, dj,  q, DIM)
                  else:
                      CODE2+="const double %s = C_p[INDEX2(%s,%s,%s)];\n"%(A_name, dj, q, DIM)
                   for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[j],x[dj]) * S[i]).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_C, EM, GLOBAL_TMP,system)
        
        if system: CODE+="}\n }\n" 
        CODE+="} else { /* constant data */\n"
   if system: 
       CODE+= """for (k=0;k<p.numEqu;k++) {
                            for (m=0;m<p.numComp;m++) {
            """ 
   DATA_C=[]
   CODE2=""
   EM = {}
   for i in range(len(S)):
       for j in range(len(S)):
          EM[(i,j)] = 0 
   for dj in range(DIM):
               A_name="C_%d"%(dj)
             A=Symbol(A_name)
             DATA_C.append(A)
             if system:
                      CODE2+="const double %s = C_p[INDEX3(k,m,%s, p.numEqu,p.numComp)];\n"%(A_name, dj)
             else:
                      CODE2+="const double %s = C_p[%s];\n"%(A_name, dj)
             for q in range(len(Q)):
                   for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[j],x[dj]) * S[i] ).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_C, EM, GLOBAL_TMP,system)    
   if system:  CODE+="}\n }\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n"  
   #DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
   CODE+="""
        /**************************************************************/
        /*   process D: */
        /**************************************************************/   
        if (NULL!=D_p) {
            add_EM_S=TRUE; 
            """
   if len(Q) > 1:
        CODE+="if (extendedD) {\n"
        if system: CODE+= """for (k=0;k<p.numEqu;k++) {
                            for (m=0;m<p.numComp;m++) {  
                         """
        DATA_D=[]
        CODE2=""
        EM = {}
        for i in range(len(S)):
            for j in range(len(S)):
                EM[(i,j)] = 0
                
        for q in range(len(Q)):
                    A_name="D_%d"%(q,)
                  A=Symbol(A_name)
                  DATA_D.append(A)
                  if system:
                      CODE2+="const double %s = D_p[INDEX3(k,m,%s, p.numEqu,p.numComp)];\n"%(A_name, q)
                  else:
                      CODE2+="const double %s = D_p[%s];\n"%(A_name, q)
                   for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * S[j] * S[i]).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_D, EM, GLOBAL_TMP,system)
        
        if system: CODE+="}\n }\n" 
        CODE+="} else { /* constant data */\n"
   if system: 
       CODE+= """for (k=0;k<p.numEqu;k++) {
                            for (m=0;m<p.numComp;m++) {
            """ 
   DATA_D=[]
   CODE2=""
   EM = {}
   for i in range(len(S)):
       for j in range(len(S)):
          EM[(i,j)] = 0 
   if 1:
               A_name="D_0"
             A=Symbol(A_name)
             DATA_D.append(A)
             if system:
                      CODE2+="const double %s = D_p[INDEX2(k,m, p.numEqu)];\n"%(A_name,)
             else:
                      CODE2+="const double %s = D_p[0];\n"%(A_name,)
             for q in range(len(Q)):
                   for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * S[j] * S[i] ).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_D, EM, GLOBAL_TMP,system)    
   if system:  CODE+="}\n }\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n"  

   #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   CODE+="""
        /**************************************************************/
        /*   process X: */
        /**************************************************************/   
        if (NULL!=X_p) {
            add_EM_F=TRUE; 
            """
   if len(Q) > 1:
        CODE+="if (extendedX) {\n"
        if system: CODE+= """for (k=0;k<p.numEqu;k++)   {
                         """
        DATA_X=[]
        CODE2=""
        EM = {}
        for i in range(len(S)):
                EM[i] = 0
                
        for q in range(len(Q)):
            for dj in range(DIM):
                    A_name="X_%d_%d"%(dj,q)
                  A=Symbol(A_name)
                  DATA_X.append(A)
                  if system:
                      CODE2+="const double %s = X_p[INDEX3(k,%s, %s, p.numEqu,%s)];\n"%(A_name, dj,  q, DIM)
                  else:
                      CODE2+="const double %s = X_p[INDEX2(%s,%s,%s)];\n"%(A_name, dj,q,DIM)
                   for j in range(len(S)):
                          EM[j] = EM[j] + (A * W[q] * diff(S[j],x[dj])).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_X, EM, GLOBAL_TMP,system)
        
        if system: CODE+="}\n" 
        CODE+="} else { /* constant data */\n"
   if system: 
       CODE+= """for (k=0;k<p.numEqu;k++) {
            """ 
   DATA_X=[]
   CODE2=""
   EM = {}
   for i in range(len(S)):
          EM[i] = 0 
   for dj in range(DIM):
               A_name="X_%d"%(dj)
             A=Symbol(A_name)
             DATA_X.append(A)
             if system:
                      CODE2+="const double %s = X_p[INDEX2(k,%s, p.numEqu)];\n"%(A_name, dj)
             else:
                      CODE2+="const double %s = X_p[%s];\n"%(A_name, dj)
             for q in range(len(Q)):
                   for j in range(len(S)):
                          EM[j] = EM[j] + (A * W[q] * diff(S[j],x[dj])).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_X, EM, GLOBAL_TMP,system)    
   if system:  CODE+="}\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n"  

   #YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYy
   CODE+="""
        /**************************************************************/
        /*   process Y: */
        /**************************************************************/   
        if (NULL!=Y_p) {
            add_EM_F=TRUE; 
            """
   if len(Q) > 1:
        CODE+="if (extendedY) {\n"
        if system: CODE+= """for (k=0;k<p.numEqu;k++) {
                         """
        DATA_Y=[]
        CODE2=""
        EM = {}
        for i in range(len(S)):
                EM[i] = 0
                
        for q in range(len(Q)):
                    A_name="Y_%d"%(q,)
                  A=Symbol(A_name)
                  DATA_Y.append(A)
                  if system:
                      CODE2+="const double %s = Y_p[INDEX3(k,%s, p.numEqu)];\n"%(A_name, q)
                  else:
                      CODE2+="const double %s = Y_p[%s];\n"%(A_name, q)
                   for i in range(len(S)):
                          EM[i] = EM[i] + (A * W[q] * S[i]).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_Y, EM, GLOBAL_TMP,system)
        
        if system: CODE+="}\n" 
        CODE+="} else { /* constant data */\n"
   if system: 
       CODE+= """for (k=0;k<p.numEqu;k++) {
            """ 
   DATA_Y=[]
   CODE2=""
   EM = {}
   for i in range(len(S)):
          EM[i] = 0 
   if 1:
               A_name="Y_0"
             A=Symbol(A_name)
             DATA_Y.append(A)
             if system:
                      CODE2+="const double %s = Y_p[k];\n"%(A_name,)
             else:
                      CODE2+="const double %s = Y_p[0];\n"%(A_name,)
             for q in range(len(Q)):
                   for i in range(len(S)):
                          EM[i] = EM[i] + (A * W[q] * S[i]).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_Y, EM, GLOBAL_TMP,system)    
   if system:  CODE+="}\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n"  
   
   for t in GLOBAL_TMP:
       PRECODE+="const double %s = %s;\n"%(GLOBAL_TMP[t],ccode(t.evalf(n=DIGITS)))
   return CODE, PRECODE

for d in [3]:
     generate(d)
