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

qD1_r=[ RealNumber(1)/2 ]
wD1_r= [RealNumber(1)]

print "1D quadrature nodes =",qD1
print "1D quadrature weights =",wD1
print "1D reduced quadrature nodes =",qD1_r
print "1D reduced quadrature weights =",wD1_r

def generate(DIM, filename):
   # variables
   h=[Symbol('h%s'%i) for i in xrange(DIM) ]
   x=[Symbol('x%s'%i) for i in xrange(DIM) ]
   #  shape functions: labeling differently from finley!!!!
   S=[ 1 ]
   GridOffset=[[]]
   GridOffsetText=[ "" ]
   idx_full=[]
   Q = [ [] ]
   W = [ 1 ]
   Q_r = [ [] ]
   W_r = [ 1 ]
   for i in xrange(DIM):
       idx_full.append(i)
       S = [ s * (1-x[i]/h[i])  for s in S] + [ s * x[i]/h[i]  for s in S]
       GridOffset = [ s + [0]  for s in GridOffset] + [ s + [1]  for s in GridOffset]
       GridOffsetText = [ s + "0"  for s in GridOffsetText] + [ s + "1"  for s in GridOffsetText]

       # generate quadrature points in element
       Qnew, Wnew =[], []
       for j in xrange(len(qD1)):
              Qnew += [ q + [qD1[j]*h[i],] for q in Q ]
              Wnew += [ w * wD1[j]*h[i] for w in W ]
       Q, W=Qnew, Wnew

       Qnew, Wnew =[], []
       for j in xrange(len(qD1_r)):
               Qnew += [ q + [qD1_r[j]*h[i],] for q in Q_r ]
               Wnew += [ w * wD1_r[j]*h[i] for w in W_r ]
       Q_r, W_r=Qnew, Wnew

   # now the same thing for the faces:
   Q_faces=[]
   W_faces=[]
   Q_r_faces=[]
   W_r_faces=[]
   idx_faces=[]
   for face in xrange(DIM):
       Q_left, W_left = [ [] ], [ 1 ]
       Q_left_r, W_left_r = [ [] ], [ 1 ]
       Q_right, W_right = [ [] ], [ 1 ]
       Q_right_r, W_right_r = [ [] ], [ 1 ]
       idx_left,idx_right=[],[]
       for i in xrange(DIM):
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
             for j in xrange(len(qD1)):
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
       for i in xrange(DIM):
          # generate quadrature points in element
          Q_left_new, W_left_new =[], []
          Q_right_new, W_right_new =[], []
          if face == i :
             Q_left_new += [ q + [0.,] for q in Q_left ]
             W_left_new += [ w * 1. for w in W_left ]
             Q_right_new += [ q + [ h[i], ] for q in Q_right ]
             W_right_new += [ w * 1. for w in W_right ]
          else:
             for j in xrange(len(qD1_r)):
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
   for i in xrange(len(GridOffsetText)):
     f_interpolation = f_interpolation + S[i] * Symbol("f_%s"%GridOffsetText[i])

   # integration to Quad points
   #CODE="\nif (out_data_type==RIPLEY_ELEMENTS) {\n"
   #CODE+=createIntegrationCode(Q, W, loopindex=idx_full)
   #CODE+="} else if (out_data_type==RIPLEY_REDUCED_ELEMENTS) {\n"
   #CODE+=createIntegrationCode(Q_r, W_r, loopindex=idx_full)
   #CODE+="} else if (out_data_type==RIPLEY_BOUNDARY_ELEMENTS) {\n"
   #for i in xrange(len(Q_faces)):
   #     CODE+="if (face_offset(%s)>-1) {\n"%i
   #     CODE+=createIntegrationCode(Q_faces[i], W_faces[i], loopindex=idx_faces[i],  gridoffset="face_offset(%s)"%i)
   #     CODE+="\n} /* end of face %s */\n"%i
   #CODE+="} else if (out_data_type==RIPLEY_REDUCED_BOUNDARY_ELEMENTS) {\n"
   #for i in xrange(len(Q_r_faces)):
   #     CODE+="if (face_offset(%s)>-1) {\n"%i
   #     CODE+=createIntegrationCode(Q_r_faces[i], W_r_faces[i], loopindex=idx_faces[i],  gridoffset="face_offset(%s)"%i)
   #     CODE+="\n} /* end of face %s */\n"%i
   #CODE+="\n} /* end of out_data_type branching */\n"
   #insertCode("Assemble_Integration_%sD.c"%DIM, { "SNIP" : CODE})
   #1/0

   # interpolation to Quad points
   #CODE=createCode(f_interpolation, x, Q, loopindex=idx_full)
   #insertCode(filename, { "SNIP_INTERPOLATE_ELEMENTS" : CODE})
   #CODE=createCode(f_interpolation, x, Q_r, loopindex=idx_full)
   #insertCode(filename, { "SNIP_INTERPOLATE_REDUCED_ELEMENTS" : CODE})
   #CODE=""
   #for i in xrange(len(Q_faces)):
   #     CODE+="if (m_faceOffset[%s] > -1) {\n"%i
   #     CODE+=createCode(f_interpolation, x, Q_faces[i], loopindex=idx_faces[i],  gridoffset="m_faceOffset[%s]"%i)
   #     CODE+="\n} /* end of face %s */\n"%i
   #insertCode(filename, { "SNIP_INTERPOLATE_FACES" : CODE})
   #CODE=""
   #for i in xrange(len(Q_r_faces)):
   #     CODE+="if (m_faceOffset[%s] > -1) {\n"%i
   #     CODE+=createCode(f_interpolation, x, Q_r_faces[i], loopindex=idx_faces[i],  gridoffset="m_faceOffset[%s]"%i)
   #     CODE+="\n} /* end of face %s */\n"%i
   #insertCode(filename, { "SNIP_INTERPOLATE_REDUCED_FACES" : CODE})

   ## gradient to Quad points
   #CODE=createGradientCode(f_interpolation, x, Q, loopindex=idx_full)
   #insertCode(filename, { "SNIP_GRAD_ELEMENTS" : CODE})
   #CODE=createGradientCode(f_interpolation, x, Q_r, loopindex=idx_full)
   #insertCode(filename, { "SNIP_GRAD_REDUCED_ELEMENTS" : CODE})
   #CODE=""
   #for i in xrange(len(Q_faces)):
   #     CODE+="if (m_faceOffset[%s] > -1) {\n"%i
   #     CODE+=createGradientCode(f_interpolation, x, Q_faces[i], loopindex=idx_faces[i],  gridoffset="m_faceOffset[%s]"%i)
   #     CODE+="\n} /* end of face %s */\n"%i
   #insertCode(filename, { "SNIP_GRAD_FACES" : CODE})
   #CODE=""
   #for i in xrange(len(Q_r_faces)):
   #     CODE+="if (m_faceOffset[%s] > -1) {\n"%i
   #     CODE+=createGradientCode(f_interpolation, x, Q_r_faces[i], loopindex=idx_faces[i],  gridoffset="m_faceOffset[%s]"%i)
   #     CODE+="\n} /* end of face %s */\n"%i
   #insertCode(filename, { "SNIP_GRAD_REDUCED_FACES" : CODE})

   #generate PDE assemblage
   CODE,PRECODE = makePDE(S, x, Q, W, DIM=DIM, system=True)
   insertCode(filename, { "SNIP_PDE_SYSTEM" : CODE , "SNIP_PDE_SYSTEM_PRE" : PRECODE } )
   CODE,PRECODE = makePDE(S, x, Q_r, W_r, DIM=DIM, system=True)
   insertCode(filename, { "SNIP_PDE_SYSTEM_REDUCED" : CODE , "SNIP_PDE_SYSTEM_REDUCED_PRE" : PRECODE })
   CODE,PRECODE = makePDE(S, x, Q, W, DIM=DIM, system=False)
   insertCode(filename, { "SNIP_PDE_SINGLE" : CODE , "SNIP_PDE_SINGLE_PRE" : PRECODE })
   CODE,PRECODE = makePDE(S, x, Q_r, W_r, DIM=DIM, system=False)
   insertCode(filename, { "SNIP_PDE_SINGLE_REDUCED" : CODE , "SNIP_PDE_SINGLE_REDUCED_PRE" : PRECODE })




def insertCode(fn, replacement):
    TOP_LINE_FMT="/* GENERATOR %s TOP"
    BOTTOM_LINE_FMT="/* GENERATOR %s BOTTOM"
    TABSIZE=4

    # read code and create back up
    bkp_stream=open(fn+".bkp",'w')
    old_code=""
    for l in open(fn, 'r').readlines():
       s=l
       if len(s)>0:
          bkp_stream.write(s)
          old_code+=s
    bkp_stream.close()

    for k in replacement.keys():
        TOP_LINE=TOP_LINE_FMT%k
        BOTTOM_LINE=BOTTOM_LINE_FMT%k
        new_code=""
        ignore=False
        for l in old_code.split("\n"):
            if l.strip().startswith(TOP_LINE):
                ignore=True
                N=len(l[:l.find(TOP_LINE)])
                new_code+=l+"\n"
                for rep_l in replacement[k].split("\n"):
                    #insert tabs:
                    if rep_l.find("/*")>-1 and rep_l.find("*/")>-1:
                        s=(rep_l[:rep_l.find("/*")] + rep_l[rep_l.find("*/")+2:]).strip()
                    else:
                        s=rep_l.strip()
                    if len(rep_l)>0:
                        if s.startswith("}"):
                            N=max(N-TABSIZE, 0)
                        if s.startswith("#"):
                            new_code+=rep_l + "\n"
                        else:
                            new_code+=" "*N + rep_l + "\n"
                        if s.endswith("{"):
                            N+=TABSIZE
            elif l.strip().startswith(BOTTOM_LINE):
                ignore=False
                new_code+=l+"\n"
            elif not ignore:
                new_code+=l+"\n"
        old_code=new_code
    new_code=old_code[:len(old_code)-1]
    open(fn,'w').write(new_code)
    print "file %s updated."%fn

def optimizeEvaluate(F, x, Q):
  F2=[]
  for f in F:
     f0=[]
     for q in xrange(len(Q)):
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
   k=[]
   N=[]
   M1=""
   TXT=""

   for i in xrange(len(Q[0])):
     if loopindex[i]>-1:
         k.append("k%s"%i)
         N.append("m_NE%s"%i)
     if len(M1)==0:
         M1+="m_N%s"%i
     else:
         if i < DIM-1: M1+=",m_N%s"%i

   if len(k)>=2:
       IDX2="INDEX%s(%s,%s)"%(len(k),",".join(k),",".join(N[:-1]))
   else:
       IDX2=k[0]

   #for i in xrange(DIM-1,-1,-1):
   #       if loopindex[i] <= -1: TXT+="const index_t k%i = 0;\n"%i
   for a,v in  consts.items():
      TXT+="const double %s = %s;\n"%(ccode(a), ccode(v.evalf(n=DIGITS)))
   TXT+="#pragma omp parallel for"

   for i in xrange(DIM-1,-1,-1):
          if loopindex[i] > -1: TXT+="\nfor (index_t k%i=0; k%i < m_NE%i; ++k%i) {"%(loopindex[i],loopindex[i],loopindex[i],loopindex[i])
   for s in args:
            k1=""
            for i in xrange(DIM):
              if s.name[2+i]=="0":
                  if loopindex[i]>-1:
                     k1+=",k%s"%i
                  elif  loopindex[i] == -1:
                     k1+=",0"
                  elif  loopindex[i] == -2:
                     k1+=",m_N%s-2"%i
              else:
                  if loopindex[i]>-1:
                     k1+=",k%s+1"%i
                  elif  loopindex[i] == -1:
                     k1+=",1"
                  elif  loopindex[i] == -2:
                     k1+=",m_N%s-1"%i
            TXT+="\nconst register double* %s = in.getSampleDataRO(INDEX%s(%s, %s));"%(s.name,DIM,k1[1:],M1)
   if len(gridoffset) > 0: IDX2=gridoffset+"+"+IDX2
   TXT+="\ndouble* o = out.getSampleDataRW(%s);"%IDX2
   TXT+="\nfor (index_t i=0; i < numComp; ++i) {\n"
   for q in xrange(len(Q)):
       for i in range(DIM):
           dFidx=dF[i][q]
           for s in args:
               dFidx=dFidx.subs(s.name, Symbol(s.name+"[i]"))
           TXT+="o[INDEX3(i,%s,%s,numComp,%s)] = %s;\n"%(i,q,DIM,ccode(dFidx))
           #TXT+="out[INDEX4(i,%s,%s,%s,numComp,%s,%s)] = %s;\n"%(i,q,IDX2,DIM,len(Q),ccode(dF[i][q]))
   TXT+="} /* end of component loop i */\n"
   for i in xrange(DIM):
           if loopindex[i]>-1 : TXT+="} /* end of k%i loop */\n"%loopindex[i]
   return TXT

def createCode(F, x, Q, gridoffset="", loopindex=(0,1,2)):
   DIM=len(loopindex)
   F, args, consts = optimizeEvaluate([F, ], x, Q)
   F=F[0]
   k=[]
   N=[]
   M1=""
   TXT=""
   for i in xrange(len(Q[0])):
     if loopindex[i]>-1:
         k.append("k%s"%i)
         N.append("m_NE%s"%i)
     if len(M1)==0:
         M1+="m_N%s"%i
     else:
         if i < DIM-1: M1+=",m_N%s"%i

   if len(k)>=2:
       IDX2="INDEX%s(%s,%s)"%(len(k),",".join(k),",".join(N[:-1]))
   else:
       IDX2=k[0]

   #for i in xrange(DIM-1,-1,-1):
   #       if loopindex[i] <= -1: TXT+="const index_t k%i = 0;\n"%i
   for a,v in  consts.items():
      TXT+="const double %s = %s;\n"%(ccode(a), ccode(v.evalf(n=DIGITS)))
   TXT+="#pragma omp parallel for"
   for i in xrange(DIM-1,-1,-1):
          if loopindex[i] > -1: TXT+="\nfor (index_t k%i=0; k%i < m_NE%i; ++k%i) {"%(loopindex[i],loopindex[i],loopindex[i],loopindex[i])
   for s in args:
        k1=""
        for i in xrange(DIM):
          if s.name[2+i]=="0":
              if loopindex[i]>-1:
                 k1+=",k%s"%i
              elif  loopindex[i] == -1:
                 k1+=",0"
              elif  loopindex[i] == -2:
                 k1+=",m_N%s-2"%i
          else:
              if loopindex[i]>-1:
                 k1+=",k%s+1"%i
              elif  loopindex[i] == -1:
                 k1+=",1"
              elif  loopindex[i] == -2:
                 k1+=",m_N%s-1"%i
        TXT+="\nconst register double* %s = in.getSampleDataRO(INDEX%s(%s, %s));"%(s.name,DIM,k1[1:],M1)
   if len(gridoffset) > 0: IDX2=gridoffset+"+"+IDX2
#        for i in xrange(DIM):
#             if loopindex[i]>-1: IDX2=gridoffset+"+k%s"%i
   TXT+="\ndouble* o = out.getSampleDataRW(%s);"%IDX2
   TXT+="\nfor (index_t i=0; i < numComp; ++i) {\n"
   #  interpolation to quadrature points
   for q in xrange(len(Q)):
       Fidx=F[q]
       for s in args:
           Fidx=Fidx.subs(s.name, Symbol(s.name+"[i]"))
       TXT+="o[INDEX2(i,numComp,%s)] = %s;\n"%(q,ccode(Fidx))
   TXT+="} /* end of component loop i */\n"
   for i in xrange(DIM):
           if loopindex[i]>-1 : TXT+="} /* end of k%i loop */\n"%loopindex[i]
   return TXT

def createIntegrationCode(Q, W, gridoffset="", loopindex=(0,1,2)):
   DIM=len(loopindex)
   GLOBAL_TMP={}
   LOCAL_TMP_E={}
   TXT_E=""
   for q in range(len(W)):
        if not GLOBAL_TMP.has_key(W[q]):
             GLOBAL_TMP[W[q]]=Symbol("w_%s"%len(GLOBAL_TMP))
        A=Symbol("f_%s"%q)
        TXT_E+="const register double %s = in[INDEX3(i,%s,e, NCOMP,%s)];\n"%(A,q,len(W))
        if not LOCAL_TMP_E.has_key(GLOBAL_TMP[W[q]]):
             LOCAL_TMP_E[GLOBAL_TMP[W[q]]]=0
        LOCAL_TMP_E[GLOBAL_TMP[W[q]]]=LOCAL_TMP_E[GLOBAL_TMP[W[q]]]+A
   for p in LOCAL_TMP_E.keys():
        TXT_E+="const register double %s = %s;\n"%( p, ccode(LOCAL_TMP_E[p]))
   print TXT_E
   print GLOBAL_TMP
   return ""
   1/0
   #=================



   k3=""
   for i in xrange(DIM-1,-1,-1):
          if loopindex[i] > -1: k3+=",k%i"%loopindex[i]
   TXT=""
   for a,v in  consts.items():
      TXT+="const double %s = %s;\n"%(ccode(a), ccode(v.evalf(n=DIGITS)))
   TXT+="#pragma omp parallel for private(i%s)"%k3
   for i in xrange(DIM-1,-1,-1):
          if loopindex[i] > -1: TXT+="\nfor (k%i =k%i_0; k%i < N%i; ++k%i) {"%(loopindex[i],loopindex[i], loopindex[i],loopindex[i],loopindex[i])
   TXT+="\nfor (i =0; i < NCOMP; ++i) {\n"
   print TXT
   1/0
   #  interpolation to quadrature points
   for q in xrange(len(Q)):
      IDX2="INDEX%s(%s,%s)"%(DIM,k,N)
      if len(gridoffset) > 0: IDX2=gridoffset+"+"+IDX2
      TXT+="out[INDEX3(i,%s,%s,NCOMP,%s)] = %s;\n"%(q,IDX2,len(Q),ccode(F[q]))
   TXT+="} /* close component loop i */\n"
   for i in xrange(DIM):
           if loopindex[i]>-1 : TXT+="} /* close k%i loop */\n"%loopindex[i]
   return TXT

def subPoint(g,x, p):
     out=g
     for i in range(len(x)): out=out.subs(x[i], p[i])
     # return out.evalf()
     return out

def generatePDECode(DATA_A, EM, GLOBAL_TMP, system=False):
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
                 OUT+="EM_F[INDEX2(k,%s,numEq)]+=%s;\n"%(p,ccode(sum_em))
              else:
                 OUT+="EM_F[%s]+=%s;\n"%(p,ccode(sum_em))
            else:
              if system:
                 OUT+="EM_S[INDEX4(k,m,%s,%s,numEq,numComp,%s)]+=%s;\n"%(p[0],p[1],max([ qq[0] for qq in EM.keys()])+1,ccode(sum_em))
              else:
                 OUT+="EM_S[INDEX2(%s,%s,%s)]+=%s;\n"%(p[0],p[1],max([ qq[0] for qq in EM.keys()])+1,ccode(sum_em))

        OUT2=""
        for p in LOCAL_TMP:
             OUT2+="const register double %s = %s;\n"%(LOCAL_TMP[p],ccode(p))
        for p in LOCAL2_TMP:
             OUT2+="const register double %s = %s;\n"%(LOCAL2_TMP[p],ccode(p))
        return OUT2+OUT

def makePDE(S, x, Q, W, DIM=2, system=False):
   GLOBAL_TMP={}
   GLOBAL_N=0
   PRECODE=""
   CODE="""
///////////////
// process A //
///////////////
if (!A.isEmpty()) {
add_EM_S=true;
const double* A_p=const_cast<escript::Data*>(&A)->getSampleDataRO(e);
"""
   if len(Q) > 1:
        CODE+="if (A.actsExpanded()) {\n"
        if system: CODE+= """for (index_t k=0; k<numEq; k++) {
for (index_t m=0; m<numComp; m++) {
"""
        DATA_A=[]
        CODE2=""
        EM = {}
        for i in range(len(S)):
            for j in range(len(S)):
                EM[(i,j)] = 0

        for q in xrange(len(Q)):
            for di in range(DIM):
               for dj in range(DIM):
                  A_name="A_%d%d_%d"%(di,dj,q)
                  A=Symbol(A_name)
                  DATA_A.append(A)
                  if system:
                      CODE2+="const register double %s = A_p[INDEX5(k,%s,m,%s,%s, numEq,%s,numComp,%s)];\n"%(A_name, di, dj, q, DIM, DIM)
                  else:
                      CODE2+="const register double %s = A_p[INDEX3(%s,%s,%s,%s,%s)];\n"%(A_name, di, dj, q, DIM, DIM)
                  for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[i],x[di]) * diff(S[j],x[dj])).subs( [ (x[jj], Q[q][jj]) for jj in xrange(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_A, EM, GLOBAL_TMP, system)

        if system: CODE+="}\n}\n"
        CODE+="} else { /* constant data */\n"
   if system:
       CODE+= """for (index_t k=0; k<numEq; k++) {
for (index_t m=0; m<numComp; m++) {
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
                      CODE2+="const register double %s = A_p[INDEX4(k,%s,m,%s, numEq,%s, numComp)];\n"%(A_name, di, dj, DIM)
             else:
                      CODE2+="const register double %s = A_p[INDEX2(%s,%s,%s)];\n"%(A_name, di, dj, DIM)
             for q in xrange(len(Q)):
                  for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[i],x[di]) * diff(S[j],x[dj])).subs( [ (x[jj], Q[q][jj]) for jj in xrange(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_A, EM, GLOBAL_TMP, system)
   if system:  CODE+="}\n}\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n"
   #BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
   CODE+="""
///////////////
// process B //
///////////////
if (!B.isEmpty()) {
add_EM_S=true;
const double* B_p=const_cast<escript::Data*>(&B)->getSampleDataRO(e);
"""
   if len(Q) > 1:
        CODE+="if (B.actsExpanded()) {\n"
        if system: CODE+= """for (index_t k=0; k<numEq; k++) {
for (index_t m=0; m<numComp; m++) {
"""
        DATA_B=[]
        CODE2=""
        EM = {}
        for i in range(len(S)):
            for j in range(len(S)):
                EM[(i,j)] = 0

        for q in xrange(len(Q)):
            for di in range(DIM):
                  A_name="B_%d_%d"%(di,q)
                  A=Symbol(A_name)
                  DATA_B.append(A)
                  if system:
                      CODE2+="const register double %s = B_p[INDEX4(k,%s,m,%s, numEq,%s,numComp)];\n"%(A_name, di,  q, DIM)
                  else:
                      CODE2+="const register double %s = B_p[INDEX2(%s,%s,%s)];\n"%(A_name, di, q, DIM)
                  for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[i],x[di]) * S[j]).subs( [ (x[jj], Q[q][jj]) for jj in xrange(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_B, EM, GLOBAL_TMP, system)

        if system: CODE+="}\n}\n"
        CODE+="} else { /* constant data */\n"
   if system:
       CODE+= """for (index_t k=0; k<numEq; k++) {
for (index_t m=0; m<numComp; m++) {
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
                      CODE2+="const register double %s = B_p[INDEX3(k,%s,m, numEq, %s)];\n"%(A_name, di,  DIM)
             else:
                      CODE2+="const register double %s = B_p[%s];\n"%(A_name, di)
             for q in xrange(len(Q)):
                  for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[i],x[di]) * S[j] ).subs( [ (x[jj], Q[q][jj]) for jj in xrange(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_B, EM, GLOBAL_TMP, system)
   if system:  CODE+="}\n}\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n"
   #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   CODE+="""
///////////////
// process C //
///////////////
if (!C.isEmpty()) {
add_EM_S=true;
const double* C_p=const_cast<escript::Data*>(&C)->getSampleDataRO(e);
"""
   if len(Q) > 1:
        CODE+="if (C.actsExpanded()) {\n"
        if system: CODE+= """for (index_t k=0; k<numEq; k++) {
for (index_t m=0; m<numComp; m++) {
"""
        DATA_C=[]
        CODE2=""
        EM = {}
        for i in range(len(S)):
            for j in range(len(S)):
                EM[(i,j)] = 0

        for q in xrange(len(Q)):
            for dj in range(DIM):
                  A_name="C_%d_%d"%(dj,q)
                  A=Symbol(A_name)
                  DATA_C.append(A)
                  if system:
                      CODE2+="const register double %s = C_p[INDEX4(k,m,%s, %s, numEq,numComp,%s)];\n"%(A_name, dj,  q, DIM)
                  else:
                      CODE2+="const register double %s = C_p[INDEX2(%s,%s,%s)];\n"%(A_name, dj, q, DIM)
                  for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[j],x[dj]) * S[i]).subs( [ (x[jj], Q[q][jj]) for jj in xrange(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_C, EM, GLOBAL_TMP, system)

        if system: CODE+="}\n}\n"
        CODE+="} else { /* constant data */\n"
   if system:
       CODE+= """for (index_t k=0; k<numEq; k++) {
for (index_t m=0; m<numComp; m++) {
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
                      CODE2+="const register double %s = C_p[INDEX3(k, m, %s, numEq, numComp)];\n"%(A_name, dj)
             else:
                      CODE2+="const register double %s = C_p[%s];\n"%(A_name, dj)
             for q in xrange(len(Q)):
                  for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[j],x[dj]) * S[i] ).subs( [ (x[jj], Q[q][jj]) for jj in xrange(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_C, EM, GLOBAL_TMP, system)
   if system: CODE+="}\n}\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n"
   #DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
   CODE+="""
///////////////
// process D //
///////////////
if (!D.isEmpty()) {
add_EM_S=true;
const double* D_p=const_cast<escript::Data*>(&D)->getSampleDataRO(e);
"""
   if len(Q) > 1:
        CODE+="if (D.actsExpanded()) {\n"
        if system: CODE+= """for (index_t k=0; k<numEq; k++) {
for (index_t m=0; m<numComp; m++) {
"""
        DATA_D=[]
        CODE2=""
        EM = {}
        for i in range(len(S)):
            for j in range(len(S)):
                EM[(i,j)] = 0

        for q in xrange(len(Q)):
                  A_name="D_%d"%(q,)
                  A=Symbol(A_name)
                  DATA_D.append(A)
                  if system:
                      CODE2+="const register double %s = D_p[INDEX3(k, m, %s, numEq, numComp)];\n"%(A_name, q)
                  else:
                      CODE2+="const register double %s = D_p[%s];\n"%(A_name, q)
                  for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * S[j] * S[i]).subs( [ (x[jj], Q[q][jj]) for jj in xrange(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_D, EM, GLOBAL_TMP, system)

        if system: CODE+="}\n }\n"
        CODE+="} else { /* constant data */\n"
   if system:
       CODE+= """for (index_t k=0; k<numEq; k++) {
for (index_t m=0; m<numComp; m++) {
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
                      CODE2+="const register double %s = D_p[INDEX2(k, m, numEq)];\n"%(A_name,)
             else:
                      CODE2+="const register double %s = D_p[0];\n"%(A_name,)
             for q in xrange(len(Q)):
                  for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * S[j] * S[i] ).subs( [ (x[jj], Q[q][jj]) for jj in xrange(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_D, EM, GLOBAL_TMP, system)
   if system: CODE+="}\n}\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n"

   #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   CODE+="""
///////////////
// process X //
///////////////
if (!X.isEmpty()) {
add_EM_F=true;
const double* X_p=const_cast<escript::Data*>(&X)->getSampleDataRO(e);
"""
   if len(Q) > 1:
        CODE+="if (X.actsExpanded()) {\n"
        if system: CODE+= "for (index_t k=0; k<numEq; k++) {\n"
        DATA_X=[]
        CODE2=""
        EM = {}
        for i in range(len(S)):
                EM[i] = 0

        for q in xrange(len(Q)):
            for dj in range(DIM):
                  A_name="X_%d_%d"%(dj,q)
                  A=Symbol(A_name)
                  DATA_X.append(A)
                  if system:
                      CODE2+="const register double %s = X_p[INDEX3(k, %s, %s, numEq, %s)];\n"%(A_name, dj,  q, DIM)
                  else:
                      CODE2+="const register double %s = X_p[INDEX2(%s,%s,%s)];\n"%(A_name, dj,q,DIM)
                  for j in range(len(S)):
                          EM[j] = EM[j] + (A * W[q] * diff(S[j],x[dj])).subs( [ (x[jj], Q[q][jj]) for jj in xrange(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_X, EM, GLOBAL_TMP, system)

        if system: CODE+="}\n"
        CODE+="} else { /* constant data */\n"
   if system:
       CODE+= "for (index_t k=0; k<numEq; k++) {\n"
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
                      CODE2+="const register double %s = X_p[INDEX2(k, %s, numEq)];\n"%(A_name, dj)
             else:
                      CODE2+="const register double %s = X_p[%s];\n"%(A_name, dj)
             for q in xrange(len(Q)):
                  for j in range(len(S)):
                          EM[j] = EM[j] + (A * W[q] * diff(S[j],x[dj])).subs( [ (x[jj], Q[q][jj]) for jj in xrange(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_X, EM, GLOBAL_TMP, system)
   if system:  CODE+="}\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n"

   #YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYy
   CODE+="""
///////////////
// process Y //
///////////////
if (!Y.isEmpty()) {
add_EM_F=true;
const double* Y_p=const_cast<escript::Data*>(&Y)->getSampleDataRO(e);
"""
   if len(Q) > 1:
        CODE+="if (Y.actsExpanded()) {\n"
        if system: CODE+= "for (index_t k=0; k<numEq; k++) {\n"
        DATA_Y=[]
        CODE2=""
        EM = {}
        for i in range(len(S)):
                EM[i] = 0

        for q in xrange(len(Q)):
                  A_name="Y_%d"%(q,)
                  A=Symbol(A_name)
                  DATA_Y.append(A)
                  if system:
                      CODE2+="const register double %s = Y_p[INDEX2(k, %s, numEq)];\n"%(A_name, q)
                  else:
                      CODE2+="const register double %s = Y_p[%s];\n"%(A_name, q)
                  for i in range(len(S)):
                          EM[i] = EM[i] + (A * W[q] * S[i]).subs( [ (x[jj], Q[q][jj]) for jj in xrange(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_Y, EM, GLOBAL_TMP, system)

        if system: CODE+="}\n"
        CODE+="} else { /* constant data */\n"
   if system:
       CODE+= "for (index_t k=0; k<numEq; k++) {\n"
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
                      CODE2+="const register double %s = Y_p[k];\n"%(A_name,)
             else:
                      CODE2+="const register double %s = Y_p[0];\n"%(A_name,)
             for q in xrange(len(Q)):
                  for i in range(len(S)):
                          EM[i] = EM[i] + (A * W[q] * S[i]).subs( [ (x[jj], Q[q][jj]) for jj in xrange(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_Y, EM, GLOBAL_TMP, system)
   if system:  CODE+="}\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n"

   import operator
   for k,v in sorted(GLOBAL_TMP.iteritems(), key=operator.itemgetter(1)):
       PRECODE+="const double %s = %s;\n"%(v,ccode(k.evalf(n=DIGITS)))
   return CODE, PRECODE

filenames={2:"Rectangle.cpp", 3:"Brick.cpp"}
for d in filenames.keys():
     generate(d, filenames[d])

