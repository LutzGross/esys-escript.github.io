"""
this script generates the assamblage routine for the ripley rectangular grid solver


"""
from sympy import *

DIGITS=20

#
#  quadrature form in [0,1]
#
q0=sqrt(RealNumber(15))/5
qD1=[ (1-q0)/2,  RealNumber(1)/2,  (1+q0)/2 ]
wD1=[RealNumber(5)/18, RealNumber(4)/ 9, RealNumber(5)/18 ]
qD1_r=[ RealNumber(1)/2 ]
wD1_r= [RealNumber(1)]
print "1D quadrature nodes =",qD1
print "1D quadrature weights =",wD1
print "1D reduced quadrature nodes =",qD1_r
print "1D reduced quadrature weights =",wD1_r

def generate(DIM):
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

   # interpolation to Quad points
   CODE="\nif (out_data_type==RIPLEY_ELEMENTS) {\n"
   CODE+=createCode(f_interpolation, x, Q, loopindex=idx_full)
   CODE+="} else if (out_data_type==RIPLEY_REDUCED_ELEMENTS) {\n"
   CODE+=createCode(f_interpolation, x, Q_r, loopindex=idx_full)
   CODE+="} else if (out_data_type==RIPLEY_BOUNDARY_ELEMENTS) {\n"
   for i in xrange(len(Q_r_faces)):
        CODE+="if (face_offset(%s)>-1) {\n"%i
        CODE+=createCode(f_interpolation, x, Q_faces[i], loopindex=idx_faces[i],  gridoffset="face_offset(%s)"%i)
        CODE+="\n} /* end of face %s */\n"%i
   CODE+="} else if (out_data_type==RIPLEY_REDUCED_BOUNDARY_ELEMENTS) {\n"
   for i in xrange(len(Q_r_faces)):
        CODE+="if (face_offset(%s)>-1) {\n"%i
        CODE+=createCode(f_interpolation, x, Q_r_faces[i], loopindex=idx_faces[i],  gridoffset="face_offset(%s)"%i)
        CODE+="\n} /* end of face %s */\n"%i
   CODE+="\n} /* end of out_data_type branching \n"
   insertCode("Assemble_Interpolation_%sD.c"%DIM, CODE)


   # interpolation to Quad points
   CODE="\nif (out_data_type==RIPLEY_ELEMENTS) {\n"
   CODE+=createGradientCode(f_interpolation, x, Q, loopindex=idx_full)
   CODE+="} else if (out_data_type==RIPLEY_REDUCED_ELEMENTS) {\n"
   CODE+=createGradientCode(f_interpolation, x, Q_r, loopindex=idx_full)
   CODE+="} else if (out_data_type==RIPLEY_BOUNDARY_ELEMENTS) {\n"
   for i in xrange(len(Q_r_faces)):
        CODE+="if (face_offset(%s)>-1) {\n"%i
        CODE+=createGradientCode(f_interpolation, x, Q_faces[i], loopindex=idx_faces[i],  gridoffset="face_offset(%s)"%i)
        CODE+="\n} /* end of face %s */\n"%i
   CODE+="} else if (out_data_type==RIPLEY_REDUCED_BOUNDARY_ELEMENTS) {\n"
   for i in xrange(len(Q_r_faces)):
        CODE+="if (face_offset(%s)>-1) {\n"%i
        CODE+=createGradientCode(f_interpolation, x, Q_r_faces[i], loopindex=idx_faces[i],  gridoffset="face_offset(%s)"%i)
        CODE+="\n} /* end of face %s */\n"%i
   CODE+="\n} /* end of out_data_type branching \n"
   insertCode("Assemble_Gradient_%sD.c"%DIM, CODE)
    

def insertCode(fn, code):
    TOP_LINE="/* GENERATOR SNIP TOP"
    BOTTOM_LINE="/* GENERATOR SNIP BOTTOM"
    TAG="   "

    bkp_stream=open(fn+".bkp",'w')
    new_code=""
    ignore=False
    for l in open(fn, 'r').readlines():
       s=l.strip()
       if len(s)>0:
          bkp_stream.write(s+"\n")
          if s.startswith(TOP_LINE): 
            ignore=True
            new_code+=s+"\n"+code
          elif s.startswith(BOTTOM_LINE): 
            ignore=False
            new_code+=s+"\n"
          elif not ignore:
            new_code+=s+"\n"
    bkp_stream.close()
    N=0
    new_code2=""
    for l in new_code.split("\n"):
       s=l.strip()
       if len(s)>0:
          if s.startswith("}"): 
             N=max(N-1, 0)
          new_code2+=TAG*N + s + "\n"
          if s.endswith("{"): 
             N+=1
    open(fn,'w').write(new_code2)

   

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
   k=""
   N=""
   M1=""
   for i in xrange(len(Q[0])):
     if len(k)==0:
         k+="k%s"%i
         N+="N%s"%i
         M1+="M%s"%i
     else:
         k+=",k%s"%i
         if i < DIM-1: N+=",N%s"%i
         if i < DIM-1: M1+=",M%s"%i
   
   k3=""
   for i in xrange(DIM-1,-1,-1): 
          if loopindex[i] > -1: k3+=",k%i"%loopindex[i]
   TXT=""
   for a,v in  consts.items():
      TXT+="const double %s = %s;\n"%(ccode(a), ccode(v.evalf(n=DIGITS)))
   TXT+="#pragma omp parallel for private(i%s)"%k3

   for i in xrange(DIM-1,-1,-1): 
          if loopindex[i] > -1: TXT+="\nfor (k%i =0; k%i < N%i; ++k%i) {"%(loopindex[i],loopindex[i],loopindex[i],loopindex[i])
   TXT+="\nfor (i =0; i < NCOMP; ++i) {\n"
   for s in args:
            k1=""
            for i in xrange(DIM):
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
            TXT+="register const double %s = in[INDEX2(i,INDEX%s(%s, %s),NCOMP)];\n"%(s.name,DIM,k1[1:],M1)
   #  interpolation to quadrature points
   for q in xrange(len(Q)):
         IDX2="INDEX%s(%s,%s)"%(DIM,k,N)
         if len(gridoffset) > 0: IDX2=gridoffset+"+"+IDX2
         for i in range(DIM):
             TXT+="out[INDEX4(i,%s,%s,%s,NCOMP,%s,%s)] = %s;\n"%(i,q,IDX2,DIM,len(Q),ccode(dF[i][q]))
   TXT+="} /* close component loop i */\n"
   for i in xrange(DIM): 
           if loopindex[i]>-1 : TXT+="} /* close k%i loop */\n"%loopindex[i]
   return TXT

def createCode(F, x, Q, gridoffset="", loopindex=(0,1,2)):
   DIM=len(loopindex)
   F, args, consts = optimizeEvaluate([F, ], x, Q)
   F=F[0]
   k=""
   N=""
   M1=""
   for i in xrange(len(Q[0])):
     if len(k)==0:
         k+="k%s"%i
         N+="N%s"%i
         M1+="M%s"%i
     else:
         k+=",k%s"%i
         if i < DIM-1: N+=",N%s"%i
         if i < DIM-1: M1+=",M%s"%i
   
   k3=""
   for i in xrange(DIM-1,-1,-1): 
          if loopindex[i] > -1: k3+=",k%i"%loopindex[i]
   TXT=""
   for a,v in  consts.items():
      TXT+="const double %s = %s;\n"%(ccode(a), ccode(v.evalf(n=DIGITS)))
   TXT+="#pragma omp parallel for private(i%s)"%k3
   for i in xrange(DIM-1,-1,-1): 
          if loopindex[i] > -1: TXT+="\nfor (k%i =0; k%i < N%i; ++k%i) {"%(loopindex[i],loopindex[i],loopindex[i],loopindex[i])
   TXT+="\nfor (i =0; i < NCOMP; ++i) {\n"
   for s in args:
            k1=""
            for i in xrange(DIM):
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
            TXT+="register const double %s = in[INDEX2(i,INDEX%s(%s, %s),NCOMP)];\n"%(s.name,DIM,k1[1:],M1)
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


for d in [2,3 ]:
     generate(d)
