
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
This script generates the assemblage routine for the ripley rectangular grid
solver.
"""


from multiprocessing import Process
from sympy import *

DIGITS=20

#
#  quadrature form in [0,1]
#
q0 = 1/sqrt(RealNumber(3))
qD1 = [ (1-q0)/2,  (1+q0)/2 ]
wD1 = [RealNumber(1)/2, RealNumber(1)/2 ]

##qD1_r = [ RealNumber(1)/2 ]
##wD1_r = [ RealNumber(1) ]

print("1D quadrature nodes = %s"%qD1)
print("1D quadrature weights = %s"%wD1)
##print("1D reduced quadrature nodes = %s"%qD1_r)
##print("1D reduced quadrature weights = %s"%wD1_r)

def generateAll(DIM, filename):
   # variables
   h = [Symbol('m_dx[%s]'%i) for i in range(DIM)]
   x = [Symbol('x%s'%i) for i in range(DIM)]
   #  shape functions: labeling differently from finley!!!!
   S = [1]
   GridOffset=[[]]
   GridOffsetText=[""]
   idx_full=[]
   Q = [ [] ]
   W = [ 1 ]
   for i in range(DIM):
        idx_full.append(i)
        S = [ s * (1-x[i]/h[i])  for s in S] + [ s * x[i]/h[i]  for s in S]
        GridOffset = [ s + [0]  for s in GridOffset] + [ s + [1]  for s in GridOffset]
        GridOffsetText = [ s + "0"  for s in GridOffsetText] + [ s + "1"  for s in GridOffsetText]

        # generate quadrature points in element
        Qnew, Wnew = [], []
        for j in range(len(qD1)):
            Qnew += [ q + [qD1[j]*h[i],] for q in Q ]
            Wnew += [ w * wD1[j]*h[i] for w in W ]
        Q, W=Qnew, Wnew

        # reduced parts are hand-optimized
        ##Qnew, Wnew = [], []
        ##for j in range(len(qD1_r)):
        ##        Qnew += [ q + [qD1_r[j]*h[i],] for q in Q_r ]
        ##        Wnew += [ w * wD1_r[j]*h[i] for w in W_r ]
        ##Q_r, W_r=Qnew, Wnew

   # now the same thing for the faces:
   Q_faces=[]
   W_faces=[]
   ##Q_r_faces=[]
   ##W_r_faces=[]
   idx_faces=[]
   for face in range(DIM):
        Q_left, W_left = [ [] ], [ 1 ]
        Q_right, W_right = [ [] ], [ 1 ]
        idx_left,idx_right=[],[]
        for i in range(DIM):
            # generate quadrature points in element
            Q_left_new, W_left_new =[], []
            Q_right_new, W_right_new =[], []

            if face == i :
                idx_left.append(-1)
                idx_right.append(-2)
                Q_left_new += [ q + [0,] for q in Q_left ]
                W_left_new += [ w * 1 for w in W_left ]
                Q_right_new += [ q + [ h[i], ] for q in Q_right ]
                W_right_new += [ w * 1 for w in W_right ]
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

        ##Q_left, W_left = [ [] ], [ 1 ]
        ##Q_right, W_right = [ [] ], [ 1 ]
        ##for i in range(DIM):
        ##   # generate quadrature points in element
        ##   Q_left_new, W_left_new =[], []
        ##   Q_right_new, W_right_new =[], []
        ##   if face == i :
        ##      Q_left_new += [ q + [0.,] for q in Q_left ]
        ##      W_left_new += [ w * 1. for w in W_left ]
        ##      Q_right_new += [ q + [ h[i], ] for q in Q_right ]
        ##      W_right_new += [ w * 1. for w in W_right ]
        ##   else:
        ##      for j in range(len(qD1_r)):
        ##          Q_left_new += [ q + [qD1_r[j]*h[i],] for q in Q_left ]
        ##          W_left_new += [ w * wD1_r[j]*h[i] for w in W_left ]
        ##          Q_right_new += [ q + [qD1_r[j]*h[i],] for q in Q_right ]
        ##          W_right_new += [ w * wD1_r[j]*h[i] for w in W_right ]
        ##   Q_left, W_left=Q_left_new, W_left_new
        ##   Q_right, W_right=Q_right_new, W_right_new
        ##Q_r_faces=Q_r_faces+[Q_left, Q_right]
        ##W_r_faces=W_r_faces+[W_left, W_right]

   #generate PDE assemblage
   #CODE,PRECODE = makePDE(S, x, Q, W, DIM=DIM, system=False)
   #insertCode(filename, { "SNIP_PDE_SINGLE" : CODE , "SNIP_PDE_SINGLE_PRE" : PRECODE })
   CODE,PRECODE = makePDE(S, x, Q, W, DIM=DIM, system=True)
   insertCode(filename, { "SNIP_PDE_SYSTEM" : CODE , "SNIP_PDE_SYSTEM_PRE" : PRECODE } )

   ##CODE,PRECODE = makePDE(S, x, Q_r, W_r, DIM=DIM, system=False)
   ##insertCode(filename, { "SNIP_PDE_SINGLE_REDUCED" : CODE , "SNIP_PDE_SINGLE_REDUCED_PRE" : PRECODE })
   ##CODE,PRECODE = makePDE(S, x, Q_r, W_r, DIM=DIM, system=True)
   ##insertCode(filename, { "SNIP_PDE_SYSTEM_REDUCED" : CODE , "SNIP_PDE_SYSTEM_REDUCED_PRE" : PRECODE })

   #generate PDEBoundary assemblage
   #CODE,PRECODE = makePDEBC(S, x, Q_faces, W_faces, DIM=DIM, system=True)
   #insertCode(filename, extendDictionary( { "SNIP_PDEBC_SYSTEM_PRE" : PRECODE }, "SNIP_PDEBC_SYSTEM", CODE))
   #CODE,PRECODE = makePDEBC(S, x, Q_faces, W_faces, DIM=DIM, system=False)
   #insertCode(filename, extendDictionary( { "SNIP_PDEBC_SINGLE_PRE" : PRECODE }, "SNIP_PDEBC_SINGLE", CODE ))
   ##CODE,PRECODE = makePDEBC(S, x, Q_r_faces, W_r_faces, DIM=DIM, system=True)
   ##insertCode(filename, extendDictionary( { "SNIP_PDEBC_SYSTEM_REDUCED_PRE" : PRECODE }, "SNIP_PDEBC_SYSTEM_REDUCED", CODE))
   ##CODE,PRECODE = makePDEBC(S, x, Q_r_faces, W_r_faces, DIM=DIM, system=False)
   ##insertCode(filename, extendDictionary( { "SNIP_PDEBC_SINGLE_REDUCED_PRE" : PRECODE }, "SNIP_PDEBC_SINGLE_REDUCED", CODE))

   # interpolation to Quad points
   #f_interpolation = 0
   #for i in range(len(GridOffsetText)):
   #  f_interpolation=f_interpolation + S[i] * Symbol("f_%s"%GridOffsetText[i])
   #CODE=createInterpolationCode(f_interpolation, x, Q, loopindex=idx_full)
   #insertCode(filename, { "SNIP_INTERPOLATE_ELEMENTS" : CODE})
   #CODE=createInterpolationCode(f_interpolation, x, Q_r, loopindex=idx_full)
   #insertCode(filename, { "SNIP_INTERPOLATE_REDUCED_ELEMENTS" : CODE})
   #CODE=""
   #for i in range(len(Q_faces)):
   #     CODE+="if (m_faceOffset[%s] > -1) {\n"%i
   #     CODE+=createInterpolationCode(f_interpolation, x, Q_faces[i], loopindex=idx_faces[i],  gridoffset="m_faceOffset[%s]"%i)
   #     CODE+="\n} /* end of face %s */\n"%i
   #insertCode(filename, { "SNIP_INTERPOLATE_FACES" : CODE})
   #CODE=""
   #for i in range(len(Q_r_faces)):
   #     CODE+="if (m_faceOffset[%s] > -1) {\n"%i
   #     CODE+=createInterpolationCode(f_interpolation, x, Q_r_faces[i], loopindex=idx_faces[i],  gridoffset="m_faceOffset[%s]"%i)
   #     CODE+="\n} /* end of face %s */\n"%i
   #insertCode(filename, { "SNIP_INTERPOLATE_REDUCED_FACES" : CODE})

   ## gradient to Quad points
   #CODE=createGradientCode(f_interpolation, x, Q, loopindex=idx_full)
   #insertCode(filename, { "SNIP_GRAD_ELEMENTS" : CODE})
   #CODE=createGradientCode(f_interpolation, x, Q_r, loopindex=idx_full)
   #insertCode(filename, { "SNIP_GRAD_REDUCED_ELEMENTS" : CODE})
   #CODE=""
   #for i in range(len(Q_faces)):
   #     CODE+="if (m_faceOffset[%s] > -1) {\n"%i
   #     CODE+=createGradientCode(f_interpolation, x, Q_faces[i], loopindex=idx_faces[i],  gridoffset="m_faceOffset[%s]"%i)
   #     CODE+="\n} /* end of face %s */\n"%i
   #insertCode(filename, { "SNIP_GRAD_FACES" : CODE})
   #CODE=""
   #for i in range(len(Q_r_faces)):
   #     CODE+="if (m_faceOffset[%s] > -1) {\n"%i
   #     CODE+=createGradientCode(f_interpolation, x, Q_r_faces[i], loopindex=idx_faces[i],  gridoffset="m_faceOffset[%s]"%i)
   #     CODE+="\n} /* end of face %s */\n"%i
   #insertCode(filename, { "SNIP_GRAD_REDUCED_FACES" : CODE})



def generatePDECode(DATA_A, EM, CONST_COEFFS, system=False):
    LOCAL_TMP={}
    for p in EM.keys():
        # first we collect the constant coefficients in the expression
        # (outside the parallel region):
        tmp={}
        for a in DATA_A:
            c=collect(EM[p], a, evaluate=False)
            if c.has_key(a):
                const_expr = c[a].simplify().subs(sqrt(3),Symbol('SQRT3'))
                if CONST_COEFFS.has_key(const_expr):
                    coeffsym = CONST_COEFFS[const_expr]
                    A_sym = a
                # cut down on constants - strategy 1
                elif CONST_COEFFS.has_key((-const_expr).simplify()):
                    coeffsym = CONST_COEFFS[(-const_expr).simplify()]
                    A_sym = -a
                # cut down on constants - strategy 2
                else:
                    A_sym = 0
                    for f in CONST_COEFFS.keys():
                        mult = f.extract_multiplicatively(const_expr)
                        if mult is not None and mult.is_Integer:
                            coeffsym = CONST_COEFFS[f]
                            A_sym = a/mult
                            break
                        mult = const_expr.extract_multiplicatively(f)
                        if mult is not None and mult.is_Integer:
                            coeffsym = CONST_COEFFS[f]
                            A_sym = a*mult
                            break
                    if A_sym == 0:
                        coeffsym = Symbol("w%s"%len(CONST_COEFFS))
                        CONST_COEFFS[const_expr]=coeffsym
                        A_sym = a
                if tmp.has_key(coeffsym):
                    tmp[coeffsym]+=A_sym
                else:
                    tmp[coeffsym]=A_sym
        # now we find temporary values which are combinations of the input
        # data:
        em_new=[]
        for coeffsym in tmp.keys():
            tt = (coeffsym * tmp[coeffsym]).simplify()
            if len(tt.free_symbols)<=2:
                em_new.append(tt)
            else:
                if not LOCAL_TMP.has_key(tt):
                    LOCAL_TMP[tt]=Symbol("tmp%s"%len(LOCAL_TMP))
                em_new.append(LOCAL_TMP[tt])
        EM[p]=em_new

    #####
    OUT=""
    tmp_key=lambda x: int(x[1].name[3:])
    for k, v in sorted(LOCAL_TMP.iteritems(), key=tmp_key):
         OUT+="const double %s = %s;\n"%(v,ccode(k))

    for p in sorted(EM.keys()):
        sum_em=0
        for t in EM[p]: sum_em=sum_em + t
        if not isinstance(sum_em, int):
            if isinstance(p, int):
                if system:
                    OUT+="EM_F[INDEX2(k,%s,numEq)]+=%s;\n"%(p, ccode(sum_em))
                else:
                    OUT+="EM_F[%s]+=%s;\n"%(p, ccode(sum_em))
            else:
                MULT=max([qq[0] for qq in EM.keys()])+1
                if system:
                    OUT+="EM_S[INDEX4(k,m,%s,%s,numEq,numComp,%s)]+=%s;\n"%(p[0],p[1], MULT, ccode(sum_em))
                else:
                    OUT+="EM_S[INDEX2(%s,%s,%s)]+=%s;\n"%(p[0], p[1], MULT, ccode(sum_em))

    return OUT

##############################################################################
##############################################################################
def makePDE(S, x, Q, W, DIM=2, system=False):
   print("Generating PDE code (%dD, %s)..."%(DIM, "System" if system else "Single"))
   CONST_COEFFS={}
   GLOBAL_N=0
   PRECODE="const double SQRT3 = %0.20f;\n"%sqrt(3)

   print("  Coefficient A...")
   CODE="""
///////////////
// process A //
///////////////
if (!A.isEmpty()) {
add_EM_S = true;
const double* A_p = A.getSampleDataRO(e);
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

        for q in range(len(Q)):
            for di in range(DIM):
               for dj in range(DIM):
                  A_name="A_%d%d_%d"%(di,dj,q)
                  A=Symbol(A_name)
                  DATA_A.append(A)
                  if system:
                      CODE2+="const double %s = A_p[INDEX5(k,%s,m,%s,%s,numEq,%s,numComp,%s)];\n"%(A_name, di, dj, q, DIM, DIM)
                  else:
                      CODE2+="const double %s = A_p[INDEX3(%s,%s,%s,%s,%s)];\n"%(A_name, di, dj, q, DIM, DIM)
                  for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[i],x[di]) * diff(S[j],x[dj])).subs( [(x[jj], Q[q][jj]) for jj in range(DIM)] )
        CODE+=CODE2+generatePDECode(DATA_A, EM, CONST_COEFFS, system)

        if system: CODE+="}\n}\n"
        CODE+="} else { // constant data\n"
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
                      CODE2+="const double %s = A_p[INDEX4(k,%s,m,%s, numEq,%s, numComp)];\n"%(A_name, di, dj, DIM)
             else:
                      CODE2+="const double %s = A_p[INDEX2(%s,%s,%s)];\n"%(A_name, di, dj, DIM)
             for q in range(len(Q)):
                  for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[i],x[di]) * diff(S[j],x[dj])).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_A, EM, CONST_COEFFS, system)
   if system:  CODE+="}\n}\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n"
   #BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
   print("  Coefficient B...")
   CODE+="""
///////////////
// process B //
///////////////
if (!B.isEmpty()) {
add_EM_S=true;
const double* B_p=B.getSampleDataRO(e);
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

        for q in range(len(Q)):
            for di in range(DIM):
                  A_name="B_%d_%d"%(di,q)
                  A=Symbol(A_name)
                  DATA_B.append(A)
                  if system:
                      CODE2+="const double %s = B_p[INDEX4(k,%s,m,%s, numEq,%s,numComp)];\n"%(A_name, di,  q, DIM)
                  else:
                      CODE2+="const double %s = B_p[INDEX2(%s,%s,%s)];\n"%(A_name, di, q, DIM)
                  for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[i],x[di]) * S[j]).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_B, EM, CONST_COEFFS, system)

        if system: CODE+="}\n}\n"
        CODE+="} else { // constant data\n"
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
                      CODE2+="const double %s = B_p[INDEX3(k,%s,m,numEq,%s)];\n"%(A_name, di,  DIM)
             else:
                      CODE2+="const double %s = B_p[%s];\n"%(A_name, di)
             for q in range(len(Q)):
                  for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[i],x[di]) * S[j] ).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_B, EM, CONST_COEFFS, system)
   if system:  CODE+="}\n}\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n"
   #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   print("  Coefficient C...")
   CODE+="""
///////////////
// process C //
///////////////
if (!C.isEmpty()) {
add_EM_S=true;
const double* C_p=C.getSampleDataRO(e);
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

        for q in range(len(Q)):
            for dj in range(DIM):
                  A_name="C_%d_%d"%(dj,q)
                  A=Symbol(A_name)
                  DATA_C.append(A)
                  if system:
                      CODE2+="const double %s = C_p[INDEX4(k,m,%s, %s, numEq,numComp,%s)];\n"%(A_name, dj,  q, DIM)
                  else:
                      CODE2+="const double %s = C_p[INDEX2(%s,%s,%s)];\n"%(A_name, dj, q, DIM)
                  for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[j],x[dj]) * S[i]).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_C, EM, CONST_COEFFS, system)

        if system: CODE+="}\n}\n"
        CODE+="} else { // constant data\n"
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
                      CODE2+="const double %s = C_p[INDEX3(k,m,%s,numEq,numComp)];\n"%(A_name, dj)
             else:
                      CODE2+="const double %s = C_p[%s];\n"%(A_name, dj)
             for q in range(len(Q)):
                  for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * diff(S[j],x[dj]) * S[i] ).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_C, EM, CONST_COEFFS, system)
   if system: CODE+="}\n}\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n"
   #DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
   print("  Coefficient D...")
   CODE+="""
///////////////
// process D //
///////////////
if (!D.isEmpty()) {
add_EM_S=true;
const double* D_p=D.getSampleDataRO(e);
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

        for q in range(len(Q)):
                  A_name="D_%d"%(q,)
                  A=Symbol(A_name)
                  DATA_D.append(A)
                  if system:
                      CODE2+="const double %s = D_p[INDEX3(k,m,%s,numEq,numComp)];\n"%(A_name, q)
                  else:
                      CODE2+="const double %s = D_p[%s];\n"%(A_name, q)
                  for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * S[j] * S[i]).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_D, EM, CONST_COEFFS, system)

        if system: CODE+="}\n }\n"
        CODE+="} else { // constant data\n"
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
                      CODE2+="const double %s = D_p[INDEX2(k, m, numEq)];\n"%(A_name,)
             else:
                      CODE2+="const double %s = D_p[0];\n"%(A_name,)
             for q in range(len(Q)):
                  for i in range(len(S)):
                      for j in range(len(S)):
                          EM[(i,j)] = EM[(i,j)] + (A * W[q] * S[j] * S[i] ).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_D, EM, CONST_COEFFS, system)
   if system: CODE+="}\n}\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n"

   #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   print("  Coefficient X...")
   CODE+="""
///////////////
// process X //
///////////////
if (!X.isEmpty()) {
add_EM_F=true;
const double* X_p=X.getSampleDataRO(e);
"""
   if len(Q) > 1:
        CODE+="if (X.actsExpanded()) {\n"
        if system: CODE+= "for (index_t k=0; k<numEq; k++) {\n"
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
                      CODE2+="const double %s = X_p[INDEX3(k,%s,%s,numEq,%s)];\n"%(A_name, dj,  q, DIM)
                  else:
                      CODE2+="const double %s = X_p[INDEX2(%s,%s,%s)];\n"%(A_name, dj,q,DIM)
                  for j in range(len(S)):
                          EM[j] = EM[j] + (A * W[q] * diff(S[j],x[dj])).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_X, EM, CONST_COEFFS, system)

        if system: CODE+="}\n"
        CODE+="} else { // constant data\n"
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
                      CODE2+="const double %s = X_p[INDEX2(k, %s, numEq)];\n"%(A_name, dj)
             else:
                      CODE2+="const double %s = X_p[%s];\n"%(A_name, dj)
             for q in range(len(Q)):
                  for j in range(len(S)):
                          EM[j] = EM[j] + (A * W[q] * diff(S[j],x[dj])).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_X, EM, CONST_COEFFS, system)
   if system:  CODE+="}\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n"

   #YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
   print("  Coefficient Y...")
   CODE+="""
///////////////
// process Y //
///////////////
if (!Y.isEmpty()) {
add_EM_F=true;
const double* Y_p=Y.getSampleDataRO(e);
"""
   if len(Q) > 1:
        CODE+="if (Y.actsExpanded()) {\n"
        if system: CODE+= "for (index_t k=0; k<numEq; k++) {\n"
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
                      CODE2+="const double %s = Y_p[INDEX2(k, %s, numEq)];\n"%(A_name, q)
                  else:
                      CODE2+="const double %s = Y_p[%s];\n"%(A_name, q)
                  for i in range(len(S)):
                          EM[i] = EM[i] + (A * W[q] * S[i]).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_Y, EM, CONST_COEFFS, system)

        if system: CODE+="}\n"
        CODE+="} else { // constant data\n"
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
                      CODE2+="const double %s = Y_p[k];\n"%(A_name,)
             else:
                      CODE2+="const double %s = Y_p[0];\n"%(A_name,)
             for q in range(len(Q)):
                  for i in range(len(S)):
                          EM[i] = EM[i] + (A * W[q] * S[i]).subs( [ (x[jj], Q[q][jj]) for jj in range(DIM) ] )
   CODE+=CODE2+generatePDECode(DATA_Y, EM, CONST_COEFFS, system)
   if system:  CODE+="}\n"
   if len(Q) > 1: CODE+="}\n"
   CODE+="}\n"

   w_key=lambda x: int(x[1].name[1:])
   # sort the constants in a special order...
   def num_oper(x):
        n=len(x[0].atoms())
        for s in x[0].free_symbols:
            if 'm_dx' in s.name:
                n+=(int(s.name.translate(None, 'm_dx[]'))+1)*10
        return n

   for k,v in sorted(CONST_COEFFS.iteritems(), key=num_oper):
       #PRECODE+="const double %s = %s;\n"%(v,ccode(k.evalf(n=DIGITS, chop=True)))
       PRECODE+="const double %s = %s;\n"%(v,ccode(k))
   return CODE, PRECODE


##############################################################################
##############################################################################
def makePDEBC(S, x, Q, W, DIM, system=True):
   print("Generating PDE BC code (%dD, %s)..."%(DIM, "System" if system else "Single"))
   CONST_COEFFS={}
   GLOBAL_N=0
   PRECODE="const double SQRT3 = %0.20f;\n"%sqrt(3)
   CODE_OUT=[]
    
   for n in range(len(Q)):
        CODE="""
///////////////
// process d //
///////////////
if (add_EM_S) {
const double* d_p=d.getSampleDataRO(e);
"""
        if len(Q[n]) > 1:
            CODE+="if (d.actsExpanded()) {\n"
            if system: CODE+= """for (index_t k=0; k<numEq; k++) {
for (index_t m=0; m<numComp; m++) {
"""
            DATA_D=[]
            CODE2=""
            EM = {}
            for i in range(len(S)):
                for j in range(len(S)):
                    EM[(i,j)] = 0

            for q in range(len(Q[n])):
                A_name="d_%d"%(q,)
                A=Symbol(A_name)
                DATA_D.append(A)
                if system:
                    CODE2+="const double %s = d_p[INDEX3(k,m,%s,numEq,numComp)];\n"%(A_name, q)
                else:
                    CODE2+="const double %s = d_p[%s];\n"%(A_name, q)
                for i in range(len(S)):
                    for j in range(len(S)):
                        EM[(i,j)] = EM[(i,j)] + (A * W[n][q] * S[j] * S[i]).subs( [ (x[jj], Q[n][q][jj]) for jj in range(DIM) ] )
            CODE+=CODE2+generatePDECode(DATA_D, EM, CONST_COEFFS, system)

            if system: CODE+="}\n }\n"
            CODE+="} else { // constant data\n"
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

        A_name="d_0"
        A=Symbol(A_name)
        DATA_D.append(A)
        if system:
            CODE2+="const double %s = d_p[INDEX2(k, m, numEq)];\n"%(A_name,)
        else:
            CODE2+="const double %s = d_p[0];\n"%(A_name,)
        for q in range(len(Q[n])):
            for i in range(len(S)):
                for j in range(len(S)):
                    EM[(i,j)] = EM[(i,j)] + (A * W[n][q] * S[j] * S[i] ).subs( [ (x[jj], Q[n][q][jj]) for jj in range(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_D, EM, CONST_COEFFS, system)
        if system: CODE+="}\n}\n"
        if len(Q[n]) > 1: CODE+="}\n"
        CODE+="}\n"

        CODE+="""
///////////////
// process y //
///////////////
if (add_EM_F) {
const double* y_p=y.getSampleDataRO(e);
"""
        if len(Q[n]) > 1:
              CODE+="if (y.actsExpanded()) {\n"
              if system: CODE+= "for (index_t k=0; k<numEq; k++) {\n"
              DATA_Y=[]
              CODE2=""
              EM = {}
              for i in range(len(S)):
                  EM[i] = 0

              for q in range(len(Q[n])):
                  A_name="y_%d"%(q,)
                  A=Symbol(A_name)
                  DATA_Y.append(A)
                  if system:
                      CODE2+="const double %s = y_p[INDEX2(k, %s, numEq)];\n"%(A_name, q)
                  else:
                      CODE2+="const double %s = y_p[%s];\n"%(A_name, q)
                  for i in range(len(S)):
                      EM[i] = EM[i] + (A * W[n][q] * S[i]).subs( [ (x[jj], Q[n][q][jj]) for jj in range(DIM) ] )
              CODE+=CODE2+generatePDECode(DATA_Y, EM, CONST_COEFFS, system)

              if system: CODE+="}\n"
              CODE+="} else { // constant data\n"
        if system:
            CODE+= "for (index_t k=0; k<numEq; k++) {\n"
        DATA_Y=[]
        CODE2=""
        EM = {}
        for i in range(len(S)):
                EM[i] = 0

        if system:
            A_name="y_p[k]"
            #CODE2+="const double %s = y_p[k];\n"%(A_name,)
        else:
            A_name="y_p[0]"
            #CODE2+="const double %s = y_p[0];\n"%(A_name,)
        A=Symbol(A_name)
        DATA_Y.append(A)
        for q in range(len(Q[n])):
            for i in range(len(S)):
                EM[i] = EM[i] + (A * W[n][q] * S[i]).subs( [ (x[jj], Q[n][q][jj]) for jj in range(DIM) ] )
        CODE+=CODE2+generatePDECode(DATA_Y, EM, CONST_COEFFS, system)
        if system:  CODE+="}\n"
        if len(Q[n]) > 1: CODE+="}\n"
        CODE+="}\n"
          
        CODE_OUT.append(CODE)

   w_key=lambda x: int(x[1].name[1:])
   # sort the constants in a special order...
   def num_oper(x):
        n=len(x[0].atoms())
        for s in x[0].free_symbols:
            if 'm_dx' in s.name:
                n+=(int(s.name.translate(None, 'm_dx[]'))+1)*10
        return n

   for k,v in sorted(CONST_COEFFS.iteritems(), key=num_oper):
       #PRECODE+="const double %s = %s;\n"%(v,ccode(k.evalf(n=DIGITS, chop=True)))
       PRECODE+="const double %s = %s;\n"%(v,ccode(k))
   return CODE_OUT, PRECODE
    
##############################################################################
##############################################################################
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
                 k0=Symbol("c%s"%cc)
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

   for i in range(len(Q[0])):
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

   #for i in range(DIM-1,-1,-1):
   #       if loopindex[i] <= -1: TXT+="const index_t k%i = 0;\n"%i
   for a,v in  consts.items():
      TXT+="const double %s = %s;\n"%(ccode(a), ccode(v.evalf(n=DIGITS)))
   TXT+="#pragma omp parallel for"

   for i in range(DIM-1,-1,-1):
          if loopindex[i] > -1: TXT+="\nfor (index_t k%i=0; k%i < m_NE%i; ++k%i) {"%(loopindex[i],loopindex[i],loopindex[i],loopindex[i])
   for s in args:
            k1=""
            for i in range(DIM):
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
            TXT+="\nconst double* %s = in.getSampleDataRO(INDEX%s(%s, %s));"%(s.name,DIM,k1[1:],M1)
   if len(gridoffset) > 0: IDX2=gridoffset+"+"+IDX2
   TXT+="\ndouble* o = out.getSampleDataRW(%s);"%IDX2
   TXT+="\nfor (index_t i=0; i < numComp; ++i) {\n"
   for q in range(len(Q)):
       for i in range(DIM):
           dFidx=dF[i][q]
           for s in args:
               dFidx=dFidx.subs(s.name, Symbol(s.name+"[i]"))
           TXT+="o[INDEX3(i,%s,%s,numComp,%s)] = %s;\n"%(i,q,DIM,ccode(dFidx))
           #TXT+="out[INDEX4(i,%s,%s,%s,numComp,%s,%s)] = %s;\n"%(i,q,IDX2,DIM,len(Q),ccode(dF[i][q]))
   TXT+="} /* end of component loop i */\n"
   for i in range(DIM):
           if loopindex[i]>-1 : TXT+="} /* end of k%i loop */\n"%loopindex[i]
   return TXT

def createInterpolationCode(F, x, Q, gridoffset="", loopindex=(0,1,2)):
   DIM=len(loopindex)
   F, args, consts = optimizeEvaluate([F, ], x, Q)
   F=F[0]
   k=[]
   N=[]
   M1=""
   TXT=""
   for i in range(len(Q[0])):
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

   #for i in range(DIM-1,-1,-1):
   #       if loopindex[i] <= -1: TXT+="const index_t k%i = 0;\n"%i
   for a,v in sorted(consts.items()):
      TXT+="const double %s = %s;\n"%(ccode(a), ccode(v.evalf(n=DIGITS)))
   TXT+="#pragma omp parallel\n{"
   for s in sorted(args):
       TXT+="\nvector<double> %s(numComp);"%s.name
   TXT+="\n#pragma omp for"
   for i in range(DIM-1,-1,-1):
          if loopindex[i] > -1: TXT+="\nfor (index_t k%i=0; k%i < m_NE%i; ++k%i) {"%(loopindex[i],loopindex[i],loopindex[i],loopindex[i])
   for s in sorted(args):
        k1=""
        for i in range(DIM):
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
        TXT+="\nmemcpy(&%s[0], in.getSampleDataRO(INDEX%s(%s, %s)), numComp*sizeof(double));"%(s.name,DIM,k1[1:],M1)
   if len(gridoffset) > 0: IDX2=gridoffset+"+"+IDX2
#        for i in range(DIM):
#             if loopindex[i]>-1: IDX2=gridoffset+"+k%s"%i
   TXT+="\ndouble* o = out.getSampleDataRW(%s);"%IDX2
   TXT+="\nfor (index_t i=0; i < numComp; ++i) {\n"
   #  interpolation to quadrature points
   for q in range(len(Q)):
       Fidx=F[q]
       for s in args:
           Fidx=Fidx.subs(s.name, Symbol(s.name+"[i]"))
       TXT+="o[INDEX2(i,numComp,%s)] = %s;\n"%(q,ccode(Fidx))
   TXT+="} /* end of component loop i */\n"
   for i in range(DIM):
           if loopindex[i]>-1 : TXT+="} /* end of k%i loop */\n"%loopindex[i]
   TXT+="} /* end of parallel section */\n"
   return TXT


def extendDictionary(this, tag, items):
    for i in range(len(items)):
         this["%s_%s"%(tag,i)]=items[i]
    return this


def insertCode(fn, replacement):
    """
    Creates a backup copy of ``fn``, then applies string replacements from the
    ``replacement`` dictionary with proper indentation (tabsize=4).
    """
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
                    elif rep_l.find("//")>-1:
                        s=(rep_l[:rep_l.find("//")]).strip()
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
        old_code=new_code[:-1]
    open(fn,'w').write(old_code)
    print("file %s backed up and updated."%fn)


if __name__ == "__main__":
    filenames={2:"Rectangle.cpp", 3:"Brick.cpp"}
    filenames={3:"Brick.cpp"}
    for d in filenames.keys():
         Process(target=generateAll, args=(d, filenames[d])).start()

