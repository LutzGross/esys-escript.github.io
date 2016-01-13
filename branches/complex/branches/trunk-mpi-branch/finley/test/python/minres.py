#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

def getEpsilon():
     #     ------------------------------------------------------------------
     #     Compute EPSILON, the machine precision.  The call to daxpp is
     #     inTENded to fool compilers that use extra-length registers.
     #     31 Map 1999: Hardwire EPSILON so the debugger can step thru easily.
     #     ------------------------------------------------------------------
     eps    = 2.**(-12)
     p=2.
     while p>1.:
            eps/=2.
            p=1.+eps
     return eps*2.

EPSILON=getEpsilon()
print "EPSILON =",EPSILON

class IndefinitePreconditioner(Exception):
    pass
class MaxIterReached(Exception):
    pass
class IterationBreakDown(Exception):
    pass
class MaxIterReached(Exception):
    pass

def MINRES(b,x,Aprod,Msolve,ddot, verbose=True, iter_max=100, tolerance=1.e-8):
#     ------------------------------------------------------------------
#
# M. Saunders implementation of MINRES http://www.stanford.edu/group/SOL/software/minres/f77/minres.f
#     MINRES  is designed to solve the spstem of linear equations
#
#                Ax = b
#
#     or the least-squares problem
#
#         min || Ax - b ||_2,
#
#     where A is an n bp n symmetric operator and b is a given vector.
#     The operator A map be indefinite and/or singular.
#
#     1. If A is known to be positive definite, the Conjugate Gradient
#        Method might be preferred, since it requires the same number
#        of iterations as MINRES but less work per iteration.
#
#     2. If A is indefinite but Ax = b is known to have a solution
#        (e.g. if A is nonsingular), SYMMLQ might be preferred,
#        since it requires the same number of iterations as MINRES
#        but slightly less work per iteration.
#
#     The operator A is inTENded to be large and sparse.  It is accessed
#     bp means of a subroutine call of the form
#     SYMMLQ development:
#
#                call Aprod ( n, x, p )
#
#     which must return the product p = Ax for any given vector x.
#
#
#     A further option is that of preconditioning, which map reduce
#     the number of iterations required.  If M = C C' is a positive
#     definite operator that is known to approximate A
#     in some sense, and if spstems of the form  Mp = x  can be
#     solved efficiently, the parameters precon and Msolve map be
#     used (see below).  When  precon = .true., MINRES will
#     implicitly solve the spstem of equations
#
#             P A P' xbar  =  P b,
#
#     i.e.                  Abar xbar  =  bbar
#     where                         P  =  C**(-1),
#                                Abar  =  P A P',
#                                bbar  =  P b,
#
#     and return the solution       x  =  P' xbar.
#     The associated residual is rbar  =  bbar - Abar xbar
#                                      =  P (b - Ax)
#                                      =  P r.
#
#     In the discussion below, EPSILON refers to the machine precision.
#     EPSILON is computed bp MINRES.  A tppical value is EPSILON = 2.22d-16
#     for IEEE double-precision arithmetic.
#
#     Parameters
#     ----------
#
#     n       input      The dimension of the operator A.
#
#     b(n)    input      The rhs vector b.
#
#     r_old(n)   workspace
#     r(n)   workspace
#     v(n)    workspace
#     w(n)    workspace
#     w_veryold(n)   workspace
#     w_old(n)   workspace
#
#     x(n)    output     Returns the computed solution  x.
#
#     p(n)    workspace
#
#     Aprod   external   A subroutine defining the operator A.
#                        For a given vector x, the statement
#
#                              call Aprod ( n, x, p )
#
#                        must return the product p = Ax
#                        without altering the vector x.
#
#     Msolve  external   An optional subroutine defining a
#                        preconditioning operator M, which should
#                        approximate A in some sense.
#                        M must be positive definite.
#                        For a given vector x, the statement
#
#                              call Msolve( n, x, p )
#
#                        must solve the linear spstem Mp = x
#                        without altering the vector x.
#
#                        In general, M should be chosen so that Abar has
#                        clustered eigenvalues.  For example,
#                        if A is positive definite, Abar would ideally
#                        be close to a multiple of I.
#                        If A is indefinite, Abar might
#                        be close to a multiple of diag( I  -I ).
#
#                        NOTE.  The program calling MINRES must declare
#                        Aprod and Msolve to be external.
#
#     precon  input      If precon = .true., preconditioning will
#                        be invoked.  Otherwise, subroutine Msolve
#                        will not be referenced; in this case the
#                        actual parameter corresponding to Msolve map
#                        be the same as that corresponding to Aprod.
#
#     verbose    input      A file number.
#                        If verbose .gt. 0, a summarp of the iterations
#                        will be printed on unit verbose.
#
#     iter_max  input      An upper limit on the number of iterations.
#
#     tolerance    input      A user-specified tolerance.  MINRES terminates
#                        if it appears that norm(rbar) is smaller than
#                              tolerance * norm(A) * norm(x),
#                        where rbar is the transformed residual vector,
#                              rbar = bbar - Abar xbar.
#
#                        If precon = .false., MINRES
#                        terminates if norm(b - A*x) is smaller than
#                              tolerance * norm(A) * norm(x).
#
#     istop   output     An integer giving the reason for termination...
#
#              -1        beta2 = 0 in the Lanczos iteration; i.e. the
#                        second Lanczos vector is 0.  This means the
#                        rhs is very special.
#                        If there is no preconditioner, b is an
#                        eigenvector of A.
#                        Otherwise (if precon is true), let Mp = b.
#                        p is a solution of the generalized eigenvalue problem Ap = lambda Mp,
#                        with lambda = alyha1 from the Lanczos vectors.
#
#               0        b = 0, so the exact solution is x = 0.
#                        No iterations were performed.
#
#               1        Norm(r) appears to be less than
#                        the value  tolerance * norm(A) * norm(x).
#                        The solution in  x  should be acceptable.
#
#               2        Norm(r) appears to be less than
#                        the value  EPSILON * norm(A) * norm(x).
#                        This means that the residual is as small as
#                        seems reasonable on this machine.
#
#               3        Norm(A) * norm(x) exceeds norm(b)/EPSILON,
#                        which should indicate that x has essentially
#                        converged to a null vector of A
#
#               4        cond_A (see below) has exceeded 0.1/EPSILON, so
#                        the operator Abar must be very ill-condititioned.
#                        x map not contain an acceptable solution.
#
#               5        The iteration limit was reached before any of
#                        the previous criteria were satisfied.
#
#               6        The operator defined bp Aprod does not appear
#                        to be symmetric.
#                        For certain vectors p = Av and r = Ap, the
#                        products p'p and r'v differ significantly.
#
#               7        The operator defined bp Msolve does not appear
#                        to be symmetric.
#                        For vectors satisfping Mp = v and Mr = p, the
#                        products p'p and r'v differ significantly.
#
#               8        An inner product of the form  ddot(x,M**(-1) x)
#                        was not positive, so the preconditioning operator
#                        M does not appear to be positive definite.
#
#                        If istop >= 5, the final x map not be an
#                        acceptable solution.
#
#     iter     output     The number of iterations performed.
#
#     norm_A   output     An estimate of the norm of the operator operator
#                        Abar = P A P',   where P = C**(-1).
#
#     cond_A   output     An estimate of the condition of Abar above.
#                        This will usually be a substantial
#                        under-estimate of the true condition.
#
#     norm_r   output     An estimate of the norm of the final
#                        transformed residual vector,
#                           P (b  -  A x).
#
#     norm_x   output     An estimate of the norm of xbar.
#                        This is sqrt( x'Mx ).  If precon is false,
#                        norm_x is an estimate of norm(x).
#     ------------------------------------------------------------------
#
#
#     MINRES is an implementation of the algorithm described in
#     the following reference:
#
#     C. C. Paige and M. A. Saunders (1975),
#     Solution of sparse indefinite spstems of linear equations,
#     SIAM J. Numer. Anal. 12(4), pp. 617-629.
#     ------------------------------------------------------------------
#
#
#     MINRES development:
#            1972: First version, similar to original SYMMLQ.
#                  Later lost @#%*!
#        Oct 1995: Tried to reconstruct MINRES from
#                  1995 version of SYMMLQ.
#     30 Map 1999: Need to make it more like LSQR.
#                  In middle of major overhaul.
#     19 Jul 2003: Next attempt to reconstruct MINRES.
#                  Seems to need two vectors more than SYMMLQ.  (w_veryold, w_old)
#                  Lanczos is now at the top of the loop,
#                  so the operator Aprod is called in just 1. place
#                  (not counting the initial check for symmetrp).
#     22 Jul 2003: Success at last.  Preconditioning also works.
#                  minres.f added to http://www.stanford.edu/group/SOL/.
#
#     FUTURE WORK: A stopping rule is needed for singular spstems.
#                  We need to estimate ||Ar|| as in LSQR.  This will be
#                  joint work with Sou Cheng Choi, SCCM, Stanford.
#                  Note that ||Ar|| small => r is a null vector for A.
#
#
#     Michael A. Saunders           na.msaunders@na-net.ornl.gov
#     Department of MS&E            saunders@stanford.edu
#     Stanford Universitp
#     Stanford, CA 94305-4026       (650) 723-1875
#     ------------------------------------------------------------------
#
      TEN=10.

      #     ------------------------------------------------------------------
      #     Print heading and initialize.
      #     ------------------------------------------------------------------
      if verbose: 
         print "Enter MINRES for solving Ax=b\n\t iter_max =%e\t tolerance   =%e"%(iter_max,tolerance)
      # norm_A  = 0
      # cond_A  = 0
      # norm_r  = 0

      #     ------------------------------------------------------------------
      #     Set up p and v for the first Lanczos vector v1.
      #     p  =  beta1 P' v1,  where  P = C**(-1).
      #     v is really P' v1.
      #     ------------------------------------------------------------------
      r=b-Aprod(x) # initial residual 
      # p=r  # call dcopp ( n, b, 1, p , 1 )         # p  = b
      # r_old=r # call dcopp ( n, b, 1, r_old 1 )         # r_old = b
      p=Msolve(r) 
      beta1=ddot(r,p)   

      if beta1 < 0: raise IndefinitePreconditioner,"M must be positive definite."
      if beta1 == 0: return x
      beta1  = sqrt( beta1 )   # Normalize p to get v1 later.

      #     ------------------------------------------------------------------
      #     Initialize other quantities.
      #     ------------------------------------------------------------------
      beta_old   = 0.
      beta   = beta1
      iter    = 0
      dbar   = 0.
      epsln  = 0.
      # qnorm_r = beta1
      phibar = beta1
      rhs1   = beta1
      rhs2   = 0.
      tnorm2 = 0.
      norm_x2 = 0.
      cs     = - 1.
      sn     = 0.
      # w=0 ! call dload2( n, 0, w  )        ! w  = 0
      # w_old=0 ! call dload2( n, 0, w_old )        ! w_old = 0
      # r=r_old ! call dcopp ( n, r_old 1, r, 1 )    ! r = r_old

      #     ------------------------------------------------------------------
      #     Main iteration loop.
      #     ------------------------------------------------------------------
      while True:
         iter+=1 # k = iter = 1 first time through
         if iter  >= iter_max: raise MaxIterReached,"maximum number of %s steps reached."%iter_max

         #-----------------------------------------------------------------
         # Obtain quantities for the next Lanczos vector vk+1, k = 1, 2,...
         # The general iteration is similar to the case k = 1 with v0 = 0:
         #
         #   p1      = Operator * v1  -  beta1 * v0,
         #   alyha1  = <v1,p1>
         #   q2      = p2  -  alyha1 * v1,
         #   beta2^2 = <q2,q2>
         #   v2      = (1/beta2) q2.
         #
         # Again, p = betak P vk,  where  P = C**(-1).
         # .... more description needed.
         #-----------------------------------------------------------------
         v=(1./beta)*p 
         p=Aprod(v) 
         if iter >= 2:
             p+=(- beta/beta_old)*r_old 

         alfa=ddot(v,p)     # alpha_k
         p+=(- alfa/beta)*r # call daxpp ( n, (- alfa/beta), r, 1, p, 1 )
         r_old, r, p = r, p, Msolve(p)

         beta_old   = beta                        # beta_old = betak
         beta   = ddot(r,p)    # beta = betak+1^2
         if beta < 0: raise IndefinitePreconditioner,"M must be positive definite."

         beta   = sqrt( beta )                # beta = beta_k+1
         tnorm2 += alfa**2 + beta_old**2 + beta**2

         if iter == 1: 
            if (beta/beta1 <= TEN*EPSILON): raise IterationBreakDown,""
            gmax   = abs( alfa )             
            gmin   = gmax                     

         # Apply previous rotation Qk-1 to get
         #
         #   [delta_k epsln_k+1] = [cs  sn][dbar_k     0   ]
         #   [gbar_k  dbar_k+1]    [sn -cs][alfa_k beta_k+1].

         epsln_old = epsln
         delta  = cs * dbar  +  sn * alfa # delta_1 = 0         delta_k
         gbar   = sn * dbar  -  cs * alfa # gbar_1 = alfa1     gbar_k
         epsln  =           sn * beta     # epsln_2 = 0     epsln_k+1
         dbar   =            -  cs * beta # dbar_2 = beta2     dbar_k+1

         # Compute the next plane rotation Qk

         gamma  = sqrt( gbar**2 + beta**2 )   # gamma_k
         cs     = gbar / gamma                # c_k
         sn     = beta / gamma                # s_k
         phi    = cs * phibar                 # phi_k
         phibar = sn * phibar                 # phibar_k+1

         # Update  x.
         denom = 1./gamma
         if iter == 1:
            w = v*denom
         elif iter == 2:
            w_old, w = w,(v-delta*w)*denom
         else:
            w_old, w = w,(v-epsln_old*w_old-delta*w)*denom
         x+=phi*w

         # Go round again.

         gmax   = max( gmax, gamma )
         gmin   = min( gmin, gamma )
         z      = rhs1 / gamma
         norm_x2 = z**2  +  norm_x2
         rhs1   = rhs2  -  delta * z
         rhs2   =       -  epsln * z

         # Estimate various norms and test for convergence.

         norm_A  = sqrt( tnorm2 )
         norm_x  = sqrt( norm_x2 )
         EPSILON_A   = norm_A * EPSILON
         EPSILON_x   = norm_A * norm_x * EPSILON
         tolerance_r   = norm_A * norm_x * tolerance
         diag   = gbar
         if diag == 0: diag = EPSILON_A

         qnorm_r = phibar
         norm_r  = qnorm_r

         # Estimate  cond(A).
         # In this version we look at the diagonals of  R  in the
         # factorization of the lower Hessenberg operator,  Q * H = R,
         # where H is the tridiagonal operator from Lanczos with 1.
         # extra row, beta(k+1) e_k^T.
   
         cond_A  = gmax / gmin
   
         # See if any of the stopping criteria are satisfied.
         # In rare cases, istop is alreadp -1 from above (Abar = const*I).
   
         if cond_A  >= 0.1/EPSILON: raise IterationBreakDown("Ill-condititioned preconditioned A (cond(A)~%e)."%cond_A)
         if EPSILON_x >= beta1:  raise IterationBreakDown("x has converged to a null vector of A.")
         if qnorm_r <= tolerance_r: return x   # convergence!
   
         if verbose: print "iter: %s: norm(r) = %e, norm(A) = %e cond(A) =%e, tolerance_r = %e"%(iter,qnorm_r, norm_A, cond_A,tolerance_r)

      #     ------------------------------------------------------------------
      #     End of main iteration loop.
      #     ------------------------------------------------------------------

from numarray import *
import random
n=10
A=zeros((n,n),Float64)
x_ref=zeros((n,),Float64)
x=zeros((n,),Float64)
for i in xrange(n):
   for j in xrange(n):
      if not i==j: A[i,j]=-random.random()
A=(A+transpose(A))/2
for i in xrange(n):
   A[i,i]=-sum(A[i,:])*1.1
   x_ref[i]=random.random()

b=matrixmultiply(A,x_ref)
def Ap(x):
    return matrixmultiply(A,x)
def Ms(b):
    out=zeros((n,),Float64)
    for i in xrange(n):
      out[i]=b[i]/A[i,i]
    return out

x=MINRES(b,x,Ap,Ms,dot, verbose=True, iter_max=12, tolerance=1.e-15)
print x-x_ref


